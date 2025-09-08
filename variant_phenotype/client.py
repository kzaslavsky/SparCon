import gzip
import io
import json
import time
import xml.etree.ElementTree as ET
from typing import Iterable, List, Optional

import requests

from .models import (
    Condition,
    PMIDCitation,
    SCVEvidence,
    VariantRecord,
)

NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
VARIATION_BASE = "https://api.ncbi.nlm.nih.gov/variation/v0"
CLINVAR_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited"


class ClinVarClient:
    """Client wrapping ClinVar/NCBI APIs for variant→phenotype pipeline."""

    def __init__(self, email: str, api_key: Optional[str] = None, tool: str = "variant_phenotype"):
        self.email = email
        self.api_key = api_key
        self.tool = tool
        self.session = requests.Session()

    # ---------------------------------------------
    # Utility helpers
    # ---------------------------------------------
    def _params(self, extra: Optional[dict] = None) -> dict:
        params = {"tool": self.tool, "email": self.email}
        if self.api_key:
            params["api_key"] = self.api_key
        if extra:
            params.update(extra)
        return params

    def _get(self, url: str, params: Optional[dict] = None, **kwargs) -> requests.Response:
        full_params = self._params(params)
        for attempt in range(3):
            resp = self.session.get(url, params=full_params, **kwargs)
            if resp.status_code == 429:
                time.sleep(1 + attempt)
                continue
            resp.raise_for_status()
            return resp
        resp.raise_for_status()

    # ---------------------------------------------
    # Step 1: HGVS → SPDI
    # ---------------------------------------------
    def normalize_hgvs(self, hgvs: str) -> VariantRecord:
        url = f"{VARIATION_BASE}/hgvs/{hgvs}/contextuals"
        resp = self._get(url)
        data = resp.json()
        record = VariantRecord(input_hgvs=hgvs)
        if "contextuals" in data:
            for contextual in data["contextuals"]:
                spdi = contextual.get("spdi")
                if spdi:
                    record.spdis.append("{seqid}:{position}:{deleted}:{inserted}".format(**spdi))
        normalized = data.get("hgvs")
        if normalized:
            record.normalized_hgvs = normalized
        return record

    # ---------------------------------------------
    # Step 2: resolve ClinVar Variation ID
    # ---------------------------------------------
    def _esearch(self, term: str) -> List[int]:
        params = {"db": "clinvar", "term": term, "retmax": 100}
        resp = self._get(f"{NCBI_EUTILS_BASE}/esearch.fcgi", params=params)
        ids = resp.json()["esearchresult"].get("idlist", [])
        return [int(i) for i in ids]

    def _spdi_to_term(self, spdi: str) -> Optional[str]:
        try:
            seq, pos, ref, alt = spdi.split(":")
            # seq is accession like NC_000013.11, convert to chr number if possible
            if seq.startswith("NC_") and len(seq) >= 15:
                chr_num = seq[7:9].lstrip("0")  # crude
                return f"{chr_num}:{int(pos)+1}:{ref}:{alt}"  # pos is 0-based in SPDI
        except Exception:
            return None
        return None

    def resolve_variation(self, record: VariantRecord) -> VariantRecord:
        ids = self._esearch(record.input_hgvs)
        if not ids and record.spdis:
            term = self._spdi_to_term(record.spdis[0])
            if term:
                ids = self._esearch(term)
        if ids:
            record.variation_id = ids[0]
            record.vcv_accession = f"VCV{ids[0]:09d}"
        return record

    # ---------------------------------------------
    # Step 3: traits via esummary
    # ---------------------------------------------
    def fetch_traits(self, record: VariantRecord) -> VariantRecord:
        if not record.variation_id:
            return record
        params = {"db": "clinvar", "id": record.variation_id, "retmode": "json"}
        resp = self._get(f"{NCBI_EUTILS_BASE}/esummary.fcgi", params=params)
        result = resp.json().get("result", {})
        v = result.get(str(record.variation_id), {})
        traits = []
        for trait in v.get("trait_set", []):
            cond = Condition(name=trait.get("trait_name"))
            ids = trait.get("xrefs", [])
            for i in ids:
                db = i.get("db")
                if db and i.get("id"):
                    cond.identifiers[db.lower()] = i["id"]
            traits.append(cond)
        record.traits = traits
        summary = v.get("variation_set" , {})
        if "allele_id" in summary:
            record.allele_ids = summary.get("allele_id", [])
        record.rsid = v.get("dbsnp_rs")
        return record

    # ---------------------------------------------
    # Step 4A: aggregated PMIDs from var_citations.txt
    # ---------------------------------------------
    def aggregated_pmids(self, record: VariantRecord) -> VariantRecord:
        if not record.variation_id:
            return record
        url = f"{CLINVAR_FTP_BASE}/var_citations.txt.gz"
        resp = self.session.get(url)
        resp.raise_for_status()
        with gzip.GzipFile(fileobj=io.BytesIO(resp.content)) as fh:
            for line in fh:
                parts = line.decode().strip().split("\t")
                if not parts or parts[0] == "#VariationID":
                    continue
                if int(parts[0]) == record.variation_id:
                    pmids = [int(p) for p in parts[1].split(",") if p]
                    record.aggregated_pmids.extend(pmids)
                    break
        return record

    # ---------------------------------------------
    # Step 4B: SCV-level PMIDs via VCV XML
    # ---------------------------------------------
    def scv_pmids(self, record: VariantRecord) -> VariantRecord:
        if not record.variation_id:
            return record
        params = {"db": "clinvar", "rettype": "vcv", "is_variationid": 1, "id": record.variation_id}
        resp = self._get(f"{NCBI_EUTILS_BASE}/efetch.fcgi", params=params)
        root = ET.fromstring(resp.content)
        evidence = []
        for assertion in root.findall('.//ClinVarAssertion'):
            scv = assertion.findtext('ID')
            submitter = assertion.findtext('SubmitterName')
            review = assertion.findtext('ReviewStatus')
            significance = assertion.findtext('ClinicalSignificance/Description')
            last_eval = assertion.findtext('ClinicalSignificance/LastEvaluated')
            cond_elem = assertion.find('.//TraitSet/Trait')
            cond = Condition()
            if cond_elem is not None:
                cond.name = cond_elem.findtext('Name/ElementValue')
                for xref in cond_elem.findall('XRef'):
                    db = xref.get('DB')
                    idx = xref.get('ID')
                    if db and idx:
                        cond.identifiers[db.lower()] = idx
            pmids = []
            for path in [
                './/ClinicalSignificance/Citation/ID',
                './/ObservedIn/ObservedData/Citation/ID',
                './/AttributeSet/Citation/ID',
            ]:
                for cite in assertion.findall(path):
                    if cite.get('Source') == 'PubMed' and cite.text and cite.text.isdigit():
                        pmids.append(PMIDCitation(pmid=int(cite.text), context=path))
            evidence.append(
                SCVEvidence(
                    scv_accession=scv,
                    submitter=submitter,
                    review_status=review,
                    significance=significance,
                    last_evaluated=last_eval,
                    condition=cond,
                    pmids=pmids,
                )
            )
        record.scv_evidence = evidence
        return record

    # ---------------------------------------------
    def variant_record(self, hgvs: str) -> VariantRecord:
        record = self.normalize_hgvs(hgvs)
        record = self.resolve_variation(record)
        record = self.fetch_traits(record)
        record = self.aggregated_pmids(record)
        record = self.scv_pmids(record)
        return record


__all__ = ["ClinVarClient"]
