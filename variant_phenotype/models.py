from dataclasses import dataclass, field
from typing import List, Optional, Dict


@dataclass
class Condition:
    name: Optional[str] = None
    identifiers: Dict[str, str] = field(default_factory=dict)


@dataclass
class PMIDCitation:
    pmid: int
    context: Optional[str] = None


@dataclass
class SCVEvidence:
    scv_accession: Optional[str] = None
    submitter: Optional[str] = None
    review_status: Optional[str] = None
    significance: Optional[str] = None
    last_evaluated: Optional[str] = None
    condition: Condition = field(default_factory=Condition)
    pmids: List[PMIDCitation] = field(default_factory=list)


@dataclass
class VariantRecord:
    input_hgvs: str
    normalized_hgvs: Optional[str] = None
    spdis: List[str] = field(default_factory=list)
    variation_id: Optional[int] = None
    vcv_accession: Optional[str] = None
    allele_ids: List[int] = field(default_factory=list)
    rsid: Optional[str] = None
    traits: List[Condition] = field(default_factory=list)
    aggregated_pmids: List[int] = field(default_factory=list)
    scv_evidence: List[SCVEvidence] = field(default_factory=list)
