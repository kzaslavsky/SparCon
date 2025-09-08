import argparse
import json

from .client import ClinVarClient


def main():
    parser = argparse.ArgumentParser(description="Fetch ClinVar phenotypes and PMIDs for HGVS variants")
    parser.add_argument("hgvs", nargs="+", help="HGVS expression(s)")
    parser.add_argument("--email", required=True, help="Email for NCBI requests")
    parser.add_argument("--api-key", help="NCBI API key")
    args = parser.parse_args()

    client = ClinVarClient(email=args.email, api_key=args.api_key)
    records = [client.variant_record(h) for h in args.hgvs]
    print(json.dumps([r.__dict__ for r in records], default=lambda o: o.__dict__, indent=2))


if __name__ == "__main__":
    main()
