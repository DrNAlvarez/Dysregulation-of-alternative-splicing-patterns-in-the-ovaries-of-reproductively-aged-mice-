#!/usr/bin/env python3
"""
fetch_ensembl_biotype.py

Fetch the full Ensembl biotypes for a given dataset and save as TSV.
"""

import argparse
from pybiomart import Server
import pandas as pd
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Fetch Ensembl transcript info and save to a TSV file."
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help="Path to the output TSV file"
    )
    parser.add_argument(
        '--host',
        default='http://www.ensembl.org',
        help="Ensembl BioMart host (default: http://www.ensembl.org)"
    )
    parser.add_argument(
        '--dataset',
        default='mmusculus_gene_ensembl',
        help="BioMart dataset name (default: mmusculus_gene_ensembl)"
    )
    args = parser.parse_args()

    # Connect to Ensembl BioMart
    try:
        server = Server(host=args.host)
        mart = server.marts['ENSEMBL_MART_ENSEMBL']
        ds = mart.datasets[args.dataset]
    except Exception as e:
        sys.exit(f"Error connecting to Ensembl BioMart: {e}")

    # Define which columns to fetch
    attrs = [
        'ensembl_transcript_id',
        'ensembl_gene_id',
        'external_gene_name',
        'transcript_biotype',
        'gene_biotype',
    ]

    # Download the table
    print(f"Downloading {len(attrs)} attributes from {args.dataset}…")
    all_tx = ds.query(attributes=attrs)
    print(f"Retrieved {len(all_tx)} records.")

    # Write to TSV
    print(f"Writing to {args.output}…")
    all_tx.to_csv(args.output, sep='\t', index=False)
    print("Done.")

if __name__ == '__main__':
    main()

