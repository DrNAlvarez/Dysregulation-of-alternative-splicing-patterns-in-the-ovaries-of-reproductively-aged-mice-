#!/usr/bin/env python3
"""
translate_cDNA.py

For each transcript in a cDNA FASTA, reverse-translate the first 5 amino acids
of its matching protein sequence and search for the correct ORF in the cDNA.
Outputs a TSV with transcript_id, gene_id, AA_sequence, start, end.
"""

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
import pandas as pd
from multiprocessing import Pool
from tqdm import tqdm

def get_possible_codons(amino_acid):
    table = CodonTable.unambiguous_dna_by_name["Standard"]
    return [codon for codon, aa in table.forward_table.items() if aa == amino_acid]

def reverse_translate(protein_seq):
    seqs = ['']
    for aa in protein_seq:
        codons = get_possible_codons(aa)
        seqs = [s + c for s in seqs for c in codons]
    return seqs

def check_cdna_translation_parallel(transcript_id, protein_seqs, cdna_seqs):
    results = []
    if transcript_id not in protein_seqs or transcript_id not in cdna_seqs:
        return results

    cdna = cdna_seqs[transcript_id]
    prot = protein_seqs[transcript_id]
    prefix = prot[:5]
    starts = reverse_translate(prefix)

    for dna_start in starts:
        pos = cdna.find(dna_start)
        if pos != -1:
            orf = Seq(cdna[pos:]).translate(to_stop=True)
            stop = pos + len(str(orf)) * 3
            gene_id = transcript_id.split("_")[1] if "_" in transcript_id else ""
            results.append({
                "id": transcript_id,
                "gene_id": gene_id,
                "transcript_id": transcript_id,
                "AA_sequence": str(orf),
                "start": pos,
                "end": stop
            })
            break
    return results

def main():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("protein_fasta", help="FASTA of protein sequences")
    p.add_argument("cdna_fasta", help="FASTA of cDNA sequences")
    p.add_argument("output_tsv", help="Where to write the results TSV")
    p.add_argument("-p", "--processes", type=int, default=4,
                   help="Number of parallel workers (default: 4)")
    args = p.parse_args()

    # Load sequences
    print(f"Loading protein sequences from {args.protein_fasta}…")
    prot_seqs = {r.id: str(r.seq) for r in SeqIO.parse(args.protein_fasta, "fasta")}
    print(f"Loaded {len(prot_seqs)} protein entries.")
    print(f"Loading cDNA sequences from {args.cdna_fasta}…")
    cdna_seqs = {r.id: str(r.seq) for r in SeqIO.parse(args.cdna_fasta, "fasta")}
    print(f"Loaded {len(cdna_seqs)} cDNA entries.")

    transcript_ids = list(cdna_seqs.keys())
    tasks = [(tid, prot_seqs, cdna_seqs) for tid in transcript_ids]

    # Parallel scan
    print(f"Searching ORFs in {len(tasks)} transcripts with {args.processes} processes…")
    with Pool(args.processes) as pool:
        raw = pool.starmap(check_cdna_translation_parallel, tqdm(tasks, desc="ORF scan"))

    # Flatten and save
    flat = [hit for sub in raw for hit in sub]
    df = pd.DataFrame(flat, columns=["id","gene_id","transcript_id","AA_sequence","start","end"])
    df.to_csv(args.output_tsv, sep="\t", index=False)
    print(f"Wrote {len(df)} results to {args.output_tsv}")

if __name__ == "__main__":
    main()

