#!/usr/bin/env python3
"""
run_diff_expr.py

Filter RNA-seq count data, compute within-group ratios, adjust to (0,1), then
run a Beta-mixture differential expression model in PyMC3.
"""

import argparse
import pandas as pd
import numpy as np
import pymc3 as pm
from Bio.Data import CodonTable  # ensure Biopython is installed if used elsewhere
from Bio.Seq import Seq          # only needed if you re-use reverse_translate
import theano.tensor as tt       # PyMC3 dependency
import arviz as az               # optional, for InferenceData conversion

def filter_rna_seq_data_V2(rna_counts, condition_a, condition_b):
    cols = condition_a + condition_b
    sub = rna_counts.iloc[:, cols].copy()
    sub.index = rna_counts['ids']
    a_ok = sub.iloc[:, :len(condition_a)].gt(0).sum(axis=1) >= 1
    b_ok = sub.iloc[:, len(condition_a):].gt(0).sum(axis=1) >= 1
    return rna_counts[a_ok & b_ok]

def calculate_column_ratios_within_promoter_group(df, group_cols, value_cols):
    df2 = df.copy()
    def ratio_calc(x): return x / x.sum() if x.sum() else 0
    for col in value_cols:
        df2[col] = df2.groupby(group_cols)[col].transform(ratio_calc)
    return df2[["ids"] + value_cols]

def adjust_data(data, lower_bound=0.001, upper_bound=0.999):
    adj = np.where(data <= 0, lower_bound, data)
    adj = np.where(adj >= 1, upper_bound, adj)
    return adj

def run_differential_expression_model_mixture_psi_v1(dataA, dataB, n_samples, n_tune, cores):
    n = dataA.shape[0]
    with pm.Model() as model:
        alpha = pm.Gamma('alpha', 5, 1.0, shape=n)[:, None]
        beta  = pm.Gamma('beta',  5, 1.0, shape=n)[:, None]
        pm.Beta('obs_a', alpha=alpha, beta=beta, observed=dataA)
        pm.Beta('obs_b', alpha=alpha, beta=beta, observed=dataB)
        trace = pm.sample(n_samples, tune=n_tune, cores=cores, random_seed=42)
    return trace

def calculate_differential_expression_mixture_v1(trace, index):
    diff = trace['beta'] - trace['alpha']
    fc   = np.exp(diff)
    hdi  = pm.stats.hdi(fc, hdi_prob=0.95)
    df = pd.DataFrame({
        'beta_a_mean': np.mean(trace['alpha'], axis=0),
        'beta_b_mean': np.mean(trace['beta'],  axis=0),
        'difference_mean': diff.mean(axis=0),
        'fold_change_mean': fc.mean(axis=0),
        'fold_change_hdi_low':  hdi[:, 0],
        'fold_change_hdi_high': hdi[:, 1],
        'p_diff_gt0': np.mean(diff > 0, axis=0),
        'p_diff_lt0': np.mean(diff < 0, axis=0),
    }, index=index)
    return df

def main():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument('--asevents',  required=True, help="TSV of ASE events (with 'ids')")
    p.add_argument('--counts',    required=True, help="TSV of count data (with 'ids')")
    p.add_argument('--condA',     required=True, nargs='+', type=int,
                   help="Zero-based column indices for condition A")
    p.add_argument('--condB',     required=True, nargs='+', type=int,
                   help="Zero-based column indices for condition B")
    p.add_argument('--groups',    required=True, nargs='+',
                   help="Grouping columns, e.g. gene_id promoter_group")
    p.add_argument('--nsamples',  type=int, default=100, help="MCMC samples")
    p.add_argument('--ntune',     type=int, default=100, help="MCMC tuning steps")
    p.add_argument('--cores',     type=int, default=4, help="Parallel cores for sampling")
    p.add_argument('--out',       required=True, help="Output TSV for DE results")
    args = p.parse_args()

    # Load
    ase = pd.read_csv(args.asevents, sep='\t').set_index('ids')
    cnt = pd.read_csv(args.counts,  sep='\t').set_index('ids')
    rna = pd.merge(cnt, ase, left_index=True, right_index=True).reset_index()

    # Filter
    filt = filter_rna_seq_data_V2(rna, args.condA, args.condB)

    # Ratios
    val_cols = [filt.columns[i] for i in (args.condA + args.condB)]
    ratios = calculate_column_ratios_within_promoter_group(filt, args.groups, val_cols)
    ratios = ratios.set_index('ids')

    # Adjust
    dataA = adjust_data(ratios.iloc[:, :len(args.condA)].values)
    dataB = adjust_data(ratios.iloc[:, len(args.condA):].values)

    # Model
    trace = run_differential_expression_model_mixture_psi_v1(
        dataA, dataB,
        n_samples=args.nsamples,
        n_tune=args.ntune,
        cores=args.cores
    )

    # Summarize
    results = calculate_differential_expression_mixture_v1(trace, ratios.index)
    results.to_csv(args.out, sep='\t')
    print(f"Wrote DE results to {args.out}")

if __name__ == '__main__':
    main()

