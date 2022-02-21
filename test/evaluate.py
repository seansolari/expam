#!/usr/bin/env python

import sys

import numpy as np
import pandas as pd


def main(result_dir, truth_dir, combine=True):
    # Classifier results.
    df = pd.read_csv(result_dir, sep="\t",
                     names=['taxid', 'c_perc', 'c_cumul', 'c_count',
                            's_perc', 's_cumul', 's_count',
                            'rank', 'lineage'])
    df['t_cumul'] = df['s_cumul'] + df['c_cumul'] if combine else df['c_cumul']

    # Truth set.
    truth_df = pd.read_csv(truth_dir, sep="\t",
                           names=['taxid', 'counts', 'abundance', 'rank', 'scientific name'])
    truth_df['taxid'] = truth_df['taxid'].apply(str)

    # Take subset of good results.
    taxa = truth_df['taxid'].tolist()
    df_subset = df[df['taxid'].isin(taxa)]

    # Total number of reads classified at species level.
    n_species = np.sum(df[df['rank'] == 'species']['t_cumul'])
    precision = np.sum(df_subset['t_cumul']) / n_species
    recall = np.sum(df_subset['t_cumul']) / 500000

    print("Precision: %f" % precision)
    print("Recall: %f" % recall)


if __name__ == "__main__":
    result_dir = sys.argv[1]
    truth_dir = sys.argv[2]

    main(result_dir, truth_dir)
