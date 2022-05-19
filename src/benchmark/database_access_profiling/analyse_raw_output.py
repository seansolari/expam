#!/usr/bin/env python3
import os
from typing import List

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


def _count_file_lines(file):
    with open(file, 'r') as f:
        for i, line in enumerate(f):
            pass

    return i


def load_histogram_counts(file):
    data = []

    with open(file, 'r') as f:
        for line in f:
            parts = line.strip().split(',')

            if parts:
                data.append(tuple(int(v) for v in parts))

    return data


def _format_compressed_form(node: str, count_str: str):
    return int(node.lstrip("p")), int(count_str)


def yield_output(file):
    n_lines = _count_file_lines(file)

    with open(file, 'r') as f:
        for line in tqdm(f, total=n_lines):
            data = line.strip().split('\t')

            if data[0] != "U" and len(data) == 5:
                yield data[0], (_format_compressed_form(*v.split(":")) for v in data[4].split(" "))


def histogram_database_access(file):
    c_counts = {}
    s_counts = {}
    max_node = 0

    for classification_type, output in yield_output(file):
        store = c_counts if classification_type == "C" else s_counts

        for node, count in output:
            max_node = max(node, max_node)

            try:
                store[node] += count
            except KeyError:
                store[node] = count

    return c_counts, s_counts


def sort_dict(node_to_counts: dict) -> List[tuple]:
    as_list = list(node_to_counts.items())
    sorted_list = sorted(as_list, key=lambda v: v[0])
    return sorted_list


def record_counts(file, c_counts: List[tuple], s_counts: List[tuple]):
    out_dir = os.path.basename(file).rstrip(".csv")
    
    try:
        os.mkdir(os.path.join('data', out_dir))
    except OSError:
        pass

    c_file = os.path.join('data', out_dir, 'c_counts.csv')
    s_file = os.path.join('data', out_dir, 's_counts.csv')

    with open(c_file, 'w') as f:
        for v in c_counts:
            f.write("%d,%d\n" % v)

    with open(s_file, 'w') as f:
        for v in s_counts:
            f.write("%d,%d\n" % v)

    return os.path.join('data', out_dir)


def weighted_percentile(x: np.ndarray, weights: np.ndarray, quantiles: list):
    quantiles = np.array(quantiles) / sum(quantiles)

    weighted_quantiles = np.cumsum(weights) - 0.5 * weights
    weighted_quantiles /= np.sum(weights)

    return np.interp(quantiles, weighted_quantiles, x)


def draw_histogram(hist_tuples: List[tuple], out_file: str):
    x, weights = [], []
    for v in hist_tuples:
        x.append(v[0])
        weights.append(v[1])

    x = np.array(x)
    weights = np.array(weights)

    q25, q75 = weighted_percentile(x, weights, [25, 75])
    bin_width = (q75 - q25) * len(x) ** (-1/3)
    bins = round((x.max() - x.min()) / bin_width)

    plt.hist(x, bins=bins, weights=weights, density=True)
    plt.grid()
    plt.savefig(out_file)
    plt.close()


def main():
    in_dir = "raw"
    files = [os.path.join(in_dir, f) for f in os.listdir(in_dir) if not f.startswith(".")]

    for file in files:
        count_dir = os.path.join('data', os.path.basename(file).rstrip(".csv"))
        c_counts_dir = os.path.join(count_dir, "c_counts.csv")
        s_counts_dir = os.path.join(count_dir, "s_counts.csv")

        if os.path.exists(c_counts_dir) and os.path.exists(s_counts_dir):
            print("Loading counts for %s..." % file)
            c_counts_list = load_histogram_counts(c_counts_dir)
            s_counts_list = load_histogram_counts(s_counts_dir)
        else:
            print("Analysing output for %s..." % file)
            c_counts, s_counts = histogram_database_access(file)

            c_counts_list = sort_dict(c_counts)
            s_counts_list = sort_dict(s_counts)

            count_dir = record_counts(file, c_counts_list, s_counts_list)

        print("Histogramming c...")
        draw_histogram(c_counts_list, os.path.join(count_dir, 'c.png'))
        print("Histogramming s...")
        draw_histogram(s_counts_list, os.path.join(count_dir, 's.png'))



if __name__ == "__main__":
    main()
