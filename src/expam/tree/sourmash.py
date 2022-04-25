from functools import partial
from multiprocessing import Pool
import os

import numpy as np                     # Pairwise distances.
from expam.sequences import get_sequence

try:
    from sourmash import signature as sm_signature
    import sourmash

except ModuleNotFoundError:
    print("Could not import sourmash! Install sourmash to use this feature.")


def make_signature(genome_url, k, s):
    print("Generating signature for %s." % genome_url)

    minhash = sourmash.MinHash(n=0, ksize=k, scaled=s)
    sequence = get_sequence(genome_url)
    minhash.add_sequence(sequence, True)

    sig = sm_signature.SourmashSignature(
        minhash, name=os.path.basename(genome_url), filename=genome_url)

    return sig


def make_signatures(n_processes, genome_dirs, sig_dir, k, s):
    print("Creating signatures!")

    if not os.path.exists(sig_dir):
        os.makedirs(sig_dir)

        print("Created directory %s." % sig_dir)

    else:
        # Check for any genomes already processed.
        sig_names = [
            os.path.basename(f).replace(".sour", "")
            for f in os.listdir(sig_dir)
        ]

        genome_dirs = [
            genome_dir
            for genome_dir in genome_dirs
            if os.path.basename(genome_dir) not in sig_names
        ]

    # Calculate chunksize.
    chunksize = 100 if (len(genome_dirs) // n_processes) > 100 else (len(genome_dirs) // n_processes)

    print("Starting process pool...")

    with Pool(n_processes) as mp_pool:
        sigs = list(mp_pool.imap_unordered(
            partial(make_signature, k=k, s=s),
            genome_dirs,
            chunksize=chunksize
        ))

    print("Process pool completed!")

    # Save sigs to disk.
    sig_data = sm_signature.save_signatures(sigs).decode('utf-8')

    with open(os.path.join(sig_dir, "mysigs.sig"), "w") as f:
        f.write(sig_data)


def make_distances(sig_dir, n_processes, matrix_dir):
    if not os.path.exists(sig_dir):
        raise Exception("No such signatures at %s!" % sig_dir)

    with open(os.path.join(sig_dir, "mysigs.sig"), "r") as f:
        sig_json = f.read()

    signatures = list(sm_signature.load_signatures(sig_json))
    distance_matrix = sourmash.compare.compare_all_pairs(
        siglist=signatures,
        ignore_abundance=True,
        downsample=False,
        n_jobs=n_processes
    )

    np.fill_diagonal(distance_matrix, 0)
    names = [
        sig.name[:-4] if sig.name[-4:] == ".fna" else sig.name
        for sig in signatures
    ]
    to_phylip_matrix(distance_matrix, names, matrix_dir)


def to_phylip_matrix(distance_matrix, names, out_dir):
    phylip_string = "%d\n" % len(distance_matrix)

    for i, row in enumerate(distance_matrix):
        phylip_string += "%s\t%s\n" % (
            names[i],
            "\t".join(
                map(str, row)
            )
        )

    with open(out_dir, "w") as f:
        f.write(phylip_string)
