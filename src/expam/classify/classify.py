from math import floor
from multiprocessing import shared_memory
import os
import re
import shutil
from typing import List, Mapping, Tuple, Union

import numpy as np
import pandas as pd

from expam.classify import ResultsPathConfig
from expam.classify.config import create_tax_results, load_results_config
from expam.classify.process import ExpamClassifierProcesses
from expam.classify.taxonomy import TaxonomyNCBI
from expam.database import EXPAM_TEMP_EXT, TIMEOUT, UNION_RATIO, FileLocationConfig
from expam.database.config import load_database_config
from expam.database.db import SharedDB
from expam.process.classify import ClassifyWorker
from expam.process.manager import ControlCenter
from expam.process.reads import ReadWorker
from expam.tree.tree import Index
from expam.utils import ls, yield_csv


def run_classifier(
    read_paths: str,
    out_dir: str,
    db_dir: str,
    k: int,
    n: int,
    phylogeny_path: str,
    keys_shape: Tuple[int],
    values_shape: Tuple[int],
    logging_dir: str,
    taxonomy: bool = False, cutoff: float = 0.0, groups: Union[None, List[tuple]] = None,
    keep_zeros: bool = False, cpm: float = 0.0, use_node_names: bool = True,
    phyla: bool = False, name_taxa: Union[Mapping[str, str], None] = None,
    colour_list: List[str] = None, paired_end: bool = False, alpha: float = 1.0,
    log_scores: bool = False, itol_mode: bool = False,
    sep: str = ","
):
    output_config: ResultsPathConfig = load_results_config(out_dir, create=True)
    database_config: FileLocationConfig = load_database_config(db_dir)

    # Load the kmer dictionary.
    print("Loading the map and phylogeny.\n")
    db_kmers: SharedDB = SharedDB(db_dir, keys_shape, values_shape)

    try:
        # Get name mapping.
        index, phylogeny_index = name_to_id(phylogeny_path)

        # Load LCA Matrix.
        lca_matrix = np.load(database_config.lca_matrix)

        # Initialise the distribution.
        d = Distribution(
            k=k,
            kmer_db=db_kmers,
            index=phylogeny_index,
            lca_matrix=lca_matrix,
            read_paths=read_paths,
            out_dir=out_dir,
            logging_dir=logging_dir,
            temp_dir=output_config.temp,
            paired_end=paired_end,
            alpha=alpha,
            sep=sep
        )

        # Run the classification with n processes.
        d.run(n)

        # Load accession ids.
        sequence_ids = list(yield_csv(database_config.accession_id))

        # Update accession ids in map.
        for (fname, accn_id, tax_id) in sequence_ids:
            if accn_id != "None":
                phylogeny_index[fname].accession_id = accn_id

            if tax_id != "None":
                phylogeny_index[fname].taxid = tax_id

        # Output results.
        results = ClassificationResults(
            index=index,
            phylogeny_index=phylogeny_index,
            results_config=output_config,
            groups=groups,
            keep_zeros=keep_zeros,
            cutoff=cutoff,
            cpm=cpm,
            use_node_names=use_node_names,
            phyla=phyla,
            name_taxa=name_taxa,
            colour_list=colour_list,
            log_scores=log_scores,
            sep=sep
        )

        if taxonomy:
            # Attempt to update taxon ids.
            tax_obj: TaxonomyNCBI = TaxonomyNCBI(database_config)
            tax_obj.accession_to_taxonomy()

            name_to_lineage, taxon_to_rank = tax_obj.load_taxonomy_map()
            results.to_taxonomy(name_to_lineage, taxon_to_rank)

        results.draw_results(itol_mode=itol_mode)
    finally:
        # Close shared memory.
        keys_shm = shared_memory.SharedMemory(
            name=db_kmers.keys_shm_name,
            create=False
        )
        values_shm = shared_memory.SharedMemory(
            name=db_kmers.values_shm_name,
            create=False
        )
        keys_shm.unlink()
        keys_shm.close()
        values_shm.unlink()
        values_shm.close()

        shutil.rmtree(output_config.temp)


def name_to_id(phylogeny_path: str):
    # Load Index from Newick file.
    with open(phylogeny_path, "r") as f:
        newick_str = f.read().strip()

    _, index_map = Index.from_newick(newick_str)

    index = {}
    for i, loc in enumerate(index_map.pool):
        index[loc.name] = i

    return index, index_map


class Distribution:
    _c_perc = "Cumulative Classified Percentage"
    _c_cumul = "Cumulative Classified Count"
    _c_count = "Raw Classified Count"
    _s_perc = "Cumulative Split Percentage"
    _s_cumul = "Cumulative Split Count"
    _s_count = "Raw Split Count"

    def __init__(self, k, kmer_db, index: Index, lca_matrix, read_paths, out_dir, temp_dir, logging_dir, alpha,
                 keep_zeros=False, cutoff=0.0, cpm=0.0, paired_end=False, sep=","):
        """

        :param k: Int
        :param kmer_db: SharedDB
        :param lca_matrix: np.ndarray
        :param url: str
        :param out_dir: str
        """
        # Set up processing references.
        self.k = k
        self.alpha = alpha
        self.kmer_db = kmer_db
        self.lca_matrix = lca_matrix
        self.out_dir = out_dir
        self.temp_dir = temp_dir
        self.logging_dir = logging_dir

        self.index: Index = index
        self.node_names = [node.name if i > 0 else "unclassified" for i, node in enumerate(index.pool)]

        self.keep_zeros = keep_zeros
        self.cutoff = cutoff
        self.cpm = cpm

        # Url containing reads.
        self.results_config: ResultsPathConfig = load_results_config(out_dir, create=True)
        self.read_paths = read_paths
        self.paired_end = paired_end

        self.SEP = sep

    def run(self, n):
        """
        Begin processing reads in self.url.
        We have some collection of processes accessing reads from storage
        and then another collection of processes that analyse the reads.

        """
        file_queue = self.prepare_queue(self.paired_end)
        n_classifiers = floor(UNION_RATIO * n)
        n_readers = n - n_classifiers

        mp_config = {
            "name": "expam",
            "phases": ["classify"],
            "layers": {
                "class": ReadWorker,
                "class_args": {
                    "paired_mode": self.paired_end,
                    "n_classifiers": n_classifiers,
                    "logging_dir": self.logging_dir,
                },
                "parent": ControlCenter,
                "child": {
                    "class": ClassifyWorker,
                    "class_args": {
                        "k": self.k,
                        "kmer_db": self.kmer_db,
                        "lca_matrix": self.lca_matrix,
                        "out_dir": self.results_config.phy,
                        "temp_dir": self.temp_dir,
                        "n_readers": n_readers,
                        "logging_dir": self.logging_dir,
                        "alpha": self.alpha
                    },
                    "parent": ReadWorker,
                    "child": None,
                    "n": n_classifiers
                },
                "n": n_readers
            },
            "timeout": TIMEOUT
        }

        process_network = ExpamClassifierProcesses.from_method_dict(self.logging_dir, mp_config)

        # Load jobs.
        for file_dir in file_queue:
            process_network.load_job(phase="classify", job=file_dir)

        process_network.run()

        # Process and combine outputs.
        output_files = self.combine_raw_read_output()
        self.create_sample_summaries(output_files)

    def combine_raw_read_output(self):
        temporary_files = ls(self.temp_dir, ext=".reads_%s" % EXPAM_TEMP_EXT)
        base_file_names = {
            re.match(r'(\S+)_\d+.reads_%s' % EXPAM_TEMP_EXT, os.path.basename(file_name)).group(1)
            for file_name in temporary_files
        }

        result_files = []
        for sample_name in base_file_names:
            results_file_name = os.path.join(self.results_config.phy_raw, sample_name + ".csv")
            raw_files = [file for file in temporary_files if re.match(r'.*{base}_\d+.reads_{ext}$'.format(base=sample_name, ext=EXPAM_TEMP_EXT), file)] 
            
            self.concatenate_files(results_file_name, *raw_files)
            result_files.append((results_file_name, sample_name))

        return result_files

    def create_sample_summaries(self, result_files):
        # Sample-wise summary of counts in kraken-like format.
        # Create matrix of counts.
        sample_names = [sample_name for _, sample_name in result_files]
        sldf = pd.DataFrame(0, index=self.node_names, columns=sample_names)
        mldf = pd.DataFrame(0, index=self.node_names, columns=sample_names)

        # Import counts into DataFrame.
        for result_file, sample_name in result_files:
            slser = sldf.loc[:, sample_name]
            mlser = mldf.loc[:, sample_name]

            with open(result_file, 'r') as f:
                for line in f:
                    parts = line.split("\t")
                    idx = int(parts[2].lstrip("p"))

                    if parts[0] == "S":
                        mlser.iloc[idx] += 1
                    else:
                        slser.iloc[idx] += 1

        # Create sample-wise summaries.
        for _, sample_name in result_files:
            df = pd.DataFrame(0, index=self.node_names, columns=[self._c_perc, self._c_cumul, self._c_count, self._s_perc, self._s_cumul, self._s_count])
            # Insert single-lineage classifications.
            df.loc[:, self._c_cumul] = sldf.loc[:, sample_name]
            df.loc[:, self._c_count] = sldf.loc[:, sample_name]
            # Insert split-lineage classifications.
            df.loc[:, self._s_cumul] = mldf.loc[:, sample_name]
            df.loc[:, self._s_count] = mldf.loc[:, sample_name]

            self.make_cumulative(df)
            total_reads = df.loc[:, self._c_count].sum() + df.loc[:, self._s_count].sum()

            # Employ cutoff.
            if self.cutoff > 0 or self.cpm > 0:
                cutoff = max((total_reads / 1e6) * self.cpm, self.cutoff)
                mask = (df.loc[:, self._c_cumul] < cutoff) & (df.loc[:, self._s_cumul])
                mask[0] = True  # Keep unclassified counts.
                df.loc[mask, :] = 0

            # Check for removal of zero rows.
            if not self.keep_zeros:
                (df, ) = self.remove_zeros(df)

            # Insert percentage form.
            self.insert_percentages(df, total_reads)

            # Write to disk.
            out_file = os.path.join(self.results_config.phy, sample_name + ".csv")
            self.save_df(df, out_file)

        # Remove null rows.
        if not self.keep_zeros:
            (sldf, mldf) = self.remove_zeros(sldf, mldf)

        # Drop unclassified from split DataFrame.
        try:
            mldf.drop(index="unclassified", inplace=True)
        except KeyError:    # If no unclassified, this will have already been removed.
            pass

        # Save counts to disk.
        self.save_df(sldf, self.results_config.phy_classified)
        self.save_df(mldf, self.results_config.phy_split)

    def save_df(self, df, out_file):
        self.fix_df_index(df)
        df.to_csv(out_file, sep=self.SEP)

    def fix_df_index(self, df):
        fixer = lambda name: 'p' + name if (name != "unclassified") and (self.index[name].type == "Branch") else name
        df.index = [fixer(n) for n in df.index]

    def make_cumulative(self, df):
        cumulative_columns = [self._c_cumul, self._s_cumul]

        for node in self.index.pool[1:][::-1]:
            df.loc[node.name, cumulative_columns] += df.loc[node.children, cumulative_columns].sum(axis=0)

    def insert_percentages(self, df, total_counts):
        df.loc[:, self._c_perc] = round(100 * df.loc[:, self._c_cumul] / total_counts, 3).map(str) + "%"
        df.loc[:, self._s_perc] = round(100 * df.loc[:, self._s_cumul] / total_counts, 3).map(str) + "%"

    @staticmethod
    def remove_zeros(*dfs):
        mask = None
        for df in dfs:
            non_zero = df.apply(lambda r: ((r > 0).any()) | (r.name == 'unclassified'), axis=1)
            mask = non_zero if mask is None else mask | non_zero

        for df in dfs:
            yield df.loc[mask, :]

    @staticmethod
    def concatenate_files(dest, *files):
        with open(dest, 'w') as out_f:
            for file in files:
                # Write data to concatenated file.
                with open(file, 'r') as in_f:
                    out_f.write(in_f.read().strip() + "\n")

                # Remove temporary file.
                os.remove(file)

    @staticmethod
    def best_match(file_name, file_list):
        def string_intersect(str_one, str_two):
            l = min(len(str_one), len(str_two))

            for i in range(l):
                if str_one[i] != str_two[i]:
                    return str_one[:i]
            else:
                return str_one[:l]

        match_ind = 0
        best_match = string_intersect(file_name, file_list[0])

        for i in range(1, len(file_list)):
            current_match = string_intersect(file_name, file_list[i])

            if len(current_match) > len(best_match):
                best_match = current_match
                match_ind = i

        return match_ind

    def prepare_queue(self, paired_end=False):
        file_queue = []

        for path in self.read_paths:
            # Create a list of all the files containing reads to be processed.
            file_list = [
                file_name
                for file_name in ls(path)
                if re.match(r"(\S+)\.(?:fa|fna|ffn|fq|fasta|fastq)(?:\.tar)*(?:\.gz)*$", file_name)
            ]

            if paired_end:  # This implies path is a folder, not a file.
                while file_list:
                    next_file = file_list.pop()

                    # Find buddy.
                    match_ind = self.best_match(next_file, file_list)
                    paired_file = file_list.pop(match_ind)

                    file_queue.append((next_file, paired_file))
            else:
                for file_name in file_list:
                    file_queue.append((file_name, ))

        return file_queue

    @staticmethod
    def _initialise_results(names, unclassified=True):
        """
        Results are returned in spectral form:
            coord --> # hits at coord.

        :param names: List - node names in the tree.
        :return: pd.Series()

        """
        index = (unclassified * ["unclassified"]) + names
        return pd.Series(
            data=np.zeros(len(names) + unclassified, dtype=int),
            index=index
        )


class ClassificationResults:
    _c_perc = "Cumulative Classified Percentage"
    _c_cumul = "Cumulative Classified Count"
    _c_count = "Raw Classified Count"
    _s_perc = "Cumulative Split Percentage"
    _s_cumul = "Cumulative Split Count"
    _s_count = "Raw Split Count"
    _rank = "Rank"
    _scientific_name = "Scientific Name"

    def __init__(self, index, phylogeny_index, results_config, groups=None, keep_zeros=False, cutoff=0.0, cpm=0.0,
                 use_node_names=True, phyla=False, name_taxa=None, colour_list=None, circle_scale=1.0,
                 log_scores=False, sep=","):
        self.index = index
        self.phylogeny_index: Index = phylogeny_index
        self.results_config: ResultsPathConfig = results_config

        self.groups = groups  # [(colour, (name1, name2, ...)), ...]
        self.keep_zeros = keep_zeros

        self.cutoff = cutoff  # Cutoff by counts.
        self.cpm = cpm  # Cutoff by counts per million.
        self.use_node_names = use_node_names

        self.phyla = phyla
        self.colour_list = colour_list
        self.name_taxa = name_taxa
        self.circle_scale = circle_scale
        self.log_scores = log_scores

        if self.phyla and self.name_taxa is None:
            raise Exception("Inclusion of phyla requires mapping of names to taxonomic lineages!")

        self.taxid_immediate_children = {"1": set()}  # Map from tax_id -> immediate children.
        self.tax_id_pool = ["1"]  # Children must appear later than parent this list.

        self.SEP = sep

    def to_taxonomy(self, name_to_lineage, taxon_to_rank):
        col_names = [self._c_perc, self._c_cumul, self._c_count, self._s_perc, self._s_cumul, self._s_count, self._rank, self._scientific_name]

        class_counts = pd.read_csv(self.results_config.phy_classified, sep=self.SEP, index_col=0, header=0)
        split_counts = pd.read_csv(self.results_config.phy_split, sep=self.SEP, index_col=0, header=0)

        # Get rid of phylogenetically printed node names.
        def fix_index(df: pd.DataFrame):
            df.index = [
                index.lstrip("p")
                if index not in self.phylogeny_index._pointers
                else index
                for index in df.index
            ]

        fix_index(class_counts)
        fix_index(split_counts)

        samples = class_counts.columns.tolist()
        phylogeny_nodes = split_counts.index.tolist()

        # Map phylogeny index to taxonomy.
        phylogeny_to_taxonomy, tax_id_data = {}, {}
        for phylogeny_node in phylogeny_nodes:
            tax_id, tax_id_map = self.map_phylogeny_node(phylogeny_node, name_to_lineage, taxon_to_rank)

            phylogeny_to_taxonomy[phylogeny_node] = tax_id
            tax_id_data.update(tax_id_map)

        # Get relevant taxon ids.
        tax_ids, ranks, lineages = ["1"], ["root"], [""]
        for tax_id in self.tax_id_pool[1:]:
            rank, lineage = tax_id_data[tax_id]

            tax_ids.append(tax_id)
            ranks.append(rank)
            lineages.append(lineage)

        # Make taxonomic output folder structure.
        create_tax_results(self.results_config)

        for sample_name in samples:
            #
            # Map summary output to taxonomy.
            #

            df = pd.DataFrame(0, index=["unclassified"] + tax_ids, columns=col_names)
            total_counts = class_counts[sample_name].sum() + split_counts[sample_name].sum()

            df.loc[tax_ids[0]:, self._rank] = ranks
            df.loc[tax_ids[0]:, self._scientific_name] = lineages

            # Manually insert unclassified counts.
            df.loc["unclassified", self._c_count] = class_counts.loc["unclassified", sample_name]
            df.loc["unclassified", self._c_cumul] = df.loc["unclassified", self._c_count]
            df.loc["unclassified", self._c_perc] = round(100 * df.loc["unclassified", self._c_cumul] / total_counts, 3)

            for col, counts_df in ((self._c_count, class_counts), (self._s_count, split_counts)):
                for node in counts_df.index[1:]:
                    df.loc[phylogeny_to_taxonomy[node], col] += counts_df.loc[node, sample_name]

            # Make counts cumulative.
            # We may not be able to use groupby, but we can iterate through the taxid
            # pool once (in hierarchical order).
            for tax_id in self.tax_id_pool[::-1]:
                for count_col, cumul_col in (
                    (self._c_count, self._c_cumul),
                    (self._s_count, self._s_cumul)
                ):
                    df.loc[tax_id, cumul_col] = df.loc[tax_id, count_col] + df.loc[list(self.taxid_immediate_children[tax_id]), cumul_col].sum(axis=0)

            # Format percentage signs.
            df.loc[:, self._c_perc] = round(100 * df.loc[:, self._c_cumul] / total_counts, 3).map(str) + "%"
            df.loc[:, self._s_perc] = round(100 * df.loc[:, self._s_cumul] / total_counts, 3).map(str) + "%"

            # Employ cutoff.
            cutoff = max(self.cutoff, (total_counts / 1e6) * self.cpm)
            df = df[(df[self._c_cumul] > cutoff) | (df[self._s_cumul] > cutoff) | (df.index == 'unclassified')]
            df.to_csv(os.path.join(self.results_config.tax, sample_name + ".csv"), sep=self.SEP, header=True)

            #
            # Map raw read output to taxonomy.
            #
            raw_counts_file = os.path.join(self.results_config.phy_raw, sample_name + ".csv")
            output_file = os.path.join(self.results_config.tax_raw, sample_name + '.csv')

            with open(raw_counts_file, 'r') as in_f, open(output_file, 'w') as out_f:
                for line in in_f:
                    read_data = line.strip().split("\t")[:-1]

                    if read_data[0] == 'U':
                        read_data[2] = 'NA'
                    else:
                        phy_id = int(read_data[2].lstrip('p'))
                        phylogeny_node = self.phylogeny_index.pool[phy_id].name
                        read_data[2] = phylogeny_to_taxonomy[phylogeny_node]

                    out_f.write("\t".join(read_data) + "\n")

    def map_phylogeny_node(self, node_name, name_to_lineage, taxon_to_rank):
        """
        Create a map from (node_name) --> (taxid, rank, taxonomic lineage).
        
        :param node_name: Name of clade.
        :param name_to_lineage: Dict[Accession ID, Tuple[Taxonomic Lineage]]
        :param taxon_to_rank: Dict[Taxonomic Name, Tuple[TaxID, Rank]]
        """
        leaves = self.get_leaves(node_name)
        lineages = [name_to_lineage[leaf] for leaf in leaves]
        common_lineage = self.tuple_intersect(*lineages)

        if not common_lineage:
            phy_tax_id, phy_rank = "1", "root"
        else:
            phy_tax_id, phy_rank = taxon_to_rank[common_lineage[-1]]

        taxid_rank_and_lineage = {}
        pool_index = 1
        for i in range(len(common_lineage)):
            # Store taxonomic data.
            parent_tax_id, parent_rank = taxon_to_rank[common_lineage[i]]
            taxid_rank_and_lineage[parent_tax_id] = [parent_rank, " ".join(common_lineage[:i + 1])]
            
            # Add top node as child of root.
            if i == 0:
                self.taxid_immediate_children["1"].add(parent_tax_id)

            # Set tax_id in pool - a list of taxids such that child taxids are necessarily further
            # in the last than any ancestors.
            try:
                pool_index = self.tax_id_pool.index(parent_tax_id) + 1
            except ValueError:
                self.tax_id_pool.insert(pool_index, parent_tax_id)
                pool_index += 1

            # Store immediate child of taxid.
            if i + 1 < len(common_lineage):
                child_tax_id, child_rank = taxon_to_rank[common_lineage[i + 1]]

                try:
                    self.taxid_immediate_children[parent_tax_id].add(child_tax_id)
                except KeyError:
                    self.taxid_immediate_children[parent_tax_id] = {child_tax_id}

            elif parent_tax_id not in self.taxid_immediate_children:
                self.taxid_immediate_children[parent_tax_id] = set()

        return phy_tax_id, taxid_rank_and_lineage

    def get_leaves(self, node_name):
        return self.phylogeny_index.get_child_leaves(node_name)

    @staticmethod
    def tuple_intersect(*tuples):
        intersection = ()

        length = min([len(tup) for tup in tuples])
        for i in range(length):
            item = tuples[0][i]

            for j in range(1, len(tuples)):
                if tuples[j][i] != item:
                    return intersection

            else:
                intersection += (item,)

        return intersection

    def draw_results(self, itol_mode=False):
        # Draw classified tree.
        self.phylogeny_index.draw_results(
            self.results_config.phy_classified,
            os.path.join(self.results_config.base, "phylotree_classified.pdf"),
            skiprows=[1],
            groups=self.groups,
            cutoff=self.cutoff,
            cpm=self.cpm,
            colour_list=self.colour_list,
            name_to_taxon=self.name_taxa,
            use_phyla=self.phyla,
            keep_zeros=self.keep_zeros,
            use_node_names=self.use_node_names,
            log_scores=self.log_scores,
            itol_mode=itol_mode,
            sep=self.SEP
        )

        # Draw unclassified tree.
        self.phylogeny_index.draw_results(
            self.results_config.phy_split,
            os.path.join(self.results_config.base, "phylotree_splits.pdf"),
            groups=self.groups,
            cutoff=self.cutoff,
            cpm=self.cpm,
            colour_list=self.colour_list,
            name_to_taxon=self.name_taxa,
            use_phyla=self.phyla,
            keep_zeros=self.keep_zeros,
            use_node_names=self.use_node_names,
            log_scores=self.log_scores,
            itol_mode=itol_mode,
            sep=self.SEP
        )
