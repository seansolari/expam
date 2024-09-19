# version: 2
from __future__ import annotations
import abc
from functools import partial
from math import floor
from multiprocess.pool import Pool   # Uses dill for serialising.
from multiprocessing import shared_memory
import os
import re
import shutil
from typing import List, Tuple, Union

import numpy as np
import pandas as pd

from expam.classify.config import create_tax_results, load_results_config
from expam.classify.process import ExpamClassifierProcesses
from expam.classify.taxonomy import TaxonomyInterface
from expam.database import EXPAM_TEMP_EXT, TIMEOUT, UNION_RATIO
from expam.database.config import load_database_config
from expam.database.db import SharedDB
from expam.process.classify import ClassifyWorker
from expam.process.manager import ControlCenter
from expam.process.reads import ReadWorker
from expam.tree.tree import Index
from expam.utils import abspath, die, ls


def run_classifier(read_paths: str, out_dir: str, db_dir: str, k: int, n: int, phylogeny_path: str,
                   keys_shape: Tuple[int], values_shape: Tuple[int], logging_dir: str,
                   taxonomy: bool = False, groups: Union[None, List[tuple]] = None,
                   keep_zeros: bool = False, cpm: float = 0.0, use_node_names: bool = True,
                   colour_list: List[str] = None, paired_end: bool = False, alpha: float = 1.0,
                   log_scores: bool = False, itol_mode: bool = False, flat_colour: bool = False,
                   sep: str = ","):
    database_config = load_database_config(db_dir)
    output_config = load_results_config(out_dir, create=True)
    # Load the kmer dictionary.
    print("Loading the map and phylogeny.\n")
    db_kmers = SharedDB(db_dir, keys_shape, values_shape)
    try:
        # Get name mapping.
        _, tree = Index.load_newick(phylogeny_path)
        # Load LCA Matrix.
        lca_matrix = np.load(database_config.lca_matrix)
        # Perform classification and phylogenetic summary.
        c = ClassificationRun(k, db_kmers, tree, lca_matrix, read_paths, out_dir,
                              output_config.temp, logging_dir, alpha,
                              keep_zeros=keep_zeros, cpm=cpm,
                              paired_end=paired_end, sep=sep, n_processes=n)
        phy = c.run()
        # Taxonomic results.
        if taxonomy:
            create_tax_results(output_config)
            my_ncbi = TaxonomyInterface(database_config, tree=tree)
            tax = TaxonomicResults.from_phylogenetic(phy, my_ncbi)
            tax.summarise(output_config.tax,
                          cpm=cpm,
                          sep=sep)
            tax.convert_raw(output_config.phy_raw, output_config.tax_raw, tree)
        # Plot phylogenetic results.
        phy.draw_results(output_config.phy, groups=groups, cpm=cpm,
                         colour_list=colour_list, keep_zeros=keep_zeros,
                         use_node_names=use_node_names, log_scores=log_scores,
                         itol_mode=itol_mode, flat_colour=flat_colour)
    finally:
        # Close shared memory.
        keys_shm = shared_memory.SharedMemory(name=db_kmers.keys_shm_name, create=False)
        values_shm = shared_memory.SharedMemory(name=db_kmers.values_shm_name, create=False)
        keys_shm.unlink()
        keys_shm.close()
        values_shm.unlink()
        values_shm.close()
        shutil.rmtree(output_config.temp)

class ClassificationRun:
    def __init__(self, k: int, kmer_db: SharedDB, tree: Index, lca_matrix: np.ndarray,
                 read_paths: List[str], out_dir: str, temp_dir: str, logging_dir: str, alpha: float,
                 keep_zeros: bool = False, cpm: float = 0.0, paired_end: bool = False,
                 sep: str = ",", n_processes: int = 1):
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

        self.tree = tree
        self.node_names: List[str] = [
            node.name if i > 0 else "unclassified"
            for i, node in enumerate(tree.pool)
        ]

        self.keep_zeros = keep_zeros
        self.cpm = cpm

        # Url containing reads.
        self.results_config = load_results_config(out_dir, create=True)
        self.read_paths = read_paths
        self.paired_end = paired_end

        self.SEP = sep
        self.n_processes = n_processes

    def run(self) -> PhylogeneticResults:
        """
        Begin processing reads in self.url.
        We have some collection of processes accessing reads from storage
        and then another collection of processes that analyse the reads.

        """
        ## sample processing
        #
        file_queue = self.prepare_queue(self.paired_end)
        if len(file_queue) == 0:
            die("Could not find any files in these directories: %s" % (", ".join(self.read_paths)))

        n_classifiers = floor(UNION_RATIO * self.n_processes)
        n_readers = self.n_processes - n_classifiers
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
        # load jobs
        for file_dir in file_queue:
            process_network.load_job(phase="classify", job=file_dir)
        process_network.run()

        ## data formatting
        #
        # concatenate process-wise output
        self.combine_raw_read_output()
        # create phylogenetic output
        phy = PhylogeneticResults.from_raw_read_output(self.results_config.phy_raw,
                                                       num_procs=self.n_processes,
                                                       tree=self.tree)
        phy.summarise(self.results_config.phy,
                      cpm=self.cpm,
                      sep=self.SEP)
        phy.to_csv(self.results_config.phy_classified,
                   self.results_config.phy_split,
                   self.keep_zeros)
        return phy

    def prepare_queue(self, paired_end: bool = False) -> List[Tuple[str, ...]]:
        file_queue = []
        for path in self.read_paths:
            # Create a list of all the files containing reads to be processed.
            file_list = [
                abspath(file_name)
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
    def best_match(file_name: str, file_list: List[str]) -> int:
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

    def combine_raw_read_output(self) -> None:
        temporary_files = ls(self.temp_dir, ext=".reads_%s" % EXPAM_TEMP_EXT)
        base_file_names = {
            re.match(r'(\S+)_\d+.reads_%s' % EXPAM_TEMP_EXT, os.path.basename(file_name)).group(1)
            for file_name in temporary_files
        }
        for sample_name in base_file_names:
            results_file_name = os.path.join(self.results_config.phy_raw, sample_name + ".csv")
            raw_files = [
                file for file in temporary_files
                if re.match(r'.*{base}_\d+.reads_{ext}$'.format(base=sample_name, ext=EXPAM_TEMP_EXT), file)
            ] 
            self.concatenate_files(results_file_name, *raw_files)
    
    @staticmethod
    def concatenate_files(dest, *files):
        with open(dest, 'w') as out_f:
            for file in files:
                # Write data to concatenated file.
                with open(file, 'r') as in_f:
                    out_f.write(in_f.read().strip() + "\n")
                # Remove temporary file.
                os.remove(file)

# base class =====================================================================

class _ClassificationOutput(metaclass=abc.ABCMeta):
    _c_perc = "Cumulative Classified Percentage"
    _c_cumul = "Cumulative Classified Count"
    _c_count = "Raw Classified Count"
    _s_perc = "Cumulative Split Percentage"
    _s_cumul = "Cumulative Split Count"
    _s_count = "Raw Split Count"
    _t_perc = "Cumulative Total Percentage"
    _t_cumul = "Cumulative Total Count"
    _t_count = "Raw Total Count"

    col_names = [_c_perc, _c_cumul, _c_count,
                 _s_perc, _s_cumul, _s_count,
                 _t_perc, _t_cumul, _t_count]
    
    def __init__(self, sl: pd.DataFrame, ml: pd.DataFrame):
        assert all(c1 == c2 for c1, c2 in zip(sl.columns, ml.columns))
        # `sl` and `ml` are both `pd.DataFrame` containing all samples as columns,
        # and every node as rows, including `unclassified`
        self.sl = sl
        self.ml = ml

    @abc.abstractmethod
    def write_sep(self, df: pd.DataFrame, out_file: str, sep: str = ",") -> None:
        """Write DataFrame to file.
        """
        df.to_csv(out_file, sep=sep, header=True)

    @staticmethod
    def remove_zeros(*dfs: pd.DataFrame):
        mask: pd.Series[bool] = None
        for df in dfs:
            non_zero = df.apply(lambda r: ((r > 0).any()) | (r.name == 'unclassified'), axis=1)
            mask = non_zero if mask is None else mask | non_zero
        for df in dfs:
            yield df.loc[mask, :]

    @abc.abstractmethod
    def accumulate(self, _df: pd.DataFrame) -> pd.DataFrame:
        """Accumulate counts within measurement units.
        """
        pass

    @staticmethod
    def get_cutoff_mask(_counts: pd.DataFrame, cpm: float) -> pd.DataFrame:
        min_counts = _counts.sum(axis=0) \
            .mul(cpm / 1e6) \
            .astype(_counts.values.dtype)
        # zero out counts below cutoff
        msk = _counts.lt(min_counts, axis=1)
        return msk
    
    @staticmethod
    def add_dfs(large: pd.DataFrame, small: pd.DataFrame) -> pd.DataFrame:
        assert all(small.columns.isin(large.columns))
        assert all(large.columns.isin(small.columns))
        assert all(small.index.isin(large.index))
        df = large.copy(deep=True)
        df.loc[small.index, small.columns] += small.values
        return df

    def summarise(self, output_base: str, cpm: float = 0.0, sep: str = ",") -> None:
        # accumulate counts
        sl_cumul = self.accumulate(self.sl)
        ml_cumul = self.accumulate(self.ml)
        # combine counts
        tot_count = self.add_dfs(self.sl, self.ml)
        tot_cumul = self.add_dfs(sl_cumul, ml_cumul)
        # apply cutoff
        min_counts = tot_count.sum(axis=0) \
            .mul(cpm / 1e6) \
            .astype(tot_cumul.values.dtype)
        msk = tot_cumul.lt(min_counts, axis=1)
        tot_cumul[msk] = 0
        # check structure of indices
        assert "unclassified" in self.sl.index
        assert self.sl.index.isin(sl_cumul.index).all()
        assert self.ml.index.isin(ml_cumul.index).all()
        assert ml_cumul.index.isin(sl_cumul.index).all()
        # write sample summary files
        for sample_name in self.sl.columns:
            df = self.summarise_sample(self.sl[sample_name], sl_cumul[sample_name],
                                       self.ml[sample_name], ml_cumul[sample_name],
                                       tot_count[sample_name], tot_cumul[sample_name],
                                       cpm=cpm)
            # save to disk
            sample_file = os.path.join(output_base, sample_name + ".csv")
            self.write_sep(df, sample_file, sep=sep)
            print("[formatter] wrote sample summary %s to %s" % (sample_name, sample_file))
        # write total cumulative to disk
        cumul_file = os.path.join(output_base, "cumulative.cpm_%.3f.csv" % cpm)
        self.write_sep(tot_cumul, cumul_file, sep=sep)
        print("[formatter] wrote cumulative counts with %.3f cutoff to %s" % (cpm, cumul_file))

    def summarise_sample(self,
                         sl: pd.Series, sl_cumul: pd.Series,
                         ml: pd.Series, ml_cumul: pd.Series,
                         tot: pd.Series, tot_cumul: pd.Series,
                         cpm: float = 0.0) -> pd.DataFrame:
        num_reads = sl.sum() + ml.sum()
        # prepare dataframe
        df = pd.DataFrame(0, index=sl_cumul.index, columns=self.col_names)
        # insert counts
        df.loc[sl.index, self._c_count] = sl
        df.loc[sl_cumul.index, self._c_cumul] = sl_cumul
        df.loc[ml.index, self._s_count] = ml
        df.loc[ml_cumul.index, self._s_cumul] = ml_cumul
        df.loc[tot.index, self._t_count] = tot
        df.loc[tot_cumul.index, self._t_cumul] = tot_cumul
        # employ cutoff
        df = df.loc[df[self._t_cumul] > 0, :]
        # format percentage signs
        df.loc[:, self._c_perc] = self.to_percentage(df.loc[:, self._c_cumul], num_reads)
        df.loc[:, self._s_perc] = self.to_percentage(df.loc[:, self._s_cumul], num_reads)
        df.loc[:, self._t_perc] = self.to_percentage(df.loc[:, self._t_cumul], num_reads)
        return df

    @staticmethod
    def to_percentage(_ser: pd.Series, total: int, decimals: int = 3) -> pd.Series:
        return _ser \
            .mul(100 / total) \
            .round(decimals) \
            .astype("string") \
            .add(" %")

# phylogenetic class =====================================================================

class PhylogeneticResults(_ClassificationOutput):
    def __init__(self, sl: pd.DataFrame, ml: pd.DataFrame, tree: Index):
        super(PhylogeneticResults, self).__init__(sl, ml)
        self.tree = tree

    @staticmethod
    def load_tree(tree: Index = None, phylogeny_path: str = None) -> Index:
        # load tree
        if tree is None and phylogeny_path is None:
            raise ValueError("Tree must be supplied!")
        elif tree is None:
            _, tree = Index.load_newick(phylogeny_path)
        return tree

    @staticmethod
    def make_index(tree: Index) -> List[str]:
        return ["unclassified" if i == 0 else node.name for i, node in enumerate(tree.pool)]
    
    @staticmethod
    def read_phy_output(csv_file: str, tree: Index, sep: str = ',') -> pd.DataFrame:
        ## read raw data
        #
        _data = pd.read_csv(csv_file, sep=sep, index_col=0, header=0)
        _data.index = [Index._format_node_name(v) for v in _data.index]
        ## target buffer
        #
        df = pd.DataFrame(0,
                          index=PhylogeneticResults.make_index(tree),
                          columns=_data.columns)
        df.loc[_data.index] = _data.values
        return df

    # read-wise parsing ---------------------------------------------------------------

    @classmethod
    def from_raw_read_output(cls, raw_read_base: str, num_procs: int, tree: Index = None, phylogeny_path: str = None):
        tree = PhylogeneticResults.load_tree(tree=tree, phylogeny_path=phylogeny_path)
        # read data
        output_files = ls(raw_read_base, ext=".csv")
        nodes = PhylogeneticResults.make_index(tree)
        ser_list = PhylogeneticResults.count_all_reads(output_files, nodes, num_procs)
        # sample-wise merge counts
        (_, sl_series_list, ml_series_list) = zip(*ser_list)
        sl = pd.concat(sl_series_list, axis=1)
        ml = pd.concat(ml_series_list, axis=1)
        return cls(sl, ml, tree)

    @staticmethod
    def count_all_reads(result_files: List[str], nodes: List[str], num_procs: int) -> List[Tuple[str, pd.Series, pd.Series]]:
        num_procs = min(num_procs, len(result_files))
        with Pool(processes=num_procs) as mp_pool:
            f = partial(PhylogeneticResults.count_sample_reads,
                        nodes=nodes)
            ser_list = mp_pool.map_async(f, result_files).get()
        return sorted(ser_list, key = lambda t: t[0])
    
    @staticmethod
    def count_sample_reads(file: str, nodes: List[str]) -> Tuple[str, pd.Series, pd.Series]:
        name = os.path.basename(file)[:-4]
        sl = pd.Series(0, index=nodes, name=name)
        ml = pd.Series(0, index=nodes, name=name)
        with open(file, 'r') as f:
            for line in f:
                data = line.split("\t")
                pid = int(data[2][1:])
                if data[0] == 'S':
                    ml.iloc[pid] += 1
                else:
                    sl.iloc[pid] += 1
        return name, sl, ml

    # preloaded -------------------------------------------------------------------------

    @classmethod
    def from_tables(cls, sl_file: str, ml_file: str, sep: str = ',',
                    tree: Index = None, phylogeny_path: str = None):
        ## load tree
        #
        tree = PhylogeneticResults.load_tree(tree=tree, phylogeny_path=phylogeny_path)
        ## data parsing
        #
        sl = PhylogeneticResults.read_phy_output(sl_file, tree, sep=sep)
        ml = PhylogeneticResults.read_phy_output(ml_file, tree, sep=sep)
        return cls(sl, ml, tree)

    # virtual methods ---------------------------------------------------------------

    def accumulate(self, _df: pd.DataFrame) -> pd.DataFrame:
        df = _df.copy()
        for node in self.tree.pool[1:][::-1]:
            df.loc[node.name, :] += df.loc[node.children, :].sum(axis=0)
        return df

    def write_sep(self, df: pd.DataFrame, out_file: str, sep: str = ","):
        df = df.copy()
        df.index = [
            'p' + pid
            if (pid != "unclassified") and (self.tree[pid].type == "Branch")
            else pid
            for pid in df.index
        ]
        _ClassificationOutput.write_sep(self, df, out_file, sep=sep)

    # helpers ---------------------------------------------------------------------------

    def to_csv(self, sl_file: str, ml_file: str, keep_zeros: bool = False):
        sl = self.sl.copy()
        ml = self.ml.copy()
        # remove null rows
        if not keep_zeros:
            (sl, ml) = self.remove_zeros(sl, ml)
        if "unclassified" in ml.index:
            ml.drop(index="unclassified", inplace=True)
        # save counts to disk
        self.write_sep(sl, sl_file)
        self.write_sep(ml, ml_file)

    # plotting ---------------------------------------------------------------------------

    def draw_results(self,
                     output_base: str,
                     groups: List[Tuple[str, Tuple[str, ...]]] = None,
                     cpm: float = 0.0,
                     colour_list: List[str] = None, keep_zeros: bool = True,
                     use_node_names: bool = True, log_scores: bool = True,
                     itol_mode: bool = False, flat_colour: bool = False):
        ## prepare dataframes
        #
        #   - SL
        sl_cumul = self.accumulate(self.sl)
        #   - ML
        ml_cumul = self.accumulate(self.ml)
        #   - total
        tot_cumul = sl_cumul.copy()
        tot_cumul.loc[ml_cumul.index, :] += ml_cumul.values

        ## render
        #
        # classified results
        phy_pdf = os.path.join(output_base, "phylotree_classified.pdf")
        self.draw_counts(sl_cumul, phy_pdf, inplace=True, groups=groups,
                         cpm=cpm, colour_list=colour_list,
                         keep_zeros=keep_zeros, use_node_names=use_node_names,
                         log_scores=log_scores, itol_mode=itol_mode, flat_colour=flat_colour)
        # split results
        split_pdf = os.path.join(output_base, "phylotree_splits.pdf")
        self.draw_counts(ml_cumul, split_pdf, inplace=True, groups=groups,
                         cpm=cpm, colour_list=colour_list,
                         keep_zeros=keep_zeros, use_node_names=use_node_names,
                         log_scores=log_scores, itol_mode=itol_mode, flat_colour=flat_colour)
        # total counts
        total_pdf = os.path.join(output_base, "phylotree_total.pdf")
        self.draw_counts(tot_cumul, total_pdf, inplace=True, groups=groups,
                         cpm=cpm, colour_list=colour_list,
                         keep_zeros=keep_zeros, use_node_names=use_node_names,
                         log_scores=log_scores, itol_mode=itol_mode, flat_colour=flat_colour)

    def draw_counts(self,
                    counts: pd.DataFrame, out_file: str, inplace: bool = True,
                    groups: List[Tuple[str, Tuple[str, ...]]] = None,
                    cpm: float = 0.0,
                    colour_list: List[str] = None, keep_zeros: bool = True,
                    use_node_names: bool = True, log_scores: bool = True,
                    itol_mode: bool = False, flat_colour: bool = False):
        if not inplace:
            counts = counts.copy()
        ## cutoff
        #
        min_count_per_sample = counts.sum(axis=0) \
            .mul(cpm / 1e6) \
            .astype(counts.values.dtype)

        # zero out counts below cutoff
        msk = counts.lt(min_count_per_sample, axis=1)
        counts[msk] = 0

        ## render
        #
        clades_with_counts = (counts > 0).any(axis=1).sum()
        if clades_with_counts > 0:
            self.tree.draw_counts(counts, out_file, groups=groups, colour_list=colour_list,
                                  keep_zeros=keep_zeros, use_node_names=use_node_names,
                                  log_scores=log_scores, itol_mode=itol_mode,
                                  flat_colour=flat_colour)
        else:
            print("Skipping plotting of %s - no samples with counts in this matrix." % out_file)

# taxonomic class =====================================================================

class TaxonomicResults(_ClassificationOutput):
    _rank = "Rank"
    _scientific_name = "Scientific Name"

    def __init__(self, sl: pd.DataFrame, ml: pd.DataFrame, taxonomy: TaxonomyInterface):
        super(TaxonomicResults, self).__init__(sl, ml)
        self.taxonomy = taxonomy
        required_taxids = {
            ptid
            for taxid in self.ml.index
            if taxid != "unclassified"
            for ptid in self.taxonomy.get_lineage(taxid)
        }
        self.taxid2rank = self.taxonomy.get_ranks(*required_taxids)
        self.taxid2name = self.taxonomy.get_names(*required_taxids)

    # hand-over -------------------------------------------------------------------------

    @classmethod
    def from_phylogenetic(cls, _data: PhylogeneticResults, taxonomy: TaxonomyInterface) :
        # convert to taxonomic counts
        t_cx_df = TaxonomicResults.gather_by_taxid(_data.sl.copy(deep=True), taxonomy)
        t_spx_df = TaxonomicResults.gather_by_taxid(_data.ml.copy(deep=True), taxonomy)
        return cls(t_cx_df, t_spx_df, taxonomy)

    # preloaded -------------------------------------------------------------------------

    @classmethod
    def from_tables(cls, sl_file: str, ml_file: str, tree: Index, taxonomy: TaxonomyInterface, sep: str = ','):
        # load phylogenetic output
        cx_df = PhylogeneticResults.read_phy_output(sl_file, tree, sep=sep)
        spx_df = PhylogeneticResults.read_phy_output(ml_file, tree, sep=sep)
        # convert to taxonomic counts
        t_cx_df = TaxonomicResults.gather_by_taxid(cx_df, taxonomy)
        t_spx_df = TaxonomicResults.gather_by_taxid(spx_df, taxonomy)
        return cls(t_cx_df, t_spx_df, taxonomy)

    @staticmethod
    def gather_by_taxid(df: pd.DataFrame, taxonomy: TaxonomyInterface) -> pd.DataFrame:
        # combine counts within taxid
        df["taxid"] = df.index \
            .to_series() \
            .map(taxonomy.clade2taxid)
        df = df \
            .groupby(["taxid"], axis=0) \
            .sum()
        assert df.index.name == "taxid"
        # sort index
        _order = {taxid: i+1 for i, taxid in enumerate(taxonomy.taxid_pool)}
        _order["unclassified"] = 0
        df.sort_index(inplace=True, key=lambda _i: _i.map(_order))
        return df
    
    # virtual methods ------------------------------------------------------------------

    def accumulate(self, _df: pd.DataFrame) -> pd.DataFrame:
        ## ensure lineage completeness
        #
        current_taxids = set(_df.index.values)
        current_taxids.remove("unclassified")
        required_taxids = {
            ptid
            for taxid in current_taxids
            for ptid in self.taxonomy.get_lineage(taxid)
        }
        # add root
        if "1" not in current_taxids:
            required_taxids.add("1")
        # append missing taxids
        missing_taxids = list(required_taxids - current_taxids)
        sec = pd.DataFrame(0, index=missing_taxids, columns=_df.columns)
        df = pd.concat([_df, sec])
        ## accumulate counts
        #
        for tid in self.taxonomy.taxid_pool[::-1]:
            if tid in df.index:
                for cid in self.taxonomy.taxid2child_taxids[tid]:
                    if cid in df.index:
                        df.loc[tid, :] += df.loc[cid, :]
        ## sort result
        #
        _order = {taxid: i+1 for i, taxid in enumerate(self.taxonomy.taxid_pool)}
        _order["unclassified"] = 0
        df.sort_index(inplace=True, key=lambda _i: _i.map(_order))
        return df
    
    def attach_taxonomic_data(self, df: pd.DataFrame) -> pd.DataFrame:
        has_taxonomy = df.index.isin(self.taxid2rank.keys()) & df.index.isin(self.taxid2name.keys())
        df.loc[has_taxonomy, self._rank] = df.index[has_taxonomy] \
            .to_series() \
            .map(self.taxid2rank)
        df.loc[has_taxonomy, self._scientific_name] = df.index[has_taxonomy] \
            .to_series() \
            .map(self.taxid2name)
        return df

    def write_sep(self, df: pd.DataFrame, out_file: str, sep: str = ","):
        df = self.attach_taxonomic_data(df)
        _ClassificationOutput.write_sep(self, df, out_file, sep=sep)

    def convert_raw(self, input_base: str, output_base: str, tree: Index):
        for read_file in ls(input_base, ext=".csv"):
            out_file = os.path.join(output_base, os.path.basename(read_file))
            with open(read_file, 'r') as infh, open(out_file, 'w') as outfh:
                for line in infh:
                    read_data = line \
                        .strip() \
                        .split("\t")[:-1]
                    if read_data[0] == 'U':
                        read_data[2] = "NA"
                    else:
                        pid = int(read_data[2][1:])
                        clade_name = tree.pool[pid].name
                        read_data[2] = self.taxonomy.clade2taxid[clade_name]
                    outfh.write("\t".join(read_data) + "\n")
            print("[taxonomy] mapped read-wise output from %s to %s" % (read_file, out_file))
