# version: 2
from __future__ import annotations
import os
from typing import List, Mapping, Set, Tuple, Type, TypeVar
from ete3 import NCBITaxa

from expam.database import FileLocationConfig
from expam.database.config import JSONConfig
from expam.tree.tree import Index
from expam.utils import work_in_directory

# NCBI -----------------------------------------------------------------------

class ncbi_taxdump:
    _remote_dump_url = "http://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
    _sqlite_file_name = "ncbi.sqlite"
    _taxdump_file_name = "ncbi.dmp.tar.gz"

    # initialisation --------------------------------------------------------

    def __init__(self, taxonomy_dir: str, force_update: bool = False) -> None:
        self.taxdump = os.path.abspath(os.path.join(taxonomy_dir, self._taxdump_file_name))
        self.sqlite = os.path.abspath(os.path.join(taxonomy_dir, self._sqlite_file_name))
        self._check_sql_db(force_update)

    def log(self, msg: str) -> None:
        print("[taxdump] %s" % msg)
        return

    @staticmethod
    def contains_taxonomic_database(taxonomy_dir: str) -> bool:
        sqlite_file = os.path.abspath(os.path.join(taxonomy_dir,
                                                   ncbi_taxdump._sqlite_file_name))
        return os.path.exists(sqlite_file)

    def _check_sql_db(self, force_update: bool) -> None:
        if not os.path.exists(self.sqlite) or force_update:
            from ete3.ncbi_taxonomy.ncbiquery import update_db
            with work_in_directory(os.path.dirname(self.taxdump)):
                self._download_taxonomy_dump()
                update_db(self.sqlite, targz_file=self.taxdump)
                os.remove(self.taxdump)
        return

    def _download_taxonomy_dump(self) -> None:
        try:
            from urllib import urlretrieve
        except ImportError:
            from urllib.request import urlretrieve
        self.log('downloading taxdump.tar.gz from NCBI FTP site (via HTTP) to %s' % self.taxdump)
        urlretrieve(self._remote_dump_url, self.taxdump)
        self.log('done. Parsing...')
        return

    # load ----------------------------------------------------------------

    def load(self) -> NCBITaxa:
        return NCBITaxa(self.sqlite)

class ncbi:
    def __init__(self,
                 taxonomy_dir: str,
                 force_update: bool = False):
        # parse taxdump
        dmp = ncbi_taxdump(taxonomy_dir, force_update=force_update)
        self._sql_handle = dmp.load()

    def get_lineage(self, taxid: str | int) -> List[str]:
        return [str(v) for v in self._sql_handle.get_lineage(taxid)]

    def get_lineages(self, *taxids: str | int) -> Mapping[str, List[str]]:
        return {taxid: self.get_lineage(taxid) for taxid in taxids}

    def get_ranks(self, *taxids: str) -> Mapping[str, str]:
        return self._sql_handle.get_rank(taxids)

    def get_names(self, *taxids: str) -> Mapping[str, str]:
        return self._sql_handle.get_taxid_translator(taxids)

# mapping ------------------------------------------------------------------------

T = TypeVar('T')

class NodeMapper:
    def __init__(self,
                 name2taxid: Mapping[str, str],
                 taxid2lineage: Mapping[str, List[str]],
                 tree: Index = None,
                 db: FileLocationConfig = None,
                 phylogeny_path: str = None) -> None:
        self.name2taxid = name2taxid
        self.taxid2lineage = taxid2lineage
        self.tree = NodeMapper.load_tree(db=db, tree=tree, phylogeny_path=phylogeny_path)
        # auxilliary data strutures
        self.taxid_pool: List[str] = []
        self.taxid2child_taxids: Mapping[str, Set[str]] = {}

    @classmethod
    def from_ncbi(cls: Type[NodeMapper],
                  ncbi_handle: ncbi,
                  tree: Index = None,
                  db: FileLocationConfig = None,
                  phylogeny_path: str = None)-> NodeMapper:
        name2taxid = NodeMapper.load_genome_taxids(db.accn_id)
        taxid2lineage = ncbi_handle.get_lineages(*set(name2taxid.values()))
        return cls(name2taxid, taxid2lineage, tree=tree, db=db, phylogeny_path=phylogeny_path)

    @staticmethod
    def load_tree(db: FileLocationConfig = None, tree: Index = None, phylogeny_path: str = None):
        # load tree
        if db is None and tree is None and phylogeny_path is None:
            raise ValueError("Must supply a tree!")
        elif tree is None:
            if phylogeny_path is None:
                conf = JSONConfig(db.conf)
                phylogeny_path = conf.get_phylogeny_path()
            _, tree = Index.load_newick(phylogeny_path)
        return tree

    @staticmethod
    def load_genome_taxids(accession_ids_file: str) -> Mapping[str, str]:
        name2taxid = {}
        missing_taxids = 0
        with open(accession_ids_file, 'r') as f:
            for line in f:
                try:
                    genome, _, taxid = line.strip() \
                        .split(",")
                except KeyError:
                    continue
                if (len(taxid) == 0) or (taxid == "None"):
                    missing_taxids += 1
                else:
                    name2taxid[genome] = taxid
        if missing_taxids > 0:
            raise ValueError("%d genomes missing taxonomic assignments. \
                             Update %s accordingly." % (missing_taxids, accession_ids_file))
        return name2taxid

    def _init_taxid_aux(self):
        self.taxid_pool = []
        self.taxid2child_taxids = {}

    def map_all_nodes(self) -> Tuple[Mapping[str, str], List[str], Mapping[str, Set[str]]]:
        nodes = [node.name for node in self.tree.pool[1:]]
        return self.map_nodes(nodes)

    def map_nodes(self, phy_nodes) -> Tuple[Mapping[str, str], List[str], Mapping[str, Set[str]]]:
        # Map phylogeny index to taxonomy.
        clade2taxid: Mapping[str, str] = {}
        for clade in phy_nodes:
            taxid = self.map_phylogeny_node(clade)
            clade2taxid[clade] = taxid
        # shift container owners
        taxid_pool__export = self.taxid_pool
        taxid2child_taxids__export = self.taxid2child_taxids
        self._init_taxid_aux()
        return clade2taxid, taxid_pool__export, taxid2child_taxids__export
    
    def map_phylogeny_node(self, node_name: str) -> str:
        # get leaves within this clade and their taxids
        taxids = [self.name2taxid[leaf] for leaf in self.tree.yield_leaves(node_name)]
        lineages = [self.taxid2lineage[taxid] for taxid in taxids]
        # find intersection of lineages of these taxids
        common_lineage = self.tuple_intersect(*lineages)
        assert len(common_lineage) >= 1
        taxid = common_lineage[-1]
        # insert lineage into taxid pool
        pool_index = 0
        for i, parent_taxid in enumerate(common_lineage):
            # insert into taxid pool
            try:
                pool_index = self.taxid_pool.index(parent_taxid) + 1
            except ValueError:
                self.taxid_pool.insert(pool_index, parent_taxid)
                pool_index += 1
            # track immediate child nodes
            if i + 1 < len(common_lineage):
                child_taxid = common_lineage[i + 1]
                try:
                    self.taxid2child_taxids[parent_taxid].add(child_taxid)
                except KeyError:
                    self.taxid2child_taxids[parent_taxid] = {child_taxid}
            elif parent_taxid not in self.taxid2child_taxids:
                self.taxid2child_taxids[parent_taxid] = set()
        return taxid
    
    @staticmethod
    def tuple_intersect(*tuples: List[T]) -> Tuple[T, ...]:
        intersection = ()
        length = min(len(tup) for tup in tuples)
        for i in range(length):
            item = tuples[0][i]
            for j in range(1, len(tuples)):
                if tuples[j][i] != item:
                    return intersection
            else:
                intersection += (item,)
        return intersection

# Expam -----------------------------------------------------------------------

class TaxonomyInterface:
    def __init__(self,
                 conf: FileLocationConfig,
                 tree: Index = None,
                 phylogeny_path: str = None,
                 force_download: bool = False,
                 force_remap: bool = False):
        self._ncbi = ncbi(conf.phylogeny, force_update=force_download)
        if force_download \
            or force_remap \
            or any(not os.path.exists(getattr(conf, attr)) \
                   for attr in ('clade_tax', 'tax_pool', 'tax_childs')):
            # create and save map files
            op = NodeMapper.from_ncbi(self._ncbi, tree=tree, db=conf, phylogeny_path=phylogeny_path)
            self.clade2taxid, self.taxid_pool, self.taxid2child_taxids = op.map_all_nodes()
            self.clade2taxid["unclassified"] = "unclassified"
            self.export(conf, op.tree)
        else:
            # load map files
            self.clade2taxid = self.load__clade2taxid(conf)
            self.taxid_pool = self.load__taxid_pool(conf)
            self.taxid2child_taxids = self.load__taxid2child_taxids(conf)

    def export(self, conf: FileLocationConfig, tree: Index) -> None:
        # CLADE_TAX
        with open(conf.clade_tax, "w") as f:
            for clade, taxid in self.clade2taxid.items():
                f.write("%s,%s\n" % (clade, taxid))
        # TAX_POOL
        with open(conf.tax_pool, "w") as f:
            f.write("\n".join(self.taxid_pool))
        # TAX_CHILDS
        with open(conf.tax_childs, "w") as f:
            for taxid, children in self.taxid2child_taxids.items():
                f.write("%s,%s\n" % (taxid, ";".join(children)))
        # CLADE_TABLE
        self.write_clade_table(conf.clade_table, tree)

    @staticmethod
    def load__clade2taxid(conf: FileLocationConfig) -> Mapping[str, str]:
        clade2tax = {}
        with open(conf.clade_tax, "r") as f:
            for line in f:
                try:
                    clade, taxid = line.strip()\
                        .split(',')
                except ValueError:
                    continue
                clade2tax[clade] = taxid
        return clade2tax
    
    @staticmethod
    def load__taxid_pool(conf: FileLocationConfig) -> List[str]:
        with open(conf.tax_pool, "r") as f:
            return [line.strip() for line in f]
        
    @staticmethod
    def load__taxid2child_taxids(conf: FileLocationConfig) -> Mapping[str, Set[str]]:
        taxid2child_taxids = {}
        with open(conf.tax_childs, "r") as f:
            for line in f:
                try:
                    taxid, childs = line.strip()\
                        .split(',')
                except ValueError:
                    continue
                taxid2child_taxids[taxid] = set(childs.split(";"))
        return taxid2child_taxids

    def get_ranks(self, *taxids: str) -> Mapping[str, str]:
        result = self._ncbi.get_ranks(*taxids)
        return {str(k): v for k, v in result.items()}

    def get_names(self, *taxids: str) -> Mapping[str, str]:
        result = self._ncbi.get_names(*taxids)
        return {str(k): v for k, v in result.items()}

    def clade_lineage(self, clade: str) -> List[str]:
        return self.get_lineage(self.clade2taxid[clade])
        
    def get_lineage(self, taxid: str) -> List[str]:
        result = self._ncbi.get_lineage(taxid)
        return [str(x) for x in result]
    
    def write_clade_table(self, file: str, tree: Index):
        # get taxonomic data for all taxids
        all_taxids = set(self.clade2taxid.values())
        all_taxids.discard("unclassified")
        taxid2name = self.get_names(*all_taxids)
        taxid2rank = self.get_ranks(*all_taxids)
        # write to disk
        with open(file, "w") as f:
            f.write("clade,taxid,taxa,rank\n")
            f.write("unclassified,unclassified,unclassified,no rank\n")
            for node in tree.pool[1:]:
                node_name = tree.give_branch_name(node.name)
                taxid = self.clade2taxid[node.name]
                f.write("%s,%s,%s,%s\n" % (node_name, taxid, taxid2name[taxid], taxid2rank[taxid]))
