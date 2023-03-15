#!/usr/bin/env python3
import click
from expam.classify.taxonomy import TaxonomyNCBI
from expam.database import FileLocationConfig
from expam.database.config import JSONConfig, load_database_config
from expam.tree.tree import Index

class NodeMapper:
    def __init__(self, db_path: str) -> None:
        # load database configuration and reference phylogeny
        database_config: FileLocationConfig = load_database_config(db_path)
        config: JSONConfig = JSONConfig(database_config.conf)
        _, _, phylogeny_path, _, _ = config.get_build_params()
        _, self.phylogeny_index = Index.load_newick(phylogeny_path)

        # load taxonomic classifications
        try:
            tax_obj = TaxonomyNCBI(database_config)
            self.name_lineage, self.taxon_rank = tax_obj.load_taxonomy_map()
        except OSError:
            if self.convert_to_taxonomy:
                click.secho("first run `download_taxonomy` to collect associated taxonomy data", fg = "red", bold = True)

        # auxilliary data structures
        self.taxid_immediate_children = {"1": set()}        # Map from tax_id -> immediate children.
        self.tax_id_pool = ["1"]                            # Children must appear later than parent this list.

    def map_all_nodes(self):
        nodes = [node.name for node in self.phylogeny_index.pool[1:]]
        return self.map_nodes(nodes)

    def map_nodes(self, phy_nodes):
        # Map phylogeny index to taxonomy.
        phylogeny_to_taxonomy, tax_id_data = {}, {}
        for node in phy_nodes:
            tax_id, tax_id_map = self.map_phylogeny_node(node)

            phylogeny_to_taxonomy[node] = tax_id
            tax_id_data.update(tax_id_map)

        return phylogeny_to_taxonomy

    def map_phylogeny_node(self, node_name):
        leaves = self.get_leaves(node_name)
        lineages = [self.name_lineage[leaf] for leaf in leaves]
        common_lineage = self.tuple_intersect(*lineages)

        if not common_lineage:
            phy_tax_id, phy_rank = "1", "root"
        else:
            phy_tax_id, phy_rank = self.taxon_rank[common_lineage[-1]]

        taxid_rank_and_lineage = {}
        pool_index = 1
        for i in range(len(common_lineage)):
            # Store taxonomic data.
            parent_tax_id, parent_rank = self.taxon_rank[common_lineage[i]]
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
                child_tax_id, child_rank = self.taxon_rank[common_lineage[i + 1]]
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

@click.command()
@click.option("db_path", "-db", "--db", type=str, help="path to expam database", required=True)
@click.option("out_file", "-o", "--out", type=str, help="where to save processed results", required=True)
def classify_clades(db_path: str, out_file: str):
    mapper = NodeMapper(db_path)
    click.secho("mapping nodes", fg = "blue", bold = True)
    data = mapper.map_all_nodes()

    with open(out_file, "w") as f:
        for k, v in data.items():
            f.write("%s,%s\n" % (k, v))
    click.secho("results saved to %s" % out_file, fg = "green", bold = True)

if __name__ == "__main__":
    classify_clades()