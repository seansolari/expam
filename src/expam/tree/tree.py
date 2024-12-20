import random
from math import floor, log
import os
import re
import sys
import traceback
from typing import List, Tuple
import pandas as pd
from expam.tree.location import Location


def propose_lca(c1, c2):
    """
    Given two binary coordinates, compute the shortest common prefix.

    :param coord: Indexable.
    :return: List.
    """
    m = min(len(c1), len(c2))
    if m == 0:
        return []

    for i in range(1, m + 1):
        if c1[-i] != c2[-i]:
            if i == 1:
                return []
            else:
                return c1[-i + 1:]
    else:
        return c1[-m:]


class Index:
    """ Phylogeny index that can load, save and manipulate Newick trees.
    """
    def __init__(self):
        self.db = None
        self.leaf_taxid = None
        self.taxid_lineage = None

        self._pointers = {}  # Pointers is a map from node name to pool location.
        self.pool = [Location()]  # Empty location to represent an 'unclassified' location.

    def __len__(self):
        return len(self.pool)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'<Phylogeny Index, length={len(self)}>'
        
    def use_db(self, db_path: str):
        if not os.path.exists(db_path):
            raise OSError("Can't find database path %s." % db_path)

        phy_path = os.path.join(db_path, 'phylogeny')
        accession_ids_path = os.path.join(phy_path, 'accession_ids.csv')
        taxid_lineage_path = os.path.join(phy_path, 'taxid_lineage.csv')

        for path in (phy_path, accession_ids_path, taxid_lineage_path):
            if not os.path.exists(path):
                raise OSError("Can't find required data %s. Have you run download_taxonomy?")

        self.load_leaf_taxid(accession_ids_path)
        self.load_taxid_lineage(taxid_lineage_path)

    def load_leaf_taxid(self, accession_ids_path: str):
        if self.leaf_taxid is None:
            self.leaf_taxid = {}
        else:
            return

        with open(accession_ids_path, 'r') as f:
            for line in f:
                name, _, taxid = line.strip().split(',')
                self.leaf_taxid[name] = taxid

    def load_taxid_lineage(self, taxid_lineage_path: str):
        if self.taxid_lineage is None:
            self.taxid_lineage = {}
        else:
            return

        with open(taxid_lineage_path, 'r') as f:
            for line in f:
                data = line.strip().split(',')
                self.taxid_lineage[data[0]] = tuple(data[1:])

    def get_node_lineages(self, node):
        taxids = [self.leaf_taxid[leaf] for leaf in self.get_child_leaves(node)]
        return [self.taxid_lineage[taxa] for taxa in taxids]

    @classmethod
    def load_newick(cls, path, keep_names=False, verbose=True):
        """load_newick Load Newick tree from file.

        :param path: path to Newick file
        :type path: str
        :raises OSError: file does not exist
        :return: name of leaves and phylogeny Index object
        :rtype: List[str], expam.tree.Index
        """
        newick_str = ""

        if not os.path.exists(path):
            raise OSError("Couldn't find path %s." % path)

        with open(path, "r") as f:
            for line in f:
                newick_str += line.strip()

        return cls.from_newick(newick_str, keep_names=keep_names, verbose=verbose)

    @classmethod
    def from_newick(cls, newick_string, keep_names=False, verbose=True):
        """from_newick Parse Newick string.

        :param newick_string: Newick string encoding tree.
        :type newick_string: str
        :return: name of leaves and phylogeny Index object
        :rtype: List[str], expam.tree.Index
        """
        stream = sys.stdout if verbose else open(os.devnull, 'w')

        # Remove whitespace.
        newick_string = newick_string.replace(" ", "").replace("\n", "")

        print("* Initialising node pool...", file=stream)
        index_pool = cls.init_pool(newick_string)

        print("* Checking for polytomies...", file=stream)
        cls.resolve_polytomies(index_pool)

        print("* Finalising index...", file=stream)
        leaf_names, index = cls.from_pool(index_pool, keep_names=keep_names)

        if stream != sys.stdout:
            stream.close()

        return leaf_names, index

    @staticmethod
    def init_pool(newick):
        pool = [Location()]

        i = 0
        parent_stack = []
        currentNode = None

        #
        # Create a new node.
        #
        def new_node(name, type):
            nonlocal currentNode

            currentNode = Location(name, type)

            pool.append(currentNode)

            # Top node has no ancestor.
            if len(parent_stack) > 0:
                parent_stack[-1].children.append(currentNode)

        #
        # Get node name.
        #
        def close_branch():
            nonlocal parent_stack, currentNode, i, newick, NEWICK_PARSE

            current_branch = parent_stack.pop()

            # Get this branch's name, if it is given.
            if newick[i + 1] not in NEWICK_PARSE:
                i += 1
                branch_name = parse_string()

                current_branch.name = branch_name

            currentNode = current_branch

        #
        # New branch in phylogeny and parent stack.
        #
        def new_branch():
            nonlocal currentNode

            new_node(name="", type="Branch")
            parent_stack.append(currentNode)

        #
        # New Leaf in phylogeny.
        #
        def new_leaf():
            nonlocal currentNode

            new_node(name="", type="Leaf")

        #
        # Parse name.
        #
        def parse_string(force_digits=False):
            nonlocal i  # Inherit current index.

            # Detect quotations.
            if newick[i] == "'":
                terminal = "'"

                i += 1

            else:
                terminal = NEWICK_PARSE

            j = i + 1
            while newick[j] not in terminal:
                j += 1

            string = newick[i:j]

            # Check formatting.
            if force_digits:
                try:
                    string = float(string)

                except ValueError:
                    raise ValueError("Invalid distance declaration %s!" % newick[i:j])

            if terminal == "'":
                i = j
            else:
                i = j - 1

            return string

        #
        # Clear current node (stack-ish).
        #
        def clear_node():
            nonlocal currentNode

            currentNode = None

        #
        # Set distance of current node.
        #
        def set_node_distance():
            nonlocal currentNode, i

            i += 1  # Skip ':' character.
            currentNode.distance = parse_string(force_digits=True)

        #
        # Do nothing.
        #
        def finisher():
            nonlocal i

            i = len(newick)

        # Parsing functions for Newick format.
        NEWICK_PARSE = {
            "(": new_branch,
            ")": close_branch,
            ",": clear_node,
            ":": set_node_distance,
            ";": finisher
        }

        # Parsing loop.
        while i < len(newick):
            if newick[i] in NEWICK_PARSE:
                NEWICK_PARSE[newick[i]]()

            else:
                new_leaf()

                currentNode.name = parse_string()

            i += 1

        # Return pool array for insertion into expam Index.
        return pool

    @staticmethod
    def resolve_polytomies(pool):
        """
        If the phylogeny contains polytomies, continually join the first two
        children with parents of distance 0 until the polytomy is resolved.

        :param pool: List.
        :return: None
        :rtype: None
        """
        # Processing loop.
        i = 1
        while i < len(pool):
            node = pool[i]

            # Continually join first two nodes.
            while len(node.children) > 2:
                print("\tPolytomy (degree=%d) detected! Resolving..."
                      % len(node.children))

                # Declare new parent of distance 0.
                newParent = Location(
                    name="",
                    type="Branch",
                    dist=0.0
                )

                # Replace children.
                for _ in range(2):
                    newParent.children.append(node.children.pop(0))

                # Insert into pool.
                node.children.append(newParent)
                pool.insert(i + 1, newParent)

            i += 1

    @classmethod
    def from_pool(cls, pool, keep_names=False, map_children_objects=True, update_nchildren=True):
        index = cls()
        index.pool = pool
        leaf_names = index.update_via_pool(keep_names=keep_names, map_children_objects=map_children_objects, update_nchildren=update_nchildren)

        return leaf_names, index

    def update_via_pool(self, keep_names=False, map_children_objects=True, update_nchildren=True):
        self._pointers = {}
        leaf_names = []

        branch_id = 1
        for i, node in enumerate(self.pool):
            if i == 0:
                continue

            if node.type == "Branch":
                if not keep_names:
                    node.name = str(branch_id)
                branch_id += 1
            else:
                # Store leaf names.
                leaf_names.append(node.name)

            # Set name mapping.
            self._pointers[node.name] = i

        if map_children_objects:
            # Replace children objects with respective names.
            # Note: This requires the name mapping to be complete...
            for node in self.pool[1:]:
                if node.type == "Branch":
                    node.children = [child.name for child in node.children]

        # Set (binary) coordinates.
        self.update_coordinates()

        if update_nchildren:
            # Update nchildren.
            self.update_nchildren()

        return leaf_names

    def append(self, item):
        # Declare an index.
        index = len(self.pool)

        # Set pointer.
        self._pointers[item.name] = index

        self.pool.append(item)

    def join(self, a, b, glue, dist=None):
        if dist is None:
            dist = [1.0, 1.0]
        i, j = self._pointers[a], self._pointers[b]
        i, j = self.smaller(i, j)  # i < j.

        one, two = j, self.neighbouring(j)
        three = self.neighbouring(i)
        self.move([one, two], three)

        self.insert(i, glue)

        # Update number of children and distances.
        glue.children.extend((a, b))

        for i, child in enumerate([a, b]):
            glue.nchildren += 1 + ((self[child].type == "Branch") * self[child].nchildren)

            self[child].distance = dist[i]

    def move(self, i, j):
        # NOTE: i > j.
        width = i[1] - i[0]

        if i[0] == j:
            return

        # Delete affected pointers.
        for n in range(j, i[1]):
            del self._pointers[n]

        # Retrieve items to be moved.
        temp = []
        for n in range(width):
            temp.append(self.pool.pop(i[1] - (n + 1)))

        # Insert these items in at j. Since they were popped in reverse order,
        # insert append them in the same index to revert back to the original
        # order.
        for obj in temp:
            self.pool.insert(j, obj)

        # Update pointers.
        for n in range(j, i[1]):
            self._pointers[n] = self.pool[n].name

    def insert(self, i, item):
        # Update pointers.
        for n in range(i, len(self.pool)):
            loc = self.pool[n]

            self._pointers[loc.name] = n + 1

        self.pool.insert(i, item)
        self._pointers[item.name] = i

    @staticmethod
    def _format_node_name(node_name):
        node_name = str(node_name)
        return node_name[1:] if (node_name[1:].isdigit() and node_name.startswith('p')) else node_name

    def give_branch_name(self, node_name):
        node = self[node_name]
        if node.type == "Branch" and node.name[0] != 'p':
            return 'p' + node_name
        else:
            return node_name

    def __setitem__(self, key, value):
        try:
            pool_index = self._pointers[key]
        except KeyError:
            pool_index = self._pointers[self._format_node_name(key)]
        self.pool[pool_index] = value

    def __getitem__(self, key):
        try:
            pool_index = self._pointers[key]
        except KeyError:
            pool_index = self._pointers[self._format_node_name(key)]
        return self.pool[pool_index]

    def __delitem__(self, key):
        self.pool.pop(self._pointers[str(key)])

    def neighbouring(self, i):
        return i + 1 + self.pool[i].nchildren

    def update_coordinates(self):
        coord, i = [-1], 1

        while i < len(self.pool):
            while coord[0] + 1 > 1:
                coord.pop(0)

            coord[0] += 1

            if self.pool[i].type == "Branch":
                self.pool[i].coordinate = list(coord[:-1])

                coord.insert(0, -1)
            else:  # It's a leaf.
                self.pool[i].coordinate = list(coord[:-1])

            i += 1

    def update_nchildren(self):
        for node in self.pool[1:][::-1]:
            node.nchildren = len(node.children) + sum([
                self[child].nchildren
                for child in node.children
            ])

    def coord(self, coordinate):
        """coord Return Location (node) at coordinate.

        :param coordinate: binary list representing path to node
        :type coordinate: list
        :return: node in tree
        :rtype: expam.tree.Location
        """
        pool_index = 1

        for step in coordinate[::-1]:
            pool_index += 1
            pool_index += int(step) * (self.pool[pool_index].nchildren + 1)

        return self.pool[pool_index]

    def print_index(self):
        for l in self.pool:
            print(l.name, l.coordinate, l.nchildren)

    @staticmethod
    def smaller(a, b):
        if a < b:
            return a, b
        return b, a

    def to_newick(self):
        """to_newick Output tree to Newick format.

        :return: Newick format tree
        :rtype: str
        """
        if len(self.pool) == 2:
            return self.pool[1].name + ";"

        NODE_FORMAT = "$node$"

        newick = "(%s,%s)%s;" % (NODE_FORMAT, NODE_FORMAT, self.pool[1].name)
        node_formats = {
            "Branch": "(%s,%s){name}:{distance}" % (NODE_FORMAT, NODE_FORMAT),
            "Leaf": "{name}:{distance}"
        }

        for node in self.pool[2:]:
            node_data = node_formats[node.type].format(name=node.name, distance=node.distance)
            newick = newick.replace(NODE_FORMAT, node_data, 1)

        return newick

    def yield_child_nodes(self, node_name):
        """yield_child_nodes Yields node and children nodes (both branches and leaves).

        :param node_name: name of node to start yielding from
        :type node_name: str
        :yield: node names at or below node_name
        :rtype: str
        """
        node_name = self._format_node_name(node_name)
        active_nodes = [node_name]

        while active_nodes:
            node = active_nodes.pop(0)

            yield node
            active_nodes.extend(self[node].children)

    def yield_leaves(self, node_name):
        """yield_leaves Yield only the leaves at or below some node.

        :param node_name: node to retrieve leaves from.
        :type node_name: str
        :yield: leaf names at or below node_name.
        :rtype: str
        """
        node_name = self._format_node_name(node_name)

        if self[node_name].type == "Leaf":
            yield node_name
            return

        branches = [node_name]

        while branches:
            branch = branches.pop(0)

            for child in self[branch].children:
                if self[child].type == "Leaf":
                    yield child
                else:
                    branches.append(child)

    def get_child_nodes(self, node_name):
        """get_child_nodes Return list of nodes at or below node_name.

        :param node_name: name of node 
        :type node_name: str
        :return: list of node names
        :rtype: List[str]
        """
        return list(self.yield_child_nodes(node_name))

    def get_child_leaves(self, node_name):
        """get_child_leaves Get list of leaves at or below node_name.

        :param node_name: name of node
        :type node_name: str
        :return: list of leaf names
        :rtype: List[str]
        """
        return list(self.yield_leaves(node_name))

    def lca(self, *names):
        """lca Return name of the lowest common ancestor of these two nodes.

        Note that coordinates are read from right-to-left.

        :param name_one: name of node
        :type name_one: str
        :param name_two: name of node
        :type name_two: str
        """
        lca_coord = self.right_intersect(*(self[name].coordinate for name in names))
        return self.coord(lca_coord).name

    def reduce(self, nodes):
        """reduce Reduce tree until all that remains are these nodes.
        
        :param nodes: Names of nodes to keep
        :type nodes: Iterable[str]
        """
        _requested_node_set = {self._format_node_name(node) for node in nodes}
        child_to_idx = {
            child: self._pointers[self._format_node_name(child)]
            for child in {
                _child
                for _node in self.pool[1:]
                for _child in _node.children
            }
        }
        def pop_nodes(*childs):
            nonlocal child_to_idx
            _objs = {}
            for child in childs:
                idx = child_to_idx[child]
                _objs[child] = self.pool.pop(idx)
                self.pool.insert(idx, None)
            return _objs

        node_to_parent = {}
        for node in self.pool[1:]:
            node_name = self._format_node_name(node.name)
            for child in node.children:
                if child not in _requested_node_set:
                    node_to_parent[child] = node_name
        def replace_parents_child(node, child):
            nonlocal node_to_parent
            try:
                parent_children_list = self[node_to_parent[node]].children
                parent_children_list[parent_children_list.index(node)] = child
            except KeyError:    # Root.
                pass

        # Reset node pool to reflect changes in tree structure.
        branch_indices = [i for i, node in enumerate(self.pool) if node.type == "Branch"]
        for branch_idx in branch_indices[::-1]:
            branch = self.pool[branch_idx]
            branch_name = self._format_node_name(branch.name)
            children = branch.children
            children_with_counts = [child for child in children if child in _requested_node_set]

            if not children_with_counts:
                # Make this branch into a leaf node.
                branch.type = "Leaf"
                pop_nodes(*children)
                branch.children.clear()
            elif len(children_with_counts) == 2:
                # Pass validity of this node - dichotomous tree structure depends on
                # this node being present.
                _requested_node_set.add(branch_name)
            elif branch_name not in _requested_node_set:
                # Only one child with counts is implied by the previous two cases.
                # Replace this node with the valid child node.
                valid_child_name, branch_dist = children_with_counts[0], branch.distance
                # Remove children from active node pool.
                removed_children = pop_nodes(*children)
                # Remove branch from node pool.
                self.pool.pop(branch_idx)
                # Insert replacement.
                removed_children[valid_child_name].distance += branch_dist
                self.pool.insert(branch_idx, removed_children[valid_child_name])
                child_to_idx[valid_child_name] = branch_idx
                # Set this node as replacement child of parent.
                replace_parents_child(branch_name, valid_child_name)
                
        # Reset data structures with the new index.
        self.pool = [node for node in self.pool if node is not None]
        self.update_via_pool(keep_names=True, map_children_objects=False)

    @staticmethod
    def right_intersect(*lists) -> list:
        lists = [l for l in lists if len(l) > 0]
        L = min(len(l) for l in lists)
        for i in range(L):
            rev_index = i + 1
            ref = lists[0][-rev_index]
            for other in lists[1:]:
                if other[-rev_index] != ref:
                    if i == 0:
                        return []
                    else:
                        return lists[0][-i:]
        else:
            return lists[0][-L:]

    def draw_counts(self,
                    counts: pd.DataFrame, out_file: str,
                    groups: List[Tuple[str, Tuple[str, ...]]] = None,
                    colour_list: List[str] = None, keep_zeros: bool = True,
                    use_node_names: bool = True, log_scores: bool = True,
                    itol_mode: bool = False, flat_colour: bool = False):
        ## specified sample groupings
        #
        from expam.sequences import format_name

        if groups is None:
            groups = [(None, (section,)) for section in counts.columns]
        else:
            # Format group names.
            groups = [
                (
                    col,
                    tuple(
                        format_name(group_member, remove_comp=True)
                        for group_member in group
                    )
                )
                for col, group in groups
            ]

        # Combine counts within the group.
        for _, group in groups:
            if len(group) > 1:
                group_name = group[0]
                counts.loc[:, group_name] = counts[list(group)].sum(axis=1)
                counts.drop(labels=list(group[1:]), axis=1, inplace=True)

        # Remove any groups that weren't specified.
        all_groups = counts.columns.tolist()
        specified_groups = set(group[0] for _, group in groups)
        unspecified_groups = [group for group in all_groups if group not in specified_groups]
        counts.drop(labels=unspecified_groups, axis=1, inplace=True)

        ## prune tree
        #
        if not keep_zeros:
            # reduce counts dataframe
            msk = (counts >= 0).any(axis=1)
            counts = counts.loc[msk, :]
            # reduce tree
            nodes_with_counts = counts.index.tolist()
            self.reduce(nodes_with_counts)

        ## configure colour generator from count scores
        #
        max_vector = tuple(counts.max())
        min_vector = tuple(int(v) for v in counts[counts > 0].min())
        colours = [colour for colour, _ in groups]
        colour_generator = ColourGenerator(
            colours,
            max_vector,
            min_vector,
            colour_list=colour_list,
            log_scores=log_scores
        )
 
        if itol_mode:
            ## plot using iTOL
            #
            itol_data = [
                "TREE_COLORS",
                "# Drag me onto the iTOL tree (in browser) to apply me!",
                "SEPARATOR TAB",
                "DATA",
            ]

            for node in counts.index:
                node_counts = counts.loc[node, :]

                if node_counts.any():
                    colour = colour_generator.generate(tuple(node_counts), flat_colour=flat_colour)
                    itol_data.append("%s\trange\t%s\t%s" % (node, colour, node.upper()))

            # Write itol_data and Newick tree to file.
            name_modifier = re.findall(r"phylotree_(\S+).pdf", os.path.basename(out_file))[0]
            itol_dir = os.path.join(os.path.dirname(out_file), 'itol_%s' % name_modifier)
            if not os.path.exists(itol_dir):
                os.mkdir(itol_dir)
                print("iTOL output directory created at %s." % itol_dir)

            # Write style file.
            itol_style_dir = os.path.join(itol_dir, 'style.txt')
            with open(itol_style_dir, 'w') as f:
                f.write('\n'.join(itol_data))

            print("iTOL style written to %s." % itol_style_dir)

            # Write Newick file.
            itol_tree_dir = os.path.join(itol_dir, 'itol_tree.nwk')
            with open(itol_tree_dir, 'w') as f:
                f.write(self.to_newick())

            print("iTOL tree written to %s. Put this in itol (in your browser)." % itol_tree_dir)

        else:
            ## plot with ete3
            #
            try:
                from ete3 import AttrFace, faces, Tree, TreeStyle, NodeStyle
            except ModuleNotFoundError as e:
                print("The ETE3 package is not installed. Either install this module, or use the\n\t--itol flag if you wish to use our native iTOL integration.")
                print("Skipping plotting of %s..." % out_file)
                return
            except ImportError as e:
                print("Could not import drawing attributes from ete3. Error raised:")
                print(traceback.format_exc())
                print("See FAQ section on GitHub for fix to this problem.\n\t(https://github.com/seansolari/expam#problems-during-installation)")
                print("Skipping plotting of %s..." % out_file)
                return

            # phylotree obj
            newick_string = self.to_newick()
            tree = Tree(newick_string, format=1)

            ## renderer
            #
            def layout(node):
                if node.is_leaf():
                    if use_node_names:
                        faces.add_face_to_node(AttrFace('name', fsize=20), node, column=0, position='aligned')

                if node.name in counts.index:
                    node_counts = counts.loc[node.name, :]

                    if node_counts.any():
                        colour = colour_generator.generate(tuple(node_counts), flat_colour=flat_colour)

                        ns = NodeStyle()
                        ns['bgcolor'] = colour
                        node.set_style(ns)

            ## execute rendering
            #
            ts = TreeStyle()
            ts.mode = "c"
            ts.show_leaf_name = False
            ts.layout_fn = layout
            ts.force_topology = False
            ts.allow_face_overlap = False
            ts.draw_guiding_lines = True
            ts.root_opening_factor = 1
            tree.render(
                out_file,
                tree_style=ts,
                dpi=300
            )
            print("Phylogenetic tree written to %s!" % out_file)

class HexColour:
    def __init__(self, hex_string):
        self.colour = hex_string.lstrip("#")

    def __add__(self, other):
        rgb_self = self.to_rgb()
        rgb_other = other.to_rgb()

        rgb_mix = rgb_self + rgb_other
        return rgb_mix.to_hex()

    def __iadd__(self, other):
        mix = self + other
        self.colour = mix.colour

    def __str__(self):
        return "#%s" % self.colour

    def to_rgb(self):
        hx = self.colour
        assert len(hx) == 6

        rgb = tuple(int(hx[2 * i:2 * (i + 1)], 16) for i in range(3))
        return RGBColour(rgb)

    def opaque(self, q):
        if q == 1:
            return self

        rgb = self.to_rgb().opaque(q)

        return rgb.to_hex()

class RGBColour:
    def __init__(self, tup):
        self.colour = tup

    def __add__(self, other):
        assert len(self.colour) == len(other.colour) == 3

        rgb = tuple(
            floor(0.5 * (self.colour[i] + other.colour[i]))
            for i in range(3)
        )
        return RGBColour(rgb)

    def __iadd__(self, other):
        mix = self + other
        self.colour = mix.colour

    @staticmethod
    def _format_rgb(num):
        f_hx = hex(num).lstrip('0x').upper()

        if len(f_hx) < 2:
            f_hx = ('0' * (2 - len(f_hx))) + f_hx

        return f_hx

    def to_hex(self):
        rgb = self.colour
        assert len(rgb) == 3

        hx = "#" + "".join(self._format_rgb(factor) for factor in rgb)
        return HexColour(hx)

    def opaque(self, q):
        return RGBColour(tuple(
            floor(255 - (q * (255 - v)))
            for v in self.colour
        ))

class RandomColour:
    alpha =  "0123456789ABCDEF"

    def __init__(self, rng=None, seed=1):
        if rng is not None:
            if not isinstance(rng, random.Random):
                raise TypeError("rng must be instance of random.Random!")

            self.rng = rng
        else:
            if seed is not None:
                self.rng = random.Random(seed)
            else:
                self.rng = random.Random()

    def __next__(self):
        return "#" + "".join(self.rng.choices(self.alpha, k=6))

class ColourGenerator:
    def __init__(self, colours, max_vector, min_vector=None, colour_list=None, log_scores=False):
        """

        :param colours:
        :param max_vector:
        :param colour_list:
        """
        # Ensure we have enough colours for plotting.
        g = RandomColour()

        # Ensure we have enough colours for each group.
        required_colours = sum(colour is None for colour in colours)
        colours = [colour for colour in colours if colour is not None]

        if required_colours > 0:
            if colour_list is None:
                colour_list = [next(g) for _ in range(required_colours)]

            else:
                colours_to_generate = required_colours - len(colour_list)
                colour_list.extend(
                    [next(g) for _ in range(colours_to_generate)]
                )

            colours.extend(colour_list)

        self.colours = [HexColour(colour) for colour in colours]

        #
        self.max_vector = max_vector
        self.min_vector = min_vector
        self.make_log_score = log_scores

        if self.make_log_score and self.min_vector is None:
            raise AttributeError("Log scores requested but not min vector supplied!")

        # Calculate normalisation factors.
        self._log_norm_factor = [log(mx / mn) for mx, mn in zip(self.max_vector, self.min_vector)] if self.make_log_score else None

    def generate(self, score_vector, flat_colour = False):
        colour = HexColour("#FFFFFF")  # White background.
        scorer = self._linear_scores if not self.make_log_score else self._log_scores

        if not flat_colour:
            for next_colour, opacity in scorer(score_vector):
                colour = colour + next_colour.opaque(opacity)
        else:
            for next_colour, _ in scorer(score_vector):
                colour = colour + next_colour

        return str(colour)

    def _linear_scores(self, vector):
        for i, score in enumerate(vector):
            if score > 0:
                yield self.colours[i], score / self.max_vector[i]

    def _log_scores(self, vector):
        for i, score in enumerate(vector):
            if score > 0:
                yield self.colours[i], log(score / self.min_vector[i]) / self._log_norm_factor[i]
