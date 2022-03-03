import random
from math import floor
import json
import os
import traceback

import numpy as np
import pandas as pd

PHYLA_COLOURS = [("Proteobacteria", "#E85343"), ("Bacteroidetes", "#30B674"),
                 ("Firmicutes", "#6585C3"), ("Actinobacteria", "#FEDB45")]


def simulate_balanced_phylogeny(names):
    """
    :param names: List - list of leaf names.
    """

    # Declare an index.
    index = Index()

    # Insert all children.
    for name in names:
        # Make child node.
        childNode = Location(
            name=name,
            type="Leaf",
            dist=1.0
        )
        # Insert into tree.
        index.append(childNode)

    # Naming string.
    _NAME = len(names) - 1

    # Repeatedly join pairs until a full tree is made.
    while len(names) > 1:
        new_names = []

        while len(names) > 0:
            try:
                # Get next two children.
                child_names = []
                i = 0
                while i < 2:
                    child_names.append(names.pop())
                    i += 1

            except IndexError:
                # Only one child left. Just pass on the name.
                new_names.append(child_names[0])
                break

            # Join these children:
            # Make a parent node.
            parentName = str(_NAME)
            _NAME -= 1
            # Make a new Location.
            parentNode = Location(
                name=parentName,
                type="Branch"
            )
            print(f'New Parent {parentName}...')
            # Apply join.
            index.join(
                child_names[0],
                child_names[1],
                parentNode
            )

            # Append parent name to new_names.
            new_names.append(parentName)

        # New lot of children.
        names = new_names

    # Convert tree to newick format.
    return index.to_newick() + ";"


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


class Location:
    """Represents a node in the phylogeny."""
    is_root = False

    def __init__(self, name="", type="", dist=0.0, coord=None, accession_id=None, taxid=None, **kwargs):
        """__init__ 

        :param name: name of node, defaults to ""
        :type name: str, optional
        :param type: Leaf or Branch, defaults to ""
        :type type: str, optional
        :param dist: distance to parent node, defaults to 0.0
        :type dist: float, optional
        :param coord: binary coordinate from root to node, defaults to None
        :type coord: list, optional
        :param accession_id: NCBI accession id, defaults to None
        :type accession_id: str, optional
        :param taxid: NCBI taxonomy id, defaults to None
        :type taxid: int, optional
        :ivar name: node name
        :ivar type: "Leaf" or "Branch"
        :ivar distance: distance to parent node
        :ivar coordinate: list of binary binary numbers representing path from root to node
        :ivar nchildren: number of children below this node
        :ivar accession_id: NCBI accession id (only valid for leaves)
        :ivar taxid: NCBI taxonomy id
        """
        self._name = name
        self._type = type
        self._distance = dist
        self._coordinate = [] if coord is None else coord
        self._nchildren = 0
        self._accession_id = accession_id
        self._taxid = taxid

        self.children = []
        self.primaryChild = None

    def __str__(self):
        return (f"Location - {self._name} ({self._type}):"
                f"\n\t+ Coordinate: {self._coordinate}"
                f"\n\t+ Children: {self._nchildren}")

    def name():
        doc = "The name property."

        def fget(self):
            return self._name

        def fset(self, value):
            self._name = value

        def fdel(self):
            del self._name

        return locals()

    name = property(**name())

    def type():
        doc = "The type property."

        def fget(self):
            return self._type

        def fset(self, value):
            self._type = value

        def fdel(self):
            del self._type

        return locals()

    type = property(**type())

    def distance():
        doc = "The distance property."

        def fget(self):
            return self._distance

        def fset(self, value):
            self._distance = value

        def fdel(self):
            del self._distance

        return locals()

    distance = property(**distance())

    def coordinate():
        doc = "The coordinate property."

        def fget(self):
            return self._coordinate

        def fset(self, value):
            self._coordinate = value

        def fdel(self):
            del self._coordinate

        return locals()

    coordinate = property(**coordinate())

    def nchildren():
        doc = "The nchildren property."

        def fget(self):
            return self._nchildren

        def fset(self, value):
            self._nchildren = value

        def fdel(self):
            del self._nchildren

        return locals()

    nchildren = property(**nchildren())

    def accession_id():
        doc = "The accession id property."

        def fget(self):
            return self._accession_id

        def fset(self, value):
            self._accession_id = value

        def fdel(self):
            del self._accession_id

        return locals()

    accession_id = property(**accession_id())

    def taxid():
        doc = "The accession id property."

        def fget(self):
            return self._taxid

        def fset(self, value):
            self._taxid = value

        def fdel(self):
            del self._taxid

        return locals()

    taxid = property(**taxid())

    @classmethod
    def make_branch(cls, branch_name):
        return Location(branch_name, "Branch", 0)

    @classmethod
    def make_leaf(cls, tree, leaf_name, dist):
        return Location(leaf_name, "Leaf", dist)

    def add_child(self, child, primary_child=False):
        self.children.__setitem__(child.name, child)

        if primary_child:
            self.set_primary(child)

    def set_primary(self, node):
        self.primaryChild = node

    def save(self, out_dir):
        fname = f'{self._name}.loc'
        with open(os.path.join(out_dir, "phylogeny", "loc", fname), 'w') as f:
            json.dump(self.__dict__, f)


class Index:
    """ Phylogeny index that can load, save and manipulate Newick trees.
    """
    def __init__(self):
        self._pointers = {}  # Pointers is a map from node name to pool location.
        self.pool = [Location()]  # Empty location to represent an 'unclassified' location.

    def __len__(self):
        return len(self.pool)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'<Phylogeny Index, length={len(self)}>'

    @classmethod
    def load_newick(cls, path):
        """load_newick Load Newick tree from file.

        :param path: path to Newick file
        :type path: str
        :raises OSError: file does not exist
        :return: name of leaves and phylogeny Index object
        :rtype: list[str], expam.tree.Index
        """
        newick_str = ""

        if not os.path.exists(path):
            raise OSError("Couldn't find path %s." % path)

        with open(path, "r") as f:
            for line in f:
                newick_str += line.strip()

        return cls.from_newick(newick_str)

    @classmethod
    def from_newick(cls, newick_string):
        """from_newick Parse Newick string.

        :param newick_string: Newick string encoding tree.
        :type newick_string: str
        :return: name of leaves and phylogeny Index object
        :rtype: list[str], expam.tree.Index
        """

        # Remove whitespace.
        newick_string = newick_string.replace(" ", "")

        print("* Initialising node pool...")
        index_pool = cls.init_pool(newick_string)

        print("* Checking for polytomies...")
        cls.resolve_polytomies(index_pool)

        print("* Finalising index...")
        leaf_names, index = cls.from_pool(index_pool)

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

            j = i + 1

            while newick[j] not in NEWICK_PARSE:
                j += 1
            string = newick[i:j]

            # Check formatting.
            if force_digits:
                try:
                    string = float(string)

                except ValueError:
                    raise ValueError("Invalid distance declaration %s!" % newick[i:j])

            i = j - 1  # Update current index position.

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
    def from_pool(cls, pool):
        leaf_names = []

        # Initialise expam Index.
        index = cls()
        index.pool = pool

        branch_id = 1
        for i, node in enumerate(pool):
            if i == 0:
                continue

            if node.type == "Branch":
                # Set branch name.
                node.name = str(branch_id)
                branch_id += 1

            else:
                # Store leaf names.
                leaf_names.append(node.name)

            # Set name mapping.
            index._pointers[node.name] = i

        # Replace children objects with respective names.
        # Note: This requires the name mapping to be complete...
        for node in pool[1:]:
            if node.type == "Branch":
                i = 0
                while i < len(node.children):
                    node.children[i] = node.children[i].name

                    i += 1

        # Set (binary) coordinates.
        index.update_coordinates()

        # Update nchildren.
        index.update_nchildren()

        return leaf_names, index

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
        return str(node_name).lstrip('p')

    def __setitem__(self, key, value):
        key = self._format_node_name(key)  # Delete node identifier from results file.
        pool_index = self._pointers[key]
        self.pool[pool_index] = value

    def __getitem__(self, key):
        key = self._format_node_name(key)
        pool_index = self._pointers[key]
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
        for j in range(len(self.pool) - 1, 0, -1):
            node = self.pool[j]

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
        :rtype: list[str]
        """
        return list(self.yield_child_nodes(node_name))

    def get_child_leaves(self, node_name):
        """get_child_leaves Get list of leaves at or below node_name.

        :param node_name: name of node
        :type node_name: str
        :return: list of leaf names
        :rtype: list[str]
        """
        return list(self.yield_leaves(node_name))

    def draw_results(self, file_path, out_dir, skiprows=None, groups=None, cutoff=None, cpm=None, colour_list=None,
                     name_to_taxon=None, use_phyla=False, keep_zeros=True, use_node_names=True):
        counts = pd.read_csv(file_path, sep='\t', index_col=0, header=0, skiprows=skiprows)
        self.draw_tree(out_dir, counts=counts, groups=groups, cutoff=cutoff, cpm=cpm, colour_list=colour_list,
                       name_to_taxon=name_to_taxon, use_phyla=use_phyla, keep_zeros=keep_zeros,
                       use_node_names=use_node_names)

    def draw_tree(self, out_dir, counts=None, groups=None, cutoff=None, cpm=None, colour_list=None, name_to_taxon=None,
                  use_phyla=False, keep_zeros=True, use_node_names=True):
        try:
            import ete3.coretype.tree
            from ete3 import AttrFace, CircleFace, faces, Tree, TreeStyle, NodeStyle, TextFace

        except ModuleNotFoundError as e:
            print("Could not import ete3 plotting modules! Error raised:")
            print(traceback.format_exc())
            print("Skipping plotting...")

            return

        """
        Phylogenetic printing of nodes.
        """
        for node in self.pool:
            if node.type == 'Branch' and node.name[0] != 'p':
                node.name = 'p' + node.name

        """
        Create PhyloTree.
        """
        newick_string = self.to_newick()
        tree = Tree(newick_string, format=1)

        """
        Take care of groupings.
        """
        if groups is None:
            groups = [(None, (section,)) for section in counts.columns]

        # *** Relies on counts being defined.
        if counts is not None:
            for _, group in groups:
                if len(group) > 1:
                    group_name = group[0]

                    # Conglomerate counts within the group.
                    counts.loc[:, group_name] = counts[list(group)].sum(axis=1)
                    counts.drop(labels=list(groups[1:]), axis=1, inplace=True)

            sections = counts.columns.tolist()
            nodes = counts.index.tolist()

            """
            Employ cutoff.
            """
            if cutoff is None and cpm is None:
                nodes_with_counts = nodes

            else:
                nodes_with_counts = set()

                for section in sections:
                    total = sum(counts[section])
                    section_cutoff = max(cutoff, (total / 1e6) * cpm)

                    for index in nodes:
                        if index not in nodes_with_counts:
                            if np.any(counts.loc[index, :] >= section_cutoff):
                                nodes_with_counts.add(index)

                nodes_with_counts = list(nodes_with_counts)

            """
            Attempt pruning of tree.
            """
            try:
                if not keep_zeros:
                    tree.prune(nodes_with_counts, preserve_branch_length=True)

            except ete3.coretype.tree.TreeError:
                print("Tree pruning failed, aborting tree drawing!")
                return

            """
            Colour generation.
            """
            max_vector = tuple(counts.max())
            colours = [colour for colour, _ in groups]
            colour_generator = ColourGenerator(colours, max_vector, colour_list=colour_list)

        """
        Function to render any given node.
        """

        def layout(node):
            nonlocal name_to_taxon

            if node.is_leaf():
                if use_node_names:
                    faces.add_face_to_node(AttrFace('name', fsize=20), node, column=0, position='aligned')

                if use_phyla:
                    for phyla, colour in PHYLA_COLOURS:
                        lineage = name_to_taxon[node.name]

                        if phyla in lineage:
                            tax_face = TextFace("   ")
                            tax_face.background.color = colour
                            faces.add_face_to_node(tax_face, node, column=1, position='aligned')

            if (counts is not None) and (node.name in counts.index):
                node_counts = tuple(counts.loc[node.name, :])

                if any(node_counts):
                    colour, _ = colour_generator.generate(node_counts)

                    ns = NodeStyle()
                    ns['bgcolor'] = colour
                    node.set_style(ns)

        """
        Render tree.
        """
        ts = TreeStyle()
        ts.mode = "c"
        ts.show_leaf_name = False
        ts.layout_fn = layout
        ts.force_topology = False
        ts.allow_face_overlap = False
        ts.draw_guiding_lines = True
        ts.root_opening_factor = 1
        tree.render(
            out_dir,
            tree_style=ts,
            dpi=300
        )
        print("Phylogenetic tree written to %s!" % out_dir)


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
    def __init__(self, rng=None, seed=1):
        self.alpha = "0123456789ABCDEF"

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
    def __init__(self, colours, max_vector, colour_list=None):
        """

        :param colours:
        :param max_vector:
        :param colour_list:
        """
        """
        Ensure we have enough colours for plotting.
        """
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

    def _calculate_colour(self, score_vector, scorer):
        colour = HexColour("#FFFFFF")  # White background.

        # Artificial scaling factor.
        non_zero_indices = [i for i, score in enumerate(score_vector) if score > 0]
        alpha = sum([score_vector[i] for i in non_zero_indices]) / sum([scorer[i] for i in non_zero_indices])

        for i, score in enumerate(score_vector):
            if score > 0:
                colour += self.colours[i].opaque(score / scorer[i])

        return str(colour), alpha

    def generate(self, vector):
        return self._calculate_colour(vector, self.max_vector)
