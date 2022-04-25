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

