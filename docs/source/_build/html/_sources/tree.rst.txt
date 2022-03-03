expam's Tree module
===================

A programmatic API to interact with phylogenetic trees, particularly those used in reference databases.

expam.tree.Location
-------------------

.. autoclass:: expam.tree.Location

.. autofunction:: expam.tree.Location.__init__

expam.tree.Index
----------------

.. autoclass:: expam.tree.Index

.. autofunction:: expam.tree.Index.load_newick

.. autofunction:: expam.tree.Index.from_newick

    Example loading an Index object from a Newick string.

    .. code-block:: python

        >>> from expam.tree import Index
        >>> tree_string = "(B:6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);"
        >>> leaves, index = Index.from_newick(tree_string)
        * Initialising node pool...
        * Checking for polytomies...
            Polytomy (degree=3) detected! Resolving...
            Polytomy (degree=3) detected! Resolving...
        * Finalising index...
        >>> leaves
        ['B', 'A', 'C', 'E', 'D']
        >>> index
        <Phylogeny Index, length=10>
        >>> index['A']
        <expam.tree.Location object at 0x109ac7970>
        >>> index['A'].name
        'A'
        >>> index['A'].coordinate
        [0, 0, 1, 0]

.. autofunction:: expam.tree.Index.resolve_polytomies

.. autofunction:: expam.tree.Index.coord

.. autofunction:: expam.tree.Index.to_newick

.. autofunction:: expam.tree.Index.yield_child_nodes

    .. code-block:: python

        >>> for node in index.yield_child_nodes('p1'):  # p1 will always be the root
        ...    print(node)
        ... 
        1
        D
        2
        B
        3
        E
        4
        A
        C

    .. note:: 

        Internal node (branch) names can start with 'p', but this may also be neglected.


.. autofunction:: expam.tree.Index.yield_leaves

.. autofunction:: expam.tree.Index.get_child_nodes

    .. code-block:: python

        >>> index.get_child_nodes('1')
        ['1', 'D', '2', 'B', '3', 'E', '4', 'A', 'C']
        >>> index.get_child_nodes('E')
        ['E']

.. autofunction:: expam.tree.Index.get_child_leaves
