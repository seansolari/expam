from expam.tree.location import Location
from expam.tree.tree import Index


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
