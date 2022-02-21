// mfil.c

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>

typedef struct {
    uint64_t    *array;
    int64_t     index;
} node_t;

typedef struct {
    node_t      *nodes;     // array of nodes.
    uint64_t    k;          // # leaves.
    uint64_t    chunks;     // wrt. node arrays.
    uint64_t    size;       // # nodes, including extra at 0.
} tree_t;

typedef struct {
    uint64_t    *array;
    int64_t     index;
} Output;

/* Forward declarations. */

// Array representation of binary tree.
static void init_tree(tree_t *tree, uint64_t k, uint64_t chunks, uint64_t **arr_ptrs, int64_t *indices);
static void set_leaves(tree_t *tree, uint64_t **arr_ptrs, int64_t *indices);
static void fill_tree(tree_t *tree);
static node_t traverse_and_compare(tree_t *tree, int64_t node_index);
static void write_to_node(tree_t *tree, int64_t dest, node_t target);
static node_t get_node(tree_t *tree, int64_t node_index);
static node_t play_game(tree_t *tree, int64_t node_index, node_t one, node_t two);
static bool compare(node_t one, node_t two, uint64_t chunks);

// K-way merging.
static void k_way_merge(uint64_t **arr_ptrs, int64_t *indices, uint64_t n, uint64_t k, uint64_t chunks);
static bool has_more(tree_t *tree);
static void write_to_out(tree_t *tree, Output *out);
static node_t get_next_player(tree_t *tree);
static int64_t get_next_index(tree_t *tree);
static void next_game(tree_t *tree, node_t next_player, int64_t node_index);

// Linear culling.
static int64_t multidisjoint(int n_arrays, int n_columns, uint64_t **arr_ptrs, int64_t *sizes);

/** Array representation of binary tree.
*   @param  k           Number of arrays to be traversed.
*   @param  arr_ptrs    Pointer on each array to be traversed.
*   @param  indices     Length of each array to be traversed.
*   @param  chunks      Breadth of each array.
*   @param  return      Binary tournament loser tree.
*
*   Given the tree is indexed using integers, the parent node of any given
*   index (i) in the tree is (i/2). There are (2k-1) nodes in the binary tree,
*   but account for extra node at top to store winner.
*/
static void init_tree(tree_t *tree, uint64_t k, uint64_t chunks, uint64_t **arr_ptrs, int64_t *indices) {
    tree->k         = k;
    tree->size      = 2 * k;
    tree->chunks    = chunks;
    tree->nodes     = malloc(sizeof(node_t) * tree->size);

    /* Set the leaves for game preparation. */
    set_leaves(tree, arr_ptrs, indices);

    /* Propagate games upward. */
    fill_tree(tree);
}

static void set_leaves(tree_t *tree, uint64_t **arr_ptrs, int64_t *indices) {
    for (uint64_t i=tree->k; i<tree->size; ++i) {

        /* Pointer on each numpy array/memoryview. */
        tree->nodes[i].array = arr_ptrs[i-tree->k];
        tree->nodes[i].index = indices[i-tree->k];
    }
}

static void fill_tree(tree_t *tree) {
    node_t  winner;

    winner = traverse_and_compare(tree, 1);
    write_to_node(tree, 0, winner);
}

static node_t traverse_and_compare(tree_t *tree, int64_t node_index) {
    node_t  one, two, winner;

    if (node_index >= tree->k)
        return get_node(tree, node_index);

    /* First fill children recursively. */
    one = traverse_and_compare(tree, 2 * node_index);
    two = traverse_and_compare(tree, (2 * node_index) + 1);

    /* Compare my children; write the loser and carry the winner. */
    winner = play_game(tree, node_index, one, two);

    return winner;
}

static void write_to_node(tree_t *tree, int64_t dest, node_t target) {
    tree->nodes[dest].array = target.array;
    tree->nodes[dest].index = target.index;
}

static node_t get_node(tree_t *tree, int64_t node_index) {
    node_t  node;

    if (node_index < tree->k)
        return tree->nodes[node_index];

    else if (tree->nodes[node_index].index < 0) {
        node.array = NULL;
        node.index = -1;

        return node;

    } else {
        node.array = tree->nodes[node_index].array + (tree->nodes[node_index].index * tree->chunks);
        node.index = node_index;

        return node;
    }
}

static node_t play_game(tree_t *tree, int64_t node_index, node_t one, node_t two) {
    /* Write the loser to that node, and return the winner. */

    if (one.index < 0) {
        write_to_node(tree, node_index, one);
        return two;

    } else if (two.index < 0) {
        write_to_node(tree, node_index, two);
        return one;

    } else if (compare(one, two, tree->chunks)) {
        write_to_node(tree, node_index, two);
        return one;

    } else {
        write_to_node(tree, node_index, one);
        return two;
    }
}

static bool compare(node_t one, node_t two, uint64_t chunks) {
    /* Look for larger kmer to win. */
    for (uint64_t x=0; x<chunks; ++x) {
        if (one.array[x] > two.array[x])
            return true;

        else if (two.array[x] > one.array[x])
            return false;

        else
            /* This chunk equal, look to next chunk. */
            continue;
    }

    return false;
}

/** K-way merge (Knuth, 1998)
*   @param  arr_ptrs    Pointer on each array to be merged.
*   @param  indices     Array of lengths of each array.
*   @param  n           Size of final array (i.e. sum(indices)).
*   @param  k           Number of arrays to be merged.
*   @param  chunks      Breadth of each array to be merged.
*   @return None        Merge done in-place on first array.
*
*   Merge k sorted arrays, whose values are distinct with respect to the
*   other arrays, into the first array. This first array has been extended
*   with 0 values to a size n, but the portion that still contains values is
*   stored in the first element of indices.
*/
static void k_way_merge(uint64_t **arr_ptrs, int64_t *indices, uint64_t n, uint64_t k, uint64_t chunks) {
    tree_t      tree;
    node_t      next_player;
    uint64_t    next_index;
    Output      out;

    /* Set output. */
    out.array = arr_ptrs[0];
    out.index = n - 1;

    /* Start competition tree. */
    init_tree(&tree, k, chunks, arr_ptrs, indices);

    /* Fill arrays. */
    while (has_more(&tree)) {
        /* Write current winner to output. */
        write_to_out(&tree, &out);

        /* Start next game. */
        next_index = get_next_index(&tree);
        next_player = get_next_player(&tree);
        next_game(&tree, next_player, next_index);
    }

    /* Free malloc tree nodes. */
    free(tree.nodes);
}

static bool has_more(tree_t *tree) {
    if (tree->nodes[0].index < 0)
        return false;

    return true;
}

static void write_to_out(tree_t *tree, Output *out) {
    for (uint64_t j=0; j<tree->chunks; ++j)
        *(out->array + (out->index * tree->chunks) + j) = tree->nodes[0].array[j];

    out->index--;
}

static node_t get_next_player(tree_t *tree) {
    int64_t leaf_index;

    leaf_index = tree->nodes[0].index;
    tree->nodes[leaf_index].index--;

    return get_node(tree, leaf_index);
}

static int64_t get_next_index(tree_t *tree) {
    return tree->nodes[0].index / 2;
}

static void next_game(tree_t *tree, node_t next_player, int64_t node_index) {
    node_t  node;

    /* Compare and take winner. */
    node = get_node(tree, node_index);
    node = play_game(tree, node_index, node, next_player);

    if (node_index == 1)
        /* Write champion. */
        write_to_node(tree, 0, node);
    else
        /* Send winner to next match. */
        next_game(tree, node, node_index / 2);
}

/** Linear culling.
*   @param  n_arrays
*   @param  n_columns
*   @param  arr_ptrs
*   @param  sizes
*
*   Remove elements from first array in place, which also appear in the other
*   arrays pointed to by arr_ptrs. Items are search for in a simultaneous
*   linear scan of all arrays.
*/
static int64_t multidisjoint(int n_arrays, int n_columns, uint64_t **arr_ptrs, int64_t *sizes) {
    bool seen = false, smaller;  // Whether kmer appears in any array.
    int64_t final_size;  // Final size of array after making disjoint.

    // Initialise place to store counters on each array.
    int64_t *upto = (int64_t *)malloc(n_arrays * sizeof(int64_t));
    for (int64_t m=0; m<n_arrays; m++) *(upto + m) = 0;

    // Disjoint loop.
    for (int64_t i=0; i<*(sizes + 0); i++) {
        seen = false;
        // Check each sub array.
        for (int j=1; j<n_arrays; j++) {
            // Compare values within array to the main array.
            smaller = true;
            while (smaller) {
                // Check it is still a valid index.
                if (*(upto + j) >= *(sizes + j)) break;

                // Compare kmers.
                for (int w=0; w<n_columns; w++) {
                    if (*(*(arr_ptrs + j) + *(upto + j)*n_columns + w) > *(*(arr_ptrs + 0) + i*n_columns + w)) {
                        smaller = false;
                        break;
                    } else if (*(*(arr_ptrs + j) + *(upto + j)*n_columns + w) < *(*(arr_ptrs + 0) + i*n_columns + w)) {
                        *(upto + j) += 1;
                        break;
                    } else {
                        if (w+1 == n_columns) {
                            seen = true;
                            smaller = false;
                        }
                    }
                }
            }
            // Check if we saw it.
            if (seen) break;
        }
        // Only write if it's not seen.
        if (!seen) {
            // Overwrite that previous value at counter.
            for (int w=0; w<n_columns; w++)
                *(*(arr_ptrs + 0) + *(upto + 0)*n_columns + w) = *(*(arr_ptrs + 0) + i*n_columns + w);
            // Increment counter.
            *(upto + 0) += 1;
        }
    }
    final_size = *(upto + 0);  // Save final size.
    free(upto);  // Free counter array.

    // Return final size.
    return final_size;
}
