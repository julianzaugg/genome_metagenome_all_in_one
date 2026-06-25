#!/usr/bin/env python

import argparse
import re

from ete3 import Tree

def find_closest_leaves(tree, query_id, exclude_pattern=None, num_neighbours=1):
    query_node = tree.search_nodes(name=query_id)
    if not query_node:
        print(f"Leaf with ID {query_id} not found")
        return []

    query_node = query_node[0]
    distances = []

    for leaf in tree.iter_leaves():
        if leaf.name == query_id:
            continue
        if exclude_pattern and re.search(exclude_pattern, leaf.name):
            continue

        distance = query_node.get_distance(leaf)
        distances.append((leaf, distance))

    distances.sort(key=lambda x: x[1])
    return [leaf for leaf, dist in distances[:num_neighbours]]

def main(tree_file, ids, exclude_pattern, num_neighbours, print_mode):
    # Load the tree from the file
    tree = Tree(tree_file,format=1, quoted_node_names=True)

    # Process each query ID
    for query_id in ids:
        closest_leaves = find_closest_leaves(tree, query_id, exclude_pattern, num_neighbours)
        if print_mode == "long":
            print('\n'.join(f"{query_id}\t{string}" for string in [leaf.name for leaf in closest_leaves]))
        else:
            if closest_leaves:
                print(f"Closest leaves to {query_id}: {[leaf.name for leaf in closest_leaves]}")
            else:
                print(f"No closest leaves found for {query_id}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process some leaf IDs for a tree.")
    parser.add_argument('-t', '--tree', metavar = 'TREE', help='Input tree in newick format.', required=True)
    parser.add_argument('ids', metavar='LEAF_IDS', type=str, help='Comma-separated leaf IDs to process')
    parser.add_argument('--exclude', metavar='PATTERN', type=str, help='Pattern to exclude from leaf names', default=None)
    parser.add_argument('-n', '--n_neighbours', help='Number of neighbours to get for each query node', required=False, default=1, type=int)
    parser.add_argument('--print_mode', metavar='PRINT_MODE', type=str, help='Print in long mode or default', default=None)

    args = parser.parse_args()
    main(args.tree, args.ids.split(','), args.exclude, args.n_neighbours, args.print_mode)
