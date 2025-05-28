import numpy as np
from algorithm import needleman_wunsch, traceback, needleman_wunsch_directions, similarity_score
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os

class Node:
    """
    A class representing a node in a phylogenetic tree.

    Attributes:
        name (str): The name of the node
        left (Node): The left child node
        right (Node): The right child node
        distance (float): The distance from this node to its parent
        height (float): The height of this node in the tree
        is_leaf (bool): Whether this node is a leaf node
    """
    def __init__(self, name, left=None, right=None, distance=0, height=0):
        self.name = name
        self.left = left
        self.right = right
        self.distance = distance
        self.height = height
        self.is_leaf = left is None and right is None

    def __str__(self):
        return self.name

def calculate_distance_matrix(sequences, match_score=1, mismatch_penalty=-1, gap_penalty=-1):
    """
    Calculate a distance matrix from a set of sequences using the Needleman-Wunsch algorithm.

    Args:
        sequences (dict): A dictionary mapping sequence names to sequences
        match_score (int): Score for matching characters
        mismatch_penalty (int): Penalty for mismatched characters
        gap_penalty (int): Penalty for introducing a gap

    Returns:
        tuple: (distance_matrix, sequence_names)
            - distance_matrix (numpy.ndarray): The distance matrix
            - sequence_names (list): The names of the sequences
    """
    sequence_names = list(sequences.keys())
    n = len(sequence_names)
    distance_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i+1, n):
            seq1 = sequences[sequence_names[i]]
            seq2 = sequences[sequence_names[j]]

            # Calculate alignment using Needleman-Wunsch
            numeric_matrix = needleman_wunsch(seq1, seq2, match_score, mismatch_penalty, gap_penalty)
            directional_matrix = needleman_wunsch_directions(numeric_matrix, match_score, mismatch_penalty, gap_penalty)
            _, aligned_seq1, aligned_seq2 = traceback(numeric_matrix, directional_matrix, seq1, seq2)

            # Calculate distance as 1 - similarity
            similarity = similarity_score(aligned_seq1, aligned_seq2)
            distance = 1 - similarity

            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance

    return distance_matrix, sequence_names

def upgma(distance_matrix, names):
    """
    Implement the UPGMA algorithm to build a phylogenetic tree.

    Args:
        distance_matrix (numpy.ndarray): The distance matrix
        names (list): The names of the sequences

    Returns:
        Node: The root node of the phylogenetic tree
    """
    n = len(names)

    # Initialize clusters with leaf nodes
    clusters = [Node(name) for name in names]

    # Make a copy of the distance matrix to modify
    dist_matrix = distance_matrix.copy()

    # Store heights of nodes
    heights = np.zeros(n)

    # Main UPGMA algorithm
    for _ in range(n-1):
        # Find the minimum distance in the matrix while ignoring the diagonal
        # Create a mask for the diagonal elements
        mask = np.eye(len(dist_matrix), dtype=bool)
        # Create a masked version of the distance matrix
        masked_dist_matrix = np.ma.array(dist_matrix, mask=mask)
        # Find the minimum value and its indices
        min_i, min_j = np.unravel_index(np.ma.argmin(masked_dist_matrix), dist_matrix.shape)

        # Create a new node
        new_height = dist_matrix[min_i, min_j] / 2
        left_distance = new_height - heights[min_i]
        right_distance = new_height - heights[min_j]

        clusters[min_i].distance = left_distance
        clusters[min_j].distance = right_distance

        new_node = Node(f"({clusters[min_i].name},{clusters[min_j].name})", 
                        left=clusters[min_i], 
                        right=clusters[min_j],
                        height=new_height)

        # Update the distance matrix
        new_dist = np.zeros(len(dist_matrix) - 1)
        for k in range(len(dist_matrix)):
            if k != min_i and k != min_j:
                idx = k if k < min_j else k - 1
                new_dist[idx] = (dist_matrix[min_i, k] + dist_matrix[min_j, k]) / 2

        # Create new distance matrix
        new_matrix = np.zeros((len(dist_matrix) - 1, len(dist_matrix) - 1))

        # Copy values from old matrix, skipping min_j
        for i in range(len(dist_matrix)):
            if i == min_j:
                continue
            i_new = i if i < min_j else i - 1
            for j in range(len(dist_matrix)):
                if j == min_j:
                    continue
                j_new = j if j < min_j else j - 1
                new_matrix[i_new, j_new] = dist_matrix[i, j]

        # Update the min_i row and column with new distances
        for k in range(len(new_dist)):
            if k != min_i and (min_i < min_j or k < min_j):
                new_matrix[min_i, k] = new_dist[k]
                new_matrix[k, min_i] = new_dist[k]

        # Update clusters list
        clusters = [clusters[i] for i in range(len(clusters)) if i != min_j]
        clusters[min_i] = new_node

        # Update heights
        heights = np.delete(heights, min_j)
        heights[min_i] = new_height

        # Update distance matrix
        dist_matrix = new_matrix

    # Return the root node
    return clusters[0]

def save_tree_to_file(root, filename="phylogenetic_tree.txt", distance_matrix=None, sequence_names=None):
    """
    Save the phylogenetic tree to a text file.

    Args:
        root (Node): The root node of the phylogenetic tree
        filename (str): The name of the file to save to
        distance_matrix (numpy.ndarray, optional): The distance matrix
        sequence_names (list, optional): The names of the sequences
    """
    with open(filename, "w", encoding="utf-8") as f:
        f.write("Phylogenetic Tree (UPGMA)\n")
        f.write("=======================\n\n")

        if distance_matrix is not None and sequence_names is not None:
            f.write("Distance Matrix:\n")
            f.write("---------------\n")
            f.write("\t" + "\t".join(sequence_names) + "\n")
            for i, name in enumerate(sequence_names):
                f.write(f"{name}\t" + "\t".join([f"{distance_matrix[i, j]:.4f}" for j in range(len(sequence_names))]) + "\n")
            f.write("\n")

        f.write("Tree in Newick format:\n")
        f.write("--------------------\n")
        f.write(get_newick_format(root) + ";\n\n")

        f.write("Tree structure:\n")
        f.write("--------------\n")
        print_tree(root, f)

    print(f"Tree saved to {filename}")

def get_newick_format(node):
    """
    Convert a tree to Newick format.

    Args:
        node (Node): The root node of the tree

    Returns:
        str: The tree in Newick format
    """
    if node.is_leaf:
        return node.name
    else:
        left = get_newick_format(node.left)
        right = get_newick_format(node.right)
        return f"({left}:{node.left.distance},{right}:{node.right.distance})"

def print_tree(node, file=None, prefix="", is_last=True):
    """
    Print a tree structure.

    Args:
        node (Node): The node to print
        file (file, optional): The file to write to
        prefix (str): The prefix to use for indentation
        is_last (bool): Whether this is the last child of its parent
    """
    if file:
        file.write(prefix + ("└── " if is_last else "├── ") + str(node) + "\n")
    else:
        print(prefix + ("└── " if is_last else "├── ") + str(node))

    prefix += "    " if is_last else "│   "

    if node.left:
        print_tree(node.left, file, prefix, node.right is None)
    if node.right:
        print_tree(node.right, file, prefix, True)

def plot_tree(root, filename="phylogenetic_tree.png"):
    """
    Plot a phylogenetic tree using matplotlib.

    Args:
        root (Node): The root node of the phylogenetic tree
        filename (str): The name of the file to save the plot to
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # Set up the plot
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Plot the tree
    leaf_count = count_leaves(root)
    leaf_positions = {}
    plot_node(root, ax, 0.1, 0.9, 0.1, 0.9, leaf_count, leaf_positions)

    # Save the plot
    plt.savefig(filename)
    plt.close()

    print(f"Tree plot saved to {filename}")

def count_leaves(node):
    """Count the number of leaf nodes in a tree."""
    if node.is_leaf:
        return 1
    return count_leaves(node.left) + count_leaves(node.right)

def plot_node(node, ax, x_min, x_max, y_min, y_max, leaf_count, leaf_positions, depth=0):
    """
    Recursively plot a node and its children.

    Args:
        node (Node): The node to plot
        ax (matplotlib.axes.Axes): The axes to plot on
        x_min, x_max, y_min, y_max (float): The boundaries of the plot area
        leaf_count (int): The total number of leaf nodes
        leaf_positions (dict): A dictionary mapping leaf names to their positions
        depth (int): The current depth in the tree
    """
    if node.is_leaf:
        # Plot a leaf node
        y = (y_min + y_max) / 2
        leaf_positions[node.name] = (x_max, y)
        ax.text(x_max + 0.01, y, node.name, va='center')
        return

    # Calculate positions for children
    left_leaves = count_leaves(node.left)
    right_leaves = count_leaves(node.right)

    left_ratio = left_leaves / leaf_count
    y_mid = y_min + (y_max - y_min) * left_ratio

    # Plot left child
    plot_node(node.left, ax, x_min, x_max - node.left.distance, y_min, y_mid, leaf_count, leaf_positions, depth + 1)

    # Plot right child
    plot_node(node.right, ax, x_min, x_max - node.right.distance, y_mid, y_max, leaf_count, leaf_positions, depth + 1)

    # Calculate node position
    x = x_max - max(node.left.distance, node.right.distance)
    y_left = (y_min + y_mid) / 2
    y_right = (y_mid + y_max) / 2

    # Draw lines to children
    left_x, left_y = leaf_positions.get(node.left.name, (x_max - node.left.distance, y_left))
    right_x, right_y = leaf_positions.get(node.right.name, (x_max - node.right.distance, y_right))

    # Horizontal lines
    ax.add_line(Line2D([x, left_x], [y_left, y_left], color='black'))
    ax.add_line(Line2D([x, right_x], [y_right, y_right], color='black'))

    # Vertical line
    ax.add_line(Line2D([x, x], [y_left, y_right], color='black'))

    # Store node position
    leaf_positions[node.name] = (x, (y_left + y_right) / 2)

def load_sequences_from_fasta(file_path):
    """
    Load sequences from a FASTA file.

    Args:
        file_path (str): The path to the FASTA file

    Returns:
        dict: A dictionary mapping sequence names to sequences
    """
    sequences = {}
    current_name = None

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_name = line[1:]  # Remove the '>' character
                sequences[current_name] = ""
            elif current_name is not None:
                sequences[current_name] += line

    return sequences

def load_distance_matrix_from_file(file_path):
    """
    Load a distance matrix from a file.

    Args:
        file_path (str): The path to the file

    Returns:
        tuple: (distance_matrix, sequence_names)
            - distance_matrix (numpy.ndarray): The distance matrix
            - sequence_names (list): The names of the sequences
    """
    with open(file_path, "r") as f:
        lines = f.readlines()

    # First line contains sequence names
    sequence_names = lines[0].strip().split()
    n = len(sequence_names)

    # Create distance matrix
    distance_matrix = np.zeros((n, n))

    # Parse the matrix
    for i in range(n):
        values = lines[i+1].strip().split()
        for j in range(n):
            distance_matrix[i, j] = float(values[j])

    return distance_matrix, sequence_names

def main_upgma(sequences=None, distance_matrix=None, sequence_names=None, 
               match_score=1, mismatch_penalty=-1, gap_penalty=-1,
               output_file="phylogenetic_tree.txt", output_image="phylogenetic_tree.png"):
    """
    Main function to run the UPGMA algorithm.

    Args:
        sequences (dict, optional): A dictionary mapping sequence names to sequences
        distance_matrix (numpy.ndarray, optional): The distance matrix
        sequence_names (list, optional): The names of the sequences
        match_score (int): Score for matching characters
        mismatch_penalty (int): Penalty for mismatched characters
        gap_penalty (int): Penalty for introducing a gap
        output_file (str): The name of the file to save the tree to
        output_image (str): The name of the file to save the tree plot to

    Returns:
        Node: The root node of the phylogenetic tree
    """
    # Calculate distance matrix if not provided
    if distance_matrix is None and sequences is not None:
        distance_matrix, sequence_names = calculate_distance_matrix(
            sequences, match_score, mismatch_penalty, gap_penalty
        )

    # Run UPGMA
    root = upgma(distance_matrix, sequence_names)

    # Save results
    save_tree_to_file(root, output_file, distance_matrix, sequence_names)
    plot_tree(root, output_image)

    return root
