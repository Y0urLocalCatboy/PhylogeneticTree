import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Node import Node



def calculate_distance_matrix(sequences):
    """
    Calculate a distance matrix from a set of sequences using a simple distance metric.

    Args:
        sequences (dict): A dictionary mapping sequence names to sequences

    Returns:
        tuple: (distance_matrix, sequence_names)
            - distance_matrix: The distance matrix
            - sequence_names: The names of the sequences
    """
    sequence_names = list(sequences.keys())
    n = len(sequence_names)
    distance_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i+1, n):
            seq1 = sequences[sequence_names[i]]
            seq2 = sequences[sequence_names[j]]

            distance = calculate_sequence_distance(seq1, seq2)

            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance

    return distance_matrix, sequence_names

def calculate_sequence_distance(seq1, seq2):
    """
    Calculate the distance between two sequences.
    For sequences of equal length, uses Hamming distance.
    For sequences of different lengths, uses a normalized edit distance.

    Args:
        seq1 (str): First sequence
        seq2 (str): Second sequence

    Returns:
        float: Distance between the sequences (0 to 1)
    """
    # If sequences are of equal length, use Hamming distance
    if len(seq1) == len(seq2):
        mismatches = sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
        return mismatches / len(seq1)


    len_diff = abs(len(seq1) - len(seq2))

    min_len = min(len(seq1), len(seq2))
    mismatches = sum(c1 != c2 for c1, c2 in zip(seq1[:min_len], seq2[:min_len]))

    total_distance = (mismatches + len_diff) / max(len(seq1), len(seq2))

    return total_distance

def upgma(distance_matrix, names):
    """
    UPGMA algorithm to build a phylogenetic tree.

    Args:
        distance_matrix: The distance matrix
        names: The names of the sequences

    Returns:
        Node: The root node of the phylogenetic tree
    """
    n = len(names)

    clusters = [Node(name) for name in names]

    dist_matrix = distance_matrix.copy()

    heights = np.zeros(n)


    while len(clusters) > 1:
        mask = np.eye(len(dist_matrix), dtype=bool)

        masked_dist_matrix = np.ma.array(dist_matrix, mask=mask)
        min_i, min_j = np.unravel_index(np.ma.argmin(masked_dist_matrix), dist_matrix.shape)

        new_height = dist_matrix[min_i, min_j] / 2
        left_distance = new_height - heights[min_i]
        right_distance = new_height - heights[min_j]

        clusters[min_i].distance = left_distance
        clusters[min_j].distance = right_distance

        new_node = Node(f"({clusters[min_i].name},{clusters[min_j].name})",
                        left=clusters[min_i],
                        right=clusters[min_j],
                        height=new_height)

        new_dist = np.zeros(len(dist_matrix) - 1)
        for k in range(len(dist_matrix)):
            if k != min_i and k != min_j:
                idx = k if k < min_j else k - 1
                new_dist[idx] = (dist_matrix[min_i, k] + dist_matrix[min_j, k]) / 2

        new_matrix = np.zeros((len(dist_matrix) - 1, len(dist_matrix) - 1))

        for i in range(len(dist_matrix)):
            if i == min_j:
                continue
            i_new = i if i < min_j else i - 1
            for j in range(len(dist_matrix)):
                if j == min_j:
                    continue
                j_new = j if j < min_j else j - 1
                new_matrix[i_new, j_new] = dist_matrix[i, j]

        for k in range(len(new_dist)):
            if k != min_i and (min_i < min_j or k < min_j):
                new_matrix[min_i, k] = new_dist[k]
                new_matrix[k, min_i] = new_dist[k]

        clusters = [clusters[i] for i in range(len(clusters)) if i != min_j]
        clusters[min_i] = new_node

        heights = np.delete(heights, min_j)
        heights[min_i] = new_height

        dist_matrix = new_matrix

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
        print_tree_structure(root, f)

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


def print_tree_structure(node, file=None, prefix="", is_last=True):  # Renamed from print_tree
    """
    Print a tree structure to console or file.

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

    children = []
    if node.left:
        children.append(node.left)
    if node.right:
        children.append(node.right)

    for i, child in enumerate(children):
        is_child_last = (i == len(children) - 1)
        new_prefix = prefix + ("    " if is_last else "│   ")
        print_tree_structure(child, file, new_prefix, is_child_last)


def plot_tree(root, filename="phylogenetic_tree.png"):
    """
    Plot a phylogenetic tree using matplotlib with improved labels and branch colors.
    The tree grows upwards, with leaf nodes at y=0 and root at maximum height.

    Args:
        root (Node): The root node of the phylogenetic tree
        filename (str): The name of the file to save the plot to
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.axis('off')

    if not root:
        print("Tree is empty, cannot plot.")
        plt.close(fig)
        return

    leaf_x_coords_map = {}
    current_x_ord = 0

    def assign_leaf_x_recursive(node):
        nonlocal current_x_ord
        if node.is_leaf:
            leaf_x_coords_map[node] = current_x_ord
            current_x_ord += 1
        else:
            if node.left: assign_leaf_x_recursive(node.left)
            if node.right: assign_leaf_x_recursive(node.right)

    assign_leaf_x_recursive(root)

    num_leaves = current_x_ord
    if num_leaves == 0:
        num_leaves = 1


    node_positions = {}

    def assign_node_positions_recursive(node):
        if node.is_leaf:
            x_pos = (leaf_x_coords_map[node] / (num_leaves - 1)) if num_leaves > 1 else 0.5
            node_positions[node] = (x_pos, node.height)
        else:
            if node.left: assign_node_positions_recursive(node.left)
            if node.right: assign_node_positions_recursive(node.right)

            x_left = node_positions[node.left][0] if node.left and node.left in node_positions else -1
            x_right = node_positions[node.right][0] if node.right and node.right in node_positions else -1

            parent_x = 0.5
            if node.left and node.left in node_positions and node.right and node.right in node_positions:
                parent_x = (x_left + x_right) / 2
            elif node.left and node.left in node_positions:  # left child
                parent_x = x_left
            elif node.right and node.right in node_positions:  # right child
                parent_x = x_right

            node_positions[node] = (parent_x, node.height)
        return node_positions[node]

    assign_node_positions_recursive(root)

    all_x = [pos[0] for pos in node_positions.values()]
    all_y = [pos[1] for pos in node_positions.values()]

    min_x_coord, max_x_coord = (min(all_x), max(all_x)) if all_x else (0, 1)
    min_y_coord, max_y_coord = (min(all_y), max(all_y)) if all_y else (0, 1)

    x_range = max_x_coord - min_x_coord if max_x_coord > min_x_coord else 1.0
    y_range = max_y_coord - min_y_coord if max_y_coord > min_y_coord else 1.0

    ax.set_xlim(min_x_coord - 0.05 * x_range, max_x_coord + 0.05 * x_range)
    y_padding_factor = 0.15 if num_leaves > 10 else 0.05
    ax.set_ylim(min_y_coord - y_padding_factor * y_range, max_y_coord + 0.05 * y_range)

    branch_colors = ["#E63946", "#457B9D", "#2A9D8F", "#F4A261", "#A8DADC", "#1D3557", "#6A0DAD", "#FFC300"]
    color_idx_counter = 0

    memo_drawn_branches = set()

    def get_branch_color():
        nonlocal color_idx_counter
        color = branch_colors[color_idx_counter % len(branch_colors)]
        color_idx_counter += 1
        return color

    def draw_connections_recursive(node, parent_color=None):
        if node not in node_positions: return

        px, py = node_positions[node]

        children_branch_color = get_branch_color() if not node.is_leaf else parent_color

        children_to_draw = []
        if node.left and node.left in node_positions: children_to_draw.append(node.left)
        if node.right and node.right in node_positions: children_to_draw.append(node.right)

        child_x_coords = []
        for child_node in children_to_draw:
            cx, cy = node_positions[child_node]
            child_x_coords.append(cx)

            branch_id = tuple(sorted((id(node), id(child_node))))
            if branch_id not in memo_drawn_branches:
                ax.add_line(Line2D([cx, cx], [cy, py], color=children_branch_color, lw=1.5))
                memo_drawn_branches.add(branch_id)

            draw_connections_recursive(child_node, children_branch_color)

        if len(child_x_coords) > 0:  # If there are children
            min_cx = min(child_x_coords) if child_x_coords else px
            max_cx = max(child_x_coords) if child_x_coords else px


            if min_cx != max_cx:
                ax.add_line(Line2D([min_cx, max_cx], [py, py], color=children_branch_color, lw=1.5))


        if node.is_leaf:
            rotation = 0
            ha = 'center'
            va = 'top'
            label_y_offset = -0.015 * y_range

            if num_leaves > 8:
                rotation = 45
                ha = 'right'
                va = 'top'
                if num_leaves > 15:
                    rotation = 60

            ax.text(px, py + label_y_offset, node.name,
                    rotation=rotation,
                    ha=ha,
                    va=va,
                    fontsize=8,
                    rotation_mode="anchor")


    color_idx_counter = 0
    draw_connections_recursive(root)

    plt.title("Phylogenetic Tree (UPGMA)")


    plt.savefig(filename, bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"Tree plot saved to {filename}")



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
                current_name = line[1:]
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

    sequence_names = lines[0].strip().split()
    n = len(sequence_names)

    distance_matrix = np.zeros((n, n))

    for i in range(n):
        values = lines[i+1].strip().split()
        for j in range(n):
            distance_matrix[i, j] = float(values[j])

    return distance_matrix, sequence_names

def main_upgma(sequences=None, distance_matrix=None, sequence_names=None,
               output_file="phylogenetic_tree.txt", output_image="phylogenetic_tree.png"):
    """
    Main function to run the UPGMA algorithm.

    Args:
        sequences (dict, optional): A dictionary mapping sequence names to sequences
        distance_matrix (numpy.ndarray, optional): The distance matrix
        sequence_names (list, optional): The names of the sequences
        output_file (str): The name of the file to save the tree to
        output_image (str): The name of the file to save the tree plot to

    Returns:
        Node: The root node of the phylogenetic tree
    """
    if distance_matrix is None and sequences is not None:
        distance_matrix, sequence_names = calculate_distance_matrix(
            sequences)

    root = upgma(distance_matrix, sequence_names)

    save_tree_to_file(root, output_file, distance_matrix, sequence_names)
    plot_tree(root, output_image)

    return root