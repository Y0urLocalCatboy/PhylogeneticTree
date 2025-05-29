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