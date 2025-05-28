import os
from upgma import load_sequences_from_fasta, load_distance_matrix_from_file, main_upgma

def test_upgma_with_sequences():
    """Test UPGMA with sequences from a FASTA file"""
    print("Testing UPGMA with sequences from FASTA file...")
    
    # Load sequences from FASTA file
    sequences = load_sequences_from_fasta("samples/sample_sequences.fasta")
    
    # Run UPGMA
    root = main_upgma(
        sequences=sequences,
        match_score=1,
        mismatch_penalty=-1,
        gap_penalty=-1,
        output_file="test_sequences_tree.txt",
        output_image="test_sequences_tree.png"
    )
    
    print(f"Tree generated successfully. Root node: {root.name}")
    print(f"Results saved to test_sequences_tree.txt and test_sequences_tree.png")
    print()

def test_upgma_with_distance_matrix():
    """Test UPGMA with a pre-calculated distance matrix"""
    print("Testing UPGMA with pre-calculated distance matrix...")
    
    # Load distance matrix from file
    distance_matrix, sequence_names = load_distance_matrix_from_file("samples/sample_distance_matrix.txt")
    
    # Run UPGMA
    root = main_upgma(
        distance_matrix=distance_matrix,
        sequence_names=sequence_names,
        output_file="test_matrix_tree.txt",
        output_image="test_matrix_tree.png"
    )
    
    print(f"Tree generated successfully. Root node: {root.name}")
    print(f"Results saved to test_matrix_tree.txt and test_matrix_tree.png")
    print()

if __name__ == "__main__":
    # Create samples directory if it doesn't exist
    if not os.path.exists("samples"):
        os.makedirs("samples")
    
    # Run tests
    test_upgma_with_sequences()
    test_upgma_with_distance_matrix()
    
    print("All tests completed successfully!")