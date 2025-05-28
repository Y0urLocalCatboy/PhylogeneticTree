# Phylogenetic Analysis Tools

This application provides tools for phylogenetic analysis, including:
1. Needleman-Wunsch pairwise sequence alignment
2. UPGMA (Unweighted Pair Group Method with Arithmetic Mean) phylogenetic tree construction

## Installation

### Requirements
- Python 3.6 or higher
- Required Python packages:
  - numpy
  - matplotlib
  - tkinter (usually comes with Python)

### Setup
1. Clone or download this repository
2. Install the required packages:
   ```
   pip install numpy matplotlib
   ```
3. Run the main application:
   ```
   python main.py
   ```

## Features

### Needleman-Wunsch Alignment
- Global pairwise sequence alignment
- Customizable scoring parameters (match score, mismatch penalty, gap penalty)
- Load sequences from FASTA files
- Save alignment results to a text file

### UPGMA Phylogenetic Tree Construction
- Build phylogenetic trees using the UPGMA method
- Two input options:
  - Multiple sequences (uses Needleman-Wunsch for distance calculation)
  - Pre-calculated distance matrix
- Visualization of the phylogenetic tree
- Save tree in text format (including Newick format) and as an image

## Usage

### Main Launcher
Run `main.py` to open the launcher, which allows you to choose between the two tools.

### Needleman-Wunsch Alignment
1. Enter two sequences directly or load them from FASTA files
2. Set the scoring parameters:
   - Match score: Score for matching characters (positive value)
   - Mismatch penalty: Penalty for mismatched characters (negative value)
   - Gap penalty: Penalty for introducing a gap (negative value)
3. Click "Calculate" to perform the alignment
4. Results are saved to "alignment_results.txt"

### UPGMA Phylogenetic Tree
1. Choose input type:
   - Sequences: Load multiple sequences from FASTA files
   - Distance Matrix: Load a pre-calculated distance matrix
2. If using sequences, set the alignment parameters for distance calculation
3. Set output file names for the text and image outputs
4. Click "Generate Phylogenetic Tree" to build the tree
5. The tree will be saved in both text and image formats, and the image will be displayed

## Input Formats

### FASTA Format
Sequences should be in standard FASTA format:
```
>Sequence_Name
ACGTACGTACGT
```

### Distance Matrix Format
The distance matrix should be in a text file with the following format:
- First line: sequence names separated by spaces
- Subsequent lines: distance values separated by spaces

Example:
```
Seq1 Seq2 Seq3
0.0 0.1 0.2
0.1 0.0 0.3
0.2 0.3 0.0
```

## Algorithms

### Needleman-Wunsch
The Needleman-Wunsch algorithm is a dynamic programming algorithm for global sequence alignment. It works by:
1. Creating a scoring matrix
2. Filling the matrix using the scoring parameters
3. Tracing back through the matrix to find the optimal alignment

### UPGMA
UPGMA (Unweighted Pair Group Method with Arithmetic Mean) is a hierarchical clustering method used to construct phylogenetic trees. It works by:
1. Starting with a distance matrix between all pairs of sequences
2. Iteratively joining the closest clusters
3. Updating the distance matrix after each join
4. Building a tree where branch lengths represent evolutionary distances

## Example

### Sample Sequences
```
>Human
ATGCGTACGTACGTAGC
>Mouse
ATGCGAACGTACGCAGC
>Rat
ATGCGAACGTACGCAGC
>Chicken
ATGCTTACGTACGTAGC
```

### Expected Output
The UPGMA algorithm will generate a tree showing the evolutionary relationships between these sequences, with Mouse and Rat likely clustering together first, then with Human, and Chicken as an outgroup.

## License
[Specify your license here]

## Author
[Your Name]

## Acknowledgments
- The Needleman-Wunsch implementation is based on the algorithm described by Needleman and Wunsch (1970)
- The UPGMA implementation follows the standard algorithm as described in phylogenetic analysis literature