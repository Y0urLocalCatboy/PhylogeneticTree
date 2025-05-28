def matrix_generation(first_sentence, second_sentence):
    """
     Generate the initial matrix.

     Args:
         first_sentence (str): The first sequence to align
         second_sentence (str): The second sequence to align

     Returns:
         list: 2D matrix with sequences as headers

     Raises:
         TypeError: If either input is not a string
         ValueError: If either input is an empty string
     """
    if isinstance(first_sentence, str) == False or isinstance(second_sentence, str) == False:
        raise TypeError("One of the inputs is not a string.")
    if len(first_sentence) == 0 or len(second_sentence) == 0:
        raise ValueError("One of the sentences are empty.")

    first_sentence = "X_" + first_sentence
    second_sentence = "X_" + second_sentence
    matrix = [[0 for x in range(len(second_sentence))] for y in range(len(first_sentence))]
    for i in range(len(first_sentence)):
        matrix[i][0] = first_sentence[i]
    for j in range(len(second_sentence)):
        matrix[0][j] = second_sentence[j]
    return matrix

def matrix_printer(matrix, file=None):
    """
    Printer of the matrix to console with writing it to a file.

    Args:
        matrix (list): 2D matrix to print
        file (file object, optional): File to write the matrix to. If None, prints to console.
    """
    for i in range(len(matrix)):
        if file:
            file.write(str(matrix[i]) + '\n')
        else:
            print(matrix[i])
    if file:
        file.write('\n')
    else:
        print("")

def edge_filler(matrix, gap_penalty):
    """
    Filler of the edges of the matrix with gap penalties.

    Args:
        matrix (list): The initial matrix with sequence headers
        gap_penalty (int): The penalty for introducing a gap

    Returns:
        list: Matrix with edge values filled according to gap penalties
    """
    for j in range(2, len(matrix[0])):
        matrix[1][j] = gap_penalty * (j-1)
    for i in range(2, len(matrix)):
        matrix[i][1] = gap_penalty * (i-1)
    return matrix

def needleman_wunsch(first_sentence, second_sentence, match_score, mismatch_penalty, gap_penalty):
    """
    Implementation of the Needleman-Wunsch algorithm for sequence alignment.

    Args:
        first_sentence (str): The first sequence to align
        second_sentence (str): The second sequence to align
        match_score (int): Score for matching characters
        mismatch_penalty (int): Penalty for mismatched characters
        gap_penalty (int): Penalty for introducing a gap

    Returns:
        list: Completed scoring matrix for the alignment
    """
    matrix = edge_filler(matrix_generation(first_sentence, second_sentence), gap_penalty)

    matrix = edge_filler(matrix, gap_penalty)

    for i in range(2, len(matrix)):
        for j in range(2, len(matrix[0])):
            if isinstance(matrix[i][0], str) and isinstance(matrix[0][j], str):
                if matrix[i][0] == matrix[0][j]:
                    score = match_score
                else:
                    score = mismatch_penalty

                match = matrix[i - 1][j - 1] + score
                delete = matrix[i - 1][j] + gap_penalty
                insert = matrix[i][j - 1] + gap_penalty

                matrix[i][j] = max(match, delete, insert)

    return matrix


def needleman_wunsch_directions(matrix, match_score, mismatch_penalty, gap_penalty):
    """
       Function createing a directional matrix, that indicates the optimal path for alignment.

       Args:
           matrix (list): The numerical matrix from needleman_wunsch
           match_score (int): Score for matching characters
           mismatch_penalty (int): Penalty for mismatched characters
           gap_penalty (int): Penalty for introducing a gap

       Returns:
           list: Matrix with directions (UP, LEFT, DIAG, etc.) indicating alignment path
       """
    directional_matrix = [[' ' for x in range(len(matrix[0]))] for y in range(len(matrix))]

    for i in range(len(matrix)):
        directional_matrix[i][0] = matrix[i][0]
    for j in range(len(matrix[0])):
        directional_matrix[0][j] = matrix[0][j]

    for j in range(2, len(matrix[0])):
        directional_matrix[1][j] = "LEFT"
    for i in range(2, len(matrix)):
        directional_matrix[i][1] = "UP"

    directional_matrix[1][1] = 0

    for i in range(2, len(matrix)):
        for j in range(2, len(matrix[0])):
            diagonal = matrix[i - 1][j - 1]
            up = matrix[i - 1][j]
            left = matrix[i][j - 1]

            if matrix[i][0] == matrix[0][j]:
                diag_score = diagonal + match_score
            else:
                diag_score = diagonal + mismatch_penalty

            up_score = up + gap_penalty
            left_score = left + gap_penalty

            max_score = max(diag_score, up_score, left_score)
            if diag_score == left_score == up_score:
                directional_matrix[i][j] = "DIAG_UP_LEFT"
            elif diag_score == up_score and diag_score > left_score:
                directional_matrix[i][j] = "DIAG_UP"
            elif diag_score == left_score and diag_score > up_score:
                directional_matrix[i][j] = "DIAG_LEFT"
            elif up_score == left_score and up_score > diag_score:
                directional_matrix[i][j] = "UP_LEFT"
            elif max_score == diag_score:
                directional_matrix[i][j] = "DIAG"
            elif max_score == up_score:
                directional_matrix[i][j] = "UP"
            else:
                directional_matrix[i][j] = "LEFT"


    return directional_matrix


def traceback(numeric_matrix, directional_matrix, first_sentence, second_sentence):
    """
    Performs the traceback to find the optimal alignment path.

    Args:
        numeric_matrix (list): The score matrix from needleman_wunsch
        directional_matrix (list): The direction matrix from needleman_wunsch_directions
        first_sentence (str): The first sequence
        second_sentence (str): The second sequence

    Returns:
        tuple: (path, aligned_seq1, aligned_seq2)
            - path (list): The sequence of directions taken
            - aligned_seq1 (str): The first sequence with gaps inserted
            - aligned_seq2 (str): The second sequence with gaps inserted
    """
    path = []
    alignment1 = []
    alignment2 = []

    i = len(directional_matrix) - 1
    j = len(directional_matrix[0]) - 1
    biggest_number = numeric_matrix[i][j]
    for counter in range(len(directional_matrix[0]) - 1, 1, -1):
        if numeric_matrix[i][counter] > biggest_number:
            biggest_number = numeric_matrix[i][j]
            j = counter
    while i > 1 or j > 1:
        current_direction = directional_matrix[i][j]

        if current_direction == 0:
            i -= 1
            j -= 1
            continue

        path.append(current_direction)

        if i <= 1 < j:
            alignment1.append('-')
            alignment2.append(second_sentence[j-2])
            j -= 1
            continue
        elif j <= 1 < i:
            alignment1.append(first_sentence[i-2])
            alignment2.append('-')
            i -= 1
            continue

        if current_direction == "DIAG" or current_direction == "DIAG_UP" or current_direction == "DIAG_LEFT" or current_direction == "DIAG_UP_LEFT":
            alignment1.append(first_sentence[i-2])
            alignment2.append(second_sentence[j-2])
            i -= 1
            j -= 1
        elif current_direction == "UP" or current_direction == "UP_LEFT":
            alignment1.append(first_sentence[i-2])
            alignment2.append('-')
            i -= 1
        elif current_direction == "LEFT":
            alignment1.append('-')
            alignment2.append(second_sentence[j-2])
            j -= 1

    path.reverse()
    alignment1.reverse()
    alignment2.reverse()

    return path, ''.join(alignment1), ''.join(alignment2)

def similarity_score(first_allingment, second_allingment):
    """
    Calculate the similarity score between two aligned sequences.

    Args:
        first_allingment (str): The first aligned sequence
        second_allingment (str): The second aligned sequence

    Returns:
        float: Similarity score as a fraction of matching positions
    """
    score = 0
    for i in range(len(first_allingment)):
        if first_allingment[i] == second_allingment[i]:
            score += 1
    return score/len(first_allingment)

def comparision_analyzer(first_allingment, second_allingment):
    """
    Create a string representation of the alignment comparison.

    Args:
        first_allingment (str): The first aligned sequence
        second_allingment (str): The second aligned sequence

    Returns:
        str: String where:
            - '*' represents a match
            - '-' represents a gap
            - '!' represents a mismatch
    """
    comparision = ""
    for i in range(len(first_allingment)):
        if first_allingment[i] == second_allingment[i]:
            comparision = comparision + "*"
        elif first_allingment[i] == "-" or second_allingment[i] == "-":
            comparision = comparision + "-"
        else:
            comparision = comparision + "!"
    return comparision


def mark_optimal_path(matrix, path_coordinates):
    """
    Marks the optimal path in a matrix with an indicator.

    Args:
        matrix (list): The matrix to mark the path in
        path_coordinates (list): List of (i, j) coordinates in the optimal path

    Returns:
        list: A copy of the matrix with the optimal path marked
    """
    marked_matrix = [row[:] for row in matrix]

    for i, j in path_coordinates:
        if isinstance(marked_matrix[i][j], str) and len(marked_matrix[i][j]) > 1:
            marked_matrix[i][j] = f"[{marked_matrix[i][j]}]"
        else:
            marked_matrix[i][j] = f"[{marked_matrix[i][j]}]"

    return marked_matrix


def main(first_sentence, second_sentence, match_score, mismatch_penalty, gap_penalty):
    """
    Executes the Needleman-Wunsch algorithm and saves the results to a file.

    Args:
        first_sentence (str): The first sequence to align
        second_sentence (str): The second sequence to align
        match_score (int): Score for matching characters
        mismatch_penalty (int): Penalty for mismatched characters
        gap_penalty (int): Penalty for introducing a gap

    Effects:
        Creates/overwrites 'alignment_results.txt' with the alignment results
    """
    with open("alignment_results.txt", "w") as output_file:
        numeric_matrix = needleman_wunsch(first_sentence, second_sentence, match_score, mismatch_penalty, gap_penalty)
        directional_matrix = needleman_wunsch_directions(numeric_matrix, match_score, mismatch_penalty, gap_penalty)

        path, aligned_seq1, aligned_seq2 = traceback(numeric_matrix, directional_matrix, first_sentence,
                                                     second_sentence)

        path_coordinates = []
        i = len(directional_matrix) - 1
        j = len(directional_matrix[0]) - 1
        path_coordinates.append((i, j))

        for direction in path:
            if "DIAG" in direction:
                i -= 1
                j -= 1
            elif direction == "UP" or "UP_LEFT" in direction:
                i -= 1
            elif direction == "LEFT":
                j -= 1
            path_coordinates.append((i, j))

        if (1, 1) not in path_coordinates:
            path_coordinates.append((1, 1))

        marked_numeric_matrix = mark_optimal_path(numeric_matrix, path_coordinates)
        marked_directional_matrix = mark_optimal_path(directional_matrix, path_coordinates)

        output_file.write("Numeric Matrix with Optimal Path Marked in:\n")
        matrix_printer(marked_numeric_matrix, output_file)

        output_file.write("Directional Matrix with Optimal Path Marked in:\n")
        matrix_printer(marked_directional_matrix, output_file)

        output_file.write(f"Missmatch = {mismatch_penalty}\n")
        output_file.write(f"Gap = {gap_penalty}\n")
        output_file.write(f"Match = {match_score}\n")
        output_file.write(f"Path: {path}\n")
        output_file.write(f"1: {aligned_seq1}\n")
        output_file.write(f"2: {aligned_seq2}\n")
        output_file.write(f"Match score: {(similarity_score(aligned_seq1, aligned_seq2) * 100).__round__(2)}%\n")
        output_file.write("* -> match\n")
        output_file.write("- -> gap\n")
        output_file.write("! -> missmatch\n")
        output_file.write(f"Comparision: {comparision_analyzer(aligned_seq1, aligned_seq2)}\n")

        print(f"Results saved to the file")