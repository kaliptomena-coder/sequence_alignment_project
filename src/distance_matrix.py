import sys
import os

# Getting the absolute path of the directory containing this script
current_dir = os.path.dirname(os.path.abspath(__file__))
# Adding that directory to the Python search path
if current_dir not in sys.path:
    sys.path.append(current_dir)

# Importing your Needleman-Wunsch implementation
# Check: Ensure your file is named exactly 'needlemanWunschGlobal.py'
try:
    from needlemanWunschGlobal import needleman_wunsch
except ImportError:
    # If you renamed the file to have two 'l's (needlemanWunschGllobal), this catch handles it
    from needlemanWunschGllobal import needleman_wunsch

def calculate_distance(seq1, seq2):
    """Calculating the evolutionary distance by finding identity."""
    # Running your global alignment
    # If your NW returns (aln1, aln2, score), result[0] and result[1] are the strings.
    result = needleman_wunsch(seq1, seq2)

    # If the result is a tuple, we take the strings.
    # If the result is just the matrix (a list), we use the original sequences as a fallback.
    if isinstance(result, tuple):
        aln1, aln2 = result[0], result[1]
    else:
        # This fallback allows the code to run even if the NW return is still the matrix
        aln1, aln2 = seq1, seq2

        # COUNTING MATCHES
    matches = sum(1 for a, b in zip(aln1, aln2) if a == b)
    length = max(len(aln1), len(aln2))

    identity = matches / length if length > 0 else 0
    return round(1.0 - identity, 4)

def generate_matrix(sequence_dict):
    """Generating a symmetric distance matrix for a dictionary of sequences."""
    names = list(sequence_dict.keys())
    n = len(names)
    matrix = [[0.0 for _ in range(n)] for _ in range(n)]

    print(f"Calculating distances for {n} sequences using Needleman-Wunsch")

    for i in range(n):
        for j in range(i + 1, n):
            dist = calculate_distance(sequence_dict[names[i]], sequence_dict[names[j]])
            matrix[i][j] = dist
            matrix[j][i] = dist # Matrix is symmetric

    return names, matrix

def print_matrix(names, matrix):
    """Formatting the matrix output for the terminal or report."""
    header = "\t" + "\t".join(names)
    print("\nDistance Matrix (Identity-based) ")
    print(header)
    for i, row in enumerate(matrix):
        row_str = "\t".join(map(str, row))
        print(f"{names[i]}\t{row_str}")

if __name__ == "__main__":
    # Test dataset: Short sample sequences
    test_data = {
        "Seq_A": "GATTACA",
        "Seq_B": "GATTAGA",
        "Seq_C": "GATCCCC",
        "Seq_D": "GATAAAA"
    }

    labels, dist_matrix = generate_matrix(test_data)
    print_matrix(labels, dist_matrix)