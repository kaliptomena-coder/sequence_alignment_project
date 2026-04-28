# UPGMA (Unweighted Pair Group Method with Arithmetic Mean) builds a rooted
# hierarchical guide tree from a pairwise distance matrix.
# The result is used to determine the order of sequence merging in progressive MSA.
#
# Algorithm:
#   1. Find the pair of clusters with the smallest distance
#   2. Merge them into a new cluster labeled as (cluster_i, cluster_j)
#   3. Recompute distances from the new cluster to all others
#      using the arithmetic mean of the two merged clusters
#   4. Repeat until only one cluster remains
#
# Time complexity: O(n³) naive (n rounds × O(n²) distance lookup and update)


def find_lowest_cell(matrix):
    """
    Finding the (i, j) coordinates of the minimum off-diagonal value in the matrix.
    This identifies the two closest clusters to be merged next.
    """
    min_dist = float('inf')
    x, y     = -1, -1
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if i != j and matrix[i][j] < min_dist:
                min_dist = matrix[i][j]
                x, y     = i, j
    return x, y

def update_labels(labels, i, j):
    """
    Merging two cluster labels into a nested tuple and removing the higher index.
    The nested tuple structure encodes the tree topology (e.g., ('A', ('B', 'C'))).
    """
    # Processing the lower index first to avoid index-shift issues
    if j < i:
        i, j = j, i

    labels[i] = (labels[i], labels[j])
    del labels[j]
    return labels

def run_upgma(labels, matrix):
    """
    Running the UPGMA clustering algorithm to produce a guide tree.

    Parameters
    ----------
    labels : list of str   - sequence names (one per row/column)
    matrix : list of list  - symmetric pairwise distance matrix

    Returns
    -------
    The root of the guide tree as a nested tuple of labels.
    Example: (('HBA_HUMAN', 'HBB_HUMAN'), 'MYG_HUMAN')
    """
    curr_matrix = [list(row) for row in matrix]
    curr_labels = list(labels)

    print("Starting UPGMA clustering...")

    while len(curr_labels) > 1:
        # Step 1: Finding the two closest clusters
        i, j = find_lowest_cell(curr_matrix)
        print(f"  Merging: {curr_labels[i]} + {curr_labels[j]}")

        # Step 2: Computing distances from the new cluster to all others
        # Using the arithmetic mean of the two merged clusters' distances
        new_row = []
        for k in range(len(curr_matrix)):
            if k != i and k != j:
                new_row.append((curr_matrix[i][k] + curr_matrix[j][k]) / 2)

        # Step 3: Removing old rows and columns (higher index first to avoid shifts)
        first, second = sorted([i, j], reverse=True)
        for row in curr_matrix:
            del row[first]
            del row[second]
        del curr_matrix[first]
        del curr_matrix[second]

        # Appending the new cluster's distances
        for idx, row in enumerate(curr_matrix):
            row.append(new_row[idx])
        new_row.append(0.0)  # Distance to itself
        curr_matrix.append(new_row)

        # Step 4: Updating labels
        curr_labels = update_labels(curr_labels, i, j)

    return curr_labels[0]

if __name__ == "__main__":
    # Quick test using the hemoglobin example
    example_labels = ["HBA_HUMAN", "HBB_HUMAN", "MYG_HUMAN"]
    example_matrix  = [
        [0.0,    0.5676, 0.7161],
        [0.5676, 0.0,    0.7451],
        [0.7161, 0.7451, 0.0   ],
    ]

    guide_tree = run_upgma(example_labels, example_matrix)
    print(f"\nGuide tree: {guide_tree}")