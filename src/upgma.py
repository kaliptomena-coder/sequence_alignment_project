import numpy as np

def find_lowest_cell(matrix):
    """Finding the coordinates of the minimum distance in the matrix."""
    min_dist = float('inf')
    x, y = -1, -1
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if i != j and matrix[i][j] < min_dist:
                min_dist = matrix[i][j]
                x, y = i, j
    return x, y

def update_labels(labels, i, j):
    """Combining labels of merged sequences into a nested tuple format."""
    # Ensuring the lower index is handled first for consistency
    if j < i:
        i, j = j, i

    # Creating a new cluster label
    new_label = (labels[i], labels[j])
    labels[i] = new_label
    del labels[j]
    return labels

def run_upgma(labels, matrix):
    """Executing the UPGMA algorithm to build a guide tree."""
    # Converting matrix to a list of lists for easier manipulation
    curr_matrix = [list(row) for row in matrix]
    curr_labels = list(labels)

    print("Starting UPGMA clustering...")

    while len(curr_labels) > 1:
        # 1. Locating the closest pair
        i, j = find_lowest_cell(curr_matrix)

        print(f"Merging: {curr_labels[i]} and {curr_labels[j]}")

        # 2. Calculating distances for the new cluster
        new_row = []
        for k in range(len(curr_matrix)):
            if k != i and k != j:
                # Computing arithmetic mean distance
                dist = (curr_matrix[i][k] + curr_matrix[j][k]) / 2
                new_row.append(dist)

        # 3. Rebuilding the matrix
        # Removing old rows/columns (removing higher index first to avoid shifts)
        first, second = sorted([i, j], reverse=True)
        for row in curr_matrix:
            del row[first]
            del row[second]
        del curr_matrix[first]
        del curr_matrix[second]

        # Adding the new cluster's distances to the matrix
        for idx, row in enumerate(curr_matrix):
            row.append(new_row[idx])
        new_row.append(0.0) # Distance to itself
        curr_matrix.append(new_row)

        # 4. Updating the label list
        curr_labels = update_labels(curr_labels, i, j)

    return curr_labels[0]

if __name__ == "__main__":
    # Testing with your actual Hemoglobin results
    example_labels = ["HBA_HUMAN", "HBB_HUMAN", "MYG_HUMAN"]
    example_matrix = [
        [0.0, 0.5676, 0.7161],
        [0.5676, 0.0, 0.7451],
        [0.7161, 0.7451, 0.0]
    ]

    guide_tree = run_upgma(example_labels, example_matrix)
    print("\nFinal Guide Tree (Cluster Order):")
    print(guide_tree)