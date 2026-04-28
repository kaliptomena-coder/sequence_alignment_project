import os

def load_fasta(filename):
    """
    Parsing a FASTA file and returning a dictionary of {sequence_id: sequence}.
    The file is expected to live in the 'data/' folder at the project root.
    """
    sequences = {}

    # Building the path to the data folder relative to this script's location
    current_dir  = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(current_dir)
    file_path    = os.path.join(project_root, 'data', filename)

    if not os.path.exists(file_path):
        print(f"Error: File not found at {file_path}")
        return {}

    with open(file_path, 'r') as f:
        current_label = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Each '>' line starts a new record; we use the first word as the ID
                current_label = line[1:].split()[0]
                sequences[current_label] = ""
            elif current_label:
                sequences[current_label] += line.replace(" ", "")

    return sequences


if __name__ == "__main__":
    # Quick test to make sure the loader works
    data = load_fasta("globins.fasta")
    if data:
        print("Loaded sequences:")
        for name, seq in data.items():
            print(f"  {name}: {len(seq)} residues")
    else:
        print("No data loaded — check that globins.fasta exists in the data/ folder.")
