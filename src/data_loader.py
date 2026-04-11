import os

def load_fasta(filename):
    """
    Parsing a FASTA file and returning a dictionary of {name: sequence}.
    """
    sequences = {}

    # Automatically finding the path to the data folder relative to this script
    # Locating the directory of the current file
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # Moving up to the project root directory
    project_root = os.path.dirname(current_dir)
    # Joining the path components to reach the target file
    file_path = os.path.join(project_root, 'data', filename)

    # Checking if the file exists before attempting to open it
    if not os.path.exists(file_path):
        print(f"Error: File {file_path} not found!")
        return {}

    # Opening the file for reading
    with open(file_path, 'r') as f:
        current_label = None
        for line in f:
            # Cleaning whitespace from the ends of the line
            line = line.strip()
            # Skipping empty lines
            if not line:
                continue

            # Identifying a new sequence header
            if line.startswith(">"):
                # Stripping the '>' character and isolating the sequence ID
                current_label = line[1:].split()[0]
                sequences[current_label] = ""
            elif current_label:
                # Appending the sequence data and removing any internal spaces
                sequences[current_label] += line.replace(" ", "")

    return sequences

if __name__ == "__main__":
    # Running a test to verify the file is being read correctly
    data = load_fasta("globins.fasta")
    if data:
        print("Data loaded:")
        for name, seq in data.items():
            # Reporting the name and the count of amino acids found
            print(f"- {name}: {len(seq)} amino acids")