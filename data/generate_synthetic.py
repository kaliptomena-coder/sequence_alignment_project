# =============================================================================
# Generate synthetic sequence datasets with known ground-truth edits:
#   data/synthetic_snps.fasta      – sequences with known SNPs
#   data/synthetic_indels.fasta    – sequences with insertions/deletions
#   data/synthetic_shuffled.fasta  – shuffled sequences (negative control)
#   data/ground_truth.txt          – exact list of every edit made
#

import random
import os

# Fix the random seed so results are reproducible
# (same seed = same sequences every time you run the script)
random.seed(42)

# Output directory
DATA_DIR = os.path.dirname(os.path.abspath(__file__))


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def random_dna(length):
    """Generate a random DNA sequence of the given length."""
    return ''.join(random.choices("ACGT", k=length))


def introduce_snps(sequence, n_snps):
    """
    Introduce exactly n_snps point substitutions (SNPs) at random positions.

    A SNP (Single Nucleotide Polymorphism) is a single-character change,
    e.g., A → T at position 42.

    Returns
    -------
    mutated_seq : str  – the sequence with SNPs applied
    edit_log    : list of str – description of each edit for ground truth
    """
    seq      = list(sequence)
    edit_log = []

    # Choose n_snps distinct positions
    positions = sorted(random.sample(range(len(seq)), n_snps))

    for pos in positions:
        original_base = seq[pos]
        # Choose a different base (not the same as original)
        new_base = random.choice([b for b in "ACGT" if b != original_base])
        seq[pos] = new_base
        edit_log.append(f"SNP  pos={pos}  {original_base}→{new_base}")

    return ''.join(seq), edit_log


def introduce_indels(sequence, n_insertions, n_deletions):
    """
    Introduce insertions and deletions (indels) into a sequence.

    An INSERTION adds a new character (simulates an extra base pair).
    A DELETION removes a character (simulates a lost base pair).

    Returns
    -------
    mutated_seq : str  – the sequence with indels applied
    edit_log    : list of str – description of each edit
    """
    seq      = list(sequence)
    edit_log = []

    # --- INSERTIONS ---
    # We apply insertions in reverse position order so earlier insertions
    # don't shift the positions of later ones.
    insert_positions = sorted(
        random.sample(range(len(seq)), n_insertions),
        reverse=True
    )
    for pos in insert_positions:
        new_base = random.choice("ACGT")
        seq.insert(pos, new_base)
        edit_log.append(f"INS  pos={pos}  inserted={new_base}")

    # --- DELETIONS ---
    # Similarly, apply in reverse to avoid position shifting.
    # after insertions, sequence is longer; sample from new length.
    delete_positions = sorted(
        random.sample(range(len(seq)), n_deletions),
        reverse=True
    )
    for pos in delete_positions:
        deleted_base = seq[pos]
        del seq[pos]
        edit_log.append(f"DEL  pos={pos}  deleted={deleted_base}")

    return ''.join(seq), edit_log


def shuffle_sequence(sequence):
    """
    Randomly permute all characters in a sequence.
    Used as a NEGATIVE CONTROL — a shuffled sequence has the same base
    composition but no meaningful similarity to the original.
    """
    seq = list(sequence)
    random.shuffle(seq)
    return ''.join(seq)


def write_fasta(filepath, sequences):
    """
    Write a dictionary of {name: sequence} to a FASTA file.

    FASTA format:
        >sequence_name
        ACGTACGTACGT...
    """
    with open(filepath, 'w') as f:
        for name, seq in sequences.items():
            f.write(f">{name}\n")
            # Write sequence in lines of 60 characters (standard FASTA format)
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")
    print(f"  Written: {filepath}  ({len(sequences)} sequences)")


def write_ground_truth(filepath, records):
    """
    Write the ground-truth edit log to a text file.

    Each record is a dict with keys: name, description, edits
    """
    with open(filepath, 'w') as f:
        f.write("GROUND TRUTH EDIT LOG\n")
        f.write("=" * 60 + "\n")
        f.write("This file describes every synthetic mutation made.\n")
        f.write("Use this to verify your alignment algorithms found the edits.\n\n")

        for record in records:
            f.write(f"\n[{record['name']}]\n")
            f.write(f"  Template: {record['template']}\n")
            f.write(f"  {record['description']}\n")
            for edit in record['edits']:
                f.write(f"  • {edit}\n")

    print(f"  Written: {filepath}")


# =============================================================================
# MAIN DATASET GENERATION
# =============================================================================

def generate_all_datasets():
    """Generate all four synthetic datasets."""

    print("\n" + "=" * 60)
    print("SYNTHETIC DATASET GENERATOR")
    print("=" * 60)

    all_ground_truth = []   # collect all edits for the ground truth file

    # -------------------------------------------------------------------------
    # DATASET 1: Sequences with known SNPs
    # -------------------------------------------------------------------------
    print("\n[1] Generating SNP dataset (10 SNPs per sequence)...")

    # The template is a 2000 bp synthetic "gene"
    template_snp = random_dna(2000)
    snp_sequences = {"TEMPLATE": template_snp}

    for i in range(1, 6):    # 5 derived sequences
        mutated, edits = introduce_snps(template_snp, n_snps=10)
        name = f"SNP_SEQ_{i:02d}"
        snp_sequences[name] = mutated
        all_ground_truth.append({
            "name": name,
            "template": "TEMPLATE",
            "description": f"10 SNPs introduced into TEMPLATE",
            "edits": edits
        })

    write_fasta(
        os.path.join(DATA_DIR, "synthetic_snps.fasta"),
        snp_sequences
    )

    # -------------------------------------------------------------------------
    # DATASET 2: Sequences with insertions and deletions (indels)
    # -------------------------------------------------------------------------
    print("[2] Generating INDEL dataset (5 insertions + 5 deletions)...")

    template_indel = random_dna(500)
    indel_sequences = {"TEMPLATE_INDEL": template_indel}

    for i in range(1, 6):
        mutated, edits = introduce_indels(template_indel,
                                          n_insertions=5, n_deletions=5)
        name = f"INDEL_SEQ_{i:02d}"
        indel_sequences[name] = mutated
        all_ground_truth.append({
            "name": name,
            "template": "TEMPLATE_INDEL",
            "description": "5 insertions + 5 deletions introduced",
            "edits": edits
        })

    write_fasta(
        os.path.join(DATA_DIR, "synthetic_indels.fasta"),
        indel_sequences
    )

    # -------------------------------------------------------------------------
    # DATASET 3: Shuffled sequences (negative control)
    # -------------------------------------------------------------------------
    print("[3] Generating SHUFFLED dataset (negative control)...")

    template_shuffle = random_dna(300)
    shuffle_sequences = {"ORIGINAL": template_shuffle}

    for i in range(1, 4):
        shuffled = shuffle_sequence(template_shuffle)
        name     = f"SHUFFLED_{i:02d}"
        shuffle_sequences[name] = shuffled
        all_ground_truth.append({
            "name": name,
            "template": "ORIGINAL",
            "description": "Complete random shuffle (negative control). "
                           "Alignment should find NO meaningful similarity.",
            "edits": ["All positions permuted randomly"]
        })

    write_fasta(
        os.path.join(DATA_DIR, "synthetic_shuffled.fasta"),
        shuffle_sequences
    )

    # -------------------------------------------------------------------------
    # DATASET 4: Protein-like sequences with known phylogeny
    # -------------------------------------------------------------------------
    print("[4] Generating PROTEIN synthetic dataset (20 amino acid alphabet)...")

    AA = "ACDEFGHIKLMNPQRSTVWY"   # 20 standard amino acids

    # Ancestor sequence
    ancestor = ''.join(random.choices(AA, k=200))
    protein_sequences = {"ANCESTOR": ancestor}

    # Two "lineages" — each diverges slightly differently
    for lineage, divergence in [("LINEAGE_A", 10), ("LINEAGE_B", 25)]:
        parent = ancestor
        for gen in range(1, 4):   # 3 generations
            mutated, edits = introduce_snps(parent, n_snps=divergence)
            # Occasionally also add an indel
            if gen == 2:
                mutated, indel_edits = introduce_indels(mutated,
                                                        n_insertions=2,
                                                        n_deletions=2)
                edits += indel_edits
            name = f"{lineage}_GEN{gen}"
            protein_sequences[name] = mutated
            all_ground_truth.append({
                "name": name,
                "template": f"ANCESTOR (via {lineage})",
                "description": f"Generation {gen}, {divergence} SNPs/gen",
                "edits": edits
            })
            parent = mutated   # next generation derives from this one

    write_fasta(
        os.path.join(DATA_DIR, "synthetic_protein_family.fasta"),
        protein_sequences
    )

    # -------------------------------------------------------------------------
    # Write the ground truth file
    # -------------------------------------------------------------------------
    print("[5] Writing ground truth log...")
    write_ground_truth(
        os.path.join(DATA_DIR, "ground_truth.txt"),
        all_ground_truth
    )

    print("\n" + "=" * 60)
    print("DONE! Files created:")
    print("  data/synthetic_snps.fasta")
    print("  data/synthetic_indels.fasta")
    print("  data/synthetic_shuffled.fasta")
    print("  data/synthetic_protein_family.fasta")
    print("  data/ground_truth.txt")
    print("\nNow you can load these with your data_loader.py:")
    print("  from data_loader import load_fasta")
    print("  seqs = load_fasta('synthetic_snps.fasta')")
    print("=" * 60)


if __name__ == "__main__":
    generate_all_datasets()
