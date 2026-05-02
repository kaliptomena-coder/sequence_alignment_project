# This script measures and compares:
#   - Runtime vs. sequence length for pairwise methods (NW, SW, Gotoh, Hirschberg, Banded)
#   - Memory usage: NW (O(nm)) vs. Hirschberg (O(n+m))
#   - BLAST-lite vs. NW runtime at increasing sequence lengths
#   - NW gap-penalty parameter sensitivity
#   - Biological dataset pairwise alignment scores (if FASTA files are present)
#
# Results are saved as CSV files in benchmarks/results/.
# Plots are generated and saved as PNG files (requires matplotlib).


import itertools
import sys
import os
import time
import random
import csv
import tracemalloc

# Resolving paths relative to this script
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC_DIR      = os.path.join(PROJECT_ROOT, 'src')
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from needlemanWunschGlobal import needleman_wunsch
from smithWatermanLocal    import smith_waterman
from gotoh                  import gotoh
from hirschberg            import hirschberg
from banded_alignment      import banded_nw
from blast_lite            import blast_lite

RESULTS_DIR = os.path.join(os.path.dirname(__file__), 'results')
os.makedirs(RESULTS_DIR, exist_ok=True)


# Utility: generating and mutating random sequences

def random_sequence(length, alphabet="ACGT"):
    """Generating a random DNA sequence of the given length."""
    return ''.join(random.choices(alphabet, k=length))


def mutate_sequence(seq, num_mutations):
    """
    Creating a related sequence by introducing point substitutions.
    Used to generate pairs that are similar but not identical.
    """
    seq      = list(seq)
    alphabet = "ACGT"
    for pos in random.sample(range(len(seq)), min(num_mutations, len(seq))):
        seq[pos] = random.choice([c for c in alphabet if c != seq[pos]])
    return ''.join(seq)

# Utility: timing and memory measurement

def time_function(func, *args, repeats=5, **kwargs):
    """
    Running a function multiple times and returning the median elapsed time in seconds.
    Using the median avoids distortion from occasional OS scheduling delays.
    """

    times = []
    for _ in range(repeats):
        start = time.perf_counter()
        func(*args, **kwargs)
        times.append(time.perf_counter() - start)
    times.sort()
    return times[len(times) // 2]


def measure_memory(func, *args):
    """
    Running a function once and returning peak memory usage in kilobytes.
    Uses Python's built-in tracemalloc module to track allocations.
    """
    tracemalloc.start()
    func(*args)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return peak / 1024  # converting bytes to KB

# Benchmark 1: Pairwise runtime vs. sequence length

def benchmark_pairwise_runtime():
    """
    Testing NW, SW, Gotoh, Hirschberg, and Banded NW on sequences of increasing length.
    Each data point is the median of 3 runs.
    """
    lengths    = [50, 100, 200, 300, 500, 750, 1000]
    algorithms = {
        "Needleman-Wunsch": lambda s1, s2: needleman_wunsch(s1, s2),
        "Smith-Waterman"  : lambda s1, s2: smith_waterman(s1, s2),
        "Gotoh"           : lambda s1, s2: gotoh(s1, s2),
        "Hirschberg"      : lambda s1, s2: hirschberg(s1, s2),
        "Banded (k=10)"   : lambda s1, s2: banded_nw(s1, s2, k=10),
    }
    results = {name: {} for name in algorithms}


    print("PAIRWISE RUNTIME BENCHMARK\n")

    print(f"{'Length':>8} | " + " | ".join(f"{n:>18}" for n in algorithms))


    for length in lengths:
        s1  = random_sequence(length)
        s2  = mutate_sequence(s1, length // 10)
        row = f"{length:>8} | "
        for name, func in algorithms.items():
            t                  = time_function(func, s1, s2, repeats=3)
            results[name][length] = t
            row               += f"{t:>18.4f} | "
        print(row)

    return results, lengths

# Benchmark 2: Memory — NW vs. Hirschberg

def benchmark_memory():
    """
    Comparing peak memory usage of NW (O(nm)) and Hirschberg (O(n+m)).
    This demonstrates the practical benefit of the divide-and-conquer approach.
    """
    lengths   = [100, 200, 400, 600, 800, 1000]
    nw_mem, hb_mem = [], []


    print("\nMEMORY BENCHMARK: NW vs Hirschberg\n")

    print(f"{'Length':>8} | {'NW (KB)':>15} | {'Hirschberg (KB)':>18}")


    for length in lengths:
        s1 = random_sequence(length)
        s2 = mutate_sequence(s1, length // 10)

        mem_nw = measure_memory(needleman_wunsch, s1, s2)
        mem_hb = measure_memory(hirschberg, s1, s2)

        nw_mem.append(mem_nw)
        hb_mem.append(mem_hb)

        print(f"{length:>8} | {mem_nw:>15.2f} | {mem_hb:>18.2f}")

    return lengths, nw_mem, hb_mem


# Benchmark 3: BLAST-lite vs. NW runtime

def benchmark_blast():
    """
    Comparing BLAST-lite (heuristic) against NW (exact) as sequence length grows.
    Expected outcome: BLAST scales roughly linearly, NW quadratically.
    """
    lengths     = [100, 250, 500, 750, 1000, 1500, 2000]
    blast_times = []
    nw_times    = []


    print("\nBLAST vs NW RUNTIME COMPARISON\n")

    print(f"{'Length':>8} | {'BLAST (s)':>12} | {'NW (s)':>12}\n")


    for length in lengths:
        s1 = random_sequence(length)
        s2 = mutate_sequence(s1, length // 5)

        t_blast = time_function(blast_lite, s1, s2, k=4, repeats=3)
        t_nw    = time_function(needleman_wunsch, s1, s2, repeats=3)

        blast_times.append(t_blast)
        nw_times.append(t_nw)

        print(f"{length:>8} | {t_blast:>12.4f} | {t_nw:>12.4f}")

    return lengths, blast_times, nw_times

# Benchmark 4: NW gap-penalty parameter sensitivity

def benchmark_gap_sensitivity():
    """
    Running NW with a range of gap penalties on the same pair of sequences.
    Shows how the penalty affects score and the number of gaps introduced.
    """
    s1 = "GATTACAACTTGCGTATGCAG"
    s2 = "GATCCAGTTCAAATGCGTATG"

    gap_values = [-1, -2, -3, -5, -7, -10, -15, -20]
    scores, gap_counts = [], []


    print("\n NW GAP PENALTY SENSITIVITY\n")

    print(f"{'Gap Penalty':>12} | {'Score':>8} | {'# Gaps':>8}\n")


    for gap in gap_values:
        a1, a2, score = needleman_wunsch(s1, s2, gap=gap)
        n_gaps        = a1.count('-') + a2.count('-')
        scores.append(score)
        gap_counts.append(n_gaps)
        print(f"{gap:>12} | {score:>8} | {n_gaps:>8}")

    return gap_values, scores, gap_counts


# Benchmark 5: Biological dataset alignment

def benchmark_all_datasets():
    """
    Running Smith-Waterman on all sequence pairs in each biological FASTA dataset.
    Skips datasets whose files are not found.
    """
    try:
        from data_loader import load_fasta
    except ImportError:
        print("Could not import data_loader — skipping biological dataset benchmark.")
        return

    datasets = {
        'Globins'         : os.path.join(PROJECT_ROOT, 'data', 'globins.fasta'),
        'Cytochrome c'    : os.path.join(PROJECT_ROOT, 'data', 'cytochrome_c.fasta'),
        'Serine Proteases': os.path.join(PROJECT_ROOT, 'data', 'serine_proteases.fasta'),
        'Synthetic SNPs'  : os.path.join(PROJECT_ROOT, 'data', 'synthetic_snps.fasta'),
    }


    print("\nBIOLOGICAL DATASETS — PAIRWISE SW SCORES\n")


    for dataset_name, filepath in datasets.items():
        if not os.path.exists(filepath):
            print(f"  Skipping {dataset_name}: file not found ({filepath})")
            continue

        seqs  = load_fasta(filepath)
        names = list(seqs.keys())
        print(f"\n=== {dataset_name} ({len(seqs)} sequences) ===")

        for name1, name2 in itertools.combinations(names, 2):
           _, _, sw_score = smith_waterman(seqs[name1], seqs[name2])
           _, _, nw_score = needleman_wunsch(seqs[name1], seqs[name2])
           _, _, gotoh_score = gotoh(seqs[name1], seqs[name2])
           print(f"  {name1} vs {name2}  |  SW score: {sw_score} |  NW score: {nw_score} | Gotoh score: {gotoh_score}")


# CSV Saving

def save_csv(filename, headers, rows):
    """Saving benchmark results to a CSV file in the results/ directory."""
    path = os.path.join(RESULTS_DIR, filename)
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)
    print(f"  Saved: {path}")


# Plotting (requires matplotlib)

def make_plots(pairwise_results, pairwise_lengths,
               mem_lengths, nw_mem, hb_mem,
               blast_lengths, blast_times, nw_times_blast,
               gap_values, gap_scores, gap_counts):
    """
    Generating four benchmark plots and saving them as PNG files.
    Falls back gracefully if matplotlib is not installed.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib not found — skipping plots. Install with: pip install matplotlib")
        return

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

    # Plot 1: Runtime vs sequence length
    fig, ax = plt.subplots(figsize=(10, 6))
    for (name, data), color in zip(pairwise_results.items(), colors):
        lengths_list = sorted(data.keys())
        ax.plot(lengths_list, [data[l] for l in lengths_list],
                marker='o', label=name, color=color, linewidth=2)
    ax.set_xlabel("Sequence Length"); ax.set_ylabel("Median Time (s)")
    ax.set_title("Runtime vs. Sequence Length — Pairwise Methods", fontweight='bold')
    ax.legend(); ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "plot_runtime_pairwise.png"), dpi=150)
    plt.close()
    print("  Saved: plot_runtime_pairwise.png")

    # Plot 2: Memory NW vs Hirschberg
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(mem_lengths, nw_mem, marker='s', label="NW O(nm)",         color='#d62728', linewidth=2)
    ax.plot(mem_lengths, hb_mem, marker='o', label="Hirschberg O(n+m)", color='#2ca02c', linewidth=2)
    ax.set_xlabel("Sequence Length"); ax.set_ylabel("Peak Memory (KB)")
    ax.set_title("Memory: NW vs. Hirschberg", fontweight='bold')
    ax.legend(); ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "plot_memory_nw_vs_hirschberg.png"), dpi=150)
    plt.close()
    print("  Saved: plot_memory_nw_vs_hirschberg.png")

    # Plot 3: BLAST vs NW
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(blast_lengths, blast_times,    marker='D', label="BLAST-lite", color='#ff7f0e', linewidth=2)
    ax.plot(blast_lengths, nw_times_blast, marker='o', label="NW (exact)",  color='#1f77b4', linewidth=2)
    ax.set_xlabel("Sequence Length"); ax.set_ylabel("Median Time (s)")
    ax.set_title("BLAST-lite vs. Needleman-Wunsch", fontweight='bold')
    ax.legend(); ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "plot_blast_vs_nw.png"), dpi=150)
    plt.close()
    print("  Saved: plot_blast_vs_nw.png")

    # Plot 4: Gap penalty sensitivity
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    ax1.plot(gap_values, gap_scores,  marker='o', color='#1f77b4', linewidth=2)
    ax1.set_xlabel("Gap Penalty"); ax1.set_ylabel("Alignment Score")
    ax1.set_title("Score vs. Gap Penalty"); ax1.grid(True, alpha=0.3)

    ax2.plot(gap_values, gap_counts, marker='s', color='#d62728', linewidth=2)
    ax2.set_xlabel("Gap Penalty"); ax2.set_ylabel("Number of Gaps")
    ax2.set_title("Gap Count vs. Gap Penalty"); ax2.grid(True, alpha=0.3)

    plt.suptitle("NW Parameter Sensitivity", fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "plot_gap_sensitivity.png"), dpi=150)
    plt.close()
    print("  Saved: plot_gap_sensitivity.png")

# Entry Point

if __name__ == "__main__":

    print("SEQUENCE ALIGNMENT BENCHMARK SUITE")
    print(f"Results will be saved to: {RESULTS_DIR}\n")


    pairwise_results, pairwise_lengths = benchmark_pairwise_runtime()
    mem_lengths, nw_mem, hb_mem        = benchmark_memory()
    blast_lengths, blast_times, nw_t   = benchmark_blast()
    gap_vals, gap_scores, gap_counts   = benchmark_gap_sensitivity()

    print("\nSaving CSV results...")
    alg_names = list(pairwise_results.keys())
    save_csv("timing_pairwise.csv",
             ["length"] + alg_names,
             [[l] + [pairwise_results[n][l] for n in alg_names] for l in pairwise_lengths])

    save_csv("memory_nw_vs_hirschberg.csv",
             ["length", "nw_kb", "hirschberg_kb"],
             list(zip(mem_lengths, nw_mem, hb_mem)))

    print("\nGenerating plots...")
    make_plots(pairwise_results, pairwise_lengths,
               mem_lengths, nw_mem, hb_mem,
               blast_lengths, blast_times, nw_t,
               gap_vals, gap_scores, gap_counts)

    benchmark_all_datasets()

    print("\nBENCHMARK COMPLETE")
    print(f"All outputs in: {RESULTS_DIR}\n")
