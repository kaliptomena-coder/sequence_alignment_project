# =============================================================================
# Measure runtime and memory for all alignment algorithms:
#   - Generates random sequences of increasing length (50 → 2000 chars)
#   - Times each algorithm on those sequences (5 repeats, takes median)
#   - Measures peak memory usage with Python's tracemalloc module
#   - Saves plots you can paste directly into your report
# =============================================================================
import itertools
import sys
import os
import time
import random
import string
import csv
import tracemalloc   # built-in Python module for memory tracking

# ---------------------------------------------------------------------------
# Add src directory to path — adjust if your structure differs
# ---------------------------------------------------------------------------
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC_DIR = os.path.join(PROJECT_ROOT, 'src')
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

# Import your algorithm functions
from needlemanWunschGlobal import needleman_wunsch
from smithWatermanLocal    import smith_waterman
from gotohAffineGap                 import gotoh_affine_gap
from hirschberg            import hirschberg
from banded_alignment      import banded_nw
from blast_lite            import blast_lite

# ---------------------------------------------------------------------------
# Output directory for results
# ---------------------------------------------------------------------------
RESULTS_DIR = os.path.join(os.path.dirname(__file__), 'results')
os.makedirs(RESULTS_DIR, exist_ok=True)


# =============================================================================
# 1. Random sequence generator
# =============================================================================

def random_sequence(length, alphabet="ACGT"):
    """
    Generate a random DNA sequence of the given length.
    For protein sequences, change alphabet to the 20 amino acid letters.
    """
    return ''.join(random.choices(alphabet, k=length))


def mutate_sequence(seq, num_mutations):
    """
    Create a similar sequence by introducing `num_mutations` point substitutions.
    Used to create pairs that are related but not identical.
    """
    seq = list(seq)
    positions = random.sample(range(len(seq)), min(num_mutations, len(seq)))
    alphabet = "ACGT"
    for pos in positions:
        # Replace with a different character
        choices = [c for c in alphabet if c != seq[pos]]
        seq[pos] = random.choice(choices)
    return ''.join(seq)


# =============================================================================
# 2. Core timing function
# =============================================================================

def time_function(func, *args, repeats=5, **kwargs):
    """
    Run a function `repeats` times and return the MEDIAN elapsed time in seconds.
    Using median instead of mean avoids distortion from occasional slow runs
    (e.g., OS scheduling interrupts).

    Parameters
    ----------
    func    : callable – the algorithm function to time
    *args   : arguments to pass to the function
    repeats : int      – number of times to run (default 5)

    Returns
    -------
    float : median elapsed time in seconds
    """
    times = []
    for _ in range(repeats):
        start = time.perf_counter()   # high-resolution timer
        func(*args, **kwargs)
        end   = time.perf_counter()
        times.append(end - start)

    times.sort()
    return times[len(times) // 2]   # return the median value


def measure_memory(func, *args):
    """
    Run a function once and return peak memory usage in kilobytes.

    Uses Python's tracemalloc module, which tracks memory allocations
    made by the Python interpreter.

    Parameters
    ----------
    func  : callable – the algorithm function to profile
    *args : arguments to pass to the function

    Returns
    -------
    float : peak memory usage in kilobytes (KB)
    """
    tracemalloc.start()              # start tracking allocations
    func(*args)                      # run the function
    _, peak = tracemalloc.get_traced_memory()   # peak is in bytes
    tracemalloc.stop()               # stop tracking
    return peak / 1024               # convert bytes → kilobytes


# =============================================================================
# 3. Benchmark pairwise runtime vs. sequence length
# =============================================================================

def benchmark_pairwise_runtime():
    """
    Benchmark NW, SW, Gotoh, and Hirschberg across increasing sequence lengths.

    Returns a dictionary:
        { algorithm_name: { length: median_time_seconds } }
    """
    # Sequence lengths to test (add/remove values as needed)
    lengths = [50, 100, 200, 300, 500, 750, 1000]

    # Which algorithms to test
    algorithms = {
        "Needleman-Wunsch": lambda s1, s2: needleman_wunsch(s1, s2),
        "Smith-Waterman"  : lambda s1, s2: smith_waterman(s1, s2),
        "Gotoh"           : lambda s1, s2: gotoh_affine_gap(s1, s2),
        "Hirschberg"      : lambda s1, s2: hirschberg(s1, s2),
        "Banded (k=10)"   : lambda s1, s2: banded_nw(s1, s2, k=10),
    }

    results = {name: {} for name in algorithms}

    print("\n" + "=" * 60)
    print("PAIRWISE RUNTIME BENCHMARK")
    print("=" * 60)
    print(f"{'Length':>8} | " + " | ".join(f"{n:>18}" for n in algorithms))
    print("-" * (8 + 3 + len(algorithms) * 21))

    for length in lengths:
        # Generate a fresh random pair for each length
        s1 = random_sequence(length)
        s2 = mutate_sequence(s1, num_mutations=length // 10)  # ~10% different

        row = f"{length:>8} | "
        for name, func in algorithms.items():
            t = time_function(func, s1, s2, repeats=3)
            results[name][length] = t
            row += f"{t:>18.4f} | "

        print(row)

    return results, lengths


# =============================================================================
# 4. Benchmark memory (NW vs Hirschberg)
# =============================================================================

def benchmark_memory():
    """
    Compare peak memory usage of standard NW (O(nm) space)
    vs Hirschberg (O(n) space).

    Returns two lists: lengths, nw_memory_kb, hirschberg_memory_kb
    """
    lengths = [100, 200, 400, 600, 800, 1000]

    nw_mem  = []
    hb_mem  = []

    print("\n" + "=" * 60)
    print("MEMORY BENCHMARK: NW vs Hirschberg")
    print("=" * 60)
    print(f"{'Length':>8} | {'NW (KB)':>15} | {'Hirschberg (KB)':>18}")
    print("-" * 50)

    for length in lengths:
        s1 = random_sequence(length)
        s2 = mutate_sequence(s1, num_mutations=length // 10)

        mem_nw = measure_memory(needleman_wunsch, s1, s2)
        mem_hb = measure_memory(hirschberg, s1, s2)

        nw_mem.append(mem_nw)
        hb_mem.append(mem_hb)

        print(f"{length:>8} | {mem_nw:>15.2f} | {mem_hb:>18.2f}")

    return lengths, nw_mem, hb_mem


# =============================================================================
# 5. Benchmark BLAST-lite vs sequence length
# =============================================================================

def benchmark_blast():
    """
    Benchmark BLAST-lite heuristic performance.
    Expected: roughly linear O(n) growth vs. O(n^2) for exact methods.
    """
    lengths = [100, 250, 500, 750, 1000, 1500, 2000]
    blast_times = []
    nw_times    = []

    print("\n" + "=" * 60)
    print("BLAST vs NW RUNTIME COMPARISON")
    print("=" * 60)
    print(f"{'Length':>8} | {'BLAST (s)':>12} | {'NW (s)':>12}")
    print("-" * 38)

    for length in lengths:
        s1 = random_sequence(length)
        s2 = mutate_sequence(s1, num_mutations=length // 5)

        t_blast = time_function(blast_lite, s1, s2, k=4, repeats=3)
        t_nw    = time_function(needleman_wunsch, s1, s2, repeats=3)

        blast_times.append(t_blast)
        nw_times.append(t_nw)

        print(f"{length:>8} | {t_blast:>12.4f} | {t_nw:>12.4f}")

    return lengths, blast_times, nw_times


# =============================================================================
# 6. Parameter sensitivity for NW gap penalty
# =============================================================================

def benchmark_gap_sensitivity():
    """
    Run NW with different gap penalties and record alignment score and gap count.
    This shows how sensitive the algorithm is to parameter choice.
    """
    s1 = "GATTACAACTTGCGTATGCAG"
    s2 = "GATCCAGTTCAAATGCGTATG"

    gap_values   = [-1, -2, -3, -5, -7, -10, -15, -20]
    scores       = []
    gap_counts   = []

    print("\n" + "=" * 60)
    print("NW GAP PENALTY SENSITIVITY")
    print("=" * 60)
    print(f"{'Gap Penalty':>12} | {'Score':>8} | {'# Gaps':>8} | Alignment")
    print("-" * 80)

    for gap in gap_values:
        a1, a2, score = needleman_wunsch(s1, s2, gap=gap)
        n_gaps = a1.count('-') + a2.count('-')
        scores.append(score)
        gap_counts.append(n_gaps)
        print(f"{gap:>12} | {score:>8} | {n_gaps:>8} | {a1}")

    return gap_values, scores, gap_counts


# =============================================================================
# 7. Save results to CSV
# =============================================================================

def save_csv(filename, headers, rows):
    """Save benchmark results to a CSV file for later analysis."""
    path = os.path.join(RESULTS_DIR, filename)
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)
    print(f"  Saved: {path}")


# =============================================================================
# 8. Plot results (requires matplotlib)
# =============================================================================

def make_plots(pairwise_results, pairwise_lengths,
               mem_lengths, nw_mem, hb_mem,
               blast_lengths, blast_times, nw_times_blast,
               gap_values, gap_scores, gap_counts):
    """
    Generate all plots and save them as PNG files.
    Each PNG can be inserted directly into your report.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')   # non-interactive backend (works without a display)
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
    except ImportError:
        print("\nWARNING: matplotlib not installed.")
        print("Install it with:  pip install matplotlib")
        print("Skipping plot generation.")
        return

    # -------------------------------------------------------------------
    # Plot 1: Runtime vs sequence length (all pairwise methods)
    # -------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))
    colors  = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

    for (name, data), color in zip(pairwise_results.items(), colors):
        lengths_list = sorted(data.keys())
        times_list   = [data[l] for l in lengths_list]
        ax.plot(lengths_list, times_list, marker='o', label=name,
                color=color, linewidth=2)

    ax.set_xlabel("Sequence Length (characters)", fontsize=12)
    ax.set_ylabel("Median Time (seconds)",        fontsize=12)
    ax.set_title("Runtime vs. Sequence Length — Pairwise Alignment Methods",
                 fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "plot_runtime_pairwise.png"), dpi=150)
    plt.close()
    print("  Saved: plot_runtime_pairwise.png")

    # -------------------------------------------------------------------
    # Plot 2: Memory usage — NW vs Hirschberg
    # -------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(mem_lengths, nw_mem, marker='s', label="Needleman-Wunsch O(nm)",
            color='#d62728', linewidth=2)
    ax.plot(mem_lengths, hb_mem, marker='o', label="Hirschberg O(n+m)",
            color='#2ca02c', linewidth=2)
    ax.set_xlabel("Sequence Length",       fontsize=12)
    ax.set_ylabel("Peak Memory (KB)",      fontsize=12)
    ax.set_title("Memory Usage: NW vs. Hirschberg", fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "plot_memory_nw_vs_hirschberg.png"), dpi=150)
    plt.close()
    print("  Saved: plot_memory_nw_vs_hirschberg.png")

    # -------------------------------------------------------------------
    # Plot 3: BLAST-lite vs NW runtime comparison
    # -------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(blast_lengths, blast_times,    marker='D', label="BLAST-lite (heuristic)",
            color='#ff7f0e', linewidth=2)
    ax.plot(blast_lengths, nw_times_blast, marker='o', label="Needleman-Wunsch (exact)",
            color='#1f77b4', linewidth=2)
    ax.set_xlabel("Sequence Length",       fontsize=12)
    ax.set_ylabel("Median Time (seconds)", fontsize=12)
    ax.set_title("Heuristic vs. Exact: BLAST-lite vs. Needleman-Wunsch",
                 fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "plot_blast_vs_nw.png"), dpi=150)
    plt.close()
    print("  Saved: plot_blast_vs_nw.png")

    # -------------------------------------------------------------------
    # Plot 4: Gap penalty sensitivity
    # -------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.plot(gap_values, gap_scores, marker='o', color='#1f77b4', linewidth=2)
    ax1.set_xlabel("Gap Penalty",  fontsize=11)
    ax1.set_ylabel("Alignment Score", fontsize=11)
    ax1.set_title("Score vs. Gap Penalty", fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)

    ax2.plot(gap_values, gap_counts, marker='s', color='#d62728', linewidth=2)
    ax2.set_xlabel("Gap Penalty", fontsize=11)
    ax2.set_ylabel("Number of Gaps", fontsize=11)
    ax2.set_title("Gap Count vs. Gap Penalty", fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)

    plt.suptitle("Needleman-Wunsch Parameter Sensitivity", fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "plot_gap_sensitivity.png"), dpi=150)
    plt.close()
    print("  Saved: plot_gap_sensitivity.png")

# =============================================================================
# 9. Biological dataset alignment experiment
# =============================================================================

def benchmark_all_datasets():
    """
    Run alignment benchmarks on real-world biological FASTA files.

    Loads pre-defined datasets, performs pairwise alignments using
    Smith-Waterman, and prints the resulting alignment scores.
    """
    from data_loader import load_fasta

    # Mapping of dataset names to their respective file paths
    datasets = {
        'Globins':          os.path.join(PROJECT_ROOT, 'data', 'globins.fasta'),
        'Cytochrome c':     os.path.join(PROJECT_ROOT, 'data', 'cytochrome_c.fasta'),
        'Serine Proteases': os.path.join(PROJECT_ROOT, 'data', 'serine_proteases.fasta'),
        'Synthetic SNPs':   os.path.join(PROJECT_ROOT, 'data', 'synthetic_snps.fasta'),
    }

    print("\n" + "=" * 60)
    print("BIOLOGICAL DATASETS ALIGNMENT EXPERIMENT")
    print("=" * 60)

    for dataset_name, filepath in datasets.items():
        # Ensure data file exists before attempting to load
        if not os.path.exists(filepath):
            print(f"  Skipping {dataset_name}: File {filepath} not found.")
            continue

        seqs = load_fasta(filepath)
        names = list(seqs.keys())
        print(f'\n=== {dataset_name} ({len(seqs)} sequences) ===')

        # Perform all-vs-all pairwise alignment within the dataset
        for name1, name2 in itertools.combinations(names, 2):
            seq1, seq2 = seqs[name1], seqs[name2]

            # Use Smith-Waterman to find the local alignment score
            _, _, score = smith_waterman(seq1, seq2)

            print(f"  {name1} vs {name2} | Score: {score}")
# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("SEQUENCE ALIGNMENT BENCHMARK SUITE")
    print("Results will be saved to:", RESULTS_DIR)
    print("=" * 60)

    # Run all benchmarks
    pairwise_results, pairwise_lengths = benchmark_pairwise_runtime()
    mem_lengths, nw_mem, hb_mem        = benchmark_memory()
    blast_lengths, blast_times, nw_t   = benchmark_blast()
    gap_vals, gap_scores, gap_counts   = benchmark_gap_sensitivity()

    # Save CSV summaries
    print("\nSaving CSV files...")
    # Pairwise timing CSV
    alg_names = list(pairwise_results.keys())
    headers   = ["length"] + alg_names
    rows      = [[l] + [pairwise_results[n][l] for n in alg_names]
                 for l in pairwise_lengths]
    save_csv("timing_pairwise.csv", headers, rows)

    # Memory CSV
    save_csv("memory_nw_vs_hirschberg.csv",
             ["length", "nw_kb", "hirschberg_kb"],
             list(zip(mem_lengths, nw_mem, hb_mem)))

    # Generate all plots
    print("\nGenerating plots...")
    make_plots(pairwise_results, pairwise_lengths,
               mem_lengths, nw_mem, hb_mem,
               blast_lengths, blast_times, nw_t,
               gap_vals, gap_scores, gap_counts)

    benchmark_all_datasets()

    print("\n" + "=" * 60)
    print("BENCHMARK COMPLETE!")
    print(f"All results saved in: {RESULTS_DIR}")
    print("=" * 60)


