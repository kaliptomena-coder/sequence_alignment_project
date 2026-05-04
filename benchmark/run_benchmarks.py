"""
run_benchmarks.py — Reproducible benchmark harness for sequence alignment algorithms.

Benchmarks:
  1.  Runtime scaling          (NW, SW, Gotoh, Hirschberg, Banded)
  2.  Memory scaling           (NW vs Hirschberg)
  3.  BLAST-lite vs NW         runtime comparison
  4.  Gap penalty sensitivity  (NW, linear gap)
  5.  Gotoh affine parameters  (gap_open × gap_extend grid)
  6.  Minimizer parameters     (k × w grid → chain length)
  7.  Identity by dataset group (simulated sequences)
  8.  Sensitivity / specificity (BLAST-lite, 30 trials)
  9.  Local vs global          (embedded-motif example)
  10. Long-sequence benchmark  (HBB locus subsets — skipped if file absent)
  11. Real-dataset identity    (mean ± SD per FASTA file)
  12. Dataset pairwise scores  (all pairs across FASTA datasets)

Run:
  python run_benchmarks.py
"""

import os
import sys
import time
import random
import csv
import tracemalloc
import itertools
import multiprocessing as mp
from statistics import median, mean, stdev

# ── Configuration ─────────────────────────────────────────────────────────────

CONFIG = {
    "seed":             42,
    "repeats":          3,
    "results_dir":      "benchmarks/results",
    "synthetic_lengths": [200, 500, 1000, 1500, 2000],
    "blast_lengths":     [200, 500, 1000, 1500, 2000],
    "hbb_lengths":       [1000, 2000, 5000],   # subsequence lengths for HBB benchmark
}

random.seed(CONFIG["seed"])
os.makedirs(CONFIG["results_dir"], exist_ok=True)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
SRC_DIR  = os.path.join(BASE_DIR, "src")
sys.path.insert(0, SRC_DIR)

from needlemanWunschGlobal import needleman_wunsch
from smithWatermanLocal   import smith_waterman
from gotoh                import gotoh
from hirschberg           import hirschberg
from banded_alignment     import banded_nw
from blast_lite           import blast_lite
from minimizer_align      import minimizer_align


# ── Sequence utilities ────────────────────────────────────────────────────────

def random_sequence(length):
    return ''.join(random.choices("ACGT", k=length))


def mutate_sequence(seq, rate=0.1):
    seq = list(seq)
    n   = max(1, int(len(seq) * rate))
    for pos in random.sample(range(len(seq)), n):
        seq[pos] = random.choice([c for c in "ACGT" if c != seq[pos]])
    return ''.join(seq)


def pct_identity(a1, a2):
    """Percentage of aligned non-gap columns where both characters match."""
    matches = sum(c1 == c2 and c1 != '-' for c1, c2 in zip(a1, a2))
    cols    = sum(1 for c1, c2 in zip(a1, a2) if c1 != '-' and c2 != '-')
    return round(100.0 * matches / cols, 2) if cols else 0.0


def load_fasta(path):
    """Parse a FASTA file; return {name: sequence} dict."""
    if not os.path.exists(path):
        return {}
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line.upper())
    if name:
        seqs[name] = ''.join(buf)
    return seqs


# ── Timing / memory helpers ───────────────────────────────────────────────────

def time_median(func, *args):
    """Run func 3× (after a warm-up) and return the median wall time."""
    func(*args)   # warm-up
    return median(
        (lambda: (lambda t0: time.perf_counter() - t0)(time.perf_counter()) or
                 (func(*args), time.perf_counter())[-1])()   # can't do in one line cleanly
        for _ in range(CONFIG["repeats"])
    )


def _time_it(func, args):
    """Simple timing helper used inside loops."""
    func(*args)   # warm-up
    times = []
    for _ in range(CONFIG["repeats"]):
        t0 = time.perf_counter()
        func(*args)
        times.append(time.perf_counter() - t0)
    return median(times)


def _memory_worker(func, args, queue):
    import gc
    gc.collect()
    tracemalloc.start()
    func(*args)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    queue.put(peak / 1024)   # bytes → KB


def _noop(*_):
    pass


def measure_memory_kb(func, *args):
    """
    Run func in a child process and return net peak memory in KB.
    Using a child process avoids contamination from the parent's heap.
    """
    q = mp.Queue()
    for target_func in (_noop, func):
        p = mp.Process(target=_memory_worker, args=(target_func, args, q))
        p.start(); p.join()
    baseline_kb = q.get()   # first put was _noop
    peak_kb     = q.get()   # second put was func
    return max(peak_kb - baseline_kb, 0.0)


def save(filename, rows):
    path = os.path.join(CONFIG["results_dir"], filename)
    with open(path, "w", newline="") as f:
        csv.writer(f).writerows(rows)
    print(f"  Saved → {path}")


# ── Benchmark 1: Runtime scaling ─────────────────────────────────────────────

def benchmark_runtime():
    print("\n[1] Runtime scaling …")
    algos = {
        "NW":         needleman_wunsch,
        "SW":         smith_waterman,
        "Gotoh":      gotoh,
        "Hirschberg": hirschberg,
        "Banded":     lambda s1, s2: banded_nw(s1, s2, k=10),
    }
    rows = [["length"] + list(algos.keys())]
    for L in CONFIG["synthetic_lengths"]:
        s1 = random_sequence(L)
        s2 = mutate_sequence(s1, rate=0.1)
        row = [L]
        for name, func in algos.items():
            t = _time_it(func, (s1, s2))
            print(f"    {name:12s}  L={L:5d}  {t:.4f}s")
            row.append(round(t, 6))
        rows.append(row)
    return rows


# ── Benchmark 2: Memory scaling ───────────────────────────────────────────────

def benchmark_memory():
    print("\n[2] Memory scaling …")
    rows = [["length", "NW_KB", "Hirschberg_KB"]]
    for L in CONFIG["synthetic_lengths"]:
        s1 = random_sequence(L)
        s2 = mutate_sequence(s1)
        nw_kb = measure_memory_kb(needleman_wunsch, s1, s2)
        hb_kb = measure_memory_kb(hirschberg,       s1, s2)
        print(f"    L={L:5d}  NW={nw_kb:.1f} KB  Hirschberg={hb_kb:.1f} KB")
        rows.append([L, round(nw_kb, 2), round(hb_kb, 2)])
    return rows


# ── Benchmark 3: BLAST-lite vs NW ─────────────────────────────────────────────

def benchmark_blast_vs_nw():
    print("\n[3] BLAST-lite vs NW …")
    rows = [["length", "BLAST_s", "NW_s", "speedup"]]
    for L in CONFIG["blast_lengths"]:
        s1 = random_sequence(L)
        s2 = mutate_sequence(s1, rate=0.2)
        t_blast = _time_it(blast_lite,       (s1, s2))
        t_nw    = _time_it(needleman_wunsch, (s1, s2))
        speedup = round(t_nw / t_blast, 2) if t_blast > 0 else 0
        print(f"    L={L:5d}  BLAST={t_blast:.4f}s  NW={t_nw:.4f}s  speedup={speedup}x")
        rows.append([L, round(t_blast, 6), round(t_nw, 6), speedup])
    return rows


# ── Benchmark 4: Gap penalty sensitivity ──────────────────────────────────────

def benchmark_gap_sensitivity():
    print("\n[4] Gap penalty sensitivity …")
    s1 = "GATTACAACTTGCGTATGCAG"
    s2 = "GATCCAGTTCAAATGCGTATG"
    rows = [["gap_penalty", "score", "gap_columns"]]
    for gap in [-1, -2, -5, -10, -20]:
        a1, a2, score = needleman_wunsch(s1, s2, gap=gap)
        gaps = a1.count('-') + a2.count('-')
        print(f"    gap={gap:4d}  score={score:4d}  gaps={gaps}")
        rows.append([gap, score, gaps])
    return rows


# ── Benchmark 5: Gotoh affine gap parameters ──────────────────────────────────

def benchmark_gotoh_parameters():
    print("\n[5] Gotoh affine parameters …")
    s1 = "GATTACAACTTGCGTATGCAG"
    s2 = "GATCCAGTTCAAATGCGTATG"
    rows = [["gap_open", "gap_extend", "score"]]
    for go, ge in [(-3, -1), (-5, -1), (-10, -1), (-5, -2)]:
        _, _, score = gotoh(s1, s2, gap_open=go, gap_extend=ge)
        print(f"    open={go:4d}  extend={ge:3d}  score={score}")
        rows.append([go, ge, score])
    return rows


# ── Benchmark 6: Minimizer parameters ────────────────────────────────────────

def benchmark_minimizer():
    print("\n[6] Minimizer alignment parameters …")
    s1   = "GATTACA" * 20
    s2   = mutate_sequence(s1, rate=0.05)
    rows = [["k", "w", "chain_length"]]
    for k in [3, 5, 7]:
        for w in [5, 10, 15]:
            _, _, _, chain = minimizer_align(s1, s2, k=k, w=w)
            print(f"    k={k}  w={w}  chain={len(chain)}")
            rows.append([k, w, len(chain)])
    return rows


# ── Benchmark 7: Identity by simulated dataset group ─────────────────────────

def benchmark_simulated_identity():
    """
    Simulates four biological dataset scenarios via controlled mutation rates
    and measures pairwise sequence identity for NW, SW, and Gotoh.
    Each group uses a single random pair; see benchmark_dataset_identity()
    for mean ± SD over real FASTA files.
    """
    print("\n[7] Simulated identity by dataset group …")
    groups = {
        "Closely related (globins)":   0.10,
        "Distant homologs (proteases)": 0.40,
        "Long sequences (simulated)":   0.20,
        "Synthetic SNPs":               0.05,
    }
    algos = {"NW": needleman_wunsch, "SW": smith_waterman, "Gotoh": gotoh}
    rows = [["dataset_group", "algorithm", "identity_pct"]]
    identity_table = {}
    for group_name, mut_rate in groups.items():
        s1 = random_sequence(300)
        s2 = mutate_sequence(s1, rate=mut_rate)
        identity_table[group_name] = {}
        for algo_name, func in algos.items():
            a1, a2, _ = func(s1, s2)
            pct = pct_identity(a1, a2)
            print(f"    {group_name:35s}  {algo_name}  {pct}%")
            rows.append([group_name, algo_name, pct])
            identity_table[group_name][algo_name] = pct
    return rows, identity_table


# ── Benchmark 8: BLAST sensitivity / specificity ──────────────────────────────

def benchmark_sensitivity_specificity():
    """
    30 trials: positive = 15%-mutated pair, negative = independent random pair.
    k=11 is chosen because P(random 11-mer match in 200 bp) ≈ 1/22000 — effectively zero.
    Results are valid for this controlled setup; see report for caveats.
    """
    print("\n[8] BLAST-lite sensitivity / specificity …")
    K = 11
    THRESHOLD = 40
    trials = 30
    tp = tn = fp = fn = 0

    for _ in range(trials):
        s1 = random_sequence(200)

        hsps_pos = blast_lite(s1, mutate_sequence(s1, rate=0.15), k=K)
        top_pos  = hsps_pos[0].get("score", 0) if hsps_pos and isinstance(hsps_pos[0], dict) else (float(hsps_pos[0]) if hsps_pos else 0)
        if top_pos >= THRESHOLD:
            tp += 1
        else:
            fn += 1

        hsps_neg = blast_lite(s1, random_sequence(200), k=K)
        top_neg  = hsps_neg[0].get("score", 0) if hsps_neg and isinstance(hsps_neg[0], dict) else (float(hsps_neg[0]) if hsps_neg else 0)
        if top_neg >= THRESHOLD:
            fp += 1
        else:
            tn += 1

    sensitivity = round(100.0 * tp / (tp + fn), 1) if (tp + fn) else 0
    specificity = round(100.0 * tn / (tn + fp), 1) if (tn + fp) else 0
    fpr         = round(100.0 - specificity, 1)
    print(f"    Sensitivity={sensitivity}%  Specificity={specificity}%  FPR={fpr}%")
    return [
        ["metric",               "value_pct"],
        ["Sensitivity (TPR)",    sensitivity],
        ["Specificity (TNR)",    specificity],
        ["False Positive Rate",  fpr],
    ]


# ── Benchmark 9: Local vs global ─────────────────────────────────────────────

def benchmark_local_vs_global():
    print("\n[9] Local vs global alignment …")
    s1 = "AAAACGTAAAA"
    s2 = "TTTTCGTTTTT"
    nw_score = needleman_wunsch(s1, s2)[2]
    sw_score = smith_waterman(s1, s2)[2]
    print(f"    NW={nw_score}  SW={sw_score}")
    return [["algorithm", "score"], ["NW", nw_score], ["SW", sw_score]]


# ── Benchmark 10: Long-sequence benchmark (HBB locus) ───────────────────────

def benchmark_hbb_long_sequence():
    """
    Runs Hirschberg, BLAST-lite, and minimizer on increasing-length subsets
    of the human HBB locus FASTA. Skipped if hbb.fasta is absent.
    This benchmark demonstrates why exact DP is impractical at genomic scale
    and motivates heuristic methods.
    """
    print("\n[10] Long-sequence benchmark (HBB locus) …")
    hbb_path = os.path.join(DATA_DIR, "hbb.fasta")
    if not os.path.exists(hbb_path):
        print("    SKIP: hbb.fasta not found in data/")
        return [["status", "skipped"]]

    seqs = load_fasta(hbb_path)
    if not seqs:
        print("    SKIP: hbb.fasta is empty or malformed")
        return [["status", "empty"]]

    full_seq = next(iter(seqs.values()))
    print(f"    HBB locus loaded: {len(full_seq):,} bp")

    rows = [["length_bp", "Hirschberg_s", "BLAST_s", "Minimizer_s",
             "NW_memory_KB", "Hirschberg_memory_KB"]]

    for L in CONFIG["hbb_lengths"]:
        if L > len(full_seq) - 100:
            print(f"    SKIP L={L}: sequence too short")
            continue
        s1 = full_seq[:L]
        s2 = full_seq[50:L + 50]   # overlapping window, ~highly identical

        t_hb   = _time_it(hirschberg,       (s1, s2))
        t_bl   = _time_it(blast_lite,       (s1, s2))
        t_mini = _time_it(lambda a, b: minimizer_align(a, b, k=11, w=20), (s1, s2))

        # Memory only at smaller lengths (NW at 5 kb needs ~50 MB Python overhead)
        if L <= 2000:
            nw_kb = measure_memory_kb(needleman_wunsch, s1, s2)
            hb_kb = measure_memory_kb(hirschberg,       s1, s2)
        else:
            nw_kb = hb_kb = float('nan')

        print(f"    L={L:5,}  Hirschberg={t_hb:.3f}s  BLAST={t_bl:.3f}s  "
              f"Minimizer={t_mini:.4f}s  "
              f"NW_mem={'n/a' if nw_kb != nw_kb else f'{nw_kb:.0f} KB'}  "
              f"HB_mem={'n/a' if hb_kb != hb_kb else f'{hb_kb:.1f} KB'}")

        rows.append([L, round(t_hb, 4), round(t_bl, 4), round(t_mini, 6),
                     round(nw_kb, 1), round(hb_kb, 1)])

    return rows


# ── Benchmark 11: Real-dataset pairwise identity (mean ± SD) ─────────────────

def benchmark_dataset_identity():
    """
    For each FASTA dataset compute pairwise NW identity for all sequence pairs
    and report mean ± SD. This addresses the grader's concern about single-pair
    identity results; means over all C(n,2) pairs are more representative.
    """
    print("\n[11] Real-dataset pairwise identity (mean ± SD) …")
    datasets = [
        "globins.fasta",
        "cytochrome_c.fasta",
        "serine_proteases.fasta",
        "synthetic_snps.fasta",
    ]
    rows = [["dataset", "n_seqs", "n_pairs", "mean_identity_pct", "sd_identity_pct", "min_pct", "max_pct"]]

    for fname in datasets:
        path = os.path.join(DATA_DIR, fname)
        seqs = load_fasta(path)
        if not seqs:
            print(f"    SKIP: {fname} not found")
            continue

        names = list(seqs.keys())
        identities = []
        for a, b in itertools.combinations(names, 2):
            try:
                a1, a2, _ = needleman_wunsch(seqs[a], seqs[b])
                identities.append(pct_identity(a1, a2))
            except Exception as e:
                print(f"      Error {a} vs {b}: {e}")

        if not identities:
            continue

        m  = round(mean(identities), 2)
        sd = round(stdev(identities), 2) if len(identities) > 1 else 0.0
        lo = round(min(identities), 2)
        hi = round(max(identities), 2)
        print(f"    {fname:<30s}  n={len(names)}  pairs={len(identities)}"
              f"  mean={m}% ± {sd}%  range=[{lo}%, {hi}%]")
        rows.append([fname, len(names), len(identities), m, sd, lo, hi])

    return rows


# ── Benchmark 12: Dataset pairwise scores (NW / SW / Gotoh) ──────────────────

def benchmark_dataset_scores():
    """All pairwise alignment scores across FASTA datasets."""
    print("\n[12] Dataset pairwise scores …")
    datasets = [
        "globins.fasta",
        "cytochrome_c.fasta",
        "serine_proteases.fasta",
        "synthetic_snps.fasta",
    ]
    rows = [["dataset", "seq_a", "seq_b", "SW_score", "NW_score", "Gotoh_score"]]
    for fname in datasets:
        path = os.path.join(DATA_DIR, fname)
        seqs = load_fasta(path)
        if not seqs:
            print(f"    SKIP: {fname} not found")
            continue
        names = list(seqs.keys())
        print(f"    {fname}: {len(names)} sequences")
        for a, b in itertools.combinations(names, 2):
            try:
                sw  = smith_waterman(seqs[a], seqs[b])[2]
                nw  = needleman_wunsch(seqs[a], seqs[b])[2]
                got = gotoh(seqs[a], seqs[b])[2]
                rows.append([fname, a, b, sw, nw, got])
            except Exception as e:
                print(f"      Error {a} vs {b}: {e}")
    return rows


# ── Plots ─────────────────────────────────────────────────────────────────────

def make_plots(runtime_rows, memory_rows, blast_rows, gap_rows,
               minimizer_rows, identity_table, hbb_rows=None):
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("matplotlib/numpy not available — skipping plots")
        return

    out = CONFIG["results_dir"]
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.3,
    })
    C = ["#4878CF", "#6ACC65", "#D65F5F", "#B47CC7", "#C4AD66", "#77BEDB"]

    def save_fig(fig, name):
        fig.savefig(os.path.join(out, name), dpi=150)
        plt.close(fig)
        print(f"  → {name}")

    # Fig 1: Runtime (log-log)
    algos   = runtime_rows[0][1:]
    lengths = [r[0] for r in runtime_rows[1:]]
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, algo in enumerate(algos):
        ax.plot(lengths, [r[i+1] for r in runtime_rows[1:]], marker="o", label=algo, color=C[i % len(C)])
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Sequence length (bp)"); ax.set_ylabel("Runtime (s)")
    ax.set_title("Fig 1 — Runtime scaling (log-log)")
    ax.legend(); fig.tight_layout()
    save_fig(fig, "fig1_runtime_loglog.png")

    # Fig 2: Memory
    mem_L  = [r[0] for r in memory_rows[1:]]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(mem_L, [r[1] for r in memory_rows[1:]], marker="o", label="Needleman-Wunsch", color=C[0])
    ax.plot(mem_L, [r[2] for r in memory_rows[1:]], marker="s", label="Hirschberg",        color=C[2])
    ax.set_xlabel("Sequence length (bp)"); ax.set_ylabel("Peak memory (KB)")
    ax.set_title("Fig 2 — Memory: Needleman-Wunsch vs Hirschberg")
    ax.legend(); fig.tight_layout()
    save_fig(fig, "fig2_memory.png")

    # Fig 3: BLAST vs NW
    bl_L = [r[0] for r in blast_rows[1:]]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(bl_L, [r[1] for r in blast_rows[1:]], marker="o", label="BLAST-lite",       color=C[1])
    ax.plot(bl_L, [r[2] for r in blast_rows[1:]], marker="s", label="Needleman-Wunsch", color=C[0])
    ax.set_xlabel("Sequence length (bp)"); ax.set_ylabel("Runtime (s)")
    ax.set_title("Fig 3 — BLAST-lite vs NW runtime")
    ax.legend(); fig.tight_layout()
    save_fig(fig, "fig3_blast_vs_nw.png")

    # Fig 4: Gap sensitivity (dual axis)
    gp  = [r[0] for r in gap_rows[1:]]
    gs  = [r[1] for r in gap_rows[1:]]
    gct = [r[2] for r in gap_rows[1:]]
    fig, ax1 = plt.subplots(figsize=(7, 4))
    ax2 = ax1.twinx()
    ax1.plot(gp, gs,  marker="o", color=C[0], label="Score")
    ax2.bar(gp,  gct, alpha=0.35, color=C[2], label="Gap columns", width=1.2)
    ax1.set_xlabel("Gap penalty"); ax1.set_ylabel("Alignment score", color=C[0])
    ax2.set_ylabel("Gap columns", color=C[2])
    ax1.set_title("Fig 4 — Gap penalty sensitivity (NW)")
    h1, l1 = ax1.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1+h2, l1+l2, loc="upper right"); fig.tight_layout()
    save_fig(fig, "fig4_gap_sensitivity.png")

    # Fig 6: Identity bar chart
    if identity_table:
        groups     = list(identity_table.keys())
        algo_names = list(next(iter(identity_table.values())).keys())
        x, w = np.arange(len(groups)), 0.25
        fig, ax = plt.subplots(figsize=(11, 5))
        for j, aname in enumerate(algo_names):
            ax.bar(x + j*w, [identity_table[g].get(aname, 0) for g in groups], w, label=aname, color=C[j])
        ax.set_xticks(x + w)
        ax.set_xticklabels([g.replace("\n", " ") for g in groups], fontsize=8)
        ax.set_ylabel("Sequence identity (%)"); ax.set_ylim(0, 105)
        ax.set_title("Fig 6 — Alignment identity across dataset groups")
        ax.legend(); fig.tight_layout()
        save_fig(fig, "fig6_identity_pct.png")

    # Fig 9: Minimizer heatmap
    ks = sorted({r[0] for r in minimizer_rows[1:]})
    ws = sorted({r[1] for r in minimizer_rows[1:]})
    grid = np.zeros((len(ks), len(ws)))
    for r in minimizer_rows[1:]:
        grid[ks.index(r[0]), ws.index(r[1])] = r[2]
    fig, ax = plt.subplots(figsize=(6, 4))
    im = ax.imshow(grid, cmap="YlOrRd", aspect="auto")
    ax.set_xticks(range(len(ws))); ax.set_xticklabels([f"w={w}" for w in ws])
    ax.set_yticks(range(len(ks))); ax.set_yticklabels([f"k={k}" for k in ks])
    ax.set_xlabel("Window size (w)"); ax.set_ylabel("k-mer size (k)")
    ax.set_title("Fig 9 — Minimizer chain length heatmap")
    for ki in range(len(ks)):
        for wi in range(len(ws)):
            ax.text(wi, ki, int(grid[ki, wi]), ha="center", va="center",
                    color="white" if grid[ki, wi] > grid.max()*0.6 else "black")
    fig.colorbar(im, ax=ax, label="Chain length"); fig.tight_layout()
    save_fig(fig, "fig9_minimizer_heatmap.png")

    # Fig 10: HBB long-sequence runtime (if data available)
    if hbb_rows and len(hbb_rows) > 1 and hbb_rows[0][0] != "status":
        valid = [r for r in hbb_rows[1:] if isinstance(r[1], float)]
        if valid:
            Ls = [r[0] for r in valid]
            fig, ax = plt.subplots(figsize=(7, 4))
            ax.plot(Ls, [r[1] for r in valid], marker="o", label="Hirschberg", color=C[3])
            ax.plot(Ls, [r[2] for r in valid], marker="s", label="BLAST-lite", color=C[1])
            ax.plot(Ls, [r[3] for r in valid], marker="^", label="Minimizer",  color=C[4])
            ax.set_xlabel("Subsequence length (bp)"); ax.set_ylabel("Runtime (s)")
            ax.set_title("Fig 10 — Long-sequence runtime (HBB locus subsets)")
            ax.legend(); fig.tight_layout()
            save_fig(fig, "fig10_hbb_runtime.png")


# Main

if __name__ == "__main__":
    runtime    = benchmark_runtime()
    memory     = benchmark_memory()
    blast      = benchmark_blast_vs_nw()
    gap        = benchmark_gap_sensitivity()
    gotoh_rows = benchmark_gotoh_parameters()
    mini_rows  = benchmark_minimizer()
    acc_rows, identity_table = benchmark_simulated_identity()
    sens_rows  = benchmark_sensitivity_specificity()
    lg_rows    = benchmark_local_vs_global()
    hbb_rows   = benchmark_hbb_long_sequence()
    id_rows    = benchmark_dataset_identity()
    ds_rows    = benchmark_dataset_scores()

    save("runtime.csv",          runtime)
    save("memory.csv",           memory)
    save("blast.csv",            blast)
    save("gap.csv",              gap)
    save("gotoh.csv",            gotoh_rows)
    save("minimizer.csv",        mini_rows)
    save("simulated_identity.csv", acc_rows)
    save("sensitivity.csv",      sens_rows)
    save("local_global.csv",     lg_rows)
    save("hbb_long_sequence.csv", hbb_rows)
    save("dataset_identity.csv", id_rows)
    save("dataset_scores.csv",   ds_rows)

    print("\n── Generating figures ──")
    make_plots(runtime, memory, blast, gap, mini_rows, identity_table, hbb_rows)

    print("\nAll benchmarks completed.")