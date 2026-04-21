# Sequence Alignment Project

**A Comprehensive Implementation, Evaluation and Comparison of Sequence Alignment Methods**

> MSc Bioinformatics 
> Repository: https://github.com/kaliptomena-coder/sequence_alignment_project

---

## Overview

This project implements **9 sequence alignment algorithms + T-Coffee consistency** entirely in Python, with a full unit-test suite (**46 tests, 0 failures**), runtime/memory benchmarks, and experiments on biological and synthetic datasets.

| Category | Algorithm | Module |
|---|---|---|
| **Pairwise exact** | Needleman–Wunsch (global) | `src/needlemanWunschGlobal.py` |
| | Smith–Waterman (local) | `src/smithWatermanLocal.py` |
| | Gotoh affine gap penalties | `src/gotohAffineGap.py` |
| | Hirschberg O(n+m) space | `src/hirschberg.py` |
| | Banded dynamic programming | `src/banded_alignment.py` |
| **Heuristic** | BLAST-style seed-and-extend | `src/blast_lite.py` |
| | Minimizer-based alignment | `src/minimizer_align.py` |
| **MSA** | Progressive ClustalW-style + UPGMA | `src/progressive_alignment.py` |
| | Iterative refinement MUSCLE-style | `src/iterative_refinement.py` |
| | Profile HMM + Viterbi | `src/profile_hmm.py` |
| **Optional (10)** | T-Coffee consistency library | `src/consistency_library.py` |

---

## Quick Start

### 1 — Clone

```bash
git clone https://github.com/kaliptomena-coder/sequence_alignment_project.git
cd sequence_alignment_project
```

### 2 — Install dependencies

**Option A: conda (recommended)**
```bash
conda env create -f environment.yml
conda activate aligner_env
```

**Option B: pip**
```bash
pip install -r requirements.txt
```

> All nine alignment algorithms use **only the Python standard library** — no numpy or biopython required. `matplotlib` is needed only for benchmark plots.

### 3 — Run unit tests

```bash
python tests/test_all.py
```

Expected final line:
```
  RESULTS: 46 passed, 0 failed, 0 errors
```

### 4 — Run benchmarks

```bash
python benchmark/run_benchmarks.py
```

Saves runtime/memory CSV files and four PNG plots to `benchmark/results/`.

### 5 — Generate synthetic datasets

```bash
python data/generate_synthetic.py
```

Creates `synthetic_snps.fasta`, `synthetic_indels.fasta`, `synthetic_shuffled.fasta`, `synthetic_protein_family.fasta`, and `ground_truth.txt` in `data/` (uses `random.seed(42)` — fully reproducible).

---

## Project Structure

```
sequence_alignment_project/
├── src/
│   ├── needlemanWunschGlobal.py    # NW global alignment
│   ├── smithWatermanLocal.py       # SW local alignment
│   ├── gotohAffineGap.py           # Gotoh affine gap alignment
│   ├── hirschberg.py               # Hirschberg divide-and-conquer (O(n+m) space)
│   ├── banded_alignment.py         # Banded NW + graceful fallback to full NW
│   ├── blast_lite.py               # Heuristic seed-and-extend (BLAST-style)
│   ├── minimizer_align.py          # Minimizer sketch → anchor chain → NW gap-fill
│   ├── progressive_alignment.py    # ClustalW-style progressive MSA
│   ├── iterative_refinement.py     # MUSCLE-style Sum-of-Pairs iterative refinement
│   ├── profile_hmm.py              # Profile HMM with Viterbi decoding
│   ├── consistency_library.py      # T-Coffee consistency (optional method 10)
│   ├── distance_matrix.py          # All-vs-all NW identity distance matrix
│   ├── upgma.py                    # UPGMA guide tree construction
│   └── data_loader.py              # FASTA file parser
├── tests/
│   └── test_all.py                 # 46 unit tests across 9 algorithm classes
├── data/
│   ├── globins.fasta               # 9 globin sequences (HBA, HBB, MYG × human/horse/mouse)
│   ├── cytochrome_c.fasta          # 3 cytochrome c sequences (human, horse, mouse)
│   ├── serine_proteases.fasta      # 5 serine protease sequences (distant homologs)
│   ├── synthetic_snps.fasta        # 2,000 bp template + 5 × 10-SNP variants
│   ├── synthetic_indels.fasta      # 500 bp template + 5 × indel variants
│   ├── synthetic_shuffled.fasta    # Negative control (shuffled sequences)
│   ├── ground_truth.txt            # Exact list of every synthetic edit made
│   └── generate_synthetic.py       # Script to regenerate synthetic datasets
├── benchmark/
│   ├── run_benchmarks.py           # Full benchmark suite (runtime, memory, sensitivity)
│   └── results/                    # CSV + PNG output (auto-generated)
├── environment.yml                 # Conda environment specification
├── requirements.txt                # pip requirements
└── README.md                       # This file
```

---

## Using Individual Algorithms

```python
import sys
sys.path.insert(0, 'src')

# Pairwise Exact

from needlemanWunschGlobal import needleman_wunsch
a1, a2, score = needleman_wunsch('GATTACA', 'GATCA')
a1, a2, score = needleman_wunsch('GATTACA', 'GATCA', match=2, mismatch=-1, gap=-2)

from smithWatermanLocal import smith_waterman
a1, a2, score = smith_waterman('TTTACGTTTTT', 'ACGT', match=2, mismatch=-1, gap=-2)

from gotohAffineGap import gotoh_affine_gap
a1, a2, score = gotoh_affine_gap('GATTACA', 'GATCA', gap_open=-5, gap_extend=-1)

from hirschberg import hirschberg
a1, a2 = hirschberg('GATTACA', 'GATCA')   # O(n+m) space

from banded_alignment import banded_nw
a1, a2, score = banded_nw('GATTACA', 'GATTACA', k=5)
# If |len(s1)-len(s2)| > k, automatically falls back to full NW

# Heuristic

from blast_lite import blast_lite
hsps = blast_lite(query, target, k=4, threshold=5)
# Returns: [{'score': int, 'query_pos': int, 'target_pos': int, 'alignment': (str,str)}, ...]

from minimizer_align import minimizer_align
a1, a2, score, chain = minimizer_align(query, target, k=4, w=8)
# chain = list of (q_start, t_start, length) anchor tuples

# Multiple Sequence Alignment

from iterative_refinement import iterative_refinement, sum_of_pairs
from needlemanWunschGlobal import needleman_wunsch

initial_msa = ['ACGT--ACGT', 'ACGTAAACGT', 'ACGT--ACGG']
refined_msa, final_score, history = iterative_refinement(
    initial_msa, needleman_wunsch, max_iterations=50
)
print(f'SP score: {sum_of_pairs(initial_msa)} -> {final_score}')

from profile_hmm import ProfileHMM
model = ProfileHMM(training_msa, gap_threshold=0.5)
path, score = model.viterbi('ACGTACGT')

# Load FASTA files

from data_loader import load_fasta
seqs = load_fasta('globins.fasta')    # returns {name: sequence}
for name, seq in seqs.items():
    print(f'{name}: {len(seq)} aa')
```

---

## Benchmark Results

All results are **fully reproducible** by running `python benchmark/run_benchmarks.py`.

### Runtime (median of 3 repeats, ~10% mutated DNA pairs)

| Length | NW (s) | SW (s) | Gotoh (s) | Hirschberg (s) | Banded k=10 (s) |
|--------|--------|--------|-----------|----------------|-----------------|
| 100    | 0.0020 | 0.0020 | 0.0046    | 0.0030         | 0.0007          |
| 300    | 0.0168 | 0.0167 | 0.0417    | 0.0220         | 0.0022          |
| 500    | 0.0485 | 0.0475 | 0.1261    | 0.0622         | 0.0047          |
| 1000   | 0.2113 | 0.2043 | 0.5343    | 0.2565         | **0.0132**      |

All exact methods scale **O(nm)**. Banded NW (k=10) is **~16x faster** than full NW at 1,000 bp.

### Memory (tracemalloc peak)

| Length | NW (KB) | Hirschberg (KB) | Ratio |
|--------|---------|-----------------|-------|
| 100    | 355     | 9               | 38x   |
| 400    | 5,703   | 48              | 120x  |
| 600    | 13,225  | 86              | 153x  |
| 1000   | 38,675  | 166             | **233x** |

Hirschberg achieves **O(n+m) space** vs O(nm) for standard NW.

### BLAST-lite speedup

Consistently **~3x faster** than NW across all tested lengths in pure Python. A compiled C implementation would show much larger speedups by revealing the true O(n) algorithmic advantage.

---

## Unit Tests

```
46 passed, 0 failed, 0 errors
```

| Class | Tests | Key Invariants Checked |
|-------|-------|------------------------|
| TestNeedlemanWunsch | 9 | Score symmetry, gap monotonicity, GATTACA alignment, equal-length output |
| TestSmithWaterman | 6 | Score >= 0, local substring detection, empty on no match |
| TestGotoh | 5 | Affine gap in output, no gap for identical seqs |
| TestHirschberg | 4 | **Score = NW score on 5 test pairs**, valid chars {ACGT-} |
| TestBandedDP | 3 | Equals full NW; narrow band triggers graceful fallback |
| TestBLASTLite | 5 | HSP fields present, sorted descending, equal alignment lengths |
| TestMinimizerAlign | 4 | Anchor 3-tuple structure, fewer anchors for different seqs |
| TestIterativeRefinement | 4 | SP score non-decreasing, equal lengths after refinement |
| TestProfileHMM | 6 | Emissions in [0,1], transitions sum=1, **Viterbi score not -inf** |

The `test_viterbi_score_not_neg_inf` test specifically validates the **Laplace pseudo-count fix** that prevents `log(0) = -inf` when unseen characters appear during Viterbi decoding.

---

## Reproducibility

| Item | Detail |
|------|--------|
| Random seed | `random.seed(42)` in all synthetic data generators |
| Timer | `time.perf_counter()`, median of 3 repeats |
| Memory | `tracemalloc` peak (bytes to KB) |
| Environment | Conda: `environment.yml` / pip: `requirements.txt` |
| Data | All FASTA files committed to `data/`; synthetic regenerated with fixed seed |
| Benchmarks | `python benchmark/run_benchmarks.py` reproduces all tables and figures |
| Tests | `python tests/test_all.py` exits 0 iff all 46 pass |

---

## Algorithm Complexity Reference

| Algorithm | Time | Space | Optimal? |
|---|---|---|---|
| Needleman-Wunsch | O(nm) | O(nm) | Yes |
| Smith-Waterman | O(nm) | O(nm) | Yes |
| Gotoh | O(nm) | O(nm) | Yes |
| Hirschberg | O(nm) | **O(n+m)** | Yes |
| Banded NW | O(nk) | O(nk) | Yes* |
| BLAST-lite | O(n) avg | O(n) | No |
| Minimizer align | O(n/w) | O(n/w) | No |
| Progressive MSA | O(N^2 nm) | O(Nm) | No |
| Iterative refinement | O(I x N^2 nm) | O(Nm) | No |
| Profile HMM Viterbi | O(T x K) | O(T x K) | Yes+ |
| T-Coffee (optional) | O(N^2 nm) | O(N^2) | No |

\* Optimal within the diagonal band.  + Optimal given the trained model.
N = number of sequences, K = model length, I = refinement iterations.

.
