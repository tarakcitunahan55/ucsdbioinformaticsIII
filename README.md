# UCSD Bioinformatics III — Comparing Genes, Proteins & Genomes

Python implementations of core dynamic programming and sequence alignment algorithms from the Comparing Genes, Proteins, and Genomes course.

## Contents

### `manhattan_problem_seq_align_week_1/` — Dynamic Programming Foundations
| File | Description |
|---|---|
| `money_change_recursive.py` / `dynamic_program_money_change.py` / `dp_change_with_coins.py` / `dp_change_space_efficient.py` | Change-making problem: recursive, DP, and space-efficient variants |
| `manhattan_recursion.py` / `manhattan_dynamic.py` | Manhattan Tourist problem via recursion and DP |
| `longest_path_dag.py` | Longest path in a Directed Acyclic Graph |
| `hanoi_recursion.py` | Tower of Hanoi (recursion warm-up) |
| `backtrack_seq_align_rec_lcs.py` / `backtrack_seq_align_rec_all_lcs.py` / `backtrack_seq_align_iterate_lcs.py` | Longest Common Subsequence via recursive and iterative backtracking |

### `global_to_local_alignment_week_2/` — Core Sequence Alignment
| File | Description |
|---|---|
| `global_align_penalty_recurse.py` / `global_align_penalty_iterate.py` / `global_all_align_penalty_recurse.py` | Global alignment with scoring/gap penalties (recursive & iterative) |
| `local_align_penalty_recurse.py` / `local_all_align_penalty_recurse.py` | Local (Smith-Waterman style) alignment |
| `fitting_alignment.py` / `fitting_alignment_all.py` | Fitting alignment (one sequence fully within another) |
| `overlap_align.py` / `overlap_align_all.py` | Overlap alignment between sequence ends |
| `edit_distance.py` | Computes edit (Levenshtein) distance between two sequences |

### `affine_gap_penalty_week_3/` — Affine Gaps & Multiple Alignment
| File | Description |
|---|---|
| `affine_gap_penalty_global.py` / `affine_gap_penalty_global_all.py` | Global alignment with affine gap penalties |
| `affine_gap_penalty_local.py` | Local alignment with affine gap penalties |
| `multiple_align_constant_gap_penalty.py` | Multiple sequence alignment, constant/linear gap penalty |
| `multiple_align_affine_gap_penalty.py` | Multiple sequence alignment, affine gap penalty |
| `multiple_align_lcs.py` | Multi-sequence LCS for conserved region detection, no gap penalty |

## Topics covered
Dynamic programming fundamentals · global/local/fitting/overlap alignment · affine gap penalties · multiple sequence alignment · edit distance

## Notes
Coursework implementations; shared for portfolio purposes. DO NOT COPY.
