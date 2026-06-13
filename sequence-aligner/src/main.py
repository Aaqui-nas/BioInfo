from sequences import Sequence, read_fasta
from affine_alignment import global_affine, semiglobal
from scoring import ScoringMatrix
from visualize import format_alignment, alignment_stats

from benchmark import time_function, generate_random_dna

from Bio.Align import PairwiseAligner, substitution_matrices


# ============================================================
# UTILS BIOPYTHON MATRIX
# ============================================================
def build_biopython_matrix(scoring_matrix: ScoringMatrix):
    alphabet = ("A", "T", "G", "C")

    matrix = substitution_matrices.Array(alphabet, dims=2)

    for a in alphabet:
        for b in alphabet:
            matrix[a, b] = scoring_matrix.score(a, b)

    return matrix


# ============================================================
# BIOPYTHON ALIGNER FACTORY (reused)
# ============================================================
def make_biopython_aligner(gap_open, gap_extend, scoring_matrix=None, semiglobal=False):
    aligner = PairwiseAligner()
    aligner.mode = "global"

    aligner.open_gap_score = gap_open + gap_extend
    aligner.extend_gap_score = gap_extend

    if semiglobal:
        aligner.left_deletion_score = 0
        aligner.right_deletion_score = 0
        aligner.left_insertion_score = 0
        aligner.right_insertion_score = 0

    if scoring_matrix is None:
        aligner.match_score = 1
        aligner.mismatch_score = -1
    else:
        aligner.substitution_matrix = build_biopython_matrix(scoring_matrix)

    return aligner


# ============================================================
# SCORE BIOPYTHON
# ============================================================
def biopython_score(seq1, seq2, gap_open=-5, gap_extend=-1, scoring_matrix=None, mode="global"):
    aligner = make_biopython_aligner(
        gap_open,
        gap_extend,
        scoring_matrix=scoring_matrix,
        mode=mode,
    )
    return aligner.score(str(seq1), str(seq2))


# ============================================================
# TIMING BIOPYTHON GLOBAL
# ============================================================
def time_biopython_global(seq1, seq2, gap_open=-5, gap_extend=-1, scoring_matrix=None):
    aligner = make_biopython_aligner(
        gap_open,
        gap_extend,
        scoring_matrix=scoring_matrix,
        mode="global",
    )

    return time_function(lambda: aligner.score(str(seq1), str(seq2)), n_repeat=5)


# ============================================================
# TIMING BIOPYTHON SEMI-GLOBAL
# ============================================================
def time_biopython_semiglobal(seq1, seq2, gap_open=-5, gap_extend=-1, scoring_matrix=None):
    aligner = make_biopython_aligner(
        gap_open,
        gap_extend,
        scoring_matrix=scoring_matrix,
        mode="semiglobal",
    )

    return time_function(lambda: aligner.score(str(seq1), str(seq2)), n_repeat=5)


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":

    # =========================
    # LOAD DATA
    # =========================
    human_seq = read_fasta("data/human.fna")[0]
    mouse_seq = read_fasta("data/mouse.fna")[0]

    scoring = ScoringMatrix.from_file("data/nuc_simple.mat")

    # ========================================================
    # GLOBAL BENCHMARK
    # ========================================================
    print("=== GLOBAL ALIGNMENT ===")

    t_ours_global = time_function(
        global_affine,
        args=(human_seq, mouse_seq),
        kwargs={"scoring_matrix": scoring},
        n_repeat=5,
    )

    t_bio_global = time_biopython_global(
        human_seq,
        mouse_seq,
        scoring_matrix=scoring,
    )

    print(f"Our global   : {t_ours_global:.6f} s")
    print(f"Bio global   : {t_bio_global:.6f} s")
    print(f"Speedup      : {t_bio_global / t_ours_global:.2f}x")

    # =========================
    # SCORE COMPARISON GLOBAL
    # =========================
    ours_score_global, h_aln, m_aln = global_affine(human_seq, mouse_seq, scoring_matrix=scoring)

    bio_score_global = biopython_score(
        human_seq,
        mouse_seq,
        scoring_matrix=scoring,
        mode="global",
    )

    print("\n=== SCORE COMPARISON (GLOBAL) ===")
    print(f"Our score : {ours_score_global}")
    print(f"Bio score : {bio_score_global}")
    print(f"Delta     : {abs(ours_score_global - bio_score_global)}")

    stats = alignment_stats(h_aln, m_aln)

    print("\n=== STATS ===")
    for k, v in stats.items():
        print(f"{k}: {v}")

    # ========================================================
    # SEMI-GLOBAL BENCHMARK
    # ========================================================
    print("\n=== SEMI-GLOBAL ALIGNMENT ===")

    t_ours_semi = time_function(
        semiglobal,
        args=(human_seq, mouse_seq),
        kwargs={"scoring_matrix": scoring},
        n_repeat=5,
    )

    t_bio_semi = time_biopython_semiglobal(
        human_seq,
        mouse_seq,
        scoring_matrix=scoring,
    )

    print(f"Our semi     : {t_ours_semi:.6f} s")
    print(f"Bio semi     : {t_bio_semi:.6f} s")
    print(f"Speedup      : {t_bio_semi / t_ours_semi:.2f}x")

    # =========================
    # SCORE COMPARISON SEMI
    # =========================
    ours_score_semi, h_aln, m_aln = semiglobal(human_seq, mouse_seq, scoring_matrix=scoring)

    bio_score_semi = biopython_score(
        human_seq,
        mouse_seq,
        scoring_matrix=scoring,
        mode="semiglobal",
    )

    print("\n=== SCORE COMPARISON (SEMI-GLOBAL) ===")
    print(f"Our score : {ours_score_semi}")
    print(f"Bio score : {bio_score_semi}")
    print(f"Delta     : {abs(ours_score_semi - bio_score_semi)}")

    # ========================================================
    # ALIGNMENT VISUALIZATION (SEMI-GLOBAL)
    # ========================================================
    print("\n=== ALIGNMENT VISUALIZATION ===")


    print(format_alignment(h_aln, m_aln, "human", "mouse", line_width=80))

    stats = alignment_stats(h_aln, m_aln)

    print("\n=== STATS ===")
    for k, v in stats.items():
        print(f"{k}: {v}")
