"""Tests pour src/affine_alignment.py (TP05)."""
import math
import pytest

from sequences import Sequence
from scoring import ScoringMatrix
from affine_alignment import (
    fill_matrices_global,
    fill_matrices_semiglobal,
    find_semiglobal_end,
    traceback_global,
    traceback_semiglobal,
    global_affine,
    semiglobal,
)


# ---------------------------------------------------------------------------
# fill_matrices_global — initialisation
# ---------------------------------------------------------------------------

class TestFillMatricesGlobalInit:
    def test_first_col(self):
        """H[i][0] = gap_open + i*gap_extend pour tout i > 0."""
        seq1 = Sequence("s1", "ACG")
        seq2 = Sequence("s2", "GT")
        H, E, F = fill_matrices_global(seq1, seq2, gap_open=-5, gap_extend=-2)
        assert H[0][0] == 0
        assert H[1][0] == -7   # -5 + 1*(-2)
        assert H[2][0] == -9   # -5 + 2*(-2)
        assert H[3][0] == -11  # -5 + 3*(-2)

    def test_first_row(self):
        """H[0][j] = gap_open + j*gap_extend pour tout j > 0."""
        seq1 = Sequence("s1", "AC")
        seq2 = Sequence("s2", "GTA")
        H, E, F = fill_matrices_global(seq1, seq2, gap_open=-5, gap_extend=-2)
        assert H[0][0] == 0
        assert H[0][1] == -7
        assert H[0][2] == -9
        assert H[0][3] == -11

    def test_e_top_row_is_neg_inf(self):
        """E[0][j] = -inf : impossible d'avoir un gap vertical quand seq1 est vide."""
        seq1 = Sequence("s1", "AC")
        seq2 = Sequence("s2", "GT")
        H, E, F = fill_matrices_global(seq1, seq2, gap_open=-5, gap_extend=-2)
        for j in range(len(seq2) + 1):
            assert math.isinf(E[0][j]) and E[0][j] < 0, \
                f"E[0][{j}] devrait être -inf, pas {E[0][j]}"

    def test_f_first_col_is_neg_inf(self):
        """F[i][0] = -inf : impossible d'avoir un gap horizontal quand seq2 est vide."""
        seq1 = Sequence("s1", "AC")
        seq2 = Sequence("s2", "GT")
        H, E, F = fill_matrices_global(seq1, seq2, gap_open=-5, gap_extend=-2)
        for i in range(len(seq1) + 1):
            assert math.isinf(F[i][0]) and F[i][0] < 0, \
                f"F[{i}][0] devrait être -inf, pas {F[i][0]}"


# ---------------------------------------------------------------------------
# fill_matrices_semiglobal — initialisation
# ---------------------------------------------------------------------------

class TestFillMatricesSemiglobalInit:
    def test_top_row_all_zeros(self):
        """H[0][j] = 0 pour tout j : pas de pénalité pour les préfixes de seq2."""
        seq1 = Sequence("s1", "ACG")
        seq2 = Sequence("s2", "GTATCG")
        H, E, F = fill_matrices_semiglobal(seq1, seq2, gap_open=-5, gap_extend=-2)
        for j in range(len(seq2) + 1):
            assert H[0][j] == 0, \
                f"H[0][{j}] devrait être 0 en semi-global, pas {H[0][j]}"

    def test_first_col_still_penalized(self):
        """H[i][0] reste pénalisé : les gaps dans seq1 (le read) coûtent toujours."""
        seq1 = Sequence("s1", "ACG")
        seq2 = Sequence("s2", "GT")
        H, E, F = fill_matrices_semiglobal(seq1, seq2, gap_open=-5, gap_extend=-2)
        assert H[1][0] == -7   # identique à global
        assert H[2][0] == -9
        assert H[3][0] == -11


# ---------------------------------------------------------------------------
# global_affine — intégration
# ---------------------------------------------------------------------------

class TestGlobalAffine:
    def test_perfect_match(self):
        """Séquences identiques : score = n * match_score, aucun gap."""
        seq1 = Sequence("s1", "ACGT")
        seq2 = Sequence("s2", "ACGT")
        score, al1, al2 = global_affine(
            seq1, seq2, gap_open=-5, gap_extend=-1, match_score=2, mismatch_score=-3
        )
        assert score == 8
        assert al1 == "ACGT"
        assert al2 == "ACGT"

    def test_single_gap(self):
        """
        seq1="ACT", seq2="ACGT" : seul alignement optimal = AC-T / ACGT.
        Score = match(A,A) + match(C,C) + gap(1) + match(T,T)
              = 2 + 2 + (-3-1) + 2 = 2.
        (mismatch=-10 rend tout mismatch pire que le gap)
        """
        seq1 = Sequence("s1", "ACT")
        seq2 = Sequence("s2", "ACGT")
        score, al1, al2 = global_affine(
            seq1, seq2, gap_open=-3, gap_extend=-1, match_score=2, mismatch_score=-10
        )
        assert score == 2
        assert al1 == "AC-T"
        assert al2 == "ACGT"

    def test_affine_gap_score_formula(self):
        """
        seq1="AAAGGG", seq2="AAACCCGGG" : un gap continu de 3 dans seq1 pour les CCC.
        Score = 3*match + (gap_open + 3*gap_extend) + 3*match.
        Vérifie que la formule affine gap(k) = gap_open + k*gap_extend est correcte.
        """
        seq1 = Sequence("s1", "AAAGGG")
        seq2 = Sequence("s2", "AAACCCGGG")
        gap_open, gap_extend, match = -10, -1, 5
        score, al1, al2 = global_affine(
            seq1, seq2, gap_open=gap_open, gap_extend=gap_extend,
            match_score=match, mismatch_score=-20
        )
        expected = 3 * match + (gap_open + 3 * gap_extend) + 3 * match
        assert score == expected

    def test_aligned_same_length(self):
        """Les deux chaînes alignées ont toujours la même longueur."""
        seq1 = Sequence("s1", "ACGTTGCA")
        seq2 = Sequence("s2", "ACGT")
        _, al1, al2 = global_affine(seq1, seq2, gap_open=-5, gap_extend=-1)
        assert len(al1) == len(al2)

    def test_with_scoring_matrix(self):
        """Intégration avec ScoringMatrix.identity_dna() : chaque match vaut +1."""
        sm = ScoringMatrix.identity_dna()
        seq1 = Sequence("s1", "ACGT")
        seq2 = Sequence("s2", "ACGT")
        score, al1, al2 = global_affine(
            seq1, seq2, gap_open=-5, gap_extend=-1, scoring_matrix=sm
        )
        assert score == 4
        assert al1 == "ACGT"
        assert al2 == "ACGT"


# ---------------------------------------------------------------------------
# semiglobal — intégration
# ---------------------------------------------------------------------------

class TestSemiglobal:
    def test_query_in_middle(self):
        """Le read s'aligne parfaitement au milieu du génome."""
        seq1 = Sequence("s1", "CGT")
        seq2 = Sequence("s2", "AAACGTCC")
        score, al1, al2 = semiglobal(
            seq1, seq2, gap_open=-5, gap_extend=-1, match_score=2, mismatch_score=-3
        )
        assert score == 6
        assert al1 == "CGT"
        assert al2 == "CGT"

    def test_query_at_start(self):
        """Le read s'aligne parfaitement au début du génome."""
        seq1 = Sequence("s1", "ACG")
        seq2 = Sequence("s2", "ACGTTT")
        score, al1, al2 = semiglobal(
            seq1, seq2, gap_open=-5, gap_extend=-1, match_score=2, mismatch_score=-3
        )
        assert score == 6
        assert al1 == "ACG"
        assert al2 == "ACG"

    def test_query_at_end(self):
        """Le read s'aligne parfaitement à la fin du génome."""
        seq1 = Sequence("s1", "CGT")
        seq2 = Sequence("s2", "AAACGT")
        score, al1, al2 = semiglobal(
            seq1, seq2, gap_open=-5, gap_extend=-1, match_score=2, mismatch_score=-3
        )
        assert score == 6
        assert al1 == "CGT"
        assert al2 == "CGT"

    def test_semiglobal_better_than_global(self):
        """
        Semi-global >> global quand le génome a de longs flancs non appariés.
        seq1="ACGT" dans seq2="TTTACGTGGG" :
          semi-global : 4*match = 8  (flancs gratuits)
          global      : 8 + 2*(gap_open + 3*gap_extend) = 8 + 2*(-8) = -8
        """
        seq1 = Sequence("s1", "ACGT")
        seq2 = Sequence("s2", "TTTACGTGGG")
        gap_open, gap_extend, match = -5, -1, 2
        score_sg, _, _ = semiglobal(
            seq1, seq2, gap_open=gap_open, gap_extend=gap_extend,
            match_score=match, mismatch_score=-10
        )
        score_gl, _, _ = global_affine(
            seq1, seq2, gap_open=gap_open, gap_extend=gap_extend,
            match_score=match, mismatch_score=-10
        )
        assert score_sg > score_gl
        assert score_sg == 4 * match

    def test_semiglobal_aligned_same_length(self):
        """Les deux chaînes retournées par semiglobal ont la même longueur."""
        seq1 = Sequence("s1", "ACGT")
        seq2 = Sequence("s2", "GGGACGTAAA")
        _, al1, al2 = semiglobal(
            seq1, seq2, gap_open=-5, gap_extend=-1, match_score=2, mismatch_score=-3
        )
        assert len(al1) == len(al2)

    def test_find_semiglobal_end(self):
        """
        find_semiglobal_end retourne la colonne et le score max de la dernière ligne.
        seq1="CGT" s'aligne sur seq2="AAACGTCC" aux positions 3–5 (indices 1-based : j=4..6).
        La dernière ligne H[3][*] doit être maximale à j=6, avec score=6.
        """
        seq1 = Sequence("s1", "CGT")
        seq2 = Sequence("s2", "AAACGTCC")
        H, _, _ = fill_matrices_semiglobal(
            seq1, seq2, gap_open=-5, gap_extend=-1, match_score=2, mismatch_score=-3
        )
        end_col, score = find_semiglobal_end(H)
        assert score == 6
        assert end_col == 6   # seq2[3..5]="CGT" → alignement finit à j=6
