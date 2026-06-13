"""Tests pour src/smith_waterman.py (TP04)."""
import pytest

from sequences import Sequence
from scoring import ScoringMatrix
from smith_waterman import fill_matrix, find_max, traceback, smith_waterman
from needleman_wunsch import needleman_wunsch


# ---------------------------------------------------------------------------
# fill_matrix
# ---------------------------------------------------------------------------

class TestFillMatrix:
    def test_border_all_zeros(self):
        """Contrairement à NW, la première ligne et colonne sont toutes à 0."""
        seq1 = Sequence("s1", "ACG")
        seq2 = Sequence("s2", "GT")
        F = fill_matrix(seq1, seq2)
        for i in range(len(seq1) + 1):
            assert F[i][0] == 0, f"F[{i}][0] doit être 0, pas {F[i][0]}"
        for j in range(len(seq2) + 1):
            assert F[0][j] == 0, f"F[0][{j}] doit être 0, pas {F[0][j]}"

    def test_floor_at_zero(self):
        """Aucune cellule ne peut être négative."""
        seq1 = Sequence("s1", "AAAA")
        seq2 = Sequence("s2", "TTTT")
        F = fill_matrix(seq1, seq2)
        for row in F:
            for val in row:
                assert val >= 0, f"Valeur négative trouvée : {val}"

    def test_known_cells(self):
        """
        seq1 = "ACGT", seq2 = "CGT", match=1, mismatch=-1, gap=-2
        F[2][1] = 1  (C↔C match depuis F[1][0]=0)
        F[3][2] = 2  (G↔G match depuis F[2][1]=1)
        F[4][3] = 3  (T↔T match depuis F[3][2]=2)
        """
        seq1 = Sequence("s1", "ACGT")
        seq2 = Sequence("s2", "CGT")
        F = fill_matrix(seq1, seq2)
        assert F[2][1] == 1
        assert F[3][2] == 2
        assert F[4][3] == 3


# ---------------------------------------------------------------------------
# find_max
# ---------------------------------------------------------------------------

class TestFindMax:
    def test_finds_maximum(self):
        matrix = [
            [0, 0, 0],
            [0, 1, 2],
            [0, 3, 1],
        ]
        row, col, score = find_max(matrix)
        assert score == 3
        assert row == 2
        assert col == 1

    def test_tie_returns_first(self):
        matrix = [
            [0, 0, 0],
            [0, 5, 5],
            [0, 5, 5],
        ]
        row, col, score = find_max(matrix)
        assert score == 5
        assert (row, col) == (1, 1)

    def test_all_zeros(self):
        matrix = [[0, 0], [0, 0]]
        row, col, score = find_max(matrix)
        assert score == 0


# ---------------------------------------------------------------------------
# smith_waterman — intégration
# ---------------------------------------------------------------------------

class TestSmithWaterman:
    def test_perfect_local_match(self):
        """seq2 est entièrement contenu dans seq1 → score = len(seq2)."""
        seq1 = Sequence("s1", "ACGT")
        seq2 = Sequence("s2", "CGT")
        score, al1, al2 = smith_waterman(seq1, seq2)
        assert score == 3
        assert al1 == "CGT"
        assert al2 == "CGT"

    def test_local_in_longer_sequence(self):
        """SW trouve le motif local même entouré de bruit."""
        seq1 = Sequence("s1", "TTTACGT")
        seq2 = Sequence("s2", "ACGT")
        score, al1, al2 = smith_waterman(seq1, seq2)
        assert score == 4
        assert al1 == "ACGT"
        assert al2 == "ACGT"

    def test_no_positive_alignment(self):
        """Quand tous les appariements sont négatifs, score = 0 et alignement vide."""
        seq1 = Sequence("s1", "AAAA")
        seq2 = Sequence("s2", "TTTT")
        score, al1, al2 = smith_waterman(seq1, seq2, match_score=1, mismatch_score=-5, gap_penalty=-5)
        assert score == 0
        assert al1 == ""
        assert al2 == ""

    def test_identical_sequences(self):
        seq1 = Sequence("s1", "AGCT")
        seq2 = Sequence("s2", "AGCT")
        score, al1, al2 = smith_waterman(seq1, seq2)
        assert score == 4
        assert al1 == "AGCT"
        assert al2 == "AGCT"

    def test_sw_score_better_than_nw_for_local_motif(self):
        """
        Quand seq2 est entouré de bruit dans seq1, SW > NW.
        seq1 = "GGGACGT", seq2 = "ACGT"
        NW aligne la totalité : GGG↔--- (-6) + ACGT↔ACGT (+4) = -2
        SW trouve uniquement ACGT : score = 4
        """
        seq1 = Sequence("s1", "GGGACGT")
        seq2 = Sequence("s2", "ACGT")
        score_sw, _, _ = smith_waterman(seq1, seq2)
        score_nw, _, _ = needleman_wunsch(seq1, seq2)
        assert score_sw > score_nw
        assert score_sw == 4

    def test_aligned_strings_same_length(self):
        seq1 = Sequence("s1", "ACGTTGCA")
        seq2 = Sequence("s2", "CGT")
        _, al1, al2 = smith_waterman(seq1, seq2)
        assert len(al1) == len(al2)

    def test_with_scoring_matrix(self):
        """Avec NUC_SIMPLE, une transversion coûte -3 → SW évite les mauvais appariements."""
        seq1 = Sequence("s1", "ACGT")
        seq2 = Sequence("s2", "ACGT")
        sm = ScoringMatrix.nuc_simple()
        score, al1, al2 = smith_waterman(seq1, seq2, scoring_matrix=sm)
        assert score == 4
        assert al1 == "ACGT"
        assert al2 == "ACGT"
