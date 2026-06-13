"""
Tests TP02 — Alignement global Needleman-Wunsch
NE PAS MODIFIER CE FICHIER.
Lancer avec : uv run pytest test_skeleton.py -v
"""
from __future__ import annotations

import pytest
from needleman_wunsch import Sequence, fill_matrix, traceback, needleman_wunsch


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def compute_score(a1: str, a2: str, match: int, mismatch: int, gap: int) -> int:
    """Calcule le score d'un alignement déjà construit — utilisé pour valider traceback."""
    assert len(a1) == len(a2)
    score = 0
    for c1, c2 in zip(a1, a2):
        if c1 == "-" or c2 == "-":
            score += gap
        elif c1 == c2:
            score += match
        else:
            score += mismatch
    return score


# ---------------------------------------------------------------------------
# fill_matrix — initialisation et valeurs connues
# ---------------------------------------------------------------------------

class TestFillMatrix:
    def test_shape(self):
        seq1 = Sequence("s1", "AG")
        seq2 = Sequence("s2", "ACG")
        m = fill_matrix(seq1, seq2)
        assert len(m) == 3        # len(seq1) + 1
        assert len(m[0]) == 4    # len(seq2) + 1

    def test_origin(self):
        seq1 = Sequence("s1", "AG")
        seq2 = Sequence("s2", "AC")
        m = fill_matrix(seq1, seq2, gap_penalty=-2)
        assert m[0][0] == 0

    def test_first_column(self):
        seq1 = Sequence("s1", "AGC")
        seq2 = Sequence("s2", "A")
        m = fill_matrix(seq1, seq2, gap_penalty=-2)
        assert m[1][0] == -2
        assert m[2][0] == -4
        assert m[3][0] == -6

    def test_first_row(self):
        seq1 = Sequence("s1", "A")
        seq2 = Sequence("s2", "AGC")
        m = fill_matrix(seq1, seq2, gap_penalty=-2)
        assert m[0][1] == -2
        assert m[0][2] == -4
        assert m[0][3] == -6

    def test_known_cell(self):
        # AGCT vs AGT — cellule (2, 2) = alignement "AG" vs "AG" = 2
        seq1 = Sequence("s1", "AGCT")
        seq2 = Sequence("s2", "AGT")
        m = fill_matrix(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-2)
        assert m[2][2] == 2

    def test_final_score(self):
        # AGCT vs AGT → score final = 1
        seq1 = Sequence("s1", "AGCT")
        seq2 = Sequence("s2", "AGT")
        m = fill_matrix(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-2)
        assert m[4][3] == 1


# ---------------------------------------------------------------------------
# needleman_wunsch — score
# ---------------------------------------------------------------------------

class TestNWScore:
    def test_identical_sequences(self):
        seq = Sequence("s", "ATGC")
        score, _, _ = needleman_wunsch(seq, seq, match_score=1, mismatch_score=-1, gap_penalty=-2)
        assert score == 4

    def test_known_example(self):
        # AGCT vs AGT → score = 1
        seq1 = Sequence("s1", "AGCT")
        seq2 = Sequence("s2", "AGT")
        score, _, _ = needleman_wunsch(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-2)
        assert score == 1

    def test_single_match(self):
        score, _, _ = needleman_wunsch(
            Sequence("s", "A"), Sequence("s", "A"),
            match_score=1, mismatch_score=-1, gap_penalty=-2,
        )
        assert score == 1

    def test_single_mismatch(self):
        score, _, _ = needleman_wunsch(
            Sequence("s", "A"), Sequence("s", "G"),
            match_score=1, mismatch_score=-1, gap_penalty=-2,
        )
        assert score == -1

    def test_empty_vs_sequence(self):
        # Aligner "" contre "ATGC" = 4 gaps → score = 4 * (-2) = -8
        seq1 = Sequence("s1", "")
        seq2 = Sequence("s2", "ATGC")
        score, _, _ = needleman_wunsch(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-2)
        assert score == -8

    def test_both_empty(self):
        seq = Sequence("s", "")
        score, a1, a2 = needleman_wunsch(seq, seq)
        assert score == 0
        assert a1 == ""
        assert a2 == ""


# ---------------------------------------------------------------------------
# needleman_wunsch — alignement
# ---------------------------------------------------------------------------

class TestNWAlignment:
    def test_aligned_strings_same_length(self):
        seq1 = Sequence("s1", "AGCT")
        seq2 = Sequence("s2", "AGT")
        _, a1, a2 = needleman_wunsch(seq1, seq2)
        assert len(a1) == len(a2)

    def test_alignment_score_matches_reported_score(self):
        seq1 = Sequence("s1", "AGCT")
        seq2 = Sequence("s2", "AGT")
        score, a1, a2 = needleman_wunsch(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-2)
        assert compute_score(a1, a2, match=1, mismatch=-1, gap=-2) == score

    def test_identical_no_gaps(self):
        seq = Sequence("s", "ATGC")
        _, a1, a2 = needleman_wunsch(seq, seq)
        assert "-" not in a1
        assert "-" not in a2

    def test_alignment_preserves_bases(self):
        # Les bases non-gap dans l'alignement doivent reproduire les séquences originales
        seq1 = Sequence("s1", "AGCT")
        seq2 = Sequence("s2", "AGT")
        _, a1, a2 = needleman_wunsch(seq1, seq2)
        assert a1.replace("-", "") == str(seq1)
        assert a2.replace("-", "") == str(seq2)

    def test_score_consistency_random(self):
        # Vérification score == score(alignement) sur un exemple différent
        seq1 = Sequence("s1", "GCATGCT")
        seq2 = Sequence("s2", "GATTACA")
        score, a1, a2 = needleman_wunsch(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-1)
        assert compute_score(a1, a2, match=1, mismatch=-1, gap=-1) == score
