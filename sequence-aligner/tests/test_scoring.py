"""Tests pour src/scoring.py (TP03)."""
from pathlib import Path

import pytest

from scoring import ScoringMatrix
from needleman_wunsch import needleman_wunsch
from sequences import Sequence

DATA_DIR = Path(__file__).parent.parent / "data"


# ---------------------------------------------------------------------------
# ScoringMatrix.__init__ et score()
# ---------------------------------------------------------------------------

def test_init_and_score_match():
    sm = ScoringMatrix({("A", "A"): 1, ("A", "T"): -1, ("T", "A"): -1, ("T", "T"): 1})
    assert sm.score("A", "A") == 1
    assert sm.score("T", "T") == 1


def test_score_mismatch():
    sm = ScoringMatrix({("A", "A"): 1, ("A", "T"): -1, ("T", "A"): -1, ("T", "T"): 1})
    assert sm.score("A", "T") == -1
    assert sm.score("T", "A") == -1


def test_score_symmetry():
    sm = ScoringMatrix.nuc_simple()
    for a in "ATGC":
        for b in "ATGC":
            assert sm.score(a, b) == sm.score(b, a), f"score({a},{b}) != score({b},{a})"


def test_score_unknown_raises():
    sm = ScoringMatrix({("A", "A"): 1})
    with pytest.raises(KeyError):
        sm.score("A", "T")


# ---------------------------------------------------------------------------
# ScoringMatrix.identity_dna()
# ---------------------------------------------------------------------------

def test_identity_dna_matches():
    sm = ScoringMatrix.identity_dna()
    for base in "ATGC":
        assert sm.score(base, base) == 1


def test_identity_dna_mismatches():
    sm = ScoringMatrix.identity_dna()
    assert sm.score("A", "T") == -1
    assert sm.score("G", "C") == -1
    assert sm.score("A", "G") == -1


# ---------------------------------------------------------------------------
# ScoringMatrix.nuc_simple()
# ---------------------------------------------------------------------------

def test_nuc_simple_matches():
    sm = ScoringMatrix.nuc_simple()
    for base in "ATGC":
        assert sm.score(base, base) == 1


def test_nuc_simple_transitions():
    # A↔G (deux purines) et C↔T (deux pyrimidines) → transitions, pénalité -1
    sm = ScoringMatrix.nuc_simple()
    assert sm.score("A", "G") == -1
    assert sm.score("G", "A") == -1
    assert sm.score("C", "T") == -1
    assert sm.score("T", "C") == -1


def test_nuc_simple_transversions():
    # A↔C, A↔T, G↔C, G↔T → transversions, pénalité -3
    sm = ScoringMatrix.nuc_simple()
    assert sm.score("A", "C") == -3
    assert sm.score("A", "T") == -3
    assert sm.score("G", "C") == -3
    assert sm.score("G", "T") == -3
    assert sm.score("C", "A") == -3
    assert sm.score("T", "A") == -3


# ---------------------------------------------------------------------------
# ScoringMatrix.from_file()
# ---------------------------------------------------------------------------

def test_from_file_loads_correctly():
    sm = ScoringMatrix.from_file(str(DATA_DIR / "nuc_simple.mat"))
    # Matches
    for base in "ATGC":
        assert sm.score(base, base) == 1
    # Transitions
    assert sm.score("A", "G") == -1
    assert sm.score("C", "T") == -1
    # Transversions
    assert sm.score("A", "C") == -3
    assert sm.score("G", "T") == -3


def test_from_file_not_found():
    with pytest.raises(FileNotFoundError):
        ScoringMatrix.from_file(str(DATA_DIR / "does_not_exist.mat"))


# ---------------------------------------------------------------------------
# Intégration : needleman_wunsch avec scoring_matrix
# ---------------------------------------------------------------------------

def test_nw_with_identity_gives_same_result_as_simple():
    """identity_dna() avec gap=-2 doit reproduire le comportement de NW par défaut."""
    seq1 = Sequence("s1", "AGCT")
    seq2 = Sequence("s2", "AGT")
    sm = ScoringMatrix.identity_dna()
    score_sm, al1_sm, al2_sm = needleman_wunsch(seq1, seq2, gap_penalty=-2, scoring_matrix=sm)
    score_def, al1_def, al2_def = needleman_wunsch(seq1, seq2, gap_penalty=-2)
    assert score_sm == score_def
    assert al1_sm == al1_def
    assert al2_sm == al2_def


def test_nw_nuc_simple_changes_score():
    """Deux transversions coûtent plus cher que deux gaps — le score doit changer."""
    # seq1 = "AT", seq2 = "CA"
    # Avec simple (match=1, mismatch=-1, gap=-2) : "AT" vs "CA", score = -1 + -1 = -2
    # Avec nuc_simple (transversion=-3, gap=-2) : -3 + -3 = -6 vs -2 + 1 + -2 = -3
    #   → l'algo préfère introduire des gaps (-AT vs CA-), score = -3
    seq1 = Sequence("s1", "AT")
    seq2 = Sequence("s2", "CA")
    sm = ScoringMatrix.nuc_simple()
    score_sm, _, _ = needleman_wunsch(seq1, seq2, gap_penalty=-2, scoring_matrix=sm)
    score_def, _, _ = needleman_wunsch(seq1, seq2, gap_penalty=-2)
    assert score_sm != score_def


def test_nw_nuc_simple_changes_alignment():
    """Vérifier que l'alignement optimal change avec nuc_simple."""
    seq1 = Sequence("s1", "AT")
    seq2 = Sequence("s2", "CA")
    sm = ScoringMatrix.nuc_simple()
    score_sm, al1, al2 = needleman_wunsch(seq1, seq2, gap_penalty=-2, scoring_matrix=sm)
    # Avec nuc_simple, deux transversions (-3 chacune) sont plus coûteuses qu'un gap (-2)
    # L'alignement optimal introduit des gaps plutôt que de subir deux transversions
    assert "-" in al1 or "-" in al2  # au moins un gap dans l'alignement
    assert score_sm == -3
