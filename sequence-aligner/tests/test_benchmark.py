"""Tests pour src/benchmark.py (TP06)."""
import time as _time

import pytest

from sequences import Sequence
from benchmark import compare_with_biopython, generate_random_dna, time_function


# ---------------------------------------------------------------------------
# generate_random_dna
# ---------------------------------------------------------------------------

class TestGenerateRandomDna:
    def test_correct_length(self):
        seq = generate_random_dna(100)
        assert len(seq) == 100

    def test_valid_alphabet(self):
        """Toutes les bases générées appartiennent à {A, C, G, T}."""
        seq = generate_random_dna(500, seed=0)
        for base in str(seq):
            assert base in "ACGT", f"Base invalide : {base!r}"

    def test_reproducible_with_seed(self):
        """Deux appels avec la même graine produisent la même séquence."""
        seq1 = generate_random_dna(50, seed=42)
        seq2 = generate_random_dna(50, seed=42)
        assert str(seq1) == str(seq2)

    def test_different_seeds_differ(self):
        """Deux graines différentes produisent (très probablement) des séquences différentes."""
        seq1 = generate_random_dna(100, seed=1)
        seq2 = generate_random_dna(100, seed=2)
        assert str(seq1) != str(seq2)


# ---------------------------------------------------------------------------
# time_function
# ---------------------------------------------------------------------------

class TestTimeFunction:
    def test_returns_non_negative_float(self):
        t = time_function(lambda: None, n_repeat=3)
        assert isinstance(t, float)
        assert t >= 0.0

    def test_slower_function_takes_longer(self):
        """Une fonction avec sleep mesurée > une fonction immédiate."""
        def fast():
            pass

        def slow():
            _time.sleep(0.02)

        t_fast = time_function(fast, n_repeat=3)
        t_slow = time_function(slow, n_repeat=3)
        assert t_slow > t_fast

    def test_args_forwarded(self):
        """Les arguments sont bien transmis à la fonction mesurée."""
        results = []

        def record(x, y):
            results.append(x + y)

        time_function(record, args=(1, 2), n_repeat=3)
        assert all(v == 3 for v in results)


# ---------------------------------------------------------------------------
# compare_with_biopython
# ---------------------------------------------------------------------------

class TestCompareWithBiopython:
    def test_perfect_match_scores_agree(self):
        """Sur un alignement parfait, nos scores et BioPython coïncident."""
        seq1 = Sequence("s1", "ACGT")
        seq2 = Sequence("s2", "ACGT")
        result = compare_with_biopython(
            seq1, seq2, gap_open=-5, gap_extend=-1, match_score=2, mismatch_score=-3
        )
        assert result["match"] is True
        assert result["our_score"] == result["biopython_score"]

    def test_result_has_required_keys(self):
        """Le dict retourné contient bien les trois clés attendues."""
        seq1 = Sequence("s1", "ACGT")
        seq2 = Sequence("s2", "ACGT")
        result = compare_with_biopython(seq1, seq2)
        assert "our_score" in result
        assert "biopython_score" in result
        assert "match" in result

    def test_with_gap_scores_agree(self):
        """Sur un alignement avec gap, les scores restent identiques après ajustement de convention."""
        seq1 = Sequence("s1", "AAAGGG")
        seq2 = Sequence("s2", "AAACCCGGG")
        result = compare_with_biopython(
            seq1, seq2, gap_open=-10, gap_extend=-1, match_score=5, mismatch_score=-20
        )
        assert result["match"] is True
