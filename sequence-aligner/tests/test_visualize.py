"""Tests pour src/visualize.py (TP06)."""
import pytest

from visualize import alignment_stats, format_alignment


# ---------------------------------------------------------------------------
# alignment_stats
# ---------------------------------------------------------------------------

class TestAlignmentStats:
    def test_perfect_match(self):
        stats = alignment_stats("ACGT", "ACGT")
        assert stats["matches"] == 4
        assert stats["mismatches"] == 0
        assert stats["gaps"] == 0
        assert stats["identity"] == 1.0

    def test_all_mismatches(self):
        stats = alignment_stats("AAAA", "TTTT")
        assert stats["matches"] == 0
        assert stats["mismatches"] == 4
        assert stats["gaps"] == 0
        assert stats["identity"] == 0.0

    def test_single_gap_in_seq1(self):
        """'-' dans seq1 → gap, même si seq2 a une base normale."""
        stats = alignment_stats("AC-T", "ACGT")
        assert stats["matches"] == 3
        assert stats["mismatches"] == 0
        assert stats["gaps"] == 1
        assert stats["length"] == 4
        assert abs(stats["identity"] - 0.75) < 1e-9

    def test_single_gap_in_seq2(self):
        """'-' dans seq2 → gap, même si seq1 a une base normale."""
        stats = alignment_stats("ACGT", "ACG-")
        assert stats["gaps"] == 1
        assert stats["matches"] == 3

    def test_mismatch_vs_gap(self):
        """Un mismatch ne compte PAS comme gap."""
        stats = alignment_stats("AACT", "AGCT")
        assert stats["matches"] == 3
        assert stats["mismatches"] == 1
        assert stats["gaps"] == 0
        assert abs(stats["identity"] - 0.75) < 1e-9

    def test_empty_sequences(self):
        stats = alignment_stats("", "")
        assert stats["length"] == 0
        assert stats["matches"] == 0
        assert stats["identity"] == 0.0

    def test_length_matches_aligned_length(self):
        stats = alignment_stats("AC-T", "ACGT")
        assert stats["length"] == len("AC-T")


# ---------------------------------------------------------------------------
# format_alignment
# ---------------------------------------------------------------------------

class TestFormatAlignment:
    def test_contains_both_sequences(self):
        """Les deux séquences apparaissent dans la sortie."""
        result = format_alignment("ACGT", "ACGT", "seq1", "seq2")
        assert "ACGT" in result
        assert "seq1" in result
        assert "seq2" in result

    def test_perfect_match_shows_pipes(self):
        """Un alignement parfait produit '||||' dans la ligne centrale."""
        result = format_alignment("ACGT", "ACGT")
        assert "||||" in result

    def test_gap_breaks_middle_line(self):
        """Un gap dans seq1 produit un espace dans la ligne centrale à cette position."""
        result = format_alignment("AC-T", "ACGT")
        assert "|| |" in result

    def test_mismatch_breaks_middle_line(self):
        """Un mismatch produit un espace dans la ligne centrale."""
        result_match = format_alignment("A", "A")
        result_mismatch = format_alignment("A", "T")
        assert "|" in result_match
        assert "|" not in result_mismatch

    def test_wraps_at_line_width(self):
        """Un alignement de 80 chars avec line_width=60 produit au moins 2 blocs."""
        al1 = "ACGT" * 20   # 80 chars
        al2 = "ACGT" * 20
        result = format_alignment(al1, al2, line_width=60)
        non_empty_lines = [l for l in result.split("\n") if l.strip()]
        assert len(non_empty_lines) >= 6   # au moins 2 blocs × 3 lignes

    def test_returns_string(self):
        result = format_alignment("ACGT", "ACGT")
        assert isinstance(result, str)

    def test_single_char(self):
        """Un alignement d'une seule base fonctionne sans erreur."""
        result = format_alignment("A", "A")
        assert "|" in result

    def test_all_gaps(self):
        """Alignment de gaps purs → aucun '|' dans la sortie."""
        result = format_alignment("----", "ACGT")
        assert "|" not in result
