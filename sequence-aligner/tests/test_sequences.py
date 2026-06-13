"""
Tests TP01 — Représentation des séquences + parsing FASTA.
NE PAS MODIFIER CE FICHIER.
Lancer avec : uv run pytest test_skeleton.py -v
"""
import pytest
from sequences import Sequence, SequenceType, read_fasta, write_fasta


# ---------------------------------------------------------------------------
# Sequence — création et validation
# ---------------------------------------------------------------------------

class TestSequenceCreation:
    def test_basic_creation(self):
        seq = Sequence("s1", "ATGC")
        assert str(seq) == "ATGC"
        assert len(seq) == 4

    def test_lowercase_normalized(self):
        seq = Sequence("s1", "atgc")
        assert str(seq) == "ATGC"

    def test_empty_sequence(self):
        seq = Sequence("s1", "")
        assert len(seq) == 0

    def test_invalid_dna_character(self):
        with pytest.raises(ValueError):
            Sequence("s1", "ATGX")

    def test_invalid_rna_character(self):
        with pytest.raises(ValueError):
            Sequence("s1", "AUGX", seq_type=SequenceType.RNA)

    def test_dna_contains_u_is_invalid(self):
        with pytest.raises(ValueError):
            Sequence("s1", "AUGC", seq_type=SequenceType.DNA)

    def test_equality(self):
        assert Sequence("s1", "ATGC") == Sequence("s1", "ATGC")
        assert Sequence("s1", "ATGC") != Sequence("s2", "ATGC")


# ---------------------------------------------------------------------------
# Complement & reverse complement
# ---------------------------------------------------------------------------

class TestComplement:
    def test_dna_complement(self):
        assert Sequence("s", "ATGC").complement() == "TACG"

    def test_dna_complement_all_bases(self):
        assert Sequence("s", "AATTGGCC").complement() == "TTAACCGG"

    def test_reverse_complement(self):
        assert Sequence("s", "ATGC").reverse_complement() == "GCAT"

    def test_palindrome_ecori(self):
        # GAATTC est un palindrome (site de restriction EcoRI)
        seq = Sequence("ecori", "GAATTC")
        assert seq.reverse_complement() == "GAATTC"

    def test_rna_complement(self):
        seq = Sequence("s", "AUGC", seq_type=SequenceType.RNA)
        assert seq.complement() == "UACG"


# ---------------------------------------------------------------------------
# GC content
# ---------------------------------------------------------------------------

class TestGCContent:
    def test_half_gc(self):
        assert Sequence("s", "ATGC").gc_content() == pytest.approx(0.5)

    def test_all_gc(self):
        assert Sequence("s", "GGCC").gc_content() == pytest.approx(1.0)

    def test_no_gc(self):
        assert Sequence("s", "AATT").gc_content() == pytest.approx(0.0)

    def test_empty_sequence(self):
        assert Sequence("s", "").gc_content() == 0.0


# ---------------------------------------------------------------------------
# Transcription
# ---------------------------------------------------------------------------

class TestTranscription:
    def test_dna_to_rna(self):
        rna = Sequence("s", "ATGC").transcribe()
        assert str(rna) == "AUGC"
        assert rna._seq_type == SequenceType.RNA

    def test_transcribe_rna_raises(self):
        with pytest.raises(TypeError):
            Sequence("s", "AUGC", seq_type=SequenceType.RNA).transcribe()

    def test_transcribe_preserves_header(self):
        assert Sequence("brca1", "ATGC").transcribe()._header == "brca1"


# ---------------------------------------------------------------------------
# Subseq
# ---------------------------------------------------------------------------

class TestSubseq:
    def test_basic_subseq(self):
        seq = Sequence("s", "ATGCATGC")
        assert str(seq.subseq(0, 4)) == "ATGC"

    def test_subseq_header_suffix(self):
        seq = Sequence("s", "ATGCATGC")
        assert "[2:6]" in seq.subseq(2, 6)._header


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------

class TestReadFasta:
    def test_single_sequence(self, tmp_path):
        f = tmp_path / "t.fasta"
        f.write_text(">seq1 desc\nATGCATGC\n")
        seqs = read_fasta(str(f))
        assert len(seqs) == 1
        assert str(seqs[0]) == "ATGCATGC"

    def test_multiline_sequence(self, tmp_path):
        f = tmp_path / "t.fasta"
        f.write_text(">seq1\nATGC\nATGC\n")
        seqs = read_fasta(str(f))
        assert str(seqs[0]) == "ATGCATGC"

    def test_multiple_sequences(self, tmp_path):
        f = tmp_path / "t.fasta"
        f.write_text(">s1\nATGC\n>s2\nGGCC\n")
        seqs = read_fasta(str(f))
        assert len(seqs) == 2
        assert str(seqs[1]) == "GGCC"

    def test_empty_file(self, tmp_path):
        f = tmp_path / "empty.fasta"
        f.write_text("")
        assert read_fasta(str(f)) == []

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            read_fasta("does_not_exist.fasta")

    def test_example_file(self):
        seqs = read_fasta("data/example.fasta")
        assert len(seqs) == 4


# ---------------------------------------------------------------------------
# FASTA write / roundtrip
# ---------------------------------------------------------------------------

class TestWriteFasta:
    def test_roundtrip(self, tmp_path):
        original = [Sequence("s1", "ATGCATGC"), Sequence("s2", "GGCCGGCC")]
        path = str(tmp_path / "out.fasta")
        write_fasta(original, path)
        recovered = read_fasta(path)
        assert [str(s) for s in recovered] == [str(s) for s in original]

    def test_line_width(self, tmp_path):
        seq = Sequence("s1", "A" * 100)
        path = str(tmp_path / "out.fasta")
        write_fasta([seq], path, line_width=60)
        lines = open(path).readlines()
        # Toutes les lignes de séquence font ≤ 61 chars (60 + \n)
        seq_lines = [l for l in lines if not l.startswith(">")]
        assert all(len(l) <= 61 for l in seq_lines)

    def test_empty_list(self, tmp_path):
        path = str(tmp_path / "empty.fasta")
        write_fasta([], path)
        assert read_fasta(path) == []
