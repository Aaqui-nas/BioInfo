import math
import pytest
from distances import DistanceMatrix, hamming_distance, jukes_cantor_distance, kimura_2p_distance


class TestHamming:
    def test_identical(self):
        assert hamming_distance("ACGT", "ACGT") == 0.0

    def test_all_different(self):
        assert hamming_distance("AAAA", "TTTT") == 1.0

    def test_half_different(self):
        assert hamming_distance("AACC", "TTCC") == 0.5

    def test_one_difference(self):
        assert hamming_distance("AATT", "ACTT") == pytest.approx(0.25)

    def test_gaps_excluded(self):
        # Position avec '-' exclue : 3 positions comparables, 1 différente
        assert hamming_distance("A-GT", "T-GT") == pytest.approx(1 / 3)

    def test_unequal_length_raises(self):
        with pytest.raises(ValueError):
            hamming_distance("ACGT", "ACG")


class TestJukesCantor:
    def test_identical(self):
        assert jukes_cantor_distance("ACGT", "ACGT") == 0.0

    def test_greater_than_hamming(self):
        jc = jukes_cantor_distance("AATT", "ACTT")
        ham = hamming_distance("AATT", "ACTT")
        assert jc > ham

    def test_saturation_returns_inf(self):
        # p = 1.0 >= 3/4 → saturation
        assert math.isinf(jukes_cantor_distance("AAAA", "TTTT"))

    def test_known_value(self):
        # p = 0.25 → d = -3/4 * ln(1 - 4/3 * 0.25) = -3/4 * ln(2/3)
        p = 0.25
        expected = -0.75 * math.log(1 - 4 / 3 * p)
        assert jukes_cantor_distance("AATT", "ACTT") == pytest.approx(expected)


class TestKimura2P:
    def test_identical(self):
        assert kimura_2p_distance("ACGT", "ACGT") == 0.0

    def test_pure_transitions_saturated(self):
        # "AATT" vs "GGTT" : A→G (ts), A→G (ts) — P=0.5, Q=0 → 1-2P=0 → inf
        assert math.isinf(kimura_2p_distance("AATT", "GGTT"))

    def test_pure_transversions_saturated(self):
        # "AATT" vs "CCTT" : A→C (tv), A→C (tv) — P=0, Q=0.5 → 1-2Q=0 → inf
        assert math.isinf(kimura_2p_distance("AATT", "CCTT"))

    def test_transition_gives_larger_distance_than_transversion(self):
        # Pour la même proportion de sites différents (1/5),
        # une transition implique une distance K2P > transversion
        # "AATTT" vs "AGTTT" : A→G = transition (P=0.2, Q=0)
        d_ts = kimura_2p_distance("AATTT", "AGTTT")
        # "AATTT" vs "ACTTT" : A→C = transversion (P=0, Q=0.2)
        d_tv = kimura_2p_distance("AATTT", "ACTTT")
        assert not math.isinf(d_ts)
        assert not math.isinf(d_tv)
        assert d_ts > d_tv

    def test_known_value_transition(self):
        # "AATT" vs "AGTT" : A→G (ts) — P=0.25, Q=0
        # d = -1/2 * ln(1 - 2*0.25 - 0) - 1/4 * ln(1 - 0)
        P, Q = 0.25, 0.0
        expected = -0.5 * math.log(1 - 2 * P - Q) - 0.25 * math.log(1 - 2 * Q)
        assert kimura_2p_distance("AATT", "AGTT") == pytest.approx(expected)

    def test_known_value_transversion(self):
        # "AATT" vs "ACTT" : A→C (tv) — P=0, Q=0.25
        # d = -1/2 * ln(1 - 0 - 0.25) - 1/4 * ln(1 - 0.5)
        P, Q = 0.0, 0.25
        expected = -0.5 * math.log(1 - 2 * P - Q) - 0.25 * math.log(1 - 2 * Q)
        assert kimura_2p_distance("AATT", "ACTT") == pytest.approx(expected)


class TestDistanceMatrix:
    def setup_method(self):
        self.seqs = {
            "A": "AATT",
            "B": "ACTT",
            "C": "AATT",
        }

    def test_symmetric(self):
        dm = DistanceMatrix.from_sequences(self.seqs)
        assert dm[("A", "B")] == pytest.approx(dm[("B", "A")])
        assert dm[("A", "C")] == pytest.approx(dm[("C", "A")])

    def test_zero_diagonal(self):
        dm = DistanceMatrix.from_sequences(self.seqs)
        for name in self.seqs:
            assert dm[(name, name)] == pytest.approx(0.0)

    def test_identical_sequences_distance_zero(self):
        dm = DistanceMatrix.from_sequences(self.seqs)
        assert dm[("A", "C")] == pytest.approx(0.0)

    def test_default_method_is_jukes_cantor(self):
        dm = DistanceMatrix.from_sequences(self.seqs)
        expected = jukes_cantor_distance("AATT", "ACTT")
        assert dm[("A", "B")] == pytest.approx(expected)

    def test_hamming_method(self):
        dm = DistanceMatrix.from_sequences(self.seqs, method="hamming")
        assert dm[("A", "B")] == pytest.approx(0.25)

    def test_kimura_method(self):
        dm = DistanceMatrix.from_sequences(self.seqs, method="kimura_2p")
        expected = kimura_2p_distance("AATT", "ACTT")
        assert dm[("A", "B")] == pytest.approx(expected)

    def test_unknown_method_raises(self):
        with pytest.raises(ValueError):
            DistanceMatrix.from_sequences(self.seqs, method="unknown")

    def test_str_contains_names_and_values(self):
        dm = DistanceMatrix.from_sequences(self.seqs, method="hamming")
        s = str(dm)
        assert "A" in s and "B" in s and "C" in s
        assert "0.25" in s or "0.2500" in s
