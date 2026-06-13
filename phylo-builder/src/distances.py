from __future__ import annotations
import math


class DistanceMatrix:
    """Matrice de distances pairwise entre N séquences."""

    def __init__(self, names: list[str], matrix: list[list[float]]) -> None:
        """
        names  : liste ordonnée de noms de séquences (N noms)
        matrix : tableau N×N, matrix[i][j] est la distance entre names[i] et names[j]
                 Invariants : symétrique, diagonale nulle.
        """
        raise NotImplementedError()

    def __getitem__(self, pair: tuple[str, str]) -> float:
        """dm[("Humain", "Souris")] retourne la distance entre Humain et Souris."""
        raise NotImplementedError()

    @classmethod
    def from_sequences(
        cls,
        sequences: dict[str, str],
        method: str = "jukes_cantor",
    ) -> DistanceMatrix:
        """
        Construit la matrice depuis un dict {nom: séquence_alignée}.

        method : "hamming" | "jukes_cantor" | "kimura_2p"
        Toutes les séquences doivent avoir la même longueur.
        Lève ValueError si method est inconnu.
        """
        raise NotImplementedError()

    def __str__(self) -> str:
        """Affichage tabulaire : noms en en-tête de ligne et de colonne, distances à 4 décimales."""
        raise NotImplementedError()


def hamming_distance(seq1: str, seq2: str) -> float:
    """
    Proportion de positions différentes entre deux séquences alignées.

    Les positions où l'un ou l'autre caractère est '-' (gap) sont exclues du calcul.
    Lève ValueError si les séquences n'ont pas la même longueur.
    Retourne 0.0 si aucune position comparable n'existe (séquences de gaps purs).
    """
    raise NotImplementedError()


def jukes_cantor_distance(seq1: str, seq2: str) -> float:
    """
    Distance de Jukes-Cantor (1969) : corrige la distance de Hamming pour les
    substitutions multiples au même site.

    Hypothèse : taux de substitution identique pour toutes les paires de nucléotides.
    Formule : d = -3/4 * ln(1 - 4/3 * p)  où p = hamming_distance(seq1, seq2)

    Retourne math.inf si p >= 3/4 (saturation — le modèle ne peut plus corriger).
    """
    raise NotImplementedError()


def kimura_2p_distance(seq1: str, seq2: str) -> float:
    """
    Distance de Kimura 2 paramètres (1980) : distingue transitions et transversions.

    Transitions  (ts) : substitutions entre nucléotides de même classe
                        A <-> G (purines) ou C <-> T (pyrimidines)
    Transversions (tv) : substitutions entre classes différentes
                        (A ou G) <-> (C ou T)

    P = proportion de sites avec une transition
    Q = proportion de sites avec une transversion
    Formule : d = -1/2 * ln(1 - 2P - Q) - 1/4 * ln(1 - 2Q)

    Retourne math.inf si l'un des arguments du logarithme est <= 0.
    """
    raise NotImplementedError()
