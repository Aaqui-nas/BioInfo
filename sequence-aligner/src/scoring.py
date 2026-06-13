"""
TP03 — Matrices de substitution
Implémenter toutes les méthodes marquées raise NotImplementedError().
Ne pas modifier les signatures ni les docstrings.
"""
from __future__ import annotations


class ScoringMatrix:
    """
    Matrice de substitution pour l'alignement de séquences.

    Stocke les scores de substitution entre paires de caractères.
    Remplace les paramètres match_score/mismatch_score dans Needleman-Wunsch
    pour capturer le fait que toutes les substitutions ne sont pas équivalentes.
    """

    def __init__(self, data: dict[tuple[str, str], int]) -> None:
        """
        Initialise la matrice à partir d'un dictionnaire de scores.

        data : {(a, b): score} — toutes les paires (a, b) et (b, a) doivent être présentes.
        Exemple minimal pour un alphabet {A, T} :
            {('A','A'): 1, ('A','T'): -1, ('T','A'): -1, ('T','T'): 1}

        La matrice est supposée symétrique : score(a, b) == score(b, a).
        Stocker data tel quel — pas de validation requise.
        """
        self._data = data

    def score(self, a: str, b: str) -> int:
        """
        Retourne le score de substitution pour la paire (a, b).

        Lever KeyError si la paire n'est pas dans la matrice.
        """
        return self._data[(a, b)]

    @classmethod
    def from_file(cls, filepath: str) -> "ScoringMatrix":
        with open(filepath) as f:
            lines = [
                line.strip()
                for line in f
                if line.strip() and not line.lstrip().startswith("#")
            ]

        columns = lines[0].split()
        scores = {}

        for line in lines[1:]:
            row, *values = line.split()

            for col, score in zip(columns, map(int, values)):
                scores[(row, col)] = score

        return cls(scores)

    @classmethod
    def identity_dna(cls) -> ScoringMatrix:
        """
        Construit la matrice identité pour l'alphabet ADN {A, T, G, C} :
        - match    (a == b) : +1
        - mismatch (a != b) : -1

        Toutes les paires (a, b) pour a, b ∈ {A, T, G, C} doivent être présentes.
        """
        nucs = ("A","T","G","C")
        res = {}
        for nuc1 in nucs:
            for nuc2 in nucs:
                res[(nuc1,nuc2)] = 1 if nuc1 == nuc2 else -1
        return cls(res)

    @classmethod
    def nuc_simple(cls) -> ScoringMatrix:
        """
        Construit la matrice NUC_SIMPLE pour l'alphabet ADN {A, T, G, C} :
        - match        (a == b)                       : +1
        - transition   (A↔G ou C↔T)                  : -1
        - transversion (A↔C, A↔T, G↔C ou G↔T)        : -3

        Rappel biologique :
          Purines    : A, G  (cycle à deux anneaux)
          Pyrimidines: C, T  (cycle à un anneau)
          Transition   = substitution dans la même classe  → moins disruptive
          Transversion = substitution entre classes         → plus disruptive

        Toutes les paires (a, b) pour a, b ∈ {A, T, G, C} doivent être présentes.
        """
        nucs = ("A","T","G","C")
        transitions = frozenset((("A","G"),("G","A"),("T","C"),("C","T")))
        res = {}
        for nuc1 in nucs:
            for nuc2 in nucs:
                res[(nuc1,nuc2)] = 1 if nuc1 == nuc2 else -1 if (nuc1,nuc2) in transitions else -3
        return cls(res)
