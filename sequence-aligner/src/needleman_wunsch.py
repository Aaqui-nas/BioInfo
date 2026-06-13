"""
TP02 — Alignement global : Needleman-Wunsch
TP03 — Intégration d'une matrice de substitution (paramètre scoring_matrix)
Ne pas modifier les signatures ni les docstrings.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from sequences import Sequence

if TYPE_CHECKING:
    from scoring import ScoringMatrix


def fill_matrix(
    seq1: Sequence,
    seq2: Sequence,
    match_score: int = 1,
    mismatch_score: int = -1,
    gap_penalty: int = -2,
    scoring_matrix: ScoringMatrix | None = None,  # TP03
) -> list[list[int]]:
    """
    Construit et remplit la matrice de programmation dynamique de Needleman-Wunsch.

    La matrice a (len(seq1) + 1) lignes et (len(seq2) + 1) colonnes.
    seq1 est indexée sur les lignes (i), seq2 sur les colonnes (j).

    - F[0][0] = 0
    - F[i][0] = i * gap_penalty  (aligner i bases de seq1 contre une séquence vide)
    - F[0][j] = j * gap_penalty  (aligner j bases de seq2 contre une séquence vide)
    - F[i][j] = max(diagonal, haut, gauche)

    Retourne la matrice complète (list de lists d'entiers).

    TP03 : si scoring_matrix est fourni, utiliser scoring_matrix.score(a, b) à la place
    de match_score / mismatch_score pour calculer le score diagonal.
    """
    m = len(seq1)
    n = len(seq2)

    F = [[0]*(n+1) for _ in range (m+1)]
    for i in range(m+1):
        F[i][0] = i*gap_penalty
    for j in range(1,n+1):
        F[0][j] = j*gap_penalty

    for i in range(1,m+1):
        for j in range(1,n+1):
            if scoring_matrix is not None:
                diag = F[i-1][j-1] + scoring_matrix.score(seq1[i-1], seq2[j-1])
            else:
                if seq1[i-1] == seq2[j-1]:
                    diag = F[i-1][j-1] + match_score
                else:
                    diag = F[i-1][j-1] + mismatch_score
            top = F[i-1][j] + gap_penalty
            left = F[i][j-1] + gap_penalty
            F[i][j] = max(diag, top, left)
    return F

def traceback(
    matrix: list[list[int]],
    seq1: Sequence,
    seq2: Sequence,
    match_score: int = 1,
    mismatch_score: int = -1,
    gap_penalty: int = -2,
    scoring_matrix: ScoringMatrix | None = None,  # TP03
) -> tuple[str, str]:
    """
    Remonte la matrice depuis (m, n) jusqu'à (0, 0) pour reconstruire l'alignement.

    Retourne un tuple (aligned_seq1, aligned_seq2) de même longueur,
    où les gaps sont représentés par '-'.

    En cas d'égalité entre plusieurs directions : priorité diagonal > haut > gauche.

    Indice d'implémentation : construire deux listes en ajoutant à la fin (append),
    puis les inverser avec [::-1] avant de joindre — plus efficace que de préfixer des strings.

    TP03 : si scoring_matrix est fourni, utiliser scoring_matrix.score() pour calculer
    diag_score à la place de match_score / mismatch_score.
    """
    i, j = len(seq1), len(seq2)
    res1, res2 = [], []
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            # TP03 : remplacer ce calcul de diag_score par scoring_matrix.score() si fourni
            if scoring_matrix is not None:
                diag_score = scoring_matrix.score(seq1[i-1], seq2[j-1])
            else:
                diag_score = match_score if seq1[i-1]==seq2[j-1] else mismatch_score
            if matrix[i][j] == matrix[i-1][j-1] + diag_score:
                res1.append(seq1[i-1])
                res2.append(seq2[j-1])
                i -= 1
                j -= 1
            elif matrix[i][j] == matrix[i-1][j] + gap_penalty:
                res1.append(seq1[i-1])
                res2.append("-")
                i -= 1
            else:
                res1.append("-")
                res2.append(seq2[j-1])
                j -= 1
        elif i > 0:
            res1.append(seq1[i-1])
            res2.append("-")
            i -= 1
        else:
            res1.append("-")
            res2.append(seq2[j-1])
            j -= 1

    return ("".join(res1[::-1]), "".join(res2[::-1]))

def needleman_wunsch(
    seq1: Sequence,
    seq2: Sequence,
    match_score: int = 1,
    mismatch_score: int = -1,
    gap_penalty: int = -2,
    scoring_matrix: ScoringMatrix | None = None,  # TP03
) -> tuple[int, str, str]:
    """
    Algorithme de Needleman-Wunsch complet.

    Orchestre fill_matrix puis traceback.

    Retourne un tuple (score, aligned_seq1, aligned_seq2) :
    - score        : valeur en bas à droite de la matrice (score optimal global)
    - aligned_seq1 : seq1 alignée, avec '-' pour les gaps
    - aligned_seq2 : seq2 alignée, avec '-' pour les gaps
    Les deux chaînes alignées ont toujours la même longueur.

    TP03 : si scoring_matrix est fourni, le propager à fill_matrix et traceback.
    """
    F = fill_matrix(seq1, seq2, match_score, mismatch_score, gap_penalty, scoring_matrix)
    al_seq1, al_seq2 = traceback(F, seq1, seq2, match_score, mismatch_score, gap_penalty, scoring_matrix)
    return F[len(seq1)][len(seq2)], al_seq1, al_seq2
