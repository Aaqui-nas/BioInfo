"""
TP04 — Alignement local : Smith-Waterman
Implémenter toutes les fonctions marquées raise NotImplementedError().
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
    gap_penalty: int = -2,
    match_score: int = 1,
    mismatch_score: int = -1,
    scoring_matrix: ScoringMatrix | None = None,
) -> list[list[int]]:
    """
    Construit et remplit la matrice de programmation dynamique de Smith-Waterman.

    Différences fondamentales vs Needleman-Wunsch :
    1. F[i][0] = 0 pour tout i  (pas de i × gap_penalty)
    2. F[0][j] = 0 pour tout j  (pas de j × gap_penalty)
    3. F[i][j] = max(0, diag, haut, gauche)  ← le 0 est la clé de l'alignement local

    Le plancher à 0 signifie : si tous les chemins mènent à un score négatif,
    on repart de zéro — on oublie ce qui précède plutôt que d'accumuler des pénalités.

    Si scoring_matrix est fourni, utiliser scoring_matrix.score(a, b) pour le diagonal.
    Sinon, utiliser match_score / mismatch_score comme en NW.

    Retourne la matrice complète (list de lists d'entiers).
    """
    m = len(seq1)
    n = len(seq2)

    F = [[0]*(n+1) for _ in range (m+1)]

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
            F[i][j] = max(0, diag, top, left)
    return F


def find_max(matrix: list[list[int]]) -> tuple[int, int, int]:
    """
    Parcourt toute la matrice et retourne (row, col, score) de la cellule maximale.

    En cas d'égalité, retourner la première rencontrée en parcourant ligne par ligne,
    de gauche à droite (ordre naturel des boucles imbriquées).

    Le score retourné est la valeur optimale de l'alignement local.
    """
    max_score, max_row, max_col = 0,0,0
    for i in range (len(matrix)):
        for j in range(len(matrix[i])):
            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_row, max_col = i, j
    return max_row, max_col, max_score


def traceback(
    matrix: list[list[int]],
    seq1: Sequence,
    seq2: Sequence,
    start_row: int,
    start_col: int,
    gap_penalty: int = -2,
    match_score: int = 1,
    mismatch_score: int = -1,
    scoring_matrix: ScoringMatrix | None = None,
) -> tuple[str, str]:
    """
    Remonte la matrice depuis (start_row, start_col) pour reconstruire l'alignement local.

    Différences fondamentales vs traceback NW :
    - Démarre à (start_row, start_col) et NON à (m, n).
    - Condition d'arrêt : matrix[i][j] == 0  (et NON i==0 ET j==0).
      Dès qu'on atteint une cellule à 0, l'alignement local est terminé.

    La logique de déplacement (diagonal / haut / gauche) est identique à NW.

    Retourne ("", "") si start_row ou start_col vaut 0, ou si matrix[start_row][start_col] == 0.
    """
    i, j = start_row, start_col
    res1, res2 = [], []
    if (
        start_row == 0
        or start_col == 0
        or matrix[start_row][start_col] == 0
    ):
        return ("", "")
    while matrix[i][j] != 0:
        if i > 0 and j > 0:
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


def smith_waterman(
    seq1: Sequence,
    seq2: Sequence,
    gap_penalty: int = -2,
    match_score: int = 1,
    mismatch_score: int = -1,
    scoring_matrix: ScoringMatrix | None = None,
) -> tuple[int, str, str]:
    """
    Algorithme de Smith-Waterman complet.

    Orchestre fill_matrix → find_max → traceback.

    Retourne un tuple (score, aligned_seq1, aligned_seq2) :
    - score        : score optimal de l'alignement local (≥ 0)
    - aligned_seq1 : sous-séquence de seq1 alignée, avec '-' pour les gaps
    - aligned_seq2 : sous-séquence de seq2 alignée, avec '-' pour les gaps
    Les deux chaînes alignées ont toujours la même longueur.

    Si le score optimal est 0 : retourner (0, "", "").
    """
    F = fill_matrix(seq1, seq2, gap_penalty, match_score, mismatch_score, scoring_matrix)
    max_row, max_col, max_score = find_max(F)
    al_seq1, al_seq2 = traceback(F, seq1, seq2, max_row, max_col,gap_penalty, match_score, mismatch_score, scoring_matrix)
    return max_score, al_seq1, al_seq2
