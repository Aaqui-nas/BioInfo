"""
TP05 — Gap penalties affines + alignement semi-global (Gotoh 1982)
Implémenter toutes les fonctions marquées raise NotImplementedError().
Ne pas modifier les signatures ni les docstrings.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from sequences import Sequence

if TYPE_CHECKING:
    from scoring import ScoringMatrix

_NEG_INF = float("-inf")

def _score(a, b, match_score, mismatch_score, scoring_matrix):
    if scoring_matrix is not None:
        return scoring_matrix.score(a, b)
    return match_score if a == b else mismatch_score


def _traceback_affine(
    H,
    E,
    F,
    seq1,
    seq2,
    start_i,
    start_j,
    stop_when_i_zero,
    gap_open,
    gap_extend,
    match_score,
    mismatch_score,
    scoring_matrix,
):
    i, j = start_i, start_j
    state = "H"

    a1 = []
    a2 = []

    while i > 0 if stop_when_i_zero else (i > 0 or j >0):

        if state == "H":

            if i > 0 and j > 0:
                diag = (
                    H[i - 1][j - 1]
                    + _score(
                        seq1[i - 1],
                        seq2[j - 1],
                        match_score,
                        mismatch_score,
                        scoring_matrix,
                    )
                )

                if H[i][j] == diag:
                    a1.append(seq1[i - 1])
                    a2.append(seq2[j - 1])
                    i -= 1
                    j -= 1
                    continue

            if H[i][j] == E[i][j]:
                state = "E"
                continue

            if H[i][j] == F[i][j]:
                state = "F"
                continue

            raise RuntimeError("Traceback H impossible")

        elif state == "E":

            a1.append(seq1[i - 1])
            a2.append("-")

            if E[i][j] == E[i - 1][j] + gap_extend:
                i -= 1
            else:
                i -= 1
                state = "H"

        else:  # state == "F"

            a1.append("-")
            a2.append(seq2[j - 1])

            if F[i][j] == F[i][j - 1] + gap_extend:
                j -= 1
            else:
                j -= 1
                state = "H"

    return "".join(reversed(a1)), "".join(reversed(a2))

def fill_matrices_global(
    seq1: Sequence,
    seq2: Sequence,
    gap_open: int = -5,
    gap_extend: int = -1,
    match_score: int = 1,
    mismatch_score: int = -1,
    scoring_matrix: "ScoringMatrix | None" = None,
) -> tuple[list[list[float]], list[list[float]], list[list[float]]]:
    """
    Construit H, E, F pour l'alignement global affin (Gotoh 1982).

    Trois matrices (toutes de taille (m+1) × (n+1)) :
    - H[i][j] : meilleur score pour aligner seq1[:i] et seq2[:j]
    - E[i][j] : meilleur score se terminant par un gap dans seq2
                (seq1[i-1] face à un gap — mouvement vertical ↓)
    - F[i][j] : meilleur score se terminant par un gap dans seq1
                (seq2[j-1] face à un gap — mouvement horizontal →)

    Initialisation (global) :
      H[0][0] = 0  ;  E[0][0] = F[0][0] = -inf
      H[i][0] = E[i][0] = gap_open + i * gap_extend  ;  F[i][0] = -inf
      H[0][j] = F[0][j] = gap_open + j * gap_extend  ;  E[0][j] = -inf

    Récurrence (i ≥ 1, j ≥ 1) :
      E[i][j] = max(H[i-1][j] + gap_open + gap_extend,  E[i-1][j] + gap_extend)
      F[i][j] = max(H[i][j-1] + gap_open + gap_extend,  F[i][j-1] + gap_extend)
      H[i][j] = max(H[i-1][j-1] + s(i,j),  E[i][j],  F[i][j])

    Si scoring_matrix fournie : s(i,j) = scoring_matrix.score(seq1[i-1], seq2[j-1]).
    Sinon : s(i,j) = match_score ou mismatch_score selon seq1[i-1] == seq2[j-1].

    Retourne (H, E, F).
    """
    m, n = len(seq1), len(seq2)

    H = [[_NEG_INF] * (n + 1) for _ in range(m + 1)]
    E = [[_NEG_INF] * (n + 1) for _ in range(m + 1)]
    F = [[_NEG_INF] * (n + 1) for _ in range(m + 1)]

    H[0][0] = 0

    for i in range(1, m + 1):
        H[i][0] = gap_open + i * gap_extend
        E[i][0] = gap_open + i * gap_extend

    for j in range(1, n + 1):
        H[0][j] = gap_open + j * gap_extend
        F[0][j] = gap_open + j * gap_extend

    for i in range(1, m + 1):
        for j in range(1, n + 1):

            E[i][j] = max(
                H[i - 1][j] + gap_open + gap_extend,
                E[i - 1][j] + gap_extend,
            )

            F[i][j] = max(
                H[i][j - 1] + gap_open + gap_extend,
                F[i][j - 1] + gap_extend,
            )

            diag = (
                H[i - 1][j - 1]
                + _score(
                    seq1[i - 1],
                    seq2[j - 1],
                    match_score,
                    mismatch_score,
                    scoring_matrix,
                )
            )

            H[i][j] = max(diag, E[i][j], F[i][j])

    return H, E, F


def fill_matrices_semiglobal(
    seq1: Sequence,
    seq2: Sequence,
    gap_open: int = -5,
    gap_extend: int = -1,
    match_score: int = 1,
    mismatch_score: int = -1,
    scoring_matrix: "ScoringMatrix | None" = None,
) -> tuple[list[list[float]], list[list[float]], list[list[float]]]:
    """
    Identique à fill_matrices_global sauf l'initialisation de la première ligne.

    Alignement semi-global : seq1 (read court) peut commencer et finir
    n'importe où dans seq2 (génome) sans pénaliser les extrémités de seq2.

    Différences avec fill_matrices_global (deux lignes seulement) :
      H[0][j] = 0   pour j = 0..n   (au lieu de gap_open + j * gap_extend)
      F[0][j] = 0   pour j = 0..n   (au lieu de gap_open + j * gap_extend)

    La première colonne et la récurrence intérieure sont identiques.

    Retourne (H, E, F).
    """
    m, n = len(seq1), len(seq2)

    H = [[_NEG_INF] * (n + 1) for _ in range(m + 1)]
    E = [[_NEG_INF] * (n + 1) for _ in range(m + 1)]
    F = [[_NEG_INF] * (n + 1) for _ in range(m + 1)]

    for j in range(n + 1):
        H[0][j] = 0
        F[0][j] = 0

    for i in range(1, m + 1):
        H[i][0] = gap_open + i * gap_extend
        E[i][0] = gap_open + i * gap_extend

    for i in range(1, m + 1):
        for j in range(1, n + 1):

            E[i][j] = max(
                H[i - 1][j] + gap_open + gap_extend,
                E[i - 1][j] + gap_extend,
            )

            F[i][j] = max(
                H[i][j - 1] + gap_open + gap_extend,
                F[i][j - 1] + gap_extend,
            )

            diag = (
                H[i - 1][j - 1]
                + _score(
                    seq1[i - 1],
                    seq2[j - 1],
                    match_score,
                    mismatch_score,
                    scoring_matrix,
                )
            )

            H[i][j] = max(diag, E[i][j], F[i][j])

    return H, E, F


def find_semiglobal_end(H: list[list[float]]) -> tuple[int, float]:
    """
    Parcourt la dernière ligne de H et retourne (end_col, score).

    end_col : indice j où H[m][j] est maximal  (m = len(H) - 1)
    score   : valeur de ce maximum

    Le traceback semi-global démarre à (m, end_col) au lieu de (m, n).
    En cas d'égalité, retourner le premier j rencontré.
    """
    m = len(H) - 1

    best_j = 0
    best_score = H[m][0]

    for j in range(1, len(H[0])):
        if H[m][j] > best_score:
            best_score = H[m][j]
            best_j = j

    return best_j, best_score


def traceback_global(
    H: list[list[float]],
    E: list[list[float]],
    F: list[list[float]],
    seq1: Sequence,
    seq2: Sequence,
    gap_open: int = -5,
    gap_extend: int = -1,
    match_score: int = 1,
    mismatch_score: int = -1,
    scoring_matrix: "ScoringMatrix | None" = None,
) -> tuple[str, str]:
    """
    Remonte les matrices depuis (m, n) pour reconstruire l'alignement global affin.

    Nécessite un suivi d'état : à chaque position, on est dans l'état H, E ou F.
    - Débuter dans l'état H (le score optimal est toujours en H[m][n]).
    - État H : vérifier diagonal (→ H), puis E (→ bascule vers E), puis F (→ bascule vers F).
    - État E : consommer seq1[i-1] face à '-', vérifier si nouveau gap (→ H) ou extension (→ E), i← i-1.
    - État F : consommer '-' face à seq2[j-1], vérifier si nouveau gap (→ H) ou extension (→ F), j← j-1.
    - Condition d'arrêt : i == 0 ET j == 0.

    Voir pseudo-code complet dans le README.

    Retourne (aligned_seq1, aligned_seq2).
    """
    return _traceback_affine(
        H,
        E,
        F,
        seq1,
        seq2,
        len(seq1),
        len(seq2),
        False,
        gap_open,
        gap_extend,
        match_score,
        mismatch_score,
        scoring_matrix,
    )


def traceback_semiglobal(
    H: list[list[float]],
    E: list[list[float]],
    F: list[list[float]],
    seq1: Sequence,
    seq2: Sequence,
    end_col: int,
    gap_open: int = -5,
    gap_extend: int = -1,
    match_score: int = 1,
    mismatch_score: int = -1,
    scoring_matrix: "ScoringMatrix | None" = None,
) -> tuple[str, str]:
    """
    Remonte les matrices depuis (m, end_col) jusqu'à i == 0.

    Deux différences vs traceback_global :
    - Point de départ : (m, end_col) au lieu de (m, n).
    - Condition d'arrêt : i == 0  (seq1 entièrement consommée).
      → on ne va pas nécessairement jusqu'à j == 0.

    La logique de suivi d'état (H, E, F) est identique à traceback_global.

    Retourne (aligned_seq1, aligned_seq2) — uniquement la région alignée,
    sans les flancs de seq2 avant/après l'alignement.
    """
    return _traceback_affine(
        H,
        E,
        F,
        seq1,
        seq2,
        len(seq1),
        end_col,
        True,
        gap_open,
        gap_extend,
        match_score,
        mismatch_score,
        scoring_matrix,
    )


def global_affine(
    seq1: Sequence,
    seq2: Sequence,
    gap_open: int = -5,
    gap_extend: int = -1,
    match_score: int = 1,
    mismatch_score: int = -1,
    scoring_matrix: "ScoringMatrix | None" = None,
) -> tuple[float, str, str]:
    """
    Alignement global avec pénalités de gap affines (Gotoh 1982).

    Orchestre : fill_matrices_global → traceback_global.

    Retourne (score, aligned_seq1, aligned_seq2) :
    - score : H[m][n], le score optimal de l'alignement global
    - Les deux chaînes alignées ont toujours la même longueur.
    """
    H, E, F = fill_matrices_global(
        seq1,
        seq2,
        gap_open,
        gap_extend,
        match_score,
        mismatch_score,
        scoring_matrix,
    )

    a1, a2 = traceback_global(
        H,
        E,
        F,
        seq1,
        seq2,
        gap_open,
        gap_extend,
        match_score,
        mismatch_score,
        scoring_matrix,
    )

    return H[len(seq1)][len(seq2)], a1, a2


def semiglobal(
    seq1: Sequence,
    seq2: Sequence,
    gap_open: int = -5,
    gap_extend: int = -1,
    match_score: int = 1,
    mismatch_score: int = -1,
    scoring_matrix: "ScoringMatrix | None" = None,
) -> tuple[float, str, str]:
    """
    Alignement semi-global : seq1 (read) flotte dans seq2 (génome).

    Orchestre : fill_matrices_semiglobal → find_semiglobal_end → traceback_semiglobal.

    Pas de pénalité pour les flancs de seq2 avant et après l'alignement.
    Pénalités normales pour les gaps à l'intérieur de seq1.

    Retourne (score, aligned_seq1, aligned_seq2) :
    - score     : max sur H[m][*] — score de la meilleure position dans seq2
    - aligned_* : uniquement la région alignée (sans les flancs de seq2)
    - Les deux chaînes alignées ont toujours la même longueur.
    """
    H, E, F = fill_matrices_semiglobal(
        seq1,
        seq2,
        gap_open,
        gap_extend,
        match_score,
        mismatch_score,
        scoring_matrix,
    )

    end_col, score = find_semiglobal_end(H)

    a1, a2 = traceback_semiglobal(
        H,
        E,
        F,
        seq1,
        seq2,
        end_col,
        gap_open,
        gap_extend,
        match_score,
        mismatch_score,
        scoring_matrix,
    )

    return score, a1, a2
