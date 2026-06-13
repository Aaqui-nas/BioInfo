"""
TP06 — Visualisation d'alignements
Implémenter toutes les fonctions marquées raise NotImplementedError().
Ne pas modifier les signatures ni les docstrings.
"""
from __future__ import annotations


def format_alignment(
    aligned1: str,
    aligned2: str,
    id1: str = "seq1",
    id2: str = "seq2",
    line_width: int = 60,
) -> str:
    """
    Affiche deux chaînes alignées en format texte (BLAST-like).

    Chaque bloc contient trois lignes :
    - Ligne seq1  : f"{id1}: {aligned1[start:end]}"
    - Ligne milieu : " " * (len(id1) + 2) + symboles de correspondance
    - Ligne seq2  : f"{id2}: {aligned2[start:end]}"

    Symboles de correspondance (position par position) :
    - "|" si les deux bases sont identiques et aucune n'est '-'
    - " " sinon (gap ou mismatch)

    Les blocs sont séparés par une ligne vide ("\n\n").

    Paramètres :
      aligned1, aligned2 : chaînes de même longueur (avec '-' pour les gaps)
      id1, id2           : identifiants des séquences
      line_width         : nombre de positions d'alignement par bloc

    Retourne la chaîne formatée.
    """
    if len(aligned1) != len(aligned2):
        raise ValueError("Les alignements doivent avoir la même longueur")

    blocks = []

    for start in range(0, len(aligned1), line_width):
        end = start + line_width

        s1 = aligned1[start:end]
        s2 = aligned2[start:end]

        match_line = "".join(
            "|" if a == b and a != "-" else " "
            for a, b in zip(s1, s2)
        )

        blocks.append(
            f"{id1}: {s1}\n"
            f"{' ' * (len(id1) + 2)}{match_line}\n"
            f"{id2}: {s2}"
        )

    return "\n\n".join(blocks)


def alignment_stats(aligned1: str, aligned2: str) -> dict:
    """
    Calcule les statistiques d'un alignement.

    Parcourt les deux chaînes position par position :
    - gap      : au moins un des deux caractères est '-'
    - match    : les deux bases sont identiques (et aucune n'est '-')
    - mismatch : les deux bases sont différentes et aucune n'est '-'

    Retourne un dict avec les clés :
      "length"     : int   — longueur de l'alignement (len(aligned1))
      "matches"    : int   — nombre de positions match
      "mismatches" : int   — nombre de positions mismatch
      "gaps"       : int   — nombre de positions contenant au moins un '-'
      "identity"   : float — matches / length  (0.0 si length == 0)

    Précondition : len(aligned1) == len(aligned2)
    """
    if len(aligned1) != len(aligned2):
        raise ValueError("Les alignements doivent avoir la même longueur")

    gaps = sum(a == "-" or b == "-" for a, b in zip(aligned1, aligned2))
    matches = sum(a == b and a != "-" for a, b in zip(aligned1, aligned2))
    mismatches = len(aligned1) - matches - gaps

    return {
        "length": len(aligned1),
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "identity": matches / len(aligned1) if aligned1 else 0.0,
    }
