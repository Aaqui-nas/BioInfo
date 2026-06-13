"""
TP06 — Benchmarking et comparaison avec BioPython
Implémenter toutes les fonctions marquées raise NotImplementedError().
Ne pas modifier les signatures ni les docstrings.
"""
from __future__ import annotations

import random
import statistics
import time
from typing import Any, Callable

from sequences import Sequence
from affine_alignment import global_affine

from Bio.Align import PairwiseAligner

def generate_random_dna(length: int, seed: int | None = None) -> Sequence:
    """
    Génère une séquence ADN aléatoire de longueur `length`.

    Si seed est fourni, initialiser random.seed(seed) en début de fonction
    pour garantir la reproductibilité de manière isolée.

    Le header de la Sequence retournée est "random".

    Utiliser random.choice("ACGT") pour chaque base.
    """
    if seed is not None:
        random.seed(seed)
    return Sequence("random", "".join(random.choice("ACGT") for _ in range(length)))


def time_function(
    func: Callable,
    args: tuple = (),
    kwargs: dict[str, Any] | None = None,
    n_repeat: int = 5,
) -> float:
    """
    Mesure le temps d'exécution médian de func(*args, **kwargs) sur n_repeat appels.

    Utiliser time.perf_counter() pour la mesure de temps.
    Utiliser statistics.median() sur les n_repeat durées.

    Retourne le temps médian en secondes (float).
    """
    if kwargs is None:
        kwargs = {}

    times = []

    for _ in range(n_repeat):
        start = time.perf_counter()
        func(*args, **kwargs)
        end = time.perf_counter()
        times.append(end - start)

    return statistics.median(times)


def compare_with_biopython(
    seq1: Sequence,
    seq2: Sequence,
    gap_open: int = -5,
    gap_extend: int = -1,
    match_score: int = 1,
    mismatch_score: int = -1,
) -> dict:
    """
    Compare le score de global_affine avec celui de Bio.Align.PairwiseAligner (mode global).

    ⚠ DIFFÉRENCE DE CONVENTION entre notre formule et BioPython :
      Notre formule : gap(k) = gap_open + k × gap_extend
      BioPython     : gap(k) = open_bio + (k-1) × extend_bio

    Pour rendre les deux équivalentes, il faut convertir :
      aligner.open_gap_score   = gap_open + gap_extend   ← correction obligatoire
      aligner.extend_gap_score = gap_extend

    Vérification :
      BioPython gap(k) = (gap_open + gap_extend) + (k-1) × gap_extend
                       = gap_open + k × gap_extend  ✓

    Retourne un dict avec les clés :
      "our_score"        : float — score de global_affine
      "biopython_score"  : float — score de Bio.Align.PairwiseAligner
      "match"            : bool  — True si |our_score - biopython_score| < 1e-6
    """
    our_score, _, _ = global_affine(
        seq1,
        seq2,
        gap_open=gap_open,
        gap_extend=gap_extend,
        match_score=match_score,
        mismatch_score=mismatch_score,
    )

    aligner = PairwiseAligner()
    aligner.mode = "global"

    aligner.open_gap_score = gap_open + gap_extend
    aligner.extend_gap_score = gap_extend

    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score

    biopython_score = aligner.score(seq1, seq2)

    return {
        "our_score": our_score,
        "biopython_score": biopython_score,
        "match": abs(our_score - biopython_score) < 1e-6,
    }
