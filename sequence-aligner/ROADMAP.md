# Roadmap — sequence-aligner

**Durée estimée totale :** 22–32 heures sur ~6 semaines
**Technologies :** Python

---

## Vue d'ensemble des TPs

| TP | Titre | Compétences | Estimé | Prérequis |Statut |
|----|-------|-------------|--------|-----------|-------|
| 01 | Représentation des séquences + parsing FASTA | BIO_BASE, FMT | 3–4 h | — | [ ] |
| 02 | Alignement global — Needleman-Wunsch | PD, ALI | 5–7 h | TP01 | [ ] |
| 03 | Matrices de substitution — BLOSUM62 | MST, STA | 3–4 h | TP02 | [ ] |
| 04 | Alignement local — Smith-Waterman | PD, ALI | 4–6 h | TP03 | [ ] |
| 05 | Gap penalties affines + alignement semi-global | ALI | 4–6 h | TP04 | [ ] |
| 06 | Benchmarking + visualisation | STA, HEU | 3–5 h | TP05 | [ ] |

`[ ]` non démarré · `[~]` en cours · `[x]` validé

---

## Progression fonctionnelle

**TP01 →** On peut lire des fichiers FASTA, représenter une séquence ADN/ARN comme objet, calculer son complément inverse et son contenu en GC.

**TP02 →** TP01 + on peut aligner deux séquences globalement avec un scoring simple (+1 match, -1 mismatch, -2 gap) et récupérer l'alignement optimal par traceback.

**TP03 →** TP02 + le scoring utilise une vraie matrice BLOSUM62 chargée depuis un fichier ; on peut comparer l'effet de la matrice sur la qualité des alignements.

**TP04 →** TP02+03 + on peut trouver les régions localement similaires entre deux séquences (Smith-Waterman) — utile pour des séquences partiellement homologues.

**TP05 →** TP04 + le modèle de gap est affine (coût d'ouverture ≠ coût d'extension) et on gère l'alignement semi-global (read court aligné sur un long génome sans pénaliser les bords).

**TP06 →** TP01–05 + on compare nos implémentations à celles de BioPython sur des séquences réelles (BRCA1 humain/souris), on mesure les performances et on visualise les alignements en mode texte.

---

## Évolutions possibles hors roadmap

- Alignement multiple (ClustalW-like) — objet de la phylogénétique
- Index k-mer pour filtrage pré-alignement (heuristique BLAST)
- Parallélisation sur matrices larges (numpy vectorisation)
