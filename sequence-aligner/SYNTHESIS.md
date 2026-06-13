# Synthesis — sequence-aligner

**Date de clôture :** 2026-06-14
**Durée totale estimée :** 22–32 h
**Durée totale réelle :** ~3h50
**TPs complétés :** 6 / 6

---

## Ce qui a été construit

Un pipeline d'alignement pairwise de séquences ADN complet, en Python pur, sans dépendance algorithmique :

| Module | Contenu |
|--------|---------|
| `src/sequences.py` | Classe `Sequence`, parsing FASTA, générateur de séquences aléatoires |
| `src/scoring.py` | `ScoringMatrix` — BLOSUM62, NUC_SIMPLE, chargement depuis fichier `.mat` |
| `src/needleman_wunsch.py` | Alignement global linéaire (NW), matrice DP + traceback |
| `src/smith_waterman.py` | Alignement local (SW), détection de la cellule maximale |
| `src/affine_alignment.py` | Gotoh 1982 — pénalités affines (H/E/F), global + semi-global |
| `src/visualize.py` | Format BLAST-like, statistiques d'alignement |
| `src/benchmark.py` | Mesure de temps médian, comparaison avec BioPython |
| `src/main.py` | Benchmark BRCA1 humain/souris + validation croisée BioPython |

**Tests :** 116 tests, tous verts, répartis en 7 fichiers `tests/test_*.py`.

---

## Ce qui a bien marché

**`_NEG_INF` comme valeur sentinelle** — initialiser toutes les matrices à `float("-inf")` et ne corriger que les cellules légitimes a rendu les bugs d'initialisation immédiatement visibles (score sentinel reconnaissable) au lieu de silencieux (0 erroné).

**`_score` helper** — extraire le choix `scoring_matrix` vs `match/mismatch` en une fonction d'une ligne a évité de dupliquer la même branche `if` dans les matrices de remplissage et dans le traceback.

**`_traceback_affine` avec `stop_when_i_zero`** — factoriser le traceback global et semi-global en une seule fonction avec un flag booléen a évité ~60 lignes de duplication sur une logique d'état H/E/F identique. Ni trop peu (duplication), ni trop (sur-abstraction).

**`make_biopython_aligner` factory** — centraliser la configuration du `PairwiseAligner` BioPython (convention de gap, substitution matrix, end-gaps semi-global) dans une seule fonction a rendu les fonctions appelantes lisibles et la convention documentée à un seul endroit.

**Séparation bibliothèque / script** — `benchmark.py` exporte des fonctions pures testables, `main.py` orchestre. Cette séparation a permis d'écrire une suite de tests unitaires sur le benchmark sans dépendre des données BRCA1.

---

## Ce qui aurait pu être mieux conçu dès le départ

**Documenter la convention de gap dès TP02** — la différence entre `gap(k) = gap_open + k × gap_extend` (notre formule) et `gap(k) = open_bio + (k-1) × extend_bio` (BioPython) n'a été découverte et prouvée qu'au TP06. Si cette preuve avait été posée au TP02 (premier TP avec gaps), la validation croisée BioPython aurait été triviale.

**Tester les valeurs limites des arguments optionnels** — le bug `if seed:` (qui ignore `seed=0`) illustre une règle générale : toute valeur sentinelle ou optionnelle doit être testée explicitement à sa valeur limite. La suite de tests couvrait `seed=42` (correct) et `seed=0` avec le mauvais test (alphabet, pas reproductibilité). Un test `generate_random_dna(10, seed=0)` appelé deux fois aurait attrapé le bug immédiatement.

**Nommer les paramètres selon leur effet réel** — le paramètre `mode` de `make_biopython_aligner` ne contrôle pas `aligner.mode` mais uniquement les end-gaps. Un nom comme `semiglobal: bool` aurait évité l'ambiguïté.

---

## Compétences acquises

| Compétence | Niveau à l'entrée | Niveau à la sortie |
|---|---|---|
| Programmation dynamique (PD) | Notion | Maîtrisé — 3 variantes (NW, SW, Gotoh) |
| Alignement de séquences (ALI) | Zéro | Maîtrisé — global, local, semi-global, affin |
| Matrices de substitution (MST) | Zéro | Couvert — BLOSUM62, NUC_SIMPLE, lecture fichier |
| Formats bio FASTA (FMT) | Zéro | Partiel — FASTA uniquement |
| Biologie fondamentale (BIO_BASE) | Zéro | Partiel — ADN, conservation évolutive, BRCA1 |
| Statistiques sur séquences (STA) | Zéro | Partiel — identité, benchmarking empirique |
