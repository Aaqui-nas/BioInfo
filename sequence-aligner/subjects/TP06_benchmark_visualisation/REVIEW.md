# Review — TP06 : Benchmarking + visualisation

**Projet :** sequence-aligner
**Date :** 2026-06-14
**Temps estimé / réel :** 3–5 h / *à remplir*

---

## 🔴 Bloquant — corriger avant le TP suivant

Aucun *(dernier TP du projet)*.

---

## 🟠 Dette — rembourser avant la clôture

**`if seed:` dans `generate_random_dna` (benchmark.py, ligne 29)** — correctness bug silencieux. `if seed:` évalue `0` comme `False`, donc `generate_random_dna(100, seed=0)` n'initialise pas le générateur et retourne une séquence non reproductible sans avertissement.

Le correctif est une seule ligne : `if seed is not None:`.

Ce bug n'est pas attrapé par la suite de tests parce que `test_valid_alphabet` utilise `seed=0` mais ne vérifie que l'alphabet (toujours vrai), pas la reproductibilité. `test_reproducible_with_seed` utilise `seed=42`, qui fonctionne correctement. Leçon : quand une fonction a un comportement spécial pour une valeur limite (ici 0), il faut un test de reproductibilité explicitement avec cette valeur limite.

---

## 🟡 Cosmétique — non urgent

- [ ] **`aligner.mode = "global"` hardcodé dans `make_biopython_aligner`** (main.py, ligne 31) : le paramètre `mode` de la fonction contrôle uniquement les pénalités de fin de gap, pas la propriété `aligner.mode` de BioPython. Le nom de paramètre est trompeur — un lecteur qui parcourt la signature s'attend à ce que `mode` soit propagé à `aligner.mode`. Renommer en `semiglobal: bool = False` ou remplacer par un commentaire explicitant que le mode BioPython est toujours `"global"` (le semi-global est un trick via end-gaps).

- [ ] **`str()` absent dans `biopython_score` et les lambdas de timing** (main.py, lignes 62, 76, 90) : `aligner.score(seq1, seq2)` est appelé avec des objets `Sequence`, pas des `str`. Ça fonctionne actuellement parce que BioPython itère implicitement les objets en mode match/mismatch et que notre `Sequence` est itérable sur des chars ACGT. Mais dès que `aligner.substitution_matrix` est défini (cas `scoring_matrix is not None`), BioPython vérifie l'alphabet strict et peut lever une `ValueError`. Le fix trivial — `aligner.score(str(seq1), str(seq2))` — est plus explicite et immunise contre ce cas.

---

## Ce qui est bien fait

**`make_biopython_aligner` factory** — excellente décision d'architecture. Sans ce helper, la configuration du `PairwiseAligner` (gap convention, mode, substitution_matrix) aurait été dupliquée dans `biopython_score`, `time_biopython_global` et `time_biopython_semiglobal`. Une factory qui centralise la configuration + des fonctions appelantes d'une ligne : c'est exactement la bonne granularité.

**Séparation `benchmark.py` (bibliothèque) / `main.py` (script)** — `benchmark.py` exporte des fonctions pures et testables (`generate_random_dna`, `time_function`, `compare_with_biopython`) sans aucune logique de présentation ni d'entrée-sortie. `main.py` orchestre. Cette séparation permet de tester les fonctions unitairement sans exécuter le benchmark complet — et c'est ce qui a rendu la suite de tests directe à écrire.

**Convention BioPython documentée avec la preuve algébrique** — ne pas juste poser `aligner.open_gap_score = gap_open + gap_extend` sans explication, mais documenter la dérivation dans le docstring : `gap(k) = (gap_open + gap_extend) + (k−1) × gap_extend = gap_open + k × gap_extend ✓`. Quand quelqu'un relira ce code dans six mois, cette preuve en trois lignes évitera une heure de débogage.

**Médiane plutôt que moyenne pour le timing** — robuste aux pics de charge OS. Une moyenne sur 5 runs est polluée par un GC ou un context switch ; la médiane les absorbe. Bon réflexe de métrologie.

**`if kwargs is None: kwargs = {}` pattern** — éviter `kwargs={}` comme valeur par défaut d'argument. En Python, les valeurs par défaut mutables sont partagées entre tous les appels. Ce pattern est la manière correcte de gérer ce cas. C'est subtil et souvent oublié ; c'est bien de l'avoir fait spontanément.

---

## Points bloquants pendant le TP

| Point bloquant | Temps perdu | Cause racine | À retenir |
|---|---|---|---|
| Score mismatch BioPython (delta ~150) | *estimation* | `scoring_matrix` non passé à `global_affine` dans le script de benchmark — convention identique, mais matrices différentes (match/mismatch vs NUC_SIMPLE) | Toujours vérifier que les deux côtés d'une comparaison utilisent exactement les mêmes paramètres avant de conclure à un bug |
| `ValueError: sequence contains letters not in the alphabet` | — | BioPython exige un `str` quand `substitution_matrix` est défini | Fixer avec `str(seq)` |

---

## Leçons pour les prochaines estimations

- TP06 terminé rapidement malgré la composante BioPython (intégration externe, convention à identifier) : la capacité à diagnostiquer des différences de convention sur des APIs tierces est maintenant démontrée.
- Le debugging de la convention gap (trouver le bon facteur d'ajustement, le vérifier algébriquement) a pris plus de temps que l'implémentation elle-même — pattern courant en bioinformatique : les outils de référence ont des conventions implicites mal documentées.

---

## Suivi des corrections

| Point | Priorité | Appliqué | Note |
|-------|----------|----------|------|
| `if seed:` → `if seed is not None:` | 🟠 | [ ] | |
| Renommer param `mode` dans `make_biopython_aligner` | 🟡 | [ ] | |
| Ajouter `str()` dans `biopython_score` et lambdas | 🟡 | [ ] | |

**Prêt pour la clôture :** après correction du 🟠 ci-dessus.
