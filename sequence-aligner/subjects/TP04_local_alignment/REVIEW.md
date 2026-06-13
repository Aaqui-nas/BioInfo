# Review — TP04 : Alignement local (Smith-Waterman)

**Projet :** sequence-aligner
**Date :** 2026-06-13
**Temps estimé / réel :** 4–6 h / ~15 min *(copier-coller de NW + modifications ciblées)*

---

## 🔴 Bloquant — corriger avant le TP suivant

| # | Problème | Impact | Correction suggérée |
|---|----------|--------|---------------------|
| 1 | `smith_waterman` retourne `F[len(seq1)][len(seq2)]` (coin bas-droite) au lieu de `max_score` | En SW, le score optimal est le **maximum de toute la matrice**, pas le coin bas-droite. Les tests passent par accident : tous les exemples ont leur alignement local qui se termine exactement à la fin de `seq1`. Contre-exemple qui casse : `seq1="ACGTTTTT"`, `seq2="ACGT"` — alignement optimal au milieu, coin bas-droite vaut bien moins que 4. | Remplacer `return F[len(seq1)][len(seq2)], ...` par `return max_score, ...` — `find_max` retourne déjà `max_score` à la ligne précédente. |

---

## 🟠 Dette — rembourser maintenant

| # | Problème | Impact | Correction suggérée |
|---|----------|--------|---------------------|
| 1 | Commentaire stale dans `traceback` : `# TP03 : remplacer ce calcul de diag_score par scoring_matrix.score() si fourni` — copié du squelette NW, mais le code gère déjà correctement `scoring_matrix`. | Un lecteur futur (ou toi dans 2 mois) voit ce commentaire et pense qu'il manque quelque chose à implémenter. Le code est correct ; le commentaire est faux. | Supprimer la ligne de commentaire. |

---

## 🟡 Cosmétique — non urgent

- [ ] `find_max` utilise `len(matrix[0])` pour la largeur de la matrice. Fonctionne car la matrice est rectangulaire, mais `len(matrix[i])` serait plus défensif (évite un IndexError si jamais la matrice était jagged).

---

## Points bloquants pendant le TP

| Point bloquant | Temps perdu | Cause racine | À retenir pour les estimations |
|---|---|---|---|
| — | — | Copier-coller NW → modifications minimes, aucun blocage | — |

---

## Leçons pour les prochaines estimations

- Le copier-coller de NW a été efficace ici parce que SW est structurellement très proche. En revanche, copier-coller peut introduire des bugs subtils si on oublie d'adapter *toutes* les différences — c'est exactement ce qui s'est passé avec le score retourné.
- **Règle à retenir :** quand tu copies une fonction, liste explicitement toutes les différences avec l'original avant de commencer. Pour SW vs NW : 3 différences (init, plancher, score retourné). La troisième a été manquée.
- Les tests de la suite de test TP04 ne détectent pas ce bug car ils ont tous été construits avec l'alignement local en fin de `seq1`. Bonne illustration de l'importance de tester des cas où le motif est **au milieu** de la séquence.

---

## Suivi des corrections

| Point | Priorité | Appliqué | Note |
|-------|----------|----------|------|
| `smith_waterman` : `F[len(seq1)][len(seq2)]` → `max_score` | 🔴 | [x] | |
| Supprimer commentaire stale `# TP03` dans `traceback` | 🟠 | [x] | |
| `find_max` : `len(matrix[0])` → `len(matrix[i])` | 🟡 | [x] | |

**Prêt pour le TP suivant :** [x]
