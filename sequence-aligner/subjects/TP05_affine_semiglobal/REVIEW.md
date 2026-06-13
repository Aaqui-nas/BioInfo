# Review — TP05 : Gap penalties affines + alignement semi-global

**Projet :** sequence-aligner
**Date :** 2026-06-13
**Temps estimé / réel :** 4–6 h / ~40 min *(écart : −5 h environ)*

---

## 🔴 Bloquant — corriger avant le TP suivant

Aucun.

---

## 🟠 Dette — rembourser maintenant

Aucune.

---

## 🟡 Cosmétique — non urgent

- [ ] **Ligne 214 redondante dans `fill_matrices_semiglobal`** : `H[0][0] = 0` est posé seul, puis écrasé par le `for j in range(n + 1): H[0][j] = 0` juste en dessous. La ligne 214 n'a aucun effet. À supprimer.

- [ ] **`stop()` closure sur variables mutables** (lignes 44–47 de `_traceback_affine`) : la condition d'arrêt est encapsulée dans une fonction imbriquée qui capture `i` et `j` — des variables réassignées dans la boucle. Ça fonctionne correctement en Python (les closures lisent la valeur courante au moment de l'appel), mais le pattern surprend les lecteurs non avertis. La condition inline est plus immédiate :
  ```python
  while i > 0 if stop_when_i_zero else (i > 0 or j > 0):
  ```

- [ ] **Double ligne vide** entre `fill_matrices_semiglobal` et `find_semiglobal_end` : une seule suffit.

---

## Ce qui est bien fait

**Abstraction `_traceback_affine` avec `stop_when_i_zero`** — c'est le point fort du code. Le traceback global et le traceback semi-global ne diffèrent que par leur point d'arrêt ; les factoriser en une seule fonction avec un flag booléen évite de dupliquer ~60 lignes de logique d'état (H/E/F). C'est un bon jugement : ni trop peu (duplication), ni trop (sur-abstraction).

**`_score` helper** — extraire le choix scoring_matrix vs match/mismatch en une ligne privée est la bonne décision. La branche if/else apparaissait 4 fois dans le code de remplissage + traceback ; un helper la réduit à 1. C'est exactement le seuil à partir duquel une factorisation vaut la peine.

**Initialisation des matrices comme `_NEG_INF` par défaut** — initialiser toute la matrice à −∞ puis ne corriger que les cellules qui le nécessitent est plus sûr que d'initialiser à 0 et de gérer les bords comme cas spéciaux. Ça garantit que toute cellule non remplie correctement produit un score sentinel reconnaissable plutôt qu'une valeur silencieusement incorrecte.

---

## Points bloquants pendant le TP

*(aucun signalé par l'étudiant)*

| Point bloquant | Temps perdu | Cause racine | À retenir |
|---|---|---|---|
| — | — | — | — |

---

## Leçons pour les prochaines estimations

- TP05 terminé en ~40 min pour 4–6 h estimés : le pattern est maintenant établi sur 5 TPs (TP01-04 entre 15 et 30 min chacun). La maîtrise de la programmation dynamique est acquise ; les prochaines estimations pour des TPs PD pur devraient être divisées par ~6 par rapport aux fourchettes standards.
- La difficulté réelle du TP06 (benchmarking + visualisation) dépendra du temps passé à manier matplotlib et BioPython, pas de l'algo — majorer cette partie.

---

## Suivi des corrections

| Point | Priorité | Appliqué | Note |
|-------|----------|----------|------|
| Ligne 214 redondante dans `fill_matrices_semiglobal` | 🟡 | [ ] | |
| `stop()` closure → condition inline | 🟡 | [ ] | |
| Double ligne vide | 🟡 | [ ] | |

**Prêt pour le TP suivant :** [x] *(aucun 🔴 ni 🟠)*
