# Review — TP02 : Alignement global Needleman-Wunsch

**Projet :** sequence-aligner
**Date :** 2026-06-13
**Temps estimé / réel :** 5–7 h / ~30 min *(écart : -4h30 à -6h30)*

---

## 🔴 Bloquant — corriger avant le TP suivant

Aucun point bloquant — tous les tests passaient dès l'implémentation.

---

## 🟠 Dette — rembourser maintenant

| # | Problème | Impact | Correction suggérée |
|---|----------|--------|---------------------|
| 1 | `np.zeros` crée des `float64` → signature `list[list[int]]` violée, score retourné en `float` | Comparaisons float/int dans `traceback` : risque silencieux sur des valeurs proches ; type de retour incorrect | Remplacé par `[[0] * (n+1) for _ in range(m+1)]` |
| 2 | `seq1[i-1]` suppose `__getitem__` sur `Sequence`, non prévu dans TP01 | `IndexError` ou `TypeError` si `Sequence` évolue sans `__getitem__` | Ajout de `__getitem__` dans TP01 (rétroactif) — bonne décision |

---

## 🟡 Cosmétique — non urgent

- [ ] `range (m+1)` → `range(m+1)` (espace parasite avant la parenthèse)
- [ ] Lignes vides manquantes entre les fonctions top-level (`def traceback`, `def needleman_wunsch`) — PEP8 demande 2 lignes
- [ ] `len(seq1),len(seq2)` et `res1, res2 = [],[]` → espaces après les virgules

---

## Points bloquants pendant le TP

| Point bloquant | Temps perdu | Cause racine | À retenir pour les estimations |
|---|---|---|---|
| Aucun signalé | — | — | — |

---

## Leçons pour les prochaines estimations

- Deuxième TP bouclé en ~30 min au lieu de 5–7 h : les estimations initiales sont très largement surestimées pour ce niveau Python. À partir de TP03, réduire les fourchettes de 70–80 % pour les TPs algorithmiques purs.
- L'ajout rétroactif de `__getitem__` à `Sequence` est une bonne réflexion : une classe qui représente une séquence devrait se comporter comme une séquence Python.
- Le choix numpy était naturel (matrice 2D) mais introduit un bug de type silencieux — réflexe à garder : numpy pour du vectoriel, listes Python pour du cell-by-cell.

---

## Suivi des corrections

| Point | Priorité | Appliqué | Note |
|-------|----------|----------|------|
| Remplacer `np.zeros` par liste de listes | 🟠 | [x] | |
| Ajouter `__getitem__` à `Sequence` (TP01) | 🟠 | [x] | Extension rétroactive justifiée |

**Prêt pour le TP suivant :** [x]
