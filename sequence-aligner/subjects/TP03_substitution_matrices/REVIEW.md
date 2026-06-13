# Review — TP03 : Matrices de substitution

**Projet :** sequence-aligner
**Date :** 2026-06-13
**Temps estimé / réel :** 3–4 h / ~30 min *(TP rapide, même pattern que TP02)*

---

## 🔴 Bloquant — corriger avant le TP suivant

| # | Problème | Impact | Correction suggérée |
|---|----------|--------|---------------------|
| 1 | `score()` utilise `.get() or exec(...)` — cassé si un score vaut `0` | En TP04 (Smith-Waterman), la matrice peut avoir des scores nuls. `get()` retourne `0` (falsy) → le `or` déclenche le `exec` → `KeyError` levée à tort pour une paire valide. Les tests actuels passent uniquement parce qu'aucun score ne vaut `0`. | `return self._data[(a, b)]` — l'indexation directe lève `KeyError` naturellement si la clé est absente, sans cas limite. |

---

## 🟠 Dette — rembourser maintenant

| # | Problème | Impact | Correction suggérée |
|---|----------|--------|---------------------|
| 1 | Dans `from_file`, la variable de boucle s'appelle `score` — même nom qu'une méthode de la classe | Pas de bug ici (c'est un classmethod sans `self`), mais toute future lecture du code fait un double-take : "est-ce la méthode ou la variable ?" Crée de la friction à la relecture. | Renommer en `val` : `for col, val in zip(columns, values): scores[(row, col)] = val` |
| 2 | Dans `nuc_simple`, la variable s'appelle `pairs` au lieu de `transitions` | `pairs` est générique — ça ne dit pas *pourquoi* ces 4 combinaisons ont un score différent. Si tu reviens dans 3 mois, tu dois te rappeler mentalement que `pairs` = transitions. Le nom devrait encoder l'intention. | Renommer : `transitions = frozenset({("A","G"), ("G","A"), ("T","C"), ("C","T")})` (le frozenset règle aussi le point cosmétique ci-dessous) |

---

## 🟡 Cosmétique — non urgent

- [ ] `nuc_simple` : `pairs` est un tuple → la vérification `(nuc1, nuc2) in pairs` est O(n). Avec 4 éléments c'est négligeable, mais `frozenset` est la structure correcte pour les tests d'appartenance — O(1) et sémantiquement plus juste (ordre non pertinent, pas de doublons).
- [ ] `from_file` : lire toutes les lignes dans `lines = []` avant de les traiter oblige à garder le fichier entier en mémoire. Pas de problème à cette taille, mais le pattern standard est de parser dans le bloc `with` sans collecte intermédiaire.

---

## Points bloquants pendant le TP

| Point bloquant | Temps perdu | Cause racine | À retenir pour les estimations |
|---|---|---|---|
| — | — | TP fluide, pas de blocage rapporté | — |

---

## Leçons pour les prochaines estimations

- TP03 a pris ~30 min comme TP02 — le pattern "classe Python + intégration dans un algo existant" est maintenant bien maîtrisé. Les estimations TP04/05 (Smith-Waterman, gap affines) peuvent rester à 4–6 h car la nouveauté algorithmique augmente.
- Le bug du `score()` est un cas classique du "ça passe les tests mais c'est fragile" : les tests ne couvrent pas les valeurs `0` parce que la matrice fournie n'en a pas. Bonne illustration de pourquoi les tests doivent couvrir les valeurs aux bords (zéro inclus), pas seulement les cas nominaux.

---

## Suivi des corrections

| Point | Priorité | Appliqué | Note |
|-------|----------|----------|------|
| `score()` → `return self._data[(a, b)]` | 🔴 | [x] | |
| `from_file` : renommer `score` → `val` en boucle | 🟠 | [x] | |
| `nuc_simple` : renommer `pairs` → `transitions`, passer en `frozenset` | 🟠 | [x] | |
| `nuc_simple` : tuple → `frozenset` (si pas déjà fait avec le point précédent) | 🟡 | [x] | |
| `from_file` : parser dans le `with` sans liste intermédiaire | 🟡 | [x] | |

**Prêt pour le TP suivant :** [x]
