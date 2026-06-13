# Review — TP01 : Représentation des séquences + parsing FASTA

**Projet :** sequence-aligner
**Date :** 2026-06-13
**Temps estimé / réel :** 3–4 h / ~25 min *(écart : -2h35 à -3h35)*

---

## 🔴 Bloquant — corriger avant le TP suivant

Aucun point bloquant — tous les tests passaient dès le départ.

---

## 🟠 Dette — rembourser maintenant

| # | Problème | Impact | Correction suggérée |
|---|----------|--------|---------------------|
| 1 | `write_fasta` accédait à `seq._header` (attribut privé depuis l'extérieur de la classe) | Encapsulation cassée — renommer `_header` aurait créé un bug silencieux | Ajouter `@property header` dans `Sequence` et utiliser `seq.header` dans `write_fasta` |
| 2 | Message d'erreur trompeur dans `transcribe` : `"Sequence must be RNA"` | Diagnostic difficile à la lecture — le message contredisait la condition | Remplacé par `"Cannot transcribe: sequence is already RNA"` |
| 3 | `[nuc for nuc in data if nuc not in alphabet]` reportait les doublons | Message d'erreur verbeux et redondant pour des séquences avec beaucoup d'invalides | Utiliser `set(...)` pour dédupliquer |

---

## 🟡 Cosmétique — non urgent

- [x] `if not(isinstance(...))` → `if not isinstance(...)` (style Python idiomatique)
- [x] `"".join([...])` → `"".join(...)` (générateur suffit, évite une liste intermédiaire)

---

## Points bloquants pendant le TP

| Point bloquant | Temps perdu | Cause racine | À retenir pour les estimations |
|---|---|---|---|
| Test `test_example_file` en échec avec `FileNotFoundError` | ~quelques minutes | Pytest lancé depuis la racine du projet au lieu du dossier TP | Toujours lancer `uv run pytest test_skeleton.py -v` depuis le dossier du TP |

---

## Leçons pour les prochaines estimations

- TP01 bouclé en ~25 min : estimation initiale (3–4 h) très largement surestimée pour ce niveau Python. Revoir les fourchettes à la baisse pour les TPs algorithmiques purs.
- Seul point de friction : le répertoire de travail de pytest — noter dans chaque README la commande exacte de lancement.

---

## Suivi des corrections

| Point | Priorité | Appliqué | Note |
|-------|----------|----------|------|
| Propriété `header` + suppression de `seq._header` dans `write_fasta` | 🟠 | [x] | |
| Message d'erreur `transcribe` | 🟠 | [x] | |
| Déduplication des caractères invalides avec `set()` | 🟠 | [x] | |
| `not isinstance` sans parenthèses | 🟡 | [x] | |
| Générateur dans `complement` | 🟡 | [x] | |

**Prêt pour le TP suivant :** [x]
