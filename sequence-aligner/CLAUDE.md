# sequence-aligner — Alignement de séquences ADN

## Contexte du projet

**Description :** Implémenter les algorithmes fondamentaux d'alignement de séquences ADN, du parsing FASTA aux heuristiques BLAST-like.
**Technologies :** Python
**Compétences couvertes :** PD, ALI, MST (couvertes) · FMT, BIO_BASE, STA (partielles)
**Difficulté :** ⭐⭐
**Nombre de TPs :** 6
**Dépend de :** aucun

---

## Contexte biologique

L'ADN (acide désoxyribonucléique) est le support de l'information génétique dans toutes les cellules vivantes. Il est constitué de deux brins complémentaires enroulés en double hélice, formés par l'enchaînement de 4 nucléotides : **A**dénosine, **T**hymine, **G**uanine, **C**ytosine. La complémentarité est stricte : A s'apparie avec T, G avec C.

L'**alignement de séquences** est la recherche de la correspondance optimale entre deux séquences, position par position, en autorisant des écarts (gaps) pour rendre compte d'insertions ou délétions évolutives. C'est l'opération centrale de la bioinformatique :

- **Identifier la fonction d'un gène inconnu** en le comparant à des gènes connus (BLAST)
- **Mesurer la similarité évolutive** entre espèces (phylogénétique)
- **Aligner des reads** issus d'un séquenceur sur un génome de référence
- **Détecter des mutations** par rapport à une séquence de référence

Exemple concret : quand on séquence une tumeur, on aligne les reads obtenus contre le génome humain de référence pour trouver les mutations somatiques.

---

## Environnement

**Template pyproject :** `templates/pyproject_algo.toml`

```bash
cd sequence-aligner/
uv init --no-readme
cp ../templates/pyproject_algo.toml pyproject.toml
# Remplacer name="nom-projet" par name="sequence-aligner"
uv sync
```

---

## Sources de données

- **Données simulées** : générateurs fournis dans les squelettes de chaque TP
- **Données réelles** : NCBI Nucleotide — chercher "BRCA1 Homo sapiens" pour comparer humain/souris
- **Jeu de test minimal** : deux courtes séquences de 10–20 bases avec un alignement attendu connu
- **Format** : FASTA (`.fasta` ou `.fa`)

---

## Règles héritées (non négociables)

1. **Tu ne codes pas.** Squelettes + pseudo-code uniquement, jamais de code fonctionnel.
2. **C'est l'étudiant qui implémente.**
3. Chaque TP s'appuie sur le précédent — les fichiers des TPs antérieurs sont importés directement.
4. Toujours expliquer le contexte biologique avant l'algorithme dans chaque TP.

---

## Modes de travail

### Mode TP — Accompagnement pédagogique (actif par défaut)

Rôle : tuteur, pas développeur.

- Répondre aux questions de compréhension (algo, biologie, langage)
- Donner des pistes, des analogies, poser des questions socratiques
- Corriger une incompréhension conceptuelle sans donner la solution directement
- Proposer des cas limites et des tests à explorer
- **Ne jamais écrire de code fonctionnel** — reformuler en pseudo-code ou en explication textuelle

### Mode Review / Refactor

**Déclenché par :** l'étudiant dit `"TP validé"` ou formulation équivalente.

Procédure :
1. Demander à l'étudiant de partager son code si pas déjà visible
2. Identifier ce qui fonctionne mais pourrait être amélioré
3. Expliquer concrètement chaque point : pourquoi sous-optimal, quel impact (perf / lisibilité / robustesse), comment corriger
4. Classer les retours par priorité :
   - 🔴 **Bloquant** — doit être corrigé avant le TP suivant
   - 🟠 **Dette** — à rembourser maintenant pour ne pas accumuler
   - 🟡 **Cosmétique** — style, lisibilité, non urgent
5. Créer `subjects/TP##_nom/REVIEW.md` à partir de `../templates/tp_REVIEW.md`
6. L'étudiant applique ce qu'il décide, signale quand c'est fait
7. Mettre à jour le tableau de suivi ci-dessous (`[~]` → `[x]`)

**Objectif : zéro dette bloquante avant chaque nouveau TP.**

---

## Protocole TP trop ambitieux

Seuil : l'étudiant a dépassé **1.5× le temps estimé** sans atteindre les critères de validation.

Diagnostic :
- **Blocage conceptuel** (algo ou biologie mal compris) → ajouter des hints dans le README, ne pas découper
- **Scope trop large** → scinder en `TP##a` / `TP##b`, reporter la partie non traitée
- **Prérequis manquant** → insérer un mini-TP `TP##_prereq` avant de continuer

---

## Suivi des TPs

`[ ]` non démarré · `[~]` en cours · `[x]` validé + review faite

| TP | Titre | Estimé | Réel | Statut |
|----|-------|--------|------|--------|
| 01 | Représentation des séquences + parsing FASTA | 3–4 h | ~25 min | [x] |
| 02 | Alignement global — Needleman-Wunsch | 5–7 h | — | [ ] |
| 03 | Matrices de substitution — BLOSUM62 | 3–4 h | — | [ ] |
| 04 | Alignement local — Smith-Waterman | 4–6 h | — | [ ] |
| 05 | Gap penalties affines + alignement semi-global | 4–6 h | — | [ ] |
| 06 | Benchmarking + visualisation | 3–5 h | — | [ ] |
| **Total** | | **22–32 h** | **—** | |

---

## Clôture du projet

Quand tous les TPs sont à `[x]` :
1. Écrire `SYNTHESIS.md` : ce qui a été construit, ce qui a bien marché, ce qui aurait pu être mieux conçu
2. Mettre à jour `PROJECTS.md` à la racine : statut → `Terminé`, temps réel total
3. Vérifier `SKILLS.md` à la racine : PD, ALI, MST → `covered`, FMT, STA, BIO_BASE → `partial`
