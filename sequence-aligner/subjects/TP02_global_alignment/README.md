 # TP02 — Alignement global : Needleman-Wunsch

**Projet :** sequence-aligner
**Prérequis :** TP01 — Représentation des séquences + parsing FASTA
**Compétences travaillées :** ALI, PD
**Durée estimée :** 5–7 h
**Durée réelle :** *à remplir*
**Fichier à implémenter :** `src/needleman_wunsch.py`

---

## Objectifs

À la fin de ce TP, tu sauras :

- Expliquer ce qu'est un alignement global et pourquoi il diffère d'un alignement local
- Implémenter l'algorithme de Needleman-Wunsch en deux phases : remplissage de la matrice + traceback
- Lire et interpréter une matrice de programmation dynamique (DP)
- Connaître la complexité temporelle et spatiale de l'algorithme

---

## Contexte biologique

### Pourquoi aligner des séquences ?

Deux espèces proches — par exemple l'humain et la souris — partagent des ancêtres communs. Au fil de l'évolution, leurs génomes ont divergé par **mutations ponctuelles** (une base remplacée par une autre) et par **indels** (insertions ou délétions de bases). Si on séquence le gène BRCA1 chez les deux espèces, on obtient deux séquences similaires mais pas identiques.

L'**alignement global** répond à la question : *quelle est la correspondance optimale, position par position, entre deux séquences homologues ?* Il suppose que les deux séquences sont comparables dans leur intégralité — du début à la fin — ce qui est typique pour des gènes orthologues de longueur similaire.

**Exemple concret :** En oncologie, on compare le gène BRCA1 d'une patiente à la séquence de référence humaine pour détecter des mutations héréditaires liées au cancer du sein. L'alignement global révèle exactement où les bases diffèrent.

### Le modèle de score

Pour quantifier la qualité d'un alignement, on attribue un score à chaque colonne :
- **Match** (même base face à face) → score positif, ex. `+1`
- **Mismatch** (bases différentes face à face) → pénalité, ex. `-1`
- **Gap** (un `-` face à une base, représentant un indel) → pénalité plus forte, ex. `-2`

L'objectif est de trouver l'alignement de score maximum — ce qui revient à choisir le meilleur compromis entre maximiser les matches et minimiser les gaps/mismatches.

---

## Problème à résoudre

**Entrée :** deux séquences ADN `seq1` et `seq2`, des paramètres de score (`match`, `mismatch`, `gap`).

**Sortie :** le score optimal, et les deux séquences alignées de même longueur (avec `-` pour les gaps).

**Exemple :**
```
seq1 = "AGCT"
seq2 = "AGT"
match=+1, mismatch=-1, gap=-2

Output :
  score    = 1
  aligned1 = "AGCT"
  aligned2 = "AG-T"

Lecture : A↔A (+1), G↔G (+1), C↔- (-2), T↔T (+1) → total = 1
```

---

## Algorithme

L'algorithme de Needleman-Wunsch (1970) repose sur la **programmation dynamique** : on divise le problème "aligner deux séquences de longueur m et n" en sous-problèmes "aligner les i premiers caractères de seq1 avec les j premiers de seq2".

Il y a deux phases bien séparées — implémente-les comme deux fonctions distinctes (`fill_matrix` et `traceback`).

### Phase 1 — Remplissage de la matrice

```
ALGORITHME fill_matrix(seq1, seq2, match_score, mismatch_score, gap_penalty)
  m ← longueur de seq1
  n ← longueur de seq2

  // Créer une matrice (m+1) lignes × (n+1) colonnes, initialisée à 0
  F ← matrice (m+1) × (n+1)

  // Ligne 0 et colonne 0 : coût d'aligner contre une séquence vide = que des gaps
  POUR i DE 0 À m : F[i][0] ← i × gap_penalty
  POUR j DE 0 À n : F[0][j] ← j × gap_penalty

  // Remplissage cellule par cellule
  POUR i DE 1 À m :
    POUR j DE 1 À n :
      SI seq1[i-1] == seq2[j-1] ALORS
        diag ← F[i-1][j-1] + match_score     // match
      SINON
        diag ← F[i-1][j-1] + mismatch_score  // mismatch
      FIN SI

      haut   ← F[i-1][j] + gap_penalty       // gap dans seq2
      gauche ← F[i][j-1] + gap_penalty       // gap dans seq1

      F[i][j] ← max(diag, haut, gauche)
  FIN POUR

  RETOURNER F

// Complexité : O(m × n) temps, O(m × n) espace
```

### Phase 2 — Traceback

Le score optimal est en `F[m][n]`. Pour reconstruire l'alignement, on **remonte** la matrice jusqu'à `F[0][0]` en suivant d'où vient chaque valeur.

```
ALGORITHME traceback(F, seq1, seq2, match_score, mismatch_score, gap_penalty)
  i ← m (= len(seq1)), j ← n (= len(seq2))
  res1 ← liste vide, res2 ← liste vide

  TANT QUE i > 0 OU j > 0 :
    SI i > 0 ET j > 0 :
      diag_score ← match_score SI seq1[i-1]==seq2[j-1], mismatch_score SINON
      SI F[i][j] == F[i-1][j-1] + diag_score ALORS
        // On vient du diagonal → match ou mismatch
        res1.ajouter_début(seq1[i-1])
        res2.ajouter_début(seq2[j-1])
        i ← i-1 ; j ← j-1
      SINON SI F[i][j] == F[i-1][j] + gap_penalty ALORS
        // On vient du haut → gap dans seq2
        res1.ajouter_début(seq1[i-1])
        res2.ajouter_début('-')
        i ← i-1
      SINON
        // On vient de gauche → gap dans seq1
        res1.ajouter_début('-')
        res2.ajouter_début(seq2[j-1])
        j ← j-1
    SINON SI i > 0 ALORS   // bord gauche : épuiser seq1
      res1.ajouter_début(seq1[i-1])
      res2.ajouter_début('-')
      i ← i-1
    SINON                   // bord haut : épuiser seq2
      res1.ajouter_début('-')
      res2.ajouter_début(seq2[j-1])
      j ← j-1
  FIN TANT QUE

  RETOURNER ("".join(res1), "".join(res2))

// En cas d'égalité entre plusieurs directions : priorité diagonal > haut > gauche
// Indice Python : construire res1 et res2 comme des listes en ajoutant à la fin,
//   puis inverser avec [::-1] avant de joindre — plus efficace que de préfixer des strings
```

---

## Fichiers

| Fichier | Rôle |
|---------|------|
| `src/needleman_wunsch.py` | `fill_matrix`, `traceback`, `needleman_wunsch` à implémenter |
| `tests/test_needleman_wunsch.py` | Tests automatisés — ne pas modifier |

---

## Critères de validation

Cocher chaque case avant de dire **"TP validé"** :

- [ ] `uv run pytest tests/test_needleman_wunsch.py -v` passe entièrement
- [ ] Je peux tracer à la main la matrice de l'exemple AGCT / AGT et retrouver le score 1
- [ ] Je sais expliquer pourquoi la première ligne et la première colonne sont initialisées avec `i × gap`
- [ ] Je sais expliquer la différence entre "venir du haut" et "venir de gauche" dans le traceback
- [ ] Je connais la complexité O(m×n) en temps et en espace

---

## Pistes si tu bloques

> Essayer dans l'ordre.

1. **Matrice mal initialisée** : trace à la main sur papier la matrice 5×4 pour AGCT/AGT avant de coder. Si ta matrice code et ta matrice papier divergent dès la 2e ligne, l'erreur est dans l'initialisation.
2. **Traceback qui boucle ou sort des bornes** : vérifie que ta condition `TANT QUE` couvre bien le cas où *i ou j* (pas seulement les deux) est encore positif. Un des deux peut atteindre 0 avant l'autre.
3. **Score correct mais alignement faux** : le traceback ne reconstruit pas depuis la bonne cellule, ou la priorité diagonal/haut/gauche n'est pas respectée. Ajoute un `print` de la cellule et de la direction choisie à chaque étape.

---

## Pour aller plus loin *(optionnel)*

- L'implémentation naïve utilise O(m×n) d'espace pour stocker toute la matrice. L'algorithme de **Hirschberg (1975)** réduit l'espace à O(min(m,n)) avec la même complexité temporelle — comment ?
- Que se passe-t-il si on remplace le `max` par un `min` et qu'on inverse les signes des scores ? (Réponse : distance d'édition de Levenshtein)
- Implémenter une visualisation ASCII de la matrice avec la direction de chaque cellule.

---

## Suivi de temps

*À remplir par l'étudiant au fil du TP.*

| | Durée |
|---|---|
| Estimée | 5–7 h |
| Réelle | |

### Points bloquants rencontrés

| Point bloquant | Temps perdu estimé | Résolution |
|---|---|---|
| | | |
