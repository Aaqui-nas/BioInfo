# TP04 — Alignement local : Smith-Waterman

**Projet :** sequence-aligner
**Prérequis :** TP02 — Needleman-Wunsch · TP03 — Matrices de substitution
**Compétences travaillées :** ALI, PD
**Durée estimée :** 4–6 h
**Durée réelle :** *à remplir*
**Fichier à implémenter :** `src/smith_waterman.py`

---

## Objectifs

À la fin de ce TP, tu sauras :

- Expliquer la différence entre alignement global et alignement local
- Identifier les deux modifications de NW qui donnent SW (initialisation + plancher)
- Implémenter `fill_matrix`, `find_max`, `traceback` et `smith_waterman`
- Montrer sur un exemple pourquoi SW est strictement meilleur que NW pour trouver des motifs locaux

---

## Contexte biologique

### Quand l'alignement global échoue

Needleman-Wunsch est conçu pour aligner deux séquences *dans leur intégralité*. C'est adapté lorsque les deux séquences sont homologues sur toute leur longueur (gènes orthologues de taille similaire).

Mais en biologie, on est souvent dans un cas différent : **chercher une région conservée dans un contexte plus large**.

Exemples concrets :
- Identifier un **domaine fonctionnel** (ex : domaine kinase, domaine de liaison à l'ADN) dans une protéine inconnue en le comparant à une base de données de domaines.
- Trouver un **exon** d'un gène dans une séquence génomique qui contient aussi des introns et des régions intergéniques.
- **BLAST** — l'outil le plus utilisé en bioinformatique — est fondé sur Smith-Waterman : il cherche des régions de haute similarité entre une requête et des millions de séquences dans une base de données.

Si on utilise NW pour aligner `GGGACGT` avec `ACGT`, l'alignement global est forcé de payer des pénalités pour les `GGG` qui n'ont pas de correspondant — le score final est négatif. SW, lui, trouve `ACGT/ACGT` et ignore le reste.

### La clé biologique : les substitutions locales fortes valent mieux qu'une couverture globale forcée

Un gène conservé entre deux espèces peut être flanqué de régions très divergentes. Forcer l'alignement de ces régions divergentes introduit du bruit et masque le signal réel. L'alignement local extrait proprement la région conservée.

---

## Problème à résoudre

**Entrée :** deux séquences ADN `seq1` et `seq2`, un `gap_penalty`, et optionnellement une `ScoringMatrix`.

**Sortie :** le score optimal de l'alignement local, et les deux sous-séquences alignées correspondantes.

**Exemple :**
```
seq1 = "TTTACGT"
seq2 = "ACGT"
match=+1, mismatch=-1, gap=-2

Output :
  score    = 4
  aligned1 = "ACGT"   ← sous-séquence de seq1, pas seq1 entière
  aligned2 = "ACGT"
```

---

## Algorithme

Smith-Waterman (1981) est une modification minime de Needleman-Wunsch. **Deux différences seulement** — mais elles changent complètement la nature du résultat.

### Différence 1 — Initialisation

```
NW : F[i][0] = i × gap_penalty,  F[0][j] = j × gap_penalty
SW : F[i][0] = 0,                 F[0][j] = 0
```

En NW, partir du coin signifie payer chaque gap pour "remonter" la séquence vide.
En SW, tout peut repartir de zéro — aucun historique forcé.

### Différence 2 — Plancher à zéro dans la récurrence

```
NW : F[i][j] = max(diag, haut, gauche)
SW : F[i][j] = max(0, diag, haut, gauche)
```

Si tous les chemins mènent à un score négatif, on préfère repartir de zéro plutôt que d'accumuler une dette. Conséquence : **aucune cellule de la matrice ne peut être négative**.

### Phase 1 — Remplissage

```
ALGORITHME fill_matrix(seq1, seq2, gap_penalty, match_score, mismatch_score, scoring_matrix)

  m ← len(seq1), n ← len(seq2)
  F ← matrice (m+1) × (n+1) initialisée à 0  ← tout à 0, pas de formule gap

  POUR i DE 1 À m :
    POUR j DE 1 À n :
      SI scoring_matrix fournie :
        s ← scoring_matrix.score(seq1[i-1], seq2[j-1])
      SINON :
        s ← match_score SI seq1[i-1] == seq2[j-1], mismatch_score SINON

      diag   ← F[i-1][j-1] + s
      haut   ← F[i-1][j]   + gap_penalty
      gauche ← F[i][j-1]   + gap_penalty

      F[i][j] ← max(0, diag, haut, gauche)   ← le plancher à 0
  FIN POUR

  RETOURNER F

// Complexité : O(m × n) temps, O(m × n) espace — identique à NW
```

### Phase 2 — Trouver le maximum

```
ALGORITHME find_max(matrix) → (row, col, score)

  max_score ← 0,  max_row ← 0,  max_col ← 0

  POUR i DE 0 À nb_lignes - 1 :
    POUR j DE 0 À nb_colonnes - 1 :
      SI matrix[i][j] > max_score :
        max_score ← matrix[i][j]
        max_row ← i,  max_col ← j
  FIN POUR

  RETOURNER (max_row, max_col, max_score)

// En cas d'égalité, la première cellule rencontrée est retenue.
```

### Phase 3 — Traceback

```
ALGORITHME traceback(matrix, seq1, seq2, start_row, start_col, ...)

  // Même logique que NW, deux différences :
  // - Point de départ : (start_row, start_col)  au lieu de (m, n)
  // - Condition d'arrêt : matrix[i][j] == 0     au lieu de (i == 0 ET j == 0)

  i ← start_row,  j ← start_col
  res1 ← [],  res2 ← []

  TANT QUE matrix[i][j] != 0 :
    SI i > 0 ET j > 0 :
      // Calculer diag_score (même logique que dans fill_matrix)
      SI F[i][j] == F[i-1][j-1] + diag_score :
        res1.ajouter(seq1[i-1]) ; res2.ajouter(seq2[j-1]) ; i← i-1 ; j← j-1
      SINON SI F[i][j] == F[i-1][j] + gap_penalty :
        res1.ajouter(seq1[i-1]) ; res2.ajouter('-') ; i← i-1
      SINON :
        res1.ajouter('-') ; res2.ajouter(seq2[j-1]) ; j← j-1
    SINON SI i > 0 :
      res1.ajouter(seq1[i-1]) ; res2.ajouter('-') ; i← i-1
    SINON :
      res1.ajouter('-') ; res2.ajouter(seq2[j-1]) ; j← j-1
  FIN TANT QUE

  RETOURNER ("".join(res1 inversé), "".join(res2 inversé))

// Si start_row == 0 ou start_col == 0 ou matrix[start_row][start_col] == 0 :
//   retourner ("", "") directement.
```

---

## Fichiers

| Fichier | Rôle |
|---------|------|
| `src/smith_waterman.py` | `fill_matrix`, `find_max`, `traceback`, `smith_waterman` à implémenter |
| `tests/test_smith_waterman.py` | 13 tests à faire passer — ne pas modifier |

---

## Critères de validation

Cocher chaque case avant de dire **"TP validé"** :

- [ ] `uv run pytest tests/test_smith_waterman.py -v` passe entièrement
- [ ] `uv run pytest tests/test_needleman_wunsch.py tests/test_scoring.py -v` passe encore (pas de régression)
- [ ] Je peux tracer à la main la matrice SW pour `"ACGT"` / `"CGT"` et identifier la cellule maximale
- [ ] Je sais expliquer pourquoi la première ligne/colonne est à 0 en SW (et pas `i × gap`)
- [ ] Je sais expliquer pourquoi le traceback s'arrête à 0 (et non en haut à gauche)
- [ ] Je connais un cas concret où SW donne un meilleur résultat que NW

---

## Pistes si tu bloques

> Essayer dans l'ordre.

1. **`fill_matrix` retourne des valeurs négatives** : tu as oublié le `max(0, ...)`. C'est la seule modification dans la récurrence par rapport à NW.

2. **`find_max` retourne (0, 0, 0) même quand la matrice a des valeurs positives** : vérifie que tu parcours bien *toutes* les cellules, y compris les cellules intérieures. Un `for i in range(len(matrix))` + `for j in range(len(matrix[i]))` suffit.

3. **`traceback` produit un alignement trop long** : ta condition d'arrêt est probablement `i > 0 or j > 0` (NW) au lieu de `matrix[i][j] != 0` (SW). Le traceback SW s'arrête dès qu'on atteint un 0, peu importe si i ou j est encore positif.

4. **Le test `test_no_positive_alignment` échoue** : quand le score maximal est 0, `smith_waterman` doit retourner `(0, "", "")` directement. Ajoute ce cas spécial avant d'appeler `traceback`.

---

## Pour aller plus loin *(optionnel)*

- SW est O(m × n) en espace. L'implémentation banded (en bande) réduit l'espace à O(k × n) quand on sait à l'avance que les deux séquences ne divergent pas de plus de k positions — comment ?
- BLAST n'implémente pas SW exactement : il utilise des *seeds* (correspondances exactes de k-mers) puis étend localement. Pourquoi ? (Indice : base de données de milliards de bases.)
- Qu'est-ce qui change si on veut trouver *tous* les alignements locaux de score ≥ seuil, pas seulement le meilleur ?

---

## Suivi de temps

*À remplir par l'étudiant au fil du TP.*

| | Durée |
|---|---|
| Estimée | 4–6 h |
| Réelle | |

### Points bloquants rencontrés

| Point bloquant | Temps perdu estimé | Résolution |
|---|---|---|
| | | |
