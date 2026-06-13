# TP05 — Gap penalties affines + alignement semi-global

**Projet :** sequence-aligner
**Prérequis :** TP04 — Smith-Waterman · TP02 — Needleman-Wunsch
**Compétences travaillées :** ALI
**Durée estimée :** 4–6 h
**Durée réelle :** *à remplir*
**Fichier(s) à implémenter :** `src/affine_alignment.py`

---

## Objectifs

À la fin de ce TP, tu sauras :

- Expliquer pourquoi les pénalités de gap linéaires sont biologiquement inexactes
- Implémenter le modèle de gap affin (Gotoh 1982) avec 3 matrices (H, E, F)
- Effectuer un traceback avec suivi d'état matriciel — plus complexe que NW/SW
- Distinguer alignement global, local et semi-global, et choisir le bon selon le contexte
- Aligner un read court sur un génome long sans pénaliser les extrémités du génome

---

## Contexte biologique

### Pourquoi les gaps linéaires sont biologiquement inexacts

Dans les TPs précédents, une insertion ou délétion (indel) de k bases coûte k × gap_penalty. Ce modèle traite chaque base du gap comme un événement indépendant — or c'est biologiquement faux.

**La réalité évolutive :** la grande majorité des indels sont des événements uniques. Une erreur de réplication peut "glisser" de plusieurs bases d'un coup (slippage), un rétrotransposon peut s'insérer d'un bloc, ou une région chromosomique peut être dupliquée ou délétée en une seule recombinaison. Une délétion de 50 bases est *un* événement, pas 50 délétions indépendantes.

**Conséquence pour l'alignement :** avec un modèle linéaire, 3 gaps séparés d'1 base (coût 3×gap) sont aussi "probables" qu'un gap continu de 3 bases (coût 3×gap). Le modèle affin corrige ça en distinguant :

| Paramètre | Notation | Rôle |
|-----------|----------|------|
| `gap_open`   | gapO < 0 | Coût de l'*événement* : ouvrir un nouveau gap |
| `gap_extend` | gapE < 0 | Coût par *base* : étendre un gap déjà ouvert |

Un gap de longueur k coûte : **`gap_open + k × gap_extend`**

En pratique : `gap_open` est beaucoup plus négatif que `gap_extend` (ex : −10 vs −1), ce qui pénalise fortement l'ouverture d'un gap mais favorise les gaps continus longs face aux gaps multiples courts.

**Exemple concret :** avec gapO=−10, gapE=−1
- Un gap de 3 continu → −10 + 3×(−1) = **−13**
- Trois gaps de 1 séparés → 3×(−10 + 1×(−1)) = **−33**

L'algorithme de Gotoh (1982) intègre ce modèle avec seulement 3 matrices au lieu d'une.

### Alignement semi-global : le cas du séquençage de lecture courte

Un séquenceur Illumina produit des *reads* de ~150 bases. On veut les aligner sur un génome de 3 milliards de bases. L'alignement **global** est inapplicable : les ~3 milliards de bases de part et d'autre du read seraient pénalisées comme autant de gaps.

L'alignement **semi-global** (ou *fitting alignment*) résout ce problème :

- **Pas de pénalité** pour les gaps aux extrémités de seq2 (le génome)
- **Pénalités normales** pour les gaps à l'intérieur de seq1 (le read)

Le read "flotte" dans le génome : on cherche la position qui maximise le score d'alignement.

Exemples concrets :
- Alignement de reads Illumina sur hg38 (bowtie2, bwa-mem)
- Recherche d'un exon (~200 bp) dans une séquence génomique contenant des introns
- Alignment d'un primer PCR (20 bp) sur une séquence de quelques kilobases

---

## Problème à résoudre

**Partie 1 — Alignement global affin :**
Entrée : seq1, seq2, gap_open, gap_extend, (optionnel) ScoringMatrix
Sortie : (score, aligned_seq1, aligned_seq2) — alignement global avec pénalités affines

**Partie 2 — Alignement semi-global :**
Entrée : seq1 (read court), seq2 (génome long), gap_open, gap_extend
Sortie : (score, aligned_seq1, aligned_seq2) — meilleure position de seq1 dans seq2

**Exemples :**
```
Global affin (gap_open=-3, gap_extend=-1, match=2, mismatch=-10) :
  seq1 = "ACT", seq2 = "ACGT"
  Output : score=2,  aligned1="AC-T",  aligned2="ACGT"
  → un gap de 1 coûte -3-1=-4  (vs mismatch=-10 au même endroit)

Semi-global (gap_open=-5, gap_extend=-1, match=2, mismatch=-3) :
  seq1 = "CGT",  seq2 = "AAACGTCC"
  Output : score=6,  aligned1="CGT",  aligned2="CGT"
  → le read s'aligne parfaitement sur positions 3–5 du génome ; les flancs AAA et CC ne sont pas pénalisés
```

---

## Algorithme

### Gotoh 1982 — Trois matrices pour les gaps affins

L'idée clé : maintenir **3 matrices** représentant l'état de l'alignement à chaque cellule.

| Matrice | Signification |
|---------|---------------|
| **H[i][j]** | Meilleur score pour aligner seq1[0..i-1] et seq2[0..j-1] — tout état confondu |
| **E[i][j]** | Meilleur score se terminant par un gap dans seq2 — seq1[i-1] face à un gap (mouvement **vertical** ↓) |
| **F[i][j]** | Meilleur score se terminant par un gap dans seq1 — seq2[j-1] face à un gap (mouvement **horizontal** →) |

### Phase 1 — Remplissage : global

```
ALGORITHME fill_matrices_global(seq1, seq2, gap_open, gap_extend, match_score, mismatch_score, scoring_matrix)

  m ← len(seq1),  n ← len(seq2)
  H, E, F ← 3 matrices (m+1) × (n+1) initialisées à 0

  // --- Initialisation ---
  H[0][0] ← 0
  E[0][0] ← −∞
  F[0][0] ← −∞

  POUR i DE 1 À m :
    H[i][0] ← gap_open + i × gap_extend    // aligner i bases de seq1 contre rien
    E[i][0] ← gap_open + i × gap_extend    // gap dans seq2, même coût
    F[i][0] ← −∞                           // impossible : seq2 est vide, pas de mouvement horizontal

  POUR j DE 1 À n :
    H[0][j] ← gap_open + j × gap_extend    // aligner j bases de seq2 contre rien
    F[0][j] ← gap_open + j × gap_extend    // gap dans seq1, même coût
    E[0][j] ← −∞                           // impossible : seq1 est vide, pas de mouvement vertical

  // --- Récurrence ---
  POUR i DE 1 À m :
    POUR j DE 1 À n :

      // Gap dans seq2 : seq1[i-1] face à un gap (mouvement vertical)
      E[i][j] ← max(
        H[i-1][j] + gap_open + gap_extend,   // ← ouvrir un nouveau gap depuis H
        E[i-1][j] + gap_extend                // ← étendre un gap déjà ouvert
      )

      // Gap dans seq1 : seq2[j-1] face à un gap (mouvement horizontal)
      F[i][j] ← max(
        H[i][j-1] + gap_open + gap_extend,   // ← ouvrir un nouveau gap depuis H
        F[i][j-1] + gap_extend                // ← étendre un gap déjà ouvert
      )

      // Score diagonal : match ou mismatch
      SI scoring_matrix fournie :
        s ← scoring_matrix.score(seq1[i-1], seq2[j-1])
      SINON :
        s ← match_score  SI seq1[i-1] == seq2[j-1],  mismatch_score  SINON

      H[i][j] ← max(
        H[i-1][j-1] + s,   // ← appariement (diagonal)
        E[i][j],             // ← gap dans seq2
        F[i][j]              // ← gap dans seq1
      )

  RETOURNER (H, E, F)
// Complexité : O(m × n) temps, O(m × n) espace
```

### Phase 2 — Traceback : global

Le traceback affin est plus délicat que NW : il faut savoir dans quel **état** on se trouve (H, E ou F) pour reconstruire correctement les gaps continus.

```
ALGORITHME traceback_global(H, E, F, seq1, seq2, gap_open, gap_extend, ...)

  i ← m = len(seq1),  j ← n = len(seq2)
  res1, res2 ← [], []
  state ← "H"   // on commence toujours dans H — le score optimal est en H[m][n]

  TANT QUE i > 0 OU j > 0 :

    SI state == "H" :
      SI i > 0 ET j > 0 :
        s ← calculer_score_diagonal(seq1, seq2, i, j, ...)

        SI H[i][j] == H[i-1][j-1] + s :
          res1 += seq1[i-1] ; res2 += seq2[j-1]
          i← i-1 ; j← j-1 ; state← "H"           // diagonal, reste dans H

        SINON SI H[i][j] == E[i][j] :
          state← "E"                                // bascule dans E, sans bouger encore

        SINON :
          state← "F"                                // bascule dans F, sans bouger encore

      SINON SI i > 0 :
        res1 += seq1[i-1] ; res2 += '-' ; i← i-1  // bord gauche : gap dans seq2
      SINON :
        res1 += '-' ; res2 += seq2[j-1] ; j← j-1  // bord haut : gap dans seq1

    SINON SI state == "E" :
      // Gap dans seq2 — on consomme seq1[i-1] face à un gap — mouvement vertical ↑
      res1 += seq1[i-1] ; res2 += '-'
      SI E[i][j] == H[i-1][j] + gap_open + gap_extend :
        state← "H"   // on a ouvert ce gap depuis H
      SINON :
        state← "E"   // on continue un gap dans E
      i← i-1

    SINON SI state == "F" :
      // Gap dans seq1 — on consomme seq2[j-1] face à un gap — mouvement horizontal ←
      res1 += '-' ; res2 += seq2[j-1]
      SI F[i][j] == H[i][j-1] + gap_open + gap_extend :
        state← "H"   // on a ouvert ce gap depuis H
      SINON :
        state← "F"   // on continue un gap dans F
      j← j-1

  FIN TANT QUE

  RETOURNER ("".join(res1 inversé), "".join(res2 inversé))

// Note : mettre state à jour AVANT de décrémenter i ou j dans les branches E et F.
```

### Phase 1 (semi-global) — Deux lignes changent

```
// fill_matrices_semiglobal = fill_matrices_global,
// sauf pour l'initialisation de la première ligne :

//   REMPLACER :
H[0][j] ← gap_open + j × gap_extend    →    H[0][j] ← 0
F[0][j] ← gap_open + j × gap_extend    →    F[0][j] ← 0

// La récurrence intérieure est identique.
// La première colonne (H[i][0], E[i][0]) reste inchangée.
```

### Trouver la fin de l'alignement semi-global

```
ALGORITHME find_semiglobal_end(H) → (end_col, score)

  m ← nombre de lignes de H − 1  // H a m+1 lignes
  max_score ← −∞,  end_col ← 0

  POUR j DE 0 À n :
    SI H[m][j] > max_score :
      max_score ← H[m][j]
      end_col ← j

  RETOURNER (end_col, max_score)

// Parcourir seulement la dernière ligne H[m][*].
// end_col indique où, dans seq2, seq1 se termine.
```

### Phase 2 (semi-global) — Traceback

```
// traceback_semiglobal = traceback_global avec deux modifications :

// 1. Point de départ : (m, end_col) au lieu de (m, n)
i← m,  j← end_col,  state← "H"

// 2. Condition d'arrêt : i == 0  (seq1 entièrement consommée)
//                       au lieu de  i == 0 ET j == 0

TANT QUE i > 0 :
  // Logique identique au traceback global — état H, E, F
FIN TANT QUE

// Les flancs de seq2 (à gauche du point de départ, à droite de end_col) ne sont
// pas inclus dans les chaînes alignées retournées.
```

---

## Fichiers

| Fichier | Rôle |
|---------|------|
| `src/affine_alignment.py` | `fill_matrices_global`, `fill_matrices_semiglobal`, `find_semiglobal_end`, `traceback_global`, `traceback_semiglobal`, `global_affine`, `semiglobal` |
| `tests/test_affine_alignment.py` | 12 tests à faire passer — ne pas modifier |

---

## Critères de validation

Cocher chaque case avant de dire **"TP validé"** :

- [ ] `uv run pytest tests/test_affine_alignment.py -v` passe entièrement
- [ ] `uv run pytest tests/ -v` passe (aucune régression sur TP01–04)
- [ ] Je sais expliquer pourquoi il faut 3 matrices (H, E, F) et non 1
- [ ] Je sais expliquer la différence entre `H[i-1][j] + gap_open + gap_extend` (ouvrir) et `E[i-1][j] + gap_extend` (étendre)
- [ ] Je sais tracer à la main H, E, F pour `"ACT"` / `"ACGT"` (3 matrices de 4×5 cases)
- [ ] Je sais expliquer pourquoi le traceback semi-global s'arrête quand `i == 0` et non quand `i == 0 ET j == 0`
- [ ] Je sais expliquer pourquoi `H[0][j] = 0` en semi-global (et `F[0][j] = 0`)

---

## Pistes si tu bloques

> Essayer dans l'ordre — chaque piste est un peu plus révélatrice que la précédente.

1. **E ou F contiennent des valeurs incohérentes sur le bord :** vérifier case par case. `E[0][j]` doit être `−∞` (seq1 vide, impossible de se trouver en état "gap dans seq2"). `E[i][0]` vaut `H[i][0]` (toutes les bases de seq1 face à des gaps). `F[i][0]` est `−∞` pour la raison symétrique.

2. **Le traceback boucle à l'infini :** la variable `state` n'est probablement pas mise à jour dans la branche E ou F. Dans `state == "E"`, il faut décider du prochain `state` **avant** de décrémenter `i`, puis décrémenter `i` à la fin.

3. **Le traceback produit un mauvais alignement avec de longs gaps :** confusion entre E et F. `state == "E"` = gap dans seq2 = `seq1` avance = on décrémente `i` et on ajoute `seq1[i-1]` + `'-'`. `state == "F"` = gap dans seq1 = `seq2` avance = on décrémente `j` et on ajoute `'-'` + `seq2[j-1]`.

4. **La semi-global donne le même résultat que la global :** vérifier que `fill_matrices_semiglobal` initialise bien `H[0][j] = 0` **et** `F[0][j] = 0` pour tout j. Si seul H est modifié mais pas F, la récurrence F[i][j] reste perturbée par F[0][j] erroné.

5. **Le traceback semi-global inclut les flancs du génome dans les chaînes retournées :** la condition d'arrêt est `i > 0` (pas `i > 0 OU j > 0`). Quand `i == 0`, le read est entièrement aligné — le reste de `j` est un flanc gratuit, à ignorer.

---

## Pour aller plus loin *(optionnel)*

- **Retourner les coordonnées dans seq2 :** `semiglobal` peut aussi renvoyer `(start_col, end_col)` — la position exacte dans le génome. Comment trouver `start_col` à partir du traceback ?
- **Complexité espace :** l'algorithme de Hirschberg réduit l'espace à O(m+n) en utilisant la récursivité divide-and-conquer. Comment adapter ce principe aux 3 matrices ?
- **Pénalités terminales asymétriques :** certains outils appliquent des pénalités de gap réduites (voire nulles) aux extrémités de seq1 aussi, pas seulement seq2. Quelle modification de l'initialisation suffit ?
- **Comparaison avec BioPython :** `Bio.Align.PairwiseAligner` supporte les gaps affins. Comparer tes scores avec les siens sur les mêmes séquences.

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
