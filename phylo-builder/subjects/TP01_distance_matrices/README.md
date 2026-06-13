# TP01 — Matrices de distance évolutive

**Projet :** phylo-builder
**Prérequis :** aucun
**Compétences travaillées :** STA, BIO_BASE
**Durée estimée :** 3–4 heures
**Durée réelle :** *à remplir*
**Fichier(s) à implémenter :** `src/distances.py`

---

## Objectifs

À la fin de ce TP, tu sauras :

- Calculer et interpréter trois mesures de distance évolutive entre séquences alignées
- Comprendre pourquoi la distance brute (Hamming) sous-estime la distance réelle, et comment les modèles de substitution corrigent cela
- Encapsuler une matrice de distances pairwise dans une structure de données utilisable par les algorithmes de construction d'arbres

---

## Contexte biologique

### Pourquoi mesurer des distances entre séquences ?

Avant de construire un arbre phylogénétique, il faut quantifier à quel point deux espèces sont "éloignées" l'une de l'autre. La première idée naïve : compter le pourcentage de positions différentes entre deux séquences alignées. C'est la **distance de Hamming** (ou p-distance).

Le problème : cette mesure **sous-estime systématiquement** la distance évolutive réelle. Pourquoi ?

Imagine deux séquences qui divergent depuis très longtemps. À la position 5, l'ancêtre commun avait un `A`. Depuis, la lignée 1 a muté `A→G`, et la lignée 2 a muté `A→T` puis `T→A`. Résultat : les deux séquences affichent un `A` à la position 5 — elles semblent identiques ! Pourtant, deux mutations ont eu lieu. Ce phénomène s'appelle les **substitutions multiples** (ou hit-back mutations).

### Les modèles de substitution

Pour corriger cette sous-estimation, on utilise des **modèles probabilistes d'évolution** qui supposent que les substitutions se produisent selon un processus de Poisson.

**Jukes-Cantor (1969)** — le modèle le plus simple : toutes les substitutions sont également probables (taux α identique entre tous les nucléotides). Si on observe une proportion p de sites différents, la distance corrigée est :

```
d_JC = -3/4 × ln(1 - 4/3 × p)
```

Pour les faibles valeurs de p, d_JC ≈ p (la correction est négligeable). Plus p est proche de 3/4, plus la correction est grande — au-delà de 3/4, le modèle est saturé (trop de substitutions multiples pour être corrigées).

**Kimura 2 paramètres / K2P (1980)** — modèle plus réaliste : les substitutions ne sont pas toutes également probables. En biologie, les **transitions** (A↔G ou C↔T — changement au sein d'une même classe chimique) sont 2 à 4 fois plus fréquentes que les **transversions** (A↔C, A↔T, G↔C, G↔T — changement entre classes). Le K2P utilise deux paramètres :
- P = proportion de sites avec une transition
- Q = proportion de sites avec une transversion

```
d_K2P = -1/2 × ln(1 - 2P - Q) - 1/4 × ln(1 - 2Q)
```

Cas limites : si (1 - 2P - Q) ≤ 0 ou (1 - 2Q) ≤ 0, le modèle est saturé → distance infinie.

### Pourquoi K2P > JC en pratique ?

Pour le même nombre de différences observées, K2P et JC donnent des distances différentes car K2P sait que les transitions sont plus probables — donc une transition observée implique plus de distance que ce que JC estimerait.

---

## Problème à résoudre

Implémenter trois fonctions de distance sur des séquences ADN **pré-alignées** (même longueur), puis une classe `DistanceMatrix` qui calcule et stocke toutes les paires.

**Exemple :**
```
Input  : seq1 = "AATT"
         seq2 = "ACTT"

Hamming : 1 position différente sur 4 → p = 0.25
JC      : -3/4 × ln(1 - 4/3 × 0.25) = -3/4 × ln(2/3) ≈ 0.3045
K2P     : A→C est une transversion → P=0, Q=0.25
          -1/2 × ln(0.75) - 1/4 × ln(0.5) ≈ 0.3171
```

---

## Algorithme

### hamming_distance

```
ALGORITHME hamming_distance(seq1 : str, seq2 : str) → float

  SI len(seq1) ≠ len(seq2) ALORS
    LEVER ValueError

  differences ← 0
  comparables ← 0

  POUR i DE 0 À len(seq1) - 1
    SI seq1[i] = '-' OU seq2[i] = '-' ALORS
      CONTINUER  // position de gap : exclue
    comparables ← comparables + 1
    SI seq1[i] ≠ seq2[i] ALORS
      differences ← differences + 1

  SI comparables = 0 ALORS
    RETOURNER 0.0

  RETOURNER differences / comparables

// Complexité : O(n) temps, O(1) espace
```

### jukes_cantor_distance

```
ALGORITHME jukes_cantor_distance(seq1 : str, seq2 : str) → float

  p ← hamming_distance(seq1, seq2)

  SI p = 0.0 ALORS RETOURNER 0.0
  SI p >= 3/4 ALORS RETOURNER +inf  // saturation

  RETOURNER -3/4 × ln(1 - 4/3 × p)

// Complexité : O(n) temps, O(1) espace
```

### kimura_2p_distance

```
ALGORITHME kimura_2p_distance(seq1 : str, seq2 : str) → float

  // Identifier si chaque substitution est transition ou transversion
  TRANSITIONS ← { (A,G), (G,A), (C,T), (T,C) }
  TRANSVERSIONS ← { (A,C), (C,A), (A,T), (T,A), (G,C), (C,G), (G,T), (T,G) }

  ts ← 0  // compteur de transitions
  tv ← 0  // compteur de transversions
  n  ← 0  // sites comparables (hors gaps)

  POUR chaque position i
    SI gap dans seq1[i] ou seq2[i] ALORS CONTINUER
    n ← n + 1
    paire ← (seq1[i], seq2[i])
    SI paire ∈ TRANSITIONS ALORS ts ← ts + 1
    SI paire ∈ TRANSVERSIONS ALORS tv ← tv + 1

  SI n = 0 ALORS RETOURNER 0.0

  P ← ts / n
  Q ← tv / n

  SI P = 0 ET Q = 0 ALORS RETOURNER 0.0

  a1 ← 1 - 2×P - Q
  a2 ← 1 - 2×Q

  SI a1 <= 0 OU a2 <= 0 ALORS RETOURNER +inf

  RETOURNER -1/2 × ln(a1) - 1/4 × ln(a2)

// Complexité : O(n) temps, O(1) espace
```

### DistanceMatrix.from_sequences

```
ALGORITHME from_sequences(sequences : dict[str, str], method : str) → DistanceMatrix

  noms ← liste des clés de sequences (ordre fixe)
  n    ← len(noms)
  distance_fn ← sélectionner la fonction selon method
                (lever ValueError si method inconnu)

  matrix ← tableau n×n initialisé à 0.0

  POUR i DE 0 À n-1
    POUR j DE i+1 À n-1
      d ← distance_fn(sequences[noms[i]], sequences[noms[j]])
      matrix[i][j] ← d
      matrix[j][i] ← d  // symétrie : d(A,B) = d(B,A)
    // matrix[i][i] reste 0.0 (distance à soi-même)

  RETOURNER DistanceMatrix(noms, matrix)

// Complexité : O(n² × L) temps, O(n²) espace  (n séquences, L longueur)
```

---

## Fichiers

| Fichier | Rôle |
|---------|------|
| `src/distances.py` | `hamming_distance`, `jukes_cantor_distance`, `kimura_2p_distance`, `DistanceMatrix` |
| `tests/test_distances.py` | suite de tests fournie — ne pas modifier |
| `data/example.fasta` | 5 séquences alignées (20 bases) pour exploration manuelle |

---

## Critères de validation

Cocher chaque case avant de dire **"TP validé"** :

- [ ] `uv run pytest tests/test_distances.py -v` passe sans erreur
- [ ] `hamming_distance("AATT", "ACTT")` retourne exactement `0.25`
- [ ] `jukes_cantor_distance("AAAA", "TTTT")` retourne `math.inf`
- [ ] `kimura_2p_distance("AATT", "AGTT")` > `kimura_2p_distance("AATT", "ACTT")` (même p, ts > tv)
- [ ] `DistanceMatrix.from_sequences(...)` produit une matrice symétrique à diagonale nulle
- [ ] Je peux expliquer pourquoi JC sous-estime moins que Hamming mais reste insuffisant pour certains cas
- [ ] Je connais la différence biologique entre transition et transversion

---

## Pistes si tu bloques

> Essayer dans l'ordre — chaque piste est un peu plus révélatrice.

1. **Hamming — les gaps** : un caractère `'-'` dans une séquence alignée représente une insertion ou délétion. Ne pas le compter comme une différence, ni même comme un site comparable. Quel type Python permet de tester facilement si un caractère est dans un ensemble donné ?

2. **JC — la condition de saturation** : si p = 0.75, que vaut `1 - 4/3 × 0.75` ? Que se passe-t-il si on passe cette valeur à `math.log` ? Quelle garde mettre avant l'appel ?

3. **K2P — classifier les paires** : avant de coder la formule, implémente d'abord une fonction auxiliaire `_classify(a, b)` qui retourne `"ts"`, `"tv"` ou `"same"`. Ça rend le corps de `kimura_2p_distance` plus lisible et testable séparément.

---

## Pour aller plus loin *(optionnel, hors critères de validation)*

- **Distance de Tamura-Nei** — généralise K2P avec des fréquences de bases non uniformes. Qu'ajoute ce troisième paramètre par rapport à K2P ?
- **Distance de LogDet** — indépendante de la composition en bases : utile quand les espèces comparées ont des GC-content très différents.
- **Visualisation de la matrice** : afficher `DistanceMatrix` comme une carte de chaleur avec `matplotlib.pyplot.imshow` — quelle propriété de la matrice garantit qu'elle sera symétrique ?

---

## Suivi de temps

*À remplir par l'étudiant au fil du TP.*

| | Durée |
|---|---|
| Estimée | 3–4 h |
| Réelle | |

### Points bloquants rencontrés

| Point bloquant | Temps perdu estimé | Résolution |
|---|---|---|
| | | |
