# TP03 — Matrices de substitution

**Projet :** sequence-aligner
**Prérequis :** TP02 — Alignement global (Needleman-Wunsch)
**Compétences travaillées :** ALI, BIO_BASE
**Durée estimée :** 3–4 h
**Durée réelle :** *à remplir*
**Fichiers à implémenter :** `src/scoring.py` · `src/needleman_wunsch.py` (extension)

---

## Objectifs

À la fin de ce TP, tu sauras :

- Expliquer pourquoi toutes les mutations ne sont pas biologiquement équivalentes
- Distinguer transitions et transversions dans l'ADN
- Implémenter une classe `ScoringMatrix` qui encapsule un système de score par paires
- Parser un fichier de matrice au format standard (BLOSUM/PAM)
- Intégrer `ScoringMatrix` dans Needleman-Wunsch sans casser les tests TP02

---

## Contexte biologique

### Toutes les substitutions ne se valent pas

Dans le modèle simple du TP02, un mismatch coûte toujours `-1`, que ce soit `A→G` ou `A→C`. Mais biologiquement, ces deux substitutions n'ont pas la même probabilité d'apparition.

**Pour l'ADN — transitions vs transversions :**

Les nucléotides se divisent en deux classes chimiques :
- **Purines** : A, G (double cycle à deux anneaux)
- **Pyrimidines** : C, T (cycle à un anneau)

Une **transition** remplace un nucléotide par un autre *de la même classe* :
- A ↔ G (purine ↔ purine)
- C ↔ T (pyrimidine ↔ pyrimidine)

Une **transversion** remplace un nucléotide par un de *l'autre classe* :
- A ↔ C, A ↔ T, G ↔ C, G ↔ T

Les transitions sont environ **2× plus fréquentes** que les transversions dans les génomes, car la structure chimique à remplacer est plus similaire. Pénaliser uniformément les deux types de substitution introduit donc un biais : l'alignement va mal interpréter les régions riches en transitions.

**Pour les protéines — matrices BLOSUM :**

BLOSUM62 est l'exemple canonique de matrice de substitution en bioinformatique. Elle a été construite empiriquement à partir de blocs de séquences protéiques conservées (base de données BLOCKS) :

1. Aligner des familles de protéines homologues (> 62 % d'identité → BLOSUM62)
2. Compter la fréquence observée de chaque paire d'acides aminés en position alignée
3. La comparer à la fréquence attendue si les substitutions étaient aléatoires
4. Score = log(observé / attendu)

Résultat : une substitution conservative (Leu → Ile, tous deux hydrophobes) a un score positif, une substitution radicale (Trp → Pro) a un score très négatif.

**Pourquoi ce TP utilise-t-il une matrice ADN et non BLOSUM62 ?**

BLOSUM62 est définie sur les 20 acides aminés, pas sur les nucléotides. Pour ce projet d'alignement ADN, on va implémenter `NUC_SIMPLE`, une matrice 4×4 {A, T, G, C} qui capture la distinction transition/transversion — l'équivalent conceptuel de BLOSUM pour l'ADN.

---

## Problème à résoudre

**Entrée :** deux séquences ADN, un `gap_penalty`, et une `ScoringMatrix` optionnelle.

**Sortie :** score optimal + alignement (comme TP02), mais le score diagonal est maintenant déterminé par la matrice de substitution plutôt que par match/mismatch.

**Exemple illustrant l'impact :**
```
seq1 = "AT",  seq2 = "CA",  gap_penalty = -2

Avec scoring simple (match=+1, mismatch=-1) :
  alignement : "AT" vs "CA"    score = -1 + -1 = -2

Avec NUC_SIMPLE (match=+1, transition=-1, transversion=-3) :
  alignement : "-AT" vs "CA-"  score = -2 + 1 + -2 = -3
  (deux gaps coûtent -4, mais on récupère +1 sur le A-A match = -3)
  (deux transversions coûteraient -3 + -3 = -6 → les gaps sont préférés)
```

La matrice change l'alignement lui-même, pas seulement le score.

---

## Algorithme

### Partie 1 — Classe `ScoringMatrix`

```
CLASSE ScoringMatrix

  __init__(data : dict[(str, str) → int]) :
    // Stocker data comme attribut interne
    // Pas de validation requise — faire confiance à l'appelant

  score(a : str, b : str) → int :
    // Retourner data[(a, b)]
    // Laisser lever KeyError si la paire est absente

  @classmethod
  from_file(filepath : str) → ScoringMatrix :
    // Lever FileNotFoundError si le fichier n'existe pas (open() le fait automatiquement)
    // Lire le fichier ligne par ligne :
    //   - ignorer les lignes vides et les lignes commençant par '#'
    //   - la PREMIÈRE ligne utile = noms des colonnes (ex: ["A", "T", "G", "C"])
    //   - chaque ligne SUIVANTE : tokens = ligne.split()
    //       tokens[0] = nom de la ligne (ex: "A")
    //       tokens[1:] = scores correspondant à chaque colonne
    // Construire le dict {(ligne, colonne): score} pour toutes les combinaisons
    // Retourner ScoringMatrix(dict)
    // Complexité : O(n²) avec n = taille de l'alphabet

  @classmethod
  identity_dna() → ScoringMatrix :
    // Construire et retourner la matrice identité pour {A, T, G, C}
    // match = +1, mismatch = -1

  @classmethod
  nuc_simple() → ScoringMatrix :
    // Construire et retourner NUC_SIMPLE pour {A, T, G, C}
    // match = +1
    // transition (A↔G, C↔T) = -1
    // transversion (A↔C, A↔T, G↔C, G↔T) = -3
    // Conseil : lister explicitement toutes les 16 paires dans un dict plutôt que
    //   de calculer — c'est plus lisible et moins sujet aux erreurs de logique
FIN CLASSE
```

### Partie 2 — Extension de `needleman_wunsch.py`

Les stubs TP03 sont déjà insérés dans `fill_matrix` et `traceback`. Il suffit de remplacer les deux `raise NotImplementedError` par la logique suivante :

```
// Dans fill_matrix — boucle interne :
SI scoring_matrix n'est pas None :
  diag ← F[i-1][j-1] + scoring_matrix.score(seq1[i-1], seq2[j-1])
SINON :
  // comportement TP02 inchangé
  SI seq1[i-1] == seq2[j-1] : diag ← F[i-1][j-1] + match_score
  SINON                      : diag ← F[i-1][j-1] + mismatch_score

// Dans traceback — même logique pour calculer diag_score :
SI scoring_matrix n'est pas None :
  diag_score ← scoring_matrix.score(seq1[i-1], seq2[j-1])
SINON :
  diag_score ← match_score SI seq1[i-1]==seq2[j-1], mismatch_score SINON
```

**Attention :** le `scoring_matrix` passé à `fill_matrix` et à `traceback` doit être le *même objet* — les deux fonctions doivent utiliser exactement la même logique de score, sinon le traceback peut ne pas retrouver le bon chemin dans la matrice.

---

## Fichiers

| Fichier | Rôle |
|---------|------|
| `src/scoring.py` | `ScoringMatrix` + méthodes `identity_dna`, `nuc_simple`, `from_file` à implémenter |
| `src/needleman_wunsch.py` | Deux stubs TP03 à remplacer dans `fill_matrix` et `traceback` |
| `tests/test_scoring.py` | Tests automatisés — ne pas modifier |
| `data/nuc_simple.mat` | Fichier de matrice fourni pour `from_file` |

---

## Critères de validation

Cocher chaque case avant de dire **"TP validé"** :

- [ ] `uv run pytest tests/test_scoring.py -v` passe entièrement
- [ ] `uv run pytest tests/test_needleman_wunsch.py -v` passe encore (pas de régression)
- [ ] Je peux expliquer la différence entre une transition et une transversion
- [ ] Je comprends pourquoi dans l'exemple "AT"/"CA", NUC_SIMPLE préfère deux gaps plutôt que deux transversions
- [ ] Je sais lire un fichier de matrice au format BLOSUM (en-tête + lignes de scores)

---

## Pistes si tu bloques

> Essayer dans l'ordre.

1. **`from_file` : la première ligne utile devient l'en-tête** — lis le fichier deux fois mentalement : une fois pour identifier quelle ligne est l'en-tête, une fois pour parser les données. Une variable booléenne `header_seen` peut aider.

2. **`nuc_simple` : tu hésites sur qui est transition et qui est transversion** — rappelle-toi A et G sont des purines (les "grands"), C et T sont des pyrimidines (les "petits"). Transition = même famille. Transversion = changement de famille.

3. **Le test d'intégration échoue même après avoir implémenté `ScoringMatrix`** — vérifie que `needleman_wunsch` propage bien `scoring_matrix` aux deux appels internes (`fill_matrix` et `traceback`). Le stub dans `needleman_wunsch` fait déjà ça, mais regarde si tu n'as pas accidentellement supprimé le paramètre.

---

## Pour aller plus loin *(optionnel)*

- La matrice **NUC44** de NCBI est l'équivalent "officiel" de BLOSUM pour l'ADN — elle inclut aussi les codes IUPAC d'ambiguïté (N, R, Y…). Comment l'adapter à notre `ScoringMatrix` ?
- Comment construire empiriquement une matrice de substitution depuis un alignement multiple existant (reproduire la démarche BLOSUM) ?
- Le paramètre `gap_penalty` est ici un gap *linéaire* (chaque position de gap coûte la même pénalité). Le TP05 introduira les **gap penalties affines** (ouverture + extension). En quoi est-ce biologiquement plus réaliste ?

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
