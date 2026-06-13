# TP06 — Benchmarking + visualisation

**Projet :** sequence-aligner
**Prérequis :** TP01–05 — tous les modules src/
**Compétences travaillées :** STA, HEU
**Durée estimée :** 3–5 h
**Durée réelle :** *à remplir*
**Fichier(s) à implémenter :** `src/visualize.py` · `src/benchmark.py`

---

## Objectifs

À la fin de ce TP, tu sauras :

- Afficher un alignement en mode texte (format BLAST-like) avec une ligne de correspondance
- Calculer les statistiques d'un alignement : identité, gaps, longueur
- Générer des séquences ADN aléatoires reproductibles pour les tests
- Mesurer le temps d'exécution de nos algorithmes en fonction de la taille des séquences
- Valider nos scores contre BioPython — et identifier une différence de convention

---

## Contexte biologique

### Pourquoi valider contre une référence ?

Implémenter un algorithme d'alignement ne suffit pas : il faut s'assurer que les résultats sont **corrects**. En bioinformatique, une erreur de score silencieuse peut conduire à de fausses conclusions biologiques — par exemple, conclure à tort qu'un gène est conservé entre deux espèces.

La validation croisée contre un outil de référence (ici BioPython, dont le code est audité par la communauté scientifique) est une pratique standard. Si nos scores divergent, il faut comprendre pourquoi : bug chez nous, ou différence de convention de scoring ?

### Pourquoi mesurer les performances ?

Needleman-Wunsch et Smith-Waterman sont O(m × n). Aligner deux gènes de 1000 bases : 10⁶ opérations — rapide. Aligner un génome de 3 Gbp contre lui-même : 9 × 10¹⁸ opérations — impossible. Comprendre la croissance empirique de nos implémentations permet de savoir jusqu'à quelle taille elles restent utilisables, et quand il faut passer à des heuristiques (BLAST).

### BRCA1 humain / souris — un exemple concret

BRCA1 est un gène suppresseur de tumeurs : des mutations germinales dans ce gène multiplient par ~10 le risque de cancer du sein et des ovaires. BRCA1 est conservé entre humain et souris depuis ~80 millions d'années d'évolution séparée.

Un alignement global affin de BRCA1 humain vs murin révèle :
- Des **régions fortement conservées** (domaines BRCT, zinc-finger) → contrainte évolutive forte → fonctionnellement critiques
- Des **régions divergentes** → moins contraintes → introns, régions régulatrices, exons moins essentiels

En pratique :
1. Télécharger BRCA1 humain sur NCBI Nucleotide (accession `NM_007294`)
2. Télécharger BRCA1 murin (`NM_009764`)
3. Les lire avec `read_fasta`, les aligner avec `global_affine` + `ScoringMatrix.nuc_simple()`
4. Visualiser avec `format_alignment` et calculer `alignment_stats`

---

## Problème à résoudre

### Partie 1 — Visualisation (`src/visualize.py`)

**`format_alignment`** : afficher deux chaînes alignées en blocs de `line_width` caractères, avec une ligne centrale montrant les correspondances.

```
Entrée : aligned1="AC-T", aligned2="ACGT", id1="seq1", id2="seq2"

Sortie :
  seq1: AC-T
        || |
  seq2: ACGT
```

Ligne centrale : `|` pour un match (deux bases identiques, aucune `-`), espace sinon.

**`alignment_stats`** : calculer les statistiques de l'alignement.

```
Entrée : aligned1="AC-T", aligned2="ACGT"

Sortie : {"length": 4, "matches": 3, "mismatches": 0, "gaps": 1, "identity": 0.75}
```

### Partie 2 — Benchmark (`src/benchmark.py`)

**`generate_random_dna`** : générer une séquence ADN aléatoire, reproductible avec une graine.

**`time_function`** : mesurer le temps médian d'une fonction sur n répétitions.

**`compare_with_biopython`** : comparer le score de `global_affine` avec celui de `Bio.Align.PairwiseAligner`.

---

## Algorithme

### Partie 1 : Visualisation

```
ALGORITHME format_alignment(aligned1, aligned2, id1, id2, line_width)

  // Construire la ligne de correspondance
  middle ← ""
  POUR i DE 0 À len(aligned1) - 1 :
    SI aligned1[i] != '-' ET aligned2[i] != '-' ET aligned1[i] == aligned2[i] :
      middle ← middle + "|"
    SINON :
      middle ← middle + " "

  // Découper en blocs
  blocs ← []
  POUR start DE 0 À len(aligned1) - 1, PAS DE line_width :
    end   ← min(start + line_width, len(aligned1))
    prefix_seq ← f"{id1}: "
    prefix_mid ← " " × len(prefix_seq)
    bloc ← prefix_seq + aligned1[start:end]  + "\n"
           + prefix_mid + middle[start:end]   + "\n"
           + f"{id2}: " + aligned2[start:end]
    blocs.ajouter(bloc)

  RETOURNER "\n\n".join(blocs)

// Complexité : O(n) où n = len(aligned1)
```

```
ALGORITHME alignment_stats(aligned1, aligned2) → dict

  length ← len(aligned1)
  matches ← 0,  mismatches ← 0,  gaps ← 0

  POUR i DE 0 À length - 1 :
    a ← aligned1[i] ;  b ← aligned2[i]
    SI a == '-' OU b == '-' :
      gaps ← gaps + 1
    SINON SI a == b :
      matches ← matches + 1
    SINON :
      mismatches ← mismatches + 1

  RETOURNER {
    "length"     : length,
    "matches"    : matches,
    "mismatches" : mismatches,
    "gaps"       : gaps,
    "identity"   : matches / length  SI length > 0  SINON  0.0
  }
```

### Partie 2 : Benchmark

```
ALGORITHME generate_random_dna(length, seed=None) → Sequence

  SI seed est fourni : random.seed(seed)
  data ← "".join(random.choice("ACGT") POUR _ DE 0 À length - 1)
  RETOURNER Sequence("random", data)

// Utiliser random.seed() pour garantir la reproductibilité.
```

```
ALGORITHME time_function(func, args=(), kwargs=None, n_repeat=5) → float

  SI kwargs est None : kwargs ← {}
  times ← []
  POUR _ DE 0 À n_repeat - 1 :
    t0 ← time.perf_counter()
    func(*args, **kwargs)
    t1 ← time.perf_counter()
    times.ajouter(t1 - t0)

  RETOURNER statistics.median(times)

// Retourner la médiane (pas la moyenne) : robuste aux pics de charge OS.
```

```
ALGORITHME compare_with_biopython(seq1, seq2, gap_open, gap_extend, match_score, mismatch_score) → dict

  // --- Notre score ---
  our_score, _, _ ← global_affine(seq1, seq2, gap_open, gap_extend, match_score, mismatch_score)

  // --- Score BioPython ---
  // ⚠ DIFFÉRENCE DE CONVENTION :
  //   Notre formule : gap(k) = gap_open + k × gap_extend
  //   BioPython     : gap(k) = open_bio + (k-1) × extend_bio
  //
  // Pour que les deux soient équivalentes, poser :
  //   open_bio   = gap_open + gap_extend
  //   extend_bio = gap_extend
  //
  // Vérification : gap(k) = (gap_open + gap_extend) + (k-1) × gap_extend
  //                       = gap_open + gap_extend + k × gap_extend − gap_extend
  //                       = gap_open + k × gap_extend  ✓

  aligner ← Bio.Align.PairwiseAligner()
  aligner.mode            ← "global"
  aligner.match_score     ← match_score
  aligner.mismatch_score  ← mismatch_score
  aligner.open_gap_score  ← gap_open + gap_extend   // ← ajustement de convention
  aligner.extend_gap_score← gap_extend

  bio_score ← aligner.score(str(seq1), str(seq2))

  RETOURNER {
    "our_score"       : our_score,
    "biopython_score" : bio_score,
    "match"           : |our_score − bio_score| < 1e-6
  }
```

---

## Fichiers

| Fichier | Rôle |
|---------|------|
| `src/visualize.py` | `format_alignment`, `alignment_stats` |
| `src/benchmark.py` | `generate_random_dna`, `time_function`, `compare_with_biopython` |
| `tests/test_visualize.py` | 11 tests — ne pas modifier |
| `tests/test_benchmark.py` | 9 tests — ne pas modifier |

---

## Critères de validation

Cocher chaque case avant de dire **"TP validé"** :

- [ ] `uv run pytest tests/test_visualize.py tests/test_benchmark.py -v` passe entièrement
- [ ] `uv run pytest tests/ -v` passe (aucune régression)
- [ ] Je peux afficher un alignement NW et un alignement SW sur les mêmes séquences et comparer visuellement
- [ ] Je sais expliquer la différence de convention BioPython vs notre formule (et le facteur d'ajustement)
- [ ] J'ai mesuré et tracé (ou imprimé) le temps d'exécution de `global_affine` pour des longueurs 50, 100, 200, 400 — et je sais dire si la croissance observée est cohérente avec O(m×n)
- [ ] J'ai aligné au moins une paire de séquences avec `compare_with_biopython` et vérifié que `result["match"] is True`

---

## Pistes si tu bloques

> Essayer dans l'ordre.

1. **`format_alignment` produit un output mal indenté** : la ligne centrale doit avoir exactement `len(id1) + 2` espaces comme préfixe (pour aligner sous `id1 + ": "`). Si id1 et id2 ont des longueurs différentes, la ligne centrale suit id1.

2. **`alignment_stats` compte mal les gaps** : vérifier que tu testes `aligned1[i] == '-' OR aligned2[i] == '-'` (pas ET). Un gap existe dès qu'un des deux est `-`, peu importe l'autre.

3. **`generate_random_dna` pas reproductible** : `random.seed(seed)` doit être appelé à l'intérieur de la fonction, pas au niveau module. Sinon d'autres appels à `random` entre les deux appels modifient l'état global.

4. **`compare_with_biopython` scores ne correspondent pas** : appliquer la correction de convention `open_bio = gap_open + gap_extend`. Voir pseudo-code ci-dessus. Sans ça, les scores diffèrent systématiquement d'un multiple de `gap_extend`.

5. **`time_function` retourne 0.0** : vérifier que tu appelles `func(*args, **kwargs)` à l'intérieur de la boucle, et que `t0` et `t1` encadrent cet appel.

---

## Pour aller plus loin *(optionnel)*

- **Tracer une courbe de performance** : utiliser `matplotlib` pour tracer le temps en fonction de la longueur de séquence pour NW, SW et global_affine sur la même figure. La courbe doit être quadratique.
- **BRCA1 réel** : télécharger `NM_007294` (humain) et `NM_009764` (souris) depuis NCBI, les lire avec `read_fasta`, les aligner avec `global_affine` + `ScoringMatrix.nuc_simple()`. Quelle est l'identité ?
- **Comparer NW avec BioPython** : BioPython `PairwiseAligner` en mode global avec `gap_score = gap_penalty` (gap linéaire). Vérifier que les scores correspondent à `needleman_wunsch`.
- **Semi-global avec BioPython** : BioPython supporte les pénalités de gap différentes aux bords (`left_gap_score`, `right_gap_score`). Comment configurer le `PairwiseAligner` pour reproduire notre mode semi-global ?

---

## Suivi de temps

*À remplir par l'étudiant au fil du TP.*

| | Durée |
|---|---|
| Estimée | 3–5 h |
| Réelle | |

### Points bloquants rencontrés

| Point bloquant | Temps perdu estimé | Résolution |
|---|---|---|
| | | |
