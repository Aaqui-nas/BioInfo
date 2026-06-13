# TP01 — Représentation des séquences + parsing FASTA

**Projet :** sequence-aligner
**Prérequis :** aucun
**Compétences travaillées :** BIO_BASE, FMT
**Durée estimée :** 3–4 heures
**Durée réelle :** *à remplir*
**Fichiers squelettes :** `skeleton.py`

---

## Objectifs

À la fin de ce TP, tu sauras :

- Représenter une séquence ADN/ARN comme un objet Python typé
- Calculer les opérations biologiques fondamentales sur une séquence (complément, GC content, transcription)
- Lire et écrire des fichiers au format FASTA
- Valider des données en entrée (alphabet biologique)

---

## Contexte biologique

### L'ADN : le support de l'information génétique

L'ADN est formé de deux brins antiparallèles enroulés en double hélice. Chaque brin est une chaîne de **nucléotides**, caractérisés par leur base azotée :

| Base | Lettre | Complément |
|------|--------|-----------|
| Adénine | A | T |
| Thymine | T | A |
| Guanine | G | C |
| Cytosine | C | G |

Les deux brins sont **complémentaires** : chaque A en face d'un T, chaque G en face d'un C. Quand on parle d'une séquence ADN, on donne conventionnellement le brin **sens** lu de 5' vers 3' (sens de la lecture par les polymérases).

Le **complément inverse** (reverse complement) est le brin opposé lu également en 5'→3'. C'est ce qu'on calcule pour trouver le brin antisens, utile pour concevoir des amorces PCR ou aligner un read qui vient de l'autre brin.

### L'ARN et la transcription

Lors de la **transcription**, l'ADN est copié en ARN messager (ARNm) : la Thymine (T) est remplacée par l'Uracile (U). L'ARNm est ensuite traduit en protéine par les ribosomes (traduction).

### Le contenu GC

Le **contenu en GC** (proportion de G et C dans une séquence) est une propriété physique importante : les paires G-C forment 3 liaisons hydrogène contre 2 pour A-T. Une séquence riche en GC a donc une température de fusion (Tm) plus élevée — information cruciale pour concevoir des amorces PCR.

### Le format FASTA

FASTA est le format texte universel pour stocker des séquences biologiques :

```
>identifiant description optionnelle
ATGCATGCATGCATGC
ATGCATGC
```

Chaque séquence commence par une ligne `>` (header), suivie de la séquence sur une ou plusieurs lignes. Un fichier FASTA peut contenir des milliers de séquences.

---

## Problème à résoudre

Créer une bibliothèque Python permettant de représenter, manipuler et stocker des séquences biologiques.

**Exemple :**
```
Input  : fichier FASTA contenant 2 séquences
Output : liste d'objets Sequence avec leurs propriétés calculables

seq = Sequence("brca1", "ATGCATGC")
seq.complement()          → "TACGTACG"
seq.reverse_complement()  → "GCATGCAT"
seq.gc_content()          → 0.5
seq.transcribe()          → Sequence("brca1", "AUGCAUGC", RNA)
```

---

## Algorithmes

### Complément d'une séquence

```
ALGORITHME complement(sequence : str, type : SequenceType) → str
  table ← {A:T, T:A, G:C, C:G}  // pour ADN
           {A:U, U:A, G:C, C:G}  // pour ARN

  résultat ← ""
  POUR chaque base dans sequence FAIRE
    résultat ← résultat + table[base]
  FIN POUR
  RETOURNER résultat
FIN ALGORITHME
// Complexité : O(n) temps, O(n) espace
```

### Complément inverse

```
ALGORITHME reverse_complement(sequence : str) → str
  // Indice : deux opérations que tu as déjà...
  RETOURNER complement(sequence) inversé
FIN ALGORITHME
// Complexité : O(n) temps, O(n) espace
```

### Contenu GC

```
ALGORITHME gc_content(sequence : str) → float
  SI longueur(sequence) = 0 ALORS
    RETOURNER 0.0
  FIN SI

  compteur ← 0
  POUR chaque base dans sequence FAIRE
    SI base ∈ {G, C} ALORS
      compteur ← compteur + 1
    FIN SI
  FIN POUR

  RETOURNER compteur / longueur(sequence)
FIN ALGORITHME
// Complexité : O(n) temps, O(1) espace
```

### Parsing FASTA

```
ALGORITHME read_fasta(filepath : str) → list[Sequence]
  SI fichier n'existe pas ALORS
    LEVER FileNotFoundError
  FIN SI

  séquences ← []
  header_courant ← None
  lignes_seq ← []

  POUR chaque ligne dans le fichier FAIRE
    ligne ← ligne.strip()
    SI ligne est vide ALORS continuer
    SI ligne commence par '>' ALORS
      SI header_courant ≠ None ALORS
        // Sauvegarder la séquence précédente
        séquences.ajouter(Sequence(header_courant, ''.join(lignes_seq)))
      header_courant ← ligne[1:]   // retirer le '>'
      lignes_seq ← []
    SINON
      lignes_seq.ajouter(ligne)
    FIN SI
  FIN POUR

  // Ne pas oublier la dernière séquence du fichier
  SI header_courant ≠ None ALORS
    séquences.ajouter(Sequence(header_courant, ''.join(lignes_seq)))

  RETOURNER séquences
FIN ALGORITHME
// Complexité : O(n) temps, O(n) espace (n = taille totale du fichier)
```

---

## Fichiers squelettes

| Fichier | Rôle |
|---------|------|
| `skeleton.py` | Classes `Sequence`, `SequenceType` + fonctions `read_fasta`, `write_fasta` |
| `test_skeleton.py` | Tests à faire passer — **ne pas modifier** |
| `data/example.fasta` | Fichier FASTA d'exemple avec 4 séquences |

---

## Critères de validation

Cocher chaque case avant de dire **"TP validé"** :

- [ ] `uv run pytest test_skeleton.py -v` passe sans erreur
- [ ] Le code gère les cas limites : fichier vide, séquence multi-lignes, fichier inexistant
- [ ] Je peux expliquer pourquoi le reverse complement se calcule en deux étapes
- [ ] Je connais la complexité temporelle et spatiale de chaque méthode
- [ ] Je peux lire `data/example.fasta` et afficher le GC content de chaque séquence

---

## Pistes si tu bloques

1. La validation de l'alphabet : pense à ce que `set("ATGC")` te donne, et à l'opérateur de différence entre ensembles.
2. Pour le parsing FASTA multi-lignes : le problème classique est de ne pas sauvegarder la *dernière* séquence du fichier. Trace ton algorithme à la main sur un fichier à 2 séquences.
3. Pour `write_fasta` avec `line_width` : `textwrap.wrap(sequence, line_width)` existe en stdlib — mais tu peux aussi faire une boucle sur des tranches.

---

## Pour aller plus loin *(optionnel)*

- Supporter les codes IUPAC ambigus (N = n'importe quelle base, R = A ou G, etc.)
- Ajouter le format FASTQ (séquence + scores de qualité Phred) — utile pour les vrais reads
- Mesurer la complexité de séquence (ratio LZ) pour détecter les répétitions

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
