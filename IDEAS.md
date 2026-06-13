# Banque d'idées de projets

Chaque idée est accompagnée d'une orientation (algo / deep learning / pipeline), d'un niveau de difficulté estimé, et des technos suggérées.

---

## Idées disponibles

### 1. Aligneur de séquences ADN
**Orientation :** algorithmique
**Difficulté :** ⭐⭐ (bonne entrée en matière)
**Technos suggérées :** Python ou C++
**Description :** Implémenter les algorithmes fondamentaux d'alignement de séquences — alignement global (Needleman-Wunsch), local (Smith-Waterman), puis semi-global et heuristique (BLAST-like). Progresser vers l'alignement multiple.
**Compétences visées :** programmation dynamique, matrices de substitution (BLOSUM/PAM), complexité algorithmique.

---

### 2. Recherche de motifs dans l'ADN
**Orientation :** algorithmique
**Difficulté :** ⭐⭐
**Technos suggérées :** C++ (performance) ou Python
**Description :** Implémenter la recherche de motifs biologiques (sites de restriction, promoteurs, répétitions) avec des structures de données avancées : suffix array, suffix tree, Aho-Corasick, puis évaluation statistique des motifs (p-value, score PWM).
**Compétences visées :** structures de données sur chaînes, statistiques sur séquences.

---

### 3. Assemblage de génomes (De Bruijn)
**Orientation :** algorithmique / graphes
**Difficulté :** ⭐⭐⭐
**Technos suggérées :** Python ou C++
**Description :** Construire un assembleur de novo simplifié à partir de reads simulés. Passer du graphe overlap-layout-consensus aux graphes de De Bruijn, gérer les erreurs de séquençage et les répétitions.
**Compétences visées :** théorie des graphes, algorithmes d'assemblage, gestion des erreurs.

---

### 4. Classificateur de séquences par deep learning
**Orientation :** deep learning
**Difficulté :** ⭐⭐⭐
**Technos suggérées :** Python (PyTorch)
**Description :** Construire un pipeline de classification de séquences ADN/ARN — encodage one-hot ou k-mer, CNN 1D, puis RNN/LSTM, et enfin Transformer. Appliquer à un problème réel (détection de promoteurs, classification taxonomique, prédiction de sites d'épissage).
**Compétences visées :** deep learning appliqué aux séquences, encodage biologique, évaluation de modèles.

---

### 5. Prédiction de structure secondaire d'ARN
**Orientation :** algorithmique + potentiellement deep learning
**Difficulté :** ⭐⭐⭐
**Technos suggérées :** Python ou C++
**Description :** Implémenter l'algorithme de Nussinov (maximisation des appariements), puis l'algorithme de Zuker (énergie minimale libre), et visualiser les structures en dot-bracket. Extension possible : comparaison avec des prédictions par réseau de neurones.
**Compétences visées :** programmation dynamique avancée, thermodynamique des acides nucléiques.

---

### 6. Mini pipeline de variant calling
**Orientation :** pipeline / algorithmique
**Difficulté :** ⭐⭐⭐⭐
**Technos suggérées :** Python + outils standards (BWA, samtools — en les appelant)
**Description :** Construire un mini-pipeline de détection de variants (SNPs, indels) à partir de reads simulés : alignement sur référence, tri/indexation BAM, pileup, appel de variants, filtrage et annotation. Introduction aux formats bioinformatiques (FASTQ, SAM/BAM, VCF).
**Compétences visées :** formats standards, pipeline bioinformatique, statistiques de variants.

---

### 7. Modèle de langage génomique (DNA-LM)
**Orientation :** deep learning avancé
**Difficulté :** ⭐⭐⭐⭐⭐
**Technos suggérées :** Python (PyTorch)
**Description :** Reproduire en miniature l'approche DNABERT/Nucleotide Transformer — tokenisation k-mer, pré-entraînement BERT sur séquences génomiques, fine-tuning sur tâche de classification. Introduction aux LLM appliqués à la génomique.
**Compétences visées :** architecture Transformer, self-supervised learning, fine-tuning, bioinformatique moderne.

---

### 8. Phylogénomique computationnelle
**Orientation :** algorithmique / statistiques
**Difficulté :** ⭐⭐⭐
**Technos suggérées :** Python
**Description :** Calculer des matrices de distance entre séquences, implémenter les algorithmes de construction d'arbres (UPGMA, Neighbor-Joining), visualiser les dendrogrammes phylogénétiques et évaluer leur robustesse (bootstrap).
**Compétences visées :** évolution moléculaire, clustering hiérarchique, visualisation.

---

### 9. Compression de séquences génomiques
**Orientation :** algorithmique / systèmes
**Difficulté :** ⭐⭐⭐
**Technos suggérées :** C++ ou Rust
**Description :** Explorer les algorithmes de compression spécialisés pour l'ADN — encodage 2 bits, run-length encoding, LZ sur bases, puis compression par référence (stockage de différences). Comparer aux compresseurs génériques.
**Compétences visées :** algorithmique bas niveau, performance, structures de données compactes.

---

### 10. Détection de gènes par HMM
**Orientation :** probabiliste / algorithmique
**Difficulté :** ⭐⭐⭐
**Technos suggérées :** Python
**Description :** Implémenter un Hidden Markov Model pour annoter les états génomiques (exon, intron, intergénique) — algorithmes Forward-Backward, Viterbi, Baum-Welch. Application à la prédiction de gènes simplifiée.
**Compétences visées :** modèles probabilistes, inférence statistique, annotation génomique.

---

### 11. BioSeqKit — Boîte à outils généraliste pour séquences ADN/ARN
**Orientation :** algorithmique + ingénierie logicielle (projet transversal)
**Difficulté :** ⭐⭐ à ⭐⭐⭐⭐ (progressive sur toute la durée)
**Technos suggérées :** Python (cœur) + C++ ou Rust (extensions performance)
**Description :** Construire de zéro une bibliothèque réutilisable et une CLI pour le traitement de séquences biologiques. Le projet est volontairement transversal : il agrège des modules issus des autres idées (parsers, alignement, motifs, stats, k-mers…) mais dans une optique de conception logicielle propre — API cohérente, modules découplés, tests, documentation.

Progression des modules :
1. **Core** — types de séquences (ADN, ARN, protéine), opérations fondamentales (complément, reverse complement, transcription, traduction), parsers FASTA/FASTQ
2. **Stats** — composition en bases, GC content, entropie, distributions de k-mers, complexité de séquence (LZ)
3. **Search** — recherche exacte de motifs, IUPAC ambiguity codes, scoring de matrices PWM
4. **Align** — alignement pairwise global/local (Needleman-Wunsch / Smith-Waterman), scoring configurable
5. **QC** — filtrage qualité de reads, trimming d'adaptateurs, détection de duplicats
6. **CLI** — interface en ligne de commande unifiée exposant tous les modules avec entrée/sortie fichiers standards
7. **Extension C++/Rust** *(optionnel)* — réimplémenter un module critique (k-mer counting, alignement) en natif et le brancher via FFI/pybind

**Ce qui rend ce projet spécial :**
- Il force à concevoir une vraie API publique (choix structurants dès le TP01)
- Les TPs suivants réutilisent obligatoirement les TPs précédents — la dette technique est immédiatement visible
- Peut servir de base aux autres projets (le classifier DL peut importer le parser FASTA de BioSeqKit)
- Introduit les bonnes pratiques d'ingénierie logicielle : typage, tests unitaires, versioning sémantique

**Compétences visées :** architecture de bibliothèque, API design, formats bio standards, algo sur séquences, interop Python/C++ ou Python/Rust, ingénierie logicielle.

---

### 12. Aligneur de séquences avancé *(suite de sequence-aligner)*
**Orientation :** algorithmique
**Difficulté :** ⭐⭐⭐
**Technos suggérées :** Python (+ C++ optionnel pour les parties performance)
**Description :** Approfondir le projet `sequence-aligner` en implémentant tout ce qui était marqué "pour aller plus loin". Ce projet suppose que NW, SW et Gotoh sont déjà maîtrisés — il travaille sur la performance, la complétude et l'heuristique.

Progression des modules :
1. **Formats étendus** — FASTQ avec scores Phred, codes IUPAC ambigus (N, R, Y, S…), matrice NUC44 officielle NCBI
2. **Complexité de séquence** — ratio LZ pour détecter les régions répétées (low-complexity masking)
3. **Optimisation espace** — algorithme de Hirschberg (O(min(m,n)) espace, même complexité temporelle que NW) ; adaptation aux 3 matrices affines (Gotoh + Hirschberg)
4. **Distance d'édition** — variante Levenshtein de la matrice NW (min + signes inversés), comparaison avec NW classique
5. **SW banded** — Smith-Waterman en bande (O(k×n) espace) quand on connaît à l'avance la divergence maximale k
6. **Tous les alignements locaux** — retourner tous les alignements locaux de score ≥ seuil, pas seulement le meilleur
7. **Heuristique BLAST-like** — seed-and-extend : trouver les k-mers exacts communs, étendre localement, filtrer par score — première approche sub-quadratique
8. **Construction de matrice de substitution** — reproduire la démarche BLOSUM depuis un alignement multiple existant (compter les paires observées vs attendues)

**Compétences visées :** optimisation mémoire (Hirschberg), distance d'édition, heuristiques d'alignement (BLAST), formats bio étendus (FASTQ, IUPAC), construction de matrices de substitution.

---

## Idées en vrac (à développer)

- Clustering de métagénomes (k-mer frequency + UMAP/t-SNE)
- Simulation de séquençage Oxford Nanopore (signal brut → base calling simplifié)
- Détection de CRISPR-Cas9 off-targets
- Analyse de ChIP-seq (pics de liaison protéine-ADN)
