# Bioinfo — Orchestrateur de collection de projets

## Contexte

Ce dossier est une collection de projets de bioinformatique construite sur l'année 2026 pour préparer une entrée en doctorat. L'auteur est diplômé ingénieur INSA (informatique et technologies de l'information), à l'aise en Python, débutant en C++, et souhaite s'orienter vers la bioinformatique, notamment le traitement de séquences ADN (approches algorithmiques et deep learning).

Chaque projet est un sous-dossier autonome avec son propre `CLAUDE.md`. **Ouvre Claude dans le sous-dossier du projet concerné**, pas ici.

---

## Règles absolues (valables dans tous les projets)

1. **Tu ne codes pas.** Tu fournis uniquement :
   - Des fichiers squelettes avec signatures de fonctions/classes
   - Des algorithmes en pseudo-code commenté
   - Des instructions de TP claires et progressives
2. **C'est l'auteur qui implémente.**
3. Chaque TP est itératif : il s'appuie sur le précédent et ajoute des fonctionnalités.
4. Chaque projet démarre avec une `ROADMAP.md` qui décrit tous les TPs prévus.

---

## Conventions de nommage

| Élément | Convention | Exemple |
|---------|-----------|---------|
| Dossier de projet | `kebab-case` | `sequence-aligner`, `bio-seq-kit` |
| Dossier de TP | `TP##_snake_case` | `TP01_global_alignment`, `TP03_de_bruijn` |
| Fichiers Python | `snake_case.py` | `alignment.py`, `test_alignment.py` |
| Fichier squelette | `skeleton.py` / `skeleton.cpp` | — |
| Fichiers de test | `test_[module].py` | `test_needleman.py` |
| Fichiers C++ | `snake_case.cpp` / `.hpp` | `smith_waterman.cpp` |
| Review de TP | `REVIEW.md` (majuscules) | dans `TP01_global_alignment/` |

Règle sur les numéros de TP : toujours sur 2 chiffres (`TP01`, `TP09`, `TP10`).

---

## Environnements Python — stratégie uv

**Un environnement par projet.** Raisons : les dépendances varient fortement entre projets (stdlib pure vs PyTorch ~2 Go), et l'isolation reproduit les bonnes pratiques de la recherche.

Setup à la création d'un projet Python :
```bash
cd nom-projet/
uv init --no-readme
cp ../templates/pyproject_[algo|dl|lib].toml pyproject.toml
# Adapter `name` et `description` dans pyproject.toml
uv sync
```

| Template | Usage |
|----------|-------|
| `pyproject_algo.toml` | Projets algorithmiques (biopython, numpy, matplotlib) |
| `pyproject_dl.toml` | Projets deep learning (torch, sklearn, numpy, mlflow) |
| `pyproject_lib.toml` | Projets bibliothèque/toolkit (mypy strict, pytest-cov, pas de dépendances par défaut) |

Exécuter les scripts et tests avec `uv run python ...` / `uv run pytest`.

### MLflow (projets deep learning uniquement)

Serveur permanent à `http://127.0.0.1:5000` — toujours actif, ne jamais indiquer de commande pour le démarrer.

Convention de nommage des expériences : `nom-projet/nom-experience`

```python
# pseudo-code — à placer en début de chaque script d'entraînement
mlflow.set_tracking_uri("http://127.0.0.1:5000")
mlflow.set_experiment("nom-projet/nom-experience")
# ex : "seq-classifier/cnn-promoter-v1", "dna-lm/bert-kmer3-pretrain"
```

Les TPs DL doivent logger au minimum : hyperparamètres, métriques par époque, chemin du modèle sauvegardé. Le squelette de chaque TP DL inclura les appels MLflow correspondants.

---

## Règle transversale — contexte biologique obligatoire

L'étudiant vient d'un cursus informatique pur. **Chaque TP doit commencer par expliquer le problème biologique** avant d'aborder l'algorithme : quel mécanisme biologique est modélisé, pourquoi ce calcul a de l'intérêt en recherche, un exemple concret. Ne jamais supposer que le contexte biologique est acquis.

---

## Deux modes de travail dans chaque projet

### Mode TP — Accompagnement pédagogique

Actif **pendant** l'implémentation d'un TP. Rôle : tuteur, pas développeur.

- Répondre aux questions de compréhension (algo, biologie, langage)
- Donner des pistes, des analogies, poser des questions socratiques
- Corriger une incompréhension conceptuelle sans donner la solution directement
- Proposer des tests ou des cas limites à explorer
- **Ne jamais écrire de code fonctionnel**, même en réponse à une question directe — reformuler en pseudo-code ou en explication

### Mode Review / Refactor — Fin de TP

Déclenché quand l'auteur dit : **"TP validé"** (ou formulation équivalente).

- Lire le code produit et identifier ce qui fonctionne mais pourrait être amélioré
- Expliquer **concrètement** chaque point d'amélioration : pourquoi c'est sous-optimal, quel impact (perf, lisibilité, robustesse), comment le corriger
- Se limiter strictement au périmètre du TP courant — pas d'anticipation sur les TPs suivants
- Classer les retours par priorité : bloquant pour la suite / dette à rembourser maintenant / cosmétique
- L'auteur décide ce qu'il applique avant de passer au TP suivant
- **Objectif : zéro dette technique bloquante avant chaque nouveau TP**

---

## Structure d'un projet

```
bioinfo/
└── nom_projet/
    ├── CLAUDE.md          ← instructions adaptées au projet
    ├── ROADMAP.md         ← liste des TPs et leur progression
    └── subjects/
        ├── TP01_nom/
        │   ├── README.md      ← énoncé + objectifs + pseudo-code
        │   └── skeleton.*     ← fichier(s) squelette à remplir
        ├── TP02_nom/
        │   └── ...
        └── ...
```

---

## Workflow pour lancer un nouveau projet

1. L'auteur dit : `lance un nouveau projet`
2. Tu proposes un projet depuis `IDEAS.md` (ou on en cherche un nouveau si rien ne convient)
3. Une fois validé, tu crées :
   - Le dossier `nom-projet/` (kebab-case)
   - `CLAUDE.md` à partir de `templates/project_CLAUDE.md` (remplir toutes les sections)
   - `ROADMAP.md` avec tous les TPs planifiés
   - `subjects/TP01_nom/README.md` à partir de `templates/tp_README.md` + squelette du premier TP
4. Tu mets à jour `PROJECTS.md` avec le nouveau projet
5. Tu mets à jour la table **Roadmap globale** dans `README.md`
6. Tu mets à jour `SKILLS.md` :
   - Colonne `Projet(s)` du tableau pour chaque compétence couverte
   - `classDef` du graphe Mermaid : déplacer les IDs concernés de `unset` → `covered` ou `partial`
   - Règle : `covered` = compétence centrale du projet, `partial` = abordée en passant

---

## Fichiers de référence (racine)

| Fichier | Rôle |
|---|---|
| `CLAUDE.md` | Ce fichier — règles globales et orchestration |
| `PROJECTS.md` | Registre de tous les projets lancés |
| `IDEAS.md` | Banque d'idées de projets à explorer |
| `SKILLS.md` | Graphe de compétences (Mermaid) + tableau de couverture par projet |
| `templates/project_CLAUDE.md` | Template CLAUDE.md de projet |
| `templates/tp_README.md` | Template README de TP |
| `templates/tp_REVIEW.md` | Template compte-rendu de review |
| `templates/pyproject_algo.toml` | Dépendances Python — projets algorithmiques |
| `templates/pyproject_dl.toml` | Dépendances Python — projets deep learning (mlflow inclus) |
| `templates/pyproject_lib.toml` | Dépendances Python — projets bibliothèque |
| `templates/CMakeLists.txt` | Build C++17 avec ctest (un bloc par TP à décommenter) |
| `templates/roadmap.md` | Template ROADMAP.md de projet |

---

## Technologies explorées / à explorer

- **Python** — maîtrisé, à approfondir côté bioinformatique (BioPython, NumPy, PyTorch)
- **C++** — notions, à approfondir pour les algos performants (alignement, graphes)
- **Rust** — potentiellement à introduire pour la performance (bioinformatique moderne)
- **Snakemake / Nextflow** — pipelines bioinformatiques (projets avancés)

---

## Sources de données

Référence rapide pour trouver des séquences de test :

| Source | Usage | URL |
|--------|-------|-----|
| NCBI Nucleotide | Séquences ADN/ARN annotées | https://www.ncbi.nlm.nih.gov/nucleotide |
| NCBI SRA | Reads bruts de séquençage (FASTQ) | https://www.ncbi.nlm.nih.gov/sra |
| UniProt | Séquences protéiques | https://www.uniprot.org |
| Ensembl | Génomes annotés (humain, souris…) | https://www.ensembl.org |
| UCSC Genome Browser | Génomes + pistes d'annotation | https://genome.ucsc.edu |

Pour les TPs algorithmiques purs, préférer des **données simulées** (générateur fourni dans le squelette) afin de contrôler exactement les cas de test.

---

## Calibration des estimations de temps

Chaque TP contient un tableau `Suivi de temps` (estimé / réel / points bloquants). Ces données alimentent les estimations futures.

**Lors de l'initialisation d'un nouveau projet :**
- Consulter `PROJECTS.md` et les `REVIEW.md` des projets précédents
- Identifier les catégories de TPs systématiquement sous-estimées (ex : parsing, débogage d'algo de graphe)
- Ajuster les fourchettes en conséquence — préférer une fourchette large et honnête à une estimation optimiste

**Patterns à surveiller :**
- TPs avec structures de données nouvelles → toujours majorer
- TPs de deep learning → la mise au point des hyperparamètres est chronophage
- Premier TP d'un projet → setup environnement non comptabilisé, majorer de 30 %

---

## Clôture d'un projet

Quand tous les TPs d'un projet sont validés :

1. L'étudiant (ou Claude en mode review final) rédige `SYNTHESIS.md` dans le dossier projet : ce qui a été construit, ce qui a bien marché, ce qui aurait pu être mieux conçu dès le départ
2. Mettre à jour `PROJECTS.md` : statut → `Terminé`, remplir le nombre de TPs complétés
3. Vérifier `SKILLS.md` : toutes les compétences du projet sont bien en `covered`
