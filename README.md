# Bioinfo — Collection de projets

Portfolio de bioinformatique construit sur l'année 2026. Objectif : préparer une entrée en doctorat en bioinformatique après un diplôme d'ingénieur INSA (informatique et technologies de l'information).

**Thème central :** traitement de séquences ADN/ARN — approches algorithmiques classiques et deep learning.

---

## Projets

| # | Nom | Technos | Thème | Statut |
|---|-----|---------|-------|--------|
| 1 | [Aligneur de séquences](sequence-aligner/) | Python | Alignement pairwise ADN (NW, SW, BLOSUM, gaps affines) | Terminé |
| 2 | [Phylogénomique computationnelle](phylo-builder/) | Python | Matrices de distance, UPGMA, NJ, bootstrap | En cours |

---

## Roadmap globale

*Mise à jour à chaque nouveau projet lancé. L'ordre reflète une progression logique des compétences.*

| Phase | Période | Projet | Orientation | Compétences clés |
|-------|---------|--------|-------------|-----------------|
| 1 | 2026 | [Aligneur de séquences](sequence-aligner/) | Algorithmique | PD, ALI, MST |
| 2 | 2026 | [Phylogénomique computationnelle](phylo-builder/) | Algorithmique / Stats | PHY, STA |
| 3 | — | *(à définir)* | Deep Learning | — |
| 4 | — | *(à définir)* | Intégrateur | — |

Détail de chaque projet : voir son `ROADMAP.md` et `CLAUDE.md`.

---

## Structure du dépôt

```
bioinfo/
├── README.md               ← ce fichier
├── PROJECTS.md             ← registre détaillé des projets (temps, TPs, dépendances)
├── IDEAS.md                ← banque d'idées de projets à venir
├── SKILLS.md               ← graphe de compétences et couverture par projet
├── CLAUDE.md               ← instructions pour l'orchestrateur IA
└── templates/              ← templates réutilisables
    ├── project_CLAUDE.md
    ├── tp_README.md
    ├── tp_REVIEW.md
    ├── pyproject_algo.toml
    ├── pyproject_dl.toml
    └── pyproject_lib.toml
```

Chaque projet occupe un sous-dossier autonome :

```
nom-projet/
├── CLAUDE.md       ← contexte + règles du projet
├── ROADMAP.md      ← plan des TPs
└── subjects/
    ├── TP01_nom/
    │   ├── README.md   ← énoncé, pseudo-code, critères de validation
    │   ├── skeleton.*  ← fichier squelette à implémenter
    │   └── REVIEW.md   ← compte-rendu de review (généré après validation)
    └── TP02_nom/
        └── ...
```

---

## Approche pédagogique

Chaque projet est une suite de TPs itératifs. L'IA joue le rôle de tuteur :
- **pendant le TP** : accompagnement, pseudo-code, pas de code fourni
- **après validation** : review de code structurée par priorité (bloquant / dette / cosmétique)

---

## Compétences couvertes

Voir [`SKILLS.md`](SKILLS.md) pour le graphe interactif (Mermaid).

Domaines : algorithmique, biologie computationnelle, probabiliste & stats, deep learning, systèmes & performance.
