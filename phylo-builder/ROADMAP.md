# Roadmap — phylo-builder

**Durée estimée totale :** 19–27 heures
**Technologies :** Python

---

## Vue d'ensemble des TPs

| TP | Titre | Compétences | Estimé | Prérequis | Statut |
|----|-------|-------------|--------|-----------|--------|
| 01 | Matrices de distance évolutive | STA, BIO_BASE | 3–4 h | — | [ ] |
| 02 | UPGMA — clustering hiérarchique | PHY, STA | 4–6 h | TP01 | [ ] |
| 03 | Neighbor-Joining | PHY | 4–6 h | TP02 | [ ] |
| 04 | Format Newick + visualisation | PHY | 3–4 h | TP03 | [ ] |
| 05 | Bootstrap — robustesse statistique | STA, PHY | 3–4 h | TP04 | [ ] |
| 06 | Application sur données réelles | PHY, BIO_BASE | 2–3 h | TP05 | [ ] |

`[ ]` non démarré · `[~]` en cours · `[x]` validé

---

## Progression fonctionnelle

**TP01 →** Étant donné un ensemble de séquences ADN alignées, calculer toutes les distances pairwise (Hamming, Jukes-Cantor, Kimura 2P) et les stocker dans une `DistanceMatrix`.

**TP02 →** TP01 + construire un arbre enraciné par UPGMA (clustering hiérarchique agglomératif sur la `DistanceMatrix`). Introduire la structure de données `Tree`.

**TP03 →** TP02 + construire un arbre non-enraciné par Neighbor-Joining — algorithme plus précis quand les taux d'évolution varient entre lignées.

**TP04 →** TP03 + sérialiser n'importe quel `Tree` au format Newick, afficher un dendrogramme ASCII dans le terminal et un dendrogramme matplotlib dans un fichier image.

**TP05 →** TP04 + évaluer la robustesse de l'arbre par bootstrap : rééchantillonner les colonnes de l'alignement N fois, reconstruire un arbre à chaque itération, reporter les valeurs de support sur les nœuds.

**TP06 →** TP05 + pipeline complet sur de vraies séquences téléchargées depuis NCBI (cytochrome b ou 16S ARNr de 6–10 espèces de vertébrés). Produire un arbre annoté avec valeurs de bootstrap, comparer UPGMA vs NJ.

---

## Évolutions possibles hors roadmap

- **Modèles d'évolution avancés** : GTR (General Time Reversible), HKY85 — paramétrer les taux de substitution par nucléotide
- **Maximum de parcimonie** : construire l'arbre qui minimise le nombre total de changements (alternative à la distance)
- **Alignement multiple** : intégrer ClustalW ou MUSCLE pour aligner les séquences avant de calculer les distances (plutôt que supposer l'alignement fourni)
- **Parseur Newick** : lire un arbre Newick depuis un fichier externe (complémentaire à la génération)
- **Correction de Felsenstein** : distance LogDet pour des compositions en bases hétérogènes entre espèces
