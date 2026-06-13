"""
TP01 — Représentation des séquences + parsing FASTA
Implémenter toutes les méthodes marquées raise NotImplementedError().
Ne pas modifier les signatures ni les docstrings.
"""
from __future__ import annotations
from enum import Enum
import os

class SequenceType(Enum):
    DNA = "DNA"
    RNA = "RNA"


class Sequence:
    """
    Représente une séquence biologique (ADN ou ARN).

    Attributes:
        header: Identifiant de la séquence (sans le '>' du format FASTA).
        data: Chaîne de nucléotides, toujours stockée en majuscules.
        seq_type: Type de séquence (DNA ou RNA).
    """

    COMPLEMENT_DNA: dict[str, str] = {"A": "T", "T": "A", "G": "C", "C": "G"}
    COMPLEMENT_RNA: dict[str, str] = {"A": "U", "U": "A", "G": "C", "C": "G"}
    ALPHABET_DNA: frozenset[str] = frozenset("ATGC")
    ALPHABET_RNA: frozenset[str] = frozenset("AUGC")

    def __init__(
        self,
        header: str,
        data: str,
        seq_type: SequenceType = SequenceType.DNA,
    ) -> None:
        """
        Initialise une séquence.

        - Normaliser data en majuscules avant toute validation.
        - Vérifier que chaque caractère appartient à l'alphabet du seq_type.
        - Lever ValueError avec un message clair si un caractère est invalide.
        """
        self._header = header
        self._seq_type = seq_type

        data = data.upper()
        alphabet = self.ALPHABET_DNA if self._seq_type == SequenceType.DNA else self.ALPHABET_RNA
        invalid = set(nuc for nuc in data if nuc not in alphabet)
        if invalid:
            raise ValueError(
                f"Invalid nucleotide(s): {', '.join(invalid)}"
            )
        self._data = data

    def __len__(self) -> int:
        """Retourne le nombre de nucléotides dans la séquence."""
        return len(self._data)

    def __str__(self) -> str:
        """Retourne la séquence sous forme de chaîne brute (ex: 'ATGC')."""
        return self._data

    def __repr__(self) -> str:
        """Retourne une représentation lisible, ex: Sequence(header='brca1', len=8, type=DNA)"""
        return (
            f"Sequence(header='{self._header}', "
            f"len={len(self)}, "
            f"type={self._seq_type.value})"
        )
    def __eq__(self, other: object) -> bool:
        """Deux séquences sont égales si header, data et seq_type sont identiques."""
        if not isinstance(other, Sequence):
            return False
        return self._header == other._header and self._data == other._data and self._seq_type == other._seq_type

    def complement(self) -> str:
        """
        Retourne le complément de la séquence (sans inverser).
        Utiliser COMPLEMENT_DNA ou COMPLEMENT_RNA selon seq_type.
        """
        complement = self.COMPLEMENT_DNA if self._seq_type == SequenceType.DNA else self.COMPLEMENT_RNA
        return "".join(complement[nuc] for nuc in self._data)


    def reverse_complement(self) -> str:
        """
        Retourne le complément inverse (brin opposé lu en 5'→3').
        Indice : deux opérations simples que tu as déjà.
        """
        return self.complement()[::-1]

    def gc_content(self) -> float:
        """
        Retourne la proportion de bases G et C (entre 0.0 et 1.0).
        Retourner 0.0 si la séquence est vide.
        """
        if len(self) == 0:
            return 0.0

        return sum(nuc in {"G", "C"} for nuc in self._data) / len(self)

    def transcribe(self) -> Sequence:
        """
        Transcrit l'ADN en ARN : remplace chaque T par U.
        Lever TypeError si seq_type n'est pas DNA.
        Retourner une nouvelle Sequence de type RNA avec le même header.
        """
        if self._seq_type == SequenceType.RNA:
            raise TypeError("Cannot transcribe: sequence is already RNA")
        return Sequence(self._header, self._data.replace("T","U"), SequenceType.RNA)

    def subseq(self, start: int, end: int) -> Sequence:
        """
        Retourne la sous-séquence data[start:end].
        Conserver le même header avec le suffixe '[start:end]'.
        Respecter la sémantique Python standard pour les indices négatifs et hors bornes.
        """
        return Sequence(
            header=f"{self._header}[{start}:{end}]",
            data=self._data[start:end],
            seq_type=self._seq_type,
        )

    @property
    def header(self) -> str:
        return self._header

    def __getitem__(self, key):
        return self._data[key]


def read_fasta(filepath: str) -> list[Sequence]:
    """
    Lit un fichier FASTA et retourne la liste des séquences dans l'ordre du fichier.

    Format FASTA attendu :
        >header description_optionnelle
        SÉQUENCE_LIGNE_1
        SÉQUENCE_LIGNE_2        ← les séquences multi-lignes doivent être concaténées

    Cas limites à gérer :
        - Fichier inexistant          → lever FileNotFoundError
        - Fichier vide                → retourner []
        - Lignes vides entre séquences → ignorer
        - Séquence sans données       → inclure avec data=""
    """
    sequences: list[Sequence] = []

    with open(filepath, "r", encoding="utf-8") as f:
        header: str | None = None
        data_parts: list[str] = []

        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    sequences.append(
                        Sequence(header, "".join(data_parts))
                    )

                header = line[1:]
                data_parts = []
            else:
                data_parts.append(line)

        if header is not None:
            sequences.append(
                Sequence(header, "".join(data_parts))
            )

    return sequences


def write_fasta(
    sequences: list[Sequence],
    filepath: str,
    line_width: int = 60,
) -> None:
    """
    Écrit une liste de séquences au format FASTA dans filepath.

    - Chaque séquence commence par '>header\n'.
    - La séquence est découpée en lignes de line_width caractères maximum.
    - Créer les répertoires parents si nécessaire (os.makedirs avec exist_ok=True).
    - Liste vide → fichier vide.
    """
    os.makedirs(os.path.dirname(filepath) or ".", exist_ok=True)

    with open(filepath, "w", encoding="utf-8") as f:
        for seq in sequences:
            f.write(f">{seq.header}\n")

            data = str(seq)
            for i in range(0, len(data), line_width):
                f.write(data[i:i + line_width] + "\n")
