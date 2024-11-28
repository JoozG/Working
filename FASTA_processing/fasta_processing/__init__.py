from file_processing import read_fasta, write_fasta
from dna_processing import dnaOperations
from protein_processing import proteinOperations
from rna2protein_converter import rnaConverter

# Eksportuj wszystkie istotne funkcje i klasy
__all__ = [
    "read_fasta",
    "write_fasta",
    "DNA_processing",
    "translate_to_protein",
    "RNAProteinConverter",
]
