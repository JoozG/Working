import sys
import os

# Dodaj katalog nadrzędny (katalog główny repozytorium) do sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fasta_processing.protein_processing import ProteinOperations
from fasta_processing.file_processing import read_fasta

proteins = read_fasta(".FASTA_processing/sampledata/example-protein.fasta")

print(proteins)
