from Protein_processing import ProteinOperations

protein_seq = "ACDEFGHIKLMNPQRSTVWY"
protein = ProteinOperations(protein_seq)

try:
    protein.validate_sequence()
    print("Amino acid counts:", protein.amino_acid_count())
    print("Molecular weight:", protein.molecular_weight(), "Da")
    print("Isoelectric point:", protein.isoelectric_point())
    print("Hydrophobicity score:", protein.hydrophobicity_score())
        
except ValueError as e:
   print(e) 