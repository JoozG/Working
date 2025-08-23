from FASTA_processing.fasta_processing.protein_processing import proteinOperations

protein_sequence = "FVNQHLCGSHLVEALYLVCGERGFFYTPKA"

protein = proteinOperations(protein_sequence)

protein.export_results(r"FASTA_processing\Examples\Results\results_protein.txt")