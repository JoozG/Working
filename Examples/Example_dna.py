from FASTA_processing.fasta_processing.dna_processing import dnaOperations

dnaStrand = dnaOperations("ATGCGTACGTAGCTAGCTTACGTAAGGCTTAGCTGAC")

dnaStrand.export_results(r"FASTA_processing/Examples/Results/results_dna.txt")