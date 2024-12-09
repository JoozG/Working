from FASTA_processing.fasta_processing.rna2protein_converter import rnaConverter

rna_sequence = "AUGGCCUUUAAUGGGAUG"
converter = rnaConverter(rna_sequence)

# Tłumaczenie RNA na białko
protein = converter.rna_to_protein()
print("Protein sequence:", "".join(protein))

# Znalezienie ORF
orfs = converter.find_open_reading_frames()
print("Open Reading Frames (ORFs):", orfs)

# Analiza użycia kodonów
codon_usage = converter.analyze_codon_usage()
print("Codon Usage:", codon_usage)

# Eksport wyników
converter.export_results(r"FASTA_processing\Examples\Results\results_rna.txt")
