from fasta_processing.rna2protein_converter import RNAConverter

rna_sequence = "AUGGCCUUUAAUGGGAUG"
converter = RNAConverter(rna_sequence)

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
converter.export_results(r"FirstONE\prodready\DNA_processing\Examples\analysis_results.txt")
