import numpy as np

from DNA_processing import DNAOperations, stringOperations, fileOperations, RNAProteinConverter

def main():
    # Load the .fasta file
    fasta_reader = fileOperations.FASTAReader(r'C:\Users\Admin\Documents\Coding\Bioinf\FirstONE\prodready\DNA_processing\sampledata\homosapiens_brca1.fasta')
    sequences = fasta_reader.sequences

    amacidCounter = stringOperations.SubstringFinder
    # Print loaded sequences
    if fasta_reader.names and sequences:
        print("Loaded Sequences:")
        for name, seq in zip(fasta_reader.names, sequences):
            print(f"{name}: {seq}")

        # Calculate and print GC content
        dna_op = DNAOperations(sequences[0])

        gc_content = dna_op.gc_content()
        print("\nGC Content:")
        print(f"{name}: {gc_content}%")

        # DNA Operations for the first sequence
        
        nucleotide_counts = dna_op.nuc_count()
        print("\nNucleotide Counts for the first sequence:")
        for nuc, count in nucleotide_counts.items():
            print(f"{nuc}: {count}")

        transcribed_rna = dna_op.transcribe_dna_to_rna()
        print(f"\nTranscribed RNA: {''.join(transcribed_rna)}")

        #complement_dna = dna_op.complement_dna()
        #print(f"\nComplement DNA: {''.join(complement_dna)}")

        # RNA to Protein Conversion for the transcribed RNA
        rna_seq = ''.join(transcribed_rna)
        rna_protein_converter = RNAProteinConverter.RNAConverter(rna_seq)
        protein_sequence = rna_protein_converter.rna_to_protein()
        print(f"\nProtein Sequence: {''.join(protein_sequence)}")

        # Hamming Distance between two sequences
        if len(sequences) > 1:
            string1 = sequences[0]
            string2 = sequences[1]
            hamming_distance = stringOperations.HammingDistance.calculate(string1, string2)
            print(f"\nHamming Distance between sequence1 and sequence2: {hamming_distance}")

            # Longest Common Substring between two sequences
            longest_common_substrings = stringOperations.SubstringFinder.longest_common_substrings(string1, string2)
            print(f"\nLongest Common Substrings between sequence1 and sequence2: {longest_common_substrings}")

        # Finding Substrings
        substring_positions = stringOperations.SubstringFinder.find_substring(sequences[0], "CGT")
        print(f"\nPositions of substring 'CGT' in sequence1: {substring_positions}")

    else:
        print("Failed to load any sequences from the file.")





   

if __name__ == "__main__":
    main()
