    
class rnaConverter:
    '''
    The dictionary 'CODON_DICT' below will be reversed shortly after being declared in the code to improve execution time.
    You may wonder why I reversed the dictionary instead of declaring it reversed in the first place.
    For the record, I did it purely for the sake of readability.
    '''
    CODON_DICT = {
        "Ala": ["GCU", "GCC", "GCA", "GCG"],
        "Ile": ["AUU", "AUC", "AUA"],
        "Arg": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "Leu": ["CUU", "CUC", "CUA", "CUG", "UUA", "UUG"],
        "Asn": ["AAU", "AAC"],
        "Lys": ["AAA", "AAG"],
        "Asp": ["GAU", "GAC"],
        "Met": ["AUG"],
        "Phe": ["UUU", "UUC"],
        "Cys": ["UGU", "UGC"],
        "Pro": ["CCU", "CCC", "CCA", "CCG"],
        "Gln": ["CAA", "CAG"],
        "Ser": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
        "Glu": ["GAA", "GAG"],
        "Thr": ["ACU", "ACC", "ACA", "ACG"],
        "Trp": ["UGG"],
        "Gly": ["GGU", "GGC", "GGA", "GGG"],
        "Tyr": ["UAU", "UAC"],
        "His": ["CAU", "CAC"],
        "Val": ["GUU", "GUC", "GUA", "GUG"],
        "STOP": ["UAA", "UGA", "UAG"]
    }

    REVERSED_CODON_DICT = {codon: key for key, codons in CODON_DICT.items() for codon in codons}

    def __init__(self, rna_sequence):
        self.rna_sequence = rna_sequence.upper()
        self.codons = self.split_into_codons()

    '''
    FUNCTION 'split_into_codons' - splitting RNA sequence into codones
    INPUT - string (self.rna_sequence)
    OUTPUT - list (codone list)
    '''
    def split_into_codons(self):

        if len(self.rna_sequence) % 3 != 0:
            raise ValueError("RNA sequence length must be a multiple of 3.")
        
        return [self.rna_sequence[i:i+3] for i in range(0, len(self.rna_sequence), 3)]

    '''
    FUNCTION 'rna_to_protein' - translates RNA to protein
    INPUT - list (self.codons)
    OUTPUT - list (protein sequence in a form of single characters)
    '''
    def rna_to_protein(self, stop_on_stop_codon=True):
        protein_sequence = []

        for codon in self.codons:
            amino_acid = self.REVERSED_CODON_DICT.get(codon, None)

            if amino_acid == "STOP" and stop_on_stop_codon:
                break
            if amino_acid:
                protein_sequence.append(amino_acid)

        return protein_sequence

    '''
    FUNCTION 'find_open_reading_frames' - identifies all open ORF (ORF)
    INPUT - list (self.codons)
    OUTPUT - list (ORF list)
    '''
    def find_open_reading_frames(self):
        orfs = []
        start_indices = [i for i, codon in enumerate(self.codons) if codon == "AUG"]

        for start_index in start_indices:
            protein_sequence = []

            for codon in self.codons[start_index:]:
                amino_acid = self.REVERSED_CODON_DICT.get(codon, None)

                if amino_acid == "STOP":
                    orfs.append(protein_sequence)
                    break

                if amino_acid:
                    protein_sequence.append(amino_acid)

        return orfs

    '''
    FUNCTION 'analyze_codon_usage' - analyzes frequency of codons in  RNA sequence
    INPUT - list (self.codons)
    OUTPUT - dictionary (dictionary with codons as keys and number of codon occurrences as values)
    '''
    def analyze_codon_usage(self):
        codon_usage = {codon: 0 for codon in self.REVERSED_CODON_DICT.keys()}

        for codon in self.codons:
            if codon in codon_usage:
                codon_usage[codon] += 1

        return codon_usage

    '''
    FUNCTION 'export_results' - exports analysis result to a .txt file
    INPUT - string (filename - path to the output file)
    OUTPUT - None
    '''
    def export_results(self, filename):

        with open(filename, "w") as f:

            f.write("Codon Usage:\n")
            codon_usage = self.analyze_codon_usage()

            for codon, count in codon_usage.items():
                f.write(f"{codon}: {count}\n")

            f.write("\nOpen Reading Frames (ORFs):\n")
            orfs = self.find_open_reading_frames()

            for i, orf in enumerate(orfs, start=1):
                f.write(f"ORF {i}: {' '.join(orf)}\n")
