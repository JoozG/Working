class RNAConverter:
    CODON_DICT = {
        "Ala": ["GCU", "GCC", "GCA", "GCG"],
        "Ile": ["AUU", "AUC", "AUA"],
        "Arg": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "Leu": ["CUU", "CUC", "CUA", "CUG", "UUA", "UUG"],
        "Asn": ["AAU", "AAC", "GAU", "GAC"],
        "Lys": ["AAA", "AAG"],
        "Asp": ["GAU", "GAC"],
        "Met": ["AUG"],
        "Phe": ["UUU", "UUC"],
        "Cys": ["UGU", "UGC"],
        "Pro": ["CCU", "CCC", "CCA", "CCG"],
        "Gln": ["CAA", "CAG", "GAA", "GAG"],
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
    
    AMINO_ACIDS_DICT = {
        "Ala": "A", "Ile": "I", "Arg": "R", "Leu": "L", "Asn": "N", 
        "Lys": "K", "Asp": "D", "Met": "M", "Phe": "F", "Cys": "C", 
        "Pro": "P", "Gln": "Q", "Ser": "S", "Glu": "E", "Thr": "T", 
        "Trp": "W", "Gly": "G", "Tyr": "Y", "His": "H", "Val": "V", 
        "STOP": "|STOP|"
    }

    def __init__(self, rna_sequence):
        self.rna_sequence = rna_sequence
        self.codons = self.split_into_codons()

    def split_into_codons(self):
        return [self.rna_sequence[i:i+3] for i in range(0, len(self.rna_sequence), 3)]

    def rna_to_protein(self):
        protein_sequence = []
        for codon in self.codons:
            for key, value in self.CODON_DICT.items():
                if codon in value:
                    protein_sequence.append(self.AMINO_ACIDS_DICT[key])
                    break
        return protein_sequence

