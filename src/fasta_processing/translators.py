    
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
    AA3_TO_AA1 = {"Ala":"A","Ile":"I","Arg":"R","Leu":"L","Asn":"N","Lys":"K","Asp":"D","Met":"M",
                  "Phe":"F","Cys":"C","Pro":"P","Gln":"Q","Ser":"S","Glu":"E","Thr":"T","Trp":"W",
                  "Gly":"G","Tyr":"Y","His":"H","Val":"V","STOP":"*"}

    def __init__(self, rna_sequence: str, frame_policy: str = "truncate"):
        self.rna_sequence = rna_sequence.upper()

        rem = len(self.rna_sequence) % 3
        if rem != 0:
            if frame_policy == "truncate":
                self.rna_sequence = self.rna_sequence[: len(self.rna_sequence) - rem]
            elif frame_policy == "pad":
                self.rna_sequence = self.rna_sequence + ("A" * (3 - rem))
            elif frame_policy == "error":
                raise ValueError("RNA sequence length must be a multiple of 3.")
            else:
                # fallback to truncate if unknown policy
                self.rna_sequence = self.rna_sequence[: len(self.rna_sequence) - rem]

        self.codons = self.split_into_codons()
    '''
    FUNCTION 'split_into_codons' - splitting RNA sequence into codons
    INPUT - string (self.rna_sequence)
    OUTPUT - list (codon list)
    '''
    def split_into_codons(self):
        return [self.rna_sequence[i:i+3] for i in range(0, len(self.rna_sequence), 3)]

    '''
    FUNCTION 'rna_to_protein' - translates RNA to protein
    INPUT - list (self.codons)
    OUTPUT - list (protein sequence in a form of single characters)
    '''
    def rna_to_protein(self, stop_on_stop_codon: bool = True, one_letter: bool = True):
        aa_list = []

        for codon in self.codons:
            aa3 = self.REVERSED_CODON_DICT.get(codon, None)

            if aa3 == "STOP":
                if stop_on_stop_codon:
                    break
                aa_list.append("*" if one_letter else "STOP")

            elif aa3:
                aa_list.append(self.AA3_TO_AA1[aa3] if one_letter else aa3)

        return aa_list

    '''
    FUNCTION 'find_open_reading_frames' - identifies all open ORFs starting with AUG and ending with a STOP codon
    INPUT - list (self.codons)
    OUTPUT - list (ORF list)
    '''
    def find_open_reading_frames(self, one_letter: bool = True, include_partial: bool = True, min_len: int = 0):

        orfs = []
        start_indices = [i for i, codon in enumerate(self.codons) if codon == "AUG"]

        for start_index in start_indices:
            aa_list = []

            for codon in self.codons[start_index:]:
                aa3 = self.REVERSED_CODON_DICT.get(codon, None)

                if aa3 == "STOP":
                    orfs.append(aa_list)
                    break

                if aa3:
                    aa_list.append(self.AA3_TO_AA1[aa3] if one_letter else aa3)
            
            # include ORF that reaches the end without STOP if requested
            if aa_list and (include_partial or (start_index + len(aa_list) < len(self.codons) and
                                               self.REVERSED_CODON_DICT.get(self.codons[start_index + len(aa_list)]) == "STOP")):
                if len(aa_list) >= min_len:
                    orfs.append("".join(aa_list))
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
