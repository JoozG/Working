import string_processing as string_processing
import rna2protein_converter as rna2protein_converter
import file_processing as file_processing

class dnaOperations:
    def __init__(self, sequence):
        self.sequence = sequence
        self.table = self.sequence_to_list()
        self.codons = self.sequence_to_codons()

    '''
    FUNCTION 'validate_sequence' - checks if the input string contains only the specified characters
    INPUT - string (variable name: sequence)
    OUTPUT - None
    '''

    def validate_sequence(self):
        valid_nucleotides = {'A', 'C', 'G', 'T'}
        for nuc in self.sequence:
            if nuc not in valid_nucleotides:
                raise ValueError(f"Invalid nucleotide {nuc} in sequence")

    '''
    FUNCTION 'to_list' - converts string to list
    INPUT - string (variable name: sequence)
    OUTPUT - list (variable name: undefined)
    '''
    def sequence_to_list(self):
        return list(self.sequence.strip().upper()) 
    
    '''
    FUNCTION 'to_codone' - converts string to list
    INPUT - list (variable name: sequence)
    OUTPUT - list (variable name: undefined)
    '''    
    def sequence_to_codons(self):
        return [self.sequence[i:i+3] for i in range(0, len(self.sequence), 3)]

    '''
    FUNCTION 'nuc_count' - counts occurences of nucleobases
    INPUT - list (variable name: table)
    OUTPUT - dictionary (variable name: counts)
    '''
    def nuc_count(self):
        counts = {
            'A': 0,
            'C': 0,
            'G': 0,
            'T': 0,
            'undefined': 0
        }
        for nuc in self.table:
            if nuc in counts:
                counts[nuc] += 1
            else:
                counts['undefined'] += 1
        return counts

    '''
    FUNCTION 'transcribe_dna_to_rna' - finds a RNA strand corresponding to declared DNA strand
    INPUT - list (variable name: table)
    OUTPUT - list (variable name: trans)
    '''
    def transcribe_dna_to_rna(self):
        trans = []
        for nuc in self.table:
            if nuc == 'T':
                trans.append('U')
            elif nuc in ('A', 'G', 'C'):
                trans.append(nuc)

            else: trans.append('-')
        return trans
    
    '''
    FUNCTION 'complement_dna' - finds a complementary DNA strand, for the one specified by the 
    INPUT - list (variable name: table)
    OUTPUT - list (variable name: revdna)
    '''
    def complement_dna(self):
        complement = {
            'A': 'T',
            'C': 'G',
            'G': 'C',
            'T': 'A'
        }
        revdna = [complement[nuc] for nuc in self.table][::-1]
        return revdna

    '''
    FUNCTION 'gc_content' - returns the gc content of a string  
    INPUT - list (variable name: table)
    OUTPUT - list (variable name: revdna)
    '''       
    def gc_content(self):
        if not self.table: #Checking if the table is not empty
            return 0
    
        gc_count = sum(1 for nuc in self.table if nuc in ['G', 'C'])

        return round(100 * gc_count / len(self.table), 6) #Rounding up to 6 decimal places