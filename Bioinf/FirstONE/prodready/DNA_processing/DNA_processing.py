import stringOperations
import RNAProteinConverter
import fileOperations

class DNAOperations:
    def __init__(self, sequence):
        self.sequence = sequence
        self.table = self.to_list()
        self.codons = self.to_codone()

    '''
    FUNCTION 'to_list' - converts string to list
    INPUT - string (variable name: sequence)
    OUTPUT - list (variable name: sequence)
    '''
    def to_list(self):
        return list(self.sequence)
    
    '''
    FUNCTION 'to_codone' - converts string to list
    INPUT - list (variable name: sequence)
    OUTPUT - list (variable name: sequence)
    '''    
    def to_codone(self):
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
    FUNCTION 'transcribe_dna_to_rna' - finds a RNA strand corresponding to declared DNA strand
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
        
    def gc_content(self):
        gc_contents = 0
        if len(self.table) == 0:
            return 0
        count_c = self.table.count('C')
        count_g = self.table.count('G')
        gc_content = round(100 * (count_c + count_g) / len(self.table), 6)
        return gc_content