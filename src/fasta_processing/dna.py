class dnaOperations:
    def __init__(self, sequence):
        self.sequence = sequence
        self.validate_sequence()
        self.table = self.sequence_to_list()
        self.codons = self.sequence_to_codons()

    '''
    FUNCTION 'validate_sequence' - checks if the input string contains only the specified characters A/C/G/T

    INPUT - string (self.sequence)
    OUTPUT - None (raises ValueError if invalid character is found)
    '''

    def validate_sequence(self):
        valid_nucleotides = {'A', 'C', 'G', 'T'}
        for nuc in self.sequence:
            if nuc not in valid_nucleotides:
                raise ValueError(f"Invalid nucleotide {nuc} in sequence")

    '''
    FUNCTION 'sequence_to_list' - converts string to list

    INPUT - string (self.equence)
    OUTPUT - list
    '''
    def sequence_to_list(self):
        return list(self.sequence.strip().upper()) 
    
    '''
    FUNCTION 'sequence_to_codons' - returns a list of codons (3-letter sequences)

    INPUT - list (self.sequence)
    OUTPUT - list 
    '''    
    def sequence_to_codons(self):
        return [self.sequence[i:i+3] for i in range(0, len(self.sequence), 3)]

    '''
    FUNCTION 'nuc_count' - counts occurences of nucleotides in a sequence

    INPUT - list (self.table)
    OUTPUT - dictionary (counts)
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
    
    INPUT - list (self.table)
    OUTPUT - list (trans)
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
    FUNCTION 'complement' - finds a complementary DNA strand, not reversed
    INPUT - list (self.table)
    OUTPUT - list 
    '''
    def complement(self):
        comp = {'A':'T','C':'G','G':'C','T':'A'}
        return [comp.get(n,'-') for n in self.table]

    '''
    FUNCTION 'reverse_complement' - finds a reverse-complementary DNA strand
    INPUT - list (self.table)
    OUTPUT - list 
    '''
    def reverse_complement(self):
        return self.complement()[::-1]

    '''
    FUNCTION 'gc_content' - calculates the GC content (up to 6 decimal places)

    INPUT - list (self.table)
    OUTPUT - float (gc content in percentage)
    '''       
    def gc_content(self):
        if not self.table: #Checking if the table is not empty
            return 0
    
        gc_count = sum(1 for nuc in self.table if nuc in ['G', 'C'])

        return round(100 * gc_count / len(self.table), 6) #Rounding up to 6 decimal places
    

    
    def export_results(self, filename):
        """
        FUNCTION: 'export_results' - exports nucleotide counts, reverse-complement,
        and GC content to a text file.

        INPUT:  string (filename)
        OUTPUT: None (creates a file with results)
        """

        with open(filename, "w") as f:
           f.write("Nucleotide usage:\n")
           
           nuc_counts = self.nuc_count()

           for nuc, count in nuc_counts.items():
                f.write(f"{nuc}: {count}\n")
           compl_dna = self.reverse_complement()

           f.write(f"\nComplementary DNA strand: 5'-{'-'.join(compl_dna)}-3'\n") # 5' to 3' 

           f.write(f"\nGC-content: {self.gc_content()} %")