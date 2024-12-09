class proteinOperations:
    # Tablica mas molekularnych aminokwasów (w daltonach)
    AMINO_ACID_MASSES = {
        'A': 89.09,  # Alanine
        'R': 174.2,  # Arginine
        'N': 132.12, # Asparagine
        'D': 133.1,  # Aspartic acid
        'C': 121.16, # Cysteine
        'E': 147.13, # Glutamic acid
        'Q': 146.15, # Glutamine
        'G': 75.07,  # Glycine
        'H': 155.16, # Histidine
        'I': 131.18, # Isoleucine
        'L': 131.18, # Leucine
        'K': 146.19, # Lysine
        'M': 149.21, # Methionine
        'F': 165.19, # Phenylalanine
        'P': 115.13, # Proline
        'S': 105.09, # Serine
        'T': 119.12, # Threonine
        'W': 204.23, # Tryptophan
        'Y': 181.19, # Tyrosine
        'V': 117.15  # Valine
    }
    
    def __init__(self, sequence):
        self.sequence = sequence.upper()
        self.amino_acid_list = self.sequence_to_list()
        self.validate_sequence()
    
    '''
    FUNCTION 'to_list' - converts string to list
    INPUT - string (self.sequence)
    OUTPUT - list (of single-character amino acid codes)
    '''
    def sequence_to_list(self):
        return list(self.sequence.strip().upper())    
    '''
    FUNCTION 'validate_sequence' - checks if the sequence contains only valid amino acid codes
    INPUT - string (self.sequence)
    OUTPUT - boolean
    '''
    def validate_sequence(self):
        valid_amino_acids = set(self.AMINO_ACID_MASSES.keys())
        invalid_chars = [aa for aa in self.sequence if aa not in valid_amino_acids]
        if invalid_chars:
            raise ValueError(f"Invalid amino acid(s) found: {invalid_chars}")
        return True
    
    '''
    FUNCTION 'amino_acid_count' - counts occurrences of each amino acid
    INPUT - list (self.amino_acid_list)
    OUTPUT - dictionary (with counts)
    '''
    def count_amino_acids(self):
        counts = {aa: 0 for aa in self.AMINO_ACID_MASSES.keys()}
        for aa in self.amino_acid_list:
            if aa in counts:
                counts[aa] += 1
        return counts

    '''
    FUNCTION 'molecular_weight' - calculates the molecular weight of the polypeptide
    INPUT - list (self.amino_acid_list)
    OUTPUT - float (molecular weight in daltons)
    '''
    def molecular_weight(self):
        total_weight = sum(self.AMINO_ACID_MASSES.get(aa, 0) for aa in self.amino_acid_list)
        return round(total_weight, 2)
    
    '''
    FUNCTION 'isoelectric_point' - simplistic isoelectric point calculation (based on average pKa)
    INPUT - list (self.amino_acid_list)
    OUTPUT - float (estimated pI)
    '''
    def isoelectric_point(self):
        acidic_residues = {'D': 3.9, 'E': 4.3}
        basic_residues = {'K': 10.5, 'R': 12.5, 'H': 6.0}
        
        total_acidic = sum(self.amino_acid_list.count(aa) * pKa for aa, pKa in acidic_residues.items())
        total_basic = sum(self.amino_acid_list.count(aa) * pKa for aa, pKa in basic_residues.items())
        
        # Simplified pI calculation (not fully accurate, for educational purposes)
        if total_acidic + total_basic == 0:
            return 7.0  # Neutral point for simplicity
        return round((total_acidic + total_basic) / (len(acidic_residues) + len(basic_residues)), 2)

    '''
    FUNCTION 'hydrophobicity_score' - calculates the hydrophobicity score of the sequence
    INPUT - list (self.amino_acid_list)
    OUTPUT - float (average hydrophobicity score)
    '''
    def hydrophobicity_score(self):
        hydrophobicity = {
            'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8, 
            'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
            'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5, 
            'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
        }
        total_score = sum(hydrophobicity.get(aa, 0) for aa in self.amino_acid_list)
        return round(total_score / len(self.amino_acid_list), 2) if self.amino_acid_list else 0

    def export_results(self, filename):

        with open(filename, "w") as f:

            f.write("Protein usage:\n")
            aminoacid_usage = self.count_amino_acids()

            for aa, count in aminoacid_usage.items():
                f.write(f"{aa}: {count}\n")

            f.write(f"\nMolecular weight: {self.molecular_weight()} Da\n")

            f.write(f"\nCalculated isoelectric point: {self.isoelectric_point()}\n")

            f.write(f"\nCalculated hydrophobicity score: {self.hydrophobicity_score()}\n")