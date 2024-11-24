#Handling fasta files
class FASTAReader:
    def __init__(self, filename):
        self.filename = filename
        self.names, self.sequences = self.load_fasta()

    def load_fasta(self):
        names = []
        sequences = []
        current_sequence = ""

        try:
            with open(self.filename, 'r') as file:
                lines = file.readlines()
        except FileNotFoundError:
            print(f"File not found: {self.filename}")
            return names, sequences

        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                names.append(line)
            else:
                current_sequence += line

        if current_sequence:
            sequences.append(current_sequence)

        if not names or not sequences:
            print("No valid sequences found in the file.")
        
        return names, sequences