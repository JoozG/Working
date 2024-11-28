"""
FUNCTION: 'read_fasta' loads a fasta file, returns a dictionary

INPUT: string(file_path)
OUTPUT: dictionary(sequences)
"""

def read_fasta(file_path):

    sequences = {} 
    with open(file_path, 'r') as file: 
        current_id = None # setting the id as none
        current_sequence = []
        for line in file: 
            line = line.strip() # stripping current line to a variable
            if line.startswith(">"): 
                if current_id is not None:
                    sequences[current_id] = ''.join(current_sequence)
                current_id = line[1:]  # updating the current_id variable without '>'
                current_sequence = [] # clearing current_sequence variable
            else:
                current_sequence.append(line) # appending the current_sequence list if the line dosent start with '>'

        if current_id is not None: # appending the output variable ('sequences')
            sequences[current_id] = ''.join(current_sequence)
    return sequences

"""
FUNCTION: 'write_fasta' saves a .fasta file in a desired path

INPUT: dictionary(sequences), string(file_path), int(line_width within .fasta file)
OUTPUT: dictionary(sequences)
"""
def write_fasta(sequences, file_path, line_width=80):

    with open(file_path, 'w') as file:
        for seq_id, sequence in sequences.items():
            file.write(f">{seq_id}\n")
            for i in range(0, len(sequence), line_width):
                file.write(sequence[i:i + line_width] + '\n')

"""
    Filtruje sekwencje z pliku FASTA i zapisuje wynik w nowym pliku.

    Parameters:
        input_file (str): Ścieżka do pliku wejściowego FASTA.
        output_file (str): Ścieżka do pliku wyjściowego FASTA.
        filter_func (function): Funkcja przyjmująca ID i sekwencję jako argumenty
                                i zwracająca True, jeśli sekwencja ma zostać uwzględniona.
"""

def filter_fasta(input_file, output_file, filter_func):

    sequences = read_fasta(input_file)
    filtered_sequences = {id: seq for id, seq in sequences.items() if filter_func(id, seq)}
    write_fasta(filtered_sequences, output_file)

"""
    Liczy liczbę sekwencji w pliku FASTA.

    Parameters:
        file_path (str): Ścieżka do pliku FASTA.

    Returns:
        int: Liczba sekwencji w pliku.
    """

def count_sequences(file_path):
    
    sequences = read_fasta(file_path)
    return len(sequences)


"""
    Dzieli plik FASTA na mniejsze pliki o określonej liczbie sekwencji.

    Parameters:
        input_file (str): Ścieżka do pliku wejściowego FASTA.
        output_dir (str): Ścieżka do katalogu wyjściowego, gdzie zapisane zostaną podzielone pliki.
        chunk_size (int): Maksymalna liczba sekwencji w jednym pliku wyjściowym.
    """
def split_fasta(input_file, output_dir, chunk_size):
    
    import os
    sequences = read_fasta(input_file)
    sequence_items = list(sequences.items())
    for i in range(0, len(sequence_items), chunk_size):
        chunk = sequence_items[i:i + chunk_size]
        chunk_dict = dict(chunk)
        output_file = os.path.join(output_dir, f"chunk_{i // chunk_size + 1}.fasta")
        write_fasta(chunk_dict, output_file)
