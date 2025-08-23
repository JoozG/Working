"""
FUNCTION: 'read_fasta' loads a fasta file, returns a dictionary

INPUT: string(file_path)
OUTPUT: dictionary(sequences)
"""

def read_fasta(file_path: str) -> dict[str, str]:

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
FUNCTION: 'read_fasta_fileobj' - loads FASTA from a file-like object with .read().

INPUT:  file-like (fobj)
OUTPUT: dictionary (id -> sequence)
"""

def read_fasta_fileobj(fobj) -> dict[str, str]:
    # fobj: ex. UploadedFile z Streamlit (ma .read() / .getvalue())
    text = fobj.read().decode('utf-8', errors='ignore') if hasattr(fobj, 'read') else str(fobj)
    return read_fasta_text(text)

"""
FUNCTION: 'read_fasta_bytes' - loads FASTA from raw bytes (e.g. Streamlit upload).

INPUT:  bytes (data)
OUTPUT: dictionary (id -> sequence)
"""
def read_fasta_bytes(data: bytes) -> dict[str, str]:
    text = data.decode('utf-8', errors='ignore')
    return read_fasta_text(text)

"""
FUNCTION: 'read_fasta_text' - loads FASTA data from a string.

INPUT:  string (text)
OUTPUT: dictionary (id -> sequence)
"""

def read_fasta_text(text: str) -> dict[str, str]:
    sequences, current_id, current_sequence = {}, None, []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                sequences[current_id] = ''.join(current_sequence)
            current_id, current_sequence = line[1:], []
        else:
            current_sequence.append(line)
    if current_id is not None:
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
    FUNCTION: 'filter_fasta' filters sequences and creates new file with filtered out sequences

    INPUT:  string (input_file), string (output_file), function (filter_func(id, seq) -> bool)
    OUTPUT: None (creates filtered file)
"""

def filter_fasta(input_file, output_file, filter_func):

    sequences = read_fasta(input_file)
    filtered_sequences = {id: seq for id, seq in sequences.items() if filter_func(id, seq)}
    write_fasta(filtered_sequences, output_file)

"""
FUNCTION: 'count_sequences' - counts the number of sequences in a FASTA file.

INPUT:  string (file_path)
OUTPUT: int (number of sequences)
"""

def count_sequences(file_path):
    
    sequences = read_fasta(file_path)
    return len(sequences)


"""
FUNCTION: 'split_fasta' - splits a FASTA file into chunks.

INPUT:  string (input_file), string (output_dir), int (chunk_size)
OUTPUT: None (creates multiple FASTA files in output_dir)
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
