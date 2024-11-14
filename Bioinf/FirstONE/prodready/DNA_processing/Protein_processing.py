import stringOperations
import fileOperations

#Protein_processing

class ProteinOperations:

    AMINO_ACIDS_DICT = {
        "Ala": {"symbol": "A", "property1": "value1", "property2": "value2"},
        "Ile": {"symbol": "I", "property1": "value1", "property2": "value2"},
        "Arg": {"symbol": "R", "property1": "value1", "property2": "value2"},
        "Leu": {"symbol": "L", "property1": "value1", "property2": "value2"},
        "Asn": {"symbol": "N", "property1": "value1", "property2": "value2"},
        "Lys": {"symbol": "K", "property1": "value1", "property2": "value2"},
        "Asp": {"symbol": "D", "property1": "value1", "property2": "value2"},
        "Met": {"symbol": "M", "property1": "value1", "property2": "value2"},
        "Phe": {"symbol": "F", "property1": "value1", "property2": "value2"},
        "Cys": {"symbol": "C", "property1": "value1", "property2": "value2"},
        "Pro": {"symbol": "P", "property1": "value1", "property2": "value2"},
        "Gln": {"symbol": "Q", "property1": "value1", "property2": "value2"},
        "Ser": {"symbol": "S", "property1": "value1", "property2": "value2"},
        "Glu": {"symbol": "E", "property1": "value1", "property2": "value2"},
        "Thr": {"symbol": "T", "property1": "value1", "property2": "value2"},
        "Trp": {"symbol": "W", "property1": "value1", "property2": "value2"},
        "Gly": {"symbol": "G", "property1": "value1", "property2": "value2"},
        "Tyr": {"symbol": "Y", "property1": "value1", "property2": "value2"},
        "His": {"symbol": "H", "property1": "value1", "property2": "value2"},
        "Val": {"symbol": "V", "property1": "value1", "property2": "value2"},
        "STOP": {"symbol": "|STOP|", "property1": "value1", "property2": "value2"}
    }

    def __init__(self, sequence):
        self.sequence = sequence
        self.table = self.tolist()


    '''
    FUNCTION 'to_list' - converts string to list
    INPUT - string (variable name: sequence)
    OUTPUT - list (variable name: sequence)
    '''
    def to_list(self):
        return list(self.sequence)


    