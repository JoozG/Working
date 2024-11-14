import numpy as np


'''
! ! ! WARNING ! ! !

The stringOperations.py script is currently redundant

'''




#Computing Hamming distance between 2 strings
class HammingDistance:
    @staticmethod
    def calculate(string, string2):
        return sum(el1 != el2 for el1, el2 in zip(string, string2))

#Finding a sequence of characters (substring) in a string
class SubstringFinder:

    def __init__(self, sequences):
        self.sequences = sequences

    @staticmethod
    def find_substring(string, substring):
        occur = []
        pointer = string.find(substring)
        while pointer != -1:
            occur.append(pointer)
            pointer = string.find(substring, pointer + 1)
        return occur if occur else 0

    @staticmethod
    def longest_common_substrings(string, string2):
        r, n = len(string), len(string2)
        queue_list = np.zeros((r + 1, n + 1), dtype=int)
        z = 0
        substr_list = []

        for i in range(1, r + 1):
            for j in range(1, n + 1):
                if string[i - 1] == string2[j - 1]:
                    if i == 1 or j == 1:
                        queue_list[i][j] = 1
                    else:
                        queue_list[i][j] = queue_list[i - 1][j - 1] + 1
                    if queue_list[i][j] > z:
                        z = queue_list[i][j]
                        substr_list = {string[i - z:i]}
                    elif queue_list[i][j] == z:
                        substr_list.add(string[i - z:i])
                else:
                    queue_list[i][j] = 0
        return substr_list

    @staticmethod
    def longest_common_substring(strings_list, common_strings_list):
        for string in strings_list:
            for common_string in common_strings_list:
                print(SubstringFinder.find_substring(string, common_string))

    def gc_content(self):
        gc_contents = []
        for seq in self.sequences:
            length = len(seq)
            if length == 0:
                continue
            count_c = seq.count('C')
            count_g = seq.count('G')
            gc_content = round(100 * (count_c + count_g) / length, 6)
            gc_contents.append(gc_content)
        return list(zip(self.names, gc_contents))

    def max_gc_content(self):
        gc_contents = self.gc_content()
        if not gc_contents:
            return None
        return max(gc_contents, key=lambda x: x[1])