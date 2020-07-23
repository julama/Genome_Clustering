import numpy as np

def ComputeDistMatrix(dict_input):
    '''This function calculates a pairwise distance matrix

    Args:
        dict_input (dict(tuple(int, int) -> tuple(string, string)): a dictionary with species paires and their sequences
    Returns:
        DistMatrix (list(list(float)): a pairwise distance matrix


    '''
    #check if input is a dict
    if not type(dict_input) == dict:
        raise Exception("malformed input")

    #l -> maximum value from the dict key tuples (+1)
    try:
        l = max([k for k in dict_input.keys()])[1] + 1
    except:
        raise Exception("malformed input")


    DistMatrix = [[0 for i in range(l)] for j in range(l)]

    for item in dict_input:

        #check if keys are tuples
        if not type(item) == tuple or not type(dict_input[item]) == tuple:
            raise Exception("malformed input")

        sequence_1 = dict_input[item][0]
        sequence_2 = dict_input[item][1]

        nucleotide_pairs = delete_indel(sequence_1, sequence_2)

        distance = compute_distance(nucleotide_pairs)

        #raise exception if key tuples entries are not int
        try:
            DistMatrix[item[0]][item[1]] = distance
            DistMatrix[item[1]][item[0]] = distance
        except:
            raise Exception("malformed input")

    return DistMatrix



def delete_indel(seq1, seq2):
    '''This function deletes indels and checks for malfromed input sequences
    Args:
        seq1 (str): nucleotide sequence
        seq2 (str): nucleotide sequence
    Returns:
        nucleotide_pairs (list(tuple)): paired nucleotides without '-'
    '''
    #check format of sequences
    if not type(seq1) == str or not type(seq2) == str:
        raise Exception("malformed input")

    if len(seq1) != len(seq2):
        raise Exception("malformed input")

    # raise exception if label contains whitespace and sequence not only contains 'ATCG'
    if len([i for i in seq1+seq2 if i not in ('ATCG-')]) > 0:
        raise Exception("malformed input")

    nucleotide_pairs = [i for i in zip(seq1, seq2)
                        if not (i[0] == '-' or i[1] == '-')]

    return nucleotide_pairs


def compute_distance(nucleotide_pairs):
    '''This function calculates the distances in the matrix
    Args:
        nucleotide_pairs (list(tuple)): paired nucleotides without '-'
    Returns:
        distance (float): the distance between the pairs
    '''
    count_differences = len([i for i, pair in enumerate(nucleotide_pairs)
                             if pair[0] != pair[1]])

    total_length = len(nucleotide_pairs)
    p = count_differences / total_length

    if p >= 3/4:
        distance = 30
    else:
        distance = - 3/4 * np.log(1 - 4/3 * p)

    return distance

