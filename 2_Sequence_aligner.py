import numpy as np

# define scores
match = 5
mismatch = -2
indel = -6


def AlignByDP(sequences):
    '''This function takes a list of tuples of species names and nucleotide sequences as inputs.
    It alignes the squences according to a score (match, mismatch and indel)

    Args:
        sequences (list(tuple(string, string))): a list of tuples tuples of species names and nucleotide sequences

    Returns:
        alignements dict(tuple(int, int) -> tuple(string, string)): a dictionary with species paires and their sequences
    '''

    #check if input is a list
    if not type(sequences) == list:
        raise Exception("malformed input")

    alignements = dict()

    for alignment in Combinations(sequences):

        sequence_1 = sequences[alignment[0]][1]
        sequence_2 = sequences[alignment[1]][1]
        name_1 = sequences[alignment[0]][0]
        name_2 = sequences[alignment[1]][0]

        # check if all elements in the list are tuples
        if not type(sequences[alignment[0]]) == tuple:
            raise Exception("malformed input")
        # check if all elements in the tuple are strings
        if not type(sequences[alignment[0]][0]) == str:
            raise Exception("malformed input")


        M = len(sequence_1)
        N = len(sequence_2)

        #create a traceback matrix
        traceback_mat = Fill_Matrices(M, N, sequence_1, sequence_2)

        traceback_index = traceback_mat[M, N]
        aligned_seq_1 = []
        aligned_seq_2 = []

        while M > 0 or N > 0:
            if traceback_index == 0:
                aligned_seq_1.insert(0, sequence_1[M - 1])
                aligned_seq_2.insert(0, sequence_2[N - 1])
                M -= 1
                N -= 1

            if traceback_index == 1:
                aligned_seq_1.insert(0, sequence_1[M - 1])
                aligned_seq_2.insert(0, '-')
                M -= 1

            if traceback_index == 2:
                aligned_seq_1.insert(0, '-')
                aligned_seq_2.insert(0, sequence_2[N - 1])
                N -= 1

            traceback_index = traceback_mat[M, N]

        alignements[alignment] = (''.join(aligned_seq_1), ''.join(aligned_seq_2))

    return alignements


def Traceback_Matrix(M, N):
    '''This function creates a traceback matrix
    Args:
        M (int): length squence 1
        N (int): length squence 2
    Returns:
        traceback_matrix (numpy array)
    '''
    traceback_matrix = np.zeros((M + 1, N + 1))
    traceback_matrix[0, 1:] = 2
    traceback_matrix[1:, 0] = 1

    return traceback_matrix

def Dynmaic_Matrix(M, N):
    '''This function creates a dynamic matrix
    Args:
        M (int): length sequence 1
        N (int): length sequence 2
    Returns:
        dp_matrix (numpy array)

    '''
    dp_matrix = np.zeros((M + 1, N + 1))
    dp_matrix[0, 1:] = [(i+1) * indel for i in range(N)]
    dp_matrix[1:, 0] = [(i+1) * indel for i in range(M)]

    return dp_matrix

def Fill_Matrices(M, N, sequence_1, sequence_2):
    '''This function calculates the scores for the dynamic matrix and the traceback matrix with the dynamic programming methos
    Args:
        M (int): length of sequence 1
        N (int): length of sequence 2
        sequence_1 (str): nucleotide sequence
        sequence_2 (str): nucleotide sequence

    Returns:
        tb_matrix (numpy array): array with
    '''
    dp_matrix = Dynmaic_Matrix(M, N)
    tb_matrix = Traceback_Matrix(M, N)

    #iterate over dynamic-matrix and traceback-matrix at once to fill the scores
    for i in range(1, M):
        for j in range(1, N):
            score = Scores(dp_matrix, i, j, sequence_1, sequence_2)
            dp_matrix[i, j] = max(score)
            tb_matrix[i, j] = np.argmax(score)

    return tb_matrix

def Scores(dp_matrix, i, j, sequence_1, sequence_2):
    '''This function calculates the scores in the matrix
    Args:
        dp_matrix (int): length of sequence 1
        i(int): length of sequence 1
        j(int): length of sequence 2
        sequence_1 (str): nucleotide sequence
        sequence_2 (str): nucleotide sequence

    Returns:
        diagonal, up, left (tuple): three values with the scores diagonal, up and left from the entry
    '''
    # raise exception if sequence not only contains 'ATCG'
    if sequence_1[i-1] not in ('ATCG') or sequence_2[j-1] not in ('ATCG'):
        raise Exception("malformed input")

    if sequence_1[i-1] == sequence_2[j-1]:
        diagonal = dp_matrix[i-1, j-1] + match
    else:
        diagonal = dp_matrix[i-1, j-1] + mismatch

    up = dp_matrix[i-1, j] + indel
    left = dp_matrix[i, j-1] + indel

    return diagonal, up, left

def Combinations(sequences):
    '''this function returns all the possible combinations from a list
    Args:
        sequences (list(tuple(string, string))): a list of tuples tuples of species names and nucleotide sequences
    Returns:
        list(tuple)

    '''
    l = range(len(sequences))
    return [(i, j) for i in l for j in l[i+1:]]
