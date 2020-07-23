
def ParseSeqFile(path):
    '''This function takes a txt-file as input and parsers its lines. Each input line should be marked with a '>'
        Each input line consist of a species name and a nucleotide sequence containing only 'ACTG'
        A cleaned up list of tuples are returned with species name and sequence.

    Args:
        path (str): absolute path to the input .txt-file

    Returns:
        output (list): the function returns a list of tuples
    '''
    output = []

    with open(path) as file:
        for line in file:
            if not line.startswith(('>', '\n', '\r\n')):
                raise Exception("malformed input")
            if line[0] == '>':
                line = line[1:].split()
                name = line[0]
                seperator = ','
                sequence = ''.join(line[1:])

                #raise exception if label contains whitespace and sequence not only contains 'ATCG'
                if len([i for i in sequence if i not in ('ATCG')]) > 0:
                    raise Exception("malformed input")

                #raise exception if sequence is empty
                if len(sequence) == 0:
                    raise Exception("malformed input")

            output.append((name, sequence))

    file.close()
    return output



