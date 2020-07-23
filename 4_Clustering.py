import numpy as np
from collections import defaultdict
#test branch update
def Cluster(pairwise_distances):
    '''This function takes a pairwise distance matrix and a list of labels as inputs.
    It clusters the labels according to the Weighted Pair Group Method with Arithmetic Mean (WPGMA) Algorithm.

    Args:
        pairwise_distances (list(list(float)), list): key-pair for the minimum value
        dist_dict (dict): 2d distance dictionary
    Returns:
        (str): the function returns a string representing the clusters in ()
    '''
    distances = pairwise_distances[0]
    labels = pairwise_distances[1]

    #check if innput format is correct
    if not len(distances) == len(labels):
        raise Exception("malformed input")
    if not type(labels) == list:
        raise Exception("malformed input")
    if not type(distances) == list:
        raise Exception("malformed input")


    dist_dict = DISTANCE_DICT(distances, labels)

    while len(dist_dict) > 2:
        min_keys = FIND_MIN(dist_dict)

        new_entry_dict = DISTANCE(min_keys, dist_dict)

        #update dictionary
        del dist_dict[min_keys[0]], dist_dict[min_keys[1]]

        for item in dist_dict.keys():
            del dist_dict[item][min_keys[0]], dist_dict[item][min_keys[1]]
            dist_dict[item][min_keys] = new_entry_dict[item]

        new_entry_dict[min_keys] = 0
        dist_dict[min_keys] = new_entry_dict

    return str(tuple([k for k in dist_dict.keys()]))



def DISTANCE(min_keys, dist_dict):
    ''' This function extracts two dictionaries from key-pair with the minimum value
        and calculates the distances by averaging their values

    Args:
        min_index (tuple): key-pair for the minimum value
        dist_dict (dict): 2d distance dictionary
    Returns:
        cluster_distances (dict): returns a dictionary with the newly clustered keys and their averaged values
    '''

    min_dict_1 = dist_dict[min_keys[0]]
    min_dict_2 = dist_dict[min_keys[1]]

    cluster_distances = [(k, (v + min_dict_2[k])/2)
                         for k, v in min_dict_1.items()
                         if v > 0 and min_dict_2[k] > 0]

    return dict(cluster_distances)

def DISTANCE_DICT(pairwise_distances, labels):
    '''This function calculates the distances for the new cluster

    Args:
        pairwise_distances (list(list(float)): a pairwise distance matrix of the species listes in lables
        labels (list(str)): a list of string representing the species
    Returns:
        dist_dict (dict): the pairwise distance matrix as a 2d dictionary
    '''

    dist_dict = defaultdict(dict)

    for i in range(len(pairwise_distances)):
        for j in range(len(pairwise_distances)):
            try:
                dist_dict[labels[i]][labels[j]] = pairwise_distances[i][j]
            except:
                raise Exception("malformed input")

    return dist_dict

def FIND_MIN(dist_dict):
    '''This function find the minimum value in the dictionary

    Args:
        dist_dict (dict): pairwise distance matrix as a 2d dictionary

    Returns:
        min_keys (tuple): a key-pair for the minimum value in the dict
    '''

    keys = list(dist_dict.keys())

    minimum = dist_dict[keys[0]][keys[1]]
    min_keys = keys[0], keys[1]

    for i in dist_dict:
        for j in dist_dict[i]:
            if 0 < dist_dict[i][j] < minimum:
                min_keys = i, j
                minimum = dist_dict[i][j]

    return min_keys
