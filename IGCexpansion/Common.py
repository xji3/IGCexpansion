# A file to store common functions that avoid iteratively import itself
# Xiang Ji
# xji3@ncsu.edu
import itertools
import numpy as np
from math import floor

def divide_configuration(configuration):
    ortho_group_to_pos = dict(extent = {}, distinct = [], loc = [])
    # extent positions that represent same paralog (by the same ortholog group number) have to be in the same state
    # distinct positions don't change states, thus only track positions
    for pos in range(len(configuration)):
        if configuration[pos][1] == 1: # extent
            ortho_group = configuration[pos][0]
            if ortho_group in ortho_group_to_pos['extent']:
                ortho_group_to_pos['extent'][ortho_group].append(pos)
            else:
                ortho_group_to_pos['extent'][ortho_group] = [pos]
                ortho_group_to_pos['loc'].append(ortho_group)
        elif configuration[pos][1] == 0: # distinct
            ortho_group_to_pos['distinct'].append(pos)

    return ortho_group_to_pos 

def get_accessible_orlg_pair(conf_list):
    accessible_orlg_pair = list() # only store pair in an upper triangular fashion
    for conf in conf_list:
        ortho_group_to_pos = divide_configuration(conf)
        accessible_orlg_pair.extend(itertools.combinations(sorted(ortho_group_to_pos['loc']), 2))

    return sorted(list(set(accessible_orlg_pair)))

def translate_two_nt_to_one_state(state):
    assert(len(state) == 2)
    assert(all([-1 < i < 4 for i in state]))
    new_state = 4 * state[0] + state[1]
    return new_state

def translate_one_state_to_two_nt(state):
    assert( -1 < state < 16)
    new_state = (int(floor(state/4)), state % 4)
    return new_state

def translate_four_nt_to_two_state(state):
    assert(len(state) == 4)
    assert(all([-1 < i < 4 for i in state]))
    return (translate_two_nt_to_one_state(state[:2]), translate_two_nt_to_one_state(state[2:]))
           
def translate_two_state_to_four_nt(state):
    assert(len(state) == 2)
    assert(all([-1 < i < 16 for i in state]))
    
    return translate_one_state_to_two_nt(state[0]) + \
           translate_one_state_to_two_nt(state[1])


def draw_from_distribution(prob, size, values):
    # convert into valid pmf first
    prob = np.array(prob)/sum(prob)
    bins = np.add.accumulate(prob)
    if size == 1:
        return values[np.digitize(np.random.random_sample(size), bins)[0]]
    else:
        return [values[i] for i in np.digitize(np.random.random_sample(size), bins)]
 
