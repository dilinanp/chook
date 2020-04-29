import numpy as np
import itertools
import random
from functools import reduce

def add_coupler(ising_coupler, hobo_set):
    """
    Converts a single Ising coupler into hobo format.
    Adds the resultant hobo couplers to hobo_set, and returns the energy
    difference dE.

    """
    ising_coupler = tuple(sorted(ising_coupler[:-1])) + (ising_coupler[-1], )
    
    k = len(ising_coupler)-1
    
    mult_fact = -2
    multiplier = -2
    
    for n in range(1, k+1):
        combs = itertools.combinations(ising_coupler[:-1], n)
        
        if n not in hobo_set.keys():
            hobo_set[n] = {'ind':[], 'w':[]}
        
        for indices in combs:
            if indices not in hobo_set[n]['ind']:
                hobo_set[n]['ind'].append(indices)
                hobo_set[n]['w'].append(multiplier*ising_coupler[-1])
            else:
                pos = hobo_set[n]['ind'].index(indices)
                hobo_set[n]['w'][pos] += multiplier*ising_coupler[-1]
                
        multiplier *= mult_fact
    
    dE = -ising_coupler[-1]
        
    return dE



def ising_to_hobo(ising_bonds):
    """
    Converts the Ising couplers in ising_bonds to hobo format.
    Returns the resultant hobo couplers and the energy difference dE.

    """
    
    hobo_bonds = []
    
    hobo_set = {}
    dE = 0
    
    for ising_coupler in ising_bonds:
        dE_local = add_coupler(ising_coupler, hobo_set)
        dE += dE_local
    
    for n in reversed(sorted(hobo_set.keys())):
        for ind, w in zip(hobo_set[n]['ind'], hobo_set[n]['w']):
            if abs(w) > 1.0e-10: # Check if the weights are non-zero
                hobo_bonds.append(ind + (w, ))

    return hobo_bonds, dE
        


def perform_gauge_transform(bonds, num_spins):
    """
    Performs a random gauge transformation on Ising systems.

    """
    elem_states = [-1, 1]
    
    C = [random.choice(elem_states) for i in range(num_spins)]

    gauged_bonds = []

    for coupler in bonds:
        gauged_bonds.append( coupler[:-1] + ( coupler[-1]*reduce(lambda i, j: C[i]*C[j], coupler[:-1]), ) )     

    return gauged_bonds

