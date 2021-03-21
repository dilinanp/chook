import numpy as np


def one_local_problem(N, convert_to_qubo=False):
    arr = [-1, 1]

    indices = np.random.choice(2, N)
    h = np.array([arr[i] for i in indices], dtype=int)

    E0 = -N

    if convert_to_qubo:
        c = -2*h
        E0 -= np.sum(h)

        bonds = list(zip(np.arange(N), c))
    else:
        bonds = list(zip(np.arange(N), h))

    return bonds, E0


class Problem:
    def __init__(self, bonds, num_spins, gs_energy, in_hobo_format=False):
        self.bonds = bonds
        self.num_spins = num_spins
        self.gs_energy = gs_energy 
        self.in_hobo_format = in_hobo_format

    def combine(self, prob):
        self.bonds.append( (-self.gs_energy, ) )
        prob.bonds.append( (-prob.gs_energy, ) )

        E = 0
        indices_list = []
        weights_list = []

        for bond1 in self.bonds:
            for bond2 in prob.bonds:
                weight = bond1[-1]*bond2[-1]
                indices = self.get_nonvanishing_spin_indices(bond1, bond2)
        
                if indices:
                    if indices in indices_list:
                        weights_list[indices_list.index(indices)] += weight
                    else:
                        indices_list.append(indices)
                        weights_list.append(weight)
                else:
                    E -= weight

        H = []

        for indices, weight in zip(indices_list, weights_list):
            if abs(weight) > 1.0e-10:
                H.append(tuple(indices) + (weight,))
                
        N = max(self.num_spins, prob.num_spins)


        return Problem(H, N, E, self.in_hobo_format)

    def get_nonvanishing_spin_indices(self, bond1, bond2):
        temp_indices = bond1[:-1]+bond2[:-1]
        
        if self.in_hobo_format:
            return list(set(temp_indices))
        else:
            existing_elems = set()
            repeating_elems = []

            for x in temp_indices:
                if x in existing_elems:
                    repeating_elems.append(x)
                existing_elems.add(x)

            nonvanishing_elems = existing_elems - set(repeating_elems)
            
            return list(nonvanishing_elems)




def build_klocal_problem(subproblems):
    result = subproblems[0]

    for subprob in subproblems[1:]:
        result = result.combine(subprob)

    return result

