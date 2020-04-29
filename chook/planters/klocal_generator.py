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

    def __init__(self, bonds, num_spins, gs_energy):
        self.bonds = bonds
        self.num_spins = num_spins
        self.gs_energy = gs_energy 

    def combine(self, prob):
        H = []

        for bond1 in self.bonds:
            for bond2 in prob.bonds:
                w = bond1[-1]*bond2[-1]

                if abs(w) > 1.0e-10:
                    H.append( bond1[:-1] + tuple(val+self.num_spins for val in bond2[:-1]) + (w, ) )

        for bond1 in self.bonds:
            w = -prob.gs_energy*bond1[-1]

            if abs(w) > 1.0e-10:
                H.append( bond1[:-1] + (w, ) )

        for bond2 in prob.bonds:
            w = -self.gs_energy*bond2[-1]

            if abs(w) > 1.0e-10:
                H.append( tuple(val+self.num_spins for val in bond2[:-1]) + (w, ) )

        E = -self.gs_energy*prob.gs_energy
                 
        N = self.num_spins + prob.num_spins

        return Problem(H, N, E)


def build_klocal_problem(subproblems):
    result = subproblems[0]

    for subprob in subproblems[1:]:
        result = result.combine(subprob)

    return result

