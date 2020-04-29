import numpy as np


def find_num_solutions(A_pass, b_pass):
    """
        Using row reduction, determines the number of solutions that satisfy 
        the linear system of equations given by A_pass.X = b_pass mod 2.
        Returns zero if no solutions exist.
 
    """    
    A = np.copy(A_pass)
    b = np.copy(b_pass)

    M, N = A.shape

    h = 0
    k = 0

    while h<M and k<N:

        max_i = h

        for i in range(h, M):
            if A[i, k] == 1:
                max_i = i
                break
        
        if A[max_i, k] == 0:
            k += 1        
        else:
            if h != max_i:
                A[[h, max_i]] = A[[max_i, h]]
                b[[h, max_i]] = b[[max_i, h]]

            for u in range((h+1), M):
                flip_val = A[u, k]
                A[u] = ( A[u] + flip_val*A[h] ) % 2
                b[u] = ( b[u] + flip_val*b[h] ) % 2

            h += 1
            k += 1

    # Find rows with all zeros
    num_all_zeros_rows = 0

    solutions_exist = True

    for i in range(M):
        if not np.any(A[i]): # All-zero row encountered

            if b[i] != 0:
                solutions_exist = False
                break

            num_all_zeros_rows += 1

    if solutions_exist:
        rank = M - num_all_zeros_rows
        num_solutions = np.power(2, N-rank)     
    else:
        num_solutions = 0

    return num_solutions



def regular_xorsat(k, n):
    """
    Generates a k-regurlar k-XORSAT instance with n variables.

    """

    indices = np.zeros((n, k), dtype=int)

    while True:
        for i in range(k):
            indices[:, i] = np.random.permutation(n)

        if all(np.unique(row).size == k for row in indices):
            break

    b = np.random.choice(2, n)

    return indices, b



def plant_regular_xorsat(k, n):
    """
    Repetitively generates k-regular k-XORSAT problems until an instance
    with non-zero solutions is obtained.
    Returns the resultant XORSAT instance formulated as an
    Ising minimization problem.
    Also returns the ground state energy and ground state degenaracy.

    """
    
    while True:
        indices, b = regular_xorsat(k, n)

        A = np.zeros( (n, n), dtype=int )

        for i, row in enumerate(indices):
            A[i, row] = 1

        num_solutions = find_num_solutions(A, b)

        if num_solutions != 0:
            break

    gs_energy = -n

    bonds = []

    for i in range(n):
        bonds.append( tuple(indices[i]) + (-np.power(-1, b[i]), ) )

    return bonds, gs_energy, num_solutions

