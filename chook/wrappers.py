import numpy as np
import math
import sys

from chook.spin_utils import perform_gauge_transform, ising_to_hobo

import chook.planters.tile_planting as tp
import chook.planters.wishart_planting as wp
import chook.planters.fcl_generator as dcl 
import chook.planters.regular_xorsat as xor
import chook.planters.klocal_generator as klocal


def tile_planting_wrapper(length, dimension, tile_params, gauge_transform=True, convert_to_hobo=False):

    if dimension == 2:
        J = tp.generate_2D_problem(length, tile_params[0], tile_params[1], tile_params[2])
    else:
        J = tp.generate_3D_problem(length, tile_params[0], tile_params[1])

    i_indices, j_indices = np.nonzero(np.triu(J))

    J_ij = J[(i_indices, j_indices)]

    J_ij = np.negative(J_ij)

    E0 = np.sum(J_ij).astype(int)

    num_spins = length**dimension

    J_ij = J_ij.astype(int)

    bonds = list(zip(i_indices, j_indices, J_ij))

    if gauge_transform:
        bonds = perform_gauge_transform(bonds, num_spins)

    if convert_to_hobo:
        bonds, dE = ising_to_hobo(bonds)
        E0 += dE

    return bonds, E0



def wishart_planting_wrapper(length, M, discretize_couplers=False, gauge_transform=True, convert_to_hobo=False):

    J = wp.generate_problem(length, M, discretize_bonds=discretize_couplers)

    if discretize_couplers:
        i_indices, j_indices = np.nonzero(np.triu(J))
    else:
        i_indices, j_indices = np.triu_indices(length, k=1) 

    J_ij = J[(i_indices, j_indices)]

    J_ij = np.negative(J_ij)

    E0 = np.sum(J_ij)

    bonds = list(zip(i_indices, j_indices, J_ij))

    if gauge_transform:
        bonds = perform_gauge_transform(bonds, length)

    if convert_to_hobo:
        bonds, dE = ising_to_hobo(bonds)
        E0 += dE

    return bonds, E0



def dcl_wrapper(Lx, Ly, M, R, scaling, convert_to_hobo=False):
    q = dcl.lattice(Lx, Ly)

    maximum_attempts = 100000

    try:
        ising = dcl.gen_fcl_ising(q, num_loops=M, min_rho=-R, max_rho=R, max_attempts=maximum_attempts)
    except:
        print('\nCould not generate instance: Number of attempts has reached max_attempsts={}.'.format(maximum_attempts))
        print('Try decreasing alpha or increasing R.\n\n')
        sys.exit()

    bonds = dcl.get_dcl_chimera(ising, lambda_par=scaling)

    if convert_to_hobo:
        bonds, __ = ising_to_hobo(bonds)
    
    return bonds    



def xorsat_wrapper(k, N, convert_to_hobo=False):
    bonds, E0, num_solutions = xor.plant_regular_xorsat(k, N)

    if convert_to_hobo:
        bonds, dE = ising_to_hobo(bonds)
        E0 += dE

    return bonds, E0, num_solutions    



def klocal_wrapper(subproblem_types, subproblem_count, subproblem_params, convert_to_hobo=False):

    subproblems = []

    for subprob_type, subprob_count, subprob_params in zip(subproblem_types, subproblem_count, subproblem_params):

        if subprob_type=='RF':
            bonds, E0 = klocal.one_local_problem(subprob_params['N'], convert_to_hobo)        

            N = subprob_params['N'] 

        elif subprob_type=='TP':
            bonds, E0 = tile_planting_wrapper(subprob_params['length'], subprob_params['dimension'], subprob_params['tile_params'], subprob_params['gauge_transform'], convert_to_hobo)       
            N = subprob_params['length']**subprob_params['dimension']

        elif subprob_type=='WP':
            bonds, E0 = wishart_planting_wrapper(subprob_params['length'], subprob_params['M'], subprob_params['discretize_couplers'], subprob_params['gauge_transform'], convert_to_hobo)
            N = subprob_params['length']

        for n in range(subprob_count):
            subproblems.append( klocal.Problem(bonds, N, E0) )
    
    result = klocal.build_klocal_problem(subproblems)
    
    return result.bonds, result.gs_energy    

