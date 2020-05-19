from  chook.wrappers import *
from chook.file_handler import *
from chook.input_handler import read_input
import sys

def run():
    params = read_input()

    dir_name = create_dir(params)
   

    gs_info = [] # Include file_name, GS energy, and GS degeneracy (if available)


    for iteration in range(params['num_instances']):

        file_name = dir_name + '_inst_' + str(iteration+1) + ('.txt' if params['file_format']=='txt' else '.json')

        if params['problem_type'] == 'TP':

            bonds, E0 =  tile_planting_wrapper(params['TP']['length'], params['TP']['dimension'], params['TP']['tile_params'], params['TP']['gauge_transform'], params['convert_to_hobo'])

            gs_info.append( (file_name, E0) )

        elif params['problem_type'] == 'WP':

            bonds, E0 = wishart_planting_wrapper(params['WP']['length'], params['WP']['M'], params['WP']['discretize_couplers'], params['WP']['gauge_transform'], params['convert_to_hobo'])

            gs_info.append( (file_name, E0) )

        elif params['problem_type'] == 'DCL':

            bonds = dcl_wrapper(params['DCL']['Lx'], params['DCL']['Ly'], params['DCL']['M'], params['DCL']['R'], params['DCL']['lambda'], params['convert_to_hobo'])

        elif params['problem_type'] == 'XORSAT':

            bonds, E0, num_solutions = xorsat_wrapper(params['XORSAT']['k'], params['XORSAT']['N'], params['convert_to_hobo'])

            gs_info.append( (file_name, E0, num_solutions) )

        elif params['problem_type'] == 'K_LOCAL':

            bonds, E0 = klocal_wrapper(params['K_LOCAL']['subproblem_types'], params['K_LOCAL']['subproblem_count'], params['K_LOCAL']['subproblem_params'], params['convert_to_hobo'])

            gs_info.append( (file_name, E0) )
        else:
            print('An unexpected error occurred.')
            sys.exit()
    

        write_instance(dir_name, file_name, bonds, params)


    if params['problem_type'] in ['TP', 'WP', 'XORSAT', 'K_LOCAL']:
        write_gs_info(dir_name, gs_info)

    print('\nSuccessfully generated', params['num_instances'], 'instances of problem type', params['problem_type'], "in directory '%s'." % dir_name, '\n');


if __name__ == '__main__':
    run()
