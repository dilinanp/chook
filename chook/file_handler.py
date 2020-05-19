import numpy as np
import os
import sys
import csv
import json
import sys

class NumpyToJsonEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NumpyToJsonEncoder, self).default(obj)


def create_dir(params):

    if  params['problem_type'] == 'TP':

        if params['TP']['dimension'] == 2:
            dir_name = 'tile_planting_{}D_L_{}_p1_{}_p2_{}_p3_{}'.format( params['TP']['dimension'], params['TP']['length'], params['TP']['tile_params'][0], \
                                                                            params['TP']['tile_params'][1], params['TP']['tile_params'][2] )
        else:
            dir_name = 'tile_planting_{}D_L_{}_p2FP_{}_p4FP_{}'.format( params['TP']['dimension'], params['TP']['length'], params['TP']['tile_params'][0], \
                                                                        params['TP']['tile_params'][1] )

    elif  params['problem_type'] == 'WP':
        dir_name = 'wishart_planting_N_{}_alpha_{:.2f}'.format( params['WP']['length'], params['WP']['alpha'] )

    elif  params['problem_type'] == 'DCL':
        dir_name = 'dcl_Lx_{}_Ly_{}_alpha_{:.2f}_R_{}_lambda_{}'.format( params['DCL']['Lx'], params['DCL']['Ly'], params['DCL']['alpha'], params['DCL']['R'], \
                                                                        params['DCL']['lambda'] )

    elif params['problem_type'] == 'XORSAT':
        dir_name = 'regular_xorsat_k_{}_N_{}'.format( params['XORSAT']['k'], params['XORSAT']['N'])

    elif params['problem_type'] == 'K_LOCAL':
        dir_name = 'combined_{}_local_N_{}'.format( params['K_LOCAL']['k'], params['K_LOCAL']['total_spins'] )
    else:
        print('An unexpected error occurred.')
        sys.exit()

    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    return dir_name



def write_instance(dir_name, file_name, bonds, params):

    if params['file_format'] == 'txt': 
        write_text_file(os.path.join(dir_name, file_name), bonds)
    else:
        # Write output in json format

        output_format = 'hobo' if params['convert_to_hobo'] else 'ising'

        output = { "problem": {"type": output_format, "terms": []} }
        
        for coupler in bonds:
            output['problem']['terms'].append( {"w":coupler[-1], "ids":coupler[:-1] } )
        
        json_output = json.dumps(output, cls=NumpyToJsonEncoder)

        with open(os.path.join(dir_name, file_name), "w") as the_file:
            the_file.write(json_output)



def write_gs_info(dir_name, gs_info):

    write_text_file(os.path.join(dir_name, 'gs_energies.txt'), gs_info)



def write_text_file(path_to_file, records):

    with open(path_to_file, "w") as the_file:
        csv.register_dialect("custom", delimiter="\t", skipinitialspace=True)
        writer = csv.writer(the_file, dialect="custom")

        for row in records:
            writer.writerow(['{:.7e}'.format(entry) if isinstance(entry, float) else entry for entry in row])

