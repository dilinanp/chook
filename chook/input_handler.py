import argparse
import configparser
import sys
import numpy as np
from collections import Counter

def positive_int(value):
    ivalue = int(value)

    if ivalue <= 0:
        raise argparse.ArgumentTypeError("invalid positive_int value: %s" % value)

    return ivalue


def parse_comm_line_params():
    
    parser = argparse.ArgumentParser(description='**************************************** CHOOK ****************************************\nA comprehensive suite for generating binary optimization problems with planted solutions', \
            epilog="The generated instances will be stored as text/JSON files in a subdirectory within the\n"
                   "current working directory. The ground state energies will be recorded in a separate text\n" 
                   "file 'gs_energies.txt'. For XORSAT problems, 'gs_energies.txt' will also contain the\n"
                   "ground state degeneracies.\n"
                   " ", \
            formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("problem_type", choices=['TP', 'WP', 'DCL', 'XORSAT', 'K_LOCAL'], metavar="problem_type", type=str.upper, 
                                help="Choose a problem type from the following options: \n"
                                     "  TP      : Tile planting \n"
                                     "  WP      : Wishart planting \n"
                                     "  DCL     : Deceptive cluster loops (DCL) \n"
                                     "  XORSAT  : Equation planting (k-regular k-XORSAT)\n"
                                     "  K_LOCAL : k-local planting \n"
                                     " ")

    parser.add_argument("config_file", help="Configuration file containing problem-type-specific parameters")

    parser.add_argument("-n", help="The number of problem instances to be generated (default: 10)", type=positive_int, dest="no_of_instances", metavar="num_instances", default=10)
    parser.add_argument("-o", help="Output format: ising/hobo (default: ising)", dest="output_format", metavar="output_format", 
                            choices=["ising", "hobo"], default="ising", type=str.lower)
    parser.add_argument("-f", help="File format: txt/json (default: txt)", dest="file_format", metavar="file_format", 
                            choices=["txt", "json"], default="txt", type=str.lower)


    args = parser.parse_args()


    params = {}

    params['num_instances']  = args.no_of_instances
    params['problem_type']  = args.problem_type.upper()
    params['config_file'] = args.config_file
    params['file_format'] = args.file_format 

    if args.output_format=='hobo':
        params['convert_to_hobo'] = True 
    else:
        params['convert_to_hobo'] = False

    return params



def read_config_file(params):

    config = configparser.ConfigParser()

    err_msg = '\nThe following error(s) occurred while parsing the configuration file {} :\n\n'.format(params['config_file'])

    try:
        with open(params['config_file']) as f:
            config.read_file(f)
    except IOError:
        print('\nUnable to read configuration file {} \n'.format(params['config_file']))
        sys.exit()
    except:
        err_msg += 'Invalid configuration file format for problem_type {}.\n'.format(params['problem_type']) 
        err_msg += 'Parameter initialization should be followed by the section header [{}].\n'.format(params['problem_type'])
        print(err_msg)
        sys.exit()

    # Begin parsing the configuration file

    # First, check if the section header exists

    try:
        if not config.has_section(params['problem_type']):
            raise KeyError
    except:
        err_msg += 'Invalid configuration file format for problem_type {}.\n'.format(params['problem_type']) 
        err_msg += 'Parameter initialization should be followed by the section header [{}].\n'.format(params['problem_type'])
        print(err_msg)
        sys.exit() 
            

    err_occurred = False

    params[ params['problem_type'] ] = {}


    if params['problem_type'] == 'TP':

        params['TP']['dimension'] = None

        try:
            params['TP']['dimension'] = int(config['TP']['dimension'])
            
            if ( params['TP']['dimension'] not in [2, 3]):
                raise ValueError
        except ValueError:
            err_msg += 'Invalid value for [TP]->dimension.\n' 
            err_msg += 'Only 2 and 3 dimensions are allowed.\n\n'
            err_occurred = True 
        except:
            err_msg += 'Unable to read required parameter [TP]->dimension.\n\n'
            err_occurred = True
        
        try:
            params['TP']['length'] = int(config['TP']['L'])
            
            if params['TP']['length']%2 != 0 or params['TP']['length']<4:
                raise ValueError
        except ValueError:
            err_msg += 'Invalid value for [TP]->L.\n'
            err_msg += 'Only even numbers greater than 2 are allowed.\n\n'
            err_occurred = True
        except:
            err_msg += 'Unable to read required parameter [TP]->dimension.\n\n'
            err_occurred = True
       

        if params['TP']['dimension'] == 2:
            try:
                params['TP']['tile_params'] = np.array([float(config['TP']['p1']), float(config['TP']['p2']), float(config['TP']['p3'])])

                if not ((params['TP']['tile_params'] >= 0.0).all() and (params['TP']['tile_params'] <= 1.0).all()):
                    raise ValueError
            
                if params['TP']['tile_params'].sum() > 1.0:
                    raise ValueError  
            except ValueError:
                err_msg += 'Invalid values for [TP]->p1, p2, and p3.\n'
                err_msg += 'Applicable conditions:  0 <= p1, p2, p3 <= 1  &  p1 + p2 + p3 <= 1.0 \n\n'
                err_occurred = True
            except:
                err_msg += 'Unable to read required parameter [TP]->p1, p2, and p3.\n\n'
                err_occurred = True


        if params['TP']['dimension'] == 3:
            try:
                params['TP']['tile_params'] = np.array([float(config['TP']['pF22']), float(config['TP']['pF42'])])

                if not ((params['TP']['tile_params'] >= 0.0).all() and (params['TP']['tile_params'] <= 1.0).all()):
                    raise ValueError
            
                if params['TP']['tile_params'].sum() > 1.0:
                    raise ValueError  
            except ValueError:
                err_msg += 'Invalid values for [TP]->pF22 and pF42.\n'
                err_msg += 'Applicable conditions:  0 <= pF22, pF42 <= 1  &  pF22 + pF42 <= 1.0 \n\n'
                err_occurred = True
            except:
                err_msg += 'Unable to read required parameter [TP]->pF22 and pF42.\n\n'
                err_occurred = True


		# Apply guage transform?
        try:
            if not config['TP']['gauge_transform']:
                raise KeyError

            params['TP']['gauge_transform'] = config['TP'].getboolean('gauge_transform')
        except ValueError:
            err_msg += 'Invalid response for [TP]->gauge_transform.\n'
            err_msg += 'Valid responses: yes/no \n\n' 
            err_occurred = True
        except:
            err_msg += 'Unable to read required parameter [TP]->gauge_transform.\n\n'
            err_occurred = True
		
        
    elif params['problem_type'] == 'WP':
        try:
            params['WP']['length'] = int(config['WP']['N'])
            
            if params['WP']['length'] <= 2:
                raise ValueError
        except ValueError:
            err_msg += 'Invalid value for [WP]->N.\n'
            err_msg += 'Must be an integer greater than 2.\n\n'
            err_occurred = True 
        except:
            err_msg += 'Unable to read required parameter [WP]->N.\n\n'
            err_occurred = True
       
        
        try:    
            params['WP']['alpha'] = float(config['WP']['alpha'])

            if not (params['WP']['alpha'] > 0):
                raise ValueError
        except ValueError:
            err_msg += 'Invalid value for [WP]->alpha.\n'
            err_msg += 'Applicable conditions: alpha > 0.0\n\n'
            err_occurred = True
        except:
            err_msg += 'Unable to read required parameter [WP]->alpha.\n\n'
            err_occurred = True


		# Discretize couplers?
        try:
            if not config['WP']['discretize_couplers']:
                raise KeyError

            params['WP']['discretize_couplers'] = config['WP'].getboolean('discretize_couplers')
        except ValueError:
            err_msg += 'Invalid response for [WP]->discretize_couplers.\n'
            err_msg += 'Valid responses: yes/no \n\n'
            err_occurred = True
        except:
            err_msg += 'Unable to read required parameter [WP]->discretize_couplers.\n\n'
            err_occurred = True


		# Apply guage transform?
        try:
            if not config['WP']['gauge_transform']:
                raise KeyError

            params['WP']['gauge_transform'] = config['WP'].getboolean('gauge_transform')
        except ValueError:
            err_msg += 'Invalid response for [WP]->gauge_transform.\n'
            err_msg += 'Valid responses: yes/no \n\n'
            err_occurred = True
        except:
            err_msg += 'Unable to read required parameter [WP]->gauge_transform.\n\n'
            err_occurred = True
        

        if not err_occurred:
            N = params['WP']['length']
            alpha_provided = params['WP']['alpha']
            M = int( round(alpha_provided*N) )

            if M == 0:
                M = 1
            
            alpha_actual = M/N 

            params['WP']['alpha'] = alpha_actual
            params['WP']['M'] = M
            
            print('\n[WP]->alpha :')
            print('Provided value = {}'.format(alpha_provided))
            print('Used value     = {:.3f}'.format(alpha_actual))

            percent_diff = abs((alpha_provided-alpha_actual))/alpha_provided*100.0

            if percent_diff > 10:
                print('Warning: Note that the difference between provided and used values differs by more than 10%')


    elif params['problem_type'] == 'DCL':

        try: 
            params['DCL']['Lx'] = int(config['DCL']['Lx'])

            if params['DCL']['Lx']<1 or params['DCL']['Lx']>16:
                raise ValueError
        except ValueError:
            err_msg += 'Invalid value for [DCL]->Lx.\n'
            err_msg += 'Must be an integer in the range [1, 16].\n\n'
            err_occurred = True
        except:
            err_msg += 'Unable to read required parameter [DCL]->Lx.\n\n'
            err_occurred = True

        try: 
            params['DCL']['Ly'] = int(config['DCL']['Ly'])

            if params['DCL']['Ly']<1 or params['DCL']['Ly']>16:
                raise ValueError
        except ValueError:
            err_msg += 'Invalid value for [DCL]->Ly.\n'
            err_msg += 'Must be an integer in the range [1, 16].\n\n'
            err_occurred = True
        except:
            err_msg += 'Unable to read required parameter [DCL]->Ly.\n\n'
            err_occurred = True


        try:    
            params['DCL']['alpha'] = float(config['DCL']['alpha'])

            if not (params['DCL']['alpha'] > 0):
                raise ValueError
        except ValueError:
            err_msg += 'Invalid value for [DCL]->alpha.\n'
            err_msg += 'Applicable conditions: alpha > 0.0\n\n'
            err_occurred = True
        except:
            err_msg += 'Unable to read required parameter [DCL]->alpha.\n\n'
            err_occurred = True


        try: 
            params['DCL']['R'] = int(config['DCL']['R'])

            if params['DCL']['R']<1:
                raise ValueError
        except ValueError:
            err_msg += 'Invalid value for [DCL]->R.\n'
            err_msg += 'Must be an integer greater than zero.\n\n'
            err_occurred = True
        except:
            err_msg += 'Unable to read required parameter [DCL]->R.\n\n'
            err_occurred = True


        try:    
            temp_lambda = float(config['DCL']['lambda'])

            params['DCL']['lambda'] = int(temp_lambda) if temp_lambda.is_integer() else temp_lambda

            if params['DCL']['lambda'] < 1:
                raise ValueError
        except ValueError:
            err_msg += 'Invalid value for [DCL]->lambda.\n'
            err_msg += 'Applicable conditions: lambda >= 1.0\n\n'
            err_occurred = True
        except:
            err_msg += 'Unable to read required parameter [DCL]->lambda.\n\n'
            err_occurred = True

        
        if not err_occurred:
            N = params['DCL']['Lx']*params['DCL']['Ly']

            alpha_provided = params['DCL']['alpha']

            M = int( round(alpha_provided*N) )

            if M == 0:
                M = 1
            
            alpha_actual = M/N 

            params['DCL']['alpha'] = alpha_actual
            params['DCL']['M'] = M
            
            print('\n[DCL]->alpha :')
            print('Provided value = {}'.format(alpha_provided))
            print('Used value     = {:.3f}'.format(alpha_actual))

            percent_diff = abs((alpha_provided-alpha_actual))/alpha_provided*100.0

            if percent_diff > 10:
                print('Warning: Note that the difference between provided and used values differs by more than 10%')


    elif params['problem_type'] == 'XORSAT':
        params['XORSAT']['k'] = None

        try:
            params['XORSAT']['k'] = int(config['XORSAT']['k'])
            
            if params['XORSAT']['k']<=1:
                raise ValueError
        except ValueError:
            err_msg += 'Invalid value for [XORSAT]]->k.\n'
            err_msg += 'Must be an integer greater than 1.\n\n'
            err_occurred = True
        except:
            err_msg += 'Unable to read required parameter [XORSAT]->k.\n\n'
            err_occurred = True

        try:
            params['XORSAT']['N'] = int(config['XORSAT']['N'])
           
            if params['XORSAT']['k'] is not None: 
                if params['XORSAT']['N'] < params['XORSAT']['k']:
                    raise ValueError
        except ValueError:
            err_msg += 'Invalid value for [XORSAT]->N.\n'
            err_msg += 'Should be an integer greater than or equal to [XORSAT]->k.\n\n'
            err_occurred = True 
        except:
            err_msg += 'Unable to read required parameter [XORSAT]->N.\n\n'
            err_occurred = True
     

    elif params['problem_type'] == 'K_LOCAL':

        try:
            params['K_LOCAL']['k'] = int(config['K_LOCAL']['k_max'])
            
            if params['K_LOCAL']['k']<=2:
                raise ValueError
        except ValueError:
            err_msg += 'Invalid value for [K_LOCAL]]->k_max.\n'
            err_msg += 'Must be an integer greater than 2.\n\n'
            err_occurred = True
        except:
            err_msg += 'Unable to read required parameter [K_LOCAL]->k_max.\n\n'
            err_occurred = True
 

        if not err_occurred:
            required_2_local_subprobs = params['K_LOCAL']['k'] // 2 
            required_1_local_subprobs = params['K_LOCAL']['k'] % 2

            try:
                params['K_LOCAL']['subproblem_list'] = [val.strip() for val in config['K_LOCAL']['subproblem_id_list'].split(',')]

                if len(params['K_LOCAL']['subproblem_list']) < (params['K_LOCAL']['k']+1)//2:
                    raise ValueError

                params['K_LOCAL']['subproblem_list'] = params['K_LOCAL']['subproblem_list'][:(params['K_LOCAL']['k']+1)//2]
            except ValueError:
                err_msg += 'Invalid values for [K_LOCAL]->subproblem_id_list.\n'
                err_msg += 'For k_max = {}, provide {} subproblem identifiers separated by commas.\n'.format(params['K_LOCAL']['k'], (params['K_LOCAL']['k']+1)//2)
                err_msg += '{} subproblems should be 2-local, and {} subproblems should be 1-local.\n\n'.format(required_2_local_subprobs, required_1_local_subprobs)

                err_occurred = True
            except:
                err_msg += 'Unable to read required parameter [K_LOCAL]->subproblem_id_list.\n\n'
                err_occurred = True


        if not err_occurred:
            for pId in params['K_LOCAL']['subproblem_list']:
                try:
                    if not config[pId]:
                        raise KeyError
                except:
                    err_msg += 'Missing subproblem specification in [K_LOCAL] for the subproblem identifier {}.\n\n'.format(pId)
                    err_occurred = True 


        # Check subproblem types

        if not err_occurred:
            params['K_LOCAL']['subproblem_ids']   = list( Counter(params['K_LOCAL']['subproblem_list']).keys() )
            params['K_LOCAL']['subproblem_count'] = list( Counter(params['K_LOCAL']['subproblem_list']).values() )

            params['K_LOCAL']['subproblem_types'] = []

            for pId in params['K_LOCAL']['subproblem_ids']:
                try:
                    params['K_LOCAL']['subproblem_types'].append( config[pId]['subproblem_type'].upper() )        

                    if not ( config[pId]['subproblem_type'] in ['TP', 'WP', 'RF'] ):
                        raise ValueError
                except ValueError:
                    err_msg += 'Invalid value for [{}]->subproblem_type.\n'.format(pId)
                    err_msg += 'Acceptable values: RF, TP, WP \n\n'
                    err_occurred = True 
                except:
                    err_msg += 'Unable to read required parameter [{}]->subproblem_type.\n\n'.format(pId)
                    err_occurred = True


        if not err_occurred:
            subproblem_localities = [ 1 if sub_type=='RF' else 2 for sub_type in params['K_LOCAL']['subproblem_types'] ]

            if np.sum(np.array(subproblem_localities)*np.array(params['K_LOCAL']['subproblem_count'])) != params['K_LOCAL']['k']:
                err_msg += 'Invalid subproblem specifications in [K_LOCAL].\n'
                err_msg += 'For k_max = {}, {} subproblems should be 2-local, and {} subproblems should be 1-local.\n\n'.format(params['K_LOCAL']['k'], required_2_local_subprobs, required_1_local_subprobs)
                err_occurred = True

        if not err_occurred:
            params['K_LOCAL']['subproblem_params'] = []

            params['K_LOCAL']['total_spins'] = 0

            system_sizes = []

            for subproblem_id, subproblem_type, subproblem_count in zip(params['K_LOCAL']['subproblem_ids'], params['K_LOCAL']['subproblem_types'], params['K_LOCAL']['subproblem_count']):
                subproblem_params = {}

                if subproblem_type == 'RF':

                    try:
                        subproblem_params['N'] = int(config[subproblem_id]['N'])
                        
                        if (subproblem_params['N'] < 1):
                            raise ValueError
                    except ValueError:
                        err_msg += 'Invalid value for [{}]->N.\n'.format(subproblem_id)
                        err_msg += 'Must be a positive integer.\n\n'
                        err_occurred = True 
                    except:
                        err_msg += 'Unable to read required parameter [{}]->N.\n\n'.format(subproblem_id)
                        err_occurred = True
            
                    if not err_occurred:
                        system_sizes.append(subproblem_params['N'])
                    
                elif subproblem_type == 'TP':
                    
                    subproblem_params['dimension'] = None

                    try:
                        subproblem_params['dimension'] = int(config[subproblem_id]['dimension'])
                        
                        if ( subproblem_params['dimension'] not in [2, 3]):
                            raise ValueError
                    except ValueError:
                        err_msg += 'Invalid value for [{}]->dimension.\n'.format(subproblem_id)
                        err_msg += 'Only 2 and 3 dimensions are allowed.\n\n'
                        err_occurred = True 
                    except:
                        err_msg += 'Unable to read required parameter [{}]->dimension.\n\n'.format(subproblem_id)
                        err_occurred = True
                    
                    try:
                        subproblem_params['length'] = int(config[subproblem_id]['L'])
                        
                        if subproblem_params['length']%2 != 0 or subproblem_params['length']<4:
                            raise ValueError
                    except ValueError:
                        err_msg += 'Invalid value for [{}]->L.\n'.format(subproblem_id)
                        err_msg += 'Only even numbers greater than 2 are allowed.\n\n'
                        err_occurred = True
                    except:
                        err_msg += 'Unable to read required parameter [{}]->length.\n\n'.format(subproblem_id)
                        err_occurred = True

                    if subproblem_params['dimension'] is not None:
                        if subproblem_params['dimension'] == 2:
                            try:
                                subproblem_params['tile_params'] = np.array([float(config[subproblem_id]['p1']), float(config[subproblem_id]['p2']), float(config[subproblem_id]['p3'])])

                                if not ((subproblem_params['tile_params'] >= 0.0).all() and (subproblem_params['tile_params'] <= 1.0).all()):
                                    raise ValueError
                            
                                if subproblem_params['tile_params'].sum() > 1.0:
                                    raise ValueError  
                            except ValueError:
                                err_msg += 'Invalid values for [{}]->p1, p2, and p3.\n'.format(subproblem_id)
                                err_msg += 'Applicable conditions:  0 <= p1, p2, p3 <= 1  &  p1 + p2 + p3 <= 1.0 \n\n'
                                err_occurred = True
                            except:
                                err_msg += 'Unable to read required parameter [{}]->p1, p2, and p3.\n\n'.format(subproblem_id)
                                err_occurred = True


                        if subproblem_params['dimension'] == 3:
                            try:
                                subproblem_params['tile_params'] = np.array([float(config[subproblem_id]['pF22']), float(config[subproblem_id]['pF42'])])

                                if not ((subproblem_params['tile_params'] >= 0.0).all() and (subproblem_params['tile_params'] <= 1.0).all()):
                                    raise ValueError
                            
                                if subproblem_params['tile_params'].sum() > 1.0:
                                    raise ValueError  
                            except ValueError:
                                err_msg += 'Invalid values for [{}]->pF22 and pF42.\n'.format(subproblem_id)
                                err_msg += 'Applicable conditions:  0 <= pF22, pF42 <= 1  &  pF22 + pF42 <= 1.0 \n\n'
                                err_occurred = True
                            except:
                                err_msg += 'Unable to read required parameter [{}]->pF22 and pF42.\n\n'.format(subproblem_id)
                                err_occurred = True

                    try:
                        if not config[subproblem_id]['gauge_transform']:
                            raise KeyError

                        subproblem_params['gauge_transform'] = config[subproblem_id].getboolean('gauge_transform')

                    except ValueError:
                        err_msg += 'Invalid response for [{}]->gauge_transform.\n'.format(subproblem_id) 
                        err_msg += 'Valid responses: yes/no \n\n'
                        err_occurred = True
                    except:
                        err_msg += 'Unable to read required parameter [{}]->gauge_transform.\n\n'.format(subproblem_id)
                        err_occurred = True

                    if not err_occurred:
                        system_sizes.append(subproblem_params['length']**subproblem_params['dimension'])
                    

                elif subproblem_type == 'WP':

                    try:
                        subproblem_params['length'] = int(config[subproblem_id]['N'])
                        
                        if subproblem_params['length'] <= 2:
                            raise ValueError
                    except ValueError:
                        err_msg += 'Invalid value for [{}]->N.\n'.format(subproblem_id)
                        err_msg += 'Must be an integer greater than 2.\n\n'
                        err_occurred = True 
                    except:
                        err_msg += 'Unable to read required parameter [{}]->N.\n\n'.format(subproblem_id)
                        err_occurred = True
                   
                    
                    try:    
                        subproblem_params['alpha'] = float(config[subproblem_id]['alpha'])

                        if not (subproblem_params['alpha'] > 0):
                            raise ValueError
                    except ValueError:
                        err_msg += 'Invalid value for [{}]->alpha.\n'.format(subproblem_id)
                        err_msg += 'Applicable conditions: alpha > 0.0\n\n'
                        err_occurred = True
                    except:
                        err_msg += 'Unable to read required parameter [{}]->alpha.\n\n'.format(subproblem_id)
                        err_occurred = True


                    # Discretize couplers?
                    try:
                        if not config[subproblem_id]['discretize_couplers']:
                            raise KeyError

                        subproblem_params['discretize_couplers'] = config[subproblem_id].getboolean('discretize_couplers')
                    except ValueError:
                        err_msg += 'Invalid response for [{}]->discretize_couplers.\n'.format(subproblem_id)
                        err_msg += 'Valid responses: yes/no \n\n'
                        err_occurred = True
                    except:
                        err_msg += 'Unable to read required parameter [{}]->discretize_couplers.\n\n'.format(subproblem_id)
                        err_occurred = True


                    # Apply guage transform?
                    try:
                        if not config[subproblem_id]['gauge_transform']:
                            raise KeyError

                        subproblem_params['gauge_transform'] = config[subproblem_id].getboolean('gauge_transform')
                    except ValueError:
                        err_msg += 'Invalid response for [{}]->gauge_transform.\n'.format(subproblem_id)
                        err_msg += 'Valid responses: yes/no \n\n'
                        err_occurred = True
                    except:
                        err_msg += 'Unable to read required parameter [{}]->gauge_transform.\n\n'.format(subproblem_id)
                        err_occurred = True
                    

                    if not err_occurred:
                        N = subproblem_params['length']
                        alpha_provided = subproblem_params['alpha']
                        M = int( round(alpha_provided*N) )

                        if M == 0:
                            M = 1
                        
                        alpha_actual = M/N 

                        subproblem_params['alpha'] = alpha_actual
                        subproblem_params['M'] = M
                        
                        print('\n[{}]->alpha :'.format(subproblem_id))
                        print('Provided value = {}'.format(alpha_provided))
                        print('Used value     = {:.3f}'.format(alpha_actual))

                        percent_diff = abs((alpha_provided-alpha_actual))/alpha_provided*100.0

                        if percent_diff > 10:
                            print('Warning: Note that the difference between provided and used values differs by more than 10%')

                    if not err_occurred:
                        system_sizes.append(subproblem_params['length'])

                if not err_occurred:
                    params['K_LOCAL']['subproblem_params'].append(subproblem_params)


        if not err_occurred:
            params['K_LOCAL']['total_spins'] = np.sum( np.array(system_sizes)*np.array(params['K_LOCAL']['subproblem_count']) )


    else: # Invalid problem type
        print('An unexpected error occurred.')
        sys.exit()


    if err_occurred:
        print(err_msg)
        sys.exit()


    return params



def read_input():
    params = parse_comm_line_params()

    read_config_file(params)

    return params
