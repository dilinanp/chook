# Chook #
#### Version 0.1.0 (Release date: May 19, 2020) ####

## License ##

This software is shared under an Apache License (attached). We ask you to
please cite the paper if you choose to use chook. Thank you. 

## Description ##

Chook is a comprehensive suite for generating binary optimization problems with planted solutions.
Chook currently supports the following problem types:
* Tile planting (2D/3D) 
* Wishart planting 
* Deceptive cluster loops (DCL)
* Equation planting (k-regular k-XORSAT)
* k-local planting

## Requirements ##

You will need Python 3.4 or later to run Chook.

## Installation ##

One can install Chook using the Python package manager [pip](https://pip.pypa.io/en/stable/):

```
pip install chook
```

Alternatively, Chook can be downloaded from GitHub and can be installed as follows:

```
cd <path_to_repo>
python setup.py install
```

## Usage ## 
```
chook problem_type config_file [-h] [-n num_instances] [-o output_format] [-f file_format] 
```
#### Required arguments ####
* **problem_type:**
Choose the type of problems to be generated. Allowed options: {`TP`, `WP`, `DCL`, `XORSAT`, `K_LOCAL`}
    * `TP`      : Tile planting (2D/3D) 
    * `WP`      : Wishart planting 
    * `DCL`     : Deceptive cluster loops (DCL)
    * `XORSAT`  : Equation planting (k-regular k-XORSAT)  
    * `K_LOCAL` : k-local planting
                     
* **config_file:**
Configuration file containing problem-type specific parameters. 
Refer the provided **params.in** file for an example.

#### Optional arguments ####
  * **-n num_instances**  
    The number of problem instances to be generated (default: 10)
  * **-o output_format**  
    Output format: ising/hobo (default: ising)
  * **-f file_format**    
    File format: txt/json (default: txt)
  * **-h, --help**        
    Show this help message and exit

The generated instances will be stored as text/JSON files in a subdirectory within the
current working directory. The ground state energies will be recorded in a separate text
file **gs_energies.txt**. For XORSAT problems, **gs_energies.txt** will also contain the
ground state degeneracies.

