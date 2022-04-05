import pandas as pd
import numpy as np
import os
import shutil
from subprocess import call, run
import sys
from glob import glob

def parse_params(fpath):
    """
    fpath : str (path to the file to parse)
    ------
    returns: {'total_porosity': value, 'absolute_permeability' : value, 'formation_factor': value} 
    
    """
    lines = []
    with open(fpath) as f:
        lines = f.readlines()
    for line in lines:
        if 'Total porosity' in line :
            total_porosity = float(line.split()[-1])
        if 'Absolute permeability' in line:
            absolute_permeability = float(line.split()[-2])
        if 'Formation factor' in line:
            formation_factor = float(line.split()[-1])
            
    return {'total_porosity': total_porosity, 'absolute_permeability' : absolute_permeability, 'formation_factor': formation_factor}



def generate_pnextract():
    path = '/home/sirius/gzprm2022/Sirius_2022/data/'
    new_path = '/home/sirius/gzprm2022/Sirius_2022/data/'
    dat_path = "/home/sirius/gzprm2022/Sirius_2022/info/software/pnflow-master_commit20210823/example/"
    mhd_files = {}
    raw_files = {}
    dirs_paths = []

    filenames = glob(path+"**/*.mhd", recursive=True)
    # for root, dirs, files in os.walk(path):
    #     for direc in dirs:
    #         if direc != 'data':
    #             new_path = os.path.join(root, direc)
    #             dirs_paths.append(new_path)
    #             files_list_mhd = []
    #             files_list_raw = []
    #             for r, d, files in os.walk(new_path):
    for file in filenames:
        # if file.endswith('.mhd'):
        name = file.replace('.mhd', '')
        name = name.split("/")[-1]
        with open(dat_path + 'input_pnflow.dat', 'r') as inp:
            llist = []
            for line in inp:
                if line.startswith('NETWORK'):
                    line = f'NETWORK  F {name};   // the base name for of the network file, without _link1.dat, _link2, _pore1'
                if line.startswith('TITLE'):
                    line = f'TITLE  {name};  // base name for the output files' + "\n"
                llist.append(line)

        with open(dat_path + 'input_pnflow.dat', 'w') as inp:
            for item in llist:
                inp.write(item)
        
            # inp.writelines(llist)
        file_dir = "/".join(file.split("/")[:-1])
        pnextract = os.system(f'cd {file_dir}; /home/sirius/gzprm2022/Sirius_2022/info/software/pnflow-master_commit20210823/bin/pnextract {name} > log.pnextract; /home/sirius/gzprm2022/Sirius_2022/info/software/pnflow-master_commit20210823/bin/pnflow input_pnflow.dat > log.pnflow')
        # files_list_mhd.append(file)
        sys.exit()

# params_dict = {}
# for root, dirs, files in os.walk(path):
#     for file in files:
#         if file.endswith('pnflow.prt') and len(file) < 15:
#             filename = file.replace('_pnflow.prt', '')
#             prt_path = os.path.join(root, file)
#             print(file)
#             params = parse_params(prt_path)
#             print(prt_path)
# #             prt_path = os.path.join(new_path, filename)
#             params_dict[prt_path] = params
#         mhd_files[direc] = files_list_mhd
#         raw_files[direc] = files_list_raw

if __name__ == "__main__":
    generate_pnextract()