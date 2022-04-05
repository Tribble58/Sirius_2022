import DRA_utils_v2 as dra
from glob import glob
import yaml
from yaml.loader import SafeLoader
import sys

# Open the file and load the file
with open('config/cfg.yaml') as f:
    data = yaml.load(f, Loader=SafeLoader)

files_dir = data["data"]["data_path"]
filenames = glob(files_dir+"/**/*.mhd", recursive=True)

for filename in filenames:
	filename_list = filename.split("/")
	file_dir = "/".join(filename_list[:-1])
	name = filename_list[-1][:-4]
	# поенять местами ось X и Z:
	dra.swapDirections(filename, file_dir+'/'+ name +"_swapXZ" + '.mhd',"XZ")

	# поенять местами ось X и Z:
	dra.swapDirections(filename, file_dir+'/'+ name +"_swapXY" + '.mhd','XY')
