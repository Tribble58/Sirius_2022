import sys  
# setting path
sys.path.insert(0, '/home/sirius/gzprm2022/Sirius_2022/src')

from src.DRA_utils_v2 import swapDirections

# поменять местами ось X и Z:
swapD('S5.mhd','S5_swapXZ.mhd','XZ')

# поменять местами ось X и Z:
swapD('S5.mhd','S5_swapXY.mhd','XY')
