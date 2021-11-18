import glob
import numpy as np
from regions import Regions
import sys

def sort_folders_glob(globlist):
    substrings = []
    for el in globlist:
        substrings.append(int(el.split('Dir')[1].split('.reg')[0]))
    substrings = np.array(substrings)
    return list(np.array(globlist)[np.argsort(substrings)])


root = sys.argv[1]
regfiles = sort_folders_glob(glob.glob(root+'/*reg'))
outname = '/'.join(regfiles[0].split('/')[:-1]) + '/combined.reg'
try:
    os.remove(outname)
except:
    pass

reglist = [Regions.read(regf,format='ds9')[0] for regf in regfiles]

newregs = Regions(reglist)
newregs.write(outname, format='ds9')
