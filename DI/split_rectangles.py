import regions
import numpy as np
import sys
import os

regionfile = sys.argv[1]
rootfolder = regionfile.rstrip('.reg')
os.mkdir(rootfolder)

list_of_regions = regions.read_ds9(regionfile)
for i, rect in enumerate(list_of_regions):
    name = f'{rootfolder}/Dir{i}.reg'
    regions.write_ds9([rect],name)
