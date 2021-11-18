import os
import sys

'''
    This code does the calibration via runwscleanauto.
    It needs the region file in the folder, along with the splitted measurement set.
    It runs after standalone_peel, but it requires a copy step where the splitted MS + boxfile
    are copied to a dedicated folder.
    Output: a lot of files, images, losoto-plots and merged.h5
'''

CHOUT = 12
STOP  =  6

peelnum = sys.argv[1]
callstring = f"python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/runwscleanLBautoR.py -b Dir{peelnum}.reg --startfromtgss --usemodeldataforsolints --usewgridder True --channelsout={CHOUT} --stop={STOP} --uvmin=60 --docircular --BLsmooth Dir{peelnum}.peel.ms"
os.system(callstring)
