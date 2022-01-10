import os
import glob,sys
import casacore.tables as pt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u


try:
    location = sys.argv[1]
except:
    location = './'

cwd = os.getcwd()
os.chdir(location)

#mslist = sorted(glob.glob('L813978_SAP000_SB033_uv.MS'))
mslist = sorted(glob.glob('L814006_SAP000_SB???_uv.MS'))
mslist = sorted(glob.glob('*_uv.MS'))


avgfreqstep = 8
avgtimestep = 4
targetsource = None #'3c196_4c' # None
#targetsource = '3c380'

c3c380= np.array([-1.44194739, 0.85078014])
c3c196= np.array([2.15374139,0.8415521])

for ms in mslist:
 mout = ms.split('.MS')[0] + '.avg.msdemix'
 print(mout)
 if not os.path.isdir(mout):
   t=pt.table(ms + '/FIELD')
   #print(t.colnames())
   adir = t.getcol('DELAY_DIR')[0][0][:]
   t.close()
   #print(adir)
   cdatta = SkyCoord(adir[0]*u.deg, adir[0]*u.deg, frame='icrs')


   if (cdatta.separation( SkyCoord(c3c380[0]*u.deg, c3c380[0]*u.deg, frame='icrs'))) < 0.2*u.deg:
     targetsource = '3c380'    
   if (cdatta.separation( SkyCoord(c3c196[0]*u.deg, c3c196[0]*u.deg, frame='icrs'))) <  0.2*u.deg:
     targetsource = '3c196_4c' 

   t=pt.table(ms)
   #print(t.colnames())
   if  np.nanmedian(t.getcol('WEIGHT_SPECTRUM')) == 1.0: 
     raw = True
     print('RAW correlator data')
   else:
     raw = False
   t.close()


   cmd = 'DPPP msin=' + ms + ' msout.storagemanager=dysco '

   instrument = mout + '/' + 'instrument' + ' '
   cmd += 'msout.writefullresflag=False msout='+ mout + ' '
  
   if raw:
     cmd += 'msin.autoweight=True '
     cmd += 'pf.type=preflagger '
     cmd += 'pf.chan=[0..nchan/32-1,31*nchan/32..nchan-1] '
     cmd += 'pf2.type=preflagger pf2.corrtype=auto '
     cmd += 'ao.strategy=/net/rijn/data2/rvweeren/software/prefactorAug29/rfistrategies/LBAdefaultnarrowband.rfis '
     cmd += 'ao.type=aoflagger ao.autocorr=False ao.keepstatistics=False '   

 
   if raw:
     cmd += 'steps=[pf,pf2,ao,demix] ' 
   else:
     cmd += 'steps=[demix] '

   cmd += 'demix.type=demix demix.baseline="[CR]S*&" '
   cmd += 'demix.corrtype=cross '
   cmd += 'demix.freqstep=' + str(avgfreqstep) + ' ' 
   cmd += 'demix.timestep='+ str(avgtimestep) + ' '
   cmd += 'demix.skymodel=demix.sourcedb demix.subtractsources=[CasA,CygA] '
   if targetsource != None:
     cmd += 'demix.targetsource=' + targetsource + ' '
   cmd += 'demix.instrumentmodel=' + instrument + ' '
   cmd += 'demix.demixtimestep=16 demix.demixfreqstep=64 '
   cmd += 'demix.ignoretarget=False'
 
   print(cmd)
   os.system(cmd)

os.chdir(cwd)
