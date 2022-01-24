import os, sys
import numpy as np
import glob
from astropy.io import fits
from astropy.wcs import WCS
import pyrap.tables as pt
import os.path
import pyregion
import argparse
import time

'''
    Standalone version of createPeelset.py
    This version is meant to be run inline - that is, it requires its own MS+all models in a separate folder
    Also, the region file needs to be present in the folder.

    The entry point is the folder containing this code, DirXX.reg, LOW_2-????-model.fits, something.MS
    The exit point is that, some rabble, and DirXX.peel.ms, which should then be copied to another folder for a good overview.
    Maybe we need to copy it to the parent folder?

    Usage:
    standalone_peel.py [DIRECTION_NAME]
'''

#https://singularity.lbl.gov/archive/docs/v2-2/docs-exec


# Not used if do_singularity = False
singularity = 'singularity exec -B /tmp,/dev/shm,/disks/paradata,/data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/botlek,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/disks/ftphome,/net/krommerijn,/net/voorrijn,/net/achterrijn,/net/tussenrijn,/net/ouderijn,/net/nieuwerijn,/net/lofar8,/net/lofar9 /net/lofar1/data1/sweijen/software/LOFAR/singularity/lofar_sksp_fedora27_ddf.sif '

def getimsize(fitsimage):
    hdu=fits.open(fitsimage)
    naxis=hdu[0].header['NAXIS1']
    hdu.close()
    return naxis

def add_dummyms(msfiles):
    '''
    Add dummy ms to create a regular freuqency grid when doing a concat with DPPP
    '''
    if len(msfiles) == 1:
      return msfiles
    keyname = 'REF_FREQUENCY'
    freqaxis = []
    newmslist  = []

    # Check for wrong REF_FREQUENCY which happens after a DPPP split in frequency
    for ms in msfiles:        
        t = pt.table(ms + '/SPECTRAL_WINDOW', readonly=True)
        freq = t.getcol('REF_FREQUENCY')[0]
        t.close()
        freqaxis.append(freq)
    freqaxis = np.sort( np.array(freqaxis))
    minfreqspacing = np.min(np.diff(freqaxis))
    if minfreqspacing == 0.0:
       keyname = 'CHAN_FREQ' 
    
    
    freqaxis = [] 
    for ms in msfiles:        
        t = pt.table(ms + '/SPECTRAL_WINDOW', readonly=True)
        if keyname == 'CHAN_FREQ':
          freq = t.getcol(keyname)[0][0]
        else:
          freq = t.getcol(keyname)[0]  
        t.close()
        freqaxis.append(freq)
    
    # put everything in order of increasing frequency
    freqaxis = np.array(freqaxis)
    idx = np.argsort(freqaxis)
    
    freqaxis = freqaxis[np.array(tuple(idx))]
    sortedmslist = list( msfiles[i] for i in idx )
    freqspacing = np.diff(freqaxis)
    minfreqspacing = np.min(np.diff(freqaxis))
 
    # insert dummies in the ms list if needed
    count = 0
    newmslist.append(sortedmslist[0]) # always start with the first ms the list
    for msnumber, ms in enumerate(sortedmslist[1::]): 
      if int(round(freqspacing[msnumber]/minfreqspacing)) > 1:
        ndummy = int(round(freqspacing[msnumber]/minfreqspacing)) - 1
 
        for dummy in range(ndummy):
          newmslist.append('dummy' + str(count) + '.ms')
          print('Added dummy:', 'dummy' + str(count) + '.ms') 
          count = count + 1
      newmslist.append(ms)
       
    print('Updated ms list with dummies inserted to create a regular frequency grid')
    print(newmslist) 
    return newmslist

def make_ms_datacol(ms, colname, do_singularity=False):
  t = pt.table(ms)
  if colname not in t.colnames():
    t.close()  
    cmd = 'DPPP msin=' + ms + ' msout=. msout.storagemanager=dysco '
    cmd+= 'msout.datacolumn=' + colname + ' steps=[]'
    print(cmd)
    if do_singularity:
      os.system(singularity + cmd)
    else:
      os.system(cmd)
  else:
    t.close()

def getregionboxcenter(regionfile):
    """
    Extract box center of a DS9 box region. 
    Input is regionfile Return NDPPP compatible string for phasecenter shifting
    """
    r = pyregion.open(regionfile)
    
    if len(r[:]) > 1:
      print('Only one region can be specified, your file contains', len(r[:]))
      sys.exit() 
    
    if r[0].name != 'box':
      print ('Only box region supported')
      sys.exit()
    
    ra  = r[0].coord_list[0]
    dec = r[0].coord_list[1]
    boxsizex = r[0].coord_list[2]
    boxsizey = r[0].coord_list[3]
    angle = r[0].coord_list[4]
    if boxsizex != boxsizey:
      print('Only a square box region supported, you have these sizes:', boxsizex, boxsizey)
      sys.exit()
    if np.abs(angle) > 1:
      print('Only normally oriented sqaure boxes are supported, your region is oriented under angle:', angle)
      sys.exit()   
    
    regioncenter =  ('{:12.8f}'.format(ra) + 'deg,' + '{:12.8f}'.format(dec) + 'deg').replace(' ', '')
    return regioncenter


def flatten(f):
    """ Flatten a fits file so that it becomes a 2D image. Return new header and data """

    naxis=f[0].header['NAXIS']
    if naxis<2:
        raise RadioError('Can\'t make map from this')
    if naxis==2:
        return fits.PrimaryHDU(header=f[0].header,data=f[0].data)

    w = WCS(f[0].header)
    wn=WCS(naxis=2)
    
    wn.wcs.crpix[0]=w.wcs.crpix[0]
    wn.wcs.crpix[1]=w.wcs.crpix[1]
    wn.wcs.cdelt=w.wcs.cdelt[0:2]
    wn.wcs.crval=w.wcs.crval[0:2]
    wn.wcs.ctype[0]=w.wcs.ctype[0]
    wn.wcs.ctype[1]=w.wcs.ctype[1]
    
    header = wn.to_header()
    header["NAXIS"]=2
    copy=('EQUINOX','EPOCH','BMAJ', 'BMIN', 'BPA', 'RESTFRQ', 'TELESCOP', 'OBSERVER')
    for k in copy:
        r=f[0].header.get(k)
        if r is not None:
            header[k]=r

    slice=[]
    for i in range(naxis,0,-1):
        if i<=2:
            slice.append(np.s_[:],)
        else:
            slice.append(0)
        
    hdu = fits.PrimaryHDU(header=header,data=f[0].data[tuple(slice)])
    return hdu

def mask_region(infilename,ds9region,outfilename):

    hdu=fits.open(infilename)
    hduflat = flatten(hdu)
    map=hdu[0].data

    r = pyregion.open(ds9region)
    manualmask = r.get_mask(hdu=hduflat)
    rmsval = np.mean(hdu[0].data[0][0][np.where(manualmask == True)])
    hdu[0].data[0][0][np.where(manualmask == True)] = 0.0
    hdu.writeto(outfilename,overwrite=True)

    return

do_singularity = False

peelregions = ['Dir1','Dir2','Dir3_larger','Dir4','Dir5','Dir6','Dir7','Dir8', 'Dir9','Dir10','Dir11']
peelregions = [f'Dir{i}' for i in range(19)]

ms          = 'corrected_L828012_concat.ms' 
globstr     = 'image_000-????-model*.fits' 
channelsout = 8
timestepavg = 4 # go to 16s
freqstepavg = 4 # go to 2 ch/sb
colname     = 'DATA_SUB'

# for name in peelregions:
ms          = sys.argv[1]
name        = sys.argv[2]
serialnum   = sys.argv[3]
boxfile     = name + '.reg'
msout       = name + '.' + serialnum + '.peel.ms'
imsize = getimsize(glob.glob(globstr)[0])
scale       = 8.0 #asec
r = pyregion.open(boxfile)
if len(r[:]) > 1:
    composite = True
else:
    print (boxfile)
    phasecenter = '[' + getregionboxcenter(boxfile) + ']'
    print (phasecenter)



if True:
    for model in glob.glob(globstr):
        modelout = 'restfield' + model
        print (model, modelout)
        mask_region(model,boxfile,modelout)
    for model in glob.glob(globstr):
        modelout = 'restfield' + model
        #print(modelout)
        if '-model-pb.fits' in modelout:
          cmdcp = 'cp ' + modelout + ' ' + modelout.replace("-pb.fits", ".fits")
          print ('WARNING, work around: ', cmdcp)
          os.system(cmdcp)
    #sys.exit()
        
time.sleep(2)


if True: # no IDG
    cmd = 'wsclean -size ' + str(imsize) + ' ' + str(imsize) + ' -reorder -parallel-reordering 4 -use-wgridder '
    cmd+= '-channels-out '+ str(channelsout)+ ' -padding 1.8 -pol i -name ' 
    cmd+= 'restfield' + glob.glob(globstr)[0].split('-')[0] + ' '
    cmd+= f'-scale {scale}arcsec  -predict ' + ms
    print(cmd)
    os.system(cmd)
    time.sleep(2)

if False: # with IDG
    #cmd = '/net/achterrijn/data2/sweijen/software/2020_11_09/wsclean/bin/wsclean '
    cmd = 'wsclean '
    cmd+= '-size '+ str(imsize) + ' ' + str(imsize) + ' -reorder -parallel-reordering 4 '
    cmd+= '-channels-out '+ str(channelsout) + ' -padding 1.2 -pol i -name '
    cmd+= 'restfield' + globstr.split('-')[0] + ' '
    cmd+= '-use-idg -idg-mode hybrid  -aterm-config aconfig.txt '
    cmd+= f'-scale {scale}arcsec -predict ' + ms
    print(cmd)
    os.system(cmd)
    time.sleep(2)
    #sys.exit()



make_ms_datacol(ms, colname, do_singularity=do_singularity)
if False: # tmp commented because of singulairty/python issues
  t = pt.table(ms)
  if 'CORRECTED_DATA' in t.colnames():
    t.close()
    print('Using CORRECTED_DATA-MODEL_DATA')
    if do_singularity:
      os.system(singularity + "taql 'update " + ms + " set DATA_SUB=CORRECTED_DATA-MODEL_DATA'")   
    else:    
      os.system("taql 'update " + ms + " set DATA_SUB=CORRECTED_DATA-MODEL_DATA'")
  else:
    t.close() 
    print('Using DATA-MODEL_DATA')
    if do_singularity:
      os.system(singularity + "taql 'update " + ms + " set DATA_SUB=DATA-MODEL_DATA'")  
    else:    
      os.system("taql 'update " + ms + " set DATA_SUB=DATA-MODEL_DATA'")  

# tmp workaround for singulairty/python issues
if do_singularity:
  os.system(singularity + "taql 'update " + ms + " set DATA_SUB=DATA-MODEL_DATA'")  
else:    
  os.system("taql 'update " + ms + " set DATA_SUB=DATA-MODEL_DATA'")  

if True:    
    cmd =  'DPPP msin="' + str(ms) + '" msout.writefullresflag=False '
    cmd += 'steps=[ps,average] '
    cmd += 'ps.type=phaseshift ps.phasecenter=' + phasecenter + ' '
    cmd += 'average.timestep=' + str(timestepavg) + ' average.freqstep=' + str(freqstepavg) + ' '   
    cmd += 'msin.weightcolumn=WEIGHT_SPECTRUM msout.storagemanager=dysco '
    cmd += 'msout=' + msout + ' '
    cmd += 'msin.datacolumn=%s '%colname
    print (cmd)
    if do_singularity:
      os.system(singularity + cmd)
    else:    
      os.system(cmd)
    #sys.exit()
# STEP 1 masks model images

# STEP 2 run wsclean to predict

# STEP 3 subtract columns

# STEP 4 phase-shift and average
