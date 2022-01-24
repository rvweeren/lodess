from astropy import units as u
import copy
import os
import gc
import glob
from astropy.io import fits
from astropy.wcs import WCS
from losoto.lib_operations import reorderAxes
from matplotlib.pyplot import figure, show
from scipy import ndimage
from scipy.interpolate import Rbf
from scipy.signal import medfilt
import matplotlib.pyplot as plt
import argparse
import casacore.tables as ct
import losoto.h5parm as h5parm
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import Angle
import sys
import scipy.ndimage as ndimage
import multiprocessing as mp


def plot_screen(img, prefix='', title='', suffix='', wcs=None):
    fig = figure(figsize=(12,12))
    if wcs is not None:
        ax = fig.add_subplot(111, projection=wcs)
    else:
        ax = fig.add_subplot(111, projection=wcs)
    
    ax.set_title(title)
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')
    i = ax.imshow(img, origin='lower', vmin=-0.1, vmax=0.1)
    cbar = fig.colorbar(i, fraction=0.046, pad=0.04)
    cbar.set_label('dTEC')
    fig.savefig(prefix + 'tecscreen_' + suffix + '.svg')
    plt.close(fig)
    del fig

def phaseref(phase, refstation = 0):
   print('Using reference station', refstation) 
   phasetmp =  np.copy(phase[:,:,refstation,:,:])
   for ant_id in range(len(phase[0,0,:,0,0])):
      #print (ant_id, phase[:,:,ant_id,:,:].shape,phasetmp.shape) 
      phase[:,:,ant_id,:,:] = phase[:,:,ant_id,:,:] - phasetmp
   return phase   

parser = argparse.ArgumentParser(description='Make a FITS IDG diaginal screen from a multi-dir H5 file using RBF interpolation, generates a FITS screen file. Important parameter to consider is boxwidth since it has to match the region covered by the WSClean IDD image')
#imaging options
parser.add_argument('--size', help='Size of the screen in pixels, default=64', default=64, type=int)
parser.add_argument('--boxwidth', help='Size of the screen in degrees, default=2.5', default=2.5, type=float)
parser.add_argument('--includeamps', help='Include the amplitude000 solutions (True/False, default=True)', type=eval, choices=[True, False], default='True')
parser.add_argument('--smoothamps', help='Create smooth amplitude screens, recommended when including amplitudes ', action='store_true')
parser.add_argument('--padding', help='Padd screen with amplitude=1 and phase=0 (True/False, default=False)', type=eval, choices=[True, False], default='False')
parser.add_argument('--stepsizepadding', help='Add padding direction every stepsizepadding pixels, default=4', default=4, type=int)
parser.add_argument('--H5file', help='H5 file', type=str, required=True)
parser.add_argument('--FITSscreen', help='Output FITS screen file', type=str, default='gainscreen_rbf.fits')
parser.add_argument('--freqdownsamplefactor', help='Downsample freq axis with this factor, needs to be a multiple on the input freq length', type=int, default=None)
parser.add_argument('--plotsfortesting', help='Only for testing.....', action='store_true')
parser.add_argument('--ms', help='Measurement set to be imaged with IDG, phasecenter location is taken from the ms', type=str,required=True)
parser.add_argument('--ncpu', help='Amount of subprocesses that are to be spawned for computation of gainscreens. This is a memory limited step, and non-unity values of ncpu might slow down this process!', type=int, default=1)
parser.add_argument('--timeblocks', help='Break up the final gainscreen in timeblocks, which will decrease the memory footprint',default=1,type=int)


if not args['includeamps']:
   if args['smoothamps']:
      print('Cannot use  --includeamps=False and --smoothamps, they are mutually exclusive...')
      sys.exit()

args = vars(parser.parse_args())
fitsfilename = args['FITSscreen']
fitsfilename_prefix = fitsfilename.split('.fits')[0]


ms = args['ms']
h5p =  args['H5file'] 
size        = args['size']   # 128 # update, not too large because you quickly run out of RAM
boxwidth    = args['boxwidth']   #1.5 # units is degr
includeamps = args['includeamps']   #True



t = ct.taql('SELECT NAME FROM {ms:s}::ANTENNA'.format(ms=ms))
names = t.getcol('NAME')
t.close()
t = ct.taql('SELECT REFERENCE_DIR FROM {ms:s}::FIELD'.format(ms=ms))
phasecenter = t.getcol('REFERENCE_DIR').squeeze()
print(phasecenter)
cdelt = np.float(boxwidth)*1.5/np.float(size)
cdelt1 = -1.*cdelt
cdelt2 = cdelt

if np.rad2deg(phasecenter[0]) < 0:
   phasecenter[0] = phasecenter[0] + (2.*np.pi)       

t.close()
# Time is stored as MJD, convert from seconds to days here as that's what FITS wants.

h5 = h5parm.h5parm(h5p)
ss = h5.getSolset('sol000')

st = ss.getSoltab('phase000')

time = st.getAxisValues('time')
stime = time[0] 
dtime = time[1] - time[0]
freq = st.getAxisValues('freq')


extension_needed = False
if args['freqdownsamplefactor'] != None:
  if (len(freq)) % (args['freqdownsamplefactor']) != 0:
    print('freqdownsamplefactor is not a multiple of',len(freq))
    freq = freq[::args['freqdownsamplefactor']]
    freq2 = np.zeros(shape=np.shape(freq)[0]+1)
    freq2[:-1] = freq
    freq2[-1] = freq2[-2] + (freq2[-2] - freq2[-3]) # Add a dummy frequency that is equally higher
    freq = freq2
    extension_needed = True
  else:
    freq = freq[::args['freqdownsamplefactor']]
print(freq)

if len(freq) <= 1:
    print('Freq axis needs to at least have length 2')
    sys.exit()

print('FREQS:', freq)

sfreq = freq[0]
dfreq = freq[1] - freq[0]





Nantenna = len(names)
Ntimes = len(time) # 60s timeslots for an 8 hour pointing
Nfreqs = len(freq)
# Set the frequency axis (that TEC doens't use) to 150 MHz as central frequency with 50 MHz of bandwidth.
header='''SIMPLE  =                    T / file does conform to FITS standard
BITPIX  =                  -32 / number of bits per data pixel
NAXIS   =                    6 / number of data axes
NAXIS1  =                 {nsizex:d} / length of RA axis
NAXIS2  =                 {nsizey:d} / length of DEC axis
NAXIS3  =                   4
NAXIS4  =                   {nant:d} / length of ANTENNA axis
NAXIS5  =                    {nfreq:d} / length of FREQ axis
NAXIS6  =                   {ntimes:d} / length of TIME axis
EXTEND  =                    T / FITS dataset may contain extensions
CTYPE1  = 'RA---SIN'           / Right ascension angle cosine
CRPIX1  =                 {crpix1:f} 
CRVAL1  =          {ra:f}
CDELT1  =              {cdelt1:f}
CUNIT1  = 'deg     '
CTYPE2  = 'DEC--SIN'           / Declination angle cosine
CRPIX2  =                  {crpix2:f} 
CRVAL2  =     {dec:f}
CDELT2  =              {cdelt2:f}
CUNIT2  = 'deg     '
CTYPE3  = 'MATRIX  '
CRPIX3  =                   1.
CDELT3  =                   1.
CTYPE4  = 'ANTENNA '
CRPIX4  =                   1.
CRVAL4  =                   0.
CTYPE5  = 'FREQ    '           / Central frequency
CRPIX5  =                   1.
CRVAL5  =     {cfreq:f}
CDELT5  =         {dfreq:f}
CUNIT5  = 'Hz      '
CTYPE6  = 'TIME    '
CRPIX6  =                   1.
CRVAL6  =                   {stime:f} / Should be an AIPS time
CDELT6  =                  {dtime:f}'''.format(nsizex=size, nsizey=size, nant=Nantenna, nfreq=len(freq),ntimes=Ntimes+1, crpix1=np.float(size)/2., ra=np.rad2deg(phasecenter[0]),cdelt1=cdelt1, crpix2=np.float(size)/2., dec=np.rad2deg(phasecenter[1]),cdelt2=cdelt2, cfreq=sfreq, dfreq=dfreq, stime=stime, dtime=dtime)
# note the  ntimes=Ntimes+1 because we append one extra time slice at the end to prevent rounding issues when solving on averaged datasets


with open('header.txt','w') as fobj:
    fobj.write(header)

# Read in h5parm.
h5 = h5parm.h5parm(h5p)
ss = h5.getSolset('sol000')
if includeamps:
  sta = ss.getSoltab('amplitude000')
stp = ss.getSoltab('phase000')
h5_stations = list(st.getAxisValues('ant'))

print(h5_stations, names)

# Find nearest pixel for a given direction.
H = fits.Header.fromstring(header, sep='\n')
wcs = WCS(H).celestial
directions = list(st.getAxisValues('dir'))
sources = ss.getSou()
RA = []
DEC = []
dirs = []
if includeamps:
  amps = reorderAxes(sta.getValues()[0], sta.getAxesNames(), ['time', 'freq', 'ant', 'dir', 'pol'])

phases = reorderAxes(stp.getValues()[0], stp.getAxesNames(), ['time', 'freq', 'ant', 'dir', 'pol'])

if False: # testing
   phases = phases[:,49:57,:,:,:]

if args['freqdownsamplefactor'] != None:
  phases = phases[:,::args['freqdownsamplefactor'],:,:,:]
  if includeamps:
    amps = amps[:,::args['freqdownsamplefactor'],:,:,:]

print(phases.shape)
# re-reference phases to first station (for international all core stations get phase=0 then)
if 'ST001' in h5_stations:
  phases = phaseref(phases, refstation=h5_stations.index('ST001'))
else:
  phases = phaseref(phases, refstation=0)  
print(phases.shape)


# Extend phases by one in the FREQ direction, if needed.
if extension_needed:
    phaseshape = list(phases.shape)
    phaseshape[1] +=1
    newshape = tuple(phaseshape)

    newphases = np.zeros(shape=newshape)
    newphases[:,:-1,:,:,:] = phases
    newphases[:,-1,:,:,:] = phases[:,-1,:,:,:]
    phases = newphases

# Check if XX and YY are the same for factor 2 speedup
scalarpol = False
if np.array_equal(phases[...,-1], phases[...,0]):
  scalarpol = True    
  if includeamps:
    if not np.array_equal(amps[...,-1], amps[...,0]): 
      scalarpol = False
if scalarpol:
  print('XX and YY are the same')      
else:
  print('XX and YY differ')   

      
if includeamps:
  gains = amps * np.exp(1j * phases)
else:
  gains = np.exp(1j * phases)  

del phases # free up RAM
if includeamps:
  del amps # free up RAM
print('min/max angles and amplitudes',np.max(np.angle(gains)) ,  np.min(np.angle(gains)), np.max(np.abs(gains)),np.min(np.abs(gains)))

# Make the RA and DEC vectors
for d in directions:
        c = sources[d]
        diridx = directions.index(d)
        dirs.append(d)
        RAd, DECd = np.rad2deg(c)
        if RAd < 0.0:
            RAd = RAd + 360.
        print('RA:', RAd, 'DEC:', DECd)
        RA.append(RAd)
        DEC.append(DECd)


#hdu.writeto('gainscreen_raw.fits', overwrite=True)
#del data

RA = np.asarray(RA)
DEC = np.asarray(DEC)
# Interpolate the grid using a nearest neighbour approach.
# https://stackoverflow.com/questions/5551286/filling-gaps-in-a-numpy-array



#def interpolate_station(antidx, interpidx, x_from, y_from, tecs, x_to, y_to):
def interpolate_station(antname,ifstep, ntimes, padding=False, includeamps=True, scalarpol=False, smoothamps=False):
    print ('Doing', antname)
    #print('Processing antenna {:s}.'.format(antname))
    if 'ST001' in h5_stations:
      refidx = h5_stations.index('ST001')
    else:
      refidx = 0  
    if ('CS' in antname) and ('ST001' in h5_stations):
        # Take ST001 solutions.
        interpidx = h5_stations.index('ST001')
    else:
        interpidx = h5_stations.index(antname)
    # These will be the new coordinates to evaluate the screen over.
    ra_max, ra_min = wcs.wcs_pix2world(0, 0, 0)[0], wcs.wcs_pix2world(size-1, size-1, 0)[0]
    dec_min, dec_max = wcs.wcs_pix2world(0, 0, 0)[1], wcs.wcs_pix2world(size-1, size-1, 0)[1]

    x = np.arange(size)
    y = np.arange(size)
    xx, yy = np.meshgrid(x, y)

  
    # Do radial basis function interpolation for each time step.
    screen = np.ones((ntimes, 1, 1, 4, size, size), dtype=np.float32)
    X, Y = np.around(wcs.wcs_world2pix(RA, DEC, 0)).astype(int)

    #X, Y = wcs.wcs_world2pix(RA, DEC, 0)
    
    for rara in RA:
        if rara > ra_max or rara< ra_min:
            #print('Info: There are directions in Right Ascension that are outside the screen area', rara)
            if padding:
                print('Info: There are directions in Right Ascension that are outside the screen area', rara)
                print('Cannot use padding in this case, turn if off')
                sys.exit()
    for decdec in DEC:
        if decdec > dec_max or decdec < dec_min:
            #print('Info: There are directions in Declination are outside the screen area', decdec)
            if padding:
                print('Info: There are directions in Declination are outside the screen area', decdec)
                print('Cannot use padding in this case, turn if off')
                sys.exit()
    #Xkeep, Ykeep = np.around(wcs.wcs_world2pix(RA, DEC, 0)).astype(int)
    #print(screen.shape)
    # data has shape (time, freq, ant, matrix, y, x)
    # gains has shape (time, freq, ant, dir, pol)
    # Iterate over all timeslots.
    for itstep in range(ntimes):
        # Interpolate the gains, not the Re/Im or Amp/Phase separately.
        gXX = gains[itstep, ifstep, interpidx, :, 0]
        gYY = gains[itstep, ifstep, interpidx, :, -1]
        
        
        #skip_because_samegain = False
        #try:
        #  if np.array_equal(gains[itstep, ifstep, interpidx, :, :], gains[itstep, ifstep, interpidx-1, :, :]):
        #     skip_because_samegain = True
        #except:
        #  pass  
        #if skip_because_samegain:
        #   return names.index(antname), screen, [0, size-1], [0, size-1], gXXcom, gYYcom, Xcom, Ycom    
             
             
        if padding: 
            # for bottom X-axis and top X-axis
            for pp in np.arange(0,size,args['stepsizepadding']):
              # bottom
              if itstep == 0:  #otherwise the X, Y arrays keep growing!!
                X = np.append(X,pp)
                Y = np.append(Y,0)
              gXX = np.append(gXX,1.0+1j*0.)
              gYY = np.append(gYY,1.0+1j*0.)
              # top
              if itstep == 0:  #otherwise the X, Y arrays keep growing!!
                X = np.append(X,pp)
                Y = np.append(Y,np.int(size)-1)
              gXX = np.append(gXX, 1.0+1j*0.)
              gYY = np.append(gYY, 1.0+1j*0.)
            # for left Y-axis and right Y-axis
            for pp in np.arange(0,size,args['stepsizepadding']):
              # left
              if itstep == 0:  #otherwise the X, Y arrays keep growing!!
                Y = np.append(Y,pp)
                X = np.append(X,0)
              gXX = np.append(gXX,1.0+1j*0.)
              gYY = np.append(gYY,1.0+1j*0.)
              # right
              if itstep == 0:  #otherwise the X, Y arrays keep growing!!
                Y = np.append(Y,pp)
                X = np.append(X,np.int(size)-1)
              gXX = np.append(gXX,1.0+1j*0.)
              gYY = np.append(gYY,1.0+1j*0.)
        
        
        rbfXX = Rbf(X, Y, gXX, smooth=1e-8)
        if scalarpol:
          rbfYY = rbfXX
        else:
          rbfYY = Rbf(X, Y, gYY, smooth=1e-8)  

        if smoothamps:
           #rbfXXphase = Rbf(X, Y, np.exp(1j*np.angle(gXX)), smooth=1e-8)
           rbfXXphase = Rbf(X, Y, gXX, smooth=1e-8) # do include amps here, makes the phases smoother at the end
           if scalarpol:
             rbfYYphase = rbfXXphase
           else:
             #rbfYYphase = Rbf(X, Y, np.exp(1j*np.angle(gYY)), smooth=1e-8)  
             rbfYYphase = Rbf(X, Y, gYY, smooth=1e-8) # do include amps here, makes the phases smoother at the end


           rbfXXamp = Rbf(X, Y, np.abs(gXX)*np.exp(1j*0.0), smooth=1e-2)
           if scalarpol:
             rbfYYamp = rbfXXamp
           else:
             rbfYYamp = Rbf(X, Y, np.abs(gYY)*np.exp(1j*0.0), smooth=1e-2)  

           tinterpXX = np.abs(rbfXXamp(xx, yy))*np.exp(1j*np.angle(rbfXXphase(xx, yy)))   
           if scalarpol:
             tinterpYY = tinterpXX
           else:
             tinterpYY = np.abs(rbfYYamp(xx, yy))*np.exp(1j*np.angle(rbfYYphase(xx, yy)))   
        else: 
           tinterpXX = rbfXX(xx, yy)
           if scalarpol:
             tinterpYY = tinterpXX
           else:
             tinterpYY = rbfYY(xx, yy)
        
        
        if not includeamps:
          # reset amps to 1.0  
          tinterpXX = np.exp(1j *np.angle(tinterpXX))
          if scalarpol:
            tinterpYY = tinterpXX  
          else:
            tinterpYY = np.exp(1j *np.angle(tinterpYY))
        
        #print ph.shape
        #print np.mean(ph), 'beg mean'
        
        #ph = ndimage.gaussian_filter(ph, sigma=(7, 7))
        #print np.mean(ph), 'end mean'
        #print ph.shape
        #sys.exit()
        #ap = np.abs(tinterpXX)
        #tinterpXX = np.exp(1j * ph)
        #tinterpYY = tinterpXX

        #data[t, f, s, 0, ind[0], ind[1]] = 1. * np.cos(pb)
        #data[t, f, s, 2, ind[0], ind[1]] = val_amp_yy * np.cos(pb)
        #data[t, f, s, 1, ind[0], ind[1]] = val_amp_xx * np.sin(pb)
        #data[t, f, s, 3, ind[0], ind[1]] = val_amp_yy * np.sin(pb)
        #print(tinterpXX.shape)
        #print(tinterpYY.shape)
        #print(Y.shape, X.shape, gXX.shape, gYY.shape)

        # remove directions outside screen area for the allclose comparison below
        Xcom=[]
        Ycom=[]
        gXXcom = []
        gYYcom = []
        #print(X)
        #sys.exit()
        for ddd in range(len(X)):
          if not (X[ddd] > (np.float(size)-1.) or X[ddd]< 0.0 or Y[ddd] > (np.float(size)-1.) or Y[ddd] < 0.0): # outside screen
            #print X[ddd], Y[ddd]
            Xcom = np.append(Xcom,X[ddd])
            Ycom = np.append(Ycom,Y[ddd])
            gXXcom = np.append(gXXcom,gXX[ddd])  
            gYYcom = np.append(gYYcom,gYY[ddd]) 
            
      
        if not np.allclose(tinterpXX[Ycom.astype(int),Xcom.astype(int)], gXXcom, rtol=1e-1, atol=1e-1):
            #print(itstep, ifstep, interpidx)
            #print('gXXcom', gXXcom)
            #print('tinterpXX', tinterpXX[Ycom.astype(int),Xcom.astype(int)])
            #print('diff', gXXcom - tinterpXX[Ycom.astype(int),Xcom.astype(int)])
            print('Interpolated screen for polarization XX does not go through nodal points.')
            #raise ValueError('Interpolated screen for polarization XX does not go through nodal points.')
        
        if not scalarpol:
          if not np.allclose(tinterpYY[Ycom.astype(int),Xcom.astype(int)], gYYcom, rtol=1e-1, atol=1e-1):
             raise ValueError('Interpolated screen for polarization YY does not go through nodal points.')
        
        del gXX, gYY #, Ycom, Xcom #, gXXcom, gYYcom # free up RAM
        matrix = np.asarray([np.real(tinterpXX), np.imag(tinterpXX), np.real(tinterpYY), np.imag(tinterpYY)])
        del tinterpXX, tinterpYY
        
        screen[itstep, 0, 0, :, :, :] = matrix

    #print( [ra_min, ra_max], [dec_min, dec_max])
    return names.index(antname), screen, [0, size-1], [0, size-1], gXXcom, gYYcom, Xcom, Ycom
    #return names.index(antname), screen, [ra_min, ra_max], [dec_min, dec_max], gXXcom, gYYcom


def single_writerun(ifreq):
    print('Processing frequency slot {:d}'.format(ifreq))
    data_int = np.ones((Ntimes+1,Nantenna,4,size,size),dtype=np.float32)
    tupples = [(station, ifreq, Ntimes, args['padding'], includeamps, scalarpol) for station in names]
    for station in names:
        antenna, screen, xxwcs, yywcs, gXXcom, gYYcom, RA_X, DEC_Y = interpolate_station(station, ifreq, Ntimes, padding=args['padding'], includeamps=includeamps, scalarpol=scalarpol)
        #print('Screen: ', screen.shape)
        #print('data_int: ', data_int.shape)
        data_int[:-1, antenna, :, :, :] = screen[:, 0, 0, :, :, :]
        data_int[-1, antenna, :, :, :] = screen[-1, 0, 0, :, :, :]
    hdu = fits.PrimaryHDU()
    hdu.data = data_int
    hdu.writeto('partfits_{}.fits'.format(ifreq),overwrite=True)
    del data_int
    del hdu.data
    del hdu
    gc.collect()

# Inspired by https://stackoverflow.com/questions/38309535/populate-numpy-array-through-concurrent-futures-multiprocessing
#print('Shape of the output SCREEN:',data.shape)


#
#data_int = np.ones(data.shape, dtype=np.float32)
# data_int = np.ones((Ntimes+1, Nfreqs, Nantenna, 4, size, size), dtype=np.float32)

ncpu = args['ncpu']

if ncpu >1:
    # Pool the single runs. Afterwards, concatenate.
    pl = mp.Pool(ncpu)
    pl.map(single_writerun,np.arange(Nfreqs))

    dims = H['NAXIS']
    shape = []

    for i in range(dims):
        shape.append(H[f'NAXIS{i+1}'])

    shape = tuple(shape[::-1]) # Byte order is inverted for fits files 

    # Prepare all the headers

    TIMEBLOCKS = args['timeblocks']

    timesteps_per_block = shape[0]//TIMEBLOCKS
    startpix = 1

    # Compute the amount of timesteps precisely

    datasizes = []
    start = 0
    starts = []
    for i in range(TIMEBLOCKS-1):
        datasizes.append(timesteps_per_block)
        starts.append(start)
        start+=timesteps_per_block
    datasizes.append(shape[0] - start) 
    starts.append(starts[-1]+timesteps_per_block)

    # Construct the headers
    headers = []
    for j in range(TIMEBLOCKS):
        Hloc = copy.deepcopy(H)
        Hloc['CRVAL6'] += Hloc['CDELT6'] * starts[j]
        headers.append(Hloc)


    if TIMEBLOCKS == 1:
        # No time averaging
        data_int = np.zeros(shape,dtype=np.float32)
        # Resort fitslist
        fitslist = glob.glob('partfits*')
        fitsnums = [int(fit.split('_')[1].split('.')[0]) for fit in fitslist]
        fitsnums = np.sort(np.array(fitsnums))
        fitslist = [f'partfits_{i}.fits' for i in fitsnums]
        for j in fitsnums:
            data_int[:,j,:,:,:,:] = fits.getdata(fitslist[j])

        hdu = fits.PrimaryHDU(header=H)
        hdu.data = data_int
        hdu.writeto(fitsfilename)

    else:
        for timeslot in range(TIMEBLOCKS):
            # Do this timeslot, one at the time.
            print('Processing timeslot '+str(timeslot))
            start = starts[timeslot]
            numtimesteps = datasizes[timeslot]
            head = headers[timeslot]

            timeslot_shape = list(shape)
            timeslot_shape[0] = numtimesteps
            data_int = np.zeros(tuple(timeslot_shape),dtype=np.float32)

            fitslist = glob.glob('partfits*')
            fitsnums = [int(fit.split('_')[1].split('.')[0]) for fit in fitslist]
            fitsnums = np.sort(np.array(fitsnums))
            fitslist = [f'partfits_{i}.fits' for i in fitsnums]
            for j in fitsnums:
                data_int[:,j,:,:,:,:] = fits.getdata(fitslist[j])[start:start+numtimesteps,:,:,:,:]

            hdu = fits.PrimaryHDU(header=head)
            hdu.data = data_int
            del data_int
            hdu.writeto(f'{fitsfilename_prefix}_timeslot{timeslot}.fits')
            del hdu.data
            del hdu
            gc.collect()

    fitslist = glob.glob('partfits*')
    fitsnums = [int(fit.split('_')[1].split('.')[0]) for fit in fitslist]
    fitsnums = np.sort(np.array(fitsnums))
    fitslist = [f'partfits_{i}.fits' for i in fitsnums]
    print('Removing temporary files...')

    for fit in fitslist:
        os.system(f'rm -rf {fit}')

if ncpu == 1:
    data_int = np.ones((Ntimes+1, Nfreqs, Nantenna, 4, size, size), dtype=np.float32)
    for ifreq in range(Nfreqs):
        print('Processing frequency slot {:d}'.format(ifreq))
        tupples = [(station, ifreq, Ntimes, args['padding'], includeamps, scalarpol) for station in names]
        for station in names:
            antenna, screen, xxwcs, yywcs, gXXcom, gYYcom, RA_X, DEC_Y = interpolate_station(station, ifreq, Ntimes, padding=args['padding'], includeamps=includeamps, scalarpol=scalarpol)
            #print('Screen: ', screen.shape)
            #print('data_int: ', data_int.shape)
            data_int[:-1, ifreq, antenna, :, :, :] = screen[:, 0, 0, :, :, :]
            data_int[-1, ifreq, antenna, :, :, :] = screen[-1, 0, 0, :, :, :]

    if TIMEBLOCKS == 1:
        hdu = fits.PrimaryHDU(header=H)
        hdu.data = data_int
        hdu.writeto(args['FITSscreen'], overwrite=True)
        print('Finished interpolating.')
    else:
        dims = H['NAXIS']
        shape = []

        for i in range(dims):
            shape.append(H[f'NAXIS{i+1}'])

        shape = tuple(shape[::-1]) # Byte order is inverted for fits files 

        # Prepare all the headers

        TIMEBLOCKS = args['timeblocks']

        timesteps_per_block = shape[0]//TIMEBLOCKS
        startpix = 1

        # Compute the amount of timesteps precisely

        datasizes = []
        start = 0
        starts = []
        for i in range(TIMEBLOCKS-1):
            datasizes.append(timesteps_per_block)
            starts.append(start)
            start+=timesteps_per_block
        datasizes.append(shape[0] - start) 
        starts.append(starts[-1]+timesteps_per_block)

        # Construct the headers
        headers = []
        for j in range(TIMEBLOCKS):
            Hloc = copy.deepcopy(H)
            Hloc['CRVAL6'] += Hloc['CDELT6'] * starts[j]
            headers.append(Hloc)

        for timeslot in range(TIMEBLOCKS):
            # Do this timeslot, one at the time.
            print('Processing timeslot '+str(timeslot))
            start = starts[timeslot]
            numtimesteps = datasizes[timeslot]
            head = headers[timeslot]

            timeslot_shape = list(shape)
            timeslot_shape[0] = numtimesteps
            data_int = np.zeros(tuple(timeslot_shape),dtype=np.float32)

            fitslist = glob.glob('partfits*')
            fitsnums = [int(fit.split('_')[1].split('.')[0]) for fit in fitslist]
            fitsnums = np.sort(np.array(fitsnums))
            fitslist = [f'partfits_{i}.fits' for i in fitsnums]
            for j in fitsnums:
                data_int[:,j,:,:,:,:] = fits.getdata(fitslist[j])[start:start+numtimesteps,:,:,:,:]

            hdu = fits.PrimaryHDU(header=head)
            hdu.data = data_int
            del data_int
            hdu.writeto(f'{fitsfilename_prefix}_timeslot{timeslot}.fits')
            del hdu.data
            del hdu
            gc.collect()


