#!/usr/bin/env python
# coding: utf-8

# In[21]:


import numpy as np
from astropy.io import fits
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from regions import write_ds9,RectangleSkyRegion
from astropy.wcs import WCS
import sys,os,glob,argparse,bdsf

parse = argparse.ArgumentParser(description='Generate regions for DD calibration')
parse.add_argument('--calibrator', help='Pointing of the calibrator')
parse.add_argument('--thres',help='Baseline flux threshold for detection', default=1,type=float)
parse.add_argument('--boxsize', help='Boxsize of the regions', default=0.3, type=float)
parse.add_argument('--sep', help='Minimal seperation between the centers of the boxes',default=1,type=float)
parse.add_argument('--largedistancebias', help='Biases larger distances, basically a first order approximation of the beam', default=0.5, type=float)

res = parse.parse_args()

# In[28]:

# Get skycoord of calibrator
cal = res.calibrator
splitted = cal.split('[')[1].split(']')[0]
ra = splitted.split('deg')[0]
dec = splitted.split('deg')[1][1:]
crd_cal = SkyCoord(ra=ra,dec=dec, unit=(u.deg,u.deg))

filename = glob.glob('*MFS-image.fits')[0]

img = bdsf.process_image(filename, rms_box = (640,160), rms_map=True, thresh='hard', thresh_isl=3.0, thresh_pix=10.0)
img.write_catalog(outfile='regions_wsclean1.fits', bbs_patches='single', catalog_type='srl', clobber=True, format='fits')

regions = fits.getdata('regions_wsclean1.fits')
reghead = fits.open('regions_wsclean1.fits')[1].header
mainfile = fits.getheader(filename)

mainpoint = SkyCoord(ra=mainfile['CRVAL1'],dec=mainfile['CRVAL2'],unit=(u.deg,u.deg))


# In[25]:





# In[3]:


RA = regions['RA']
DEC = regions['DEC']
tflx = regions['Total_flux']

sort_indices = np.argsort(tflx)[::-1]
RA_sorted = RA[sort_indices]
DEC_sorted = DEC[sort_indices]
tflx_sorted = tflx[sort_indices]


# In[4]:


def testInBox(ra_1,dec_1,ra_2,dec_2,sep=1):
    sep_ra = np.abs(ra_2-ra_1)
    sep_dec = np.abs(dec_2-dec_1)
    if sep_ra > sep or sep_dec > sep:
        return False
    else:
        return True
    
def testInRadius(ra_1,dec_1,ra_2,dec_2,sep=1):
    coord1 = SkyCoord(ra_1,dec_1,frame='fk5',unit='deg')
    coord2 = SkyCoord(ra_2,dec_2,frame='fk5',unit='deg')
    separation = coord1.separation(coord2).degree
    if separation>sep:
        return False
    else:
        return True


# In[33]:


directions_list = []
thres = res.thres # Jy
boxsize = res.boxsize # degree
sep = res.sep # degree
largedistancebias = res.largedistancebias # factor per degree which makes faint edge sources more prominent
separate_files = False

for i in range(len(RA)):
    if len(directions_list) == 0:
        directions_list.append((RA_sorted[i],DEC_sorted[i],tflx_sorted[i]))
    else:
        inbox = False
        for j in range(len(directions_list)):
            boxRA,boxDEC,boxflx = directions_list[j]
            # inbox = testInBox(RA_sorted[i],DEC_sorted[i],boxRA,boxDEC,sep)
            inbox = testInRadius(RA_sorted[i],DEC_sorted[i],boxRA,boxDEC,sep)
            if inbox:
                break
        if inbox:
            # Match found
            new_tuple = (boxRA,boxDEC,boxflx+tflx_sorted[i])
            directions_list[j] = new_tuple
        else:
            # No match found
            directions_list.append((RA_sorted[i],DEC_sorted[i],tflx_sorted[i]))

significant_directions = []
for ii in range(len(directions_list)):
    direc = directions_list[ii]
    pointing = SkyCoord(ra=direc[0],dec=direc[1],unit=(u.deg,u.deg))
    dist = pointing.separation(mainpoint).value
    caldist = pointing.separation(crd_cal).value
    if caldist > 0.2:
        if direc[2] > (thres-largedistancebias*dist):
            significant_directions.append(direc)
    if len(significant_directions) > 50:
        break


# In[34]:


# Here be some pyregion stuff
centers = [SkyCoord(ii[0],ii[1], frame='fk5',unit='deg') for ii in significant_directions]
shapelist = [RectangleSkyRegion(center,2*boxsize*u.deg,2*boxsize*u.deg) for center in centers]
foldname = 'regions_ws1'

if separate_files:
    try:
        os.mkdir(foldname)
    except FileExistsError:
        os.system(f'rm -rf {foldname}')
        os.mkdir(foldname)
    for i,shape in enumerate(shapelist):
        write_ds9([shapelist[i]],f'{foldname}/Dir{i}.reg',overwrite=True)
else:
    parser = write_ds9(shapelist,f'{foldname}.reg',overwrite=True)
