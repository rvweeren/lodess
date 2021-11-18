#!/usr/bin/env python
# coding: utf-8



import matplotlib.pyplot as plt
import numpy as np
from astroquery.vizier import Vizier
from astroquery.ned import Ned
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u
import argparse as argp

parser = argp.ArgumentParser()
parser.add_argument('crdstring', help='Coordinates to look at.')
parser.add_argument('--maxdistance',default=3,type=int,help='Maximum distance for a given source')

res = parser.parse_args()

Vizier.ROW_LIMIT = -1


# In[2]:


# CRD = "337.6670deg 82.2650deg"
# maxrad = 3

CRD = res.crdstring
maxrad = res.maxdistance

coord = SkyCoord(CRD, unit=(u.deg,u.deg))


# In[3]:


sources = Vizier.query_region(coord,radius=maxrad*u.deg,catalog='VIII/97')[0]
maxfluxpos = np.argmax(sources['Sp'])
target = sources[maxfluxpos]
targetpos = SkyCoord(target['RAJ2000']+target['DEJ2000'],unit=(u.hourangle,u.deg))
print(targetpos)

simb = Simbad.query_region(targetpos,radius=30*u.arcsec)
objname = simb['MAIN_ID'][0]


# In[16]:


res = Ned.get_table(objname,table='photometry')
freq = res['Frequency']
flxd = res['Flux Density']
lowfreqs = freq < 5e9

freq = freq[lowfreqs]
flxd = flxd[lowfreqs]


# In[20]:


print(target)
print(f'\n\nDistance to    {objname}:  ', round(targetpos.separation(coord).to('deg').value,4),' degrees')


plt.scatter(freq,flxd,marker='x',color='black')
plt.loglog()
plt.title(objname)
plt.ylabel('Flux density (Jy)')
plt.xlabel('Frequency (Hz)')
plt.show()
