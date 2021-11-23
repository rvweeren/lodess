import numpy as np
import os,sys
from astropy.coordinates import SkyCoord
import astropy.units as u

targets = [ {'name' : 'CasA', 'ra' : 6.123487680622104,  'dec' : 1.0265153995604648},
            {'name' : 'CygA', 'ra' : 5.233686575770755,  'dec' : 0.7109409582180791},
            {'name' : 'TauA', 'ra' : 1.4596748493730913, 'dec' : 0.38422502335921294},
            {'name' : 'HerA', 'ra' : 4.4119087330382163, 'dec' : 0.087135562905816893},
            {'name' : 'VirA', 'ra' : 3.276086511413598,  'dec' : 0.21626589533567378},
            {'name' : 'HydraA', 'ra' : 2.4352, 'dec' : -0.21099}]


point = SkyCoord(sys.argv[1], unit=(u.deg,u.deg))
for target in targets:
    name = target['name']
    ra = target['ra']
    dec = target['dec']
    target_coord = SkyCoord(ra=ra,dec=dec,unit=(u.rad,u.rad))
    print(f'Separation to {name} : {target_coord.separation(point)}')


pointgal = point.transform_to('galactic')
print(f'Galactic Coordinates: l: {pointgal.l.value:.4f} deg     b: {pointgal.b.value:.4f} deg')
