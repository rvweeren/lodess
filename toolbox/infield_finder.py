import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
import argparse
from astroquery.vizier import Vizier
from astroquery.ned import Ned
from astroquery.simbad import Simbad
from matplotlib.patches import Circle

Vizier.ROW_LIMIT = -1
targets = [ {'name' : 'CasA', 'ra' : 6.123487680622104,  'dec' : 1.0265153995604648},
            {'name' : 'CygA', 'ra' : 5.233686575770755,  'dec' : 0.7109409582180791},
            {'name' : 'TauA', 'ra' : 1.4596748493730913, 'dec' : 0.38422502335921294},
            {'name' : 'HerA', 'ra' : 4.4119087330382163, 'dec' : 0.087135562905816893},
            {'name' : 'VirA', 'ra' : 3.276086511413598,  'dec' : 0.21626589533567378},
            {'name' : 'HydraA', 'ra' : 2.4352, 'dec' : -0.21099}]

def detect_ateam(point):
    for target in targets:
        name = target['name']
        ra = target['ra']
        dec = target['dec']
        target_coord = SkyCoord(ra=ra,dec=dec,unit=(u.rad,u.rad))
        print(f'Separation to {name} : {target_coord.separation(point)}')

    pointgal = point.transform_to('galactic')
    print(f'Galactic Coordinates: l: {pointgal.l.value:.4f} deg     b: {pointgal.b.value:.4f} deg')


def fix_name(name):
    # I really hate that this is necessary. But well
    clist = ['8C','4C','3C']
    for c in clist:
        if c in name:
            cataloglist = name.split(c)[1].split(' ')
            name = c + next(s for s in cataloglist if s)
    if '2C' in name:
        res = Simbad.query_objectids(name)['ID']
        for i in res:
            if '3C ' in i:
                name = i
    return name


def on_pick(event):
    artist = event.artist
    data = artist.get_offsets()
    # find brightest object near click
    selfluxes = np.array(fluxes)[event.ind]
    selected_index = event.ind[np.argmax(selfluxes)]

    coord = w.pixel_to_world(*data[selected_index])
    print(f'"({coord.ra.value}, {coord.dec.value})"')
    process_pointing(coord)

def process_pointing(coord):
    spectrum.clear()
    simb = Simbad.query_region(coord,radius=2*u.arcmin)
    namelist = simb['MAIN_ID']
    # Order of precedence in the naming convention
    name = [namelist[0]]
    name += list(filter(lambda x: "NVSS" in x,namelist))
    name += list(filter(lambda x: "8C" in x,namelist))
    name += list(filter(lambda x: "4C" in x,namelist))
    name += list(filter(lambda x: "3C" in x,namelist))
    name = name[-1]
    name = fix_name(name)
    spectrum.set_title(name)

    res = Ned.get_table(name, table='photometry')
    freq = res['Frequency']
    flxd = res['Flux Density']
    lowfreqs = freq < 5e9
    freq = freq[lowfreqs]
    flxd = flxd[lowfreqs]
    spectrum.scatter(freq,flxd,color='black')
    spectrum.loglog()
    spectrum.set_ylabel('Flux (Jy)')
    spectrum.set_xlabel('Frequency (Hz)')
    spectrum.axvspan(15e6,30e6,color='lightgreen',alpha=0.4)

    fig.canvas.draw()

def find_sources(pointing,radius):
    sources = Vizier.query_region(pointing,radius=radius,catalog='VIII/97',column_filters={'Sp':'>3'})[0]
    ra = sources['RAJ2000']
    dec = sources['DEJ2000']
    flux = list(sources['Sp'])
    coords = [SkyCoord(rasel,desel,unit=(u.hourangle,u.deg)) for rasel,desel in zip(ra,dec)]
    return coords,flux

if __name__ == '__main__':
    args = argparse.ArgumentParser()
    args.add_argument('coord', help='Coordinates of the pointing')
    args.add_argument('--fov', help='Size of the final image in degrees',default = 18,type=float)
    args.add_argument('--circlerad', help='Size of the "guiding circle"',default=5, type=float)

    parser = args.parse_args()
    crd = SkyCoord(parser.coord,unit=(u.deg,u.deg))
    detect_ateam(crd)

    w = WCS(naxis=2)
    w.wcs.crpix = [1024,1024] # for resolution of 2048pix
    w.wcs.crval = parser.coord.split('+')
    w.wcs.cdelt = [parser.fov/2048,parser.fov/2048]
    w.wcs.ctype = ["RA---SIN","DEC--SIN"]

    fig = plt.figure(figsize=(18,7))
    axes = fig.add_subplot(121,projection=w)
    axes.set_transform(w)
    axes.coords[0].set_ticks(number=10)
    axes.coords[1].set_ticks(number=10)
    axes.grid(True)
    axes.set_ylim(0,2048)
    axes.set_xlim(0,2048)

    coords,fluxes = find_sources(crd,0.5*parser.fov*u.deg)
    pixels_positions = np.array([w.world_to_pixel(coordsel) for coordsel in coords]).T
    plt.scatter(pixels_positions[0],pixels_positions[1],c=fluxes,picker=10,label = 'VLSSR points')
    cb = plt.colorbar()
    cb.set_label("Flux from VLSSR (Jy, @74MHz)")

    axes.scatter([1024],[1024],marker='X', color='orangered',s=60,picker=10,label='Pointing center')
    axes.add_artist(Circle((1024,1024),parser.circlerad/parser.fov*2048,fill=False,edgecolor='orangered',linestyle='--',label=f'Circle (radius: {parser.circlerad} deg)'))

    fig.canvas.mpl_connect('pick_event', on_pick)
    
    spectrum = fig.add_subplot(122)

    plt.show()
