import sys
import shutil
import os
import glob
from h5_merger import merge_h5

'''
    This code automatically extracts all directions from the run files.
'''

def elstrip(s,prefix):
    if s.startswith(prefix):
        s = s[len(prefix):]
    return s

def process_direction(path, dirnum):
    print(f'direction{dirnum}...')
    os.system(f'cp -r {path}/image_005-MFS-image.fits RESULTS/fits/direction{dirnum}.fits')
    plotlosotodir = glob.glob(f'{path}/plotlosoto*')[0]
    tecandphaseplots = sorted(glob.glob(f'{plotlosotodir}/tecandphase0*5*'))[0]
    os.system(f'cp -r {tecandphaseplots} RESULTS/tecandphase/direction{dirnum}.png')
    h5file = sorted(glob.glob(f'{path}/merged*h5'))[-2]
    os.system(f'cp -r {h5file} RESULTS/h5files/direction{dirnum}.h5')
    images = sorted(glob.glob(f'{path}/*png'))
    os.mkdir(f'RESULTS/images/direction{dirnum}')
    for img in images:
        os.system(f'cp -r {img} RESULTS/images/direction{dirnum}/')

separate_runs = glob.glob('run*')

try:
    os.mkdir ('RESULTS')
except:
    shutil.rmtree('RESULTS')
    os.mkdir('RESULTS')

os.mkdir('RESULTS/fits')
os.mkdir('RESULTS/images')
os.mkdir('RESULTS/tecandphase/')
os.mkdir('RESULTS/h5files/')
os.system(f'cp -r ../LOW_3-MFS-image.fits RESULTS/fits/DI.fits')

dirnums = []
for run in separate_runs:
    directions = glob.glob(f'{run}/direction*')
    for dirry in directions:
        dirnum = elstrip(dirry,f'{run}/direction')
        process_direction(f'{dirry}',dirnum)
        dirnums.append(dirnum)

available_regions = [elstrip(reg,'rectangles/Dir').rstrip('.reg') for reg in glob.glob('rectangles/*')]

for i in dirnums:
    if i not in available_regions:
        print(f'Direction {i} is not calibrated!!')

# Now merge all the H5 files
h5files = sorted(glob.glob('RESULTS/h5files/direction*'))
print('Merging directions...')

example_dir = glob.glob(f'{separate_runs[0]}/direction*')[-1]
example_ms = glob.glob(f'{example_dir}/Dir*peel.ms')[-1]
print(example_ms)
# merge_h5("RESULTS/merged_all_directions_DD.h5", h5_tables=h5files, ms_files=example_ms,convert_tec=True,add_directions=[5.629128,1.465345])
