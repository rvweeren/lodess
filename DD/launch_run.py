import os
import sys
import glob

# Test if this is the first run
try:
    os.mkdir('PROGRESS')
    directions = glob.glob('rectangles/*')
    directionnames = [dirry.lstrip('rectangles/').rstrip('.reg') for dirry in directions]
    for directionname in directionnames:
        os.system(f'touch PROGRESS/{directionname}')
except:
    pass

# Initialize this run

make_success = False
try:
    runname = sys.argv[1]
except:
    runname = None
while make_success == False:
    if runname == None:
        runname = input('Name this run: ')
    else:
        make_success = True
        print(runname)
    try:
        os.mkdir(f'run_{runname}')
        make_success = True
    except:
        print('Very funny. Now name it something unique')
        runname = None

modelfiles = glob.glob('*fits')
for model in modelfiles:
    os.system(f'cp -r {model} run_{runname}/{model}')
os.system(f'cp -r *.py run_{runname}/')
msnames = glob.glob('*.ms')
os.system(f'cp -r {msnames} run_{runname}')

# iterate through each direction
while len(os.listdir('PROGRESS/')) != 0:
    dirs_avail_nums = [int(num.lstrip('PROGRESS/Dir').rstrip('.reg')) for num in glob.glob('PROGRESS/*')]
    dirs_avail_nums.sort()
    chosen_direction = 'Dir'+str(dirs_avail_nums[0])+'.reg'
    chosen_dir = 'Dir'+str(dirs_avail_nums[0])
    i = str(dirs_avail_nums[0])
    print(i)

    os.system(f'rm -rf PROGRESS/{chosen_dir}')
    os.system(f'cp -r rectangles/{chosen_direction} run_{runname}/{chosen_direction}')
    os.chdir(f'run_{runname}')

    for j,msname in enumerate(msnames):
        os.system(f'python3 standalone_peel.py {msname} {chosen_dir} {j}')

    # Now, copy files to dedicated location

    os.mkdir(f'direction{i}')
    os.system(f'cp -r {chosen_direction} direction{i}')
    os.system(f'cp -r calibrate.py direction{i}')
    os.system(f'cp -r {chosen_dir}.*.peel.ms direction{i}')
    os.chdir(f'direction{i}')
    os.system(f'python3 calibrate.py {i}')
    os.chdir('../../')

    os.system(f'rm -rf run_{runname}/{chosen_direction}')
