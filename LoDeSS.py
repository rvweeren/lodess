#! /usr/bin/python3
import sys
import time
import threading
import argparse
import multiprocessing as mp
import subprocess
import numpy as np
import pyrap.tables as pt
import os
from regions import RectangleSkyRegion
import losoto.h5parm as h5parm
from astropy.coordinates import SkyCoord
import astropy.units as u
import glob
import bdsf


'''
    Input: *msdemix files
           cal_solutions.h5
           skymodel, with the name of the skymodel prepended to it (e.g. 3C380-alsjfoieaj.skymodel)
    Output: stuff
            fitsfiles
'''

freqstep = 1

c3c380= np.array([-1.44194739, 0.85078014])
c3c196= np.array([2.15374139,0.8415521])

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

def _run_sing(runname):
    # sing_extension = 'singularity exec -B /data1,/net/rijn1/data2,/net/rijn2/data2/,/net/rijn3/data2,/net/rijn4/data2,/net/rijn5/data2,/net/rijn6/data2/net/rijn7/data2/,/net/rijn/data2,/net/hoendiep/data2,/net/bovenrijn/data1/,/net/ouderijn/data2,/net/tussenrijn/data2 /net/rijn/data2/groeneveld/LoDeSS_files/pill-latestSep.simg'
    cmd = f'python launch_run.py {runname}'
    print(cmd)
    os.system(cmd)

def find_skymodel():
    measurementset = glob.glob('*msdemix')[2]
    tab = pt.table(measurementset+'::FIELD')
    adir = tab.getcol('DELAY_DIR')[0][0][:]
    tab.close()

    cdatta = SkyCoord(adir[0]*u.deg, adir[0]*u.deg, frame='icrs')

    targetsource = None

    if (cdatta.separation( SkyCoord(c3c380[0]*u.deg, c3c380[0]*u.deg, frame='icrs'))) < 0.2*u.deg:
        targetsource = '3c380'    
    if (cdatta.separation( SkyCoord(c3c196[0]*u.deg, c3c196[0]*u.deg, frame='icrs'))) <  0.2*u.deg:
        targetsource = '3c196' 
    return targetsource

def find_missing_stations():
    measurementset = glob.glob('*msdemix')[0]
    h5 = 'Band_PA.h5'

    # First look at the h5

    H5 = h5parm.h5parm(h5)
    a = H5.getSolset('calibrator')
    antlist_h5 = a.getSoltab('bandpass').getValues()[1]['ant']
    H5.close()

    # Now look at ms

    tab = pt.table(measurementset+'::ANTENNA')
    antlist_ms = tab.getcol('NAME')

    # compare them

    missinglist = []
    for ant in antlist_ms:
        if ant not in antlist_h5:
            missinglist.append(ant)
    return missinglist

def generate_boxfile(direction):
    width = 8*512*u.arcsec
    height = width
    parsed_dir = direction[1:-1].replace(',','+')
    coord = SkyCoord(parsed_dir,unit=(u.deg,u.deg))
    region = RectangleSkyRegion(center=coord,width=width,height = height)
    region.write('boxfile.reg',format='ds9')
    
def find_rms(thres = 0.002):
    '''
        Run this in the DD_cal directory. Iterates through all the run_X folders,
        finds the measurement sets and computes the 
    '''
    import runwsclean as runw
    msses = glob.glob('run*/direction*/Dir*ms')
    msnames = [ms.split('/')[-1].split('.')[0] for ms in msses]
    msnums = [int(ms.split('Dir')[-1]) for ms in msnames]
    sorting = np.argsort(msnums)

    msses = np.array(msses)[sorting]
    msnames = np.array(msnames)[sorting]
    msnums = np.array(msnums)[sorting]

    snrs = []
    for ms,msname in zip(msses,msnames):
        try:
            noise,flux,_,__ = runw.getmsmodelinfo(ms,'MODEL_DATA',fastrms=True)
            snr = flux/noise
            snrs.append(snr)
        except:
            snrs.append(np.nan)
    
    snrs = np.array(snrs)
    toreject = np.where(snrs < thres)[0]
    os.chdir('RESULTS')
    os.mkdir('rejected')
    for torej in toreject:
        num = msnums[torej]
        os.system(f'mv h5files/direction{num}.h5 rejected/direction{num}.h5')
    os.chdir('../')

def run(comb):
    '''
        Processes demixed datasets and outputs a 'raw' dataset, where only bandpass+polaligns
        are applied
    '''
    ms,missinglist = comb[0],comb[1]
    # First filter out missing stations
    if len(missinglist) > 0:
        msout = ms.split('.msdemix')[0] + '.split.ms'
        cmd =  'DPPP numthreads=80 msin=' + ms + ' msout.storagemanager=dysco '
        cmd += 'msout='+msout + ' '
        cmd += 'msout.writefullresflag=False '
        cmd += 'steps=[f2] '

        cmd += 'f2.type=filter f2.remove=True f2.baseline="'
        for missing in missinglist:
            cmd += f'!{missing}&&*;'
        cmd = cmd.rstrip(';')
        cmd += '"'
        print(cmd)
        subprocess.call(cmd, shell=True)
        msin = ms.split('.msdemix')[0] + '.split.ms'
    else:
        msin = ms
    msout = ms.split('.msdemix')[0] + '.corr.ms'
    msout = msout.split('archive/')[-1]
    cmd =  'DPPP numthreads=80 msin=' + msin + ' msout.storagemanager=dysco '
    cmd += 'msout='+msout + ' '
    cmd += 'msout.writefullresflag=False '
    cmd += 'steps=[applyPA,applyBandpass,applyBeam,avg] '

    cmd += 'applyPA.type=applycal applyPA.correction=polalign '
    cmd += 'applyPA.parmdb=Band_PA.h5 '

    cmd += 'applyBandpass.type=applycal applyBandpass.correction=bandpass '
    cmd += 'applyBandpass.parmdb=Band_PA.h5 '
    cmd += 'applyBandpass.updateweights=True '

    cmd += 'applyBeam.type=applybeam '
    cmd += 'applyBeam.updateweights=True '
    cmd += 'applyBeam.usechannelfreq=False ' # because we work on a single SB

    cmd += f'avg.type=averager avg.freqstep={freqstep} '

    print(cmd)
    subprocess.call(cmd, shell=True)

def initrun(LnumLoc):
    Lnum = LnumLoc.split('/')[-2]
    os.mkdir(Lnum)
    os.system(f'cp -r {glob.glob("*py")[0]} {Lnum}')
    os.system(f'cp -r /net/rijn/data2/groeneveld/largefiles/Band_PA.h5 {Lnum}')
    os.chdir(Lnum)
    if LnumLoc[0] == '/': #Absolute path
        tocopy = glob.glob(LnumLoc + '*msdemix')
    else:
        tocopy = glob.glob('../'+LnumLoc+'*msdemix')
    for cop in sorted(tocopy):
        print(f'Copying SB{cop.split("SB")[1][:3]}')
        os.system(f'cp -r {cop} .')
    
    target_source = find_skymodel()
    if target_source == '3c196':
        os.system(f'cp -r /net/bovenrijn/data1/groeneveld/software/prefactor/skymodels/3C196-pandey.skymodel .')
    elif target_source == '3c380':
        os.system(f'cp -r /net/bovenrijn/data1/groeneveld/software/prefactor/skymodels/3C380_8h_SH.skymodel 3C380-SH.skymodel')
    print(target_source)

def extract_directions(calibrator):
    filename = glob.glob('*MFS-image.fits')[0]
    # img = bdsf.process_image(filename, rms_box = (640,160), rms_map=True, thresh='hard', thresh_isl=10.0, thresh_pix=25.0)
    # img.write_catalog(outfile='regions_wsclean1.fits', bbs_patches='single', catalog_type='srl', clobber=True, format='fits')

    cmd = f'''python extract.py'''
    print(cmd)
    os.system(cmd)

def _run_demix(location):
    os.system(f'python3 {location}averageandemix.py {location}')

def pre_init(location):
    ncpu = 4 # Be patient...
    os.system(f'cp -r /net/rijn/data2/groeneveld/LoDeSS_files/prerun/*py {location}')
    os.system(f'cp -r /net/rijn/data2/groeneveld/LoDeSS_files/prerun/demix.sourcedb {location}')
    demix_pool = []

    for i in range(ncpu):
        proc = mp.Process(target=_run_demix, args=(location,))
        proc.start()
        demix_pool.append(proc)
        time.sleep(3*60) # Wait 3 minutes for the SB to load
    
    for proc in demix_pool:
        proc.join()

def generate_ddf_bashfile(msname):
    base_cmd = f'''export NUMEXPR_MAX_THREADS=96
echo $NUMEXPR_MAX_THREADS
DDF.py --Data-ChunkHours=0.5 --Debug-Pdb=never --Parallel-NCPU=32 --Cache-Dir ./ --Data-MS {msname} --Data-ColName DATA --Data-Sort 1 --Output-Mode Clean --Deconv-CycleFactor 0 --Deconv-MaxMinorIter 1000000 --Deconv-RMSFactor 2.0 --Deconv-FluxThreshold 0.0 --Deconv-Mode HMP --HMP-AllowResidIncrease 1.0 --Weight-Robust -0.5 --Image-NPix 8192 --CF-wmax 50000 --CF-Nw 100 --Beam-CenterNorm 1 --Beam-Smooth 1 --Beam-Model LOFAR --Beam-LOFARBeamMode A --Beam-NBand 1 --Beam-DtBeamMin 5 --Output-Also onNeds --Image-Cell 8.0 --Freq-NDegridBand 7 --Freq-NBand 7 --Mask-Auto 1 --Mask-SigTh 2.0 --GAClean-MinSizeInit 10 --GAClean-MaxMinorIterInitHMP 100000 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Weight-ColName WEIGHT_SPECTRUM --Output-Name run1 --DDESolutions-DDModeGrid AP --DDESolutions-DDModeDeGrid AP --RIME-ForwardMode BDA-degrid --Output-RestoringBeam 45.0 --DDESolutions-DDSols merged.h5:sol000/phase000+amplitude000 --Deconv-MaxMajorIter 8 --Deconv-PeakFactor 0.005 --Cache-Reset 1 #>> ddfacet-c0.log 2>&'''
    with open('cmd1.sh','w') as handle:
        handle.write(base_cmd)
    
    second_cmd = f'''export NUMEXPR_MAX_THREADS=96
echo $NUMEXPR_MAX_THREADS
DDF.py --Data-ChunkHours=0.5 --Debug-Pdb=never --Parallel-NCPU=32 --Cache-Dir ./ --Mask-External=run1.int.restored.fits.mask.fits --Predict-InitDicoModel=run1.01.DicoModel --Data-MS {msname} --Data-ColName DATA --Data-Sort 1 --Output-Mode Clean --Deconv-CycleFactor 0 --Deconv-MaxMinorIter 1000000 --Deconv-RMSFactor 2.0 --Deconv-FluxThreshold 0.0 --Deconv-Mode HMP --HMP-AllowResidIncrease 1.0 --Weight-Robust -0.5 --Image-NPix 8192 --CF-wmax 50000 --CF-Nw 100 --Beam-CenterNorm 1 --Beam-Smooth 1 --Beam-Model LOFAR --Beam-LOFARBeamMode A --Beam-NBand 1 --Beam-DtBeamMin 5 --Output-Also onNeds --Image-Cell 8.0 --Freq-NDegridBand 7 --Freq-NBand 7 --Mask-Auto 1 --Mask-SigTh 2.0 --GAClean-MinSizeInit 10 --GAClean-MaxMinorIterInitHMP 100000 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Weight-ColName WEIGHT_SPECTRUM --Output-Name run2 --DDESolutions-DDModeGrid AP --DDESolutions-DDModeDeGrid AP --RIME-ForwardMode BDA-degrid --Output-RestoringBeam 45.0 --DDESolutions-DDSols merged.h5:sol000/phase000+amplitude000 --Deconv-MaxMajorIter 8 --Deconv-PeakFactor 0.005 --Cache-Reset 1 #>> ddfacet-c1.log 2>&'''
    with open('cmd2.sh','w') as handle:
        handle.write(second_cmd)

    
def calibrator():
    missinglist = find_missing_stations()

    mslist = sorted(glob.glob('*msdemix'))
    comblist = [(ms,missinglist) for ms in mslist]
    for c in comblist:
        run(c)

    msnames = glob.glob('*corr*')
    msname = msnames[2].split('SB')[0]
    corrected_msnames = add_dummyms(msnames)
    outname = msname + 'concat.ms'
    retstr = '['
    for j in corrected_msnames:
        retstr += f'{j},'
    retstr = retstr[:-1] + ']'
    cmd = f'DPPP numthreads=80 msin={retstr} msout={outname} msout.storagemanager=dysco msout.writefullresflag=false msin.missingdata=true msin.orderms=false steps=[]'
    print(cmd)
    os.system(cmd)

    input_concat = glob.glob('*concat.ms')[0]
    skymodel = glob.glob('*skymodel')[0]
    sourcename = skymodel.split('-')[0]
    if sourcename == '3C380':
        # I am so sorry for this line
        sourcename = '3c380'
    
    cmd = f'''python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/runwscleanLBautoR.py --BLsmooth --ionfactor 0.02 --docircular --no-beamcor --skymodel={skymodel} --skymodelsource={sourcename} --soltype-list="['scalarphasediff','scalarphase','complexgain']" --solint-list="[4,1,8]" --nchan-list="[1,1,1]" --smoothnessconstraint-list="[0.6,0.3,1]" --imsize=4096 --uvmin=300 --stopafterskysolve --channelsout=24 --fitspectralpol=False --soltypecycles-list="[0,0,0]" --normamps=False --stop=1 --smoothnessreffrequency-list="[30.,20.,0.]" --doflagging=True --doflagslowphases=False --flagslowamprms=25 {input_concat}'''
    print(cmd)
    os.system(cmd)   

def target(calfile,target):
    '''
        To do:
        Find missing stations, apply Band+PA and concatenate (similar to calibrator)
        Next, apply calfile in three applycal steps ( do this in circular, and go back )
        Next, phaseshift + average to target (new MS)
        runwscleanLBAuto
        DPPP apply to original ms, and use wsclean to image

        phaseshift format: [xxx.xxdeg, yy.yydeg]
    '''
    os.system(f'cp -r {calfile} calibrator.h5')
    missinglist = find_missing_stations()

    mslist = sorted(glob.glob('*msdemix'))
    comblist = [(ms,missinglist) for ms in mslist]
    for c in comblist:
        run(c)

    msnames = glob.glob('*corr*')
    msname = msnames[2].split('SB')[0]
    corrected_msnames = add_dummyms(msnames)
    outname = msname + 'concat.ms'
    retstr = '['
    for j in corrected_msnames:
        retstr += f'{j},'
    retstr = retstr[:-1] + ']'
    cmd = f'DPPP numthreads=80 msin={retstr} msout={outname} msout.storagemanager=dysco msin.missingdata=true msin.orderms=false msout.writefullresflag=false steps=[]'
    print(cmd)
    os.system(cmd)

    # Go to circular ...
    os.system('cp -r /net/rijn/data2/groeneveld/LoDeSS_files/lin2circ.py .')
    os.system(f'python lin2circ.py -i {outname} -c DATA -o DATA_CIRC')

    # Now, apply the calfile

    cmd = f'DPPP msin={outname} msout=. steps=[ac1,ac2] msout.datacolumn=CALCORRECT_DATA_CIRC msin.datacolumn=DATA_CIRC '
    cmd += f'ac1.type=applycal ac1.parmdb=calibrator.h5 ac1.solset=sol000 ac1.correction=phase000 '
    cmd += f'ac2.type=applycal ac2.parmdb=calibrator.h5 ac2.solset=sol000 ac2.correction=amplitude000 '
    print(cmd)
    os.system(cmd)

    # ... and go back to linear

    os.system(f'python lin2circ.py -i {outname} -c CALCORRECT_DATA_CIRC -b -l CALCORRECT_DATA')

    # Phaseshift to target+average

    cmd = f'DPPP msin={outname} msin.datacolumn=CALCORRECT_DATA msout=phaseshifted_{outname} msout.storagemanager=dysco steps=[phaseshift,averager] '
    cmd += f'phaseshift.phasecenter={target} averager.freqstep=4 msout.writefullresflag=false '
    print(cmd)
    os.system(cmd)

    os.mkdir('target_cal')
    os.system(f'cp -r phaseshifted_{outname} target_cal/')
    os.chdir('target_cal')

    # Copy + calibrate
    os.system(f'cp -r phaseshifted_{outname} backup_phaseshifted_{outname}')
    generate_boxfile(target)
    cmd = f'''python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/runwscleanLBautoR.py --pixelscale 8 -b boxfile.reg --antennaconstraint="['core',None]" --BLsmooth --ionfactor 0.02 --docircular --startfromtgss --soltype-list="['scalarphasediffFR','tecandphase']" --solint-list="[6,1]" --nchan-list="[1,1]" --smoothnessconstraint-list="[1.0,0.0]" --uvmin=300 --channelsout=24 --fitspectralpol=False --soltypecycles-list="[0,0]" --normamps=False --stop=5 --smoothnessreffrequency-list="[30.,0]" --doflagging=True --doflagslowphases=False --flagslowamprms=25 phaseshifted_{outname}'''
    print(cmd)
    os.system(cmd)   

    # Make a direction independent image of the whole field
    os.chdir('..')
    os.mkdir('DI_image')
    os.system('cp -r target_cal/merged_selfcalcyle004* merged_target.h5')
    cmd = f'DPPP msin={outname} msout=DI_image/corrected_{outname} msin.datacolumn=CALCORRECT_DATA steps=[ac1,ac2] msout.writefullresflag=false msout.storagemanager=dysco '
    cmd += f'ac1.type=applycal ac1.parmdb=merged_target.h5 ac1.solset=sol000 ac1.correction=phase000 '
    cmd += f'ac2.type=applycal ac2.parmdb=merged_target.h5 ac2.solset=sol000 ac2.correction=amplitude000 '
    print(cmd)
    os.system(cmd)
    os.chdir('DI_image')
    
    wscleancmd = f'wsclean -no-update-model-required -minuv-l 80.0 -size 8192 8192 -reorder -parallel-deconvolution 2048 -weight briggs -0.5 -weighting-rank-filter 3 -clean-border 1 -parallel-reordering 4 -mgain 0.8 -fit-beam -data-column DATA -padding 1.4 -join-channels -channels-out 8 -auto-mask 2.5 -auto-threshold 0.5 -pol i -baseline-averaging 2.396844981071314 -use-wgridder -name image_000 -scale 8.0arcsec -niter 150000 corrected_{outname}'
    print(wscleancmd)
    os.system(wscleancmd)
    os.chdir('..')

    # Make a guesstimate of the regions
    os.mkdir('extract_directions')
    os.chdir('extract_directions')
    os.system(f'cp -r ../DI_image/image_000-MFS-image.fits .')
    os.system(f'cp -r /net/rijn/data2/groeneveld/LoDeSS_files/DI/extract.py .')
    os.system(f'cp -r /net/rijn/data2/groeneveld/LoDeSS_files/DI/split_rectangles.py .')
    extract_directions(target)
    os.system(f'python split_rectangles.py regions_ws1.reg')
    # I don't intend on making this part of a 'full' pipeline,
    # as the input of the user is critical in this step. Maybe
    # once we have found that the choice of the directions
    # does not significantly impacts the DD performance, we might

def dd_pipeline(location,boxes,nthreads,target):
    '''
        This pipeline requires boxes to be pre-determined, as this is 
        a difficult step to automize. Maybe in the future...
    '''
    boxes = os.path.abspath(boxes)
    os.chdir(location)
    os.mkdir('DD_cal')
    os.chdir('DD_cal')
    os.system(f'cp -r {boxes} ./rectangles')
    os.system(f'cp -r /net/rijn/data2/groeneveld/LoDeSS_files/DD/* .')
    os.system(f'cp -r ../DI_image/image_000-????-model.fits .')
    os.system(f'cp -r ../DI_image/*ms .')

    spawn_delay = 3600 # == 1 hr
    threadlist = []
    for ii in range(nthreads):
        t = threading.Thread(target=_run_sing,args=str(ii))
        t.daemon = True
        t.start()
        threadlist.append(t)
        time.sleep(spawn_delay)
    for t in threadlist:
        t.join()
    # Go back to the root directory
    os.chdir('../../')


def DDF_pipeline(location,direction):
    '''
        This pipeline starts off where the DD pipeline stops:
        it checks what the noise is for 
    '''
    os.system(f'cp -r /net/rijn/data2/groeneveld/LoDeSS_files/runwsclean.py .')
    os.chdir(location)
    if not os.path.isdir('DD_cal'):
        print("You need to perform DD calibration before running the facet-imaging pipeline")
        return 1
    else:
        pass
    os.chdir('DD_cal')
    os.system('python extract_results.py')
    find_rms()

    cmd = f'python h5_merger.py -out merged.h5 -in RESULTS/h5files/* --ms run_0/direction0/Dir0.peel.ms '
    if direction != None:
        crdlist = direction.lstrip('[').rstrip(']').split(',')
        crdlist = [crd.split('deg')[0] for crd in crdlist]
        crd = SkyCoord(*crdlist,unit = (u.deg,u.deg))
        radiancoord = str([crd.ra.to(u.radian).value,crd.dec.to(u.radian).value]).replace(' ','')
        cmd += f'--add_direction {radiancoord}'
    print(cmd)
    os.system(cmd)

    # now make facet imaging folder and generate the shell files
    msname = glob.glob('*ms')[0]
    os.chdir('../')
    os.mkdir('facet_imaging')
    os.chdir('facet_imaging')
    os.system(f'cp -r ../DD_cal/merged.h5 .')
    os.system(f'cp -r ../DD_cal/{msname} ./')
    generate_ddf_bashfile(msname)
    
    # And now run the two DDF runs
    os.system('bash cmd1.sh')
    os.system('MakeMask.py --RestoredIm run1.int.restored.fits --Th 3')
    os.system('bash cmd2.sh')


if __name__ == "__main__":
    parse = argparse.ArgumentParser(description='LoDeSS calibrator+target pipeline')
    parse.add_argument('location',help='Location of the downloaded+demixed data. For now, it is important that the final folder begins with L??????.',type=str)
    parse.add_argument('--cal_H5',help='H5 file from the calibrator source. This is used to make an initial correction', default=None)
    parse.add_argument('--direction',help='Direction to go to when using the target pipeline. Format: "[xxx.xxdeg,yyy.yydeg]"', default=None,type=str)
    parse.add_argument('--boxes', help='Folder with boxes, called DirXX. Needed for direction dependent calibration')
    parse.add_argument('--nthreads', default=6, help='Amount of threads to be spawned by DD calibration. 5 will basically fill up a 96 core node (~100 load avg)')
    parse.add_argument('--prerun', action = 'store_true', help='Do this if the folder contains raw .tar files instead of demixed folders. Untarring has to happen on the node itself - so from a performance POV this might not be a good choice.')
    parse.add_argument('--pipeline', help='Pipeline of choice', choices=['DD','DI_target','DI_calibrator','DDF'])

    res = parse.parse_args()

    location = res.location
    call = ' '.join(sys.argv)
    if not os.path.isfile('calls.log'):
        os.system('touch calls.log')
    with open('calls.log','a') as handle:
        handle.write('\n')
        handle.write(call)

    if res.prerun:
        pre_init(location)

    if res.pipeline=='DI_calibrator':
        initrun(location)
        calibrator()
    elif res.pipeline=='DD':
        # This step doesn't necessarily need a target
        dd_pipeline(location,res.boxes,res.nthreads,res.direction)
    elif res.pipeline=='DI_target':
        # This step absolutely needs a target
        calfile_abs = os.path.abspath(res.cal_H5)
        initrun(location)
        target(calfile_abs,res.direction)
    elif res.pipeline=='DDF':
        DDF_pipeline(location,res.direction)
