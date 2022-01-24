# The LOFAR Decameter Sky Survey Pipeline (DeSSPi)

The LOFAR Decameter Sky Survey (LoDeSS) is a survey of the northern sky, performed at frequencies between 14 and 30 MHz. Due to the unique nature of this 
frequency range, we created a new pipeline for calibration. In addition, we include some other programs that help with finding the right parameters.

## Survey Setup

LoDeSS consists of around 360 overlapping pointings covering the entire northern sky at declinations above 20°. Each pointing is observed for 5 hours,
which do not necessarily need to be in one night: especially for objects with a low declination, pointings can be split into several nights.
Each night, five pointings and one calibrator are recorded. The data is available on the long-term archive (LTA).

## Rough overview of calibration

Calibrating a pointing consists of roughly four steps, each accounting for a specific set of effects:

| | *Calibration step* | *Effect* |
|----|----|----|
|1. | Apply pre-determined bandpass response and polarization alignment corrections | Effects that are completely time-independent (only change on the order of ~years)|
|2. | Apply corrections from the primary calibrator | Corrections on (1.), and effects that do not or only slightly depend on the direction of your pointing (especially clock effects)|
|3. | Calibrate on a bright source in the pointing (in-field calibrator) | Effects that depend on time and the direction, but give a good "starting point" for the next step. In particular, we correct for Faraday Rotation and TEC in the ionosphere|
|4. | Break the field up in facets, and calibrate for each facet individually | Effects that strongly change in direction, in this case we only calibrate for TEC effects through the field|

These steps are always performed in this specific sequence: generally this means that we begin with correcting for the most "broad" or "general" effects, and 
consequently move to more specific effects.

## Requirements
a lot of patience

## Workflow

In this section, I describe the workflow that I use. This is not necessarily the "correct" workflow, but just something that seems to work for me - your mileage
may vary.

- Select the pointings that you want to reduce, and find the corresponding L-numbers and the L-numbers of the calibrators from the LTA. Create a stage request
- Use [preprocessor.py](prerun/preprocessor.py) to quickly download and untar many different objects using parallel downloads from the LTA
- Use [infield_finder.py](toolbox/infield_finder.py) to get an overview of each pointing that you want to reduce. Make sure that each pointing is sufficiently far away from A-team sources (CasA, CygA and to a lesser extend TauA or HerA). Identify a good, bright in-field calibrator by clicking on it and verify that it is not turning over.
- Run: ` [LoDeSS.py](LoDeSS.py) --pipeline DI_calibrator --prerun ` to calibrate the calibrator (`--prerun` will also demix the data, you probably should leave this enabled)
- Check the FITS image and the waterfall plots to see if there are no major problems (severe scintillation, A-team interference or particularly strong ionospheric effects)
- Run: `LoDeSS.py --pipeline DI_target --direction "(xx.xxxx,yy.yyyy)" --cal_H5 /path/to/cal/h5/1 /path/to/cal/h5/2 --prerun /path/to/folder/with/ms1 /path/to/folder/with/ms2` . This will start the direction independent pipeline. Note that you need to give Lodess.py the calibrator h5s from the previous step
- Check if the in-field calibrator has been succesfully calibrated (and not diverged)
- Check if the DI-image is succesful, and off-calibrator sources are still recoverable (good way to see if FR corrections were properly applied)
- Check the boxes placed by extract_directions whether or not they are properly placed
- STILL WORK IN PROGRESS! Run: `LoDeSS.py --pipeline DD --rectangles /path/to/boxes --nthreads 6 /path/to/ms1 /path/to/ms2` to start the direction dependent pipeline - `--nthreads` should be drastically lowered on machines with less cpu/memory availability (this number works well on the Leiden LOFAR nodes, which contain 96 cpus and ~500G of RAM)
- STILL WORK IN PROGRESS! Run: `LoDeSS.py --pipeline DDF /path/to/ms1 /path/to/ms2` to use DDF to make an image.
- enjoy your nice picture
