# lodess

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
- Use [preprocessor.py](prerun/preprocessor.py) 
