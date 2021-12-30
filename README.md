# emit-sds-l1b

Welcome to the EMIT Level 1b science data system repository.  To understand how this repository is linked to the rest of the emit-sds repositories, please see [the repository guide](https://github.jpl.nasa.gov/emit-sds/emit-main/wiki/Repository-Guide).

The top level calibration utility is emitrdn.py.  To calibrate a raw EMIT data file to spectrally resolved radiance units, simply type:

> python emitrdn.py --dark_file [path_to_dark_file.raw]  [input_dn_file.raw]  [output_radiance.img]

This utility is configured at runtime by the JSON configuration files in the config/ subdirectory.  These files provide options for varying different components of the radiometric calibration process, including the spectral calibration, radiometric calibration, optical and electronic corrections.  It selects a default configuration file if none is specified.

The data/ subdirectory contains all calibration data files.

The utils/ subdirectory contains library code for executing the steps of the calibration process.  Each can be done independently by using the following files as executables:

- utils/emit2dark.py - average a dark data sequence to create a "dark frame" for calibration
- utils/darksubtract.py - subtract the dark level indicated by a dark frame from the reference observation
- utils/pedestal.py - fix electronic pedestal shift.
- utils/fixlinearity - correct radiometric response for nonlinearity.
- utils/fixbad - infer values for bad detector elements.
- utils/fixscatter - correct blurring caused by extended spatial and spectral response functions.
- utils/fixghost - fix the symetric spatial Dyson ghost artifact.

To see the calling syntax for these utilities, along with example calibration data files for each, see script/quickcorrect.py.  This script chains the corrections together, emulating the action of emitrdn.py
