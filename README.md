<h1 align="center"> emit-sds-l1b </h1>

_NOTE - at this time the EMIT repositories are not supporting Pull Requests from members outside of the EMIT-SDS Team.  This is expected to change in March, and guidance on branches will be provided at that time. At present, public migration of this repository is a work-in-progress, and not all features are yet fully up to date.  See the **develop** branch - set as default - for the latest code._

Welcome to the EMIT Level 1b science data system repository.  To understand how this repository is linked to the rest of the emit-sds repositories, please see [the repository guide](https://github.jpl.nasa.gov/emit-sds/emit-main/wiki/Repository-Guide).

## Calibrating an image

The top level calibration utility is emitrdn.py.  To calibrate a raw EMIT data file to spectrally resolved radiance units, simply type:

```
python emitrdn.py [path_to_dark_file.raw] [path_to_configuration_file.json] [input_dn_file.raw]  [output_radiance.img]
```

This utility is configured at runtime by the JSON configuration files in the config/ subdirectory.  These files provide options for varying different components of the radiometric calibration process, including the spectral calibration, radiometric calibration, optical and electronic corrections.  It selects a default configuration file if none is specified.

Often instruments will have more than one operating mode corresponding to different integration time settings.  This changes the radiometric calibration of the sensor.  For example, the EMIT top level configuration file has two blocks associated with normal and half-integration-time settings.  At run-time, half integration time observations (also known as "gypsum mode" observations) are specified using the "--mode half" flag:

```
python emitrdn.py --mode half [path_to_dark_file.raw] ...
```

Other command line parameters control verbosity and parallelism.

## Contents

The data/ subdirectory contains all calibration data files.  Files are organized by sensor; each sensor has its own subdirectory. 

The utils/ subdirectory contains library code for executing the steps of the calibration process.  Each can be done independently by using the following files as executables:

- utils/emit2dark.py - average a dark data sequence to create a "dark frame" for calibration
- utils/darksubtract.py - subtract the dark level indicated by a dark frame from the reference observation
- utils/pedestal.py - fix electronic pedestal shift.
- utils/fixlinearity - correct radiometric response for nonlinearity.
- utils/fixbad - infer values for bad detector elements.
- utils/fixscatter - correct blurring caused by extended spatial and spectral response functions.
- utils/fixghost - fix the symetric spatial Dyson ghost artifact.

The scripts/ subdirectory contains all of the scripts used to generate calibration files from laboratory data. Files are organized in subdirectories, one per sensor. These contain hard-coded filepaths that make them non-general; they are provided for transparency and as a template for any user who wants to adapt these analyses to their own datasets.

## Dependencies

The code uses numpy/scipy, matplotlib, scikit learn/image, and a variety of other public packages.  It uses the ray library for multicore parallelism.

