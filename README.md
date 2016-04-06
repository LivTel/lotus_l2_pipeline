lotus\_l2\_pipeline
=================

# Overview

This package provides a reduction pipeline to reduce data taken with the LOTUS spectrograph. Its operation is similar 
to that of **frodo-l2-pipeline**, except extraction and curvature correction operations have been adjusted/added to cater 
for the nature of long-slit rather than fibre-fed data. 

The core procedures (`src/`) have been written in C, and are mostly ported from the FRODOSpec codebase. The wrapper for the 
pipeline, `scripts/L2_exec.py`, is written in Python2.

# Installation

1. Clone the repository

2. Edit `scripts/L2_setup` and set the `L2_BASE_DIR` variable to the root directory of the repository.

3. Set up the environment by sourcing the `scripts/L2_setup` file

4. Edit `src/Makefile` and set the `LIBS` and `INCLUDES` parameters accordingly. These paths need to include locations of the GSL and CFITSIO headers/libraries. 

5. Make the binaries and library: `src/make all`

n.b. for LT operations, `L2_BASE_DIR`, `LIBS` and `INCLUDES` should already have a commented option for `lt-qc`.

# Configuration Files

Configuration file lookup tables (`config.tab`) are kept in `config/lookup_tables`. The file provides a record of: config file location (relative to `config/configs`), binning mode, date active from, time active from, date active to and time active to. The last two fields can be replaced by a single value, "now", signifying that the file can be used up to the current date/time, e.g.

`010715/config_4x4.ini	4x4	01/07/15	12:00:00	now`

This table should be appended to for any change in configuration.

The actual configuration files themselves are normally kept in `config/configs/[**DATE**]/` but obviously could be kept anywhere, as long as the location is correctly specified in the corresponding lookup table. The bulk of the config parameters may be understood by looking at the routine subheading's `man/` file. 

Other non-routine options include [**max\_curvature\_post\_cor**] and [**operations**]. The former sets the maximum curvature allowed post correction, in pixels. The latter defines the extension order of the resulting output file. Options are **LSS\_NONSS**, **SPEC\_NONSS** and **SPEC\_SS**.


# Arc Calibration Solution

The arc solution for LOTUS is a 2nd order hard-coded polynomial. Its coefficients can be found in `config/ARC_SOLUTION.dat`.

# Invoking the Pipeline

The parameters available at run-time can be found by invoking the pipeline with the help (`--h`) flag:

`[rmb@rmb-tower scripts]$ python L2_exec.py --h`

The `--f` flag specifies the target frame to be reduced, the `--r` flag specifies a reference file that is used to 
remove spectral curvature; it should be of a bright target taken at roughly the same alt/az as the target observation 
(reference and target frames must suffer the same flexure, otherwise this correction doesn't make sense).
The `--c` flag is only used if the spectrum position is to be automatically found, rather than using hardcoded limits - 
see `man/LOTUS_RED_CLIP` for more details about the [**force\_***] parameters.

The pipeline can perform a check the validity of the reference file for spectral curvature correction without doing a 
full reduction by using the `--rc` flag. Clobbering existing files can be allowed with `--o`.

# Known Issues

Artifacts are introduced at spectral curvature correction stage when the bias is non-zero. This was a problem in earlier 
data, and RJS has effectively disabled it for now (ln 655-657, L2_exec.py).

As there is only a single static solution for the arc calibration, if different binning modes are used, this stage will 
not work.

