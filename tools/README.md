# Biogeochemical Flux Model Tools

A set of tools is here provided to perform different processing of BFM output and restart data:

- **bnmerge**: merge the output files produced by MPI parallel `BFM_NEMO` simulations when XIOS is not used
- **chlsat**: compute satellite-like Chlorophyll-a fields from 3D BFM outputs
- **standalone_diag**: create plots and HTML index for STANDALONE experiments

## bnmerge
This tool merges the output files produced by `BFM_NEMO` when run in parallel without the use of XIOS library.
The 1D BFM output files are remapped into NEMO 3D grid and all variables are copied in output file along with NEMO `tmask`.

**Setup**

Change the variable `COMPILER` in the `Makefile` with the name of the included compiler file from the Compilers directory. Run gmake. 

The tool requires NetCDF-4.3+ libraries. To debug (OMP parallel or serial) export variable DEBUG equal to yes, namely `$> export DEBUG=yes`

**Usage**

Data processing is controlled using the namelist `bnmerge.nml` (Check the example provided for explanations of the input parameters):

`$>./bnmerge.x -f bnmerge.nml`

A batch execution script (`bnmerge.sh`) is provided as an example to generate and use the reconstruction tool on CMCC clusters.

## chlsat
This code computes the chlorophyll concentration as seen by satellite considering:
1. the optical depth and a tolerance level (as described in Eq. 2 of Vichi et al., 2007)
2. the 1% light level
3. the 0.1% light level

Input files are the chlorophyll concentration (variable `Chla`) and the attenuation coefficient
(variable `xEPS`), both with the same number of time stamps, and the mask file.
It also allows to compute the attenuation coefficient using the BFM formula from Chl concentration, background attenuation and chl specific absorption but neglecting the contribution from inorganic suspended matter and detritus.

This tool also derives the integrated (gross and net) primary production down to 1% and 0.1% light level by setting the flag `compute_intpp` within the input namelist and providing the path to the files containing BFM output fields `ruPPYc` and `resPPYc`.

**Setup**

Change the variable `COMPILER` in the `Makefile` with the name of the included compiler file from the Compilers directory. Run gmake. The tool requires NetCDF-4.3+ libraries. 

**Usage**

Data processing is controlled using the namelist `chlsat.nml` (Check the example provided for explanations of the input parameters):

`$>./chlsat.x -f chlsat.nml`


## standalone_diag
Create plots of BFM output state variables and nutriets ratios for STANDALONE simualtions with HTML index file.

**Setup**

A working python environment for conda `bfm_standalone_diag` can be created using the provided `environment.yml` file with

`$> conda env create -f environment.yml`

**Usage**

Processing is done by `standalone_diag.py` that requires BFM STANDALONE output file as input:

`(bfm_standalone_diag)$> python standalone_diag.py <bfm_output.nc>`

The script will generate a new folder `standalone_html` in the same location of the output file that contains plots produced on the base of available output fields (`standalone_html/images` directory) and an HTML index file to visualize the figure. 


