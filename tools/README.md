## Biogeochemical Flux Model Tools

A set of tools is here provided to perform different processing of BFM output and restart data

- bnmerge: merge the output files produced by MPI parallel `BFM_NEMO` simulations when XIOS is not used
- chlsat: compute satellite-like Chlorophyll-a fields from 3D BFM outputs
- standalone_diag: create plots and HTML index for STANDALONE experiments

### bnmerge
This tool merges the output files produced by `BFM_NEMO` when run in parallel.
The 1D BFM output files are remapped into NEMO 3D grid.
All variables are copied and the output file contains also the Tmask and lat lon arrays.

**Usage**
bnmerge is controlled by the namelist bnmerge.nml (see bnmerge.nml example)


**Compilation**
- To debug (OMP parallel or serial) export variable DEBUG equal to yes: export DEBUG=yes
- To execute load netcdf-4.3+ libraries

Change the variabl`COMPILER` in the `Makefile` with the name of the included compiler file from the Compilers directory. Run gmake.

**Excution Script (bnmerge.sh)**
Script will help to generate a bacth runscript for CMCC clusters

### chlsat
This code computes the chlorophyll concentration as seen by satellite considering:
1. the optical depth and a tolerance level, as described in Eq. 2 of Vichi et al. (2007)
2. the 1% light level
3. the 0.1% light level

Input files are the chlorophyll concentration (variable `Chla`) and the attenuation coefficient
(variable `xEPS`), both with the same number of time stamps, and the mask file.
It also allows to compute the attenuation coefficient using the BFM formula from Chl
concentration, background attenuation and chl specific absorption but neglecting
the contribution from inorganic suspended matter and detritus.
Check the ORCA2 mask to see in the history how to generate the mask file from the `mesh_mask`.
The standard `mesh_mask.nc` file can also be used directly.
Check the namelist `chlsat.nml` for explanations of the input parameters

This tool also computes the integrated primary production (gross and net) down to 1% and
0.1% light level by setting the flag `compute_intpp` and providing the paths to the files
containing BFM diagnostics `ruPPYc` and `resPPYc`.

### standalone_diag
Create plots of BFM output state variables and nutriets ratios for STANDALONE experiments.
