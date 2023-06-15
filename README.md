## What is the BFM

The Biogeochemical Flux Model (BFM) is numerical model to simulate the dynamics of major biogeochemical properties in lower trophic levels of marine ecosystems.

The BFM is open source software freely available under the GNU Public License. 
The model description and documentation is available on the BFM website http://www.bfm-community.eu.

In the `doc/` folder of the repository are also available the manuals and logos.


#### Members of the BFM agreement

The development of the BFM is supported by a system team through a formal agreement (see https://bfm-community.github.io/www.bfm-community.eu/bfm-consortium/). 
The agreement establishes a program for model development that is discussed and prepared every year by the participants. Model developments are inserted by the system team into the different releases of the model. 
The agreement is open to entrance of more contributors and more information can be addressed to the **BFM system Team** at bfm_st(at)cmcc.it

Current members of the consortium are: 
**CMCC, OGS, UNIBO, UCT, SYKE**

##  Installation Requirements

Supported architectures:
- Linux
- Mac OSX (Darwin)

Software requirements:
- PERL (version 5.8.2 and above).
- FORTRAN 90/95 compiler. Under linux and Mac OSX the model can be currently compiled with gfortran (version 4.6 or higher) and ifort.
- NetCDF Libraries (https://www.unidata.ucar.edu/software/netcdf). It is mandatory that the library has been compiled with the same compiler used for the model compilation (as F90 netcdf module is used).
- GNU make. Makefile only works with GNU make, therefore substitute your system make or use an alias to ensure that the right one is set in case you are not on a linux machine.

## Get started with STANDALONE presets
Configuration and deployment of the model is done automatically by the script `bfm_configure.sh` found in the `build` directory.
```bash
$> cd $BFMDIR/build
[$BFMDIR/build]> ./bfm_configure.sh -h
```

The user-dependent options are set either through the command line of the script or by adjusting (or adding) an architecture file in directory `$BFMDIR/compilers`. Default file is `gfortran.inc` that automatically detects the path of NetCDF libraries. The standard GNU gmake variables are used for compiler and archiver names.

If preset is not specified, the **STANDALONE** zero-dimensional configuration is compiled by default with the command
```bash
[$BFMDIR/build]> ./bfm_configure.sh -gcd
```
(see the script help for details on the options) and if successful will produce the following output
```bash
.................................. Makefile is ready.
STANDALONE generation done!
Starting STANDALONE compilation...
STANDALONE compilation done!
Go to $BFMDIR/run/standalone and execute command:
./bfm_standalone.x
```
The configuration script creates in `$BFMDIR/run/standalone` the executable link, a copy of the namelists, and the running environment. 
To launch a simulation go to the experiment run folder and execute the following 
```bash
$> cd $BFMDIR/run/standalone
[$BFMDIR/run/standalone]> ./bfm_standalone.x
```
**TIP**: use `./bfm_standalone.x &> outputfile` to redirect the screen output messages to a file in bash.

Finally, realize an organized set of plots for BFM STANDALONE preset using the `standalone_diag` python script from the execution folder (usage details in tools README) :
```bash
[$BFMDIR/run/standalone]> python  ../../tools/standalone_diag/standalone_diag.py BFM_Standalone.nc`
```

A list of the available presets can be obtained using the command
```bash
[$BFMDIR/build]> ./bfm_configure.sh -P 
```
Refer to the model documentation for a description of the available presets.

To compile a specific preset use the following
```bash
[$BFMDIR/build]> ./bfm_configure.sh -gcd -p <preset_name> 
```

## Acknowledgment
BFM System Team gratefully acknowledge the work of the development team of the original ERSEM versions I and II under the lead of J. W. Baretta providing the ancestor model that the original version of the BFM has emerged from.
