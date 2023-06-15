#NCDF_DIR     netcdf path
#NCDF_INC     netcdf include file
#NCDF_LIB     netcdf library (also hdf5 with nc4)
#FC           Fortran compiler command
#FCFLAGS      Fortran compiler flags
#FFLAGS       Fortran 77 flags
#LD           Linker
#LDFLAGS      Linker flags
#AR           Assembler
#ARFLAGS      Assembler flags
#OTHERFLAGS   User defined additional flags
#NOTE: in athena you have to load module: NETCDF/netcdf-4.3

NCDF_LIB= $(shell nc-config --libs) $(shell nf-config --flibs )
NCDF_INC= -I$(shell nf-config --includedir )

FC=ifort
LD=ifort

AR=ar
ARFLAGS=-r

CPPFLAGS=
# production
FFLAGS= -O3 -r8 -xHost -fp-model source -traceback -qopt-zmm-usage=high -no-fma ${NCDF_INC}

# debug
# FFLAGS= -g -O0 -r8 -xHost -fp-model source -traceback -fp-stack-check -check all -fpe0 ${NCDF_INC}

LDFLAGS= ${NCDF_INC} ${NCDF_LIB}

