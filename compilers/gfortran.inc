# Include file for compilation with gfortran
# Note that the debugging option -g is default and can be removed
# you may also need the following netcdf libraries 
# -L/usr/local/lib -lhdf5_hl -lhdf5 -lz

NCDF_LIB= $(shell nc-config --libs) $(shell nf-config --flibs )
NCDF_INC= -I$(shell nf-config --includedir )

FC=gfortran
LD=gfortran

AR=ar
ARFLAGS=-r

CPPFLAGS=

# production
FFLAGS= -O3 -fdefault-real-8 -funroll-all-loops -fcray-pointer -ftree-vectorize -ffree-line-length-256 ${NCDF_INC}

# debug
# FFLAGS= -g -fdefault-real-8 -funroll-all-loops -fcray-pointer -ftree-vectorize -ffree-line-length-256 -ffpe-trap=invalid,zero,overflow -fbounds-check -Wall -Wuninitialized ${NCDF_INC}
  
LDFLAGS= ${NCDF_INC} ${NCDF_LIB} 
