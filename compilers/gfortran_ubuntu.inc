# Include file for compilation with gfortran 4.2-4.6
# Note that the debugging option -g is default and can be removed
# you may also need the following netcdf libraries 
# -L/usr/local/lib -lhdf5_hl -lhdf5 -lz

  NCDF_LIB= -L/usr/lib/x86_64-linux-gnu -lnetcdff -lnetcdf 
  NCDF_INC= -I/usr/include

  FC=gfortran
  LD=gfortran

  AR=ar
  ARFLAGS=-r

  CPPFLAGS=
  FFLAGS= -O3 -fdefault-real-8 -funroll-all-loops -fcray-pointer -ftree-vectorize -ffree-line-length-256 ${NCDF_INC}
  LDFLAGS= ${NCDF_INC} ${NCDF_LIB} 
