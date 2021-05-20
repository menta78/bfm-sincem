# Include file for compilation with gfortran 4.2-4.6
# Note that the debugging option -g is default and can be removed
# you may also need the following netcdf libraries 
# -L/usr/local/lib -lhdf5_hl -lhdf5 -lz
  NCLIB=  -L${NETCDFF_LIB} -lnetcdff -L${NETCDF_LIB} -lnetcdf

  FC=gfortran
  LD=gfortran

  AR=ar
  ARFLAGS=-r

  CPPFLAGS=
  FFLAGS= -cpp -g -fdefault-real-8 -funroll-all-loops -fcray-pointer -ftree-vectorize ${NETCDFF_INC}
  LDFLAGS= ${NCLIB}
