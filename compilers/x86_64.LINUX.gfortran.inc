# Include file for compilation with gfortran 4.2-4.6
# Note that the debugging option -g is default and can be removed
# you may also need the following netcdf libraries 
# -L/usr/local/lib -lhdf5_hl -lhdf5 -lz
  NCLIB=  -L${NETCDFF_LIB} -lnetcdff -L${NETCDF_LIB} -lnetcdf

  FC=mpif90
  LD=mpif90

  AR=ar
  ARFLAGS=-r

  CPPFLAGS=
  FFLAGS= -O2 -cpp -g -I${NETCDF_INC} -I${NETCDFF_INC} --free-line-length-160
  LDFLAGS= ${NCLIB}
