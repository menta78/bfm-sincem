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

  NCDF_DIR=${NETCDF}
  NCDF_LIB= -L${NCDF_DIR}/lib -lnetcdff -lnetcdf
  NCDF_INC= -I${NCDF_DIR}/include

  FC=ifort
  LD=ifort

  AR=ar
  ARFLAGS=-r

  CPPFLAGS=
  FFLAGS= -O3 -r8 -fp-model source ${NCDF_INC} -xSKYLAKE-AVX512 -qopt-zmm-usage=high -no-fma
  LDFLAGS= ${NCDF_INC} ${NCDF_LIB}
