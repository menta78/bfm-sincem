# Include file for compilation with gfortran

  NCDF_LIB= -L${NC_LIB}/lib -lnetcdff -lnetcdf 
  NCDF_INC= -I${NC_INC}/include

  FC=gfortran
  LD=gfortran

  AR=ar
  ARFLAGS=-r

  CPPFLAGS=

  # Production
  FFLAGS= -O3 -fdefault-real-8 -funroll-all-loops -fcray-pointer -ftree-vectorize -ffree-line-length-256 ${NCDF_INC}

  # Debug
  # FFLAGS= -g -fdefault-real-8 -funroll-all-loops -fcray-pointer -ftree-vectorize -ffree-line-length-256 -ffpe-trap=invalid,zero,overflow -fbounds-check -Wall -Wuninitialized ${NCDF_INC}
  
LDFLAGS= ${NCDF_INC} ${NCDF_LIB} 
