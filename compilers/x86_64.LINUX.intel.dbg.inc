        NCLIB= -L${NETCDF_LIB} -lnetcdff -lnetcdf

        FC=mpif90 $(OPENMP_FLAG)
        LD=mpif90 $(OPENMP_FLAG)
	AR=ar
	ARFLAGS=-r
	CPPFLAGS=
        FFLAGS_DBG=  -O0 -g -traceback -fp-stack-check -check bounds -fpe0
	FFLAGS= ${FFLAGS_DBG}  -I${NETCDF_INC}
	LDFLAGS=${NCLIB}
