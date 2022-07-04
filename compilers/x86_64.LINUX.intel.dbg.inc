        NCLIB=  -L${NETCDFF_LIB} -lnetcdff -L${NETCDF_LIB} -lnetcdf

        FC=mpiifort $(OPENMP_FLAG)
        LD=mpiifort $(OPENMP_FLAG)
	AR=ar
	ARFLAGS=-r
	CPPFLAGS=
        FFLAGS_DBG=  -O0 -g -traceback -fp-stack-check -check all -fpe0 
	FFLAGS= ${FFLAGS_DBG}  -I${NETCDFF_INC}
	LDFLAGS=${NCLIB}
