        NCLIB= -L${NETCDF_LIB} -lnetcdff -lnetcdf

        FC=mpif90 $(OPENMP_FLAG)
        LD=mpif90 $(OPENMP_FLAG)
	AR=ar
	ARFLAGS=-r
	CPPFLAGS=
        FFLAGS_OPT= -fno-math-errno -unroll=3 -opt-subscript-in-range -align all  -cpp -heap-arrays
	FFLAGS= ${FFLAGS_OPT}  -I${NETCDF_INC}
	LDFLAGS=${NCLIB}
