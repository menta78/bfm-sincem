        NCLIB=  -L${NETCDFF_LIB} -lnetcdff -L${NETCDF_LIB} -lnetcdf

        FC=mpiifort $(OPENMP_FLAG)
        LD=mpiifort $(OPENMP_FLAG)
	AR=ar
	ARFLAGS=-r
	CPPFLAGS=
        FFLAGS_OPT= -fno-math-errno -unroll=3 -qopt-subscript-in-range -align all  -cpp -heap-arrays
	FFLAGS= ${FFLAGS_OPT}  -I${NETCDFF_INC}
	LDFLAGS=${NCLIB}
