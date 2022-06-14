MODULE BFM_PARAM3D

   USE dom_oce
   USE par_oce
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O library
   USE fldread

   USE mem, ONLY : ESS
   USE api_bfm
   USE init_var_bfm_local

   IMPLICIT NONE
   
   PRIVATE
   TYPE(FLD), ALLOCATABLE, SAVE, DIMENSION(:) :: sf_param3d  ! structure of input fields (file informations, fields read)
   LOGICAL    :: moduleInitialized = .FALSE.
   INTEGER  , SAVE      ::   jfld   = 1
   INTEGER  , SAVE      ::   jf_ess  = 1    ! index of suspended sediment
   LOGICAL              ::   ln_ess3d = .FALSE.
   
   PUBLIC param3d_dyn

CONTAINS

   SUBROUTINE param3d_dyn( kt )
      INTEGER, INTENT(IN) :: kt

      IF (.NOT. moduleInitialized) THEN
         CALL param3d_dyn_init
         moduleInitialized = .TRUE.
      END IF  
 
      IF (ln_ess3d) THEN
         WRITE(numout,*) "  BFM: loading dynamic ess data"
         CALL fld_read( kt, 1, sf_param3d )      !=  read data at kt time step   ==!
         ESS = pack(sf_param3d(jf_ess)%fnow(:,:,:)*tmask(:,:,:), SEAmask)
      END IF
   END SUBROUTINE



   SUBROUTINE param3d_dyn_init
      INTEGER NMLUNIT, ifpr, idv, idimv, ierr, ierr0, ierr1
      TYPE(FLD_N), DIMENSION(jfld) ::  slf_d         ! array of namelist informations on the fields to read

      CHARACTER(len=100)            ::  cn_dir        !   Root directory for location of core files
      TYPE(FLD_N)                   ::  sn_ess        !   informations about the fields to be read
      NAMELIST/bfm_param3d/cn_dir, ln_ess3d, sn_ess 

      NMLUNIT=GetLun()
      open(NMLUNIT,file=bfm_init_fname,action='read',status='old')
      cn_dir = './'
      ln_ess3d = .FALSE.
      rewind(NMLUNIT)
      read(NMLUNIT,nml=bfm_param3d,iostat=ierr)
      close(NMLUNIT)

      IF (ln_ess3d) THEN
         slf_d(jf_ess) = sn_ess
         ALLOCATE( sf_param3d(jfld), STAT=ierr )         ! set sf structure
         IF( ierr > 0 )  THEN
            CALL ctl_stop( 'param3d_dyn_init: unable to allocate sf structure' )
            RETURN
         ENDIF
         CALL fld_fill( sf_param3d, slf_d, cn_dir, 'param3d_dyn_init', 'Data in file', 'param3d' )
   
         !
         ! Open file for each variable to get his number of dimension
         DO ifpr = 1, jfld
            CALL fld_clopn( sf_param3d(ifpr), nyear, nmonth, nday )
            idv   = iom_varid( sf_param3d(ifpr)%num , slf_d(ifpr)%clvar )        ! id of the variable sdjf%clvar
            idimv = iom_file ( sf_param3d(ifpr)%num )%ndims(idv)                 ! number of dimension for variable sdjf%clvar
            IF( sf_param3d(ifpr)%num /= 0 )   CALL iom_close( sf_param3d(ifpr)%num ) ! close file if already open
            ierr1=0
            IF( idimv == 3 ) THEN    ! 2D variable
                                         ALLOCATE( sf_param3d(ifpr)%fnow(jpi,jpj,1), STAT=ierr0 )
               IF( slf_d(ifpr)%ln_tint ) ALLOCATE( sf_param3d(ifpr)%fdta(jpi,jpj,1,2), STAT=ierr1 )
            ELSE                     ! 3D variable
                                         ALLOCATE( sf_param3d(ifpr)%fnow(jpi,jpj,jpk), STAT=ierr0 )
               IF( slf_d(ifpr)%ln_tint ) ALLOCATE( sf_param3d(ifpr)%fdta(jpi,jpj,jpk,2), STAT=ierr1 )
            ENDIF
            IF( ierr0 + ierr1 > 0 ) THEN
               CALL ctl_stop( 'param3d_dyn_init : unable to allocate sf_param3d array structure' )   ;   RETURN
            ENDIF
         END DO
      END IF

   END SUBROUTINE

END MODULE BFM_PARAM3D
