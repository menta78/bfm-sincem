!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! DESCRIPTION
!   NetCDF I/O interface
!   
!    ncchlid   : Chl input file identifier
!    ncepsid   : attenuation coeff input file identifier
!    ncmaskid  : mesh_mask input file identifier
!    ncid      : output file identifier
!
! COPYING
!
!   Copyright (C) 2023 BFM System Team (bfm_st@cmcc.it)
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation.
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
subroutine read_create

  use netcdf
  use mod_chlsat
  implicit none
  integer           :: ncid, ncchlid, ncepsid, ncmaskid, errst
  integer           :: p,n,d,t
  integer           :: ndims, nVars, nGlobalAtts
  integer           :: IDtime,IDtimetmp,IDunlimdim,IDvartime
  integer           :: IDx, IDy, IDz, IDvar, IDtarget, IDatt
  integer           :: IDmask, IDe3t, IDchla, IDeps
  integer           :: IDgpp, IDrsp, ncrspid, ncgppid
  real, allocatable, dimension(:) :: time
  integer           :: npoints
  character(len = NF90_MAX_NAME) :: fname1,fname2,attname,timename,coorlabel
  integer,parameter    :: namlst=10
  namelist /chlsat_nml/ inp_dir,out_dir,out_fname,chla_fname,eps_fname,      &
  &                     chla_name,eps_name,mask_fname,tolerance,compute_eps, &
  &                     p_eps0,p_epsChla,gpp_fname,rsp_fname,gpp_name,rsp_name, &
  &                     compute_intpp,compute_chlsat,dimsname,coorname,opt_chlsat

     ! Reading directory names and file name specifications
     open(namlst,file=trim(cf_nml),action='read',status='old',err=99)
     read(namlst,nml=chlsat_nml,err=98)
     close(namlst)
   
     ! check chlsat option
     if (compute_chlsat) then
        if (trim(opt_chlsat) == "MEAN" .OR. trim(opt_chlsat) == "V07" .OR. trim(opt_chlsat) == "V07mod" ) then
           write(*,*) "Compute chlsat with "//trim(opt_chlsat)//" algorithm."
           write(*,*)
        else 
           write(*,*) "Unrecognized option for chlsat algorithm (opt_chlsat). Available options are: "
           write(*,*) " MEAN   - compute mean value over optical extinction depth "
           write(*,*) " V07    - use Eq. 2 from Vichi et al. 2007 (doi:10.1016/j.jmarsys.2006.03.014)"
           write(*,*) " V07mod - modified version of V07 using cumulative optical extinction "
           stop
        endif
     endif

     ! open input files, read dimensions and get data 
     fname1 = trim(inp_dir)//"/"//trim(chla_fname)
     write(*,*) "Open Chla input file: ",trim(fname1)

     call handle_err( nf90_open(path = fname1, mode = NF90_WRITE, ncid = ncchlid), &
     &                errstring="Error opening "//trim(chla_fname))
     ! EPS file
     fname2 = trim(inp_dir)//"/"//trim(eps_fname)
     if ((fname1 .eq. fname2) .or. compute_eps) then
        ncepsid = ncchlid
     else
        write(*,*) "Open eps input file: ",trim(fname2)
        call handle_err(nf90_open(path = fname2, mode = NF90_WRITE, ncid = ncepsid), &
        &               "Error opening "//trim(eps_fname))
     end if

     call handle_err( nf90_inquire(ncchlid, nDims, nVars, nGlobalAtts, IDunlimdim), "0" )
     call handle_err( nf90_inquire_dimension(ncchlid, IDunlimdim, name=timename, len=ntime), "1" )
     call handle_err( nf90_inq_dimid(ncchlid, TRIM(dimsname(1)), IDx), "2" )
     call handle_err( nf90_inq_dimid(ncchlid, TRIM(dimsname(2)), IDy), "3" )
     call handle_err( nf90_inq_dimid(ncchlid, TRIM(dimsname(3)), IDz), "4" )
     ! get dimensions value
     call handle_err( nf90_inquire_dimension(ncchlid, IDunlimdim, len = ntime), "1" )
     call handle_err( nf90_inquire_dimension(ncchlid, IDx, len = jpi), "5" )
     call handle_err( nf90_inquire_dimension(ncchlid, IDy, len = jpj), "6" )
     call handle_err( nf90_inquire_dimension(ncchlid, IDz, len = jpk), "7" )

     write(*,*) " Input file dimensions :"
     write(*,*) " X dim - name: ",   TRIM(dimsname(1)), ", size:",jpi
     write(*,*) " Y dim - name: ",   TRIM(dimsname(2)), ", size:",jpj
     write(*,*) " Z dim - name: ",   TRIM(dimsname(3)), ", size:",jpk
     write(*,*) " Unlim dim - name: ", TRIM(dimsname(4)), ", size:",ntime

     ! open mesh_mask file, get mask and vertical length scale
     fname1 = trim(mask_fname)
     call handle_err( nf90_open(path = fname1, mode = NF90_NOWRITE, ncid = ncmaskid), &
     &                errstring="Error opening file "//trim(mask_fname) )
     allocate(mask(jpi,jpj,jpk))
     call handle_err(nf90_inq_varid(ncmaskid, "tmask", IDmask), &
     &               errstring="Error inquiring tmask values")
     call handle_err(nf90_get_var(ncmaskid, IDmask,mask), &
     &               errstring="Error in getting tmask values")
     allocate(e3t(jpi,jpj,jpk))
     call handle_err(nf90_inq_varid(ncmaskid, "e3t_0", IDe3t), &
     &               errstring="Error inquiring e3t_0 values")
     call handle_err(nf90_get_var(ncmaskid, IDe3t,e3t), &
     &               errstring="Error in getting e3t_0 values")
     call handle_err( nf90_close(ncmaskid) )

     write(*,*) ""
     write(*,*) "Read I/O coordinates: ", TRIM(coorname(1))//" , ", TRIM(coorname(2))//" , " &
                                        , TRIM(coorname(3))//" , ", TRIM(coorname(4))
     coorlabel=TRIM(coorname(1))//" "//TRIM(coorname(2))
     ! get coordinates
     allocate(lon(jpi,jpj))
     call handle_err( nf90_inq_varid(ncchlid, TRIM(coorname(1)), IDvar) )
     call handle_err( nf90_get_var(ncchlid, IDvar, lon),errstring="variable: "//TRIM(coorname(1)) )
     allocate(lat(jpi,jpj))
     call handle_err( nf90_inq_varid(ncchlid, TRIM(coorname(2)), IDvar) )
     call handle_err( nf90_get_var(ncchlid, IDvar, lat), errstring="variable: "//TRIM(coorname(2)) )
     allocate(depth(jpk))
     call handle_err( nf90_inq_varid(ncchlid, TRIM(coorname(3)), IDvar) )
     call handle_err( nf90_get_var(ncchlid, IDvar, depth),errstring="variable: "//TRIM(coorname(3)) )

     ! get chl data
     allocate(chla(jpi,jpj,jpk,ntime))
     call handle_err(nf90_inq_varid(ncchlid, chla_name, IDchla), &
     &               errstring="Error inquiring chl values")
     call handle_err(nf90_get_var(ncchlid, IDchla,chla), &
     &               errstring="Error in getting chl values")

     ! get extinction coefficient
     allocate(eps(jpi,jpj,jpk,ntime))
     if (compute_eps) then
        ! approximate computation from chlorophyll using 
        ! background attenuation and specific Chla absorption
        ! It neglects suspended inorganic matter and detritus
        eps(:,:,:,:) = p_eps0 + p_epsChla*chla(:,:,:,:)
     else
        call handle_err(nf90_inq_varid(ncepsid, eps_name, IDeps), &
        &               errstring="Error inquiring xeps values")
        call handle_err(nf90_get_var(ncepsid, IDeps,eps), &
        &               errstring="Error in getting xeps values")
     end if

     ! compute the satellite-like chl
     ! this is always done because optical depths are computed here
     allocate(ezd_od(jpi,jpj,ntime))
     allocate(chlsat_od(jpi,jpj,ntime))
     call chlcalc(ONE+tolerance,ezd_od,chlsat_od)
     allocate(ezd_01(jpi,jpj,ntime))
     allocate(chlsat_01(jpi,jpj,ntime))
     call chlcalc(-log(0.01_RLEN),ezd_01,chlsat_01)
     allocate(ezd_001(jpi,jpj,ntime))
     allocate(chlsat_001(jpi,jpj,ntime))
     call chlcalc(-log(0.001_RLEN),ezd_001,chlsat_001)

     if (compute_intpp) then
        fname1 = trim(inp_dir)//"/"//trim(gpp_fname)
        call handle_err( nf90_open(path = fname1, mode = NF90_WRITE, ncid = ncgppid), &
        &                errstring="Error opening "//trim(gpp_fname))
        fname1 = trim(inp_dir)//"/"//trim(rsp_fname)
        call handle_err( nf90_open(path = fname1, mode = NF90_WRITE, ncid = ncrspid), &
        &                errstring="Error opening "//trim(rsp_fname))
        ! get GPP and RSP data
        allocate(gpp(jpi,jpj,jpk,ntime))
        call handle_err(nf90_inq_varid(ncgppid, gpp_name, IDgpp), &
        &               errstring="Error inquiring GPP values")
        call handle_err(nf90_get_var(ncgppid, IDgpp,gpp), &
        &               errstring="Error in getting GPP values")
        allocate(rsp(jpi,jpj,jpk,ntime))
        call handle_err(nf90_inq_varid(ncrspid, rsp_name, IDrsp), &
        &               errstring="Error inquiring RSP values")
        call handle_err(nf90_get_var(ncrspid, IDrsp,rsp), &
        &               errstring="Error in getting RSP values")
        ! compute the integrated gross and net PP
        allocate(gpp_01(jpi,jpj,ntime))
        allocate(npp_01(jpi,jpj,ntime))
        call intppcalc(-log(0.01_RLEN),gpp_01,npp_01)
        allocate(gpp_001(jpi,jpj,ntime))
        allocate(npp_001(jpi,jpj,ntime))
        call intppcalc(-log(0.001_RLEN),gpp_001,npp_001)
     end if

     ! create the output file
     call handle_err( nf90_create(trim(out_dir)//"/"//trim(out_fname), NF90_NOCLOBBER, ncid), &
     &                errstring="A file named "//trim(out_fname)//" already exists!" )
     ! Define the dimensions
     call handle_err( nf90_def_dim(ncid, TRIM(dimsname(1)) , jpi, IDx) )
     call handle_err( nf90_def_dim(ncid, TRIM(dimsname(2)) , jpj, IDy) )
     call handle_err( nf90_def_dim(ncid, TRIM(dimsname(3)) , jpk, IDz) )
     call handle_err( nf90_def_dim(ncid, TRIM(dimsname(4)) , NF90_UNLIMITED, IDtime) )
     ! define geographic variables
     call handle_err (nf90_def_var(ncid, TRIM(coorname(2)), NF90_REAL, (/ IDx, IDy /), IDtarget))
     call handle_err (nf90_put_att(ncid, IDtarget, "units", "degrees_north"))
     call handle_err (nf90_def_var(ncid, TRIM(coorname(1)), NF90_REAL, (/ IDx, IDy /), IDtarget))
     call handle_err (nf90_put_att(ncid, IDtarget, "units", "degrees_east"))
     call handle_err (nf90_def_var(ncid, TRIM(coorname(3)), NF90_REAL, (/ IDz /), IDtarget))
     call handle_err (nf90_put_att(ncid, IDtarget, "long_name", "depth_below_sea"))
     call handle_err (nf90_put_att(ncid, IDtarget, "units", "m"))
     call handle_err (nf90_put_att(ncid, IDtarget, "positive", "down"))
     call handle_err (nf90_put_att(ncid, IDtarget, "axis", "Z"))
     ! copy global attributes
     write(*,*) "Creating file:",trim(out_fname)//".nc"," containing",ntime,"time frames"
#ifdef DEBUG
        write(*,*) "Start copying global attributes ..."
#endif
     do IDatt=1,nGlobalAtts
        errst=nf90_inq_attname(ncchlid, NF90_GLOBAL, IDatt, name=attname)
        call handle_err( nf90_copy_att(ncchlid, NF90_GLOBAL, trim(attname), ncid, NF90_GLOBAL), &
        &                errstring="copying attribute "//trim(attname))
     end do
     ! copy time variable and attributes
     call handle_err( nf90_inq_varid(ncchlid, TRIM(timename), IDtimetmp), errstring="inquiring time var in "//fname1)
     allocate(time(ntime))
     call handle_err( nf90_get_var(ncchlid, IDtimetmp, time, start = (/ 1 /), count = (/ ntime /)) )
     call handle_err( nf90_def_var(ncid, TRIM(coorname(4)) , NF90_DOUBLE, (/ IDtime /), IDvartime) )
     call handle_err( nf90_copy_att(ncchlid, IDtimetmp, "units", ncid, IDvartime) , errstring='copying time units from source')
     errst = nf90_inquire_attribute(ncchlid, IDtimetmp,"calendar")
     if ( errst == nf90_noerr) &
        call handle_err( nf90_copy_att(ncchlid, IDtimetmp, "calendar", ncid, IDvartime), errstring='missinf calendar attr in time')

     if (compute_chlsat) then
     ! define the output variable: chlsat
     call handle_err( nf90_def_var(ncid, TRIM(chla_name)//"sat_od", NF90_REAL, (/ IDx, IDy, IDtime /), IDvar) )
     ! add attributes
     call handle_err( nf90_put_att(ncid, IDvar, "units", "mg Chla/m3") )
     call handle_err( nf90_put_att(ncid, IDvar, "long_name", "Satellite-like Chl (1st optical depth)") )
     call handle_err( nf90_put_att(ncid, IDvar, "optical_depth_tolerance", tolerance) )
     ! Add fill value
     call handle_err( nf90_put_att(ncid, IDvar, "_FillValue", NF90_FILL_REAL) )
     ! Add reference coordinate names (needed with CDO)
     call handle_err( nf90_put_att(ncid, IDvar, "coordinates", TRIM(coorlabel) ) )
     ! define the output variable: depth
     call handle_err( nf90_def_var(ncid, "PZdepth_od", NF90_REAL, (/ IDx, IDy, IDtime /), IDvar) )
     ! add attributes
     call handle_err( nf90_put_att(ncid, IDvar, "units", "m") )
     call handle_err( nf90_put_att(ncid, IDvar, "long_name", "Photic zone depth (1st optical depth)") )
     ! Add fill value
     call handle_err( nf90_put_att(ncid, IDvar, "_FillValue", NF90_FILL_REAL) )
     ! Add reference coordinate names (needed with CDO)
     call handle_err( nf90_put_att(ncid, IDvar, "coordinates", TRIM(coorlabel) ) )

     ! define the output variable: chlsat
     call handle_err( nf90_def_var(ncid, TRIM(chla_name)//"sat_01", NF90_REAL, (/ IDx, IDy, IDtime /), IDvar) )
     ! add attributes
     call handle_err( nf90_put_att(ncid, IDvar, "units", "mg Chla/m3") )
     call handle_err( nf90_put_att(ncid, IDvar, "long_name", "Satellite-like Chl (1% light)") )
     ! Add fill value
     call handle_err( nf90_put_att(ncid, IDvar, "_FillValue", NF90_FILL_REAL) )
     ! Add reference coordinate names (needed with CDO)
     call handle_err( nf90_put_att(ncid, IDvar, "coordinates", TRIM(coorlabel) ) )
     ! define the output variable: chlsat
     call handle_err( nf90_def_var(ncid, TRIM(chla_name)//"sat_001", NF90_REAL, (/ IDx, IDy, IDtime /), IDvar) )
     ! add attributes
     call handle_err( nf90_put_att(ncid, IDvar, "units", "mg Chla/m3") )
     call handle_err( nf90_put_att(ncid, IDvar, "long_name", "Satellite-like Chl (0.1% light)") )
     ! Add fill value
     call handle_err( nf90_put_att(ncid, IDvar, "_FillValue", NF90_FILL_REAL) )
     ! Add reference coordinate names (needed with CDO)
     call handle_err( nf90_put_att(ncid, IDvar, "coordinates", TRIM(coorlabel) ) )
     end if

     ! define the output variable: depth
     call handle_err( nf90_def_var(ncid, "PZdepth_01", NF90_REAL, (/ IDx, IDy, IDtime /), IDvar) )
     ! add attributes
     call handle_err( nf90_put_att(ncid, IDvar, "units", "m") )
     call handle_err( nf90_put_att(ncid, IDvar, "long_name", "Photic zone depth (1% light)") )
     ! Add fill value
     call handle_err( nf90_put_att(ncid, IDvar, "_FillValue", NF90_FILL_REAL) )
     ! Add reference coordinate names (needed with CDO)
     call handle_err( nf90_put_att(ncid, IDvar, "coordinates", TRIM(coorlabel) ) )
     ! define the output variable: depth
     call handle_err( nf90_def_var(ncid, "PZdepth_001", NF90_REAL, (/ IDx, IDy, IDtime /), IDvar) )
     ! add attributes
     call handle_err( nf90_put_att(ncid, IDvar, "units", "m") )
     call handle_err( nf90_put_att(ncid, IDvar, "long_name", "Photic zone depth (0.1% light)") )
     ! Add fill value
     call handle_err( nf90_put_att(ncid, IDvar, "_FillValue", NF90_FILL_REAL) )
     ! Add reference coordinate names (needed with CDO)
     call handle_err( nf90_put_att(ncid, IDvar, "coordinates", TRIM(coorlabel) ) )

     if (compute_intpp) then
        ! define the output variable: gpp
        call handle_err( nf90_def_var(ncid, "gpp_01", NF90_REAL, (/ IDx, IDy, IDtime /), IDvar) )
        ! add attributes
        call handle_err( nf90_put_att(ncid, IDvar, "units", "mg C/m2/day") )
        call handle_err( nf90_put_att(ncid, IDvar, "long_name", "Integrated Gross Primary Production (1% light)") )
        ! Add fill value
        call handle_err( nf90_put_att(ncid, IDvar, "_FillValue", NF90_FILL_REAL) )
        ! Add reference coordinate names (needed with CDO)
        call handle_err( nf90_put_att(ncid, IDvar, "coordinates", TRIM(coorlabel) ) )
        !
        call handle_err( nf90_def_var(ncid, "gpp_001", NF90_REAL, (/ IDx, IDy, IDtime /), IDvar) )
        ! add attributes
        call handle_err( nf90_put_att(ncid, IDvar, "units", "mg C/m2/day") )
        call handle_err( nf90_put_att(ncid, IDvar, "long_name", "Integrated Gross Primary Production (0.1% light)") )
        ! Add fill value
        call handle_err( nf90_put_att(ncid, IDvar, "_FillValue", NF90_FILL_REAL) )
        ! Add reference coordinate names (needed with CDO)
        call handle_err( nf90_put_att(ncid, IDvar, "coordinates", TRIM(coorlabel) ) )
   
        ! define the output variable: npp
        call handle_err( nf90_def_var(ncid, "npp_01", NF90_REAL, (/ IDx, IDy, IDtime /), IDvar) )
        ! add attributes
        call handle_err( nf90_put_att(ncid, IDvar, "units", "mg C/m2/day") )
        call handle_err( nf90_put_att(ncid, IDvar, "long_name", "Integrated Net Primary Production (1% light)") )
        ! Add fill value
        call handle_err( nf90_put_att(ncid, IDvar, "_FillValue", NF90_FILL_REAL) )
        ! Add reference coordinate names (needed with CDO)
        call handle_err( nf90_put_att(ncid, IDvar, "coordinates", TRIM(coorlabel) ) )
        !
        call handle_err( nf90_def_var(ncid, "npp_001", NF90_REAL, (/ IDx, IDy, IDtime /), IDvar) )
        ! add attributes
        call handle_err( nf90_put_att(ncid, IDvar, "units", "mg C/m2/day") )
        call handle_err( nf90_put_att(ncid, IDvar, "long_name", "Integrated Net Primary Production (0.1% light)") )
        ! Add fill value
        call handle_err( nf90_put_att(ncid, IDvar, "_FillValue", NF90_FILL_REAL) )
        ! Add reference coordinate names (needed with CDO)
        call handle_err( nf90_put_att(ncid, IDvar, "coordinates", TRIM(coorlabel) ) )
     end if

     ! close the input netcdf files
     call handle_err( nf90_close(ncchlid) )
     if (ncepsid /= ncchlid) call handle_err( nf90_close(ncepsid) )

     ! exit definition mode
     call handle_err( nf90_enddef(ncid), errstring="while exiting definition mode")

     ! write time values
     call handle_err( nf90_put_var(ncid, IDvartime, time, &
     &               start = (/ 1 /), count = (/ ntime /)), errstring="variable: time")
     ! Latitude
     call handle_err( nf90_inq_varid(ncid, TRIM(coorname(2)), IDvar), &
     &                       errstring="inquiring variable: lat")
     call handle_err( nf90_put_var(ncid, IDvar, real(lat,4), start = (/ 1, 1 /),     &
     &                      count = (/ jpi, jpj/)), errstring="Writing: lat")
     ! Longitude
     call handle_err( nf90_inq_varid(ncid, TRIM(coorname(1)), IDvar), &
     &                       errstring="inquiring variable: lon")
     call handle_err( nf90_put_var(ncid, IDvar, real(lon,4), start = (/ 1, 1 /),     &
     &                      count = (/ jpi, jpj/)), errstring="Writing: lon")
     ! Depth levels
     call handle_err(nf90_inq_varid(ncid, TRIM(coorname(3)), IDvar))
     call handle_err(nf90_put_var(ncid, IDvar, real(depth,4)),errstring="Writing: depth")

     if (compute_chlsat) then
        ! Write chl
        call handle_err( nf90_inq_varid(ncid, TRIM(chla_name)//"sat_od", IDvar), &
        &                       errstring="inquiring variable: chlsat")
        call handle_err( nf90_put_var(ncid, IDvar, real(chlsat_od,4), start = (/ 1, 1, 1 /),     &
        &                      count = (/ jpi, jpj, ntime/)), errstring="Writing: chlsat")
        call handle_err( nf90_inq_varid(ncid, TRIM(chla_name)//"sat_01", IDvar), &
        &                       errstring="inquiring variable: chlsat")
        call handle_err( nf90_put_var(ncid, IDvar, real(chlsat_01,4), start = (/ 1, 1, 1 /),     &
        &                      count = (/ jpi, jpj, ntime/)), errstring="Writing: chlsat")
        call handle_err( nf90_inq_varid(ncid, TRIM(chla_name)//"sat_001", IDvar), &
        &                       errstring="inquiring variable: chlsat")
        call handle_err( nf90_put_var(ncid, IDvar, real(chlsat_001,4), start = (/ 1, 1, 1 /),     &
        &                      count = (/ jpi, jpj, ntime/)), errstring="Writing: chlsat")
        ! Write depth
        call handle_err( nf90_inq_varid(ncid, "PZdepth_od", IDvar), &
        &                       errstring="inquiring variable: chlsat")
        call handle_err( nf90_put_var(ncid, IDvar, real(ezd_od,4), start = (/ 1, 1, 1 /),     &
        &                      count = (/ jpi, jpj, ntime/)), errstring="Writing: chlsat")
     end if

     call handle_err( nf90_inq_varid(ncid, "PZdepth_01", IDvar), &
     &                       errstring="inquiring variable: chlsat")
     call handle_err( nf90_put_var(ncid, IDvar, real(ezd_01,4), start = (/ 1, 1, 1 /),     &
     &                      count = (/ jpi, jpj, ntime/)), errstring="Writing: chlsat")
     call handle_err( nf90_inq_varid(ncid, "PZdepth_001", IDvar), &
     &                       errstring="inquiring variable: chlsat")
     call handle_err( nf90_put_var(ncid, IDvar, real(ezd_001,4), start = (/ 1, 1, 1 /),     &
     &                      count = (/ jpi, jpj, ntime/)), errstring="Writing: chlsat")

     if (compute_intpp) then
        ! Write gpp
        call handle_err( nf90_inq_varid(ncid, "gpp_01", IDvar), &
        &                       errstring="inquiring variable: gpp_01")
        call handle_err( nf90_put_var(ncid, IDvar, real(gpp_01,4), start = (/ 1, 1, 1 /),     &
        &                      count = (/ jpi, jpj, ntime/)), errstring="Writing: gpp_01")
        call handle_err( nf90_inq_varid(ncid, "gpp_001", IDvar), &
        &                       errstring="inquiring variable: gpp_001")
        call handle_err( nf90_put_var(ncid, IDvar, real(gpp_001,4), start = (/ 1, 1, 1 /),     &
        &                      count = (/ jpi, jpj, ntime/)), errstring="Writing: gpp_001")
        ! Write npp
        call handle_err( nf90_inq_varid(ncid, "npp_01", IDvar), &
        &                       errstring="inquiring variable: npp_01")
        call handle_err( nf90_put_var(ncid, IDvar, real(gpp_01,4), start = (/ 1, 1, 1 /),     &
        &                      count = (/ jpi, jpj, ntime/)), errstring="Writing: npp_01")
        call handle_err( nf90_inq_varid(ncid, "npp_001", IDvar), &
        &                       errstring="inquiring variable: npp_001")
        call handle_err( nf90_put_var(ncid, IDvar, real(npp_001,4), start = (/ 1, 1, 1 /),     &
        &                      count = (/ jpi, jpj, ntime/)), errstring="Writing: npp_001")
     end if
     ! close
     call handle_err( nf90_close(ncid) )

#ifdef DEBUG
        write(*,*) "Output file created with ",ntime,"time frames"
#endif


     return

99   write (*,*) 'I could not open the namelist: ', trim(cf_nml); write (*,*)
     stop
98   write (*,*) 'I could not read the namelist: ', trim(cf_nml); write (*,*)
     stop


end subroutine read_create
