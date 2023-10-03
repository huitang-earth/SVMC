MODULE io_mod

  use readctrl_mod
  use readvegpara_mod
  use netcdf

  implicit none

contains

  subroutine ncdf_handle_error(status)
    use netcdf    
    integer status 

    write(*,*) "******* Fatal NETCDF error *******"
    write(*,*) trim(NF90_strerror(status))
    write(*,*) "*************** RIP **************" 
    stop
  end subroutine ncdf_handle_error

  subroutine check(status)
    use netcdf
    integer, intent ( in) :: status

    !print*, 'Check standard '
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  

  subroutine netCDF_prepareOUTPUT(filename, lon, lat, ntim)
    !Initialize a netCDF-file for emission fields
    use readctrl_mod
    use netcdf
    implicit none

    character(*), intent(in) :: filename
    character(80):: str_time
    character    :: adate*8,atime*6,timeunit*32
    integer     :: nc_id, status
    integer :: londim_id, latdim_id, timedim_id,lonvar_id, latvar_id, timevar_id, &
               pftdim_id, pftvar_id, gppvar_id, neevar_id, nppvar_id, trvar_id, arvar_id, &
               hrvar_id, srvar_id, laivar_id, scvar_id, stvar_id, evvar_id, travar_id, &
               smvar_id, smpvar_id, cyvar_id, abvar_id, tbvar_id, lcvar_id, rcvar_id, fpvar_id, &
               namedim_id, jmvar_id, vcvar_id, dpsivar_id, chivar_id, provar_id
    integer :: nx_lon=1, ny_lat=1, ntim
    integer :: yyyy,mm,dd,hh,mi,ss

    !Some vars for standard netcdf example
    real(8), dimension(1) :: lon
    real(8), dimension(1) :: lat

    write(adate,'(i8.8)') start_date_day
    write(atime,'(i6.6)') start_date_hour


    call check(nf90_create(filename, cmode = NF90_HDF5, ncid = nc_id) )
     
    !Define dimensions
    call check(NF90_DEF_DIM(nc_id, "lon", nx_lon, londim_id))
    call check(NF90_DEF_DIM(nc_id, "lat", ny_lat, latdim_id))
    call check(NF90_DEF_DIM(nc_id, "time", ntim, timedim_id))
    call check(NF90_DEF_DIM(nc_id, "pft", num_pft, pftdim_id))
    call check(NF90_DEF_DIM(nc_id, "namelen", 20, namedim_id))

    !Define variables
    !   status=NF90_DEF_VAR(nc_id,"lon",NF90_FLOAT, (/ londim_id /), lonvar_id)
    call check(nf90_def_var(nc_id, "lon", nf90_float, (/ londim_id /), lonvar_id))
    call check(nf90_def_var(nc_id, "lat", nf90_float, (/ latdim_id /), latvar_id))
    call check(nf90_def_var(nc_id, "time", nf90_float, (/timedim_id/), timevar_id))
    call check(nf90_def_var(nc_id, "pftname", nf90_char, (/namedim_id, pftdim_id/), pftvar_id))
    call check(nf90_def_var(nc_id, "GPP", nf90_float, (/londim_id,latdim_id,timedim_id/), gppvar_id))
    call check(nf90_def_var(nc_id, "NEE", nf90_float, (/londim_id,latdim_id,timedim_id/), neevar_id))
    call check(nf90_def_var(nc_id, "NPP", nf90_float, (/londim_id,latdim_id,timedim_id/), nppvar_id))
    call check(nf90_def_var(nc_id, "TotalResp", nf90_float, (/londim_id,latdim_id,timedim_id/), trvar_id))
    call check(nf90_def_var(nc_id, "AutoResp", nf90_float, (/londim_id,latdim_id,timedim_id/), arvar_id))
    call check(nf90_def_var(nc_id, "HeteroResp", nf90_float, (/londim_id,latdim_id,timedim_id/), hrvar_id))
    call check(nf90_def_var(nc_id, "SoilResp", nf90_float, (/londim_id,latdim_id,timedim_id/), srvar_id))
    call check(nf90_def_var(nc_id, "LAI", nf90_float, (/londim_id,latdim_id,timedim_id/), laivar_id))
    call check(nf90_def_var(nc_id, "soil_carbon_content", nf90_float, (/londim_id,latdim_id,timedim_id/), scvar_id))    
    call check(nf90_def_var(nc_id, "stomatal_conductance", nf90_float, (/londim_id,latdim_id,timedim_id/), stvar_id))   
    call check(nf90_def_var(nc_id, "Evap", nf90_float, (/londim_id,latdim_id,timedim_id/), evvar_id))   
    call check(nf90_def_var(nc_id, "Transp", nf90_float, (/londim_id,latdim_id,timedim_id/), travar_id))
    call check(nf90_def_var(nc_id, "SoilMoist", nf90_float, (/londim_id,latdim_id,timedim_id/), smvar_id))
    call check(nf90_def_var(nc_id, "SoilMoistPot", nf90_double, (/londim_id,latdim_id,timedim_id/), smpvar_id))
    call check(nf90_def_var(nc_id, "CropYield", nf90_float, (/londim_id,latdim_id,timedim_id/), cyvar_id))
    call check(nf90_def_var(nc_id, "AGB", nf90_float, (/londim_id,latdim_id,timedim_id/), abvar_id))
    call check(nf90_def_var(nc_id, "TotLivBiom", nf90_float, (/londim_id,latdim_id,timedim_id/), tbvar_id))
    call check(nf90_def_var(nc_id, "leaf_carbon_content", nf90_float, (/londim_id,latdim_id,timedim_id/), lcvar_id))
    call check(nf90_def_var(nc_id, "root_carbon_content", nf90_float, (/londim_id,latdim_id,timedim_id/), rcvar_id))
    call check(nf90_def_var(nc_id, "fPAR", nf90_float, (/londim_id,latdim_id,timedim_id/), fpvar_id))

    call check(nf90_def_var(nc_id, "Jmax", nf90_float, (/londim_id,latdim_id,timedim_id/), jmvar_id))
    call check(nf90_def_var(nc_id, "Vcmax", nf90_float, (/londim_id,latdim_id,timedim_id/), vcvar_id))
    call check(nf90_def_var(nc_id, "Chi", nf90_float, (/londim_id,latdim_id,timedim_id/), chivar_id))
    call check(nf90_def_var(nc_id, "Dpsi", nf90_float, (/londim_id,latdim_id,timedim_id/), dpsivar_id))
    call check(nf90_def_var(nc_id, "Profit", nf90_float, (/londim_id,latdim_id,timedim_id/), provar_id))

    !Attributes
    call check(NF90_PUT_ATT(nc_id, lonvar_id, "units", "degrees_east"))
    call check(NF90_PUT_ATT(nc_id, lonvar_id, "standard_name", "longitude"))
    call check(NF90_PUT_ATT(nc_id, lonvar_id, "long_name", "Longitude"))
    call check(NF90_PUT_ATT(nc_id, latvar_id, "units", "degrees_north"))
    call check(NF90_PUT_ATT(nc_id, latvar_id, "standard_name", "latitude"))
    call check(NF90_PUT_ATT(nc_id, latvar_id, "long_name", "Latitude"))
    
    !write(str_time, *) "days since ", yyyy,"-",mm,"-",dd," ",hh,":",mi,":",ss
    !print *, str_time
    str_time = 'days since '//adate(1:4)//'-'//adate(5:6)// &
     '-'//adate(7:8)//' '//atime(1:2)//':'//atime(3:4)//':'//atime(5:6)
    print *, str_time

    call check(NF90_PUT_ATT(nc_id, timevar_id, "units", trim(str_time)))
    call check(NF90_PUT_ATT(nc_id, timevar_id, "standard_name", "time"))
    call check(NF90_PUT_ATT(nc_id, timevar_id, "long_name", "Time middle averaging period"))
    
    call check(NF90_PUT_ATT(nc_id, pftvar_id, "units", "-"))
    call check(NF90_PUT_ATT(nc_id, pftvar_id, "standard_name", "pft"))
    call check(NF90_PUT_ATT(nc_id, pftvar_id, "long_name", "Plant Functional Type"))    

    call check(NF90_PUT_ATT(nc_id, gppvar_id, "units", "kg C m-2 s-1"))
    call check(NF90_PUT_ATT(nc_id, gppvar_id, "standard_name", "GPP"))
    call check(NF90_PUT_ATT(nc_id, gppvar_id, "long_name", "Gross Primary Productivity"))
    call check(NF90_PUT_ATT(nc_id, neevar_id, "units", "kg C m-2 s-1"))
    call check(NF90_PUT_ATT(nc_id, neevar_id, "standard_name", "NEE"))
    call check(NF90_PUT_ATT(nc_id, neevar_id, "long_name", "Net Ecosystem Exchange"))
    call check(NF90_PUT_ATT(nc_id, nppvar_id, "units", "kg C m-2 s-1"))
    call check(NF90_PUT_ATT(nc_id, nppvar_id, "standard_name", "NPP"))
    call check(NF90_PUT_ATT(nc_id, nppvar_id, "long_name", "Net Primary Productivity"))
    call check(NF90_PUT_ATT(nc_id, trvar_id,  "units", "kg C m-2 s-1"))
    call check(NF90_PUT_ATT(nc_id, trvar_id,  "standard_name", "TotalResp"))
    call check(NF90_PUT_ATT(nc_id, trvar_id,  "long_name", "Total Respiration"))
    call check(NF90_PUT_ATT(nc_id, arvar_id, "units", "kg C m-2 s-1"))
    call check(NF90_PUT_ATT(nc_id, arvar_id, "standard_name", "AutoResp"))
    call check(NF90_PUT_ATT(nc_id, arvar_id, "long_name", "Autotrophic Respiration"))
    call check(NF90_PUT_ATT(nc_id, hrvar_id, "units", "kg C m-2 s-1"))
    call check(NF90_PUT_ATT(nc_id, hrvar_id, "standard_name", "HeteroResp"))
    call check(NF90_PUT_ATT(nc_id, hrvar_id, "long_name", "Heterotrophic Respiration"))
    call check(NF90_PUT_ATT(nc_id, srvar_id, "units", "kg C m-2 s-1"))
    call check(NF90_PUT_ATT(nc_id, srvar_id, "standard_name", "SoilResp"))
    call check(NF90_PUT_ATT(nc_id, srvar_id, "long_name", "Soil Respiration"))

    call check(NF90_PUT_ATT(nc_id, laivar_id, "units", "m2 m-2"))
    call check(NF90_PUT_ATT(nc_id, laivar_id, "standard_name", "LAI"))
    call check(NF90_PUT_ATT(nc_id, laivar_id, "long_name", "Leaf Area Index"))
    call check(NF90_PUT_ATT(nc_id, scvar_id, "units", "kg C m-2"))
    call check(NF90_PUT_ATT(nc_id, scvar_id, "standard_name", "	soil_carbon_content_of_soil_layer"))
    call check(NF90_PUT_ATT(nc_id, scvar_id, "long_name", "Soil Carbon Content by Layer"))
    call check(NF90_PUT_ATT(nc_id, stvar_id, "units", "mol CO2 m-2 s-1"))     ! Unit not consistent with pecan
    call check(NF90_PUT_ATT(nc_id, stvar_id, "standard_name", "stomatal_conductance"))
    call check(NF90_PUT_ATT(nc_id, stvar_id, "long_name", "Stomatal Conductance"))
    call check(NF90_PUT_ATT(nc_id, evvar_id, "units", "kg m-2 s-1"))
    call check(NF90_PUT_ATT(nc_id, evvar_id, "standard_name", "Evaporation"))
    call check(NF90_PUT_ATT(nc_id, evvar_id, "long_name", "Total Evaporation"))
    call check(NF90_PUT_ATT(nc_id, travar_id, "units", "kg m-2 s-1"))
    call check(NF90_PUT_ATT(nc_id, travar_id, "standard_name", "Transpiration"))
    call check(NF90_PUT_ATT(nc_id, travar_id, "long_name", "Total transpiration"))
    call check(NF90_PUT_ATT(nc_id, smvar_id, "units", "m3 m-3"))
    call check(NF90_PUT_ATT(nc_id, smvar_id, "standard_name", "Volumetric soil moisture"))
    call check(NF90_PUT_ATT(nc_id, smvar_id, "long_name", "Average Layer Soil Moisture"))
    call check(NF90_PUT_ATT(nc_id, smpvar_id, "units", "MPa"))
    call check(NF90_PUT_ATT(nc_id, smpvar_id, "standard_name", "soil water potential"))
    call check(NF90_PUT_ATT(nc_id, smpvar_id, "long_name", "Average Layer Soil water potential"))    

    call check(NF90_PUT_ATT(nc_id, cyvar_id, "units", "kg m-2"))
    call check(NF90_PUT_ATT(nc_id, cyvar_id, "standard_name", "CropYield"))
    call check(NF90_PUT_ATT(nc_id, cyvar_id, "long_name", "CropYield"))
    call check(NF90_PUT_ATT(nc_id, abvar_id, "units", "kg C m-2"))
    call check(NF90_PUT_ATT(nc_id, abvar_id, "standard_name", "AGB"))
    call check(NF90_PUT_ATT(nc_id, abvar_id, "long_name", "Total aboveground biomass"))
    call check(NF90_PUT_ATT(nc_id, tbvar_id, "units", "kg C m-2"))
    call check(NF90_PUT_ATT(nc_id, tbvar_id, "standard_name", "TotLivBiom"))
    call check(NF90_PUT_ATT(nc_id, tbvar_id, "long_name", "Total living biomass"))
    call check(NF90_PUT_ATT(nc_id, lcvar_id, "units", "kg C m-2"))
    call check(NF90_PUT_ATT(nc_id, lcvar_id, "standard_name", "leaf_carbon_content"))
    call check(NF90_PUT_ATT(nc_id, lcvar_id, "long_name", "Leaf Carbon Content"))
    call check(NF90_PUT_ATT(nc_id, rcvar_id, "units", "kg C m-2"))
    call check(NF90_PUT_ATT(nc_id, rcvar_id, "standard_name", "root_carbon_content_of_size_class"))
    call check(NF90_PUT_ATT(nc_id, rcvar_id, "long_name", "Root Carbon Content"))
    call check(NF90_PUT_ATT(nc_id, fpvar_id, "units", "-"))
    call check(NF90_PUT_ATT(nc_id, fpvar_id, "standard_name", "fPAR"))
    call check(NF90_PUT_ATT(nc_id, fpvar_id, "long_name", "Absorbed fraction incoming PAR"))
    call check(NF90_PUT_ATT(nc_id, jmvar_id, "units", "umol C m-2 s-1"))
    call check(NF90_PUT_ATT(nc_id, jmvar_id, "standard_name", "Jmax"))
    call check(NF90_PUT_ATT(nc_id, jmvar_id, "long_name", "Maximum electron-transport capacity of leaves, under light saturation"))
    call check(NF90_PUT_ATT(nc_id, vcvar_id, "units", "umol C m-2 s-1"))
    call check(NF90_PUT_ATT(nc_id, vcvar_id, "standard_name", "Vcmax"))
    call check(NF90_PUT_ATT(nc_id, vcvar_id, "long_name", "Maximum carboxylation capacity of leaves"))
    call check(NF90_PUT_ATT(nc_id, dpsivar_id, "units", "MPa"))
    call check(NF90_PUT_ATT(nc_id, dpsivar_id, "standard_name", "Dpsi"))
    call check(NF90_PUT_ATT(nc_id, dpsivar_id, "long_name", "Soil-to-leaf water potential difference"))
    call check(NF90_PUT_ATT(nc_id, chivar_id, "units", "-"))
    call check(NF90_PUT_ATT(nc_id, chivar_id, "standard_name", "Leaf internal-to-external CO2 ratio"))
    call check(NF90_PUT_ATT(nc_id, chivar_id, "long_name", "Leaf internal-to-external CO2 ratio"))
    call check(NF90_PUT_ATT(nc_id, provar_id, "units", "umol m-2 s-1"))
    call check(NF90_PUT_ATT(nc_id, provar_id, "standard_name", "Profit"))
    call check(NF90_PUT_ATT(nc_id, provar_id, "long_name", "Optimized profit of photosynthesis"))


    !Finished defining
    call check( nf90_enddef(nc_id) )

    !Save data to variables
    !print *, "lon=", lon
    call check( nf90_put_var(nc_id, lonvar_id, lon))
    !print *, "lat=", lat
    call check( nf90_put_var(nc_id, latvar_id, lat))
    !print *, "pft=", pft_type
    call check( nf90_put_var(nc_id, pftvar_id, trim(pft_type),(/1,1/)))

    print*, "Prepared NetCDF output file"
    call check( nf90_close(nc_id) )
  end subroutine netCDF_prepareOUTPUT


  subroutine netCDF_writeOUTPUT(filename, var_name, var_data, time, i)

     use netcdf
     implicit none

     character(*) :: filename, var_name
     integer :: i, ncid, VarId_date, VarId_var
     real(8)    :: time
     real(8), dimension(:,:,:) :: var_data

     call check (nf90_open(filename, nf90_Write, ncid))
     !Get id number
     call check (nf90_inq_varid(ncid, "time", VarId_date))
     call check (nf90_inq_varid(ncid, var_name, VarId_var))
     call check( nf90_put_var(ncid, VarId_date, time, (/i+1/)))
     call check( nf90_put_var(ncid, VarId_var,var_data, (/1,1,i+1/)))
     call check( nf90_close(ncid) )

  end subroutine netCDF_writeOUTPUT


  subroutine netCDF_readvar(filename, var_name, var_data, ntim)
    ! read  one time step in the file at once?
    use netcdf
    implicit none

    character(*), intent(in) :: filename, var_name
    integer, intent(in)      :: ntim
    real(8),dimension(:,:,:), intent (out)       :: var_data

    integer      :: i, ncid, varid, ndimsi
    character    :: var_name_tmp
    character    :: dim_name_tmp(4)
    integer      :: xtype
    integer      :: dimids(4)                     ! netCDF dimension ids
    integer      :: bego(4),leno(4)               ! netCDF bounds
    integer      :: begi(4),leni(4)               ! netCDF bounds 
    

    call check (nf90_open(filename, nf90_nowrite, ncid))
    !Get id number
    call check (nf90_inq_varid(ncid, var_name, varid))
    !Get dimension id of variables, Get dimension number of variables 
    !call check (nf90_inq_vardimid(ncid, varid, dimids))
    !call check (nf90_inq_varndims(ncid, varid, ndimsi))
    call check (nf90_inquire_variable(ncid, varid, var_name_tmp, xtype, ndimsi, dimids))

    ! Get the length of each dimension
    ! the order of dimension: lon, lat, lev, time
    ! Assuming reading whole lon, lat, lev and time range?

    if (ndimsi ==4) then
      begi(1) = 1
      begi(2) = 1
      begi(3) = 1
      begi(4) = ntim
      call check (nf90_inquire_dimension(ncid, dimids(1), dim_name_tmp(1), leni(1)))
      call check (nf90_inquire_dimension(ncid, dimids(2), dim_name_tmp(2), leni(2)))
      call check (nf90_inquire_dimension(ncid, dimids(3), dim_name_tmp(3), leni(3)))
      leni(4) = 1
    else if (ndimsi== 3) then
      begi(1) = 1
      begi(2) = 1
      begi(3) = ntim
      call check (nf90_inquire_dimension(ncid, dimids(1), dim_name_tmp(1), leni(1)))
      call check (nf90_inquire_dimension(ncid, dimids(2), dim_name_tmp(2), leni(2)))
      leni(3) = 1
    end if

    !do m = 1, ntim
      !if (ndimsi == 4) begi(4)=m
      !if (ndimsi == 3) begi(3)=m
        
    ! Get values of variables
    call check (nf90_get_var(ncid, varid, var_data, start=begi(1:ndimsi), count=leni(1:ndimsi)))

    ! Close netcdf file
    call check( nf90_close(ncid))

  end subroutine netCDF_readvar


  subroutine netCDF_readTime(filename, ntim, start_input_time, end_input_time)
    ! Read time info included in netcdf file

    use netcdf
    implicit none

    character(*), intent(in) :: filename
    integer, intent(out)     :: ntim
    real(8), dimension(1), intent(out)     :: start_input_time, end_input_time
 
    integer  :: ncid, dimid, varid
    character :: dim_name_tmp

    call check (nf90_open(filename, nf90_nowrite, ncid))
    ! Get time dimension id number
    call check (nf90_inq_dimid(ncid, 'time', dimid))
    ! Get time dimension length 
    call check (nf90_inquire_dimension(ncid, dimid, dim_name_tmp, ntim))

    ! Get values of variables
    call check (nf90_inq_varid(ncid, 'time', varid))
    call check (nf90_get_var(ncid, varid, start_input_time, start=(/1/), count=(/1/)))
    call check (nf90_get_var(ncid, varid, end_input_time, start=(/ntim/), count=(/1/)))

    ! Close netcdf file
    call check( nf90_close(ncid))

  end subroutine

 
  subroutine netCDF_readlonlat(filename, nsites, lat, lon)
    ! Read time info included in netcdf file

    use netcdf
    implicit none

    character(*), intent(in) :: filename
    integer, intent(out)     :: nsites
    real(8), dimension(:), intent(out)     :: lat, lon
 
    integer  :: ncid, dimid, varid1,varid2
    character :: dim_name_tmp

    call check (nf90_open(filename, nf90_nowrite, ncid))
    ! Get time dimension id number
    call check (nf90_inq_dimid(ncid, 'latitude', dimid))
    ! Get time dimension length 
    call check (nf90_inquire_dimension(ncid, dimid, dim_name_tmp, nsites))

    ! Get values of variables
    call check (nf90_inq_varid(ncid, 'latitude', varid1))
    call check (nf90_inq_varid(ncid, 'longitude', varid2))
    call check (nf90_get_var(ncid, varid1, lat, start=(/1/), count=(/nsites/)))
    call check (nf90_get_var(ncid, varid2, lon, start=(/1/), count=(/nsites/)))

    ! Close netcdf file
    call check( nf90_close(ncid))

  end subroutine
  

  subroutine netCDF_readClim(filename,temp, ppfd, Rn, prec, sh, rh, vpd, pres, co2, wind, i)
    ! Read meteorological forcing variables
    ! temp.... other forcing variables need to be defined
  
    use netcdf
    implicit none

    character(*), intent(in) :: filename
    integer, intent(in)      :: i
    real(8),dimension(:,:,:), intent(out) :: temp, ppfd, Rn, prec, sh, rh, vpd, pres, co2, wind

    call netCDF_readvar(filename, "air_temperature", temp, i)
    call netCDF_readvar(filename, "surface_downwelling_photosynthetic_photon_flux_in_air", ppfd,i)
    call netCDF_readvar(filename, "surface_downwelling_shortwave_flux_in_air", Rn ,i)
    call netCDF_readvar(filename, "precipitation_flux", prec, i)
    call netCDF_readvar(filename, "specific_humidity", sh, i)
    call netCDF_readvar(filename, "relative_humidity", rh, i)
    call netCDF_readvar(filename, "water_vapor_saturation_deficit", vpd, i)
    call netCDF_readvar(filename, "air_pressure", pres, i)
    call netCDF_readvar(filename, "mole_fraction_of_carbon_dioxide_in_air", co2, i)
    call netCDF_readvar(filename, "wind_speed", wind, i)
    
  end subroutine netCDF_readClim 

  subroutine netCDF_readlai(filename, lai, i)

     use netcdf
     implicit none

     character(*), intent(in) :: filename
     integer, intent(in)      :: i
     real(8), dimension(:,:,:), intent (out)       :: lai

     call netCDF_readvar(filename, "LAI", lai, i)

  end subroutine netCDF_readlai

  subroutine netCDF_readsoilmoist(filename, soilmoist, i)

     use netcdf
     implicit none


     character(*), intent(in) :: filename
     integer, intent(in)      :: i
     real(8), dimension(:,:,:), intent (out)       :: soilmoist

     call netCDF_readvar(filename, "SoilMoist", soilmoist, i)

  end subroutine netCDF_readsoilmoist

end module io_mod