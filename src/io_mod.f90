MODULE write_output_mod

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

  subroutine netCDF_prepareOUTPUT(output_filename, lons, lats, ntim)
    !Initialize a netCDF-file for emission fields
    use readctrl_mod
    use netcdf
    implicit none

    character(*) :: input_filename
    character(80):: str_time
    integer     :: nc_id, status
    integer :: londim_id, latdim_id, timedim_id,lonvar_id, latvar_id, timevar_id,emitvar_id
    ! integer :: soil_id, area_id,tot_em_id,time_s_dim_id
    integer :: soil_id,area_id,clay_id,sand_id,tot_em_id,time_s_dim_id
    integer :: singdim_id, hourvar_id,dayvar_id
    !Some vars for standard netcdf example
    real, dimension(0:nx_lon_out-1) :: lons
    real, dimension(0:ny_lat_out-1) :: lats
    !real, dimension( int(releaseDays*24/time_step)) :: dates

    call check(nf90_create(trim(input_filename), cmode = NF90_HDF5, ncid = nc_id) )
     
    !Define dimensions
    call check(NF90_DEF_DIM(nc_id, "lon", nx_lon, londim_id))
    call check(NF90_DEF_DIM(nc_id, "lat", ny_lat, latdim_id))
    call check(NF90_DEF_DIM(nc_id, "time", ntim, timedim_id))
    call check(NF90_DEF_DIM(nc_id, "pft", num_pft, pftdim_id))

    !Define variables
    !   status=NF90_DEF_VAR(nc_id,"lon",NF90_FLOAT, (/ londim_id /), lonvar_id)
    call check(nf90_def_var(nc_id, "lon", nf90_float, (/ londim_id /), lonvar_id))
    call check(nf90_def_var(nc_id, "lat", nf90_float, (/ latdim_id /), latvar_id))
    call check(nf90_def_var(nc_id, "time", nf90_int, (/timedim_id/), timevar_id))
    call check(nf90_def_var(nc_id, "pftname", nf90_char, (/pftdim_id/), pftvar_id))
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
    call check(nf90_def_var(nc_id, "Transp", nf90_float, (/londim_id,latdim_id,timedim_id/), trvar_id))
    call check(nf90_def_var(nc_id, "SoilMoist", nf90_float, (/londim_id,latdim_id,timedim_id/), smvar_id))
    call check(nf90_def_var(nc_id, "SoilMoistPot", nf90_float, (/londim_id,latdim_id,timedim_id/), smpvar_id))
    call check(nf90_def_var(nc_id, "CropYield", nf90_float, (/londim_id,latdim_id,timedim_id/), cyvar_id))
    call check(nf90_def_var(nc_id, "AGB", nf90_float, (/londim_id,latdim_id,timedim_id/), abvar_id))
    call check(nf90_def_var(nc_id, "TotLivBiom", nf90_float, (/londim_id,latdim_id,timedim_id/), tbvar_id))
    call check(nf90_def_var(nc_id, "leaf_carbon_content", nf90_float, (/londim_id,latdim_id,timedim_id/), lcvar_id))
    call check(nf90_def_var(nc_id, "root_carbon_content", nf90_float, (/londim_id,latdim_id,timedim_id/), rcvar_id))
    call check(nf90_def_var(nc_id, "fPAR", nf90_float, (/londim_id,latdim_id,timedim_id/), fpvar_id))

    !Attributes
    call check(NF90_PUT_ATT(nc_id, lonvar_id, "units", "degrees_east"))
    call check(NF90_PUT_ATT(nc_id, lonvar_id, "standard_name", "longitude"))
    call check(NF90_PUT_ATT(nc_id, lonvar_id, "long_name", "Longitude"))
    call check(NF90_PUT_ATT(nc_id, latvar_id, "units", "degrees_north"))
    call check(NF90_PUT_ATT(nc_id, latvar_id, "standard_name", "latitude"))
    call check(NF90_PUT_ATT(nc_id, latvar_id, "long_name", "Latitude"))
    
    write(str_time, *) "days since ", start_date_day,"-",start_date_hour," UTC"
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
    call check(NF90_PUT_ATT(nc_id, stvar_id, "units", "kg m-2 s-1"))
    call check(NF90_PUT_ATT(nc_id, stvar_id, "standard_name", "stomatal_conductance"))
    call check(NF90_PUT_ATT(nc_id, stvar_id, "long_name", "Stomatal Conductance"))
    call check(NF90_PUT_ATT(nc_id, evvar_id, "units", "kg m-2 s-1"))
    call check(NF90_PUT_ATT(nc_id, evvar_id, "standard_name", "Evaporation"))
    call check(NF90_PUT_ATT(nc_id, evvar_id, "long_name", "Total Evaporation"))
    call check(NF90_PUT_ATT(nc_id, trvar_id, "units", "kg m-2 s-1"))
    call check(NF90_PUT_ATT(nc_id, trvar_id, "standard_name", "Transpiration"))
    call check(NF90_PUT_ATT(nc_id, trvar_id, "long_name", "Total transpiration"))
    call check(NF90_PUT_ATT(nc_id, smvar_id, "units", "kg m-2"))
    call check(NF90_PUT_ATT(nc_id, smvar_id, "standard_name", "Soil moisture"))
    call check(NF90_PUT_ATT(nc_id, smvar_id, "long_name", "Average Layer Soil Moisture"))
    call check(NF90_PUT_ATT(nc_id, smpvar_id, "units", "mm"))
    call check(NF90_PUT_ATT(nc_id, smpvar_id, "standard_name", "soil suction, negative [mm]"))
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

    !Finished defining
    call check( nf90_enddef(nc_id) )

    !Save data to variables
    call check( nf90_put_var(nc_id, lonvar_id, lons) )
    call check( nf90_put_var(nc_id, latvar_id, lats) )
    call check( nf90_put_var(nc_id, timevar_id, time) )
    call check( nf90_put_var(nc_id, pftvar_id, pftname) )

    print*, "Prepared NetCDF output file"
    call check( nf90_close(nc_id) )
  end subroutine netCDF_prepareOUTPUT

  subroutine netCDF_writeOUTPUT(output_filename, var_name, var_data, time, i)

     use netcdf
     implicit none

     character(*) :: output_filename, var_name
     integer :: i, ncid, VarId_date, VarId_var
     integer, dimension(1) :: time
     real, dimension(0:nx_lon_out-1, 0:ny_lat_out-1) :: var_data

     call check (nf90_open(output_filename, nf90_Write, ncid))
     !Get id number
     call check (nf90_inq_varid(ncid, "time", VarId_date))
     call check (nf90_inq_varid(ncid, var_name, VarId_var))
     call check( nf90_put_var(ncid, VarId_date, time, (/i+1/)))
     call check( nf90_put_var(ncid, VarId_var,var_data, (/1,1,i+1/)))
     call check( nf90_close(ncid) )

  end subroutine netCDF_writeOUTPUT

  subroutine netCDF_readvar(input_filename, var_name, var_dat, ntim)
    ! Need to read all the time steps available in the file at once?

    use netcdf
    implicit none

    character(*) :: input_filename, var_name
    integer :: i, ncid, VarId_date, VarId_var
    integer, dimension(0) :: time
    real, dimension(0:nx_lon_out-1, 0:ny_lat_out-1), intent (out) :: var_data

    call check (nf90_open(input_filename, nf90_nowrite, ncid))
    !Get id number
    call check (nf90_inq_varid(ncid, var_name, varid))
    !Get dimension id of variables
    call check (nf90_inq_vardimid(ncidi, varid, dimids))
    !Get dimension number of variables 
    call check (nf90_inq_varndims(ncidi, varid, ndimsi))
    ! Get the length of each dimension
    ! the order of dimension: lon, lat, lev, time
    ! Assuming reading whole lon, lat, lev and time range?
    if (ndimsi ==4) then
      begi(1) = 1
      begi(2) = 1
      begi(3) = 1
      begi(4) = ntim
      call check (nf90_inq_dimlen(ncidi, dimids(1), leni(1)))
      call check (nf90_inq_dimlen(ncidi, dimids(2), leni(2)))
      call check (nf90_inq_dimlen(ncidi, dimids(3), leni(3)))
      leni(4) = 1
    else if (ndimsi== 3) then
      begi(1) = 1
      begi(2) = 1
      begi(3) = ntim
      call check (nf90_inq_dimlen(ncidi, dimids(1), leni(1)))
      call check (nf90_inq_dimlen(ncidi, dimids(2), leni(2)))
      leni(3) = 1
    end if

    !do m = 1, ntim
      !if (ndimsi == 4) begi(4)=m
      !if (ndimsi == 3) begi(3)=m
        
    ! Get values of variables
    call check(nf90_get_var (ncidi, varid, var_data, begi(1:ndimsi), leni(1:ndimsi)))

    ! Close netcdf file
    call check( nf90_close(ncid))

  end subroutine netCDF_readvar

  subroutine netCDF_readTime(input_filename, ntim, start_input_time, end_input_time)
    ! Read time info included in netcdf file

    call check (nf90_open(input_filename, nf90_nowrite, ncid))
    ! Get time dimension id number
    call check (nf90_inq_dimid(ncidi, 'time', dimid))
    ! Get time dimension length 
    call check (nf909_inq_dimlen(ncidi, dimid, ntim))

    ! Get values of variables
    call check (nf90_inq_varid(ncid, 'time', varid))
    call check (nf90_get_var (ncidi, varid, start_input_time, 1, 1))
    call check (nf90_get_var (ncidi, varid, end_input_time, ntim, 1))

    ! Close netcdf file
    call check( nf90_close(ncid))

  end subroutine

  subroutine netCDF_readClim(input_filename,temp, ppfd, prec, sh, rh, vpd, pres, co2, i)
    ! Read meteorological forcing variables
    ! temp.... other forcing variables need to be defined
    ! start_input_time, end_input_time need to be set

    use netcdf
    implicit none

    character(*) :: input_filename, var_name
    integer :: i, ncid, VarId_date, VarId_var
    integer, dimension(0) :: time
    real, dimension(0:nx_lon_out-1, 0:ny_lat_out-1), intent (out) :: var_data
    real, intent (out) :: start_input_time, end_input_time

    call netCDF_readvar(input_filename, "air_temperature", temp, i)
    call netCDF_readvar(input_filename, "surface_downwelling_photosynthetic_photon_flux_in_air", ppfd,i)
    call netCDF_readvar(input_filename, "precipitation_flux", prec, i)
    call netCDF_readvar(input_filename, "specific_humidity", sh, i)
    call netCDF_readvar(input_filename, "relative_humidity", rh, i)
    call netCDF_readvar(input_filename, "water_vapor_saturation_deficit", vpd, i)
    call netCDF_readvar(input_filename, "air_pressure", pres, i)
    call netCDF_readvar(input_filename, "mole_fraction_of_carbon_dioxide_in_air", co2, i)
    
  end subroutine netCDF_readClim 

  subroutine netCDF_readlai(input_filename, lai, i)

     use netcdf
     implicit none

     character(*) :: input_filename, var_name
     integer :: i, ncid, VarId_date, VarId_var
     integer, dimension(0) :: time
     real, dimension(0:nx_lon_out-1, 0:ny_lat_out-1), intent (out) :: var_data

     call netCDF_readvar(input_filename, "LAI", lai, i)

  end subroutine netCDF_readlai

end module write_output_mod