!**********************************************************************
! Copyright 2023                                                      *
! FMI                                                                 *
!                                                                     *                                    *
! ***.****@fmi.fi                                                     *                                                              *
!                                                                     *
! This file is part of FMI-SVM.                                       *
!                                                                     *
! FMI-SVM is free software: you can redistribute it and/or modify     *
! it under the terms of the GNU General Public License as published by*
! the Free Software Foundation, either version 3 of the License, or   *
! (at your option) any later version.                                 *
!                                                                     *
! You should have received a copy of the GNU General Public License   *
! along with FMI-SVM.  If not, see <http://www.gnu.org/licenses/>.    *
!**********************************************************************

program SVMC

  use netcdf                    ! library for processing netcdf files
  use initialization_mod        ! initialize svm model
  use readctrl_mod              ! module for reading control parameters from namelist file
  use readvegpara_mod           ! module for reading vegetation parameter from namelist file
  use readsoilpara_mod          ! module for reading soil parameter from namelist file 
  use readclim_mod              ! module for reading reading meteorological forcing data
  use io_mod                    ! manage input/output of the model

  use phydro_mod                ! module for p-hydro
  use spafhy_mod                ! module for soil water bucket model, which will provide psi_soil for p-hydro
  !use alloc_mod                ! module for carbon allocation and yield
  use yasso                     ! module for soil decomposition model, which will provide heterogeneous respiration (hr)   

  implicit none

  ! Loop variables
  !***********************************
  integer         :: i, j, k, t
  integer         :: step_nc_hr, step_nc_day, step_clim, step_lai
  real(kind=dp)   :: tot_hour, tot_hour_end, juldate, start_date, end_date
  real            :: start_clim_time, end_clim_time 
  real            :: start_lai_time, end_lai_time
  
  integer         :: ntim_clim, ntim_lai, ntim_out_hr, ntime_out_day
 
  !***********************************
  !Model variables
  !***********************************
  integer   ::

  ! Input variables for p-hydro (the definition is from rpmodel.R)
  real     ::    tc        ! Air temperature (tc), degrees C
  real     ::    ppfd      ! Photosynthetic photon flux density (mol m-2 d-1) (incoming solar radiation from forcing data?)
  real     ::    vpd       ! Vapour pressure deficit (Pa) (will be calculated using pressure & humidity)
  real     ::    co2       ! Atmospheric CO2 concentration (ppm)
  real     ::    elv       ! Elevation above sea-level (m.a.s.l.) (not needed if we have surface pressure!)
  real     ::    fapar     ! Fraction of absorbed photosynthetically active radiation (unitless) (will be calculated using LAI)
  real     ::    kphio     ! Apparent quantum yield efficiency (unitless).
  real     ::    psi_soil  ! soil water potential (Mpa)
  real     ::    rdark = 0 !
  character(len=200)  :: opt_hypothesis=''   ! character, Either "Lc" or "PM"

  type(par_plant_type)          :: par_plant           ! A list of plant hydraulic parameters (will be defined in readpara_mod.f90).
  type(par_cost_type)           :: par_cost            ! A list of cost parameters (will be defined in readpara_mod.f90).

  character(len=200)  ::
  real,dimension(:,:), allocatable
  logical,dimension(:,:), allocatable ::
  !*********************************** 
 
  ! Initialize the model, including
  ! (1) reading spatial and temporal configurations, e.g., longitude, latitude, beginning and ending years.
  ! (2) read and initialize model parameters
  ! (3) read model inputdata including: (1) initial status file if availabe? (2) meteorological forcing, (3) soil properties
  ! (4) Different modes: cold start (everything bareground); 
  !                      initial file (soil carbon and moisture status provided either from obs or previous experiment)
  !                      restart file (reproduce the status from previous experiment)
  
  !***********************************
  call initialization
   
  psi_soil=0


  !***********************************
  ! Set time control parameters
  !***********************************
  
  start_date=juldate(start_date_day, start_date_hour)
  end_date  =juldate(end_date_day, end_date_hour)
 
  ! Calculate total hours of the simulation
  tot_hour_end= (juldate(start_date_day, start_date_hour) - juldate(end_date_day, end_date_hour))*24  
  ntim_out_hr = mod(tot_hour_end, time_step_output)
  ntim_out_day= mod(tot_hour_end, 24)

  ! Run the model
  tot_hour=0.0
  step_nc_hr=0
  step_nc_day=0

  ! Read time series of input data
  call netCDF_readTime(input_climfile, ntim_clim, start_clim_time, end_clim_time)
  call netCDF_readTime(input_laifile, ntim_lai, start_lai_time, end_lai_time)
  
  ! To simplify the time management, specify the julian start date of the inputdata by hand.
  start_clim_juldate=juldate(20201231,233000)
  start_lai_juldate =juldate(20210101,000000)
  
  ! step the starting time steps for reading input files
  step_clim = floor((start_date-start_clim_juldate)*24)+1
  step_lai  = floor(start_date-start_lai_juldate)+1
   
  ! Loop over time, and locations
  do while (tot_hour .le. tot_hour_end)  ! in hour or 30 minutes, time loop 
    
    do i=1,num_sites  ! site loop (we do not use lon-lat box to allow the flexibility to run sites or regional/global simulations)
      
      do m=1,num_pft   ! crop/pft loop
                       ! Do we really need it? We can run individual experiments to represent different crops 
                       ! This may help for scaling up?

        !if (is_pheno_on()) then ! Only when crop/vegetation is present, can be turned off with prescribed LAI
        ! need a bit input preparation for runnnig phydro
        
          ! call pmodel_input_prep()

          ! Determine whether to read new lai data
          ! Only update LAI daily
          if (mod(tot_hour,24) .eq. 0) then
            call netCDF_readlai(input_laifile, lai, step_lai)
            step_lai=step_lai+1
            fapar= 1-exp(-k*LAI)
          end if

          ! Determine whether to read new climate variables
          ! if () then
          call netCDF_readClim(input_climfile, temp, ppfd, prec, sh, rh, vpd, pres, co2, step_clim)
          step_clim=step_clim+1
          ! end if

          ! run phydro to estimate photosynthetic rate (a) and stomatal conductance (gs)
          ! At what time scale the optimization should work need to be tested!!!!
          
          call pmodel_hydraulics_numerical(temp-273.15, ppfd*3600*24, vpd, co2*1000000, pres, fapar, &
                                 kphio, psi_soil, 0,                                                 &
                                 conductivity, psi50, b, alpha, gamma,                               &
                                 opt_hypothesis = "PM",                                              &
                                 jmax, dpsi, gs, aj, ci, chi, vcmax, profit, chi_jmax_lim            &
                                 )

          ! Carbon allocation: update gpp, npp, ar ....
          ! call carbon_allocation_hr(a,....)          
        
          ! Update transpiration, canopy evaporation...
          ! CanopyGrid in spafhy...
          call canopy_water_flux(gs, ppfd*fapar*, temp, prec, Rg, Par, VPD, U=2.0, co2, Rew=1.0, beta=1.0, P=101300.0)       
        

        !end if  !phenology
        
        ! Calculation of NEE & LATENT heat flux
        nee= gpp-ar-hr
        et = tr+evap_can(temp, rad)+evap_soil(temp, rad)

        ! Write hourly output at output time step frequency

        if ( mod((tot_hour,time_step_output).eq.0. ) then
            
          !Write hourly output for this time step
          !************************************************************************
          if(step_nc_hr.eq.0)then
            !Initialize netcdf file
            call netCDF_prepareOUTPUT(output_filename_hr, lon_sites, lat_sites, ntim_out_hr)
          end if

          call netCDF_writeOUTPUT(output_filename_hr, "GPP", gpp, tot_hour/24, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "Evap", tot_evap, tot_hour/24, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "Transp", tr, tot_hour/24, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "SoilMoist", soilwater_state%WatSto, tot_hour/24, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "SoilMoistPot", smp, tot_hour/24, step_nc_hr)
          step_nc_hr= step_nc_hr+1 
        endif
        
        ! Comulative parts can be a separate module in the future
        ! Cumulative GPP or NPP which will be used for carbon allocation on daily or yearly scale. 
        gpp_sum_day=gpp_sum_day + gpp_hr*3600*...           
        npp_sum_day=npp_sum_day + npp_hr*3600*...
        gpp_sum_year=gpp_sum_year + gpp_hr*3600*...                   
        npp_sum_year=npp_sum_year + npp_hr*3600*...
          
        ! Comulative temperature which will be used for GDD calculation and phenology
        temp_day=temp_day+temp_hr

        ! Cumulative transpiration and canopy evaporation
        tr_day=tr_day + tr
        evap_can_day=evap_can_day+evap_can

        ! Cumulative solar radiation (energy)

        
        if ((mod(tot_hour,24) .eq. 0).and.(tot_hour .gt. 0)) then    ! here assume the start time is always the beginning of the day!
          !Update GDD which is the criteria for phenology stages         
          gdd_sum= gdd_sum+temp_day/24.....
          temp_day=0
          
          ! run Topmodel
          ! catchment average ground water recharge [m per unit area]
          call topmodel(gs or tr_veg, qd(t-1), qr(t)...) 

          ! run CanopyGrid (spafhy) (move this function in hourly cycle)
          ! call canopy_water_flux(gs, Rn, Ta, Prec, Rg, Par, VPD, U=2.0, CO2=380.0, Rew=1.0, beta=1.0, P=101300.0)

          ! update soil water with the bucket model
          ! run BucketGrid water balance: watbal
          call soil_water(gs or tr_veg, qr(t), qd(t) ...)  

          ! calculate soil water potential based on volumetric soil moisture (P-V curve)
          call soil_water_retention_curve(soil_water_v) 

          ! Update phenological stages according to growing degree day     
          call phenology(gdd_sum.....)                                    
          
          if (is_pheno_on()) then ! Only when crop/vegetation is present
            ! Update daily variables: LAI, Aboveground & Belowground biomass according to phenological stages
            ! also litter input if running yasso on hourly or daily step
            call carbon_allocation_day(pheno_stage, npp_sum_day.....)     
            
            ! Yasso: split input c into various yasso fractions
            call inputs_to_fractions(leaf, root, soluble, compost, fract)
          
          end if

          ! Yasso: create average meteorological forcings for yasso
          call average_met(met_daily, met_rolling, aver_size, met_state, met_ind)
          
          ! Yasso: update soil respiration, should we update soil respiration hourly or daily? 
          ! resp should be an output variable
          call decompose(param, timestep_days, c_input_awenh_day, nitr_input_day, tempr_c, &
                                      precip_day, cstate, nstate, ctend, ntend) 

          ! We could also put yasso here can do the daily calculation for NEE
          ! call yasso20_day()
          ! nee_day= .......

          ! Write output for daily variables
          if(step_nc_day.eq.0)then
            !Initialize netcdf file
            call netCDF_prepareOUTPUT(output_filename_day, lon_sites, lat_sites, ntim_out_day)
          end if

          call netCDF_writeOUTPUT(output_filename_day, "HeteroResp", hr, tot_hour/24, step_nc_day)           
          step_nc_day= step_nc_day+1                

          gpp_sum_day =0
          npp_sum_day =0
          tr_day      =0
          evap_can_day=0
          
        end if

        if (is_end_curr_year() .or. is_pheno_harvest() ) then

          call carbon_allocation_yr(pheno_stage, gpp_sum_year, .....)   ! calculate yield, harvest biomass, 
                                                                   ! and litter input if running yasso  
          
          ! We could also put yasso here to do yearly calculation for NEE
          ! call yasso20_yr()

          call write_output_day()                   ! Write output for yearly variables

          gpp_sum_year=0
          npp_sum_year=0
          gdd_sum=0  

        end if

      end do ! m
    end do  !i

    ! advance time step
    tot_hour=tot_hour + time_step
  end do ! t

  call write_restart()

end program SVMC
