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
!  use initialization_mod        ! initialize svm model
  use readctrl_mod              ! module for reading control parameters from namelist file
  use readvegpara_mod           ! module for reading vegetation parameter from namelist file
  use readsoilpara_mod          ! module for reading soil parameter from namelist file 
  use io_mod                    ! manage input/output of the model

  use phydro_mod                ! module for p-hydro
  use spafhy_mod                ! module for soil water bucket model, which will provide psi_soil for p-hydro
  !use alloc_mod                ! module for carbon allocation and yield
 ! use yasso                     ! module for soil decomposition model, which will provide heterogeneous respiration (hr)   

  implicit none

  ! Loop variables
  !***********************************
  integer         :: i, m
  integer         :: step_nc_hr, step_nc_day, step_clim, step_lai
  real(kind=dp)   :: tot_hour, tot_hour_end, juldate, start_date, end_date
  real(8), dimension(1)     :: start_clim_time, end_clim_time
  real(8)            :: start_clim_juldate,start_lai_juldate    
  real(8), dimension(1)    :: start_lai_time, end_lai_time 
  integer         :: ntim_clim, ntim_lai, ntim_out_hr, ntim_out_day
 
  !***********************************
  !Model variables
  !***********************************

  ! p-hydro variabless
  ! Input variables for p-hydro (the definition is from rpmodel.R)
  real(8)     ::    temp        ! Air temperature (tc), degrees C
  real(8)     ::    ppfd      ! Photosynthetic photon flux density (mol m-2 d-1) (incoming solar radiation from forcing data?)
  real(8)     ::    vpd       ! Vapour pressure deficit (Pa) (will be calculated using pressure & humidity)
  real(8)     ::    co2       ! Atmospheric CO2 concentration (ppm)
  real(8)     ::    elv       ! Elevation above sea-level (m.a.s.l.) (not needed if we have surface pressure!)
  real(8)     ::    fapar     ! Fraction of absorbed photosynthetically active radiation (unitless) (will be calculated using LAI)
  real(8)     ::    prec      ! 
  real(8)     ::    pres
  real(8)     ::    sh
  real(8)     ::    rh
  real(8)     ::    wind
  real(8)     ::    rg, rn
  real(8)     ::    psi_soil, psi_soil_spafhy  ! soil water potential (Mpa)
  real(8)     ::    soilmoist
  real(8)     ::    rdark     

  real(8), dimension(1,1,1)  :: lai_matrix, soilmoist_matrix
  real(8), dimension(1,1,1)  :: temp_matrix, ppfd_matrix, rg_matrix, prec_matrix, &
                                sh_matrix, rh_matrix, vpd_matrix, wind_matrix, &
                                pres_matrix, co2_matrix, gpp_matrix, &
                                jmax_matrix, vcmax_matrix, dpsi_matrix, &
                                chi_matrix, profit_matrix, gs_matrix, &
                                evap_matrix, psi_soil_matrix, soilmoist1_matrix

  real(8)     ::    lai

  real(8)    :: jmax       !  The maximum rate of RuBP regeneration (umol/m2/s) at growth temperature (argument\code{tc}), calculated using
                                                ! \deqn{A_J = A_C} 
                                                !  Electron transport capacity (umol/m2/s)
  real(8)    :: dpsi       ! soil-to-leaf water potential difference (\eqn{\psi_s-\psi_l}), Mpa
  real(8)    :: gs         ! Stomatal conductance (gs, in mol m-2 s-1)
  real(8)    :: aj         ! electron-transport limited assimilation rate (umol/m2/s)
  real(8)    :: ci         !  leaf-internal CO2 concentration, converted to partial pressure (Pa)
  real(8)    :: chi        ! Optimal ratio of leaf internal to ambient CO2 (unitless).
  real(8)    :: vcmax      !   Carboxylation capacity (umol/m2/s)
  real(8)    :: profit                  ! Net assimilation rate after accounting for costs
  real(8)    :: chi_jmax_lim      ! Analytical chi in the case of strong Jmax limitation
  real(8)    :: gpp
  real(8)    :: tr_phydro   ! Transpiration estimated by p-hydro model
  real(8)    :: LE          ! latent heat flux (Wm-2)

  ! spafhy variables
  real(8)    :: tr_spafhy   ! Transpiration estimated by spafhy model
  real(8)    :: retflow     ! return flow from ground water [m]
  real(8)    :: LE          ! latent heat flux [Wm-2]
  real(8)    :: LatentHeat  ! latent heat of vaporization [J kg-1]

  ! For soil water retention curve
 ! real(8) :: watsat      ! v/v saturate moisture
  !real(8) :: watres      ! v/v, residual soil moisture for Van Genuchten
  !real(8) :: vol_ice     ! v/v, volumetric ice in soil bucket 
  !real(8) :: vol_liq     ! v/v, volumetric of liq in soil bucket     
  !real(8) :: satfrac     ! parameter for Van Genuchten
  !real(8) :: n1, m1, alpha_van           ! (-), pore-size-distribution parameter for Van Genuchten 1.07
  !real(8) :: eff_porosity! v/v, volume of ice

  type(soilwater_type)          :: soilwater_state  
  type(canopywater_type)        :: canopywater_state
  type(snowwater_type)        :: snowwater_state

  !character(len=200)  ::
  !real,dimension(:,:), allocatable
  !logical,dimension(:,:), allocatable ::
  !*********************************** 
 
  ! Initialize the model, including
  ! (1) reading spatial and temporal configurations, e.g., longitude, latitude, beginning and ending years.
  ! (2) read and initialize model parameters
  ! (3) read model inputdata including: (1) initial status file if availabe? (2) meteorological forcing, (3) soil properties
  ! (4) Different modes: cold start (everything bareground); 
  !                      initial file (soil carbon and moisture status provided either from obs or previous experiment)
  !                      restart file (reproduce the status from previous experiment)
  
  !***********************************
  ! call initialization
   
  ! psi_soil=0
  call readctrl_namelist
  call readvegpara_namelist
  call readsoilhydro_namelist

  !call set_soilwaterState(soilwater_state, canopywater_state)
  call initialization_spafhy(canopywater_state, snowwater_state, soilwater_state)


  !***********************************
  ! Set time control parameters
  !***********************************
  
  print *, start_date_day

  start_date=juldate(start_date_day, start_date_hour)
  end_date  =juldate(end_date_day, end_date_hour)
 
  ! Calculate total hours of the simulation
  tot_hour_end= (juldate(end_date_day, end_date_hour) - juldate(start_date_day, start_date_hour))*24  
  ntim_out_hr = tot_hour_end/time_step_output
  ntim_out_day= tot_hour_end/24.0

  ! Run the model
  tot_hour=0.0
  step_nc_hr=0
  step_nc_day=0

  print *, tot_hour_end, ntim_out_hr, ntim_out_day
  ! Read time series of input data
  call netCDF_readTime(input_climfile, ntim_clim, start_clim_time, end_clim_time)
  call netCDF_readTime(input_laifile, ntim_lai, start_lai_time, end_lai_time)
  call netCDF_readlonlat(input_climfile, num_sites, lat_sites, lon_sites)

  ! To simplify the time management, specify the julian start date of the inputdata by hand.
  start_clim_juldate=juldate(20201231,233000)
  start_lai_juldate =juldate(20210101,000000)
  
  ! step the starting time steps for reading input files
  step_clim = floor((start_date-start_clim_juldate)*24)+1
  step_lai  = floor(start_date-start_lai_juldate)+1

  print *, num_sites, lat_sites, lon_sites, tot_hour, num_pft   
  ! Loop over time, and locations
  do while (tot_hour .lt. tot_hour_end)  ! in hour or 30 minutes, time loop 
    
    do i=1,num_sites  ! site loop (we do not use lon-lat box to allow the flexibility to run sites or regional/global simulations)
      
      do m=1,num_pft   ! crop/pft loop
                       ! Do we really need it? We can run individual experiments to represent different crops 
                       ! This may help for scaling up?

        !if (is_pheno_on()) then ! Only when crop/vegetation is present, can be turned off with prescribed LAI
        ! need a bit input preparation for runnnig phydro
        
          ! call pmodel_input_prep()

          ! Determine whether to read new lai data
          ! Only update LAI daily

          if (mod(tot_hour,24.0) .eq. 0) then
            call netCDF_readlai(input_laifile, lai_matrix, step_lai)
            call netCDF_readsoilmoist('../data/FieldObs_Qvidja.2021.soilmoist.nc', soilmoist_matrix, step_lai)
            
            lai=lai_matrix(1,1,1)
            soilmoist=soilmoist_matrix(1,1,1)             
            
            call soil_water_retention_curve(soilmoist, psi_soil)
            ! add soil rentention curve here to test soil water potential calculation
            !n1=1.07       !Launiainen et al. 2022: C1-5: 1.12, 1.14, 1.07, 1.27, 1.18  
            !m1=1.0/n1  
            !watres=0.0   !Launiainen et al. 2022: C1-5: 0.0
            !alpha_van=2.02   !Launiainen et al. 2022: C1-5: 4.45, 5.92, 2.02, 4.49, 3.35
            !watsat=0.46  !Launiainen et al. 2022: C1-5: 0.75, 0.68, 0.46, 0.47, 0.54  

            !vol_liq=soilmoist
            !vol_ice = 0.0   
            !eff_porosity = max(0.01, watsat-vol_ice)
    
            !satfrac  = (vol_liq-watres)/(eff_porosity-watres)
            !print *, "satfrac=", satfrac 
            !psi_soil = -(1.0/alpha_van)*(satfrac**(1.0/(m1-1.0)) - 1.0 )**m1  ! psi_soil in kPa

            ! calculate fapar:
            step_lai=step_lai+1
            fapar= 1-exp(-k*lai)
          end if

          ! Determine whether to read new climate variables
          ! if () then
          call netCDF_readClim(input_climfile, temp_matrix, ppfd_matrix, rg_matrix, prec_matrix, &
                       sh_matrix, rh_matrix, vpd_matrix, pres_matrix, &
                       co2_matrix, wind_matrix, step_clim)

          temp=temp_matrix(1,1,1) ! deg C or K???
          ppfd=ppfd_matrix(1,1,1)
          rg  = rg_matrix(1,1,1)
          prec=prec_matrix(1,1,1)
          sh=sh_matrix(1,1,1)
          rh=rh_matrix(1,1,1)
          vpd=vpd_matrix(1,1,1)
          pres=pres_matrix(1,1,1)
          co2=co2_matrix(1,1,1) 
          wind=wind_matrix(1,1,1)
          
          step_clim=step_clim+1
          ! end if

          ! run phydro to estimate photosynthetic rate (a) and stomatal conductance (gs)
          ! At what time scale the optimization should work need to be tested!!!!
          
          rdark=0.0

          print *, "temp =", temp-273.15                  ! unit should be C
          print *, "ppfd =", ppfd*1000000.0/lai          ! umol/m2/s, current unit is wrong
                                                         ! Averaging light absorption to each unit area of leaf (multi-layer leaf), 
                                                         !     multiply lai when calculating GPP
                                                         ! Alternatively, use total light absorption (one big leaf)
                                                         !     no need to multiply lai when calculating GPP  
          print *, "vpd =", vpd                          ! pa
          print *, "co2 =", co2*1000000                  ! ppm
          print *, "pres =", pres                        ! pa
          print *, "fapar =", fapar                      ! frac
          !print *, "vol_liq =", vol_liq
          print *, "psi_soil =", psi_soil*0.001          ! convert from Kpa to MPa
          
          ! for coupling with SpaFHy: psi_soil = soilwater_state%Psi
          call pmodel_hydraulics_numerical(temp-273.15, ppfd*1000000.0/lai, vpd, co2*1000000, pres, fapar, &
                                 psi_soil*0.001, rdark,                                                 &
                                 jmax, dpsi, gs, aj, ci, chi, vcmax, profit, chi_jmax_lim            &
                                 )
          
          ! HT: need to multiply lai or not? probably not as fapar has considered the effect of lai.
          gpp= aj * c_molmass * 1e-6 * 1e-3 * lai       ! aj in umol/m2/s, gpp kg C/m2/s multi-layer hypothesis
          !gpp= aj * c_molmass * 1e-6 * 1e-3            ! big leaf hypothesis
          print *, "gpp=", gpp, aj, c_molmass, lai
          ! Carbon allocation: update gpp, npp, ar ....
          ! call carbon_allocation_hr(a,....)          
        
          ! Solve plant canopy and soil water budget
          
          ! Transpiration derived from P-hydro
          ! HUI - check units of gs and conversion to tr_phydro. We want it to be [mm s-1 = kg H2O m-2 s-1] 
          ! [tr_phydro] = [1] * [??] * [Pa Pa-1] * [kg mol-1] / [kg m-3] 
          tr_phydro = 1.6*gs*(vpd/pres)*h2o_molmass/density_h2o(temp-273.15, pres)
          
          ! net radiation of the whole canopy-soil system [W m-2]
          ! Samuli will revise later!
          rn= rg * 0.7
          !rn = max(2.57*lai/(2.57*lai+0.57)-0.2, 0.55)*rg  ! Launiainen et al. 2016 GCB, fit to Fig 2a

          !call canopy_water_flux(gs, rn, temp-273.15, prec, ppfd, vpd,  &
          !                        wind, co2*1000000, soilwater_state%Rew, & 
          !                        pres, lai, canopywater_state, snowwater_state)       
          
          ! call SpaFHy code to compute new canopywater_state and snowwater_state
          ! returns water fluxes integrated over time_step in units [mm = kg H2O m-2]
          call canopy_water_flux(rn, temp-273.15, prec, vpd, wind, press, fapar, lai, &
                                  canopywater_state, snowwater_state, soilwater_state)

          ! ET [mm]
          canopywater_state%ET =  tr_phydro * (time_step*3600.0) +  canopywater_state%GroundEvap + &
                                      canopywater_state%CanopyEvap
          
          LatentHeat = 1.0e3 * (3147.5 - 2.37 * (temp))
          LE = canopywater_state%ET / (time_step * 3600.0) * LatentHeat ! Wm-2

          ! Solve soil water balance

          !tr_spafhy=tr_phydro*(time_step*3600.0*1.0e-3)
          retflow=0.0
          ! water fluxes must be in units [m]. Updates soilwater_state, including soilwater_state%Psi.
          call soil_water(soilwater_state, snowwater_state%PotInf*1.0e-3, &
                          tr_phydro*(time_step*3600.0*1.0e-3),  &
                          canopywater_state%GroundEvapfloor*1.0e-3, retflow)  

          !call soil_water_retention_curve(soilwater_state%Wliq, psi_soil_spafhy) 

        !end if  !phenology
        
        ! Calculation of NEE & LATENT heat flux
        ! nee= gpp-ar-hr

        ! Write hourly output at output time step frequency

        if ( mod(tot_hour,time_step_output) .eq. 0.0 ) then
            
          !Write hourly output for this time step
          !************************************************************************
          if(step_nc_hr.eq.0)then
            !Initialize netcdf file
            print *, "ntim_out_hr=", ntim_out_hr
            call netCDF_prepareOUTPUT(output_filename_hr, lon_sites, lat_sites, ntim_out_hr)
          end if

          gpp_matrix(1,1,1)   =gpp
          gs_matrix(1,1,1)    =gs
          jmax_matrix(1,1,1)  =jmax
          vcmax_matrix(1,1,1) =vcmax
          chi_matrix(1,1,1)   =chi
          dpsi_matrix(1,1,1)  =dpsi
          profit_matrix(1,1,1)=profit
          soilmoist1_matrix(1,1,1)=soilwater_state%Wliq
          !psi_soil_matrix(1,1,1)=psi_soil_spafhy
          psi_soil_matrix(1,1,1)=soilwater_state%Psi ! MPa
          evap_matrix(1,1,1)=canopywater_state%ET
          print *, "step_nc_hr=", step_nc_hr

          call netCDF_writeOUTPUT(output_filename_hr, "GPP", gpp_matrix, tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "stomatal_conductance", gs_matrix, tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "Jmax", jmax_matrix, tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "Vcmax", vcmax_matrix, tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "Chi", chi_matrix, tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "Dpsi", dpsi_matrix, tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "Profit", profit_matrix, tot_hour/24.0, step_nc_hr)

          call netCDF_writeOUTPUT(output_filename_hr, "Evap", evap_matrix, tot_hour/24.0, step_nc_hr)
          !call netCDF_writeOUTPUT(output_filename_hr, "Transp", tr, tot_hour/24, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "SoilMoist", soilmoist1_matrix, tot_hour/24.0, step_nc_hr)
          print *, "psi=", psi_soil_matrix, soilmoist1_matrix
          call netCDF_writeOUTPUT(output_filename_hr, "SoilMoistPot", psi_soil_matrix, tot_hour/24.0, step_nc_hr)
          
          step_nc_hr= step_nc_hr+1 
        endif
        
        ! Comulative parts can be a separate module in the future
        ! Cumulative GPP or NPP which will be used for carbon allocation on daily or yearly scale. 
        !gpp_sum_day=gpp_sum_day + gpp_hr*3600*...           
        !npp_sum_day=npp_sum_day + npp_hr*3600*...
        !gpp_sum_year=gpp_sum_year + gpp_hr*3600*...                   
        !npp_sum_year=npp_sum_year + npp_hr*3600*...
          
        ! Comulative temperature which will be used for GDD calculation and phenology
        !temp_day=temp_day+temp_hr

        ! Cumulative transpiration and canopy evaporation
        !tr_day=tr_day + tr
        !evap_can_day=evap_can_day+evap_can

        ! Cumulative solar radiation (energy)

        
        if ((mod(tot_hour,24.0) .eq. 0).and.(tot_hour .gt. 0)) then    ! here assume the start time is always the beginning of the day!
          !Update GDD which is the criteria for phenology stages         
        !  gdd_sum= gdd_sum+temp_day/24.....
        !  temp_day=0
          
          ! run Topmodel
          ! catchment average ground water recharge [m per unit area]
        !  call topmodel(gs or tr_veg, qd(t-1), qr(t)...) 

          ! run CanopyGrid (spafhy) (move this function in hourly cycle)
          ! call canopy_water_flux(gs, Rn, Ta, Prec, Rg, Par, VPD, U=2.0, CO2=380.0, Rew=1.0, beta=1.0, P=101300.0)

          ! update soil water with the bucket model
          ! run BucketGrid water balance: watbal
        !  call soil_water(gs or tr_veg, qr(t), qd(t) ...)  

          ! calculate soil water potential based on volumetric soil moisture (P-V curve)
        !  call soil_water_retention_curve(soil_water_v) 

          ! Update phenological stages according to growing degree day     
        !  call phenology(gdd_sum.....)                                    
          
       !   if (is_pheno_on()) then ! Only when crop/vegetation is present
            ! Update daily variables: LAI, Aboveground & Belowground biomass according to phenological stages
            ! also litter input if running yasso on hourly or daily step
      !    call carbon_allocation_day(pheno_stage, npp_sum_day.....)     
            
            ! Yasso: split input c into various yasso fractions
      !      call inputs_to_fractions(leaf, root, soluble, compost, fract)
          
      !    end if

          ! Yasso: create average meteorological forcings for yasso
      !    call average_met(met_daily, met_rolling, aver_size, met_state, met_ind)
          
          ! Yasso: update soil respiration, should we update soil respiration hourly or daily? 
          ! resp should be an output variable
      !    call decompose(param, timestep_days, c_input_awenh_day, nitr_input_day, tempr_c, &
      !                                precip_day, cstate, nstate, ctend, ntend) 

          ! We could also put yasso here can do the daily calculation for NEE
          ! call yasso20_day()
          ! nee_day= .......

          ! Write output for daily variables
          if(step_nc_day.eq.0)then
            !Initialize netcdf file
            call netCDF_prepareOUTPUT(output_filename_day, lon_sites, lat_sites, ntim_out_day)
          end if

         ! call netCDF_writeOUTPUT(output_filename_day, "HeteroResp", hr, tot_hour/24, step_nc_day)           
         ! step_nc_day= step_nc_day+1                

      !    gpp_sum_day =0
      !    npp_sum_day =0
      !    tr_day      =0
      !    evap_can_day=0
          
      !  end if

      !  if (is_end_curr_year() .or. is_pheno_harvest() ) then

      !    call carbon_allocation_yr(pheno_stage, gpp_sum_year, .....)   ! calculate yield, harvest biomass, 
                                                                   ! and litter input if running yasso  
          
          ! We could also put yasso here to do yearly calculation for NEE
          ! call yasso20_yr()

      !    call write_output_day()                   ! Write output for yearly variables

      !    gpp_sum_year=0
      !    npp_sum_year=0
      !    gdd_sum=0  

        end if

      end do ! m
    end do  !i

    ! advance time step
    tot_hour=tot_hour + time_step
  end do ! t

  !call write_restart()

end program SVMC


