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
  use water_mod                 ! module for canopy and soil water budget
  !use spafhy_mod                ! module for soil water bucket model, which will provide psi_soil for p-hydro
  
  !use alloc_mod                ! module for carbon allocation and yield
  use wrapper_yasso             ! module for soil decomposition model, which will provide heterogeneous respiration (hr)   
  use yasso
  use allocation 

  implicit none

  ! Loop variables
  !***********************************
  integer         :: i, m
  integer         :: step_nc_hr, step_nc_day, step_clim, step_lai, step_soilmoist, step_management
  real(kind=dp)   :: tot_hour, tot_hour_end, juldate, start_date, end_date
  real(8), dimension(1)     :: start_clim_time, end_clim_time
  real(8)            :: start_clim_juldate,start_lai_juldate,start_soilmoist_juldate, start_manage_juldate
                            
  real(8), dimension(1)    :: start_lai_time, end_lai_time 
  integer         :: ntim_clim, ntim_lai, ntim_out_hr, ntim_out_day
 
  ! epsilon
  real(8) :: eps = 1e-16

  !***********************************
  !Model variables
  !***********************************

  ! p-hydro variabless
  ! Input variables for p-hydro (the definition is from rpmodel.R)
  real(8)     ::    temp, temp_day, temp_year         ! Air temperature (tc), K
  real(8), dimension(12) :: temp_mon, month, mon_daynum
  real(8)     ::    ppfd      ! Photosynthetic photon flux density (mol m-2 d-1) (incoming solar radiation from forcing data?)
  real(8)     ::    vpd       ! Vapour pressure deficit (Pa) (will be calculated using pressure & humidity)
  real(8)     ::    co2       ! Atmospheric CO2 concentration (ppm)
  real(8)     ::    elv       ! Elevation above sea-level (m.a.s.l.) (not needed if we have surface pressure!)
  real(8)     ::    fapar     ! Fraction of absorbed photosynthetically active radiation (unitless) (will be calculated using LAI)
  real(8)     ::    prec, precip_day, precip_year       ! 
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
                                evap_matrix, psi_soil_matrix, soilmoist1_matrix, tmp_matrix

  real(8)     ::    lai, lai_alloc

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
  real(8)    :: gpp, gpp_day, npp_day, nee_day
  real(8)    :: tr_phydro   ! Transpiration estimated by p-hydro model

  ! spafhy variables
  real(8)    :: tr_spafhy   ! Transpiration estimated by spafhy model
  real(8)    :: latflow     ! lateral flow to ditches (placeholder) [m]
  real(8)    :: LE          ! latent heat flux [Wm-2]
  real(8)    :: LatentHeat  ! latent heat of vaporization [J kg-1]
  
  ! yasso variables
  real(8)    :: HeteroResp, AutoResp,TotalResp
  real(8)    :: leaf_litter_c, root_litter_c, soluble=0.0, compost=0.0
  real(8)    :: leaf_litter_c_year, root_litter_c_year, soluble_year=0.0, compost_year=0.0
  real(8)    :: metyasso_roll(2), metyasso(2)
  real(8)    :: metyasso_state(2,24*30+1)
  integer    :: metyasso_ind
  real(8)    :: tmp_input_cf1, tmp_input_cf2, tmp_input_cf3, tmp_input_cf4, tmp_input_cf5

  ! alloc variables
  real(8)    :: croot=0.0, cleaf=0.0, cstem=0.0
  real(8)    :: above_biomass, below_biomass, yield
  integer    :: pheno_stage=1, num_gpp_day=0

  ! For soil water retention curve
 ! real(8) :: watsat      ! v/v saturate moisture
  !real(8) :: watres      ! v/v, residual soil moisture for Van Genuchten
  !real(8) :: vol_ice     ! v/v, volumetric ice in soil bucket 
  !real(8) :: vol_liq     ! v/v, volumetric of liq in soil bucket     
  !real(8) :: satfrac     ! parameter for Van Genuchten
  !real(8) :: n1, m1, alpha_van           ! (-), pore-size-distribution parameter for Van Genuchten 1.07
  !real(8) :: eff_porosity! v/v, volume of ice

  type(soilwater_state_type)          :: soilwater_state
  type(soilwater_flux_type)           :: soilwater_flux  
  type(canopywater_state_type)        :: canopywater_state
  type(canopywater_flux_type)         :: canopywater_flux
  type(spafhy_para_type)              :: spafhy_para
  
  type(soilcn_state_type)       :: soilcn_state, soilcn_state0
  type(soilcn_flux_type)        :: soilcn_flux, soilcn_flux0
  type(yasso_para_type)         :: yasso_para
  type(alloc_para_type)         :: alloc_para
  type(management_data_type)    :: manage_data

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
   
  call readctrl_namelist
  call readvegpara_namelist

  call readsoilhydro_namelist(spafhy_para)
  call readsoilyasso_namelist(yasso_para)
  
  !call set_soilwaterState(soilwater_state, canopywater_state)
  call initialization_spafhy(canopywater_state, soilwater_state, spafhy_para)
  
  ! initialize yasso model: need temperature & precipitation input
  call wrapper_yasso_initialize_totc(soilcn_state, yasso_para)
  soilcn_state0=soilcn_state

  ! call initialize_yasso_totc(param, totc, cn_input, fract_root_input, fract_legacy_soc, &
  !    tempr_c, precip_day, tempr_ampl, cstate, nstate)

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
  metyasso_ind=1
  temp_day=0.0
  precip_day=0.0
  gpp_day=0.0

  temp_year=0.0
  temp_mon(:)=0.0
  month=(/31,59,90,120,151,181,212,243,273,304,334,365/)
  mon_daynum=(/31,28,31,30,31,30,31,31,30,31,30,31/)

  precip_year=0.0
  root_litter_c_year=0.0
  leaf_litter_c_year=0.0

  print *, tot_hour_end, ntim_out_hr, ntim_out_day
  ! Read time series of input data
  call netCDF_readTime(input_climfile, ntim_clim, start_clim_time, end_clim_time)
  call netCDF_readTime(input_laifile, ntim_lai, start_lai_time, end_lai_time)
  call netCDF_readlonlat(input_climfile, num_sites, lat_sites, lon_sites)

  ! To simplify the time management, specify the julian start date of the inputdata by hand.
  start_clim_juldate=juldate(20201231,233000)
  start_lai_juldate =juldate(20210101,000000)
  start_soilmoist_juldate =juldate(20210101,000000)
  start_manage_juldate =juldate(20210101,000000)

  ! step the starting time steps for reading input files
  step_clim = floor((start_date-start_clim_juldate)*24)+1
  step_lai  = floor(start_date-start_lai_juldate)+1
  step_soilmoist = floor(start_date-start_soilmoist_juldate)+1
  step_management = floor(start_date-start_manage_juldate)+1

  !********* open file for writing SpaFHy test outputs
  open(99, file = 'logbook.txt', status = 'old')
  !write(99,*) "prec,T,Wliq,WliqTop,PsiS,Mbe,tr,ground_evap,infil,drain,roff,pondsto,swe,&
  !    &swe_l,swe_i,canopy_evap,canopy_mbe,CanopyStorage,Trfall,PotInf, LE, tr_phydro, gs, WatSto, WatStoTop"
  write(99, *) "prec,T,Wliq,PsiS,soil_mbe,soil_ET,soil_Infiltration,soil_Drainage,soil_Runoff,soil_PondSto,&
                &canopy_SWE,canopy_swe_l,canopy_swe_i,canopy_CanopyEvap,canopy_SoilEvap,tr_spafhy,&
                &canopy_Interception,canopy_Throughfall,canopy_PotInfiltration,canopy_mbe,canopy_CanopyStorage,&
                &LE,tr_phydro,gs,soil_Kh,Melt,Freeze,beta_SoilEvap,canopy_Unloading"
                
  open(999, file = 'yassodebug.txt', status = 'old')
  write(999,*) "cstate1,cstate2, cstate3, cstate4, cstate5, nstate, &
              input_cfract1, input_cfract2, input_cfract3, input_cfract4, input_cfract5, input_nfract, &
              ctend1, ctend2, ctend3, ctend4, ctend5, ntend"
  !*******************
  
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
            
            if (obs_lai) then
              call netCDF_readlai(input_laifile, lai_matrix, step_lai)
              lai=lai_matrix(1,1,1)
              step_lai=step_lai+1
              ! calculate fapar:
              fapar= 1-exp(-k*lai)
            end if

            if (obs_soilmoist) then
              call netCDF_readsoilmoist('../data/FieldObs_Qvidja.2021.soilmoist.nc', soilmoist_matrix, step_soilmoist)
              soilmoist=soilmoist_matrix(1,1,1)               
              call soil_water_retention_curve(soilmoist, spafhy_para, psi_soil)
              step_soilmoist=step_soilmoist+1
            end if
            
            ! Read management information
            call netCDF_readmanagement('../data/management_qvidja_2021.nc', manage_data%management_type, & 
                            manage_data%management_c_input, manage_data%management_c_output, & 
                            manage_data%management_n_input, manage_data%management_n_output, step_management)
            step_management=step_management+1              
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
          if(.not. obs_soilmoist) then
            !psi_soil = -1.0
            psi_soil = soilwater_state%Psi !root zone, MPa
            !psi_soil = min(-eps, max(psi_soil, -3.0))  ! ensures psi_soil <0 and >-3.0 MPa

          end if

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
          print *, "psi_soil =", psi_soil         ! MPa

          
          ! for coupling with SpaFHy: psi_soil = soilwater_state%Psi
          call pmodel_hydraulics_numerical(temp-273.15, ppfd*1000000.0/lai, vpd, co2*1000000, pres, fapar, &
                                 psi_soil, rdark,                                                 &
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
          ! HUI - check units of gs and conversion to tr_phydro. We want it to be ≈ kg H2O m-2 s-1] 
          ! [tr_phydro] [mm s-1] = [1] * [mol/m2/s] * [Pa Pa-1] * [g mol-1] / [kg m-3]  
          ! Density of water is used here to more accurately convert the unit to mm -s
          if (ISNAN(gs)) then
            tr_phydro = 0.0
          else
             tr_phydro = 1.6*gs*(vpd/pres)*h2o_molmass/density_h2o(temp-273.15, pres) * lai
          end if
          
          ! net radiation of the whole canopy-soil system [W m-2]
          ! Samuli will revise later!
          rn= rg * 0.7
          !rn = max(2.57*lai/(2.57*lai+0.57)-0.2, 0.55)*rg  ! Launiainen et al. 2016 GCB, fit to Fig 2a     
          
          ! reset water fluxes to zero
          call reset_spafhy_flux(canopywater_flux, soilwater_flux)

          ! call SpaFHy code to compute new canopywater_state and snowwater_state
          ! returns water fluxes integrated over time_step in units [mm = kg H2O m-2]
          call canopy_water_flux(rn, temp-273.15, prec, vpd, wind, pres, fapar, lai, &
                                  canopywater_state, canopywater_flux, soilwater_state, spafhy_para)

          ! ET [mm]
          canopywater_flux%ET =  tr_phydro * (time_step*3600.0) +  canopywater_flux%SoilEvap + &
                                      canopywater_flux%CanopyEvap
          
          LatentHeat = 1.0e3 * (3147.5 - 2.37 * (temp))
          LE = canopywater_flux%ET / (time_step * 3600.0) * LatentHeat ! Wm-2

          ! Solve soil water balance
          tr_spafhy = tr_phydro*(time_step*3600.0) ! mm
          latflow=0.0
          ! water fluxes must be in units [m]. Updates soilwater_state, including soilwater_state%Psi.
          call soil_water(soilwater_state, soilwater_flux, spafhy_para, canopywater_flux%PotInfiltration, &
                          tr_spafhy, canopywater_flux%SoilEvap, latflow)

          !tr_spafhy=tr_phydro*(time_step*3600.0*1.0e-3)   ! This variable has to be used to be modified in soil_water
          !retflow=0.0
          ! water fluxes must be in units [m]. Updates soilwater_state, including soilwater_state%Psi.
          !call soil_water(soilwater_state, soilwater_flux, canopywater_flux%PotInf*1.0e-3, &
          !                tr_spafhy,  &
          !                canopywater_flux%GroundEvap*1.0e-3, retflow, spafhy_para)  

          !call soil_water_retention_curve(soilwater_state%Wliq, psi_soil_spafhy) 

        !end if  !phenology
        
        ! Calculation of NEE & LATENT heat flux
        ! nee= gpp-ar-hr

        ! Write hourly output at output time step frequency
        write(99,'(*(G0.6,:,","))') & 
        prec*time_step, temp - 273.15, soilwater_state%Wliq, soilwater_state%Psi, soilwater_flux%mbe, &
        soilwater_flux%ET, soilwater_flux%Infiltration, soilwater_flux%Drainage, soilwater_flux%Runoff, &
        soilwater_state%PondSto, canopywater_state%SWE, canopywater_state%swe_l, canopywater_state%swe_i, &
        canopywater_flux%CanopyEvap, canopywater_flux%SoilEvap, tr_spafhy, canopywater_flux%Interception, &
        canopywater_flux%Throughfall, canopywater_flux%PotInfiltration, canopywater_flux%mbe, & 
        canopywater_state%CanopyStorage, LE, tr_phydro, gs, soilwater_state%Kh, canopywater_flux%Melt, canopywater_flux%Freeze, &
        soilwater_state%beta,canopywater_flux%Unloading

        ! *** test output for de-bugging
        !write(99,'(*(G0.6,:,","))') & 
        !prec, temp - 273.15, soilwater_state%Wliq, soilwater_state%Wliq_top, &
        !soilwater_state%Psi, soilwater_state%mbe, tr_spafhy, &
        !canopywater_flux%GroundEvap*1.0e-3, soilwater_flux%Inflow, soilwater_flux%Drain, soilwater_flux%Roff, &
        !soilwater_state%PondSto, canopywater_state%swe, canopywater_state%SWEl, canopywater_state%SWEi, &
        !canopywater_flux%CanopyEvap*1e-3, &
        !canopywater_state%MBE, canopywater_state%CanopyStorage, canopywater_flux%Trfall, canopywater_flux%PotInf, &
        !LE, tr_phydro, gs, soilwater_state%WatSto, soilwater_state%WatStoTop

        if ( mod(tot_hour,time_step_output) .eq. 0.0 ) then
            
          !Write hourly output for this time step
          !************************************************************************
          if(step_nc_hr.eq.0)then
            !Initialize netcdf file
            print *, "ntim_out_hr=", ntim_out_hr
            call netCDF_prepareOUTPUT(output_filename_hr, lon_sites, lat_sites, ntim_out_hr)
          end if

          print *, "step_nc_hr=", step_nc_hr

          call netCDF_writeOUTPUT(output_filename_hr, "GPP", gpp, tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "stomatal_conductance", gs, tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "Jmax", jmax, tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "Vcmax", vcmax, tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "Chi", chi, tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "Dpsi", dpsi, tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "Profit", profit, tot_hour/24.0, step_nc_hr)

          !Use common unit [mm s-1] ≈ [kg H2O m-2 s-1] for flux
          call netCDF_writeOUTPUT(output_filename_hr, "Qle", LE, tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "Evap", canopywater_flux%ET/(time_step*3600.0), & 
                                      tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "Transp", tr_spafhy/(time_step*3600.0), &
                                      tot_hour/24, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "CanopyEvap", canopywater_flux%CanopyEvap/(time_step*3600.0), &
                                      tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "GroundEvap", canopywater_flux%SoilEvap/(time_step*3600.0), & 
                                      tot_hour/24, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "Trfall", canopywater_flux%Throughfall/(time_step*3600.0), &
                                      tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "CanopyInterc", canopywater_flux%Interception/(time_step*3600.0), &
                                      tot_hour/24, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "PotInf", canopywater_flux%PotInfiltration/(time_step*3600.0), & 
                                      tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "Unload", canopywater_flux%Unloading/(time_step*3600.0), & 
                                      tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "Roff", soilwater_flux%Runoff /(time_step*3600.0), & 
                                      tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "Drain", soilwater_flux%Drainage/(time_step*3600.0), & 
                                      tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "Inflow", soilwater_flux%LateralFlow/(time_step*3600.0), & 
                                      tot_hour/24.0, step_nc_hr)
          !call netCDF_writeOUTPUT(output_filename_hr, "TopSoilInterc", soilwater_flux%Interc/(time_step*3600.0), & 
          !                           tot_hour/24.0, step_nc_hr)                    
         
          !Use the units used by original spafhy: [m] for soil water storage, [mm] for canopy water storage
          call netCDF_writeOUTPUT(output_filename_hr, "WatSto", soilwater_state%WatSto, tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "PondSto", soilwater_state%PondSto, tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "WatStoTop", soilwater_state%WatSto, tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "SoilMoist", soilwater_state%Wliq, tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "SoilMoistPot", soilwater_state%Psi, tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "CanopyStorage", canopywater_state%CanopyStorage, &
                                      tot_hour/24.0, step_nc_hr)
          call netCDF_writeOUTPUT(output_filename_hr, "swe", canopywater_state%swe, tot_hour/24.0, step_nc_hr)
          
          step_nc_hr= step_nc_hr+1 
        endif
        
        ! Comulative parts can be a separate module in the future
        ! Cumulative GPP or NPP which will be used for carbon allocation on daily or yearly scale. 
        !gpp_sum_day=gpp_sum_day + gpp_hr*3600*...           
        !npp_sum_day=npp_sum_day + npp_hr*3600*...
        !gpp_sum_year=gpp_sum_year + gpp_hr*3600*...                   
        !npp_sum_year=npp_sum_year + npp_hr*3600*...
          
        ! Comulative temperature which will be used for GDD calculation and phenology
        !temp_day=temp_day+temp
        !precip_day=precip_day+prec*time_step*3600

        ! Cumulative transpiration and canopy evaporation
        !tr_day=tr_day + tr
        !evap_can_day=evap_can_day+evap_can

        ! Cumulative solar radiation (energy)

        ! Yasso: create average meteorological forcings for yasso
        ! HT: Currently, we use monthly rolling average for each hour calculation
        ! JV: compare annual balance with daily calculation
        metyasso(1)=temp-273.15
        metyasso(2)=prec


      
        ! 30-day moving averaging 
        ! call average_met(metyasso, metyasso_roll, 24*30, &
        !                           metyasso_state, metyasso_ind)

        ! Exponential smoothing
        call exponential_smooth_met(metyasso, metyasso_roll, metyasso_ind)                            
        
        temp_day=temp_day+metyasso_roll(1)
        precip_day=precip_day+metyasso_roll(2)*time_step*3600
        

        temp_year=temp_year + temp-273.15
        precip_year= precip_year + prec*time_step*3600


        write(999,'(*(G0.6,:,","))') & 
          metyasso(1), metyasso(2), metyasso_roll(1), metyasso_roll(2)

        if (ISNAN(gpp) .or. (gpp .lt. 0.0)) then
          gpp = 0.0
        else 
          num_gpp_day= num_gpp_day+1
        end if
        gpp_day=gpp_day + gpp
        
        if ((mod(tot_hour,24.0) .eq. 0).and.(tot_hour .gt. 0)) then    ! here assume the start time is always the beginning of the day!
          !Update GDD which is the criteria for phenology stages         
        !  gdd_sum= gdd_sum+temp_day/24.....

          temp_day=temp_day/24

          if (num_gpp_day .eq. 0) then
            gpp_day = 0.0
          else
            gpp_day = gpp_day/num_gpp_day
          end if

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
            
          !call litter_input(npp_sum, pheno_stage, management, leaf_litter_c, root_litter_c, root_exudent_c)  
          !HT: Some sensitivity test with harvest residues input
          !leaf_litter_c= gpp_day * 0.5 * 0.5 * 0.5 * 3600 * 24
          !root_litter_c= gpp_day * 0.5 * 0.5 * 0.5 * 3600 * 24
        
          !call alloc_hypothesis_1(gpp_day, npp_day,  leaf_litter_c, root_litter_c, alloc_para)
          !AutoResp=gpp_day*0.5*3600*24
          call alloc_hypothesis_2(gpp_day, npp_day, AutoResp, croot, cleaf, cstem, leaf_litter_c, root_litter_c, &
                                  compost, above_biomass, below_biomass, yield, &
                                  lai_alloc, alloc_para, manage_data, pheno_stage)
          
          leaf_litter_c_year = leaf_litter_c_year + leaf_litter_c
          root_litter_c_year = root_litter_c_year + root_litter_c
          compost_year       = compost_year + compost

          ! Yasso: split input c into various yasso fractions
      !   end if
          
          ! Yasso: update soil respiration, should we update soil respiration hourly or daily? 
          ! resp should be an output variable


          write(999,'(*(G0.6,:,","))') & 
             soilcn_state%cstate(1), soilcn_state%cstate(2), soilcn_state%cstate(3), soilcn_state%cstate(4), & 
             soilcn_state%cstate(5), soilcn_state%nstate, soilcn_flux%input_cfract(1), soilcn_flux%input_cfract(2), & 
             soilcn_flux%input_cfract(3), soilcn_flux%input_cfract(4), soilcn_flux%input_cfract(5), &
             soilcn_flux%input_nfract, soilcn_flux%ctend(1), soilcn_flux%ctend(2), &
             soilcn_flux%ctend(3), soilcn_flux%ctend(4), soilcn_flux%ctend(5), soilcn_flux%ntend, & 
             leaf_litter_c, root_litter_c, metyasso_roll(1),metyasso_roll(2)

          call wrapper_yasso_initialize_flux(soilcn_flux)
          call inputs_to_fractions(leaf_litter_c, root_litter_c, soluble, compost, soilcn_flux%input_cfract)

          write(999,'(*(G0.6,:,","))') & 
             soilcn_state%cstate(1), soilcn_state%cstate(2), soilcn_state%cstate(3), soilcn_state%cstate(4), & 
             soilcn_state%cstate(5), soilcn_state%nstate, soilcn_flux%input_cfract(1), soilcn_flux%input_cfract(2), & 
             soilcn_flux%input_cfract(3), soilcn_flux%input_cfract(4), soilcn_flux%input_cfract(5), &
             soilcn_flux%input_nfract, soilcn_flux%ctend(1), soilcn_flux%ctend(2), &
             soilcn_flux%ctend(3), soilcn_flux%ctend(4), soilcn_flux%ctend(5), soilcn_flux%ntend

          call wrapper_yasso_decompose(soilcn_state, soilcn_flux, yasso_para, 1.0, temp_day, precip_day)

          write(999,'(*(G0.6,:,","))') & 
             soilcn_state%cstate(1), soilcn_state%cstate(2), soilcn_state%cstate(3), soilcn_state%cstate(4), & 
             soilcn_state%cstate(5), soilcn_state%nstate, soilcn_flux%input_cfract(1), soilcn_flux%input_cfract(2), & 
             soilcn_flux%input_cfract(3), soilcn_flux%input_cfract(4), soilcn_flux%input_cfract(5), &
             soilcn_flux%input_nfract, soilcn_flux%ctend(1), soilcn_flux%ctend(2), &
             soilcn_flux%ctend(3), soilcn_flux%ctend(4), soilcn_flux%ctend(5), soilcn_flux%ntend

          tmp_input_cf1=tmp_input_cf1+soilcn_flux%input_cfract(1)
          tmp_input_cf2=tmp_input_cf2+soilcn_flux%input_cfract(2)
          tmp_input_cf3=tmp_input_cf3+soilcn_flux%input_cfract(3)
          tmp_input_cf4=tmp_input_cf4+soilcn_flux%input_cfract(4)
          tmp_input_cf5=tmp_input_cf5+soilcn_flux%input_cfract(5)


          HeteroResp= sum(-soilcn_flux%ctend)/24/3600  
          TotalResp=HeteroResp+AutoResp/24/3600
          
          ! We could also put yasso here can do the daily calculation for NEE
          ! call yasso20_day()
          nee_day=TotalResp - gpp_day

          ! Write output for daily variables
          if(step_nc_day.eq.0)then
            !Initialize netcdf file
            call netCDF_prepareOUTPUT(output_filename_day, lon_sites, lat_sites, ntim_out_day-1)
          end if
          
          call netCDF_writeOUTPUT(output_filename_day, "HeteroResp", HeteroResp, tot_hour/24.0, step_nc_day) 
          call netCDF_writeOUTPUT(output_filename_day, "AutoResp", AutoResp, tot_hour/24.0, step_nc_day)         
          call netCDF_writeOUTPUT(output_filename_day, "TotalResp", TotalResp, tot_hour/24.0, step_nc_day) 
          call netCDF_writeOUTPUT(output_filename_day, "GPP", gpp_day, tot_hour/24.0, step_nc_day)
          call netCDF_writeOUTPUT(output_filename_day, "NEE", nee_day, tot_hour/24.0, step_nc_day)
          call netCDF_writeOUTPUT(output_filename_day, "LAI", lai_alloc, tot_hour/24.0, step_nc_day) 
          call netCDF_writeOUTPUT(output_filename_day, "TotLivBiom", above_biomass + below_biomass, tot_hour/24.0, step_nc_day)
          call netCDF_writeOUTPUT(output_filename_day, "leaf_carbon_content", cleaf, tot_hour/24.0, step_nc_day)  
          call netCDF_writeOUTPUT(output_filename_day, "root_carbon_content", croot, tot_hour/24.0, step_nc_day) 
          call netCDF_writeOUTPUT(output_filename_day, "soil_carbon_content", sum(soilcn_state%cstate), tot_hour/24.0, step_nc_day) 
          call netCDF_writeOUTPUT(output_filename_day, "temperature_yasso", metyasso_roll(1), tot_hour/24.0, step_nc_day) 
          call netCDF_writeOUTPUT(output_filename_day, "precipitation_yasso", metyasso_roll(2), tot_hour/24.0, step_nc_day) 

          step_nc_day= step_nc_day+1
          
          if(step_nc_day <= month(1)) then
            temp_mon(1) = temp_mon(1) + temp_day
          else if (step_nc_day <= month(2)) then
            temp_mon(2) = temp_mon(2) + temp_day
          else if (step_nc_day <= month(3)) then
            temp_mon(3) = temp_mon(3) + temp_day
          else if (step_nc_day <= month(4)) then
            temp_mon(4) = temp_mon(4) + temp_day
          else if (step_nc_day < month(5)) then
            temp_mon(5) = temp_mon(5) + temp_day
          else if (step_nc_day < month(6)) then
            temp_mon(6) = temp_mon(6) + temp_day
          else if (step_nc_day < month(7)) then
            temp_mon(7) = temp_mon(7) + temp_day
          else if (step_nc_day < month(8)) then
            temp_mon(8) = temp_mon(8) + temp_day
          else if (step_nc_day < month(9)) then
            temp_mon(9) = temp_mon(9) + temp_day
          else if (step_nc_day < month(10)) then
            temp_mon(10) = temp_mon(10) + temp_day
          else if (step_nc_day < month(11)) then
            temp_mon(11) = temp_mon(11) + temp_day
          else if (step_nc_day < month(12)) then
            temp_mon(12) = temp_mon(12) + temp_day
          end if

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
          
          gpp_day=0.0
          temp_day=0.0
          precip_day=0.0
          num_gpp_day=0

        end if

        if ((mod(tot_hour,24*363.0) .eq. 0).and.(tot_hour .gt. 0)) then
          
          temp_mon(:)=temp_mon(:)/mon_daynum(:)
          
          call wrapper_yasso_initialize_flux(soilcn_flux0)
          call inputs_to_fractions(leaf_litter_c_year, root_litter_c_year, soluble_year, compost_year, soilcn_flux0%input_cfract)
          call wrapper_yasso_annual(soilcn_state0, soilcn_flux0, yasso_para, real(step_nc_day+1), &
                                                        temp_mon, precip_year)

          write(999,'(*(G0.6,:,","))') & 
            soilcn_state0%cstate(1), soilcn_state0%cstate(2), soilcn_state0%cstate(3), soilcn_state0%cstate(4), & 
            soilcn_state0%cstate(5), leaf_litter_c_year, root_litter_c_year, step_nc_day, soilcn_flux0%input_cfract, &
            temp_mon, precip_year, sum(soilcn_flux0%ctend)/(step_nc_day-1)/24/3600

          write(999,'(*(G0.6,:,","))') & 
            tmp_input_cf1, tmp_input_cf2, tmp_input_cf3, tmp_input_cf4, tmp_input_cf5

          temp_mon(:)=0.0
          !soilcn_flux0%ctend/(step_nc_day-1)/24/3600

        end if

      end do ! m
    end do  !i

    ! advance time step
    tot_hour=tot_hour + time_step
  end do ! t

  !call write_restart()


  !call wrapper_yasso_initialize_flux(soilcn_flux)
  !call inputs_to_fractions(leaf_litter_c_year, root_litter_c_year, soluble_year, compost_year, soilcn_flux%input_cfract)
  !call wrapper_yasso_decompose(soilcn_state0, soilcn_flux, yasso_para, real(step_nc_day+1), &
  !                                                   temp_year/(tot_hour), precip_year)

  !write(99,'(*(G0.6,:,","))') & 
  !     soilcn_state0%cstate(1), soilcn_state0%cstate(2), soilcn_state0%cstate(3), soilcn_state0%cstate(4), & 
  !     soilcn_state0%cstate(5), soilcn_state0%nstate, soilcn_flux%input_cfract(1), soilcn_flux%input_cfract(2), & 
  !     soilcn_flux%input_cfract(3), soilcn_flux%input_cfract(4), soilcn_flux%input_cfract(5), &
  !     soilcn_flux%input_nfract, soilcn_flux%ctend(1), soilcn_flux%ctend(2), &
  !     soilcn_flux%ctend(3), soilcn_flux%ctend(4), soilcn_flux%ctend(5), soilcn_flux%ntend

end program SVMC


