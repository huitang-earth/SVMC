MODULE spafhy_mod

!--------------------------------------------------------------
! References:
! Launiainen, S., Guan, M., Salmivaara, A., Kieloaho, A.-J., 2019. 
! Modeling boreal forest evapotranspiration and water balance at stand and catchment scales: a spatial approach. 
! Hydrology and Earth System Sciences 23, 3457–3480. https://doi.org/10.5194/hess-23-3457-2019
!--------------------------------------------------------------

! Modules
  use netcdf             ! library for processing netcdf files
  use readctrl_mod       ! module for reading control parameters
  use readvegpara_mod    ! module for reading vegetation parameters
  use readsoilpara_mod   ! module for reading soil properties (shared with yasso)

  implicit none

  !Public member functions:
  public :: initialization_spafhy
  public :: canopy_water_flux   ! require inputdata for p-hydro module
  !public :: topmodel   !
  public :: soil_water
  public :: soil_water_retention_curve

  private :: aerodynamics        ! aerodynamic conductances
  private :: canopy_water_snow   ! canopy interception, evaporation and snowpack
  private :: ground_evaporation  ! ground evaporation from top soil layer

  !private :: dry_canopy_et       ! Computes ET from 2-layer canopy in absense of intercepted precipitiation
  !                               ! This will be replaced by p-hydro when p-hydro is turned on
  private :: penman_monteith      !
  private :: e_sat                !
  private :: set_soilwaterState   !
  !private :: hydrCond             ! 
  !private :: relative_evaporation !

contains

  SUBROUTINE initialization_spafhy(canopywater_state, soilwater_state, spafhy_para)
    type(soilwater_state_type), intent(inout)    :: soilwater_state
    type(canopywater_state_type), intent(inout)  :: canopywater_state
    type(spafhy_para_type), intent(in)    :: spafhy_para
        
    soilwater_state%MaxWatSto = spafhy_para%soil_depth * spafhy_para%max_poros
    soilwater_state%MaxStoTop = spafhy_para%org_depth * spafhy_para%org_fc
    ! initial state
    soilwater_state%WatSto  = 0.9 * soilwater_state%MaxWatSto
    soilwater_state%WatStoTop = 0.9 * soilwater_state%MaxStoTop
    soilwater_state%PondSto = 0.0 ! no ponding

    call set_soilwaterState(soilwater_state, spafhy_para)

    ! canopywater
    canopywater_state%CanopyStorage=0.0

    canopywater_state%swe=0.0
    canopywater_state%SWEi=0.0
    canopywater_state%SWEl=0.0

  END SUBROUTINE initialization_spafhy
  
  SUBROUTINE soil_water(soilwater_state, soilwater_flux, rr, tr, evap, retflow, spafhy_para)
  ! ---------------------------------------------------------------
  !      Computes 2-layer bucket model water balance for one timestep dt
  !      Top layer is interception storage and contributes only to evap.
  !      Lower layer is rootzone and contributes only tr and creates drainage.
  !      Capillary interaction between layers is neglected and connection from bottom up
  !      is only in case of excess returnflow.
  !      Pond storage can exist above top layer.   
  !      IN:
  !          soilwater_state
  !          rr = potential infiltration [m]
  !          tr = transpiration from root zone [m]
  !          evap = evaporation from top layer [m]
  !          retflow = return flow from ground water [m]
  ! ----------------------------------------------------------------

  ! !ARGUMENTS:
    type(soilwater_state_type), intent(inout)    :: soilwater_state    
    type(soilwater_flux_type), intent(inout)    :: soilwater_flux    
    type(spafhy_para_type), intent(in)    :: spafhy_para
    real(8)      , intent(in)    :: rr      ! potential infiltration [m]
    real(8)      , intent(inout) :: tr      ! transpiration from root zone [m]
    real(8)      , intent(in)    :: evap    ! evaporation from top layer [m]
    real(8)      , intent(in)    :: retflow ! return flow from ground water [m]


  ! ! Local variables 
    real(8)     ::    rr0, rr1,evap1, PondSto0, WatSto0, WatStoTop0
    real(8)     ::    eps=1e-16
    real(8)     ::    Qin, dSto, exfil, to_pond, to_top_layer
    
    ! storages are in [m] and input fluxes in [m during timestep]
    ! m H2O = 1e-3 mm H2O = 1 g H2O m-2
  
    ! old state
    WatSto0 = soilwater_state%WatSto
    WatStoTop0 = soilwater_state%WatStoTop
    PondSto0 = soilwater_state%PondSto

    ! add current Pond storage to rr & update storage
    rr1 = rr + soilwater_state%PondSto
    soilwater_state%PondSto = 0.0
        
    ! top layer:
    ! interception
    soilwater_flux%Interc = max(0.0, (soilwater_state%MaxStoTop - soilwater_state%WatStoTop)) &
                    * (1.0 - exp(-(rr1 / soilwater_state%MaxStoTop)))

    soilwater_state%WatStoTop = max(0.0, soilwater_state%WatStoTop + soilwater_flux%Interc)  
    ! evaporation 
    evap1 = min(evap, soilwater_state%WatStoTop)
    soilwater_state%WatStoTop = soilwater_state%WatStoTop - evap1
      
    ! infiltration to rootzone
    rr1 = rr1 - soilwater_flux%Interc
                
    ! root one
    ! transpiration
    tr = min(tr, soilwater_state%WatSto - eps)
    soilwater_state%WatSto = soilwater_state%WatSto - tr
        
    ! drainage: if retflow from groundwater storage, set drain to zero.

    if (retflow .gt. 0.0) then 
     soilwater_flux%Drain = 0.0
    else
      soilwater_flux%Drain = min(soilwater_state%Kh * (time_step * 3600.0), & ! conductivity
                               max(0.0, (soilwater_state%Wliq - spafhy_para%fc))* spafhy_para%soil_depth) ! available water
    end if 
    
    ! inflow to root zone: restricted by potential inflow or available pore space
    Qin = (retflow + rr1)         ! m, pot. inflow
    soilwater_flux%Inflow = min(Qin, soilwater_state%MaxWatSto - soilwater_state%WatSto + soilwater_flux%Drain)    
    dSto = (soilwater_flux%Inflow - soilwater_flux%Drain)
    soilwater_state%WatSto = min(soilwater_state%MaxWatSto, max(soilwater_state%WatSto + dSto, eps))
                
    ! if inflow excess after filling rootzone, update first top layer storage
    exfil = Qin - soilwater_flux%Inflow
    to_top_layer = min(exfil, soilwater_state%MaxStoTop - soilwater_state%WatStoTop - eps)
    ! self.WatStoTop = self.WatStoTop + to_top_layer
    soilwater_state%WatStoTop = soilwater_state%WatStoTop + to_top_layer
        
    ! ... and then pond storage ...
    to_pond = min(exfil - to_top_layer, spafhy_para%maxpond - soilwater_state%PondSto - eps)
    soilwater_state%PondSto = soilwater_state%PondSto + to_pond
 
    ! ... and route remaining to surface runoff
    soilwater_flux%Roff = exfil - to_top_layer - to_pond

    ! update soil water state variables
    call set_soilwaterState(soilwater_state, spafhy_para)

    ! mass balance error [m]
    soilwater_state%mbe = (soilwater_state%WatSto - WatSto0)  &
                     + (soilwater_state%WatStoTop - WatStoTop0) &
                     + (soilwater_state%PondSto - PondSto0) &
                     - (rr + retflow - tr - evap1 - soilwater_flux%Drain - soilwater_flux%Roff)
  
  END SUBROUTINE soil_water

  SUBROUTINE soil_water_retention_curve(vol_liq, smp)
  ! converts vol. water content to soil water potential (in MPa)
  ! To do: 
  !      - fc, wp and kh need to be derived in this soubroutine

    implicit none
    real(8), intent(in) :: vol_liq        ! v/v, volumetric of liq in soil bucket
    real(8), intent(out):: smp            !soil suction, negative, MPa
    !real(r8), intent(out)            :: hk      !hydraulic conductivity [mm/s]
    !real(8), optional, intent(out)  :: dsmpds   !d[smp]/ds, [mm]
    !real(r8), optional, intent(out)  :: dhkds    !d[hk]/ds   [mm/s]
  !
  ! !LOCAL VARIABLES:

    real(8) :: vol_ice     ! v/v, volumetric ice in soil bucket 
    real(8) :: satfrac     ! parameter for Van Genuchten
    real(8) :: n1, m1, alpha_van, watsat,watres         ! (-), pore-size-distribution parameter for Van Genuchten 1.07
    real(8) :: eff_porosity! v/v, volume of ice
    
  !--- Van Genuchten scheme based on FatesHydro (parameters from Launiainen et al. Forests, 2022)
  ! MOVE TO PARAMETERS - reading alpha_sat from namelist does not work!!!
    n1= 1.07 !n_van       !Launiainen et al. 2022: C1-5: 1.12, 1.14, 1.07, 1.27, 1.18  
    m1=1.0/n1  
    alpha_van=2.02   !Launiainen et al. 2022: C1-5: 4.45, 5.92, 2.02, 4.49, 3.35
    watsat=0.46  !Launiainen et al. 2022: C1-5: 0.75, 0.68, 0.46, 0.47, 0.54  
    watres=0.0
  
    vol_ice = 0.0  
    eff_porosity = max(0.01, watsat - vol_ice)
    
    satfrac = (vol_liq-watres)/(eff_porosity-watres)
    smp = -(1.0/alpha_van)*(satfrac**(1.0/(m1-1.0)) - 1.0 )**m1 !kPa
    smp = smp * 0.001 !MPa

  END SUBROUTINE soil_water_retention_curve


  SUBROUTINE canopy_water_flux(Rn, Ta, Prec, VPD, U, P, fapar, LAI,            & 
                             canopywater_state, canopywater_flux, soilwater_state, spafhy_para)
    !
    ! Computes plant canopy interception, throughfall, ground evaporation, snowpack dynamics
    ! NOTE: Re-think canopy snow interception when snow depth > canopy height 
  
    ! !ARGUMENTS:
    real(8)      , intent(in)    :: Rn     ! Net solar radiation obsorbed by canopy & soil (W m-2)
    real(8)      , intent(in)    :: Ta     ! Air temperature (tc), (degrees C)
    real(8)      , intent(in)    :: Prec   ! precipitatation rate [mm/s] = kg m-2 s-1
    ! real(8)      , intent(in)    :: Par    ! photos. act. radiation [Wm-2]
    real(8)      , intent(in)    :: VPD    ! Vapour pressure deficit (Pa) (will be calculated using pressure & humidity)
    real(8)      , intent(in)    :: U      ! mean wind speed at ref. height above canopy top [ms-1] 
    real(8)      , intent(in)    :: P      ! pressure [Pa], scalar or matrix
    real(8)      , intent(in)    :: fapar  ! fraction of canopy absorbed PAR [-]
    real(8)      , intent(in)    :: LAI    ! leaf area index

    type(canopywater_state_type), intent(inout)   :: canopywater_state  ! canopy state
    type(canopywater_flux_type), intent(inout)    :: canopywater_flux    ! snowpack state
    type(soilwater_state_type), intent(in)        :: soilwater_state  ! soilwater state (for ground evaporation) 
    type(spafhy_para_type), intent(in)    :: spafhy_para

    ! ! Local variables 
    real(8)     ::    Ra, Rb, Ras, ustar, Uh, Ug, fPheno, AE
  
    ! Calculate aerodynamic resistances for canopy and soil layers
    call aerodynamics(LAI, U, Ra, Rb, Ras, ustar, Uh, Ug, spafhy_para)

    ! Calculate canopy interception, canopy evaporation, snowpack dynamics, ground evaporation
    AE = Rn * fapar
    call canopy_water_snow(canopywater_state, canopywater_flux, spafhy_para, Ta, Prec, AE, VPD, Ra, U, LAI, P)

    ! Calculate soil evaporation rate
    AE = Rn * (1 - fapar)
    print *, "fapar=", fapar

    call ground_evaporation(canopywater_state, canopywater_flux, soilwater_state, spafhy_para, Ta, AE, VPD, Ras, P)

  END SUBROUTINE canopy_water_flux


  SUBROUTINE ground_evaporation(canopywater_state, canopywater_flux, soilwater_state, spafhy_para, T, AE, VPD, Ras, P)
    !
    ! Calculates evaporation from top soil layer [mm]

    ! !ARGUMENTS:
    real(8)      , intent(in)    :: T      ! air temperature (degC)
    real(8)      , intent(in)    :: AE     ! available energy (~net radiation) (Wm-2)
    real(8)      , intent(in)    :: VPD    ! vapor pressure deficit (Pa)
    real(8)      , intent(in)    :: Ras     ! ground aerodynamic resistance (s m-1)
    real(8)      , intent(in)    :: P      ! pressure [Pa], scalar or matrix
    type(canopywater_flux_type), intent(inout)   :: canopywater_flux
    type(soilwater_state_type), intent(in)          :: soilwater_state  ! soilwater state (for ground evaporation)
    type(canopywater_state_type), intent(in)   :: canopywater_state
    type(spafhy_para_type), intent(in)    :: spafhy_para

    ! ! Local variables 
    real(8)     ::    Lv, erate, Gas, eps=1E-16
    Lv = 1.0e3 * (3147.5 - 2.37 * (T + 273.15))
    Gas = 1 / Ras
    ! gsoil is soil surface conductance when fully wet, defined in readvegpara

    erate = (time_step * 3600) * soilwater_state%beta * penman_monteith(AE, VPD, T, spafhy_para%gsoil, Gas, P) / Lv ! mm

    ! maximum equals available water 
    canopywater_flux%GroundEvap = min(1e3*soilwater_state%WatStoTop, erate)

    if (canopywater_state%swe>eps) then
      canopywater_flux%GroundEvap = 0.0  ! no evaporation from floor if snow on ground
    end if
  END SUBROUTINE ground_evaporation

  SUBROUTINE canopy_water_snow(canopywater_state, canopywater_flux, spafhy_para, T, Pre, AE, D, Ra, U, LAI, P)
    !
    ! Calculates canopy interception, throughfall and snowpack change during timestep dt
    ! Updates canopy and snow storages

    ! !ARGUMENTS:
    real(8)      , intent(in)    :: T      ! air temperature [deg C]
    real(8)      , intent(in)    :: Pre    ! precipitation rate during [mm s-1]
    real(8)      , intent(in)    :: AE     ! available energy (~net radiation) [W m-2]
    real(8)      , intent(in)    :: D      ! vapor pressure deficit [Pa]
    real(8)      , intent(in)    :: Ra     ! canopy aerodynamic resistance [s m-1]
    real(8)      , intent(in)    :: U      ! mean wind speed at ref. height above canopy top [ms-1]
    real(8)      , intent(in)    :: LAI    ! leaf area index [m2 m-2]
    real(8)      , intent(in)    :: P      ! pressure [Pa]
    type(canopywater_state_type), intent(inout)   :: canopywater_state
    type(canopywater_flux_type), intent(inout)    :: canopywater_flux
    type(spafhy_para_type), intent(in)    :: spafhy_para
    
    ! Local variables.
    real(8)  :: fW, fS, Tmin, Tmax, Tmelt, wmax_tot, wmaxsnow_tot
    real(8)  :: Ga, Ce, Sh, gi, erate, gs, Sice, Sliq
    real(8)  :: Melt, Freeze, Lv, Ls, SWEo, Wo, Prec
    real(8)  :: eps = 1e-16

    ! quality of precipitation depends on temperature
    Tmin = 0.0  ! 'C, below all is snow
    Tmax = 1.0  ! 'C, above all is water
    Tmelt = 0.0  ! 'C, T when melting starts

    ! state of precipitation [as water (fW) or as snow(fS)]
    if ( T <= Tmin) then
      fS = 1.0
    else if (T >= Tmax) then
      fW = 1.0
    else if ((T > Tmin) .and. (T < Tmax)) then
      fW = (T - Tmin) / (Tmax - Tmin)
      fS = 1.0 - fW
    end if

    !canopy storage capacities [mm]
    wmax_tot     = spafhy_para%wmax * LAI
    wmaxsnow_tot = spafhy_para%wmaxsnow * LAI

    ! latent heat of vaporization (Lv) and sublimation (Ls) J kg-1
    Lv = 1.0e3 * (3147.5 - 2.37 * (T + 273.15))
    Ls = Lv + 3.3e5

    Prec = Pre * time_step * 3600  ! mm during timestep
  
    Ga = 1.0 / Ra ! aerodynamic conductance

    ! resistance for snow sublimation adopted from:
    ! Pomeroy et al. 1998 Hydrol proc; Essery et al. 2003 J. Climate;
    ! Best et al. 2011 Geosci. Mod. Dev.

    Ce = 0.01*((canopywater_state%CanopyStorage + eps) / wmaxsnow_tot)**(-0.4)  ! exposure coeff (-)
    Sh = (1.79 + 3.0*U**0.5)                      ! Sherwood numbner (-)
    gi = Sh*canopywater_state%CanopyStorage * Ce / 7.68 + eps                ! m s-1

    if ((Prec == 0) .and. (T <= Tmin)) then
      ! sublimation
      erate =  (time_step * 3600) / Ls * penman_monteith(AE, D, T, gi, Ga, P) ! mm in timestep
    
    else if ((Prec == 0) .and. (T > Tmin)) then
      ! evaporation
      gs = 1e6   ! set to large number for free evaporation from wet surface
      erate =  (time_step * 3600) / Lv * penman_monteith(AE, D, T, gs, Ga, P)  ! mm in timestep
    end if 

    !----- Initial conditions for calculating mass balance error
    Wo = canopywater_state%CanopyStorage     ! canopy storage, mm
    SWEo = canopywater_state%swe    ! Snow water equivalent mm

    !----- Canopy water storage change
    ! snow unloading from canopy, ensures also that seasonal LAI development does not mess up computations
    ! HT: (1) unloading first at each timestep? (2) why no snow unload? do we need it here?
    ! SL: below describes snow unloading. If T > Tmin, maximum storage is  that of liquid water.
    if (T >= Tmin) then
      canopywater_flux%Unload = max(canopywater_state%CanopyStorage - wmax_tot, 0.0)
      canopywater_state%CanopyStorage = canopywater_state%CanopyStorage - canopywater_flux%Unload
    end if

    !----- Interception of rain or snow: asymptotic approach of saturation.
    !      based on: Hedstrom & Pomeroy 1998. Hydrol. Proc 12, 1611-1625;
    !                Koivusalo & Kokkonen 2002 J.Hydrol. 262, 145-164.
    if (T < Tmin) then
      canopywater_flux%Interc = (wmaxsnow_tot- canopywater_state%CanopyStorage) &
                * (1.0 - exp(-Prec/wmaxsnow_tot))
    end if
        
    ! Above Tmin, interception capacity equals that of liquid precip
    if (T >= Tmin) then
      canopywater_flux%Interc = max(0.0, (wmax_tot - canopywater_state%CanopyStorage)) &
                * (1.0 - exp(-Prec/wmax_tot))
    end if

    ! update canopy storage after interception
    canopywater_state%CanopyStorage =  canopywater_state%CanopyStorage + canopywater_flux%Interc  ! new canopy storage, mm
    canopywater_flux%Trfall = Prec + canopywater_flux%Unload - canopywater_flux%Interc  ! Throughfall to field layer or snowpack

    ! evaporate from canopy and update storage
    canopywater_flux%CanopyEvap = min(erate,  canopywater_state%CanopyStorage + eps)  ! mm
    canopywater_state%CanopyStorage = canopywater_state%CanopyStorage - canopywater_flux%CanopyEvap

    !---- Snowpack (in case no snow, all Trfall routed to floor) """
    if (T >= Tmelt) then
      Melt = min(canopywater_state%SWEi, spafhy_para%kmelt * (time_step*3600) * (T - Tmelt))  ! mm
      Freeze = 0.0
    else if (T < Tmelt) then
      Freeze = min(canopywater_state%SWEl, spafhy_para%kfreeze * (time_step*3600) * (Tmelt - T))  ! mm
      Melt = 0.0
    end if

    !---- amount of water as ice and liquid in snowpack
    Sice = max(0.0, canopywater_state%SWEi + fS * canopywater_flux%Trfall + Freeze - Melt)
    Sliq = max(0.0, canopywater_state%SWEl + fW * canopywater_flux%Trfall - Freeze + Melt)

    ! The water that can not be hold by snow will penetrate to soil
    canopywater_flux%PotInf = max(0.0, Sliq - Sice * spafhy_para%frac_snowliq)  ! mm,
    Sliq   = max(0.0, Sliq - canopywater_flux%PotInf)  ! mm, liquid water in snow

    ! update Snowpack state variables
    canopywater_state%SWEl = Sliq
    canopywater_state%SWEi = Sice
    canopywater_state%swe  = canopywater_state%SWEl + canopywater_state%SWEi
        
    ! mass-balance error mm
    canopywater_state%MBE = (canopywater_state%CanopyStorage + canopywater_state%swe) - & 
                               (Wo + SWEo) - (Prec - canopywater_flux%CanopyEvap - & 
                               canopywater_flux%PotInf)

  END SUBROUTINE canopy_water_snow


  SUBROUTINE aerodynamics(LAI, Uo, ra, rb, ras, ustar, Uh, Ug, spafhy_para)
    !
    ! computes wind speed at ground and canopy + boundary layer conductances
    ! Computes wind speed at ground height assuming logarithmic profile above and
    ! exponential within canopy
    ! SOURCE:
    !   Cammalleri et al. 2010 Hydrol. Earth Syst. Sci
    !   Massman 1987, BLM 40, 179 - 197.
    !   Magnani et al. 1998 Plant Cell Env.

    ! !ARGUMENTS:
    real(8)      , intent(in)    :: LAI    ! one-sided leaf-area /plant area index (m2m-2)
    real(8)      , intent(in)    :: Uo     ! mean wind speed at reference height zm (ms-1)
    type(spafhy_para_type), intent(in) :: spafhy_para

    real(8)      , intent(out)   :: ra         ! canopy aerodynamic resistance (sm-1)
    real(8)      , intent(out)   :: rb         ! canopy boundary layer resistance (sm-1)
    real(8)      , intent(out)   :: ras        ! forest floor aerod. resistance (sm-1)
    real(8)      , intent(out)   :: ustar      ! friction velocity (ms-1)
    real(8)      , intent(out)   :: Uh         ! wind speed at hc (ms-1)
    real(8)      , intent(out)   :: Ug         ! wind speed at zg (ms-1)

    ! local variables
    real(8) :: zm1, zg1, alpha1, d, zom, zov, zosv, zn
    real(8) :: kv = 0.4  ! von Karman constant (-)
    real(8) :: beta_aero=285.0   ! s/m, from Campbell & Norman eq. (7.33) x 42.0 molm-3

    zm1 = spafhy_para%hc + spafhy_para%zmeas  ! m
    zg1 = min(spafhy_para%zground, 0.1 * spafhy_para%hc)
    alpha1 = LAI / 2.0  ! wind attenuation coeff (Yi, 2008 eq. 23)
    d = 0.66*spafhy_para%hc     ! displacement height [m]
    zom = 0.123*spafhy_para%hc  ! roughness lenght for momentum [m]
    zov = 0.1*zom   ! scalar roughness length [m]
    zosv = 0.1*spafhy_para%zo_ground ! soil scalar roughness length [m]

    ! solve ustar and U(hc) from log-profile above canopy
    ustar = Uo * kv / log((zm1 - d) / zom) 
    Uh = ustar / kv * log((spafhy_para%hc - d) / zom)
    
    ! U(zg) from exponential wind profile
    zn = min(zg1 / spafhy_para%hc, 1.0)  ! zground can't be above canopy top
    Ug = Uh * exp(alpha1*(zn - 1.0))

    ! canopy aerodynamic & boundary-layer resistances (sm-1). Magnani et al. 1998 PCE eq. B1 & B5
    !ra = 1. / (kv*ustar) * log((zm - d) / zom)
    ra = 1./(kv**2.0 * Uo) * log((zm1-d)/zom) * log((zm1-d)/zov)    
    rb = 1./LAI * beta_aero * ((spafhy_para%w_leaf / Uh)*(alpha1/(1.0-exp(-alpha1/2.0))))**0.5

    ! soil aerodynamic resistance (sm-1)
    ras = 1.0/(kv**2.0*Ug) * (log(spafhy_para%zground/spafhy_para%zo_ground))*log(spafhy_para%zground/(zosv))
    
    ra = ra + rb

  END SUBROUTINE aerodynamics


  real(8) FUNCTION penman_monteith(AE, D, T, Gs, Ga, P)

    !------------------------------------------------------
    !  Computes latent heat flux LE (Wm-2) using Penman-Monteith equation
    !  INPUT:
    !     AE - available energy [Wm-2]
    !     D  - vapor pressure deficit [Pa]
    !     T  - ambient air temperature [degC]
    !     Gs - surface (or stomatal) conductance [ms-1]
    !     Ga - aerodynamic conductance [ms-1]
    !     P  - ambient pressure [Pa]
    !  OUTPUT:
    !     x - latent heat flux LE [W m-2]
    !------------------------------------------------------ 
  
    ! !ARGUMENTS
    real(8), intent(in) :: AE                        ! available energy [Wm-2]
    real(8), intent(in) :: D                       ! vapor pressure deficit [Pa]
    real(8), intent(in) :: T                         ! ambient air temperature [degC]
    real(8), intent(in) :: Gs                        ! surface conductance [ms-1]
    real(8), intent(in) :: Ga                        ! aerodynamic conductance [ms-1]
    real(8), intent(in) :: P                         ! ambient pressure [Pa] 
    !character (*), intent(in) :: units          ! W (Wm-2), mm (mms-1=kg m-2 s-1), mol (mol m-2 s-1)

    ! !LOCAL VARIABLES:
    real(8)             :: cp, rho, Mw, s, g, esat, L !P

    ! --- constants
    cp = 1004.67  ! J kg-1 K-1
    rho= 1.25    ! kg m-3
    Mw = 18e-3    ! kg mol-1
    !P = 10130.0  ! standard sea-level pressure (Pa)
  
    call e_sat(T, P, s, g, esat)  ! slope of sat. vapor pressure, psycrom const

    L = 1e3 * (3147.5 - 2.37 * (T + 273.15))

    penman_monteith = (s * AE + rho * cp * Ga * D) / (s + g * (1.0 + Ga / Gs))  ! Wm-2

    ! if (units == 'mm') then
    !   penman_monteith = penman_monteith/L     ! kgm-2s-1 = mms-1
    ! else if (units == 'mol') then
    !   penman_monteith = penman_monteith/L/Mw  ! mol m-2 s-1
    ! end if

    penman_monteith = max(penman_monteith, 0.0)

    return

  END FUNCTION penman_monteith

  
  SUBROUTINE e_sat(T, P, s, g, esat)
  !-------------------------------------------
  !  Computes saturation vapor pressure (Pa), slope of vapor pressure curve
  !  [Pa K-1]  and psychrometric constant [Pa K-1]
  !  IN:
  !      T - air temperature (degC)
  !      P - ambient pressure (Pa)
  !  OUT:
  !      esa - saturation vapor pressure in Pa
  !      s - slope of saturation vapor pressure curve (Pa K-1)
  !      g - psychrometric constant (Pa K-1)
  !-------------------------------------------

  ! !ARGUMENTS
    real(8), intent(in) :: T                         ! air temperature (degC)
    real(8), intent(in) :: P                         ! ambient pressure (Pa)
    real(8), intent(out):: esat                       ! saturation vapor pressure in Pa
    real(8), intent(out):: s                         ! slope of saturation vapor pressure curve (Pa K-1)
    real(8), intent(out):: g                         ! psychrometric constant (Pa K-1)
    
  ! !LOCAL VARIABLES:
    real(8)             ::  NT = 273.15
    real(8)             ::  cp = 1004.67             ! J/kg/K
    real(8)             ::  Lambda
    real(8)             ::  eps = 1e-16

    Lambda = 1.0e3 * (3147.5 - 2.37 * (T + NT))         ! lat heat of vapor [J/kg]
    esat = 1.0e3 * (0.6112 * exp((17.67 * T) / (T + 273.16 - 29.66)))  ! Pa

    s = 17.502 * 240.97 * esat / ((240.97 + T) ** 2)
    g = P * cp / (0.622 * Lambda)

  END SUBROUTINE e_sat


  SUBROUTINE set_soilwaterState(soilwater_state, spafhy_para)
  !-------------------------------------------
  ! Updates soil water state variables
  !-------------------------------------------
    ! !ARGUMENTS
    type(soilwater_state_type), intent(inout) :: soilwater_state                         ! 
    type(spafhy_para_type), intent(in)  :: spafhy_para
  ! !LOCAL VARIABLES:  
    real(8)   :: eps=1e-16
    real(8)   :: psis

  ! root zone
    soilwater_state%Wliq =  spafhy_para%max_poros * min(1.0, (soilwater_state%WatSto / soilwater_state%MaxWatSto))
    soilwater_state%Sat  = soilwater_state%Wliq / spafhy_para%max_poros
    soilwater_state%Kh =  spafhy_para%ksat ! [m s-1] REVISE and replace by vanGenuchten -formulation

    call soil_water_retention_curve(soilwater_state%Wliq, psis) ! soil water potential [MPa]
    soilwater_state%Psi = psis

   ! organic top layer; maximum that can be hold is Fc
    soilwater_state%Wliq_top =  spafhy_para%org_fc * min(1.0, (soilwater_state%WatStoTop / soilwater_state%MaxStoTop))
    soilwater_state%beta = min(1.0, soilwater_state%Wliq_top / spafhy_para%org_fc) ! modifier for soil evaporation [-]

   END SUBROUTINE set_soilwaterState
        
  ! real(8) FUNCTION hydrCond(soilwater_state)
  
  ! ! ------- hydrCond
  ! ! returns hydraulic conductivity [ms-1] based on Campbell -formulation
  ! ! assumes drainage from root zone is vertical
  ! ! DEV: for consistency, predict hydr. conductivity dependency on Sat using VanGenuchten -approach
  ! ! -------
  ! ! !ARGUMENTS
  !   type(soilwater_type), intent(in) :: soilwater_state                         ! 

  !   hydrCond = ksat * soilwater_state%Sat**(2*beta + 3.0)
  !   return

  ! END FUNCTION hydrCond


  !real(r8) FUNCTION relative_evaporation(soilwater_state)
  ! ---------------------
  ! returns relative evaporation rate from the organic top layer; loosely
  ! based on Launiainen et al. 2015 Ecol. Mod. Moss-module
  ! Returns:
  !  f - [-], array or grid of 
  ! ---------------------
  !  relative_evaporation = max(0.0, min(0.98*soilwater_state%Wliq_top/soilwater_state%rw_top, 1.0))
  !  return
  !END FUNCTION relative_evaporation
 
  !   SUBROUTINE dry_canopy_et(canopywater_state, snowwater_state, D, Qp, AE, Ta, Ra, Ras, CO2, Rew, P, LAI, fPheno)
!   !
!   !  Computes ET from 2-layer canopy in absense of intercepted precipitiation,
!   !      i.e. in dry-canopy conditions
!   !      IN:
!   !         canopy_type - object
!   !         D  - vpd in kPa
!   !         Qp - PAR in Wm-2
!   !         AE - available energy in Wm-2
!   !         Ta - air temperature degC
!   !         Ra - aerodynamic resistance (s/m)
!   !         Ras - soil aerodynamic resistance (s/m)
!   !         CO2 - atm. CO2 mixing ratio (ppm)
!   !         Rew - relative extractable water [-]
!   !         beta - relative soil conductance for evaporation [-]
!   !         fPheno - phenology modifier [-]
!   !      Args:
!   !         Tr - transpiration rate (mm s-1)
!   !         Efloor - forest floor evaporation rate (mm s-1)
!   !         Gc - canopy conductance (integrated stomatal conductance)  (m s-1)
!   !      SOURCES:
!   !      Launiainen et al. (2016). Do the energy fluxes and surface conductance
!   !      of boreal coniferous forests in Europe scale with leaf area?
!   !      Global Change Biol.
!   !      Modified from: Leuning et al. 2008. A Simple surface conductance model
!   !      to estimate regional evaporation using MODIS leaf area index and the
!   !      Penman-Montheith equation. Water. Resources. Res., 44, W10419
!   !      Original idea Kelliher et al. (1995). Maximum conductances for
!   !      evaporation from global vegetation types. Agric. For. Met 85, 135-147
!   !      Samuli Launiainen, Luke
!   !      Last edit: 13.6.2018: TESTING UPSCALING
!   !

!   !     ---Amax and g1 as LAI -weighted average of conifers and decid.

!     ! !USES:

!     ! use soilwater_mod   ! Placeholder for soil water bucket model, which will provide psi_soil?

!     ! !ARGUMENTS:        
! !    type(canopy_type), intent(in)    :: canopy_type  ! canopy state object
!     real(8)         , intent(in)    :: D            ! vpd in Pa
!     real(8)         , intent(in)    :: Qp           ! PAR in Wm-2
!     real(8)         , intent(in)    :: AE           ! available energy in Wm-2
!     real(8)         , intent(in)    :: Ta           ! air temperature degC
!     real(8)         , intent(in)    :: Ra           ! aerodynamic resistance (s/m)
!     real(8)         , intent(in)    :: Ras          ! soil aerodynamic resistance (s/m)
!     real(8)         , intent(in)    :: CO2          ! atm. CO2 mixing ratio (ppm)
!     real(8)         , intent(in)    :: Rew          ! relative extractable water [-]
!     real(8)         , intent(in)    :: P            ! air pressure
!     real(8)         , intent(in)    :: fPheno       ! phenology modifier [-]
!     real(8)      , intent(in)    :: LAI    ! leaf area index

!     type(canopywater_type), intent(inout)   :: canopywater_state
!     type(snowwater_type), intent(inout)   :: snowwater_state

!     real(8)   :: Gcs, rhoa, tau

!     !real(r8)         , intent(out)   :: Tr           ! transpiration rate (mm s-1), not needed when P-hydro is used
!     !real(r8)         , intent(out)   :: Efloor       ! forest floor evaporation rate (mm s-1), not needed when Efloor is used
!     !real(r8)         , intent(out)   :: Gc           ! canopy conductance (integrated stomatal conductance)(m s-1)

!     !rhoa = 101300.0 / (8.31 * (Ta + 273.15)) ! mol m-3
        
!     !Amax => canopy_type%Amax
!     !g1   => canopy_type%g1
!     !kp   => canopy_type%kp      ! (-) attenuation coefficient for PAR
!     !q50  => canopy_type%q50     ! Wm-2, half-sat. of leaf light response
!     !rw   => canopy_type%rw      ! rew parameter
!     !rwmin => canopy_type%rwmin  ! rew parameter

!     tau = exp(-k * LAI)    ! fraction of Qp at ground relative to canopy top
!                                         ! Hui: repetition of P-hydro model?
    
!     !--- canopy conductance Gc (integrated stomatal conductance)----- 
!     !fQ: Saugier & Katerji, 1991 Agric. For. Met., eq. 4. Leaf light response = Qp / (Qp + q50)
!     !fQ = 1./ kp * log((kp*Qp + q50) / (kp*Qp*exp(-kp * canopy_type%LAI) + q50 + eps) )

!     ! The next formulation is from Leuning et al., 2008 WRR for daily Gc; they refer to 
!     ! Kelliher et al. 1995 AFM but the resulting equation is not exact integral of K95.        
!     ! fQ = 1./ kp * np.log((Qp + q50) / (Qp*np.exp(-kp*self.LAI) + q50))

!     ! soil moisture response: Lagergren & Lindroth, xxxx"""
!     !fRew = min(1.0, max(Rew / rw, rwmin))
!       !fRew = 1.0

!     ! CO2 -response of canopy conductance, derived from APES-simulations
!     ! (Launiainen et al. 2016, Global Change Biology). relative to 380 ppm
!     !fCO2 = 1.0 - 0.387 * log(CO2 / 380.0)
        
!     ! Hui: This part will be delt by P-hydro, can be used as an alternative to compare with P-hydro.
!     ! leaf level light-saturated gs (m/s)
!     !gs = 1.6*(1.0 + g1 / sqrt(D)) * Amax / CO2 / rhoa
        
!     ! canopy conductance
!     !Gc = gs * fQ * fRew * fCO2 * fPheno
    
!     ! Hui: not sure if this is really needed!
!     !Gc[np.isnan(Gc)] = eps

!     ! --- transpiration rate ---
!     !Tr= penman_monteith((1.-tau)*AE, 1e3*D, Ta, Gc, 1./Ra, units='mm')) 

!     ! --- forest floor evaporation rate---
!     ! soil conductance is function of relative water availability
!     ! gcs = 1. / self.soilrp * beta**2.0
!     ! beta = Wliq / FC; Best et al., 2011 Geosci. Model. Dev. JULES
!     Gcs = gsoil
        
!     canopywater_state%Efloor = beta * penman_monteith(tau * AE, D, Ta, Gcs, 1.0/Ras, P, units='mm')

!     if (snowwater_state%swe>0) then
!       canopywater_state%Efloor = 0.0  ! no evaporation from floor if snow on ground or beta == 0
!     end if

!   END SUBROUTINE dry_canopy_et

 ! SUBROUTINE topmodel(soilwater_state, R)
  ! -----------------------------
  ! runs a timestep, updates saturation deficit and returns fluxes
  !      Args:
  !          R - recharge [m per unit catchment area] during timestep
  !      OUT:
  !          Qb - baseflow [m per unit area]
  !          Qr - returnflow [m per unit area]
  !          qr - distributed returnflow [m]
  !          fsat - saturated area fraction [-]
  !      Note: 
  !          R is the mean drainage [m] from bucketgrid.
  ! ------------------------------

  ! !ARGUMENTS:
  !  type(soilwater_type), intent(inout)    :: soilwater_state    
  !  real(8)      , intent(in)    :: R       ! recharge [m per unit catchment area] during timestep
  !  real(8)      , intent(out)   :: Qb      ! baseflow [m per unit area]
  !  real(8)      , intent(out)   :: Qr      ! returnflow [m per unit area]
  !  real(8)      , intent(out)   :: qr      ! distributed returnflow [m]
  !  real(8)      , intent(out)   :: fsat    ! saturated area fraction [-]

  ! ! Local variables 
  !  real     ::    S0, s, Qb, S

  ! initial conditions
  !  So = soilwater_state%S
  !  s  = soilwater_state%local_s(So)

  ! subsurface flow, based on initial state
  !  Qb = soilwater_state%subsurfaceflow()

  ! update storage deficit and check where we have returnflow 
  !  S = So + Qb - R
  !  s = soilwater_state%local_s(S)

  ! returnflow grid
  !  qr = -s
        
  !  if (qr < 0) then
  !    qr = 0.0            ! returnflow grid, m
  !  end if
        
  ! average returnflow per unit area
  !  Qr = qr*soilwater_state%CellArea/soilwater_state%CatchmentArea

  ! now all saturation excess is in Qr so update s and S. 
  ! Deficit increases when Qr is removed 
  !  S = S + Qr
  !  soilwater_state%S = S

  ! saturated area fraction
  !  if (s <= 0) then
    ! fsat = np.max(np.shape(ix))*self.CellArea/self.CatchmentArea
  !    fsat = s*soilwater_state%CellArea/soilwater_state%CatchmentArea
  !  end if

    ! check mass balance
  !  dS = (So - soilwater_state%S)
  !  dF = R - Qb - Qr
  !  mbe = dS - dF

  ! append results
   ! if hasattr(self, 'results'):
   !       self.results['R'].append(R)
   !       self.results['S'].append(self.S)
   !       self.results['Qb'].append(Qb)
   !       self.results['Qr'].append(Qr)
   !       self.results['qr'].append(qr)
   !       self.results['fsat'].append(fsat)
   !       self.results['Mbe'].append(mbe)

   !     return Qb, Qr, qr, fsat
  
  !END SUBROUTINE topmodel

end module spafhy_mod
