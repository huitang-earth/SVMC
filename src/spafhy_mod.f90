MODULE spafhy_mod

!--------------------------------------------------------------
! References:
! Launiainen, S., Guan, M., Salmivaara, A., Kieloaho, A.-J., 2019. 
! Modeling boreal forest evapotranspiration and water balance at stand and catchment scales: a spatial approach. 
! Hydrology and Earth System Sciences 23, 3457–3480. https://doi.org/10.5194/hess-23-3457-2019
!--------------------------------------------------------------

! Modules
  use netcdf             ! library for processing netcdf files
  use readpara_mod       ! module for reading parameter files in ASCII
  use readclim_mod       ! module for reading reading meteorological forcing data
  use readsoil_mod       ! module for reading soil properties (shared with yasso?)

  implicit none

  !Public member functions:
  public :: canopy_water_flux   ! require inputdata for p-hydro module
  public :: topmodel   !
  public :: soil_water
  public :: soil_water_retention_curve

  private :: aerodynamics        ! aerodynamic conductances
  private :: canopy_water_snow   ! interception, evaporation and snowpack
  private :: dry_canopy_et       ! Computes ET from 2-layer canopy in absense of intercepted precipitiation
                                 ! This will be replaced by p-hydro when p-hydro is turned on
  private :: penman_monteith     !
  private :: e_sat               !
  private :: set_soilwaterState   !
  private :: hydrCond             ! 
  private :: relative_evaporation !

contains
  
  SUBROUTINE topmodel(soilwater_state, R)
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
    type(soilwater_state), intent(inout)    :: soilwater_state    
    real(r8)      , intent(in)    :: R       ! recharge [m per unit catchment area] during timestep
    real(r8)      , intent(out)   :: Qb      ! baseflow [m per unit area]
    real(r8)      , intent(out)   :: Qr      ! returnflow [m per unit area]
    real(r8)      , intent(out)   :: qr      ! distributed returnflow [m]
    real(r8)      , intent(out)   :: fsat    ! saturated area fraction [-]

  ! ! Local variables 
    real     ::    S0, s, Qb, S

  ! initial conditions
    So = soilwater_state%S
    s  = soilwater_state%local_s(So)

  ! subsurface flow, based on initial state
    Qb = soilwater_state%subsurfaceflow()

  ! update storage deficit and check where we have returnflow 
    S = So + Qb - R
    s = soilwater_state%local_s(S)

  ! returnflow grid
    qr = -s
        
    if (qr < 0) then
      qr = 0.0            ! returnflow grid, m
    end if
        
  ! average returnflow per unit area
    Qr = qr*soilwater_state%CellArea/soilwater_state%CatchmentArea

  ! now all saturation excess is in Qr so update s and S. 
  ! Deficit increases when Qr is removed 
    S = S + Qr
    soilwater_state%S = S

  ! saturated area fraction
    if (s <= 0) then
    ! fsat = np.max(np.shape(ix))*self.CellArea/self.CatchmentArea
      fsat = s*soilwater_state%CellArea/soilwater_state%CatchmentArea
    end if

    ! check mass balance
    dS = (So - soilwater_state%S)
    dF = R - Qb - Qr
    mbe = dS - dF

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
  
  END SUBROUTINE topmodel

  SUBROUTINE soil_water(soilwater_state, dt, rr, tr, evap, retflow)
  ! ---------------------------------------------------------------
  !      Computes 2-layer bucket model water balance for one timestep dt
  !      Top layer is interception storage and contributes only to evap.
  !      Lower layer is rootzone and contributes only tr and creates drainage.
  !      Capillary interaction between layers is neglected and connection from bottom up
  !      is only in case of excess returnflow.
  !      Pond storage can exist above top layer.   
  !      IN:
  !          dt [s]
  !          rr = potential infiltration [m]
  !          tr = transpiration from root zone [m]
  !          evap = evaporation from top layer [m]
  !          retflow = return flow from ground water [m]
  !      OUT: dict with 
  !          inflow [m] - total inflow to root zone
  !          roff [m] - surface runoff
  !          drain [m] - drainage from root zone
  !          tr [m] - transpiration from root zone
  !          mbe [m] - mass balance error
  ! ----------------------------------------------------------------

  ! !ARGUMENTS:
    type(soilwater_state), intent(inout)    :: soilwater_state    
    real(r8)      , intent(in)    :: dt      ! [s]
    real(r8)      , intent(in)    :: rr      ! potential infiltration [m]
    real(r8)      , intent(inout) :: tr      ! transpiration from root zone [m]
    real(r8)      , intent(in)    :: evap    ! evaporation from top layer [m]
    real(r8)      , intent(in)    :: retflow ! return flow from ground water [m]

    real(r8)      , intent(out)   :: inflow  ! [m] - total inflow to root zone
    real(r8)      , intent(out)   :: roff    ! [m] - surface runoff
    real(r8)      , intent(out)   :: drain   ! [m] - drainage from root zone
    real(r8)      , intent(out)   :: tr      ! [m] - transpiration from root zone
    real(r8)      , intent(out)   :: mbe     ! [m] - mass balance error

  ! ! Local variables 
    real     ::    rr0, PondSto0, WatSto0, interc, eps=???
    real     ::    Qin, dSto, exfil, to_pond, to_top_layer
        
    rr0 = rr
       
    ! add current Pond storage to rr & update storage
    PondSto0 = soilwater_state%PondSto
    rr = rr+soilwater_state%PondSto
    soilwater_state%PondSto = 0.0

    WatSto0 = soilwater_state%WatSto
    WatStoTop0 = soilwater_state%WatStoTop
        
    !top layer interception & water balance
    interc = max(0.0, (soilwater_state%MaxStoTop - soilwater_state%WatStoTop)) &
                    * (1.0 - exp(-(rr / soilwater_state%MaxStoTop)))

    soilwater_state%WatStoTop = max(0.0, soilwater_state%WatStoTop + interc)  
    evap = min(evap, soilwater_state%WatStoTop)
    soilwater_state%WatStoTop = soilwater_state%WatStoTop - evap
      
    ! infiltration to rootzone
    rr = rr - interc
                
    ! ********* compute bottom layer (root zone) water balance ***********
    ! transpiration removes water from rootzone
    tr = min(tr, soilwater_state%WatSto - eps)
    soilwater_state%WatSto = soilwater_state%WatSto - tr
        
    ! drainage: at gridcells where retflow > 0, set drain to zero.
    ! This delays drying of cells which receive water from topmodel storage
    ! ... and removes oscillation of water content at those cells.
    if (retflow .gt. 0.0) then 
      drain = 0.0
    else
      drain = min(soilwater_state%hydrCond * dt, max(0.0, (soilwater_state%Wliq - soilwater_state%Fc))*soilwater_state%D)
    end if 
    
    ! inflow to root zone: restricted by potential inflow or available pore space
    Qin = (retflow + rr)          ! m, pot. inflow
    inflow = min(Qin, soilwater_state%MaxWatSto - soilwater_state%WatSto + drain)    
    dSto = (inflow - drain)
    soilwater_state%WatSto = min(soilwater_state%MaxWatSto, max(soilwater_state%WatSto + dSto, eps))
                
    ! if inflow excess after filling rootzone, update first top layer storage
    exfil = Qin - inflow
    to_top_layer = min(exfil, soilwater_state%MaxStoTop - soilwater_state%WatStoTop - eps)
    ! self.WatStoTop = self.WatStoTop + to_top_layer
    soilwater_state%WatStoTop = soilwater_state%WatStoTop + to_top_layer
        
    ! ... and then pond storage ...
    to_pond = min(exfil - to_top_layer, soilwater_state%MaxPond - soilwater_state%PondSto - eps)
    soilwater_state%PondSto = soilwater_state%PondSto + to_pond
 
    ! ... and route remaining to surface runoff
    roff = exfil - to_top_layer - to_pond

    ! compute diagnostic state variables at root zone:
    call set_soilwaterState(soilwater_state)

    ! mass balance error [m]
    mbe = (soilwater_state%WatSto - WatSto0)  + (soilwater_state%WatStoTop - WatStoTop0) + (soilwater_state%PondSto - PondSto0) &
            - (rr0 + retflow - tr - evap - drain - roff)

    ! update grid total drainage to ground water [m]
    !    self._drainage_to_gw = np.nansum(drain)
           
    ! append results to lists; use only for testing small grids!
    !    if hasattr(self, 'results'):
    !        self.results['Infil'].append(inflow - retflow)   # infiltration through top boundary
    !        self.results['Retflow'].append(retflow)         # return flow from below 
    !        self.results['Roff'].append(roff)       # surface runoff
    !        self.results['Drain'].append(drain)     # drainage
    !        self.results['ET'].append(tr + evap)            
    !        self.results['Mbe'].append(mbe)
    !        self.results['Wliq'].append(self.Wliq)
    !        self.results['PondSto'].append(self.PondSto)
    !        self.results['Wliq_top'].append(self.Wliq_top)
    !        self.results['Ree'].append(self.Ree)
        
    !    return inflow, roff, drain, tr, evap, mbe
  
  END SUBROUTINE soil_water

  SUBROUTINE soil_water_retention_curve(sand, clay, om_frac, soil_depth, soilwater_state, smp, dsmpds)
  !
  ! !DESCRIPTION:
  ! Two steps: 
  !    1- compute hydraulic properties based on functions derived 
  !       from Table 5 in cosby et al, 1984 (Based on CTSM: )
  !    2- Compute soil suction potential from soil water storage
  !       Implementation of soil_water_retention_curve_type using different approaches
  !        -  Clapp-Hornberg 1978 parameterizations (Based on CTSM: )
  !        -  Van Genuchten 1980 parameterizations (Based on CTSM-FATES:)
  !             -  Use parameters derived from Launiainen et al. 2022 
  ! !ARGUMENTS:

    implicit none
    real(r8), intent(in) :: sand        !% sand
    real(r8), intent(in) :: clay        !% clay
    real(r8), intent(in) :: om_frac     !%, organic matter fraction
    real(r8), intent(in) :: soil_depth  !m, depth of soil bucket           
    type(soilwater_state_type), intent(in) :: soilwater_state
    real(r8), intent(out)            :: smp      !soil suction, negative, [mm]
    !real(r8), intent(out)            :: hk      !hydraulic conductivity [mm/s]
    real(r8), optional, intent(out)  :: dsmpds   !d[smp]/ds, [mm]
    !real(r8), optional, intent(out)  :: dhkds    !d[hk]/ds   [mm/s]
  !
  ! !LOCAL VARIABLES:
    real(r8) :: bsw         ! b shape parameter 
    real(r8) :: min_bsw     ! b shape parameter for mineral soil
    real(r8) :: om_bsw      ! b shape parameter for organic matter
    real(r8) :: watsat      ! v/v saturate moisture
    real(r8) :: min_watsat  ! v/v saturate moisture for mineral soil 
    real(r8) :: om_watsat   ! v/v saturate moisture for organic matter 
    real(r8) :: sucsat  ! mm, soil matric potential at saturation point
    real(r8) :: min_sucsat  ! mm, soil matric potential at saturation point, mineral soil
    real(r8) :: om_sucsat   ! mm, soil matric potential at saturation point, organic soil
    real(r8) :: hksat       ! mm/s, saturated hydraulic conductivity
    real(r8) :: min_xksat   ! mm/s, saturated hydraulic conductivity, mineral soil
    real(r8) :: om_xksat    ! mm/s, saturated hydraulic conductivity, organic soil
    real(r8) :: vol_ice     ! v/v, volumetric ice in soil bucket  
    real(r8) :: vol_liq     ! v/v, volumetric of liq in soil bucket
    real(r8) :: alpha       ! parameter for Van Genuchten
    real(r8) :: watres      ! v/v, residual soil moisture for Van Genuchten
    real(r8) :: n           ! (-), pore-size-distribution parameter for Van Genuchten
    real(r8) :: satfrac     ! parameter for Van Genuchten
    real(r8) :: eff_porosity! v/v, volume of ice
 
  ! Cosby et al. Table 5     
  !  min_watsat = 0.489_r8 - 0.00126_r8*sand
  !  min_bsw    = 2.91 + 0.159*clay
  !  min_sucsat = 10._r8 * ( 10._r8**(1.88_r8-0.0131_r8*sand) )            
    !min_xksat  = 0.0070556 *( 10.**(-0.884+0.0153*sand) ) ! mm/s, from table 5 

  !  om_watsat  = 0.93_r8 
  !  om_bsw     = 2.7_r8
  !  om_sucsat  = 10.3_r8    
    ! om_xksat = 0.28_r8

  !  watsat     = (1._r8 - om_frac) * min_watsat + om_watsat*om_frac
  !  bsw        = (1._r8 - om_frac) * min_bsw + om_bsw*om_frac  
  !  sucsat     = (1._r8 - om_frac) * min_sucsat + om_sucsat*om_frac  
    ! hksat    = ?

    ! Hui: there is no soil ice in spaphy???? 
    !vol_ice = min(watsat, h2osoi_ice_col/(dz(c,j)*denice))
  !  vol_ice = 0.0_r8   
  !  eff_porosity = max(0.01_r8, watsat-vol_ice)

    !soil_depth=soilwater_state%MaxWatSto/denh2o/watsat
    !soil depth is parameterized in spaphy (40 cm by default)
  !  vol_liq = min(eff_porosity, soilwater_state%WatSto/(soil_depth*denh2o))

    ! relative saturation, [0, 1]
  !  s = max(vol_liq/eff_porosity,0.01_r8)
      
    !compute soil suction potential, negative
  !  smp = -sucsat*s**(-bsw)
    !hk=imped*hksat*s**(2._r8*bsw+3._r8)  ! Hui: This is similar to Spaphy

    !compute derivative
  !  if(present(dsmpds))then
  !     dsmpds=-bsw*smp/s
       !dhkds=(2._r8*bsw+3._r8)*hk/s
  !  endif

  !--- Van Genuchten scheme based on FatesHydro (parameters from Launiainen et al. Forests, 2022)
    
    n=1.07       !Launiainen et al. 2022: C1-5: 1.12, 1.14, 1.07, 1.27, 1.18  
    m=1._r8/n  
    watres=0.0   !Launiainen et al. 2022: C1-5: 0.0
    alpha=2.02   !Launiainen et al. 2022: C1-5: 4.45, 5.92, 2.02, 4.49, 3.35
    watsat=0.46  !Launiainen et al. 2022: C1-5: 0.75, 0.68, 0.46, 0.47, 0.54  

    vol_ice = 0.0_r8   
    eff_porosity = max(0.01_r8, watsat-vol_ice)
    
    satfrac = (vol_liq-watres)/(eff_porosity-watres)
    smp = -(1._r8/alpha)*(satfrac**(1._r8/(m-1._r8)) - 1._r8 )**m

  END SUBROUTINE soil_water_retention_curve


  SUBROUTINE canopy_water_flux(gs, Rn, Ta, Prec, Rg, Par, VPD, U=2.0, CO2=380.0, Rew=1.0, beta=1.0, P=101300.0)
    !
    ! !DESCRIPTION:
    ! Adapted from Gridded canopy and snow hydrology model for SpaFHy.
    ! Based on simple schemes for computing canopy water transpiration and evaporation
    ! Can be both daily or sub-daily timesteps.
    !
    ! !USES:

    ! use soilwater_mod   ! Placeholder for soil water bucket model, which will provide psi_soil?

    ! !ARGUMENTS:
    real(r8)      , intent(in)    :: gs     ! Stomatal conductance (mol C m-2 Pa-1)
    real(r8)      , intent(in)    :: Rn     ! Net solar radiation obsorbed by canopy (W m-2)
    real(r8)      , intent(in)    :: Ta     ! Air temperature (tc), (degrees C)
    real(r8)      , intent(in)    :: Par    ! photos. act. radiation [Wm-2]
    real(r8)      , intent(in)    :: VPD    ! Vapour pressure deficit (Pa) (will be calculated using pressure & humidity)
    real(r8)      , intent(in)    :: CO2    ! Atmospheric CO2 concentration (ppm)
    real(r8)      , intent(in)    :: prec   ! precipitatation rate [mm/s]
    real(r8)      , intent(in)    :: U      ! mean wind speed at ref. height above canopy top [ms-1]
    real(r8)      , intent(in)    :: P      ! pressure [Pa], scalar or matrix
    real(r8)      , intent(in)    :: beta   ! term for soil evaporation resistance (Wliq/FC) [-]

    real(r8)      , intent(out)   :: PotInf     ! potential infiltration to soil profile (mm)
    real(r8)      , intent(out)   :: Trfall     ! throughfall to snow / soil surface (mm)
    real(r8)      , intent(out)   :: Interc     ! interception of canopy (mm)
    real(r8)      , intent(out)   :: Evap       ! evaporation / sublimation from canopy store (mm)
    real(r8)      , intent(out)   :: ET         ! total evapo-transpiration
    real(r8)      , intent(out)   :: Transpi    ! transpiration rate (mm s-1)
    real(r8)      , intent(out)   :: Efloor     ! forest floor evaporation rate (mm s-1)
    real(r8)      , intent(out)   :: MBE        ! mass balance error (mm)

    ! ! Local variables 

    real     ::    fn_profit
    real     ::    lj_dps
    real     ::    aj
    real     ::    
    real     ::    
    real     ::   
    integer  :: i, j, k, t
  

    ! Calculate aerodynamic conductances
    Call aerodynamics(LAI, hc, U, w=0.01, zm=zmeas, zg=zground, zos=zo_ground, Ra, Rb, Ras, ustar, Uh, Ug)

    ! Calculate interception, evaporation and snowpack
    call canopy_water_snow(dt, Ta, Prec, Rn, VPD, Ra=Ra, PotInf, Trfall, Evap, Interc, MBE, unload)

    ! Dry-canopy evapotranspiration [mm s-1]
    call dry_canopy_et(VPD, Par, Rn, Ta, Ra=Ra, Ras=Ras, CO2=CO2, Rew=Rew, beta=beta, fPheno=self.fPheno, Transpi, Efloor, Gc)
        
    Transpi = Transpi * dt
    Efloor = Efloor * dt
    ET = Transpi + Efloor + Evap

  END SUBROUTINE canopy_water_flux

  SUBROUTINE canopy_water_snow(canopy_type, dt, T, Prec, AE, D, Ra, U, snowtype)
    !
    ! Calculates canopy water interception and SWE during timestep dt
    ! !USES:

    ! use soilwater_mod   ! Placeholder for soil water bucket model, which will provide psi_soil?

    ! !ARGUMENTS:
    real(r8)      , intent(in)    :: dt     ! timestep [s]
    real(r8)      , intent(in)    :: T      ! air temperature (degC)
    real(r8)      , intent(in)    :: Prec   ! precipitation rate during (mm d-1)
    real(r8)      , intent(in)    :: AE     ! available energy (~net radiation) (Wm-2) CO2    ! Atmospheric CO2 concentration (ppm)
    real(r8)      , intent(in)    :: D      ! vapor pressure deficit (kPa)
    real(r8)      , intent(in)    :: U      ! mean wind speed at ref. height above canopy top [ms-1]
    real(r8)      , intent(in)    :: Ra     ! canopy aerodynamic resistance (s m-1)

    real(r8)      , intent(out)   :: PotInf     ! potential infiltration to soil profile (mm)
    real(r8)      , intent(out)   :: Trfall     ! throughfall to snow / soil surface (mm)
    real(r8)      , intent(out)   :: Interc     ! interception of canopy (mm)
    real(r8)      , intent(out)   :: Evap       ! evaporation / sublimation from canopy store (mm)
    real(r8)      , intent(out)   :: Unload     ! undloading from canopy storage (mm)
    real(r8)      , intent(out)   :: MBE        ! mass balance error (mm)       
    !self - updated state W, Wf, SWE, SWEi, SWEl
    
    ! Local variables.
    real(r8)  :: fW, fS, Tmin, Tmax, Tmelt, Wmax, Wmaxsnow, Kmelt, Kfreeze


    ! quality of precipitation
    Tmin = 0.0  ! 'C, below all is snow
    Tmax = 1.0  ! 'C, above all is water
    Tmelt = 0.0  ! 'C, T when melting starts

    !storage capacities mm
    Wmax = wmax * canopytype(?)%LAI
    Wmaxsnow = wmaxsnow * canopytype(?)%LAI

    ! melting/freezing coefficients mm/s
    Kmelt = Kmelt - 1.64 * cf / dt  ! Kuusisto E, 'Lumi Suomessa'
    Kfreeze = Kfreeze

    tau = exp(-kp*canopytype(?)%LAI)  ! fraction of Rn at ground

    ! inputs to arrays, needed for indexing later in the code
    !   gridshape = np.shape(self.LAI)  # rows, cols
    !    if np.shape(T) != gridshape:
    !        T = np.ones(gridshape) * T
    !       Prec = np.ones(gridshape) * Prec
    !        AE = np.ones(gridshape) * AE
    !        D = np.ones(gridshape) * D
    !        Ra = np.ones(gridshape) * Ra

    Prec = Prec * dt  ! mm
    
    ! latent heat of vaporization (Lv) and sublimation (Ls) J kg-1
    Lv = 1e3 * (3147.5 - 2.37 * (T + 273.15))
    Ls = Lv + 3.3e5

    ! compute 'potential' evaporation / sublimation rates for each grid cell
    Ga = 1. / Ra ! aerodynamic conductance
    Ce = 0.01*((canopy_type%W + eps) / Wmaxsnow)**(-0.4)  ! exposure coeff (-)
    Sh = (1.79 + 3.0*U**0.5)                      ! Sherwood numbner (-)
    gi = Sh*canopy_type%W*Ce / 7.68 + eps                ! m s-1

    ! resistance for snow sublimation adopted from:
    ! Pomeroy et al. 1998 Hydrol proc; Essery et al. 2003 J. Climate;
    ! Best et al. 2011 Geosci. Mod. Dev.

    if ((Prec == 0) .and. (T <= Tmin)) then
      erate = dt / Ls * penman_monteith((1.0 - tau)*AE, 1e3*D, T, gi, Ga, units='W')
    else if ((Prec == 0) .and. (T > Tmin)) then
      ! evaporation of intercepted water, mm
      gs = 1e6
      erate = dt / Lv * penman_monteith((1.0 - tau)*AE,1e3*D, T, gs, Ga, units='W')
    end if 

    ! state of precipitation [as water (fW) or as snow(fS)]
    if ( T <= Tmin)) then
      fS = 1.0
    else if (T >= Tmax) then
      fW = 1.0
    else if ((T > Tmin) .and. (T < Tmax)) then
      fW = (T - Tmin) / (Tmax - Tmin)
      fS = 1.0 - fW
    end if

    !----- Initial conditions for calculating mass balance error
    Wo = canopy_type%W      ! canopy storage
    SWEo = canopy_type%SWE  ! Snow water equivalent mm

    !----- Canopy water storage change
    ! snow unloading from canopy, ensures also that seasonal LAI development does
    ! not mess up computations
    ! HT: (1) unloading first at each timestep? (2) why no snow unload? do we need it here? 
    if (T >= Tmax) then
      Unload = max(canopy_type%W - Wmax, 0.0)
      canopy_type%W = canopy_type%W - Unload
      ! dW = self.W - Wo
    end if

    !----- Interception of rain or snow: asymptotic approach of saturation.
    !      based on: Hedstrom & Pomeroy 1998. Hydrol. Proc 12, 1611-1625;
    !                Koivusalo & Kokkonen 2002 J.Hydrol. 262, 145-164.
    if (T < Tmin) then
      Interc = (Wmaxsnow - canopy_type%W) &
                * (1.0 - exp(-(canopy_type%cf / Wmaxsnow) * Prec))
    end if
        
    ! Above Tmin, interception capacity equals that of liquid precip
    if (T >= Tmin) then
      Interc = max(0.0, (Wmax - canopy_type%W)) &
                * (1.0 - exp(-(canopy_type%cf/Wmax) * Prec))
    end if

    canopy_type%W = canopy_type%W + Interc  ! new canopy storage, mm
    Trfall = Prec + Unload - Interc  ! Throughfall to field layer or snowpack

    ! evaporate from canopy and update storage
    Evap = min(erate, canopy_type%W)  ! mm
    canopy_type%W = canopy_type%W - Evap

    !---- Snowpack (in case no snow, all Trfall routed to floor) """
    if (T >= Tmelt) then
      Melt = min(canopy_type%SWEi, Kmelt * dt * (T - Tmelt))  ! mm
    else if (T < Tmelt) then
      Freeze = min(canopy_type%SWEl, Kfreeze * dt * (Tmelt - T))  ! mm
    end if

    !---- amount of water as ice and liquid in snowpack
    Sice = max(0.0, canopy_type%SWEi + fS * Trfall + Freeze - Melt)
    Sliq = max(0.0, canopy_type%SWEl + fW * Trfall - Freeze + Melt)

    PotInf = max(0.0, Sliq - Sice * canopy_type%R)  ! mm
    Sliq   = max(0.0, Sliq - PotInf)  ! mm, liquid water in snow

    ! update Snowpack state variables
    canopy_type%SWEl = Sliq
    canopy_type%SWEi = Sice
    canopy_type%SWE  = canopy_type%SWEl + canopy_type%SWEi
        
    ! mass-balance error mm
    MBE = (canopy_type%W + canopy_type%SWE) - (Wo + SWEo) - (Prec - Evap - PotInf)

  END SUBROUTINE canopy_water_snow


  SUBROUTINE aerodynamics(LAI, hc, Uo, w=0.01, zm=2.0, zg=0.5, zos=0.01, ra, rb, ras, ustar, Uh, Ug)
    !
    ! computes wind speed at ground and canopy + boundary layer conductances
    ! Computes wind speed at ground height assuming logarithmic profile above and
    ! exponential within canopy
    ! SOURCE:
    !   Cammalleri et al. 2010 Hydrol. Earth Syst. Sci
    !   Massman 1987, BLM 40, 179 - 197.
    !   Magnani et al. 1998 Plant Cell Env.
    
    ! !USES:

    ! use soilwater_mod   ! Placeholder for soil water bucket model, which will provide psi_soil?

    ! !ARGUMENTS:
    real(r8)      , intent(in)    :: LAI    ! one-sided leaf-area /plant area index (m2m-2)
    real(r8)      , intent(in)    :: hc     ! canopy height (m)
    real(r8)      , intent(in)    :: Uo     ! mean wind speed at height zm (ms-1)
    real(r8)      , intent(in)    :: w      ! leaf length scale (m)
    real(r8)      , intent(in)    :: zm     ! wind speed measurement height above canopy (m)
    real(r8)      , intent(in)    :: zg     ! height above ground where Ug is computed (m)
    real(r8)      , intent(in)    :: zos    ! forest floor roughness length, ~ 0.1*roughness element height (m)

    real(r8)      , intent(out)   :: ra         ! canopy aerodynamic resistance (sm-1)
    real(r8)      , intent(out)   :: rb         ! canopy boundary layer resistance (sm-1)
    real(r8)      , intent(out)   :: ras        ! forest floor aerod. resistance (sm-1)
    real(r8)      , intent(out)   :: ustar      ! friction velocity (ms-1)
    real(r8)      , intent(out)   :: Uh         ! wind speed at hc (ms-1)
    real(r8)      , intent(out)   :: Ug         ! wind speed at zg (ms-1)

    ! local variables:
    real(r8)                      :: beta = 285.0   ! s/m, from Campbell & Norman eq. (7.33) x 42.0 molm-3
    kv = 0.4  ! von Karman constant (-)


    zm = hc + zm  ! m
    zg = min(zg, 0.1 * hc)
    alpha = LAI / 2.0  ! wind attenuation coeff (Yi, 2008 eq. 23)
    d = 0.66*hc     !? These coefficients need to be moved readpara_mod.f90?
    zom = 0.123*hc  ! 
    zov = 0.1*zom
    zosv = 0.1*zos

    ! solve ustar and U(hc) from log-profile above canopy
    ustar = Uo * kv / log((zm - d) / zom) 
    Uh = ustar / kv * log((hc - d) / zom)
    
    ! U(zg) from exponential wind profile
    zn = min(zg / hc, 1.0)  ! zground can't be above canopy top
    Ug = Uh * exp(alpha*(zn - 1.0))

    ! canopy aerodynamic & boundary-layer resistances (sm-1). Magnani et al. 1998 PCE eq. B1 & B5
    !ra = 1. / (kv*ustar) * log((zm - d) / zom)
    ra = 1./(kv**2.0 * Uo) * log((zm - d) / zom) * log((zm - d) / zov)    
    rb = 1. / LAI * beta * ((w / Uh)*(alpha / (1.0 - exp(-alpha / 2.0))))**0.5

    ! soil aerodynamic resistance (sm-1)
    ras = 1. / (kv**2.0*Ug) * (log(zg / zos))*log(zg / (zosv))
    
    !print('ra', ra, 'rb', rb)
    ra = ra + rb

  END SUBROUTINE aerodynamics

  SUBROUTINE dry_canopy_et(canopy_type, D, Qp, AE, Ta, Ra=25.0, Ras=250.0, CO2=380.0, Rew=1.0, beta=1.0, fPheno=1.0)
  !
  !  Computes ET from 2-layer canopy in absense of intercepted precipitiation,
  !      i.e. in dry-canopy conditions
  !      IN:
  !         canopy_type - object
  !         D  - vpd in kPa
  !         Qp - PAR in Wm-2
  !         AE - available energy in Wm-2
  !         Ta - air temperature degC
  !         Ra - aerodynamic resistance (s/m)
  !         Ras - soil aerodynamic resistance (s/m)
  !         CO2 - atm. CO2 mixing ratio (ppm)
  !         Rew - relative extractable water [-]
  !         beta - relative soil conductance for evaporation [-]
  !         fPheno - phenology modifier [-]
  !      Args:
  !         Tr - transpiration rate (mm s-1)
  !         Efloor - forest floor evaporation rate (mm s-1)
  !         Gc - canopy conductance (integrated stomatal conductance)  (m s-1)
  !      SOURCES:
  !      Launiainen et al. (2016). Do the energy fluxes and surface conductance
  !      of boreal coniferous forests in Europe scale with leaf area?
  !      Global Change Biol.
  !      Modified from: Leuning et al. 2008. A Simple surface conductance model
  !      to estimate regional evaporation using MODIS leaf area index and the
  !      Penman-Montheith equation. Water. Resources. Res., 44, W10419
  !      Original idea Kelliher et al. (1995). Maximum conductances for
  !      evaporation from global vegetation types. Agric. For. Met 85, 135-147
  !      Samuli Launiainen, Luke
  !      Last edit: 13.6.2018: TESTING UPSCALING
  !

  !     ---Amax and g1 as LAI -weighted average of conifers and decid.

    ! !USES:

    ! use soilwater_mod   ! Placeholder for soil water bucket model, which will provide psi_soil?

    ! !ARGUMENTS:        
    type(canopy_type), intent(in)    :: canopy_type  ! canopy state object
    real(r8)         , intent(in)    :: D            ! vpd in kPa
    real(r8)         , intent(in)    :: Qp           ! PAR in Wm-2
    real(r8)         , intent(in)    :: AE           ! available energy in Wm-2
    real(r8)         , intent(in)    :: Ta           ! air temperature degC
    real(r8)         , intent(in)    :: Ra           ! aerodynamic resistance (s/m)
    real(r8)         , intent(in)    :: Ras          ! soil aerodynamic resistance (s/m)
    real(r8)         , intent(in)    :: CO2          ! atm. CO2 mixing ratio (ppm)
    real(r8)         , intent(in)    :: Rew          ! relative extractable water [-]
    real(r8)         , intent(in)    :: beta         ! relative soil conductance for evaporation [-]
    real(r8)         , intent(in)    :: fPheno       ! phenology modifier [-]

    real(r8)         , intent(out)   :: Tr           ! transpiration rate (mm s-1), not needed when P-hydro is used
    real(r8)         , intent(out)   :: Efloor       ! forest floor evaporation rate (mm s-1), not needed when Efloor is used
    real(r8)         , intent(out)   :: Gc           ! canopy conductance (integrated stomatal conductance)(m s-1)

    rhoa = 101300.0 / (8.31 * (Ta + 273.15)) ! mol m-3
        
    Amax => canopy_type%Amax
    g1   => canopy_type%g1
    kp   => canopy_type%kp      ! (-) attenuation coefficient for PAR
    q50  => canopy_type%q50     ! Wm-2, half-sat. of leaf light response
    rw   => canopy_type%rw      ! rew parameter
    rwmin => canopy_type%rwmin  ! rew parameter

    tau = exp(-kp * canopy_type%LAI)    ! fraction of Qp at ground relative to canopy top
                                        ! Hui: repetition of P-hydro model?
    
    !--- canopy conductance Gc (integrated stomatal conductance)----- 
    !fQ: Saugier & Katerji, 1991 Agric. For. Met., eq. 4. Leaf light response = Qp / (Qp + q50)
    fQ = 1./ kp * log((kp*Qp + q50) / (kp*Qp*exp(-kp * canopy_type%LAI) + q50 + eps) )

    ! The next formulation is from Leuning et al., 2008 WRR for daily Gc; they refer to 
    ! Kelliher et al. 1995 AFM but the resulting equation is not exact integral of K95.        
    ! fQ = 1./ kp * np.log((Qp + q50) / (Qp*np.exp(-kp*self.LAI) + q50))

    ! soil moisture response: Lagergren & Lindroth, xxxx"""
    fRew = min(1.0, max(Rew / rw, rwmin))
      !fRew = 1.0

    ! CO2 -response of canopy conductance, derived from APES-simulations
    ! (Launiainen et al. 2016, Global Change Biology). relative to 380 ppm
    fCO2 = 1.0 - 0.387 * log(CO2 / 380.0)
        
    ! Hui: This part will be delt by P-hydro, can be used as an alternative to compare with P-hydro.
    ! leaf level light-saturated gs (m/s)
    gs = 1.6*(1.0 + g1 / sqrt(D)) * Amax / CO2 / rhoa
        
    ! canopy conductance
    Gc = gs * fQ * fRew * fCO2 * fPheno
    
    ! Hui: not sure if this is really needed!
    !Gc[np.isnan(Gc)] = eps

    ! --- transpiration rate ---
    Tr= penman_monteith((1.-tau)*AE, 1e3*D, Ta, Gc, 1./Ra, units='mm')) 

    ! --- forest floor evaporation rate---
    ! soil conductance is function of relative water availability
    ! gcs = 1. / self.soilrp * beta**2.0
    ! beta = Wliq / FC; Best et al., 2011 Geosci. Model. Dev. JULES
    Gcs = canopy_type%gsoil
        
    Efloor = beta * penman_monteith(tau * AE, 1e3*D, Ta, Gcs, 1./Ras, P, units='mm')

    if (canopy_state%SWE>0) then
      Efloor = 0.0  ! no evaporation from floor if snow on ground or beta == 0
    end if

    return Tr, Efloor, Gc
  END SUBROUTINE dry_canopy_et

  real(r8) FUNCTION penman_monteith(AE, D, T, Gs, Ga, P, units='W')

    !------------------------------------------------------
    !  Computes latent heat flux LE (Wm-2) i.e evapotranspiration rate ET (mm/s)
    !  from Penman-Monteith equation
    !  INPUT:
    !     AE - available energy [Wm-2]
    !     VPD- vapor pressure deficit [Pa]
    !     T  - ambient air temperature [degC]
    !     Gs - surface conductance [ms-1]
    !     Ga - aerodynamic conductance [ms-1]
    !     P  - ambient pressure [Pa] 
    !        - HT: should be from atm forcing not a fixed value (10130 Pa)?
    !     units - W (Wm-2), mm (mms-1=kg m-2 s-1), mol (mol m-2 s-1)
    !  OUTPUT:
    !     x - evaporation rate in 'units'
    !------------------------------------------------------ 
  
    ! !ARGUMENTS
    real(r8), intent(in) :: AE                        ! available energy [Wm-2]
    real(r8), intent(in) :: VPD                       ! vapor pressure deficit [Pa]
    real(r8), intent(in) :: T                         ! ambient air temperature [degC]
    real(r8), intent(in) :: Gs                        ! surface conductance [ms-1]
    real(r8), intent(in) :: Ga                        ! aerodynamic conductance [ms-1]
    real(r8), intent(in) :: P                         ! ambient pressure [Pa] 
    character (len=200), intent(in) :: units          ! W (Wm-2), mm (mms-1=kg m-2 s-1), mol (mol m-2 s-1)

    ! !LOCAL VARIABLES:
    real(r8)             :: cp, rho, Mw, s, g, esat

    ! --- constants
    cp = 1004.67  ! J kg-1 K-1
    rho = 1.25    ! kg m-3
    Mw = 18e-3    ! kg mol-1
    call e_sat(T, P, s, g, esat)  ! slope of sat. vapor pressure, psycrom const
    L = 1e3 * (3147.5 - 2.37 * (T + 273.15))

    x = (s * AE + rho * cp * Ga * D) / (s + g * (1.0 + Ga / Gs))  ! Wm-2

    if (units == 'mm') then
      penman_monteith = x/L     ! kgm-2s-1 = mms-1
    else if (units == 'mol') then
      penman_monteith = x/L/Mw  ! mol m-2 s-1
    end if

    penman_monteith=max(penman_monteith, 0.0)
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
    real(r8), intent(in) :: T                         ! air temperature (degC)
    real(r8), intent(in) :: P                         ! ambient pressure (Pa)
    real(r8), intent(out):: esa                       ! saturation vapor pressure in Pa
    real(r8), intent(out):: s                         ! slope of saturation vapor pressure curve (Pa K-1)
    real(r8), intent(out):: g                         ! psychrometric constant (Pa K-1)
    
  ! !LOCAL VARIABLES:
    real(r8)             ::  NT = 273.15
    real(r8)             ::  cp = 1004.67             ! J/kg/K
            
    Lambda = 1e3 * (3147.5 - 2.37 * (T + NT))         ! lat heat of vapor [J/kg]
    esa = 1e3 * (0.6112 * exp((17.67 * T) / (T + 273.16 - 29.66)))  ! Pa

    s = 17.502 * 240.97 * esa / ((240.97 + T) ** 2)
    g = P * cp / (0.622 * Lambda)

  END SUBROUTINE e_sat

  SUBROUTINE set_soilwaterState(soilwater_state)
  ! !ARGUMENTS
    type(soilwater_type), intent(inout) :: soilwater_state                         ! 

  ! ------- updates state variables
  ! root zone
  soilwater_state%Wliq = soilwater_state%poros*soilwater_state%WatSto/soilwater_state%MaxWatSto
  soilwater_state%Wair = soilwater_state%poros-soilwater_state%Wliq
  soilwater_state%Sat  = soilwater_state%Wliq / soilwater_state%poros
  soilwater_state%Rew  = max(0.0, min((soilwater_state%Wliq - soilwater_state%Wp) / (soilwater_state%Fc - soilwater_state%Wp + eps), 1.0))
        
  ! organic top layer; maximum that can be hold is Fc
  soilwater_state%Wliq_top = soilwater_state%Fc_top * soilwater_state%WatStoTop / soilwater_state%MaxStoTop
  soilwater_state%Ree      = relative_evaporation(soilwater_state)
  END SUBROUTINE set_soilwaterState
        
  real(r8) FUNCTION hydrCond(soilwater_state)
  
  ! ------- hydrCond
  ! returns hydraulic conductivity [ms-1] based on Campbell -formulation
  ! -------
  ! !ARGUMENTS
    type(soilwater_type), intent(inout) :: soilwater_state                         ! 

    hydrCond = soilwater_state%Ksat*soilwater_state%Sat**(2*soilwater_state%beta + 3.0)
    return

  END FUNCTION hydrCond


  real(r8) FUNCTION relative_evaporation(soilwater_state)
  ! ---------------------
  ! returns relative evaporation rate from the organic top layer; loosely
  ! based on Launiainen et al. 2015 Ecol. Mod. Moss-module
  ! Returns:
  !  f - [-], array or grid of 
  ! ---------------------
    relative_evaporation = max(0.0, min(0.98*soilwater_state%Wliq_top/soilwater_state%rw_top, 1.0))
    return
  END FUNCTION relative_evaporation
  

end module spafhy_mod
