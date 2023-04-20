MODULE spafpy_mod

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
  !use yasso20

  implicit none

  !Public member functions:
  public :: canopy_water_flux   ! p-hydro module
  public :: topmodel   !
  public :: soil_water
  public :: soil_water_p

  private :: aerodynamics        ! aerodynamic conductances
  private :: canopy_water_snow   ! interception, evaporation and snowpack


contains

  SUBROUTINE canopy_water_flux(gs, Ta, Prec, Rg, Par, VPD, U=2.0, CO2=380.0, Rew=1.0, beta=1.0, P=101300.0)
    !
    ! !DESCRIPTION:
    ! Adapted from Gridded canopy and snow hydrology model for SpaFHy.
    ! Based on simple schemes for computing canopy water transpiration and evaporation
    ! Can be both daily or sub-daily timesteps.
    !
    ! !USES:

    ! use soilwater_mod   ! Placeholder for soil water bucket model, which will provide psi_soil?

    ! !ARGUMENTS:
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

  SUBROUTINE canopy_water_snow(dt, T, Prec, AE, D, Ra, U, snowtype):
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

    ! quality of precipitation
    Tmin = 0.0  ! 'C, below all is snow
    Tmax = 1.0  ! 'C, above all is water
    Tmelt = 0.0  ! 'C, T when melting starts

    !storage capacities mm
    Wmax = wmax * self.LAI
    Wmaxsnow = self.wmaxsnow * self.LAI

    ! melting/freezing coefficients mm/s
    Kmelt = self.Kmelt - 1.64 * self.cf / dt  ! Kuusisto E, 'Lumi Suomessa'
    Kfreeze = self.Kfreeze

    kp = self.physpara['kp']
    tau = exp(-kp*self.LAI)  ! fraction of Rn at ground

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
    ! erate = np.zeros(gridshape)
    ! ixs = np.where((Prec == 0) & (T <= Tmin))
    ! ixr = np.where((Prec == 0) & (T > Tmin))
    ! Ga = 1. / Ra  # aerodynamic conductance

    ! resistance for snow sublimation adopted from:
    ! Pomeroy et al. 1998 Hydrol proc; Essery et al. 2003 J. Climate;
    ! Best et al. 2011 Geosci. Mod. Dev.

    Ce = 0.01*((self.W + eps) / Wmaxsnow)**(-0.4)  ! exposure coeff (-)
    Sh = (1.79 + 3.0*U**0.5)                      ! Sherwood numbner (-)
    gi = Sh*self.W*Ce / 7.68 + eps                   ! m s-1

    erate[ixs] = dt / Ls[ixs] * penman_monteith((1.0 - tau[ixs])*AE[ixs],
                                                     1e3*D[ixs], T[ixs], gi[ixs],
                                                     Ga[ixs], units='W')

        # evaporation of intercepted water, mm
        gs = 1e6
        erate[ixr] = dt / Lv[ixr] * penman_monteith((1.0 - tau[ixr])*AE[ixr],
                                                     1e3*D[ixr], T[ixr], gs,
                                                     Ga[ixr], units='W')

        # ---state of precipitation [as water (fW) or as snow(fS)]
        fW = np.zeros(gridshape)
        fS = np.zeros(gridshape)

        fW[T >= Tmax] = 1.0
        fS[T <= Tmin] = 1.0

        ix = np.where((T > Tmin) & (T < Tmax))
        fW[ix] = (T[ix] - Tmin) / (Tmax - Tmin)
        fS[ix] = 1.0 - fW[ix]
        del ix

        # --- Local fluxes (mm)
        Unload = np.zeros(gridshape)  # snow unloading
        Interc = np.zeros(gridshape)  # interception
        Melt = np.zeros(gridshape)   # melting
        Freeze = np.zeros(gridshape)  # freezing
        Evap = np.zeros(gridshape)

        """ --- initial conditions for calculating mass balance error --"""
        Wo = self.W  # canopy storage
        SWEo = self.SWE  # Snow water equivalent mm

        """ --------- Canopy water storage change -----"""
        # snow unloading from canopy, ensures also that seasonal LAI development does
        # not mess up computations
        ix = (T >= Tmax)
        Unload[ix] = np.maximum(self.W[ix] - Wmax[ix], 0.0)
        self.W = self.W - Unload
        del ix
        # dW = self.W - Wo

        # Interception of rain or snow: asymptotic approach of saturation.
        # Hedstrom & Pomeroy 1998. Hydrol. Proc 12, 1611-1625;
        # Koivusalo & Kokkonen 2002 J.Hydrol. 262, 145-164.
        ix = (T < Tmin)
        Interc[ix] = (Wmaxsnow[ix] - self.W[ix]) \
                    * (1.0 - np.exp(-(self.cf[ix] / Wmaxsnow[ix]) * Prec[ix]))
        del ix
        
        # above Tmin, interception capacity equals that of liquid precip
        ix = (T >= Tmin)
        Interc[ix] = np.maximum(0.0, (Wmax[ix] - self.W[ix]))\
                    * (1.0 - np.exp(-(self.cf[ix] / Wmax[ix]) * Prec[ix]))
        del ix
        self.W = self.W + Interc  # new canopy storage, mm

        Trfall = Prec + Unload - Interc  # Throughfall to field layer or snowpack

        # evaporate from canopy and update storage
        Evap = np.minimum(erate, self.W)  # mm
        self.W = self.W - Evap

        """ Snowpack (in case no snow, all Trfall routed to floor) """
        ix = np.where(T >= Tmelt)
        Melt[ix] = np.minimum(self.SWEi[ix], Kmelt[ix] * dt * (T[ix] - Tmelt))  # mm
        del ix
        ix = np.where(T < Tmelt)
        Freeze[ix] = np.minimum(self.SWEl[ix], Kfreeze * dt * (Tmelt - T[ix]))  # mm
        del ix

        # amount of water as ice and liquid in snowpack
        Sice = np.maximum(0.0, self.SWEi + fS * Trfall + Freeze - Melt)
        Sliq = np.maximum(0.0, self.SWEl + fW * Trfall - Freeze + Melt)

        PotInf = np.maximum(0.0, Sliq - Sice * self.R)  # mm
        Sliq = np.maximum(0.0, Sliq - PotInf)  # mm, liquid water in snow

        # update Snowpack state variables
        self.SWEl = Sliq
        self.SWEi = Sice
        self.SWE = self.SWEl + self.SWEi
        
        # mass-balance error mm
        MBE = (self.W + self.SWE) - (Wo + SWEo) - (Prec - Evap - PotInf)

        return PotInf, Trfall, Evap, Interc, MBE, Unload


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

    zm = hc + zm  ! m
    zg = min(zg, 0.1 * hc)
    kv = 0.4  ! von Karman constant (-)
    beta = 285.0  ! s/m, from Campbell & Norman eq. (7.33) x 42.0 molm-3
    alpha = LAI / 2.0  ! wind attenuation coeff (Yi, 2008 eq. 23)
    d = 0.66*hc  ! m
    zom = 0.123*hc  ! m
    zov = 0.1*zom
    zosv = 0.1*zos

    ! solve ustar and U(hc) from log-profile above canopy
    ustar = Uo * kv / np.log((zm - d) / zom) 
    Uh = ustar / kv * np.log((hc - d) / zom)
    
    ! U(zg) from exponential wind profile
    zn = min(zg / hc, 1.0)  ! zground can't be above canopy top
    Ug = Uh * exp(alpha*(zn - 1.0))

    ! canopy aerodynamic & boundary-layer resistances (sm-1). Magnani et al. 1998 PCE eq. B1 & B5
    !ra = 1. / (kv*ustar) * np.log((zm - d) / zom)
    ra = 1./(kv**2.0 * Uo) * log((zm - d) / zom) * log((zm - d) / zov)    
    rb = 1. / LAI * beta * ((w / Uh)*(alpha / (1.0 - exp(-alpha / 2.0))))**0.5

    ! soil aerodynamic resistance (sm-1)
    ras = 1. / (kv**2.0*Ug) * (log(zg / zos))*log(zg / (zosv))
    
    !print('ra', ra, 'rb', rb)
    ra = ra + rb

  END SUBROUTINE aerodynamics

end module spafpy_mod
