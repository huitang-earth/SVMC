MODULE phydro_mod

! Modules
  use netcdf             ! library for processing netcdf files
  use readpara_mod       ! module for reading parameter files in ASCII
  use readclim_mod       ! module for reading reading meteorological forcing data
  use readsoil_mod       ! module for reading soil properties (shared with yasso?)
  use phydro_mod         ! module for p-hydro
  use alloc_mod          ! module for carbon allocation and yield
  !use yasso20

  implicit none

  !Public member functions:
  public :: pmodel_hydraulics_numerical   ! p-hydro module

contains
  
  !---------------------------------------------------------
  SUBROUTINE pmodel_hydraulics_numerical(tc, ppfd, vpd, co2, sp, fapar, kphio, psi_soil, rdark = 0, par_plant, par_cost = NULL, opt_hypothesis = "PM")
    !
    ! !DESCRIPTION:
    ! The hydraulic p-model (HyPE), numerical version
    ! Calculates the carboxylation capacity, as coordinated to a given electron-transport limited assimilation rate.
    !
    ! !USES:

    ! use soilwater_mod   ! Placeholder for soil water bucket model, which will provide psi_soil?

    ! !ARGUMENTS:
    real(r8)      , intent(in)    :: tc     ! Air temperature (tc), degrees C
    real(r8)      , intent(in)    :: ppfd   ! Photosynthetic photon flux density (mol m-2 d-1) (incoming solar radiation from forcing data?)
    real(r8)      , intent(in)    :: vpd    ! Vapour pressure deficit (Pa) (will be calculated using pressure & humidity)
    real(r8)      , intent(in)    :: co2    ! Atmospheric CO2 concentration (ppm)
    real(r8)      , intent(in)    :: sp     ! Surface pressure (pa)
    real(r8)      , intent(in)    :: fapar  ! Fraction of absorbed photosynthetically active radiation (unitless) (will be calculated using LAI) 
    real(r8)      , intent(in)    :: kphio  ! Apparent quantum yield efficiency (unitless).
    real(r8)      , intent(in)    :: psi_soil  ! soil water potential (Mpa)
    real(r8)      , intent(in)    :: rdark = 0 !
    type(par_plant_type), intent(in)  :: par_plant           ! A list of plant hydraulic parameters (will be defined in readpara_mod.f90).
    type(par_cost_type) , intent(in)  :: par_cost            ! A list of cost parameters.
    character(len=200)  , intent(in)  :: opt_hypothesis=''   ! character, Either "Lc" or "PM"

    real(r8)      , intent(out)   :: jmax     
    real(r8)      , intent(out)   :: dpsi
    real(r8)      , intent(out)   :: gs
    real(r8)      , intent(out)   :: a
    real(r8)      , intent(out)   :: ci
    real(r8)      , intent(out)   :: chi = ci/par_photosynth_now$ca,
    real(r8)      , intent(out)   :: vcmax
    real(r8)      , intent(out)   :: profit,    
    real(r8)      , intent(out)   :: chi_jmax_lim = 0

    ! ! Local variables 

  real     ::    
  real     ::    
  real     ::    
  real     ::    
  real     ::    
  real     ::   


     !
  !Loop variables
  !***********************************
  integer   :: i, j, k, t
  
    par_photosynth_now <- list(
    kmm = rpmodel::calc_kmm(tc, p),  # Why does this use std. atm pressure, and not p(z)?
    gammastar = rpmodel::calc_gammastar(tc, p),
    phi0 = kphio*rpmodel::calc_ftemp_kphio(tc),
    Iabs = ppfd*fapar,
    ca = co2*p*1e-6,  # Convert to partial pressure
    patm = p,
    delta = rdark
  )
  
  par_env_now = list(
    viscosity_water = rpmodel::calc_viscosity_h2o(tc, p),  # Needs to be imported from rpmodel.R
    density_water = rpmodel::calc_density_h2o(tc, p),  # Needs to be imported from rpmodel.R
    patm = p,
    tc = tc,
    vpd = vpd
  )
  
  par_plant_now = par_plant
  
  if (!is.null(par_cost)){
    par_cost_now = par_cost
  }
  else{
    if (opt_hypothesis == "PM"){
      par_cost_now = list(
        alpha = 0.1,       # cost of Jmax
        gamma = 1          # cost of hydraulic repair
      )
    } else if (opt_hypothesis == "LC"){
      par_cost_now = list(
        alpha = .1,        # cost of Jmax
        gamma = 0.5          # cost of hydraulic repair
      )
    }
  }
  
  lj_dps = optimise_midterm_multi(fn_profit, psi_soil = psi_soil, par_cost  = par_cost_now, par_photosynth = par_photosynth_now, par_plant = par_plant_now, par_env = par_env_now, opt_hypothesis = opt_hypothesis)
  
  profit = fn_profit(par = lj_dps, psi_soil = psi_soil, par_cost  = par_cost_now, par_photosynth = par_photosynth_now, par_plant = par_plant_now, par_env = par_env_now, opt_hypothesis = opt_hypothesis)
  
  jmax = exp(lj_dps[1])
  dpsi = lj_dps[2]
  
  gs = calc_gs(dpsi=dpsi, psi_soil=psi_soil, par_plant = par_plant_now, par_env = par_env_now)
  
  a_j = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth_now)
  a = a_j$a
  ci = a_j$ci
  
  vcmax = calc_vcmax_coordinated_numerical(a,ci, par_photosynth_now)
  
  return(list(
    jmax=jmax,
    dpsi=dpsi,
    gs=gs,
    a=a,
    ci=ci,
    chi = ci/par_photosynth_now$ca,
    vcmax=vcmax,
    profit = profit,
    chi_jmax_lim = 0
  ))

  END SUBROUTINE




end module phydro_mod
