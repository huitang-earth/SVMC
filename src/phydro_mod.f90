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
    real(r8)      , intent(in)    :: tc     ! Air temperature (tc), (degrees C)
    real(r8)      , intent(in)    :: ppfd   ! Photosynthetic photon flux density (mol m-2 d-1) (incoming solar radiation from forcing data?)
    real(r8)      , intent(in)    :: vpd    ! Vapour pressure deficit (Pa) (will be calculated using pressure & humidity)
    real(r8)      , intent(in)    :: co2    ! Atmospheric CO2 concentration (ppm)
    real(r8)      , intent(in)    :: sp     ! Surface pressure (pa)
    real(r8)      , intent(in)    :: fapar  ! Fraction of absorbed photosynthetically active radiation (unitless) (will be calculated using LAI) 
    real(r8)      , intent(in)    :: kphio  ! Apparent quantum yield efficiency (unitless).
    real(r8)      , intent(in)    :: psi_soil  ! soil water potential (Mpa)
    real(r8)      , intent(in)    :: rdark = 0 ! Dark respiration \eqn{Rd} (mol C m-2)
    type(par_plant_type), intent(in)  :: par_plant           ! A list of plant hydraulic parameters (will be defined in readpara_mod.f90).
    type(par_cost_type) , intent(in)  :: par_cost            ! A list of cost parameters (will be defined in readpara_mod.f90)..
    character(len=200)  , intent(in)  :: opt_hypothesis='PM'   ! character, Either "Lc" or "PM"
                                                             ! "Lc": Least Cost (see: rpmodel_hydraulics_numerical.R)
                                                             ! "PM": Profit Maximisation (see: rpmodel_hydraulics_numerical.R)
    real(r8)      , intent(out)   :: jmax       !  The maximum rate of RuBP regeneration (umol/m2/s) at growth temperature (argument\code{tc}), calculated using
                                                ! \deqn{A_J = A_C} 
                                                !  Electron transport capacity (umol/m2/s)
    real(r8)      , intent(out)   :: dpsi       ! soil-to-leaf water potential difference (\eqn{\psi_s-\psi_l}), Mpa
    real(r8)      , intent(out)   :: gs         !  Stomatal conductance (gs, in mol C m-2 Pa-1), calculated as
                                                !  \deqn{ gs = A / (ca (1-\chi)) } where \eqn{A} is \code{gpp}\eqn{/Mc}.
    real(r8)      , intent(out)   :: a          !  electron-transport limited assimilation rate (umol/m2/s)
    real(r8)      , intent(out)   :: ci         !  leaf-internal CO2 concentration, converted to partial pressure (Pa)
    real(r8)      , intent(out)   :: chi = ci/par_photosynth_now$ca  ! Optimal ratio of leaf internal to ambient CO2 (unitless).
    real(r8)      , intent(out)   :: vcmax                         !   Carboxylation capacity (umol/m2/s)
    real(r8)      , intent(out)   :: profit                   ! Net assimilation rate after accounting for costs
    real(r8)      , intent(out)   :: chi_jmax_lim = 0      ! Analytical chi in the case of strong Jmax limitation

    ! ! Local variables 
    real     ::    
    real     ::    
    real     ::    
    real     ::    
    real     ::    
    real     ::   
    integer  :: i, j, k, t
  
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
