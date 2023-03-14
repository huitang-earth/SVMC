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
  
    call calc_kmm(tc, p, kmm)  !Why does this use std. atm pressure, and not p(z)?
    
    
    par_photosynth_now <- list(
    calc_kmm(tc, p, kmm)
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

  END SUBROUTINE pmodel_hydraulics_numerical


  calc_kmm <- function( tc, patm ) {
  

  
  return(kmm)
}

  SUBROUTINE calc_kmm(tc, patm, kmm)
    !
    ! !DESCRIPTION:
    ! Calculates the Michaelis Menten coefficient for Rubisco-limited photosynthesis
    ! Calculates the Michaelis Menten coefficient of Rubisco-limited assimilation as a function of temperature and atmospheric pressure.
    ! From function "calc_kmm" in rpmodel
    ! References: 
    !    Farquhar,  G.  D.,  von  Caemmerer,  S.,  and  Berry,  J.  A.:
    !    A  biochemical  model  of photosynthetic CO2 assimilation in leaves of
    !    C 3 species, Planta, 149, 78–90, 1980.
    !
    !    Bernacchi,  C.  J.,  Singsaas,  E.  L.,  Pimentel,  C.,  Portis,  A.
    !    R.  J.,  and  Long,  S.  P.:Improved temperature response functions
    !    for models of Rubisco-limited photosyn-thesis, Plant, Cell and
    !    Environment, 24, 253–259, 2001
    
    ! !USES:
    use readpara_mod   

    ! !ARGUMENTS:
    real(r8)      , intent(in)    :: tc     ! Air temperature (tc), degrees C
    real(r8)      , intent(in)    :: patm   ! Atmospheric pressure (Pa)
    real(r8)      , intent(out)   :: kmm    ! Michaelis-Menten coefficient at specific temperature and pressure (in Pa)    

    ! ! Local variables 
    real(r8)     ::    dhac   = 79430      ! (J/mol) Activation energy, Bernacchi et al. (2001)
    real(r8)     ::    dhao   = 36380      ! (J/mol) Activation energy, Bernacchi et al. (2001)
    real(r8)     ::    kco    = 2.09476e5  ! (ppm) O2 partial pressure, Standard Atmosphere
    ! k25 parameters are not dependent on atmospheric pressure
    real(r8)     ::    kc25   = 39.97   ! Pa, value based on Bernacchi et al. (2001), converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.
    real(r8)     ::    ko25   = 27480   ! Pa, value based on Bernacchi et al. (2001), converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.
    real(r8)     ::    tk, kc, ko, po 
  
    ! conversion to Kelvin
    tk     = tc + 273.15
  
    kc     = kc25 * ftemp_arrh( tk, dha=dhac )
    ko     = ko25 * ftemp_arrh( tk, dha=dhao )
  
    po     = kco * (1e-6) * patm         ! O2 partial pressure
    kmm    = kc * (1.0 + po/ko)

  END SUBROUTINE calc_kmm

  real(r8) FUNCTION ftemp_arrh(tk, dha, tkref)

    !---------------------------------------------------------------
    ! Calculates the Arrhenius-type temperature response
    ! Given a kinetic rate at a reference temperature (argument \code{tkref})
    ! this function calculates its temperature-scaling factor
    ! following Arrhenius kinetics.
    !---------------------------------------------------------------

    ! !ARGUMENTS
    real(r8), intent(in) :: tk                        ! Air temperature (Kelvin)
    real(r8), intent(in) :: dha                       ! Activation energy (J mol-1)
    real(r8), intent(in) :: tkref = 298.15            ! tkref Reference temperature (Kelvin)
    
    ! !LOCAL VARIABLES:
    real(r8)             :: kR=8.3145                 ! Universal gas constant, J/mol/K
  
    ! Note that the following forms are equivalent:
    ! ftemp_arrh = exp( dha * (tk - 298.15) / (298.15 * kR * tk) )
    ! ftemp_arrh = exp( dha * (tc - 25.0)/(298.15 * kR * (tc + 273.15)) )
    ! ftemp_arrh = exp( (dha/kR) * (1/298.15 - 1/tk) )
    ftemp_arrh = exp(dha * (tk - tkref) / (tkref * kR * tk))
  
    return
  
  END FUNCTION ftemp_arrh

  real(r8) FUNCTION gammastar(tk, dha, tkref)

   Calculates the CO2 compensation point
#'
#' Calculates the photorespiratory CO2 compensation point in absence of dark
#' respiration, \eqn{\Gamma*} (Farquhar, 1980).
#'
#' @param tc Temperature, relevant for photosynthesis (degrees Celsius)
#' @param patm Atmospheric pressure (Pa)
#'
#' @details The temperature and pressure-dependent photorespiratory
#' compensation point in absence of dark respiration \eqn{\Gamma* (T,p)}
#' is calculated from its value at standard temperature (\eqn{T0 = 25 }deg C)
#' and atmospheric pressure (\eqn{p0 = 101325} Pa), referred to as \eqn{\Gamma*0},
#' quantified by Bernacchi et al. (2001) to 4.332 Pa (their value in molar
#' concentration units is multiplied here with 101325 Pa to yield 4.332 Pa).
#' \eqn{\Gamma*0} is modified by temperature following an Arrhenius-type temperature
#' response function \eqn{f(T, \Delta Ha)} (implemented by \link{ftemp_arrh})
#' with activation energy \eqn{\Delta Ha = 37830} J mol-1  and is corrected for
#' atmospheric pressure \eqn{p(z)} (see \link{calc_patm}) at elevation \eqn{z}.
#' \deqn{
#'       \Gamma* = \Gamma*0 f(T, \Delta Ha) p(z) / p_0
#' }
#' \eqn{p(z)} is given by argument \code{patm}.
#'
#' @references Farquhar,  G.  D.,  von  Caemmerer,  S.,  and  Berry,  J.  A.:
#'             A  biochemical  model  of photosynthetic CO2 assimilation in leaves of
#'             C 3 species, Planta, 149, 78–90, 1980.
#'
#'             Bernacchi,  C.  J.,  Singsaas,  E.  L.,  Pimentel,  C.,  Portis,  A.
#'             R.  J.,  and  Long,  S.  P.:Improved temperature response functions
#'             for models of Rubisco-limited photosyn-thesis, Plant, Cell and
#'             Environment, 24, 253–259, 2001
#'
#' @return A numeric value for \eqn{\Gamma*} (in Pa)
#'
#' @examples print("CO2 compensation point at 20 degrees Celsius and standard atmosphere (in Pa):")
#' print(calc_gammastar(20, 101325))
#'
#' @export
#'
calc_gammastar <- function( tc, patm ) {
  
  # (J/mol) Activation energy, Bernacchi et al. (2001)
  dha    <- 37830
  
  # Pa, value based on Bernacchi et al. (2001), 
  # converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.
  gs25_0 <- 4.332
  
  gammastar <- gs25_0 * patm / calc_patm(0.0) * ftemp_arrh( (tc + 273.15), dha=dha )
  
  return( gammastar )
}
  
  END FUNCTION gammastar
  
  real(r8) FUNCTION ftemp_kphio(tk, dha, tkref)

    #' Calculates the temperature dependence of the quantum yield efficiency
#'
#' Calculates the temperature dependence of the quantum yield efficiency
#' following the temperature dependence of the maximum quantum yield of photosystem II
#' in light-adapted tobacco leaves, determined by Bernacchi et al. (2003)
#'
#' @param tc Temperature, relevant for photosynthesis (degrees Celsius)
#' @param c4 Boolean specifying whether fitted temperature response for C4 plants
#' is used. Defaults to \code{FALSE} (C3 photoynthesis temperature resposne following
#' Bernacchi et al., 2003 is used).
#'
#' @details The temperature factor for C3 photosynthesis (argument \code{c4 = FALSE}) is calculated
#' based on Bernacchi et al. (2003) as
#' 			\deqn{
#' 				\phi(T) = 0.352 + 0.022 T - 0.00034 T^2
#'       }
#'
#' The temperature factor for C4 (argument \code{c4 = TRUE}) photosynthesis is calculated based on
#' pers. comm. by David Orme, correcting values provided in Cai & Prentice (2020). Corrected 
#' parametrisation is:
#' 			\deqn{
#' 				\phi(T) = -0.064 + 0.03 T - 0.000464 T^2 
#'       }
#'
#' The factor \eqn{\phi(T)} is to be multiplied with leaf absorptance and the fraction
#' of absorbed light that reaches photosystem II. In the P-model these additional factors
#' are lumped into a single apparent quantum yield efficiency parameter (argument \code{kphio}
#' to function \link{rpmodel}).
#'
#' @return A numeric value for \eqn{\phi(T)}
#'
#' @examples
#' ## Relative change in the quantum yield efficiency
#' ## between 5 and 25 degrees celsius (percent change):
#' print(paste((ftemp_kphio(25.0)/ftemp_kphio(5.0)-1)*100 ))
#'
#' @references  
#' Bernacchi, C. J., Pimentel, C., and Long, S. P.:  In vivo temperature
#' 				response func-tions  of  parameters required  to  model  RuBP-limited
#' 				photosynthesis,  Plant  Cell Environ., 26, 1419–1430, 2003
#' Cai, W., and Prentice, I. C.: Recent trends in gross primary production 
#'        and their drivers: analysis and modelling at flux-site and global scales,
#'        Environ. Res. Lett. 15 124050 https://doi.org/10.1088/1748-9326/abc64e, 2020
#'
#' @export
#'
ftemp_kphio <- function( tc, c4 = FALSE ){
  
  if (c4){
    ftemp = -0.064 + 0.03 * tc - 0.000464 * tc^2     # correcting erroneous values provided in Cai & Prentice, 2020, according to D. Orme (issue #19) 
    # XXX THIS IS NOT CORRECT: ftemp = -0.008 + 0.00375 * tc - 0.58e-4 * tc^2   # Based on calibrated values by Shirley
  } else {
    ftemp <- 0.352 + 0.022 * tc - 3.4e-4 * tc^2
  }
  
  ## avoid negative values
  ftemp <- ifelse(ftemp < 0.0, 0.0, ftemp)
  
  return(ftemp)
}
  
  END FUNCTION ftemp_kphio
  
  real(r8) FUNCTION viscosity_h2o(tk, dha, tkref)

#' Viscosity of water
#'
#' Calculates the viscosity of water as a function of temperature and atmospheric
#' pressure.
#'
#' @param tc numeric, air temperature (tc), degrees C
#' @param p numeric, atmospheric pressure (p), Pa
#'
#' @return numeric, viscosity of water (mu), Pa s
#'
#' @examples print("Density of water at 20 degrees C and standard atmospheric pressure:")
#' print(density_h2o(20, 101325))
#'
#' @references  Huber, M. L., R. A. Perkins, A. Laesecke, D. G. Friend, J. V.
#' Sengers, M. J. Assael, ..., K. Miyagawa (2009) New
#' international formulation for the viscosity of H2O, J. Phys.
#' Chem. Ref. Data, Vol. 38(2), pp. 101-125.
#'
#' @export
#'
viscosity_h2o <- function(tc, p) {
  
  # Define reference temperature, density, and pressure values:
  tk_ast  <- 647.096    # Kelvin
  rho_ast <- 322.0      # kg/m^3
  mu_ast  <- 1e-6       # Pa s
  
  # Get the density of water, kg/m^3
  rho <- density_h2o(tc, p)
  
  # Calculate dimensionless parameters:
  tbar  <- (tc + 273.15)/tk_ast
  tbarx <- tbar^(0.5)
  tbar2 <- tbar^2
  tbar3 <- tbar^3
  rbar  <- rho/rho_ast
  
  # Calculate mu0 (Eq. 11 & Table 2, Huber et al., 2009):
  mu0 <- 1.67752 + 2.20462/tbar + 0.6366564/tbar2 - 0.241605/tbar3
  mu0 <- 1e2*tbarx/mu0
  
  # Create Table 3, Huber et al. (2009):
  h_array <- array(0.0, dim=c(7,6))
  h_array[1,] <- c(0.520094, 0.0850895, -1.08374, -0.289555, 0.0, 0.0)  # hj0
  h_array[2,] <- c(0.222531, 0.999115, 1.88797, 1.26613, 0.0, 0.120573) # hj1
  h_array[3,] <- c(-0.281378, -0.906851, -0.772479, -0.489837, -0.257040, 0.0) # hj2
  h_array[4,] <- c(0.161913,  0.257399, 0.0, 0.0, 0.0, 0.0) # hj3
  h_array[5,] <- c(-0.0325372, 0.0, 0.0, 0.0698452, 0.0, 0.0) # hj4
  h_array[6,] <- c(0.0, 0.0, 0.0, 0.0, 0.00872102, 0.0) # hj5
  h_array[7,] <- c(0.0, 0.0, 0.0, -0.00435673, 0.0, -0.000593264) # hj6
  
  # Calculate mu1 (Eq. 12 & Table 3, Huber et al., 2009):
  mu1 <- 0.0
  ctbar <- (1.0/tbar) - 1.0
  # print(paste("ctbar",ctbar))
  # for i in xrange(6):
  for (i in 1:6){
    coef1 <- ctbar^(i-1)
    # print(paste("i, coef1", i, coef1))
    coef2 <- 0.0
    for (j in 1:7){
      coef2 <- coef2 + h_array[j,i] * (rbar - 1.0)^(j-1)
    }
    mu1 <- mu1 + coef1 * coef2
  }
  mu1 <- exp( rbar * mu1 )
  # print(paste("mu1",mu1))
  
  # Calculate mu_bar (Eq. 2, Huber et al., 2009)
  #   assumes mu2 = 1
  mu_bar <- mu0 * mu1
  
  # Calculate mu (Eq. 1, Huber et al., 2009)
  mu <- mu_bar * mu_ast    # Pa s
  
  return( mu )
}

  
    return
  
  END FUNCTION viscosity_h2o
  
  
  #' Density of water
#'
#' Calculates the density of water as a function of temperature and atmospheric
#' pressure, using the Tumlirz Equation.
#'
#' @param tc numeric, air temperature (tc), degrees C
#' @param p numeric, atmospheric pressure (p), Pa
#'
#' @return numeric, density of water, kg/m^3
#'
#' @examples 
#'  # Density of water at 20 degrees C and standard atmospheric pressure
#'  print(density_h2o(20, 101325))
#'
#' @references  F.H. Fisher and O.E Dial, Jr. (1975) Equation of state of
#' pure water and sea water, Tech. Rept., Marine Physical
#' Laboratory, San Diego, CA.
#'
#' @export
#'
density_h2o <- function(tc, p){
  
  # Calculate lambda, (bar cm^3)/g:
  my_lambda <- 1788.316 +
    21.55053*tc +
    -0.4695911*tc*tc +
    (3.096363e-3)*tc*tc*tc +
    -(7.341182e-6)*tc*tc*tc*tc
  
  # Calculate po, bar
  po <- 5918.499 +
    58.05267*tc +
    -1.1253317*tc*tc +
    (6.6123869e-3)*tc*tc*tc +
    -(1.4661625e-5)*tc*tc*tc*tc
  
  # Calculate vinf, cm^3/g
  vinf <- 0.6980547 +
    -(7.435626e-4)*tc +
    (3.704258e-5)*tc*tc +
    -(6.315724e-7)*tc*tc*tc +
    (9.829576e-9)*tc*tc*tc*tc +
    -(1.197269e-10)*tc*tc*tc*tc*tc +
    (1.005461e-12)*tc*tc*tc*tc*tc*tc +
    -(5.437898e-15)*tc*tc*tc*tc*tc*tc*tc +
    (1.69946e-17)*tc*tc*tc*tc*tc*tc*tc*tc +
    -(2.295063e-20)*tc*tc*tc*tc*tc*tc*tc*tc*tc
  
  # Convert pressure to bars (1 bar <- 100000 Pa)
  pbar <- (1e-5)*p
  
  # Calculate the specific volume (cm^3 g^-1):
  v <- vinf + my_lambda/(po + pbar)
  
  # Convert to density (g cm^-3) -> 1000 g/kg; 1000000 cm^3/m^3 -> kg/m^3:
  rho <- (1e3/v)
  
  return(rho)
}

  
end module phydro_mod
