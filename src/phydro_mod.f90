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

  type, public :: par_env_type
    real(r8) :: viscosity_water     
    real(r8) :: density_water             
    real(r8) :: patm
    real(r8) :: tc
    real(r8) :: vpd                   
  end type par_env_type

  type, public :: par_photosynth_type
    real(r8) :: kmm  
    real(r8) :: gammastar             
    real(r8) :: phi0
    real(r8) :: Iabs
    real(r8) :: ca  
    real(r8) :: patm
    real(r8) :: delta      
  end par_photosynth_type

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
    type(par_cost_type) , intent(in)  :: par_cost            ! A list of cost parameters (will be defined in readpara_mod.f90).

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

    type(par_env_type)            :: par_env_now           ! 
    type(par_photosynth_type)     :: par_cost              ! 
    type(par_plant_type)          :: par_plant_now         ! 
    type(par_cost_type)           :: par_cost_now          !

    real     ::    fn_profit
    real     ::    lj_dps
    real     ::    
    real     ::    
    real     ::    
    real     ::   
    integer  :: i, j, k, t
  
    call calc_kmm(tc, p, par_photosynth_now%kmm)  !Why does this use std. atm pressure, and not p(z)?
    par_photosynth_now%gammastar = gammastar(tc, sp)
    par_photosynth_now%phi0 = kphio*ftemp_kphio(tc)
    par_photosynth_now%Iabs = ppfd*fapar
    par_photosynth_now%ca = co2*sp*1e-6             ! Convert to partial pressure
    par_photosynth_now%patm = sp,
    par_photosynth_now%delta = rdark
  
    par_env_now%viscosity_water = viscosity_h2o(tc, sp)
    par_env_now%density_water = density_h2o(tc, sp)
    par_env_now%patm = sp
    par_env_now%tc = tc
    par_env_now%vpd = vpd
  
    par_plant_now = par_plant
    par_cost_now = par_cost
    
    !! if par_cost is empty, use pre-defined parameter
    !if (opt_hypothesis == "PM") then
    !  par_cost_now%alpha = 0.1          ! cost of Jmax
    !  par_cost_now%gamma = 1.0          ! cost of hydraulic repair
    !else if (opt_hypothesis == "LC") then
    !  par_cost_now%alpha = 0.1          ! cost of Jmax
    !  par_cost_now%gamma = 0.5          ! cost of hydraulic repair
    !end if
  
     call optimise_midterm_multi(fn_profit, psi_soil, par_cost_now, par_photosynth_now, par_plant_now, par_env_now, opt_hypothesis)

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


  SUBROUTINE optimise_midterm_multi(fn_profit, psi_soil, par_cost, par_photosynth, par_plant, par_env, return_all = FALSE, opt_hypothesis)
    


  END SUBROUTINE optimise_midterm_multi

  optimise_midterm_multi <- function(){
  
  out_optim <- optimr::optimr(
    par       = c(logjmax=0, dpsi=1),  
    lower     = c(-10, .0001),
    upper     = c(10, 1e6),
    fn        = fn_profit,
    psi_soil  = psi_soil,
    par_cost  = par_cost,
    par_photosynth = par_photosynth,
    par_plant = par_plant,
    par_env   = par_env,
    do_optim  = TRUE,
    opt_hypothesis = opt_hypothesis,
    method    = "L-BFGS-B",
    control   = list( maxit = 500, maximize = TRUE, fnscale=1e4 )
  )
  
  out_optim$value <- -out_optim$value
  
  if (return_all){
    out_optim
  } else {
    return(out_optim$par)
  }
}


  SUBROUTINE calc_kmm(tc, patm, kmm)
    !-------------------------------------------------------------------------
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
    !--------------------------------------------------------------------------

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

  real(r8) FUNCTION gammastar(tc, patm)
    
    !---------------------------------------------------------------
    ! Calculates the CO2 compensation point
    ! Calculates the photorespiratory CO2 compensation point in absence of dark
    ! respiration, \eqn{\Gamma*} (Farquhar, 1980).
    ! Temperature and pressure-dependent
    !---------------------------------------------------------------

    ! !ARGUMENTS
    real(r8), intent(in) :: tc                        ! Temperature, relevant for photosynthesis (degrees Celsius)
    real(r8), intent(in) :: patm                      ! Atmospheric pressure (Pa)
    
    ! !LOCAL VARIABLES:
    real(r8)             ::  dha    = 37830           ! Activation energy (J/mol), Bernacchi et al. (2001)
    real(r8)             ::  gs25_0 = 4.332           ! Photorespiratory CO2 compensation point at standard temperature 
                                                      ! (T = 25 degC, p0 = 101325) (Pa)
                                                      ! Quantified by Bernacchi et al. (2001) to 4.332 Pa 
                                                      ! Their value in molar concentration units is multiplied with 101325 Pa to yield 4.332 Pa                                                 
    real(r8)             ::  patm0 = 101325           ! Atmospheric pressure at sea level (Pa), defaults to 101325 Pa.

    gammastar = gs25_0 * patm / patm0 * ftemp_arrh((tc + 273.15), dha=dha)
  
    return
  
  END FUNCTION gammastar
  
  real(r8) FUNCTION ftemp_kphio(tc, c4 = FALSE )
    
    !---------------------------------------------------------------
    ! Calculates the temperature dependence of the quantum yield efficiency.
    ! following the temperature dependence of the maximum quantum yield of photosystem II
    ! in light-adapted tobacco leaves, determined by Bernacchi et al. (2003).
    ! The factor is to be multiplied with leaf absorptance and the fraction
    ! of absorbed light that reaches photosystem II.
    ! From function "ftemp_kphio" in p-model
    !---------------------------------------------------------------

    ! !ARGUMENTS
    real(r8), intent(in) :: tc                        ! Temperature, relevant for photosynthesis (degrees Celsius)
    logical,  intent(in) :: c4                        ! Boolean specifying whether fitted temperature response for C4 plants
                                                      ! Defaults to FALSE (C3 photoynthesis temperature resposne following
                                                      ! Bernacchi et al., 2003 is used

    if (c4) then
      ! correcting erroneous values provided in Cai & Prentice, 2020, according to D. Orme (issue #19) 
      ! XXX THIS IS NOT CORRECT: ftemp = -0.008 + 0.00375 * tc - 0.58e-4 * tc^2   # Based on calibrated values by Shirley
      ftemp_kphio = -0.064 + 0.03 * tc - 0.000464 * tc^2     
    else 
      ! The temperature factor for C3 photosynthesis (c4 = FALSE) is calculated based on Bernacchi et al. (2003)
      ftemp_kphio = 0.352 + 0.022 * tc - 3.4e-4 * tc^2
    
    end if
  
    ! Avoid negative values
    if (ftemp_kphio<0.0) then
      ftemp_kphio=0.0
    end if

    return

  END FUNCTION ftemp_kphio


  real(r8) FUNCTION viscosity_h2o(tc, patm)

    !---------------------------------------------------------------
    ! Calculates the viscosity of water ((mu), Pa s) as a function of temperature and atmospheric pressure.
    ! From function "viscosity_h2o" in p-model
    ! References:
    ! Huber, M. L., R. A. Perkins, A. Laesecke, D. G. Friend, J. V. Sengers, M. J. Assael, ..., K. Miyagawa (2009) 
    !    New international formulation for the viscosity of H2O, J. Phys. Chem. Ref. Data, Vol. 38(2), pp. 101-125.
    !---------------------------------------------------------------

    ! !ARGUMENTS
    real(r8), intent(in) :: tc                        ! Air temperature, degrees C
    real(r8), intent(in) :: patm                      ! Atmospheric pressure (Pa)

    ! ! Local variables 
    ! Define reference temperature, density, and pressure values
    real(r8)     ::    tk_ast  = 647.096    ! Kelvin
    real(r8)     ::    rho_ast = 322.0      ! kg/m^3
    real(r8)     ::    mu_ast  = 1e-6       ! Pa s
    real(r8)     ::    rho                  ! density of water, kg/m^3
    real(r8)     ::    tbar, tbarx, tbar2, tbar3, rbar, ctbar     ! dimensionless parameters
    real(r8)     ::    coef1,
    real(r8)     ::    mu0, mu1            !
    
    real(r8), dimension (1:7,1:6)     ::    h_array     ! Create Table 3, Huber et al. (2009)

    integer      ::    i, j

    ! Get the density of water, kg/m^3
    rho = density_h2o(tc, patm)
  
    ! Calculate dimensionless parameters:
    tbar  = (tc + 273.15)/tk_ast
    tbarx = tbar^(0.5)
    tbar2 = tbar^2
    tbar3 = tbar^3
    rbar  = rho/rho_ast
  
    ! Calculate mu0 (Eq. 11 & Table 2, Huber et al., 2009):
    mu0 = 1.67752 + 2.20462/tbar + 0.6366564/tbar2 - 0.241605/tbar3
    mu0 = 1e2*tbarx/mu0
  
    ! Create Table 3, Huber et al. (2009):

    h_array(1,:) = (/0.520094, 0.0850895, -1.08374, -0.289555, 0.0, 0.0/)  ! hj0
    h_array(2,:) = (/0.222531, 0.999115, 1.88797, 1.26613, 0.0, 0.120573/) ! hj1
    h_array(3,:) = (/-0.281378, -0.906851, -0.772479, -0.489837, -0.257040, 0.0/) ! hj2
    h_array(4,:) = (/0.161913,  0.257399, 0.0, 0.0, 0.0, 0.0/) ! hj3
    h_array(5,:) = (/-0.0325372, 0.0, 0.0, 0.0698452, 0.0, 0.0/) ! hj4
    h_array(6,:) = (/0.0, 0.0, 0.0, 0.0, 0.00872102, 0.0/) ! hj5
    h_array(7,:) = (/0.0, 0.0, 0.0, -0.00435673, 0.0, -0.000593264/) ! hj6
  
    ! Calculate mu1 (Eq. 12 & Table 3, Huber et al., 2009):
    mu1 = 0.0
    ctbar = (1.0/tbar) - 1.0

    do i=1,6
      coef1 = ctbar^(i-1)
      coef2 = 0.0
      do j=1,7
        coef2 = coef2 + h_array(j,i) * (rbar - 1.0)^(j-1)
      end do
      mu1 = mu1 + coef1 * coef2
    end do

    mu1 = exp(rbar*mu1)
  
    ! Calculate mu_bar (Eq. 2, Huber et al., 2009)
    !  assumes mu2 = 1
    mu_bar = mu0 * mu1
  
    ! Calculate mu (Eq. 1, Huber et al., 2009)
    viscosity_h2o = mu_bar * mu_ast    ! Pa s
  
    return
  
  END FUNCTION viscosity_h2o
  
  real(r8) FUNCTION density_h2o(tc, patm)

    !---------------------------------------------------------------
    ! Calculates the density of water (kg/m^3) as a function of temperature and atmospheric pressure, using the Tumlirz Equation.
    ! From function "viscosity_h2o" in p-model
    ! References:
    ! F.H. Fisher and O.E Dial, Jr. (1975) Equation of state of pure water and sea water, Tech. Rept., Marine Physical Laboratory, San Diego, CA.
    !---------------------------------------------------------------

    ! !ARGUMENTS
    real(r8), intent(in) :: tc                        ! Air temperature, degrees C
    real(r8), intent(in) :: patm                      ! Atmospheric pressure (Pa)

    ! ! Local variables 
    real(r8)     ::    my_lambda    ! lambda, (bar cm^3)/g
    real(r8)     ::    po           ! bar
    real(r8)     ::    vinf, v      ! cm^3/g
    real(r8)     ::    pbar         ! pressure to bars
  
    ! Calculate lambda, (bar cm^3)/g
    my_lambda = 1788.316 + 21.55053*tc +  &
                -0.4695911*tc*tc + (3.096363e-3)*tc*tc*tc + &
                -(7.341182e-6)*tc*tc*tc*tc
  
    ! Calculate po, bar
    po = 5918.499 +  &
          58.05267*tc + &
          -1.1253317*tc*tc + &
          (6.6123869e-3)*tc*tc*tc + &
          -(1.4661625e-5)*tc*tc*tc*tc
  
    ! Calculate vinf, cm^3/g
    vinf = 0.6980547 + &
           -(7.435626e-4)*tc + &
          (3.704258e-5)*tc*tc + &
          -(6.315724e-7)*tc*tc*tc + &
          (9.829576e-9)*tc*tc*tc*tc + &
          -(1.197269e-10)*tc*tc*tc*tc*tc + &
          (1.005461e-12)*tc*tc*tc*tc*tc*tc + &
          -(5.437898e-15)*tc*tc*tc*tc*tc*tc*tc + &
          (1.69946e-17)*tc*tc*tc*tc*tc*tc*tc*tc + &
          -(2.295063e-20)*tc*tc*tc*tc*tc*tc*tc*tc*tc
  
    ! Convert pressure to bars (1 bar <- 100000 Pa)
    pbar = (1e-5)*p
  
    ! Calculate the specific volume (cm^3 g^-1):
    v = vinf + my_lambda/(po + pbar)
  
    ! Convert to density (g cm^-3) -> 1000 g/kg; 1000000 cm^3/m^3 -> kg/m^3:
    density_h2o = (1e3/v)

    return

  END FUNCTION density_h2o
  
end module phydro_mod
