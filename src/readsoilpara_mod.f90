MODULE readsoilpara_mod

implicit none
  type, public :: soilwater_type
    real(8)  :: PondSto ! [m] - pond storage
    real(8)  :: Inflow  ! [m] - total inflow to root zone during timestep
    real(8)  :: Roff    ! [m] - surface runoff -"-
    real(8)  :: Drain   ! [m] - drainage from root zone -"-
    real(8)  :: Interc  ! [m] - interception of top layer -"-
    real(8)  :: WatSto  ! [m] - root zone storage
    real(8)  :: MaxWatSto ! [m] - root zone storage capacity
    real(8)  :: WatStoTop ! [m] - top layer storage
    real(8)  :: MaxStoTop ! [m] - top layer storage capacity
    real(8)  :: Wliq      ! [m3 m-3] - root zone water content
    real(8)  :: Wliq_top  ! [m3 m-3] - top layer water content
    real(8)  :: Psi  ! water potential, root zone (MPa)
    real(8)  :: PsiTop  ! water potential, organic top layer (MPa)
    real(8)  :: Sat   ! saturation ratio (-), root zone
    real(8)  :: Kh   ! Hydraulic conductivity at Sat [m s-1]
    real(8)  :: beta  ! modifier for soil evaporation rate, WliqTop/FCtop

    real(8)  :: mbe     ! [m] - mass balance error
  end type soilwater_type

  type, public :: canopywater_type
    real(8):: CanopyStorage  ! canopy water storage (mm = kg m-2(ground)) 
    real(8):: Trfall     ! throughfall to snow / soil surface (mm, during timestep)
    real(8):: Interc     ! interception of canopy (mm)
    real(8):: CanopyEvap ! evaporation / sublimation from canopy store (mm)
    real(8):: Unload     ! undloading from canopy storage (mm)    
    real(8):: GroundEvap ! evaporation from ground (mm)
    real(8):: MBE        ! mass balance error (mm)      
    real(8):: ET         ! total evapo-transpiration (mm)
    real(8):: Transpi    ! transpiration rate (mm)  
  end type canopywater_type

  type, public :: snowwater_type
    real(8):: swe        ! [mm] snow water equivalent
    real(8):: SWEi      ! [mm] snow water equivalent as ice
    real(8):: SWEl     ! [mm] snow water equivalent as liquid
    real(8):: PotInf     ! potential infiltration to soil profile (mm)            
  end type snowwater_type

! Soil hydraulic properties for spafhy
  real(8) :: soil_depth    ! root zone depth (m)
  real(8) :: max_poros     ! [m3 m-3], porosity
  real(8) :: fc            ! [m3 m-3], field capacity. For consistency, must be computed from soil water retention curve at Psi= xx KPa
  real(8) :: wp            ! [m3 m-3], wilting point. For consistency, must be computed from soil water retention curve at Psi=xx KPa
  real(8) :: ksat          ! [m s-1], saturated hyd, conductivity
  real(8) :: beta          ! [-], term for soil evaporation resistance (Wliq/FC). THIS IS in soilwater_state%beta

  ! organic (moss) layer
  real(8) :: org_depth     ! depth of organic top layer (m)
  real(8) :: org_poros      ! porosity (-)
  real(8) :: org_fc         ! field capacity (-)
  !real(8) :: org_rw        ! critical vol. moisture content (-) for decreasing phase in Ef
  real(8) :: maxpond        ! max ponding allowed (m)
  
  ! initial states: rootzone and toplayer soil saturation ratio [-] and pond storage [m]
  !real(8) :: rootzone_sat   ! root zone saturation ratio (-)
  real(8) :: org_sat        ! organic top layer saturation ratio (-)
  !real(8) :: pond_sto       ! pond storage
  !real(8) :: soilcode    = -1   ! site-specific values

  ! soil retention curve (Van Genuchten) parameters
  !--- parameters from Launiainen et al. Forests, 2022  
  real(8) :: n_van       !Launiainen et al. 2022: C1-5: 1.12, 1.14, 1.07, 1.27, 1.18    
  real(8) :: watres  !Launiainen et al. 2022: C1-5: 0.0
  real(8) :: alpha_van   !Launiainen et al. 2022: C1-5: 4.45, 5.92, 2.02, 4.49, 3.35
  real(8) :: watsat  !Launiainen et al. 2022: C1-5: 0.75, 0.68, 0.46, 0.47, 0.54  


contains

  subroutine readsoilhydro_namelist
  
    implicit none

    logical :: old
    integer :: readerror
    integer,parameter :: unitsoilhydro=3

    namelist /soilhydro_namelist/ &
      soil_depth, &
      max_poros, &
      fc, &
      wp, &
      ksat, &
      beta, &
      org_depth, & 
      org_poros, &
      org_fc, & 
      !org_rw, &
      maxpond, &
      !rootzone_sat, &
      org_sat, & 
      !pond_sto, &  
      n_van, &     
      watres, &
      alpha_van, &
      watsat

    old=.false.
  ! Presetting namelist command
    soil_depth=0.4
    max_poros=0.46           ! should be equivalent to watsat here.
    fc=0.36                   ! based on C3 in Launiainen et al. 2022 ! Must be computed from water-retention curve
    wp=0.22     ! Must be computed from water-retention curve
    ksat=2.0e-6
    !beta=4.7            ! default 
    org_depth=0.04
    org_poros=0.9
    org_fc=0.3        
    !org_rw=0.24       
    maxpond=0.0
    !rootzone_sat= 0.6 
    org_sat     = 1.0
    !pond_sto    = 0.0
    n_van=1.07            !Launiainen et al. 2022: C1-5: 1.12, 1.14, 1.07, 1.27, 1.18 
    watres=0.0            !Launiainen et al. 2022: C1-5: 0.0
    alpha_van=2.02          !Launiainen et al. 2022: C1-5: 4.45, 5.92, 2.02, 4.49, 3.35
    watsat=0.46          !Launiainen et al. 2022: C1-5: 0.75, 0.68, 0.46, 0.47, 0.54  

  ! Reading namelist
  !  open(unitsoilhydro, file='./soilhydro_namelist', status='old', form='formatted', err=999)
  !  read(unitsoilhydro, soilhydro_namelist, iostat=readerror)
  !  close(unitsoilhydro)

!999 write(*,*) ' #### MODEL ERROR! FILE "soilhydro_namelist"    #### '
!    write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
!    stop

  end subroutine readsoilhydro_namelist


!  subroutine readsoilyasso_namelist
  
!    implicit none

!    logical :: old
!    integer :: readerror
!    integer,parameter :: unitsoilyasso=4

!    old=.false.

!    namelist /soilyasso_namelist/ &
!       param_y20_map
!       nc_mb
!       cue_min
!       nc_h_max
!       awenh_fineroot
!       awenh_leaf
!       awenh_soluble
!       awenh_compost
    
  ! Presetting namelist command
!    param_y20_map(num_params_y20) = (/ &
!     0.51, &
!     5.19, &
!     0.13, &
!     0.1, &
!     0.5, &
!     0., &
!     1., &
!     1., &
!     0.99, &
!     0., &
!     0., &
!     0., &
!     0., &
!     0., &
!     0.163, &
!     0., &
!     -0., &
!     0., &
!     0., &
!     0., &
!     0., &
!     0.158, &
!     -0.002, &
!     0.17, &
!     -0.005, &
!     0.067, &
!     -0., &
!     -1.44, &
!     -2.0, &
!     -6.9, &
!     0.0042, &
!     0.0015, &
!     -2.55, &
!     1.24, &
!     0.25/)
!    nc_mb = 0.1
!    cue_min = 0.1
!    nc_h_max = 0.1
!    awenh_fineroot(statesize_yasso) = (/0.46, 0.32, 0.04, 0.18, 0.0/)
!    awenh_leaf(statesize_yasso) = (/0.46, 0.32, 0.04, 0.18, 0.0/)
!    awenh_soluble(statesize_yasso) = (/0.0, 1.0, 0.0, 0.0, 0.0/)
!    awenh_compost(statesize_yasso) = (/0.69, 0.09, 0.02, 0.20, 0.0/)

  ! Reading namelist
!    open(unitsoilyasso, file='./soilyasso_namelist', status='old', form='formatted', err=999)
!    read(unitsoilyasso, soilyasso_namelist, iostat=readerror)
!    close(unitsoilyasso)

!999 write(*,*) ' #### MODEL ERROR! FILE "soilyasso_namelist"    #### '
!    write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
!    stop
!  end subroutine readsoilyasso_namelist

END MODULE readsoilpara_mod

!----------------------------
  ! Soil properties for yasso
  !-----------------------------
  !real(r8), parameter :: days_yr = 365.0
  !integer, parameter, public :: statesize_yasso = 5

! The yasso parameter vector:
! 1-16 matrix A entries: 4*alpha, 12*p
! 17-21 Leaching parameters: w1,...,w5 IGNORED IN THIS FUNCTION
! 22-23 Temperature-dependence parameters for AWE fractions: beta_1, beta_2
! 24-25 Temperature-dependence parameters for N fraction: beta_N1, beta_N2
! 26-27 Temperature-dependence parameters for H fraction: beta_H1, beta_H2
! 28-30 Precipitation-dependence parameters for AWE, N and H fraction: gamma, gamma_N, gamma_H
! 31-32 Humus decomposition parameters: p_H, alpha_H (Note the order!)
! 33-35 Woody parameters: theta_1, theta_2, r 

! The Yasso20 maximum a posteriori parameters:
  !integer, public :: num_params_y20
  !real, public :: param_y20_map(num_params_y20)

  ! Nitrogen-specific parameters
  !real, public :: nc_mb ! N-C ratio of the microbial biomass 
  !real, public :: cue_min  ! minimum microbial carbon use efficiency
  !real, public :: nc_h_max ! N-C ratio of the H pool

  ! AWENH composition from Palosuo et al. (2015), for grasses. For now, we'll use the same
  ! composition for both above and below ground inputs. The last values (H) are always 0.
  !real :: awenh_fineroot(statesize_yasso)
  !real :: awenh_leaf(statesize_yasso)
  ! A soil amendment consisting of soluble carbon (and nitrogen)
  !real :: awenh_soluble(statesize_yasso)
  ! From Heikkinen et al 2021, composted horse manure with straw litter
  !real :: awenh_compost(statesize_yasso)

  !integer, parameter, public :: met_ind_init = 1