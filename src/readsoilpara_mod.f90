MODULE readsoilpara_mod

implicit none
  type, public :: soilwater_state_type
    real(8)  :: PondSto ! [m] - pond storage
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
    real(8)  :: beta    ! modifier for soil evaporation rate, WliqTop/FCtop
    real(8)  :: mbe     ! [m] - mass balance error
  end type soilwater_state_type

  type, public :: soilwater_flux_type
    real(8)  :: Inflow  ! [m] - total inflow to root zone during timestep
    real(8)  :: Roff    ! [m] - surface runoff -"-
    real(8)  :: Drain   ! [m] - drainage from root zone -"-
    real(8)  :: Interc  ! [m] - interception of top layer -"-
  end type soilwater_flux_type

  type, public :: canopywater_state_type
    real(8) :: CanopyStorage ! canopy water storage (mm = kg m-2(ground)) 
    real(8) :: swe        ! [mm] snow water equivalent
    real(8) :: SWEi      ! [mm] snow water equivalent as ice
    real(8) :: SWEl     ! [mm] snow water equivalent as liquid
    real(8) :: MBE       ! mass balance error (mm)      
  end type canopywater_state_type

  type, public :: canopywater_flux_type
    real(8) :: Trfall     ! throughfall to snow / soil surface (mm, during timestep)
    real(8) :: Interc     ! interception of canopy (mm)
    real(8) :: CanopyEvap ! evaporation / sublimation from canopy store (mm)
    real(8) :: Unload     ! undloading from canopy storage (mm)    
    real(8) :: GroundEvap ! evaporation from ground (mm)   
    real(8) :: ET         ! total evapo-transpiration (mm)
    real(8) :: Transpi    ! transpiration rate (mm)  
    real(8) :: PotInf     ! potential infiltration to soil profile (mm)    
  end type canopywater_flux_type

! Hydraulic properties (organic layer & soil) for spafhy
  type, public :: spafhy_para_type
    ! soil related parameters
    real(8) :: soil_depth    ! root zone depth (m)
    real(8) :: max_poros     ! [m3 m-3], porosity
    real(8) :: fc            ! [m3 m-3], field capacity. For consistency, must be computed from soil water retention curve at Psi= xx KPa
    real(8) :: wp            ! [m3 m-3], wilting point. For consistency, must be computed from soil water retention curve at Psi=xx KPa
    real(8) :: ksat          ! [m s-1], saturated hyd, conductivity
 
    ! soil retention curve (Van Genuchten) parameters
    !--- parameters from Launiainen et al. Forests, 2022  
    real(8) :: n_van       !Launiainen et al. 2022: C1-5: 1.12, 1.14, 1.07, 1.27, 1.18    
    real(8) :: watres      !Launiainen et al. 2022: C1-5: 0.0
    real(8) :: alpha_van   !Launiainen et al. 2022: C1-5: 4.45, 5.92, 2.02, 4.49, 3.35
    real(8) :: watsat      !Launiainen et al. 2022: C1-5: 0.75, 0.68, 0.46, 0.47, 0.54  

  
    ! organic (moss) layer
    real(8) :: org_depth     ! depth of organic top layer (m)
    real(8) :: org_poros      ! porosity (-), redundent.
    real(8) :: org_fc         ! field capacity (-)
    real(8) :: org_sat        ! organic top layer saturation ratio (-)
    !real(8) :: org_rw        ! critical vol. moisture content (-) for decreasing phase in Ef
    real(8) :: maxpond        ! max ponding allowed (m)
   
    ! Canopy water parameters
    real(8) :: wmax     ! storage capacity for rain (mm/LAI), Hui: this is too much compared to CTSM
    real(8) :: wmaxsnow ! storage capacity for snow (mm/LAI), Hui: this is reasonable
    real(8) :: hc       ! canopy height (m)
    real(8) :: cf       ! canopy closure fraction (-)
    real(8) :: w_leaf   ! leaf length scale (m)
    real(8) :: rw       ! critical value for REW (-),
    real(8) :: rwmin    ! minimum relative conductance (-)
    real(8) :: gsoil    ! soil surface conductance if soil is fully wet (m/s)
    
    ! degree-day snow model
    real(8) :: kmelt         ! melt coefficient in open (mm/s)
    real(8) :: kfreeze       ! freezing coefficient (mm/s)
    real(8) :: frac_snowliq              ! r, maximum fraction of liquid in snow (-)

    ! flow field
    real(8) :: zmeas    
    real(8) :: zground   
    real(8) :: zo_ground 
  
  end type spafhy_para_type
  
  ! canopy interception
  ! real(8) :: kv = 0.4  ! von Karman constant (-) ! moved to spafhy_mod
  ! real(8) :: beta_aero=285.0   ! s/m, from Campbell & Norman eq. (7.33) x 42.0 molm-3 ! moved to spafhy_mod

  ! LAI is annual maximum LAI and for gridded simulations are input from GisData!
  ! keys must be 'LAI_ + key in spec_para
  !real(r8)  :: LAI_conif
  !real(r8)  :: LAI_decid

  ! canopy conductance                     
  ! real(r8), parameter :: kp =0.6         ! Hui: This overlaps with p-hydro

  ! soil evaporation
  ! CFT parameters (not needed in spafhy, replaced by p-hydro)
  ! 'amax': 10.0, # maximum photosynthetic rate (umolm-2(leaf)s-1)
  ! 'g1': 2.1, # stomatal parameter
  ! 'q50': 50.0, # light response parameter (Wm-2)
  ! 'lai_cycle': False,

  ! phenology (not needed, will be replaced by penology module in svm)

  ! 'smax': 18.5, # degC
  ! 'tau': 13.0,  # days
  ! 'xo': -4.0, # degC
  ! 'fmin': 0.05, # minimum photosynthetic capacity in winter (-)
                           
  ! annual cycle of leaf-area in deciduous trees
  !real(r8) :: lai_decid_min = 0.1     ! minimum relative LAI (-)
  !real(r8) :: ddo= 45.0               ! degree-days for bud-burst (5degC threshold)
  !real(r8) :: ddur= 23.0              ! duration of leaf development (days)
  !real(r8) :: sdl= 9.0                ! daylength for senescence start (h)
  !real(r8) :: sdur= 30.0              ! duration of leaf senescence (days),     


  ! initial states: rootzone and toplayer soil saturation ratio [-] and pond storage [m]
  !real(8) :: rootzone_sat   ! root zone saturation ratio (-)
  !real(8) :: soilcode    = -1   ! site-specific values

contains

  subroutine readsoilhydro_namelist(spafhy_para)
  
    type(spafhy_para_type), intent(inout)    :: spafhy_para

    ! local parameters
    real(8) :: soil_depth    ! root zone depth (m)
    real(8) :: max_poros     ! [m3 m-3], porosity
    real(8) :: fc            ! [m3 m-3], field capacity. For consistency, must be computed from soil water retention curve at Psi= xx KPa
    real(8) :: wp            ! [m3 m-3], wilting point. For consistency, must be computed from soil water retention curve at Psi=xx KPa
    real(8) :: ksat          ! [m s-1], saturated hyd, conductivity
    real(8) :: n_van       !Launiainen et al. 2022: C1-5: 1.12, 1.14, 1.07, 1.27, 1.18    
    real(8) :: watres      !Launiainen et al. 2022: C1-5: 0.0
    real(8) :: alpha_van   !Launiainen et al. 2022: C1-5: 4.45, 5.92, 2.02, 4.49, 3.35
    real(8) :: watsat      !Launiainen et al. 2022: C1-5: 0.75, 0.68, 0.46, 0.47, 0.54  
    real(8) :: org_depth     ! depth of organic top layer (m)
    real(8) :: org_poros      ! porosity (-)
    real(8) :: org_fc         ! field capacity (-)
    real(8) :: org_sat        ! organic top layer saturation ratio (-)
    !real(8) :: org_rw        ! critical vol. moisture content (-) for decreasing phase in Ef
    real(8) :: maxpond        ! max ponding allowed (m)
    
    ! Canopy water parameters
    real(8) :: wmax     ! storage capacity for rain (mm/LAI), Hui: this is too much compared to CTSM
    real(8) :: wmaxsnow ! storage capacity for snow (mm/LAI), Hui: this is reasonable
    real(8) :: hc       ! canopy height (m)
    real(8) :: cf       ! canopy closure fraction (-)
    real(8) :: w_leaf   ! leaf length scale (m)
    real(8) :: rw       ! critical value for REW (-),
    real(8) :: rwmin    ! minimum relative conductance (-)
    real(8) :: gsoil    ! soil surface conductance if soil is fully wet (m/s)
    real(8) :: kmelt         ! melt coefficient in open (mm/s)
    real(8) :: kfreeze       ! freezing coefficient (mm/s)
    real(8) :: frac_snowliq              ! r, maximum fraction of liquid in snow (-)
    real(8) :: zmeas    
    real(8) :: zground   
    real(8) :: zo_ground 

    logical :: old
    integer :: readerror
    integer,parameter :: unitsoilhydro=3

    namelist /soilhydro_namelist/ &
      soil_depth, &
      max_poros, &
      fc, &
      wp, &
      ksat, &
      org_depth, & 
      org_poros, &
      org_fc, & 
      maxpond, &
      org_sat, &   
      n_van, &     
      watres, &
      alpha_van, &
      watsat, &
      wmax, &
      wmaxsnow, &
      hc ,&
      cf, &
      w_leaf, &
      rw, &
      rwmin, &
      gsoil, &
      kmelt, &
      kfreeze, &
      zmeas, &
      zground, &
      zo_ground


    old=.false.
  ! Presetting namelist command
    soil_depth=0.4
    max_poros=0.46           ! should be equivalent to watsat here.
    fc=0.36                   ! based on C3 in Launiainen et al. 2022 ! Must be computed from water-retention curve
    wp=0.22                   ! Must be computed from water-retention curve
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
  
    wmax = 0.5      ! storage capacity for rain (mm/LAI), default: 1.5 too high? Hui: this is too much compared to CTSM
    wmaxsnow = 4.5  ! storage capacity for snow (mm/LAI), Hui: this is reasonable
    hc = 0.6         ! canopy height (m)
    !cf = 0.6          ! canopy closure fraction (-) LET's omit for simplicity!!
    w_leaf=0.01       !leaf length scale (m)
    ! canopy conductance                     
    ! real(r8), parameter :: kp =0.6         ! canopy light attenuation parameter (-) Hui: This overlaps with p-hydro
    rw =0.20          ! critical value for REW (-),
    rwmin=0.02        ! minimum relative conductance (-)
    ! soil evaporation
    gsoil=1e-2              ! soil surface conductance if soil is fully wet (m/s)
    ! degree-day snow model
    kmelt   = 2.8934e-05    ! melt coefficient in open (mm/s)
    kfreeze = 5.79e-6       ! freezing coefficient (mm/s)
    frac_snowliq = 0.05          ! maximum fraction of liquid in snow (-)
    ! flow field
    zmeas     = 2.0
    zground   = 0.5
    zo_ground = 0.01 

  ! Reading namelist
  !  open(unitsoilhydro, file='./soilhydro_namelist', status='old', form='formatted', err=999)
  !  read(unitsoilhydro, soilhydro_namelist, iostat=readerror)
  !  close(unitsoilhydro)

!999 write(*,*) ' #### MODEL ERROR! FILE "soilhydro_namelist"    #### '
!    write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
!    stop

     spafhy_para%soil_depth=soil_depth
     spafhy_para%max_poros=max_poros
     spafhy_para%fc=fc
     spafhy_para%wp=wp
     spafhy_para%ksat=ksat
     spafhy_para%org_depth=org_depth
     spafhy_para%org_poros=org_poros
     spafhy_para%org_fc=org_fc
     spafhy_para%org_sat=org_sat
     spafhy_para%maxpond=maxpond
     spafhy_para%n_van=n_van
     spafhy_para%watres=watres
     spafhy_para%alpha_van=alpha_van
     spafhy_para%watsat=watsat
     
     spafhy_para%wmax=wmax
     spafhy_para%wmaxsnow=wmaxsnow
     spafhy_para%hc=hc
     spafhy_para%w_leaf=w_leaf
     spafhy_para%rw=rw
     spafhy_para%rwmin=rwmin
     spafhy_para%gsoil=gsoil
     spafhy_para%kmelt=kmelt
     spafhy_para%kfreeze=kfreeze
     spafhy_para%frac_snowliq=frac_snowliq
     spafhy_para%zmeas=zmeas
     spafhy_para%zground=zground
     spafhy_para%zo_ground=zo_ground

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