MODULE readsoilpara_mod

implicit none
  type, public :: soilwater_type
    real(8)  :: PondSto
    real(8)  :: inflow  ! [m] - total inflow to root zone
    real(8)  :: roff    ! [m] - surface runoff
    real(8)  :: drain   ! [m] - drainage from root zone
    real(8)  :: Interc     ! interception of top layer (mm)
    real(8)  :: WatSto
    real(8)  :: MaxWatSto
    real(8)  :: WatStoTop
    real(8)  :: MaxStoTop
    real(8)  :: Wliq
    real(8)  :: Wliq_top
    real(8)  :: Sat
    real(8)  :: Rew
    real(8)  :: mbe     ! [m] - mass balance error
  end type soilwater_type

  type, public :: canopywater_type
    real(8):: h2o
    real(8):: Trfall     ! throughfall to snow / soil surface (mm)
    real(8):: Interc     ! interception of canopy (mm)
    real(8):: Evap       ! evaporation / sublimation from canopy store (mm)
    real(8):: Unload     ! undloading from canopy storage (mm)    
    real(8):: Efloor     ! forest floor evaporation rate (mm s-1)
    real(8):: MBE        ! mass balance error (mm)      
    real(8):: ET         ! total evapo-transpiration
    !real(8):: Transpi    ! transpiration rate (mm s-1)  
  end type canopywater_type

  type, public :: snowwater_type
    real(8):: swe    
    real(8):: SWEi    
    real(8):: SWEl   
    real(8):: PotInf     ! potential infiltration to soil profile (mm)            
  end type snowwater_type

! Soil hydraulic properties for spafhy
  real(8) :: soil_depth    ! root zone depth (m)
  real(8) :: max_poros     ! porosity (-)
  real(8) :: fc            ! field capacity (-)
  real(8) :: wp            ! wilting point (-)
  real(8) :: ksat        ! conductivity at saturation point
  real(8) :: beta           !  term for soil evaporation resistance (Wliq/FC) [-]

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
    max_poros =0.46           ! should be equivalent to watsat here.
    fc=0.36                   ! based on C3 in Launiainen et al. 2022
    wp=0.22
    ksat=2.0e-6
    beta=4.7            ! default 
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


  

END MODULE readsoilpara_mod

