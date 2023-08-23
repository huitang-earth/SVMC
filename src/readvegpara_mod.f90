MODULE readvegpara_mod

  implicit none

  ! !PUBLIC TYPES:
  ! Define type for plant hydraulic parameters
  
  !? Do we need types or just have variables separately
  type, public :: par_plant_type
    !!**** Parameters for P-hydro model
    ! Plant hydraulic parameters
    real(r8) :: conductivity     ! Leaf conductivity (m) (for stem, this could be Ks*HV/Height)
    real(r8) :: psi50             ! Leaf P50 (Mpa)
    real(r8) :: b                     ! Slope of leaf vulnerability curve 
  end type par_plant_type

  type, public :: par_cost_type
    ! A list of cost parameters
    real(r8) :: alpha  !cost of Jmax
    real(r8) :: gamma    !cost of hydraulic repair
  end type par_cost_type

  character(len=256), public  :: opt_hypothesis        ! character, Either "Lc" or "PM"      

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

  type, public :: optimizer_type
    real(r8) :: logjmax  
    real(r8) :: dpsi                
  end optimizer_type

  ! vegetation status (p-hydro)
  real(r8), parameter :: kv = 0.4  ! von Karman constant (-)
  real(r8)            ::  beta   ! s/m, from Campbell & Norman eq. (7.33) x 42.0 molm-3

  real(r8)            :: kphio=0.087182   ! Apparent quantum yield efficiency (unitless).
  real(r8)            :: k=0.5            ! dimensionless constant, assigned a generic value of 0.5, Beer's law.

  ! Canopy water parameters (Spafhy)
  ! canopy interception
  real(r8) :: wmax     ! storage capacity for rain (mm/LAI), Hui: this is too much compared to CTSM
  real(r8) :: wmaxsnow ! storage capacity for snow (mm/LAI), Hui: this is reasonable

  ! LAI is annual maximum LAI and for gridded simulations are input from GisData!
  ! keys must be 'LAI_ + key in spec_para
  !real(r8)  :: LAI_conif
  !real(r8)  :: LAI_decid
  real(r8)  :: hc         ! canopy height (m)
  real(r8)  :: cf         ! canopy closure fraction (-)

  ! canopy conductance                     
  ! real(r8), parameter :: kp =0.6         ! canopy light attenuation parameter (-) Hui: This overlaps with p-hydro
  real(r8) :: rw         ! critical value for REW (-),
  real(r8) :: rwmin        ! minimum relative conductance (-)
  ! soil evaporation
  real(r8) :: gsoil              ! soil surface conductance if soil is fully wet (m/s)

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

  ! degree-day snow model
  real(r8) :: kmelt         ! melt coefficient in open (mm/s)
  real(r8) :: kfreeze       ! freezing coefficient (mm/s)
  real(r8) :: r             ! maximum fraction of liquid in snow (-)

  ! flow field
  real(r8) :: zmeas    
  real(r8) :: zground   
  real(r8) :: zo_ground 

  contains

  subroutine par_init(this)

    class(canopystate_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

    if ( this%leaf_mr_vcm == spval ) then
       call endrun(msg="ERROR canopystate Init called before ReadNML"//errmsg(sourcefile, __LINE__))
    end if

  end subroutine Init

  subroutine readvegpara_namelist
    implicit none
    
    logical :: old
    integer :: readerror
    integer,parameter :: unitvegpara=2

    old=.false.

    namelist /veg_namelist/ &
      num_pft,
      pft_type,
      conductivity, &
      psi50, &
      b, &
      alpha, &
      gamma, &
      opt_hypothesis, &
      beta, &
      wmax, &
      wmaxsnow, &
      hc ,&
      cf, &
      rw, &
      rwmin, &
      gsoil, &
      kmelt, &
      kfreeze, &
      r, &
      zmeas, &
      zground, &
      zo_ground
 
    ! Presetting namelist command
    !--- Parameters for P-hydro model
    ! Plant hydraulic parameters
    conductivity=3e-17     ! Leaf conductivity (m) (for stem, this could be Ks*HV/Height)
    psi50 = -2             ! Leaf P50 (Mpa)
    b=2                     ! Slope of leaf vulnerability curve 

    ! A list of cost parameters
    alpha=0.1  !cost of Jmax
    gamma=1     !cost of hydraulic repair

    opt_hypothesis= 'Lc'          ! character, Either "Lc" or "PM"      

    ! vegetation status (p-hydro)
    beta = 285.0   ! s/m, from Campbell & Norman eq. (7.33) x 42.0 molm-3

    ! Canopy water parameters (Spafhy)
    ! canopy interception
    wmax = 1.5      ! storage capacity for rain (mm/LAI), Hui: this is too much compared to CTSM
    wmaxsnow = 4.5, ! storage capacity for snow (mm/LAI), Hui: this is reasonable

    ! LAI is annual maximum LAI and for gridded simulations are input from GisData!
    ! keys must be 'LAI_ + key in spec_para
    hc = 0.6         ! canopy height (m)
    cf = 0.6          ! canopy closure fraction (-)

    ! canopy conductance                     
    ! real(r8), parameter :: kp =0.6         ! canopy light attenuation parameter (-) Hui: This overlaps with p-hydro
    rw =0.20          ! critical value for REW (-),
    rwmin=0.02        ! minimum relative conductance (-)
  
    ! soil evaporation
    gsoil=1e-2              ! soil surface conductance if soil is fully wet (m/s)

    ! degree-day snow model
    kmelt   = 2.8934e-05    ! melt coefficient in open (mm/s)
    kfreeze = 5.79e-6       ! freezing coefficient (mm/s)
    r       = 0.05          ! maximum fraction of liquid in snow (-)

    ! flow field
    zmeas     = 2.0
    zground   = 0.5
    zo_ground = 0.01 
  
    ! Reading namelist
    open(unitvegpara, file='./veg_namelist', status='old', form='formatted', err=999)
    read(unitvegpara,veg_namelist,iostat=readerror)
    close(unitvegpara)

999   write(*,*) ' #### MODEL ERROR! FILE "veg_namelist"    #### '
    write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
    stop

  end subroutine readvegpara_namelist

                                     
end module readvegpara_mod
