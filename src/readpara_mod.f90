MODULE readpara_mod

  implicit none

  ! !PUBLIC TYPES:
  ! Define type for plant hydraulic parameters
  
  !? Do we need types or just have variables separately
  type, public :: par_plant_type

     !!**** Parameters for P-hydro model
     ! Plant hydraulic parameters
     real(r8) :: conductivity=3e-17,     ! Leaf conductivity (m) (for stem, this could be Ks*HV/Height)
     real(r8) :: psi50 = -2,             ! Leaf P50 (Mpa)
     real(r8) :: b=2                     ! Slope of leaf vulnerability curve 

  end type par_plant_type

  type, public :: par_cost_type

    ! A list of cost parameters
    real(r8) :: alpha=0.1,  !cost of Jmax
    real(r8) :: gamma=1     !cost of hydraulic repair

  end type par_cost_type    

  character(len=256), public  :: opt_hypothesis= ''          ! character, Either "Lc" or "PM"      

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


end module readpara_mod
