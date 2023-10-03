MODULE readctrl_mod

  !*****************************************************************************
  !                                                                            *
  !     This routine reads the user specifications for the current model run.  *
  !                                                                            *
  !     Author: H. Tang                                                        *
  !*****************************************************************************
  !                                                                            *
  !                                                                            *
  ! Variables:                                                                 *
  ! unitcommand          unit connected to file COMMAND                        *
  ! start_date_day       start date of simulation                              *
  ! start_date_hour      start time of simulation                              *
  ! time_step            output time step                                      *
  ! end_date_day         end date of simulation                                *
  ! end_date_hour        end time of simulation                                *
  ! output_directory     directory of output                                   *
  ! output_filename      file name of output                                   *
  !                                                                            *
  !*****************************************************************************

  implicit none

  public :: readctrl_namelist
  
  !Input files/settings

  integer,parameter :: dp=selected_real_kind(P=15)
  
  integer :: start_date_day, start_date_hour, &
             end_date_day, end_date_hour
  integer :: num_sites              ! This can be set by reading input file
  real(8), dimension(1)    :: lat_sites, lon_sites   ! This can be set by reading input file 
  real    :: time_step, time_step_output
  character(len=256)  :: output_filename_day, output_filename_hr, input_climfile, input_laifile
  logical :: obs_lai, obs_soilmoist

contains
  !------------------------------------------------------
  subroutine readctrl_namelist
    logical :: old
    integer :: readerror
    integer,parameter :: unitcommand=1

    namelist /ctrl_namelist/ &
    start_date_day, &
    start_date_hour, &
    time_step, &
    end_date_day, &
    end_date_hour, &
    num_sites, & 
    input_climfile, &
    input_laifile, &
    time_step_output, &
    output_filename_hr, &
    obs_lai, &
    obs_soilmoist

    old=.false.
    
    ! Presetting namelist command
    start_date_day  =20210601 !20210101
    start_date_hour =000000
    time_step       =1               ! hours
    end_date_day    =20210615 !20211231
    end_date_hour   =000000
    num_sites       =1                 ! number of sites, should also be read from input file?
    ! lon_sites     = (/ /)            ! longitude of sites (not needed, can be well defined input file)
    ! lat_sites     = (/ /)            ! latitude of sites (not needed, can be well defined by input file) 
    input_climfile      ='../data/FieldObs_Qvidja.2021.hr.nc'         
    input_laifile       ='../data/FieldObs_Qvidja.2021.lai.nc'
    time_step_output=1.0
    output_filename_day ='../data/test.nc' 
    output_filename_hr ='../data/test_hr_2021_nolai_default_alpha0.08_gs_spafhy.nc' 
    obs_lai=.true.
    obs_soilmoist=.true.

    ! Reading namelist
    !open(unitcommand, file='./ctrl_namelist', status='old', form='formatted', err=999)
    !read(unitcommand, ctrl_namelist, iostat=readerror)
    !close(unitcommand)

    if (time_step.eq.0) then
      write(*,*) ' #### MODEL ERROR! TIME STEP MUST    #### '
      write(*,*) ' #### NOT BE ZERO                             #### '
      write(*,*) ' #### CHANGE INPUT IN FILE COMMAND.           #### '
      stop
    endif

    return

!999   write(*,*) ' #### MODEL ERROR! FILE "ctrl_namelist"    #### '
!    write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
!    stop

  end subroutine readctrl_namelist
                   
end module readctrl_mod
