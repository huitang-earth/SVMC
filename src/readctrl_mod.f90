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

  !Input files/settings
  
  integer :: start_date_day, start_date_hour, time_step, &
             end_date_day, end_date_hour
  real    :: time_step
  character(len=256)  :: output_directory, output_filename


  subroutine readctrl_namelist

    use par_mod
    use com_mod
    use dust_mod   
    implicit none

    logical :: old
    integer :: readerror
    integer,parameter :: unitcommand=1

    old=.false.
:
    namelist /ctrl_namelist/ &
    start_date_day, &
    start_date_hour, &
    time_step, &
    end_date_day, &
    end_date_hour, &
    output_directory, &
    output_filename
    
    ! Presetting namelist command
    start_date_day  =20190301
    start_date_hour =000000
    time_step       =1               ! hour
    end_date_day    =20191001
    end_date_hour   =000000
    output_directory='/cluster/work/users/ovewh/'
    output_filename ='test.nc' 

    ! Reading namelist
    open(unitcommand, file='./ctrl_namelist', status='old', form='formatted', err=999)
    read(unitcommand, ctrl_namelist, iostat=readerror)
    close(unitcommand)

    if (time_step.eq.0) then
      write(*,*) ' #### MODEL ERROR! TIME STEP MUST    #### '
      write(*,*) ' #### NOT BE ZERO                             #### '
      write(*,*) ' #### CHANGE INPUT IN FILE COMMAND.           #### '
      stop
    endif

    return

999   write(*,*) ' #### MODEL ERROR! FILE "ctrl_namelist"    #### '
    write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
    stop

  end subroutine readctrl_namelist
                   
end module readctrl_mod
