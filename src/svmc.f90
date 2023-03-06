!**********************************************************************
! Copyright 2023                                                      *
! FMI                                                                 *
!                                                                     *                                    *
! ***.****@fmi.fi                                                     *                                                              *
!                                                                     *
! This file is part of FMI-SVM.                                       *
!                                                                     *
! FMI-SVM is free software: you can redistribute it and/or modify     *
! it under the terms of the GNU General Public License as published by*
! the Free Software Foundation, either version 3 of the License, or   *
! (at your option) any later version.                                 *
!                                                                     *
! You should have received a copy of the GNU General Public License   *
! along with FMI-SVM.  If not, see <http://www.gnu.org/licenses/>.    *
!**********************************************************************

program SVMC

  use readpara_mod       ! module for reading parameter files in ASCII
  use readclim_mod       ! module for reading reading meteorological forcing data
  use readsoil_mod       ! module for reading soil properties (shared with yasso?)
  use phydro_mod         ! module for p-hydro
  use alloc_mod          ! module for carbon allocation and yield
  !use yasso20
  !use soilwater_mod

  implicit none

  !Loop variables
  !***********************************
  integer   :: i, j, k, t
  real     ::    
  !***********************************

  !Model variables
  !***********************************
  integer   ::
  real     ::    
  character(len=200)  ::
  real,dimension(:,:), allocatable
  logical,dimension(:,:), allocatable ::
  !*********************************** 
 
  !Call

  !Call photosynthetic module from p-hydro
  call ***()  

  !Call Yasso model
  !call mod5c20(theta,time,temp,prec,init,b,d,leac,xt,steadystate_pred)
  

end program SVMC
