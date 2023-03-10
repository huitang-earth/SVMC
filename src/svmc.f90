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

  use netcdf             ! library for processing netcdf files
  use readpara_mod       ! module for reading parameter files in ASCII
  use readclim_mod       ! module for reading reading meteorological forcing data
  use readsoil_mod       ! module for reading soil properties (shared with yasso?)
  use phydro_mod         ! module for p-hydro
  use alloc_mod          ! module for carbon allocation and yield
  !use yasso20           ! placeholder for soil decomposition model, which will provide hr   
  !use soilwater_mod     ! placeholder for soil water bucket model, which will provide psi_soil for p-hydro

  implicit none

  !Loop variables
  !***********************************
  integer   :: i, j, k, t

  !***********************************

  !Model variables
  !***********************************
  integer   ::

  ! Input variables for p-hydro (the definition is from rpmodel.R)
  real     ::    tc        ! Air temperature (tc), degrees C
  real     ::    ppfd      ! Photosynthetic photon flux density (mol m-2 d-1) (incoming solar radiation from forcing data?)
  real     ::    vpd       ! Vapour pressure deficit (Pa) (will be calculated using pressure & humidity)
  real     ::    co2       ! Atmospheric CO2 concentration (ppm)
  real     ::    elv       ! Elevation above sea-level (m.a.s.l.) (not needed if we have surface pressure!)
  real     ::    fapar     ! Fraction of absorbed photosynthetically active radiation (unitless) (will be calculated using LAI)
  real     ::    kphio     ! Apparent quantum yield efficiency (unitless).
  real     ::    psi_soil  ! soil water potential (Mpa)
  real     ::    rdark = 0 !
  real     ::    par_plant !  A list of plant hydraulic parameters (will be defined in readpara_mod.f90)
  real     ::   par_cost = NULL  ! A list of cost parameters
  character(len=200)  :: opt_hypothesis=''   ! character, Either "Lc" or "PM"


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
