MODULE readpara_mod

implicit none

!!**** Parameters for P-hydro model
! Plant hydraulic parameters
  conductivity=3e-17,     ! Leaf conductivity (m) (for stem, this could be Ks*HV/Height)
  psi50 = -2,             ! Leaf P50 (Mpa)
  b=2                     ! Slope of leaf vulnerability curve 

! A list of cost parameters
  alpha=0.1,  !cost of Jmax
  gamma=1     !cost of hydraulic repair



contains


