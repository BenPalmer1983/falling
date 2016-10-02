! -------------------------------------------------
!  Include File:   Classical Mechanics
!
! -------------------------------------------------

  Function F_forceMA (mass, acceleration) RESULT (force)
! F = MA  SI units         N = kg m s^-2
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal) :: mass, acceleration
! Vars:  Out
    Real(kind=DoubleReal) :: force
! Calculate
    force = mass * acceleration
  End Function F_forceMA


!---------------------------
! MD units eV, ang, ps, amu
!---------------------------

  Function F_forceMA_MD (mass, acceleration) RESULT (force)
! F = MA  MD units
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal) :: mass, acceleration
! Vars:  Out
    Real(kind=DoubleReal) :: force
! Calculate
    force = (mass * acceleration) * 1.03643D-4
  End Function F_forceMA_MD

  Function F_accelerationFM_MD (force, mass) RESULT (acceleration)
! A = F/M  MD units
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal) :: mass, acceleration
! Vars:  Out
    Real(kind=DoubleReal) :: force
! Calculate
    !force = mass * acceleration
    acceleration = (force/mass)*9.64853D3
  End Function F_accelerationFM_MD














! -------------------------------------------------
