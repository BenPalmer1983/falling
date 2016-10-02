! -------------------------------------------------
!  Include File:   Materials Functions
!
! -------------------------------------------------


  Function F_BirchMurn(volume,coefficients) RESULT (energy)
! Calculate energy from volume using Murnaghan EoS
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal) :: volume, energy, eta
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients    ! E0, V0, B0, B'0
    Real(kind=DoubleReal) :: E, V, B, BP
! Shorten var names
    E = coefficients(1)
    V = coefficients(2)
    B = coefficients(3)
    BP = coefficients(4)
! Calculate energy
    eta = ((1.0D0*volume)/(1.0D0*V))**(1.0D0/3.0D0)
    energy = E+(9.0D0/16.0D0)*(B*V)*&
    ((eta**2-1.0D0)**2)*(6.0D0+BP*(eta**2-1.0D0)-4.0D0*eta**2)
! Rearranged:
! energy = E+(9.0D0/16.0D0)*(B)*(&
!  V**(-1.0D0)*volume**2.0D0*(BP-4.0D0)&
!  +V**(-1.0D0/3.0D0)*volume**(4.0D0/3.0D0)*(14.0D0-3.0D0*BP)&
!  +V**(1.0D0/3.0D0)*volume**(2.0D0/3.0D0)*(3.0D0*BP-16.0D0)&
!  +V*(6.0D0-BP))
  End Function F_BirchMurn














! -------------------------------------------------
