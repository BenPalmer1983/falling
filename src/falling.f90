PROGRAM falling
! University of Birmingham
! Ben Palmer
!!
! Internal units:
! - energy, eV
! - length, Angstrom
! - forces, ev/Angstrom
! - stress, ev/ang^3
! - mass, atomic mass A
! - acceleration, ang/ns^2
! - time, ns  (nano seconds)
!
! Printed units:
! - energy, eV
! - length, Angstrom
! - forces, ev/Ang
! - stress, GPa
!
! Calculations:
! F = m a
! m = F (1/a)
! a = F/m        a [ang/ns^2] = 9.6485E9 * (F[ev/ang]/m[amu])
!
! Setup Modules
  Use kinds          ! libBP  kinds.f90
  Use programMod     ! libBP  programMod.f90
! Program Selection
  Use fallingMod

! Force declaration of all variables
  Implicit None
! Include MPI header
  Include 'mpif.h'
! Variables
  Integer(kind=StandardInteger) :: error
!======================================
! Init MPI
  Call MPI_Init(error)
!======================================
!--------------------------------------------------------------------------------------------------------
! Program Starts !
!----------------!
  Call ProgramInit()  ! libBP  programMod.f90
  Call ProgramSelect()  ! eampa.selection.f90

!======================================
! Finalise MPI
  Call MPI_Finalize(error)
!======================================

  Contains

! ------------------------------------------------------------
!               Sub Program Selection
! ------------------------------------------------------------

  Include "falling.selection.f90"



! --------------------------------------------------------------


End Program falling
