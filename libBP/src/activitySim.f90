! --------------------------------------------------------------!
! title
! moduleTypes, module
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Module Description
!
! ----------------------------------------
! Updated:
! ----------------------------------------

Module moduleTypes
! Setup Modules
  Use kinds
! Force declaration of all variables
  Implicit None
! Vars:  Module Parameters
!  Integer(kind=StandardInteger), Parameter :: p_ = 3
! Make private
  Private
! Public Variables and Parameters
!  Public :: p_
! Public derived types
!  Public ::

!  Type :: oType
!    Character(Len=64) :: c
!    Integer(kind=StandardInteger) :: i
!    Logical :: l
!    Real(kind=DoubleReal) :: d
!  End Type oType


End Module moduleTypes


Module module
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use mpi
  Use kinds
! Force declaration of all variables
  Implicit None
  Private
! ---- Variables
!  Public :: nl
! ---- Subroutines
!  Public ::
! ---- Functions
!  Public ::
! Interfaces
!  Interface iName
!    Module Procedure subA, subB
!  End Interface iName


!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! -----------------------------------------------
!        Module Subroutines
!
! -----------------------------------------------


!  Subroutine subA()
! Subroutine description
!    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
! Vars:  Private
!  End Subroutine startPotentialFit





! -----------------------------------------------
!        Module Functions
!
! -----------------------------------------------


!  Function subA()
! Function description
!    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
! Vars:  Private
!  End Function description


End Module module

























!-----------------------------------
