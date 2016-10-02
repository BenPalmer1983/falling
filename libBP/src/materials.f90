! --------------------------------------------------------------!
! Materials module
! materialsTypes, materials
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
!
!
! ----------------------------------------
! Updated: 26th June 2016
! ----------------------------------------


Module materialsTypes
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
  Public :: oMaterial
!------------------------------------------------------------------------------
! Defined Types
!------------------------------------------------------------------------------
  Type :: oM_Isotope
    Integer(kind=StandardInteger) :: isotopeP
    Integer(kind=StandardInteger) :: isotopeN
    Real(kind=DoubleReal) :: ratioByMass
    Real(kind=DoubleReal) :: percentageByMass
    Real(kind=DoubleReal) :: percentageByNumber
  End Type oM_Isotope

  Type :: oMaterial
    Integer(kind=StandardInteger) :: isotopeCount
    Type(oM_Isotope), Allocatable, Dimension(:) :: isotope
  End Type oMaterial

End Module materialsTypes


Module materials
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
! Force declaration of all variables
  Implicit None
! Make private
  Private
! Public
! Interfaces
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! ------------------------------------------------------------------------!
! Material Subroutines
! ------------------------------------------------------------------------!







End Module materials

!-----------------------------------------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------------------------------------
