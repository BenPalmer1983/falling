! --------------------------------------------------------------!
! Science Functions module
! scienceFunctions
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
!
!
! ----------------------------------------
! Updated: 6th July 2016
! ----------------------------------------

!Module _Types
! Setup Modules
!  Use kinds
! Force declaration of all variables
!  Implicit None
! Vars:  Module Parameters
!  Integer(kind=StandardInteger), Parameter :: p_ = 32
! Make private
!  Private
! Public Variables and Parameters
!  Public :: p_
! Public derived types
!  Public ::
!End Module _Types




Module scienceFunctions
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Setup Modules
  Use mpi
  Use kinds
! Force declaration of all variables
  Implicit None
! Vars:  Module scope parameters
! Vars:  Module scope variables
! Make private
  Private
! Public Variables and Parameters
! Public Subroutines and Functions
! Classical Mechanics
  Public :: F_forceMA
  Public :: F_forceMA_MD
  Public :: F_accelerationFM_MD
! Potential Functions
  Public :: F_ZBL
  Public :: F_ZblFull
  Public :: F_Morse
  Public :: F_MorseFull
! Materials Functions
  Public :: F_BirchMurn
! Misc Functions
  Public :: F_LoopTemp

!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! -----------------------------------------------
!        Module Functions
!
! -----------------------------------------------


!------------------------------------------
! Classical Mechanics
!------------------------------------------
  Include "scienceFunctions.classicalMechanics.f90"

!------------------------------------------
! Potential Functions
!------------------------------------------
  Include "scienceFunctions.potentialFunctions.f90"

!------------------------------------------
! Potential Functions
!------------------------------------------
  Include "scienceFunctions.materialsFunctions.f90"

!------------------------------------------
! Misc Functions
!------------------------------------------
  Include "scienceFunctions.miscFunctions.f90"




End Module scienceFunctions

!
