! --------------------------------------------------------------!
! Static Caculations module
! staticCalcsTypes, staticCalcs
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Calculates the forces between atoms
! Calculates the energies of atoms/total energy
! Calculates stresses
!
! ----------------------------------------
! Updated: 21st May 2016
! ----------------------------------------

! ------------------------------------------------------------
!     DYNAMICCALCS: Types
! ------------------------------------------------------------
  Include "dynamicCalcs._types.f90"
! ------------------------------------------------------------

Module dynamicCalcs
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use mpi
  Use kinds
  Use rng
  Use scienceFunctions
  Use strings
  Use constants
  Use units
  Use printModTypes
  Use printMod
  Use matrix
  Use mpiSubsTypes
  Use mpiSubs
  Use basicMaths
  Use keysMod
  Use linearAlgebra
  Use coordFunctions
  Use fitting
  Use splinesFitting
  Use simulatedAnnealing
  Use geomTypes
  Use geom
  Use potentialsTypes
  Use potentials
  Use staticCalcsTypes
  Use staticCalcs
  Use dynamicCalcsTypes
! Force declaration of all variables
  Implicit None
  Private
! ---- Variables
!  Public :: nl
! ---- Subroutines
  Public :: runDynamics
  Public :: make1D_coords
  Public :: make1D_nl
  Public :: calcEF_Nodes
  Public :: printEF_Nodes
  Public :: run1D_SaOpt
  Public :: run1D_GD
  Public :: run1D_NewtonOpt
  Public :: run1D_CG



!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! -----------------------------------------------
!        Module Subroutines
!
! -----------------------------------------------


! ------------------------------------------------------------
!      Dynamics subroutines for 3D Molecular Dynamics
! ------------------------------------------------------------

  Include "dynamicCalcs.md_3d.f90"


! ------------------------------------------------------------
!      Dynamics subroutines for 3D Molecular Dynamics
! ------------------------------------------------------------

  Include "dynamicCalcs.nodes_1d.f90"



End Module dynamicCalcs


















!---------------------------------
