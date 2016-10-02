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

Module staticCalcsTypes
! Setup Modules
  Use kinds
! Force declaration of all variables
  Implicit None
! Vars:  Module Parameters
  Integer(kind=StandardInteger), Parameter :: p_bpPoints = 5  ! Must be odd
! Make private
  Private
! Public Variables and Parameters
  Public :: p_bpPoints
! Public derived data types
  Public :: oCalcSettings, oBulkProperty

  Type :: oCalcSettings
    Integer(kind=StandardInteger) :: cKey = 0
    Integer(kind=StandardInteger) :: atomKey = 0   ! for calcEF
    Logical :: calcEnergy = .true.
    Logical :: calcForces = .true.
    Logical :: calcStress = .true.
    Logical :: resetForces = .false.
    Logical :: useRD_Moved = .false.
    Logical :: recalculateED = .false.
  End Type oCalcSettings

  Type :: oBulkProperty
! BCC
    Real(kind=DoubleReal) :: bccAlat_Estimate
    Real(kind=DoubleReal) :: bccAlat
    Real(kind=DoubleReal) :: bccV0
    Real(kind=DoubleReal) :: bccE0
    Real(kind=DoubleReal) :: bccB0
    Real(kind=DoubleReal) :: bccBp0
    Real(kind=DoubleReal) :: bccB0_GPA
    Real(kind=DoubleReal) :: bccC11
    Real(kind=DoubleReal) :: bccC11_GPA
    Real(kind=DoubleReal) :: bccC12
    Real(kind=DoubleReal) :: bccC12_GPA
    Real(kind=DoubleReal) :: bccC44
    Real(kind=DoubleReal) :: bccC44_GPA
    Real(kind=DoubleReal), Dimension(1:15,1:2) :: bccEosPoints
! FCC
    Real(kind=DoubleReal) :: fccAlat_Estimate
    Real(kind=DoubleReal) :: fccAlat
    Real(kind=DoubleReal) :: fccV0
    Real(kind=DoubleReal) :: fccE0
    Real(kind=DoubleReal) :: fccB0
    Real(kind=DoubleReal) :: fccBp0
    Real(kind=DoubleReal) :: fccB0_GPA
    Real(kind=DoubleReal) :: fccC11
    Real(kind=DoubleReal) :: fccC11_GPA
    Real(kind=DoubleReal) :: fccC12
    Real(kind=DoubleReal) :: fccC12_GPA
    Real(kind=DoubleReal) :: fccC44
    Real(kind=DoubleReal) :: fccC44_GPA
    Real(kind=DoubleReal), Dimension(1:15,1:2) :: fccEosPoints
  End Type oBulkProperty



End Module staticCalcsTypes


Module staticCalcs
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use mpi
  Use kinds
  Use time
  Use strings
  Use constants
  Use scienceFunctions
  Use units
  Use printModTypes
  Use printMod
  Use matrix
  Use mpiSubsTypes
  Use mpiSubs
  Use basicMaths
  Use keysMod
  Use isotopesTypes
  Use isotopes
  Use rng
  Use linearAlgebra
  Use coordFunctions
  Use fitting
  Use splinesFitting
  Use simulatedAnnealingTypes
  Use simulatedAnnealing
  Use geom
  Use geomTypes
  Use potentialsTypes
  Use potentials
  Use staticCalcsTypes
! Force declaration of all variables
  Implicit None
  Private
! ---- Variables
!  Public :: nl
! ---- Subroutines
  Public :: atomLabelIDs
  Public :: printAtomLabelIDs
  Public :: calcEFS
  Public :: calcEF
  Public :: calcE
  Public :: optGeom
  Public :: calcBP
  Public :: nlPotentialKeys
  Public :: nlPotentialKeys_opt


!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! -----------------------------------------------
!        Module Subroutines
!
! -----------------------------------------------

! ----------------
! Prep
! ----------------


!------------------------------------------
! Atom Labels
!------------------------------------------
  Include "staticCalcs.atomLabels.f90"

!------------------------------------------
! Key arrays
!------------------------------------------
  Include "staticCalcs.keys.f90"

!------------------------------------------
! Energy. Force, Stress calculations
!------------------------------------------
  Include "staticCalcs.calcEFS.f90"
! calcEFS, calcEFS_MPI, calcEFS_Action

!------------------------------------------
! Energy. Force, Stress calculations
!------------------------------------------
Include "staticCalcs.calcEF.f90"
! calcEFS, calcEFS_MPI, calcEFS_Action

!------------------------------------------
! Energy only calculations
!------------------------------------------
  Include "staticCalcs.calcE.f90"
! calcE, calcE_MPI, calcE_Action

!------------------------------------------
! Bulk Properties
!------------------------------------------
  Include "staticCalcs.calcBP.f90"

!------------------------------------------
! Opt Geom
!------------------------------------------
  Include "staticCalcs.optGeom.f90"

!------------------------------------------
! Print
!------------------------------------------
  Include "staticCalcs.print.f90"

!------------------------------------------
! Elastic Constants
!------------------------------------------
!  Include "staticCalcs.calcElasticConstants.f90"



! -----------------------------------------------
!        Module Functions
!
! -----------------------------------------------













End Module staticCalcs
