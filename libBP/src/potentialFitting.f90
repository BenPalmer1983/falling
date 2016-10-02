! --------------------------------------------------------------!
! Module to fit potentials toe DFT Data
! potentialFittingTypes, potentialFitting
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Module Description
!
! ----------------------------------------
! Updated:
! ----------------------------------------

Module potentialFittingTypes
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
  Public :: oPotentialFitting
  Public :: oPfInputFiles
!  Public :: oFittingReference
!------------------------------------------------------------------------------
! Defined Types
!------------------------------------------------------------------------------
  Type :: oPfInputFiles
    Character(Len=128) :: filePath
    Integer :: fileType
  End Type oPfInputFiles

  Type :: oPfCoords
    Character(Len=16) :: label
    Real(kind=DoubleReal), Dimension(1:3) :: xyz
    Real(kind=DoubleReal), Dimension(1:3) :: forces
  End Type oPfCoords

  Type :: oPfConfig
    Character(Len=128) :: sourceFile
    Real(kind=DoubleReal) :: aLat     ! Alat of entire config
    Real(kind=DoubleReal) :: energy
    Real(kind=DoubleReal) :: epa
    Integer(kind=StandardInteger) :: coordLength
    Type(oPfCoords), Allocatable, Dimension(:) :: coords
  End Type oPfConfig

  Type :: oPotentialFitting
! Directories
    Character(Len=128) :: outputDirectory
! Input file, record of coords file
    Character(Len=128) :: inputFile
    Character(Len=128) :: recordFile
! Potential
    Character(Len=128) :: potentialFile
! Config Count
    Integer(kind=StandardInteger) :: configCount
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: configMap   ! link configs array to n=1,configCount
! Config files/data
    Integer(kind=StandardInteger) :: cFileCount
    Integer(kind=StandardInteger) :: dftFileCount
    Integer(kind=StandardInteger) :: fileCount  ! File count
    Type(oPfInputFiles), Allocatable, Dimension(:) :: fileList
    Type(oPfConfig), Allocatable, Dimension(:) :: configs
! Calc info
    Integer(kind=StandardInteger) :: calcCount = 0
    Integer(kind=StandardInteger) :: calcCount_TotalEFS = 0
    Logical :: calcFile = .true.
  End Type oPotentialFitting






! Old, delete?
  Type :: oFittingReference
    Real(kind=DoubleReal) :: totalEnergy
    Real(kind=DoubleReal) :: aLat
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: forces
  End Type oFittingReference


End Module potentialFittingTypes


Module potentialFitting
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use mpi
  Use kinds
  Use strings
  Use generalTypes
  Use general
  Use units
  Use mpiSubsTypes
  Use mpiSubs
  Use env
  Use units
  Use printMod
  Use geomTypes
  Use geom
  Use potentialsTypes
  Use potentials
  Use staticCalcsTypes
  Use staticCalcs
  Use qeTypes
  Use qe
  Use potentialFittingTypes
! Force declaration of all variables
  Implicit None
  Private
! ---- Variables
!  Public :: nl
! ---- Subroutines
  Public :: startPotentialFit
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

  Subroutine startPotentialFit()
! Start potential fitting
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
! Vars:  Private
    Type(oPotentialFitting) :: pfObj
    Type(coordsType), Allocatable, Dimension(:) :: pfCoords
    Type(nlType), Allocatable, Dimension(:) :: pfNl
    Type(potentialType) :: pfPotential




!-------------------------------------------------------------------------
! Init Page
    Call initPage(mainPage)
! Title
    Call addLinePage("Potential Fitting","T")

!---------------------------------------------
! Load files
    Call pfLoadFiles(pfObj)
! Make configs
    Call pfMakeNL(pfObj, pfCoords, pfNl)
! Load potential
    Call pfLoadPotential(pfObj, pfPotential)
! Assign atom IDs
    Call addLinePage("Assign IDs")
    Call atomLabelIDs(pfPotential, pfCoords)

! Print potential summary
    Call printPotentialSummary(pfPotential)      ! potentials.f90

! Run calculations
    Call pfRunCalc(pfObj, pfCoords, pfNl, pfPotential)







!Call addLinePage("Input file: "//trim(cmdArgs%args(1)))

! Count configurations
    !Call countInputConfigs(trim(cmdArgs%args(1)), configCount)
! Allocate arrays
    !Allocate(coords(1:configCount))
    !Allocate(refData(1:configCount))
! Load input files
    !Call loadInputFiles(trim(cmdArgs%args(1)),coords,refData)

    !Call printCoords(coords, 1)









! Deallocate arrays
    Deallocate(pfNl)
! Output Page
    Call printPage(mainPage)
!----
  End Subroutine startPotentialFit



! ------------------------------------------------------------
!               Load Files
! ------------------------------------------------------------

  Include "potentialFitting.loadFiles.f90"

! ------------------------------------------------------------
!               Make Configs
! ------------------------------------------------------------

  Include "potentialFitting.makeConfigs.f90"

! ------------------------------------------------------------
!               Calculations
! ------------------------------------------------------------

  Include "potentialFitting.calculations.f90"


! -----------------------------------------------
!        Module Functions
!
! -----------------------------------------------




End Module potentialFitting

























!-----------------------------------
