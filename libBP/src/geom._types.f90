! ------------------------------------------------------------
!               GEOM: Types
! ------------------------------------------------------------

Module geomTypes
! Setup Modules
  Use kinds
! Force declaration of all variables
  Implicit None
! Vars:  Module Parameters
  Integer(kind=StandardInteger), Parameter :: p_confs = 10
  Integer(kind=StandardInteger), Parameter :: p_cMax = 4096        ! Coords max length
  Integer(kind=StandardInteger), Parameter :: p_cuMax = 512        ! Coords unit type max length
  Integer(kind=StandardInteger), Parameter :: p_nlMax = 50000      ! Neighbour list length - no longer used
! Make private
  Private
! Public Variables and Parameters
  Public :: p_confs, p_cuMax, p_cMax, p_nlMax
! Public derived types
  Public :: oCoord
  Public :: coordsUnitType, coordsType, oSubCell, nlType, oNLTempRD, nlType_Opt
  Public :: o1D_Nodes, o1D_NodesNL

! oCoord (label, labelID, xyz, force, velocity)
  Type :: oCoord
    Character(Len=16) :: label
    Integer(kind=StandardInteger) :: labelID
    Real(kind=DoubleReal), Dimension(1:3) :: xyz
    Real(kind=DoubleReal), Dimension(1:3) :: force
    Real(kind=DoubleReal), Dimension(1:3) :: velocity
  End Type oCoord

  Type :: oParticle
    Integer(kind=StandardInteger) :: charge   ! e
    Integer(kind=StandardInteger) :: mass     ! amu
  End Type oParticle

  Type :: coordsUnitType
    Integer(kind=StandardInteger) :: points = 0
    Integer(kind=StandardInteger) :: xCopy = 1
    Integer(kind=StandardInteger) :: yCopy = 1
    Integer(kind=StandardInteger) :: zCopy = 1
    Real(kind=DoubleReal) :: aLat = 1.0D0
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: unitCell
    Character(Len=16), Dimension(1:p_cuMax) :: label           ! Unit labels
    Integer(kind=StandardInteger), Dimension(1:p_cuMax) :: labelID
    Real(kind=DoubleReal), Dimension(1:p_cuMax,1:3) :: unitCoords   ! Unit coords - fractional
    Real(kind=DoubleReal), Dimension(1:p_cuMax,1:3) :: unitForces   ! Unit coords - fractional
! particle details
    Integer(kind=StandardInteger), Dimension(1:p_cuMax) :: charge = 1   ! e
    Real(kind=DoubleReal), Dimension(1:p_cuMax) :: mass = 1.0D0     ! amu
  End Type coordsUnitType

  Type :: coordsType
    Integer(kind=StandardInteger) :: length = 0
    Real(kind=DoubleReal) :: aLat = 0.0D0
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: unitCell
    Character(Len=16), Dimension(1:p_cMax) :: label           ! Unit labels
    Integer(kind=StandardInteger), Dimension(1:p_cMax) :: labelID
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: coords                ! coords (fractional)          Doesn't change
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: forces = 0.0D0        ! initial forces     ev/ang    Doesn't change
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: velocity = 0.0D0      ! initial velocity   ang/ps    Doesn't change
    Integer(kind=StandardInteger), Dimension(1:p_cMax) :: coordsFixed = 0           ! coords fixed
    Integer(kind=StandardInteger) :: IDcount
    Character(Len=16), Dimension(1:p_cuMax) :: atomIDs
    Integer(kind=StandardInteger) :: atomID_Count
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: electronDensity   ! Moved to nl object
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: atomEnergy        ! Moved to nl object
    Real(kind=DoubleReal) :: pairEnergy
    Real(kind=DoubleReal) :: embeddingEnergy
    Real(kind=DoubleReal) :: totalEnergy
! particle details
    Integer(kind=StandardInteger), Dimension(1:p_cMax) :: charge = 1   ! e
    Real(kind=DoubleReal), Dimension(1:p_cMax) :: mass = 1.0D0     ! amu
  End Type coordsType


  Type :: oSubCell
    Integer(kind=StandardInteger) :: width
    Integer(kind=StandardInteger) :: count
    Real(kind=DoubleReal) :: fracAlat
    Real(kind=DoubleReal) :: aLat
    Integer(kind=StandardInteger) :: maxAtoms
    Integer(kind=StandardInteger) :: cellCombinations
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: atomCount
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: keyArr
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:,:) :: atomIDs
    Real(kind=DoubleReal), Allocatable, Dimension(:,:,:) :: atomCoords
    Real(kind=DoubleReal), Dimension(1:27,1:3) :: ghostCellShift
  End Type oSubCell


  Type :: nlType
! 6.2MB per neighbour list
    Integer(kind=StandardInteger) :: length = 0
    Integer(kind=StandardInteger) :: arrayLength = 0
! Coord details - unit cell
    Integer(kind=StandardInteger) :: coordsLength
    Real(kind=DoubleReal) :: aLat
    Real(kind=DoubleReal) :: volume     ! Volume of supercell
    Real(kind=DoubleReal) :: volumeDistorted     ! Volume of supercell
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: unitCell
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: unitCellScaled   ! aLat * (UnitCell)
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: unitCellScaledInverse   ! (aLat * (UnitCell))^-1
    Integer(kind=StandardInteger) :: subCellWidth
! Unique atom IDs count & list
    Integer(kind=StandardInteger) :: atomID_Count
    Character(Len=16), Dimension(1:p_cuMax) :: atomIDs
! Coord details - coords, velocity etc
    Character(Len=16), Allocatable, Dimension(:) :: label                ! Unit labels
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: labelID
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: coords    ! Initial Coords (Unit coords [0,1]
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: coordsMD  ! MD coords for Update NL list - transformed coords
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: coordsFixed      ! Fix some atoms in MD
! Coord details - subcell
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: subCellKey
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: cell
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: ghostCell  !x,y,z,uKey
! Coord details - forces/velocity/acceleration
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: forces        ! Initial Coords (Unit coords [0,1]
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: velocity      ! Initial Coords (Unit coords [0,1]
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: velocityH     ! MD step velocity half time step
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: acceleration  ! Initial Coords (Unit coords [0,1]
! Initial force/velocity
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: forcesInitial
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: velocityInitial
! Coord details: particle details
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: charge   ! e
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: mass             ! amu
! Potential Keys
    Logical :: keysSet = .false.
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: pairKeyArray  ! Total potentials for each combination pair stored in "0"
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:,:) :: densityKeyArray  ! Total potentials for each combination pair stored in "0"
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:,:) :: embeddingKeyArray  ! Total potentials for each combination pair stored in "0"
! Atom details - EFS Calculation
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: electronDensity   ! Moved to nl object
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: atomEnergy        ! Moved to nl object
    Real(kind=DoubleReal) :: pairEnergy
    Real(kind=DoubleReal) :: embeddingEnergy
    Real(kind=DoubleReal) :: totalEnergy
! Neighbour List Details
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: nlOverKey   ! Used to isolate individual atoms and their local neighbour list
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: atomA_ID
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: atomB_ID
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: atomA_Type
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: atomB_Type
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: atomPairKey
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: inCell
! Position Details
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: rD
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: rD_original
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: vecAB
! Subcell
    Integer(kind=StandardInteger) :: scCount = 0
    Real(kind=DoubleReal) :: scALat
    Integer(kind=StandardInteger) :: maxAtomsPerSC = 0
! Misc
    Real(kind=DoubleReal) :: rMin
    Real(kind=DoubleReal) :: rMax
    Real(kind=DoubleReal) :: rVerlet
    Real(kind=DoubleReal) :: rVerletSq
    Real(kind=DoubleReal) :: totalRD
    Real(kind=DoubleReal) :: totalRDSq
! NL Update Parameters
    Real(kind=DoubleReal) :: aLat_Original    ! Lattice parameter when NL built
! MD vars
    Real(kind=DoubleReal) :: approxSeparation
    Real(kind=DoubleReal) :: approxTimeStep


! Remove these
! MD forces/coord - updated as they change
!    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: forcesMD         ! MD step forces
!    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: velocityMD       ! MD step velocity
!    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: accelerationMD   ! MD step acceleration
!    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: subCell
  End Type nlType


  Type :: oNLTempRD
! Temporary storage array for atom separations
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: rD
  End Type oNLTempRD




  Type :: nlType_Opt
! 6.2MB per neighbour list
    Integer(kind=StandardInteger) :: length = 0
    Integer(kind=StandardInteger) :: arrayLength = 0
! Coord details - unit cell
    Integer(kind=StandardInteger) :: coordsLength
    Real(kind=DoubleReal) :: aLat
    Real(kind=DoubleReal) :: volume     ! Volume of supercell
    Real(kind=DoubleReal) :: volumeDistorted     ! Volume of supercell
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: unitCell
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: unitCellScaled   ! aLat * (UnitCell)
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: unitCellScaledInverse   ! (aLat * (UnitCell))^-1
    Integer(kind=StandardInteger) :: subCellWidth
! Unique atom IDs count & list
    Integer(kind=StandardInteger) :: atomID_Count
    Character(Len=16), Dimension(1:p_cuMax) :: atomIDs
! Coord details - coords, velocity etc
    Character(Len=16), Allocatable, Dimension(:) :: label                ! Unit labels
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: labelID
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: coords    ! Initial Coords (Unit coords [0,1]
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: coordsMD  ! MD coords for Update NL list - transformed coords
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: coords_Moved    !
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: coordsFixed      ! Fix some atoms in MD
! Coord details - subcell
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: subCellKey
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: cell
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: ghostCell  !x,y,z,uKey
! Coord details - forces/velocity/acceleration
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: forces        ! Initial Coords (Unit coords [0,1]
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: velocity      ! Initial Coords (Unit coords [0,1]
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: velocityH     ! MD step velocity half time step
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: acceleration  ! Initial Coords (Unit coords [0,1]
! Initial force/velocity
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: forcesInitial
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: velocityInitial
! Coord details: particle details
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: charge   ! e
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: mass             ! amu
! Potential Keys
    Logical :: keysSet = .false.
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: pairKeyArray  ! Total potentials for each combination pair stored in "0"
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:,:) :: densityKeyArray  ! Total potentials for each combination pair stored in "0"
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:,:) :: embeddingKeyArray  ! Total potentials for each combination pair stored in "0"
! Atom details - EFS Calculation
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: electronDensity   ! Moved to nl object
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: atomEnergy        ! Moved to nl object
    Real(kind=DoubleReal) :: pairEnergy
    Real(kind=DoubleReal) :: embeddingEnergy
    Real(kind=DoubleReal) :: totalEnergy
! Neighbour List Details
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: nlOverKey   ! Used to isolate individual atoms and their local neighbour list
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: atomA_ID
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: atomB_ID
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: atomA_Type
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: atomB_Type
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: atomPairKey
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: inCell
! Position Details
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: rD
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: rD_original
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: vecAB
! Subcell
    Integer(kind=StandardInteger) :: scCount = 0
    Real(kind=DoubleReal) :: scALat
    Integer(kind=StandardInteger) :: maxAtomsPerSC = 0
! Misc
    Real(kind=DoubleReal) :: rMin
    Real(kind=DoubleReal) :: rMax
    Real(kind=DoubleReal) :: rVerlet
    Real(kind=DoubleReal) :: rVerletSq
    Real(kind=DoubleReal) :: totalRD
    Real(kind=DoubleReal) :: totalRDSq
! NL Update Parameters
    Real(kind=DoubleReal) :: aLat_Original    ! Lattice parameter when NL built
! MD vars
    Real(kind=DoubleReal) :: approxSeparation
    Real(kind=DoubleReal) :: approxTimeStep
! Opt variables
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: electronDensityFixed    ! Initial Coords (Unit coords [0,1]
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: rD_moved   ! used when moving atoms "nlMoveAtom_Opt"
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: atomEnergyFixed
    Real(kind=DoubleReal) :: totalEnergyFixed
    Logical :: fixedDensityCalculated = .false.
    Logical :: fixedEnergyCalculated = .false.
  End Type nlType_Opt


  Type :: o1D_Nodes
    Integer(kind=StandardInteger) :: count
    Real(kind=DoubleReal) :: aLat
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: nodeCoords
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: nodeMass
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: nodeFixed
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: force
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: velocityH
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: velocity
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: acceleration
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: optStep
    Real(kind=DoubleReal) :: potentialEnergy
  End Type o1D_Nodes

  Type :: o1D_NodesNL
    Integer(kind=StandardInteger) :: length
    Real(kind=DoubleReal) :: rCutoff = 6.5D0
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: nodeA
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: nodeB
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: vecAB
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: RD
  End Type o1D_NodesNL


End Module geomTypes




!------------------------------------------------------------------------------------------------------------------------------------------
