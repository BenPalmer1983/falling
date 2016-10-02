! ------------------------------------------------------------
!           DYNAMICSCALCS: Types
! ------------------------------------------------------------


Module dynamicCalcsTypes
! Setup Modules
  Use kinds
! Force declaration of all variables
  Implicit None
! Vars:  Module Parameters
!  Integer(kind=StandardInteger), Parameter :: p_
! Make private
  Private
! Public Variables and Parameters
!  Public :: p_
! Public derived data types
  Public :: oDynamics

  Type :: oDynamics
    Integer(kind=StandardInteger) :: timeSteps=250
    Real(kind=DoubleReal) :: stepSize = 0.005D0 ! picoseconds
    Logical :: adaptiveTimestep = .false.
    Real(kind=DoubleReal) :: stepSizeMin = 0.001D0
    Real(kind=DoubleReal) :: stepSizeMax = 0.010D0
    Logical :: damp = .false.
    Integer(kind=StandardInteger) :: rebuildSteps=10
    Integer(kind=StandardInteger) :: rebuildCount=0
    Integer(kind=StandardInteger) :: mdStep = 0
    Real(kind=DoubleReal) :: mdTime = 0.0D0
  End Type oDynamics

End Module dynamicCalcsTypes














!------------------------------------------------------------------------------------------------------------------------------------------
