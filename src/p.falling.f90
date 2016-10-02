! --------------------------------------------------------------!
! Crystal module
! crystalTypes, crystal
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
!
!
! ----------------------------------------
! Updated: 26th June 2016
! ----------------------------------------


Module fallingModTypes
! Setup Modules
  Use kinds
  Use geomTypes
! Force declaration of all variables
  Implicit None
! Vars:  Module Parameters
!  Integer(kind=StandardInteger), Parameter :: p_ = 3
! Make private
  Private
! Public Variables and Parameters
!  Public :: p_
! Public derived types
  Public :: oFalling
!------------------------------------------------------------------------------
! Defined Types
!------------------------------------------------------------------------------
! oCoord (label, labelID, xyz, force, velocity)

  Type :: oFallingData
    Integer :: pointCount
    Real(kind=DoubleReal), Dimension(:,:), Allocatable :: dataPoints
  End Type oFallingData

  Type :: oMassFunction
    Integer :: pointCount
    Real(kind=DoubleReal), Dimension(:,:), Allocatable :: dataPoints
  End Type oMassFunction

  Type :: oDensityFunction
    Integer :: pointCount
    Real(kind=DoubleReal), Dimension(:,:), Allocatable :: dataPoints
  End Type oDensityFunction

  Type :: oFalling
    Character(len=128) :: outputDirectory
    Character(len=128) :: densityFile
    Integer(kind=StandardInteger) :: steps
    Type(oDensityFunction) :: densityFunction
    Type(oMassFunction) :: massFunction
    Type(oFallingData) :: fallingData
    Real(kind=DoubleReal) :: maxR
    Real(kind=DoubleReal) :: maxMass
    Real(kind=DoubleReal) :: totalTime
    Real(kind=DoubleReal) :: initialV
    Real(kind=DoubleReal) :: initialX
  End Type oFalling

End Module fallingModTypes


Module fallingMod
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use rng
  Use constants
  Use strings
  Use general
  Use sortMod
  Use units
  Use interpolation
  Use programMod
  Use fallingModTypes
! Force declaration of all variables
  Implicit None
! Make private
  Private
! Public
  Public :: runFalling
! Interfaces
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! ------------------------------------------------------------------------!
! Crystal Subroutines
! ------------------------------------------------------------------------!

  Subroutine runFalling()
! Run crystal subroutine
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
! Vars:  Private
    Type(oFalling) :: fallingObj
! Read input file
    Call readInput(fallingObj)
! Make mass table from density input
    Call makeMassTable(fallingObj)
! Run
    Call runFallingCalc(fallingObj)
! Output
    Call outputLog(fallingObj)
  End Subroutine runFalling


  Subroutine readInput(fallingObj)
! Read input file
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oFalling) :: fallingObj
! Vars:  Private
    Integer(kind=StandardInteger) :: i, n
    Character(len=128), Dimension(1:10) :: fieldArray
    Integer(kind=StandardInteger) :: fieldCount
    Integer(kind=StandardInteger) :: densityFileRows
    Character(len=128), Dimension(1:4096) :: densityFileData
!-----------------------------------
! Init
!-----------------------------------
    fallingObj%densityFile = BlankString(fallingObj%densityFile)
!-----------------------------------
! Loop through file rows
!-----------------------------------
    Do i=1,programObj%programFiles(1)%rowsP
! Load keywords
      If(programObj%programFiles(1)%dataP_UC(i)(1:10).eq."#OUTPUTDIR")Then
        Call readFieldsCharacter(programObj%programFiles(1)%dataP(i),fieldArray,fieldCount)
        fallingObj%outputDirectory = BlankString(fallingObj%outputDirectory)
        fallingObj%outputDirectory = trim(fieldArray(2))
        fallingObj%outputDirectory = fillPath(fallingObj%outputDirectory)
      End If
      If(programObj%programFiles(1)%dataP_UC(i)(1:8).eq."#DENSITY")Then
        Call readFieldsCharacter(programObj%programFiles(1)%dataP(i),fieldArray,fieldCount)
        fallingObj%densityFile = trim(fieldArray(2))
        fallingObj%densityFile = fillPath(fallingObj%densityFile)
      End If
      If(programObj%programFiles(1)%dataP_UC(i)(1:6).eq."#STEPS")Then
        Call readFieldsCharacter(programObj%programFiles(1)%dataP_UC(i),fieldArray,fieldCount)
        fallingObj%steps = StrToInt(fieldArray(2))
      End If
      If(programObj%programFiles(1)%dataP_UC(i)(1:7).eq."#STARTV")Then
        Call readFieldsCharacter(programObj%programFiles(1)%dataP_UC(i),fieldArray,fieldCount)
        fallingObj%initialV = StrToDp(fieldArray(2))
        fallingObj%initialV = UnitConvert(fallingObj%initialV,fieldArray(3),"MS-1")
      End If
      If(programObj%programFiles(1)%dataP_UC(i)(1:7).eq."#STARTX")Then
        Call readFieldsCharacter(programObj%programFiles(1)%dataP_UC(i),fieldArray,fieldCount)
        fallingObj%initialX = StrToDp(fieldArray(2))
        fallingObj%initialX = UnitConvert(fallingObj%initialX,fieldArray(3),"M")
      End If
    End Do
    print *,"DENSITY FILE:         ",trim(fallingObj%densityFile)
    print *,"STEPS:                ",fallingObj%steps
    print *,"Initial X (m):        ",fallingObj%initialX
    print *,"Initial V (ms-1):     ",fallingObj%initialV
!-----------------------------------
! Read Density File
!-----------------------------------
    Call readFile(fallingObj%densityFile, densityFileData, densityFileRows)
! Allocate rows
    Allocate(fallingObj%densityFunction%dataPoints(1:densityFileRows,1:2))
! Loop through rows
    fallingObj%densityFunction%pointCount = densityFileRows
    Do n=1,densityFileRows
      Call readFieldsCharacter(densityFileData(n),fieldArray,fieldCount,",")
      fallingObj%densityFunction%dataPoints(n,1) = StrToDp(fieldArray(1)(1:32))
      fallingObj%densityFunction%dataPoints(n,2) = StrToDp(fieldArray(2)(1:32))
    End Do
! Sort array
    Call sortArray(fallingObj%densityFunction%dataPoints)
! Set maximum R
    fallingObj%maxR = fallingObj%densityFunction%dataPoints(fallingObj%densityFunction%pointCount,1)
  End Subroutine readInput

  Subroutine makeMassTable(fallingObj)
! Read input file
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oFalling) :: fallingObj
! Vars:  Private
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal) :: radius, radius3, radius3Last, mass, totalMass
! Allocate array
    Allocate(fallingObj%massFunction%dataPoints(1:fallingObj%densityFunction%pointCount,1:2))
! Calculate total mass below each "step"
    radius3Last = 0.0D0
    totalMass = 0.0D0
    Do i=1,fallingObj%densityFunction%pointCount
      radius = fallingObj%densityFunction%dataPoints(i,1)
      radius3 = radius**3
      If(radius.le.0.0D0)Then
        mass = 0.0D0
      Else
        mass = spherepi*(radius3-radius3Last)*fallingObj%densityFunction%dataPoints(i,2)
      End If
      totalMass = totalMass + mass
! Store
      fallingObj%massFunction%dataPoints(i,1) = radius
      fallingObj%massFunction%dataPoints(i,2) = totalMass
      radius3Last = radius3
    End Do
    fallingObj%maxMass = totalMass
  End Subroutine makeMassTable

  Subroutine runFallingCalc(fallingObj)
! Read input file
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oFalling) :: fallingObj
! Vars:  Private
    Integer(kind=StandardInteger) :: i, n, k
    Integer(kind=StandardInteger) :: mins, intPoints
    Real(kind=DoubleReal) :: x, s, sStart, sEnd, sInc, sInc_K
    Real(kind=DoubleReal) :: velocity, mass, accel, u, v, t, time, seconds, accelAvg
    Real(kind=DoubleReal) :: velocityLast, timeLast
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Allocate falling data array
    Allocate(fallingObj%fallingData%dataPoints(1:fallingObj%steps,1:7))
! start position - height above ground/m
    x = fallingObj%initialX
    velocity = fallingObj%initialV
    time = 0.0D0
    intPoints = 50
! distance from center
    s = fallingObj%maxR+x
    sStart = s
! increment
    sInc = (s/fallingObj%steps)
    Do n=1,fallingObj%steps
      velocityLast = velocity
      timeLast = time
      sEnd = sStart-sInc
      sInc_K = sInc/(1.0D0*intPoints)
      s = sStart+0.5D0*sInc_K
      accelAvg = 0.0D0
      Do k=1,intPoints
        s = s-sInc_K
        If(s.ge.fallingObj%maxR)Then
          mass = fallingObj%maxMass
        Else
          yArray = PointInterp(fallingObj%massFunction%dataPoints,s,4)
          mass = yArray(1)
        End If
        If(mass.lt.0.0D0)Then
          mass = 0.0D0
        End If
! Calculate acceleration
        accel = (gravitationalConstant * mass)/(s*s)
        accelAvg = accelAvg + accel
        u = velocity
        v = sqrt(u*u+2.0D0*accel*sInc_K)
        If(accel.gt.0)Then
          t = (v-u)/accel
        Else
          t = 0.0D0
        End If
        velocity = v        ! Update velocity
        time = time + t     !
      End Do
      fallingObj%fallingData%dataPoints(n,1) = sStart
      fallingObj%fallingData%dataPoints(n,2) = sEnd
      fallingObj%fallingData%dataPoints(n,3) = velocityLast
      fallingObj%fallingData%dataPoints(n,4) = velocity
      fallingObj%fallingData%dataPoints(n,5) = accelAvg/(1.0D0*intPoints)
      fallingObj%fallingData%dataPoints(n,6) = timeLast
      fallingObj%fallingData%dataPoints(n,7) = time
      sStart = sEnd
    End Do
! Store total time through
    fallingObj%totalTime = 2.0D0 * time
! Output
    mins = Floor(fallingObj%totalTime/60.0D0)
    seconds = fallingObj%totalTime-60.0D0*mins
    print "(A12,I3,A6,F10.6,A8)","Total time: ",mins," mins ",seconds," seconds"
  End Subroutine runFallingCalc

  Subroutine outputLog(fallingObj)
! Read input file
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oFalling) :: fallingObj
! Vars:  Private
    Integer(kind=StandardInteger) :: n, m
    Character(len=128) :: filePath
    Real(kind=DoubleReal) :: time, dt
! Make dir
    Call makeDir(fallingObj%outputDirectory)
! Save to file
    filePath = JoinStr(fallingObj%outputDirectory,"/falling.log")
    open(unit=999,file=trim(filePath))
    write(999,"(A16,A16,A16,A16,A16,A16,A16)") &
    "r start (m)","r end (m)","v start (ms-1)","v end (ms-1)",&
    "accel (ms-2)","t start (s)","t end (s)"
    Do n = 1,fallingObj%steps
      write(999,"(F16.5,F16.5,F16.5,F16.5,F16.5,F16.5,F16.5)") &
      fallingObj%fallingData%dataPoints(n,1), &
      fallingObj%fallingData%dataPoints(n,2), &
      fallingObj%fallingData%dataPoints(n,3), &
      fallingObj%fallingData%dataPoints(n,4), &
      fallingObj%fallingData%dataPoints(n,5), &
      fallingObj%fallingData%dataPoints(n,6), &
      fallingObj%fallingData%dataPoints(n,7)
    End Do
    time = fallingObj%fallingData%dataPoints(fallingObj%steps,7)
    Do m = 1,fallingObj%steps
      n = fallingObj%steps+1-m
      dt = fallingObj%fallingData%dataPoints(n,7)-fallingObj%fallingData%dataPoints(n,6)
      write(999,"(F16.5,F16.5,F16.5,F16.5,F16.5,F16.5,F16.5)") &
      (-1.0D0*fallingObj%fallingData%dataPoints(n,2)), &
      (-1.0D0*fallingObj%fallingData%dataPoints(n,1)), &
      fallingObj%fallingData%dataPoints(n,4), &
      fallingObj%fallingData%dataPoints(n,3), &
      fallingObj%fallingData%dataPoints(n,5), &
      time, &
      (time+dt)
      time = time + dt
    End Do
    close(999)
  End Subroutine outputLog

End Module fallingMod

!-----------------------------------------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------------------------------------
