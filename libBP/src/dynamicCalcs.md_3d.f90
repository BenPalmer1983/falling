! ------------------------------------------------------------
!           DYNAMICSCALCS: 3D Molecular Dynamics
! ------------------------------------------------------------


  Subroutine runDynamics(dynamicsSettings, coords, potential, cKey)
! Run Molecular Dynamics Subroutine
! Makes neighbour list
! Sorts out the potential keys
! Runs Verlet scheme
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oDynamics) :: dynamicsSettings
    Type(coordsType), Dimension(:) :: coords
    Type(potentialType) :: potential
    Integer(kind=StandardInteger) :: cKey
! Vars:  Private
    Integer(kind=StandardInteger) :: i, configCount
    Type(nlType), Allocatable, Dimension(:) :: nl
    Type(coordsType), Dimension(1:size(coords,1)) :: coordsMD

! Init vars
    configCount = size(coords,1)
    coordsMD = coords
!----------------------------------
! Allocate arrays on HEAP
!----------------------------------
    Allocate(nl(1:configCount))

    print *,"Run Dynamics"
! Make initial neighbour list
    Call makeNL(nl, coordsMD, 6.5D0)       ! geom.f90
! Build potential keys arrays
    Call nlPotentialKeys(nl, potential)
! First EFS calculation
    Call calcEFS(nl, potential, cKey)
! Calc starting acceleration
    Call calcAcceleration(nl, cKey)

    Call estimateTimeStep(nl, cKey)
    dynamicsSettings%mdTime = 0.0D0
    dynamicsSettings%mdStep = 0
    print *,nl(cKey)%approxTimeStep
!
    !print *,"0"
    !print *,"C",nl(ckey)%coordsMD(1,1),nl(ckey)%coordsMD(1,2),nl(ckey)%coordsMD(1,3)
    !print *,"V",nl(ckey)%velocity(1,1),nl(ckey)%velocity(1,2),nl(ckey)%velocity(1,3)
    !print *,"H",nl(ckey)%velocityH(1,1),nl(ckey)%velocityH(1,2),nl(ckey)%velocityH(1,3)
    !print *,"F",nl(ckey)%forces(1,1),nl(ckey)%forces(1,2),nl(ckey)%forces(1,3)
    !print *,"A",nl(ckey)%acceleration(1,1),nl(ckey)%acceleration(1,2),nl(ckey)%acceleration(1,3)
    !print *,""
    print *,"0",nl(ckey)%coordsMD(1,1),nl(ckey)%acceleration(1,1)

    Do i=1,dynamicsSettings%timeSteps
      Call verletStep(dynamicsSettings, nl, potential, cKey)
      !print *,i
      !print *,"C",nl(ckey)%coordsMD(1,1),nl(ckey)%coordsMD(1,2),nl(ckey)%coordsMD(1,3)
      !print *,"V",nl(ckey)%velocity(1,1),nl(ckey)%velocity(1,2),nl(ckey)%velocity(1,3)
      !print *,"H",nl(ckey)%velocityH(1,1),nl(ckey)%velocityH(1,2),nl(ckey)%velocityH(1,3)
      !print *,"F",nl(ckey)%forces(1,1),nl(ckey)%forces(1,2),nl(ckey)%forces(1,3)
      !print *,"A",nl(ckey)%acceleration(1,1),nl(ckey)%acceleration(1,2),nl(ckey)%acceleration(1,3)
      !print *,""
      print *,i,nl(ckey)%coordsMD(1,1),nl(ckey)%acceleration(1,1)
    End Do
!----------------------------------
! Deallocate arrays on HEAP
!----------------------------------
    Deallocate(nl)
  End Subroutine runDynamics




  Subroutine estimateTimeStep(nl, cKey)
! Estimate time step
! Uses initial velocity and first calculated force
! !! Only working for initial guess !!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType), Dimension(:) :: nl
    Integer(kind=StandardInteger) :: cKey
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j
    Real(kind=DoubleReal) :: maxForce, testForce, minMass, maxAcceleration, maxVelocity
! Calc approximate atom separation
    nl(cKey)%approxSeparation = (nl(cKey)%volume/(1.0D0*nl(cKey)%coordsLength))**(1/3.0D0)
! Get max force
    maxForce = 0.0D0
    maxVelocity = 0.0D0
    minMass = nl(cKey)%mass(1)
    Do i=1,nl(cKey)%coordsLength
      Do j=1,3
        If(nl(cKey)%forces(i,j).gt.maxForce)Then
          maxForce = nl(cKey)%forces(i,j)
        Else
          testForce = -1.0D0 * nl(cKey)%forces(i,j)
          If(testForce.gt.maxForce)Then
            maxForce = testForce
          End If
        End If
        If(nl(cKey)%velocity(i,j).gt.maxVelocity)Then
          maxVelocity = nl(cKey)%velocity(i,j)
        End If
      End Do
      If(nl(cKey)%mass(i).lt.minMass)Then
        minMass = nl(cKey)%mass(i)
      End If
    End Do
    maxAcceleration = F_accelerationFM_MD(maxForce, minMass)
! Estimate for a 1/10 estimated atom separation distance step
! (A/2)t^2 + v0 t - x = 0
    nl(cKey)%approxTimeStep = Quadratic(maxAcceleration/2.0D0,maxVelocity,(-1.0D0*nl(cKey)%approxSeparation),1)
  End Subroutine estimateTimeStep


  Subroutine verletStep(dynamicsSettings, nl, potential, cKey)
! Calculate new position and velocity of each atom
! using Verlet step scheme
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oDynamics) :: dynamicsSettings
    Type(nlType), Dimension(:) :: nl
    Type(potentialType) :: potential
    Integer(kind=StandardInteger) :: cKey
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j
    Real(kind=DoubleReal) :: dt, dtHalf, dtsq
! Init vars
    If(dynamicsSettings%adaptiveTimestep)Then
      Call estimateTimeStep(nl, cKey)
      dt = nl(cKey)%approxTimeStep
      If(dt.lt.dynamicsSettings%stepSizeMin)Then
        dt = dynamicsSettings%stepSizeMin
      End If
      If(dt.gt.dynamicsSettings%stepSizeMax)Then
        dt = dynamicsSettings%stepSizeMax
      End If
    Else
      dt = dynamicsSettings%stepSize
    End If
    dtHalf = 0.5D0*dt
    dtsq = dt**2
! Increase time step
    dynamicsSettings%mdTime = dynamicsSettings%mdTime + dt
    dynamicsSettings%mdStep = dynamicsSettings%mdStep + 1
!------------------------------------------------------------------------------------------
! Step 1 + 2
    Do i=1,nl(cKey)%coordsLength
      If(nl(cKey)%coordsFixed(i).eq.0)Then  ! Only move/calc velocity for unfixed points
! Step 1 - Calculate Half dt Velocity vec(v)(t+0.5dt)
        Do j=1,3
          nl(cKey)%velocityH(i,j) = &
            nl(cKey)%velocity(i,j) + dtHalf * nl(cKey)%acceleration(i,j)
        End Do
! Step 2 - Calculate vec(r)(t+dt)
        Call moveAtom_V(nl, cKey, i, dt)
      End If
    End Do
! Step 3 - Refresh or Rebuild Neighbour List
    dynamicsSettings%rebuildCount = dynamicsSettings%rebuildCount + 1
    If(dynamicsSettings%rebuildCount.eq.dynamicsSettings%rebuildSteps)Then
      !Call rebuildNL(nl, 6.5D0, cKey)   ! geom.nl.f90
      Call refreshNL(nl, cKey)         ! geom.nl.f90
      dynamicsSettings%rebuildCount = 0
    Else
      Call refreshNL(nl, cKey)         ! geom.nl.f90
    End If
! Step 4a - Recalculate forces on atoms in new position
    Call calcEFS(nl, potential, cKey, .false.)   ! staticCalcs.calcEFS.f90  ! Fills in forcesMD array
! Step 4b - Calculate acceleration
    Call calcAcceleration(nl, 1)                                   ! Fills in accelerationMD array
! Step 5 - Calculate final velocity
    Do i=1,nl(cKey)%coordsLength
      Do j=1,3
        nl(cKey)%velocity(i,j) = &
          nl(cKey)%velocityH(i,j) + dtHalf * nl(cKey)%acceleration(i,j)
      End Do
    End Do
  End Subroutine verletStep

  Subroutine calcAcceleration(nl, cKey)
! Calculate acceleration ang/ps^2
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType), Dimension(:) :: nl
    Integer(kind=StandardInteger) :: cKey
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j
! Loop through coords
    Do i=1,nl(cKey)%coordsLength
      Do j=1,3
        nl(cKey)%acceleration(i,j) = &
          F_accelerationFM_MD(nl(cKey)%forces(i,j), nl(cKey)%mass(i))  ! scienceFunctions.classicalMechanics.f90 in MD units
      End Do
    End Do
  End Subroutine calcAcceleration

  Subroutine moveAtom_V(nl, cKey, atomID, dt)
! Move an atom due to it's dt/2 velocity and a period of time dt
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType), Dimension(:) :: nl
    Integer(kind=StandardInteger) :: cKey, atomID
    Real(kind=DoubleReal) :: dt
! Vars:  Private
    Integer(kind=StandardInteger) :: j
! Move atom
    Do j=1,3
      nl(cKey)%coordsMD(atomID,j) = &
        nl(cKey)%coordsMD(atomID,j) + nl(cKey)%velocityH(atomID,j) * dt
    End Do
  End Subroutine moveAtom_V
