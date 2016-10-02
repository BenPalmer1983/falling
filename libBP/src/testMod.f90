Module testMod
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
!
! Test functions for Library
! --------------------------------------------------------------!
  Use kinds
  Use time
  Use logicalMod
  Use strings
  Use general
  Use envTypes
  Use env
  Use printModTypes
  Use printMod
  Use activityFunctionsTypes
  Use activityFunctions
  Use basicMaths
  Use MatrixCalculus
  Use keysMod
  Use newtonMethods
  Use lmaM
  Use conjugateGradients
  Use coordFunctions
  Use laplaceTransforms
  Use interpolation
  Use regression
  Use splinesFitting
  Use geomTypes
  Use geom
  Use rng
  Use rngDistTypes
  Use rngDist
  Use potentialsTypes
  Use potentials
  Use staticCalcsTypes
  Use staticCalcs
  Use dynamicCalcsTypes
  Use dynamicCalcs
  Use isotopesTypes
  Use isotopes
  Use potentialFitting
! Force declaration of all variables
  Implicit None
! Public variables
! Make private
  Private
! Public
! --Subroutines--!
  Public :: testActivity
  Public :: testCalcIsotopeAmount
  Public :: testGaverStehfest
  Public :: testDecayChain
  Public :: testNeighbourList
  Public :: testNewton
  Public :: testGeomOptimisation
  Public :: testStaticCalc
  Public :: testNodesOpt
  Public :: testPotentialFitting
  Public :: testCombinations
  Public :: testLinearRegression
  Public :: testMorseFit
  Public :: testRandomDistribution

!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! ---------------------------------------------------------
! MODULE SUBROUTINES
! ---------------------------------------------------------


! ------------------------------------------------------------------------------------------
!  Test subroutines for activityFunctions.f90
! ------------------------------------------------------------------------------------------

  Subroutine testActivity()
! Tests activity function
    Implicit None ! Force declaration of all variables
! Vars: In/Out
!   None
! Vars: Private
    Print *,"Testing activity function"
! Polonium-218
    Call ActivityCompareGS(0.002236D0, 2000, 100000.0D0, 1.0D6, 1000.0D0)
  End Subroutine testActivity


  Subroutine testCalcIsotopeAmount()
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: n
    Real(kind=DoubleReal) :: w, t
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: decayDataArray
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: isotopeChange

    t = 60.0D0
    w = 0.0D0
    Allocate(decayDataArray(1:5,1:6))
    decayDataArray(1,1) = 1.0D0
    decayDataArray(1,2) = 1.0D6
    decayDataArray(1,3) = 30.0D0
    decayDataArray(1,4) = 1.0D0
    decayDataArray(1,5) = 1.0D0
    decayDataArray(1,6) = 1.0D0
    decayDataArray(2,1) = 2.0D0
    decayDataArray(2,2) = 0.0D0
    decayDataArray(2,3) = 122.0D0
    decayDataArray(2,4) = 0.9998D0
    decayDataArray(2,5) = 1.0D0
    decayDataArray(2,6) = 1.0D0
    decayDataArray(3,1) = 2.0D0
    decayDataArray(3,2) = 0.0D0
    decayDataArray(3,3) = 1196.0D0
    decayDataArray(3,4) = 1.0D0
    decayDataArray(3,5) = 1.0D0
    decayDataArray(3,6) = 1.0D0
    decayDataArray(4,1) = 2.0D0
    decayDataArray(4,2) = 0.0D0
    decayDataArray(4,3) = 1648.0D0
    decayDataArray(4,4) = 0.9998D0
    decayDataArray(4,5) = 1.0D0
    decayDataArray(4,6) = 1.0D0
    decayDataArray(5,1) = 2.0D0
    decayDataArray(5,2) = 0.0D0
    decayDataArray(5,3) = 7006.0D0
    decayDataArray(5,4) = 1.0D0
    decayDataArray(5,5) = 1.0D0
    decayDataArray(5,6) = 1.0D0
    Allocate(isotopeChange(1:5,1:12))
    print *,""
    isotopeChange = CalcIsotopeAmount(w,decayDataArray,t,1)
    Do n = 1,5
      Print *,n,isotopeChange(n,3),isotopeChange(n,4)
    End Do
    print *,""
    isotopeChange = CalcIsotopeAmount(w,decayDataArray,t,2)
    Do n = 1,5
      Print *,n,isotopeChange(n,3),isotopeChange(n,4)
    End Do
  End Subroutine testCalcIsotopeAmount


  Subroutine testGaverStehfest()
! Tests activity function
    Implicit None ! Force declaration of all variables
! Vars: In/Out
!   None
! Vars: Private
    Real(kind=DoubleReal) :: t, ft, exact
    Real(kind=DoubleReal), Dimension(1:10) :: p
    Print *,"Testing G-S Inversion"
! Init
    t = 0.024D0
    p(1) = 0.731D0
    ft = GaverStehfest(ltExp, t, p, 7)
    exact = exp(p(1)*t)
    print *,"m=7 (double) ",t,exact,ft
    ft = GaverStehfest(ltExp, t, p, 8)
    exact = exp(p(1)*t)
    print *,"m=8 (double) ",t,exact,ft
  End Subroutine testGaverStehfest

  Subroutine testDecayChain()
! Tests activity function
    Implicit None ! Force declaration of all variables
! Vars: In/Out
!   None
! Vars: Private
    Type(decayChainObj) :: decayChain
    Type(activityTimeObj), Dimension(1:100) :: activityTime
    Real(kind=DoubleReal) :: endTime, zeroProductionTime
!-------------------------------------------
    Print *,"Testing decay chain"
! Polonium-218 decay chain
    !decayChain%productionRate = 0.0D0
    decayChain%time = 60.0D0
! Parent
    decayChain%label(1) = "Po-218"
    decayChain%productionRate(1) = 5.0D4
    decayChain%branchFactor(1) = 1.0D0  ! Not used
    decayChain%halfLife(1) = 185.88D0
    decayChain%amountStart(1) = 1.0D6
! Daughter 1
    decayChain%label(2) = "Pb-214"
    decayChain%productionRate(2) = 0.0D0
    decayChain%branchFactor(2) = 0.99981D0   ! bf from 1 to 2
    decayChain%halfLife(2) = 1608.0D0
    decayChain%amountStart(2) = 1.0D5
! Daughter 2
    decayChain%label(3) = "Bi-214"
    decayChain%productionRate(3) = 1.0D4
    decayChain%branchFactor(3) = 1.0D0  ! bf from 2 to 3
    decayChain%halfLife(3) = 1194.0D0
    decayChain%amountStart(3) = 500.0D0
! Daughter 3
    decayChain%label(4) = "Po-214"
    decayChain%productionRate(4) = 0.0D0
    decayChain%branchFactor(4) = 0.99979D0  ! bf from 3 to 4
    decayChain%halfLife(4) = 0.0001637D0
    decayChain%amountStart(4) = 0.0D0
! Daughter 4
    decayChain%label(5) = "Pb-210"
    decayChain%productionRate(5) = 0.0D0
    decayChain%branchFactor(5) = 1.0D0  ! bf from 3 to 4
    decayChain%halfLife(5) = 6.9930D8
    decayChain%amountStart(5) = 0.0D0
! Daughter 5
    decayChain%label(6) = "Bi-210"
    decayChain%productionRate(6) = 0.0D0
    decayChain%branchFactor(6) = 1.0D0  ! bf from 4 to 5
    decayChain%halfLife(6) = 1.600665D-6
    decayChain%amountStart(6) = 0.0D0
! Daughter 6
    decayChain%label(7) = "Tl-206"
    decayChain%productionRate(7) = 0.0D0
    decayChain%branchFactor(7) = 1.0D0  ! bf from 5 to 6
    decayChain%halfLife(7) = 2.5212D2
    decayChain%amountStart(7) = 0.0D0
! Daughter 7
    decayChain%label(8) = "Pb-206"
    decayChain%productionRate(8) = 0.0D0
    decayChain%branchFactor(8) = 1.0D0  ! bf from 6 to 7
    decayChain%halfLife(8) = -1.0D0
    decayChain%amountStart(8) = 0.0D0
! Calculate
    endTime = 10000.0D0
    zeroProductionTime = 1000.0D0
!   print *,endTime,zeroProductionTime
    Call CalcActivities(decayChain,activityTime,endTime,zeroProductionTime)
    Call CalcActivitiesPrint(decayChain,activityTime,.false.)

  End Subroutine testDecayChain



  Subroutine testNeighbourList()
! Uses inverse laplace transform to calculate isotope amounts at time t (after time = 0)
! t time in seconds after t=0
! w production rate of parent isotope
! isotope chain data
    Implicit None ! Force declaration of all variables
! Vars Private
    Type(coordsUnitType), Dimension(1:4) :: coordsUnit
    Type(coordsType), Dimension(1:4) :: coords
    Type(nlType), Dimension(1:4) :: nl
    Integer(kind=StandardInteger) :: i, cKey
! Start test
    print *,"Neighbour List Testing"
    print *,"==========================================================="
    Call initUnitCoords(coordsUnit)
    !Call initCoords(coords)
    coordsUnit%aLat = 4.04D0
    coordsUnit%xCopy = 4
    coordsUnit%yCopy = 4
    coordsUnit%zCopy = 4
    !Call standardCoords("FCC", coordsUnit)
    !Call expandUnitCoords(coordsUnit, coords)
    Do cKey=1,size(coords%length)
      Do i=1,coords(cKey)%length
        print *,coords(cKey)%label(i),&
        coords(cKey)%coords(i,1),coords(cKey)%coords(i,2),coords(cKey)%coords(i,3)
      End Do
    End Do
    print *,coords%length
    Call makeNL(nl, coords, 6.5D0)
  End Subroutine testNeighbourList






  Subroutine testNewton()
! Test static calculation interatomic potentials
    Implicit None ! Force declaration of all variables
! Vars Private
    Real(kind=DoubleReal), Dimension(1:2) :: parametersIn, parametersOut
    Real(kind=DoubleReal), Dimension(1:3) :: pQuadIn, pQuadOut
    Real(kind=DoubleReal), Dimension(1:4,1:2) :: points
    Real(kind=DoubleReal), Dimension(1:4,1:4) :: Ap
    Real(kind=DoubleReal), Dimension(1:4) :: bp, xp, x0p
    Real(kind=DoubleReal), Dimension(1:2) :: xVec, gradVec
    Real(kind=DoubleReal), Dimension(1:3) :: pVec
    Real(kind=DoubleReal), Dimension(1:2,1:2) :: hessian


    parametersIn(1) = -1.0D0
    parametersIn(2) = 1.0D0

    parametersOut = NewtonSolve(NewtonSolve_Sample, parametersIn)
    print *,"x",parametersOut(1)
    print *,"y",parametersOut(2)



    points(1,1) = 0.0D0
    points(2,1) = 0.5D0
    points(3,1) = 1.0D0
    points(4,1) = 1.5D0

    points(1,2) = 2.70D0
    points(2,2) = 1.36D0
    points(3,2) = 1.22D0
    points(4,2) = 2.28D0

    pQuadIn(1) = 1.0D0
    pQuadIn(2) = 1.0D0
    pQuadIn(3) = 1.0D0
    pQuadOut = NewtonGauss(NewtonGauss_Quadratic, pQuadIn, points)
    print *,"NG ",pQuadOut(1),pQuadOut(2),pQuadOut(3)


    pQuadOut = LMA(points, LMA_Quadratic, pQuadIn)
    print *,"LMA ",pQuadOut(1),pQuadOut(2),pQuadOut(3)


    Ap(1,1) = 2.4D0
    Ap(1,2) = 5.2D0
    Ap(1,3) = 13.0D0
    Ap(1,4) = 1.2D0

    Ap(2,1) = 5.2D0
    Ap(2,2) = 1.35D0
    Ap(2,3) = 4.25D0
    Ap(2,4) = 1.1D0

    Ap(3,1) = 13.0D0
    Ap(3,2) = 4.25D0
    Ap(3,3) = 6.2D0
    Ap(3,4) = 7.6D0

    Ap(4,1) = 1.2D0
    Ap(4,2) = 1.1D0
    Ap(4,3) = 7.6D0
    Ap(4,4) = 0.15D0


    bp(1) = 2.4D0
    bp(2) = 9.2D0
    bp(3) = 1.4D0
    bp(4) = 3.7D0

    x0p(1) = 1
    x0p(2) = 1
    x0p(3) = 1
    x0p(4) = 1

    xp = cgLinear(Ap, bp, x0p)
    print *,xp(1),xp(2),xp(3),xp(4)


    xVec(1) = 0.2D0
    xVec(2) = 0.7D0
    pVec(1) = 0.3D0
    pVec(2) = 2.1D0
    pVec(3) = 1.3D0
    gradVec = matrixGradient(MC_SampleA,xVec,pVec)
    print *,"grad",gradVec(1),gradVec(2)

    xVec(1) = 0.2D0
    xVec(2) = 0.7D0
    gradVec = matrixGradient(MC_SampleB,xVec,pVec)
    print *,"grad",gradVec(1),gradVec(2)


    xVec(1) = 0.2D0
    xVec(2) = 0.7D0
    hessian = matrixHessian(MC_SampleB,xVec,pVec)
    print *,"H ",hessian(1,1),hessian(1,2)
    print *,"H ",hessian(2,1),hessian(2,2)
    print *,""

    Call matrixGH(MC_SampleB,xVec,pVec,gradVec,hessian)
    print *,"grad",gradVec(1),gradVec(2)
    print *,"H ",hessian(1,1),hessian(1,2)
    print *,"H ",hessian(2,1),hessian(2,2)
    xVec(1) = 0
    xVec(2) = 0
    xVec = NewtonOptimise(NewtonOpt_SampleA, xVec)
    print *,xVec(1),xVec(2)
    xVec(1) = 0
    xVec(2) = 0
    xVec = NewtonOptimise(NewtonOpt_SampleB, xVec)
    print *,xVec(1),xVec(2)
  End Subroutine testNewton







  Subroutine testGeomOptimisation()
! Test static calculation interatomic potentials
    Implicit None ! Force declaration of all variables
! Vars Private
    Type(coordsUnitType), Allocatable, Dimension(:) :: coordsUnit
    Type(coordsType), Allocatable, Dimension(:) :: coords
    Type(nlType), Allocatable, Dimension(:) :: nl
    Type(potentialType) :: potential
    Character(Len=2), Dimension(1:4) :: atomLabels
    Integer(kind=StandardInteger) :: cKey
! Start time
    Call setStartTime()
!----------------------------------
! Allocate arrays on HEAP
!----------------------------------
    Allocate(coordsUnit(1:1))
    Allocate(coords(1:1))
    Allocate(nl(1:1))
! Set up coords
    Call initUnitCoords(coordsUnit)     ! geom.f90
    Call initCoords(coords)             ! geom.f90
    cKey = 1
    coordsUnit%aLat = 4.04D0
    coordsUnit%xCopy = 4
    coordsUnit%yCopy = 4
    coordsUnit%zCopy = 4
    atomLabels = "AL"
    Call standardCoords("FCC", coordsUnit,cKey,atomLabels)  ! geom.f90
    Call expandUnitCoords(coordsUnit, coords) ! geom.f90
! Init + load potential
    Call initPotential(potential)              ! potentials.f90
    Call loadPotential("test.pot", potential)  ! potentials.f90
! Assign atom IDs
    Call addLinePage("Assign IDs")
    Call atomLabelIDs(potential, coords)
! Print potential summary
    Call printPotentialSummary(potential)      ! potentials.f90
! Print
    Call printAtomLabelIDs(coords)
!
! Starting/opt config
    Call makeNL(nl, coords, 6.5D0, cKey)       ! geom.f90
    Call calcEFS(nl, potential, cKey)
    print *,"Start calc A ",nl(cKey)%totalEnergy,(nl(cKey)%totalEnergy/nl(cKey)%coordsLength)
    Call calcE(nl, potential, cKey)
    print *,"Start calc B ",nl(cKey)%totalEnergy,(nl(cKey)%totalEnergy/nl(cKey)%coordsLength)

! Heat up coords
    Call HeatCoords(coords, 0.01D0, cKey)
! Make neighbour list
    Call makeNL(nl, coords, 6.5D0, cKey)       ! geom.f90
    Call calcEFS(nl, potential, cKey)
    print *,"Heated calc A ",nl(cKey)%totalEnergy,(nl(cKey)%totalEnergy/nl(cKey)%coordsLength)
    Call calcE(nl, potential, cKey)
    print *,"Heated calc B ",nl(cKey)%totalEnergy,(nl(cKey)%totalEnergy/nl(cKey)%coordsLength)

! Optimise
    Call optGeom(coords, potential, .true., .false.)

    Call printCoords(coords, cKey)


    Call expandUnitCoords(coordsUnit, coords) ! geom.f90
    coords(cKey)%coords(1,1) = 0.05
    coords(cKey)%coords(1,2) = 0.05
    coords(cKey)%coords(1,3) = 0.05


End Subroutine testGeomOptimisation

  Subroutine testNodesOpt()
! Test static calculation interatomic potentials
    Implicit None ! Force declaration of all variables
! Vars Private
    Type(o1D_Nodes) :: nodes
    Type(o1D_NodesNL) :: nl
    Type(o1D_Nodes) :: nodesTemp
    Type(o1D_NodesNL) :: nlTemp
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal) :: optEnergy

    Call make1D_coords(nodes)
    Call make1D_nl(nodes, nl)

    nodesTemp = nodes
    nlTemp = nl

    print *,""
    print *,""
    print *,"coords start"
    Call calcEF_Nodes(nodes, nl)
    Call printEF_Nodes(nodes)
    print *, "Unoptimised ",nodes%potentialEnergy

    print *,""
    print *,""
    print *,"Simulated Annealing"
    Call run1D_SaOpt(nodes, optEnergy)

    print *,"SA result ",optEnergy
    Do i=1,nodes%count
      print *,"   ",i,nodes%nodeCoords(i)
    End Do

    nodes = nodesTemp
    nl = nlTemp
    print *,""
    print *,""
    print *,"Gradient Descent"
    Call run1D_GD(nodes, optEnergy)
    print *,"GD result ",optEnergy
    Do i=1,nodes%count
      print *,"   ",i,nodes%nodeCoords(i)
    End Do

    nodes = nodesTemp
    nl = nlTemp
    print *,""
    print *,""
    print *,"Newton Opt"
    Call run1D_NewtonOpt(nodes, optEnergy)
    print *,"NO result ",optEnergy
    Do i=1,nodes%count
      print *,"   ",i,nodes%nodeCoords(i)
    End Do


    nodes = nodesTemp
    nl = nlTemp
    print *,""
    print *,""
    print *,"coords CG opt"
    Call run1D_CG(nodes)
    Do i=1,nodes%count
      print *,"   ",i,nodes%nodeCoords(i)
    End Do



    nodes%nodeCoords(1) = 0.0D0
    nodes%nodeCoords(2) = 0.1D0
    nodes%nodeCoords(3) = 0.2D0
    nodes%nodeCoords(4) = 0.3D0
    nodes%nodeCoords(5) = 0.4D0
    nodes%nodeCoords(6) = 0.5D0
    nodes%nodeCoords(7) = 0.6D0
    nodes%nodeCoords(8) = 0.7D0
    nodes%nodeCoords(9) = 0.8D0
    nodes%nodeCoords(10) = 0.9D0
    Call make1D_nl(nodes, nl)
    Call calcEF_Nodes(nodes, nl)
    print *,""
    print *,""
    print *, "Real opt ",nodes%potentialEnergy

  End Subroutine testNodesOpt



  Subroutine testStaticCalc()
! Test static calculation interatomic potentials
    Implicit None ! Force declaration of all variables
! Vars Private
    Type(coordsUnitType), Allocatable, Dimension(:) :: coordsUnit
    Type(coordsType), Allocatable, Dimension(:) :: coords
    Type(nlType), Allocatable, Dimension(:) :: nl, nlB
    Type(potentialType) :: potential
    Type(oBulkProperty), Dimension(1:32) :: bpObj
    Real(kind=DoubleReal), Dimension(1:3) :: nodeA, nodeB
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: distortion
    Real(kind=DoubleReal) :: tempDp
    Real(kind=DoubleReal) :: rD
    Character(Len=16) :: inputStr
    Character(Len=2), Dimension(1:4) :: atomLabels
    Integer(kind=StandardInteger) :: i, j, cKey, atomKey
    Type(oDynamics) :: dynamicsSettings
! Start time
    Call setStartTime()
!----------------------------------
! Allocate arrays on HEAP
!----------------------------------
    Allocate(coordsUnit(1:p_confs))
    Allocate(coords(1:p_confs))
    Allocate(nl(1:p_confs))
    Allocate(nlB(1:p_confs))
! Load Vars
    Call loadVars()
    print *,envVars%cwd
    print *,envVars%user
    print *,envVars%home
!    Type(isotopesObj) :: isotopesList
!    Type(isotopeObj) :: isotopeA
! Load isotope data
    !Call loadIsotopes(isotopesList)
    !isotopeA = SearchIsotopes(26,28,isotopesList)
    !Call printIsotopeDetails(isotopeA)
    inputStr = "0.0D0"
    tempDp = StrToDp(inputStr)
    inputStr = "1.3D0"
    tempDp = StrToDp(inputStr)
! Init Page
    Call initPage(mainPage)
! Start test
    Call addLinePage("Static Calculation Testing","T")
! Set up coords
    Call initUnitCoords(coordsUnit)     ! geom.f90
    Call initCoords(coords)             ! geom.f90
    coordsUnit%aLat = 4.04D0
    coordsUnit%xCopy = 4
    coordsUnit%yCopy = 4
    coordsUnit%zCopy = 4
    atomLabels = "AL"
    Call standardCoords("FCC", coordsUnit,1,atomLabels)  ! geom.f90
    Call standardCoords("FCC", coordsUnit,2,atomLabels)  ! geom.f90
    Call expandUnitCoords(coordsUnit, coords) ! geom.f90
! Init potential
    Call addLinePage("Init potential")
    Call initPotential(potential)        ! potentials.f90
! Load potential
    Call addLinePage("Load potential")
    Call loadPotential("test.pot", potential)  ! potentials.f90
! Assign atom IDs
    Call addLinePage("Assign IDs")
    Call atomLabelIDs(potential, coords)       ! staticCalcs.f90
! Print potential summary
    Call printPotentialSummary(potential)      ! potentials.f90
! Print
    Call printAtomLabelIDs(coords)

    j = 1


    If(j.eq.1)Then
      cKey = 1
      Call makeNL(nl, coords, 6.5D0, cKey)       ! geom.f90
      Call calcEFS(nl, potential, cKey)
      print *,nl(1)%totalEnergy,(nl(1)%totalEnergy/nl(1)%coordsLength)

      Call calcBP(bpObj,potential)
    End If



    If(j.eq.2)Then
      cKey = 1
      atomKey = 1
      Call fitStandardPotentials(potential)
      Call HeatCoords(coords, 0.01D0, cKey)
      !coords(cKey)%coords(atomKey,1) = 0.001D0
      coords(cKey)%coordsFixed = 0
      !coords(cKey)%coordsFixed = 1
      !coords(cKey)%coordsFixed(1) = 0
      dynamicsSettings%adaptiveTimestep = .false.
      dynamicsSettings%damp = .true.
      Call RunDynamics(dynamicsSettings, coords, potential, cKey)


    End If


    If(j.eq.3)Then
! Distorion
      distortion(1,1) = 1.1D0
      distortion(1,2) = 0.02D0
      distortion(1,3) = 0.007D0
      distortion(2,1) = 0.003D0
      distortion(2,2) = 0.98D0
      distortion(2,3) = -0.001D0
      distortion(3,1) = -0.0079D0
      distortion(3,2) = 0.0021D0
      distortion(3,3) = 0.99D0
! Node A
      nodeA(1) = 1.6D0
      nodeA(2) = 3.2D0
      nodeA(3) = 3.2D0
! Node B
      nodeB(1) = 11.2D0
      nodeB(2) = 6.4D0
      nodeB(3) = 4.8D0
! RD
      rD = RdCoords(nodeA, nodeB, distortion)
      print *,rD
    End If






    If(j.eq.12)Then

    Call fitStandardPotentials(potential)
! Output Potential
    !Call outputPotential(potential,envVars%cwd)
! Build neighbour list and make keys
    Call addLinePage("Build neighbour list")
    Call makeNL(nl, coords, 6.5D0)       ! geom.f90
    Call nlPotentialKeys(nl, potential)

    !Call calcBP(bpObj,potential)

    Call HeatCoords(coords, 0.00D0, 1)
    Call makeNL(nl, coords, 6.5D0, 1)       ! geom.f90
    Call calcEFS(nl, potential, 1)
    print *,"EFS ",nl(1)%totalEnergy,nl(1)%pairEnergy,nl(1)%embeddingEnergy
    Do i=1,5
      print *,nl(1)%electronDensity(i,1),nl(1)%atomEnergy(i,1),nl(1)%atomEnergy(i,2),nl(1)%atomEnergy(i,3)
    End Do
    print *,""

    Call HeatCoords(coords, 0.01D0, 1)
    Call makeNL(nlB, coords, 6.5D0, 1)       ! geom.f90
    Call calcEFS(nlB, potential, 1)
    print *,"EFS ",nlB(1)%totalEnergy,nlB(1)%pairEnergy,nlB(1)%embeddingEnergy
    Do i=1,5
      print *,nlB(1)%electronDensity(i,1),nlB(1)%atomEnergy(i,1),nlB(1)%atomEnergy(i,2),nlB(1)%atomEnergy(i,3)
    End Do
    print *,""

    !Call getTime(theTime)
    !Print *,theTime
    Do i=1,30
      !Call calcEFS(nlB, potential, 1)
    End Do
    !Call getTime(theTime)
    !Print *,theTime

    print *,"--------------"
    !Call optGeom(coords, potential)

    End If



!----------------------------------
! Deallocate arrays on HEAP
!----------------------------------
    Deallocate(nl)
    Deallocate(nlB)
    Deallocate(coords)
    Deallocate(coordsUnit)

  End Subroutine testStaticCalc


  Subroutine testPotentialFitting()
! test potential fitting subroutine
    Implicit None ! Force declaration of all variables
! Vars:  Private
    Call startPotentialFit()
  End Subroutine testPotentialFitting



  Subroutine testCombinations()
    Implicit None ! Force declaration of all variables
  ! Vars Private
    Call combinationsPrint(4,3)
    print *,""
    print *,""
    Call combinationsPrint(5,4)
    print *,""
    print *,""
  End Subroutine testCombinations



  Subroutine testLinearRegression()
! Test static calculation interatomic potentials
    Implicit None ! Force declaration of all variables
! Vars Private
    Real(kind=DoubleReal), Dimension(1:4,1:2) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:2) :: parameters

    dataPoints(1,1) = 1.0D0
    dataPoints(1,2) = 13.8D0
    dataPoints(2,1) = 2.0D0
    dataPoints(2,2) = 18.9D0
    dataPoints(3,1) = 3.0D0
    dataPoints(3,2) = 22.6781D0
    dataPoints(4,1) = 4.0D0
    dataPoints(4,2) = 25.2D0

    parameters = BestFitLine(dataPoints)
    print *,"Best fit test: "
    print *,parameters(1),parameters(2)

    dataPoints(1,1) = 1.0D0
    dataPoints(1,2) = 1.58907D0
    dataPoints(2,1) = 2.0D0
    dataPoints(2,2) = 0.78911D0
    dataPoints(3,1) = 3.0D0
    dataPoints(3,2) = 0.39186D0
    dataPoints(4,1) = 8.0D0
    dataPoints(4,2) = 0.011833D0
    parameters = BestFitExp(dataPoints)
    print *,"Best fit exp test: "
    print *,parameters(1),parameters(2)

    dataPoints(1,1) = 1.0D0
    dataPoints(1,2) = 0.537325D0
    dataPoints(2,1) = 2.0D0
    dataPoints(2,2) = 0.10423D0
    dataPoints(3,1) = 3.0D0
    dataPoints(3,2) = 0.0202186D0
    dataPoints(4,1) = 8.0D0
    dataPoints(4,2) = 5.5531D-6
    parameters = BestFitExp(dataPoints)
    print *,"Best fit exp^2 test: "
    print *,parameters(1),parameters(2)


  End Subroutine testLinearRegression

  Subroutine testMorseFit()
! Test static calculation interatomic potentials
    Implicit None ! Force declaration of all variables
! Vars Private
    Real(kind=DoubleReal), Dimension(1:8,1:2) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:8,1:2) :: dataPointsDYDX
    Real(kind=DoubleReal), Dimension(1:3) :: parameters
    Real(kind=DoubleReal), Dimension(1:2) :: parametersLJ

    dataPoints(1,1) = 1.0D0
    dataPoints(1,2) = 43.9370D0
    dataPoints(2,1) = 1.6D0
    dataPoints(2,2) = 8.9981D0
    dataPoints(3,1) = 2.0D0
    dataPoints(3,2) = 2.6781D0
    dataPoints(4,1) = 2.6D0
    dataPoints(4,2) = 0.080557D0
    dataPoints(5,1) = 3.0D0
    dataPoints(5,2) = -0.23856D0
    dataPoints(6,1) = 3.4D0
    dataPoints(6,2) = -0.26361D0
    dataPoints(7,1) = 5.0D0
    dataPoints(7,2) = -0.066057D0
    dataPoints(8,1) = 5.6D0
    dataPoints(8,2) = -0.0339987D0

    dataPointsDYDX(1,1) = 1.0D0
    dataPointsDYDX(1,2) = -111.0190D0
    dataPointsDYDX(2,1) = 1.6D0
    dataPointsDYDX(2,2) = -25.2747D0
    dataPointsDYDX(3,1) = 2.0D0
    dataPointsDYDX(3,2) = -8.94672D0
    dataPointsDYDX(4,1) = 2.6D0
    dataPointsDYDX(4,2) = -1.53451D0
    dataPointsDYDX(5,1) = 3.0D0
    dataPointsDYDX(5,2) = -0.28964D0
    dataPointsDYDX(6,1) = 3.4D0
    dataPointsDYDX(6,2) = 0.083474D0
    dataPointsDYDX(7,1) = 5.0D0
    dataPointsDYDX(7,2) = 0.07155D0
    dataPointsDYDX(8,1) = 5.6D0
    dataPointsDYDX(8,2) = 0.038265D0

    parameters = MorseFit(dataPoints)
    print *,parameters(1),parameters(2),parameters(3)


    parameters = MorseFit_Extended(dataPoints, dataPointsDYDX)
    print *,parameters(1),parameters(2),parameters(3)


    parametersLJ = LJFit(dataPoints)
    print *,parametersLJ(1),parametersLJ(2)



  End Subroutine testMorseFit




  Subroutine testRandomDistribution()
! Test static calculation interatomic potentials
    Implicit None ! Force declaration of all variables
! Vars Private
    Integer(kind=LongInteger), Dimension(1:50,1:2) :: pointsBin
    Real(kind=DoubleReal) :: randResult = 0.0D0
    Real(kind=DoubleReal) :: difference
    Real(kind=DoubleReal) :: x, y
    Integer(kind=LongInteger) :: i, n, k, t, choice
    Type(oRandDist) :: randObj
! Init bin
    pointsBin = 0
! Set choice of distribution
    choice = 6
    n = 1000000
! bin size
    t = size(pointsBin,1)

    If(choice.eq.1)Then
! Flat distribution
! Set object
      randObj%distType = "FLAT"
      randObj%lowerBound = 0.0D0
      randObj%upperBound = 1.0D0
! Set difference
      difference = randObj%upperBound - randObj%lowerBound
! Loop through
      Do i=1,n
        randResult = RandDist(randObj)
        If(randResult.eq.randObj%lowerBound)Then
          k = 1
        Else
          k = ceiling((randResult-randObj%lowerBound)*(t/difference))
        End If
        If(k.lt.1)Then
          k = 1
        End If
        If(k.gt.100)Then
          k = 100
        End If
        pointsBin(k,1) = k
        pointsBin(k,2) = pointsBin(k,2) + 1
      End Do
    End If

    If(choice.eq.2)Then
! Set object
      randObj%distType = "GAUSS"
      randObj%lowerBound = -4.0D0
      randObj%upperBound = 4.0D0
      randObj%parameters(1) = 1.0D0  ! Sigma
      randObj%parameters(2) = 0.0D0  ! Mu
! Set difference
      difference = randObj%upperBound - randObj%lowerBound
! Loop through
      Do i=1,n
        randResult = RandDist(randObj)
        If(randResult.eq.randObj%lowerBound)Then
          k = 1
        Else
          k = ceiling((randResult-randObj%lowerBound)*(t/difference))
        End If
        If(k.lt.1)Then
          k = 1
        End If
        If(k.gt.100)Then
          k = 100
        End If
        pointsBin(k,1) = k
        pointsBin(k,2) = pointsBin(k,2) + 1
      End Do
    End If

    If(choice.eq.3)Then
! Set object
      randObj%distType = "SAMPLE"
      randObj%lowerBound = 0.0D0
      randObj%upperBound = 4.0D0
      print *,"sample"
! Set difference
      difference = randObj%upperBound - randObj%lowerBound
! Loop through
      Do i=1,n
        randResult = RandDist(randObj)
        If(randResult.eq.randObj%lowerBound)Then
          k = 1
        Else
          k = ceiling((randResult-randObj%lowerBound)*(t/difference))
        End If
        If(k.lt.1)Then
          k = 1
        End If
        If(k.gt.100)Then
          k = 100
        End If
        pointsBin(k,1) = k
        pointsBin(k,2) = pointsBin(k,2) + 1
      End Do
    End If

    If(choice.eq.4)Then
! Set object
      randObj%distType = "MAXBOLTZ"
      randObj%lowerBound = 0.0D0
      randObj%upperBound = 4.0D0
      randObj%parameters(1) = 1.0D0  ! a
      print *,"sample"
! Set difference
      difference = randObj%upperBound - randObj%lowerBound
! Loop through
      Do i=1,n
        randResult = RandDist(randObj)
        If(randResult.eq.randObj%lowerBound)Then
          k = 1
        Else
          k = ceiling((randResult-randObj%lowerBound)*(t/difference))
        End If
        If(k.lt.1)Then
          k = 1
        End If
        If(k.gt.100)Then
          k = 100
        End If
        pointsBin(k,1) = k
        pointsBin(k,2) = pointsBin(k,2) + 1
      End Do
    End If

    If(choice.eq.5)Then
! Set object
      randObj%distType = "DBGAUSS"
      randObj%lowerBound = 0.0D0
      randObj%upperBound = 4.0D0
      randObj%parameters(1) = 1.0D0  ! Factor
      randObj%parameters(2) = 1.0D0  ! Sigma
      randObj%parameters(3) = 0.0D0  ! Mu
      randObj%parameters(4) = 1.5D0  ! Factor
      randObj%parameters(5) = 1.0D0  ! Sigma
      randObj%parameters(6) = 3.0D0  ! Mu
! Set difference
      difference = randObj%upperBound - randObj%lowerBound
! Loop through
      Do i=1,n
        randResult = RandDist(randObj)
        If(randResult.eq.randObj%lowerBound)Then
          k = 1
        Else
          k = ceiling((randResult-randObj%lowerBound)*(t/difference))
        End If
        If(k.lt.1)Then
          k = 1
        End If
        If(k.gt.100)Then
          k = 100
        End If
        pointsBin(k,1) = k
        pointsBin(k,2) = pointsBin(k,2) + 1
      End Do
    End If

    If(choice.eq.6)Then
! Set object
      randObj%distType = "GHEAT"
      randObj%lowerBound = -1.0D0
      randObj%upperBound = 1.0D0
! Set difference
      difference = randObj%upperBound - randObj%lowerBound
! Loop through
      Do i=1,n
        randResult = RandDist(randObj)
        If(randResult.eq.randObj%lowerBound)Then
          k = 1
        Else
          k = ceiling((randResult-randObj%lowerBound)*(t/difference))
        End If
        If(k.lt.1)Then
          k = 1
        End If
        If(k.gt.100)Then
          k = 100
        End If
        pointsBin(k,1) = k
        pointsBin(k,2) = pointsBin(k,2) + 1
      End Do
    End If




! output
    Do i=1,t
      x = randObj%lowerBound+(i/(t*1.0D0))*&
        (randObj%upperBound-randObj%lowerBound)
      y = (pointsBin(i,2))/(1.0D0*n)
      print *,x,y
    End Do

    !randResult = RandDistF_Fx(RandDistF_Gaussian, randObj)

  End Subroutine testRandomDistribution


End Module testMod






















!
