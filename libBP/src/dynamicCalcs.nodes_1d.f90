! ------------------------------------------------------------
!           DYNAMICSCALCS: 1D Nodes
! ------------------------------------------------------------

! No practical use really, just a very simple model to use for testing and applying ideas to 3D MD
!
!


  Subroutine make1D_coords(nodes)
! Make 1D Node coords
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(o1D_Nodes) :: nodes
! Vars:  Private
! Number of nodes
    nodes%count = 10
! Allocate node arrays
    If(.not.Allocated(nodes%nodeCoords))Then
      Allocate(nodes%nodeCoords(1:nodes%count))
      Allocate(nodes%nodeMass(1:nodes%count))
      Allocate(nodes%nodeFixed(1:nodes%count))
      Allocate(nodes%force(1:nodes%count))
      Allocate(nodes%velocityH(1:nodes%count))
      Allocate(nodes%velocity(1:nodes%count))
      Allocate(nodes%acceleration(1:nodes%count))
      Allocate(nodes%optStep(1:nodes%count))
    End If
! Set alat
    !nodes%aLat = 16.265D0
    nodes%aLat = 32.53D0
! Set coords
    !nodes%nodeCoords(1) = 0.0D0
    !nodes%nodeCoords(2) = 0.22D0
    !nodes%nodeCoords(3) = 0.39D0
    !nodes%nodeCoords(4) = 0.615D0
    !nodes%nodeCoords(5) = 0.8D0
    nodes%nodeCoords(1) = 0.0211D0
    nodes%nodeCoords(2) = 0.0982D0
    nodes%nodeCoords(3) = 0.2081D0
    nodes%nodeCoords(4) = 0.3122D0
    nodes%nodeCoords(5) = 0.3983D0
    nodes%nodeCoords(6) = 0.4988D0
    nodes%nodeCoords(7) = 0.6122D0
    nodes%nodeCoords(8) = 0.6895D0
    nodes%nodeCoords(9) = 0.7995D0
    nodes%nodeCoords(10) = 0.9022D0
! Set Mass
    nodes%nodeMass = 1
! fixed nodes (none)
    nodes%nodeFixed = 0
  End Subroutine make1D_coords


  Subroutine make1D_nl(nodes, nl)
! Make 1D Node coords
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(o1D_Nodes) :: nodes
    Type(o1D_NodesNL) :: nl
    Integer(kind=StandardInteger) :: i, nodeA, nodeB, nlKey
    Real(kind=DoubleReal) :: rD, aPos, bPos
! Allocate nl arrays
    If(.not.Allocated(nl%nodeA))Then
      Allocate(nl%nodeA(1:100))
      Allocate(nl%nodeB(1:100))
      Allocate(nl%vecAB(1:100))
      Allocate(nl%RD(1:100))
    End If
! Keep within 0,1
    Do nodeA = 1,nodes%count
      nodes%nodeCoords(nodeA) = Modulus(nodes%nodeCoords(nodeA),1.0D0)
    End Do
! Build Neighbour List
    nlKey = 0
    nl%rCutoff = 7.50D0
    Do nodeA = 1,nodes%count
      Do nodeB = nodeA,nodes%count
        Do i=-1,1  ! ghost "cell"
          aPos = nodes%aLat*nodes%nodeCoords(nodeA)
          bPos = nodes%aLat*(nodes%nodeCoords(nodeB)+i)
          rD = sqrt((aPos-bPos)**2)
          If(rD.le.nl%rCutoff)Then
            If(nodeA.lt.nodeB)Then
              nlKey = nlKey + 1
              nl%nodeA(nlKey) = nodeA
              nl%nodeB(nlKey) = nodeB
              nl%RD(nlKey) = rD
              nl%vecAB(nlKey) = (aPos-bPos)/rD
            End If
          End If
        End Do
      End Do
    End Do
    nl%length = nlKey
  End Subroutine make1D_nl


  Subroutine calcEF_Nodes(nodes, nl)
! 1D Nodes Example
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(o1D_Nodes) :: nodes
    Type(o1D_NodesNL) :: nl
! Vars:  Private
    Integer(kind=StandardInteger) :: nlKey, cKey
    Real(kind=DoubleReal), Dimension(1:3) :: params, yArray
    Real(kind=DoubleReal) :: force
! Calculate forces
    nodes%force = 0.0D0
    nodes%potentialEnergy = 0.0D0
    Do nlKey = 1,nl%length
      params(1) = 0.2703D0 ! De
      params(2) = 1.1646D0 ! Alpha
      params(3) = 3.253D0  ! rc
      yArray = F_MorseFull(params,nl%RD(nlKey))    ! Use a Morse potential as an example
      nodes%potentialEnergy = nodes%potentialEnergy + yArray(1)
      force = nl%vecAB(nlKey)*yArray(2)
      nodes%force(nl%nodeA(nlKey)) = nodes%force(nl%nodeA(nlKey)) + force
      nodes%force(nl%nodeB(nlKey)) = nodes%force(nl%nodeB(nlKey)) - force
    End Do
    Do cKey = 1,nodes%count
      nodes%acceleration(cKey) = nodes%force(cKey)/nodes%nodeMass(cKey)
    End Do
  End Subroutine calcEF_Nodes



  Subroutine printEF_Nodes(nodes)
! 1D Nodes Example
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(o1D_Nodes) :: nodes
! Vars:  Private
    Integer(kind=StandardInteger) :: cKey
! Calculate forces
    print *,"Energy:  ",nodes%potentialEnergy
    print *,"Forces: "
    Do cKey = 1,nodes%count
      print *,cKey,nodes%force(cKey)
    End Do
  End Subroutine printEF_Nodes


  Subroutine run1D_SaOpt(nodes, fx_opt)
! 1D simulated annealing
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(o1D_Nodes) :: nodes
    Real(kind=DoubleReal) :: fx_opt
! Vars:  Private
    Type(o1D_NodesNL) :: nl
    Type(o1D_Nodes) :: nodesTemp, nodesOpt
    Integer(kind=StandardInteger) :: saLoops, i, j
    Real(kind=DoubleReal) :: optEnergy, energy, testEnergy, saTemp, varyMax, alpha, bestAlpha
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: direction
!    Real(kind=DoubleReal) ::
    Integer(kind=StandardInteger) :: accept
    Logical :: improved
! Number of nodes
    nodesTemp%count = nodes%count
! Allocate node arrays
    Allocate(nodesTemp%nodeCoords(1:nodes%count))
    Allocate(nodesTemp%nodeMass(1:nodes%count))
    Allocate(nodesTemp%nodeFixed(1:nodes%count))
    Allocate(nodesTemp%force(1:nodes%count))
    Allocate(nodesTemp%velocityH(1:nodes%count))
    Allocate(nodesTemp%velocity(1:nodes%count))
    Allocate(nodesTemp%acceleration(1:nodes%count))
    Allocate(nodesTemp%optStep(1:nodes%count))
!
    Allocate(direction(1:nodes%count))
! Initialise objects
    saLoops = 50
    Call make1D_nl(nodes, nl)
    Call calcEF_Nodes(nodes, nl)
    energy = nodes%potentialEnergy     ! set starting energy
    optEnergy = nodes%potentialEnergy  ! set opt energy as the starting energy
    nodesOpt = nodes
    print *,"SA start ",optEnergy
! Find maximum step size
    Do i=1,nodes%count
      direction(i) = nodes%force(i)/abs(nodes%force(i))
    End Do
    alpha = 1.0D0/(2.0D0*nodes%count)
    bestAlpha = alpha
    nodesTemp = nodes
    improved = .false.
    Do i=1,8
      Do j=1,nodes%count
        nodesTemp%nodeCoords(j) = nodes%nodeCoords(j)-alpha*direction(j)
        nodesTemp%nodeCoords(j) = Modulus(nodesTemp%nodeCoords(j),1.0D0)
      End Do
      Call make1D_nl(nodesTemp, nl)
      Call calcEF_Nodes(nodesTemp, nl)
      testEnergy = nodesTemp%potentialEnergy
      If(testEnergy.lt.optEnergy)Then
        bestAlpha = alpha
        optEnergy = testEnergy
        improved = .true.
      Else
        If(improved)Then
          Exit
        End If
      End If
      alpha = 0.2D0*alpha
    End Do
! Update nodes
    Do j=1,nodes%count
      nodes%nodeCoords(j) = nodes%nodeCoords(j)-bestAlpha*direction(j)
    End Do
    alpha = 2.0D0 * bestAlpha
    Do i=1,saLoops
! Calc
      nodesTemp = nodes
      varyMax = saVaryMax(i,saLoops,alpha,0.01D0*alpha)
      Call run1D_MoveNode(nodesTemp, varyMax, .true.)
      Call make1D_nl(nodesTemp, nl)
      Call calcEF_Nodes(nodesTemp, nl)
      testEnergy = nodesTemp%potentialEnergy
      saTemp = saLoopTemp(i,saLoops,0.001D0,0.00001D0)
      accept = saAccept(energy, testEnergy, saTemp)
      If(accept.gt.0)Then
        energy = testEnergy ! update energy
        nodes = nodesTemp   ! update node
        If(accept.eq.1)Then  ! If it was a better result than the optimum, save (only where accepted as a better result)
          If(energy.lt.optEnergy)Then
            optEnergy = energy
            nodesOpt = nodes
          End If
        End If
      End If
    End Do
    nodes = nodesOpt
    fx_opt = optEnergy
  End Subroutine run1D_SaOpt

  Subroutine run1D_MoveNode(nodes, maxVar, moveWithForce_In)
! 1D simulated annealing
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(o1D_Nodes) :: nodes
    Real(kind=DoubleReal) :: randNumber, maxVar
    Logical, Optional :: moveWithForce_In
! Vars:  Private
    Integer(kind=StandardInteger) :: i
    Logical :: moveWithForce
! Optional argument
    moveWithForce = .true. ! make the node mode in the direction of the force on it
    If(Present(moveWithForce_In))Then
      moveWithForce = moveWithForce_In
    End If
! Adjust
    Do i=1,nodes%count
      randNumber = RandomLCG()
      If(moveWithForce)Then
        nodes%nodeCoords(i) = nodes%nodeCoords(i) - &
          randNumber * maxVar * abs(nodes%force(i))/nodes%force(i)
      Else
        nodes%nodeCoords(i) = nodes%nodeCoords(i) + &
          (randNumber-0.5D0) * 2.0D0 * maxVar
      End If
      nodes%nodeCoords(i) = Modulus(nodes%nodeCoords(i),1.0D0)
    End Do
  End Subroutine run1D_MoveNode


  Subroutine run1D_GD(nodes, fx_opt)
! 1D simulated annealing
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(o1D_Nodes) :: nodes
    Real(kind=DoubleReal) :: fx_opt
! Vars:  Private
    Type(o1D_NodesNL) :: nl
    Type(o1D_Nodes) :: nodesTemp
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: gradVec
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: weightVec
    Real(kind=DoubleReal) :: alpha, bestAlpha
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: direction
    Integer(kind=StandardInteger) :: i, n, j
    Real(kind=DoubleReal) :: h, hsq, fx0, fx0_t, fx_min, maxGrad, forceMag
    Logical :: weights
    Logical :: improved
! Allocate opt arrays
    Allocate(gradVec(1:1*nodes%count))
    Allocate(weightVec(1:1*nodes%count))
    Allocate(direction(1:1*nodes%count))
! Init vars
    h = 1.0D-6
    hsq = 1.0D-12
    alpha = 0.1D0
    weights = .false.
! Iterate and find a good value for alpha
    Do n=1,8  ! Optimisation loops
      Call make1D_nl(nodes, nl)
      Call calcEF_Nodes(nodes, nl)
      fx0 = nodes%potentialEnergy
      nodesTemp = nodes
      maxGrad = 0.0D0
      Do i=1,nodes%count
        forceMag = abs(nodes%force(i))
        direction(i) = nodes%force(i)/forceMag
        If(forceMag.gt.maxGrad)Then
          maxGrad = forceMag
        End If
        weightVec(i) = forceMag
      End Do
      Do i=1,nodes%count
        weightVec(i) = weightVec(i) / maxGrad
      End Do
      alpha = 10.0D0 * alpha
      If(alpha.lt.1.0D-4)Then
        weights = .true.
        alpha = 1.0D0
      End If
! Line search
      improved = .false.
      Do i=1,8
        nodesTemp = nodes
        Do j=1,nodes%count
          If(weights)Then
            nodesTemp%nodeCoords(j) = nodes%nodeCoords(j)-alpha*weightVec(j)*direction(j)
          Else
            nodesTemp%nodeCoords(j) = nodes%nodeCoords(j)-alpha*direction(j)
          End If
        End Do
        If(i.eq.1)Then
          bestAlpha = alpha
          fx_min = fx0
        End If
        Call make1D_nl(nodesTemp, nl)
        Call calcEF_Nodes(nodesTemp, nl)
        fx0_t = nodesTemp%potentialEnergy
        If(fx0_t.lt.fx_min)Then
          bestAlpha = alpha
          fx_min = fx0_t
          improved = .true.
        Else
          If(improved)Then
            Exit
          End If
        End If
        alpha = alpha / 10.0D0
      End Do
      Do j=1,nodes%count
        nodes%nodeCoords(j) = nodes%nodeCoords(j)-bestAlpha*direction(j)
      End Do
    End Do
    fx_opt = fx_min
  End Subroutine run1D_GD


  Subroutine run1D_NewtonOpt(nodes, fx_opt)
! 1D simulated annealing
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(o1D_Nodes) :: nodes
    Real(kind=DoubleReal) :: fx_opt
! Vars:  Private
    Type(o1D_NodesNL) :: nl
    Type(o1D_Nodes) :: nodesTemp
    Integer(kind=StandardInteger) :: optLoop, i, j
    Real(kind=DoubleReal) :: optEnergy
    Real(kind=DoubleReal) :: fx0
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: xVec
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: eVec
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: gradVec
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: hessian
! Number of nodes
    nodesTemp%count = nodes%count
! Allocate temp node arrays
    Allocate(nodesTemp%nodeCoords(1:nodes%count))
    Allocate(nodesTemp%nodeMass(1:nodes%count))
    Allocate(nodesTemp%nodeFixed(1:nodes%count))
    Allocate(nodesTemp%force(1:nodes%count))
    Allocate(nodesTemp%velocityH(1:nodes%count))
    Allocate(nodesTemp%velocity(1:nodes%count))
    Allocate(nodesTemp%acceleration(1:nodes%count))
    Allocate(nodesTemp%optStep(1:nodes%count))
! Allocate opt arrays
    Allocate(xVec(1:1*nodes%count))
    Allocate(eVec(1:1*nodes%count))
    Allocate(gradVec(1:1*nodes%count))
    Allocate(hessian(1:1*nodes%count,1:1*nodes%count))
! Initialise objects
    Call make1D_nl(nodes, nl)
    Call calcEF_Nodes(nodes, nl)
    optEnergy = nodes%potentialEnergy
    Do optLoop=1,50
! Evaluate current node positions
      Call make1D_nl(nodes, nl)
      Call calcEF_Nodes(nodes, nl)
      fx0 = nodes%potentialEnergy
! Gradient
      Do i=1,nodes%count
        Do j=1,1  ! 1,3 for 3D
          gradVec(i) = nodes%force(i)
        End Do
      End Do
! Build Hessian
      Call run1D_Hessian(nodes,hessian)
! Calculate change vector
      eVec = SolveLinearSet(hessian, gradVec)
      eVec = -1.0D0 * eVec
      Do i=1,nodes%count
        nodes%nodeCoords(i) = nodes%nodeCoords(i)+eVec(i)
      End Do
    End Do
    Call make1D_nl(nodes, nl)
    Call calcEF_Nodes(nodes, nl)
    fx_opt = nodes%potentialEnergy
  End Subroutine run1D_NewtonOpt


  Subroutine run1D_Gradient(nodes,gradient)
! 1D simulated annealing
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(o1D_Nodes) :: nodes
    Real(kind=DoubleReal), Dimension(:) :: gradient
! Vars:  Private
    Type(o1D_NodesNL) :: nl
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal) :: fx0
    Real(kind=DoubleReal) :: h, fx_f
! Init
    h = 1.0D-6
! Unperturbed
    Call make1D_nl(nodes, nl)
    Call calcEF_Nodes(nodes, nl)
    fx0 = nodes%potentialEnergy
! Use finite difference
    Do i=1,nodes%count
      nodes%nodeCoords(i) = nodes%nodeCoords(i) + h
      Call make1D_nl(nodes, nl)
      Call calcEF_Nodes(nodes, nl)
      nodes%nodeCoords(i) = nodes%nodeCoords(i) - h
      fx_f = nodes%potentialEnergy
      gradient(i) = (fx_f-fx0)/h
    End Do
  End Subroutine run1D_Gradient


  Subroutine run1D_Hessian(nodes,hessian)
! 1D simulated annealing
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(o1D_Nodes) :: nodes
    Real(kind=DoubleReal), Dimension(:,:) :: hessian
! Vars:  Private
    Type(o1D_Nodes) :: nodesTemp
    Type(o1D_NodesNL) :: nl
    Integer(kind=StandardInteger) :: i, j
    Real(kind=DoubleReal) :: fx0
    Real(kind=DoubleReal) :: h, hsq, fx_ff, fx_fb, fx_bf, fx_bb, coordI, coordJ
! Init
    h = 1.0D-6
    hsq = 1.0D-12
! Unperturbed
    Call make1D_nl(nodes, nl)
    Call calcEF_Nodes(nodes, nl)
    fx0 = nodes%potentialEnergy
    nodesTemp = nodes
! Use finite difference - diagonal of hessian
    Do i=1,nodes%count
      coordI = nodesTemp%nodeCoords(i)
! calc fx_bb
      nodesTemp%nodeCoords(i) = nodesTemp%nodeCoords(i) - h
      Call make1D_nl(nodesTemp, nl)
      Call calcEF_Nodes(nodesTemp, nl)
      fx_bb = nodesTemp%potentialEnergy
! calc fx_ff
      nodesTemp%nodeCoords(i) = nodesTemp%nodeCoords(i) + h + h
      Call make1D_nl(nodesTemp, nl)
      Call calcEF_Nodes(nodesTemp, nl)
      fx_ff = nodesTemp%potentialEnergy
! reset
      nodesTemp%nodeCoords(i) = coordI
! Store
      hessian(i,i) = (fx_ff+fx_bb-fx0-fx0)/hsq
    End Do
! Upper/lower
    Do i=1,nodes%count-1
      Do j=i,nodes%count
        coordI = nodesTemp%nodeCoords(i)
        coordJ = nodesTemp%nodeCoords(j)
! FF
        nodesTemp%nodeCoords(i) = nodesTemp%nodeCoords(i) + h
        nodesTemp%nodeCoords(j) = nodesTemp%nodeCoords(j) + h
        Call make1D_nl(nodesTemp, nl)
        Call calcEF_Nodes(nodesTemp, nl)
        fx_ff = nodesTemp%potentialEnergy
! FB
        nodesTemp%nodeCoords(j) = nodesTemp%nodeCoords(j) - h - h
        Call make1D_nl(nodesTemp, nl)
        Call calcEF_Nodes(nodesTemp, nl)
        fx_fb = nodesTemp%potentialEnergy
! BB
        nodesTemp%nodeCoords(i) = nodesTemp%nodeCoords(i) - h - h
        Call make1D_nl(nodesTemp, nl)
        Call calcEF_Nodes(nodesTemp, nl)
        fx_bb = nodesTemp%potentialEnergy
! BF
        nodesTemp%nodeCoords(j) = nodesTemp%nodeCoords(j) + h + h
        Call make1D_nl(nodesTemp, nl)
        Call calcEF_Nodes(nodesTemp, nl)
        fx_bf = nodesTemp%potentialEnergy
  ! reset
        nodesTemp%nodeCoords(i) = coordI
        nodesTemp%nodeCoords(j) = coordJ
! Store
        hessian(i,j) = (fx_ff+fx_bb-fx_fb-fx_bf)/(4.0D0*hsq)
        hessian(j,i) = hessian(i,j)
      End Do
    End Do
  End Subroutine run1D_Hessian


  Subroutine run1D_CG(nodes)
! 1D simulated annealing
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(o1D_Nodes) :: nodes
! Vars:  Private
    Type(o1D_NodesNL) :: nl
    Type(o1D_Nodes) :: nodesTemp
    Integer(kind=StandardInteger) :: optLoop, i, j
    Real(kind=DoubleReal) :: optEnergy
    Real(kind=DoubleReal) :: fx0
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: xVec
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: eVec
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: gradVec
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: hessian
! Number of nodes
    nodesTemp%count = nodes%count
! Allocate temp node arrays
    Allocate(nodesTemp%nodeCoords(1:nodes%count))
    Allocate(nodesTemp%nodeMass(1:nodes%count))
    Allocate(nodesTemp%nodeFixed(1:nodes%count))
    Allocate(nodesTemp%force(1:nodes%count))
    Allocate(nodesTemp%velocityH(1:nodes%count))
    Allocate(nodesTemp%velocity(1:nodes%count))
    Allocate(nodesTemp%acceleration(1:nodes%count))
    Allocate(nodesTemp%optStep(1:nodes%count))
! Allocate opt arrays
    Allocate(xVec(1:1*nodes%count))
    Allocate(eVec(1:1*nodes%count))
    Allocate(gradVec(1:1*nodes%count))
    Allocate(hessian(1:1*nodes%count,1:1*nodes%count))
! Initialise objects

    Call make1D_nl(nodes, nl)
    Call calcEF_Nodes(nodes, nl)
    optEnergy = nodes%potentialEnergy
! Build gradients
    Call run1D_Gradient(nodes,gradVec)
! Build Hessian
    Call run1D_Hessian(nodes,hessian)





    print *,optLoop,"Opt Energy ",optEnergy
    Do optLoop=1,50
! Evaluate current node positions
      Call make1D_nl(nodes, nl)
      Call calcEF_Nodes(nodes, nl)
      fx0 = nodes%potentialEnergy
      !print *,optLoop,nodes%potentialEnergy
! Gradient
      Do i=1,nodes%count
        Do j=1,1  ! 1,3 for 3D
          gradVec(i) = nodes%force(i)
        End Do
      End Do
! Calculate change vector
      eVec = SolveLinearSet(hessian, gradVec)
      eVec = -1.0D0 * eVec
      Do i=1,nodes%count
        nodes%nodeCoords(i) = nodes%nodeCoords(i)+eVec(i)
      End Do
    End Do
    Call make1D_nl(nodes, nl)
    Call calcEF_Nodes(nodes, nl)
    print *,nodes%potentialEnergy
  End Subroutine run1D_CG



! Functions

  Function NodeU(r) Result (ur)
! 1D Nodes Example
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal) :: r
! Vars:  Out
    Real(kind=DoubleReal) :: ur
! Vars:  Private
! Calculate
    ur = ((r-2.0D0)**4-10*r**2)*exp(-1*r**2)
  End Function

  Function NodeDU(r) Result (dudr)
! 1D Nodes Example
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal) :: r
! Vars:  Out
    Real(kind=DoubleReal) :: dudr
! Vars:  Private
    Real(kind=DoubleReal) :: urA, urB
! Calculate
    urA = NodeU(r)
    urB = NodeU(r+1.0D-5)
    dudr = (urB-urA)/(1.0D-5)
  End Function
