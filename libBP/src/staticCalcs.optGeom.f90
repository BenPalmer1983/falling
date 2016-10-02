! -------------------------------------------------
!  Include File:   Optimise geometry
!
! -------------------------------------------------
  Subroutine optGeom(coordsGeom, potential, followGrad_In, gradWeighting_In, cKey_In)
! Optimise Geometry
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(coordsType), Allocatable, Dimension(:) :: coordsGeom
    Type(potentialType) :: potential
    Logical, Optional :: followGrad_In
    Logical, Optional :: gradWeighting_In
    Integer(kind=StandardInteger), Optional :: cKey_In
! Vars Private
    Type(nlType), Allocatable, Dimension(:) :: nlGeom
    Integer(kind=StandardInteger) :: cKey, cKey_s, cKey_e
    Integer(kind=StandardInteger) :: saLoops, saLoopsW, i, j, n
    Real(kind=DoubleReal) :: energy, optEnergy, testEnergy, saTemp, tStart, tEnd
    Real(kind=DoubleReal) :: alpha, bestAlpha, varyMax
    Real(kind=DoubleReal) :: randFloat, maxForce, absForce, lastImprovement
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: direction
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: gradWeight
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: coordsTemp
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: coordsOpt
    Real(kind=DoubleReal) :: rD, xD, yD, zD, rdMin
    Integer(kind=StandardInteger) :: accept
    Logical :: improved
    Logical :: followGrad
    Logical :: gradWeighting
! Optional Arguments
    cKey_s = 1
    cKey_e = size(coordsGeom,1)
    If(Present(cKey_In))Then
      cKey_s = cKey_In
      cKey_e = cKey_In
    End If
    followGrad = .false.
    If(Present(followGrad_In))Then
      followGrad = followGrad_In
    End If
    gradWeighting = .false.
    If(Present(gradWeighting_In))Then
      gradWeighting = gradWeighting_In
    End If
!----------------------------------
! Allocate arrays on HEAP
!----------------------------------
! NL
    Allocate(nlGeom(1:size(coordsGeom,1)))
! Loop through configs
    Do cKey = cKey_s,cKey_e
! Allocate direction array
      Allocate(direction(1:coordsGeom(cKey)%length,1:3))
      Allocate(gradWeight(1:coordsGeom(cKey)%length,1:3))
      Allocate(coordsTemp(1:coordsGeom(cKey)%length,1:3))
      Allocate(coordsOpt(1:coordsGeom(cKey)%length,1:3))
!
! 1: Find step size
!----------------------------------------------------------
! First evaluation
      Call makeNL(nlGeom, coordsGeom, 6.5D0, cKey)
      Call calcEFS(nlGeom, potential, cKey)
! Calc direction vector
      Do i=1,coordsGeom(cKey)%length
        Do j=1,3
          absForce = abs(nlGeom(cKey)%forces(i,j))
          direction(i,j) = nlGeom(cKey)%forces(i,j)/absForce
        End Do
      End Do
      optEnergy = nlGeom(cKey)%totalEnergy
! Init vars
      improved = .false.
      alpha = 0.2D0/((coordsGeom(cKey)%length)**(1.0D0/3.0D0))
      bestAlpha = alpha
      lastImprovement = 0.0D0
! Search for best alpha
      Do n=1,8
! Update positions
        Do i=1,coordsGeom(cKey)%length
          Do j=1,3
            coordsTemp(i,j) = coordsGeom(cKey)%coords(i,j)
            coordsGeom(cKey)%coords(i,j) = coordsGeom(cKey)%coords(i,j)+alpha*direction(i,j)
          End Do
        End Do
! Perfom calc
        Call makeNL(nlGeom, coordsGeom, 6.5D0, cKey)
        Call calcEFS(nlGeom, potential, cKey)
! Reset coords
        Do i=1,coordsGeom(cKey)%length
          Do j=1,3
            coordsGeom(cKey)%coords(i,j) = coordsTemp(i,j)
          End Do
        End Do
! Check result
        testEnergy = nlGeom(cKey)%totalEnergy
        If(testEnergy.lt.optEnergy)Then
          lastImprovement = optEnergy - testEnergy
          improved = .true.
          bestAlpha = alpha
          optEnergy = testEnergy
          maxForce = 0.0D0
        Else
          If(improved)Then
            Exit
          End If
        End If
        alpha = 0.2D0 * alpha
      End Do
! Update with best alpha
      Do i=1,coordsGeom(cKey)%length
        Do j=1,3
          coordsGeom(cKey)%coords(i,j) = coordsGeom(cKey)%coords(i,j)+bestAlpha*direction(i,j)
          coordsOpt(i,j) = coordsGeom(cKey)%coords(i,j)
        End Do
      End Do
! Prepare new staring calc for SA
      Call makeNL(nlGeom, coordsGeom, 6.5D0, cKey)
      Call calcEFS(nlGeom, potential, cKey)
      maxForce = 0.0D0
      Do i=1,coordsGeom(cKey)%length
        Do j=1,3
          absForce = abs(nlGeom(cKey)%forces(i,j))
          direction(i,j) = nlGeom(cKey)%forces(i,j)/absForce
          If(absForce.gt.maxForce)Then
            maxForce = absForce
          End If
        End Do
      End Do
      If(gradWeighting)Then
        Do i=1,coordsGeom(cKey)%length
          Do j=1,3
            gradWeight(i,j) = abs(direction(i,j)/maxForce)
          End Do
        End Do
      End If
! Estimate start/end temperatures
      tStart = (-1.0D0*(lastImprovement))/log(0.3D0)  ! accept based on last error with 30% chance
      tEnd = (-1.0D0*(lastImprovement))/log(0.001D0)  ! accept based on last error with 0.1% chance
!----------------------------------------------------------
! Start SA
!----------------------------------------------------------
! Set alpha for SA algorithm
      alpha = 2.0D0*bestAlpha
!--------------------------------
! Initialise
!--------------------------------
      saLoops = 10
      saLoopsW = 1
      energy = optEnergy
      Do n=1,saLoops
        saTemp = saLoopTemp(n,saLoops,tStart,tEnd)
        varyMax = saVaryMax(n,saLoops,alpha,0.1D0*alpha)
! Vary coords
        Do i=1,coordsGeom(cKey)%length
          Do j=1,3
            randFloat = RandomLCG()
            coordsTemp(i,j) = coordsGeom(cKey)%coords(i,j)
            If(followGrad)Then
              If(gradWeighting.and.(n.ge.saLoopsW))Then
                coordsGeom(cKey)%coords(i,j) = coordsGeom(cKey)%coords(i,j)+&
                randFloat*varyMax*direction(i,j)*2.0D0*gradWeight(i,j)
              Else
                coordsGeom(cKey)%coords(i,j) = coordsGeom(cKey)%coords(i,j)+randFloat*varyMax*direction(i,j)
              End If
            Else
              coordsGeom(cKey)%coords(i,j) = coordsGeom(cKey)%coords(i,j)+(randFloat-0.5D0)*2.0D0*varyMax
            End If
          End Do
        End Do
        Call makeNL(nlGeom, coordsGeom, 6.5D0, cKey)
        Call calcEFS(nlGeom, potential, cKey)
        testEnergy = nlGeom(cKey)%totalEnergy
        accept = saAccept(energy, testEnergy, saTemp)
        If(accept.gt.0)Then
! update energy
          energy = testEnergy
! update direction vector
          maxForce = 0.0D0
          Do i=1,coordsGeom(cKey)%length
            Do j=1,3
              absForce = abs(nlGeom(cKey)%forces(i,j))
              direction(i,j) = nlGeom(cKey)%forces(i,j)/absForce
              If(absForce.gt.maxForce)Then
                maxForce = absForce
              End If
            End Do
          End Do
          If(gradWeighting)Then
            Do i=1,coordsGeom(cKey)%length
              Do j=1,3
                gradWeight(i,j) = abs(direction(i,j)/maxForce)
              End Do
            End Do
          End If
! If it was a better result than the optimum, save (only where accepted as a better result)
          If(accept.eq.1)Then
            If(energy.lt.optEnergy)Then
              optEnergy = energy
              Do i=1,coordsGeom(cKey)%length
                Do j=1,3
                  coordsOpt(i,j) = coordsGeom(cKey)%coords(i,j)
                End Do
              End Do
            End If
          End If
        Else
! Reverse step
          coordsGeom(cKey)%coords(i,j) = coordsTemp(i,j)
        End If
      End Do
! Align with origin and save optimum coords in coordGeom object
! find distance from 0,0,0 or 1,1,1
      n = 1
      rdMin = 3.0D0
      Do i=1,coordsGeom(cKey)%length
        rD = 0.0D0
        Do j=1,3
          If(coordsOpt(i,j).lt.0.5D0)Then
            rD = rD + coordsOpt(i,j)*coordsOpt(i,j)
          Else
            rD = rD + (1.0D0-coordsOpt(i,j))*(1.0D0-coordsOpt(i,j))
          End If
        End Do
        If(rD.lt.rdMin)Then
          n = i
          rdMin = rD
        End If
      End Do
      xD = coordsOpt(n,1)
      yD = coordsOpt(n,2)
      zD = coordsOpt(n,3)
      Do i=1,coordsGeom(cKey)%length
        coordsGeom(cKey)%coords(i,1) = Modulus(coordsOpt(i,1)-xD,1.0D0)
        coordsGeom(cKey)%coords(i,2) = Modulus(coordsOpt(i,2)-yD,1.0D0)
        coordsGeom(cKey)%coords(i,3) = Modulus(coordsOpt(i,3)-zD,1.0D0)
      End Do
! Deallocate direction array
      Deallocate(direction)
      Deallocate(coordsTemp)
    End Do
!----------------------------------
! Deallocate arrays on HEAP
!----------------------------------
    Deallocate(nlGeom)

  End Subroutine optGeom





  Subroutine optGeom_BuildMatrix(coordsGeom, potential, cKey)
! Optimise Geometry
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(coordsType), Allocatable, Dimension(:) :: coordsGeom
    Type(potentialType) :: potential
    Integer(kind=StandardInteger) :: cKey
! Vars Private
    Type(nlType_Opt), Allocatable, Dimension(:) :: nlGeom
    Integer(kind=StandardInteger) :: i, atomKey
    Real(kind=DoubleReal), Dimension(1:3) :: coordChange
    Real(kind=DoubleReal) :: theTime
    Type(oCalcSettings) :: calcSettings
!----------------------------------
! Allocate arrays on HEAP
!----------------------------------
    Allocate(nlGeom(1:1))
! Make neighbour list
    Call makeNL_opt(nlGeom, coordsGeom, 6.5D0, cKey)   ! geom.nl_opt.f90   - add option to set rVerlet, maybe in a data type

    atomKey = 1
    calcSettings%cKey = cKey
    calcSettings%atomKey = atomKey
    calcSettings%calcEnergy = .true.
    calcSettings%calcForces = .false.
    calcSettings%useRD_Moved = .false.

    print *,"Approx energy 0:  "
    Call calcE_Approx(nlGeom, potential, calcSettings)
    print *,nlGeom(1)%totalEnergyFixed,nlGeom(1)%totalEnergy
    print *,""
    print *,""

    Call getTime(theTime)
    Print *,theTime
    Do i=1,30
      atomKey = i
      calcSettings%cKey = cKey
      calcSettings%atomKey = atomKey
      calcSettings%calcEnergy = .true.
      calcSettings%calcForces = .false.
      calcSettings%useRD_Moved = .true.
      coordChange(1) = 0.005D0
      coordChange(2) = 0.0D0
      coordChange(3) = 0.0D0
      Call nlMoveAtom_Opt(nlGeom, cKey, atomKey, coordChange)
      Call calcE_Approx(nlGeom, potential, calcSettings)
    End Do
    Call getTime(theTime)
    Print *,theTime


! Loop through each atom and vary position slightly



!----------------------------------
! Dellocate arrays on HEAP
!----------------------------------
    Deallocate(nlGeom)
  End Subroutine optGeom_BuildMatrix








!
