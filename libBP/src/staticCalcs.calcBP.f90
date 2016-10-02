! -------------------------------------------------
!  Include File:   Calc bulk properties
!
! -------------------------------------------------
  Subroutine calcBP(bpObj,potential)
! Calculate bulk properties
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oBulkProperty), Dimension(:) :: bpObj
    Type(potentialType) :: potential
! Vars:  Private
! Run Calculations
!
! Equation of State
    Call calcEOS(bpObj,potential)
! Elastic Constants
    Call calcEC(bpObj,potential)
  End Subroutine calcBP



  Subroutine calcEOS(bpObj,potential,structure_In)
! Calculate bulk properties
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oBulkProperty), Dimension(:) :: bpObj
    Type(potentialType) :: potential
    Character(*), Optional :: structure_In
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, n, m, cKey
    Character(Len=3) :: structure
    Type(coordsUnitType), Allocatable, Dimension(:) :: coordsUnitBP
    Type(coordsType), Allocatable, Dimension(:) :: coordsBP
    Type(nlType), Allocatable, Dimension(:) :: nlBP
    Type(oElements) :: elementsList
    Type(oElement) :: atomData
    Character(Len=16), Dimension(1:4) :: atomLabels
    Integer(kind=StandardInteger), Dimension(1:4) :: atomIDs
! EOS
    Integer(kind=StandardInteger) :: loops, loopsMidPoint
    Real(kind=DoubleReal), Dimension(1:p_bpPoints,1:2) :: points
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients
    Real(kind=DoubleReal) :: startAlat, dAlat, estAlat, aLat, optAlat
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: distortionArray
    Real(kind=DoubleReal) :: V, E
! Optional Arguments
    structure = "   "
    If(Present(structure_In))Then
      structure = structure_In(1:3)
      structure = StrToUpper(structure)
    End If
!----------------------------------
! Allocate arrays on HEAP
!----------------------------------
    Allocate(coordsUnitBP(1:4))
    Allocate(coordsBP(1:4))
    Allocate(nlBP(1:4))
! Initialise objects
    Call initUnitCoords(coordsUnitBP)     ! geom.f90
    Call initCoords(coordsBP)             ! geom.f90
! Load a list of elements data
    Call loadElements(elementsList)
! Set loops
    loops = p_bpPoints
    loopsMidPoint = (loops + 1)/2
! Loop through atom IDs
    cKey = 0
    n = 0
    m = 0
!----------------------------------
! FCC Structure
!----------------------------------
    If((structure.eq."   ").or.(structure.eq."FCC"))Then
      Do i=1,potential%atomID_Count
        atomLabels = potential%atomIDs(i)  ! Set atomLabels array and atomIDs array
        atomIDs = i
! Look up atom details (if the label matches a value element symbol)
        atomData = SearchElements(potential%atomIDs(i), elementsList)
        dAlat = 0.05D0
!----------------------------------
! FCC Equation of State
!----------------------------------
        estAlat = LatticeParameter(atomData,4)
        bpObj(i)%fccAlat_Estimate = estAlat
        startAlat = estAlat * (1.0D0 + dAlat)
        cKey = cKey + 1
        coordsUnitBP%aLat = startAlat
        coordsUnitBP%xCopy = 4
        coordsUnitBP%yCopy = 4
        coordsUnitBP%zCopy = 4
        Call standardCoords("FCC", coordsUnitBP, cKey, atomLabels, atomIDs)  ! geom.f90
        Call expandUnitCoords(coordsUnitBP, coordsBP)   ! geom.f90
! Find e0 and aLat
!----------------------------------
        Call makeNL(nlBP, coordsBP, 6.5D0, cKey)       ! geom.f90
        Call nlPotentialKeys(nlBP, potential)    ! staticCalcs.keys.f90
        Call calcE(nlBP, potential, cKey)
        points(loops,1) = nlBP(cKey)%volume/nlBP(cKey)%coordsLength
        points(loops,2) = (nlBP(cKey)%totalEnergy/nlBP(cKey)%coordsLength)
        Do j=1,loops-1
          aLat = (estAlat*(1.0D0-(j-loopsMidPoint+1)*dAlat))
          distortionArray = IdentityMatrix(distortionArray)*(aLat/startAlat)
          Call distortNL(nlBP, distortionArray, cKey)       ! geom.f90
          Call calcE(nlBP, potential, cKey)                 ! staticCalcs.calcE.f90
          points(loops-j,1) = nlBP(cKey)%volumeDistorted/nlBP(cKey)%coordsLength
          points(loops-j,2) = (nlBP(cKey)%totalEnergy/nlBP(cKey)%coordsLength)
        End Do
        coefficients = BirchMurnFit(points)
        optAlat = (coefficients(2)*nlBP(cKey)%coordsLength)**(1.0D0/3.0D0)
! Store data
        bpObj(i)%fccAlat = (coefficients(2)*4)**(1.0D0/3.0D0)
        bpObj(i)%fccV0 = (coefficients(2))
        bpObj(i)%fccE0 = coefficients(1)
        bpObj(i)%fccB0 = coefficients(3)
        bpObj(i)%fccB0_GPA = UnitConvert(bpObj(i)%fccB0, "EVAN3", "GPA")   ! Convert to GPA
        bpObj(i)%fccBp0 = coefficients(4)
! Make EoS points
        Do j=1,15
          V = bpObj(i)%fccV0*(1.0D0-(8-j)*0.01D0)
          E = F_BirchMurn(V,coefficients)  ! scienceFunctions.materialsFunctions.f90
          bpObj(i)%fccEosPoints(j,1) = V
          bpObj(i)%fccEosPoints(j,2) = E
        End Do
      End Do
    End If
!----------------------------------
! BCC Structure
!----------------------------------
    If((structure.eq."   ").or.(structure.eq."BCC"))Then
      Do i=1,potential%atomID_Count
        atomLabels = potential%atomIDs(i)  ! Set atomLabels array and atomIDs array
        atomIDs = i
! Look up atom details (if the label matches a value element symbol)
        atomData = SearchElements(potential%atomIDs(i), elementsList)
        dAlat = 0.05D0
!----------------------------------
! BCC Equation of State
!----------------------------------
        estAlat = LatticeParameter(atomData,2)
        bpObj(i)%bccAlat_Estimate = estAlat
        startAlat = estAlat * (1.0D0 + dAlat)
        cKey = cKey + 1
        coordsUnitBP%aLat = startAlat
        coordsUnitBP%xCopy = 4
        coordsUnitBP%yCopy = 4
        coordsUnitBP%zCopy = 4
        Call standardCoords("BCC", coordsUnitBP, cKey, atomLabels, atomIDs)  ! geom.f90
        Call expandUnitCoords(coordsUnitBP, coordsBP)   ! geom.f90
! Find e0 and aLat
!----------------------------------
        Call makeNL(nlBP, coordsBP, 6.5D0, cKey)       ! geom.f90
        Call nlPotentialKeys(nlBP, potential)    ! staticCalcs.keys.f90
        Call calcE(nlBP, potential, cKey)        ! staticCalcs.calcE.f90
        points(loops,1) = nlBP(cKey)%volume/nlBP(cKey)%coordsLength
        points(loops,2) = (nlBP(cKey)%totalEnergy/nlBP(cKey)%coordsLength)
        Do j=1,loops-1
          aLat = (estAlat*(1.0D0-(j-loopsMidPoint+1)*dAlat))
          distortionArray = IdentityMatrix(distortionArray)*(aLat/startAlat)
          Call distortNL(nlBP, distortionArray, cKey)       ! geom.f90
          Call calcE(nlBP, potential, cKey)                 ! staticCalcs.calcE.f90
          points(loops-j,1) = nlBP(cKey)%volumeDistorted/nlBP(cKey)%coordsLength
          points(loops-j,2) = (nlBP(cKey)%totalEnergy/nlBP(cKey)%coordsLength)
        End Do
        coefficients = BirchMurnFit(points)
        optAlat = (coefficients(2)*nlBP(cKey)%coordsLength)**(1.0D0/3.0D0)
! Store data
        bpObj(i)%bccAlat = (coefficients(2)*2)**(1.0D0/3.0D0)
        bpObj(i)%bccV0 = (coefficients(2))
        bpObj(i)%bccE0 = coefficients(1)
        bpObj(i)%bccB0 = coefficients(3)
        bpObj(i)%bccB0_GPA = UnitConvert(bpObj(i)%bccB0, "EVAN3", "GPA")   ! Convert to GPA
        bpObj(i)%bccBp0 = coefficients(4)
! Make EoS points
        Do j=1,15
          V = bpObj(i)%bccV0*(1.0D0-(8-j)*0.01D0)
          E = F_BirchMurn(V,coefficients)  ! scienceFunctions.materialsFunctions.f90
          bpObj(i)%bccEosPoints(j,1) = V
          bpObj(i)%bccEosPoints(j,2) = E
        End Do
      End Do
    End If
!----------------------------------
! Deallocate arrays on HEAP
!----------------------------------
    Deallocate(coordsUnitBP)
    Deallocate(coordsBP)
    Deallocate(nlBP)
  End Subroutine calcEOS

  Subroutine calcEC(bpObj,potential,structure_In)
! Calculate bulk properties
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oBulkProperty), Dimension(:) :: bpObj
    Type(potentialType) :: potential
    Character(*), Optional :: structure_In
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, n, m, cKey
    Character(Len=3) :: structure
    Type(coordsUnitType), Allocatable, Dimension(:) :: coordsUnitBP
    Type(coordsType), Allocatable, Dimension(:) :: coordsBP
    Type(nlType), Allocatable, Dimension(:) :: nlBP
    Character(Len=16), Dimension(1:4) :: atomLabels
    Integer(kind=StandardInteger), Dimension(1:4) :: atomIDs
! EOS
    Integer(kind=StandardInteger) :: loops, loopsMidPoint
    Real(kind=DoubleReal), Dimension(1:p_bpPoints,1:2) :: points
    Real(kind=DoubleReal), Dimension(1:3) :: polyCoefficients
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: distortionArray
! Elastic Constants
    Real(kind=DoubleReal) :: sigma
! Optional Arguments
    structure = "   "
    If(Present(structure_In))Then
      structure = structure_In(1:3)
      structure = StrToUpper(structure)
    End If
!----------------------------------
! Allocate arrays on HEAP
!----------------------------------
    Allocate(coordsUnitBP(1:4))
    Allocate(coordsBP(1:4))
    Allocate(nlBP(1:4))
! Initialise objects
    Call initUnitCoords(coordsUnitBP)     ! geom.f90
! Set loops
    loops = p_bpPoints
    loopsMidPoint = (loops + 1)/2
! Loop through atom IDs
    cKey = 0
    n = 0
    m = 0
!----------------------------------
! FCC Structure
!----------------------------------
    If((structure.eq."   ").or.(structure.eq."FCC"))Then
      Do i=1,potential%atomID_Count
        atomLabels = potential%atomIDs(i)  ! Set atomLabels array and atomIDs array
        atomIDs = i
!----------------------------------
! FCC Build Coords
!----------------------------------
        cKey = cKey + 1
        coordsUnitBP%aLat = bpObj(i)%fccAlat
        coordsUnitBP%xCopy = 4
        coordsUnitBP%yCopy = 4
        coordsUnitBP%zCopy = 4
        Call standardCoords("FCC", coordsUnitBP, cKey, atomLabels, atomIDs)  ! geom.f90
        Call expandUnitCoords(coordsUnitBP, coordsBP)   ! geom.f90
!----------------------------------
! Elastic constants - rebuild neighbour list with new alat
!----------------------------------
        Call makeNL(nlBP, coordsBP, 6.5D0, cKey)
!----------------------------------
! C11, C12
!----------------------------------
        points(1,1) = 0.0D0
        points(1,2) = bpObj(i)%fccE0
        sigma = 0.025D0
        Do j=1,loops-1
          distortionArray = Tensor_Orthorhombic(j*sigma)
          Call distortNL(nlBP, distortionArray, cKey)       ! geom.f90
          Call calcE(nlBP, potential, cKey)
          points(j+1,1) = j*sigma
          points(j+1,2) = (nlBP(cKey)%totalEnergy/nlBP(cKey)%coordsLength)
        End Do
        polyCoefficients = PolyFit(points,2)
        bpObj(i)%fccC11 = (2.0D0*polyCoefficients(3))/(3.0D0*bpObj(i)%fccV0)+bpObj(i)%fccB0
        bpObj(i)%fccC11_GPA = UnitConvert(bpObj(i)%fccC11, "EVAN3", "GPA")   ! Convert to GPA
        bpObj(i)%fccC12 = (3.0D0*bpObj(i)%fccB0-bpObj(i)%fccC11)/2.0D0
        bpObj(i)%fccC12_GPA = UnitConvert(bpObj(i)%fccC12, "EVAN3", "GPA")   ! Convert to GPA
!----------------------------------
! C44
!----------------------------------
        points(1,1) = 0.0D0
        points(1,2) = bpObj(i)%fccE0
        sigma = 0.025D0
        Do j=1,loops-1
          distortionArray = Tensor_Tetragonal(j*sigma)
          Call distortNL(nlBP, distortionArray, cKey)       ! geom.f90
          Call calcE(nlBP, potential, cKey)
          points(j+1,1) = j*sigma
          points(j+1,2) = (nlBP(cKey)%totalEnergy/nlBP(cKey)%coordsLength)
        End Do
        polyCoefficients = PolyFit(points,2)
        bpObj(i)%fccC44 = (2.0D0*polyCoefficients(3))/(bpObj(i)%fccV0)
        bpObj(i)%fccC44_GPA = UnitConvert(bpObj(i)%fccC44, "EVAN3", "GPA")   ! Convert to GPA
      End Do
    End If
!----------------------------------
! BCC Structure
!----------------------------------
    If((structure.eq."   ").or.(structure.eq."BCC"))Then
      Do i=1,potential%atomID_Count
        atomLabels = potential%atomIDs(i)  ! Set atomLabels array and atomIDs array
        atomIDs = i
!----------------------------------
! BCC Build Coords
!----------------------------------
        cKey = cKey + 1
        coordsUnitBP%aLat = bpObj(i)%bccAlat
        coordsUnitBP%xCopy = 4
        coordsUnitBP%yCopy = 4
        coordsUnitBP%zCopy = 4
        Call standardCoords("BCC", coordsUnitBP, cKey, atomLabels, atomIDs)  ! geom.f90
        Call expandUnitCoords(coordsUnitBP, coordsBP)   ! geom.f90
!----------------------------------
! Elastic constants - rebuild neighbour list with new alat
!----------------------------------
        Call makeNL(nlBP, coordsBP, 6.5D0, cKey)
        Call calcE(nlBP, potential, cKey)
!----------------------------------
! C11, C12
!----------------------------------
        points(1,1) = 0.0D0
        points(1,2) = bpObj(i)%bccE0
        sigma = 0.025D0
        Do j=1,loops-1
          distortionArray = Tensor_Orthorhombic(j*sigma)
          Call distortNL(nlBP, distortionArray, cKey)       ! geom.f90
          Call calcE(nlBP, potential, cKey)
          points(j+1,1) = j*sigma
          points(j+1,2) = (nlBP(cKey)%totalEnergy/nlBP(cKey)%coordsLength)
        End Do
        polyCoefficients = PolyFit(points,2)
        bpObj(i)%bccC11 = (2.0D0*polyCoefficients(3))/(3.0D0*bpObj(i)%bccV0)+bpObj(i)%bccB0
        bpObj(i)%bccC11_GPA = UnitConvert(bpObj(i)%bccC11, "EVAN3", "GPA")   ! Convert to GPA
        bpObj(i)%bccC12 = (3.0D0*bpObj(i)%bccB0-bpObj(i)%bccC11)/2.0D0
        bpObj(i)%bccC12_GPA = UnitConvert(bpObj(i)%bccC12, "EVAN3", "GPA")   ! Convert to GPA
!----------------------------------
! C44
!----------------------------------
        points(1,1) = 0.0D0
        points(1,2) = bpObj(i)%bccE0
        sigma = 0.025D0
        Do j=1,loops-1
          distortionArray = Tensor_Tetragonal(j*sigma)
          Call distortNL(nlBP, distortionArray, cKey)       ! geom.f90
          Call calcE(nlBP, potential, cKey)
          points(j+1,1) = j*sigma
          points(j+1,2) = (nlBP(cKey)%totalEnergy/nlBP(cKey)%coordsLength)
        End Do
        polyCoefficients = PolyFit(points,2)
        bpObj(i)%bccC44 = (2.0D0*polyCoefficients(3))/(bpObj(i)%bccV0)
        bpObj(i)%bccC44_GPA = UnitConvert(bpObj(i)%bccC44, "EVAN3", "GPA")   ! Convert to GPA
      End Do
    End If
!----------------------------------
! Deallocate arrays on HEAP
!----------------------------------
    Deallocate(coordsUnitBP)
    Deallocate(coordsBP)
    Deallocate(nlBP)
  End Subroutine calcEC


  Subroutine initBP(bpObj)
! Initialise bulk properties object
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oBulkProperty) :: bpObj
! Vars:  Private
! Initialise
    bpObj%bccAlat = 0.0D0
    bpObj%bccE0 = 0.0D0

    bpObj%fccAlat = 0.0D0
    bpObj%fccE0 = 0.0D0



  End Subroutine initBP
