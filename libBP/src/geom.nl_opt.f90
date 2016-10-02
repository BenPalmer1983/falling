! ------------------------------------------------------------
!               GEOM: Neighbour List
!                     Structure Optimisation
! ------------------------------------------------------------

  Subroutine makeNL_opt(nl, coords, rVerlet, cKey_NL_in)
! Make neighbour list for atoms
! The input coords and alat must be large enough so the same atom does not interact
! with a copy in a periodic cell surrounding the original
! e.g. if the rVerlet cutoff is 5, the alat must be greater than 5
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType_Opt), Dimension(:) :: nl
    Type(coordsType), Dimension(:)  :: coords
    Real(kind=DoubleReal) :: rVerlet
    Integer(kind=StandardInteger), Optional :: cKey_NL_in
! Vars:  Private
    Integer(kind=StandardInteger) :: cKey_NL
    Integer(kind=StandardInteger) :: cKey
    Integer(kind=StandardInteger) :: i, j, k
    Integer(kind=StandardInteger) :: l, m ,n
    Integer(kind=StandardInteger) :: scX, scY, scZ
    Integer(kind=StandardInteger) :: atomA, atomB, atomA_ID, atomB_ID
    Type(oSubCell) :: subCell
    Integer(kind=StandardInteger) :: scKey
    Integer(kind=StandardInteger) :: scKeyA, scKeyB
    Real(kind=DoubleReal) :: aLat
    Real(kind=DoubleReal) :: rVerletSQ
    Integer(kind=StandardInteger) :: nlKey
    Integer(kind=StandardInteger) :: xKey, yKey, zKey
    Integer(kind=StandardInteger) :: xGhostKey, yGhostKey, zGhostKey
    Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB
    Real(kind=DoubleReal) :: xD, yD, zD, rD
    Real(kind=DoubleReal) :: xDsq, yDsq, zDsq, rDsq
    Real(kind=DoubleReal) :: xShift, yShift, zShift
    Real(kind=DoubleReal), Dimension(1:3) :: position, positionT
    Integer(kind=StandardInteger), Dimension(1:3) :: shiftArr
    Integer(kind=StandardInteger) :: inCell
    Integer(kind=StandardInteger) :: ghostCellKey
    !Real(kind=DoubleReal), Dimension(1:27,1:3) :: ghostCellShift
! Vars: Neighbour list estimate
    Integer(kind=StandardInteger) :: nlArrayLength, nlArrayOverKey
!----------------------------------------------
! Optional arguments
!----------------------------------------------
    cKey_NL = 0
    If(Present(cKey_NL_in))Then
      cKey_NL = cKey_NL_in
    End If
!
!----------------------------------------------
! Array Management
!----------------------------------------------
!
! Deallocate if allocated
    Do cKey=1,size(nl,1)
      If((cKey.le.size(nl,1)).and.(coords(cKey)%length.gt.0))Then  ! Check that there is space in the nl object and there are coordinates
        If((cKey_NL.eq.0).or.(cKey_NL.eq.cKey))Then
!-------------------
    If(Allocated(nl(cKey)%label))Then
      Deallocate(nl(cKey)%label)
    End If
    If(Allocated(nl(cKey)%labelID))Then
      Deallocate(nl(cKey)%labelID)
    End If
    If(Allocated(nl(cKey)%coords))Then
      Deallocate(nl(cKey)%coords)
    End If
    If(Allocated(nl(cKey)%coordsMD))Then
      Deallocate(nl(cKey)%coordsMD)
    End If
    If(Allocated(nl(cKey)%coords_Moved))Then
      Deallocate(nl(cKey)%coords_Moved)
    End If
    If(Allocated(nl(cKey)%coordsFixed))Then
      Deallocate(nl(cKey)%coordsFixed)
    End If
    If(Allocated(nl(cKey)%subCellKey))Then
      Deallocate(nl(cKey)%subCellKey)
    End If
    If(Allocated(nl(cKey)%forces))Then
      Deallocate(nl(cKey)%forces)
    End If
    If(Allocated(nl(cKey)%velocity))Then
      Deallocate(nl(cKey)%velocity)
    End If
    If(Allocated(nl(cKey)%velocityH))Then
      Deallocate(nl(cKey)%velocityH)
    End If
    If(Allocated(nl(cKey)%acceleration))Then
      Deallocate(nl(cKey)%acceleration)
    End If
    If(Allocated(nl(cKey)%forcesInitial))Then
      Deallocate(nl(cKey)%forcesInitial)
    End If
    If(Allocated(nl(cKey)%velocityInitial))Then
      Deallocate(nl(cKey)%velocityInitial)
    End If
    If(Allocated(nl(cKey)%charge))Then
      Deallocate(nl(cKey)%charge)
    End If
    If(Allocated(nl(cKey)%mass))Then
      Deallocate(nl(cKey)%mass)
    End If
    If(Allocated(nl(cKey)%nlOverKey))Then
      Deallocate(nl(cKey)%nlOverKey)
    End If
    If(Allocated(nl(cKey)%atomA_ID))Then
      Deallocate(nl(cKey)%atomA_ID)
    End If
    If(Allocated(nl(cKey)%atomB_ID))Then
      Deallocate(nl(cKey)%atomB_ID)
    End If
    If(Allocated(nl(cKey)%atomA_Type))Then
      Deallocate(nl(cKey)%atomA_Type)
    End If
    If(Allocated(nl(cKey)%atomB_Type))Then
      Deallocate(nl(cKey)%atomB_Type)
    End If
    If(Allocated(nl(cKey)%atomPairKey))Then
      Deallocate(nl(cKey)%atomPairKey)
    End If
    If(Allocated(nl(cKey)%inCell))Then
      Deallocate(nl(cKey)%inCell)
    End If
    If(Allocated(nl(cKey)%rD))Then
      Deallocate(nl(cKey)%rD)
    End If
    If(Allocated(nl(cKey)%rD_original))Then
      Deallocate(nl(cKey)%rD_original)
    End If
    If(Allocated(nl(cKey)%vecAB))Then
      Deallocate(nl(cKey)%vecAB)
    End If
    If(Allocated(nl(cKey)%ghostCell))Then
      Deallocate(nl(cKey)%ghostCell)
    End If
    If(Allocated(nl(cKey)%cell))Then
      Deallocate(nl(cKey)%cell)
    End If
    If(Allocated(nl(cKey)%electronDensityFixed))Then
      Deallocate(nl(cKey)%electronDensityFixed)
    End If
    If(Allocated(nl(cKey)%rD_moved))Then
      Deallocate(nl(cKey)%rD_moved)
    End If
    If(Allocated(nl(cKey)%atomEnergyFixed))Then
      Deallocate(nl(cKey)%atomEnergyFixed)
    End If
!-------------------
        End If
      End If
    End Do
! Allocate if not allocated
    Do cKey=1,size(nl,1)
      If((cKey.le.size(nl,1)).and.(coords(cKey)%length.gt.0))Then  ! Check that there is space in the nl object and there are coordinates
        If((cKey_NL.eq.0).or.(cKey_NL.eq.cKey))Then
!-------------------
! Estimate NL size
      nlArrayLength = 10000+2000+&
        ceiling(0.35D0*(coords(cKey)%length)**2*(8.0D0*rVerlet**3)/((coords(cKey)%aLat)**3))
      nl(cKey)%arrayLength = nlArrayLength
      nlArrayOverKey = ceiling((3.0D0*nlArrayLength)/coords(cKey)%length)
! Coord details
      Allocate(nl(cKey)%label(1:coords(cKey)%length))
      Allocate(nl(cKey)%labelID(1:coords(cKey)%length))
      Allocate(nl(cKey)%coords(1:coords(cKey)%length,1:3))
      Allocate(nl(cKey)%coordsMD(1:coords(cKey)%length,1:3))
      Allocate(nl(cKey)%coords_Moved(1:coords(cKey)%length,1:3))
      Allocate(nl(cKey)%coordsFixed(1:coords(cKey)%length))
      Allocate(nl(cKey)%subCellKey(1:coords(cKey)%length))
      Allocate(nl(cKey)%forces(1:coords(cKey)%length,1:3))
      Allocate(nl(cKey)%velocity(1:coords(cKey)%length,1:3))
      Allocate(nl(cKey)%velocityH(1:coords(cKey)%length,1:3))
      Allocate(nl(cKey)%acceleration(1:coords(cKey)%length,1:3))
      Allocate(nl(cKey)%forcesInitial(1:coords(cKey)%length,1:3))
      Allocate(nl(cKey)%velocityInitial(1:coords(cKey)%length,1:3))
      Allocate(nl(cKey)%charge(1:coords(cKey)%length))
      Allocate(nl(cKey)%mass(1:coords(cKey)%length))
! Atom details
      Allocate(nl(cKey)%nlOverKey(1:coords(cKey)%length,0:nlArrayOverKey))
      Allocate(nl(cKey)%atomA_ID(1:nlArrayLength))
      Allocate(nl(cKey)%atomB_ID(1:nlArrayLength))
      Allocate(nl(cKey)%atomA_Type(1:nlArrayLength))
      Allocate(nl(cKey)%atomB_Type(1:nlArrayLength))
      Allocate(nl(cKey)%atomPairKey(1:nlArrayLength))
      Allocate(nl(cKey)%inCell(1:nlArrayLength))
! Position Details
      Allocate(nl(cKey)%rD(1:nlArrayLength))
      Allocate(nl(cKey)%rD_original(1:nlArrayLength))
      Allocate(nl(cKey)%vecAB(1:nlArrayLength,1:3))
! Subcell
      Allocate(nl(cKey)%ghostCell(1:nlArrayLength,1:4))
      Allocate(nl(cKey)%cell(1:nlArrayLength,1:3))
! Vars EFS calculation
      Allocate(nl(cKey)%electronDensity(1:coords(cKey)%length,1:3))
      Allocate(nl(cKey)%atomEnergy(1:coords(cKey)%length,1:3))
! Vars used for opt calcs
      Allocate(nl(cKey)%electronDensityFixed(1:coords(cKey)%length,1:3))
      Allocate(nl(cKey)%rD_moved(1:coords(cKey)%length))
      Allocate(nl(cKey)%atomEnergyFixed(1:coords(cKey)%length))
!-------------------
        End If
      End If
    End Do
!
!----------------------------------------------
! Build NL
!----------------------------------------------
!
! Loop through configs
!
    Do cKey=1,size(coords,1)
      If((cKey.le.size(nl,1)).and.(coords(cKey)%length.gt.0))Then  ! Check that there is space in the nl object and there are coordinates
        If((cKey_NL.eq.0).or.(cKey_NL.eq.cKey))Then
!----------------------------------------------
! Store from coords
    nl(cKey)%coordsLength = coords(cKey)%length
    nl(cKey)%aLat = coords(cKey)%aLat
    nl(cKey)%aLat_Original = nl(cKey)%aLat
    nl(cKey)%atomID_Count = coords(cKey)%atomID_Count  ! number of unique atom ids
    nl(cKey)%atomIDs = coords(cKey)%atomIDs            ! list of unique atom ids
! Unit cell from coords
    nl(cKey)%unitCell = coords(cKey)%unitCell
    nl(cKey)%unitCellScaled = nl(cKey)%aLat*nl(cKey)%unitCell
    nl(cKey)%unitCellScaledInverse = InvertMatrix(nl(cKey)%unitCellScaled)
    nl(cKey)%volume = TripleProductSq(nl(cKey)%unitCellScaled)
! Load coords
    Do i=1,nl(cKey)%coordsLength
      nl(cKey)%label(i) = coords(cKey)%label(i)
      nl(cKey)%labelID(i) = coords(cKey)%labelID(i)
      nl(cKey)%charge(i) = coords(cKey)%charge(i)
      nl(cKey)%mass(i) = coords(cKey)%mass(i)
      nl(cKey)%coordsFixed(i) = coords(cKey)%coordsFixed(i)
      Do j=1,3
        position(j) = Modulus(coords(cKey)%coords(i,j),1.0D0)    ! x/y/z coord, keep within a 1x1x1 box
      End Do
      positionT = MatMul(nl(cKey)%unitCellScaled,position)  ! Transform
      Do j=1,3
        nl(cKey)%coords(i,j) = position(j)
        nl(cKey)%coordsMD(i,j) = positionT(j)
        nl(cKey)%forces(i,j) = coords(cKey)%forces(i,j)      ! Load initial forces as the starting MD forces (will be recalculated at each MD step anyway)
        nl(cKey)%velocity(i,j) = coords(cKey)%velocity(i,j)  ! Load initial velocity
        nl(cKey)%acceleration(i,j) = 0.0D0
        nl(cKey)%forcesInitial(i,j) = coords(cKey)%forces(i,j)
        nl(cKey)%velocityInitial(i,j) = coords(cKey)%velocity(i,j)
      End Do
    End Do
!----------------------------------------------
! Init config specific variables
    rVerletSQ = rVerlet**2
    aLat = coords(cKey)%aLat
    nl(cKey)%totalRD = 0.0D0
    nl(cKey)%totalRDSq = 0.0D0
    nl(cKey)%rVerlet = rVerlet
    nl(cKey)%rVerletSQ = rVerletSQ
    nl(cKey)%rMin = rVerlet
    nl(cKey)%rMax = 0.0D0
    nlKey = 0
!----------------------------------------------
! Init overkey count
    Do i = 1,nl(cKey)%coordsLength
      nl(cKey)%nlOverKey(i,0) = 0
    End Do
!----------------------------------------------
! calculate sub cell parameters
    subCell%width = floor(aLat/(1.0D0*rVerlet))   ! Number of subcells (1D)
    If(subCell%width.lt.1)Then
      subCell%width = 1
    End If
    subCell%count = subCell%width**3                                        ! Total number of subcells in 3D
    subCell%fracAlat = 1.0D0/(1.0D0*subCell%width)
    subCell%aLat = aLat/(1.0D0*subCell%width)                               ! Subcell alat size
    subCell%maxAtoms = 5*ceiling(nl(cKey)%coordsLength/(1.0D0*subCell%count))  ! Estimate max atoms per subcell
    subCell%cellCombinations = DoubleKey(subCell%count,subCell%count)
!----------------------------------------------
! Allocate arrays
!----------------------------------------------
    Allocate(subCell%atomCount(1:subCell%count))
    Allocate(subCell%keyArr(1:subCell%count,1:3))
    Allocate(subCell%atomIDs(1:subCell%count,1:subCell%maxAtoms,1:2))
    Allocate(subCell%atomCoords(1:subCell%count,1:subCell%maxAtoms,1:3))
! Store in NL
    nl(cKey)%scCount = subCell%count
    nl(cKey)%scAlat = subCell%aLat
    nl(cKey)%maxAtomsPerSC = subCell%maxAtoms
    nl(cKey)%subCellWidth = subCell%width
! Zero atom count array
    subCell%atomCount = 0
! Fill in sub cell keys
    Do k=1,subCell%width
      Do j=1,subCell%width
        Do i=1,subCell%width
          scKey = i+subCell%width*(j-1)+subCell%width**2*(k-1)
          subCell%keyArr(scKey,1) = i
          subCell%keyArr(scKey,2) = j
          subCell%keyArr(scKey,3) = k
        End Do
      End Do
    End Do
! Load atom ids and coords into arrays
! Loop through atoms
    Do i=1,nl(cKey)%coordsLength
! Make sub cell key
      scKey = SubCellKey(subCell%fracAlat,subCell%width,nl(cKey)%coords(i,1),nl(cKey)%coords(i,2),nl(cKey)%coords(i,3))
      nl(cKey)%subCellKey(i) = scKey
      subCell%atomCount(scKey) = subCell%atomCount(scKey) + 1
! Store atom IDs + coords in sub cell arrays
      subCell%atomIDs(scKey,subCell%atomCount(scKey),1) = i                          ! unique atom id
      subCell%atomIDs(scKey,subCell%atomCount(scKey),2) = nl(cKey)%labelID(i)        ! atom type
      Do j=1,3
        subCell%atomCoords(scKey,subCell%atomCount(scKey),j) = nl(cKey)%coords(i,j)     ! atom coords
      End Do
    End Do
! Transform 1x1x1 coords
    Do i=1,subCell%count
      Do j=1,subCell%atomCount(i)
        Do k=1,3
          position(k) = subCell%atomCoords(i,j,k)
        End Do
        positionT = MatMul(nl(cKey)%unitCellScaled,position)
        Do k=1,3
          subCell%atomCoords(i,j,k) = positionT(k)
        End Do
      End Do
    End Do
! Calculate ghostCellShift
    Do l=-1,1
      Do m=-1,1
        Do n=-1,1
          ghostCellKey = (l+2)+3*(m+1)+9*(n+1)
          position(1) = l
          position(2) = m
          position(3) = n
          positionT = MatMul(nl(cKey)%unitCellScaled,position)
          subCell%ghostCellShift(ghostCellKey,1) = positionT(1)
          subCell%ghostCellShift(ghostCellKey,2) = positionT(2)
          subCell%ghostCellShift(ghostCellKey,3) = positionT(3)
        End Do
      End Do
    End Do
! Loop through subcells
    Do scKeyA=1,subCell%count
! loop through surrounding and image sub cells
      Do l=-1,1
        Do m=-1,1
          Do n=-1,1
! Transform (l,m,n) array
            position(1) = l
            position(2) = m
            position(3) = m
! Atom B subcell key
            xKey = subCell%keyArr(scKeyA,1)+l
            yKey = subCell%keyArr(scKeyA,2)+m
            zKey = subCell%keyArr(scKeyA,3)+n
! PBC subcell key
            scX = Modulus(xKey-1,subCell%width)+1
            scY = Modulus(yKey-1,subCell%width)+1
            scZ = Modulus(zKey-1,subCell%width)+1
            scKeyB = scX+subCell%width*(scY-1)+subCell%width**2*(scZ-1)
! GhostKey settings
            shiftArr = 0
            inCell = 1
            xGhostKey = 0
            yGhostKey = 0
            zGhostKey = 0
            If(xKey.lt.1)Then
              xGhostKey = -1
              shiftArr(1) = -1
              inCell = 0
            End If
            If(xKey.gt.subCell%width)Then
              xGhostKey = 1
              shiftArr(1) = 1
              inCell = 0
            End If
            If(yKey.lt.1)Then
              yGhostKey = -1
              shiftArr(2) = -1
              inCell = 0
            End If
            If(yKey.gt.subCell%width)Then
              yGhostKey = 1
              shiftArr(2) = 1
              inCell = 0
            End If
            If(zKey.lt.1)Then
              zGhostKey = -1
              shiftArr(3) = -1
              inCell = 0
            End If
            If(zKey.gt.subCell%width)Then
              zGhostKey = 1
              shiftArr(3) = 1
              inCell = 0
            End If
            ghostCellKey = (xGhostKey+2)+(subCell%width+1)*(yGhostKey+1)+(subCell%width+1)**2*(zGhostKey+1)
            xShift = subCell%ghostCellShift(ghostCellKey,1)
            yShift = subCell%ghostCellShift(ghostCellKey,2)
            zShift = subCell%ghostCellShift(ghostCellKey,3)
! loop through A atoms
            Do atomA = 1,subCell%atomCount(scKeyA)
! loop through B atoms
              Do atomB = 1,subCell%atomCount(scKeyB)
                atomA_ID = subCell%atomIDs(scKeyA,atomA,1)
                atomB_ID = subCell%atomIDs(scKeyB,atomB,1)
! Only count if A is lower than B (half length neighbour list)
                If(atomA_ID.lt.atomB_ID)Then
                  xA = subCell%atomCoords(scKeyA,atomA,1)
                  xB = subCell%atomCoords(scKeyB,atomB,1)+xShift
                  xD = xA-xB
                  xDsq = xd**2
                  If(xDsq.le.rVerletSQ)Then
                    yA = subCell%atomCoords(scKeyA,atomA,2)
                    yB = subCell%atomCoords(scKeyB,atomB,2)+yShift
                    yD = yA-yB
                    yDsq = yd**2
                    If(yDsq.le.rVerletSQ)Then
                      zA = subCell%atomCoords(scKeyA,atomA,3)
                      zB = subCell%atomCoords(scKeyB,atomB,3)+zShift
                      zD = zA-zB
                      zDsq = zd**2
                      If(zDsq.le.rVerletSQ)Then
                        rdSq = xDsq + yDsq + zDsq
                        If(rdSq.le.rVerletSq)Then
                          rD = sqrt(rdSq)
                          ! Key
                          nlKey = nlKey + 1
                          ! Overkey
                          nl(cKey)%nlOverKey(atomA_ID,0) = nl(cKey)%nlOverKey(atomA_ID,0) + 1
                          nl(cKey)%nlOverKey(atomA_ID,nl(cKey)%nlOverKey(atomA_ID,0)) = nlKey
                          nl(cKey)%nlOverKey(atomB_ID,0) = nl(cKey)%nlOverKey(atomB_ID,0) + 1
                          nl(cKey)%nlOverKey(atomB_ID,nl(cKey)%nlOverKey(atomB_ID,0)) = nlKey
                          ! Key/Type
                          nl(cKey)%atomA_ID(nlKey) = atomA_ID                           ! A ID
                          nl(cKey)%atomB_ID(nlKey) = atomB_ID                           ! B ID
                          nl(cKey)%atomA_Type(nlKey) = subCell%atomIDs(scKeyA,atomA,2)  ! A type
                          nl(cKey)%atomB_Type(nlKey) = subCell%atomIDs(scKeyB,atomB,2)  ! B type
                          nl(cKey)%inCell(nlKey) = inCell                               ! 1 if in cell, 0 if out of cell
                          nl(cKey)%atomPairKey = 0
                          ! Displacement/Direction
                          nl(cKey)%rD(nlKey) = sqrt(rdSq)
                          nl(cKey)%rD_original(nlKey) = nl(cKey)%rD(nlKey)
                          nl(cKey)%vecAB(nlKey,1) = xD/nl(cKey)%rD(nlKey)               ! Vector from B to A (x)
                          nl(cKey)%vecAB(nlKey,2) = yD/nl(cKey)%rD(nlKey)               ! Vector from B to A (y)
                          nl(cKey)%vecAB(nlKey,3) = zD/nl(cKey)%rD(nlKey)               ! Vector from B to A (z)
                          ! super cell coords of atom B
                          nl(cKey)%ghostCell(nlKey,1) = shiftArr(1)
                          nl(cKey)%ghostCell(nlKey,2) = shiftArr(2)
                          nl(cKey)%ghostCell(nlKey,3) = shiftArr(3)
                          nl(cKey)%ghostCell(nlKey,4) = ghostCellKey
                          If(rD.gt.nl(cKey)%rMax)Then
                            nl(cKey)%rMax = rD
                          End If
                          If(rD.lt.nl(cKey)%rMin)Then
                            nl(cKey)%rMin = rD
                          End If
                        End If
                      End If
                    End If
                  End If
                End If
              End Do ! End loop through B atoms
            End Do ! End loop through A atoms
          End Do
        End Do
      End Do
    End Do ! End scKeyA loop
    nl(cKey)%length = nlKey
!----------------------------------------------
! Dellocate arrays
!----------------------------------------------
    Deallocate(subCell%atomCount)
    Deallocate(subCell%keyArr)
    Deallocate(subCell%atomIDs)
    Deallocate(subCell%atomCoords)
!----------------------------------------------
        End If
      End If
    End Do
!----------------------------------------------
  End Subroutine makeNL_Opt



  Subroutine nlMoveAtom_Opt(nlGeom, cKey, atomKey, coordChange)
! Moves 1 atom, stores moved coords and moved rD in coords_Moved and rD_moved
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType_Opt), Dimension(:) :: nlGeom
    Integer(kind=StandardInteger) :: cKey
    Integer(kind=StandardInteger) :: atomKey, atomKeyBs
    Real(kind=DoubleReal), Dimension(1:3) :: coordChange
    Real(kind=DoubleReal), Dimension(1:3) :: aVec, bVec, aPVec
! Vars:  Private
    Integer(kind=StandardInteger) :: i, nKey
! Init vars
    aVec(1) = nlGeom(cKey)%coords(atomKey,1)
    aVec(2) = nlGeom(cKey)%coords(atomKey,2)
    aVec(3) = nlGeom(cKey)%coords(atomKey,3)
    aPVec(1) = nlGeom(cKey)%coords(atomKey,1)+coordChange(1)
    aPVec(2) = nlGeom(cKey)%coords(atomKey,2)+coordChange(2)
    aPVec(3) = nlGeom(cKey)%coords(atomKey,3)+coordChange(3)
! Store moved atom coord
    nlGeom(cKey)%coords_Moved(atomKey,1) = aPVec(1)
    nlGeom(cKey)%coords_Moved(atomKey,2) = aPVec(2)
    nlGeom(cKey)%coords_Moved(atomKey,3) = aPVec(3)
! Loop through nl
    Do i=1,nlGeom(cKey)%nlOverKey(atomKey,0)
      nKey = nlGeom(cKey)%nlOverKey(atomKey,i)
      atomKeyBs = nlGeom(cKey)%atomB_ID(nKey)
      !bVec(1) = nlGeom(cKey)%coords(atomKeyBs,1)+nlGeom(cKey)%subCell(nKey,1)
      !bVec(2) = nlGeom(cKey)%coords(atomKeyBs,2)+nlGeom(cKey)%subCell(nKey,2)
      !bVec(3) = nlGeom(cKey)%coords(atomKeyBs,3)+nlGeom(cKey)%subCell(nKey,3)
      nlGeom(cKey)%rD_moved(nKey) = nlGeom(cKey)%aLat*sqrt(&
      (aPVec(1)-bVec(1))**2+(aPVec(2)-bVec(2))**2+(aPVec(3)-bVec(3))**2)
    End Do
  End Subroutine nlMoveAtom_Opt

  Subroutine nlResetAtom_Opt(nlGeom, cKey, atomKey)
! Moves atom, updates RD only
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType_Opt), Dimension(:) :: nlGeom
    Integer(kind=StandardInteger) :: cKey
    Integer(kind=StandardInteger) :: atomKey
! Vars:  Private
    Integer(kind=StandardInteger) :: i, nKey
! Restore original coords
    nlGeom(cKey)%coords_Moved(atomKey,1) = nlGeom(cKey)%coords(atomKey,1)
    nlGeom(cKey)%coords_Moved(atomKey,2) = nlGeom(cKey)%coords(atomKey,2)
    nlGeom(cKey)%coords_Moved(atomKey,3) = nlGeom(cKey)%coords(atomKey,3)
! Loop through nl and restore original coords
    Do i=1,nlGeom(cKey)%nlOverKey(atomKey,0)
      nKey = 1
      nlGeom(cKey)%rD_moved(nKey) = nlGeom(cKey)%rD(nKey)
    End Do
  End Subroutine nlResetAtom_Opt




!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------




  Function SubCellKey_Opt (scAlat, scW, x, y, z) Result (key)
! Return sub cell key
    Implicit None   ! Force declaration of all variables
! In
    Real(kind=DoubleReal) :: scAlat, x, y, z
    Integer(kind=StandardInteger) :: scW
! Out
    Integer(kind=StandardInteger) :: key
! Private
    Integer(kind=StandardInteger) :: i, j, k
    i = floor(x/(1.0D0*scAlat))+1
    j = floor(y/(1.0D0*scAlat))+1
    k = floor(z/(1.0D0*scAlat))+1
    key = i+scW*(j-1)+scW**2*(k-1)
  End Function SubCellKey_Opt
