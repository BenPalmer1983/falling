! ------------------------------------------------------------
!               GEOM: Neighbour List
! ------------------------------------------------------------


  Subroutine makeNL(nl, coords, rVerlet, cKey_NL_in)
! Make neighbour list for atoms
! The input coords and alat must be large enough so the same atom does not interact
! with a copy in a periodic cell surrounding the original
! e.g. if the rVerlet cutoff is 5, the alat must be greater than 5
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType), Dimension(:) :: nl
    Type(coordsType), Dimension(:)  :: coords
    Real(kind=DoubleReal) :: rVerlet
    Integer(kind=StandardInteger), Optional :: cKey_NL_in
! Vars:  Private
    Integer(kind=StandardInteger) :: cKey_NL
    Integer(kind=StandardInteger) :: cKey
    Integer(kind=StandardInteger) :: i, j, k
    Integer(kind=StandardInteger) :: l, m ,n
    Integer(kind=StandardInteger) :: scX, scY, scZ
    Integer(kind=StandardInteger) :: atomA, atomB, atomA_ID, atomB_ID, loopStartAtom
    Type(oSubCell) :: subCell
    Integer(kind=StandardInteger) :: coordLength, scW, scCount, scKey, maxAtomsPerSC
    Integer(kind=StandardInteger) :: scKeyA, scKeyB, scKeyA_Loop, scKeyB_Loop, scKeyB_P, scKeyPair
    Real(kind=DoubleReal) :: aLat
    Real(kind=DoubleReal) :: rVerletSQ, scAlat
    Integer(kind=StandardInteger) :: nlKey, uKey
    Integer(kind=StandardInteger) :: xKey, yKey, zKey
    Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB
    Real(kind=DoubleReal) :: xD, yD, zD, rD
    Real(kind=DoubleReal) :: xDsq, yDsq, zDsq, rDsq
    Real(kind=DoubleReal) :: xShift, yShift, zShift
    Real(kind=DoubleReal), Dimension(1:3) :: position, positionT
    Integer(kind=StandardInteger), Dimension(1:3) :: shiftArr
    Integer(kind=StandardInteger) :: inCell
    Logical :: ghostCell, processCellPair
! Vars: Neighbour list estimate
    Integer(kind=StandardInteger) :: nlArrayLength, nlArrayOverKey
! Vars:  Allocatable arrays
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: scAtomCount     ! Number of atoms per sub cell
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: scKeyArr        ! array of scKey and x,y,z position
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:,:) :: scCoordsI
    Real(kind=DoubleReal), Allocatable, Dimension(:,:,:) :: scCoordsR
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: scKeyMap
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: uKeyArr
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
    If(Allocated(nl(cKey)%subCellKey))Then
      Deallocate(nl(cKey)%subCellKey)
    End If
    If(Allocated(nl(cKey)%forces))Then
      Deallocate(nl(cKey)%forces)
    End If
    If(Allocated(nl(cKey)%velocity))Then
      Deallocate(nl(cKey)%velocity)
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
    If(Allocated(nl(cKey)%vecAB))Then
      Deallocate(nl(cKey)%vecAB)
    End If
    If(Allocated(nl(cKey)%subCell))Then
      Deallocate(nl(cKey)%subCell)
    End If
    If(Allocated(nl(cKey)%cell))Then
      Deallocate(nl(cKey)%cell)
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
      Allocate(nl(cKey)%subCellKey(1:coords(cKey)%length))
      Allocate(nl(cKey)%forces(1:coords(cKey)%length,1:3))
      Allocate(nl(cKey)%velocity(1:coords(cKey)%length,1:3))
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
      Allocate(nl(cKey)%vecAB(1:nlArrayLength,1:3))
! Subcell
      Allocate(nl(cKey)%subCell(1:nlArrayLength,1:3))
      Allocate(nl(cKey)%cell(1:nlArrayLength,1:3))
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
! Load coords
    Do i=1,nl(cKey)%coordsLength
      nl(cKey)%label(i) = coords(cKey)%label(i)
      nl(cKey)%labelID(i) = coords(cKey)%labelID(i)
      nl(cKey)%charge(i) = coords(cKey)%charge(i)
      nl(cKey)%mass(i) = coords(cKey)%mass(i)
      Do j=1,3
        position(j) = coords(cKey)%coords(i,j)   ! x/y/z coord
        position(j) = Modulus(position(j),1.0D0) ! Keep within a 1x1x1 box
      End Do
      positionT = MatMul(nl(cKey)%unitCellScaled,position)  ! Transform
      Do j=1,3
        nl(cKey)%coords(i,j) = positionT(j)
        nl(cKey)%forces(i,j) = 0.0D0
        nl(cKey)%velocity(i,j) = 0.0D0
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
    subCell%aLat = aLat/(1.0D0*subCell%width)                               ! Subcell width
    subCell%maxAtoms = 5*ceiling(nl(cKey)%coordsLength/(1.0D0*subCell%count))  ! Estimate max atoms per subcell
    subCell%cellCombinations = DoubleKey(subCell%count,subCell%count)
    print *,subCell%cellCombinations
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
    nl(cKey)%maxAtomsPerSC = maxAtomsPerSC
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
      scKey = SubCellKey(subCell%aLat,subCell%width,nl(cKey)%coords(i,1),nl(cKey)%coords(i,2),nl(cKey)%coords(i,3))
      nl(cKey)%subCellKey(i) = scKey
      subCell%atomCount(scKey) = subCell%atomCount(scKey) + 1
! Store atom IDs + coords in sub cell arrays
      subCell%atomIDs(scKey,subCell%atomCount(scKey),1) = i                          ! unique atom id
      subCell%atomIDs(scKey,subCell%atomCount(scKey),2) = nl(cKey)%labelID(i)        ! atom type
      Do j=1,3
        subCell%atomCoords(scKey,subCell%atomCount(scKey),j) = nl(cKey)%coords(i,j)     ! atom coords
      End Do
    End Do
! Loop through subcells
    Do scKeyA=1,subCell%count
! loop through surrounding and image sub cells
      Do l=-1,1
        Do m=-1,1
          Do n=-1,1
            ghostCell = .false.
! Atom B subcell key
            xKey = subCell%keyArr(scKeyA,1)+l
            yKey = subCell%keyArr(scKeyA,2)+m
            zKey = subCell%keyArr(scKeyA,3)+n
! PBC subcell key
            scX = Modulus(xKey-1,subCell%width)+1
            scY = Modulus(yKey-1,subCell%width)+1
            scZ = Modulus(zKey-1,subCell%width)+1
            scKeyB = scX+subCell%width*(scY-1)+subCell%width**2*(scZ-1)
! set default coord shift
            xShift = 0.0D0
            yShift = 0.0D0
            zShift = 0.0D0
            shiftArr = 0
! adjust shift
            inCell = 1
            If(xKey.lt.1)Then
              xShift = -1.0D0 * aLat
              shiftArr(1) = -1
              inCell = 0
            End If
            If(yKey.lt.1)Then
              yShift = -1.0D0 * aLat
              shiftArr(2) = -1
              inCell = 0
            End If
            If(zKey.lt.1)Then
              zShift = -1.0D0 * aLat
              shiftArr(3) = -1
              inCell = 0
            End If
            If(xKey.gt.subCell%width)Then
              xShift = 1.0D0 * aLat
              shiftArr(1) = 1
              inCell = 0
            End If
            If(yKey.gt.subCell%width)Then
              yShift = 1.0D0 * aLat
              shiftArr(2) = 1
              inCell = 0
            End If
            If(zKey.gt.subCell%width)Then
              zShift = 1.0D0 * aLat
              shiftArr(3) = 1
              inCell = 0
            End If
! loop through A atoms
            Do atomA = 1,subCell%atomCount(scKeyA)
! loop through B atoms
              Do atomB = 1,subCell%atomCount(scKeyB)

                atomA_ID = subCell%atomIDs(scKeyA,atomA,1)
                atomB_ID = subCell%atomIDs(scKeyB,atomB,1)

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
                        nl(cKey)%vecAB(nlKey,1) = xD/nl(cKey)%rD(nlKey)               ! Vector from B to A (x)
                        nl(cKey)%vecAB(nlKey,2) = yD/nl(cKey)%rD(nlKey)               ! Vector from B to A (y)
                        nl(cKey)%vecAB(nlKey,3) = zD/nl(cKey)%rD(nlKey)               ! Vector from B to A (z)
                        ! super cell coords of atom B
                        nl(cKey)%subCell(nlKey,1) = shiftArr(1)
                        nl(cKey)%subCell(nlKey,2) = shiftArr(2)
                        nl(cKey)%subCell(nlKey,3) = shiftArr(3)
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
    End Do
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
  End Subroutine makeNL







  Subroutine makeNL_Broke(nl, coords, rVerlet, cKey_NL_in)
! Make neighbour list for atoms
! The input coords and alat must be large enough so the same atom does not interact
! with a copy in a periodic cell surrounding the original
! e.g. if the rVerlet cutoff is 5, the alat must be greater than 5
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType), Dimension(:) :: nl
    Type(coordsType), Dimension(:)  :: coords
    Real(kind=DoubleReal) :: rVerlet
    Integer(kind=StandardInteger), Optional :: cKey_NL_in
! Vars:  Private
    Integer(kind=StandardInteger) :: cKey_NL
    Integer(kind=StandardInteger) :: cKey
    Integer(kind=StandardInteger) :: i, j, k
    Integer(kind=StandardInteger) :: l, m ,n
    Integer(kind=StandardInteger) :: scX, scY, scZ
    Integer(kind=StandardInteger) :: atomA, atomB, atomA_ID, atomB_ID, loopStartAtom
    Integer(kind=StandardInteger) :: coordLength, scW, scCount, scKey, maxAtomsPerSC
    Integer(kind=StandardInteger) :: scKeyA, scKeyB, scKeyA_Loop, scKeyB_Loop, scKeyB_P
    Real(kind=DoubleReal) :: aLat
    Real(kind=DoubleReal) :: rVerletSQ, scAlat
    Integer(kind=StandardInteger) :: nlKey, uKey
    Integer(kind=StandardInteger) :: xKey, yKey, zKey
    Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB
    Real(kind=DoubleReal) :: xD, yD, zD, rD
    Real(kind=DoubleReal) :: xDsq, yDsq, zDsq, rDsq
    Real(kind=DoubleReal) :: xShift, yShift, zShift
    Real(kind=DoubleReal), Dimension(1:3) :: position, positionT
    Integer(kind=StandardInteger), Dimension(1:3) :: shiftArr
    Integer(kind=StandardInteger) :: inCell
! Vars: Neighbour list estimate
    Integer(kind=StandardInteger) :: nlArrayLength, nlArrayOverKey
! Vars:  Allocatable arrays
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: scAtomCount     ! Number of atoms per sub cell
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: scKeyArr        ! array of scKey and x,y,z position
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:,:) :: scCoordsI
    Real(kind=DoubleReal), Allocatable, Dimension(:,:,:) :: scCoordsR
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: scKeyMap
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: uKeyArr
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
    If(Allocated(nl(cKey)%subCellKey))Then
      Deallocate(nl(cKey)%subCellKey)
    End If
    If(Allocated(nl(cKey)%forces))Then
      Deallocate(nl(cKey)%forces)
    End If
    If(Allocated(nl(cKey)%velocity))Then
      Deallocate(nl(cKey)%velocity)
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
    If(Allocated(nl(cKey)%vecAB))Then
      Deallocate(nl(cKey)%vecAB)
    End If
    If(Allocated(nl(cKey)%subCell))Then
      Deallocate(nl(cKey)%subCell)
    End If
    If(Allocated(nl(cKey)%cell))Then
      Deallocate(nl(cKey)%cell)
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
      nlArrayLength = 2000+&
        ceiling(0.35D0*(coords(cKey)%length)**2*(8.0D0*rVerlet**3)/((coords(cKey)%aLat)**3))
      nl(cKey)%arrayLength = nlArrayLength
!-------------------
! Estimate NL size
      nlArrayLength = 2000+&
        ceiling(0.35D0*(coords(cKey)%length)**2*(8.0D0*rVerlet**3)/((coords(cKey)%aLat)**3))
      nl(cKey)%arrayLength = nlArrayLength
      nlArrayOverKey = ceiling((3.0D0*nlArrayLength)/coords(cKey)%length)
! Coord details
      Allocate(nl(cKey)%label(1:coords(cKey)%length))
      Allocate(nl(cKey)%labelID(1:coords(cKey)%length))
      Allocate(nl(cKey)%coords(1:coords(cKey)%length,1:3))
      Allocate(nl(cKey)%subCellKey(1:coords(cKey)%length))
      Allocate(nl(cKey)%forces(1:coords(cKey)%length,1:3))
      Allocate(nl(cKey)%velocity(1:coords(cKey)%length,1:3))
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
      Allocate(nl(cKey)%vecAB(1:nlArrayLength,1:3))
! Subcell
      Allocate(nl(cKey)%subCell(1:nlArrayLength,1:3))
      Allocate(nl(cKey)%cell(1:nlArrayLength,1:3))
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
! Load coords
    Do i=1,nl(cKey)%coordsLength
      nl(cKey)%label(i) = coords(cKey)%label(i)
      nl(cKey)%labelID(i) = coords(cKey)%labelID(i)
      nl(cKey)%charge(i) = coords(cKey)%charge(i)
      nl(cKey)%mass(i) = coords(cKey)%mass(i)
      Do j=1,3
        position(j) = coords(cKey)%coords(i,j)   ! x/y/z coord
        position(j) = Modulus(position(j),1.0D0) ! Keep within a 1x1x1 box
      End Do
      positionT = MatMul(nl(cKey)%unitCellScaled,position)  ! Transform
      Do j=1,3
        nl(cKey)%coords(i,j) = positionT(j)
        nl(cKey)%forces(i,j) = 0.0D0
        nl(cKey)%velocity(i,j) = 0.0D0
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
    scW = floor(aLat/(1.0D0*rVerlet))   ! Number of subcells (1D)
    If(scW.lt.1)Then
      scW = 1
    End If
    scCount = scW**3                                        ! Total number of subcells in 3D
    scAlat = aLat/(1.0D0*scW)                               ! Subcell width
    maxAtomsPerSC = 5*ceiling(nl(cKey)%coordsLength/(1.0D0*scCount))  ! Estimate max atoms per subcell
    nl(cKey)%scCount = scCount
    nl(cKey)%scAlat = scAlat
    nl(cKey)%maxAtomsPerSC = maxAtomsPerSC
! Allocate arrays
    Allocate(scAtomCount(1:scCount))
    Allocate(scKeyArr(1:scCount,1:3))
    Allocate(scCoordsI(1:scCount,1:maxAtomsPerSC,1:2))
    Allocate(scCoordsR(1:scCount,1:maxAtomsPerSC,1:3))
    Allocate(scKeyMap(1:scCount,1:2))
    Allocate(uKeyArr(1:10000))
    uKeyArr = 0
! Init subcell arrays
    scAtomCount = 0
    Do k=1,scW
      Do j=1,scW
        Do i=1,scW
          scKey = i+scW*(j-1)+scW**2*(k-1)
          scKeyArr(scKey,1) = i
          scKeyArr(scKey,2) = j
          scKeyArr(scKey,3) = k
        End Do
      End Do
    End Do
! Loop through atoms
    Do i=1,nl(cKey)%coordsLength
! Make sub cell key
      scKey = SubCellKey(scAlat,scW,nl(cKey)%coords(i,1),nl(cKey)%coords(i,2),nl(cKey)%coords(i,3))
      nl(cKey)%subCellKey(i) = scKey
      scAtomCount(scKey) = scAtomCount(scKey) + 1
! Store in sub cell arrays
      scCoordsI(scKey,scAtomCount(scKey),1) = i                          ! unique atom id
      scCoordsI(scKey,scAtomCount(scKey),2) = nl(cKey)%labelID(i)        ! atom type
      Do j=1,3
        scCoordsR(scKey,scAtomCount(scKey),j) = nl(cKey)%coords(i,j)     ! atom coords
      End Do
    End Do
! Build key map (probably NOT needed)
    !Do i=1,scCount
    !  scKeyMap(i,1) = scAtomCount(i)
    !  scKeyMap(i,2) = i
    !End Do
    !Call sortArray(scKeyMap,"A")


    !Do i=1,scCount
  !    print *,i,scKeyMap(i,1),scKeyMap(i,2)
!    End Do

!    Do k=1,scW
!      Do j=1,scW
        !Do i=1,scW
!          scKey = i+scW*(j-1)+scW**2*(k-1)
!          print *,scAtomCount(scKey)
!        End Do
!      End Do
!    End Do



! Make neighbour list
! Loop through subcells
    Do scKeyA=1,scCount
! loop through surrounding and image sub cells
      Do l=-1,1
        Do m=-1,1
          Do n=-1,1
! Atom B subcell key
            xKey = scKeyArr(scKeyA,1)+l
            yKey = scKeyArr(scKeyA,2)+m
            zKey = scKeyArr(scKeyA,3)+n
! PBC subcell key
            scX = Modulus(xKey-1,scW)+1
            scY = Modulus(yKey-1,scW)+1
            scZ = Modulus(zKey-1,scW)+1
            scKeyB = scX+scW*(scY-1)+scW**2*(scZ-1)
            !scKeyB_P = (scX+l+1)+(scW+2)*(scY+m)+(scW+2)**2*(scZ+n)

! set default coord shift
            xShift = 0.0D0
            yShift = 0.0D0
            zShift = 0.0D0
            shiftArr = 0
! adjust shift
            inCell = 1
            If(xKey.lt.1)Then
              xShift = -1.0D0 * aLat
              shiftArr(1) = -1
              inCell = 0
            End If
            If(yKey.lt.1)Then
              yShift = -1.0D0 * aLat
              shiftArr(2) = -1
              inCell = 0
            End If
            If(zKey.lt.1)Then
              zShift = -1.0D0 * aLat
              shiftArr(3) = -1
              inCell = 0
            End If
            If(xKey.gt.scW)Then
              xShift = 1.0D0 * aLat
              shiftArr(1) = 1
              inCell = 0
            End If
            If(yKey.gt.scW)Then
              yShift = 1.0D0 * aLat
              shiftArr(2) = 1
              inCell = 0
            End If
            If(zKey.gt.scW)Then
              zShift = 1.0D0 * aLat
              shiftArr(3) = 1
              inCell = 0
            End If
! loop through A atoms
            Do atomA = 1, scAtomCount(scKeyA)
! loop through B atoms
              Do atomB = atomA, scAtomCount(scKeyB)
                atomA_ID = scCoordsI(scKeyA,atomA,1)
                atomB_ID = scCoordsI(scKeyB,atomB,1)
                If(atomA_ID.ne.atomB_ID)Then
                xA = scCoordsR(scKeyA,atomA,1)
                xB = scCoordsR(scKeyB,atomB,1)+xShift
                xD = xA-xB
                xDsq = xd**2
                If(xDsq.le.rVerletSQ)Then
                  yA = scCoordsR(scKeyA,atomA,2)
                  yB = scCoordsR(scKeyB,atomB,2)+yShift
                  yD = yA-yB
                  yDsq = yd**2
                  If(yDsq.le.rVerletSQ)Then
                    zA = scCoordsR(scKeyA,atomA,3)
                    zB = scCoordsR(scKeyB,atomB,3)+zShift
                    zD = zA-zB
                    zDsq = zd**2
                    If(zDsq.le.rVerletSQ)Then
                      rdSq = xDsq + yDsq + zDsq
                      If(rdSq.le.rVerletSq)Then
                        rd = sqrt(rdSq)
                        ! Key
                        nlKey = nlKey + 1
                        ! Overkey
                        nl(cKey)%nlOverKey(atomA_ID,0) = nl(cKey)%nlOverKey(atomA_ID,0) + 1
                        nl(cKey)%nlOverKey(atomA_ID,nl(cKey)%nlOverKey(atomA_ID,0)) = nlKey
                        nl(cKey)%nlOverKey(atomB_ID,0) = nl(cKey)%nlOverKey(atomB_ID,0) + 1
                        nl(cKey)%nlOverKey(atomB_ID,nl(cKey)%nlOverKey(atomB_ID,0)) = nlKey
                        ! Key/Type
                        nl(cKey)%atomA_ID(nlKey) = atomA_ID                   ! A ID
                        nl(cKey)%atomB_ID(nlKey) = atomB_ID                   ! B ID
                        nl(cKey)%atomA_Type(nlKey) = scCoordsI(scKeyA,atomA,2)  ! A type
                        nl(cKey)%atomB_Type(nlKey) = scCoordsI(scKeyB,atomB,2)  ! B type
                        nl(cKey)%inCell(nlKey) = inCell                     ! 1 if in cell, 0 if out of cell
                        nl(cKey)%atomPairKey = DoubleKey(scCoordsI(scKeyA,atomA,2),scCoordsI(scKeyB,atomB,2))
                        ! Displacement/Direction
                        nl(cKey)%rD(nlKey) = rD
                        nl(cKey)%vecAB(nlKey,1) = xD/rD                ! Vector from B to A (x)
                        nl(cKey)%vecAB(nlKey,2) = yD/rD                ! Vector from B to A (y)
                        nl(cKey)%vecAB(nlKey,3) = zD/rD                ! Vector from B to A (z)
                        ! super cell coords of atom B
                        nl(cKey)%subCell(nlKey,1) = shiftArr(1)
                        nl(cKey)%subCell(nlKey,2) = shiftArr(2)
                        nl(cKey)%subCell(nlKey,3) = shiftArr(3)
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
              End Do
            End Do
          End Do
        End Do
      End Do
    End Do
    nl(cKey)%length = nlKey
!----------------------------------------------
! Deallocate subcell arrays
    Deallocate(scAtomCount)
    Deallocate(scKeyArr)
    Deallocate(scCoordsI)
    Deallocate(scCoordsR)
    Deallocate(scKeyMap)
    Deallocate(uKeyArr)
!----------------------------------------------
        End If
      End If
    End Do
!----------------------------------------------
  End Subroutine makeNL_Broke




  Subroutine makeNL_Plain(nl, coords, rVerlet, cKey_NL_in)
! Make neighbour list for atoms
! The input coords and alat must be large enough so the same atom does not interact
! with a copy in a periodic cell surrounding the original
! e.g. if the rVerlet cutoff is 5, the alat must be greater than 5
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType), Dimension(:) :: nl
    Type(coordsType), Dimension(:)  :: coords
    Real(kind=DoubleReal) :: rVerlet
    Integer(kind=StandardInteger), Optional :: cKey_NL_in
! Vars:  Private
    Integer(kind=StandardInteger) :: cKey_NL
    Integer(kind=StandardInteger) :: cKey
    Integer(kind=StandardInteger) :: i, j, k
    Integer(kind=StandardInteger) :: l, m ,n
    Integer(kind=StandardInteger) :: scX, scY, scZ
    Integer(kind=StandardInteger) :: atomA, atomB, atomA_ID, atomB_ID, loopStartAtom
    Integer(kind=StandardInteger) :: coordLength, scW, scCount, scKey, maxAtomsPerSC
    Integer(kind=StandardInteger) :: scKeyA, scKeyB, scKeyA_Loop, scKeyB_Loop, scKeyB_P
    Real(kind=DoubleReal) :: aLat
    Real(kind=DoubleReal) :: rVerletSQ, scAlat
    Integer(kind=StandardInteger) :: nlKey, uKey
    Integer(kind=StandardInteger) :: xKey, yKey, zKey
    Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB
    Real(kind=DoubleReal) :: xD, yD, zD, rD
    Real(kind=DoubleReal) :: xDsq, yDsq, zDsq, rDsq
    Real(kind=DoubleReal) :: xShift, yShift, zShift
    Real(kind=DoubleReal), Dimension(1:3) :: position, positionT
    Integer(kind=StandardInteger), Dimension(1:3) :: shiftArr
    Integer(kind=StandardInteger) :: inCell
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
    If(Allocated(nl(cKey)%subCellKey))Then
      Deallocate(nl(cKey)%subCellKey)
    End If
    If(Allocated(nl(cKey)%forces))Then
      Deallocate(nl(cKey)%forces)
    End If
    If(Allocated(nl(cKey)%velocity))Then
      Deallocate(nl(cKey)%velocity)
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
    If(Allocated(nl(cKey)%vecAB))Then
      Deallocate(nl(cKey)%vecAB)
    End If
    If(Allocated(nl(cKey)%subCell))Then
      Deallocate(nl(cKey)%subCell)
    End If
    If(Allocated(nl(cKey)%cell))Then
      Deallocate(nl(cKey)%cell)
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
      nlArrayLength = 2000+&
        ceiling(0.35D0*(coords(cKey)%length)**2*(8.0D0*rVerlet**3)/((coords(cKey)%aLat)**3))
      nl(cKey)%arrayLength = nlArrayLength
!-------------------
! Estimate NL size
      nlArrayLength = 2000+&
        ceiling(0.35D0*(coords(cKey)%length)**2*(8.0D0*rVerlet**3)/((coords(cKey)%aLat)**3))
      nl(cKey)%arrayLength = nlArrayLength
      nlArrayOverKey = ceiling((3.0D0*nlArrayLength)/coords(cKey)%length)
! Coord details
      Allocate(nl(cKey)%label(1:coords(cKey)%length))
      Allocate(nl(cKey)%labelID(1:coords(cKey)%length))
      Allocate(nl(cKey)%coords(1:coords(cKey)%length,1:3))
      Allocate(nl(cKey)%subCellKey(1:coords(cKey)%length))
      Allocate(nl(cKey)%forces(1:coords(cKey)%length,1:3))
      Allocate(nl(cKey)%velocity(1:coords(cKey)%length,1:3))
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
      Allocate(nl(cKey)%vecAB(1:nlArrayLength,1:3))
! Subcell
      Allocate(nl(cKey)%subCell(1:nlArrayLength,1:3))
      Allocate(nl(cKey)%cell(1:nlArrayLength,1:3))
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
! Load coords
    Do i=1,nl(cKey)%coordsLength
      nl(cKey)%label(i) = coords(cKey)%label(i)
      nl(cKey)%labelID(i) = coords(cKey)%labelID(i)
      nl(cKey)%charge(i) = coords(cKey)%charge(i)
      nl(cKey)%mass(i) = coords(cKey)%mass(i)
      Do j=1,3
        position(j) = coords(cKey)%coords(i,j)   ! x/y/z coord
        position(j) = Modulus(position(j),1.0D0) ! Keep within a 1x1x1 box
      End Do
      positionT = MatMul(nl(cKey)%unitCellScaled,position)  ! Transform
      Do j=1,3
        nl(cKey)%coords(i,j) = positionT(j)
        nl(cKey)%forces(i,j) = 0.0D0
        nl(cKey)%velocity(i,j) = 0.0D0
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
    nl(cKey)%rMin = rVerlet
    nl(cKey)%rMax = 0.0D0
    nlKey = 0
!----------------------------------------------
! Init overkey count
    Do i = 1,nl(cKey)%coordsLength
      nl(cKey)%nlOverKey(i,0) = 0
    End Do

! Loop through atoms
    Do atomA_ID=1,nl(cKey)%coordsLength
      Do atomB_ID=atomA_ID,nl(cKey)%coordsLength
        Do l=-1,1
          Do m=-1,1
            Do n=-1,1
! set default coord shift
              xShift = 1.0D0*l*aLat
              yShift = 1.0D0*m*aLat
              zShift = 1.0D0*n*aLat
              If(atomA_ID.ne.atomB_ID)Then
                xA = nl(cKey)%coords(atomA_ID,1)
                xA = nl(cKey)%coords(atomB_ID,1)+xShift
                xD = xA-xB
                xDsq = xd**2
                If(xDsq.le.rVerletSQ)Then
                  yA = nl(cKey)%coords(atomA_ID,2)
                  yB = nl(cKey)%coords(atomB_ID,2)+yShift
                  yD = yA-yB
                  yDsq = yd**2
                  If(yDsq.le.rVerletSQ)Then
                    zA = nl(cKey)%coords(atomA_ID,3)
                    zB = nl(cKey)%coords(atomB_ID,3)+zShift
                    zD = zA-zB
                    zDsq = zd**2
                    If(zDsq.le.rVerletSQ)Then
                      rdSq = xDsq + yDsq + zDsq
                      If(rdSq.le.rVerletSq)Then
                        rd = sqrt(rdSq)
                        ! Key
                        nlKey = nlKey + 1
                        ! Overkey
                        nl(cKey)%nlOverKey(atomA_ID,0) = nl(cKey)%nlOverKey(atomA_ID,0) + 1
                        nl(cKey)%nlOverKey(atomA_ID,nl(cKey)%nlOverKey(atomA_ID,0)) = nlKey
                        nl(cKey)%nlOverKey(atomB_ID,0) = nl(cKey)%nlOverKey(atomB_ID,0) + 1
                        nl(cKey)%nlOverKey(atomB_ID,nl(cKey)%nlOverKey(atomB_ID,0)) = nlKey
                        ! Key/Type
                        nl(cKey)%atomA_ID(nlKey) = atomA_ID                   ! A ID
                        nl(cKey)%atomB_ID(nlKey) = atomB_ID                   ! B ID
                        nl(cKey)%atomA_Type(nlKey) = nl(cKey)%labelID(atomA_ID)  ! A type
                        nl(cKey)%atomB_Type(nlKey) = nl(cKey)%labelID(atomB_ID)  ! B type
                        nl(cKey)%inCell(nlKey) = inCell                     ! 1 if in cell, 0 if out of cell
                        nl(cKey)%atomPairKey = DoubleKey(atomA_ID,atomB_ID)
                        ! Displacement/Direction
                        nl(cKey)%rD(nlKey) = rD
                        nl(cKey)%vecAB(nlKey,1) = xD/rD                ! Vector from B to A (x)
                        nl(cKey)%vecAB(nlKey,2) = yD/rD                ! Vector from B to A (y)
                        nl(cKey)%vecAB(nlKey,3) = zD/rD                ! Vector from B to A (z)
                        ! Store rmin/rmax
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
            End Do
          End Do
        End Do
      End Do
    End Do
    nl(cKey)%length = nlKey
!----------------------------------------------
        End If
      End If
    End Do
!----------------------------------------------
  End Subroutine makeNL_Plain






  Subroutine updateNL(nl, cKey_NL_in, distortion_in)
! Update neighbour list without rebuilding
! after applying a distortion to all coordinates
! Uses fractional coords - static calculations
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType), Dimension(:) :: nl
    Integer(kind=StandardInteger), Optional :: cKey_NL_in
    Real(kind=DoubleReal), Dimension(1:3,1:3), Optional :: distortion_in
! Vars:  Private
    Integer(kind=StandardInteger) :: n, cKey, cKey_NL
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: distortion
    Real(kind=DoubleReal) :: aLat
    Real(kind=DoubleReal) :: xD, yD, zD
    Real(kind=DoubleReal), Dimension(1:3) :: aVec, bVec
    Logical :: cellDistortion
!----------------------------------------------
! Optional arguments
!----------------------------------------------
    cKey_NL = 0
    If(Present(cKey_NL_in))Then
      cKey_NL = cKey_NL_in
    End If
    distortion = IdentityMatrix(distortion)
    cellDistortion = .false.
    If(Present(distortion_in))Then
      distortion = distortion_in
      cellDistortion = .true.
    End If
!----------------------------------------------
! Loop through neighbour lists
!----------------------------------------------
    Do cKey=1,size(nl,1)
      If((cKey_NL.eq.0).or.(cKey_NL.eq.cKey))Then
!------------------------
    aLat = nl(cKey)%aLat
! Loop through pairs
    Do n=1,nl(cKey)%length
      aVec(1) = nl(cKey)%coords(nl(cKey)%atomA_ID(n),1)
      aVec(2) = nl(cKey)%coords(nl(cKey)%atomA_ID(n),2)
      aVec(3) = nl(cKey)%coords(nl(cKey)%atomA_ID(n),3)
      bVec(1) = (nl(cKey)%coords(nl(cKey)%atomB_ID(n),1)+nl(cKey)%subCell(n,1))
      bVec(2) = (nl(cKey)%coords(nl(cKey)%atomB_ID(n),2)+nl(cKey)%subCell(n,2))
      bVec(3) = (nl(cKey)%coords(nl(cKey)%atomB_ID(n),3)+nl(cKey)%subCell(n,3))
! Distort position vectors
      If(cellDistortion)Then
        aVec = MatMul(distortion, aVec)
        bVec = MatMul(distortion, bVec)
      End If
! Update atom seperation
      xD = aLat*(aVec(1)-bVec(1)) ! xA - xB
      yD = aLat*(aVec(2)-bVec(2)) ! yA - yA
      zD = aLat*(aVec(3)-bVec(3)) ! zA - zB
      nl(cKey)%rD(n) = sqrt(xD**2+yD**2+zD**2)
! Update atom A to atom B vector
      nl(cKey)%vecAB(n,1) = xD/nl(cKey)%rD(n)
      nl(cKey)%vecAB(n,2) = yD/nl(cKey)%rD(n)
      nl(cKey)%vecAB(n,3) = zD/nl(cKey)%rD(n)
    End Do
!------------------------
      End If
    End Do
  End Subroutine updateNL


  Subroutine refreshNL(nl, cKey_NL_in)
! Refresh neighbour list using coodsMD
! without rebuilding entire list
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType), Dimension(:) :: nl
    Integer(kind=StandardInteger), Optional :: cKey_NL_in
! Vars:  Private
    Integer(kind=StandardInteger) :: n, cKey, cKey_NL
    Real(kind=DoubleReal) :: xA, yA, zA
    Real(kind=DoubleReal) :: xB, yB, zB
    Real(kind=DoubleReal) :: xD, yD, zD
!----------------------------------------------
! Optional arguments
!----------------------------------------------
    cKey_NL = 0
    If(Present(cKey_NL_in))Then
      cKey_NL = cKey_NL_in
    End If
!----------------------------------------------
! Loop through neighbour lists
!----------------------------------------------
    Do cKey=1,size(nl,1)
      If((cKey_NL.eq.0).or.(cKey_NL.eq.cKey))Then
!------------------------
! Loop through neighbour list
    Do n=1,nl(cKey)%length
! Calculate atom seperation
! atom positions
      xA = nl(ckey)%coordsMD(nl(ckey)%atomA_ID(n),1)
      yA = nl(ckey)%coordsMD(nl(ckey)%atomA_ID(n),2)
      zA = nl(ckey)%coordsMD(nl(ckey)%atomA_ID(n),3)
      xB = nl(ckey)%coordsMD(nl(ckey)%atomB_ID(n),1)+nl(cKey)%aLat*nl(cKey)%subCell(n,1)
      yB = nl(ckey)%coordsMD(nl(ckey)%atomB_ID(n),2)+nl(cKey)%aLat*nl(cKey)%subCell(n,2)
      zB = nl(ckey)%coordsMD(nl(ckey)%atomB_ID(n),3)+nl(cKey)%aLat*nl(cKey)%subCell(n,3)
! atom displacements
      xD = xA - xB
      yD = yA - yB
      zD = zA - zB
      nl(ckey)%rD(n) = sqrt(xD**2+yD**2+zD**2)
      !If(n.eq.10)Then
      !  print *,"Alat ",nl(cKey)%aLat
      !  print *,"A ",xA,yA,zA
      !  print *,"B ",xB,yB,zB
      !  print *,"SC B ",nl(cKey)%subCell(n,1)
      !  print *,"Bx new ",nl(ckey)%coordsMD(nl(ckey)%atomB_ID(n),1)+nl(cKey)%aLat*nl(cKey)%subCell(n,1)
      !  print *,"X Y Z  ",xD,yD,zD
      !  print *,"RD ",nl(ckey)%rD(n)
      !End If
! Vectors
      nl(ckey)%vecAB(n,1) = xD/nl(ckey)%rD(n)
      nl(ckey)%vecAB(n,2) = yD/nl(ckey)%rD(n)
      nl(ckey)%vecAB(n,3) = zD/nl(ckey)%rD(n)
    End Do
!------------------------
      End If
    End Do
  End Subroutine refreshNL






! ------------------------------------------------------------
!    GEOM.NL Functions
! ------------------------------------------------------------

  Function SubCellKey (scAlat, scW, x, y, z) Result (key)
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
  End Function SubCellKey
