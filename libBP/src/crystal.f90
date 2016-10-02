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


Module crystalTypes
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
  Public :: oCrystal, oCrystals
!------------------------------------------------------------------------------
! Defined Types
!------------------------------------------------------------------------------
! oCoord (label, labelID, xyz, force, velocity)


  Type :: oCrystal
    Integer(kind=StandardInteger) :: batch = 0    ! key for the batch
    Integer(kind=StandardInteger) :: batchID = 0  ! key for crystal in batch
    Integer(kind=StandardInteger) :: batchCount = 0
    Integer(kind=StandardInteger) :: coordInCount = 0
    Integer(kind=StandardInteger) :: coordCount = 0
    Character(len=128) :: outputDir
    Real(kind=DoubleReal) :: aLat
    Real(kind=DoubleReal) :: aLatInverse
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: unitCell
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: aLat_unitCell
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: unitCellInverse
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: aLat_unitCellInverse
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: transformToExpanded
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: transformFromExpanded
    Integer(kind=StandardInteger), Dimension(1:3) :: copy  ! x,y,z
    Real(kind=DoubleReal), Dimension(1:2) :: dRange = -1.0D0
    Integer(kind=StandardInteger) :: seed = 0
    Real(kind=DoubleReal) :: randVariation
    Type(oCoord), Allocatable, Dimension(:) :: coordsIn
    Type(oCoord), Allocatable, Dimension(:) :: coords
    Type(oCoord), Allocatable, Dimension(:) :: coordsFractional
  End Type oCrystal

  Type :: oCrystals
    Integer(kind=StandardInteger) :: crystalCount
    Type(oCrystal), Allocatable, Dimension(:) :: crystal
  End Type oCrystals

End Module crystalTypes


Module crystal
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use rng
  Use constants
  Use strings
  Use general
  Use basicMaths, Only: Modulus
  Use units
  Use matrix, Only: InvertMatrix
  Use rngDistTypes, Only: oRandDist
  Use programMod
  Use geom, Only: expandCoords_UnitCell, geomTransformCoords, perturbCoords
  Use crystalTypes
! Force declaration of all variables
  Implicit None
! Make private
  Private
! Public
  Public :: runCrystal
! Interfaces
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! ------------------------------------------------------------------------!
! Crystal Subroutines
! ------------------------------------------------------------------------!

  Subroutine runCrystal()
! Run crystal subroutine
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, k, configs, totalConfigs
    Integer(kind=StandardInteger) :: configID, configID_L, coordID, batch
    Type(oCrystals) :: cData
    Character(len=32), Dimension(1:10) :: fieldArray
    Integer(kind=StandardInteger) :: fieldCount, tempInt
    Real(kind=DoubleReal) :: tempDP, maxPerturbation, randNum
    Type(oRandDist) :: randDistObj
    Real(kind=DoubleReal), Dimension(1:3) :: modPeriod
    Integer(kind=LongInteger) :: randSeed
!
    print *,"CRYSTAL"
! Init
    modPeriod(1) = 1.0D0
    modPeriod(2) = 1.0D0
    modPeriod(3) = 1.0D0

!-----------------------------------
! Count crystal configs
!-----------------------------------
    configs = 0
    totalConfigs = 0
    Do i=1,programObj%programFiles(1)%rowsP
! Read row
      Call readFieldsCharacter(programObj%programFiles(1)%dataP_UC(i),fieldArray,fieldCount)
! Count
      If(fieldArray(1)(1:4).eq."#NEW")Then
        configs = configs + 1
        totalConfigs = totalConfigs + 1
      End If
      If(fieldArray(1)(1:2).eq."#G")Then
        tempInt = StrToInt(fieldArray(2))
        totalConfigs = totalConfigs -1 + tempInt
      End If
    End Do
    print *,configs,totalConfigs
! Store count + allocate array
    cData%crystalCount = totalConfigs
    Allocate(cData%crystal(1:cData%crystalCount))
!-----------------------------------
! Load Data
!-----------------------------------
    configID = 0
    batch = 0
! Read Input File - File 1
    Do i=1,programObj%programFiles(1)%rowsP
! Read row
      Call readFieldsCharacter(programObj%programFiles(1)%dataP_UC(i),fieldArray,fieldCount)
!---------------------
! File Keywords
!---------------------
      If(fieldArray(1)(1:4).eq."#NEW")Then
        batch = batch + 1
        configID = configID + 1
        cData%crystal(configID)%batch = batch
        cData%crystal(configID)%batchID = 1
        cData%crystal(configID)%coordInCount = 0
      End If
      If(fieldArray(1)(1:7).eq."#OUTDIR")Then
        cData%crystal(configID)%outputDir = BlankString(cData%crystal(configID)%outputDir)
        cData%crystal(configID)%outputDir = programObj%programFiles(1)%dataP(i)(8:128)
        cData%crystal(configID)%outputDir = trim(adjustl(cData%crystal(configID)%outputDir))
        cData%crystal(configID)%outputDir = fillPath(cData%crystal(configID)%outputDir)  ! general.f90
      End If
      If(fieldArray(1)(1:3).eq."#LP")Then
        tempDP = StrToDP(fieldArray(2))
        cData%crystal(configID)%aLat = tempDP
      End If
      If(fieldArray(1)(1:2).eq."#X")Then
        tempDP = StrToDP(fieldArray(2))
        cData%crystal(configID)%unitCell(1,1) = tempDP
        tempDP = StrToDP(fieldArray(3))
        cData%crystal(configID)%unitCell(1,2) = tempDP
        tempDP = StrToDP(fieldArray(4))
        cData%crystal(configID)%unitCell(1,3) = tempDP
      End If
      If(fieldArray(1)(1:2).eq."#Y")Then
        tempDP = StrToDP(fieldArray(2))
        cData%crystal(configID)%unitCell(2,1) = tempDP
        tempDP = StrToDP(fieldArray(3))
        cData%crystal(configID)%unitCell(2,2) = tempDP
        tempDP = StrToDP(fieldArray(4))
        cData%crystal(configID)%unitCell(2,3) = tempDP
      End If
      If(fieldArray(1)(1:2).eq."#Z")Then
        tempDP = StrToDP(fieldArray(2))
        cData%crystal(configID)%unitCell(3,1) = tempDP
        tempDP = StrToDP(fieldArray(3))
        cData%crystal(configID)%unitCell(3,2) = tempDP
        tempDP = StrToDP(fieldArray(4))
        cData%crystal(configID)%unitCell(3,3) = tempDP
      End If
      If(fieldArray(1)(1:3).eq."#CC")Then
        tempInt = StrToInt(fieldArray(2))
        cData%crystal(configID)%copy(1) = tempInt
        tempInt = StrToInt(fieldArray(3))
        cData%crystal(configID)%copy(2) = tempInt
        tempInt = StrToInt(fieldArray(4))
        cData%crystal(configID)%copy(3) = tempInt
      End If
      If(fieldArray(1)(1:2).eq."#D")Then
        cData%crystal(configID)%dRange(1) = -1.0D0
        cData%crystal(configID)%dRange(2) = -1.0D0
        tempDP = StrToDP(fieldArray(2))
        cData%crystal(configID)%dRange(1) = UnitConvert(tempDP,fieldArray(3),"Ang")
        If(fieldCount.eq.5)Then
          tempDP = StrToDP(fieldArray(4))
          cData%crystal(configID)%dRange(2) = UnitConvert(tempDP,fieldArray(5),"Ang")
        End If
      End If
      If(fieldArray(1)(1:2).eq."#S")Then
        tempInt = StrToInt(fieldArray(2))
        cData%crystal(configID)%seed = tempInt
      End If
      If(fieldArray(1)(1:2).eq."#G")Then
        tempInt = StrToInt(fieldArray(2))
        cData%crystal(configID)%batchCount = tempInt
      End If
!-----
      If(fieldArray(1)(1:1).ne."#".and.fieldCount.eq.4)Then
        cData%crystal(configID)%coordInCount = cData%crystal(configID)%coordInCount + 1
      End If
!-----
      If(fieldArray(1)(1:4).eq."#END")Then
! Calculate alat/unit vars
        cData%crystal(configID)%aLatInverse = 1.0D0/cData%crystal(configID)%aLat
        cData%crystal(configID)%unitCellInverse = InvertMatrix(cData%crystal(configID)%unitCell)
        cData%crystal(configID)%aLat_unitCell = cData%crystal(configID)%aLat*&
          cData%crystal(configID)%unitCell
        cData%crystal(configID)%aLat_unitCellInverse = cData%crystal(configID)%aLatInverse*&
          cData%crystal(configID)%unitCellInverse
        Do j=1,3
          Do k=1,3
            cData%crystal(configID)%transformToExpanded(j,k) = &
              cData%crystal(configID)%aLat_unitCell(j,k) * cData%crystal(configID)%copy(k)
            cData%crystal(configID)%transformFromExpanded(j,k) = &
              cData%crystal(configID)%aLat_unitCellInverse(j,k) * (1.0D0/cData%crystal(configID)%copy(k))
          End Do
        End Do
! Allocate coords array
        Allocate(cData%crystal(configID)%coordsIn(1:cData%crystal(configID)%coordInCount))
! Loop through batch
        If(cData%crystal(configID)%batchCount.gt.1)Then
          Do j=1,cData%crystal(configID)%batchCount-1
            configID_L = configID
            configID = configID + 1
            cData%crystal(configID)%batch = cData%crystal(configID_L)%batch
            cData%crystal(configID)%batchID = cData%crystal(configID_L)%batchID+1
            cData%crystal(configID)%batchCount = cData%crystal(configID_L)%batchCount
            cData%crystal(configID)%coordInCount = cData%crystal(configID_L)%coordInCount
! Copy from previous
            cData%crystal(configID)%outputDir = cData%crystal(configID_L)%outputDir
            cData%crystal(configID)%aLat = cData%crystal(configID_L)%aLat
            cData%crystal(configID)%aLatInverse = cData%crystal(configID_L)%aLatInverse
            cData%crystal(configID)%unitCell = cData%crystal(configID_L)%unitCell
            cData%crystal(configID)%aLat_unitCell = cData%crystal(configID_L)%aLat_unitCell
            cData%crystal(configID)%unitCellInverse = cData%crystal(configID_L)%unitCellInverse
            cData%crystal(configID)%aLat_unitCellInverse = cData%crystal(configID_L)%aLat_unitCellInverse
            cData%crystal(configID)%transformToExpanded = cData%crystal(configID_L)%transformToExpanded
            cData%crystal(configID)%transformFromExpanded = cData%crystal(configID_L)%transformFromExpanded
            cData%crystal(configID)%copy = cData%crystal(configID_L)%copy
            cData%crystal(configID)%dRange = cData%crystal(configID_L)%dRange
            cData%crystal(configID)%seed = cData%crystal(configID_L)%seed
! Allocate coords array
            Allocate(cData%crystal(configID)%coordsIn(1:cData%crystal(configID)%coordInCount))
          End Do
        End If
      End If
    End Do
!-----------------------------------
! Load Coords
!-----------------------------------
    configID = 0
! Read Input File - File 1
    Do i=1,programObj%programFiles(1)%rowsP
! Read row
      Call readFieldsCharacter(programObj%programFiles(1)%dataP_UC(i),fieldArray,fieldCount)
!-----
      If(fieldArray(1)(1:4).eq."#NEW")Then
        configID = configID + 1
        coordID = 0
      End If
!-----
      If(fieldArray(1)(1:1).ne."#".and.fieldCount.eq.4)Then
        coordID = coordID + 1
        Call readFieldsCharacter(programObj%programFiles(1)%dataP(i),fieldArray,fieldCount)
        cData%crystal(configID)%coordsIn(coordID)%label = StrReplace(fieldArray(1),char(0),char(32))
        !cData%crystal(configID)%coordsIn(coordID)%label = &
        !  StrReplace(cData%crystal(configID)%coordsIn(coordID)%label),char(0),char(32))
        tempDP = StrToDP(fieldArray(2))
        cData%crystal(configID)%coordsIn(coordID)%xyz(1) = tempDP
        tempDP = StrToDP(fieldArray(3))
        cData%crystal(configID)%coordsIn(coordID)%xyz(2) = tempDP
        tempDP = StrToDP(fieldArray(4))
        cData%crystal(configID)%coordsIn(coordID)%xyz(3) = tempDP
      End If
!-----
      If(fieldArray(1)(1:4).eq."#END")Then
! Allocate expanded "coords" array
        cData%crystal(configID)%coordCount = cData%crystal(configID)%copy(1)*&
          cData%crystal(configID)%copy(2)*cData%crystal(configID)%copy(3)*&
          cData%crystal(configID)%coordInCount
        Allocate(cData%crystal(configID)%coords(1:cData%crystal(configID)%coordCount))
        Allocate(cData%crystal(configID)%coordsFractional(1:cData%crystal(configID)%coordCount))
        If(cData%crystal(configID)%batchCount.gt.1)Then
          Do j=1,cData%crystal(configID)%batchCount-1
            configID_L = configID
            configID = configID + 1
! Allocate expanded "coords" array
            cData%crystal(configID)%coordCount = cData%crystal(configID)%copy(1)*&
              cData%crystal(configID)%copy(2)*cData%crystal(configID)%copy(3)*&
              cData%crystal(configID)%coordInCount
            Allocate(cData%crystal(configID)%coords(1:cData%crystal(configID)%coordCount))
            Allocate(cData%crystal(configID)%coordsFractional(1:cData%crystal(configID)%coordCount))
! Loop through coords
            Do coordID = 1,cData%crystal(configID)%coordInCount
              cData%crystal(configID)%coordsIn(coordID)%label = &
                cData%crystal(configID_L)%coordsIn(coordID)%label
              cData%crystal(configID)%coordsIn(coordID)%xyz(1) = &
                cData%crystal(configID_L)%coordsIn(coordID)%xyz(1)
              cData%crystal(configID)%coordsIn(coordID)%xyz(2) = &
                cData%crystal(configID_L)%coordsIn(coordID)%xyz(2)
              cData%crystal(configID)%coordsIn(coordID)%xyz(3) = &
                cData%crystal(configID_L)%coordsIn(coordID)%xyz(3)
            End Do
          End Do
        End If
      End If
    End Do
!-----------------------------------
! Expand, Transform, Randomise
!-----------------------------------
! Expand Coords
    Do configID=1,cData%crystalCount
      Call expandCoords_UnitCell(cData%crystal(configID)%coordsIn, cData%crystal(configID)%coords, cData%crystal(configID)%copy)
    End Do
! Transform coords
    Do configID=1,cData%crystalCount
      Call geomTransformCoords(cData%crystal(configID)%coords, cData%crystal(configID)%unitCell, &
        cData%crystal(configID)%aLat, cData%crystal(configID)%copy)
    End Do
! Randomise
    randDistObj%distType = "GHEAT"
    randDistObj%lowerBound = -1.0D0
    randDistObj%upperBound = 1.0D0
    Do configID=1,cData%crystalCount
      cData%crystal(configID)%randVariation = 0.0D0
      If(cData%crystal(configID)%dRange(1).ge.0.0D0)Then ! If random variation selected
! Set seed (use configID as seed)
        randSeed = 0
        If(cData%crystal(configID)%seed.eq.0)Then   ! If 0, use configID as the seed
          randSeed = configID
        Else If(cData%crystal(configID)%seed.lt.0)Then  ! -1 or lower, random seed
          randSeed = RandomSeed()
        Else
          randSeed = cData%crystal(configID)%seed
        End If
        randNum = RandomLCG(randSeed)
        If(cData%crystal(configID)%dRange(2).ge.0.0D0)Then
          maxPerturbation = cData%crystal(configID)%dRange(1) + &
           randNum * (cData%crystal(configID)%dRange(2)-cData%crystal(configID)%dRange(1))
        Else
          maxPerturbation = cData%crystal(configID)%dRange(1)
        End If
        cData%crystal(configID)%randVariation = maxPerturbation
        Call perturbCoords(cData%crystal(configID)%coords, randDistObj, 0.0D0, maxPerturbation)
      End If
    End Do
! Modulus
    Do configID=1,cData%crystalCount
      Do i=1,3
        modPeriod(i) = cData%crystal(configID)%aLat*cData%crystal(configID)%copy(i)
      End Do
      Do coordID = 1,cData%crystal(configID)%coordCount
        Do i=1,3
          cData%crystal(configID)%coords(coordID)%xyz(i) = &
            Modulus(cData%crystal(configID)%coords(coordID)%xyz(i),modPeriod(i))
        End Do
      End Do
    End Do
! Make fractional
    Do configID=1,cData%crystalCount
      Do coordID = 1,cData%crystal(configID)%coordCount
        cData%crystal(configID)%coordsFractional(coordID)%xyz = &
          MatMul(cData%crystal(configID)%transformFromExpanded,cData%crystal(configID)%coords(coordID)%xyz)
      End Do
    End Do
!-----------------------------------
! Output
!-----------------------------------

    Call outputSummary(cData)
    Call outputCrystalFiles(cData)
    !Call outputCrystalPrint(cData)


  End Subroutine runCrystal

  Subroutine outputSummary(cData)
! Run crystal subroutine
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oCrystals) :: cData
! Vars:  Private
    Integer(kind=StandardInteger) :: configID
! Loop through
    Do configID=1,cData%crystalCount
      print *,configID,cData%crystal(configID)%randVariation
    End Do
  End Subroutine outputSummary


  Subroutine outputCrystalFiles(cData)
! Run crystal subroutine
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oCrystals) :: cData
! Vars:  Private
    Integer(kind=StandardInteger) :: configID, coordID
    Character(len=2) :: batchID_Str
    Character(len=4) :: batchCount_Str
    Character(len=16) :: fileName
    Character(len=12) :: batchName
    Character(len=128) :: filePath
! Make fractional coord files
    Do configID=1,cData%crystalCount
! Make dir
      Call makeDir(cData%crystal(configID)%outputDir)
! Make file name
      batchID_Str = IntToStrL(cData%crystal(configID)%batch, 2)
      batchCount_Str = IntToStrL(cData%crystal(configID)%batchID, 4)
      fileName = BlankString(fileName)
      fileName(1:5) = "Frac_"
      fileName(6:7) = batchID_Str
      fileName(8:8) = "_"
      fileName(9:12) = batchCount_Str
      fileName(13:16) = ".dat"
      filePath = trim(cData%crystal(configID)%outputDir)//"/"//fileName
      open(unit=999,file=trim(filePath))
      Do coordID = 1,cData%crystal(configID)%coordCount
        write(999,"(A16,F12.5,F12.5,F12.5)") cData%crystal(configID)%coords(coordID)%label,&
        cData%crystal(configID)%coordsFractional(coordID)%xyz(1),&
        cData%crystal(configID)%coordsFractional(coordID)%xyz(2),&
        cData%crystal(configID)%coordsFractional(coordID)%xyz(3)
      End Do
      close(999)
    End Do
! Make summary files
    Do configID=1,cData%crystalCount
      If(cData%crystal(configID)%batchID.eq.1)Then
        batchID_Str = IntToStrL(cData%crystal(configID)%batch, 2)
        batchName(1:6) = "Batch_"
        batchName(7:8) = batchID_Str
        batchName(9:12) = ".dat"
        filePath = trim(cData%crystal(configID)%outputDir)//"/"//batchName
        open(unit=999,file=trim(filePath))
        write(999,"(I8,F12.5)") configID,cData%crystal(configID)%randVariation
      Else
        write(999,"(I8,F12.5)") configID,cData%crystal(configID)%randVariation
      End If
      If(cData%crystal(configID)%batchID.eq.cData%crystal(configID)%batchCount)Then
        close(999)
      End If
    End Do

  End Subroutine outputCrystalFiles


  Subroutine outputCrystalPrint(cData)
! Run crystal subroutine
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oCrystals) :: cData
! Vars:  Private
    Integer(kind=StandardInteger) :: configID, coordID
! Print
    Do configID=1,cData%crystalCount
      print *,"Crystal  ",configID
      print *,"...Batch ID             ",cData%crystal(configID)%batchID," of ",cData%crystal(configID)%batchCount
      print *,"...aLat                 ",cData%crystal(configID)%batchID," of ",cData%crystal(configID)%batchCount
      print *,"...coord count in       ",cData%crystal(configID)%coordInCount
      print *,"...coord count expanded ",cData%crystal(configID)%coordCount
      Do coordID = 1,cData%crystal(configID)%coordInCount
        print *,coordID,cData%crystal(configID)%coordsIn(coordID)%label,&
        cData%crystal(configID)%coordsIn(coordID)%xyz(1),&
        cData%crystal(configID)%coordsIn(coordID)%xyz(2),&
        cData%crystal(configID)%coordsIn(coordID)%xyz(3)
      End Do
      print *,"-----------------------------------------"
      Do coordID = 1,cData%crystal(configID)%coordCount
        print *,coordID,cData%crystal(configID)%coords(coordID)%label,&
        cData%crystal(configID)%coords(coordID)%xyz(1),&
        cData%crystal(configID)%coords(coordID)%xyz(2),&
        cData%crystal(configID)%coords(coordID)%xyz(3),"            ",&
        cData%crystal(configID)%coordsFractional(coordID)%xyz(1),&
        cData%crystal(configID)%coordsFractional(coordID)%xyz(2),&
        cData%crystal(configID)%coordsFractional(coordID)%xyz(3)
      End Do
      print *,""
    End Do
  End Subroutine outputCrystalPrint


End Module crystal

!-----------------------------------------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------------------------------------
