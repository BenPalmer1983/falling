! ------------------------------------------------------------
!               POTENTIALFITTING: Load Files
!                   Load fitting files
! ------------------------------------------------------------



  Subroutine pfLoadFiles(pfObj)
! Load potential fitting files:
!   main file
!   config files
!   DFT files
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oPotentialFitting) :: pfObj
! Vars:  Private
! Load Arguments - input file name/path
    Call loadCmdArgs()
    pfObj%inputFile = trim(cmdArgs%args(1))
    Call getcwd(pfObj%recordFile)
    pfObj%recordFile = trim(pfObj%recordFile)//"/potentialFit.configs"
! Count configs
    Call countInputConfigs(pfObj)
! Load input file
    Call loadInputFile(pfObj)
! Load config files
    Call loadConfigFiles(pfObj)
! Make dirs
    Call pfMakeDirs(pfObj)
! Count configs and make config make (on the offchance some have zero coord length)
    Call pfConfigCount(pfObj)
! Output File
    Call pfConfigToFile(pfObj)
!--------------------------------------------------------
  End Subroutine pfLoadFiles



  Subroutine countInputConfigs(pfObj)
! Start potential fitting
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oPotentialFitting) :: pfObj
! Vars:  Private
    Integer(kind=StandardInteger) :: i, n
    Integer(kind=StandardInteger) :: fileLength
    Character(Len=64), Dimension(1:256) :: fileArray
    Character(Len=64) :: fileRow, nextRow
!--------------------------------------------------------
! Read potential input file
    Call readFile(pfObj%inputFile, fileArray, fileLength)
!--------------------------------------------------------
! Loop through input file rows and count configs
    n = 0
    pfObj%cFileCount = 0
    pfObj%dftFileCount = 0
    Do While(n.le.fileLength)
      n = n + 1
      fileRow = StrToUpper(fileArray(n))
! config files
      If(fileRow(1:7).eq."#CFILES")Then
        Do i=1,fileLength
          nextRow = Trim(fileArray(n+1))
          If((n+1.gt.fileLength).or.(nextRow(1:1).eq."#"))Then
            Exit
          Else
            n = n + 1
            fileRow = Trim(fileArray(n))
            pfObj%cFileCount = pfObj%cFileCount + 1
          End If
        End Do
      End If
! DFT files
      If(fileRow(1:9).eq."#DFTFILES")Then
        Do i=1,fileLength
          nextRow = Trim(fileArray(n+1))
          If((n+1.gt.fileLength).or.(nextRow(1:1).eq."#"))Then
            Exit
          Else
            n = n + 1
            fileRow = Trim(fileArray(n))
            pfObj%dftFileCount = pfObj%dftFileCount + 1
          End If
        End Do
      End If
    End Do
    pfObj%fileCount = pfObj%cFileCount + pfObj%dftFileCount
!--------------------------------------------------------
  End Subroutine countInputConfigs


  Subroutine loadInputFile(pfObj)
! Start potential fitting
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oPotentialFitting) :: pfObj
! Vars:  Private
    Integer(kind=StandardInteger) :: i, n, fileCount
    Integer(kind=StandardInteger) :: fileLength
    Character(Len=64), Dimension(1:256) :: fileArray
    Character(Len=64) :: fileRow, nextRow
!--------------------------------------------------------
! Read potential input file
    Call readFile(pfObj%inputFile, fileArray, fileLength)
!--------------------------------------------------------
! Allocate
    If(Allocated(pfObj%fileList))Then
      Deallocate(pfObj%fileList)
    End If
    Allocate(pfObj%fileList(1:pfObj%fileCount))
    If(Allocated(pfObj%configs))Then
      Deallocate(pfObj%configs)
    End If
    Allocate(pfObj%configs(1:pfObj%fileCount))
!--------------------------------------------------------
! Loop through input file rows and load data
!--------------------------------------------------------
! Loop through input file rows and count configs
    n = 0
    fileCount = 0
    Do While(n.le.fileLength)
      n = n + 1
      fileRow = StrToUpper(fileArray(n))
! potential file
      If(fileRow(1:10).eq."#OUTPUTDIR")Then
        n = n + 1
        pfObj%outputDirectory = fileArray(n)
      End If
! potential file
      If(fileRow(1:8).eq."#POTFILE")Then
        n = n + 1
        pfObj%potentialFile = fileArray(n)
      End If
! config files
      If(fileRow(1:7).eq."#CFILES")Then
        Do i=1,fileLength
          nextRow = Trim(fileArray(n+1))
          If((n+1.gt.fileLength).or.(nextRow(1:1).eq."#"))Then
            Exit
          Else
            n = n + 1
            fileRow = Trim(fileArray(n))
            fileCount = fileCount + 1
            pfObj%fileList(fileCount)%filePath = Trim(fileRow)
            pfObj%fileList(fileCount)%fileType = 1
          End If
        End Do
      End If
! DFT files
      If(fileRow(1:9).eq."#DFTFILES")Then
        Do i=1,fileLength
          nextRow = Trim(fileArray(n+1))
          If((n+1.gt.fileLength).or.(nextRow(1:1).eq."#"))Then
            Exit
          Else
            n = n + 1
            fileRow = Trim(fileArray(n))
            fileCount = fileCount + 1
            pfObj%fileList(fileCount)%filePath = Trim(fileRow)
            pfObj%fileList(fileCount)%fileType = 2
          End If
        End Do
      End If
    End Do
    pfObj%fileCount = pfObj%cFileCount + pfObj%dftFileCount
  End Subroutine loadInputFile


  Subroutine loadConfigFiles(pfObj)
! Loop through and load config files
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oPotentialFitting) :: pfObj
! Vars:  Private
    Type(oPwscfOut) :: pwscfOut
    Integer(kind=StandardInteger) :: n
! Loop through
    Do n=1,pfObj%fileCount
      If(pfObj%fileList(n)%fileType.eq.1)Then  ! config file
        Call loadConfig(pfObj, n)
      End If
      If(pfObj%fileList(n)%fileType.eq.2)Then  ! dft qe config file
        Call readPwscfOutput(trim(pfObj%fileList(n)%filePath), pwscfOut)
        !Call printPwscfOutput(pwscfOut)
        Call pfPwscfToConfig(n,pfObj,pwscfOut)
      End If
    End Do
  End Subroutine loadConfigFiles


  Subroutine loadConfig(pfObj, configKey)
! Load from a config file
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oPotentialFitting) :: pfObj
    Integer(kind=StandardInteger) :: configKey
! Vars:  Private
    Character(Len=128), Dimension(1:1024) :: fileArray
    Integer(kind=StandardInteger) :: n, coordCount, coordKey, fileLength
    Character(Len=128) :: fileRow
    Character(Len=64), Dimension(1:10) :: fieldArray
    Integer(kind=StandardInteger) :: fieldCount
    Character(Len=16) :: label, unitStr
    Real(kind=DoubleReal) :: aLat, energy
!--------------------------------------------------------
! Read potential input file
    Call readFile(pfObj%fileList(configKey)%filePath, fileArray, fileLength)
!--------------------------------------------------------
! Count number of coords
    n = 0
    coordCount = 0
    Do While(n.le.fileLength)
      n = n + 1
      fileRow = StrToUpper(fileArray(n))
      If(fileRow(1:1).ne."#")Then
        Call readFieldsCharacter(fileRow,fieldArray,fieldCount)
        If((fieldCount.eq.4).or.(fieldCount.eq.7))Then
          coordCount = coordCount + 1
        End If
      End If
    End Do
    pfObj%configs(configKey)%coordLength = coordCount
!--------------------------------------------------------
! Allocate arrays and initialise
    Call initPFCoords(pfObj, configKey)
!--------------------------------------------------------
! Load data
    n = 0
    coordKey = 0
    Do While(n.le.fileLength)
      n = n + 1
      fileRow = StrToUpper(fileArray(n))
      If(fileRow(1:1).eq."#")Then
        If(fileRow(1:5).eq."#ALAT")Then
          Call readFieldsCharacter(fileRow,fieldArray,fieldCount)
          aLat = StrToDp(trim(adjustl(fieldArray(2))))
          unitStr = "                "
          unitStr = StrMerge(trim(adjustl(fieldArray(3))),unitStr)
          pfObj%configs(configKey)%aLat = UnitConvert(aLat, unitStr, "ANGS")
        End If
        If(fileRow(1:2).eq."#E")Then
          Call readFieldsCharacter(fileRow,fieldArray,fieldCount)
          energy = StrToDp(trim(adjustl(fieldArray(2))))
          unitStr = "                "
          unitStr = StrMerge(trim(adjustl(fieldArray(3))),unitStr)
          pfObj%configs(configKey)%energy = UnitConvert(energy, unitStr, "EV")
          pfObj%configs(configKey)%epa = pfObj%configs(configKey)%energy/pfObj%configs(configKey)%coordLength
        End If
      Else
        Call readFieldsCharacter(fileRow,fieldArray,fieldCount)
        If((fieldCount.eq.4).or.(fieldCount.eq.7))Then
          coordKey = coordKey + 1
          label = "                "
          label = StrMerge(trim(adjustl(fieldArray(1))),label)
          pfObj%configs(configKey)%coords(coordKey)%label = label
          pfObj%configs(configKey)%coords(coordKey)%xyz(1) = StrToDp(trim(adjustl(fieldArray(2))))
          pfObj%configs(configKey)%coords(coordKey)%xyz(2) = StrToDp(trim(adjustl(fieldArray(3))))
          pfObj%configs(configKey)%coords(coordKey)%xyz(3) = StrToDp(trim(adjustl(fieldArray(4))))
        End If
      End If
    End Do
  End Subroutine loadConfig


  Subroutine initPFCoords(pfObj, configKey)
! Start potential fitting
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oPotentialFitting) :: pfObj
    Integer(kind=StandardInteger) :: configKey
! Vars:  Private
    Integer(kind=StandardInteger) :: n, i
! Allocate array
    If(Allocated(pfObj%configs(configKey)%coords))Then
      Deallocate(pfObj%configs(configKey)%coords)
    End If
    Allocate(pfObj%configs(configKey)%coords(1:pfObj%configs(configKey)%coordLength))
! Zero coords/forces
    Do n=1,pfObj%configs(configKey)%coordLength
      pfObj%configs(configKey)%coords(n)%label = "                "
      Do i=1,3
        pfObj%configs(configKey)%coords(n)%xyz(i) = 0.0D0
        pfObj%configs(configKey)%coords(n)%forces(i) = 0.0D0
      End Do
    End Do
  End Subroutine initPFCoords


  Subroutine pfPwscfToConfig(n,pfObj,pwscfOut)
! Start potential fitting
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Integer(kind=StandardInteger) :: n
    Type(oPotentialFitting) :: pfObj
    Type(oPwscfOut) :: pwscfOut
! Vars: Private
    Integer(kind=StandardInteger) :: i, j
! Copy pwscf to config
    pfObj%configs(n)%sourceFile = pwscfOut%filePath
    pfObj%configs(n)%aLat = pwscfOut%aLat
    pfObj%configs(n)%energy = pwscfOut%totalEnergy
    pfObj%configs(n)%coordLength = pwscfOut%atomCount
    pfObj%configs(n)%epa = pwscfOut%totalEnergy/pwscfOut%atomCount
    Call initPFCoords(pfObj, n)
    Do i=1,pwscfOut%atomCount
      Do j=1,3
        pfObj%configs(n)%coords(i)%label = pwscfOut%atoms(i)%label
        pfObj%configs(n)%coords(i)%xyz(j) = pwscfOut%atoms(i)%xyz(j)
        pfObj%configs(n)%coords(i)%forces(j) = pwscfOut%atoms(i)%forces(j)
      End Do
    End Do
  End Subroutine pfPwscfToConfig


  Subroutine pfConfigToFile(pfObj)
! Start potential fitting
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oPotentialFitting) :: pfObj
! Vars:  Private
    Integer(kind=StandardInteger) :: i, n
! Write output file
    Open(UNIT=422,FILE=Trim(pfObj%recordFile))
! Title
    write(422,"(A25)") "Potential Fitting Configs"
    write(422,"(A50)") "--------------------------------------------------"
    write(422,"(A1)") " "
! File Count
    write(422,"(A12,I4)") "File Count: ",pfObj%fileCount
    write(422,"(A1)") " "
! File List
    write(422,"(A16)") "Input File List:"
    Do n=1,pfObj%fileCount
      write(422,"(A1,A)") " ",AdjustL(trim(pfObj%fileList(n)%filePath))
    End Do
    write(422,"(A1)") " "
! Coords
    Do n=1,pfObj%fileCount
      write(422,"(A24)") "========================"
      write(422,"(A7,I8)") "Config ",n
      write(422,"(A24)") "========================"
      write(422,"(A16,E16.10)") "aLat    (ang)   ",pfObj%configs(n)%aLat
      write(422,"(A16,E16.10)") "epa     (ev)    ",pfObj%configs(n)%epa
      write(422,"(A16,E16.10)") "energy  (ev)    ",pfObj%configs(n)%energy
      write(422,"(A16,I8)") "coord length    ",pfObj%configs(n)%coordLength
      Do i=1,pfObj%configs(n)%coordLength
        write(422,"(I4,A16,F12.7,F12.7,F12.7,F12.7,F12.7,F12.7)") &
        i,trim(pfObj%configs(n)%coords(i)%label),&
        pfObj%configs(n)%coords(i)%xyz(1),&
        pfObj%configs(n)%coords(i)%xyz(2),&
        pfObj%configs(n)%coords(i)%xyz(3),&
        pfObj%configs(n)%coords(i)%forces(1),&
        pfObj%configs(n)%coords(i)%forces(2),&
        pfObj%configs(n)%coords(i)%forces(3)
      End Do
      write(422,"(A1)") " "
    End Do
! Close file
    Close(422)
  End Subroutine pfConfigToFile


  Subroutine pfConfigCount(pfObj)
! count configs (where coords gt 0) and create a map
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oPotentialFitting) :: pfObj
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j
! Loop through
    j = 0
    Do i=1,size(pfObj%configs,1)
      If(pfObj%configs(i)%coordLength.gt.0)Then
        j = j + 1
      End If
    End Do
    pfObj%configCount = j
    If(Allocated(pfObj%configMap))Then
      Deallocate(pfObj%configMap)
    End If
    Allocate(pfObj%configMap(1:pfObj%configCount))
    j = 0
    Do i=1,size(pfObj%configs,1)
      If(pfObj%configs(i)%coordLength.gt.0)Then
        j = j + 1
        pfObj%configMap(j) = i
      End If
    End Do
  End Subroutine pfConfigCount


  Subroutine pfMakeDirs(pfObj)
! count configs (where coords gt 0) and create a map
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oPotentialFitting) :: pfObj
! Vars:  Private
! Complete path
    pfObj%outputDirectory = fillPath(pfObj%outputDirectory)  ! general.f90
! Make output directory
    Call makeDir(pfObj%outputDirectory)
  End Subroutine pfMakeDirs











!---------------------------------------------------------------------------------------
