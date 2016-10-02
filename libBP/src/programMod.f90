! --------------------------------------------------------------!
! Materials module
! materialsTypes, materials
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
!
!
! ----------------------------------------
! Updated: 26th June 2016
! ----------------------------------------


Module programModTypes
! Setup Modules
  Use kinds
  Use generalTypes
! Force declaration of all variables
  Implicit None
! Vars:  Module Parameters
!  Integer(kind=StandardInteger), Parameter :: p_ = 3
! Make private
  Private
! Public Variables and Parameters
!  Public :: p_
! Public derived types
  Public :: oProgram
!------------------------------------------------------------------------------
! Defined Types
!------------------------------------------------------------------------------

  Type :: oProgramFile
    Type(oFile) :: file
    Character(Len=128), Dimension(1:8192) :: data
    Integer(kind=StandardInteger) :: rows
    Character(Len=128), Dimension(1:8192) :: dataP
    Integer(kind=StandardInteger) :: rowsP
    Character(Len=128), Dimension(1:8192) :: dataP_UC
  End Type oProgramFile

  Type :: oProgram
    Real(kind=DoubleReal) :: startTime
    Character(Len=64) :: program
    Integer(kind=StandardInteger) :: args
    Character(Len=64), Dimension(:), Allocatable :: argList
    Integer(kind=StandardInteger) :: fileCount
    Type(oFile), Dimension(:), Allocatable :: filesList
    Type(oProgramFile), Dimension(:), Allocatable :: programFiles
  End Type oProgram
End Module programModTypes


Module programMod
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
  Use strings
  Use generalTypes
  Use general
  Use programModTypes
! Force declaration of all variables
  Implicit None
! Public variables
  Type(oProgram) :: programObj
! Make private
  Private
! Public:  Vars
  Public :: programObj
! Public:  Functions/Subroutines
  Public :: ProgramInit
! Interfaces
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! ------------------------------------------------------------------------!
! programMod
! ------------------------------------------------------------------------!

  Subroutine ProgramInit
! Run
    Implicit None ! Force declaration of all variables
! Vars:  In/Out
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, n
    Type(oFile), Dimension(:), Allocatable :: filesListTemp
    Character(Len=2), Dimension(1:5) :: commentChars
    Character(Len=128) :: tempRow
    Logical :: breakLoop
! Program Start Time
    Call Cpu_Time(programObj%startTime)
! Input Arguments
    programObj%args = command_argument_count()
! Allocate array
    Allocate(programObj%argList(1:programObj%args))
    Allocate(filesListTemp(1:programObj%args))
! Store arguments
    If(programObj%args.gt.0)Then
      Do n=1,programObj%args
        call get_command_argument(n, programObj%argList(n))
        programObj%argList(n) = trim(adjustl(programObj%argList(n)))
      End Do
    End If
! Check if arguments are valid files
    j = 0
    Do i=1,programObj%args
      filesListTemp(i) = fileInfo(programObj%argList(i))
      If(filesListTemp(i)%exists)Then
        j = j + 1
      End If
    End Do
    n = j
    programObj%fileCount = n
! Allocate
    Allocate(programObj%filesList(1:n))
    Allocate(programObj%programFiles(1:n))
! Save
    j = 0
    Do i=1,n
      If(filesListTemp(i)%exists)Then
        j = j + 1
        programObj%filesList(j) = filesListTemp(i)
      End If
    End Do
! Comment Chars
    commentChars = BlankStringArray(commentChars)
    commentChars(1) = "!"
    commentChars(2) = "//"
! Load files
    Do i=1,n
      Call readFileRaw(programObj%filesList(i)%path, &
      programObj%programFiles(i)%data, programObj%programFiles(i)%rows)
      Call prepFileData(programObj%programFiles(i)%data, programObj%programFiles(i)%rows, &
      commentChars, programObj%programFiles(i)%dataP, programObj%programFiles(i)%rowsP)
      Do j=1,programObj%programFiles(i)%rowsP
        programObj%programFiles(i)%dataP_UC(j) = StrToUpper(programObj%programFiles(i)%dataP(j))
      End Do
    End Do
! Check for progam
    breakLoop = .false.
    Do i=1,programObj%fileCount
      Do j=1,programObj%programFiles(i)%rowsP
        tempRow = StrToUpper(programObj%programFiles(i)%dataP(j))
        If(tempRow(1:9).eq."##PROGRAM")Then
          programObj%program = trim(adjustl(tempRow(10:64)))
        End If
      End Do
      If(breakLoop)Then
        Exit
      End If
    End Do
  End Subroutine ProgramInit


End Module programMod

!-----------------------------------------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------------------------------------
