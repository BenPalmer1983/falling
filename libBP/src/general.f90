! --------------------------------------------------------------!
! General module
! generalTypes, general
! collection of general functions
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
!
!
! ----------------------------------------
! Updated: 26th June 2016
! ----------------------------------------


Module generalTypes
! Setup Modules
  Use kinds
! Force declaration of all variables
  Implicit None
! Vars:  Module Parameters
!  Integer(kind=StandardInteger), Parameter :: p_ = 3
! Make private
  Private
! Public Variables and Parameters
!  Public :: p_
! Public derived types
  Public :: oFile
!------------------------------------------------------------------------------
! Defined Types
!------------------------------------------------------------------------------
  Type :: oFile
    Logical :: exists
    Character(Len=128) :: input
    Character(Len=128) :: path  ! dir + name
    Character(Len=128) :: name
    Character(Len=128) :: dir
    Character(Len=8) ::   ext
    Integer(kind=StandardInteger) :: dirLen
    Integer(kind=StandardInteger) :: nameLen
    Integer(kind=StandardInteger) :: extLen
  End Type oFile
End Module generalTypes


Module general
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use strings
  Use generalTypes
! Force declaration of all variables
  Implicit None
! Public variables
! Make private
  Private
! Public
! ---- Variables
! ---- Functions
  Public :: GetClockTime
  Public :: dpToString
  Public :: intToString
  Public :: fileInfo
  Public :: fillPath
! ---- Subroutines
  Public :: makeDir
  Public :: rmFile
  Public :: rmDir
  Public :: correctFilePath
  Public :: printFileInfo
  Public :: readFile
  Public :: readFileRaw
  Public :: prepFileData
  Public :: readCSV
  Public :: readFieldsCharacter
  Public :: extractArrayColumnDP
  Public :: extractArrayColumnInt
  Public :: swapArrayRows1D
  Public :: swapArrayRows2D
  Public :: strToIntArr
  Public :: strToDPArr
  Public :: strToStrArr
  Public :: timeAcc
  Public :: completePath
! Functions
  Public :: FileExists
  Public :: CountRowFields

! ---- Subroutines
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------


! ---------------------------------------------------------
! MODULE FUNCTIONS
! ---------------------------------------------------------

! TIMES

  Function GetClockTime () RESULT (outputTime)
! Force declaration of all variables
    Implicit None
    Real(kind=DoubleReal) :: outputTime
    Call cpu_time(outputTime)
  End Function GetClockTime


! TYPE CONVERSION

  Function dpToString(inputDP) RESULT (outputString)
! Force declaration of all variables
    Implicit None
! declare private variables
    Real(kind=DoubleReal) :: inputDP
    Character(len=32) :: outputString
! Read dp to string
    inputDP = 1.0D0 * inputDP
    Write(outputString,"(ES16.8E3)") inputDP
  End Function dpToString

  Function intToString(inputInt) RESULT (outputString)
! force declaration of all variables
    Implicit None
! declare private variables
    Integer(kind=StandardInteger) :: inputInt
    Character(len=32) :: outputString
! Read int to string
    Write(outputString,"(I16)") inputInt
  End Function intToString



! FILE INFO
  Function fileInfo(fileInput) Result (fileObj)
! Create file info object
    Implicit None  !Force declaration of all variables
! Vars:  In
    Character(*) :: fileInput
    Character(Len=128) :: fileInput_w
! Vars:  Out
    Type(oFile) :: fileObj
! Vars:  Private
    Character(Len=128) :: cwd
    Integer(kind=StandardInteger) :: i, j, k, n, levels
    Logical :: levelsFlag
! Store input
    fileObj%input = fileInput
    fileInput_w = TrimStr(fileInput)
! Strip out weird characters that are not allowed
! Check if relative or absolute
    levels = 0
    If(fileInput_w(1:1).eq."/")Then ! Absolute
      fileObj%path = fileInput_w
    Else
      Call GetCWD(cwd)
      ! Check levels "../../" etc
      levelsFlag = .true.
      Do While(levelsFlag)
        If(fileInput_w(1:3).eq."../")Then
          fileInput_w = StrRmChunk(fileInput_w,1,3," ")
          levels = levels + 1
        Else
          levelsFlag = .false.
        End If
      End Do
      Do i=1,levels
        Do j=1,len(cwd)
          k = len(cwd)+1-j
          If(cwd(k:k).eq."/")Then
            Exit
          End If
        End Do
        Do j=k,len(cwd)
          cwd(j:j) = " "
        End Do
      End Do
      fileObj%path = trim(cwd)//"/"//fileInput_w
    End If
! split up directory and file name
    Do i=1,len(fileObj%path)
      k = len(fileObj%path)+1-i
      If(fileObj%path(k:k).eq."/")Then
        Exit
      End If
    End Do
    fileObj%dir = BlankString(fileObj%dir)
    fileObj%name = BlankString(fileObj%name)
    Do i=1,k-1
      fileObj%dir(i:i) = fileObj%path(i:i)
    End Do
    j = 0
    Do i=k+1,len(fileObj%path)
      If((fileObj%path(i:i).eq." ").or.(fileObj%path(i:i).eq.Char(0)))Then
        exit
      End If
      j = j + 1
      fileObj%name(j:j) = fileObj%path(i:i)
    End Do
    fileObj%dirLen = k-1
    fileObj%nameLen = j
    fileObj%extLen = 0
! Look for extension
    fileObj%ext = "        "
    Do i=1,fileObj%nameLen
      k = fileObj%nameLen+1-i
      If(fileObj%name(k:k).eq.".")Then
        fileObj%extLen = fileObj%nameLen-k
        n = 0
        Do j=k+1,fileObj%nameLen
          n = n + 1
          If(n.le.len(fileObj%ext))Then
            fileObj%ext(n:n) = fileObj%name(j:j)
          End If
        End Do
        Exit
      End If
    End Do
! Check if file exists
    Inquire(File=trim(fileObj%path),Exist=fileObj%exists)
  End Function fileInfo

  Function fillPath(pathIn) Result (pathOut)
! Complete the full file path to a file
    Implicit None  !Force declaration of all variables
! Vars:  In
    Character(*) :: pathIn
! Vars:  Out
    Character(Len=128) :: pathOut
! Vars:  Private
    Character(Len=128) :: pathTemp
    Character(Len=128) :: cwd
    Integer(kind=StandardInteger) :: i, j, k, levels
    Logical :: levelsFlag
! Blank strings
    pathTemp = BlankString(pathTemp)
    pathOut = BlankString(pathOut)
    cwd = BlankString(cwd)
! Load input
    Call StrIntoStr(pathIn, pathTemp)
    pathTemp = Adjustl(pathTemp)
! Check if relative or absolute
    levels = 0
    If(pathTemp(1:1).eq."/")Then ! Absolute
      pathOut = trim(pathTemp)
    Else
      Call GetCWD(cwd)
! Check levels "../../" etc
      levelsFlag = .true.
      Do While(levelsFlag)
        If(pathTemp(1:3).eq."../")Then
          pathTemp = StrRmChunk(pathTemp,1,3," ")
          levels = levels + 1
        Else
          levelsFlag = .false.
        End If
      End Do
      Do i=1,levels
        Do j=1,len(cwd)
          k = len(cwd)+1-j
          If(cwd(k:k).eq."/")Then
            Exit
          End If
        End Do
        Do j=k,len(cwd)
          cwd(j:j) = " "
        End Do
      End Do
      pathOut = trim(cwd)//"/"//pathTemp
    End If
  End Function fillPath




!---------------------------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------------------



! ---------------------------------------------------------
! MODULE SUBROUTINES
! ---------------------------------------------------------

! FILES AND DIRECTORIES

  Subroutine makeDir(path)
! Make dir if it doesn't exist
    Implicit None ! Force declaration of all variables
! Vars:  In/Out
    Character(*) :: path
    Logical :: dirExists
! Check if exists
    Inquire(file=path, exist=dirExists)
    If(.not.dirExists)Then
      Call system("mkdir -p "//trim(path))
      print *,"made dir"
    End If
  End Subroutine makeDir

  Subroutine rmFile(path)
! Swap rows of square dp matrix
! force declaration of all variables
    Implicit None
! declare private variables
    Character(*) :: path
    Call system("rm -f "//trim(path))
  End Subroutine rmFile

  Subroutine rmDir(path)
! Swap rows of square dp matrix
! force declaration of all variables
    Implicit None
! declare private variables
    Character(*) :: path
    Call system("rm -fR "//trim(path))
  End Subroutine rmDir


  Subroutine correctFilePath (filePath)
! Correct file path
    Implicit None  !Force declaration of all variables
! Declare variables
    Character(*), Intent(INOUT) :: filePath
    Character(Len(filePath)) :: tempFilePath
    Integer(kind=StandardInteger) :: i, n
! Blank string
    tempFilePath = BlankString(tempFilePath)
    n = 0
    Do i=1,Len(filePath)
      If((iachar(filePath(i:i)).ge.32.and.iachar(filePath(i:i)).le.126))Then
        If(filePath(i:i).ne."?".and.filePath(i:i).ne."*".and.&
          filePath(i:i).ne."%".and.filePath(i:i).ne."+")Then
          n = n + 1
          tempFilePath(n:n) = filePath(i:i)
        End If
      End If
    End Do
    filePath = tempFilePath
  End Subroutine correctFilePath

  Subroutine printFileInfo(fileIn)
! Print file info
    Implicit None ! Force declaration of all variables
! Vars:  In
    Type(oFile) :: fileIn
! Print
    print *,trim(fileIn%path)
    If(fileIn%exists)Then
      print *,"Exists"
    Else
      print *,"Does not Exist"
    End If
    print *,"Directory:   ",trim(fileIn%dir)," (",fileIn%dirLen,")"
    print *,"File name:   ",trim(fileIn%name)," (",fileIn%nameLen,")"
    print *,"Extension:   ",trim(fileIn%ext)," (",fileIn%extLen,")"
    print *,""
  End Subroutine printFileInfo

! READING FILES

  Subroutine readFile(inputFilePath, fileArray, n, commentsOff_In, blankLinesOff_In)
! Subroutine to read file into an array
! Removes comments !.....
! Removes blank lines
! Removes leading spaces
    Implicit None ! Force declaration of all variables
! Vars:  In/Out
    Character(*) :: inputFilePath
    Character(*), Dimension(:) :: fileArray
    Integer(kind=StandardInteger) :: n
    Logical, Optional :: commentsOff_In
    Logical, Optional :: blankLinesOff_In
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, ios
    Character(len=255) :: fileRow, fileRowTemp
    Logical :: commentsFlag, commentsOff, blankLinesOff
! Optional arguments
    commentsOff = .true.    ! strip comments by default
    blankLinesOff = .true.  ! strip blank lines by default
    If(Present(commentsOff_In))Then
      commentsOff = commentsOff_In
    End If
    If(Present(blankLinesOff_In))Then
      blankLinesOff = blankLinesOff_In
    End If
! Exit if file does not exist
    If(.not.FileExists(inputFilePath))Then
      print *,"File Does Not Exist"
      stop
    End If
! open file
    Open(UNIT=9999,FILE=trim(inputFilePath),status='old',action='read')
    n = 0
    Do i=1,size(fileArray,1)
      Read(9999,"(A255)",IOSTAT=ios) fileRow
      If(ios /= 0)Then
        EXIT
      End If
! remove comments
      !If(commentsOff)Then
        commentsFlag = .false.
        Do j=1,255
          If(fileRow(j:j).eq."!")Then
            commentsFlag = .true.
          End If
          If(commentsFlag)Then
            fileRow(j:j) = char(0)
          End If
        End Do
      !End If
! remove blank lines
      !If(blankLinesOff)Then
        fileRowTemp = trim(adjustl(fileRow))
        If(fileRowTemp(1:1).ne." ")Then
          n = n + 1
          fileArray(n) = fileRowTemp
        End If
      !End If
    End Do
! Close file
    close(9999)
  End Subroutine readFile

  Subroutine readFileRaw(inputFilePath, fileArray, n)
! Subroutine to read file into an array
! Leaves comments, leading spaces, blank rows
    Implicit None ! Force declaration of all variables
! Vars:  In/Out
    Character(*) :: inputFilePath
    Character(*), Dimension(:) :: fileArray
    Integer(kind=StandardInteger) :: n
! Vars:  Private
    Integer(kind=StandardInteger) :: i, ios
    Character(len=255) :: fileRow
! Exit if file does not exist
    If(.not.FileExists(inputFilePath))Then
      print *,"File Does Not Exist"
      stop
    End If
! open file
    Open(UNIT=9999,FILE=trim(inputFilePath),status='old',action='read')
    n = 0
    Do i=1,size(fileArray,1)
      Read(9999,"(A255)",IOSTAT=ios) fileRow
      If(ios /= 0)Then
        EXIT
      End If
      n = n + 1
      fileArray(n) = fileRow
    End Do
! Close file
    close(9999)
  End Subroutine readFileRaw



  Subroutine prepFileData(inputArray, inRows, commentChars, outputArray, outRows)
    Implicit None ! Force declaration of all variables
! Vars:  In/Out
    Character(*), Dimension(:) :: inputArray
    Integer(kind=StandardInteger) :: inRows
    Character(*), Dimension(:) :: commentChars
    Character(*), Dimension(:) :: outputArray
    Integer(kind=StandardInteger) :: outRows
! Vars:  Private
    Integer(kind=StandardInteger) :: i, ii, j, n, nn, k, commentArrLen
    Character(Len(inputArray)) :: tempLine
    Integer(kind=StandardInteger), Dimension(1:size(commentChars,1)) :: commentLength
    Logical :: loopBreak, storeRow
! Count comment lengths
    Do j=1,size(commentChars,1)
      commentLength(j) = len(trim(commentChars(j)))
      If(IsBlank(commentChars(j)))Then
        commentArrLen = j - 1
        Exit
      End If
    End Do
    nn = 0
    Do n=1,inRows ! Loop through rows
      tempLine = BlankString(tempLine)
      ii = 0
      loopBreak = .false.
      storeRow = .false.
      Do i=1,Len(inputArray(n)) ! Loop through characters in row
        If((.not.storeRow).and.(inputArray(n)(i:i).ne." "))Then
          storeRow = .true.
        End If
        If(storeRow)Then
          Do j=1,commentArrLen
            If(.not.IsBlank(commentChars(j)))Then
              k = i-1+commentLength(j)
              If(k.le.Len(inputArray(n)))Then
                If(inputArray(n)(i:i+commentLength(j)-1).eq.trim(commentChars(j)))Then
                  loopBreak = .true.
                End If
              End If
            End If
          End Do
          If(loopBreak)Then
            Exit
          End If
          ii = ii + 1
          tempLine(ii:ii) = inputArray(n)(i:i)
        End If
      End Do
      If(.not.IsBlank(tempLine))Then
        nn = nn + 1
        outputArray(nn) = tempLine
      End If
    End Do
    outRows = nn
  End Subroutine prepFileData


  Subroutine readCSV(inputFilePath, fieldSeparator, csvArray, rows, columns)
! Subroutine to read csv file into an array
! Removes comments !.....
! Removes blank lines
! Removes leading spaces
    Implicit None ! Force declaration of all variables
! Vars: In/Out
    Character(*) :: inputFilePath
    Character(*) :: fieldSeparator
    Real(kind=DoubleReal), Dimension(:,:) :: csvArray
    Integer(kind=StandardInteger) :: rows, columns
! Vars: Private
    Integer(kind=StandardInteger) :: i, j, k, n, ios
    Integer(kind=StandardInteger) :: row, col
    Character(len=255), Dimension(1:32000) :: fileArray
    Character(len=255) :: fileRow, fileRowTemp
    Character(len=255) :: fieldTemp
    Logical :: commentsFlag
!
! Set vars
    rows = 0
    columns = 0
!
! Part 1: read file into array, strip out any comments and blank lines
! open file
    Open(UNIT=9999,FILE=trim(inputFilePath),status='old',action='read')
    n = 0
    Do i=1,size(fileArray,1)
      Read(9999,"(A255)",IOSTAT=ios) fileRow
    If(ios /= 0)Then
        EXIT
      End If
! remove comments
      commentsFlag = .false.
      Do j=1,255
        If(fileRow(j:j).eq."!")Then
          commentsFlag = .true.
        End If
        If(commentsFlag)Then
          fileRow(j:j) = " "
        End If
      End Do
! remove blank lines
      fileRowTemp = trim(adjustl(fileRow))
      If(fileRowTemp(1:1).ne." ")Then
        n = n + 1
        fileArray(n) = fileRowTemp
      End If
    End Do
! Close file
    close(9999)
!
! Part 2: read into array
    Do i=1,n    ! Loop through file rows
      fileRow = fileArray(i)
      row = i
      col = 0
      fieldTemp = BlankString(fieldTemp)
      k = 0
      Do j=1,255
        If(fileRow(j:j).eq.fieldSeparator.or.j.eq.255)Then
          col = col + 1
          fieldTemp = trim(fieldTemp)
          Read(fieldTemp,*) csvArray(row,col)
          fieldTemp = BlankString(fieldTemp)
          k = 0
        Else
          k = k + 1
          fieldTemp(k:k) = fileRow(j:j)
        End if
      End Do
      If(col.gt.columns)Then
        columns = col
      End If
    End Do
    rows = n
  End Subroutine readCSV

  Subroutine readFieldsCharacter(inputRow,outputArray,fieldCount,splitCharIn)
! Subroutine to read input row into an array [CHARACTER]
    Implicit None ! Force declaration of all variables
! Vars:  In/Out
    Character(*) :: inputRow
    Character(*), Dimension(:) :: outputArray
    Integer(kind=StandardInteger) :: fieldCount
    Character(Len=1), Optional :: splitCharIn
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, n, charNum, charNum_Temp, endChar, splitCharNum
    Character(Len(inputRow)) :: inputRowTemp
    Character(Len(outputArray)) :: fieldTemp
    Logical :: inQuotes
    Character(Len=1) :: splitChar
! Optional Arguments
    splitChar = " "
    If(Present(splitCharIn))Then
      splitChar = splitCharIn
    End If
    splitCharNum = iChar(splitChar)
! Init
    outputArray = WipeStringArray(outputArray)
    inputRowTemp = WipeString(inputRowTemp)
    inQuotes = .false.
! Filter out multiple spaces
    j = 0
    Do i=1,Len(inputRow)
      charNum = iChar(inputRow(i:i))
      If(charNum.eq.34)Then
        If(inQuotes)Then
          inQuotes = .false.
        Else
          inQuotes = .true.
        End If
      End If
      If(j.eq.0)Then
        j = 1
        inputRowTemp(j:j) = inputRow(i:i)
      Else
        If(inQuotes)Then
          j = j + 1
          inputRowTemp(j:j) = inputRow(i:i)
        Else
          charNum_Temp = iChar(inputRowTemp(j:j))
          If(charNum_Temp.eq.32)Then
            If(charNum.ne.32)Then
              j = j + 1
              inputRowTemp(j:j) = inputRow(i:i)
            End If
          Else
            j = j + 1
            inputRowTemp(j:j) = inputRow(i:i)
          End If
        End If
      End If
    End Do
    endChar = j
! Remove trailing space
    charNum = iChar(inputRowTemp(endChar:endChar))
    If(charNum.eq.32)Then
      endChar = endChar - 1
    End If
    Do i=(endChar+1),Len(inputRow)
      inputRowTemp(i:i) = char(0)
    End Do
! Store in array
    inQuotes = .false.
    j = 0
    n = 0
! Loop through characters
    fieldTemp = WipeString(fieldTemp)
    Do i=1,endChar
! Check if in quotes
      charNum = iChar(inputRowTemp(i:i))
      If(charNum.eq.34)Then
        If(inQuotes)Then
          inQuotes = .false.
        Else
          inQuotes = .true.
        End If
      End If
! Store
      If(inQuotes)Then
        j = j + 1
        If(j.le.len(fieldTemp))Then
          fieldTemp(j:j) = inputRowTemp(i:i)
        End If
      Else
        If(charNum.eq.splitCharNum)Then
          n = n + 1
          outputArray(n) = fieldTemp
          fieldTemp = WipeString(fieldTemp)
          j = 0
        Else
          j = j + 1
          If(j.le.len(fieldTemp))Then
            fieldTemp(j:j) = inputRowTemp(i:i)
          End If
        End If
      End If
    End Do
! Store final
    n = n + 1
    outputArray(n) = fieldTemp
    fieldCount = n
! ----------------------
  End Subroutine readFieldsCharacter




! ARRAYS

  Subroutine PrintMatrix(xMatrix)
! Prints out a 2D matrix
    Implicit None  ! Force declaration of all variables
! Declare private variables
    Real(kind=DoubleReal), Dimension(:,:) :: xMatrix
    Integer(kind=StandardInteger) :: i, j, k, n
    Character(len=16) :: tempStr
    Character(len=4096) :: printRow
! Print
    Do i=1,size(xMatrix,1)
      Do j=1,4096
        printRow(j:j) = " "
      End Do
      Do j=1,size(xMatrix,2)
        Do k=1,16
          tempStr(k:k) = " "
        End Do
        write(tempStr,"(E10.4)") xMatrix(i,j)
        Do k=1,16
          n = (j-1)*16+k
          printRow(n:n) = tempStr(k:k)
        End Do
      End Do
      print *,trim(printRow)
    End Do
  End Subroutine PrintMatrix

 Subroutine extractArrayColumnDP(inputArray,outputArray,column)
! Extract one column of a 2D dp array array(row,col)
! force declaration of all variables
    Implicit None
! declare private variables
    Integer(kind=StandardInteger) :: i, column
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputArray
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: outputArray
! Allocate output array
    Allocate(outputArray(1:size(inputArray,1)))
! Copy column
    Do i=1,size(inputArray,1)
      outputArray(i) = inputArray(i,column)
    End Do
  End Subroutine extractArrayColumnDP

  Subroutine extractArrayColumnInt(inputArray,outputArray)
! Subroutine extractArrayColumnInt(inputArray,outputArray,column)
! Extract one column of a 2D int array array(row,col)
! force declaration of all variables
    Implicit None
! declare private variables
    Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: inputArray
    Integer(kind=StandardInteger), Dimension( : ), Allocatable :: outputArray
! Allocate output array
    Allocate(outputArray(1:size(inputArray,1)))
! Copy column
! Do i=1,size(inputArray,1)
!  outputArray(i) = inputArray(i,column)
! End Do
  End Subroutine extractArrayColumnInt

  Subroutine swapArrayRows1D(matrix,rowA,rowB)
! Swap rows of square dp matrix
! force declaration of all variables
    Implicit None
! declare private variables
    Integer(kind=StandardInteger) :: i, rowA, rowB, matH, matW
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: matrix
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: rowAArr
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: rowBArr
! Set variables
    matH = size(matrix,1)
    matW = 1
! Only do if rows are in the matrix
    If(rowA.ge.1.and.rowA.le.matH.and.rowB.ge.1.and.rowB.le.matH)Then
! Allocate arrays
      Allocate(rowAArr(1:matW))
      Allocate(rowBArr(1:matW))
! Swap rows
      Do i=1,matW
        rowAArr(i) = matrix(rowA)
        rowBArr(i) = matrix(rowB)
      End Do
      Do i=1,matW
        matrix(rowA) = rowBArr(i)
        matrix(rowB) = rowAArr(i)
      End Do
    End If
  End Subroutine swapArrayRows1D

  Subroutine swapArrayRows2D(matrix,rowA,rowB)
! Swap rows of square dp matrix
! force declaration of all variables
    Implicit None
! declare private variables
    Integer(kind=StandardInteger) :: i, rowA, rowB, matH, matW
    Real(kind=DoubleReal), Dimension( : , :), Allocatable :: matrix
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: rowAArr
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: rowBArr
! Set variables
    matH = size(matrix,1)
    matW = size(matrix,2)
! Only do if rows are in the matrix
    If(rowA.ge.1.and.rowA.le.matH.and.rowB.ge.1.and.rowB.le.matH)Then
! Allocate arrays
      Allocate(rowAArr(1:matW))
      Allocate(rowBArr(1:matW))
! Swap rows
      Do i=1,matW
        rowAArr(i) = matrix(rowA,i)
        rowBArr(i) = matrix(rowB,i)
      End Do
      Do i=1,matW
        matrix(rowA,i) = rowBArr(i)
        matrix(rowB,i) = rowAArr(i)
      End Do
    End If
  End Subroutine swapArrayRows2D

  ! Integer, 1D:
  Subroutine swapRows_Integer_1D(matrix,rowA,rowB)
    Integer(kind=StandardInteger) :: matrix(:)
    Integer(kind=StandardInteger) :: temp,rowA,rowB
! Swap rows
    temp = matrix(rowA)
    matrix(rowA) = matrix(rowB)
    matrix(rowB) = temp
  End Subroutine swapRows_Integer_1D

! Integer, 2D:
  Subroutine swapRows_Integer_2D(matrix,rowA,rowB)
    Integer(kind=StandardInteger) :: matrix(:,:)
    Integer(kind=StandardInteger) :: temp, i, rowA, rowB
! Swap rows
    Do i=1,size(matrix,2)  !Loop through columns
      temp = matrix(rowA,i)
      matrix(rowA,i) = matrix(rowB,i)
      matrix(rowB,i) = temp
    End Do
  End Subroutine swapRows_Integer_2D

! Single, 1D:
  Subroutine swapRows_Single_1D(matrix,rowA,rowB)
    Real(kind=SingleReal) :: matrix(:)
    Integer(kind=StandardInteger) :: rowA, rowB
    Real(kind=SingleReal) :: temp
! Swap rows
    temp = matrix(rowA)
    matrix(rowA) = matrix(rowB)
    matrix(rowB) = temp
  End Subroutine swapRows_Single_1D

! Single, 2D:
  Subroutine swapRows_Single_2D(matrix,rowA,rowB)
    Real(kind=SingleReal) :: matrix(:,:)
    Integer(kind=StandardInteger) :: rowA, rowB
    Real(kind=SingleReal) :: temp
    Integer(kind=StandardInteger) :: i
! Swap rows
    Do i=1,size(matrix,2)  !Loop through columns
      temp = matrix(rowA,i)
      matrix(rowA,i) = matrix(rowB,i)
      matrix(rowB,i) = temp
    End Do
  End Subroutine swapRows_Single_2D

! Double, 1D:
  Subroutine swapRows_Double_1D(matrix,rowA,rowB)
    Real(kind=DoubleReal) :: matrix(:)
    Integer(kind=StandardInteger) :: rowA, rowB
    Real(kind=DoubleReal) :: temp
! Swap rows
    temp = matrix(rowA)
    matrix(rowA) = matrix(rowB)
    matrix(rowB) = temp
  End Subroutine swapRows_Double_1D

! Double, 2D:
  Subroutine swapRows_Double_2D(matrix,rowA,rowB)
    Real(kind=DoubleReal) :: matrix(:,:)
    Integer(kind=StandardInteger) :: rowA, rowB
    Real(kind=DoubleReal) :: temp
    Integer(kind=StandardInteger) :: i
! Swap rows
    Do i=1,size(matrix,2)  !Loop through columns
      temp = matrix(rowA,i)
      matrix(rowA,i) = matrix(rowB,i)
      matrix(rowB,i) = temp
    End Do
  End Subroutine swapRows_Double_2D

! Integer, 1D:
  Subroutine sort_Integer_1D(list)
    Integer(kind=StandardInteger) :: list(:)
    Integer(kind=StandardInteger) :: i, sortComplete
! Sort list
    sortComplete = 0
    Do While(sortComplete.eq.0)
      sortComplete = 1
      Do i=2,size(list,1)
        If(list(i-1).gt.list(i))Then
          Call swapRows_Integer_1D(list,i-1,i)
          sortComplete = 0
        End If
      End Do
    End Do
  End Subroutine sort_Integer_1D

! Integer, 2D:
  Subroutine sort_Integer_2D(list, sortRow)
    Integer(kind=StandardInteger) :: list(:,:)
    Integer(kind=StandardInteger) :: i, sortRow, sortComplete
! Sort list
    sortComplete = 0
    Do While(sortComplete.eq.0)
      sortComplete = 1
      Do i=2,size(list,1)
        If(list(i-1,sortRow).gt.list(i,sortRow))Then
          Call swapRows_Integer_2D(list,i-1,i)
          sortComplete = 0
        End If
      End Do
    End Do
  End Subroutine sort_Integer_2D


! TYPE CONVERSION

  Subroutine strToIntArr(stringIn,intArr)
! Take space separated integers and convert to array
    Implicit None   ! Force declaration of all variables
! Declare private variables
    Integer(kind=StandardInteger) :: i, j, k
    Integer(kind=StandardInteger), Dimension(:) :: intArr
    Character(*) :: stringIn
    Character(len(stringIn)) :: stringPrep
    Character(32) :: intString
! Prepare input string
    stringIn = Trim(Adjustl(stringIn))
    stringPrep = BlankString(stringPrep)
! One space only
    j = 0
    Do i=1,len(stringIn)
      If(i.eq.1)Then
        j = j + 1
        stringPrep(j:j) = stringIn(i:i)
      Else
        If(ichar(stringIn(i:i)).eq.32.and.ichar(stringIn(i-1:i-1)).eq.32)Then
! Do not add
        Else
          j = j + 1
          stringPrep(j:j) = stringIn(i:i)
        End If
      End If
    End Do
    intString = BlankString(intString)
    j = 0
    k = 0
    Do i=1,len(stringPrep)
      If(ichar(stringPrep(i:i)).eq.32)Then  !Space
        j = 0
        k = k + 1
        Read(intString,*) intArr(k)
        intString = BlankString(intString)
        If(ichar(stringPrep(i+1:i+1)).eq.32)Then
          Exit
        End If
      Else
        j = j + 1
        intString(j:j) = stringPrep(i:i)
      End If
    End Do
  End Subroutine strToIntArr

  Subroutine strToDPArr(stringIn,dpArr)
! Take space separated integers and convert to array
    Implicit None   ! Force declaration of all variables
! Declare private variables
    Integer(kind=StandardInteger) :: i, j, k
    Real(kind=DoubleReal), Dimension(:) :: dpArr
    Character(*) :: stringIn
    Character(len(stringIn)) :: stringPrep
    Character(32) :: intString
! Prepare input string
    stringIn = Trim(Adjustl(stringIn))
    stringPrep = BlankString(stringPrep)
! One space only
    j = 0
    Do i=1,len(stringIn)
      If(i.eq.1)Then
        j = j + 1
        stringPrep(j:j) = stringIn(i:i)
      Else
        If(ichar(stringIn(i:i)).eq.32.and.ichar(stringIn(i-1:i-1)).eq.32)Then
! Do not add
        Else
          j = j + 1
          stringPrep(j:j) = stringIn(i:i)
        End If
      End If
    End Do
    intString = BlankString(intString)
    j = 0
    k = 0
    Do i=1,len(stringPrep)
      If(ichar(stringPrep(i:i)).eq.32)Then  !Space
        j = 0
        k = k + 1
        Read(intString,*) dpArr(k)
        intString = BlankString(intString)
        If(ichar(stringPrep(i+1:i+1)).eq.32)Then
          Exit
        End If
      Else
        j = j + 1
        intString(j:j) = stringPrep(i:i)
      End If
    End Do
  End Subroutine strToDPArr

  Subroutine strToStrArr(stringIn,strArr)
! Take space separated string and convert to array
    Implicit None   ! Force declaration of all variables
! Declare private variables
    Integer(kind=StandardInteger) :: i, j, k
    Character(*), Dimension(:) :: strArr
    Character(*) :: stringIn
    Character(len(stringIn)) :: stringPrep
    Character(Len(strArr)) :: tempString
! Prepare input string
    stringIn = Trim(Adjustl(stringIn))
    stringPrep = BlankString(stringPrep)
! One space only
    j = 0
    Do i=1,len(stringIn)
      If(i.eq.1)Then
        j = j + 1
        stringPrep(j:j) = stringIn(i:i)
      Else
        If(ichar(stringIn(i:i)).eq.32.and.ichar(stringIn(i-1:i-1)).eq.32)Then
! Do not add
        Else
          j = j + 1
          stringPrep(j:j) = stringIn(i:i)
        End If
      End If
    End Do
    tempString = BlankString(tempString)
    j = 0
    k = 0
    Do i=1,len(stringPrep)
      If(ichar(stringPrep(i:i)).eq.32)Then  !Space
        j = 0
        k = k + 1
        strArr(k) = tempString
        tempString = BlankString(tempString)
        If(ichar(stringPrep(i+1:i+1)).eq.32)Then
          Exit
        End If
      Else
        j = j + 1
        tempString(j:j) = stringPrep(i:i)
      End If
    End Do
  End Subroutine strToStrArr

! Time Accumulator subroutine
  Subroutine timeAcc(time,timeStart,timeEnd)
! Take space separated integers and convert to array
    Implicit None   ! Force declaration of all variables
! Declare private variables
    Real(kind=DoubleReal) :: time,timeStart,timeEnd
    time = time + timeEnd - timeStart
  End Subroutine timeAcc

! Complete path subroutine
  Subroutine completePath(pathIn)
! Complete path using the current working directory
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Character(*) :: pathIn
! Vars:  Private
    Character(Len(pathIn)) :: pathCWD
! Check if absolute i.e. starts with /
    If(pathIn(1:1).eq."/")Then
! Do nothing, already absolute path
    Else
      Call getcwd(pathCWD)
      pathCWD = trim(pathCWD)//"/"
      pathIn = ConcatStr(pathCWD,pathIn)
    End If
  End Subroutine completePath



! --------------------------------------------------------------------------------------------------
!    Functions
! --------------------------------------------------------------------------------------------------


  Function FileExists(filePath) Result (boolOut)
    Implicit None   ! Force declaration of all variables
! Private variables
    Character(*) :: filePath
    Logical :: boolOut
! Inquire
    INQUIRE(FILE=trim(filePath), EXIST=boolOut)
  End Function FileExists

  Function CountRowFields(fileRow) Result (fieldCount)
    Implicit None   ! Force declaration of all variables
! In
    Character(*) :: fileRow
! Out
    Integer(kind=StandardInteger) :: fieldCount
! Private variables
    Integer(kind=StandardInteger) :: i
    Logical :: inQuotes
!    Character(Len(fileRow)) :: fileRowTemp

! Blank string
!   fileRowTemp = BlankString(fileRowTemp)
    fieldCount = 1
    inQuotes = .false.
    Do i=1,Len(trim(fileRow))
      If(fileRow(i:i).eq.'"')Then
        If(inQuotes)Then
          inQuotes = .false.
        Else
          inQuotes = .true.
        End If
      End If
      If(.Not.inQuotes)Then
        If(i.gt.1)Then
          If(fileRow(i:i).eq." ".and.fileRow(i-1:i-1).ne." ")Then
            fieldCount = fieldCount + 1
          End If
        End If
      End If
    End Do
  End Function CountRowFields






!---------------------------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------------------
End Module general
