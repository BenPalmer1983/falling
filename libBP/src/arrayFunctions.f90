Module arrayFunctions
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use strings
! Force declaration of all variables
  Implicit None
! Make private
  Private
! Public
! ---- Functions
! ---- Subroutines
  Public :: PrintMatrix
!  Public :: SortArray
  Public :: swapRows
! Interfaces
!  Interface SortArray
!    Module Procedure SortArray_1D_DP
!  End Interface SortArray
  Interface PrintMatrix
    Module Procedure PrintMatrix_1D, PrintMatrix_2D
  End Interface PrintMatrix
  Interface swapRows
    Module Procedure swapRows_Integer_1D, swapRows_Integer_2D, swapRows_Single_1D,&
    swapRows_Single_2D, swapRows_Double_1D, swapRows_Double_2D
  End Interface swapRows
! ---- Subroutines
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! ---------------------------------------------------------
! MODULE FUNCTIONS
! ---------------------------------------------------------

! ---------------------------------------------------------
! MODULE SUBROUTINES
! ---------------------------------------------------------

! ARRAYS

!  Subroutine SortArray_1D_DP(arrayIn)
! Sort array
!    Implicit None  ! Force declaration of all variables
! Vars:  In/Out
!    Real(kind=DoubleReal), Dimension(:) :: arrayIn
! Vars:  Private
!    Integer(kind=StandardInteger) :: i
! Sort array
!  End Subroutine SortArray_1D_DP

  Subroutine PrintMatrix_1D(xMatrix)
! Prints out a 2D matrix
    Implicit None  ! Force declaration of all variables
! Declare private variables
    Real(kind=DoubleReal), Dimension(:) :: xMatrix
    Integer(kind=StandardInteger) :: i
    Character(len=16) :: tempStr
! Print
    Do i=1,size(xMatrix,1)
      write(tempStr,"(E10.4)") xMatrix(i)
      print *,trim(tempStr)
    End Do
  End Subroutine PrintMatrix_1D

  Subroutine PrintMatrix_2D(xMatrix)
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
  End Subroutine PrintMatrix_2D

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

! ---------------
! Swap Rows
! ---------------

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


  Subroutine swapRows_Char_1D(matrix,rowA,rowB)
! Swap rows, character, 1D array:
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Character(*), Dimension(:) :: matrix
    Integer(kind=StandardInteger) :: rowA, rowB
    Character(Len(matrix)) :: temp
! Swap rows
    temp = matrix(rowA)
    matrix(rowA) = matrix(rowB)
    matrix(rowB) = temp
  End Subroutine swapRows_Char_1D


  Subroutine swapRows_Char_2D(matrix,rowA,rowB)
! Swap rows, character, 2D array:
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Character(*), Dimension(:,:) :: matrix
    Integer(kind=StandardInteger) :: rowA, rowB
! Vars:  Private
    Character(Len(matrix)) :: temp
    Integer(kind=StandardInteger) :: i
! Swap rows
    Do i=1,size(matrix,2)
      temp = matrix(rowA,i)
      matrix(rowA,i) = matrix(rowB,i)
      matrix(rowB,i) = temp
    End Do
  End Subroutine swapRows_Char_2D


















!---------------------------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------------------
End Module arrayFunctions
