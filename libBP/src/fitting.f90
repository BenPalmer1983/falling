Module fitting
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
  Use matrix
  Use linearAlgebra
  Use lmaM
! Force declaration of all variables
  Implicit None
! Public variables
! Make private
  Private
! Public
! --variables--!
! --functions--!
  Public :: PolyFit
! Interfaces
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

  Function PolyFit(points,order) RESULT (coefficients)
! Fits a polynomial of order to the points input
! Uses Vandermonde Matrix
    Implicit None  !Force declaration of all variables
! Declare variables
! In
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Integer(kind=StandardInteger) :: order
! Out
    Real(kind=DoubleReal), Dimension(1:(order+1)) :: coefficients
! Private
    Integer(kind=StandardInteger) :: k,col,row,exponentValue
    Real(kind=DoubleReal), Dimension(1:(order+1),1:(order+1)) :: xMatrix
    Real(kind=DoubleReal), Dimension(1:(order+1)) :: yMatrix
! Build Least Squares Fitting Vandermonde matrix
    Do row=1,(order+1)
      Do col=1,(order+1)
        exponentValue = row+col-2
        xMatrix(row,col) = 0.0D0
        Do k=1,size(points,1)
          xMatrix(row,col) = 1.0D0*xMatrix(row,col)+1.0D0*points(k,1)&
          **exponentValue
        End Do
      End Do
    End Do
    Do row=1,(order+1)
      exponentValue = row-1
      yMatrix(row) = 0.0D0
      Do k=1,size(points,1)
        yMatrix(row) = 1.0D0*yMatrix(row)+1.0D0*points(k,2)*&
        points(k,1)**exponentValue
      End Do
    End Do
! invert xMatrix
    xMatrix = InvertMatrix(xMatrix)         ! matrix.f90
! multiply inverse by y to get coefficients
    coefficients = matMul(xMatrix,yMatrix)
  End Function PolyFit






End Module fitting
