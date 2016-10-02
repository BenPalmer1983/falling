! --------------------------------------------------------------!
! Conjugate Gradients module
! conjugateGradients
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Functions/Subroutines for a number of Conjugate Gradient based optimisation algorithms
!
! cgSolve Ax = b
!
! ----------------------------------------
! Updated: 4th November 2015
! ----------------------------------------

Module conjugateGradients
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
  Use matrix
  Use linearAlgebra
! Force declaration of all variables
  Implicit None
! Public variables
! Make private
  Private
! Public
! --variables--!
! --functions--!
  Public :: cgLinear
! Interfaces



!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------
! Use Newton method to solve multivariable f(x) = 0
!-----------------------------------------------------


  Function cgLinear(A, b, x0_in, convThresh_In) RESULT (x)
! Solve Ax = b using conjugate gradients method
! A must be symmertic, positive definite
    Implicit None  !Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal), Dimension(:,:) :: A
    Real(kind=DoubleReal), Dimension(:) :: b
    Real(kind=DoubleReal), Optional, Dimension(1:size(b,1)) :: x0_in
    Real(kind=DoubleReal), Optional :: convThresh_In
! Vars:  Out
    Real(kind=DoubleReal), Dimension(1:size(b,1)) :: x
! Vars:  Private
    Integer(kind=StandardInteger) :: i, n
    Real(kind=DoubleReal), Dimension(1:size(b,1)) :: xk
    Real(kind=DoubleReal), Dimension(1:size(b,1)) :: Ax
    Real(kind=DoubleReal), Dimension(1:size(b,1)) :: R    ! residuals
    Real(kind=DoubleReal), Dimension(1:size(b,1)) :: P    ! direction vector
    Real(kind=DoubleReal), Dimension(1:size(b,1)) :: AP    !
    Real(kind=DoubleReal) :: ptr, ptAp, aptr
    Real(kind=DoubleReal) :: alpha, beta, rSq
    Real(kind=DoubleReal) :: convThresh
    Logical :: iterate
! Optional arguments
    convThresh = 1.0D-8
    If(Present(convThresh_In))Then
      convThresh = convThresh_In
    End If
    xk = 0.0D0
    If(Present(convThresh_In))Then
      xk = x0_in
    End If
! Calc Ax
    Ax = matmul(A,xk)
! set up residual and direction vectors
    Do i=1,size(x,1)
      R = b - Ax
      P = R
    End Do
!Run loop
    iterate = .true.
    n = 0
    Do while (iterate)
      n = n + 1
! Prepare matrices
      AP = Matmul(A,P)
! alpha
      ptr = MatMult1D(P,R)
      ptAp = MatMult1D(P,AP)
      alpha = ptr/ptAp
! Update x
      xk = xk + alpha*P
! recalc residuals
      R = R - alpha * AP
! Breakout
      rSq = 1.0D0*MatMult1D(R,R)
      If(rSq.le.convThresh)Then
        iterate = .false.
      End If
      If(n.ge.(2*size(x,1)))Then
        ! Stop infinite loop if non-symmetric non-positive definite A used
        iterate = .false.
      End If
! calc beta
      aptr = MatMult1D(AP,R)
      beta = aptr/ptAp
! Update direction vector
      P = R - beta * P
    End Do
    x = xk  ! set x to final iteration
  End Function cgLinear







  Function cgLinearBroken(A, b, x0_in) RESULT (x)
! Solve Ax = b using conjugate gradients method
! A must be symmertic, positive definite
    Implicit None  !Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal), Dimension(:,:) :: A
    Real(kind=DoubleReal), Dimension(:) :: b
    Real(kind=DoubleReal), Dimension(1:size(b,1)) :: x0_in
! Vars:  Out
    Real(kind=DoubleReal), Dimension(1:size(b,1)) :: x
! Vars:  Private
    Integer(kind=StandardInteger) :: i, n
    Real(kind=DoubleReal), Dimension(1:size(b,1)) :: xk
    Real(kind=DoubleReal), Dimension(1:size(b,1)) :: Ax
    Real(kind=DoubleReal), Dimension(1:size(b,1)) :: R    ! residuals
    Real(kind=DoubleReal), Dimension(1:size(b,1)) :: P    ! direction vector
    Real(kind=DoubleReal), Dimension(1:size(b,1)) :: AP    !
    Real(kind=DoubleReal) :: ptAp
    Real(kind=DoubleReal) :: alpha, beta, r_sq, r_sq_new
    Logical :: iterate
! Set 0
    xk = x0_in
    Ax = matmul(A,xk)
! set up residual and direction vectors
    Do i=1,size(x,1)
      R = b - Ax
      P = R
    End Do
    beta = 0.0D0
!Run loop
    iterate = .true.
    n = 0
    Do while (iterate)
      n = n + 1
! Residual square
      r_sq = MatMult1D(R,R)
! Search Direction
      P = R - beta * P
! Prepare matrices
      AP = Matmul(A,P)
      ptAp = MatMult1D(P,AP)
! alpha - Line search
      alpha = r_sq/ptAp
! Update x
      xk = xk + alpha*P
! recalc residuals
      R = R - alpha * AP
! Breakout
      If(n.eq.(size(b,1)+1))Then
        iterate = .false.
      End If
! calc beta
      r_sq_new = MatMult1D(R,R)
      beta = r_sq_new/r_sq

    End Do
    x = xk  ! set x to final iteration
  End Function cgLinearBroken




End Module conjugateGradients






!
