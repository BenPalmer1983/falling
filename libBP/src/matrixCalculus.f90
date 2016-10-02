! --------------------------------------------------------------!
! Matrix Calculus Module
! matrixCalculusTypes, matrixCalculus
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Module Description
!
! ----------------------------------------
! Updated:
! ----------------------------------------

Module matrixCalculusTypes
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
!  Public ::

!  Type :: oType
!    Character(Len=64) :: c
!    Integer(kind=StandardInteger) :: i
!    Logical :: l
!    Real(kind=DoubleReal) :: d
!  End Type oType


End Module matrixCalculusTypes


Module matrixCalculus
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use mpi
  Use kinds
! Force declaration of all variables
  Implicit None
  Private
! ---- Variables
!  Public :: nl
! ---- Subroutines
  Public :: matrixGH
! ---- Functions
  Public :: matrixGradient
  Public :: matrixHessian
  Public :: MC_SampleA
  Public :: MC_SampleB
! Interfaces
!  Interface iName
!    Module Procedure subA, subB
!  End Interface iName


!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------



! -----------------------------------------------
!        Module Subroutines + Functions
!
! -----------------------------------------------


  Function matrixGradient(MGrad_Call,xVec,pVec) Result (gradVec)
! Function description
    Implicit None   ! Force declaration of all variables
! Vars:  In - Function call interface
    Interface
      Function MGrad_Call(xVec, pVec)
        Import DoubleReal
        Real(kind=DoubleReal), Dimension(:) :: xVec ! x vector
        Real(kind=DoubleReal), Dimension(:) :: pVec
        Real(kind=DoubleReal) :: MGrad_Call ! result
      End function MGrad_Call
    End Interface
! Vars:  In
    Real(kind=DoubleReal), Dimension(:) :: xVec   ! vector being varied
    Real(kind=DoubleReal), Dimension(:) :: pVec   ! fixed parameters if required by the function e.g. A B C for the quadratic equation
! Vars:  Out
    Real(kind=DoubleReal), Dimension(1:size(xVec,1)) :: gradVec
! Vars:  Private
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal) :: fx, fxp
    Real(kind=DoubleReal) :: h
    Real(kind=DoubleReal), Dimension(1:size(xVec,1)) :: xVecP
! Init
    gradVec = 0.0D0
    h = 1.0D-7
! Calc unperturbed
    fx = MGrad_Call(xVec, pVec)
! Use finite difference
    Do i=1,size(xVec,1)
      xVecP = xVec
      xVecP(i) = xVecP(i) + h
      fxp = MGrad_Call(xVecP, pVec)
      gradVec(i) = (fxp-fx)/h
    End Do
  End Function matrixGradient




  Function matrixHessian(MCFunc_Call,xVec,pVec) Result (hessian)
! Function description
    Implicit None   ! Force declaration of all variables
! Vars:  In - Function call interface
    Interface
      Function MCFunc_Call(xVec, pVec)
        Import DoubleReal
        Real(kind=DoubleReal), Dimension(:) :: xVec ! x vector
        Real(kind=DoubleReal), Dimension(:) :: pVec
        Real(kind=DoubleReal) :: MCFunc_Call ! result
      End function MCFunc_Call
    End Interface
! Vars:  In
    Real(kind=DoubleReal), Dimension(:) :: xVec   ! vector being varied
    Real(kind=DoubleReal), Dimension(:) :: pVec   ! fixed parameters if required by the function e.g. A B C for the quadratic equation
! Vars:  Out
    Real(kind=DoubleReal), Dimension(1:size(xVec,1),1:size(xVec,1)) :: hessian
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j
    Real(kind=DoubleReal) :: fx, fx_f, fx_b, fx_ff, fx_fb, fx_bf, fx_bb
    Real(kind=DoubleReal) :: h, hsq
    Real(kind=DoubleReal), Dimension(1:size(xVec,1)) :: xVecP
! Check matrix
    If(size(xVec,1).gt.1)Then
! Init
      hessian = 0.0D0
      h = 1.0D-5
      hsq = 1.0D-10
! Calc unperturbed
      fx = MCFunc_Call(xVec, pVec)
! Use finite difference - diagonal of hessian
      Do i=1,size(xVec,1)
        xVecP = xVec
        xVecP(i) = xVecP(i) - h
        fx_b = MCFunc_Call(xVecP, pVec)
        xVecP(i) = xVecP(i) + h + h
        fx_f = MCFunc_Call(xVecP, pVec)
! Store
        hessian(i,i) = (fx_f+fx_b-fx-fx)/hsq
      End Do
! Upper/lower
      Do i=1,size(xVec,1)-1
        Do j=i,size(xVec,1)
          xVecP = xVec
! FF
          xVecP(i) = xVecP(i) + h
          xVecP(j) = xVecP(j) + h
          fx_ff = MCFunc_Call(xVecP, pVec)
! FB
          xVecP(j) = xVecP(j) - h - h
          fx_fb = MCFunc_Call(xVecP, pVec)
! BB
          xVecP(i) = xVecP(i) - h - h
          fx_bb = MCFunc_Call(xVecP, pVec)
! BB
          xVecP(j) = xVecP(j) + h + h
          fx_bf = MCFunc_Call(xVecP, pVec)
! Store
          hessian(i,j) = (fx_ff+fx_bb-fx_fb-fx_bf)/(4.0D0*hsq)
          hessian(j,i) = hessian(i,j)
        End Do
      End Do
    End If
  End Function matrixHessian


  Subroutine matrixGH(MCFunc_Call,xVec,pVec,gradVec,hessian)
! Function description
    Implicit None   ! Force declaration of all variables
! Vars:  In - Function call interface
    Interface
      Function MCFunc_Call(xVec, pVec)
        Import DoubleReal
        Real(kind=DoubleReal), Dimension(:) :: xVec ! x vector
        Real(kind=DoubleReal), Dimension(:) :: pVec
        Real(kind=DoubleReal) :: MCFunc_Call ! result
      End function MCFunc_Call
    End Interface
! Vars:  In
    Real(kind=DoubleReal), Dimension(:) :: xVec   ! vector being varied
    Real(kind=DoubleReal), Dimension(:) :: pVec   ! fixed parameters if required by the function e.g. A B C for the quadratic equation
! Vars:  Out
    Real(kind=DoubleReal), Dimension(1:size(xVec,1)) :: gradVec
    Real(kind=DoubleReal), Dimension(1:size(xVec,1),1:size(xVec,1)) :: hessian
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j
    Real(kind=DoubleReal) :: fx, fx_f, fx_b, fx_ff, fx_fb, fx_bf, fx_bb
    Real(kind=DoubleReal) :: h, hsq
    Real(kind=DoubleReal), Dimension(1:size(xVec,1)) :: xVecP
! Check matrix
    If(size(xVec,1).gt.1)Then
! Init
      gradVec = 0.0D0
      hessian = 0.0D0
      h = 1.0D-5
      hsq = 1.0D-10
! Calc unperturbed
      fx = MCFunc_Call(xVec, pVec)
! Use finite difference - diagonal of hessian
      Do i=1,size(xVec,1)
        xVecP = xVec
        xVecP(i) = xVecP(i) - h
        fx_b = MCFunc_Call(xVecP, pVec)
        xVecP(i) = xVecP(i) + h + h
        fx_f = MCFunc_Call(xVecP, pVec)
! Store
        gradVec(i) = (fx_f-fx_b)/(h+h)
        hessian(i,i) = (fx_f+fx_b-fx-fx)/hsq
      End Do
! Upper/lower
      Do i=1,size(xVec,1)-1
        Do j=i,size(xVec,1)
          xVecP = xVec
! FF
          xVecP(i) = xVecP(i) + h
          xVecP(j) = xVecP(j) + h
          fx_ff = MCFunc_Call(xVecP, pVec)
! FB
          xVecP(j) = xVecP(j) - h - h
          fx_fb = MCFunc_Call(xVecP, pVec)
! BB
          xVecP(i) = xVecP(i) - h - h
          fx_bb = MCFunc_Call(xVecP, pVec)
! BB
          xVecP(j) = xVecP(j) + h + h
          fx_bf = MCFunc_Call(xVecP, pVec)
! Store
          hessian(i,j) = (fx_ff+fx_bb-fx_fb-fx_bf)/(4.0D0*hsq)
          hessian(j,i) = hessian(i,j)
        End Do
      End Do
    End If
  End Subroutine matrixGH




  Function MC_SampleA(xVec, pVec) RESULT (y)
! Double Quadratic function
! Ax^2+Bx+C+Ay^2+By+C
! P(1) = A, P(2) = B, P(3) = C
    Implicit None  !Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal), Dimension(:) :: xVec
    Real(kind=DoubleReal), Dimension(:) :: pVec
! Vars:  Out
    Real(kind=DoubleReal) :: y
! Evaluate function
    y = pVec(1)*xVec(1)**2+pVec(2)*xVec(1)+pVec(3)+pVec(1)*xVec(2)**2+pVec(2)*xVec(2)+pVec(3)
  End Function MC_SampleA


  Function MC_SampleB(xVec, pVec) RESULT (y)
! Double Quadratic function
! Ax^2+Bx+C+Ay^2+By+C
! P(1) = A, P(2) = B, P(3) = C
    Implicit None  !Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal), Dimension(:) :: xVec
    Real(kind=DoubleReal), Dimension(:) :: pVec
! Vars:  Out
    Real(kind=DoubleReal) :: y
! Evaluate function
    pVec = 0.0D0
    y = 5*(xVec(1)+2)**4*(xVec(2)-1)**3+3*(xVec(1)-3)**3*(0.1D0*xVec(2)+6)**4
  End Function MC_SampleB




End Module matrixCalculus

























!-----------------------------------
