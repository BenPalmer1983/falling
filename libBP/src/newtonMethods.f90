! --------------------------------------------------------------!
! Newton Methods module
! newtonMethods
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Functions/Subroutines for a number of Newton based optimisation algorithms
!
! NewtonSolve - f(x) = 0 where x is a multivariable vector
!
! ----------------------------------------
! Updated: 4th November 2015
! ----------------------------------------

Module newtonMethods
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
  Use matrix
  Use matrixCalculus
  Use linearAlgebra
! Force declaration of all variables
  Implicit None
! Public variables
! Make private
  Private
! Public
! --variables--!
! --functions--!
  Public :: NewtonSolve
  Public :: NewtonSolve_Sample
  Public :: NewtonGauss
  Public :: NewtonGauss_Quadratic
  Public :: NewtonOptimise
  Public :: NewtonOpt_SampleA
  Public :: NewtonOpt_SampleB
! Interfaces



!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------
! Use Newton method to solve multivariable f(x) = 0
!-----------------------------------------------------


  Function NewtonSolve(NewtonSolve_Call, parametersIn, convThr_In) RESULT (parametersOut)
! Uses Newton method to solve F(x) = 0 where x is a vector
    Implicit None  !Force declaration of all variables
    Interface
      Function NewtonSolve_Call(parameters)
        Import DoubleReal
        Real(kind=DoubleReal), Dimension(:), Intent(IN) :: parameters ! x vector
        Real(kind=DoubleReal), Dimension(1:size(parameters,1)) :: NewtonSolve_Call ! result
      End function NewtonSolve_Call
    End Interface
! Vars:  In
    Real(kind=DoubleReal), Dimension(:) :: parametersIn
    Real(kind=DoubleReal), Optional :: convThr_In
! Vars:  Out
    Real(kind=DoubleReal), Dimension(1:size(parametersIn,1)) :: parametersOut
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, pSize, n
    Real(kind=DoubleReal), Dimension(1:size(parametersIn,1)) :: fx, fxB, fxF
    Real(kind=DoubleReal), Dimension(1:size(parametersIn,1),1:size(parametersIn,1)) :: Jacobian
    Real(kind=DoubleReal) :: h, convThr, errA, errB, errConv
    Real(kind=DoubleReal), Dimension(1:size(parametersIn,1)) :: parametersTrial
    Real(kind=DoubleReal), Dimension(1:size(parametersIn,1)) :: dX
    Logical :: iterate
! Optional Arguments
    convThr = 1.0D-8
    If(Present(convThr_In))Then
      convThr = convThr_In
    End If
! Init vars
    pSize = size(parametersIn,1)
    h = 0.000001D0
    parametersOut = parametersIn
! Iterate
    iterate = .true.
    n = 0
    Do while (iterate)
      n = n + 1
! Calc fx
      fx = NewtonSolve_Call(parametersOut)
      errA = NewtonSolve_Err(fx)
! Calc Jacobian
      Do i =1,pSize ! loop through parameters
! Back difference
        parametersTrial(i) = parametersOut(i) - h
        fxB = NewtonSolve_Call(parametersTrial)
! Forward Difference
        parametersTrial(i) = parametersOut(i) + h
        fxF = NewtonSolve_Call(parametersTrial)
! Store central difference
        Do j=1,pSize ! Loop through functions - functions rows, parameters columns
          Jacobian(j,i) = (fxF(j)-fxB(j))/(2.0D0*h)
        End Do
      End Do
! Invert Jacobian
      Jacobian = InvertMatrix(Jacobian)
! Calculate update matrix
      dX = matmul(Jacobian,fx)
! Update variables vector
      parametersOut = parametersOut - dX
      fx = NewtonSolve_Call(parametersOut)
      errB = NewtonSolve_Err(fx)
      errConv = abs(errA - errB)
! Break out if converged
      If(errConv.lt.convThr)Then
        iterate = .false.
      End If
! Break out whatever the answer after 100 iterations
      If(n.eq.100)Then
        iterate = .false.
      End If
    End Do
  End Function NewtonSolve

  Function NewtonSolve_Err(fx) RESULT (rss)
! rss error from f(x) = 0
    Implicit None  !Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal), Dimension(:) :: fx
! Vars:  Out
    Real(kind=DoubleReal) :: rss
! Vars:  Private
    Integer(kind=StandardInteger) :: i
! Calc error
    rss = 0.0D0
    Do i=1,size(fx,1)
      rss = rss + fx(i)**2
    End Do
  End Function NewtonSolve_Err

  Function NewtonSolve_Sample(parameters) RESULT (fx)
! Sample function to solve with two functions and two unknowns
    Implicit None  !Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal), Dimension(:), Intent(IN) :: parameters ! x vector
! Vars:  Out
    Real(kind=DoubleReal), Dimension(1:size(parameters,1)) :: fx ! result
! Two functions, two variable vector
    fx(1) = 4.0D0*parameters(1)**2-parameters(2)**3+28    !  aim, fx(1) = 0
    fx(2) = 3.0D0*parameters(1)**3-parameters(2)**2-145   !  aim, fx(2) = 0
  End Function NewtonSolve_Sample


!-----------------------------------------------------
! Use Newton method to minimise least squares between data points and a function f(x)
!-----------------------------------------------------

  Function NewtonGauss(NewtonGauss_Call, parametersIn, pointsIn, convThr_In) RESULT (parametersOut)
! Uses Newton method to solve F(x) = 0 where x is a vector
    Implicit None  !Force declaration of all variables
! Vars:  In - Function call interface
    Interface
      Function NewtonGauss_Call(parameters, x)
        Import DoubleReal
        Real(kind=DoubleReal), Dimension(:), Intent(IN) :: parameters ! x vector
        Real(kind=DoubleReal) :: x
        Real(kind=DoubleReal) :: NewtonGauss_Call ! result
      End function NewtonGauss_Call
    End Interface
! Vars:  In
    Real(kind=DoubleReal), Dimension(:,:) :: pointsIn
    Real(kind=DoubleReal), Dimension(:) :: parametersIn
    Real(kind=DoubleReal), Optional :: convThr_In
! Vars:  Out
    Real(kind=DoubleReal), Dimension(1:size(parametersIn,1)) :: parametersOut
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, pSize, n, nPoints
    Real(kind=DoubleReal), Dimension(1:size(pointsIn,1),1:size(parametersIn,1)) :: Jacobian
    Real(kind=DoubleReal), Dimension(1:size(parametersIn,1),1:size(pointsIn,1)) :: JacobianT
    Real(kind=DoubleReal), Dimension(1:size(parametersIn,1),1:size(parametersIn,1)) :: JTJ
    Real(kind=DoubleReal), Dimension(1:size(parametersIn,1)) :: JTR, B, P
    Real(kind=DoubleReal), Dimension(1:size(pointsIn,1)) :: Residuals
    Real(kind=DoubleReal) :: h, convThr
    Real(kind=DoubleReal), Dimension(1:size(parametersIn,1)) :: parametersTrial
    Real(kind=DoubleReal) :: fx, fxP, rss, rssA
    Logical :: iterate
! Optional Arguments
    rssA = 0.0D0
    convThr = 1.0D-8
    If(Present(convThr_In))Then
      convThr = convThr_In
    End If
! Init vars
    pSize = size(parametersIn,1)
    nPoints = size(pointsIn,1)
    h = 0.000001D0
    parametersOut = parametersIn
! Iterate
    iterate = .true.
    n = 0
    Do while (iterate)
      n = n + 1
! Cal residuals + rss
      rss = 0.0D0
      Do i=1,nPoints
        fx = NewtonGauss_Call(parametersOut,pointsIn(i,1))
        Residuals(i) = fx - pointsIn(i,2)
        rss = rss + Residuals(i)**2
! Calc Jacobian
        Do j=1,pSize
          parametersTrial = parametersOut
          parametersTrial(j) = parametersOut(j) + h
          fxP = NewtonGauss_Call(parametersTrial,pointsIn(i,1))
          Jacobian(i,j) = (fxP-fx)/h
        End Do
      End Do
! Calc new parameters
      JacobianT = TransposeMatrix(Jacobian)
      JTJ = matmul(JacobianT,Jacobian)
      JTR = matmul(JacobianT,Residuals)
      B = -1.0D0*JTR
      P = SolveLinearSet(JTJ, B)
! save changes
      Do i=1,pSize
        parametersOut(i) = parametersOut(i) + p(i)
      End Do
! RSS check threshold
      If(n.gt.1)Then
        If(abs(rss-rssA).lt.convThr)Then
          iterate = .false.
        End If
      End If
      rssA = rss
! Break out whatever the answer after 100 iterations
      If(n.eq.100)Then
        iterate = .false.
      End If
    End Do
  End Function NewtonGauss

  Function NewtonGauss_RSS(NewtonGauss_Call, parametersIn, pointsIn) RESULT (rss)
! Uses Newton method to solve F(x) = 0 where x is a vector
    Implicit None  !Force declaration of all variables
! Vars:  In - Function call interface
    Interface
      Function NewtonGauss_Call(parameters, x)
        Import DoubleReal
        Real(kind=DoubleReal), Dimension(:), Intent(IN) :: parameters ! x vector
        Real(kind=DoubleReal) :: x
        Real(kind=DoubleReal) :: NewtonGauss_Call ! result
      End function NewtonGauss_Call
    End Interface
! Vars:  In
    Real(kind=DoubleReal), Dimension(:,:) :: pointsIn
    Real(kind=DoubleReal), Dimension(:) :: parametersIn
! Vars:  Out
    Real(kind=DoubleReal) :: rss
! Vars:  Private
    Integer(kind=StandardInteger) :: i
! Calc
    rss = 0.0D0
    Do i=1,size(pointsIn,1)
      rss = rss + NewtonGauss_Call(parametersIn,pointsIn(i,1))
    End Do
  End Function NewtonGauss_RSS

  Function NewtonGauss_Quadratic(parameters,x) RESULT (fx)
! Sample function to solve with two functions and two unknowns
    Implicit None  !Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal), Dimension(:), Intent(IN) :: parameters ! parameter vector
    Real(kind=DoubleReal) :: x                                          ! x val
! Vars:  Out
    Real(kind=DoubleReal) :: fx ! result
! Quadratic function
    fx = parameters(1)*x**2+parameters(2)*x+parameters(3)
  End Function NewtonGauss_Quadratic



!-----------------------------------------------------
! Use Newton method to optimise f(x) = 0
! i.e. find f'(x) = 0
!-----------------------------------------------------

  Function NewtonOptimise(NewtonOpt_Call, xIn, convThr_In) RESULT (xOut)
! Uses Newton method to solve F(x) = 0 where x is a vector
    Implicit None  !Force declaration of all variables
    Interface
      Function NewtonOpt_Call(xVec, pVec)
        Import DoubleReal
        Real(kind=DoubleReal), Dimension(:) :: xVec ! x vector
        Real(kind=DoubleReal), Dimension(:) :: pVec ! parameter vector
        Real(kind=DoubleReal) :: NewtonOpt_Call ! result
      End function NewtonOpt_Call
    End Interface
! Vars:  In
    Real(kind=DoubleReal), Dimension(:) :: xIn
    Real(kind=DoubleReal), Optional :: convThr_In
! Vars:  Out
    Real(kind=DoubleReal), Dimension(1:size(xIn,1)) :: xOut
! Vars:  Private
    Integer(kind=StandardInteger) :: j, n
    Real(kind=DoubleReal) :: convThr
    Real(kind=DoubleReal), Dimension(1:2) :: pVec  ! currently unused
    Real(kind=DoubleReal), Dimension(1:size(xIn,1)) :: eVec
    Real(kind=DoubleReal), Dimension(1:size(xIn,1)) :: xVec
    Real(kind=DoubleReal), Dimension(1:size(xIn,1)) :: gradVec
    Real(kind=DoubleReal), Dimension(1:size(xIn,1),1:size(xIn,1)) :: hessian
    Real(kind=DoubleReal) :: gradSum
    Logical :: iterate
! Set threshold
    convThr = 1.0D-8
    If(Present(convThr_In))Then
      convThr = convThr_In
    End If
! Init
    xVec = xIn
    iterate = .true.
    n = 0
    Do While(iterate)
      n = n + 1
      Call matrixGH(NewtonOpt_Call,xVec,pVec,gradVec,hessian)
      eVec = SolveLinearSet(hessian, gradVec)
      eVec = -1.0D0 * eVec
      xVec = MatAdd(xVec,eVec)
      gradSum = 0.0D0
      Do j=1,size(xIn,1)
        gradSum = gradSum + gradVec(j)
      End Do
      If(gradSum.lt.convThr)Then
        iterate = .false.
      End If
      If(n.eq.100)Then
        iterate = .false.
      End If
    End Do
    xOut = xVec
  End Function NewtonOptimise


  Function NewtonOpt_SampleA(xVec, pVec) RESULT (y)
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
    y = (xVec(1)+4)**2+(xVec(2)-3)**2
  End Function NewtonOpt_SampleA

  Function NewtonOpt_SampleB(xVec, pVec) RESULT (y)
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
    y = (xVec(1)+4)**2+(xVec(2)+xVec(1)-3)**2
  End Function NewtonOpt_SampleB



End Module newtonMethods








!------------------
