! --------------------------------------------------------------!
! RNG Distribution module
! rngDistTypes, rngDist
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
!
!
! ----------------------------------------
! Updated: 20th September 2016
! ----------------------------------------


Module rngDistTypes
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
  Public :: oRandDist, oRandFunction
!------------------------------------------------------------------------------
! Defined Types
!------------------------------------------------------------------------------
  Type :: oRandDist
    Real(kind=DoubleReal) :: lowerBound = 0.0D0
    Real(kind=DoubleReal) :: upperBound = 1.0D0
    Character(Len=8) :: distType = "        "
    Real(kind=DoubleReal), Dimension(1:10) :: parameters = 0.0D0
  End Type oRandDist

  Type :: oRandFunction
    Real(kind=DoubleReal), Dimension(1:1024,1:7) :: table
    Integer(kind=StandardInteger) :: tableSize
    Logical :: tableSet = .false.
    Integer(kind=StandardInteger) :: targetSize = 256
  End Type oRandFunction
End Module rngDistTypes


Module rngDist
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
  Use strings
  Use rng
  Use calcFunctions
  Use fitting
  Use regression
  Use interpolation
  Use rngDistTypes
! Force declaration of all variables
  Implicit None
! Public variables
  Real(kind=DoubleReal) :: randomDist_randNumberG = -1.0D0
  Real(kind=DoubleReal) :: randomDist_randNumberH = -1.0D0
  Real(kind=DoubleReal), Dimension(0:100,1:2) :: randomDist_inverseInt = 0.0D0
  Type(oRandFunction) :: functionTable
! Make private
  Private
! Public
! --functions--!
  Public :: RandomDist
  Public :: RandomDist_GP
  Public :: RandomVaryPoint
  Public :: RandDist
  Public :: RandDistF_Fx
  Public :: RandDistF_Gaussian
  Public :: RandDistF_Sample
  Public :: functionTable
! Interfaces
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------


  Function RandomDist(distTypeIn,setupDistIn,lowerIn,upperIn,sigmaIn) RESULT (output)
! Random float
! F flat distribution
! S square root (diagonal)
! G Gaussian - Box Muller - fits mu=2.5 sigma=0.5 multiplied by 0.05
! H Half Gaussian - Box Muller - fits mu=0.0 sigma=1.0 multiplied by 0.1
! P Test distribution
! M Maxwell-Boltzman distribution
    Implicit None ! Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal) :: lower, upper, sigma
    Character(len=1) :: distType, setupDist  !F flat, G Inverse gaussian distribution
    Real(kind=DoubleReal) :: output
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal) :: randNumber, randNumberB, r, theta, x, y
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Real(kind=DoubleReal), optional :: lowerIn, upperIn, sigmaIn
    Character(len=1), optional :: distTypeIn, setupDistIn
    Real(kind=DoubleReal), Dimension(1:6,1:2) :: fitPoints
    Real(kind=DoubleReal), Dimension(1:21,1:2) :: mbPoints
! Optional arguments
    distType = "F"
    setupDist = "N"
    sigma = 0.4D0
    lower = 0.0D0
    upper = 1.0D0
    If(Present(distTypeIn))Then
      distType = distTypeIn
    End If
    If(Present(setupDistIn))Then
      setupDist = setupDistIn
    End If
    If(Present(sigmaIn))Then
      sigma = sigmaIn
    End If
    If(Present(lowerIn))Then
      lower = lowerIn
    End If
    If(Present(upperIn))Then
      upper = upperIn
    End If
! Get random number 0<= x <=1
    randNumber = RandomLCG()  ! maths.f90
! If Square Root Type
    If(distTypeIn.eq."S")Then
      randNumber = Sqrt(randNumber)
    End If
! If Gaussian Type - full curve, centered on 0.5
! Box Muller method - fits mu=2.5 sigma=0.5 multiplied by 0.05
    If(distTypeIn.eq."G")Then
      If(randomDist_randNumberG.gt.-1.0D0)Then
        randNumber = randomDist_randNumberG
        randomDist_randNumberG = -1.0D0
      Else
! Get second random number
        randNumberB = RandomLCG()
! Calculate x and y
        r = sqrt(-2.0D0*log(randNumber))
        theta = 2.0D0*pi*randNumberB
        x = r*cos(theta)
        y = r*sin(theta)
! Adjust values
        x = (x+5.0D0)/10.0D0
        y = (y+5.0D0)/10.0D0
! store second random number
        randNumber = x
        randomDist_randNumberG = y
      End If
      If(randNumber.lt.lower)Then
        randNumber = lower
      End If
      If(randNumber.gt.upper)Then
        randNumber = upper
      End If
    End If
! If Gaussian Type - half curve, centered on 0.0
! Box Muller method - fits mu=0.0 sigma=1.0 multiplied by 0.1
    If(distTypeIn.eq."H")Then
      If(randomDist_randNumberH.gt.-1.0D0)Then
        randNumber = randomDist_randNumberH
        randomDist_randNumberH = -1.0D0
      Else
! Get second random number
        randNumberB = RandomLCG()  ! maths.f90
! Calculate x and y
        r = sqrt(-2.0D0*log(randNumber))
        theta = 2.0D0*pi*randNumberB
        x = r*cos(theta)
        y = r*sin(theta)
! Adjust values
        x = abs(x)/5.0D0
        y = abs(y)/5.0D0
! store second random number
        randNumber = x
        randomDist_randNumberH = y
      End If
      If(randNumber.lt.lower)Then
        randNumber = lower
      End If
      If(randNumber.gt.upper)Then
        randNumber = upper
      End If
    End If
! Testing dist type - points
    If(distTypeIn.eq."P")Then
      If(setupDist.eq."Y")Then
! Example set of points
        fitPoints(1,1) = 0.0D0
        fitPoints(2,1) = 0.2D0
        fitPoints(3,1) = 0.4D0
        fitPoints(4,1) = 0.6D0
        fitPoints(5,1) = 0.8D0
        fitPoints(6,1) = 1.0D0
        fitPoints(1,2) = 0.0D0
        fitPoints(2,2) = 1.2D0
        fitPoints(3,2) = 0.8D0
        fitPoints(4,2) = 1.8D0
        fitPoints(5,2) = 1.9D0
        fitPoints(6,2) = 0.5D0
        randomDist_inverseInt = RandomDist_GP(fitPoints,"T")  ! maths.f90
      End If
      yArray = PointInterp(randomDist_inverseInt,randNumber,4)  ! maths.f90
      randNumber = yArray(1)
    End If
! Maxwell-Boltzman Distribution
! P(x) = srqt(2/pi)*(x^2*exp(...))
    If(distTypeIn.eq."M")Then
      If(setupDist.eq."Y")Then
! a = 0.25
        Do i=1,21
          mbPoints(i,1) = (i-0)/20.0D0
          mbPoints(i,2) = MaxwellBoltzman(mbPoints(i,1),0.25D0)
        End Do
        randomDist_inverseInt = RandomDist_GP(mbPoints,"T")
      End If
      yArray = PointInterp(randomDist_inverseInt,randNumber,4)  ! maths.f90
      randNumber = yArray(1)
    End If
! Output (adjust to fall in range)
    output = lower + randNumber*(upper-lower)
  End Function RandomDist

  Function RandomDist_GP(inputPoints, integratorIn) RESULT (outputPoints)
! General purpose distribution function
    Implicit None ! Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: inputPoints ! 6+ pairs
    Integer(kind=LongInteger) :: i
    Real(kind=DoubleReal), Dimension(0:100,1:2) :: outputPoints
    Real(kind=DoubleReal), Dimension(1:6) :: coefficients
    Real(kind=DoubleReal), Dimension(1:7) :: coefficientsI
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Real(kind=DoubleReal) :: maxY, dy, y, x, xLast, iVal, iValMax
    Character(len=1), optional :: integratorIn
    Character(len=1) :: integrator
! Init
    integrator = "P"  ! P Polynomial Fit  T Trapezoidal
    xLast = 0.0D0
    iVal = 0.0D0
! Optional arguments
    If(Present(integratorIn))Then
      integrator = integratorIn
    End If
! Polynomial fit + integration
    If(integrator.eq."P")Then
! Fit polynomial to the data points
      coefficients = PolyFit(inputPoints,5)  ! maths.f90
! Integral coefficients
      coefficientsI(1) = 0.0D0
      Do i=1,6
        coefficientsI(i+1) = (1.0D0/(1.0D0*i))*coefficients(i)
      End Do
! Output points
      maxY = CalcPolynomial(coefficientsI, 1.0D0)  ! maths.f90 ! f(x) should be always positive, Int(f(x)) always increasing
      Do i=0,100
        outputPoints(i,1) = CalcPolynomial(coefficientsI, (i/100.0D0)) / maxY
        outputPoints(i,2) = (i/100.0D0)
      End Do
    End If
! Trapezoidal integrate
    If(integrator.eq."T")Then
      dy = 1.0D0/100.0D0
      Do i=0,100
        y = (i/100.0D0)
        outputPoints(i,2) = y
        If(i.eq.0)Then
          iVal = 0.0D0
          outputPoints(i,1) = iVal
          x = inputPoints(1,1)
        Else
          yArray = PointInterp(inputPoints,y,3)  ! maths.f90
          x = yArray(1)
          If(x.lt.0.0D0)Then
            x = 0.0D0
          End If
          iVal = iVal + (0.5D0*dy*(x+xLast))
          outputPoints(i,1) = iVal
        End If
        xLast = x
      End Do
      iValMax = iVal
! Normalize so total integral = 1
      Do i=0,100
        outputPoints(i,1) = outputPoints(i,1)/iValMax
      End Do
    End If
  End Function RandomDist_GP

  Function RandomVaryPoint(pointValue, maxVariation, sigma) RESULT (output)
! Vary a point +/- a max amount using inverse Gaussian distribution
    Implicit None ! Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal) :: pointValue, maxVariation, sigma
    Real(kind=DoubleReal) :: randNumber, variation, output
! Initialise variables
    output = 0.0D0
! Make varied point
    variation = RandomDist("G","N",0.0D0,maxVariation,sigma)
    Call RANDOM_NUMBER(randNumber)
    If(randNumber.gt.0.5D0)Then
      variation = -1.0D0 * variation
    End If
    output = pointValue + variation
  End Function RandomVaryPoint



!-------------------------------------------------------------------------------
!  Random Distribution
!-------------------------------------------------------------------------------

  Function RandDist(input) RESULT (output)
! Random Distribution Function
    Implicit None ! Force declaration of all variables
! Vars:  In
    Type(oRandDist) :: input
! Vars:  Out
    Real(kind=DoubleReal) :: output
! Vars:  Private
! Init
    output = 0.0D0
    input%distType = StrToUpper(input%distType)
! Choose function
    If(input%distType(1:4).eq."FLAT")Then
      output = RandDist_Flat()
    End If
    If(input%distType(1:6).eq."SQROOT")Then
      output = RandDist_SqRoot(input)
    End If
    If(input%distType(1:5).eq."GHEAT")Then
      input%lowerBound = -1.0D0
      input%upperBound = 1.0D0
      input%parameters(1) = 1.0D0 ! sigma
      input%parameters(2) = 0.0D0 ! mu
      input%parameters(3) = 1.0D0 ! factor
      input%parameters(4) = 4.0D0 ! factor
      output = RandDistF_Fx(RandDistF_Gaussian, input)
    End If
    If(input%distType(1:5).eq."GAUSS")Then
      input%parameters(3) = 1.0D0
      input%parameters(4) = 1.0D0 ! factor
      output = RandDistF_Fx(RandDistF_Gaussian, input)
    End If
    If(input%distType(1:7).eq."DBGAUSS")Then
      output = RandDistF_Fx(RandDistF_DoubleGaussian, input)
    End If
    If(input%distType(1:8).eq."MAXBOLTZ")Then
      output = RandDistF_Fx(RandDistF_MaxwellBoltzman, input)
    End If
    If(input%distType(1:6).eq."SAMPLE")Then
      output = RandDistF_Fx(RandDistF_Sample, input)
    End If
  End Function RandDist

! FLAT
  Function RandDist_Flat() RESULT (output)
! Random Distribution - Flat
    Implicit None ! Force declaration of all variables
! Vars:  In
! Vars:  Out
    Real(kind=DoubleReal) :: output
! Return random number
    output = RandomLCG_B()
  End Function RandDist_Flat

! SQUARE ROOT
  Function RandDist_SqRoot(input) RESULT (output)
! Random Distribution - Flat
    Implicit None ! Force declaration of all variables
! Vars:  In
    Type(oRandDist) :: input
! Vars:  Out
    Real(kind=DoubleReal) :: output
! Return random number
    output = RandomFloat(input%lowerBound,input%upperBound)
    output = sqrt(output)
  End Function RandDist_SqRoot


  Subroutine RandDistF_Init(functionTable)
! Reset table
    Implicit None ! Force declaration of all variables
! Vars:  In/Out - Function call interface
    Type(oRandFunction) :: functionTable
! Reset object
    functionTable%table = 0.0D0
    functionTable%tableSet = .false.
    functionTable%tableSize = 0
  End Subroutine


  Function RandDistF_Fx(RandDistF_Call, input) RESULT (randNumber)
! RNG - "Any function"
    Implicit None ! Force declaration of all variables
! Vars:  In - Function call interface
    Interface
      Function RandDistF_Call(pSize, parameters, x)   ! Define the format of the call function
        Import StandardInteger, DoubleReal
        Integer(kind=StandardInteger), Intent(IN) ::  pSize
        Real(kind=DoubleReal), Dimension(1:pSize), Intent(IN) :: parameters ! parameters vector
        Real(kind=DoubleReal), Intent(IN) :: x
        Real(kind=DoubleReal) :: RandDistF_Call ! result
      End function RandDistF_Call
    End Interface
! Vars:  In
    Type(oRandDist) :: input
! Vars:  Out
    Real(kind=DoubleReal) :: randNumber
! Vars:  Private
    Real(kind=DoubleReal) :: uA, uB, uC
    Real(kind=DoubleReal) :: range, area, blockArea, y, yA, yB, yMax, x, xA, x_m, xB, xInc, fx
    Integer(kind=StandardInteger) :: i, randBlock
    Real(kind=DoubleReal), Dimension(1:3,1:2) :: pointsArr
    Logical :: runLoop, loopTrials, findArea
! Build function
    If(.not.functionTable%tableSet)Then
      functionTable%tableSet = .true.
      functionTable%targetSize = 256
! Integrate - estimate the total area
      range = input%upperBound-input%lowerBound
      xInc = range/(1.0D0*(functionTable%targetSize-1))
      xA = input%lowerBound
      xB = xA + xInc
      area = 0.0D0
      Do i=1,functionTable%targetSize
        yA = RandDistF_Call(size(input%parameters,1), input%parameters, xA)
        yB = RandDistF_Call(size(input%parameters,1), input%parameters, xB)
        area = area + 0.5D0*xInc*(yA+yB)
        xA = xB
        xB = xA + xInc
      End Do
      !print *,"Area: ",area
      blockArea = area/(1.0D0*functionTable%targetSize)
      xA = input%lowerBound
      runLoop = .true.
      i = 0
      Do while(runLoop)
        i = i + 1
        yA = RandDistF_Call(size(input%parameters,1), input%parameters, xA)
! Find xB
        xInc = range * 1.0D-4
        xB = xA
        area = 0.0D0
        yMax = yA
        findArea = .true.
        Do while(findArea)
          xB = xB + xInc
          yB = RandDistF_Call(size(input%parameters,1), input%parameters, xB)
          If(yB.gt.yA)Then
            yMax = yB
          End If
          area = area + 0.5D0*xInc*yMax
          If(area.gt.blockArea)Then
            findArea = .false.
          End If
        End Do
        xB = xA + blockArea/(1.0D0*yMax)
        If(xB.gt.input%upperBound)Then
          xB = input%upperBound
          yMax = area/(xB-xA)
          functionTable%tableSize = i
          runLoop = .false.
        End If
! Store
        x_m = 0.5D0*(xA+xB)
        functionTable%table(i,1) = yMax
        functionTable%table(i,2) = xA
        functionTable%table(i,3) = x_m
        functionTable%table(i,4) = xB
        functionTable%table(i,5) = RandDistF_Call(size(input%parameters,1), input%parameters, xA)   ! yA
        functionTable%table(i,6) = RandDistF_Call(size(input%parameters,1), input%parameters, x_m)  ! y_m
        functionTable%table(i,7) = RandDistF_Call(size(input%parameters,1), input%parameters, xB)   ! yB
! Set next xA to xB
        xA = xB
      End Do
    End If
    loopTrials = .true.
    Do While(loopTrials)
      uA = RandomLCG_A()  ! Block val
      uB = RandomLCG_A()  ! x val
      uC = RandomLCG_A()  ! random y val
! Select random block
      If(uA.eq.0.0D0)Then
        randBlock = 1
      Else
        randBlock = Ceiling(uA*functionTable%tableSize)
      End If
! Select random value of x in block
      x = functionTable%table(randBlock,2)+uB*(functionTable%table(randBlock,4)-functionTable%table(randBlock,2))
! Select random val of y from 0 to yMax, to check if within the function
      y = functionTable%table(randBlock,1)*uC
! 3 point interp
      pointsArr(1,1) = functionTable%table(randBlock,2)  ! xA
      pointsArr(2,1) = functionTable%table(randBlock,3)  ! x_m
      pointsArr(3,1) = functionTable%table(randBlock,4)  ! xB
      pointsArr(1,2) = functionTable%table(randBlock,5)  ! yA
      pointsArr(2,2) = functionTable%table(randBlock,6)  ! y_m
      pointsArr(3,2) = functionTable%table(randBlock,7)  ! yB
! Calculate 3 point interp of value of function
      fx = Lagrange_FX(x, pointsArr)
      If(y.gt.fx)Then
        ! Discard
      Else
        ! OK, break loop
        loopTrials = .false.
      End If
      randNumber = x
    End Do
  End Function RandDistF_Fx


  Function RandDistF_Gaussian(pSize, parameters, x) RESULT (fx)
! Gaussian function
    Implicit None ! Force declaration of all variables
! Vars:  In
    Integer(kind=StandardInteger), Intent(IN) ::  pSize
    Real(kind=DoubleReal), Dimension(1:pSize), Intent(IN) :: parameters ! parameters vector
    Real(kind=DoubleReal), Intent(IN) :: x
! Vars:  Out
    Real(kind=DoubleReal) :: fx
! Calculate
    fx = parameters(3)*(0.398942D0/parameters(1))*&
      exp((-0.5D0/(parameters(1)*parameters(1)))*((parameters(4)*(x-parameters(2)))**2))
  End Function RandDistF_Gaussian

  Function RandDistF_DoubleGaussian(pSize, parameters, x) RESULT (fx)
! Gaussian function
    Implicit None ! Force declaration of all variables
! Vars:  In
    Integer(kind=StandardInteger), Intent(IN) ::  pSize
    Real(kind=DoubleReal), Dimension(1:pSize), Intent(IN) :: parameters ! parameters vector
    Real(kind=DoubleReal), Intent(IN) :: x
! Vars:  Out
    Real(kind=DoubleReal) :: fx
! Calculate
    fx = parameters(1)*(0.398942D0/parameters(2))*exp((-0.5D0/(parameters(2)*parameters(2)))*((x-parameters(3))**2))
    fx = fx + parameters(4)*(0.398942D0/parameters(5))*exp((-0.5D0/(parameters(5)*parameters(5)))*((x-parameters(6))**2))
  End Function RandDistF_DoubleGaussian


  Function RandDistF_MaxwellBoltzman(pSize, parameters, x) RESULT (fx)
! Gaussian function
    Implicit None ! Force declaration of all variables
! Vars:  In
    Integer(kind=StandardInteger), Intent(IN) ::  pSize
    Real(kind=DoubleReal), Dimension(1:pSize), Intent(IN) :: parameters ! parameters vector
    Real(kind=DoubleReal), Intent(IN) :: x
! Vars:  Out
    Real(kind=DoubleReal) :: fx
! Calculate
    fx = sqrtTwoPi*((x**2*exp(-1.0D0*((x*x)/(2.0D0*parameters(1)*parameters(1)))))/(parameters(1)**3))
  End Function RandDistF_MaxwellBoltzman


  Function RandDistF_Sample(pSize, parameters, x) RESULT (fx)
! Gaussian function
    Implicit None ! Force declaration of all variables
! Vars:  In
    Integer(kind=StandardInteger), Intent(IN) ::  pSize
    Real(kind=DoubleReal), Dimension(1:pSize), Intent(IN) :: parameters ! parameters vector
    Real(kind=DoubleReal), Intent(IN) :: x
! Vars:  Out
    Real(kind=DoubleReal) :: fx, p
! Calculate
    p = parameters(1) ! fix to stop unused dummy var error 
    fx = 0.2D0+2.0D0*x**2-0.5D0*x**3
  End Function RandDistF_Sample


















End Module rngDist












!
