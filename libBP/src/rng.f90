Module rng
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
! Force declaration of all variables
  Implicit None
! Public variables
! RandomLCG
  Integer(kind=LongInteger) :: randomLCG_n=0
  Integer(kind=LongInteger) :: randomLCG_xn
  Integer(kind=LongInteger) :: randomLCG_R_n=0
  Integer(kind=LongInteger) :: randomLCG_R_xn
! RandomLCG_A
  Integer(kind=LongInteger) :: RandomLCG_A_n=0
  Integer(kind=LongInteger) :: randomLCG_A_xn
  Integer(kind=LongInteger) :: RandomLCG_A_Seed=0
! RandomLCG_A
  Integer(kind=LongInteger) :: RandomLCG_B_n=0
  Integer(kind=LongInteger) :: randomLCG_B_xn
  Integer(kind=LongInteger) :: RandomLCG_B_Seed=0
! Make private
  Private
! Public
! ---- Variables
! RandomLCG
  Public :: randomLCG_n
  Public :: randomLCG_xn
  Public :: randomLCG_R_n
  Public :: randomLCG_R_xn
! RandomLCG_A
  Public :: RandomLCG_A_n
  Public :: randomLCG_A_xn
  Public :: RandomLCG_A_Seed
! RandomLCG_B
  Public :: randomLCG_B_n
  Public :: randomLCG_B_xn
  Public :: RandomLCG_B_Seed
! ---- Functions
  Public :: RandomLCG
  Public :: RandomLCG_R
  Public :: RandomInteger
  Public :: RandomFloat
  Public :: RandomLCG_A
  Public :: RandomLCG_B
  Public :: RandomSeed
! ---- Subroutines
! Interfaces
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

  Function RandomLCG(seedIn) RESULT (output)
! Random number - linear congruential generator
    Implicit None ! Force declaration of all variables
! Declare variables
    Integer(kind=LongInteger) :: m, a, c, clockTime
    Integer(kind=LongInteger) :: seed
    Integer(kind=LongInteger), Optional :: seedIn
    Real(kind=DoubleReal) :: output
! init
    seed = 0
    output = 0.0D0
    m = 4294967296_LongInteger
    a = 1103515245_LongInteger
    c = 12345_LongInteger
! Read input, reset counter
    If(Present(seedIn))Then
      seed = seedIn
      If(seed.lt.0)Then ! If seed -1 (for example) get random seed
        Call SYSTEM_CLOCK(clockTime) ! "nano seconds" - well, an estimate
        seed = mod(clockTime,m) ! fit in m
      End If
      randomLCG_n = 0
    End If
! If first iteration
    If(randomLCG_n.eq.0)Then
      If(seed.eq.0)Then
        seed = 12791244+45778951 ! Use default seed
      End If
      randomLCG_n = 0
      randomLCG_xn = seed
    End If
! Increment counter
    randomLCG_n = randomLCG_n + 1
! calculate
    randomLCG_xn = mod((a*randomLCG_xn+c),m)
    output = (1.0D0*randomLCG_xn)/(1.0D0*m)
  End Function RandomLCG

  Function RandomLCG_R() RESULT (output)
! Random number - linear congruential generator
! This function starts with a random seed
    Implicit None ! Force declaration of all variables
! Declare variables
    Integer(kind=LongInteger) :: m, a, c, clockTime
    Integer(kind=LongInteger) :: seed
    Real(kind=DoubleReal) :: output
! init
    seed = 0
    output = 0.0D0
    m = 4294967296_LongInteger
    a = 1103515245_LongInteger
    c = 12345_LongInteger
! If first iteration
    If(randomLCG_R_n.eq.0)Then
! Make "random" seed
      Call SYSTEM_CLOCK(clockTime) ! "nano seconds" - well, an estimate
      seed = RandomSeed()
      randomLCG_R_n = 0
      randomLCG_R_xn = seed
    End If
! Increment counter
    randomLCG_R_n = randomLCG_R_n + 1
! calculate
    randomLCG_R_xn = mod((a*randomLCG_R_xn+c),m)
    output = (1.0D0*randomLCG_R_xn)/(1.0D0*m)
  End Function RandomLCG_R

  Function RandomInteger(lower,upper) RESULT (randInt)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: lower, upper
    Integer(kind=StandardInteger) :: randInt
    Integer(kind=StandardInteger) :: diff, tempInt
    Real(kind=DoubleReal) :: randDouble
! Make random integer
    If(lower.gt.upper)Then
      tempInt = lower
      lower = upper
      upper = tempInt
    End If
    diff = (upper - lower) + 1
! Call RANDOM_NUMBER(randDouble)
    randDouble = RandomLCG()
    randInt = lower+floor(1.0D0*diff*randDouble)
  End Function RandomInteger

  Function RandomFloat(lower,upper) RESULT (randFloat)
! force declaration of all variables
    Implicit None
! declare variables
    Real(kind=DoubleReal) :: lower, upper
    Real(kind=DoubleReal) :: randFloat
    Real(kind=DoubleReal) :: diff, tempFloat
    Real(kind=DoubleReal) :: randDouble
! Make random float
    If(lower.eq.upper)Then
      lower = 0.0D0
    Else If(lower.gt.upper)Then
      tempFloat = lower
      lower = upper
      upper = tempFloat
    End If
    diff = 1.0D0*(upper - lower)
    randDouble = RandomLCG()
    randFloat = lower+1.0D0*diff*randDouble
  End Function RandomFloat


!-----------------------------------------
! Two independent LCG RNGs
!

  Function RandomLCG_A() RESULT (output)
! Random number - linear congruential generator
    Implicit None ! Force declaration of all variables
! Declare variables
    Integer(kind=LongInteger) :: m, a, c
    Integer(kind=LongInteger) :: seed
    Real(kind=DoubleReal) :: output
! init
    seed = 0_LongInteger
    output = 0.0D0
    m = 4294967296_LongInteger
    a = 1103515245_LongInteger
    c = 12345_LongInteger
! Read input, reset counter
! If first iteration
    If(randomLCG_A_n.eq.0)Then
      If(seed.eq.0)Then
        seed = 16778951_LongInteger ! Use default seed
      End If
      randomLCG_A_xn = seed
    End If
! Increment counter
    randomLCG_A_n = randomLCG_A_n + 1
! calculate
    randomLCG_A_xn = mod((a*randomLCG_A_xn+c),m)
    output = (1.0D0*randomLCG_A_xn)/(1.0D0*m)
  End Function RandomLCG_A


  Function RandomLCG_B() RESULT (output)
! Random number - linear congruential generator
    Implicit None ! Force declaration of all variables
! Declare variables
    Integer(kind=LongInteger) :: m, a, c
    Integer(kind=LongInteger) :: seed
    Real(kind=DoubleReal) :: output
! init
    seed = 0_LongInteger
    output = 0.0D0
    m = 4294967296_LongInteger
    a = 1103515245_LongInteger
    c = 12345_LongInteger
! Read input, reset counter
! If first iteration
    If(randomLCG_B_n.eq.0)Then
      If(seed.eq.0)Then
        seed = 37573826_LongInteger ! Use default seed
      End If
      randomLCG_B_xn = seed
    End If
! Increment counter
    randomLCG_B_n = randomLCG_B_n + 1
! calculate
    randomLCG_B_xn = mod((a*randomLCG_B_xn+c),m)
    output = (1.0D0*randomLCG_B_xn)/(1.0D0*m)
  End Function RandomLCG_B

  Function RandomSeed() RESULT (seed)
! Pick random seed
! Typical modulus 4294967296_LongInteger
    Implicit None ! Force declaration of all variables
! Vars:  In
! Vars:  Out
    Integer(kind=LongInteger) :: seed
    Integer(kind=LongInteger), Dimension(1:5) :: seedArr
    Integer(kind=LongInteger) :: exponent
    Integer(kind=LongInteger) :: m, a, c
    Integer(kind=LongInteger) :: clockTime
    Integer(kind=StandardInteger) :: i, j, k
! Vars:  Private

! Init
    m = 4294967296_LongInteger
    a = 1103515245_LongInteger
    c = 12345_LongInteger

    seed = 0_LongInteger
      Do k=1,4
      seedArr = 0_LongInteger
      Do j=1,5
        Do i=1,4
          Call SYSTEM_CLOCK(clockTime)
          seedArr(j) = Mod(seedArr(j)+clockTime,1000)
        End Do
      End Do
      exponent = 0_LongInteger
      Do j=1,5
        exponent = 10_LongInteger**(3_LongInteger*(j-1_LongInteger))
        seed = seed+seedArr(j)*exponent
      End Do
      seed = Mod(seed,m)
    End Do
  End Function RandomSeed



End Module rng
























!
