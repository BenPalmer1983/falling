! --------------------------------------------------------------!
! Simulated Annealing module
! simulatedAnnealingTypes, simulatedAnnealing
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Optimises using simulated annealing
!
! ----------------------------------------
! Updated: 21st May 2016
! ----------------------------------------

Module simulatedAnnealingTypes
! Setup Modules
  Use kinds
! Force declaration of all variables
  Implicit None
! Vars:  Module Parameters
!  Integer(kind=StandardInteger), Parameter :: p_
! Make private
  Private
! Public Variables and Parameters
!  Public :: p_
! Public derived data types
!  Public :: o

End Module simulatedAnnealingTypes

Module simulatedAnnealing
! Setup Modules
  Use kinds
  Use rng
  Use constants
  Use matrix
  Use calcFunctions
! Force declaration of all variables
  Implicit None
! Public variables
! Make private
  Private
! Public
! --variables--!
! --functions--!
  Public :: saLoopTemp
  Public :: saAccept
  Public :: saVaryMax


! Interfaces
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------


! ----------------------------------------
! Module Subroutines
! ----------------------------------------



! ----------------------------------------
! Module Functions
! ----------------------------------------

  Function saLoopTemp (loop, loops, startTemp, endTemp) RESULT (loopTemp)
! Simulated annealing loop temperature
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Integer(kind=StandardInteger) :: loop, loops
    Real(kind=DoubleReal) :: startTemp, endTemp
! Vars:  Out
    Real(kind=DoubleReal) :: loopTemp
! Vars:  Private
!
    loopTemp = startTemp/((startTemp/endTemp)**(1.0D0*((loop-1)/(1.0D0*(loops-1)))))
  End Function saLoopTemp

  Function saAccept (optVal, testVal, temp) RESULT (accept)
! Simulated annealing loop temperature
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal) :: optVal, testVal, temp
! Vars:  Out
    Integer(kind=StandardInteger) :: accept
! Vars:  Private
    Real(kind=DoubleReal) :: randNumber, probThreshold
! Check
    If(testVal.lt.optVal)Then
      accept = 1  ! Accept, better value
    Else
      accept = 0
      randNumber = RandomLCG()
      probThreshold = exp(-1.0D0*((testVal-optVal)/temp))
      If(probThreshold.le.randNumber)Then
        accept = 2 ! Accept, worse value, but give it a chance
      End If
    End If
  End Function saAccept

  Function saVaryMax (loop, loops, startVaryMax, endVaryMax) RESULT (loopTemp)
! Simulated annealing loop temperature
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Integer(kind=StandardInteger) :: loop, loops
    Real(kind=DoubleReal) :: startVaryMax, endVaryMax
! Vars:  Out
    Real(kind=DoubleReal) :: loopTemp
! Vars:  Private
! Calculate
    loopTemp = startVaryMax/((startVaryMax/endVaryMax)**(1.0D0*((loop-1)/(1.0D0*(loops-1)))))
  End Function saVaryMax



End Module simulatedAnnealing











!-----------------------------------------------
