! -------------------------------------------------
!  Include File:   Misc Functions
!
! -------------------------------------------------


  Function F_LoopTemp (loop, loops, startTemp, endTemp) RESULT (loopTemp)
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
  End Function F_LoopTemp













! -------------------------------------------------
