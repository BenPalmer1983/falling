Module time
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
! Force declaration of all variables
  Implicit None
! Declare public variables
  Real(kind=DoubleReal) :: g_startTime
! Make private
  Private
! Public
! --variables--!
  Public :: g_startTime
! --Subroutines--!
  Public :: setStartTime
  Public :: getTime
! Interfaces
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

  Subroutine setStartTime()
! Set start time for g_startTime global variable
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
! Vars:  Private
    Real(kind=DoubleReal) :: cpuTime
! Get CPU time
    CALL CPU_TIME(cpuTime)
! Set start time
    g_startTime = cpuTime
  End Subroutine setStartTime


  Subroutine getTime(cpuTime)
! Init the unit coords data type
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Real(kind=DoubleReal) :: cpuTime
! Vars:  Private
! Set time
    CALL CPU_TIME(cpuTime)
    cpuTime = cpuTime-g_startTime
  End Subroutine getTime

End Module time
