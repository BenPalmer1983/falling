! ------------------------------------------------------------
!               POTENTIALFITTING: Calculations
!             Run calculations to calculate RSS
! ------------------------------------------------------------




  Subroutine pfRunCalc(pfObj, pfCoords, pfNl, pfPotential)
! count configs (where coords gt 0) and create a map
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oPotentialFitting) :: pfObj
    Type(coordsType), Allocatable, Dimension(:) :: pfCoords
    Type(nlType), Allocatable, Dimension(:) :: pfNl
    Type(potentialType) :: pfPotential
! Vars:  Private
    Integer(kind=StandardInteger) :: i, cKey
    Type(oMpi) :: mpiObj
    Character(Len=128) :: calcOutFile
! Init mpi obj (ID, count)
    Call m_initMpi(mpiObj)
! Init output file name
    calcOutFile = trim(pfObj%outputDirectory)//"/calcOutput.out"
! Count calculation
    pfObj%calcCount = pfObj%calcCount + 1
! Make output file
    If((pfObj%calcFile).and.(pfObj%calcCount.eq.1))Then
      If(mpiObj%ID.eq.0)Then
        open(unit=999,file=trim(calcOutFile))
        write(999,"(A1)") " "
        close(999)
      End If
    End If


! Loop through configs
    Do i=1,pfObj%configCount
      pfObj%calcCount_TotalEFS = pfObj%calcCount_TotalEFS + 1
      cKey = pfObj%configMap(i)
!
      Call makeNL(pfNl, pfCoords, 6.5D0,cKey)       ! geom.f90
      Call calcEFS(pfNl, pfPotential,cKey)

      If(mpiObj%ID.eq.0)Then
        Open(UNIT=999,FILE=Trim(calcOutFile),&
        status="old",position="append",action="write")

        write(999,"(A16,I8)") "Calculation:    ",pfObj%calcCount
        write(999,"(A16,I8)") "Config:         ",i
        write(999,"(A30,A30)") "Reference","Calculated"
        write(999,"(A16,F14.6,A16,F14.6)") &
          "Energy:         ",pfObj%configs(cKey)%energy,"              ",pfNl(cKey)%totalEnergy

      ! Close file
        Close(999)
      End If
    End Do

  End Subroutine pfRunCalc
















!---------------------------------------------------------------------------------------
