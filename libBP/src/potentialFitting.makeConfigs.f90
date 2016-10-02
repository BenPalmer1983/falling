! ------------------------------------------------------------
!               POTENTIALFITTING: Make Configs
!                Make config from input files
! ------------------------------------------------------------

  Subroutine pfMakeNL(pfObj, pfCoords, pfNl)
! Make config files
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oPotentialFitting) :: pfObj
    Type(coordsType), Allocatable, Dimension(:) :: pfCoords
    Type(nlType), Allocatable, Dimension(:) :: pfNl
! Vars:  Private
    Integer(kind=StandardInteger) :: cKey, i, j
! Allocate arrays
    Allocate(pfCoords(1:pfObj%fileCount))
    Allocate(pfNl(1:pfObj%fileCount))

! Copy coords
    Do cKey=1,pfObj%fileCount
      pfCoords(cKey)%length = pfObj%configs(cKey)%coordLength
      pfCoords(cKey)%aLat = pfObj%configs(cKey)%aLat
      pfCoords(cKey)%unitCell = 0.0D0
      Do i=1,3
        pfCoords(cKey)%unitCell(i,i) = 1.0D0
      End Do
      Do i=1,pfObj%configs(cKey)%coordLength
        pfCoords(cKey)%label(i) = pfObj%configs(cKey)%coords(i)%label
        Do j=1,3
          pfCoords(cKey)%coords(i,j) = pfObj%configs(cKey)%coords(i)%xyz(j)
        End Do
      End Do
    End Do



  End Subroutine pfMakeNL



  Subroutine pfLoadPotential(pfObj, potential)
! Make config files
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oPotentialFitting) :: pfObj
    Type(potentialType) :: potential
! Vars:  Private
! Init potential
    Call addLinePage("Init potential")
    Call initPotential(potential)        ! potentials.f90
! Load potential
    Call addLinePage("Load potential")
    Call loadPotential(pfObj%potentialFile, potential)  ! potentials.f90
  End Subroutine pfLoadPotential
