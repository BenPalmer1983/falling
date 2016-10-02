! ------------------------------------------------------------
!                  STATICCALCS: Print
!                      Print Data
! ------------------------------------------------------------

  Subroutine scPrintEFS(nl)
! Calculate bulk properties
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType), Allocatable, Dimension(:) :: nl
! Vars:  Private
    Character(Len=128) :: tempLine
! Add line to page
    Write(tempLine,*) "Calculation Results: Energy"
    Call addLinePage(tempLine,"T")
    Write(tempLine,*) "NL Size: ",nl%length
    Call addLinePage(tempLine)


  End Subroutine scPrintEFS
