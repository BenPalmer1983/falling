! ------------------------------------------------------------------------------
!               EAMPA: select program
! ------------------------------------------------------------------------------


Subroutine ProgramSelect()
! Select Program
  Implicit None   ! Force declaration of all variables
! Select and run program
! programObj%program
  print *,programObj%program
  Select Case (programObj%program)
    Case("FALLING")
      Call runFalling()
    Case DEFAULT
      Print *,"No program selected"
  End Select

End Subroutine ProgramSelect







































! ------------------------------------------------------------------------------
