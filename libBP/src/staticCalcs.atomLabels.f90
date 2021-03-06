! -------------------------------------------------
!  Include File:   Atom Labels
!
! -------------------------------------------------

  Subroutine atomLabelIDs(potential, coords)
! Look through the atom labels from the potential and coordinates and assign unique
! atom IDs - case insensitive
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Type(potentialType) :: potential
    Type(coordsType), Dimension(:) :: coords
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, cKey, k
    Character(Len=16) :: tempLabel
    Character(Len=16), Dimension(1:64) :: tempLabels
! Blank tempLabels array
    tempLabels = BlankStringArray(tempLabels)
! Loop through potential functions
    Do i=1,potential%fCount
      tempLabel = StrToUpper(potential%atomLabel_A(i))
      If(IsBlank(tempLabel))Then
        potential%atomID_A(i) = 0
      Else
        Do j=1,Size(tempLabels,1)
          If(StrMatch(tempLabel,tempLabels(j)))Then
            potential%atomID_A(i) = j
            Exit
          End If
          If(IfStringEmpty(tempLabels(j)))Then
            tempLabels(j) = tempLabel
            potential%atomID_A(i) = j
            Exit
          End If
        End Do
      End If
      tempLabel = StrToUpper(potential%atomLabel_B(i))
      If(IsBlank(tempLabel))Then
        potential%atomID_B(i) = 0
      Else
        Do j=1,Size(tempLabels,1)
          If(StrMatch(tempLabel,tempLabels(j)))Then
            potential%atomID_B(i) = j
            Exit
          End If
          If(IfStringEmpty(tempLabels(j)))Then
            tempLabels(j) = tempLabel
            potential%atomID_B(i) = j
            Exit
          End If
        End Do
      End If
    End Do
! Loop through coords
    Do cKey=1,Size(coords,1)
      Do i=1,coords(cKey)%length
        tempLabel = StrToUpper(coords(cKey)%label(i))
        Do j=1,Size(tempLabels,1)
          If(StrMatch(tempLabel,tempLabels(j)))Then
            coords(cKey)%labelID(i) = j
            Exit
          End If
          If(IfStringEmpty(tempLabels(j)))Then
            tempLabels(j) = tempLabel
            coords(cKey)%labelID(i) = j
            Exit
          End If
        End Do
      End Do
    End Do
! Store IDs
    k = 0
    Do j=1,64
      If(tempLabels(j)(1:1).eq." ")Then
        Exit
      End If
      k = k + 1
      potential%atomIDs(j) = tempLabels(j)
    End Do
    potential%atomID_Count = k
    Do cKey=1,Size(coords,1)
      Do j=1,64
        coords(cKey)%atomIDs(j) = tempLabels(j)
      End Do
      coords(cKey)%atomID_Count = k
    End Do
! Fill in uKey for potentials
    Do i=1,potential%fCount
      If(potential%atomID_B(i).eq.0)Then
        potential%uKey(i) = potential%atomID_A(i)
      Else
        potential%uKey(i) = DoubleKey(potential%atomID_A(i),potential%atomID_B(i))
      End If
    End Do
  End Subroutine atomLabelIDs


  Subroutine printAtomLabelIDs(coords)
! Print Label
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Type(coordsType), Dimension(:) :: coords
! Vars:  Private
    Integer(kind=StandardInteger) :: i
    Character(Len=64) :: tempLine
! Print out atom ids
    Call addLinePage("Atom Labels and IDs","T")
    Do i=1,32
      If(IfStringEmpty(coords(1)%atomIDs(i)))Then
        Exit
      Else
        Write(tempLine,*) i,coords(1)%atomIDs(i)
        Call addLinePage(tempLine)
      End If
    End Do
  End Subroutine printAtomLabelIDs




!-------------------------------
