! ------------------------------------------------------------
!               GEOM: Coords
! ------------------------------------------------------------


!----------------------------------------------
! Make Neighbour List
!----------------------------------------------

  Subroutine expandCoords_UnitCell(coordsIn, coordsOut, expandArr)
! expand a unit cell
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oCoord), Allocatable, Dimension(:) :: coordsIn  ! oCoord (label, labelID, xyz, force, velocity)
    Type(oCoord), Allocatable, Dimension(:) :: coordsOut
    Integer(kind=StandardInteger), Dimension(1:3) :: expandArr
! Vars:  Private
    Integer(kind=StandardInteger) :: x, y, z
    Integer(kind=StandardInteger) :: coordI, coordJ

! Loop through coords
    coordJ = 0
    Do coordI=1,size(coordsIn,1)
      Do x=1,expandArr(1)
        Do y=1,expandArr(2)
          Do z=1,expandArr(3)
            coordJ = coordJ + 1
            If(coordJ.gt.size(coordsOut,1))Then
              Exit
            End If
            coordsOut(coordJ)%label = coordsIn(coordI)%label
            coordsOut(coordJ)%labelID = coordsIn(coordI)%labelID
            coordsOut(coordJ)%xyz(1) = (x-1.0D0+coordsIn(coordI)%xyz(1))/expandArr(1)  ! divide by x count for all
            coordsOut(coordJ)%xyz(2) = (y-1.0D0+coordsIn(coordI)%xyz(2))/expandArr(1)
            coordsOut(coordJ)%xyz(3) = (z-1.0D0+coordsIn(coordI)%xyz(3))/expandArr(1)
            coordsOut(coordJ)%force = coordsIn(coordI)%force
            coordsOut(coordJ)%velocity = coordsIn(coordI)%velocity
          End Do
          If(coordJ.gt.size(coordsOut,1))Then
            Exit
          End If
        End Do
        If(coordJ.gt.size(coordsOut,1))Then
          Exit
        End If
      End Do
      If(coordJ.gt.size(coordsOut,1))Then
        Exit
      End If
    End Do
!----------------------------------------------
  End Subroutine expandCoords_UnitCell



  Subroutine perturbCoords(coordsIn, randDistObj, minPerturbation, maxPerturbation)
! Randomise
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oCoord), Allocatable, Dimension(:) :: coordsIn
    Type(oRandDist) :: randDistObj
    Real(kind=DoubleReal) :: randResult, perturbAmount, minPerturbation, maxPerturbation
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j
!----------------------------------------------
    Do i=1,size(coordsIn,1)
      Do j=1,3
        randResult = RandDist(randDistObj)
        perturbAmount = randResult*(minPerturbation+(maxPerturbation-minPerturbation))
        coordsIn(i)%xyz(j) = coordsIn(i)%xyz(j)+perturbAmount
      End Do
    End Do
!----------------------------------------------
  End Subroutine perturbCoords


  Subroutine geomTransformCoords(coordsIn, transformVecIn, aLat, copyArr)
! Randomise
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: transformVecIn
    Type(oCoord), Allocatable, Dimension(:) :: coordsIn
    Real(kind=DoubleReal), Optional :: aLat
    Integer(kind=StandardInteger), Dimension(1:3), Optional :: copyArr
! Vars:  Out
! Vars:  Private
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: transformVec
    Integer(kind=StandardInteger) :: i, j
! Optional
! Multiply transform vector by aLat
    If(Present(aLat))Then
      Do i=1,3
        Do j=1,3
          transformVec(i,j) = 1.0D0 * aLat * transformVecIn(i,j)
        End Do
      End Do
    Else
      Do i=1,3
        Do j=1,3
          transformVec(i,j) = 1.0D0 * transformVecIn(i,j)
        End Do
      End Do
    End If
! Multiply by "copyArr" if supplied
    If(Present(copyArr))Then
      Do i=1,3
        Do j=1,3
          transformVec(i,j) = transformVec(i,j)*copyArr(i)
        End Do
      End Do
    End If
! Multiply
    Do i=1,size(coordsIn,1)
      coordsIn(i)%xyz = Matmul(transformVec,coordsIn(i)%xyz)
    End Do
!----------------------------------------------
  End Subroutine geomTransformCoords







!
