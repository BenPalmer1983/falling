! --------------------------------------------------------------!
! Materials module
! materialsTypes, materials
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
!
!
! ----------------------------------------
! Updated: 26th June 2016
! ----------------------------------------


Module isotopesTypes
! Setup Modules
  Use kinds
! Force declaration of all variables
  Implicit None
! Vars:  Module Parameters
  Integer(kind=StandardInteger), Parameter :: p_maxIsotopes = 4000
  Integer(kind=StandardInteger), Parameter :: p_maxElements = 150
! Make private
  Private
! Public Variables and Parameters
  Public :: p_maxIsotopes, p_maxElements
! Public derived data types
  Public :: oIsotopes, oIsotope, oElements, oElement

  Type :: oIsotopes  ! Many
    Integer(kind=StandardInteger) :: isotopeCount
    Character(Len=3), Dimension(1:p_maxIsotopes) :: atomicSymbol
    Integer(kind=StandardInteger), Dimension(1:p_maxIsotopes) :: atomicNumber
    Integer(kind=StandardInteger), Dimension(1:p_maxIsotopes) :: massNumber
    Real(kind=DoubleReal), Dimension(1:p_maxIsotopes) :: isotopicComposition
    Real(kind=DoubleReal), Dimension(1:p_maxIsotopes) :: relativeAtomicMass    ! Mass of this isotope
    Real(kind=DoubleReal), Dimension(1:p_maxIsotopes) :: standardAtomicMass    ! Mass of natuarlly occuring element
    Integer(kind=StandardInteger), Dimension(1:57830) :: keyTable = 0
  End Type oIsotopes

  Type :: oIsotope   ! Singular
    Character(Len=3) :: atomicSymbol
    Integer(kind=StandardInteger) :: atomicNumber
    Integer(kind=StandardInteger) :: massNumber
    Real(kind=DoubleReal) :: isotopicComposition
    Real(kind=DoubleReal) :: relativeAtomicMass    ! Mass of this isotope
    Real(kind=DoubleReal) :: standardAtomicMass    ! Mass of natuarlly occuring element
  End Type oIsotope

  Type :: oElements  ! Many
    Integer(kind=StandardInteger) :: elementCount
    Integer(kind=StandardInteger), Dimension(1:p_maxElements) :: atomicNumber
    Character(Len=3), Dimension(1:p_maxElements) :: atomicSymbol
    Character(Len=32), Dimension(1:p_maxElements) :: elementName
    Integer(kind=StandardInteger), Dimension(1:p_maxElements) :: group
    Integer(kind=StandardInteger), Dimension(1:p_maxElements) :: period
    Real(kind=DoubleReal), Dimension(1:p_maxElements) :: standardAtomicMass
    Real(kind=DoubleReal), Dimension(1:p_maxElements) :: density
  End Type oElements

  Type :: oElement  ! Singular
    Integer(kind=StandardInteger) :: atomicNumber
    Character(Len=3) :: atomicSymbol
    Character(Len=32) :: elementName
    Integer(kind=StandardInteger) :: group
    Integer(kind=StandardInteger) :: period
    Real(kind=DoubleReal) :: standardAtomicMass
    Real(kind=DoubleReal) :: density
  End Type oElement

End Module isotopesTypes


Module isotopes
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use strings
  Use constants
  Use keysMod
  Use isotopesTypes
! Force declaration of all variables
  Implicit None
! Make private
  Private
! Public
  Public :: loadIsotopes          ! loads isotope list into object
  Public :: printIsotopeDetails   ! prints out details for one isitope
  Public :: SearchIsotopes        ! Searches isotope list/db by proton/neutron number
  Public :: loadElements
  Public :: SearchElements
  Public :: LatticeParameter
! Interfaces
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! ------------------------------------------------------------------------!
! Vector Functions
! ------------------------------------------------------------------------!

  Subroutine loadIsotopes(isotopesList)
! Calculates unique key for two keys (order of keys NOT important)
! (A,B) = (B,A)
    Implicit None ! Force declaration of all variables
! Vars:  In/Out
    Type(oIsotopes) :: isotopesList
! Vars:  Private
! Load data
    Call loadIsotopesData(isotopesList)
! Tie keys to data
    Call loadIsotopesKeys(isotopesList)
  End Subroutine loadIsotopes

  Subroutine loadIsotopesKeys(isotopesList)
! Calculates unique key for two keys (order of keys NOT important)
! (A,B) = (B,A)
    Implicit None ! Force declaration of all variables
! Vars:  In/Out
    Type(oIsotopes) :: isotopesList
! Vars:  Private
    Integer(kind=StandardInteger) :: i, keyU
    Integer(kind=StandardInteger) :: protons, neutrons
! Init
! Loop
    Do i=1,isotopesList%isotopeCount
      protons = isotopesList%atomicNumber(i)
      neutrons = isotopesList%massNumber(i) - protons
      keyU = IsotopeKey(protons,neutrons)
      isotopesList%keyTable(keyU) = i
    End Do
  End Subroutine loadIsotopesKeys

  Subroutine printIsotopeDetails(isotopeDetails)
! Calculates unique key for two keys (order of keys NOT important)
! (A,B) = (B,A)
    Implicit None ! Force declaration of all variables
! Vars:  In/Out
    Type(oIsotope) :: isotopeDetails
! Print
    print *,"Symbol:                ",isotopeDetails%atomicSymbol
    print *,"Z:                     ",isotopeDetails%atomicNumber
    print *,"A:                     ",isotopeDetails%massNumber
    print *,"Fraction Composition:  ",isotopeDetails%isotopicComposition
    print *,"Relative atomic mass:  ",isotopeDetails%relativeAtomicMass
    print *,"Standard atomic mass:  ",isotopeDetails%standardAtomicMass
  End Subroutine printIsotopeDetails

  Function SearchIsotopes(protons, neutrons, isotopesList) Result (isotopeDetails)
! Search isotope data
    Implicit None ! Force declaration of all variables
! Vars:  In
    Integer(kind=StandardInteger) :: protons, neutrons
    Type(oIsotopes) :: isotopesList
! Vars:  Out
    Type(oIsotope) :: isotopeDetails
! Vars:  Private
    Integer(kind=StandardInteger) :: keyP, keyN, keyU, keyI
! Store keys
    keyP = protons
    keyN = neutrons
! Get key
    keyU = IsotopeKey(keyP,keyN)
    keyI = isotopesList%keyTable(keyU)
! Store details
    isotopeDetails%atomicSymbol = isotopesList%atomicSymbol(keyI)
    isotopeDetails%atomicNumber = isotopesList%atomicNumber(keyI)
    isotopeDetails%massNumber = isotopesList%massNumber(keyI)
    isotopeDetails%isotopicComposition = isotopesList%isotopicComposition(keyI)
    isotopeDetails%relativeAtomicMass = isotopesList%relativeAtomicMass(keyI)
    isotopeDetails%standardAtomicMass = isotopesList%standardAtomicMass(keyI)
  End Function SearchIsotopes


  Subroutine loadElements(elementsList)
! Calculates unique key for two keys (order of keys NOT important)
! (A,B) = (B,A)
    Implicit None ! Force declaration of all variables
! Vars:  In/Out
    Type(oElements) :: elementsList
! Vars:  Private
! Load data
    Call loadElementsData(elementsList)
  End Subroutine loadElements

  Function SearchElements(elementSymbol, elementsList) Result (elementDetails)
! Search isotope data
    Implicit None ! Force declaration of all variables
! Vars:  In
    Character(*) :: elementSymbol
    Type(oElements) :: elementsList
! Vars:  Out
    Type(oElement) :: elementDetails
! Vars: Private
    Character(Len=3) :: strIn, strCheck
    Integer(kind=StandardInteger) :: i, eN
! Input
    strIn = "   "
    strIn = StrToUpper(elementSymbol)
! Loop through list
    Do i=1,p_maxElements
      strCheck = StrToUpper(elementsList%atomicSymbol(i))
      If(StrMatch(strIn,strCheck))Then
        eN = i
        Exit
      End If
    End Do
! save details
    If(eN.lt.150)Then
      elementDetails%atomicNumber = elementsList%atomicNumber(eN)
      elementDetails%atomicSymbol = elementsList%atomicSymbol(eN)
      elementDetails%elementName = elementsList%elementName(eN)
      elementDetails%group = elementsList%group(eN)
      elementDetails%period = elementsList%period(eN)
      elementDetails%standardAtomicMass = elementsList%standardAtomicMass(eN)
      elementDetails%density = elementsList%density(eN)
    End If
  End Function SearchElements

  Function LatticeParameter(elementDetails, unitSize) Result (aLat)
! Estimate lattice parameter
    Implicit None ! Force declaration of all variables
! Vars:  In
    Type(oElement) :: elementDetails
    Integer(kind=StandardInteger) :: unitSize
! Vars:  Private
    Real(kind=DoubleReal) :: aLat
! Calculation
    aLat = (1.66058D0*((elementDetails%standardAtomicMass*unitSize)/elementDetails%density))**(1/3.0D0)
  End Function LatticeParameter


!---------------------------------------------------------------------------------------------------------------------------------------
! Isotope Data
!---------------------------------------------------------------------------------------------------------------------------------------

  Subroutine loadIsotopesData(isotopesList)
  ! Calculates unique key for two keys (order of keys NOT important)
  ! (A,B) = (B,A)
      Implicit None ! Force declaration of all variables
  ! Vars:  In/Out
      Type(oIsotopes) :: isotopesList
  ! Vars:  Private
      Include "isotopes.isotopeData.f90"
  End Subroutine loadIsotopesData

!---------------------------------------------------------------------------------------------------------------------------------------
! Element Data
!---------------------------------------------------------------------------------------------------------------------------------------

  Subroutine loadElementsData(elementsList)
  ! Loads element data
      Implicit None ! Force declaration of all variables
  ! Vars:  In/Out
      Type(oElements) :: elementsList
  ! Vars:  Private
      Include "isotopes.elementData.f90"
  End Subroutine loadElementsData
End Module isotopes

!-----------------------------------------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------------------------------------
