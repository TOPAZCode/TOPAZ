MODULE ModKinematics
implicit none



type :: Histogram
    integer :: NBins
    real(8) :: BinSize
    real(8) :: LowVal
    real(8) :: SetScale
    real(8),allocatable :: Value(:)
    real(8),allocatable :: Value2(:)
    integer,allocatable :: Hits(:)
    character :: Info*(50)
end type



contains



SUBROUTINE EvalPhasespace(N,EHat,Masses,xRndPS,Mom,PSWgt)
use ModParameters
implicit none
real(8) :: EHat,Masses(1:N),xRndPS(1:3*N-4)
real(8) :: Mom(1:4,1:N),PSWgt
real(8) :: PiWgt2
integer :: N,NVar


   NVar=3*N-4
   PiWgt2 = (2d0*Pi)**(-NVar) * (4d0*Pi)**(N-1)

!  generate PS: massless + massless --> massive(anti-top) + massive(top)
   call genps(N,Ehat,xRndPS(1:NVar),Masses,Mom(1:4,1:N),PSWgt)
   PSWgt = PSWgt*PiWgt2


return
END SUBROUTINE






END MODULE
