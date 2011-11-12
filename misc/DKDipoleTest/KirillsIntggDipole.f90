PROGRAM CallAmps
    implicit none
integer, parameter :: dp = selected_real_kind(15)
real(dp) :: res,dec_fi_gg


 res = dec_fi_gg(1,1)

 print*, res

END PROGRAM



   double precision function  dec_fi_gg(ij,s)
   implicit none
   integer, parameter :: dp = selected_real_kind(15)
!   real(dp), intent(in)  :: p(4,5)
   integer, intent(in) :: s,ij
   real(dp) :: L,L2,epinv,epinv2,alpha_if,dilog
   real(dp) :: prec(4)
   real(dp):: r,r2,r4,r6,r8,r10
   real(dp) :: x1,x2,x3,x4,x5,x6,x7,x8,x9
   real(dp), parameter :: Pi2 = 9.86960440108935861883449


!-- calculate the invariant mass of the recoiling system / assuming that the recoiling
!-- mass  // we assume that top quark is labeled by the position ``s'' and that
!--- the position of the ``merged'' parton is given by ``ij''


   prec(:) = -p(:,s) -p(:,ij)  ! assumes all are outgoing

!--- input parameters to write the formula




   r2 = scr(prec,prec)/mt**2
   r = dsqrt(r2)

!      r = 0.6d0
!      r2=r**2
!      L=2d0
!      L2 = L**2
!      epinv = 1d-2
!      epinv2 = 1d-2
!      alpha_if = 0.8d0




    r4=r2**2
    r6=r4*r2
    r8=r6*r2
    r10=r8*r2

    x1=alpha_if
    x2=x1*x1
    x3=x2*x1
    x4=x3*x1
    x5=x4*x1
    x6=x5*x1
    x7=x6*x1
    x8=x7*x1
    x9=x8*x1


      L = dlog(mt**2/MuRen**2)


      dec_fi_gg = 0.5d0*epinv*epinv2  + epinv*(17d0/12d0-dlog(1d0-r2)-L/2d0)
      dec_fi_gg =  dec_fi_gg +    &
       + dilog(1d0-r2)-5d0/12d0*Pi2  &
       +1d0/6d0*r2*(12d0*r8*x1+12d0*x1-60d0*x1*r2+84d0*r4*x1-48d0*r6*x1 &
       +6d0*r4*x3-15d0*x2*r4+59d0*r4+6d0*x3*r2-3d0*x2*r2-27d0*r2+11d0*r8 &
       -43d0*r6+12d0)/(-1d0+r2)**5*dlog(r) &
       +dlog(x1)**2-1d0/12d0*(2d0*x1+1d0)*(-1d0+x1)**2*dlog(x1) &
       +dlog(1d0-r2)**2+(L-17d0/6d0)*dlog(1d0-r2) &
       -1d0/4d0*r2*(-1d0+x1)*(2d0*x2*r4+25d0*r4-3d0*r4*x1 &
       -19d0*r2+2d0*x2*r2+x1*r2+4d0*r8+4d0-16d0*r6)/(-1+r2)**5*dlog(1d0-x1+x1*r2) &
       -1d0/240d0*(-892d0+393d0*r8*x3-67d0*x2-40d0*r10*x2-272d0*x2*r4+358d0*x2*r2 &
       +205d0*x3-283d0*x4+241d0*x5+300d0*L2*r8*x1-141d0*x6-852d0*r8-602d0*r6*x2 &
       -115d0*x7*r10+42d0*x8*r10-8d0*x9*r10-852d0*r10*x1+2040d0*L*r4 &
       +263d0*r8*x2+240d0*L2*r2-60d0*L2-300d0*L2*x1*r2+340d0*L*r10*x1-60d0*L2*r10*x1 &
       -5112d0*r4+340d0*L+3708d0*r2+892d0*x1-60d0*L2*r8-1360d0*L*r6+3388d0*r6-97d0*r10*x3 &
       +340d0*L*r8-4680d0*x1*r2+8520d0*r4*x1-8200d0*r6*x1+4200d0*r8*x1-3400d0*L*r4*x1 &
       +802d0*x7*r4-348d0*x8*r4+80d0*x9*r4+3400d0*L*r6*x1-1700d0*L*r8*x1-80d0*x9*r6 &
       -898d0*x7*r6+507d0*x7*r8+372d0*x8*r6+40d0*x9*r8-198d0*x8*r8-363d0*x7*r2 &
       +162d0*x8*r2-40d0*x9*r2-340d0*L*x1+1700d0*L*x1*r2-360d0*L2*r4+730d0*x6*r2 &
       +1630d0*x6*r6-1520d0*x6*r4-30d0*x8+8d0*x9+67d0*x7-1360d0*L*r2+216d0*x6*r10 &
       -915d0*x6*r8+60d0*L2*x1-305d0*x5*r10+1291d0*x5*r8-2524d0*x5*r6 &
       +2516d0*x5*r4-1219d0*x5*r2+1360d0*x4*r2-3050d0*x4*r4+600d0*L2*r4*x1 &
       -600d0*L2*r6*x1-1195d0*x4*r8+2650d0*x4*r6+278d0*x4*r10-917d0*x3*r2+2078d0*r4*x3 &
       -1062d0*r6*x3+240d0*L2*r6)/(1d0-x1+x1*r2)/(1d0-r2)**4




   end function dec_fi_gg




FUNCTION dilog(xIn)
IMPLICIT NONE
double precision xIn,x,z,z2
double precision dilog, Li2tmp
double precision Const,Fact
double precision Pi26,Pi23
parameter (Pi26=1.6449340668482264d0)
parameter (Pi23=3.2898681336964528d0)
double precision B1,B2,B4,B6,B8,B10,B12,B14,B16,B18
parameter (B1= -0.25d0)
parameter (B2=  2.7777777777777778d-02)
parameter (B4= -2.7777777777777778d-04)
parameter (B6=  4.7241118669690098d-06)
parameter (B8= -9.1857730746619635d-08)
parameter (B10= 1.8978869988970999d-09)
parameter (B12=-4.0647616451442255d-11)
parameter (B14= 8.9216910204564525d-13)
parameter (B16=-1.9939295860721075d-14)
parameter (B18= 4.5189800296199181d-16)

 Const=0d0
 Fact =1d0

 x = xIn
 if ( x.gt.1d0 ) then
   Fact = -1d0
   Const = Pi23-0.5d0*dlog(x)**2
   x = 1d0/x
 elseif ( x.lt.-1d0 ) then
   Fact = -1d0
   Const =-Pi26-0.5d0*dlog(-x)**2
   x = 1d0/x
 elseif ( x.eq.1d0) then
   dilog = Pi26
   return
 endif

 if ( x.gt.0.5d0 ) then
   Fact = -1d0*Fact
   Const = Const + Pi26 - dlog(x)*dlog(1d0-x)
   x = 1d0-x
 endif

 z = -dlog(1d0-x)
 Li2tmp = z
 z2 = z*z
 Li2tmp = Li2tmp + B1  * z2
 z = z*z2
 Li2tmp = Li2tmp + B2  * z
 z = z*z2
 Li2tmp = Li2tmp + B4  * z
 z = z*z2
 Li2tmp = Li2tmp + B6  * z
 z = z*z2
 Li2tmp = Li2tmp + B8  * z
 z = z*z2
 Li2tmp = Li2tmp + B10 * z
 z = z*z2
 Li2tmp = Li2tmp + B12 * z
 z = z*z2
 Li2tmp = Li2tmp + B14 * z
 z = z*z2
 Li2tmp = Li2tmp + B16 * z
 z = z*z2
 Li2tmp = Li2tmp + B18 * z

 dilog = Fact*Li2tmp + Const

return
END FUNCTION



















