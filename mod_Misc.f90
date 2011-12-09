MODULE ModMisc
implicit none


INTERFACE OPERATOR (.dot.)
   MODULE PROCEDURE MinkowskyProduct
   MODULE PROCEDURE MinkowskyProduct_128
   MODULE PROCEDURE MinkowskyProductC
   MODULE PROCEDURE MinkowskyProductC_128
   MODULE PROCEDURE MinkowskyProductRC
   MODULE PROCEDURE MinkowskyProductRC_128
   MODULE PROCEDURE MinkowskyProductCR
   MODULE PROCEDURE MinkowskyProductCR_128
END INTERFACE OPERATOR (.dot.)

INTERFACE Dot
   MODULE PROCEDURE MinkowskyProduct
   MODULE PROCEDURE MinkowskyProduct_128
   MODULE PROCEDURE MinkowskyProductC
   MODULE PROCEDURE MinkowskyProductC_128
   MODULE PROCEDURE MinkowskyProductRC
   MODULE PROCEDURE MinkowskyProductRC_128
   MODULE PROCEDURE MinkowskyProductCR
   MODULE PROCEDURE MinkowskyProductCR_128
END INTERFACE Dot

INTERFACE OPERATOR (.cross.)
   MODULE PROCEDURE VectorCross
END INTERFACE OPERATOR (.cross.)

INTERFACE sc_
    module procedure sc__64
    module procedure sc__128
END INTERFACE sc_

INTERFACE sc
    module procedure sc_64
    module procedure sc_128
END INTERFACE sc


INTERFACE go_Gauss
    module procedure go_Gauss_64
    module procedure go_Gauss_128
END INTERFACE go_Gauss



INTERFACE swapMom
    module procedure swapMomR
    module procedure swapMomC
END INTERFACE swapMom



INTERFACE ASpiStr! works only in Weyl conventions
   module procedure ASpiStr0
   module procedure ASpiStr1
   module procedure ASpiStr2
   module procedure ASpiStr3
END INTERFACE ASpiStr

INTERFACE SpiStr! works only in Weyl conventions
   module procedure SpiStr0
   module procedure SpiStr1
   module procedure SpiStr2
   module procedure SpiStr3
END INTERFACE SpiStr




    complex(8),private :: u(8,16), ux(8,16), uy(8,16), uz(8,16),uyc(8,16)
        complex(8),private :: one,cone,mcone,mone,sqrt2,msqrt2,csqrt2,cmsqrt2
        parameter(one=1d0)
        parameter(mone=-1d0)
        parameter(sqrt2 = 0.70710678118654752440084436210484903929d0)
        parameter(msqrt2=-0.70710678118654752440084436210484903929d0)
        parameter(cone = (0d0,1d0))
        parameter(csqrt2=cone*sqrt2)
        parameter(cmsqrt2=cone*msqrt2)
        parameter(mcone = (0d0,-1d0))

        data u(1,1)/ one/, u(2,2)/ one/, u(3,5)/ one/,u(4,6)/ one/, u(5,9)/ one/, u(6,10)/ one/, u(7,13)/ one/, u(8,14)/ one/
        data uz(1,1)/ one/, uz(1,3)/ one/,uz(2,2)/ one/, uz(2,4)/ mone/,uz(3,5)/ one/, uz(3,7)/ one/,uz(4,6)/ one/, uz(4,8)/ mone/,  &
                            uz(5,9)/ one/, uz(5,11)/ one/,uz(6,10)/ one/, uz(6,12)/ mone/,uz(7,13)/ one/, uz(7,15)/ one/,uz(8,14)/ one/, uz(8,16)/ mone/
        data ux(1,1)/ sqrt2/, ux(1,2)/ sqrt2/,ux(1,3)/ sqrt2/, ux(1,4)/ sqrt2/,ux(2,1)/ msqrt2/, ux(2,2)/ sqrt2/,ux(2,3)/ sqrt2/,  ux(2,4)/ msqrt2/,  &
                    ux(3,5)/ sqrt2/,  ux(3,6)/ sqrt2/,ux(3,7)/ sqrt2/,  ux(3,8)/ sqrt2/,ux(4,5)/ msqrt2/, ux(4,6)/ sqrt2/,ux(4,7)/ sqrt2/,  ux(4,8)/ msqrt2/,  &
                    ux(5,9)/ sqrt2/, ux(5,10)/ sqrt2/,ux(5,11)/ sqrt2/, ux(5,12)/ sqrt2/,ux(6,9)/ msqrt2/, ux(6,10)/ sqrt2/,ux(6,11)/ sqrt2/, ux(6,12)/ msqrt2/,  &
                    ux(7,13)/ sqrt2/, ux(7,14)/ sqrt2/,ux(7,15)/ sqrt2/, ux(7,16)/ sqrt2/,ux(8,13)/ msqrt2/, ux(8,14)/ sqrt2/,ux(8,15)/ sqrt2/, ux(8,16)/ msqrt2/

     data uy(1,1)/ sqrt2/,  uy(1,2)/ sqrt2/,uy(1,3)/ cmsqrt2/, uy(1,4)/ csqrt2/,    &
          uy(2,1)/ msqrt2/, uy(2,2)/ sqrt2/,uy(2,3)/ cmsqrt2/, uy(2,4)/ cmsqrt2/,   &
          uy(3,5)/ sqrt2/,  uy(3,6)/ sqrt2/,uy(3,7)/ cmsqrt2/, uy(3,8)/ csqrt2/,    &
          uy(4,5)/ msqrt2/, uy(4,6)/ sqrt2/,uy(4,7)/ cmsqrt2/, uy(4,8)/ cmsqrt2/,   &
          uy(5,9)/ sqrt2/,  uy(5,10)/ sqrt2/,uy(5,11)/ cmsqrt2/, uy(5,12)/ csqrt2/,   &
          uy(6,9)/ msqrt2/, uy(6,10)/ sqrt2/,uy(6,11)/ cmsqrt2/, uy(6,12)/ cmsqrt2/,  &
          uy(7,13)/ sqrt2/, uy(7,14)/ sqrt2/,uy(7,15)/ cmsqrt2/, uy(7,16)/ csqrt2/,   &
          uy(8,13)/ msqrt2/,uy(8,14)/ sqrt2/,uy(8,15)/ cmsqrt2/, uy(8,16)/ cmsqrt2/

     data uyc(1,1)/ sqrt2/,  uyc(1,2)/ sqrt2/,uyc(1,3)/ csqrt2/, uyc(1,4)/ cmsqrt2/,    &
          uyc(2,1)/ msqrt2/, uyc(2,2)/ sqrt2/,uyc(2,3)/ csqrt2/, uyc(2,4)/ csqrt2/,   &
          uyc(3,5)/ sqrt2/,  uyc(3,6)/ sqrt2/,uyc(3,7)/ csqrt2/, uyc(3,8)/ cmsqrt2/,    &
          uyc(4,5)/ msqrt2/, uyc(4,6)/ sqrt2/,uyc(4,7)/ csqrt2/, uyc(4,8)/ csqrt2/,   &
          uyc(5,9)/ sqrt2/,  uyc(5,10)/ sqrt2/,uyc(5,11)/ csqrt2/, uyc(5,12)/ cmsqrt2/,   &
          uyc(6,9)/ msqrt2/, uyc(6,10)/ sqrt2/,uyc(6,11)/ csqrt2/, uyc(6,12)/ csqrt2/,  &
          uyc(7,13)/ sqrt2/, uyc(7,14)/ sqrt2/,uyc(7,15)/ csqrt2/, uyc(7,16)/ cmsqrt2/,   &
          uyc(8,13)/ msqrt2/,uyc(8,14)/ sqrt2/,uyc(8,15)/ csqrt2/, uyc(8,16)/ csqrt2/


contains





FUNCTION ASpiStr0(ASpi)
implicit none
complex(8) :: ASpiStr0(1:4)
complex(8) :: ASpi(1:4)

   ASpiStr0(1:4) = ASpi(1:4)

END FUNCTION


FUNCTION ASpiStr1(ASpi,Sl1)! works only in Weyl conventions
implicit none
complex(8) :: ASpiStr1(1:4)
complex(8) :: ASpi(1:4),Sl1(1:4)


   ASpiStr1(1:4) = spb2_Weyl(ASpi(1:4),Sl1(1:4))

! print *, "check"
! print *, ASpiStr1
!
!      Tmp(1:4) = vbqq_Weyl(4,ASpi(1:4),Spi(1:4))
!      ASpiStr1 = psp1_(Tmp(1:4),Sl1(1:4))
! print *, ASpiStr1
! pause

END FUNCTION


FUNCTION ASpiStr2(ASpi,Sl1,Sl2)! works only in Weyl conventions
implicit none
complex(8) :: ASpiStr2(1:4)
complex(8) :: ASpi(1:4),Sl1(1:4),Sl2(1:4)


   ASpiStr2(1:4) = spb2_Weyl(ASpi(1:4),Sl1(1:4))
   ASpiStr2(1:4) = spb2_Weyl( ASpiStr2(1:4),Sl2(1:4))

END FUNCTION


FUNCTION ASpiStr3(ASpi,Sl1,Sl2,Sl3)! works only in Weyl conventions
implicit none
complex(8) :: ASpiStr3(1:4)
complex(8) :: ASpi(1:4),Sl1(1:4),Sl2(1:4),Sl3(1:4)

   ASpiStr3(1:4) = spb2_Weyl(ASpi(1:4),Sl1(1:4))
   ASpiStr3(1:4) = spb2_Weyl( ASpiStr3(1:4),Sl2(1:4))
   ASpiStr3(1:4) = spb2_Weyl( ASpiStr3(1:4),Sl3(1:4))

END FUNCTION



FUNCTION SpiStr0(Spi)
implicit none
complex(8) :: SpiStr0(1:4)
complex(8) :: Spi(1:4)

   SpiStr0(1:4) = Spi(1:4)

END FUNCTION


FUNCTION SpiStr1(Sl1,Spi)! works only in Weyl conventions
implicit none
complex(8) :: SpiStr1(1:4)
complex(8) :: Spi(1:4),Sl1(1:4)

   SpiStr1(1:4) = spi2_Weyl(Sl1(1:4),Spi(1:4))

END FUNCTION


FUNCTION SpiStr2(Sl2,Sl1,Spi)! works only in Weyl conventions
implicit none
complex(8) :: SpiStr2(1:4)
complex(8) :: Spi(1:4),Sl1(1:4),Sl2(1:4)


   SpiStr2(1:4) = spi2_Weyl(Sl1(1:4),Spi(1:4))
   SpiStr2(1:4) = spi2_Weyl(Sl2(1:4),SpiStr2(1:4))

END FUNCTION


FUNCTION SpiStr3(Sl3,Sl2,Sl1,Spi)! works only in Weyl conventions
implicit none
complex(8) :: SpiStr3(1:4)
complex(8) :: Spi(1:4),Sl1(1:4),Sl2(1:4),Sl3(1:4)

   SpiStr3(1:4) = spi2_Weyl(Sl1(1:4),Spi(1:4))
   SpiStr3(1:4) = spi2_Weyl(Sl2(1:4),SpiStr3(1:4))
   SpiStr3(1:4) = spi2_Weyl(Sl3(1:4),SpiStr3(1:4))

END FUNCTION



FUNCTION VectorProd(p1,p2)
implicit none
real(8), intent(in) :: p1(1:3),p2(1:3)
real(8)             :: VectorProd

    VectorProd  = p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)
END FUNCTION VectorProd



FUNCTION VectorCross(p1,p2)
implicit none
real(8), intent(in) :: p1(1:3),p2(1:3)
real(8)             :: VectorCross(1:3)

    VectorCross(1)  = p1(2)*p2(3) - p1(3)*p2(2)
    VectorCross(2)  = p1(3)*p2(1) - p1(1)*p2(3)
    VectorCross(3)  = p1(1)*p2(2) - p1(2)*p2(1)
END FUNCTION VectorCross



FUNCTION MinkowskyProduct(p1,p2)
implicit none
real(8), intent(in) :: p1(1:4),p2(1:4)
real(8)             :: MinkowskyProduct

     MinkowskyProduct = p1(1)*p2(1)  &
                      - p1(2)*p2(2)  &
                      - p1(3)*p2(3)  &
                      -p1(4)*p2(4)
END FUNCTION MinkowskyProduct


FUNCTION MinkowskyProduct_128(p1,p2)
implicit none
real(16), intent(in) :: p1(1:4),p2(1:4)
real(16)             :: MinkowskyProduct_128

     MinkowskyProduct_128 = p1(1)*p2(1)  &
                      - p1(2)*p2(2)  &
                      - p1(3)*p2(3)  &
                       -p1(4)*p2(4)
END FUNCTION MinkowskyProduct_128


FUNCTION MinkowskyProductC(p1,p2)
implicit none
complex(8), intent(in) :: p1(1:4),p2(1:4)
complex(8)             :: MinkowskyProductC

     MinkowskyProductC = p1(1)*p2(1)  &
                       - p1(2)*p2(2)  &
                       - p1(3)*p2(3)  &
                       - p1(4)*p2(4)
END FUNCTION MinkowskyProductC



FUNCTION MinkowskyProductC_128(p1,p2)
implicit none
complex(16), intent(in) :: p1(1:4),p2(1:4)
complex(16)             :: MinkowskyProductC_128

     MinkowskyProductC_128 = p1(1)*p2(1)  &
                       - p1(2)*p2(2)  &
                       - p1(3)*p2(3)  &
                       - p1(4)*p2(4)
END FUNCTION MinkowskyProductC_128


FUNCTION MinkowskyProductRC(p1,p2)
implicit none
real(8), intent(in) :: p1(1:4)
complex(8), intent(in) :: p2(1:4)
complex(8)             :: MinkowskyProductRC

     MinkowskyProductRC = dcmplx(p1(1))*p2(1)  &
                       - dcmplx(p1(2))*p2(2)  &
                       - dcmplx(p1(3))*p2(3)  &
                       - dcmplx(p1(4))*p2(4)
END FUNCTION MinkowskyProductRC


FUNCTION MinkowskyProductRC_128(p1,p2)
implicit none
real(16), intent(in) :: p1(1:4)
complex(16), intent(in) :: p2(1:4)
complex(16)             :: MinkowskyProductRC_128

     MinkowskyProductRC_128 = dcmplx(p1(1))*p2(1)  &
                       - dcmplx(p1(2))*p2(2)  &
                       - dcmplx(p1(3))*p2(3)  &
                       - dcmplx(p1(4))*p2(4)
END FUNCTION MinkowskyProductRC_128


FUNCTION MinkowskyProductCR(p1,p2)
implicit none
real(8), intent(in) :: p2(1:4)
complex(8), intent(in) :: p1(1:4)
complex(8)             :: MinkowskyProductCR

     MinkowskyProductCR = dcmplx(p2(1))*p1(1)  &
                       - dcmplx(p2(2))*p1(2)  &
                       - dcmplx(p2(3))*p1(3)  &
                       - dcmplx(p2(4))*p1(4)
END FUNCTION MinkowskyProductCR



FUNCTION MinkowskyProductCR_128(p1,p2)
implicit none
real(16), intent(in) :: p2(1:4)
complex(16), intent(in) :: p1(1:4)
complex(16)             :: MinkowskyProductCR_128

     MinkowskyProductCR_128 = dcmplx(p2(1))*p1(1)  &
                       - dcmplx(p2(2))*p1(2)  &
                       - dcmplx(p2(3))*p1(3)  &
                       - dcmplx(p2(4))*p1(4)
END FUNCTION MinkowskyProductCR_128



FUNCTION Fac(N)
implicit none
integer :: N, Fac,i

   if    ( N.le.1 ) then
      Fac=1
   elseif( N.eq.2 ) then
      Fac=2
   elseif( N.eq.3 ) then
      Fac=6
   elseif( N.eq.4 ) then
      Fac=24
   elseif( N.eq.5 ) then
      Fac=120
   elseif( N.eq.6 ) then
      Fac=720
   elseif( N.ge.7 ) then
      Fac=720
      do i=7,N
         Fac = Fac*i
      enddo
   endif
END FUNCTION



FUNCTION Binomial(N,K)
implicit none
integer :: N,K, Binomial

   Binomial = Fac(N)/Fac(K)/Fac(N-K)
END FUNCTION





FUNCTION Lambda(x,y,z)
implicit none
real(8) :: Lambda,x,y,z

    Lambda = x**2 + y**2 + z**2 - 2d0*x*y - 2d0*x*z - 2d0*y*z

END FUNCTION



FUNCTION SqrtLambda(x,y,z)
implicit none
real(8) :: SqrtLambda,x,y,z

    SqrtLambda = dsqrt(x**2 + y**2 + z**2 - 2d0*x*y - 2d0*x*z - 2d0*y*z)

END FUNCTION






function DLi2(x)
implicit none
double precision DLi2,x,y,t,s,a,pi3,pi6,zero,one,half,malf,mone,mtwo
double precision c(0:18),h,alfa,b0,b1,b2
integer i
data zero /0.0d0/, one /1.0d0/
data half /0.5d0/, malf /-0.5d0/, mone /-1.0d0/, mtwo /-2.0d0/
data pi3 /3.289868133696453d0/, pi6 /1.644934066848226d0/
data c( 0) / 0.4299669356081370d0/
data c( 1) / 0.4097598753307711d0/
data c( 2) /-0.0185884366501460d0/
data c( 3) / 0.0014575108406227d0/
data c( 4) /-0.0001430418444234d0/
data c( 5) / 0.0000158841554188d0/
data c( 6) /-0.0000019078495939d0/
data c( 7) / 0.0000002419518085d0/
data c( 8) /-0.0000000319334127d0/
data c( 9) / 0.0000000043454506d0/
data c(10) /-0.0000000006057848d0/
data c(11) / 0.0000000000861210d0/
data c(12) /-0.0000000000124433d0/
data c(13) / 0.0000000000018226d0/
data c(14) /-0.0000000000002701d0/
data c(15) / 0.0000000000000404d0/
data c(16) /-0.0000000000000061d0/
data c(17) / 0.0000000000000009d0/
data c(18) /-0.0000000000000001d0/

      if(x .eq. one) then
       DLi2=pi6
       return
      else if(x .eq. mone) then
       DLi2=malf*pi6
       return
      end if
      t=-x
      if(t .le. mtwo) then
       y=mone/(one+t)
       s=one
       a=-pi3+half*(log(-t)**2-log(one+one/t)**2)
      else if(t .lt. mone) then
       y=mone-t
       s=mone
       a=log(-t)
       a=-pi6+a*(a+log(one+one/t))
      else if(t .le. malf) then
       y=(mone-t)/t
       s=one
       a=log(-t)
       a=-pi6+a*(malf*a+log(one+t))
      else if(t .lt. zero) then
       y=-t/(one+t)
       s=mone
       a=half*log(one+t)**2
      else if(t .le. one) then
       y=t
       s=one
       a=zero
      else
       y=one/t
       s=mone
       a=pi6+half*log(t)**2
      end if

      h=y+y-one
      alfa=h+h
      b1=zero
      b2=zero
      do i = 18,0,-1
          b0=c(i)+alfa*b1-b2
          b2=b1
          b1=b0
      enddo
      DLi2=-(s*(b0-h*b2)+a)


return
END FUNCTION




! FUNCTION DLi2(xIn)! something is wrong if x is between +1d0 and +2d0
! implicit none
! double precision xIn,x,z,z2
! double precision DLi2,Li2tmp
! double precision Const,Fact
! double precision Pi26,Pi23
! parameter (Pi26=1.6449340668482264d0)
! parameter (Pi23=3.2898681336964528d0)
! double precision B1,B2,B4,B6,B8,B10,B12,B14,B16,B18
! parameter (B1= -0.25d0)
! parameter (B2=  2.7777777777777778d-02)
! parameter (B4= -2.7777777777777778d-04)
! parameter (B6=  4.7241118669690098d-06)
! parameter (B8= -9.1857730746619635d-08)
! parameter (B10= 1.8978869988970999d-09)
! parameter (B12=-4.0647616451442255d-11)
! parameter (B14= 8.9216910204564525d-13)
! parameter (B16=-1.9939295860721075d-14)
! parameter (B18= 4.5189800296199181d-16)
!
!   Const=0d0
!   Fact =1d0
!
! if(xIn.gt.1d0 .and. xIn.lt.2d0) call Error("")
!
!   x = xIn
!   if ( x.gt.1d0 ) then
!     Fact = -1d0
!     Const = Pi23-0.5d0*dlog(x)**2
!     x = 1d0/x
!   elseif ( x.lt.-1d0 ) then
!     Fact = -1d0
!     Const =-Pi26-0.5d0*dlog(-x)**2
!     x = 1d0/x
!   elseif ( x.eq.1d0) then
!     DLi2 = Pi26
!     return
!   endif
!
!   if ( x.gt.0.5d0 ) then
!     Fact = -1d0*Fact
!     Const = Const + Pi26 - dlog(x)*dlog(1d0-x)
!     x = 1d0-x
!   endif
!
!   z = -dlog(1d0-x)
!   Li2tmp = z
!   z2 = z*z
!   Li2tmp = Li2tmp + B1  * z2
!   z = z*z2
!   Li2tmp = Li2tmp + B2  * z
!   z = z*z2
!   Li2tmp = Li2tmp + B4  * z
!   z = z*z2
!   Li2tmp = Li2tmp + B6  * z
!   z = z*z2
!   Li2tmp = Li2tmp + B8  * z
!   z = z*z2
!   Li2tmp = Li2tmp + B10 * z
!   z = z*z2
!   Li2tmp = Li2tmp + B12 * z
!   z = z*z2
!   Li2tmp = Li2tmp + B14 * z
!   z = z*z2
!   Li2tmp = Li2tmp + B16 * z
!   z = z*z2
!   Li2tmp = Li2tmp + B18 * z
!
!   DLi2 = Fact*Li2tmp + Const
!
! return
! END FUNCTION





!    subroutine CLi2(zxdilo,ipi12,x,ieps)! complex dilog from LoopTools
!    implicit none
!    integer ipi12,ieps,ier
!    double precision x
!    double complex zxdilo,zlog
!    integer jsgn
!    double precision fact,u,u2,dfflo1,ffbnd,a,xdilo,xprec,bdn02,bdn05,bdn10,bdn15
!    double complex cy,cfact
!    double precision :: xlg2 = .6931471805599453094172321214581d0
!    double precision :: pi = 3.1415926535897932384626433832795028842d0,bf(20)
!    double precision,parameter :: precx=4d-16
!    double precision,parameter :: xloss=1d0/8d0,xalogm=2.d-30, xalog2=xalogm**2
!    external ffbnd,dfflo1
!    save xprec,bdn02,bdn05,bdn10,bdn15
!    data xprec /-1d0/
!
!    bf(1) = - 1.d+0/4.D+0
!    bf(2) = + 1.D+0/36.D+0
!    bf(3) = - 1.D+0/36.D+2
!    bf(4) = + 1.D+0/21168.D+1
!    bf(5) = - 1.D+0/108864.D+2
!    bf(6) = + 1.D+0/52690176.D+1
!    bf(7) = - 691.D+0/16999766784.D+3
!    bf(8) = + 1.D+0/1120863744.D+3
!    bf(9) = - 3617.D+0/18140058832896.D+4
!    bf(10) = + 43867.D+0/97072790126247936.D+3
!    bf(11) = - 174611.D+0/168600109166641152.D+5
!    bf(12) = + 77683.D+0/32432530090601152512.D+4
!    bf(13) = - 236364091.D+0/4234560341829359173632.D+7
!    bf(14) = + 657931.D+0/5025632054039239458816.D+6
!    bf(15) = - 3392780147.D+0/109890470493622010006470656.D+7
!    bf(16)=+172.3168255201D+0/2355349904102724211909.3102313472D+6
!    bf(17)=-770.9321041217D+0/4428491985594062112714.2791446528D+8
!    bf(18)=( 0.4157635644614046176D-28)
!    bf(19)=(-0.9962148488284986022D-30)
!    bf(20)=( 0.2394034424896265390D-31)
!
!
!    if ( xprec .ne. precx ) then
!        xprec = precx
!        bdn02 = ffbnd(1,2,bf)
!        bdn05 = ffbnd(1,5,bf)
!        bdn10 = ffbnd(1,10,bf)
!        bdn15 = ffbnd(1,15,bf)
!    endif
!    if ( x .eq. 1) then
!        zxdilo = 0
!        zlog = -99999
!        ipi12 = 2
!        return
!    elseif (x .eq. -1) then
!        zxdilo = 0
!        zlog = xlg2
!        ipi12 = -1
!        return
!    elseif (x .eq. .5D0) then
!        zxdilo = - xlg2**2/2
!        zlog = -xlg2
!        ipi12 = 1
!        return
!    elseif ( abs(x) .lt. precx ) then
!        zxdilo = x
!        zlog = -x
!        ipi12 = 0
!        return
!    endif
!    if (x .lt. -1) then
!        fact = log(-x)
!        cy = - fact**2/2
!        ipi12 = -2
!        if ( -x*xloss .gt. 1 ) then
!       u = -dfflo1(1/x,ier)
!        else
!       u = -log(1-1/x)
!        endif
!        zlog = log(1-x)
!        jsgn = -1
!    elseif ( x .lt. .5D0) then
!        cy = 0
!        ipi12 = 0
!        if ( abs(x) .lt. xloss ) then
!       zlog = dfflo1(x,ier)
!        else
!       zlog = log(1-x)
!        endif
!        u = -DBLE(zlog)
!        jsgn = 1
!    elseif ( x .le. 2 ) then
!        u = -log(x)
!        if ( abs(1-x) .lt. xalogm ) then
!       cy = 0
!        elseif ( x .lt. 1 ) then
!       zlog = log(1-x)
!       cy = DBLE(u)*zlog
!        elseif ( ieps .gt. 0 ) then
!       zlog = DCMPLX(log(x-1),-pi)
!       cy = DBLE(u)*zlog
!        else
!       zlog = DCMPLX(log(x-1),+pi)
!       cy = DBLE(u)*zlog
!        endif
!        ipi12 = 2
!        jsgn = -1
!    else
!        if ( ieps .gt. 0 ) then
!       cfact = DCMPLX(log(x),-pi)
!       zlog = DCMPLX(log(x-1),-pi)
!        else
!       cfact = DCMPLX(log(x),+pi)
!       zlog = DCMPLX(log(x-1),+pi)
!        endif
!        cy = - cfact**2/2
!        ipi12 = -2
!        if ( x*xloss .gt. 1 ) then
!       u = -dfflo1(1/x,ier)
!        else
!       u = -log(1-1/x)
!        endif
!        jsgn = -1
!    endif
!    if ( abs(u) .lt. xalog2 ) then
!        xdilo = u
!    else
!    u2 = u**2
!    a = abs(u2)
!    if ( a .gt. bdn15 ) then
!        xdilo = u2*(bf(16) + u2*(bf(17) + u2*(bf(18) +u2*(bf(19) + u2*bf(20) ))))
!    else
!        xdilo = 0
!    endif
!    if ( a .gt. bdn10 ) then
!        xdilo = u2*(bf(11) + u2*(bf(12) + u2*(bf(13) + u2*(bf(14) + u2*(bf(15) + xdilo)))))
!    endif
!    if ( a .gt. bdn05 ) then
!        xdilo = u2*(bf(6) + u2*(bf(7) + u2*(bf(8) + u2*(bf(9) + u2*(bf(10) + xdilo)))))
!    endif
!    if ( a .gt. bdn02 ) then
!        xdilo = u2*(bf(3) + u2*(bf(4) + u2*(bf(5) + xdilo)))
!    endif
!    xdilo = u + u2*(bf(1) + u*(bf(2) + xdilo))
!    endif
!    if(jsgn.eq.1)then
!        zxdilo =  DBLE(xdilo) + cy
!    else
!        zxdilo = -DBLE(xdilo) + cy
!    endif
!    end subroutine









       double precision function vl(x1,x2,x3)
       real(8), intent(in) ::  x1,x2,x3
       vl = x1**2+x2**2+x3**2-2d0*x1*x2-2d0*x1*x3-2d0*x2*x3
       end function vl


       double precision function vel(p1,p2)
       real(8), intent(in) ::  p1(4),p2(4)
       real(8) :: p12(4),x,y,z
       p12 = p1+p2
       x = scr(p12,p12)
       y = scr(p1,p1)
       z = scr(p2,p2)
       vel = sqrt(vl(x,y,z))/(x-y-z)
       end function vel


       double precision function scr(p1,p2)
       real(8), intent(in) :: p1(4), p2(4)
       scr = p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)-p1(4)*p2(4)
       end function scr

       double complex function scc(cp1,p2)
       complex(8), intent(in) :: cp1(4)
       real(8), intent(in) :: p2(4)
       scc = cp1(1)*p2(1)-cp1(2)*p2(2)-cp1(3)*p2(3)-cp1(4)*p2(4)
       end function scc


 


! split two sets into INTERsection and COMPLement
! e.g.: A=(1,2,3,4,5), B=(2,4,5)
!       INTER=(2,4,5), COMPL=(1,3), numMatch=3
! healthy requirements: size(A).ge.size(B), size(INTER)=size(B), size(COMPL)=size(A)
FUNCTION MatchSets(A,B,INTER,COMPL)
implicit none
integer :: MatchSets
integer, intent(in) :: A(:),B(:)
integer, intent(out):: INTER(:),COMPL(:)
integer :: n,numMismatch,numMatch

  numMatch=0; numMismatch=0;
  do n=1,size(A,1)
    if( any(B(:).eq.A(n)) ) then
        numMatch=numMatch+1
        INTER(numMatch) = A(n)
    else
        numMismatch=numMismatch+1
        COMPL(numMismatch) = A(n)
    endif
  enddo

  MatchSets = numMatch

RETURN
END FUNCTION





FUNCTION IsNan(x)
implicit none
logical IsNan
real(8) :: x

   if( .not.x.le.0d0 .and. .not.x.gt.0d0 ) then
       IsNaN=.true.
   else
      IsNaN=.false.
   endif
END FUNCTION





FUNCTION StepFunc(x)
implicit none
real(8) :: StepFunc,x

   StepFunc = 0d0
   if( x.ge.0d0 ) then
       StepFunc = 1d0
   endif
END FUNCTION




FUNCTION IsAScalar(Type)
implicit none
logical IsAScalar
integer :: Type

   if( abs(Type).eq.15 .or. abs(Type).eq.16 ) then
      IsAScalar = .true.
   else
      IsAScalar = .false.
   endif
END FUNCTION



FUNCTION IsAQuark(Type)
implicit none
logical IsAQuark
integer :: Type

   if( abs(Type).le.6 .and. Type.ne.0 ) then
      IsAQuark = .true.
   else
      IsAQuark = .false.
   endif
END FUNCTION




FUNCTION IsABoson(Type)
implicit none
logical IsABoson
integer :: Type

   if( abs(Type).eq.11 .or. abs(Type).eq.12 .or. abs(Type).eq.13 ) then
      IsABoson = .true.
   else
      IsABoson = .false.
   endif
END FUNCTION




FUNCTION QuarkToAntiQuark(Type)
implicit none
integer :: QuarkToAntiQuark,Type

   if( IsAQuark(Type) .and. Type.ge.1 ) then
     QuarkToAntiQuark = -Type
   else
     QuarkToAntiQuark = Type
   endif
END FUNCTION



FUNCTION ChargeConj(Part)
implicit none
integer :: ChargeConj,Part

   if( IsAQuark(Part) .or. IsAScalar(Part) ) then
     ChargeConj =-Part
   else
     ChargeConj = Part
   endif
END FUNCTION



FUNCTION Get_ET(Mom)
implicit none
real(8) ::Mom(1:4),Get_ET,sinTheta

   sinTheta = dsqrt(1d0-Mom(4)**2/(Mom(2)**2+Mom(3)**2+Mom(4)**2)) ! = pT/|pVec|
   Get_ET = Mom(1) * sinTheta

RETURN
END FUNCTION



FUNCTION Get_PT(Mom)
implicit none
real(8) ::Mom(1:4),Get_PT

   Get_PT = dsqrt( Mom(2)**2 + Mom(3)**2 )

RETURN
END FUNCTION


FUNCTION Get_PT2(Mom)
implicit none
real(8) ::Mom(1:4),Get_PT2

   Get_PT2 = Mom(2)**2 + Mom(3)**2

RETURN
END FUNCTION


FUNCTION Get_MInv(Mom)
implicit none
real(8) ::Mom(1:4),Get_MInv

   Get_MInv = dsqrt( Mom(1:4).dot.Mom(1:4) )

RETURN
END FUNCTION


FUNCTION Get_MInv2(Mom)
implicit none
real(8) ::Mom(1:4),Get_MInv2

   Get_MInv2 = Mom(1:4).dot.Mom(1:4)

RETURN
END FUNCTION


FUNCTION Get_ETA(Mom)
implicit none
real(8) ::Mom(1:4),Get_ETA

   Get_ETA = 0.5d0*dlog( (Mom(1)+Mom(4))/(Mom(1)-Mom(4)) )

RETURN
END FUNCTION


FUNCTION Get_PseudoETA(Mom)
implicit none
real(8) ::Mom(1:4),Get_PseudoETA

   Get_PseudoETA = 0.5d0*dlog( (dsqrt(Mom(2)**2+Mom(3)**2+Mom(4)**2)+Mom(4))/(dsqrt(Mom(2)**2+Mom(3)**2+Mom(4)**2)-Mom(4)) )

RETURN
END FUNCTION



FUNCTION Get_PHI(Mom)
implicit none
real(8) ::Mom(1:4),Get_PHI

   Get_PHI = datan2( Mom(3),Mom(2) )

RETURN
END FUNCTION


FUNCTION Get_CosAlpha(Mom1,Mom2)
implicit none
real(8) ::Mom1(1:4),Mom2(1:4),Get_CosAlpha

    Get_CosAlpha = (Mom1(2)*Mom2(2)+Mom1(3)*Mom2(3)+Mom1(4)*Mom2(4))/dsqrt(Mom1(2)**2+Mom1(3)**2+Mom1(4)**2)/dsqrt(Mom2(2)**2+Mom2(3)**2+Mom2(4)**2)

RETURN
END FUNCTION



FUNCTION Get_R(Mom1,Mom2)
implicit none
real(8),parameter :: Pi=3.141592653589793d0
real(8) :: Mom1(1:4),Mom2(1:4),Get_R
real(8) :: eta1,eta2,phi1,phi2,DeltaPhi,r2,delphi

   eta1 = Get_ETA(Mom1(1:4))
   eta2 = Get_ETA(Mom2(1:4))
   phi1 = Get_PHI(Mom1(1:4))
   phi2 = Get_PHI(Mom2(1:4))
   DeltaPhi = dabs(phi1-phi2)
   if( DeltaPhi.gt.Pi ) DeltaPhi=2d0*Pi-DeltaPhi
   Get_R = dsqrt((eta1-eta2)**2 + DeltaPhi**2)

!       r2= (Mom1(2)*Mom2(2)+Mom1(3)*Mom2(3)) /dsqrt((Mom1(2)**2+Mom1(3)**2)*(Mom2(2)**2+Mom2(3)**2))
!       delphi=dacos(r2)
!       print *, delphi-DeltaPhi
RETURN
END FUNCTION






FUNCTION Get_MT(Mom12,MomMiss)
implicit none
real(8) :: Get_MT,Get_MTsq
real(8) :: Mom12(1:4),MomMiss(1:4)
real(8) :: pT_12,Minv_12sq,ET_12,pT_miss,pT_x,pT_y

    pT_12 = get_PT(Mom12(1:4))
    Minv_12sq = Mom12(1:4).dot.Mom12(1:4)

    ET_12 = dsqrt(Mom12(1)**2-Mom12(4)**2) !=dsqrt(pT_12**2 + Minv_12sq)
    pT_miss = get_PT(MomMiss(1:4))
    pT_x = Mom12(2)+MomMiss(2)
    pT_y = Mom12(3)+MomMiss(3)
    Get_MTsq = (ET_12+pT_miss)**2 - (pT_x**2+pT_y**2)  ! there exist different definitions with pT_miss and ET_miss

    Get_MT = dsqrt(Get_MTsq)
RETURN
END FUNCTION




FUNCTION Vertex_WPW(MomP,MomW,PolP,PolW)! W-W-Photon vertex, MomP,MomW outgoing, missing prefactor: +i*e*Q_(W^sig), W^sig --> W^sig + gamma
implicit none
complex(8) :: Vertex_WPW(1:4)
complex(8) :: PolP(1:4),PolW(1:4)
real(8) :: MomP(1:4),MomW(1:4),Mom1(1:4),Mom2(1:4),Mom3(1:4)

  Mom1(1:4) = MomP(1:4) - MomW(1:4)
  Mom2(1:4) = MomP(1:4) + 2d0*MomW(1:4)
  Mom3(1:4) = 2d0*MomP(1:4) + MomW(1:4)

  Vertex_WPW(1:4) = - ( Mom1(1:4) * (PolP(1:4).dot.PolW(1:4)) &
                      + PolW(1:4) * (Mom2(1:4).dot.PolP(1:4)) &
                      - PolP(1:4) * (Mom3(1:4).dot.PolW(1:4)) )

RETURN
END FUNCTION
!    eval_TripVert(1:Dv) = IOverSqrt2 * ( (k1(1:Dv)-k2(1:Dv))*(v1.Ndot.v2)  &
!                            - 2d0*v1(1:Dv)*(k1.Ndot.v2)  &
!                            + 2d0*v2(1:Dv)*(k2.Ndot.v1) )


FUNCTION WritePartType( ParticleList )
implicit none
integer :: NumPart
integer :: PartType,i,ParticleList(:)
character(9)   :: Temp(1:9)
character(100) :: WritePartType

   NumPart = size(ParticleList,dim=1)

   if(NumPart.gt.9) call Error("NumPart too large.")

   do i=1,NumPart
      PartType = ParticleList(i)
      if    ( PartType .eq. 1 ) then
         Temp(i) = " up "
      elseif( PartType .eq. 2 ) then
         Temp(i) = " dn "
      elseif( PartType .eq. 3 ) then
         Temp(i) = " chm "
      elseif( PartType .eq. 4 ) then
         Temp(i) = " str "
      elseif( PartType .eq. 5 ) then
         Temp(i) = " top "
      elseif( PartType .eq. 6 ) then
         Temp(i) = " bot "
      elseif( PartType .eq. 10) then
         Temp(i) = " glu "
      elseif( PartType .eq. 11) then
         Temp(i) = " pho "
      elseif( PartType .eq. 12) then
         Temp(i) = " Z "
      elseif( PartType .eq. 13) then
         Temp(i) = " Wp "
      elseif( PartType .eq. 14) then
         Temp(i) = " HTop "
      elseif( PartType .eq. 15) then
         Temp(i) = " Stop "
      elseif( PartType .eq. 16) then
         Temp(i) = " Sbot "
      elseif( PartType .eq. 0 ) then
         Temp(i) = " discd "
      elseif ( PartType .eq. -1 ) then
         Temp(i) = " Aup "
      elseif( PartType .eq. -2 ) then
         Temp(i) = " Adn "
      elseif( PartType .eq. -3 ) then
         Temp(i) = " Achm "
      elseif( PartType .eq. -4 ) then
         Temp(i) = " Astr "
      elseif( PartType .eq. -5 ) then
         Temp(i) = " Atop "
      elseif( PartType .eq. -6 ) then
         Temp(i) = " Abot "
      elseif( PartType .eq. -13) then
         Temp(i) = " Wm "
      elseif( PartType .eq. -14) then
         Temp(i) = " AHTop "
      elseif( PartType .eq. -15) then
         Temp(i) = " AStop "
      elseif( PartType .eq. -16) then
         Temp(i) = " ASbot "
      else
         Temp(i) = "error"
      endif
   enddo

   if ( NumPart.eq.1 ) WritePartType = Temp(1)
   if ( NumPart.eq.2 ) WritePartType = Temp(1)//Temp(2)
   if ( NumPart.eq.3 ) WritePartType = Temp(1)//Temp(2)//Temp(3)
   if ( NumPart.eq.4 ) WritePartType = Temp(1)//Temp(2)//Temp(3)//Temp(4)
   if ( NumPart.eq.5 ) WritePartType = Temp(1)//Temp(2)//Temp(3)//Temp(4)//Temp(5)
   if ( NumPart.eq.6 ) WritePartType = Temp(1)//Temp(2)//Temp(3)//Temp(4)//Temp(5)//Temp(6)
   if ( NumPart.eq.7 ) WritePartType = Temp(1)//Temp(2)//Temp(3)//Temp(4)//Temp(5)//Temp(6)//Temp(7)
   if ( NumPart.eq.8 ) WritePartType = Temp(1)//Temp(2)//Temp(3)//Temp(4)//Temp(5)//Temp(6)//Temp(7)//Temp(8)
   if ( NumPart.eq.9 ) WritePartType = Temp(1)//Temp(2)//Temp(3)//Temp(4)//Temp(5)//Temp(6)//Temp(7)//Temp(8)//Temp(9)
   return
END FUNCTION




FUNCTION check_MatSol(N,matrix,sol,vec)
implicit none
complex(8) :: matrix(1:N,1:N),sol(1:N),vec(1:N)
complex(8) :: check_MatSol(1:N)
integer, intent(in) :: N
integer :: i,j

  do i=1,N
    check_MatSol(i) = (0d0,0d0)
    do j=1,N
        check_MatSol(i) = check_MatSol(i) + matrix(i,j)*sol(j)
    enddo
    check_MatSol(i) = check_MatSol(i) - vec(i)
  enddo

return
END FUNCTION




FUNCTION go_Gauss_64(N,matrix)
implicit none
complex(8), intent(in) :: matrix(1:N,1:N+1)
integer, intent(in) :: N
complex(8) :: go_Gauss_64(1:N)
complex(8) :: work(1:N,1:N+1)
integer :: i,j
complex(8) :: mult, tsum

    work(1:N,1:N+1) = matrix(1:N,1:N+1)
    do i = 1,N
      if (cdabs(work(i,i)) <= 1.0e-6) then
!          print *, i,N
!          print *, "go_Gauss, Zero pivot element: ",cdabs(work(i,i))
!          call Error("go_Gauss, Zero pivot element")
      endif
      do j = i+1,N
        mult = work(j,i)/work(i,i)
        work(j,1:N+1) = work(j,1:N+1) - mult*work(i,1:N+1)
      end do
    end do
    do i=N,1,-1
      tsum = work(i,N+1)
      do j=i+1,N
        tsum = tsum - work(i,j)*work(j,N+1)
      end do
      work(i,N+1) = tsum/work(i,i)
    end do
    go_Gauss_64(1:N) = work(1:N,N+1)
    return
END FUNCTION





FUNCTION go_Gauss_128(N,matrix)
implicit none
complex(16), intent(in) :: matrix(1:N,1:N+1)
integer, intent(in) :: N
complex(16) :: go_Gauss_128(1:N)
complex(16) :: work(1:N,1:N+1)
integer :: i,j
complex(16) :: mult, tsum

    work(1:N,1:N+1) = matrix(1:N,1:N+1)
    do i = 1,N
      if (cqabs(work(i,i)) <= 1.0e-6) then
!          print *, i,N
!          print *, "go_Gauss, Zero pivot element: ",cdabs(work(i,i))
!          call Error("go_Gauss, Zero pivot element")
      endif
      do j = i+1,N
        mult = work(j,i)/work(i,i)
        work(j,1:N+1) = work(j,1:N+1) - mult*work(i,1:N+1)
      end do
    end do
    do i=N,1,-1
      tsum = work(i,N+1)
      do j=i+1,N
        tsum = tsum - work(i,j)*work(j,N+1)
      end do
      work(i,N+1) = tsum/work(i,i)
    end do
    go_Gauss_128(1:N) = work(1:N,N+1)
    return
END FUNCTION



SUBROUTINE printYRnd(y)
implicit none
real(8) :: y(:)
integer :: n


do n=1,size(y)
  write (*,"(A,I2,A,1F20.16,A)") "yRnd(",n,")=",y(n),"d0"
enddo

END SUBROUTINE



SUBROUTINE SwitchEnergyComponent(p)
implicit none
real(8) :: p(:,:),px,py,pz,E
integer :: NMom,i

  NMom = size(p,2)
  do i=1,NMom
      E =p(1,i)
      px=p(2,i)
      py=p(3,i)
      pz=p(4,i)
      p(1,i) = px
      p(2,i) = py
      p(3,i) = pz
      p(4,i) = E
  enddo

END SUBROUTINE



SUBROUTINE SwitchEnergyComponentBack(p)
implicit none
real(8) :: p(:,:),px,py,pz,E
integer :: NMom,i

  NMom = size(p,2)
  do i=1,NMom
      E =p(4,i)
      px=p(1,i)
      py=p(2,i)
      pz=p(3,i)
      p(1,i) = E
      p(2,i) = px
      p(3,i) = py
      p(4,i) = pz
  enddo

END SUBROUTINE





SUBROUTINE swapMomR(p1,p2)
implicit none
real(8) :: p1(1:4),p2(1:4),ptmp(1:4)

  ptmp(1:4) = p2(1:4)
  p2(1:4)   = p1(1:4)
  p1(1:4)   = ptmp(1:4)

END SUBROUTINE



SUBROUTINE swapMomC(p1,p2)
implicit none
complex(8) :: p1(1:4),p2(1:4),ptmp(1:4)

  ptmp(1:4) = p2(1:4)
  p2(1:4)   = p1(1:4)
  p1(1:4)   = ptmp(1:4)

END SUBROUTINE


SUBROUTINE swap(r1,r2)
implicit none
real(8) :: r1,r2,rtmp

  rtmp = r2
  r2   = r1
  r1   = rtmp

END SUBROUTINE



SUBROUTINE pT_order(N,Mom)
implicit none
integer :: N
real(8) :: Mom(1:4,1:N),Mom_Tmp(1:4,1:N),pTList(1:N)
integer :: i,MomOrder(1:N)


    if(N.lt.1) return
    do i=1,N
      pTList(i) = get_PT(Mom(1:4,i))
      MomOrder(i) = i
    enddo

    call BubleSort(N,pTList(1:N),MomOrder(1:N))

    Mom_Tmp(1:4,1:N) = Mom(1:4,1:N)
    do i=1,N
        Mom(1:4,i) = Mom_Tmp(1:4,MomOrder(i))
    enddo

END SUBROUTINE




SUBROUTINE BubleSort(N,X, IY)
IMPLICIT NONE
integer n
real(8) x(1:n)
integer iy(1:n)
real(8) temp
integer i, j, jmax, itemp

      jmax=n-1
      do i=1,n-1
         temp=1d38
         do j=1,jmax
            if(x(j).gt.x(j+1)) cycle
              temp=x(j)
              x(j)=x(j+1)
              x(j+1)=temp
              itemp=iy(j)
              iy(j)=iy(j+1)
              iy(j+1)=itemp
         enddo
         if(temp.eq.1d38) return
         jmax=jmax-1
       enddo

! check the routine
! real(8) :: x(1:10)
! integer :: iy(1:10)
!     x(1:10) = (/3d0,62d0,2d0,78d0,32d0,87d0,1d0,199d0,4d0,73d0/)
!     iy(1:10) = (/1,2,3,4,5,6,7,8,9,10/)
!     print *, x(:)
!     call  BubleSort(10,X(1:10), IY)
!     print *, x(:)
!     print *, iy(:)
!     stop

return
END SUBROUTINE



FUNCTION go_Gauss2(N,matrix)
! same routine but destroyes the original matrix <--> no matrix duplication necessary
implicit none
complex(8) :: matrix(1:N,1:N+1)
integer, intent(in) :: N
complex(8) :: go_Gauss2(1:N)
integer :: i, j
complex(8) :: mult, tsum


    do i = 1,N
      if (cdabs(matrix(i,i)) <= 1.0e-6) then
        call Error("go_Gauss2, Zero pivot element")
      end if
      do j = i+1,N
        mult = matrix(j,i)/matrix(i,i)
        matrix(j,1:N+1) = matrix(j,1:N+1) - mult*matrix(i,1:N+1)
      end do
    end do
    do i=N,1,-1
      tsum=matrix(i,N+1)
      do j=i+1,N
        tsum = tsum - matrix(i,j)*matrix(j,N+1)
      end do
      matrix(i,N+1) = tsum/matrix(i,i)
    end do
    go_Gauss2(1:N) = matrix(1:N,N+1)
    return
END FUNCTION





FUNCTION go_GaussLU(n,A)
implicit none
integer n, perm(1:n)
complex(8) A(1:n,1:n+1)
complex(8) go_GaussLU(n)

      call XLUDecomp(A(1:n,1:n), n, perm(1:n))
      call XLUBackSubst(A(1:n,1:n), n, perm(1:n), A(1:n,n+1))
      go_GaussLU(1:n) = A(1:n,n+1)

END FUNCTION


! * Solution of the linear equation A.x = B by Gaussian elimination
! * with partial pivoting
! * this file is part of LoopTools
! * last modified 24 Jan 06 th
! * Author: Michael Rauch, 7 Dec 2004
! * Reference: Folkmar Bornemann, course notes to
! * Numerische Mathematik 1, Technische Universitaet, Munich, Germany
!
! *#include "defs.h"
!
! *#define MAXDIM 8
!
! ************************************************************************
! * LUDecomp computes the LU decomposition of the n-by-n matrix A
! * by Gaussian Elimination with partial pivoting;
! * compact (in situ) storage scheme
! * Input:
! *   A: n-by-n matrix to LU-decompose
! *   n: dimension of A
! * Output:
! *   A: mangled LU decomposition of A in the form
! *     ( y11 y12 ... y1n )
! *     ( x21 y22 ... y2n )
! *     ( x31 x32 ... y3n )
! *     ( ............... )
! *     ( xn1 xn2 ... ynn )
! *   where
! *     (   1   0 ...   0 )  ( y11 y12 ... y1n )
! *     ( x21   1 ...   0 )  (   0 y22 ... y2n )
! *     ( x31 x32 ...   0 )  (   0   0 ... y3n )  =  Permutation(A)
! *     ( ............... )  ( ............... )
! *     ( xn1 xn2 ...   1 )  (   0   0 ... ynn )
! *   perm: permutation vector


   SUBROUTINE XLUDecomp(A, n, perm)
   implicit none
   integer n, perm(n)
   complex(8) A(n,n)

   integer i, j, k, imax
   complex(8) tmp,czip
   real(8) Amax
   parameter(czip=(0d0,0d0))
   do j = 1, n
! do U part (minus diagonal one)
     do i = 1, j - 1
       do k = 1, i - 1
         A(i,j) = A(i,j) - A(i,k)*A(k,j)
       enddo
     enddo

! do L part (plus diagonal from U case)
     Amax = 0
     do i = j, n
       tmp = czip
       do k = 1, j - 1
         tmp = tmp + A(i,k)*A(k,j)
       enddo
       A(i,j) = A(i,j) - tmp

! do partial pivoting
! find the pivot
       if( abs(A(i,j)) .gt. Amax ) then
         Amax = abs(A(i,j))
         imax = i
       endif
     enddo

! exchange rows
     perm(j) = imax
     do k = 1, n
       tmp = A(j,k)
       A(j,k) = A(imax,k)
       A(imax,k) = tmp
     enddo

! division by the pivot element
     if( A(j,j) .eq. czip ) then
       tmp = dcmplx(1D123)
     else
       tmp = 1/A(j,j)
     endif
     do i = j + 1, n
       A(i,j) = A(i,j)*tmp
     enddo
   enddo
   END SUBROUTINE

! ************************************************************************
! * LUBackSubst computes the x in A.x = b from the LU-decomposed A.
! * Input:
! *   A: LU-decomposed n-by-n matrix A
! *   b: input vector b in A.x = b
! *   n: dimension of A
! *   p: permutation vector from LU decomposition
! * Output:
! *   b: solution vector x in A.x = b

   SUBROUTINE XLUBackSubst(A, n, p, b)
   implicit none
   integer n, p(n)
   complex(8) A(n,n)
   complex(8) b(*)

   integer i, j
   complex(8) tmp

! permute b
   do i = 1, n
     tmp = b(i)
     b(i) = b(p(i))
     b(p(i)) = tmp
   enddo

! forward substitution L.Y = B
   do i = 1, n
     do j = 1, i - 1
       b(i) = b(i) - A(i,j)*b(j)
     enddo
   enddo

! backward substitution U.X = Y
   do i = n, 1, -1
     do j = i + 1, n
       b(i) = b(i) - A(i,j)*b(j)
     enddo
     b(i) = b(i)/A(i,i)
   enddo
   END SUBROUTINE




SUBROUTINE Error(Message,ErrNum)
implicit none
character(*) :: Message
integer,optional :: ErrNum

   if( present(ErrNum) ) then
      print *, "ERROR: ",Message,ErrNum
   else
      print *, "ERROR: ",Message
   endif
   stop
END SUBROUTINE






! Lorentz transformation for the dipole matrix element in top decays (taken from MCFM)
SUBROUTINE WTransform(MomDK,MomDKTd,pbDpg,ptDpg,ptDpb)
implicit none
real(8) :: MomDK(1:4,1:4),MomDKTd(1:4,1:3),pw(4),pt(4),lDt(3:4),lDw(3:4)
real(8) :: root,hsin,hcos,a,b
real(8) :: ptDpt,pwDpw,ptDpw,ptDpg,pbDpg,ptDpb,pWDpl,ptDpl,pWDpn,ptDpn


    pw(1:4) = MomDK(1:4,2) + MomDK(1:4,3)
    pt(1:4) = pw(1:4) + MomDK(1:4,1) + MomDK(1:4,4)

    pbDpg = MomDK(1:4,1).dot.MomDK(1:4,4)
    ptDpg = pt(1:4).dot.MomDK(1:4,4)
    ptDpb = pt(1:4).dot.MomDK(1:4,1)
    ptDpw = pt(1:4).dot.pw(1:4)
    ptDpt = pt(1:4).dot.pt(1:4)
    pwDpw = pw(1:4).dot.pw(1:4)
    pWDpl = pw(1:4).dot.MomDK(1:4,2)
    ptDpl = pt(1:4).dot.MomDK(1:4,2)
    pWDpn = pw(1:4).dot.MomDK(1:4,3)
    ptDpn = pt(1:4).dot.MomDK(1:4,3)

    root=dsqrt(ptDpw**2-ptDpt*pwDpw)
    hsin=0.5d0/(ptDpt*pwDpw)*(-(ptDpt-pwDpw)*ptDpw+(ptDpt+pwDpw)*root)
    hcos=0.5d0/(ptDpt*pwDpw)*(+(ptDpt+pwDpw)*ptDpw-(ptDpt-pwDpw)*root)
    a=hsin/root
    b=(hcos-1d0)/root**2

    MomDKTd(1:4,2) = MomDK(1:4,2) + a*( pt(1:4)*pwDpl-pw(1:4)*ptDpl )  &
                                  + b*( ptDpw*(pt(1:4)*pwDpl+pw(1:4)*ptDpl) - pt(1:4)*pwDpw*ptDpl - pw(1:4)*ptDpt*pwDpl )
    MomDKTd(1:4,3) = MomDK(1:4,3) + a*( pt(1:4)*pwDpn-pw(1:4)*ptDpn )  &
                                  + b*( ptDpw*(pt(1:4)*pwDpn+pw(1:4)*ptDpn) - pt(1:4)*pwDpw*ptDpn - pw(1:4)*ptDpt*pwDpn )
    MomDKTd(1:4,1) = pt(1:4) - MomDKTd(1:4,2) - MomDKTd(1:4,3)

RETURN
END SUBROUTINE




! Lorentz transformation for the dipole matrix element in top decays with additional photon in w decay
SUBROUTINE WTransform2(MomDK,MomDKTd,pbDpg,ptDpg,ptDpb)
implicit none
real(8) :: MomDK(1:4,1:5),MomDKTd(1:4,1:4),pw(4),pt(4),lDt(3:4),lDw(3:4)
real(8) :: root,hsin,hcos,a,b
real(8) :: ptDpt,pwDpw,ptDpw,ptDpg,pbDpg,ptDpb,pWDpl,ptDpl,pWDpn,ptDpn,pWDpp,ptDpp


    pw(1:4) = MomDK(1:4,2) + MomDK(1:4,3) + MomDK(1:4,4)
    pt(1:4) = pw(1:4) + MomDK(1:4,1) + MomDK(1:4,5)

    pbDpg = MomDK(1:4,1).dot.MomDK(1:4,5)
    ptDpg = pt(1:4).dot.MomDK(1:4,5)
    ptDpb = pt(1:4).dot.MomDK(1:4,1)
    ptDpw = pt(1:4).dot.pw(1:4)
    ptDpt = pt(1:4).dot.pt(1:4)
    pwDpw = pw(1:4).dot.pw(1:4)
    pWDpl = pw(1:4).dot.MomDK(1:4,2)
    ptDpl = pt(1:4).dot.MomDK(1:4,2)
    pWDpn = pw(1:4).dot.MomDK(1:4,3)
    ptDpn = pt(1:4).dot.MomDK(1:4,3)
    pWDpp = pw(1:4).dot.MomDK(1:4,4)
    ptDpp = pt(1:4).dot.MomDK(1:4,4)


    root=dsqrt(ptDpw**2-ptDpt*pwDpw)
    hsin=0.5d0/(ptDpt*pwDpw)*(-(ptDpt-pwDpw)*ptDpw+(ptDpt+pwDpw)*root)
    hcos=0.5d0/(ptDpt*pwDpw)*(+(ptDpt+pwDpw)*ptDpw-(ptDpt-pwDpw)*root)
    a=hsin/root
    b=(hcos-1d0)/root**2

    MomDKTd(1:4,2) = MomDK(1:4,2) + a*( pt(1:4)*pwDpl-pw(1:4)*ptDpl )  &
                                  + b*( ptDpw*(pt(1:4)*pwDpl+pw(1:4)*ptDpl) - pt(1:4)*pwDpw*ptDpl - pw(1:4)*ptDpt*pwDpl )
    MomDKTd(1:4,3) = MomDK(1:4,3) + a*( pt(1:4)*pwDpn-pw(1:4)*ptDpn )  &
                                  + b*( ptDpw*(pt(1:4)*pwDpn+pw(1:4)*ptDpn) - pt(1:4)*pwDpw*ptDpn - pw(1:4)*ptDpt*pwDpn )
    MomDKTd(1:4,4) = MomDK(1:4,4) + a*( pt(1:4)*pwDpp-pw(1:4)*ptDpp )  &
                                  + b*( ptDpw*(pt(1:4)*pwDpp+pw(1:4)*ptDpp) - pt(1:4)*pwDpw*ptDpp - pw(1:4)*ptDpt*pwDpp )

    MomDKTd(1:4,1) = pt(1:4) - MomDKTd(1:4,2) - MomDKTd(1:4,3) - MomDKTd(1:4,4)

RETURN
END SUBROUTINE




! Lorentz transformation for the dipole matrix element in top decays with additional photon in w decay
SUBROUTINE DipoleTrafo(MomT,MomK,MomIn,MomOut)
implicit none
real(8) :: MomT(1:4),MomK(1:4),MomOut(1:4),MomIn(1:4)
real(8) :: xh,sinhyp,coshyp,a,b
real(8) :: mt2,mK2,s_tP,s_KP,s_tK


    mt2  = MomT(1:4).dot.MomT(1:4)
    mK2  = MomK(1:4).dot.MomK(1:4)
    s_tP = MomT(1:4).dot.MomIn(1:4)
    s_KP = MomK(1:4).dot.MomIn(1:4)
    s_tK = MomT(1:4).dot.MomK(1:4)

    xh = s_tK**2 - mt2 * mK2
    sinhyp = 0.5d0/(mt2 * mK2)*( -(mt2-mK2)*s_tK + (mt2+mK2)*dsqrt(xh) )
    coshyp = 0.5d0/(mt2 * mK2)*( +(mt2+mK2)*s_tK - (mt2-mK2)*dsqrt(xh) )

    a = sinhyp/dsqrt(xh)
    b = (coshyp-1d0)/xh

    MomOut(1:4) = MomIn(1:4) - a*( s_tP*MomK(1:4) - s_KP*MomT(1:4) )  &
                  + b*( s_tK*( s_tP*MomK(1:4) + s_KP*MomT(1:4) )  -  mK2*s_tP*MomT(1:4) - mt2*s_KP*MomK(1:4))


RETURN
END SUBROUTINE






         subroutine sc_64(n,x,y,r)
         implicit none
         integer i,n
         complex(8) x(*),y(*)
         complex(8) r

            r = x(1)*y(1)
            do i=2, n
              r = r - x(i)*y(i)
            enddo

         return
         end subroutine


          function sc__64(p1,p2)
          implicit none
          complex(8) :: p1(:),p2(:)
          complex(8) :: sc__64
          integer :: sizemin

              sizemin=min(size(p1),size(p2))
              call sc_64(sizemin,p1,p2,sc__64)

          return
          end function


         subroutine sc_128(n,x,y,r)
         implicit none
         integer i,n
         complex(16) x(*),y(*)
         complex(16) r

            r = x(1)*y(1)
            do i=2, n
              r = r - x(i)*y(i)
            enddo

         return
         end subroutine


          function sc__128(p1,p2)
          implicit none
          complex(16) :: p1(:),p2(:)
          complex(16) :: sc__128
          integer :: sizemin

              sizemin=min(size(p1),size(p2))
              call sc_128(sizemin,p1,p2,sc__128)

          return
          end function




         function spi2_(v,sp)
         implicit none
         double complex, intent(in) :: sp(:),v(:)
         double complex :: spi2_(size(sp))
         integer :: Dv,Ds

          Ds = size(sp)
          if (Ds == 4) Dv = 4
          if (Ds == 8) Dv = 6
          if (Ds == 16) Dv = 8
          call spi2(Dv,Ds,v,sp,spi2_)
        return
        end function


        function spb2_(sp,v)
        implicit none
        double complex, intent(in) :: sp(:),v(:)
        double complex :: spb2_(size(sp))
        integer :: Dv,Ds

          Ds = size(sp)
          if (Ds == 4) Dv = 4
          if (Ds == 8) Dv = 6
          if (Ds == 16) Dv = 8
          call spb2(Dv,Ds,sp,v,spb2_)
        return
        end function





      function WeylToDirac(sp)   ! unitary transformation U to convert a Weyl spinor into the Dirac representation
      implicit none              ! sp can be spinor or bar-spinor, i.e. U^dagger.sp = barsp.U
      double complex :: sp(1:4)
      double complex :: WeylToDirac(1:4)
      double precision,parameter :: SqrtFac=1d0/dsqrt(2d0)

          WeylToDirac(1) = SqrtFac*(sp(1)+sp(3))
          WeylToDirac(2) = SqrtFac*(sp(2)+sp(4))
          WeylToDirac(3) = SqrtFac*(sp(1)-sp(3))
          WeylToDirac(4) = SqrtFac*(sp(2)-sp(4))
      return
      end function




      function SpiVL(sp,v)   ! SpiVL=sp.(v_mu*gamma^mu) =spb2(4,4,...)
      implicit none
      double complex :: sp(1:4),v(1:4)
      double complex :: SpiVL(1:4)

        SpiVL(1) = sp(1)*v(1) + sp(4)*(v(2) + (0d0,1d0)*v(3)) + sp(3)*v(4)
        SpiVL(2) = sp(2)*v(1) + sp(3)*(v(2) - (0d0,1d0)*v(3)) - sp(4)*v(4)
        SpiVL(3) =-sp(3)*v(1) - sp(2)*(v(2) + (0d0,1d0)*v(3)) - sp(1)*v(4)
        SpiVL(4) =-sp(4)*v(1) - sp(1)*(v(2) - (0d0,1d0)*v(3)) + sp(2)*v(4)
      return
      end function






!-----------modified procedurepittau.f:            subroutine psp1(Ds,sp1,sp2,r)
!---------- Dv is the dimensionality of the vector space
!---------- Ds is the dimensionality of the spinorial representation
         subroutine spb2(Dv,Ds,sp,v,f)
         implicit none
         integer i,i1,i2,i3,Dv,Ds,imax
         double complex sp(Ds),v(Dv),f(Ds)
         double complex x0(4,4),xx(4,4),xy(4,4)
         double complex xz(4,4),x5(4,4)
         double complex y1,y2,y3,y4,bp,bm,cp,cm

           imax = Ds/4

           do i=1,imax
           i1= 1+4*(i-1)
           i2=i1+3

           y1=sp(i1)
           y2=sp(i1+1)
           y3=sp(i1+2)
           y4=sp(i1+3)

           x0(1,i)=y1
           x0(2,i)=y2
           x0(3,i)=-y3
           x0(4,i)=-y4

           xx(1,i) = -y4
           xx(2,i) = -y3
           xx(3,i) = y2
           xx(4,i) = y1

           xy(1,i)=dcmplx(0d0,-1d0)*y4
           xy(2,i)=dcmplx(0d0,1d0)*y3
           xy(3,i)=dcmplx(0d0,1d0)*y2
           xy(4,i)=dcmplx(0d0,-1d0)*y1

           xz(1,i)=-y3
           xz(2,i)=y4
           xz(3,i)=y1
           xz(4,i)=-y2

           x5(1,i)=y3
           x5(2,i)=y4
           x5(3,i)=y1
           x5(4,i)=y2

           enddo

           if (Dv.eq.4) then

           do i=1,4

           f(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)

           enddo

           endif

           if (Dv.eq.6) then
           bp = (v(5)+dcmplx(0d0,1d0)*v(6))
           bm=(v(5)-dcmplx(0d0,1d0)*v(6))

           do i=1,4

           f(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)+bm*x5(i,2)

            i1 = i+4

            f(i1)= v(1)*x0(i,2)-v(2)*xx(i,2)-v(3)*xy(i,2)-v(4)*xz(i,2)-bp*x5(i,1)


            enddo

           endif

           if (Dv.eq.8) then
           bp=(v(5)+dcmplx(0d0,1d0)*v(6))
           bm=(v(5)-dcmplx(0d0,1d0)*v(6))
           cp=(v(7)+dcmplx(0d0,1d0)*v(8))
           cm=(v(7)-dcmplx(0d0,1d0)*v(8))

           do i=1,4

           f(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)+bm*x5(i,2)-cm*x5(i,3)

            i1 = i+4

            f(i1)= v(1)*x0(i,2)-v(2)*xx(i,2)-v(3)*xy(i,2)-v(4)*xz(i,2)-bp*x5(i,1)+cm*x5(i,4)

             i2 = i1+4

             f(i2)=v(1)*x0(i,3)-v(2)*xx(i,3)-v(3)*xy(i,3)-v(4)*xz(i,3)+bm*x5(i,4)+cp*x5(i,1)

              i3=i2+4

              f(i3)=v(1)*x0(i,4)-v(2)*xx(i,4)-v(3)*xy(i,4)-v(4)*xz(i,4)-bp*x5(i,3)-cp*x5(i,2)

              enddo

              endif

               return
               end SUBROUTINE











!-----------modified procedure
!---------- Dv is the dimensionality of the vector space
!---------- Ds is the dimensionality of the spinorial representation
         subroutine spi2(Dv,Ds,v,sp,f)
         implicit none
         integer i,i1,i2,i3,imax,Dv,Ds
         double complex sp(Ds),v(Dv),f(Ds)
         double complex x0(4,4),xx(4,4),xy(4,4)
         double complex xz(4,4),x5(4,4)
         double complex y1,y2,y3,y4,bp,bm,cp,cm

         imax = Ds/4

           do i=1,imax
           i1= 1+4*(i-1)
           i2=i1+3

           y1=sp(i1)
           y2=sp(i1+1)
           y3=sp(i1+2)
           y4=sp(i1+3)

           x0(1,i)=y1
           x0(2,i)=y2
           x0(3,i)=-y3
           x0(4,i)=-y4


           xx(1,i) = y4
           xx(2,i) = y3
           xx(3,i) = -y2
           xx(4,i) = -y1


           xy(1,i)=dcmplx(0d0,-1d0)*y4
           xy(2,i)=dcmplx(0d0,1d0)*y3
           xy(3,i)=dcmplx(0d0,1d0)*y2
           xy(4,i)=dcmplx(0d0,-1d0)*y1

           xz(1,i)=y3
           xz(2,i)=-y4
           xz(3,i)=-y1
           xz(4,i)=y2

           x5(1,i)=y3
           x5(2,i)=y4
           x5(3,i)=y1
           x5(4,i)=y2

           enddo

           if(Dv.eq.4) then

           do i=1,4

           f(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)
           enddo

           endif


          if (Dv.eq.6) then
           bp = (v(5)+dcmplx(0d0,1d0)*v(6))
           bm=(v(5)-dcmplx(0d0,1d0)*v(6))


           do i=1,4

           f(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)-bp*x5(i,2)

            i1=i+4

            f(i1)=v(1)*x0(i,2)-v(2)*xx(i,2)-v(3)*xy(i,2)-v(4)*xz(i,2)+bm*x5(i,1)

            enddo

          endif

          if (Dv.eq.8) then

           bp = (v(5)+dcmplx(0d0,1d0)*v(6))
           bm=(v(5)-dcmplx(0d0,1d0)*v(6))
           cp=(v(7)+dcmplx(0d0,1d0)*v(8))
           cm=(v(7)-dcmplx(0d0,1d0)*v(8))



           do i=1,4

           f(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)-bp*x5(i,2)+ cp*x5(i,3)

            i1=i+4

            f(i1)=v(1)*x0(i,2)-v(2)*xx(i,2)-v(3)*xy(i,2)-v(4)*xz(i,2)+bm*x5(i,1)-cp*x5(i,4)

             i2=i1+4

            f(i2)=v(1)*x0(i,3)-v(2)*xx(i,3)-v(3)*xy(i,3)-v(4)*xz(i,3)-bp*x5(i,4)-cm*x5(i,1)

             i3=i2+4

             f(i3)=v(1)*x0(i,4)-v(2)*xx(i,4)-v(3)*xy(i,4)-v(4)*xz(i,4)+bm*x5(i,3)+cm*x5(i,2)


            enddo

            endif

           return
           end subroutine





!-----------``scalar product'' for the two spinors
            subroutine psp1(Ds,sp1,sp2,r)
            implicit none
            integer i,Ds
            double complex sp1(Ds),sp2(Ds)
            double complex r

            r=dcmplx(0d0,0d0)

               do i=1,Ds
                 r = r+sp1(i)*sp2(i)
               enddo
            return
            end subroutine




          function  psp1_(sp1,sp2)
          implicit none
          double complex, intent(in) :: sp1(:)
          double complex, intent(in) :: sp2(:)
          double complex :: psp1_

            psp1_ = sum(sp1(1:)*sp2(1:))

            end function psp1_









      function vggg(e1,k1,e2,k2)
      implicit none
      complex(8), intent(in) :: e1(:), e2(:)
      complex(8), intent(in) :: k1(:), k2(:)
      complex(8)             :: vggg(size(e1))
      complex(8):: sk1e2,se1e2,sk2e1,xx
      real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0


          sk1e2=sc_(k1,e2)
          sk2e1=sc_(k2,e1)
          se1e2=sc_(e1,e2)
          xx=(0.0d0,1.0d0)*sqrt2
          vggg = xx*(-sk1e2*e1+sk2e1*e2+se1e2/2d0*(k1-k2))

       end function vggg


       function  vgggg(e1,e2,e3)
       implicit none
       complex(8), intent(in) :: e1(:),e2(:),e3(:)
       complex(8)             :: vgggg(size(e1))
       complex(8):: se1e3,se2e3,se1e2

          se1e3=sc_(e1,e3)
          se2e3=sc_(e2,e3)
          se1e2=sc_(e1,e2)
          vgggg = (0.0d0,1.0d0)*(e2*se1e3-0.5d0*(e1*se2e3+ e3*se1e2))

       end function vgggg


      function vqg(sp,e1)
      implicit none
      complex(8), intent(in) :: e1(:)
      complex(8), intent(in) :: sp(:)
      complex(8) :: vqg(size(sp))
      real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

          vqg = (0d0,1d0)/sqrt2*spb2_(sp,e1)

      end function vqg


      function vqg_Weyl(sp,e1)
      implicit none
      complex(8), intent(in) :: e1(:)
      complex(8), intent(in) :: sp(:)
      complex(8) :: vqg_Weyl(size(sp))
      real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

          vqg_Weyl = (0d0,1d0)/sqrt2*spb2_Weyl(sp,e1)

      end function vqg_Weyl





      function vgq(e1,sp)
      implicit none
      complex(8), intent(in) :: e1(:)
      complex(8), intent(in) :: sp(:)
      complex(8) :: vgq(size(sp))
      real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

          vgq = (0d0,-1d0)/sqrt2*spb2_(sp,e1)

      end function vgq





      function vbqg(sp,e1)
      implicit none
      complex(8), intent(in) :: e1(:)
      complex(8), intent(in) :: sp(:)
      complex(8) :: vbqg(size(sp))
      real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

          vbqg = (0d0,-1d0)/sqrt2*spi2_(e1,sp)

      end function vbqg



      function vgbq(e1,sp)
      implicit none
      complex(8), intent(in) :: e1(:)
      complex(8), intent(in) :: sp(:)
      complex(8) :: vgbq(size(sp))
      real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

          vgbq = (0d0,1d0)/sqrt2*spi2_(e1,sp)

      end function vgbq



      function vgbq_Weyl(e1,sp)
      implicit none
      complex(8), intent(in) :: e1(:)
      complex(8), intent(in) :: sp(:)
      complex(8) :: vgbq_Weyl(size(sp))
      real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

          vgbq_Weyl = (0d0,1d0)/sqrt2*spi2_Weyl(e1,sp)

      end function vgbq_Weyl


      function vbqq(Dv,sp1,sp2)
      implicit none
      complex(8), intent(in) :: sp1(:), sp2(:)
      integer, intent(in) ::  Dv
      integer :: i
      complex(8) :: vbqq(Dv)
      complex(8) :: rr, va(Dv),sp1a(size(sp1))
      real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

          va=(0d0,0d0)
          vbqq=(0d0,0d0)

          do i=1,Dv
             if (i.eq.1) then
               va(1)=(1d0,0d0)
             else
               va(i)=(-1d0,0d0)
             endif
             sp1a=spb2_(sp1,va)

             rr=(0d0,-1d0)/sqrt2*psp1_(sp1a,sp2)
             if (i.eq.1) then
                  vbqq = vbqq + rr*va
              else
                  vbqq = vbqq - rr*va
             endif
             va(i)=(0d0,0d0)
          enddo

      end function vbqq





      function vbqq2(sp1,sp2)! this is my own simpler 4-dim version,   without -i/sqrt(2) norm!!
      implicit none
      complex(8), intent(in) :: sp1(:), sp2(:)
      integer :: i
      complex(8) :: vbqq2(4)
      complex(8) :: sp1a(4)
      real(8) :: va(1:4,1:4)

         va(1,1:4)=(/+1d0,0d0,0d0,0d0/)      
         va(2,1:4)=(/0d0,-1d0,0d0,0d0/)      
         va(3,1:4)=(/0d0,0d0,-1d0,0d0/)      
         va(4,1:4)=(/0d0,0d0,0d0,-1d0/)

          do i=1,4
             sp1a=spb2_(sp1,dcmplx(va(i,1:4)))
             vbqq2(i) = psp1_(sp1a,sp2)
          enddo

      end function vbqq2




      function vvbqq(sp1,sp2)
      implicit none
      complex(8), intent(in) :: sp1(:), sp2(:)
      integer :: i,j
      complex(8) :: vvbqq(4,4)
      complex(8) :: sp1a(4)
      real(8) :: va(1:4,1:4)

         va(1,1:4)=(/+1d0,0d0,0d0,0d0/)      
         va(2,1:4)=(/0d0,-1d0,0d0,0d0/)      
         va(3,1:4)=(/0d0,0d0,-1d0,0d0/)      
         va(4,1:4)=(/0d0,0d0,0d0,-1d0/)

         do i=1,4
          do j=1,4
             sp1a=spb2_(sp1, dcmplx(va(i,1:4)))
             sp1a=spb2_(sp1a,dcmplx(va(j,1:4)))
             vvbqq(i,j) = psp1_(sp1a,sp2)
          enddo
         enddo

      end function vvbqq




      function vvvbqq(sp1,sp2)
      implicit none
      complex(8), intent(in) :: sp1(:), sp2(:)
      integer :: i,j,k
      complex(8) :: vvvbqq(4,4,4)
      complex(8) :: sp1a(4)
      real(8) :: va(1:4,1:4)

         va(1,1:4)=(/+1d0,0d0,0d0,0d0/)      
         va(2,1:4)=(/0d0,-1d0,0d0,0d0/)      
         va(3,1:4)=(/0d0,0d0,-1d0,0d0/)      
         va(4,1:4)=(/0d0,0d0,0d0,-1d0/)

         do i=1,4
          do j=1,4
           do k=1,4
              sp1a=spb2_(sp1, dcmplx(va(i,1:4)))
              sp1a=spb2_(sp1a,dcmplx(va(j,1:4)))
              sp1a=spb2_(sp1a,dcmplx(va(k,1:4)))
              vvvbqq(i,j,k) = psp1_(sp1a,sp2)
           enddo
          enddo
         enddo

      end function vvvbqq







      function vbqW(sp,e1)
      implicit none
      complex(8), intent(in) :: e1(:)
      complex(8), intent(in) :: sp(:)
      complex(8) :: vbqW(size(sp))

            vbqW = -(0d0,1d0)*Chir(.false.,spb2_(sp,e1) )

      end function vbqW


      function vWq(e1,sp)
      implicit none
      complex(8), intent(in) :: e1(:)
      complex(8), intent(in) :: sp(:)
      complex(8) :: vWq(size(sp))

            vWq = -(0d0,1d0)*Chir(.true., spi2_(e1,sp) )

      end function vWq




      function vbqW_Weyl(sp,e1)
      implicit none
      complex(8), intent(in) :: e1(:)
      complex(8), intent(in) :: sp(:)
      complex(8) :: vbqW_Weyl(size(sp))

            vbqW_Weyl = -(0d0,1d0)*spb2_Weyl(sp,e1)

      end function vbqW_Weyl


      function vWq_Weyl(e1,sp)
      implicit none
      complex(8), intent(in) :: e1(:)
      complex(8), intent(in) :: sp(:)
      complex(8) :: vWq_Weyl(size(sp))

            vWq_Weyl = -(0d0,1d0)*spi2_Weyl(e1,sp)

      end function vWq_Weyl






!-------------- color charged scalar couplings

      function csg(e1,k1,p)!! are all momenta really outgoing??
      implicit none
      complex(8), intent(in) :: e1(:)
      complex(8), intent(in) :: k1(:)
      complex(8), intent(in) ::  p(:)
      complex(8) :: csg
      real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

          csg = (0d0,+1d0)/sqrt2*( sc_(k1,e1)+2d0*sc_(p,e1) )!  minus sign corrected

      end function csg

      function cgs(e1,k1,p)
      implicit none
      complex(8), intent(in) :: e1(:)
      complex(8), intent(in) :: k1(:)
      complex(8), intent(in) ::  p(:)
      complex(8) :: cgs
      real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

          cgs = (0d0,-1d0)/sqrt2*( sc_(k1,e1)+2d0*sc_(p,e1) )

      end function cgs





      function cbsg(e1,k1,p)!! are all momenta really outgoing??
      implicit none
      complex(8), intent(in) :: e1(:)
      complex(8), intent(in) :: k1(:)
      complex(8), intent(in) ::  p(:)
      complex(8) :: cbsg
      real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

          cbsg = (0d0,1d0)/sqrt2*( sc_(k1,e1)+2d0*sc_(p,e1) )

      end function cbsg


      function cgbs(e1,k1,p)
      implicit none
      complex(8), intent(in) :: e1(:)
      complex(8), intent(in) :: k1(:)
      complex(8), intent(in) ::  p(:)
      complex(8) :: cgbs
      real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0
 
          cgbs = (0d0,-1d0)/sqrt2*( sc_(k1,e1)+2d0*sc_(p,e1) )

      end function cgbs


      function cggs(e1,e2)
      implicit none
      complex(8), intent(in) :: e1(:),e2(:)
      complex(8) :: cggs

          cggs = (0d0,1d0)/2d0 * sc_(e1,e2)!  i/2 corrected

      end function cggs


      function csgg(e1,e2)
      implicit none
      complex(8), intent(in) :: e1(:),e2(:)
      complex(8) :: csgg

          csgg =  (0d0,1d0)/2d0 * sc_(e1,e2)!  i/2 corrected

      end function csgg


      function cgsg(e1,e2)
      implicit none
      complex(8), intent(in) :: e1(:),e2(:)
      complex(8) :: cgsg

          cgsg = -(0d0,1d0)/2d0 * sc_(e1,e2)

      end function cgsg





      function vbss(Dv,p1,p2)
      implicit none
      complex(8), intent(in) :: p1(:),p2(:)
      integer :: Dv
      complex(8) :: vbss(1:Dv)

          vbss(:) = -(0d0,1d0)/dsqrt(2d0) * ( p1(:)-p2(:) )! corrected: sqrt(2)

      end function vbss


      function vsbs(Dv,p1,p2)
      implicit none
      complex(8), intent(in) :: p1(:),p2(:)
      integer :: Dv
      complex(8) :: vsbs(1:Dv)

          vsbs(:) = +(0d0,1d0)/dsqrt(2d0) * ( p2(:)-p1(:) )! corrected: sqrt(2)

      end function vsbs



      function vggss(Dv,e1)
      implicit none
      complex(8), intent(in) :: e1(:)
      integer :: Dv
      complex(8) :: vggss(1:Dv)

          vggss(:) = +(0d0,1d0)/2d0 * e1(:)

      end function vggss


!-------------- END: color charged scalar couplings







          subroutine ubarSpi(p,m,i,f)
          implicit none
          integer i
          real(8) m
          complex(8) p(4)
          complex(8) f(4),fc
          real(8)  p0,px,py,pz,fc2

          p0=dreal(p(1))
          px=dreal(p(2))
          py=dreal(p(3))
          pz=dreal(p(4))

          fc2=p0+m
          fc=cdsqrt( dcmplx(fc2))
!           fc=dsqrt(fc2)

          if (i.eq.1) then
            f(1)=fc
            f(2)=dcmplx(0d0,0d0)
            f(3)=-1d0*pz*fc/fc2
            f(4)=-(px-(0d0,1d0)*py)*fc/fc2
          elseif (i.eq.-1) then
            f(1)=dcmplx(0d0,0d0)
            f(2)=fc
            f(3)=-(px+(0d0,1d0)*py)*fc/fc2
            f(4)=pz*fc/fc2
          else
              call Error("wrong helicity setting in ubarSpi")
          endif

          return
          end subroutine





          subroutine uSpi(p,m,i,f)
          implicit none
          integer i
          real(8) m
          complex(8) p(4)
          complex(8) f(4),fc
          real(8)  p0,px,py,pz,fc2

          p0=dreal(p(1))
          px=dreal(p(2))
          py=dreal(p(3))
          pz=dreal(p(4))

          fc2=p0+m
          fc=cdsqrt( dcmplx(fc2))
!           fc=dsqrt(fc2)

          if (i.eq.1) then
            f(1)=fc
            f(2)=dcmplx(0d0,0d0)
            f(3)=pz*fc/fc2
            f(4)=(px+(0d0,1d0)*py)*fc/fc2
          elseif (i.eq.-1) then
            f(1)=dcmplx(0d0,0d0)
            f(2)=fc
            f(3)=(px-(0d0,1d0)*py)*fc/fc2
            f(4)=-pz*fc/fc2
          else
              call Error("wrong helicity setting in uSpi")
          endif

          return
          end subroutine






          subroutine vSpi(p,m,i,f)
          implicit none
          integer i
          real(8) m
          complex(8) p(4)
          complex(8) f(4),fc
          real(8) p0,px,py,pz,fc2

          p0=dreal(p(1))
          px=dreal(p(2))
          py=dreal(p(3))
          pz=dreal(p(4))

          fc2 = p0+m
          fc=cdsqrt(dcmplx(fc2))
!           fc=dsqrt(fc2)

          if (i.eq.1) then
            f(1)=pz*fc/fc2
            f(2)=(px+(0d0,1d0)*py)*fc/fc2
            f(3)=fc
            f(4)=dcmplx(0d0,0d0)
          elseif (i.eq.-1) then
            f(1)=(px-(0d0,1d0)*py)*fc/fc2
            f(2)=-pz*fc/fc2
            f(3)=dcmplx(0d0,0d0)
            f(4)=fc
          else
              call Error("wrong helicity setting in vSpi")
          endif

          return
          end SUBROUTINE





          subroutine vbarSpi(p,m,i,f)
          implicit none
          integer i
          real(8) m
          complex(8) p(4)
          complex(8) f(4),fc
          real(8) p0,px,py,pz,fc2

          p0=dreal(p(1))
          px=dreal(p(2))
          py=dreal(p(3))
          pz=dreal(p(4))

          fc2 = p0+m
          fc=cdsqrt(dcmplx(fc2))
!           fc=dsqrt(fc2)


          if (i.eq.1) then
            f(1)=pz*fc/fc2
            f(2)=(px-(0d0,1d0)*py)*fc/fc2
            f(3)=-fc
            f(4)=dcmplx(0d0,0d0)
          elseif (i.eq.-1) then
            f(1)=(px+(0d0,1d0)*py)*fc/fc2
            f(2)=-pz*fc/fc2
            f(3)=dcmplx(0d0,0d0)
            f(4)=-fc
          else
              call Error("wrong helicity setting in vbarSpi")
          endif
          return
          end SUBROUTINE






!-------massive vector boson polarization routine

      function pol_mass(p,m,i)
      implicit none
      integer, intent(in) :: i
      integer :: pol
      complex(8), intent(in) :: p(4)
      complex(8) :: pol_mass(4)
      real(8),  intent(in) :: m
      real(8) :: p0,px,py,pz, pv
      real(8) :: ct,st,cphi,sphi

          p0=dreal(p(1))
          px=dreal(p(2))
          py=dreal(p(3))
          pz=dreal(p(4))

          pv= dsqrt(dabs(p0**2 - m**2))
          if(pv/m.lt.1d-8) then
                if(i.eq.0) then
                    pol_mass(1:4)=(0d0,0d0)
                    return
                endif
                ct = 1d0; st=0d0
          else
                ct= pz/pv
                st= dsqrt(dabs(1.0d0-ct**2))
          endif


          if (st .lt. 1D-15) then
              cphi=1.0d0
              sphi=0.0d0
          else
              cphi= px/pv/st
              sphi= py/pv/st
          endif


!         i=0 is longitudinal polarization
!         the following ifstatement distinguishes between
!         positive and negative energies
          if ( p0 .gt. 0.0d0) then
          pol=i
          else
          pol=-i
          endif

          if(pol .eq. -1.or.pol .eq. 1) then
              pol_mass(1)=dcmplx(0.0d0,0.0d0)
              pol_mass(2)=dcmplx(ct*cphi/dsqrt(2d0),-pol*sphi/dsqrt(2d0))
              pol_mass(3)=dcmplx(ct*sphi/dsqrt(2d0), pol*cphi/dsqrt(2d0))
              pol_mass(4)=dcmplx(-st/dsqrt(2d0),0.0d0)
          elseif (pol .eq. 0) then
              pol_mass(1)= dcmplx(pv/m,0.0d0)
              pol_mass(2)= dcmplx(p0/m/pv*px,0.0d0)
              pol_mass(3)= dcmplx(p0/m/pv*py,0.0d0)
              pol_mass(4)= dcmplx(p0/m/pv*pz,0.0d0)
          else
              call Error("wrong helicity setting in pol_mass")
          endif

        end function pol_mass




      subroutine pol_massSR(p,m,i,pol_mass)
      implicit none
      integer, intent(in) :: i
      integer :: pol
      complex(8), intent(in) :: p(4)
      complex(8) :: pol_mass(4)
      real(8),  intent(in) :: m
      real(8) :: p0,px,py,pz, pv
      real(8) :: ct,st,cphi,sphi

          p0=dreal(p(1))
          px=dreal(p(2))
          py=dreal(p(3))
          pz=dreal(p(4))

          pv= dsqrt(dabs(p0**2 - m**2))
          ct= pz/pv
          st= dsqrt(dabs(1.0d0-ct**2))

          if (st .lt. 1D-15) then
              cphi=1.0d0
              sphi=0.0d0
          else
              cphi= px/pv/st
              sphi= py/pv/st
          endif


!         i=0 is longitudinal polarization
!         the following ifstatement distinguishes between
!         positive and negative energies
          if ( p0 .gt. 0.0d0) then
          pol=i
          else
          pol=-i
          endif

      if(pol .eq. -1.or.pol .eq. 1) then
            pol_mass(1)=dcmplx(0.0d0,0.0d0)
            pol_mass(2)=dcmplx(ct*cphi/dsqrt(2d0),-pol*sphi/dsqrt(2d0))
            pol_mass(3)=dcmplx(ct*sphi/dsqrt(2d0), pol*cphi/dsqrt(2d0))
            pol_mass(4)=dcmplx(-st/dsqrt(2d0),0.0d0)
      elseif (pol == 0) then
            pol_mass(1)= dcmplx(pv/m,0.0d0)
            pol_mass(2)= dcmplx(p0/m/pv*px,0.0d0)
            pol_mass(3)= dcmplx(p0/m/pv*py,0.0d0)
            pol_mass(4)= dcmplx(p0/m/pv*pz,0.0d0)
      else
              call Error("wrong helicity setting in pol_mass")
      endif

        end subroutine pol_massSR






!         massless vector polarization  subroutine
          subroutine pol_mless(p,i,f)
          implicit none
          integer i,pol
          complex(8) p(4)
          real(8) p0,px,py,pz,pv,ct,st,cphi,sphi
          complex(8) f(4),phase

          p0=dreal(p(1))
          px=dreal(p(2))
          py=dreal(p(3))
          pz=dreal(p(4))

!!        abs(3-vec)
          pv=dsqrt(dabs(p0**2))
!!        cos,sin theta
          ct=pz/pv
          st=dsqrt(dabs(1d0-ct**2))

!!        cos,sin phi   (! only pos. values !)
          if (st.lt.1d-8) then
               cphi=1d0
               sphi=0d0
          else
               cphi= px/pv/st
               sphi= py/pv/st
          endif


!         the following ifstatement distinguishes between
!         positive and negative energies
          if ( p0.gt.0d0) then
            pol=i
          else
            pol=-i
          endif

!           phase = -pol*(cphi + (0d0,1d0)*sphi)

          f(1)=dcmplx(0d0,0d0)
          f(2)=dcmplx(ct*cphi/dsqrt(2d0),-pol*sphi/dsqrt(2d0)) !*phase
          f(3)=dcmplx(ct*sphi/dsqrt(2d0), pol*cphi/dsqrt(2d0)) !*phase
          f(4)=dcmplx(-st/dsqrt(2d0),0d0)                      !*phase

            if( abs(pol).ne.1) call Error("wrong helicity setting in pol_mless")

          return
          END SUBROUTINE








! !        procedure to get parts of a 16-dim aray
!          subroutine getparts(x,i1,i2,y)
!          implicit none
!          integer i1,i2,i3,i
!          complex(8) x(16),y(4)
!
!            do i=i1,i2
!            i3=i-i1+1
!            y(i3) = x(i)
!            enddo
!
!           return
!           end SUBROUTINE








!           Spinorial eigenstates for on-shell momenta
             subroutine give_usp_old(Nf,Dv,Ds,p,m,r)
             implicit none
             integer i,j,Nf,Dv,Ds
             real(8) m
             complex(8) p(Dv),u1(Ds),r(Ds,8)
             complex(8) f(Ds)
             complex(8) fc
!              include './misc/ucom_import'

!------------massive case
             if (dabs(m).gt.1D-8) then
                  fc=1D0/cdsqrt(p(1)+dcmplx(m,0d0))
                  do i=1,Nf
                      do j=1,Ds
                          u1(j)=u(i,j)
                      enddo
                      call spi2(Dv,Ds,p,u1,f)
                      do j = 1,Ds
                          r(j,i)=(f(j)+dcmplx(m,0d0)*u1(j))*fc
                      enddo
                  enddo
             endif

!-------------massless case
              if (dabs(m).lt.1D-8) then

                  do i=1,Nf
      !               if (zabs(p(1)-p(4)).gt.1D-5) then
                      if (zabs(p(1)-p(4)).gt.zabs(p(1)-p(2))) then
                          fc=1D0/cdsqrt(2d0*(p(1)-p(4)))
                          do j=1,Ds
                              u1(j)=uz(i,j)
                          enddo
                      else
      !                   if (zabs(p(1)-p(2)).gt.1D-5) then
                          fc=1D0/zsqrt(2d0*(p(1)-p(2)))
                            do j=1,Ds
                                u1(j)=ux(i,j)
                            enddo
      !                   else
                              if (zabs(p(1)-p(2)).lt.1D-5) then
                                   print *, 'MASSLESS FERMION FAILED in give_usp',zabs(p(1)-p(4)),zabs(p(1)-p(2)),p(:)
!                                     pause
                              endif
      !                   endif
                      endif
                      call spi2(Dv,Ds,p,u1,f)
                      do j = 1,Ds
                          r(j,i)=f(j)*fc
                      enddo
                  enddo
              endif


             return
             end subroutine






             subroutine give_barusp_old(Nf,Dv,Ds,p,m,r)
             implicit none
             integer i,j,Nf,Dv,Ds
             real(8) m
             real(8) hc(16)
             complex(8) p(Dv),u1(Ds),r(Ds,8)
             complex(8) f(Ds)
             complex(8) fc
             data hc(1)/ 1d0/, hc(2)/ 1d0/,hc(3)/ -1d0/, hc(4)/ -1d0/, hc(5)/ 1d0/, hc(6)/ 1d0/, &
                  hc(7)/ -1d0/, hc(8)/ -1d0/, hc(9)/ 1d0/, hc(10)/ 1d0/,hc(11)/ -1d0/, &
                  hc(12)/ -1d0/, hc(13)/ 1d0/, hc(14)/ 1d0/,hc(15)/ -1d0/, hc(16)/ -1d0/
!              include './misc/ucom_import'


!-----------massive case
             if (dabs(m).gt.1D-8) then
                fc=1d0/zsqrt(p(1)+dcmplx(m,0d0))
                do i=1,Nf
                    do j=1,Ds
                        u1(j)=u(i,j)
                    enddo
                    call spb2(Dv,Ds,u1,p,f)
                    do j = 1,Ds
                        r(j,i)=(f(j)+dcmplx(m,0d0)*u1(j))*fc
                    enddo
                enddo
             endif

!-----------massless case
             if (dabs(m).lt.1D-8) then

                do i=1,Nf
!                         if (zabs(p(1)-p(4)).gt.1D-5) then
                        if (zabs(p(1)-p(4)).gt.zabs(p(1)-p(2))) then
                            fc=1D0/zsqrt(2d0*(p(1)-p(4)))
                            do j=1,Ds
                                u1(j)=hc(j)*uz(i,j)
                            enddo
                        else
!                             if (zabs(p(1)-p(2)).gt.1D-5) then
                                  fc=1D0/zsqrt(2d0*(p(1)-p(2)))
                                  do j=1,Ds
                                      u1(j)=hc(j)*ux(i,j)
                                  enddo
!                             else
                              if (zabs(p(1)-p(2)).lt.1D-5)  then
                                    print *, 'MASSLESS FERMION FAILED in give_barusp',zabs(p(1)-p(4)),zabs(p(1)-p(2)),p(:)
!                                     pause
                              endif
!                             endif
                        endif

                        call spb2(Dv,Ds,u1,p,f)
                        do j = 1,Ds
                            r(j,i)=f(j)*fc
                        enddo
                enddo

             endif

             return
             end SUBROUTINE



!           Spinorial eigenstates for on-shell momenta
             subroutine give_usp(Nf,Dv,Ds,p,m,r)
             implicit none
             integer i,j,Nf,Dv,Ds, case
             real(8) m, rmax
             complex(8) p(Dv),u1(Ds),r(Ds,8)
             complex(8) f(Ds)
             complex(8) fc
!              include './misc/ucom_import'



!------------massive case
             if (dabs(m).gt.1D-8) then
                  fc=1D0/cdsqrt(p(1)+dcmplx(m,0d0))
                  do i=1,Nf
                      do j=1,Ds
                          u1(j)=u(i,j)
                      enddo
                      call spi2(Dv,Ds,p,u1,f)
                      do j = 1,Ds
                          r(j,i)=(f(j)+dcmplx(m,0d0)*u1(j))*fc
                      enddo
                  enddo
             endif

!-------------massless case
              if (dabs(m).lt.1D-8) then


            if (zabs(p(1)-p(4)).gt.zabs(p(1)-p(2))) then
                  case = 1
                  rmax = zabs(p(1)-p(4))
            else
                  rmax = zabs(p(1)-p(2))
                  case = 2
            endif

            if (rmax.gt.zabs(p(1)-p(3))) then
               case = case
            else
                 case = 3
            endif


                  do i=1,Nf
                      if (case.eq.1) then
                          fc=1D0/cdsqrt(2d0*(p(1)-p(4)))
                          do j=1,Ds
                              u1(j)=uz(i,j)
                          enddo
                      elseif(case.eq.2) then
                          fc=1D0/zsqrt(2d0*(p(1)-p(2)))
                            do j=1,Ds
                                u1(j)=ux(i,j)
                             enddo
                      elseif(case.eq.3) then
                          fc=1D0/zsqrt(2d0*(p(1)-p(3)))
                            do j=1,Ds
                                u1(j)=uy(i,j)
                             enddo
                      endif

                      call spi2(Dv,Ds,p,u1,f)
                      do j = 1,Ds
                          r(j,i)=f(j)*fc
                      enddo
                  enddo


              endif


             return
             end subroutine






             subroutine give_barusp(Nf,Dv,Ds,p,m,r)
             implicit none
             integer i,j,Nf,Dv,Ds, case
             real(8) m
             real(8) hc(16), rmax
             complex(8) p(Dv),u1(Ds),r(Ds,8)
             complex(8) f(Ds)
             complex(8) fc
             data hc(1)/ 1d0/, hc(2)/ 1d0/,hc(3)/ -1d0/, hc(4)/ -1d0/, hc(5)/ 1d0/, hc(6)/ 1d0/, &
                  hc(7)/ -1d0/, hc(8)/ -1d0/, hc(9)/ 1d0/, hc(10)/ 1d0/,hc(11)/ -1d0/, &
                  hc(12)/ -1d0/, hc(13)/ 1d0/, hc(14)/ 1d0/,hc(15)/ -1d0/, hc(16)/ -1d0/
!              include './misc/ucom_import'


!-----------massive case
             if (dabs(m).gt.1D-8) then
                fc=1d0/zsqrt(p(1)+dcmplx(m,0d0))
                do i=1,Nf
                    do j=1,Ds
                        u1(j)=u(i,j)
                    enddo
                    call spb2(Dv,Ds,u1,p,f)
                    do j = 1,Ds
                        r(j,i)=(f(j)+dcmplx(m,0d0)*u1(j))*fc
                    enddo
                enddo
             endif

!-----------massless case

             if (dabs(m).lt.1D-8) then

            if (zabs(p(1)-p(4)).gt.zabs(p(1)-p(2))) then
                  case = 1
                  rmax = zabs(p(1)-p(4))
            else
                  rmax = zabs(p(1)-p(2))
                  case = 2
            endif

            if (rmax.gt.zabs(p(1)-p(3))) then
               case = case
            else
                 case = 3
            endif

                do i=1,Nf
                        if (case.eq.1) then
                            fc=1D0/zsqrt(2d0*(p(1)-p(4)))
                            do j=1,Ds
                                u1(j)=hc(j)*uz(i,j)
                            enddo
                        elseif(case.eq.2) then
                                  fc=1D0/zsqrt(2d0*(p(1)-p(2)))
                                  do j=1,Ds
                                      u1(j)=hc(j)*ux(i,j)
                                  enddo
                        elseif(case.eq.3) then
                                  fc=1D0/zsqrt(2d0*(p(1)-p(3)))
                                  do j=1,Ds
                                      u1(j)=hc(j)*uyc(i,j)
                                  enddo
                        endif
                        call spb2(Dv,Ds,u1,p,f)
                        do j = 1,Ds
                            r(j,i)=f(j)*fc
                        enddo
                enddo


             endif


             return
             end SUBROUTINE





      function Chir(sign,sp)   ! Chir = sp.omega_sign = omega_sign.sp
      implicit none
      logical :: sign
      double complex :: sp(1:4)
      double complex :: Chir(1:4)

        if(sign) then !omega_+
          Chir(1) = 0.5d0*(sp(1)+sp(3))
          Chir(2) = 0.5d0*(sp(2)+sp(4))
          Chir(3) = Chir(1)
          Chir(4) = Chir(2)
        else !omega_-
          Chir(1) = 0.5d0*(sp(1)-sp(3))
          Chir(2) = 0.5d0*(sp(2)-sp(4))
          Chir(3) =-Chir(1)
          Chir(4) =-Chir(2)
        endif
      return
      end function



!--------------------Weyl routines-----------------------

          SUBROUTINE ubarSpi_Weyl(p,i,ubarSpi)  ! i=+1 is ES to Chir_Weyl(.false.), i=-1 is ES to Chir_Weyl(.true.)
          implicit none
          integer, intent(in):: i
          complex(8), intent(in) :: p(4)
          complex(8) :: ubarSpi(4)
          complex(8) :: ephi
          real(8) :: p0,px,py,pz
          complex(8) :: fc, fc2

          p0=dreal(p(1))
          px=dreal(p(2))
          py=dreal(p(3))
          pz=dreal(p(4))

         fc2 = p0 + pz
         fc=cdsqrt(fc2)

         if (cdabs(fc2).gt.1D-15) then

         if (i.eq.1) then
            ubarSpi(1)=(0d0,0d0)
            ubarSpi(2)=(0d0,0d0)
            ubarSpi(3)=fc
            ubarSpi(4)=(px-(0d0,1d0)*py)/fc
         elseif (i.eq.-1) then
            ubarSpi(1)=(px+(0d0,1d0)*py)/fc
            ubarSpi(2)=-fc
            ubarSpi(3)=(0d0,0d0)
            ubarSpi(4)=(0d0,0d0)
         else
          call Error("wrong helicity setting in ubarSpi_Weyl")
         endif

         else

         if (i.eq.1) then
            ubarSpi(1) = (0d0,0d0)
            ubarSpi(2) = (0d0,0d0)
            ubarSpi(3) = (0d0,0d0)
            ubarSpi(4) = dsqrt(2d0*p0)
         elseif (i.eq.-1) then
            ubarSpi(1) = dsqrt(2d0*p0)
            ubarSpi(2) = (0d0,0d0)
            ubarSpi(3) = (0d0,0d0)
            ubarSpi(4) = (0d0,0d0)
         else
            call Error("wrong helicity setting in ubarSpi_Weyl")
         endif

         endif
        return
        END SUBROUTINE



          SUBROUTINE vSpi_Weyl(p,i,vSpi)  ! i=+1 is ES to Chir_Weyl(.false.), i=-1 is ES to Chir_Weyl(.true.)
          implicit none
          integer, intent(in):: i
          complex(8), intent(in) :: p(4)
          complex(8) :: vSpi(4)
          complex(8) :: ephi
          real(8) :: p0,px,py,pz
          real(8) :: nx,ny,nz,theta,phi
          real(8) :: ct,ct2,st,st2,cphi,sphi
          complex(8) :: fc2, fc

          p0=dreal(p(1))
          px=dreal(p(2))
          py=dreal(p(3))
          pz=dreal(p(4))

         fc2 = p0 + pz
         fc=cdsqrt(fc2)

         if (cdabs(fc2).gt.1D-15) then

         if (i.eq.1) then
            vSpi(1)=(0d0,0d0)
            vSpi(2)=(0d0,0d0)
            vSpi(3)=(px-(0d0,1d0)*py)/fc
            vSpi(4)=-fc
         elseif (i.eq.-1) then
            vSpi(1)=fc
            vSpi(2)=(px+(0d0,1d0)*py)/fc
            vSpi(3)=(0d0,0d0)
            vSpi(4)=(0d0,0d0)
         else
            call Error("wrong helicity setting in vSpi_Weyl")
         endif

         else

         if (i.eq.1) then
            vSpi(1)=(0d0,0d0)
            vSpi(2)=(0d0,0d0)
            vSpi(3)=dsqrt(2d0*p0)
            vSpi(4)=(0d0,0d0)
         elseif (i.eq.-1) then
            vSpi(1)=(0d0,0d0)
            vSpi(2)=dsqrt(2d0*p0)
            vSpi(3)=(0d0,0d0)
            vSpi(4)=(0d0,0d0)
         else
            call Error("wrong helicity setting in vSpi_Weyl")
         endif

         endif

         RETURN
         END SUBROUTINE



          SUBROUTINE uSpi_Weyl(p,i,uSpi)!  u = gamma0.ubar^+   ! NOT CHECKED
          implicit none
          integer, intent(in):: i
          complex(8), intent(in) :: p(4)
          complex(8) :: uSpi(4),ubarSpi(1:4)

               call ubarSpi_Weyl(p,i,ubarSpi)
               uSpi(1:4) = dconjg( (/ ubarSpi(3),ubarSpi(4),ubarSpi(1),ubarSpi(2)/)  )

          RETURN 
          END SUBROUTINE



          SUBROUTINE vbarSpi_Weyl(p,i,vbarSpi)!  vbar = v^+.gamma0  ! NOT CHECKED
          implicit none
          integer, intent(in):: i
          complex(8), intent(in) :: p(4)
          complex(8) :: vSpi(4),vbarSpi(1:4)

               call vSpi_Weyl(p,i,vSpi)
               vbarSpi(1:4) = dconjg( (/ vSpi(3),vSpi(4),vSpi(1),vSpi(2)/)  )

          RETURN 
          END SUBROUTINE








      function Chir_Weyl(sign,sp)   ! omega_+ * uSpi_+ = uSpi_+   and     ubarSpi_+ * omega_- = ubarSpi_+
      implicit none                 ! and of course uSpi_lambda=vSpi_-lambda  and  ubarSpi_lambda=vbarSpi_-lambda
      logical :: sign
      double complex :: sp(1:4)
      double complex :: Chir_Weyl(1:4)

        if(sign) then !omega_+
          Chir_Weyl(1) = sp(1)
          Chir_Weyl(2) = sp(2)
          Chir_Weyl(3) = 0d0
          Chir_Weyl(4) = 0d0
        else !omega_-
          Chir_Weyl(1) = 0d0
          Chir_Weyl(2) = 0d0
          Chir_Weyl(3) = sp(3)
          Chir_Weyl(4) = sp(4)
        endif
      return
      end function




        FUNCTION vbqg_Weyl(sp,e1)
        implicit none
        complex(8), intent(in) :: e1(:)
        complex(8), intent(in) :: sp(:)
        complex(8) :: vbqg_Weyl(size(sp))
        real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

            vbqg_Weyl = (0d0,-1d0)/sqrt2*spi2_Weyl(e1,sp)

        END FUNCTION




        FUNCTION vbqq_Weyl(Dv,sp1,sp2)
        implicit none
        complex(8), intent(in) :: sp1(:), sp2(:)
        integer, intent(in) ::  Dv
        integer :: i
        complex(8) :: vbqq_Weyl(Dv)
        complex(8) :: rr, va(Dv),sp1a(size(sp1))
        real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

            va=(0d0,0d0)
            vbqq_Weyl=(0d0,0d0)

            do i=1,Dv
              if (i.eq.1) then
                va(1)=(1d0,0d0)
              else
                va(i)=(-1d0,0d0)
              endif
              sp1a=spb2_Weyl(sp1,va)

              rr=(0d0,-1d0)/sqrt2*psp1_(sp1a,sp2)
              if (i.eq.1) then
                    vbqq_Weyl = vbqq_Weyl + rr*va
                else
                    vbqq_Weyl = vbqq_Weyl - rr*va
              endif
              va(i)=(0d0,0d0)
            enddo

        END FUNCTION





      function vvbqq_Weyl(sp1,sp2)
      implicit none
      complex(8), intent(in) :: sp1(:), sp2(:)
      integer :: i,j
      complex(8) :: vvbqq_Weyl(4,4)
      complex(8) :: sp1a(4)
      real(8) :: va(1:4,1:4)

         va(1,1:4)=(/+1d0,0d0,0d0,0d0/)      
         va(2,1:4)=(/0d0,-1d0,0d0,0d0/)      
         va(3,1:4)=(/0d0,0d0,-1d0,0d0/)      
         va(4,1:4)=(/0d0,0d0,0d0,-1d0/)

         do i=1,4
          do j=1,4
             sp1a=spb2_Weyl(sp1, dcmplx(va(i,1:4)))
             sp1a=spb2_Weyl(sp1a,dcmplx(va(j,1:4)))
             vvbqq_Weyl(i,j) = psp1_(sp1a,sp2)
          enddo
         enddo

      end function vvbqq_Weyl




      function vvvbqq_Weyl(sp1,sp2)
      implicit none
      complex(8), intent(in) :: sp1(:), sp2(:)
      integer :: i,j,k
      complex(8) :: vvvbqq_Weyl(4,4,4)
      complex(8) :: sp1a(4)
      real(8) :: va(1:4,1:4)

         va(1,1:4)=(/+1d0,0d0,0d0,0d0/)      
         va(2,1:4)=(/0d0,-1d0,0d0,0d0/)      
         va(3,1:4)=(/0d0,0d0,-1d0,0d0/)      
         va(4,1:4)=(/0d0,0d0,0d0,-1d0/)

         do i=1,4
          do j=1,4
           do k=1,4
              sp1a=spb2_Weyl(sp1, dcmplx(va(i,1:4)))
              sp1a=spb2_Weyl(sp1a,dcmplx(va(j,1:4)))
              sp1a=spb2_Weyl(sp1a,dcmplx(va(k,1:4)))
              vvvbqq_Weyl(i,j,k) = psp1_(sp1a,sp2)
           enddo
          enddo
         enddo

      end function vvvbqq_Weyl






        function vgq_Weyl(e1,sp)
        implicit none
        complex(8), intent(in) :: e1(:)
        complex(8), intent(in) :: sp(:)
        complex(8) :: vgq_Weyl(size(sp))
        real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

            vgq_Weyl = (0d0,-1d0)/sqrt2*spb2_Weyl(sp,e1)

        end function





      function spb2_Weyl(sp,v)
      implicit none
      complex(8), intent(in) :: sp(:),v(:)
      complex (8) :: spb2_Weyl(size(sp))
      complex(8) :: x0(4,4),xx(4,4),xy(4,4)
      complex(8) :: xz(4,4),x5(4,4)
      complex(8) :: y1,y2,y3,y4,bp,bm,cp,cm
      integer :: i,i1,i2,i3,Dv,Ds,imax

      Ds = size(sp)

      if (Ds == 4) Dv = 4
      if (Ds == 8) Dv = 6
      if (Ds == 16) Dv = 8

      imax = Ds/4

           do i=1,imax
           i1= 1+4*(i-1)
           i2=i1+3

           y1=sp(i1)
           y2=sp(i1+1)
           y3=sp(i1+2)
           y4=sp(i1+3)

           x0(1,i)=y3
           x0(2,i)=y4
           x0(3,i)=y1
           x0(4,i)=y2

           xx(1,i) = y4
           xx(2,i) = y3
           xx(3,i) = -y2
           xx(4,i) = -y1

           xy(1,i)=(0d0,1d0)*y4
           xy(2,i)=-(0d0,1d0)*y3
           xy(3,i)=-(0d0,1d0)*y2
           xy(4,i)=(0d0,1d0)*y1

           xz(1,i)=y3
           xz(2,i)=-y4
           xz(3,i)=-y1
           xz(4,i)=y2

           x5(1,i)=y1
           x5(2,i)=y2
           x5(3,i)=-y3
           x5(4,i)=-y4

           enddo

           if (Dv.eq.4) then

           do i=1,4

           spb2_Weyl(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)

           enddo

           endif

           if (Dv.eq.6) then
           bp = (v(5)+(0d0,1d0)*v(6))
           bm=(v(5)-(0d0,1d0)*v(6))

           do i=1,4

           spb2_Weyl(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)+bm*x5(i,2)

            i1 = i+4

            spb2_Weyl(i1)= v(1)*x0(i,2)-v(2)*xx(i,2)-v(3)*xy(i,2)-v(4)*xz(i,2) -bp*x5(i,1)


            enddo

           endif

           if (Dv.eq.8) then
           bp=(v(5)+(0d0,1d0)*v(6))
           bm=(v(5)-(0d0,1d0)*v(6))
           cp=(v(7)+(0d0,1d0)*v(8))
           cm=(v(7)-(0d0,1d0)*v(8))

           do i=1,4

       spb2_Weyl(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)+bm*x5(i,2)-cm*x5(i,3)

            i1 = i+4

       spb2_Weyl(i1) = v(1)*x0(i,2)-v(2)*xx(i,2) -v(3)*xy(i,2)-v(4)*xz(i,2) -bp*x5(i,1)+cm*x5(i,4)

             i2 = i1+4

       spb2_Weyl(i2)=v(1)*x0(i,3)-v(2)*xx(i,3)  -v(3)*xy(i,3)-v(4)*xz(i,3)+bm*x5(i,4)+cp*x5(i,1)

              i3=i2+4

       spb2_Weyl(i3)=v(1)*x0(i,4)-v(2)*xx(i,4)  -v(3)*xy(i,4)-v(4)*xz(i,4) -bp*x5(i,3)-cp*x5(i,2)

              enddo

              endif

               end function



         function spi2_Weyl(v,sp)
         implicit none
         complex(8), intent(in) :: sp(:),v(:)
         complex(8) :: spi2_Weyl(size(sp))
         complex(8) :: x0(4,4),xx(4,4),xy(4,4)
         complex(8) :: xz(4,4),x5(4,4)
         complex(8) ::  y1,y2,y3,y4,bp,bm,cp,cm
         integer :: i,i1,i2,i3,imax,Dv,Ds

         Ds = size(sp)

         if (Ds == 4) Dv = 4
         if (Ds == 8) Dv = 6
         if (Ds == 16) Dv = 8

         imax = Ds/4

           do i=1,imax
           i1= 1+4*(i-1)
           i2=i1+3

           y1=sp(i1)
           y2=sp(i1+1)
           y3=sp(i1+2)
           y4=sp(i1+3)

           x0(1,i)=y3
           x0(2,i)=y4
           x0(3,i)=y1
           x0(4,i)=y2


           xx(1,i) = -y4
           xx(2,i) = -y3
           xx(3,i) = y2
           xx(4,i) = y1


           xy(1,i)=(0d0,1d0)*y4
           xy(2,i)=-(0d0,1d0)*y3
           xy(3,i)=-(0d0,1d0)*y2
           xy(4,i)=(0d0,1d0)*y1

           xz(1,i)=-y3
           xz(2,i)=y4
           xz(3,i)=y1
           xz(4,i)=-y2

           x5(1,i)=y1
           x5(2,i)=y2
           x5(3,i)=-y3
           x5(4,i)=-y4

           enddo

           if(Dv.eq.4) then

           do i=1,4

           spi2_Weyl(i)=v(1)*x0(i,1)-v(2)*xx(i,1) -v(3)*xy(i,1)-v(4)*xz(i,1)
           enddo

           endif


          if (Dv.eq.6) then
           bp = (v(5)+(0d0,1d0)*v(6))
           bm=(v(5)-(0d0,1d0)*v(6))


           do i=1,4

           spi2_Weyl(i)=v(1)*x0(i,1)-v(2)*xx(i,1) -v(3)*xy(i,1)-v(4)*xz(i,1)-bp*x5(i,2)

            i1=i+4

            spi2_Weyl(i1)=v(1)*x0(i,2)-v(2)*xx(i,2) -v(3)*xy(i,2)-v(4)*xz(i,2)  +bm*x5(i,1)

            enddo

          endif

          if (Dv.eq.8) then

           bp = (v(5)+(0d0,1d0)*v(6))
           bm=(v(5)-(0d0,1d0)*v(6))
           cp=(v(7)+(0d0,1d0)*v(8))
           cm=(v(7)-(0d0,1d0)*v(8))



           do i=1,4

           spi2_Weyl(i)=v(1)*x0(i,1)-v(2)*xx(i,1) -v(3)*xy(i,1)-v(4)*xz(i,1)-bp*x5(i,2)+ cp*x5(i,3)

            i1=i+4

            spi2_Weyl(i1)=v(1)*x0(i,2)-v(2)*xx(i,2)-v(3)*xy(i,2)-v(4)*xz(i,2)+bm*x5(i,1)-cp*x5(i,4)

             i2=i1+4

            spi2_Weyl(i2)=v(1)*x0(i,3)-v(2)*xx(i,3)   -v(3)*xy(i,3)-v(4)*xz(i,3) -bp*x5(i,4)-cm*x5(i,1)

             i3=i2+4

             spi2_Weyl(i3)=v(1)*x0(i,4)-v(2)*xx(i,4) -v(3)*xy(i,4)-v(4)*xz(i,4) +bm*x5(i,3)+cm*x5(i,2)


            enddo

            endif

           end function




!--------------------Majorana spinors in Dirac convention-----------------------

          SUBROUTINE uMajoSpi(p,m,i,uSpi)!  u_Majo = C.vbar^T    where C=i*gamma2.gamma0 and vbar are in Dirac convention
          implicit none
          integer, intent(in):: i
          complex(8), intent(in) :: p(4)
          real(8), intent(in) :: m
          complex(8) :: uSpi(4),vbar(1:4)

             call vbarSpi(p,m,i,vbar)
             uSpi(1:4) = (/ -vbar(4), vbar(3), -vbar(2), vbar(1) /)

          END SUBROUTINE



          SUBROUTINE vMajoSpi(p,m,i,vSpi)!  v_Majo = C.ubar^T    where C=i*gamma2.gamma0 and ubar are in Dirac convention
          implicit none
          integer, intent(in):: i
          complex(8), intent(in) :: p(4)
          complex(8) :: vSpi(4), ubar(1:4)
          real(8), intent(in) :: m

             call ubarSpi(p,m,i,ubar)
             vSpi(1:4) = (/ -ubar(4), ubar(3), -ubar(2), ubar(1) /)

          END SUBROUTINE



          SUBROUTINE ubarMajoSpi(p,m,i,ubarSpi)!  ubar_Maja = u_Maja^+.gamma0
          implicit none
          integer, intent(in):: i
          complex(8), intent(in) :: p(4)
          real(8), intent(in) :: m
          complex(8) :: uSpi(4),ubarSpi(1:4)

             call uMajoSpi(p,m,i,uSpi)
             ubarSpi(1:4) = dconjg( (/ uSpi(1), uSpi(2), -uSpi(3), -uSpi(4) /) )

          END SUBROUTINE



          SUBROUTINE vbarMajoSpi(p,m,i,vbarSpi)!  vbar_Maja = v_Maja^+.gamma0
          implicit none
          integer, intent(in):: i
          complex(8), intent(in) :: p(4)
          real(8), intent(in) :: m
          complex(8) :: vSpi(4),vbarSpi(1:4)

             call vMajoSpi(p,m,i,vSpi)
             vbarSpi(1:4) = dconjg( (/ vSpi(1), vSpi(2), -vSpi(3), -vSpi(4) /) )

          END SUBROUTINE



!--------------------Majorana spinors in Weyl convention-----------------------



          SUBROUTINE uMajoSpi_Weyl(p,i,uSpi)!  u_Majo = C.vbar^T    where C=i*gamma2.gamma0 and vbar are in Weyl convention
          implicit none
          integer, intent(in):: i
          complex(8), intent(in) :: p(4)
          complex(8) :: uSpi(4),vbar(1:4)

             call vbarSpi_Weyl(p,i,vbar)
             uSpi(1:4) = (/ -vbar(2), vbar(1), vbar(4), -vbar(3) /)

          END SUBROUTINE



          SUBROUTINE vMajoSpi_Weyl(p,i,vSpi)!  v_Majo = C.ubar^T    where C=i*gamma2.gamma0 and ubar are in Weyl convention
          implicit none
          integer, intent(in):: i
          complex(8), intent(in) :: p(4)
          complex(8) :: vSpi(4), ubar(1:4)

             call ubarSpi_Weyl(p,i,ubar)
             vSpi(1:4) = (/ -ubar(2), ubar(1), ubar(4), -ubar(3) /)

          END SUBROUTINE



          SUBROUTINE ubarMajoSpi_Weyl(p,i,ubarSpi)!  ubar_Maja = u_Maja^+.gamma0
          implicit none
          integer, intent(in):: i
          complex(8), intent(in) :: p(4)
          complex(8) :: uSpi(4),ubarSpi(1:4)

             call uMajoSpi_Weyl(p,i,uSpi)
             ubarSpi(1:4) = dconjg( (/ uSpi(3), uSpi(4), uSpi(1), uSpi(2) /) )

          END SUBROUTINE



          SUBROUTINE vbarMajoSpi_Weyl(p,i,vbarSpi)!  vbar_Maja = v_Maja^+.gamma0
          implicit none
          integer, intent(in):: i
          complex(8), intent(in) :: p(4)
          complex(8) :: vSpi(4),vbarSpi(1:4)

             call vMajoSpi_Weyl(p,i,vSpi)
             vbarSpi(1:4) = dconjg( (/ vSpi(3), vSpi(4), vSpi(1), vSpi(2) /) )

          END SUBROUTINE




END MODULE

