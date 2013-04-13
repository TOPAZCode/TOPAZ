MODULE ModTTBJ_NLODKW
use ModMisc
use ModParameters
implicit none
!! VIRTUAL VARS ONLY
complex(8),private :: Integrals(1:13)
complex(8),private :: SME(1:9,1:4), CoeffsWm(1:9), CoeffsWmLO(1:9), CoeffsWmCTgs(1:9), CoeffsWmCTgluon(1:9)
complex(8),private :: SP4(1:4)
real(8),private :: ThetaEps
!! VIRTUAL + REAL VARS
complex(8),private :: p1w(1:4),k1ub(1:4),k2d(1:4),k3g(1:4),k4g(1:4)
complex(8),private :: Spi2d(1:4),EpsS3g(1:4),Spi1ub(1:4)
real(8),private,parameter ::  Nc =3.d0
real(8),private,parameter ::  DenNc =1.d0/3.d0
!real(8),private,parameter :: N_f=5.d0, GW=1d0,EEL=1d0
real(8),parameter :: g2_weakh = 4d0*dsqrt(2d0)*m_W**2*GF
real(8),private,parameter :: N_f=Nf_light,EEL=dsqrt(alpha*4.d0*DblPi), swq = EEL**2/g2_weakh, GW=1.d0/dsqrt(2.d0*swq)
real(8), private :: GS, alpha_S_scharf

!! REAL VARS ONLY
complex(8),private :: ggSP4(1:15), Dirac(1:34,1:4), CoeffsR(1:34), qqbSP4(1:5)
complex(8),private :: EpsS4g(1:4)
complex(8),private :: Spi3c(1:4), Spi4cb(1:4),CoeffsqqbR(1:8)
complex(8),private :: k3c(1:4), k4cb(1:4), Diracqqb(1:8,1:4)
complex(8),private :: uubSP41(1:5), uubSP42(1:5)
complex(8),private :: Spi3u(1:4), Spi4ub(1:4)
complex(8),private :: k3u(1:4),k4ub(1:4)
complex(8),private :: Diracuub1(1:8,1:4), Diracuub2(1:8,1:4), Coeffsuub1R(1:8), Coeffsuub2R(1:8)
complex(8),private :: Spi3d(1:4), Spi4db(1:4)
complex(8),private :: k3d(1:4),k4db(1:4)
complex(8),private :: ddbSP41(1:5), ddbSP42(1:5)
complex(8),private :: Diracddb1(1:8,1:4), Diracddb2(1:8,1:4), Coeffsddb1R(1:8), Coeffsddb2R(1:8)

contains


subroutine Wto3Jet_init
implicit none
alpha_S_scharf = alpha_S*RunAlphaS(NLOParam,MuRen)
 GS= dsqrt(alpha_S_scharf*4.d0*Pi)
end subroutine Wto3Jet_init

subroutine phase_space_point_virtual(momset,pol_g,LO_NLO_sw,xe,muR,EpsW)
implicit none
complex(8) ::  EpsW(1:4),amp_weyl_lo(1:4), CT(1:4), Wpol(1:4)
real(8) ::  momset(1:4,1:3)
real(8) :: iop, muR, iop2
integer :: xe, pol_g, LO_NLO_sw, polg
p1w(1:4) = dcmplx(momset(1:4,1) + momset(1:4,2) + momset(1:4,3),0.d0)
k1ub(1:4) =  dcmplx(momset(1:4,1),0.d0)
k2d(1:4) =  dcmplx(momset(1:4,2),0.d0)
k3g(1:4) =  dcmplx(momset(1:4,3),0.d0)
!xe =-1
polg =pol_g
!polg=1
!polg =2
call SPandPolWm(polg)
EpsW(1:4) = czero
if(LO_NLO_sw .eq.  0) then
   call CalcSMEWm
   call LeadingOrderWm
   EpsW(1:4) = (CoeffsWmLO(1)*SME(1,1:4) + CoeffsWmLO(2)*SME(2,1:4) + CoeffsWmLO(3)*SME(3,1:4) + CoeffsWmLO(4)*SME(4,1:4) &
        + CoeffsWmLO(5)*SME(5,1:4) + CoeffsWmLO(6)*SME(6,1:4) + CoeffsWmLO(7)*SME(7,1:4)+ CoeffsWmLO(8)*SME(8,1:4) + CoeffsWmLO(9)*SME(9,1:4))*2.d0 ! * Sqrt(Color_LO)
elseif(LO_NLO_sw .eq. 1) then
   call WmIntegrals(xe,muR)
   call CalcSMEWm
   call CalcCTsWm(xe,muR,CT)
   call LeadingOrderWm
   call WmDecSMCoeff
   iop =0.d0
   call iopwm(xe,muR, iop)

    amp_weyl_lo(1:4) = (CoeffsWmLO(1)*SME(1,1:4) + CoeffsWmLO(2)*SME(2,1:4) + CoeffsWmLO(3)*SME(3,1:4) + CoeffsWmLO(4)*SME(4,1:4) &
        + CoeffsWmLO(5)*SME(5,1:4) + CoeffsWmLO(6)*SME(6,1:4) + CoeffsWmLO(7)*SME(7,1:4)+ CoeffsWmLO(8)*SME(8,1:4) + CoeffsWmLO(9)*SME(9,1:4))

   EpsW(1:4) = (2.d0*((CoeffsWm(1)*SME(1,1:4) + CoeffsWm(2)*SME(2,1:4) + CoeffsWm(3)*SME(3,1:4) + CoeffsWm(4)*SME(4,1:4) &
        + CoeffsWm(5)*SME(5,1:4) + CoeffsWm(6)*SME(6,1:4) + CoeffsWm(7)*SME(7,1:4) + CoeffsWm(8)*SME(8,1:4) + CoeffsWm(9)*SME(9,1:4))/(4.d0*Pi)**2 &
        +CT(1:4)) + iop*amp_weyl_lo(1:4)) *2.d0 ! * Sqrt(Color_LO)

!wpol(1:4) = pol_mass(p1w,m_w,1)
!write(*,*)  "! NLO.dot.Wpol(1)", ((2.d0*((CoeffsWm(1)*SME(1,1:4) + CoeffsWm(2)*SME(2,1:4) + CoeffsWm(3)*SME(3,1:4) + CoeffsWm(4)*SME(4,1:4) &
!        + CoeffsWm(5)*SME(5,1:4) + CoeffsWm(6)*SME(6,1:4) + CoeffsWm(7)*SME(7,1:4) + CoeffsWm(8)*SME(8,1:4) + CoeffsWm(9)*SME(9,1:4))/(4.d0*Pi)**2 &
!        +CT(1:4))).dot.Wpol)*dconjg(amp_weyl_lo(1:4).dot.Wpol)
!write(*,*)  "!iop*born(1:4).dot.wpol(1)",((iop*amp_weyl_lo(1:4)).dot.Wpol)*dconjg(amp_weyl_lo(1:4).dot.Wpol)
endif

return
end subroutine phase_space_point_virtual



subroutine SPandPOLWm(pol_g)
implicit none
integer :: pol_g

if(abs(pol_g) .gt. 1) then
   EpsS3g(1:4)  = k3g(1:4)
else
   call pol_mless(k3g,pol_g,EpsS3g)
endif

call ubarSpi_Weyl(k2d,-1,Spi2d)
call vSpi_Weyl(k1ub,+1,Spi1ub)
SP4(1) = EpsS3g.dot.k1ub
SP4(2) = EpsS3g.dot.k2d
SP4(3) = k1ub.dot.p1w
SP4(4) = k2d.dot.p1w

end subroutine SPandPOLWm


subroutine CalcSMEWm
implicit none
complex(8) :: help, help4(1:4)
real(8) ::   sqrt2
integer i, Dv
Dv = 4
sqrt2 = 1.414213562373095d0

SME(1,1:4) = vbqq_weyl(Dv,Spi2d,Spi1ub)/(-ne)*sqrt2


help = weyl_sme1(Spi2d,EpsS3g,Spi1ub)
do i=1,4
   SME(2,i) = help*k1ub(i)
enddo


do i=1,4
   SME(3,i) = help*k2d(i)
enddo


do i=1,4
   SME(5,i) = help*p1w(i)
enddo


help = weyl_sme1(Spi2d,p1w,Spi1ub)
do i=1,4
   SME(4,i) = help*EpsS3g(i)
enddo

do i=1,4
   SME(6,i) = help*p1w(i)
enddo


do i=1,4
   SME(7,i) = help*k1ub(i)
enddo


do i=1,4
   SME(8,i) = help*k2d(i)
enddo

help4 = weyl_spin2(Spi2d,p1w,EpsS3g)

SME(9,1:4) = vbqq_weyl(Dv,help4,Spi1ub)/(-ne)*sqrt2


end subroutine CalcSMEWm

subroutine WmIntegrals(xe,muRF)
implicit none
real(8) ::Mw, mt,mzero,muRF, p1q,p2q,p3q, p4q, m1q,m2q,m3q,m4q, p12q,p23q, muq
complex(8) :: qlI1,qlI2, qlI3, qlI4
integer :: xe,pri
muq =muRF**2
mzero = 0.d0
mt = m_top
mw = m_w
p1q = dble(-Mw**2 + 2*SP4(3) + 2*SP4(4))
       m1q = mzero
       m2q = mzero
       Integrals(1) = qlI2(p1q,m1q,m2q,muq,xe)

p1q = dble(Mw**2 - 2*SP4(4))
       m1q = mzero
       m2q = mzero
       Integrals(2) = qlI2(p1q,m1q,m2q,muq,xe)

p1q = dble(Mw**2 - 2*SP4(3))
       m1q = mzero
       m2q = mzero
       Integrals(3) = qlI2(p1q,m1q,m2q,muq,xe)

p1q = dble(Mw**2)
       m1q = mzero
       m2q = mzero
       Integrals(4) = qlI2(p1q,m1q,m2q,muq,xe)

p1q = dble(mzero)
       p2q = dble(mzero)
       p3q = dble(-Mw**2 + 2*SP4(3) + 2*SP4(4))
       m1q = mzero
       m2q = mzero
       m3q = mzero
       Integrals(5) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(mzero)
       p2q = dble(Mw**2 - 2*SP4(3))
       p3q = dble(Mw**2)
       m1q = mzero
       m2q = mzero
       m3q = mzero
       Integrals(6) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(mzero)
       p2q = dble(mzero)
       p3q = dble(Mw**2 - 2*SP4(4))
       m1q = mzero
       m2q = mzero
       m3q = mzero
       Integrals(7) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(-Mw**2 + 2*SP4(3) + 2*SP4(4))
       p2q = dble(mzero)
       p3q = dble(Mw**2)
       m1q = mzero
       m2q = mzero
       m3q = mzero
       Integrals(8) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(Mw**2 - 2*SP4(4))
       p2q = dble(mzero)
       p3q = dble(Mw**2)
       m1q = mzero
       m2q = mzero
       m3q = mzero
       Integrals(9) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(mzero)
       p2q = dble(mzero)
       p3q = dble(Mw**2 - 2*SP4(3))
       m1q = mzero
       m2q = mzero
       m3q = mzero
       Integrals(10) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(mzero)
       p2q = dble(mzero)
       p3q = dble(mzero)
       p4q = dble(Mw**2)
       p12q = dble(-Mw**2 + 2*SP4(3) + 2*SP4(4))
       p23q = dble(Mw**2 - 2*SP4(3))

       m1q = mzero
       m2q = mzero
       m3q = mzero
       m4q = mzero

       Integrals(11) = qlI4(p1q,p2q,p3q,p4q,p12q,p23q,m1q,m2q,m3q,m4q,muq,xe)

p1q = dble(mzero)
       p2q = dble(mzero)
       p3q = dble(mzero)
       p4q = dble(Mw**2)
       p12q = dble(Mw**2 - 2*SP4(4))
       p23q = dble(Mw**2 - 2*SP4(3))

       m1q = mzero
       m2q = mzero
       m3q = mzero
       m4q = mzero

       Integrals(12) = qlI4(p1q,p2q,p3q,p4q,p12q,p23q,m1q,m2q,m3q,m4q,muq,xe)

p1q = dble(mzero)
       p2q = dble(mzero)
       p3q = dble(mzero)
       p4q = dble(Mw**2)
       p12q = dble(-Mw**2 + 2*SP4(3) + 2*SP4(4))
       p23q = dble(Mw**2 - 2*SP4(4))

       m1q = mzero
       m2q = mzero
       m3q = mzero
       m4q = mzero

       Integrals(13) = qlI4(p1q,p2q,p3q,p4q,p12q,p23q,m1q,m2q,m3q,m4q,muq,xe)


ThetaEps = 0.d0
If(xe .eq. 0) then
   ThetaEps =1.d0
endif
end subroutine WmIntegrals


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! REAL CORRECTIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SPandPOLWmR(pol_g3,pol_g4)
implicit none
integer :: pol_g3, pol_g4, Dv

if(abs(pol_g3) .gt. 1) then
   EpsS3g(1:4)  = k3g(1:4)
else
   call pol_mless(k3g,pol_g3,EpsS3g)
endif

if(abs(pol_g4) .gt. 1) then
   EpsS4g(1:4)  = k4g(1:4)
else
   call pol_mless(k4g,pol_g4,EpsS4g)
endif

call ubarSpi_Weyl(k2d,-1,Spi2d)
call vSpi_Weyl(k1ub,+1,Spi1ub)
ggSP4(1) = EpsS3g.dot.EpsS4g
ggSP4(2) = EpsS3g.dot.k1ub
ggSP4(3) = EpsS3g.dot.k2d
ggSP4(4) = EpsS3g.dot.p1w
ggSP4(5) = EpsS4g.dot.k1ub
ggSP4(6) = EpsS4g.dot.k2d
ggSP4(7) = EpsS4g.dot.k3g
ggSP4(8) = k1ub.dot.k3g
ggSP4(9) = k1ub.dot.k4g
ggSP4(10) = k1ub.dot.p1w
ggSP4(11) = k2d.dot.k3g
ggSP4(12) = k2d.dot.k4g
ggSP4(13) = k2d.dot.p1w
ggSP4(14) = k3g.dot.k4g
ggSP4(15) = k3g.dot.p1w

end subroutine SPandPOLWmR

subroutine CalcDiracWm
implicit none
complex(8) :: help, help4(1:4)
real(8) :: sqrt2
integer :: Dv
Dv = 4
sqrt2 = 1.4142135623730950d0


          help = weyl_sme1(Spi2d,EpsS3g,Spi1ub)
          Dirac(1,1:4) = help*EpsS4g(1:4)
          Dirac(2,1:4) = help*k1ub(1:4)
          Dirac(3,1:4) = help*k2d(1:4)
          Dirac(4,1:4) = help*k3g(1:4)
          Dirac(5,1:4) = help*p1w(1:4)
          help = weyl_sme1(Spi2d,EpsS4g,Spi1ub)
          Dirac(6,1:4) = help*EpsS3g(1:4)
          Dirac(7,1:4) = help*k1ub(1:4)
          Dirac(8,1:4) = help*k2d(1:4)
          Dirac(9,1:4) = help*k3g(1:4)
          Dirac(10,1:4) = help*p1w(1:4)
          help = weyl_sme1(Spi2d,k3g,Spi1ub)
          Dirac(11,1:4) = help*EpsS3g(1:4)
          Dirac(12,1:4) = help*EpsS4g(1:4)
          Dirac(13,1:4) = help*k1ub(1:4)
          Dirac(14,1:4) = help*k2d(1:4)
          Dirac(15,1:4) = help*k3g(1:4)
          Dirac(16,1:4) = help*p1w(1:4)
          help = weyl_sme1(Spi2d,p1w,Spi1ub)
          Dirac(17,1:4) = help*EpsS3g(1:4)
          Dirac(18,1:4) = help*EpsS4g(1:4)
          Dirac(19,1:4) = help*k3g(1:4)
          help = weyl_sme3(Spi2d,k3g,EpsS4g,EpsS3g,Spi1ub)
          Dirac(20,1:4) = help*k1ub(1:4)
          Dirac(21,1:4) = help*k2d(1:4)
          Dirac(22,1:4) = help*k3g(1:4)
          Dirac(23,1:4) = help*p1w(1:4)
          help = weyl_sme3(Spi2d,p1w,EpsS4g,EpsS3g,Spi1ub)
          Dirac(24,1:4) = help*k3g(1:4)
          help = weyl_sme3(Spi2d,p1w,k3g,EpsS3g,Spi1ub)
          Dirac(25,1:4) = help*EpsS4g(1:4)
          help = weyl_sme3(Spi2d,p1w,k3g,EpsS4g,Spi1ub)
          Dirac(26,1:4) = help*EpsS3g(1:4)

          help4(1:4) = weyl_spin0(Spi2d)
          Dirac(27,1:4) = vbqq_weyl(Dv,help4,Spi1ub)/(-ne)*sqrt2

          help4(1:4) = weyl_spin2(Spi2d,EpsS4g,EpsS3g)
          Dirac(28,1:4) = vbqq_weyl(Dv,help4,Spi1ub)/(-ne)*sqrt2

          help4(1:4) = weyl_spin2(Spi2d,k3g,EpsS3g)
          Dirac(29,1:4) = vbqq_weyl(Dv,help4,Spi1ub)/(-ne)*sqrt2

          help4(1:4) = weyl_spin2(Spi2d,k3g,EpsS4g)
          Dirac(30,1:4) = vbqq_weyl(Dv,help4,Spi1ub)/(-ne)*sqrt2

          help4(1:4) = weyl_spin2(Spi2d,p1w,EpsS3g)
          Dirac(31,1:4) = vbqq_weyl(Dv,help4,Spi1ub)/(-ne)*sqrt2

          help4(1:4) = weyl_spin2(Spi2d,p1w,EpsS4g)
          Dirac(32,1:4) = vbqq_weyl(Dv,help4,Spi1ub)/(-ne)*sqrt2

          help4(1:4) = weyl_spin2(Spi2d,p1w,k3g)
          Dirac(33,1:4) = vbqq_weyl(Dv,help4,Spi1ub)/(-ne)*sqrt2

          help4(1:4) = weyl_spin4(Spi2d,p1w,k3g,EpsS4g,EpsS3g)
          Dirac(34,1:4) = vbqq_weyl(Dv,help4,Spi1ub)/(-ne)*sqrt2


end subroutine CalcDiracWm

subroutine CoeffsWmRealAB
implicit none
real(8) :: Mw
complex(8) ::  s1,s2,s3,s4,s5,s6
Mw =M_W
      CoeffsR(1) = -4.d0/ggSP4(9)+2.d0*Mw**2.d0/ggSP4(9)/(ggSP4(8)+ggSP4&
      (9)+ggSP4(14))-4.d0*ggSP4(13)/ggSP4(9)/(ggSP4(8)+ggSP4&
      (9)+ggSP4(14))+4.d0*ggSP4(15)/ggSP4(9)/ggSP4(11)
      CoeffsR(2) = -4.d0*ggSP4(7)/ggSP4(14)/(ggSP4(11)+ggSP4(12&
       )+ggSP4(14))
      CoeffsR(3) = 4.d0*ggSP4(7)/ggSP4(14)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))
      CoeffsR(4) = czero
      CoeffsR(5) = -4.d0*ggSP4(7)/ggSP4(14)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))
      CoeffsR(6) = -2.d0*Mw**2.d0/ggSP4(9)/(ggSP4(8)+ggSP4(9)+ggSP4(14&
       ))+4.d0*ggSP4(13)/ggSP4(9)/(ggSP4(8)+ggSP4(9)+ggSP4(14&
       ))
      CoeffsR(7) = -4.d0*ggSP4(3)/ggSP4(11)/(ggSP4(11)+ggSP4(12&
       )+ggSP4(14))-4.d0*ggSP4(2)/ggSP4(14)/(ggSP4(11)+ggSP4(12&
       )+ggSP4(14))-4.d0*ggSP4(3)/ggSP4(14)/(ggSP4(11)+ggSP4(12&
       )+ggSP4(14))+4.d0*ggSP4(4)/ggSP4(14)/(ggSP4(11)+ggSP4&
      (12)+ggSP4(14))
      CoeffsR(8) = 4.d0*ggSP4(3)/ggSP4(9)/ggSP4(11)+4.d0*ggSP4(3)/&
      ggSP4(9)/(ggSP4(8)+ggSP4(9)+ggSP4(14))-4.d0*ggSP4(4)/&
      ggSP4(9)/(ggSP4(8)+ggSP4(9)+ggSP4(14))+4.d0*ggSP4(2)/&
      ggSP4(14)/(ggSP4(8)+ggSP4(9)+ggSP4(14))+4.d0*ggSP4(3)/&
      ggSP4(14)/(ggSP4(8)+ggSP4(9)+ggSP4(14))-4.d0*ggSP4(4)/&
      ggSP4(14)/(ggSP4(8)+ggSP4(9)+ggSP4(14))
      CoeffsR(9) = 4.d0*ggSP4(3)/ggSP4(9)/ggSP4(11)
      CoeffsR(10) = -4.d0*ggSP4(3)/ggSP4(9)/ggSP4(11)-4.d0*ggSP4(3&
       )/ggSP4(9)/(ggSP4(8)+ggSP4(9)+ggSP4(14))+4.d0*ggSP4(4)&
      /ggSP4(9)/(ggSP4(8)+ggSP4(9)+ggSP4(14))-4.d0*ggSP4(2)/&
      ggSP4(14)/(ggSP4(8)+ggSP4(9)+ggSP4(14))-4.d0*ggSP4(3)/&
      ggSP4(14)/(ggSP4(8)+ggSP4(9)+ggSP4(14))+4.d0*ggSP4(4)/&
      ggSP4(14)/(ggSP4(8)+ggSP4(9)+ggSP4(14))
      CoeffsR(11) = czero
      CoeffsR(12) = -4.d0*ggSP4(4)/ggSP4(9)/ggSP4(11)
      CoeffsR(13) = 4.d0*ggSP4(1)/ggSP4(11)/(ggSP4(11)+ggSP4(12&
       )+ggSP4(14))+4.d0*ggSP4(1)/ggSP4(14)/(ggSP4(11)+ggSP4(12&
       )+ggSP4(14))
      CoeffsR(14) = -4.d0*ggSP4(1)/ggSP4(9)/ggSP4(11)-4.d0*ggSP4(1&
       )/ggSP4(9)/(ggSP4(8)+ggSP4(9)+ggSP4(14))-4.d0*ggSP4(1)&
      /ggSP4(14)/(ggSP4(8)+ggSP4(9)+ggSP4(14))
      CoeffsR(15) = -4.d0*ggSP4(1)/ggSP4(9)/ggSP4(11)
      CoeffsR(16) = 4.d0*ggSP4(1)/ggSP4(9)/ggSP4(11)+4.d0*ggSP4(1)&
      /ggSP4(9)/(ggSP4(8)+ggSP4(9)+ggSP4(14))+4.d0*ggSP4(1)/&
      ggSP4(14)/(ggSP4(8)+ggSP4(9)+ggSP4(14))
      CoeffsR(17) = 4.d0*ggSP4(7)/ggSP4(14)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))
      CoeffsR(18) = 4.d0*ggSP4(3)/ggSP4(9)/ggSP4(11)+4.d0*ggSP4(3)&
      /ggSP4(9)/(ggSP4(8)+ggSP4(9)+ggSP4(14))-4.d0*ggSP4(4)/&
      ggSP4(9)/(ggSP4(8)+ggSP4(9)+ggSP4(14))+4.d0*ggSP4(2)/&
      ggSP4(14)/(ggSP4(8)+ggSP4(9)+ggSP4(14))+4.d0*ggSP4(3)/&
      ggSP4(14)/(ggSP4(8)+ggSP4(9)+ggSP4(14))-4.d0*ggSP4(4)/&
      ggSP4(14)/(ggSP4(8)+ggSP4(9)+ggSP4(14))
      CoeffsR(19) = -4.d0*ggSP4(1)/ggSP4(9)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))-4.d0*ggSP4(1)/ggSP4(14)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))
      CoeffsR(20) = -2.d0/ggSP4(11)/(ggSP4(11)+ggSP4(12)+ggSP4(14&
       ))
      CoeffsR(21) = 2.d0/ggSP4(9)/ggSP4(11)+2.d0/ggSP4(9)/(ggSP4(8&
       )+ggSP4(9)+ggSP4(14))
      CoeffsR(22) = 2.d0/ggSP4(9)/ggSP4(11)
      CoeffsR(23) = -2.d0/ggSP4(9)/ggSP4(11)-2.d0/ggSP4(9)/(ggSP4(8&
       )+ggSP4(9)+ggSP4(14))
      CoeffsR(24) = 2.d0/ggSP4(9)/(ggSP4(8)+ggSP4(9)+ggSP4(14))
      CoeffsR(25) = -2.d0/ggSP4(9)/ggSP4(11)-2.d0/ggSP4(9)/(ggSP4(8&
       )+ggSP4(9)+ggSP4(14))
      CoeffsR(26) = 2.d0/ggSP4(9)/(ggSP4(8)+ggSP4(9)+ggSP4(14))
      s1 = 4.d0*ggSP4(1)/ggSP4(9)-2.d0*Mw**2*ggSP4(1)/ggSP4(9)/(&
      ggSP4(8)+ggSP4(9)+ggSP4(14))+4.d0*ggSP4(1)*ggSP4(13)/&
      ggSP4(9)/(ggSP4(8)+ggSP4(9)+ggSP4(14))-Mw**2*ggSP4(1)/&
      ggSP4(14)/(ggSP4(8)+ggSP4(9)+ggSP4(14))+2.d0*ggSP4(1)*&
      ggSP4(13)/ggSP4(14)/(ggSP4(8)+ggSP4(9)+ggSP4(14))+4.d0&
      *ggSP4(3)*ggSP4(5)/ggSP4(11)/(ggSP4(11)+ggSP4(12)+&
      ggSP4(14))+4.d0*ggSP4(3)*ggSP4(6)/ggSP4(11)/(ggSP4(11)+&
      ggSP4(12)+ggSP4(14))+4.d0*ggSP4(3)*ggSP4(7)/ggSP4(11)/&
      (ggSP4(11)+ggSP4(12)+ggSP4(14))+Mw**2*ggSP4(1)/ggSP4(14&
       )/(ggSP4(11)+ggSP4(12)+ggSP4(14))+4.d0*ggSP4(2)*ggSP4&
      (5)/ggSP4(14)/(ggSP4(11)+ggSP4(12)+ggSP4(14))
      s2 = s1+4.d0*ggSP4(3)*ggSP4(5)/ggSP4(14)/(ggSP4(11)+ggSP4&
      (12)+ggSP4(14))-4.d0*ggSP4(4)*ggSP4(5)/ggSP4(14)/(ggSP4&
       (11)+ggSP4(12)+ggSP4(14))+4.d0*ggSP4(2)*ggSP4(6)/ggSP4&
       (14)/(ggSP4(11)+ggSP4(12)+ggSP4(14))+4.d0*ggSP4(3)*&
      ggSP4(6)/ggSP4(14)/(ggSP4(11)+ggSP4(12)+ggSP4(14))-4.d0*&
      ggSP4(4)*ggSP4(6)/ggSP4(14)/(ggSP4(11)+ggSP4(12)+&
      ggSP4(14))
      CoeffsR(27) = s2+4.d0*ggSP4(2)*ggSP4(7)/ggSP4(14)/(ggSP4(11&
       )+ggSP4(12)+ggSP4(14))+4.d0*ggSP4(3)*ggSP4(7)/ggSP4(14&
       )/(ggSP4(11)+ggSP4(12)+ggSP4(14))-2.d0*ggSP4(1)*ggSP4(10&
       )/ggSP4(14)/(ggSP4(11)+ggSP4(12)+ggSP4(14))-4.d0*ggSP4&
       (1)*ggSP4(15)/ggSP4(9)/ggSP4(11)-4.d0*ggSP4(1)*ggSP4(15&
       )/ggSP4(11)/(ggSP4(11)+ggSP4(12)+ggSP4(14))-4.d0*&
      ggSP4(1)*ggSP4(15)/ggSP4(14)/(ggSP4(11)+ggSP4(12)+ggSP4&
       (14))
      CoeffsR(28) = -2.d0/ggSP4(9)+Mw**2.d0/ggSP4(9)/(ggSP4(8)+ggSP4(9&
       )+ggSP4(14))-2.d0*ggSP4(13)/ggSP4(9)/(ggSP4(8)+ggSP4(9&
       )+ggSP4(14))+2.d0*ggSP4(15)/ggSP4(9)/ggSP4(11)+2.d0*ggSP4&
       (15)/ggSP4(11)/(ggSP4(11)+ggSP4(12)+ggSP4(14))
      CoeffsR(29) = -2.d0*ggSP4(5)/ggSP4(11)/(ggSP4(11)+ggSP4(12&
       )+ggSP4(14))-2.d0*ggSP4(6)/ggSP4(11)/(ggSP4(11)+ggSP4(12&
       )+ggSP4(14))-2.d0*ggSP4(7)/ggSP4(11)/(ggSP4(11)+ggSP4(12&
       )+ggSP4(14))
      CoeffsR(30) = 2.d0*ggSP4(4)/ggSP4(9)/ggSP4(11)+2.d0*ggSP4(4)&
      /ggSP4(11)/(ggSP4(11)+ggSP4(12)+ggSP4(14))
      CoeffsR(31) = -2.d0*ggSP4(7)/ggSP4(14)/(ggSP4(8)+ggSP4(9)&
      +ggSP4(14))-2.d0*ggSP4(7)/ggSP4(14)/(ggSP4(11)+ggSP4(12&
       )+ggSP4(14))
      CoeffsR(32) = -2.d0*ggSP4(3)/ggSP4(9)/ggSP4(11)-2.d0*ggSP4(3&
       )/ggSP4(9)/(ggSP4(8)+ggSP4(9)+ggSP4(14))+2.d0*ggSP4(4)&
      /ggSP4(9)/(ggSP4(8)+ggSP4(9)+ggSP4(14))-2.d0*ggSP4(2)/&
      ggSP4(14)/(ggSP4(8)+ggSP4(9)+ggSP4(14))-2.d0*ggSP4(3)/&
      ggSP4(14)/(ggSP4(8)+ggSP4(9)+ggSP4(14))+2.d0*ggSP4(4)/&
      ggSP4(14)/(ggSP4(8)+ggSP4(9)+ggSP4(14))-2.d0*ggSP4(3)/&
      ggSP4(11)/(ggSP4(11)+ggSP4(12)+ggSP4(14))-2.d0*ggSP4(2&
       )/ggSP4(14)/(ggSP4(11)+ggSP4(12)+ggSP4(14))-2.d0*ggSP4(3&
       )/ggSP4(14)/(ggSP4(11)+ggSP4(12)+ggSP4(14))+2.d0*ggSP4(4&
       )/ggSP4(14)/(ggSP4(11)+ggSP4(12)+ggSP4(14))
      CoeffsR(33) = 2.d0*ggSP4(1)/ggSP4(9)/ggSP4(11)+2.d0*ggSP4(1)&
      /ggSP4(9)/(ggSP4(8)+ggSP4(9)+ggSP4(14))+2.d0*ggSP4(1)/&
      ggSP4(14)/(ggSP4(8)+ggSP4(9)+ggSP4(14))+2.d0*ggSP4(1)/&
      ggSP4(11)/(ggSP4(11)+ggSP4(12)+ggSP4(14))+2.d0*ggSP4(1&
       )/ggSP4(14)/(ggSP4(11)+ggSP4(12)+ggSP4(14))
      CoeffsR(34) = -1.d0/ggSP4(9)/ggSP4(11)-1.d0/ggSP4(9)/(ggSP4(8&
       )+ggSP4(9)+ggSP4(14))-1.d0/ggSP4(11)/(ggSP4(11)+ggSP4(12&
       )+ggSP4(14))

end subroutine CoeffsWmRealAB

subroutine CoeffsWmRealBA
implicit none
real(8) :: Mw
complex(8) ::  s1,s2,s3,s4,s5,s6
Mw =M_W
      CoeffsR(1) = czero
      CoeffsR(2) = 4.d0*ggSP4(7)/ggSP4(8)/ggSP4(12)+4.d0*ggSP4(7)/&
      ggSP4(12)/(ggSP4(11)+ggSP4(12)+ggSP4(14))+4.d0*ggSP4(7&
       )/ggSP4(14)/(ggSP4(11)+ggSP4(12)+ggSP4(14))
      CoeffsR(3) = -4.d0*ggSP4(7)/ggSP4(8)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))-4.d0*ggSP4(7)/ggSP4(14)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))
      CoeffsR(4) = -4.d0*ggSP4(5)/ggSP4(8)/ggSP4(12)-4.d0*ggSP4(6)&
      /ggSP4(8)/ggSP4(12)
      CoeffsR(5) = 4.d0*ggSP4(7)/ggSP4(8)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))+4.d0*ggSP4(7)/ggSP4(14)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))
      CoeffsR(6) = -4.d0/ggSP4(12)
      CoeffsR(7) = 4.d0*ggSP4(2)/ggSP4(8)/ggSP4(12)+4.d0*ggSP4(2)/&
      ggSP4(12)/(ggSP4(11)+ggSP4(12)+ggSP4(14))-4.d0*ggSP4(4&
       )/ggSP4(12)/(ggSP4(11)+ggSP4(12)+ggSP4(14))+4.d0*ggSP4(2&
       )/ggSP4(14)/(ggSP4(11)+ggSP4(12)+ggSP4(14))+4.d0*ggSP4(3&
       )/ggSP4(14)/(ggSP4(11)+ggSP4(12)+ggSP4(14))-4.d0*ggSP4&
       (4)/ggSP4(14)/(ggSP4(11)+ggSP4(12)+ggSP4(14))
      CoeffsR(8) = -4.d0*ggSP4(2)/ggSP4(8)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))-4.d0*ggSP4(2)/ggSP4(14)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))-4.d0*ggSP4(3)/ggSP4(14)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))+4.d0*ggSP4(4)/ggSP4(14)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))
      CoeffsR(9) = 4.d0*ggSP4(2)/ggSP4(8)/ggSP4(12)
      CoeffsR(10) = 4.d0*ggSP4(2)/ggSP4(8)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))+4.d0*ggSP4(2)/ggSP4(14)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))+4.d0*ggSP4(3)/ggSP4(14)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))-4.d0*ggSP4(4)/ggSP4(14)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))
      CoeffsR(11) = 4.d0*ggSP4(5)/ggSP4(8)/ggSP4(12)+4.d0*ggSP4(6)&
      /ggSP4(8)/ggSP4(12)+4.d0*ggSP4(7)/ggSP4(8)/ggSP4(12)
      CoeffsR(12) = czero
      CoeffsR(13) = -4.d0*ggSP4(1)/ggSP4(14)/(ggSP4(11)+ggSP4(12&
       )+ggSP4(14))
      CoeffsR(14) = 4.d0*ggSP4(1)/ggSP4(14)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))
      CoeffsR(15) = czero
      CoeffsR(16) = -4.d0*ggSP4(1)/ggSP4(14)/(ggSP4(8)+ggSP4(9)&
      +ggSP4(14))
      CoeffsR(17) = -4.d0*ggSP4(7)/ggSP4(8)/ggSP4(12)-4.d0*ggSP4(7&
       )/ggSP4(8)/(ggSP4(8)+ggSP4(9)+ggSP4(14))-4.d0*ggSP4(7)&
      /ggSP4(14)/(ggSP4(8)+ggSP4(9)+ggSP4(14))
      CoeffsR(18) = -4.d0*ggSP4(2)/ggSP4(8)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))-4.d0*ggSP4(2)/ggSP4(14)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))-4.d0*ggSP4(3)/ggSP4(14)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))+4.d0*ggSP4(4)/ggSP4(14)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))
      CoeffsR(19) = 4.d0*ggSP4(1)/ggSP4(14)/(ggSP4(8)+ggSP4(9)+&
      ggSP4(14))
      CoeffsR(20) = -2.d0/ggSP4(8)/ggSP4(12)-2.d0/ggSP4(12)/(ggSP4(11&
       )+ggSP4(12)+ggSP4(14))
      CoeffsR(21) = 2.d0/ggSP4(8)/(ggSP4(8)+ggSP4(9)+ggSP4(14))
      CoeffsR(22) = -2.d0/ggSP4(8)/ggSP4(12)
      CoeffsR(23) = -2.d0/ggSP4(8)/(ggSP4(8)+ggSP4(9)+ggSP4(14)&
      )
      CoeffsR(24) = 2.d0/ggSP4(8)/ggSP4(12)+2.d0/ggSP4(8)/(ggSP4(8&
       )+ggSP4(9)+ggSP4(14))
      CoeffsR(25) = -2.d0/ggSP4(8)/(ggSP4(8)+ggSP4(9)+ggSP4(14)&
      )
      CoeffsR(26) = 2.d0/ggSP4(8)/ggSP4(12)+2.d0/ggSP4(8)/(ggSP4(8&
       )+ggSP4(9)+ggSP4(14))
      s1 = -4.d0*ggSP4(2)*ggSP4(5)/ggSP4(8)/ggSP4(12)-4.d0*ggSP4(2&
       )*ggSP4(6)/ggSP4(8)/ggSP4(12)-4.d0*ggSP4(2)*ggSP4(7&
       )/ggSP4(8)/ggSP4(12)+Mw**2*ggSP4(1)/ggSP4(14)/(ggSP4(8&
       )+ggSP4(9)+ggSP4(14))-2.d0*ggSP4(1)*ggSP4(13)/ggSP4(14&
       )/(ggSP4(8)+ggSP4(9)+ggSP4(14))-4.d0*ggSP4(2)*ggSP4(5&
       )/ggSP4(12)/(ggSP4(11)+ggSP4(12)+ggSP4(14))+4.d0*ggSP4(4&
       )*ggSP4(5)/ggSP4(12)/(ggSP4(11)+ggSP4(12)+ggSP4(14&
       ))-4.d0*ggSP4(2)*ggSP4(6)/ggSP4(12)/(ggSP4(11)+ggSP4(12&
       )+ggSP4(14))+4.d0*ggSP4(4)*ggSP4(6)/ggSP4(12)/(ggSP4(11&
       )+ggSP4(12)+ggSP4(14))-4.d0*ggSP4(2)*ggSP4(7)/ggSP4(12&
       )/(ggSP4(11)+ggSP4(12)+ggSP4(14))
      s2 = s1-Mw**2*ggSP4(1)/ggSP4(14)/(ggSP4(11)+ggSP4(12)+&
      ggSP4(14))-4.d0*ggSP4(2)*ggSP4(5)/ggSP4(14)/(ggSP4(11)&
      +ggSP4(12)+ggSP4(14))-4.d0*ggSP4(3)*ggSP4(5)/ggSP4(14)&
      /(ggSP4(11)+ggSP4(12)+ggSP4(14))+4.d0*ggSP4(4)*ggSP4(5&
       )/ggSP4(14)/(ggSP4(11)+ggSP4(12)+ggSP4(14))-4.d0*ggSP4(2&
       )*ggSP4(6)/ggSP4(14)/(ggSP4(11)+ggSP4(12)+ggSP4(14&
       ))
      CoeffsR(27) = s2-4.d0*ggSP4(3)*ggSP4(6)/ggSP4(14)/(ggSP4(11&
       )+ggSP4(12)+ggSP4(14))+4.d0*ggSP4(4)*ggSP4(6)/ggSP4(14&
       )/(ggSP4(11)+ggSP4(12)+ggSP4(14))-4.d0*ggSP4(2)*ggSP4(7&
       )/ggSP4(14)/(ggSP4(11)+ggSP4(12)+ggSP4(14))-4.d0*ggSP4&
      (3)*ggSP4(7)/ggSP4(14)/(ggSP4(11)+ggSP4(12)+ggSP4(14&
       ))+2.d0*ggSP4(1)*ggSP4(10)/ggSP4(14)/(ggSP4(11)+ggSP4&
      (12)+ggSP4(14))+4.d0*ggSP4(1)*ggSP4(15)/ggSP4(14)/(&
      ggSP4(11)+ggSP4(12)+ggSP4(14))
      CoeffsR(28) = 2.d0/ggSP4(12)-Mw**2.d0/ggSP4(12)/(ggSP4(11)+ggSP4&
       (12)+ggSP4(14))+2.d0*ggSP4(10)/ggSP4(12)/(ggSP4(11)+&
      ggSP4(12)+ggSP4(14))+2.d0*ggSP4(15)/ggSP4(12)/(ggSP4(11&
       )+ggSP4(12)+ggSP4(14))
      CoeffsR(29) = -2.d0*ggSP4(5)/ggSP4(8)/ggSP4(12)-2.d0*ggSP4(6&
       )/ggSP4(8)/ggSP4(12)-2.d0*ggSP4(7)/ggSP4(8)/ggSP4(12)-&
      2.d0*ggSP4(5)/ggSP4(12)/(ggSP4(11)+ggSP4(12)+ggSP4(14)&
      )-2.d0*ggSP4(6)/ggSP4(12)/(ggSP4(11)+ggSP4(12)+ggSP4(14&
       ))-2.d0*ggSP4(7)/ggSP4(12)/(ggSP4(11)+ggSP4(12)+ggSP4(14&
       ))
      CoeffsR(30) = 2.d0*ggSP4(4)/ggSP4(12)/(ggSP4(11)+ggSP4(12&
       )+ggSP4(14))
      CoeffsR(31) = 2.d0*ggSP4(7)/ggSP4(8)/ggSP4(12)+2.d0*ggSP4(7)&
      /ggSP4(8)/(ggSP4(8)+ggSP4(9)+ggSP4(14))+2.d0*ggSP4(7)/&
      ggSP4(14)/(ggSP4(8)+ggSP4(9)+ggSP4(14))+2.d0*ggSP4(7)/&
      ggSP4(12)/(ggSP4(11)+ggSP4(12)+ggSP4(14))+2.d0*ggSP4(7&
       )/ggSP4(14)/(ggSP4(11)+ggSP4(12)+ggSP4(14))
      CoeffsR(32) = 2.d0*ggSP4(2)/ggSP4(8)/ggSP4(12)+2.d0*ggSP4(2)&
      /ggSP4(8)/(ggSP4(8)+ggSP4(9)+ggSP4(14))+2.d0*ggSP4(2)/&
      ggSP4(14)/(ggSP4(8)+ggSP4(9)+ggSP4(14))+2.d0*ggSP4(3)/&
      ggSP4(14)/(ggSP4(8)+ggSP4(9)+ggSP4(14))-2.d0*ggSP4(4)/&
      ggSP4(14)/(ggSP4(8)+ggSP4(9)+ggSP4(14))+2.d0*ggSP4(2)/&
      ggSP4(12)/(ggSP4(11)+ggSP4(12)+ggSP4(14))-2.d0*ggSP4(4&
       )/ggSP4(12)/(ggSP4(11)+ggSP4(12)+ggSP4(14))+2.d0*ggSP4(2&
       )/ggSP4(14)/(ggSP4(11)+ggSP4(12)+ggSP4(14))+2.d0*ggSP4(3&
       )/ggSP4(14)/(ggSP4(11)+ggSP4(12)+ggSP4(14))-2.d0*ggSP4&
       (4)/ggSP4(14)/(ggSP4(11)+ggSP4(12)+ggSP4(14))
      CoeffsR(33) = -2.d0*ggSP4(1)/ggSP4(14)/(ggSP4(8)+ggSP4(9)&
      +ggSP4(14))-2.d0*ggSP4(1)/ggSP4(14)/(ggSP4(11)+ggSP4(12&
       )+ggSP4(14))
      CoeffsR(34) = -1.d0/ggSP4(8)/ggSP4(12)-1.d0/ggSP4(8)/(ggSP4(8&
       )+ggSP4(9)+ggSP4(14))-1.d0/ggSP4(12)/(ggSP4(11)+ggSP4(12&
       )+ggSP4(14))

end subroutine CoeffsWmRealBA

subroutine CalcRealWm(momset,pol3g,pol4g,color,res)
implicit none
integer :: pol3g,pol4g, color, pol_3g, pol_4g
complex(8) :: res(1:4)
real(8) :: momset(1:4,1:4)
p1w(1:4) = dcmplx(momset(1:4,1) +momset(1:4,2) +momset(1:4,3) +momset(1:4,4),0.d0)
k1ub(1:4) = dcmplx(momset(1:4,1),0.d0)
k2d(1:4) = dcmplx(momset(1:4,2),0.d0)
k3g(1:4) = dcmplx(momset(1:4,3),0.d0)
k4g(1:4) = dcmplx(momset(1:4,4),0.d0)
pol_3g = pol3g
pol_4g = pol4g
! Check gauge invariance
! pol_4g =2
res(1:4) = czero
call SPandPOLWmR(pol_3g,pol_4g)
call CalcDiracWm
if(color .eq. 1) then
  call CoeffsWmRealAB
elseif(color .eq. 2) then
  call CoeffsWmRealBA
endif

res(1:4) =  1.d0/4.d0* EEL* GS**2* GW* ne*(&
+ CoeffsR(1)*Dirac(1,1:4) + CoeffsR(2)*Dirac(2,1:4) + CoeffsR(3)*Dirac(3,1:4) + CoeffsR(4)*Dirac(4,1:4) + CoeffsR(5)*Dirac(5,1:4) &
+ CoeffsR(6)*Dirac(6,1:4) + CoeffsR(7)*Dirac(7,1:4) + CoeffsR(8)*Dirac(8,1:4) + CoeffsR(9)*Dirac(9,1:4) + CoeffsR(10)*Dirac(10,1:4) &
+ CoeffsR(11)*Dirac(11,1:4) + CoeffsR(12)*Dirac(12,1:4) + CoeffsR(13)*Dirac(13,1:4) + CoeffsR(14)*Dirac(14,1:4) + CoeffsR(15)*Dirac(15,1:4) &
+ CoeffsR(16)*Dirac(16,1:4) + CoeffsR(17)*Dirac(17,1:4) + CoeffsR(18)*Dirac(18,1:4) + CoeffsR(19)*Dirac(19,1:4) + CoeffsR(20)*Dirac(20,1:4) &
+ CoeffsR(21)*Dirac(21,1:4) + CoeffsR(22)*Dirac(22,1:4) + CoeffsR(23)*Dirac(23,1:4) + CoeffsR(24)*Dirac(24,1:4) + CoeffsR(25)*Dirac(25,1:4) &
+ CoeffsR(26)*Dirac(26,1:4) + CoeffsR(27)*Dirac(27,1:4) + CoeffsR(28)*Dirac(28,1:4) + CoeffsR(29)*Dirac(29,1:4) + CoeffsR(30)*Dirac(30,1:4) &
+ CoeffsR(31)*Dirac(31,1:4) + CoeffsR(32)*Dirac(32,1:4) + CoeffsR(33)*Dirac(33,1:4) + CoeffsR(34)*Dirac(34,1:4) )/dsqrt(2.d0) ! identical particles

!write(*,*) "GS,alphas, EEL, GW, swq", GS, alpha_S,EEL, GW, swq, dsqrt(2.d0)/8.d0*g2_weakh/m_W**2
!write(*,*) "", k1ub.dot.k1ub,k2d.dot.k2d,k3g.dot.k3g, p1w.dot.p1w, k4g.dot.k4g

!write(*,*) "res", res(1:4)
end subroutine CalcRealWm

subroutine SPandPOLWmqqbR(pol_c3,pol_cb4)
implicit none
integer :: pol_c3, pol_cb4, Dv


call ubarSpi_Weyl(k2d,-1,Spi2d)
call vSpi_Weyl(k1ub,+1,Spi1ub)

call ubarSpi_Weyl(k3c,pol_c3,Spi3c)
call vSpi_Weyl(k4cb,pol_cb4,Spi4cb)

qqbSP4(1) = k1ub.dot.k3c
qqbSP4(2) = k1ub.dot.k4cb
qqbSP4(3) = k2d.dot.k3c
qqbSP4(4) = k2d.dot.k4cb
qqbSP4(5) = k3c.dot.k4cb

end subroutine SPandPOLWmqqbR

subroutine CalcDiracWmqqb
implicit none
complex(8) :: help, help4(1:4), help_one4(1:4,1:4)
complex(8) :: qqb_vec1(1:4,1:4),qqb_vec2(1:4,1:4),qqb_vec3(1:4)
real(8) , parameter :: sqrt2 =1.4142135623730950d0
integer :: Dv, i
Dv = 4
          help = weyl_sme1(Spi3c,k1ub,Spi4cb)
          Diracqqb(1,1:4) = help*vbqq_weyl(Dv,Spi2d,Spi1ub)/(-ne)*sqrt2

          help = weyl_sme1(Spi3c,k2d,Spi4cb)
          Diracqqb(2,1:4) = help*vbqq_weyl(Dv,Spi2d,Spi1ub)/(-ne)*sqrt2

          help = weyl_sme1(Spi2d,k3c,Spi1ub)
          Diracqqb(3,1:4) = help*vbqq_weyl(Dv,Spi3c,Spi4cb)/(-ne)*sqrt2

          help = weyl_sme1(Spi2d,k4cb,Spi1ub)
          Diracqqb(4,1:4) = help*vbqq_weyl(Dv,Spi3c,Spi4cb)/(-ne)*sqrt2


           help_one4(1:4,1:4) = czero
           help4(1:4) = weyl_spin1(Spi2d,k3c)

           do i=1,4
              help_one4(i,i) = dcmplx(1.d0,0.d0)
              qqb_vec1(1:4,i) =  vbqq_weyl(Dv,help4,help_one4(1:4,i))/(-ne)*sqrt2
              qqb_vec2(1:4,i) =  vbqq_weyl(Dv,help_one4(1:4,i),Spi1ub)/(-ne)*sqrt2
           enddo
           qqb_vec3(1:4) = vbqq_weyl(Dv,Spi3c,Spi4cb)/(-ne)*sqrt2

           Diracqqb(5,1:4) = (qqb_vec1(1:4,1).dot.qqb_vec3(1:4))*qqb_vec2(1:4,1) + (qqb_vec1(1:4,2).dot.qqb_vec3(1:4))*qqb_vec2(1:4,2) &
                                          +  (qqb_vec1(1:4,3).dot.qqb_vec3(1:4))*qqb_vec2(1:4,3) + (qqb_vec1(1:4,4).dot.qqb_vec3(1:4))*qqb_vec2(1:4,4)


           help_one4(1:4,1:4) = czero
           help4(1:4) = weyl_spin1(Spi2d,k4cb)

           do i=1,4
              help_one4(i,i) = dcmplx(1.d0,0.d0)
              qqb_vec1(1:4,i) =  vbqq_weyl(Dv,help4,help_one4(1:4,i))/(-ne)*sqrt2
              qqb_vec2(1:4,i) =  vbqq_weyl(Dv,help_one4(1:4,i),Spi1ub)/(-ne)*sqrt2
           enddo
           qqb_vec3(1:4) = vbqq_weyl(Dv,Spi3c,Spi4cb)/(-ne)*sqrt2

           Diracqqb(6,1:4) = (qqb_vec1(1:4,1).dot.qqb_vec3(1:4))*qqb_vec2(1:4,1) + (qqb_vec1(1:4,2).dot.qqb_vec3(1:4))*qqb_vec2(1:4,2) &
                                          +  (qqb_vec1(1:4,3).dot.qqb_vec3(1:4))*qqb_vec2(1:4,3) + (qqb_vec1(1:4,4).dot.qqb_vec3(1:4))*qqb_vec2(1:4,4)


           qqb_vec3(1:4) = vbqq_weyl(Dv,Spi3c,Spi4cb)/(-ne)*sqrt2
           help4(1:4) = vbqq_weyl(Dv,Spi2d,Spi1ub)/(-ne)*sqrt2

           Diracqqb(7,1:4) = (qqb_vec3.dot.help4)*k3c(1:4)


           qqb_vec3(1:4) = vbqq_weyl(Dv,Spi3c,Spi4cb)/(-ne)*sqrt2
           help4(1:4) = vbqq_weyl(Dv,Spi2d,Spi1ub)/(-ne)*sqrt2

           Diracqqb(8,1:4) = (qqb_vec3.dot.help4)*k4cb(1:4)


end subroutine CalcDiracWmqqb

subroutine CoeffsWmRealqqb
implicit none
real(8) :: Mw
Mw =M_W
      CoeffsqqbR(1) = -1.d0/qqbSP4(5)/(qqbSP4(1)+qqbSP4(2)+qqbSP4(5&
       ))/2.d0
      CoeffsqqbR(2) = 1.d0/qqbSP4(5)/(qqbSP4(3)+qqbSP4(4)+qqbSP4(5&
       ))/2.d0
      CoeffsqqbR(3) = 1.d0/qqbSP4(5)/(qqbSP4(1)+qqbSP4(2)+qqbSP4(5&
       ))/2.d0
      CoeffsqqbR(4) = 1.d0/qqbSP4(5)/(qqbSP4(1)+qqbSP4(2)+qqbSP4(5&
       ))/2.d0
      CoeffsqqbR(5) = -(qqbSP4(1)+qqbSP4(2)+qqbSP4(3)+qqbSP4(4&
       )+2.d0*qqbSP4(5))/qqbSP4(5)/(qqbSP4(1)+qqbSP4(2)+qqbSP4&
      (5))/(qqbSP4(3)+qqbSP4(4)+qqbSP4(5))/4.d0
      CoeffsqqbR(6) = -(qqbSP4(1)+qqbSP4(2)+qqbSP4(3)+qqbSP4(4&
       )+2.d0*qqbSP4(5))/qqbSP4(5)/(qqbSP4(1)+qqbSP4(2)+qqbSP4&
      (5))/(qqbSP4(3)+qqbSP4(4)+qqbSP4(5))/4.d0
      CoeffsqqbR(7) = -1.d0/qqbSP4(5)/(qqbSP4(1)+qqbSP4(2)+qqbSP4(5&
       ))/2.d0
      CoeffsqqbR(8) = -1.d0/qqbSP4(5)/(qqbSP4(1)+qqbSP4(2)+qqbSP4(5&
       ))/2.d0

end subroutine CoeffsWmRealqqb

subroutine CalcRealWmqqb(momset,pol3c,pol4cb,res)
implicit none
integer :: pol3c,pol4cb, color
complex(8) :: res(1:4,2)
real(8) :: momset(1:4,1:4)

p1w(1:4) = dcmplx(momset(1:4,1) + momset(1:4,2) + momset(1:4,3) + momset(1:4,4),0.d0)
k1ub(1:4) = dcmplx(momset(1:4,1),0.d0)
k2d(1:4) = dcmplx(momset(1:4,2),0.d0)
k3c(1:4) = dcmplx(momset(1:4,3),0.d0)
k4cb(1:4) = dcmplx(momset(1:4,4),0.d0)

res(1:4,1) = czero
res(1:4,2) = czero
call SPandPOLWmqqbR(pol3c,pol4cb)
call CalcDiracWmqqb
call CoeffsWmRealqqb

res(1:4,1) =  1.d0/2.d0* EEL* GS**2* GW* ne*(&
+ CoeffsqqbR(1)*Diracqqb(1,1:4) + CoeffsqqbR(2)*Diracqqb(2,1:4) + CoeffsqqbR(3)*Diracqqb(3,1:4) + CoeffsqqbR(4)*Diracqqb(4,1:4) + CoeffsqqbR(5)*Diracqqb(5,1:4) &
+ CoeffsqqbR(6)*Diracqqb(6,1:4) + CoeffsqqbR(7)*Diracqqb(7,1:4) + CoeffsqqbR(8)*Diracqqb(8,1:4) )

res(1:4,2) = -DenNc*res(1:4,1)
end subroutine CalcRealWmqqb

subroutine SPandPOLWmuub1R(pol_ub1,pol_u3,pol_ub4)
implicit none
integer :: pol_u3, pol_ub4, Dv, pol_ub1


call ubarSpi_Weyl(k2d,-1,Spi2d)
call vSpi_Weyl(k1ub,pol_ub1,Spi1ub)

call ubarSpi_Weyl(k3u,pol_u3,Spi3u)
call vSpi_Weyl(k4ub,pol_ub4,Spi4ub)

uubSP41(1) = k1ub.dot.k3u
uubSP41(2) = k1ub.dot.k4ub
uubSP41(3) = k2d.dot.k3u
uubSP41(4) = k2d.dot.k4ub
uubSP41(5) = k3u.dot.k4ub

end subroutine SPandPOLWmuub1R

subroutine CalcDiracWmuub1
implicit none
complex(8) :: help, help4(1:4), help_one4(1:4,1:4)
complex(8) :: uub1_vec1(1:4,1:4),uub1_vec2(1:4,1:4),uub1_vec3(1:4)
real(8) , parameter :: sqrt2 =1.4142135623730950d0
integer :: Dv, i
Dv = 4
          help = weyl_sme1(Spi3u,k1ub,Spi4ub)
          Diracuub1(1,1:4) = help*vbqq_weyl(Dv,Spi2d,Spi1ub)/(-ne)*sqrt2

          help = weyl_sme1(Spi3u,k2d,Spi4ub)
          Diracuub1(2,1:4) = help*vbqq_weyl(Dv,Spi2d,Spi1ub)/(-ne)*sqrt2

          help = weyl_sme1(Spi2d,k3u,Spi1ub)
          Diracuub1(3,1:4) = help*vbqq_weyl(Dv,Spi3u,Spi4ub)/(-ne)*sqrt2

          help = weyl_sme1(Spi2d,k4ub,Spi1ub)
          Diracuub1(4,1:4) = help*vbqq_weyl(Dv,Spi3u,Spi4ub)/(-ne)*sqrt2


           help_one4(1:4,1:4) = czero
           help4(1:4) = weyl_spin1(Spi2d,k3u)

           do i=1,4
              help_one4(i,i) = dcmplx(1.d0,0.d0)
              uub1_vec1(1:4,i) =  vbqq_weyl(Dv,help4,help_one4(1:4,i))/(-ne)*sqrt2
              uub1_vec2(1:4,i) =  vbqq_weyl(Dv,help_one4(1:4,i),Spi1ub)/(-ne)*sqrt2
           enddo
           uub1_vec3(1:4) = vbqq_weyl(Dv,Spi3u,Spi4ub)/(-ne)*sqrt2

           Diracuub1(5,1:4) = (uub1_vec1(1:4,1).dot.uub1_vec3(1:4))*uub1_vec2(1:4,1) + (uub1_vec1(1:4,2).dot.uub1_vec3(1:4))*uub1_vec2(1:4,2) &
                                          +  (uub1_vec1(1:4,3).dot.uub1_vec3(1:4))*uub1_vec2(1:4,3) + (uub1_vec1(1:4,4).dot.uub1_vec3(1:4))*uub1_vec2(1:4,4)


           help_one4(1:4,1:4) = czero
           help4(1:4) = weyl_spin1(Spi2d,k4ub)

           do i=1,4
              help_one4(i,i) = dcmplx(1.d0,0.d0)
              uub1_vec1(1:4,i) =  vbqq_weyl(Dv,help4,help_one4(1:4,i))/(-ne)*sqrt2
              uub1_vec2(1:4,i) =  vbqq_weyl(Dv,help_one4(1:4,i),Spi1ub)/(-ne)*sqrt2
           enddo
           uub1_vec3(1:4) = vbqq_weyl(Dv,Spi3u,Spi4ub)/(-ne)*sqrt2

           Diracuub1(6,1:4) = (uub1_vec1(1:4,1).dot.uub1_vec3(1:4))*uub1_vec2(1:4,1) + (uub1_vec1(1:4,2).dot.uub1_vec3(1:4))*uub1_vec2(1:4,2) &
                                          +  (uub1_vec1(1:4,3).dot.uub1_vec3(1:4))*uub1_vec2(1:4,3) + (uub1_vec1(1:4,4).dot.uub1_vec3(1:4))*uub1_vec2(1:4,4)


           uub1_vec3(1:4) = vbqq_weyl(Dv,Spi3u,Spi4ub)/(-ne)*sqrt2
           help4(1:4) = vbqq_weyl(Dv,Spi2d,Spi1ub)/(-ne)*sqrt2

           Diracuub1(7,1:4) = (uub1_vec3.dot.help4)*k3u(1:4)


           uub1_vec3(1:4) = vbqq_weyl(Dv,Spi3u,Spi4ub)/(-ne)*sqrt2
           help4(1:4) = vbqq_weyl(Dv,Spi2d,Spi1ub)/(-ne)*sqrt2

           Diracuub1(8,1:4) = (uub1_vec3.dot.help4)*k4ub(1:4)


end subroutine CalcDiracWmuub1

subroutine CoeffsWmRealuub1
implicit none
real(8) :: Mw
Mw =M_W
      Coeffsuub1R(1) = -1.d0/uubSP41(5)/(uubSP41(1)+uubSP41(2)+&
      uubSP41(5))/2.d0
      Coeffsuub1R(2) = 1.d0/uubSP41(5)/(uubSP41(3)+uubSP41(4)+&
      uubSP41(5))/2.d0
      Coeffsuub1R(3) = 1.d0/uubSP41(5)/(uubSP41(1)+uubSP41(2)+&
      uubSP41(5))/2.d0
      Coeffsuub1R(4) = 1.d0/uubSP41(5)/(uubSP41(1)+uubSP41(2)+&
      uubSP41(5))/2.d0
      Coeffsuub1R(5) = -(uubSP41(1)+uubSP41(2)+uubSP41(3)+uubSP41&
       (4)+2.d0*uubSP41(5))/uubSP41(5)/(uubSP41(1)+uubSP41(2&
       )+uubSP41(5))/(uubSP41(3)+uubSP41(4)+uubSP41(5))/4.d0
      Coeffsuub1R(6) = -(uubSP41(1)+uubSP41(2)+uubSP41(3)+uubSP41&
       (4)+2.d0*uubSP41(5))/uubSP41(5)/(uubSP41(1)+uubSP41(2&
       )+uubSP41(5))/(uubSP41(3)+uubSP41(4)+uubSP41(5))/4.d0
      Coeffsuub1R(7) = -1.d0/uubSP41(5)/(uubSP41(1)+uubSP41(2)+&
      uubSP41(5))/2.d0
      Coeffsuub1R(8) = -1.d0/uubSP41(5)/(uubSP41(1)+uubSP41(2)+&
      uubSP41(5))/2.d0

end subroutine CoeffsWmRealuub1

subroutine CalcRealWmuub1(momset,pol1ub,pol3u,pol4ub,res)
implicit none
integer :: pol1ub,pol3u,pol4ub, color
complex(8) :: res(1:4,2)
real(8) :: momset(1:4,1:4)

p1w(1:4) = dcmplx(momset(1:4,1) + momset(1:4,2) + momset(1:4,3) + momset(1:4,4),0.d0)
k1ub(1:4) = dcmplx( momset(1:4,1),0.d0)
k2d(1:4) =  dcmplx(momset(1:4,2),0.d0)
k3u(1:4) =  dcmplx(momset(1:4,3),0.d0)
k4ub(1:4) =  dcmplx(momset(1:4,4),0.d0)

res(1:4,1) = czero
res(1:4,2) = czero
call SPandPOLWmuub1R(pol1ub,pol3u,pol4ub)
call CalcDiracWmuub1
call CoeffsWmRealuub1

res(1:4,1) =  1.d0/2.d0* EEL* GS**2* GW* ne*(&
+ Coeffsuub1R(1)*Diracuub1(1,1:4) + Coeffsuub1R(2)*Diracuub1(2,1:4) + Coeffsuub1R(3)*Diracuub1(3,1:4) + Coeffsuub1R(4)*Diracuub1(4,1:4) + Coeffsuub1R(5)*Diracuub1(5,1:4) &
+ Coeffsuub1R(6)*Diracuub1(6,1:4) + Coeffsuub1R(7)*Diracuub1(7,1:4) + Coeffsuub1R(8)*Diracuub1(8,1:4) )/dsqrt(2.d0)! identical particles

res(1:4,2) = -DenNc*res(1:4,1)
end subroutine CalcRealWmuub1


subroutine SPandPOLWmuub2R(pol_ub1,pol_u3,pol_ub4)
implicit none
integer :: pol_u3, pol_ub4, Dv, pol_ub1


call ubarSpi_Weyl(k2d,-1,Spi2d)
call vSpi_Weyl(k4ub,pol_ub4,Spi4ub)

call ubarSpi_Weyl(k3u,pol_u3,Spi3u)
call vSpi_Weyl(k1ub,pol_ub1,Spi1ub)

uubSP42(1) = k1ub.dot.k2d
uubSP42(2) = k1ub.dot.k3u
uubSP42(3) = k1ub.dot.k4ub
uubSP42(4) = k2d.dot.k3u
uubSP42(5) = k3u.dot.k4ub

end subroutine SPandPOLWmuub2R

subroutine CalcDiracWmuub2
implicit none
complex(8) :: help, help4(1:4), help_one4(1:4,1:4)
complex(8) :: uub2_vec1(1:4,1:4),uub2_vec2(1:4,1:4),uub2_vec3(1:4)
real(8) , parameter :: sqrt2 =1.4142135623730950d0
integer :: Dv, i
Dv = 4
          help = weyl_sme1(Spi3u,k2d,Spi1ub)
          Diracuub2(1,1:4) = help*vbqq_weyl(Dv,Spi2d,Spi4ub)/(-ne)*sqrt2

          help = weyl_sme1(Spi3u,k4ub,Spi1ub)
          Diracuub2(2,1:4) = help*vbqq_weyl(Dv,Spi2d,Spi4ub)/(-ne)*sqrt2

          help = weyl_sme1(Spi2d,k1ub,Spi4ub)
          Diracuub2(3,1:4) = help*vbqq_weyl(Dv,Spi3u,Spi1ub)/(-ne)*sqrt2

          help = weyl_sme1(Spi2d,k3u,Spi4ub)
          Diracuub2(4,1:4) = help*vbqq_weyl(Dv,Spi3u,Spi1ub)/(-ne)*sqrt2


           help_one4(1:4,1:4) = czero
           help4(1:4) = weyl_spin1(Spi2d,k1ub)

           do i=1,4
              help_one4(i,i) = dcmplx(1.d0,0.d0)
              uub2_vec1(1:4,i) =  vbqq_weyl(Dv,help4,help_one4(1:4,i))/(-ne)*sqrt2
              uub2_vec2(1:4,i) =  vbqq_weyl(Dv,help_one4(1:4,i),Spi4ub)/(-ne)*sqrt2
           enddo
           uub2_vec3(1:4) = vbqq_weyl(Dv,Spi3u,Spi1ub)/(-ne)*sqrt2

           Diracuub2(5,1:4) = (uub2_vec1(1:4,1).dot.uub2_vec3(1:4))*uub2_vec2(1:4,1) + (uub2_vec1(1:4,2).dot.uub2_vec3(1:4))*uub2_vec2(1:4,2) &
                                          +  (uub2_vec1(1:4,3).dot.uub2_vec3(1:4))*uub2_vec2(1:4,3) + (uub2_vec1(1:4,4).dot.uub2_vec3(1:4))*uub2_vec2(1:4,4)


           help_one4(1:4,1:4) = czero
           help4(1:4) = weyl_spin1(Spi2d,k3u)

           do i=1,4
              help_one4(i,i) = dcmplx(1.d0,0.d0)
              uub2_vec1(1:4,i) =  vbqq_weyl(Dv,help4,help_one4(1:4,i))/(-ne)*sqrt2
              uub2_vec2(1:4,i) =  vbqq_weyl(Dv,help_one4(1:4,i),Spi4ub)/(-ne)*sqrt2
           enddo
           uub2_vec3(1:4) = vbqq_weyl(Dv,Spi3u,Spi1ub)/(-ne)*sqrt2

           Diracuub2(6,1:4) = (uub2_vec1(1:4,1).dot.uub2_vec3(1:4))*uub2_vec2(1:4,1) + (uub2_vec1(1:4,2).dot.uub2_vec3(1:4))*uub2_vec2(1:4,2) &
                                          +  (uub2_vec1(1:4,3).dot.uub2_vec3(1:4))*uub2_vec2(1:4,3) + (uub2_vec1(1:4,4).dot.uub2_vec3(1:4))*uub2_vec2(1:4,4)


           uub2_vec3(1:4) = vbqq_weyl(Dv,Spi3u,Spi1ub)/(-ne)*sqrt2
           help4(1:4) = vbqq_weyl(Dv,Spi2d,Spi4ub)/(-ne)*sqrt2

           Diracuub2(7,1:4) = (uub2_vec3.dot.help4)*k1ub(1:4)


           uub2_vec3(1:4) = vbqq_weyl(Dv,Spi3u,Spi1ub)/(-ne)*sqrt2
           help4(1:4) = vbqq_weyl(Dv,Spi2d,Spi4ub)/(-ne)*sqrt2

           Diracuub2(8,1:4) = (uub2_vec3.dot.help4)*k3u(1:4)


end subroutine CalcDiracWmuub2



subroutine CoeffsWmRealuub2
implicit none
real(8) :: Mw

Mw =M_W
      Coeffsuub2R(1) = DenNc/uubSP42(2)/(uubSP42(1)+uubSP42(2)+&
      uubSP42(4))/2.d0
      Coeffsuub2R(2) = -DenNc/uubSP42(2)/(uubSP42(2)+uubSP42(3)&
      +uubSP42(5))/2.d0
      Coeffsuub2R(3) = DenNc/uubSP42(2)/(uubSP42(2)+uubSP42(3)+&
      uubSP42(5))/2.d0
      Coeffsuub2R(4) = DenNc/uubSP42(2)/(uubSP42(2)+uubSP42(3)+&
      uubSP42(5))/2.d0
      Coeffsuub2R(5) = -DenNc*(uubSP42(1)+2.d0*uubSP42(2)+uubSP42(3&
       )+uubSP42(4)+uubSP42(5))/uubSP42(2)/(uubSP42(1)+uubSP42&
       (2)+uubSP42(4))/(uubSP42(2)+uubSP42(3)+uubSP42(5)&
      )/4.d0
      Coeffsuub2R(6) = -DenNc*(uubSP42(1)+2.d0*uubSP42(2)+uubSP42(3&
       )+uubSP42(4)+uubSP42(5))/uubSP42(2)/(uubSP42(1)+uubSP42&
       (2)+uubSP42(4))/(uubSP42(2)+uubSP42(3)+uubSP42(5)&
      )/4.d0
      Coeffsuub2R(7) = -DenNc/uubSP42(2)/(uubSP42(2)+uubSP42(3)&
      +uubSP42(5))/2.d0
      Coeffsuub2R(8) = -DenNc/uubSP42(2)/(uubSP42(2)+uubSP42(3)&
      +uubSP42(5))/2.d0

end subroutine CoeffsWmRealuub2

subroutine CalcRealWmuub2(momset,pol1ub,pol3u,pol4ub,res)
implicit none
integer :: pol3u,pol4ub, color, pol1ub
complex(8) :: res(1:4,2)
real(8) :: momset(1:4,1:4)

p1w(1:4) =  dcmplx(momset(1:4,1) + momset(1:4,2) + momset(1:4,3) + momset(1:4,4),0.d0)
k1ub(1:4) =  dcmplx(momset(1:4,1),0.d0)
k2d(1:4) =  dcmplx(momset(1:4,2),0.d0)
k3u(1:4) =  dcmplx(momset(1:4,3),0.d0)
k4ub(1:4) =  dcmplx(momset(1:4,4),0.d0)

res(1:4,1) = czero
res(1:4,2) = czero
call SPandPOLWmuub2R(pol1ub,pol3u,pol4ub)
call CalcDiracWmuub2
call CoeffsWmRealuub2

res(1:4,1) =  1.d0/2.d0* EEL* GS**2* GW* ne*(&
+ Coeffsuub2R(1)*Diracuub2(1,1:4) + Coeffsuub2R(2)*Diracuub2(2,1:4) + Coeffsuub2R(3)*Diracuub2(3,1:4) + Coeffsuub2R(4)*Diracuub2(4,1:4) + Coeffsuub2R(5)*Diracuub2(5,1:4) &
+ Coeffsuub2R(6)*Diracuub2(6,1:4) + Coeffsuub2R(7)*Diracuub2(7,1:4) + Coeffsuub2R(8)*Diracuub2(8,1:4) )/dsqrt(2.d0)! identical particles

res(1:4,2) = -1.d0/DenNc*res(1:4,1)
end subroutine CalcRealWmuub2

subroutine SPandPOLWmddb1R(pol_d2,pol_d3,pol_db4)
implicit none
integer :: pol_d3, pol_db4, Dv, pol_d2


call ubarSpi_Weyl(k2d,pol_d2,Spi2d)
call vSpi_Weyl(k1ub,+1,Spi1ub)

call ubarSpi_Weyl(k3d,pol_d3,Spi3d)
call vSpi_Weyl(k4db,pol_db4,Spi4db)

ddbSP41(1) = k1ub.dot.k3d
ddbSP41(2) = k1ub.dot.k4db
ddbSP41(3) = k2d.dot.k3d
ddbSP41(4) = k2d.dot.k4db
ddbSP41(5) = k3d.dot.k4db

end subroutine SPandPOLWmddb1R

subroutine CalcDiracWmddb1
implicit none
complex(8) :: help, help4(1:4), help_one4(1:4,1:4)
complex(8) :: ddb1_vec1(1:4,1:4),ddb1_vec2(1:4,1:4),ddb1_vec3(1:4)
real(8) , parameter :: sqrt2 =1.4142135623730950d0
integer :: Dv, i
Dv = 4
          help = weyl_sme1(Spi3d,k1ub,Spi4db)
          Diracddb1(1,1:4) = help*vbqq_weyl(Dv,Spi2d,Spi1ub)/(-ne)*sqrt2

          help = weyl_sme1(Spi3d,k2d,Spi4db)
          Diracddb1(2,1:4) = help*vbqq_weyl(Dv,Spi2d,Spi1ub)/(-ne)*sqrt2

          help = weyl_sme1(Spi2d,k3d,Spi1ub)
          Diracddb1(3,1:4) = help*vbqq_weyl(Dv,Spi3d,Spi4db)/(-ne)*sqrt2

          help = weyl_sme1(Spi2d,k4db,Spi1ub)
          Diracddb1(4,1:4) = help*vbqq_weyl(Dv,Spi3d,Spi4db)/(-ne)*sqrt2


           help_one4(1:4,1:4) = czero
           help4(1:4) = weyl_spin1(Spi2d,k3d)

           do i=1,4
              help_one4(i,i) = dcmplx(1.d0,0.d0)
              ddb1_vec1(1:4,i) =  vbqq_weyl(Dv,help4,help_one4(1:4,i))/(-ne)*sqrt2
              ddb1_vec2(1:4,i) =  vbqq_weyl(Dv,help_one4(1:4,i),Spi1ub)/(-ne)*sqrt2
           enddo
           ddb1_vec3(1:4) = vbqq_weyl(Dv,Spi3d,Spi4db)/(-ne)*sqrt2

           Diracddb1(5,1:4) = (ddb1_vec1(1:4,1).dot.ddb1_vec3(1:4))*ddb1_vec2(1:4,1) + (ddb1_vec1(1:4,2).dot.ddb1_vec3(1:4))*ddb1_vec2(1:4,2) &
                                          +  (ddb1_vec1(1:4,3).dot.ddb1_vec3(1:4))*ddb1_vec2(1:4,3) + (ddb1_vec1(1:4,4).dot.ddb1_vec3(1:4))*ddb1_vec2(1:4,4)


           help_one4(1:4,1:4) = czero
           help4(1:4) = weyl_spin1(Spi2d,k4db)

           do i=1,4
              help_one4(i,i) = dcmplx(1.d0,0.d0)
              ddb1_vec1(1:4,i) =  vbqq_weyl(Dv,help4,help_one4(1:4,i))/(-ne)*sqrt2
              ddb1_vec2(1:4,i) =  vbqq_weyl(Dv,help_one4(1:4,i),Spi1ub)/(-ne)*sqrt2
           enddo
           ddb1_vec3(1:4) = vbqq_weyl(Dv,Spi3d,Spi4db)/(-ne)*sqrt2

           Diracddb1(6,1:4) = (ddb1_vec1(1:4,1).dot.ddb1_vec3(1:4))*ddb1_vec2(1:4,1) + (ddb1_vec1(1:4,2).dot.ddb1_vec3(1:4))*ddb1_vec2(1:4,2) &
                                          +  (ddb1_vec1(1:4,3).dot.ddb1_vec3(1:4))*ddb1_vec2(1:4,3) + (ddb1_vec1(1:4,4).dot.ddb1_vec3(1:4))*ddb1_vec2(1:4,4)


           ddb1_vec3(1:4) = vbqq_weyl(Dv,Spi3d,Spi4db)/(-ne)*sqrt2
           help4(1:4) = vbqq_weyl(Dv,Spi2d,Spi1ub)/(-ne)*sqrt2

           Diracddb1(7,1:4) = (ddb1_vec3.dot.help4)*k3d(1:4)


           ddb1_vec3(1:4) = vbqq_weyl(Dv,Spi3d,Spi4db)/(-ne)*sqrt2
           help4(1:4) = vbqq_weyl(Dv,Spi2d,Spi1ub)/(-ne)*sqrt2

           Diracddb1(8,1:4) = (ddb1_vec3.dot.help4)*k4db(1:4)


end subroutine CalcDiracWmddb1

subroutine CoeffsWmRealddb1
implicit none
real(8) :: Mw
Mw =M_W
      Coeffsddb1R(1) = -1.d0/ddbSP41(5)/(ddbSP41(1)+ddbSP41(2)+&
      ddbSP41(5))/2.d0
      Coeffsddb1R(2) = 1.d0/ddbSP41(5)/(ddbSP41(3)+ddbSP41(4)+&
      ddbSP41(5))/2.d0
      Coeffsddb1R(3) = 1.d0/ddbSP41(5)/(ddbSP41(1)+ddbSP41(2)+&
      ddbSP41(5))/2.d0
      Coeffsddb1R(4) = 1.d0/ddbSP41(5)/(ddbSP41(1)+ddbSP41(2)+&
      ddbSP41(5))/2.d0
      Coeffsddb1R(5) = -(ddbSP41(1)+ddbSP41(2)+ddbSP41(3)+ddbSP41&
       (4)+2.d0*ddbSP41(5))/ddbSP41(5)/(ddbSP41(1)+ddbSP41(2&
       )+ddbSP41(5))/(ddbSP41(3)+ddbSP41(4)+ddbSP41(5))/4.d0
      Coeffsddb1R(6) = -(ddbSP41(1)+ddbSP41(2)+ddbSP41(3)+ddbSP41&
       (4)+2.d0*ddbSP41(5))/ddbSP41(5)/(ddbSP41(1)+ddbSP41(2&
       )+ddbSP41(5))/(ddbSP41(3)+ddbSP41(4)+ddbSP41(5))/4.d0
      Coeffsddb1R(7) = -1.d0/ddbSP41(5)/(ddbSP41(1)+ddbSP41(2)+&
      ddbSP41(5))/2.d0
      Coeffsddb1R(8) = -1.d0/ddbSP41(5)/(ddbSP41(1)+ddbSP41(2)+&
      ddbSP41(5))/2.d0

end subroutine CoeffsWmRealddb1

subroutine CalcRealWmddb1(momset,pol2d,pol3d,pol4db,res)
implicit none
integer :: pol2d,pol3d,pol4db, color
complex(8) :: res(1:4,2)
real(8) :: momset(1:4,1:4)

p1w(1:4) = dcmplx(momset(1:4,1) + momset(1:4,2) + momset(1:4,3) + momset(1:4,4),0.d0)
k1ub(1:4) = dcmplx(momset(1:4,1),0.d0)
k2d(1:4) = dcmplx(momset(1:4,2),0.d0)
k3d(1:4) = dcmplx(momset(1:4,3),0.d0)
k4db(1:4) = dcmplx(momset(1:4,4),0.d0)

res(1:4,1) = czero
res(1:4,2) = czero
call SPandPOLWmddb1R(pol2d,pol3d,pol4db)
call CalcDiracWmddb1
call CoeffsWmRealddb1

res(1:4,1) =  1.d0/2.d0* EEL* GS**2* GW* ne*(&
+ Coeffsddb1R(1)*Diracddb1(1,1:4) + Coeffsddb1R(2)*Diracddb1(2,1:4) + Coeffsddb1R(3)*Diracddb1(3,1:4) + Coeffsddb1R(4)*Diracddb1(4,1:4) + Coeffsddb1R(5)*Diracddb1(5,1:4) &
+ Coeffsddb1R(6)*Diracddb1(6,1:4) + Coeffsddb1R(7)*Diracddb1(7,1:4) + Coeffsddb1R(8)*Diracddb1(8,1:4) )/dsqrt(2.d0)! identical particles

res(1:4,2) = -DenNc*res(1:4,1)
end subroutine CalcRealWmddb1


subroutine SPandPOLWmddb2R(pol_d2,pol_d3,pol_db4)
implicit none
integer :: pol_d3, pol_db4, Dv, pol_d2


call ubarSpi_Weyl(k2d,pol_d2,Spi2d)
call vSpi_Weyl(k4db,pol_db4,Spi4db)

call ubarSpi_Weyl(k3d,pol_d3,Spi3d)
call vSpi_Weyl(k1ub,+1,Spi1ub)

ddbSP42(1) = k1ub.dot.k2d
ddbSP42(2) = k1ub.dot.k4db
ddbSP42(3) = k2d.dot.k3d
ddbSP42(4) = k2d.dot.k4db
ddbSP42(5) = k3d.dot.k4db

end subroutine SPandPOLWmddb2R

subroutine CalcDiracWmddb2
implicit none
complex(8) :: help, help4(1:4), help_one4(1:4,1:4)
complex(8) :: ddb2_vec1(1:4,1:4),ddb2_vec2(1:4,1:4),ddb2_vec3(1:4)
real(8) , parameter :: sqrt2 =1.4142135623730950d0
integer :: Dv, i
Dv = 4
          help = weyl_sme1(Spi2d,k1ub,Spi4db)
          Diracddb2(1,1:4) = help*vbqq_weyl(Dv,Spi3d,Spi1ub)/(-ne)*sqrt2

          help = weyl_sme1(Spi2d,k3d,Spi4db)
          Diracddb2(2,1:4) = help*vbqq_weyl(Dv,Spi3d,Spi1ub)/(-ne)*sqrt2

          help = weyl_sme1(Spi3d,k2d,Spi1ub)
          Diracddb2(3,1:4) = help*vbqq_weyl(Dv,Spi2d,Spi4db)/(-ne)*sqrt2

          help = weyl_sme1(Spi3d,k4db,Spi1ub)
          Diracddb2(4,1:4) = help*vbqq_weyl(Dv,Spi2d,Spi4db)/(-ne)*sqrt2


           help_one4(1:4,1:4) = czero
           help4(1:4) = weyl_spin1(Spi3d,k2d)

           do i=1,4
              help_one4(i,i) = dcmplx(1.d0,0.d0)
              ddb2_vec1(1:4,i) =  vbqq_weyl(Dv,help4,help_one4(1:4,i))/(-ne)*sqrt2
              ddb2_vec2(1:4,i) =  vbqq_weyl(Dv,help_one4(1:4,i),Spi1ub)/(-ne)*sqrt2
           enddo
           ddb2_vec3(1:4) = vbqq_weyl(Dv,Spi2d,Spi4db)/(-ne)*sqrt2

           Diracddb2(5,1:4) = (ddb2_vec1(1:4,1).dot.ddb2_vec3(1:4))*ddb2_vec2(1:4,1) + (ddb2_vec1(1:4,2).dot.ddb2_vec3(1:4))*ddb2_vec2(1:4,2) &
                                          +  (ddb2_vec1(1:4,3).dot.ddb2_vec3(1:4))*ddb2_vec2(1:4,3) + (ddb2_vec1(1:4,4).dot.ddb2_vec3(1:4))*ddb2_vec2(1:4,4)


           help_one4(1:4,1:4) = czero
           help4(1:4) = weyl_spin1(Spi3d,k4db)

           do i=1,4
              help_one4(i,i) = dcmplx(1.d0,0.d0)
              ddb2_vec1(1:4,i) =  vbqq_weyl(Dv,help4,help_one4(1:4,i))/(-ne)*sqrt2
              ddb2_vec2(1:4,i) =  vbqq_weyl(Dv,help_one4(1:4,i),Spi1ub)/(-ne)*sqrt2
           enddo
           ddb2_vec3(1:4) = vbqq_weyl(Dv,Spi2d,Spi4db)/(-ne)*sqrt2

           Diracddb2(6,1:4) = (ddb2_vec1(1:4,1).dot.ddb2_vec3(1:4))*ddb2_vec2(1:4,1) + (ddb2_vec1(1:4,2).dot.ddb2_vec3(1:4))*ddb2_vec2(1:4,2) &
                                          +  (ddb2_vec1(1:4,3).dot.ddb2_vec3(1:4))*ddb2_vec2(1:4,3) + (ddb2_vec1(1:4,4).dot.ddb2_vec3(1:4))*ddb2_vec2(1:4,4)


           ddb2_vec3(1:4) = vbqq_weyl(Dv,Spi3d,Spi1ub)/(-ne)*sqrt2
           help4(1:4) = vbqq_weyl(Dv,Spi2d,Spi4db)/(-ne)*sqrt2

           Diracddb2(7,1:4) = (ddb2_vec3.dot.help4)*k2d(1:4)


           ddb2_vec3(1:4) = vbqq_weyl(Dv,Spi3d,Spi1ub)/(-ne)*sqrt2
           help4(1:4) = vbqq_weyl(Dv,Spi2d,Spi4db)/(-ne)*sqrt2

           Diracddb2(8,1:4) = (ddb2_vec3.dot.help4)*k4db(1:4)


end subroutine CalcDiracWmddb2



subroutine CoeffsWmRealddb2
implicit none
real(8) :: Mw

Mw =M_W
      Coeffsddb2R(1) = -DenNc/ddbSP42(4)/(ddbSP42(1)+ddbSP42(2)&
      +ddbSP42(4))/2.d0
      Coeffsddb2R(2) = DenNc/ddbSP42(4)/(ddbSP42(3)+ddbSP42(4)+&
      ddbSP42(5))/2.d0
      Coeffsddb2R(3) = DenNc/ddbSP42(4)/(ddbSP42(1)+ddbSP42(2)+&
      ddbSP42(4))/2.d0
      Coeffsddb2R(4) = DenNc/ddbSP42(4)/(ddbSP42(1)+ddbSP42(2)+&
      ddbSP42(4))/2.d0
      Coeffsddb2R(5) = -DenNc*(ddbSP42(1)+ddbSP42(2)+ddbSP42(3)&
      +2.d0*ddbSP42(4)+ddbSP42(5))/ddbSP42(4)/(ddbSP42(1)+ddbSP42&
       (2)+ddbSP42(4))/(ddbSP42(3)+ddbSP42(4)+ddbSP42(5)&
      )/4.d0
      Coeffsddb2R(6) = -DenNc*(ddbSP42(1)+ddbSP42(2)+ddbSP42(3)&
      +2.d0*ddbSP42(4)+ddbSP42(5))/ddbSP42(4)/(ddbSP42(1)+ddbSP42&
       (2)+ddbSP42(4))/(ddbSP42(3)+ddbSP42(4)+ddbSP42(5)&
      )/4.d0
      Coeffsddb2R(7) = -DenNc/ddbSP42(4)/(ddbSP42(1)+ddbSP42(2)&
      +ddbSP42(4))/2.d0
      Coeffsddb2R(8) = -DenNc/ddbSP42(4)/(ddbSP42(1)+ddbSP42(2)&
      +ddbSP42(4))/2.d0

end subroutine CoeffsWmRealddb2

subroutine CalcRealWmddb2(momset,pol2d,pol3d,pol4db,res)
implicit none
integer :: pol3d,pol4db, color, pol2d
complex(8) :: res(1:4,2)
real(8) :: momset(1:4,1:4)

p1w(1:4) = dcmplx(momset(1:4,1) + momset(1:4,2) + momset(1:4,3) + momset(1:4,4),0.d0)
k1ub(1:4) = dcmplx(momset(1:4,1),0.d0)
k2d(1:4) = dcmplx(momset(1:4,2),0.d0)
k3d(1:4) = dcmplx(momset(1:4,3),0.d0)
k4db(1:4) = dcmplx(momset(1:4,4),0.d0)

res(1:4,1) = czero
res(1:4,2) = czero
call SPandPOLWmddb2R(pol2d,pol3d,pol4db)
call CalcDiracWmddb2
call CoeffsWmRealddb2

res(1:4,1) =  1.d0/2.d0* EEL* GS**2* GW* ne*(&
+ Coeffsddb2R(1)*Diracddb2(1,1:4) + Coeffsddb2R(2)*Diracddb2(2,1:4) + Coeffsddb2R(3)*Diracddb2(3,1:4) + Coeffsddb2R(4)*Diracddb2(4,1:4) + Coeffsddb2R(5)*Diracddb2(5,1:4) &
+ Coeffsddb2R(6)*Diracddb2(6,1:4) + Coeffsddb2R(7)*Diracddb2(7,1:4) + Coeffsddb2R(8)*Diracddb2(8,1:4) )/dsqrt(2.d0)! identical particles

res(1:4,2) = -1.d0/DenNc*res(1:4,1)
end subroutine CalcRealWmddb2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!DIPOLES
subroutine WmDipolesQQB(ndip,momset,dipole,pol1,pol2,momsetT)
implicit none
complex(8) :: p1w(1:4),k1ub(1:4),k2d(1:4),k3c(1:4),k4cb(1:4)
complex(8) :: y, dip, dipole
integer :: ndip,pol1,pol2
complex(8) :: k1ubT(1:4), k2dT(1:4), kgT(1:4)
real(8) ::momsetT(1:4,1:3), momset(1:4,1:4)
real(8), parameter :: CF = 4.d0/3.d0
real(8), parameter :: CA = 3.d0

k1ub(1:4) = dcmplx(momset(1:4,1),0.d0)
k2d(1:4) =  dcmplx(momset(1:4,2),0.d0)
k3c(1:4) = dcmplx(momset(1:4,3),0.d0)
k4cb(1:4) = dcmplx(momset(1:4,4),0.d0)
k1ubT(1:4) = dcmplx(momset(1:4,1),0.d0)
k2dT(1:4) =  dcmplx(momset(1:4,2),0.d0)



if(ndip .eq. 1) then
   ! 34_1
   call qqb_fin_fin(k3c,k4cb,k1ub,pol1,pol2,y,dip)
   call finfinKinematics(k3c,k4cb,k1ub,y,kgT,k1ubT)
elseif(ndip .eq. 2) then
   ! 34_2
   call qqb_fin_fin(k3c,k4cb,k2d,pol1,pol2,y,dip)
   call finfinKinematics(k3c,k4cb,k2d,y,kgT,k2dT)
endif

momsetT(1:4,1) = dble(k1ubT(1:4))
momsetT(1:4,2) = dble(k2dT(1:4))
momsetT(1:4,3) = dble(kgT(1:4))

dipole = dip*(-CA/2.d0/CA)
return
end subroutine WmDipolesQQB



subroutine WmDipolesUUB(ndip,momset,dipole,pol1,pol2,momsetT)
implicit none
complex(8) :: p1w(1:4),k1ub(1:4),k2d(1:4),k3c(1:4),k4cb(1:4)
complex(8) :: y, dip, dipole
integer :: ndip,pol1,pol2
complex(8) :: k1ubT(1:4), k2dT(1:4), kgT(1:4)
real(8) ::momsetT(1:4,1:3), momset(1:4,1:4)
real(8), parameter :: CF = 4.d0/3.d0
real(8), parameter :: CA = 3.d0

k1ub(1:4) = dcmplx(momset(1:4,1),0.d0)
k2d(1:4) =  dcmplx(momset(1:4,2),0.d0)
k3c(1:4) = dcmplx(momset(1:4,3),0.d0)
k4cb(1:4) = dcmplx(momset(1:4,4),0.d0)

k1ubT(1:4) = dcmplx(momset(1:4,1),0.d0)
k2dT(1:4) =  dcmplx(momset(1:4,2),0.d0)



if(ndip .eq. 1) then
   ! 34_1
   call qqb_fin_fin(k3c,k4cb,k1ub,pol1,pol2,y,dip)
   call finfinKinematics(k3c,k4cb,k1ub,y,kgT,k1ubT)
elseif(ndip .eq. 2) then
   ! 34_2
   call qqb_fin_fin(k3c,k4cb,k2d,pol1,pol2,y,dip)
   call finfinKinematics(k3c,k4cb,k2d,y,kgT,k2dT)
elseif(ndip .eq. 3) then
   ! 13_2
   call qqb_fin_fin(k3c,k1ub,k2d,pol1,pol2,y,dip)
   call finfinKinematics(k3c,k1ub,k2d,y,kgT,k2dT)
   k1ubT(1:4) = k4cb(1:4)
elseif(ndip .eq. 4) then
   ! 13_4
   call qqb_fin_fin(k3c,k1ub,k4cb,pol1,pol2,y,dip)
   call finfinKinematics(k1ub,k3c,k4cb,y,kgT,k1ubT)
endif

momsetT(1:4,1) = dble(k1ubT(1:4))
momsetT(1:4,2) = dble(k2dT(1:4))
momsetT(1:4,3) = dble(kgT(1:4))

dipole = dip*(-CA/2.d0/CA)/2.d0! identical particles
return
end subroutine WmDipolesUUB


subroutine WmDipolesDDB(ndip,momset,dipole,pol1,pol2,momsetT)
implicit none
complex(8) :: p1w(1:4),k1ub(1:4),k2d(1:4),k3c(1:4),k4cb(1:4)
complex(8) :: y, dip, dipole
integer :: ndip,pol1,pol2
complex(8) :: k1ubT(1:4), k2dT(1:4), kgT(1:4)
real(8) :: momsetT(1:4,1:3), momset(1:4,1:4)
real(8), parameter :: CF = 4.d0/3.d0
real(8), parameter :: CA = 3.d0

k1ub(1:4) = dcmplx(momset(1:4,1),0.d0)
k2d(1:4) =  dcmplx(momset(1:4,2),0.d0)
k3c(1:4) = dcmplx(momset(1:4,3),0.d0)
k4cb(1:4) = dcmplx(momset(1:4,4),0.d0)

k1ubT(1:4) = dcmplx(momset(1:4,1),0.d0)
k2dT(1:4) =  dcmplx(momset(1:4,2),0.d0)



if(ndip .eq. 1) then
   ! 34_1
   call qqb_fin_fin(k3c,k4cb,k1ub,pol1,pol2,y,dip)
   call finfinKinematics(k3c,k4cb,k1ub,y,kgT,k1ubT)
elseif(ndip .eq. 2) then
   ! 34_2
   call qqb_fin_fin(k3c,k4cb,k2d,pol1,pol2,y,dip)
   call finfinKinematics(k3c,k4cb,k2d,y,kgT,k2dT)
elseif(ndip .eq. 3) then
   ! 24_3
   call qqb_fin_fin(k2d,k4cb,k3c,pol1,pol2,y,dip)
   call finfinKinematics(k2d,k4cb,k3c,y,kgT,k2dT)

elseif(ndip .eq. 4) then
   ! 24_1
   call qqb_fin_fin(k2d,k4cb,k1ub,pol1,pol2,y,dip)
   call finfinKinematics(k2d,k4cb,k1ub,y,kgT,k1ubT)
   k2dT(1:4) = k3c(1:4)
endif

momsetT(1:4,1) = dble(k1ubT(1:4))
momsetT(1:4,2) = dble(k2dT(1:4))
momsetT(1:4,3) = dble(kgT(1:4))

dipole = dip*(-CA/2.d0/CA)/2.d0! identical particles
return
end subroutine WmDipolesDDB




subroutine WmDipolesGG(ndip,momset,dipole,pol1,pol2,momsetT)
implicit none
complex(8) :: p1w(1:4),k1ub(1:4),k2d(1:4),k3g(1:4),k4g(1:4)
complex(8) :: y, dip, dipole
integer :: ndip, pol1,pol2
complex(8) :: k1ubT(1:4), k2dT(1:4), kgT(1:4)
real(8) :: momsetT(1:4,1:3), momset(1:4,1:4)
real(8), parameter :: CF = 4.d0/3.d0
real(8), parameter :: CA = 3.d0

k1ub(1:4) = dcmplx(momset(1:4,1),0.d0)
k2d(1:4) =  dcmplx(momset(1:4,2),0.d0)
k3g(1:4) = dcmplx(momset(1:4,3),0.d0)
k4g(1:4) = dcmplx(momset(1:4,4),0.d0)

k1ubT(1:4) = dcmplx(momset(1:4,1),0.d0)
k2dT(1:4) = dcmplx(momset(1:4,2),0.d0)

dipole =(0.d0,0.d0)

if(ndip .eq. 9) then
   ! 34_1
   call gg_fin_fin(k3g,k4g,k1ub,pol1,pol2,y,dip)
   call finfinKinematics(k3g,k4g,k1ub,y,kgT,k1ubT)
elseif(ndip .eq. 10) then
   ! 34_2
   call gg_fin_fin(k3g,k4g,k2d,pol1,pol2,y,dip)
   call finfinKinematics(k3g,k4g,k2d,y,kgT,k2dT)
endif

   momsetT(1:4,1) = dble(k1ubT(1:4))
   momsetT(1:4,2) = dble(k2dT(1:4))
   momsetT(1:4,3) = dble(kgT(1:4))

dipole = dip*(-CA/2.d0)/2.d0! identical particles
!write(*,*) ((k1ubT+k2dT+kgT).dot.(k1ubT+k2dT+kgT))/m_W**2
return
end subroutine WmDipolesGG




subroutine WmDipolesQG(ndip,momset,dipole,momsetT)
implicit none
complex(8) :: p1w(1:4),k1ub(1:4),k2d(1:4),k3g(1:4),k4g(1:4)
complex(8) :: y, dip, dipole
integer :: pol_g,polg, ndip
complex(8) :: k1ubT(1:4), k2dT(1:4), kgT(1:4)
real(8) :: momsetT(1:4,1:3), momset(1:4,1:4)
real(8), parameter :: CF = 4.d0/3.d0
real(8), parameter :: CA = 3.d0

k1ub(1:4) = dcmplx(momset(1:4,1),0.d0)
k2d(1:4) =  dcmplx(momset(1:4,2),0.d0)
k3g(1:4) = dcmplx(momset(1:4,3),0.d0)
k4g(1:4) = dcmplx(momset(1:4,4),0.d0)

k1ubT(1:4) = dcmplx(momset(1:4,1),0.d0)
k2dT(1:4) = dcmplx(momset(1:4,2),0.d0)

dip =(0.d0,0.d0)
if(ndip .eq. 1) then
   ! 13_2

   call qg_fin_fin(k1ub,k3g,k2d,y,dip)
   call finfinKinematics(k1ub,k3g,k2d,y,k1ubT,k2dT)

   kgT(1:4) = k4g(1:4)
   dip = dip*(-CF+CA/2.d0)
elseif(ndip .eq. 2) then
   ! 13_4
   call qg_fin_fin(k1ub,k3g,k4g,y,dip)
   call finfinKinematics(k1ub,k3g,k4g,y,k1ubT,kgT)
   dip = dip*(-CA/2.d0)

elseif(ndip .eq. 3) then
   ! 14_2
   call qg_fin_fin(k1ub,k4g,k2d,y,dip)
   call finfinKinematics(k1ub,k4g,k2d,y,k1ubT,k2dT)
   kgT(1:4) = k3g(1:4)
   dip = dip*(-CF+CA/2.d0)

elseif(ndip .eq. 4) then
   ! 14_3
   call qg_fin_fin(k1ub,k4g,k3g,y,dip)
   call finfinKinematics(k1ub,k4g,k3g,y,k1ubT,kgT)
   dip = dip*(-CA/2.d0)

elseif(ndip .eq. 5) then
   ! 23_1
   call qg_fin_fin(k2d,k3g,k1ub,y,dip)
   call finfinKinematics(k2d,k3g,k1ub,y,k2dT,k1ubT)
   kgT(1:4) = k4g(1:4)
   dip = dip*(-CF+CA/2.d0)

elseif(ndip .eq. 6) then
   ! 23_4
   call qg_fin_fin(k2d,k3g,k4g,y,dip)
   call finfinKinematics(k2d,k3g,k4g,y,k2dT,kgT)
   dip = dip*(-CA/2.d0)

elseif(ndip .eq. 7) then
   ! 24_1
   call qg_fin_fin(k2d,k4g,k1ub,y,dip)
   call finfinKinematics(k2d,k4g,k1ub,y,k2dT,k1ubT)
   kgT(1:4) = k3g(1:4)
   dip = dip*(-CF+CA/2.d0)

elseif(ndip .eq. 8) then
   ! 24_3
   call qg_fin_fin(k2d,k4g,k3g,y,dip)
   call finfinKinematics(k2d,k4g,k3g,y,k2dT,kgT)
   dip = dip*(-CA/2.d0)

endif


momsetT(1:4,1) = dble(k1ubT(1:4))
momsetT(1:4,2) = dble(k2dT(1:4))
momsetT(1:4,3) = dble(kgT(1:4))

dipole = dip/2.d0! identical particles

return
end subroutine WmDipolesQG





subroutine qg_fin_fin(emit,mited,spec,y,dip)
implicit none
complex(8) :: spec(1:4), emit(1:4), mited(1:4)
complex(8) :: y,z,dip
y = (emit.dot.mited)/((emit.dot.mited) + (emit.dot.spec) + (mited.dot.spec))
z = (emit.dot.spec)/((mited.dot.spec) + (emit.dot.spec))
dip = czero

If(dble(y) .lt. alpha_DKWff) then
   dip = -1.d0/2.d0/(emit.dot.mited)      * 8.d0*Pi*alpha_S_scharf*( 2.d0/(1.d0-z*(1.d0-y)) - (1.d0 + z))
endif
return
end subroutine qg_fin_fin





subroutine gg_fin_fin(emit,mited,spec,pol1,pol2,y,dip)
implicit none
integer:: pol1,pol2
complex(8) :: spec(1:4), emit(1:4), mited(1:4), Epsg1(1:4), Epsg2(1:4), temit(1:4), tspec(1:4)
complex(8) :: y,zi, zj, dip
y = (emit.dot.mited)/((emit.dot.mited) + (emit.dot.spec) + (mited.dot.spec))
zi = (emit.dot.spec)/((mited.dot.spec) + (emit.dot.spec))

!zj = (mited.dot.spec)/((mited.dot.spec) + (emit.dot.spec))
zj = dcmplx(1.d0,0.d0)-zi

call finfinKinematics(emit,mited,spec,y,temit,tspec)
call pol_mless(temit,pol1,Epsg1)
call pol_mless(temit,pol2,Epsg2)
Epsg1(1:4) = dconjg(Epsg1(1:4))

dip = czero
If(dble(y) .lt. alpha_DKWff) then
   dip = -1.d0/2.d0/(emit.dot.mited)  *  16.d0*Pi*alpha_S_scharf*( -(Epsg1.dot.Epsg2)*(1.d0/(1.d0-zi*(1.d0-y))+ 1.d0/(1.d0-zj*(1.d0-y)) - 2.d0 ) &
+ 1.d0/(emit.dot.mited)*(zi*(emit.dot.Epsg1)-zj*(mited.dot.Epsg1))*(zi*(emit.dot.Epsg2)-zj*(mited.dot.Epsg2)))
endif
!write(*,*) "check spin-corr-tensor 1",(zi*(emit.dot.temit)-zj*(mited.dot.temit))/M_W
return
end subroutine gg_fin_fin


subroutine qqb_fin_fin(emit,mited,spec,pol1,pol2,y,dip)
implicit none
integer:: pol1,pol2
complex(8) :: spec(1:4), emit(1:4), mited(1:4), Epsg1(1:4), Epsg2(1:4), temit(1:4), tspec(1:4)
complex(8) :: y,zi,dip, zj
y = (emit.dot.mited)/((emit.dot.mited) + (emit.dot.spec) + (mited.dot.spec))
zi = (emit.dot.spec)/((mited.dot.spec) + (emit.dot.spec))
zj =  dcmplx(1.d0,0.d0)-zi
call finfinKinematics(emit,mited,spec,y,temit,tspec)
call pol_mless(temit,pol1,Epsg1)
call pol_mless(temit,pol2,Epsg2)
Epsg1(1:4) = dconjg(Epsg1(1:4))

dip = czero
If(dble(y) .lt. alpha_DKWff) then
   dip = -1.d0/2.d0/(emit.dot.mited)  *  8.d0*Pi*alpha_S_scharf/2.d0*( -(Epsg1.dot.Epsg2) &
- 2.d0/(emit.dot.mited)*(zi*(emit.dot.Epsg1)-zj*(mited.dot.Epsg1))*(zi*(emit.dot.Epsg2)-zj*(mited.dot.Epsg2)))
endif
return
end subroutine qqb_fin_fin


subroutine finfinKinematics(emit,mited,spec,y,temit,tspec)
implicit none
complex(8) :: emit(1:4),spec(1:4),mited(1:4),temit(1:4),tspec(1:4)
complex(8) :: y

temit(1:4) = emit(1:4) + mited(1:4) - y/(1.d0-y)*spec(1:4)
tspec(1:4) = 1.d0/(1.d0-y)*spec(1:4)
return
end subroutine finfinKinematics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Integrated Dipoles

subroutine IopWm(xe,mu, iop)
use ModDKIntDipoles
use Modparameters
use Modmisc
implicit none
!complex(8) :: k1ub(1:4),k2d(1:4),k3g(1:4)
real(8), parameter :: CF = 4.d0/3.d0
real(8), parameter :: CA = 3.d0
real(8) :: s12, s13, s23
real(8) :: int_dip, dip(1:20), muq, Iop, mu
integer :: xe

muq =mu**2
int_dip =0.d0
s12 = 2.d0*dble(k1ub.dot.k2d)
s13 = 2.d0*dble(k1ub.dot.k3g)
s23 = 2.d0*dble(k2d.dot.k3g)

! q -> qg-Splitting (8 contributions)
! EMITTER IS ALWAYS a quark -> 1/T_i^2 = 1/C_F cancells C_F from mod_DKDipoles
   ! 13_2
dip(1) = DKff_qq(xe,muq,s12)
dip(1) = dip(1)*(-CF+CA/2.d0)/2.d0 !ColorCorrelation* Symmetry-Factor

   ! 13_4
dip(2) = DKff_qq(xe,muq,s13)
dip(2) = dip(2)*(-CA/2.d0)/2.d0 !ColorCorrelation* Symmetry-Factor

   ! 14_2
dip(3) = DKff_qq(xe,muq,s12)
dip(3) = dip(3)*(-CF+CA/2.d0)/2.d0 !ColorCorrelation* Symmetry-Factor

   ! 14_3
dip(4) = DKff_qq(xe,muq,s13)
dip(4) = dip(4)*(-CA/2.d0)/2.d0 !ColorCorrelation* Symmetry-Factor

   ! 23_1
dip(5) = DKff_qq(xe,muq,s12)
dip(5) = dip(5)*(-CF+CA/2.d0)/2.d0 !ColorCorrelation* Symmetry-Factor

   ! 23_4
dip(6) = DKff_qq(xe,muq,s23)
dip(6) = dip(6)*(-CA/2.d0)/2.d0 !ColorCorrelation* Symmetry-Factor

   ! 24_1
dip(7) = DKff_qq(xe,muq,s12)
dip(7) = dip(7)*(-CF+CA/2.d0)/2.d0 !ColorCorrelation* Symmetry-Factor

   ! 24_3
dip(8) = DKff_qq(xe,muq,s23)
dip(8) = dip(8)*(-CA/2.d0)/2.d0 !ColorCorrelation* Symmetry-Factor


! g -> gg-Splitting (2 contributions)
! EMITTER IS ALWAYS a gluon -> 1/T_i^2 = 1/C_A cancells C_A from mod_DKDipoles
   ! 34_1
dip(9) = DKff_gg(xe,muq,s13)
dip(9) = dip(9)*(-CA/2.d0) /2.d0 !ColorCorrelation* Symmetry-Factor

   ! 34_2
dip(10) = DKff_gg(xe,muq,s23)
dip(10) = dip(10)*(-CA/2.d0) /2.d0 !ColorCorrelation* Symmetry-Factor


! EMITTER IS ALWAYS a gluon -> 1/T_i^2 = 1/C_A DOES NOT CANCELL, since DKff_qg ~ T_R only
! g-> qqb-Splitting (2 contributions for not-identical quarks)

dip(11) = DKff_qg(xe,muq,s13)
dip(11) = dip(11)*(-CA/2.d0/CA)/2.d0 !*ColorCorrelation/ T_i^2 * T_R

dip(12) = DKff_qg(xe,muq,s23)
dip(12) = dip(12)*(-CA/2.d0/CA)/2.d0 !*ColorCorrelation/ T_i^2 * T_R

! g-> uub-Splitting (4 contributions for identical quarks)
  ! 34_2
dip(13) = DKff_qg(xe,muq,s23)
dip(13) = dip(13)*(-CA/2.d0/CA)/2.d0/2.d0 !*ColorCorrelation/ T_i^2 * T_R * Symmetry-Factor
  ! 34_1
dip(14) = DKff_qg(xe,muq,s13)
dip(14) = dip(14)*(-CA/2.d0/CA)/2.d0/2.d0 !*ColorCorrelation/ T_i^2 * T_R * Symmetry-Factor
  ! 13_2
dip(15) = DKff_qg(xe,muq,s23)
dip(15) = dip(15)*(-CA/2.d0/CA)/2.d0/2.d0 !*ColorCorrelation/ T_i^2 * T_R * Symmetry-Factor
  ! 13_4
dip(16) = DKff_qg(xe,muq,s13)
dip(16) = dip(16)*(-CA/2.d0/CA)/2.d0/2.d0 !*ColorCorrelation/ T_i^2 * T_R * Symmetry-Factor

! g-> ddb-Splitting (4 contributions for identical quarks)
  ! 34_2
dip(17) = DKff_qg(xe,muq,s23)
dip(17) = dip(17)*(-CA/2.d0/CA)/2.d0/2.d0 !*ColorCorrelation/ T_i^2 * T_R * Symmetry-Factor
  ! 34_1
dip(18) = DKff_qg(xe,muq,s13)
dip(18) = dip(18)*(-CA/2.d0/CA)/2.d0/2.d0 !*ColorCorrelation/ T_i^2 * T_R * Symmetry-Factor
  ! 24_3
dip(19) = DKff_qg(xe,muq,s23)
dip(19) = dip(19)*(-CA/2.d0/CA)/2.d0/2.d0 !*ColorCorrelation/ T_i^2 * T_R * Symmetry-Factor
  ! 24_1
dip(20) = DKff_qg(xe,muq,s13)
dip(20) = dip(20)*(-CA/2.d0/CA)/2.d0/2.d0 !*ColorCorrelation/ T_i^2 * T_R * Symmetry-Factor


int_dip =  dip(1) + dip(2) + dip(3) + dip(4) + dip(5) + dip(6) + dip(7) + dip(8) + dip(9) + dip(10) + 3.d0*(dip(11) + dip(12)) &
 + dip(13) + dip(14) + dip(15) + dip(16) + dip(17) + dip(18) + dip(19) + dip(20)
!dip(1) + dip(2) + dip(3) + dip(4) + dip(5) + dip(6) + dip(7) + dip(8) + dip(9) + dip(10)
!+ dip(11) + dip(12)
!dip(17) + dip(18) + dip(19) + dip(20)
!int_dip = dip(1)
! Looks suspicious, these n_f terms at the last two dipoles
iop = -alpha_S_scharf/2.d0/Pi* int_dip

end subroutine IopWm




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SME-COEFFICIENTS for VIRTUALS

      subroutine WmDecSMCoeff
      implicit none
            complex(8) t1
            complex(8) t3
            complex(8) t4
            complex(8) t5
            complex(8) t6
            complex(8) t7
            complex(8) t10
            complex(8) t11
            complex(8) t12
            complex(8) t13
            complex(8) t14
            complex(8) t17
            complex(8) t20
            complex(8) t21
            complex(8) t22
            complex(8) t24
            complex(8) t25
            complex(8) t26
            complex(8) t27
            complex(8) t30
            complex(8) t31
            complex(8) t33
            complex(8) t34
            complex(8) t38
            complex(8) t40
            complex(8) t41
            complex(8) t42
            complex(8) t43
            complex(8) t47
            complex(8) t48
            complex(8) t49
            complex(8) t50
            complex(8) t51
            complex(8) t53
            complex(8) t54
            complex(8) t57
            complex(8) t61
            complex(8) t62
            complex(8) t63
            complex(8) t64
            complex(8) t69
            complex(8) t70
            complex(8) t74
            complex(8) t75
            complex(8) t76
            complex(8) t77
            complex(8) t78
            complex(8) t85
            complex(8) t86
            complex(8) t87
            complex(8) t88
            complex(8) t91
            complex(8) t92
            complex(8) t93
            complex(8) t96
            complex(8) t98
            complex(8) t99
            complex(8) t101
            complex(8) t102
            complex(8) t105
            complex(8) t106
            complex(8) t109
            complex(8) t110
            complex(8) t113
            complex(8) t114
            complex(8) t115
            complex(8) t117
            complex(8) t118
            complex(8) t119
            complex(8) t121
            complex(8) t122
            complex(8) t123
            complex(8) t124
            complex(8) t126
            complex(8) t133
            complex(8) t134
            complex(8) t135
            complex(8) t136
            complex(8) t137
            complex(8) t139
            complex(8) t144
            complex(8) t145
            complex(8) t149
            complex(8) t152
            complex(8) t155
            complex(8) t158
            complex(8) t159
            complex(8) t160
            complex(8) t161
            complex(8) t163
            complex(8) t171
            complex(8) t175
            complex(8) t176
            complex(8) t177
            complex(8) t178
            complex(8) t179
            complex(8) t180
            complex(8) t183
            complex(8) t185
            complex(8) t186
            complex(8) t194
            complex(8) t195
            complex(8) t205
            complex(8) t206
            complex(8) t207
            complex(8) t208
            complex(8) t209
            complex(8) t210
            complex(8) t213
            complex(8) t214
            complex(8) t223
            complex(8) t224
            complex(8) t230
            complex(8) t239
            complex(8) t245
            complex(8) t248
            complex(8) t254
            complex(8) t259
            complex(8) t260
            complex(8) t261
            complex(8) t267
            complex(8) t272
            complex(8) t281
            complex(8) t284
            complex(8) t287
            complex(8) t288
            complex(8) t289
            complex(8) t290
            complex(8) t291
            complex(8) t294
            complex(8) t299
            complex(8) t303
            complex(8) t315
            complex(8) t316
            complex(8) t321
            complex(8) t322
            complex(8) t324
            complex(8) t325
            complex(8) t335
            complex(8) t338
            complex(8) t341
            complex(8) t343
            complex(8) t345
            complex(8) t347
            complex(8) t358
            complex(8) t359
            complex(8) t367
            complex(8) t369
            complex(8) t371
            complex(8) t372
            complex(8) t378
            complex(8) t379
            complex(8) t381
            complex(8) t384
            complex(8) t388
            complex(8) t397
            complex(8) t398
            complex(8) t401
            complex(8) t426
            complex(8) t428
            complex(8) t431
            complex(8) t432
            complex(8) t463
            complex(8) t466
            complex(8) t475
            complex(8) t492
            complex(8) t493
            complex(8) t494
            complex(8) t495
            complex(8) t511
            complex(8) t545
            complex(8) t546
            complex(8) t547
            complex(8) t549
            complex(8) t550
            complex(8) t556
            complex(8) t557
            complex(8) t559
            complex(8) t560
            complex(8) t564
            complex(8) t568
            complex(8) t571
            complex(8) t574
            complex(8) t588
            complex(8) t589
            complex(8) t599
            complex(8) t601
            complex(8) t605
            complex(8) t606
            complex(8) t620
            complex(8) t627
            complex(8) t630
            complex(8) t635
            complex(8) t640
            complex(8) t663
            complex(8) t668
            complex(8) t671
            complex(8) t678
            complex(8) t683
            complex(8) t684
            complex(8) t702
            complex(8) t703
            complex(8) t705
            complex(8) t706
            complex(8) t725
            complex(8) t727
            complex(8) t728
            complex(8) t732
            complex(8) t736
            complex(8) t738
            complex(8) t744
            complex(8) t745
            complex(8) t747
            complex(8) t749
            complex(8) t753
            complex(8) t758
            complex(8) t765
            complex(8) t767
            complex(8) t768
            complex(8) t769
            complex(8) t771
            complex(8) t772
            complex(8) t776
            complex(8) t778
            complex(8) t788
            complex(8) t794
            complex(8) t795
            complex(8) t800
            complex(8) t803
            complex(8) t807
            complex(8) t825
            complex(8) t830
            complex(8) t831
            complex(8) t832
            complex(8) t836
            complex(8) t868
            complex(8) t882
            complex(8) t910
            complex(8) t912
            complex(8) t919
            complex(8) t923
            complex(8) t926
            complex(8) t930
            complex(8) t931
            complex(8) t934
            complex(8) t937
            complex(8) t941
            complex(8) t942
            complex(8) t950
            complex(8) t954
            complex(8) t956
            complex(8) t965
            complex(8) t966
            complex(8) t982
            complex(8) t986
            complex(8) t990
            complex(8) t993
            complex(8) t997
            complex(8) t1000
            complex(8) t1010
            complex(8) t1011
            complex(8) t1026
            complex(8) t1030
            complex(8) t1036
            complex(8) t1048
            complex(8) t1068
            complex(8) t1071
            complex(8) t1074
            complex(8) t1098
            complex(8) t1102
            complex(8) t1103
            complex(8) t1105
            complex(8) t1112
            complex(8) t1113
            complex(8) t1133
            complex(8) t1140
            complex(8) t1146
            complex(8) t1149
            complex(8) t1164
            complex(8) t1192
            complex(8) t1199
            complex(8) t1202
            complex(8) t1224
            complex(8) t1225
            complex(8) t1234
            complex(8) t1277
            complex(8) t1307
            complex(8) t1308
            complex(8) t1355
            complex(8) t1465
            complex(8) t1506
      complex(8) :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
      complex(8) :: s11,s12,s13,s14,s15,s16,s17,s18,s19
      real(8) :: Mw
      Mw = M_W

      t1 = GS**2
      t3 = EEL*t1*GS
      t4 = t3*GW
      t5 = Integrals(3)
      t6 = SP4(1)
      t7 = SP4(2)
      t10 = Mw**2
      t11 = SP4(3)
      t12 = 2.d0*t11
      t13 = t10-t12
      t14 = 1.d0/t13
      t17 = t6+t7
      t20 = t4*t5*t17*t14
      t21 = 3.d0*t10
      t22 = -t21+t12
      t24 = Integrals(4)
      t25 = 3.d0*t24
      t26 = Integrals(6)
      t27 = t26*t11
      t30 = t10*(t25+4.d0*t27)
      t31 = t5*t22+t30
      t33 = 1.d0/t11
      t34 = t14*t33
      t38 = Integrals(2)
      t40 = SP4(4)
      t41 = 2.d0*t40
      t42 = t10-t41
      t43 = 1.d0/t42
      t47 = t10*t6
      t48 = 2.d0*t47
      t49 = t10*t7
      t50 = t6*t11
      t51 = t6*t40
      t53 = t7*t40
      t54 = 2.d0*t53
      t57 = t14*t43
      t61 = Integrals(1)
      t62 = t61*t13
      t63 = t6*t13
      t64 = 2.d0*t10
      t69 = t10-t11-t40
      t70 = 1.d0/t69
      t74 = t38*t13
      t75 = t42**2
      t76 = t7*t75
      t77 = 6.d0*t11
      t78 = t21-t77-t41
      t85 = t10-t12-t41
      t86 = 1.d0/t85
      t87 = t43*t86
      t88 = 1.d0/t40
      t91 = Integrals(7)
      t92 = t91*t42
      t93 = 4.d0*t11
      t96 = t63+t7*(-t21+t93+t41)
      t98 = Integrals(5)
      t99 = t98*t85
      t101 = Integrals(13)
      t102 = t101*t42
      t105 = Integrals(8)
      t106 = t105*t69
      t109 = Integrals(9)
      t110 = t109*t40
      t113 = t10*t24
      t114 = t113*t13
      t115 = t63*t40
      t117 = t10**2
      t118 = 3.d0*t117
      t119 = t11**2
      t121 = t10*t40
      t122 = 7.d0*t121
      t123 = t40**2
      t124 = 2.d0*t123
      t126 = 12.d0*t40
      t133 = 2.d0*t119
      t134 = 3.d0*t121
      t135 = 4.d0*t40
      t136 = -t21+t135
      t137 = t11*t136
      t139 = 1.d0/(t117+t133-t134+t124+t137)
      t144 = t13**2
      t145 = 1.d0/t144
      t149 = 3.d0*t40
      t152 = t7*t42
      t155 = t42*t70
      t158 = t5*t42
      t159 = t6*t144
      t160 = 6.d0*t40
      t161 = t21-t12-t160
      t163 = t7*t11
      t171 = Integrals(10)
      t175 = 4.d0*t119
      t176 = t10-t40
      t177 = t11*t176
      t178 = 8.d0*t177
      t179 = 4.d0*t121
      t180 = 8.d0*t123
      t183 = -t7*(t21-t12-t135)*t42+t6*(t117+t175-t178+t179-t180)
      t185 = Integrals(11)
      t186 = t185*t85
      t194 = 5.d0*t117
      t195 = 12.d0*t121
      t205 = t119*t11
      t206 = 4.d0*t205
      t207 = t119*t176
      t208 = 8.d0*t207
      t209 = t10*Mw
      t210 = Mw*t40
      t213 = (t209-2.d0*t210)**2
      t214 = 2.d0*t213
      t223 = t113*t42
      t224 = t163*t42
      t230 = 6.d0*t123
      t239 = 1.d0/t75
      t245 = GW*t91
      t248 = GW*t171
      t254 = -t5*t22-t30
      t259 = 2.d0*t49
      t260 = 2.d0*t50
      t261 = 2.d0*t51
      t267 = t38*t85
      t272 = t113*t85
      t281 = t11*t42
      t284 = t6*(t117+t175-4.d0*t281-t179-t180)
      t287 = Integrals(12)
      t288 = t287*t42
      t289 = 8.d0*t119
      t290 = 2.d0*t281
      t291 = 2.d0*t121
      t294 = t284+t7*(t117-t289+t290+t291-t180)
      t299 = t5*t85
      t303 = 5.d0*t10
      t315 = 4.d0*t119*(t10+t41)
      t316 = 7.d0*t117
      t321 = 8.d0*t205
      t322 = t303-t41
      t324 = 2.d0*t119*t322
      t325 = 18.d0*t121
      t335 = 14.d0*t121
      t338 = t117*t10
      t341 = t117*t40
      t343 = t10*t123
      t345 = t123*t40
      t347 = 6.d0*t121
      t358 = t85**2
      t359 = 1.d0/t358
      t367 = t3*t248
      t369 = GW*t5
      t371 = t3*t369*t14
      t372 = t10*ThetaEps
      t378 = 4.d0*t117*t26*t11
      t379 = t10*t11
      t381 = t194-4.d0*t379
      t384 = -t194+8.d0*t379-t175
      t388 = 1.d0/t119
      t397 = 4.d0*t123
      t398 = 4.d0*t10
      t401 = t117+t175-t121+t397+t11*(-t398+t41)
      t426 = t119**2
      t428 = 8.d0*t10
      t431 = t117*t75
      t432 = 2.d0*t341
      t463 = t78*t86
      t466 = t10-t12-t40
      t475 = t62*t78*t70
      t492 = 2.d0*t117
      t493 = 2.d0*t177
      t494 = 5.d0*t121
      t495 = t492-t493-t494+t397
      t511 = 17.d0*t121
      t545 = t3*t245
      t546 = 2.d0*t545
      t547 = GW*t38
      t549 = t3*t547*t43
      t550 = t43*t88
      t556 = 4.d0*t117*t109*t40
      t557 = -t194+t179
      t559 = 8.d0*t121
      t560 = t194-t559+t397
      t564 = 1.d0/t123
      t568 = t175-t281+t75
      t571 = t287*t13
      t574 = t171*t13
      t588 = 2.d0*t338
      t589 = 7.d0*t341
      t599 = t117-t124
      t601 = 4.d0*t119*t599
      t605 = 2.d0*t343
      t606 = 4.d0*t345
      t620 = 2.d0*t549
      t627 = t42*t86
      t630 = t10-t11-t41
      t635 = t185*t13
      t640 = t61*t161*t155
      t663 = t13*t70
      t668 = t175+2.d0*t10*t176-t11*t322
      t671 = t101*t85
      t678 = t11*(t21+t41)
      t683 = 5.d0*t341
      t684 = 8.d0*t117
      t702 = t117**2
      t703 = t117*t123
      t705 = t10*t345
      t706 = 4.d0*t705
      t725 = GW*ThetaEps
      t727 = t3*t725*t43
      t728 = t21-t41
      t732 = t10*(t25+4.d0*t110)
      t736 = t4*(t38*t728-t732)*t43*t88
      t738 = t38*t43
      t744 = t574*(-t93+t303-10.d0*t40)*t43*t86
      t745 = t10+t93-t41
      t747 = t91*t745*t86
      t749 = t571*t745*t86
      t753 = 2.d0*t27*t745*t43*t86
      t758 = 2.d0*t109*(t492-4.d0*t177-t134-t124)*t87
      t765 = 3.d0*t727
      t767 = t574*t85
      t768 = t98*t358
      t769 = t635*t358
      t771 = t61*t42
      t772 = t85*t70
      t776 = 2.d0*t105*t85*t69
      t778 = -2.d0*t27*t85
      t788 = t4*(-t38*t728+t732)*t43*t88
      t794 = t117+t289+t179-t397-2.d0*t678
      t795 = t794*t43
      t800 = t118+t289-t179+t397+t11*(-10.d0*t10+t135)
      t803 = t800*t43
      t807 = 4.d0*t343
      t825 = t5-t24
      t830 = t38*t176
      t831 = 2.d0*t830
      t832 = t64-t40
      t836 = 2.d0*t10*t109*t40
      t868 = t24*t42
      t882 = t24*t13
      t910 = t38*t7
      t912 = t5*t6
      t919 = t33*t43
      t923 = t3*t725
      t926 = t43*t70
      t930 = t91*t7
      t931 = t42*t85
      t934 = t98*t7
      t937 = t101*t7
      t941 = t105*t7
      t942 = t85*t69
      t950 = t109*t7
      t954 = t21-t12-t41
      t956 = 12.d0*t119
      t965 = t69**2
      t966 = 1.d0/t965
      t982 = 1.d0/t144/t13
      t986 = t171*t6
      t990 = t98*t6
      t993 = t185*t6
      t997 = t105*t6
      t1000 = t26*t6
      t1010 = 16.d0*t121
      t1011 = 12.d0*t123
      t1026 = t728**2
      t1030 = (t209-t210)**2
      t1036 = 8.d0*t345
      t1048 = 1.d0/t75/t42
      t1068 = t13*t42
      t1071 = t85*t14
      t1074 = 6.d0*t224
      t1098 = t49*t11
      t1102 = t47*t40
      t1103 = t50*t40
      t1105 = t163*t40
      t1112 = t21-t77-t135
      t1113 = t1112*t42
      t1133 = 9.d0*t117
      t1140 = t86*t966
      t1146 = (t209-2.d0*Mw*t11)**2
      t1149 = 22.d0*t10
      t1164 = t48-t152
      t1192 = 22.d0*t121
      t1199 = t64-t149
      t1202 = t176**2
      t1224 = (t117*Mw-3.d0*t209*t40+2.d0*Mw*t123)**2
      t1225 = 6.d0*t338
      t1234 = t123**2
      t1277 = 6.d0*t115
      t1307 = t21-t93-t160
      t1308 = t13*t1307
      t1355 = -t259+t63
      t1465 = t4*ThetaEps*t69*t57
      t1506 = t57*t86
      s1 = ne
      s4 = DenNc
      s6 = t4*t5*(t6-t7)*t14+t20+t4*t17*t31*t34/2.d0
      s5 = s6+2.d0*t4*t38*t6*t43+2.d0*t4*ThetaEps*(t48+t49-t50-3.d0*t51-t54)*t57+&
      t3*GW*(-2.d0*t62*(t63+t7*(t64-3.d0*t11-t40))*t70-t74*(t76*t78-8.d0*t50*t13*&
      t40)*t87*t88+t92*t96-t99*t96+t102*t85*t96-2.d0*t106*t96-2.d0*t110*t96+&
      t114*(-2.d0*t115+t7*(t118+6.d0*t119-t122+t124+t11*(-9.d0*t10+t126)))*t88*t139&
      )*t145/2.d0+t3*GW*(2.d0*t61*(t6*(t64-t11-t149)+t152)*t155+t158*(t159*&
      t161-8.d0*t163*t40*t42)*t34*t86+t171*t183+t186*t183-2.d0*t106*t183*t14-t99*&
      (t7*(t10+t12-t135)*t42+t6*(t194+t175-t178-t195+t180))*t14-2.d0*t26*(&
      t152*(t133+2.d0*t10*t42+t137)+t6*(t206-t208+t214+t11*(t117+t179-t180))&
      )*t14-t223*(-2.d0*t224+t6*(t133+t11*(-7.d0*t10+t126)+t118-9.d0*t121+t230))*&
      t33*t139)*t239/2.d0
      s3 = s4*s5
      s4 = Nc*(-t3*t245*t6-t3*t248*(2.d0*t6+t7)-t20+t4*t17*t254*t34/2.d0-t4*&
      ThetaEps*(t259+t260-t261-4.d0*t53)*t57+t3*GW*(3.d0*t267*(t152-t261)*t88-3.d0*&
      t272*(t163-t51)*t33*t88-t171*(3.d0*t7*t13*t42+t284)+t288*t294+2.d0*t110*&
      t294*t14-t299*(4.d0*t163*(-t10+t12-t40)+t6*(t118+t289-2.d0*t11*(t303+t41&
      )))*t14*t33+2.d0*t26*(t6*(t206+t315+t214+t11*(-t316+t195-t180))+t7*(-&
      t321+t324+t214+t11*(-t316+t325-t180)))*t14-t91*(t152*(t194+t289-14.d0&
      *t281-t335+t180)+t6*(3.d0*t338+16.d0*t205-22.d0*t341+40.d0*t343-16.d0*t345-t315-8.d0&
      *t11*(t117-t347+t180)))*t14)*t359/2.d0)
      s2 = s3+s4
      CoeffsWm(1) = s1*s2
      s1 = ne
      s3 = Nc*(2.d0*t367+t371+t4*t372*t34/2.d0+t4*(t378+t24*t381+t5*t384)*t14*&
      t388/4.d0+t3*GW*(6.d0*t267-t113*(t64-t11-t135)*t85*t388-2.d0*t171*t401+2.d0*&
      t288*t401+4.d0*t110*t401*t14+2.d0*t92*(t117+t175-t122+t397+t11*(-t398+14.d0*&
      t40))*t14+t299*(2.d0*t117*t42+6.d0*t119*t42+t11*(-t316+t195))*t14*t388-2.d0*&
      t26*(-8.d0*t426+t205*(t428-t135)+t431-4.d0*t11*(t338-t432)+2.d0*t119*(t117+&
      t121-t397))*t34)*t359/2.d0)
      s4 = DenNc*(-2.d0*t371-t4*ThetaEps*(t10+t93)*t34/2.d0+t4*(-t378-t24*t381&
      -t5*t384)*t14*t388/4.d0+t3*GW*(-2.d0*t74*t463+2.d0*t92*t466-2.d0*t99*t466+2.d0*&
      t102*t85*t466+t475+t114*t463*t70-4.d0*t105*t466*t69-4.d0*t109*t466*t40)*&
      t145/2.d0+t3*GW*(-t61*(t10-t135)*t155+2.d0*t171*t495+2.d0*t186*t495-4.d0*t106*&
      t495*t14+2.d0*t99*(t493+t40*t136)*t14+t223*(2.d0*t207+2.d0*t75*t176+t11*(-&
      t194+t511-14.d0*t123))*t388*t139-t158*(2.d0*t10*t119+t214+t11*(-t194+t335-&
      t180))*t388*t86+2.d0*t26*(t431+4.d0*t205*t176-2.d0*t119*(t492-t494+t397))*&
      t34)*t239/2.d0)
      s2 = s3+s4
      CoeffsWm(2) = s1*s2
      s1 = ne
      s3 = Nc*(-t546-t549-t4*t372*t550/2.d0+t4*(-t556+t24*t557+t38*t560)*&
      t43*t564/4.d0+t3*GW*(2.d0*t91*t568-2.d0*t571*t568-2.d0*t574*(t175-7.d0*t281+t75)*&
      t43-4.d0*t27*t568*t43-6.d0*t299+t113*t85*(t64-t93-t40)*t564-t267*(t588-&
      t589+6.d0*t343-4.d0*t11*(t117-t134+3.d0*t123))*t43*t564+2.d0*t109*(t601+t75*t599+&
      t11*(-4.d0*t338+8.d0*t341+t605-t606))*t550)*t359/2.d0)
      s4 = DenNc*(t620+t4*ThetaEps*(t10+t135)*t550/2.d0+t3*GW*(2.d0*t5*t161*&
      t627-2.d0*t574*t630+2.d0*t99*t630-2.d0*t635*t85*t630-t640-t113*t161*t627*t70+&
      4.d0*t105*t630*t69+4.d0*t27*t630)*t239/2.d0+t4*(t556-t24*t557-t38*t560)*t43&
      *t564/4.d0+t3*GW*(t61*(t10-t93)*t663-2.d0*t91*t668-2.d0*t671*t668+4.d0*t106*&
      t668*t43+2.d0*t99*(-t175-t291+t678)*t43+t74*(t588+t208-t683+t605+t11*(-&
      t684+t335))*t86*t564-t114*(t588-t321+2.d0*t119*(t428-7.d0*t40)-t683+t605&
      +t11*(-10.d0*t117+t511-t124))*t564*t139-2.d0*t109*(t702-4.d0*t703+t706+t601&
      -2.d0*t11*(t588-5.d0*t343+2.d0*t345))*t550)*t145/2.d0)
      s2 = s3+s4
      CoeffsWm(3) = s1*s2
      CoeffsWm(4) = ne*(Nc*(-t546-t727-t549+t736/2.d0+t3*GW*(-2.d0*t738-t744-&
      t747+t749+t753+t758)/2.d0)+DenNc*(t765+t620+t3*GW*(-2.d0*t158+t767-t768+&
      t769+t113*t155+t771*t772-t776+t778)*t239/2.d0+t788/2.d0+t3*GW*(2.d0*t74-t113&
      *t663-t62*t772-t99*t795-t91*t800-t671*t800+2.d0*t106*t803-2.d0*t109*(&
      t588+t208-3.d0*t341+t807-t606-2.d0*t11*(4.d0*t117-t494+t124))*t43)*t145/2.d0))
      CoeffsWm(5) = ne*(Nc*(t546+t3*GW*t825*t33/2.d0+t727+t549+t4*(-t831+&
      t24*t832+t836)*t43*t88+t3*GW*(-2.d0*t5*t33+t744+t747-t749-t753+2.d0*t830*&
      t550-2.d0*t24*(t11-t40)*t33*t88-t758)/2.d0)+DenNc*(-t3*GW*t825*t33/2.d0-t765&
      -t620+t4*(t831-t836-t24*t832)*t43*t88+t3*GW*(t767-t768+t769+t640-&
      t776+2.d0*t158*t630*t33+t778+t868*(-t492+t379+t347-t397)*t33/t69)*t239&
      /2.d0+t3*GW*(-t475-2.d0*t74*t466*t88+t882*(t492-6.d0*t379+t175-t121)*t70*&
      t88+t91*t794+t671*t794-2.d0*t106*t795+t99*t803+2.d0*t109*(t588+t208-t341-&
      t807+t606+t11*(-t684+t347+t397))*t43)*t145/2.d0))
      s1 = ne
      s3 = t3*GW*Nc*(-t910*t281+t912*t13*t40+t24*(t224-t115))*t14*t919*&
      t88
      s4 = DenNc*(-2.d0*t923*(t47-t49-t260+t54)*t14*t926+t3*GW*(4.d0*t930*t931&
      -4.d0*t934*t358+4.d0*t937*t42*t358-8.d0*t941*t942+2.d0*t910*t13*(t10-t12-t135)&
      *t88-8.d0*t950*t85*t40+t62*(t159*t954+t152*(t194+t956-4.d0*t11*(t398-5.d0*&
      t40)-t335+t180))*t43*t966-t882*(t159*t954*t40+t152*(t588-t206+t324-&
      t589+t807+t11*(-t684+t195)))*t550*t966)*t982/2.d0+t3*GW*(-4.d0*t986*t13*&
      t85+4.d0*t990*t358-4.d0*t993*t13*t358+8.d0*t997*t942+8.d0*t1000*t11*t85+2.d0*t912&
      *t42*(-t10+t93+t41)*t33-t771*(t76*t954+t63*(t194+t289-t1010+t1011+&
      t11*(-14.d0*t10+20.d0*t40)))*t14*t966+t868*(t163*t75*t954+t6*(-8.d0*t10*&
      t205+2.d0*t119*t1026+2.d0*t42*t1030+t11*(-11.d0*t338+28.d0*t341-24.d0*t343+t1036)))*&
      t34*t966)*t1048/2.d0)
      s2 = s3+s4
      CoeffsWm(6) = s1*s2
      s1 = ne
      s3 = Nc*(-t4*ThetaEps*t6*t919+t3*GW*(6.d0*t171*t7*t13-12.d0*t26*t7*t11+6.d0&
      *t930*t42-6.d0*t287*t7*t1068-12.d0*t910*t1071+t5*(-t1074+t47*t85)*t85*&
      t388*t43-t113*(-t1074+t63*t85)*t1071*t388*t43-12.d0*t950*t40)*t359/2.d0)
      s5 = DenNc
      s7 = t923*(t117*t6-t47*t11-2.d0*t1098-2.d0*t6*t119-t1102+2.d0*t1103+4.d0*t1105&
      )*t14*t919*t70
      s9 = t3*GW*(-2.d0*t930*t1113+2.d0*t934*t1112*t85-2.d0*t937*t1113*t85+16.d0*&
      t910*t13*t466*t86+4.d0*t941*t1112*t69+4.d0*t950*t1112*t40+t114*(t159*t85-&
      t152*(t1133+20.d0*t119-28.d0*t177-t325+t180))*t43*t1140-t62*(t6*t1146+t152&
      *(t316+16.d0*t119-t1010+t180+t11*(-t1149+24.d0*t40)))*t43*t966)*t982/2.d0
      s10 = t3*GW*(-2.d0*t574*t1164+4.d0*t27*t1164+2.d0*t98*t1164*t85-2.d0*t635*&
      t1164*t85+4.d0*t105*t1164*t69+t158*(-t47*(-t289+t290+t75)+2.d0*t163*t42*t161&
      )*t388*t86-t771*(t47*t13*(t303-t93-t160)-t152*(t1133+t175+4.d0*t11*&
      t136-t1192+t1011))*t14*t966+t113*t42*(-t163*t42*(-t206+8.d0*t119*t1199+&
      6.d0*t42*t1202+t11*(-19.d0*t117+52.d0*t121-32.d0*t123))+t6*(-8.d0*t10*t426+16.d0*&
      t205*(t117-t121)+t1224+t119*(-t1225+4.d0*t341+12.d0*t343-t1036)-2.d0*t11*(t702&
      -7.d0*t338*t40+16.d0*t703-14.d0*t705+4.d0*t1234)))*t14*t388*t1140)*t1048/2.d0
      s8 = s9+s10
      s6 = s7+s8
      s4 = s5*s6
      s2 = s3+s4
      CoeffsWm(7) = s1*s2
      s1 = ne
      s3 = Nc*(t4*ThetaEps*t7*t14*t88+t3*GW*(-6.d0*t986*t13+12.d0*t1000*t11-6.d0*&
      t91*t6*t42+6.d0*t287*t6*t1068+12.d0*t912*t85*t43+12.d0*t109*t6*t40-t267*(&
      t49*t85-t1277)*t14*t564+t272*(t152*t85-t1277)*t14*t43*t564)*t359/2.d0)
      s5 = DenNc
      s7 = -t923*(t117*t7-t1098-2.d0*t1102-t49*t40+4.d0*t1103+2.d0*t1105-2.d0*t7*&
      t123)*t14*t926*t88
      s9 = t3*GW*(2.d0*t986*t1308-2.d0*t990*t1307*t85+2.d0*t993*t1308*t85-16.d0*t912&
      *t42*t630*t86-4.d0*t997*t1307*t69-4.d0*t1000*t11*t1307+t771*(t7*t213+t63&
      *(t316+t289-8.d0*t11*t1199-t1192+16.d0*t123))*t14*t966+t223*(-t76*t85+&
      t63*(t1133+t289-28.d0*t121+20.d0*t123+t11*(-18.d0*t10+28.d0*t40)))*t14*t1140)*&
      t1048/2.d0
      s10 = t3*GW*(-2.d0*t91*t1355*t42+2.d0*t98*t1355*t85-2.d0*t101*t1355*t931+4.d0*&
      t105*t1355*t69+4.d0*t109*t1355*t40+t74*(-2.d0*t63*t78*t40+t49*(t117+t175&
      +t291-t180-4.d0*t11*(t10+t40)))*t86*t564+t62*(t49*(t303-t77-t135)*t42&
      -t63*(t956+t1026+t11*(-t1149+16.d0*t40)))*t43*t966-t113*t13*(-t63*t40&
      *(t1225-12.d0*t205+t119*(30.d0*t10-32.d0*t40)-19.d0*t341+16.d0*t343-t606-4.d0*t11*(6&
      *t117-13.d0*t121+t230))+t152*(t702+4.d0*t426-6.d0*t703+t706+4.d0*t205*(-t21+&
      t40)+t119*(13.d0*t117-t347)+t11*(-t1225+t432+8.d0*t343)))*t87*t564*t966)*&
      t982/2.d0
      s8 = s9+s10
      s6 = s7+s8
      s4 = s5*s6
      s2 = s3+s4
      CoeffsWm(8) = s1*s2
      CoeffsWm(9) = ne*(DenNc*(-t371+t4*t254*t14*t33/4.d0+t3*GW*(t171+2.d0*t26&
      *(t10-t11)*t14+t186+t99*t14-2.d0*t106*t14)/2.d0+t3*GW*(t91+t671+t99*t43+&
      2.d0*t109*t176*t43-2.d0*t106*t43)/2.d0-t549-3.d0*t1465+t736/4.d0)+Nc*(t545+t367+&
      t3*t369/(t64-t93)+t4*t31*t14*t33/4.d0+t3*t547/(t64-t135)+t1465+t788/4.d0+&
      t3*GW*(t5*t14+t738-t287*(t117-t175+4.d0*t11*t40-t397)*t86+t91*(t118-&
      t175-4.d0*t11*(t10-t149)-t559+t397)*t14*t86+t171*(t118+t175-t179-t397+&
      t11*(-t428+t126))*t87-2.d0*t26*(-t206-t11*t75+4.d0*t119*t40+t213)*t1506-&
      2.d0*t109*(t338+4.d0*t207-t341-t606+4.d0*t11*(-t117+t121+t123))*t1506)/2.d0))

      end subroutine WmDecSMCoeff

      subroutine LeadingOrderWm
      implicit none
      complex(8) :: t31,t32,t33,t34,t35,t36,t37,t38,t39,t30
      complex(8) :: t21,t22,t23,t24,t25,t26,t27,t28,t29,t20
      complex(8) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
      complex(8) :: t11,t12,t13,t14,t15,t16,t17,t18,t19
      real(8) :: Mw
      Mw = M_W
            t1 = EEL*GS
      t2 = t1*GW
      t3 = SP4(1)
      t4 = SP4(2)
      t7 = Mw**2
      t8 = SP4(3)
      t11 = 1.d0/(t7-2.d0*t8)
      t15 = GW*ne
      t19 = SP4(4)
      t22 = 1.d0/(t7-2.d0*t19)
      t25 = 2.d0*t1*t15*t22
      CoeffsWmLO(1) = 2.d0*t2*ne*(t3+t4)*t11
      CoeffsWmLO(2) = -2.d0*t1*t15*t11
      CoeffsWmLO(3) = t25
      CoeffsWmLO(4) = t25
      CoeffsWmLO(5) = -t25
      CoeffsWmLO(6) = czero
      CoeffsWmLO(7) = czero
      CoeffsWmLO(8) = czero
      CoeffsWmLO(9) = -2.d0*t2*ne*(t7-t8-t19)*t11*t22

      end subroutine LeadingOrderWm

      subroutine CTgsWm
      implicit none
      complex(8) :: t31,t32,t33,t34,t35,t36,t37,t38,t39,t30
      complex(8) :: t21,t22,t23,t24,t25,t26,t27,t28,t29,t20
      complex(8) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
      complex(8) :: t11,t12,t13,t14,t15,t16,t17,t18,t19
      real(8) ::  Mw
      Mw = M_W
      t1 = EEL*GS
      t2 = t1*GW
      t3 = SP4(1)
      t4 = SP4(2)
      t7 = Mw**2
      t8 = SP4(3)
      t11 = 1.d0/(t7-2.d0*t8)
      t15 = GW*ne
      t19 = SP4(4)
      t22 = 1.d0/(t7-2.d0*t19)
      t25 = 2.d0*t1*t15*t22
      CoeffsWmCTgs(1) = 2.d0*t2*ne*(t3+t4)*t11
      CoeffsWmCTgs(2) = -2.d0*t1*t15*t11
      CoeffsWmCTgs(3) = t25
      CoeffsWmCTgs(4) = t25
      CoeffsWmCTgs(5) = -t25
      CoeffsWmCTgs(6) = czero
      CoeffsWmCTgs(7) = czero
      CoeffsWmCTgs(8) = czero
      CoeffsWmCTgs(9) = -2.d0*t2*ne*(t7-t8-t19)*t11*t22

      end subroutine CTgsWm

      subroutine CTgluonWm
      implicit none
      complex(8) :: t31,t32,t33,t34,t35,t36,t37,t38,t39,t30
      complex(8) :: t21,t22,t23,t24,t25,t26,t27,t28,t29,t20
      complex(8) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
      complex(8) :: t11,t12,t13,t14,t15,t16,t17,t18,t19
      real(8) ::  Mw
      Mw = M_W
      t1 = EEL*GS
      t2 = t1*GW
      t3 = SP4(1)
      t4 = SP4(2)
      t7 = Mw**2
      t8 = SP4(3)
      t11 = 1.d0/(t7-2.d0*t8)
      t14 = GW*ne
      t17 = SP4(4)
      t20 = 1.d0/(t7-2.d0*t17)
      t22 = t1*t14*t20
      CoeffsWmCTgluon(1) = t2*ne*(t3+t4)*t11
      CoeffsWmCTgluon(2) = -t1*t14*t11
      CoeffsWmCTgluon(3) = t22
      CoeffsWmCTgluon(4) = t22
      CoeffsWmCTgluon(5) = -t22
      CoeffsWmCTgluon(6) = czero
      CoeffsWmCTgluon(7) = czero
      CoeffsWmCTgluon(8) = czero
      CoeffsWmCTgluon(9) = -t2*ne*(t7-t8-t17)*t11*t20

      end subroutine CTgluonWm

subroutine CalcCTsWm(xe,mu, CT)
implicit none
integer :: xe
real(8) :: dZMwm,dZgsWm,dZgluonWm, mu, Ca, Cf
complex(8) :: amp_weyl_ctgsWm(1:4),amp_weyl_ctgluonWm(1:4), CT(1:4)
dZgsWm = 0.d0
dZgluonWm = 0.d0


Ca = Nc
Cf = (Nc**2-1.d0)/2.d0/Nc

CT(1:4) = czero

  if(xe .eq. 0 .or. xe .eq. -1) then
     call CTgluonWm
     call CTgsWm
     if(xe .eq. 0 ) then
        dZgsWm = alpha_S_scharf/12.d0/Pi*dlog(mu**2/m_w**2)
        dZgluonWm = -alpha_S_scharf/6.d0/Pi*dlog(mu**2/m_w**2)
     elseif(xe .eq. -1 ) then
        dZgsWm    = alpha_S_scharf/4.d0/Pi*(N_f/3.d0-11.d0/2.d0) + alpha_S_scharf/12.d0/Pi
        dZgluonWm = -alpha_S_scharf/6.d0/Pi
     endif

     amp_weyl_ctgsWm(1:4) = dZgsWm* &
          ( CoeffswmCtgs(1)*SME(1,1:4) + CoeffswmCtgs(2)*SME(2,1:4) + CoeffswmCtgs(3)*SME(3,1:4) + CoeffswmCtgs(4)*SME(4,1:4) &
          + CoeffswmCtgs(5)*SME(5,1:4) + CoeffswmCtgs(6)*SME(6,1:4) + CoeffswmCtgs(7)*SME(7,1:4) + CoeffswmCtgs(8)*SME(8,1:4) + CoeffswmCtgs(9)*SME(9,1:4))

     amp_weyl_ctgluonWm(1:4) = dZgluonWm* &
          ( CoeffswmCTgluon(1)*SME(1,1:4) + CoeffswmCTgluon(2)*SME(2,1:4) + CoeffswmCTgluon(3)*SME(3,1:4) + CoeffswmCTgluon(4)*SME(4,1:4) &
          + CoeffswmCTgluon(5)*SME(5,1:4) + CoeffswmCTgluon(6)*SME(6,1:4) + CoeffswmCTgluon(7)*SME(7,1:4) + CoeffswmCTgluon(8)*SME(8,1:4)+ CoeffswmCtgluon(9)*SME(9,1:4))

     CT(1:4) = +amp_weyl_ctgluonWm(1:4) + amp_weyl_ctgsWm(1:4)
  end if


! write(*,*) N_F
! stop
end subroutine CalcCTswm



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! STANDARD-MATRIX-ELEMENTS


function weyl_sme0(sp1,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4)
complex(8) :: weyl_sme0
weyl_sme0 = psp1_(sp1,sp2)
return
end function weyl_sme0


function weyl_sme1(sp1,v1,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4)
complex(8) :: weyl_sme1
sp9(1:4) =  spb2_weyl(sp1,v1)
weyl_sme1 = psp1_(sp9,sp2)
return
end function weyl_sme1


function weyl_sme2(sp1,v1,v2,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4)
complex(8) :: weyl_sme2
sp9(1:4) =  spb2_weyl(sp1,v1)
sp7(1:4) =  spb2_weyl(sp9,v2)
weyl_sme2 = psp1_(sp7,sp2)
return
end function weyl_sme2


function weyl_sme3(sp1,v1,v2,v3,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4), v3(1:4)
complex(8) :: weyl_sme3
sp9(1:4) =  spb2_weyl(sp1,v1)
sp7(1:4) =  spb2_weyl(sp9,v2)
sp9(1:4) =  spb2_weyl(sp7,v3)
weyl_sme3 = psp1_(sp9,sp2)
return
end function weyl_sme3



function weyl_sme4(sp1,v1,v2,v3,v4,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4), v3(1:4), v4(1:4)
complex(8) :: weyl_sme4
sp9(1:4) =  spb2_weyl(sp1,v1)
sp7(1:4) =  spb2_weyl(sp9,v2)
sp9(1:4) =  spb2_weyl(sp7,v3)
sp7(1:4) =  spb2_weyl(sp9,v4)
weyl_sme4 = psp1_(sp7,sp2)
return
end function weyl_sme4


function weyl_sme5(sp1,v1,v2,v3,v4,v5,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4), v3(1:4), v4(1:4), v5(1:4)
complex(8) :: weyl_sme5
sp9(1:4) =  spb2_weyl(sp1,v1)
sp7(1:4) =  spb2_weyl(sp9,v2)
sp9(1:4) =  spb2_weyl(sp7,v3)
sp7(1:4) =  spb2_weyl(sp9,v4)
sp9(1:4) =  spb2_weyl(sp7,v5)
weyl_sme5 = psp1_(sp9,sp2)
return
end function weyl_sme5


function weyl_sme6(sp1,v1,v2,v3,v4,v5,v6,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4), v3(1:4), v4(1:4), v5(1:4), v6(1:4)
complex(8) :: weyl_sme6
sp9(1:4) =  spb2_weyl(sp1,v1)
sp7(1:4) =  spb2_weyl(sp9,v2)
sp9(1:4) =  spb2_weyl(sp7,v3)
sp7(1:4) =  spb2_weyl(sp9,v4)
sp9(1:4) =  spb2_weyl(sp7,v5)
sp7(1:4) =  spb2_weyl(sp9,v6)
weyl_sme6 = psp1_(sp7,sp2)
return
end function weyl_sme6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


function weyl_ime0(sp1,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4)
complex(8) :: weyl_ime0
weyl_ime0 = psp1_(sp1,sp2)
return
end function weyl_ime0


function weyl_ime1(sp1,v1,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4)
complex(8) :: weyl_ime1
sp9(1:4) =  spi2_weyl(v1,sp2)
weyl_ime1 = psp1_(sp9,sp1)
return
end function weyl_ime1


function weyl_ime2(sp1,v1,v2,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4)
complex(8) :: weyl_ime2
sp9(1:4) =  spi2_weyl(v1,sp2)
sp7(1:4) =  spi2_weyl(v2,sp9)
weyl_ime2 = psp1_(sp7,sp1)
return
end function weyl_ime2


function weyl_ime3(sp1,v1,v2,v3,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4), v3(1:4)
complex(8) :: weyl_ime3
sp9(1:4) =  spi2_weyl(v1,sp2)
sp7(1:4) =  spi2_weyl(v2,sp9)
sp9(1:4) =  spi2_weyl(v3,sp7)
weyl_ime3 = psp1_(sp9,sp1)
return
end function weyl_ime3



function weyl_ime4(sp1,v1,v2,v3,v4,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4), v3(1:4), v4(1:4)
complex(8) :: weyl_ime4
sp9(1:4) =  spi2_weyl(v1,sp2)
sp7(1:4) =  spi2_weyl(v2,sp9)
sp9(1:4) =  spi2_weyl(v3,sp7)
sp7(1:4) =  spi2_weyl(v4,sp9)
weyl_ime4 = psp1_(sp7,sp1)
return
end function weyl_ime4


function weyl_ime5(sp1,v1,v2,v3,v4,v5,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4), v3(1:4), v4(1:4), v5(1:4)
complex(8) :: weyl_ime5
sp9(1:4) =  spi2_weyl(v1,sp2)
sp7(1:4) =  spi2_weyl(v2,sp9)
sp9(1:4) =  spi2_weyl(v3,sp7)
sp7(1:4) =  spi2_weyl(v4,sp9)
sp9(1:4) =  spi2_weyl(v5,sp7)
weyl_ime5 = psp1_(sp9,sp1)
return
end function weyl_ime5


function weyl_ime6(sp1,v1,v2,v3,v4,v5,v6,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4), v3(1:4), v4(1:4), v5(1:4), v6(1:4)
complex(8) :: weyl_ime6
sp9(1:4) =  spi2_weyl(v1,sp2)
sp7(1:4) =  spi2_weyl(v2,sp9)
sp9(1:4) =  spi2_weyl(v3,sp7)
sp7(1:4) =  spi2_weyl(v4,sp9)
sp9(1:4) =  spi2_weyl(v5,sp7)
sp7(1:4) =  spi2_weyl(v6,sp9)
weyl_ime6 = psp1_(sp7,sp1)
end function weyl_ime6


!!!SPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINOR!!! SPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINOR!!!
!!!SPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINOR!!! SPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINOR!!!
!!!SPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINOR!!! SPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINOR!!!
!!!SPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINOR!!! SPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINORSPINOR!!!



function weyl_spin0(sp1)
implicit none
complex(8) :: sp1(1:4)
complex(8) :: weyl_spin0(1:4)
weyl_spin0(1:4)= sp1(1:4)
return
end function weyl_spin0


function weyl_spin1(sp1,v1)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4)
complex(8) :: weyl_spin1(1:4)
sp9(1:4) =  spb2_weyl(sp1,v1)
weyl_spin1(1:4)= sp9(1:4)
return
end function weyl_spin1


function weyl_spin2(sp1,v1,v2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4)
complex(8) :: weyl_spin2(1:4)
sp9(1:4) =  spb2_weyl(sp1,v1)
sp7(1:4) =  spb2_weyl(sp9,v2)
weyl_spin2(1:4)= sp7(1:4)
return
end function weyl_spin2


function weyl_spin3(sp1,v1,v2,v3)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4), v3(1:4)
complex(8) :: weyl_spin3(1:4)
sp9(1:4) =  spb2_weyl(sp1,v1)
sp7(1:4) =  spb2_weyl(sp9,v2)
sp9(1:4) =  spb2_weyl(sp7,v3)
weyl_spin3(1:4)= sp9(1:4)
return
end function weyl_spin3



function weyl_spin4(sp1,v1,v2,v3,v4)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4), v3(1:4), v4(1:4)
complex(8) :: weyl_spin4(1:4)
sp9(1:4) =  spb2_weyl(sp1,v1)
sp7(1:4) =  spb2_weyl(sp9,v2)
sp9(1:4) =  spb2_weyl(sp7,v3)
sp7(1:4) =  spb2_weyl(sp9,v4)
weyl_spin4(1:4)= sp7(1:4)
return
end function weyl_spin4


function weyl_spin5(sp1,v1,v2,v3,v4,v5)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4), v3(1:4), v4(1:4), v5(1:4)
complex(8) :: weyl_spin5(1:4)
sp9(1:4) =  spb2_weyl(sp1,v1)
sp7(1:4) =  spb2_weyl(sp9,v2)
sp9(1:4) =  spb2_weyl(sp7,v3)
sp7(1:4) =  spb2_weyl(sp9,v4)
sp9(1:4) =  spb2_weyl(sp7,v5)
weyl_spin5(1:4)= sp9(1:4)
return
end function weyl_spin5


function weyl_spin6(sp1,v1,v2,v3,v4,v5,v6)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4), v3(1:4), v4(1:4), v5(1:4), v6(1:4)
complex(8) :: weyl_spin6(1:4)
sp9(1:4) =  spb2_weyl(sp1,v1)
sp7(1:4) =  spb2_weyl(sp9,v2)
sp9(1:4) =  spb2_weyl(sp7,v3)
sp7(1:4) =  spb2_weyl(sp9,v4)
sp9(1:4) =  spb2_weyl(sp7,v5)
sp7(1:4) =  spb2_weyl(sp9,v6)
weyl_spin6(1:4)= sp7(1:4)
return
end function weyl_spin6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




function weyl_ispin1(v1,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4)
complex(8) :: weyl_ispin1(1:4)
sp9(1:4) =  spi2_weyl(v1,sp2)
weyl_ispin1(1:4) = sp9(1:4)
return
end function weyl_ispin1


function weyl_ispin2(v1,v2,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4)
complex(8) :: weyl_ispin2(1:4)
sp9(1:4) =  spi2_weyl(v2,sp2)
sp7(1:4) =  spi2_weyl(v1,sp9)
weyl_ispin2(1:4) = sp7(1:4)
return
end function weyl_ispin2


function weyl_ispin3(v1,v2,v3,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4), v3(1:4)
complex(8) :: weyl_ispin3(1:4)
sp9(1:4) =  spi2_weyl(v3,sp2)
sp7(1:4) =  spi2_weyl(v2,sp9)
sp9(1:4) =  spi2_weyl(v1,sp7)
weyl_ispin3(1:4) = sp9(1:4)
return
end function weyl_ispin3



function weyl_ispin4(v1,v2,v3,v4,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4), v3(1:4), v4(1:4)
complex(8) :: weyl_ispin4(1:4)
sp9(1:4) =  spi2_weyl(v4,sp2)
sp7(1:4) =  spi2_weyl(v3,sp9)
sp9(1:4) =  spi2_weyl(v2,sp7)
sp7(1:4) =  spi2_weyl(v1,sp9)
weyl_ispin4(1:4) = sp7(1:4)
return
end function weyl_ispin4


function weyl_ispin5(v1,v2,v3,v4,v5,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4), v3(1:4), v4(1:4), v5(1:4)
complex(8) :: weyl_ispin5(1:4)
sp9(1:4) =  spi2_weyl(v5,sp2)
sp7(1:4) =  spi2_weyl(v4,sp9)
sp9(1:4) =  spi2_weyl(v3,sp7)
sp7(1:4) =  spi2_weyl(v2,sp9)
sp9(1:4) =  spi2_weyl(v1,sp7)
weyl_ispin5(1:4) = sp9(1:4)
return
end function weyl_ispin5


function weyl_ispin6(v1,v2,v3,v4,v5,v6,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4), v3(1:4), v4(1:4), v5(1:4), v6(1:4)
complex(8) :: weyl_ispin6(1:4)
sp9(1:4) =  spi2_weyl(v1,sp2)
sp7(1:4) =  spi2_weyl(v2,sp9)
sp9(1:4) =  spi2_weyl(v3,sp7)
sp7(1:4) =  spi2_weyl(v4,sp9)
sp9(1:4) =  spi2_weyl(v5,sp7)
sp7(1:4) =  spi2_weyl(v6,sp9)
weyl_ispin6(1:4) = sp7(1:4)
end function weyl_ispin6

END MODULE
