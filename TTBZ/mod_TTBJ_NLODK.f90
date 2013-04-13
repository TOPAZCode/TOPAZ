Module ModTTBJ_NLODK
use ModMisc
use ModParameters

complex(8),private :: p1t(1:4),k1b(1:4),k2w(1:4),k3g(1:4), p1tb(1:4),k1bb(1:4)
complex(8),private :: EpsS3g(1:4), EpsS2w(1:4), Spi1b(1:4), Spi1t(1:4), Spi1bb(1:4), Spi1tb(1:4)
complex(8),private :: Integrals(1:18)
complex(8),private :: SME(1:8,1:4), CoeffsTopLO(1:8), CoeffsTopCTmass(1:8), CoeffsTopCTwave(1:8), CoeffsTopCTgs(1:8), CoeffsTopCTgluon(1:8)
complex(8),private :: SMEA(1:8,1:4), CoeffsAtopLO(1:8), CoeffsAtopCTmass(1:8), CoeffsAtopCTwave(1:8), CoeffsAtopCTgs(1:8), CoeffsAtopCTgluon(1:8)
complex(8),private :: SP4(1:7)
real(8),private ::  Nc, DenNc
real(8),private :: Mw
real(8),private :: mt,GS
real(8),private,parameter :: N_f=Nf_light, GW=1d0,EEL=1d0


contains


! changes: module name, mt mw declarations, privates, couplings +nf, arg. order neu--lep, w pol vec., removed CF, xe if


subroutine calc_topdecay(pt,kb,klep,kneu,kg,pol_g,pol_w,xe,mu,lo_nlo,gau,Spinor)
implicit none
real(8):: pt(1:4),kb(1:4),klep(1:4),kneu(1:4),kg(1:4), mu,ThetaEps, gram_det
complex(8) :: Spinor(1:4), CT(1:4), CTMass(1:4), weyl(1:4), CoeffsTopNLO(1:8),pol_w(1:4)
integer :: pol_g, lo_nlo, gau, xe, i
Spinor(1:4) = czero
Nc =3.d0
DenNc =1.d0/3.d0
call qlinit
call ffini


    Mw = M_w
    mt = M_top
    GS=dsqrt(4d0*DblPi*alpha_s*RunAlphaS(NLOParam,MuRen))


ThetaEps =0.d0
if(xe .eq. 0) then
   ThetaEps = 1.d0
endif

call SPandPoltop(pt,kb,kneu,klep,kg,pol_g,pol_w, gau)
gram_det = dabs(dble(((k1b.dot.k2w)*(m_top**2 - 2.d0*(k1b.dot.p1t)) -  M_w**2*(k1b.dot.p1t))/m_top**4))
if(LO_NLO.eq.0 .and. xe.eq.0) then
   call CalcSMEtop
   call LeadingOrdertop
   weyl(1:4) = CoeffsTopLO(1)*SME(1,1:4) + CoeffsTopLO(2)*SME(2,1:4) + CoeffsTopLO(3)*SME(3,1:4) + CoeffsTopLO(4)*SME(4,1:4) &
             + CoeffsTopLO(5)*SME(5,1:4) + CoeffsTopLO(6)*SME(6,1:4) + CoeffsTopLO(7)*SME(7,1:4) + CoeffsTopLO(8)*SME(8,1:4)
   Spinor(1:4) = (spb2_(weyltodirac(weyl(1:4)),p1t(1:4))+m_top*weyltodirac(weyl(1:4)))*ne  *dsqrt(4.d0/3.d0) ! COLOR-FACTOR
elseif(LO_NLO .eq.  1 .and. gram_det .gt. 1.d-4) then
   call TopIntegralsTTBJ(xe,mu)
   call CalcSMETop
   call TopDecSMCoeff(SP4,mt,Mw,Nc,DenNc,ThetaEps,Integrals,CoeffsTopNLO)
   call CalcCTsTop(xe, mu, CT, CTmass)
   weyl(1:4) = (CoeffsTopNLO(1)*SME(1,1:4) + CoeffsTopNLO(2)*SME(2,1:4) + CoeffsTopNLO(3)*SME(3,1:4) + CoeffsTopNLO(4)*SME(4,1:4) &
              + CoeffsTopNLO(5)*SME(5,1:4) + CoeffsTopNLO(6)*SME(6,1:4) + CoeffsTopNLO(7)*SME(7,1:4) + CoeffsTopNLO(8)*SME(8,1:4))/(4.d0*Pi)**2*ne*EEL*GW*GS**3 &
       + CTmass(1:4) + Ct(1:4)

   Spinor(1:4) = (spb2_(weyltodirac(weyl(1:4)),p1t(1:4))+m_top*weyltodirac(weyl(1:4)))*ne  *dsqrt(4.d0/3.d0)! COLOR-FACTOR
else
   Spinor(1:4) = (0d0,0d0)
endif

return
end subroutine calc_topdecay


subroutine calc_atopdecay(pt,kb,klep,kneu,kg,pol_g,pol_w,xe,mu,lo_nlo,gau,Spinor)
implicit none
real(8):: pt(1:4),kb(1:4),klep(1:4),kneu(1:4),kg(1:4), mu,ThetaEps, gram_det
complex(8) :: Spinor(1:4), CT(1:4), CTMass(1:4), weyl(1:4), CoeffsAtopNLO(1:8),pol_w(1:4)
integer :: pol_g, lo_nlo, gau, xe, i
Spinor(1:4) = czero
Nc =3.d0
DenNc =1.d0/3.d0
call qlinit
call ffini

mw=m_w
mt=m_top
GS=dsqrt(4d0*DblPi*alpha_s*RunAlphaS(NLOParam,MuRen))

ThetaEps =0.d0
if(xe .eq. 0) then
   ThetaEps = 1.d0
endif

call SPandPolAtop(pt,kb,kneu,klep,kg,pol_g,pol_w, gau)
gram_det = dabs(dble(((k1bb.dot.k2w)*(m_top**2 - 2.d0*(k1bb.dot.p1tb)) -  M_w**2*(k1bb.dot.p1tb))/m_top**4))
if(LO_NLO.eq.0 .and. xe.eq.0) then
   call CalcSMEAtop
   call LeadingOrderAtop
   weyl(1:4) = CoeffsAtopLO(1)*SMEA(1,1:4) + CoeffsAtopLO(2)*SMEA(2,1:4) + CoeffsAtopLO(3)*SMEA(3,1:4) + CoeffsAtopLO(4)*SMEA(4,1:4) &
             + CoeffsAtopLO(5)*SMEA(5,1:4) + CoeffsAtopLO(6)*SMEA(6,1:4) + CoeffsAtopLO(7)*SMEA(7,1:4) + CoeffsAtopLO(8)*SMEA(8,1:4)
   Spinor(1:4) = (-spi2_(p1tb(1:4),weyltodirac(weyl(1:4)))+m_top*weyltodirac(weyl(1:4)))*ne  *dsqrt(4.d0/3.d0)! COLOR-FACTOR
elseif(LO_NLO .eq.  1 .and. gram_det .gt. 1.d-4) then
   call AtopIntegralsTTBJ(xe,mu)
   call CalcSMEAtop
   call AtopDecSMCoeff(SP4,mt,Mw,Nc,DenNc,ThetaEps,Integrals,CoeffsATopNLO)
   call CalcCTsAtop(xe, mu, CT, CTmass)
   weyl(1:4) = (CoeffsAtopNLO(1)*SMEA(1,1:4) + CoeffsAtopNLO(2)*SMEA(2,1:4) + CoeffsAtopNLO(3)*SMEA(3,1:4) + CoeffsAtopNLO(4)*SMEA(4,1:4) &
              + CoeffsAtopNLO(5)*SMEA(5,1:4) + CoeffsAtopNLO(6)*SMEA(6,1:4) + CoeffsAtopNLO(7)*SMEA(7,1:4) + CoeffsAtopNLO(8)*SMEA(8,1:4))/(4.d0*Pi)**2*ne*EEL*GW*GS**3 &
       + CTmass(1:4) + Ct(1:4)
   Spinor(1:4) = (-spi2_(p1tb(1:4),weyltodirac(weyl(1:4)))+m_top*weyltodirac(weyl(1:4)))*ne  *dsqrt(4.d0/3.d0)! COLOR-FACTOR
else
   Spinor(1:4) = (0d0,0d0)
endif
return
end subroutine calc_atopdecay





! subroutine calc_topdecay(pt,kb,klep,kneu,kg,pol_g,pol_w,xe,mu,lo_nlo,gau,Spinor)
! implicit none
! real(8):: pt(1:4),kb(1:4),klep(1:4),kneu(1:4),kg(1:4), mu,ThetaEps
! complex(8) :: Spinor(1:4), CT(1:4), CTMass(1:4), weyl(1:4), CoeffsTopNLO(1:8),pol_w(1:4)
! integer :: pol_g, lo_nlo, gau, xe, i
! Spinor(1:4) = zero
! Nc =3.d0
! DenNc =1.d0/3.d0
! call qlinit
! call ffini
!
!
!     Mw = M_w
!     mt = M_top
!     GS=dsqrt(4d0*DblPi*alpha_s*RunAlphaS(NLOParam,MuRen))
!
! ThetaEps =0.d0
! if(xe .eq. 0) then
!    ThetaEps = 1.d0
! endif
!
! call SPandPoltop(pt,kb,kneu,klep,kg,pol_g,pol_w, gau)
! if(LO_NLO.eq.0 .and. xe.eq.0) then
!    call CalcSMEtop
!    call LeadingOrdertop
!    weyl(1:4) = CoeffsTopLO(1)*SME(1,1:4) + CoeffsTopLO(2)*SME(2,1:4) + CoeffsTopLO(3)*SME(3,1:4) + CoeffsTopLO(4)*SME(4,1:4) &
!              + CoeffsTopLO(5)*SME(5,1:4) + CoeffsTopLO(6)*SME(6,1:4) + CoeffsTopLO(7)*SME(7,1:4) + CoeffsTopLO(8)*SME(8,1:4)
!    Spinor(1:4) = (spb2_(weyltodirac(weyl(1:4)),p1t(1:4))+m_top*weyltodirac(weyl(1:4)))*ne  !*dsqrt(4.d0/3.d0)
! elseif(LO_NLO .eq.  1) then
!    call TopIntegralsTTBJ(xe,mu)
!    call CalcSMETop
!    call TopDecSMCoeff(SP4,mt,Mw,Nc,DenNc,ThetaEps,Integrals,CoeffsTopNLO)
!    call CalcCTsTop(xe, mu, CT, CTmass)
!    weyl(1:4) = (CoeffsTopNLO(1)*SME(1,1:4) + CoeffsTopNLO(2)*SME(2,1:4) + CoeffsTopNLO(3)*SME(3,1:4) + CoeffsTopNLO(4)*SME(4,1:4) &
!               + CoeffsTopNLO(5)*SME(5,1:4) + CoeffsTopNLO(6)*SME(6,1:4) + CoeffsTopNLO(7)*SME(7,1:4) + CoeffsTopNLO(8)*SME(8,1:4))/(4.d0*Pi)**2*ne*EEL*GW*GS**3 &
!        + CTmass(1:4) + Ct(1:4)
!
!    Spinor(1:4) = (spb2_(weyltodirac(weyl(1:4)),p1t(1:4))+m_top*weyltodirac(weyl(1:4)))*ne  !*dsqrt(4.d0/3.d0)
! else
!    Spinor(1:4) = (0d0,0d0)
! endif
! return
! end subroutine calc_topdecay
!
!
! subroutine calc_atopdecay(pt,kb,klep,kneu,kg,pol_g,pol_w,xe,mu,lo_nlo,gau,Spinor)
! implicit none
! real(8):: pt(1:4),kb(1:4),klep(1:4),kneu(1:4),kg(1:4), mu,ThetaEps
! complex(8) :: Spinor(1:4), CT(1:4), CTMass(1:4), weyl(1:4), CoeffsAtopNLO(1:8),pol_w(1:4)
! integer :: pol_g, lo_nlo, gau, xe, i
! Spinor(1:4) = zero
! Nc =3.d0
! DenNc =1.d0/3.d0
! call qlinit
! call ffini
!
! mw=m_w
! mt=m_top
! GS=dsqrt(4d0*DblPi*alpha_s*RunAlphaS(NLOParam,MuRen))
!
! ThetaEps =0.d0
! if(xe .eq. 0) then
!    ThetaEps = 1.d0
! endif
!
! call SPandPolAtop(pt,kb,kneu,klep,kg,pol_g,pol_w, gau)
! if(LO_NLO.eq.0 .and. xe.eq.0) then
!    call CalcSMEAtop
!    call LeadingOrderAtop
!    weyl(1:4) = CoeffsAtopLO(1)*SMEA(1,1:4) + CoeffsAtopLO(2)*SMEA(2,1:4) + CoeffsAtopLO(3)*SMEA(3,1:4) + CoeffsAtopLO(4)*SMEA(4,1:4) &
!              + CoeffsAtopLO(5)*SMEA(5,1:4) + CoeffsAtopLO(6)*SMEA(6,1:4) + CoeffsAtopLO(7)*SMEA(7,1:4) + CoeffsAtopLO(8)*SMEA(8,1:4)
!    Spinor(1:4) = (-spi2_(p1tb(1:4),weyltodirac(weyl(1:4)))+m_top*weyltodirac(weyl(1:4)))*ne  !*dsqrt(4.d0/3.d0)
! elseif(LO_NLO .eq.  1) then
!    call AtopIntegralsTTBJ(xe,mu)
!    call CalcSMEAtop
!    call AtopDecSMCoeff(SP4,mt,Mw,Nc,DenNc,ThetaEps,Integrals,CoeffsATopNLO)
!    call CalcCTsAtop(xe, mu, CT, CTmass)
!    weyl(1:4) = (CoeffsAtopNLO(1)*SMEA(1,1:4) + CoeffsAtopNLO(2)*SMEA(2,1:4) + CoeffsAtopNLO(3)*SMEA(3,1:4) + CoeffsAtopNLO(4)*SMEA(4,1:4) &
!               + CoeffsAtopNLO(5)*SMEA(5,1:4) + CoeffsAtopNLO(6)*SMEA(6,1:4) + CoeffsAtopNLO(7)*SMEA(7,1:4) + CoeffsAtopNLO(8)*SMEA(8,1:4))/(4.d0*Pi)**2*ne*EEL*GW*GS**3 &
!        + CTmass(1:4) + Ct(1:4)
!    Spinor(1:4) = (-spi2_(p1tb(1:4),weyltodirac(weyl(1:4)))+m_top*weyltodirac(weyl(1:4)))*ne  !*dsqrt(4.d0/3.d0)
! else
!    Spinor(1:4) = (0d0,0d0)
! endif
! return
! end subroutine calc_atopdecay




subroutine SPandPOLAtop(pt,kb,kneu,klep,kg,pol_g,pol_w,gau)
implicit none
integer :: gau,pol_g, i
real(8) ::  pt(1:4),kb(1:4),klep(1:4),kneu(1:4),kg(1:4)
complex(8) :: pol_w(1:4)
do i =1,4
   p1tb(i) = dcmplx(pt(i), 0.d0)
   k1bb(i) = dcmplx(kb(i),0.d0)
   k2w(i) = dcmplx(klep(i) + kneu(i),0.d0)
   k3g(i) = dcmplx(kg(i),0.d0)
enddo

if(gau .eq. 1) then
   EpsS3g(1:4)  = k3g(1:4)
else
   call pol_mless(k3g,pol_g,EpsS3g)
endif
!WPOL
! pol_w =-1
! EpsS2w(1:4) = pol_mass(k2w,M_w,pol_w)
EpsS2w(1:4) = pol_w(1:4)
call vSpi_Weyl(k1bb,+1,Spi1bb)

SP4(1) = EpsS2w.dot.k1bb
SP4(2) = EpsS2w.dot.p1tb
SP4(3) = EpsS3g.dot.EpsS2w
SP4(4) = EpsS3g.dot.k1bb
SP4(5) = EpsS3g.dot.k2w
SP4(6) = k1bb.dot.k2w
SP4(7) = k1bb.dot.p1tb

end subroutine SPandPOLAtop


subroutine CalcSMEAtop
implicit none
complex(8) :: help(1:4)
integer i

do i=1,4
   SMEA(1,i) = Spi1bb(i)
enddo

help =  weyl_ispin1(EpsS2w,Spi1bb)
do i=1,4
   SMEA(2,i) = help(i)
enddo

help = weyl_ispin1(k2w,Spi1bb) - m_top*Spi1bb(1:4)
do i=1,4
   SMEA(3,i) = help(i)
enddo

help = weyl_ispin1(EpsS3g,Spi1bb)
do i=1,4
   SMEA(4,i) = help(i)
enddo

help = weyl_ispin2(k2w,EpsS2w,Spi1bb)
do i=1,4
   SMEA(5,i) = help(i)
enddo

help = weyl_ispin2(EpsS3g,EpsS2w,Spi1bb)
do i=1,4
   SMEA(6,i) = help(i)
enddo

help = weyl_ispin2(k2w,EpsS3g,Spi1bb)
do i=1,4
   SMEA(7,i) = help(i)
enddo

help = weyl_ispin3(k2w,EpsS3g,EpsS2w,Spi1bb) - m_top*SMEA(6,1:4)
do i=1,4
   SMEA(8,i) = help(i)
enddo

end subroutine CalcSMEAtop

subroutine AtopIntegralsTTBJ(xe,muRF)
implicit none
real(8) :: mzero,muRF, p1q,p2q,p3q, p4q, m1q,m2q,m3q,m4q, p12q,p23q, muq
complex(8) :: qlI1,qlI2, qlI3, qlI4
integer :: xe
muq =muRF**2
mzero = 0.d0
m1q = mt**2
       Integrals(1) = qlI1(m1q,muq,xe)

p1q = dble(Mw**2 + 2*SP4(6))
       m1q = mt**2
       m2q = mzero
       Integrals(2) = qlI2(p1q,m1q,m2q,muq,xe)

p1q = dble(Mw**2)
       m1q = mt**2
       m2q = mzero
       Integrals(3) = qlI2(p1q,m1q,m2q,muq,xe)

p1q = dble(-2*(SP4(6) - SP4(7)))
       m1q = mzero
       m2q = mzero
       Integrals(4) = qlI2(p1q,m1q,m2q,muq,xe)

p1q = dble(mt**2 - 2*SP4(7))
       m1q = mt**2
       m2q = mzero
       Integrals(5) = qlI2(p1q,m1q,m2q,muq,xe)

p1q = dble(mt**2)
       m1q = mt**2
       m2q = mzero
       Integrals(6) = qlI2(p1q,m1q,m2q,muq,xe)

p1q = dble(mzero)
       p2q = dble(mt**2 - 2*SP4(7))
       p3q = dble(mt**2)
       m1q = mzero
       m2q = mzero
       m3q = mt**2
       Integrals(7) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(mzero)
       p2q = dble(mzero)
       p3q = dble(-2*(SP4(6) - SP4(7)))
       m1q = mzero
       m2q = mzero
       m3q = mzero
       Integrals(8) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(Mw**2 + 2*SP4(6))
       p2q = dble(mzero)
       p3q = dble(mt**2)
       m1q = mzero
       m2q = mt**2
       m3q = mt**2
       Integrals(9) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(Mw**2 + 2*SP4(6))
       p2q = dble(mzero)
       p3q = dble(mt**2)
       m1q = mt**2
       m2q = mzero
       m3q = mzero
       Integrals(10) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(Mw**2)
       p2q = dble(mzero)
       p3q = dble(Mw**2 + 2*SP4(6))
       m1q = mt**2
       m2q = mzero
       m3q = mzero
       Integrals(11) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(Mw**2)
       p2q = dble(-2*(SP4(6) - SP4(7)))
       p3q = dble(mt**2)
       m1q = mt**2
       m2q = mzero
       m3q = mzero
       Integrals(12) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(Mw**2)
       p2q = dble(mzero)
       p3q = dble(mt**2 - 2*SP4(7))
       m1q = mzero
       m2q = mt**2
       m3q = mt**2
       Integrals(13) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(mzero)
       p2q = dble(Mw**2 + 2*SP4(6))
       p3q = dble(mt**2)
       m1q = mt**2
       m2q = mt**2
       m3q = mzero
       Integrals(14) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(mzero)
       p2q = dble(Mw**2)
       p3q = dble(mt**2 - 2*SP4(7))
       m1q = mzero
       m2q = mzero
       m3q = mt**2
       Integrals(15) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(mzero)
       p2q = dble(Mw**2)
       p3q = dble(mzero)
       p4q = dble(mt**2)
       p12q = dble(Mw**2 + 2*SP4(6))
       p23q = dble(mt**2 - 2*SP4(7))

       m1q = mzero
       m2q = mzero
       m3q = mt**2
       m4q = mt**2

       Integrals(16) = qlI4(p1q,p2q,p3q,p4q,p12q,p23q,m1q,m2q,m3q,m4q,muq,xe)

p1q = dble(mzero)
       p2q = dble(mzero)
       p3q = dble(Mw**2)
       p4q = dble(mt**2)
       p12q = dble(-2*(SP4(6) - SP4(7)))
       p23q = dble(mt**2 - 2*SP4(7))

       m1q = mzero
       m2q = mzero
       m3q = mzero
       m4q = mt**2

       Integrals(17) = qlI4(p1q,p2q,p3q,p4q,p12q,p23q,m1q,m2q,m3q,m4q,muq,xe)

p1q = dble(Mw**2)
       p2q = dble(mzero)
       p3q = dble(mzero)
       p4q = dble(mt**2)
       p12q = dble(Mw**2 + 2*SP4(6))
       p23q = dble(-2*(SP4(6) - SP4(7)))

       m1q = mt**2
       m2q = mzero
       m3q = mzero
       m4q = mzero

       Integrals(18) = qlI4(p1q,p2q,p3q,p4q,p12q,p23q,m1q,m2q,m3q,m4q,muq,xe)


end subroutine AtopIntegralsTTBJ


subroutine SPandPOLTop(pt,kb,kneu,klep,kg,pol_g,pol_w,gau)
implicit none
integer :: pol_g, gau, i
real(8) ::  pt(1:4),kb(1:4),klep(1:4),kneu(1:4),kg(1:4)
complex(8) :: pol_w(1:4)
do i =1,4
   p1t(i) = dcmplx(pt(i), 0.d0)
   k1b(i) = dcmplx(kb(i),0.d0)
   k2w(i) = dcmplx(klep(i) + kneu(i),0.d0)
   k3g(i) = dcmplx(kg(i),0.d0)
enddo

if(gau .eq. 1) then
   EpsS3g(1:4)  = k3g(1:4)
   write(*,*) "Checking gauge invariance"
else
   call pol_mless(k3g,pol_g,EpsS3g)
endif
!WPOL1
! pol_W = -1
! EpsS2w(1:4) = pol_mass(k2w,M_w,pol_w)
EpsS2w(1:4) = pol_w(1:4)
call ubarSpi_Weyl(k1b,-1,Spi1b)


SP4(1) = EpsS2w.dot.k1b
SP4(2) = EpsS2w.dot.p1t
SP4(3) = EpsS3g.dot.EpsS2w
SP4(4) = EpsS3g.dot.k1b
SP4(5) = EpsS3g.dot.k2w
SP4(6) = k1b.dot.k2w
SP4(7) = k1b.dot.p1t

end subroutine SPandPOLTop


subroutine CalcSMETop
implicit none
complex(8) :: help(1:4)
!complex(8) :: weyl_spin1(1:4),weyl_spin2(1:4),weyl_spin3(1:4)
integer i

do i=1,4
   SME(1,i) = Spi1b(i)
enddo

help =  weyl_spin1(Spi1b,EpsS2w)
do i=1,4
   SME(2,i) = help(i)
enddo

help = weyl_spin1(Spi1b,k2w) - m_top*Spi1b(1:4)
do i=1,4
   SME(3,i) = help(i)
enddo

help = weyl_spin1(Spi1b,EpsS3g)
do i=1,4
   SME(4,i) = help(i)
enddo

help = weyl_spin2(Spi1b,k2w,EpsS2w)
do i=1,4
   SME(5,i) = help(i)
enddo

help = weyl_spin2(Spi1b,EpsS3g,EpsS2w)
do i=1,4
   SME(6,i) = help(i)
enddo

help = weyl_spin2(Spi1b,k2w,EpsS3g)
do i=1,4
   SME(7,i) = help(i)
enddo

help = weyl_spin3(Spi1b,k2w,EpsS3g,EpsS2w) - m_top*SME(6,1:4)
do i=1,4
   SME(8,i) = help(i)
enddo

end subroutine CalcSMETop

subroutine TopIntegralsTTBJ(xe,muRF)
implicit none
real(8) ::Mw, mt,mzero,muRF, p1q,p2q,p3q, p4q, m1q,m2q,m3q,m4q, p12q,p23q, muq
complex(8) :: qlI1,qlI2, qlI3, qlI4
integer :: xe
muq =muRF**2
mzero = 0.d0
mt = m_top
mw = m_w
m1q = mt**2
       Integrals(1) = qlI1(m1q,muq,xe)
p1q = dble(Mw**2 + 2*SP4(6))
       m1q = mt**2
       m2q = mzero
       Integrals(2) = qlI2(p1q,m1q,m2q,muq,xe)

p1q = dble(Mw**2)
       m1q = mt**2
       m2q = mzero
       Integrals(3) = qlI2(p1q,m1q,m2q,muq,xe)

p1q = dble(-2*(SP4(6) - SP4(7)))
       m1q = mzero
       m2q = mzero
       Integrals(4) = qlI2(p1q,m1q,m2q,muq,xe)

p1q = dble(mt**2 - 2*SP4(7))
       m1q = mt**2
       m2q = mzero
       Integrals(5) = qlI2(p1q,m1q,m2q,muq,xe)

p1q = dble(mt**2)
       m1q = mt**2
       m2q = mzero
       Integrals(6) = qlI2(p1q,m1q,m2q,muq,xe)

p1q = dble(mzero)
       p2q = dble(mt**2 - 2*SP4(7))
              p3q = dble(mt**2)
       m1q = mzero
       m2q = mzero
       m3q = mt**2
       Integrals(7) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(mzero)
       p2q = dble(mzero)
              p3q = dble(-2*(SP4(6) - SP4(7)))
       m1q = mzero
       m2q = mzero
       m3q = mzero
       Integrals(8) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(Mw**2 + 2*SP4(6))
       p2q = dble(mzero)
              p3q = dble(mt**2)
       m1q = mt**2
       m2q = mzero
       m3q = mzero
       Integrals(9) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(Mw**2)
       p2q = dble(mzero)
              p3q = dble(Mw**2 + 2*SP4(6))
       m1q = mt**2
       m2q = mzero
       m3q = mzero
       Integrals(10) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(Mw**2)
       p2q = dble(-2*(SP4(6) - SP4(7)))
              p3q = dble(mt**2)
       m1q = mt**2
       m2q = mzero
       m3q = mzero
       Integrals(11) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(Mw**2)
       p2q = dble(mzero)
              p3q = dble(mt**2 - 2*SP4(7))
       m1q = mzero
       m2q = mt**2
       m3q = mt**2
       Integrals(12) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(mzero)
       p2q = dble(Mw**2 + 2*SP4(6))
              p3q = dble(mt**2)
       m1q = mt**2
       m2q = mt**2
       m3q = mzero
       Integrals(13) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(mzero)
       p2q = dble(Mw**2)
              p3q = dble(mt**2 - 2*SP4(7))
       m1q = mzero
       m2q = mzero
       m3q = mt**2
       Integrals(14) = qlI3(p1q,p2q,p3q,m1q,m2q,m3q,muq,xe)

p1q = dble(mzero)
       p2q = dble(Mw**2)
       p3q = dble(mzero)
       p4q = dble(mt**2)
       p12q = dble(Mw**2 + 2*SP4(6))
       p23q = dble(mt**2 - 2*SP4(7))

       m1q = mzero
       m2q = mzero
       m3q = mt**2
       m4q = mt**2

       Integrals(15) = qlI4(p1q,p2q,p3q,p4q,p12q,p23q,m1q,m2q,m3q,m4q,muq,xe)

p1q = dble(mzero)
       p2q = dble(mzero)
       p3q = dble(Mw**2)
       p4q = dble(mt**2)
       p12q = dble(-2*(SP4(6) - SP4(7)))
       p23q = dble(mt**2 - 2*SP4(7))

       m1q = mzero
       m2q = mzero
       m3q = mzero
       m4q = mt**2

       Integrals(16) = qlI4(p1q,p2q,p3q,p4q,p12q,p23q,m1q,m2q,m3q,m4q,muq,xe)

p1q = dble(Mw**2)
       p2q = dble(mzero)
       p3q = dble(mzero)
       p4q = dble(mt**2)
       p12q = dble(Mw**2 + 2*SP4(6))
       p23q = dble(-2*(SP4(6) - SP4(7)))

       m1q = mt**2
       m2q = mzero
       m3q = mzero
       m4q = mzero

       Integrals(17) = qlI4(p1q,p2q,p3q,p4q,p12q,p23q,m1q,m2q,m3q,m4q,muq,xe)


end subroutine TopIntegralsTTBJ

subroutine CalcCTsAtop(xe,mu, CT, CTmass)
implicit none
integer :: xe
real(8) :: dZMtop,dZgsTop,dZgluonTop,dZwave_top, mu, Ca, Cf
complex(8) :: amp_weyl_ctgsTop(1:4), amp_weyl_ctmass(1:4),amp_weyl_ctgluonTop(1:4), amp_weyl_ctwave(1:4), CT(1:4), CTmass(1:4)
dZMtop = 0.d0
dZgsTop = 0.d0
dZgluonTop = 0.d0
dZwave_top = 0.d0

Ca = Nc
Cf = (Nc**2-1.d0)/2.d0/Nc

CTmass(1:4) = zero
CT(1:4) = zero

  if(xe .eq. 0 .or. xe .eq. -1) then
     call CTmassAtop
     call CTwaveAtop
     call CTgluonATop
     call CTgsATop
     if(xe .eq. 0 ) then
        dZMtop = 1.d0/2.d0*(DenNc-Nc)*m_top*alpha_S*RunAlphaS(NLOParam,MuRen)/4.d0/Pi*(3.d0* dlog(mu**2/m_top**2) + 4.d0)
        dZgsTop = alpha_S*RunAlphaS(NLOParam,MuRen)/12.d0/Pi*dlog(mu**2/m_top**2)
        dZgluonTop = -alpha_S*RunAlphaS(NLOParam,MuRen)/6.d0/Pi*dlog(mu**2/m_top**2)
        dZwave_top = GS**2/(4.d0*Pi)**2* 1.d0/2.d0*(Nc-DenNc)*(-4.d0 -3.d0*dlog(mu**2/m_top**2))
     elseif(xe .eq. -1 ) then
        dZMtop = 1.d0/2.d0*(DenNc-Nc)*m_top*alpha_S*RunAlphaS(NLOParam,MuRen)/4.d0/Pi*3.d0
        dZgsTop    = alpha_S*RunAlphaS(NLOParam,MuRen)/4.d0/Pi*(N_f/3.d0-11.d0/2.d0) + alpha_S*RunAlphaS(NLOParam,MuRen)/12.d0/Pi
        dZgluonTop = -alpha_S*RunAlphaS(NLOParam,MuRen)/6.d0/Pi
        dZwave_top = GS**2/(4.d0*Pi)**2*3.d0/2.d0*(DenNc-Nc)
     endif

     amp_weyl_ctmass(1:4) = dZMtop* &
          ( CoeffsAtopCTmass(1)*SMEA(1,1:4) + CoeffsAtopCTmass(2)*SMEA(2,1:4) + CoeffsAtopCTmass(3)*SMEA(3,1:4) + CoeffsAtopCTmass(4)*SMEA(4,1:4) &
          + CoeffsAtopCTmass(5)*SMEA(5,1:4) + CoeffsAtopCTmass(6)*SMEA(6,1:4) + CoeffsAtopCTmass(7)*SMEA(7,1:4) + CoeffsAtopCTmass(8)*SMEA(8,1:4) )

     amp_weyl_ctgsTop(1:4) = dZgsTop* &
          ( CoeffsAtopCtgs(1)*SMEA(1,1:4) + CoeffsAtopCtgs(2)*SMEA(2,1:4) + CoeffsAtopCtgs(3)*SMEA(3,1:4) + CoeffsAtopCtgs(4)*SMEA(4,1:4) &
          + CoeffsAtopCtgs(5)*SMEA(5,1:4) + CoeffsAtopCtgs(6)*SMEA(6,1:4) + CoeffsAtopCtgs(7)*SMEA(7,1:4) + CoeffsAtopCtgs(8)*SMEA(8,1:4) )

     amp_weyl_ctgluonTop(1:4) = dZgluonTop* &
          ( CoeffsAtopCTgluon(1)*SMEA(1,1:4) + CoeffsAtopCTgluon(2)*SMEA(2,1:4) + CoeffsAtopCTgluon(3)*SMEA(3,1:4) + CoeffsAtopCTgluon(4)*SMEA(4,1:4) &
          + CoeffsAtopCTgluon(5)*SMEA(5,1:4) + CoeffsAtopCTgluon(6)*SMEA(6,1:4) + CoeffsAtopCTgluon(7)*SMEA(7,1:4) + CoeffsAtopCTgluon(8)*SMEA(8,1:4))

     amp_weyl_ctwave(1:4) = dZwave_top* &
          ( CoeffsAtopCTwave(1)*SMEA(1,1:4) + CoeffsAtopCTwave(2)*SMEA(2,1:4) + CoeffsAtopCTwave(3)*SMEA(3,1:4) + CoeffsAtopCTwave(4)*SMEA(4,1:4) &
          + CoeffsAtopCTwave(5)*SMEA(5,1:4) + CoeffsAtopCTwave(6)*SMEA(6,1:4) + CoeffsAtopCTwave(7)*SMEA(7,1:4) + CoeffsAtopCTwave(8)*SMEA(8,1:4) )
  endif

  CTmass(1:4) = amp_weyl_ctmass(1:4)
  CT(1:4) = amp_weyl_ctwave(1:4) + amp_weyl_ctgsTop(1:4)+amp_weyl_ctgluonTop(1:4)

end subroutine CalcCTsAtop


subroutine CalcCTsTop(xe,mu, CT, CTmass)
implicit none
integer :: xe
real(8) :: dZMtop,dZgstop,dZgluonTop,dZwave_top, mu, Ca, CF
complex(8) :: amp_weyl_ctgsTop(1:4), amp_weyl_ctmass(1:4),amp_weyl_ctgluonTop(1:4), amp_weyl_ctwave(1:4), CT(1:4), CTmass(1:4)
dZMtop = 0.d0
dZgsTop = 0.d0
dZgluonTop = 0.d0
dZwave_top = 0.d0

Ca = Nc
Cf = (Nc-DenNc)/2.d0

CTmass(1:4) = zero
CT(1:4) = zero

  if(xe .eq. 0 .or. xe .eq. -1) then
     call CTmassTop
     call CTwaveTop
     call CTgluonTop
     call CTgsTop
     if(xe .eq. 0 ) then
        dZMtop = 1.d0/2.d0*(DenNc-Nc)*m_top*alpha_S*RunAlphaS(NLOParam,MuRen)/4.d0/Pi*(3.d0* dlog(mu**2/m_top**2) + 4.d0)
        dZgsTop = alpha_S*RunAlphaS(NLOParam,MuRen)/12.d0/Pi*dlog(mu**2/m_top**2)
        dZgluonTop = -alpha_S*RunAlphaS(NLOParam,MuRen)/6.d0/Pi*dlog(mu**2/m_top**2)
        dZwave_top = GS**2/(4.d0*Pi)**2* 1.d0/2.d0*(Nc-DenNc)*(-4.d0 -3.d0*dlog(mu**2/m_top**2))
     elseif(xe .eq. -1 ) then
        dZMtop = 1.d0/2.d0*(DenNc-Nc)*m_top*alpha_S*RunAlphaS(NLOParam,MuRen)/4.d0/Pi*3.d0
        dZgsTop    = alpha_S*RunAlphaS(NLOParam,MuRen)/4.d0/Pi*(N_f/3.d0-11.d0/2.d0) + alpha_S*RunAlphaS(NLOParam,MuRen)/12.d0/Pi
        dZgluonTop = -alpha_S*RunAlphaS(NLOParam,MuRen)/6.d0/Pi
        dZwave_top = GS**2/(4.d0*Pi)**2*3.d0/2.d0*(DenNc-Nc)
     endif

     amp_weyl_ctmass(1:4) = dZMtop* &
          ( CoeffsTopCTmass(1)*SME(1,1:4) + CoeffsTopCTmass(2)*SME(2,1:4) + CoeffsTopCTmass(3)*SME(3,1:4) + CoeffsTopCTmass(4)*SME(4,1:4) &
          + CoeffsTopCTmass(5)*SME(5,1:4) + CoeffsTopCTmass(6)*SME(6,1:4) + CoeffsTopCTmass(7)*SME(7,1:4) + CoeffsTopCTmass(8)*SME(8,1:4) )

     amp_weyl_ctgsTop(1:4) = dZgsTop* &
          ( CoeffsTopCTgs(1)*SME(1,1:4) + CoeffsTopCTgs(2)*SME(2,1:4) + CoeffsTopCTgs(3)*SME(3,1:4) + CoeffsTopCTgs(4)*SME(4,1:4) &
          + CoeffsTopCTgs(5)*SME(5,1:4) + CoeffsTopCTgs(6)*SME(6,1:4) + CoeffsTopCTgs(7)*SME(7,1:4) + CoeffsTopCTgs(8)*SME(8,1:4) )

     amp_weyl_ctgluonTop(1:4) = dZgluonTop* &
          ( CoeffsTopCTgluon(1)*SME(1,1:4) + CoeffsTopCTgluon(2)*SME(2,1:4) + CoeffsTopCTgluon(3)*SME(3,1:4) + CoeffsTopCTgluon(4)*SME(4,1:4) &
          + CoeffsTopCTgluon(5)*SME(5,1:4) + CoeffsTopCTgluon(6)*SME(6,1:4) + CoeffsTopCTgluon(7)*SME(7,1:4) + CoeffsTopCTgluon(8)*SME(8,1:4) )

     amp_weyl_ctwave(1:4) = dZwave_top* &
          ( CoeffsTopCTwave(1)*SME(1,1:4) + CoeffsTopCTwave(2)*SME(2,1:4) + CoeffsTopCTwave(3)*SME(3,1:4) + CoeffsTopCTwave(4)*SME(4,1:4) &
          + CoeffsTopCTwave(5)*SME(5,1:4) + CoeffsTopCTwave(6)*SME(6,1:4) + CoeffsTopCTwave(7)*SME(7,1:4) + CoeffsTopCTwave(8)*SME(8,1:4) )
  end if

  CTmass(1:4) = amp_weyl_ctmass(1:4)
 CT(1:4) =  amp_weyl_ctwave(1:4) + amp_weyl_ctgsTop(1:4)+amp_weyl_ctgluonTop(1:4)

end subroutine CalcCTsTop


      subroutine LeadingOrderTop
      implicit none
      complex(8) :: t31,t32,t33,t34,t35,t36,t37,t38,t39,t30
      complex(8) :: t21,t22,t23,t24,t25,t26,t27,t28,t29,t20
      complex(8) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
      complex(8) :: t11,t12,t13,t14,t15,t16,t17,t18,t19
      real(8) :: mt, Mw
      mt = m_Top
      Mw = M_W
            t2 = EEL*GS*GW
      t3 = SP4(5)
      t5 = SP4(6)
      t6 = SP4(7)
      t8 = 1.d0/(t5-t6)
      t11 = SP4(3)
      t13 = mt**2
      t14 = Mw**2
      t17 = 1.d0/(t13-t14-2.d0*t5)
      t21 = SP4(2)
      t24 = SP4(1)
      t33 = t17*t8
      CoeffsTopLO(1) = czero
      CoeffsTopLO(2) = t2*ne*t3*t8
      CoeffsTopLO(3) = 2.d0*t2*ne*t11*t17
      CoeffsTopLO(4) = -t2*ne*(t13*t21-t14*t21+2.d0*t24*t5-2.d0*t21*t5-2.d0*t24*&
      t6)*t33
      CoeffsTopLO(5) = czero
      CoeffsTopLO(6) = czero
      CoeffsTopLO(7) = czero
      CoeffsTopLO(8) = -t2*ne*(t13-t14-2.d0*t6)*t33/2.d0

      end subroutine LeadingOrderTop

      subroutine CTmassTop
      implicit none
      complex(8) :: t31,t32,t33,t34,t35,t36,t37,t38,t39,t30
      complex(8) :: t21,t22,t23,t24,t25,t26,t27,t28,t29,t20
      complex(8) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
      complex(8) :: t11,t12,t13,t14,t15,t16,t17,t18,t19
      real(8) :: mt, Mw
      mt = m_Top
      Mw = M_W
      t1 = EEL*GS
      t2 = t1*GW
      t3 = SP4(3)
      t5 = mt**2
      t6 = Mw**2
      t7 = SP4(6)
      t9 = t5-t6-2.d0*t7
      t10 = 1.d0/t9
      t14 = mt*ne
      t15 = t9**2
      t16 = 1.d0/t15
      t21 = SP4(1)
      CoeffsTopCTmass(1) = -2.d0*t2*ne*t3*t10
      CoeffsTopCTmass(2) = czero
      CoeffsTopCTmass(3) = -4.d0*t2*t14*t3*t16
      CoeffsTopCTmass(4) = 4.d0*t2*t14*t21*t16
      CoeffsTopCTmass(5) = czero
      CoeffsTopCTmass(6) = t1*GW*ne*t10
      CoeffsTopCTmass(7) = czero
      CoeffsTopCTmass(8) = 2.d0*t2*t14*t16

      end subroutine CTmassTop

      subroutine CTwaveTop
      implicit none
      complex(8) :: t31,t32,t33,t34,t35,t36,t37,t38,t39,t30
      complex(8) :: t21,t22,t23,t24,t25,t26,t27,t28,t29,t20
      complex(8) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
      complex(8) :: t11,t12,t13,t14,t15,t16,t17,t18,t19
      real(8) :: mt, Mw
      mt = m_Top
      Mw = M_W
      t2 = EEL*GS*GW
      t3 = SP4(5)
      t5 = SP4(6)
      t6 = SP4(7)
      t7 = t5-t6
      t11 = SP4(3)
      t13 = mt**2
      t14 = Mw**2
      t17 = 1.d0/(t13-t14-2.d0*t5)
      t20 = SP4(2)
      t23 = SP4(1)
      t33 = t17/t7
      CoeffsTopCTwave(1) = czero
      CoeffsTopCTwave(2) = t2*ne*t3/t7/2.d0
      CoeffsTopCTwave(3) = t2*ne*t11*t17
      CoeffsTopCTwave(4) = -t2*ne*(t13*t20-t14*t20+2.d0*t23*t5-2.d0*t20*t5-2.d0*&
      t23*t6)*t33/2.d0
      CoeffsTopCTwave(5) = czero
      CoeffsTopCTwave(6) = czero
      CoeffsTopCTwave(7) = czero
      CoeffsTopCTwave(8) = -t2*ne*(t13-t14-2.d0*t6)*t33/4.d0

      end subroutine CTwaveTop

      subroutine CTgsTop
      implicit none
      complex(8) :: t31,t32,t33,t34,t35,t36,t37,t38,t39,t30
      complex(8) :: t21,t22,t23,t24,t25,t26,t27,t28,t29,t20
      complex(8) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
      complex(8) :: t11,t12,t13,t14,t15,t16,t17,t18,t19
      real(8) :: mt, Mw
      mt = m_Top
      Mw = M_W
      t2 = EEL*GS*GW
      t3 = SP4(5)
      t5 = SP4(6)
      t6 = SP4(7)
      t8 = 1.d0/(t5-t6)
      t11 = SP4(3)
      t13 = mt**2
      t14 = Mw**2
      t17 = 1.d0/(t13-t14-2.d0*t5)
      t21 = SP4(2)
      t24 = SP4(1)
      t33 = t17*t8
      CoeffsTopCTgs(1) = czero
      CoeffsTopCTgs(2) = t2*ne*t3*t8
      CoeffsTopCTgs(3) = 2.d0*t2*ne*t11*t17
      CoeffsTopCTgs(4) = -t2*ne*(t13*t21-t14*t21+2.d0*t24*t5-2.d0*t21*t5-2.d0*t24&
      *t6)*t33
      CoeffsTopCTgs(5) = czero
      CoeffsTopCTgs(6) = czero
      CoeffsTopCTgs(7) = czero
      CoeffsTopCTgs(8) = -t2*ne*(t13-t14-2.d0*t6)*t33/2.d0

      end subroutine CTgsTop

      subroutine CTgluonTop
      implicit none
      complex(8) :: t31,t32,t33,t34,t35,t36,t37,t38,t39,t30
      complex(8) :: t21,t22,t23,t24,t25,t26,t27,t28,t29,t20
      complex(8) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
      complex(8) :: t11,t12,t13,t14,t15,t16,t17,t18,t19
      real(8) :: mt, Mw
      mt = m_Top
      Mw = M_W
      t2 = EEL*GS*GW
      t3 = SP4(5)
      t5 = SP4(6)
      t6 = SP4(7)
      t7 = t5-t6
      t11 = SP4(3)
      t13 = mt**2
      t14 = Mw**2
      t17 = 1.d0/(t13-t14-2.d0*t5)
      t20 = SP4(2)
      t23 = SP4(1)
      t33 = t17/t7
      CoeffsTopCTgluon(1) = czero
      CoeffsTopCTgluon(2) = t2*ne*t3/t7/2.d0
      CoeffsTopCTgluon(3) = t2*ne*t11*t17
      CoeffsTopCTgluon(4) = -t2*ne*(t13*t20-t14*t20+2.d0*t23*t5-2.d0*t20*t5-2.d0*&
      t23*t6)*t33/2.d0
      CoeffsTopCTgluon(5) = czero
      CoeffsTopCTgluon(6) = czero
      CoeffsTopCTgluon(7) = czero
      CoeffsTopCTgluon(8) = -t2*ne*(t13-t14-2.d0*t6)*t33/4.d0

      end subroutine CTgluonTop


      subroutine LeadingOrderAtop
      implicit none
      complex(8) :: t31,t32,t33,t34,t35,t36,t37,t38,t39,t30,t46
      complex(8) :: t21,t22,t23,t24,t25,t26,t27,t28,t29,t20
      complex(8) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
      complex(8) :: t11,t12,t13,t14,t15,t16,t17,t18,t19
      real(8) :: mt, Mw
      mt = m_Top
      Mw = M_W
      t1 = EEL*GS
      t2 = t1*GW
      t4 = SP4(3)
      t5 = SP4(6)
      t6 = SP4(7)
      t8 = 1.d0/(t5-t6)
      t13 = SP4(5)
      t15 = mt**2
      t16 = Mw**2
      t19 = 1.d0/(t15-t16-2.d0*t5)
      t26 = SP4(2)
      t29 = SP4(1)
      t38 = t19*t8
      t46 = ne*(t15-t16-2.d0*t6)*t38
      CoeffsAtopLO(1) = -2.d0*t2*mt*ne*t4*t8
      CoeffsAtopLO(2) = -2.d0*t2*ne*t13*t19
      CoeffsAtopLO(3) = -t2*ne*t4*t8
      CoeffsAtopLO(4) = -t2*ne*(t15*t26-t16*t26+2.d0*t29*t5-2.d0*t26*t5-2.d0*t29*&
      t6)*t38
      CoeffsAtopLO(5) = czero
      CoeffsAtopLO(6) = t1*GW*mt*t46
      CoeffsAtopLO(7) = czero
      CoeffsAtopLO(8) = t2*t46/2.d0

      end subroutine LeadingOrderAtop

      subroutine CTmassAtop
      implicit none
      complex(8) :: t31,t32,t33,t34,t35,t36,t37,t38,t39,t30,t46
      complex(8) :: t21,t22,t23,t24,t25,t26,t27,t28,t29,t20
      complex(8) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
      complex(8) :: t11,t12,t13,t14,t15,t16,t17,t18,t19
      real(8) :: mt, Mw
      mt = m_Top
      Mw = M_W
      t2 = EEL*GS*GW
      t3 = mt*ne
      t4 = SP4(5)
      t5 = mt**2
      t6 = Mw**2
      t7 = SP4(6)
      t8 = 2.d0*t7
      t10 = (t5-t6-t8)**2
      t11 = 1.d0/t10
      t16 = SP4(1)
      CoeffsAtopCTmass(1) = czero
      CoeffsAtopCTmass(2) = 4.d0*t2*t3*t4*t11
      CoeffsAtopCTmass(3) = czero
      CoeffsAtopCTmass(4) = 4.d0*t2*t3*t16*t11
      CoeffsAtopCTmass(5) = czero
      CoeffsAtopCTmass(6) = -t2*ne*(3.d0*t5+t6+t8)*t11
      CoeffsAtopCTmass(7) = czero
      CoeffsAtopCTmass(8) = -2.d0*t2*t3*t11

      end subroutine CTmassAtop

      subroutine CTwaveAtop
      implicit none
      complex(8) :: t31,t32,t33,t34,t35,t36,t37,t38,t39,t30,t46
      complex(8) :: t21,t22,t23,t24,t25,t26,t27,t28,t29,t20
      complex(8) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
      complex(8) :: t11,t12,t13,t14,t15,t16,t17,t18,t19
      real(8) :: mt, Mw
      mt = m_Top
      Mw = M_W
      t1 = EEL*GS
      t2 = t1*GW
      t4 = SP4(3)
      t5 = SP4(6)
      t6 = SP4(7)
      t7 = t5-t6
      t8 = 1.d0/t7
      t12 = SP4(5)
      t14 = mt**2
      t15 = Mw**2
      t18 = 1.d0/(t14-t15-2.d0*t5)
      t25 = SP4(2)
      t28 = SP4(1)
      t37 = t18*t8
      t46 = ne*(t14-t15-2.d0*t6)*t37
      CoeffsAtopCTwave(1) = -t2*mt*ne*t4*t8
      CoeffsAtopCTwave(2) = -t2*ne*t12*t18
      CoeffsAtopCTwave(3) = -t2*ne*t4/t7/2.d0
      CoeffsAtopCTwave(4) = -t2*ne*(t14*t25-t15*t25+2.d0*t28*t5-2.d0*t25*t5-2.d0*&
      t28*t6)*t37/2.d0
      CoeffsAtopCTwave(5) = czero
      CoeffsAtopCTwave(6) = t1*GW*mt*t46/2.d0
      CoeffsAtopCTwave(7) = czero
      CoeffsAtopCTwave(8) = t2*t46/4.d0

      end subroutine CTwaveAtop

      subroutine CTgsAtop
      implicit none
      complex(8) :: t31,t32,t33,t34,t35,t36,t37,t38,t39,t30,t46
      complex(8) :: t21,t22,t23,t24,t25,t26,t27,t28,t29,t20
      complex(8) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
      complex(8) :: t11,t12,t13,t14,t15,t16,t17,t18,t19
      real(8) :: mt, Mw
      mt = m_Top
      Mw = M_W
      t1 = EEL*GS
      t2 = t1*GW
      t4 = SP4(3)
      t5 = SP4(6)
      t6 = SP4(7)
      t8 = 1.d0/(t5-t6)
      t13 = SP4(5)
      t15 = mt**2
      t16 = Mw**2
      t19 = 1.d0/(t15-t16-2.d0*t5)
      t26 = SP4(2)
      t29 = SP4(1)
      t38 = t19*t8
      t46 = ne*(t15-t16-2.d0*t6)*t38
      CoeffsAtopCTgs(1) = -2.d0*t2*mt*ne*t4*t8
      CoeffsAtopCTgs(2) = -2.d0*t2*ne*t13*t19
      CoeffsAtopCTgs(3) = -t2*ne*t4*t8
      CoeffsAtopCTgs(4) = -t2*ne*(t15*t26-t16*t26+2.d0*t29*t5-2.d0*t26*t5-2.d0*&
      t29*t6)*t38
      CoeffsAtopCTgs(5) = czero
      CoeffsAtopCTgs(6) = t1*GW*mt*t46
      CoeffsAtopCTgs(7) = czero
      CoeffsAtopCTgs(8) = t2*t46/2.d0

      end subroutine CTgsAtop

      subroutine CTgluonAtop
      implicit none
      complex(8) :: t31,t32,t33,t34,t35,t36,t37,t38,t39,t30,t46
      complex(8) :: t21,t22,t23,t24,t25,t26,t27,t28,t29,t20
      complex(8) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
      complex(8) :: t11,t12,t13,t14,t15,t16,t17,t18,t19
      real(8) :: mt, Mw
      mt = m_Top
      Mw = M_W
      t1 = EEL*GS
      t2 = t1*GW
      t4 = SP4(3)
      t5 = SP4(6)
      t6 = SP4(7)
      t7 = t5-t6
      t8 = 1.d0/t7
      t12 = SP4(5)
      t14 = mt**2
      t15 = Mw**2
      t18 = 1.d0/(t14-t15-2.d0*t5)
      t25 = SP4(2)
      t28 = SP4(1)
      t37 = t18*t8
      t46 = ne*(t14-t15-2.d0*t6)*t37
      CoeffsAtopCTgluon(1) = -t2*mt*ne*t4*t8
      CoeffsAtopCTgluon(2) = -t2*ne*t12*t18
      CoeffsAtopCTgluon(3) = -t2*ne*t4/t7/2.d0
      CoeffsAtopCTgluon(4) = -t2*ne*(t14*t25-t15*t25+2.d0*t28*t5-2.d0*t25*t5-2.d0&
      *t28*t6)*t37/2.d0
      CoeffsAtopCTgluon(5) = czero
      CoeffsAtopCTgluon(6) = t1*GW*mt*t46/2.d0
      CoeffsAtopCTgluon(7) = czero
      CoeffsAtopCTgluon(8) = t2*t46/4.d0

      end subroutine CTgluonAtop

!!! NLO





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!STANDARD MATRIX ELEMENTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!STANDARD MATRIX ELEMENTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!STANDARD MATRIX ELEMENTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Dirac_sme0(sp1,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4)
complex(8) :: Dirac_sme0
Dirac_sme0 = psp1_(sp1,sp2)
return
end function Dirac_sme0


function Dirac_sme1(sp1,v1,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4)
complex(8) :: Dirac_sme1
sp9(1:4) =  spb2_(sp1,v1)
Dirac_sme1 = psp1_(sp9,sp2)
return
end function Dirac_sme1


function Dirac_sme2(sp1,v1,v2,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4)
complex(8) :: Dirac_sme2
sp9(1:4) =  spb2_(sp1,v1)
sp7(1:4) =  spb2_(sp9,v2)
Dirac_sme2 = psp1_(sp7,sp2)
return
end function Dirac_sme2


function Dirac_sme3(sp1,v1,v2,v3,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4), v3(1:4)
complex(8) :: Dirac_sme3
sp9(1:4) =  spb2_(sp1,v1)
sp7(1:4) =  spb2_(sp9,v2)
sp9(1:4) =  spb2_(sp7,v3)
Dirac_sme3 = psp1_(sp9,sp2)
return
end function Dirac_sme3



function Dirac_sme4(sp1,v1,v2,v3,v4,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4), v3(1:4), v4(1:4)
complex(8) :: Dirac_sme4
sp9(1:4) =  spb2_(sp1,v1)
sp7(1:4) =  spb2_(sp9,v2)
sp9(1:4) =  spb2_(sp7,v3)
sp7(1:4) =  spb2_(sp9,v4)
Dirac_sme4 = psp1_(sp7,sp2)
return
end function Dirac_sme4


function Dirac_sme5(sp1,v1,v2,v3,v4,v5,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4), v3(1:4), v4(1:4), v5(1:4)
complex(8) :: Dirac_sme5
sp9(1:4) =  spb2_(sp1,v1)
sp7(1:4) =  spb2_(sp9,v2)
sp9(1:4) =  spb2_(sp7,v3)
sp7(1:4) =  spb2_(sp9,v4)
sp9(1:4) =  spb2_(sp7,v5)
Dirac_sme5 = psp1_(sp9,sp2)
return
end function Dirac_sme5


function Dirac_sme6(sp1,v1,v2,v3,v4,v5,v6,sp2)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp9(1:4), sp7(1:4)
complex(8) :: v1(1:4),v2(1:4), v3(1:4), v4(1:4), v5(1:4), v6(1:4)
complex(8) :: Dirac_sme6
sp9(1:4) =  spb2_(sp1,v1)
sp7(1:4) =  spb2_(sp9,v2)
sp9(1:4) =  spb2_(sp7,v3)
sp7(1:4) =  spb2_(sp9,v4)
sp9(1:4) =  spb2_(sp7,v5)
sp7(1:4) =  spb2_(sp9,v6)
Dirac_sme6 = psp1_(sp7,sp2)
return
end function Dirac_sme6



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!MULTIPLY SPINOR-CHAINS
function Spi0Spi0(sp1,sp2,sp3,sp4)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp3(1:4), sp4(1:4)
complex(8) :: vsp1(1:4),vsp2(1:4)
complex(8) :: Spi0Spi0, res
integer :: Dv
Dv = 4
vsp1(1:4) = vbqq(Dv,sp1,sp2)
vsp2(1:4) = vbqq(Dv,sp3,sp4)
res = vsp1.dot.vsp2

Spi0Spi0 =res
return
end function Spi0Spi0


function Spi1Spi0(sp1,v1,sp2,sp3,sp4)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp3(1:4), sp4(1:4), sp7(1:4), sp9(1:4)
complex(8) :: vsp1(1:4),vsp2(1:4)
complex(8) :: v1(1:4)
complex(8) :: Spi1Spi0, res
integer :: Dv
Dv = 4

sp9(1:4) =  spb2_(sp1,v1)
vsp1(1:4) = vbqq(Dv,sp9,sp2)
vsp2(1:4) = vbqq(Dv,sp3,sp4)
res = vsp1.dot.vsp2

Spi1Spi0 =res
return
end function Spi1Spi0


function Spi0Spi1(sp1,sp2,sp3,v1,sp4)
implicit none
complex(8) :: sp1(1:4),sp2(1:4), sp3(1:4), sp4(1:4), sp7(1:4), sp9(1:4)
complex(8) :: vsp1(1:4),vsp2(1:4)
complex(8) :: v1(1:4)
complex(8) :: Spi0Spi1, res
integer :: Dv
Dv = 4

vsp1(1:4) = vbqq(Dv,sp1,sp2)
sp9(1:4) = spb2_(sp3,v1)
vsp2(1:4) = vbqq(Dv,sp9,sp4)
res = vsp1.dot.vsp2

Spi0Spi1 =res
return
end function Spi0Spi1

end Module
