MODULE ModAmplitudes_Zprime
  use ModMisc
  implicit none

contains

  subroutine Tree_Zprime_tbtqbq(LO_Res_Pol_Left, LO_Res_Pol_Right)
    use ModProcess
    use ModParameters
    complex(8), intent(out) :: LO_Res_Pol_Left, LO_Res_Pol_Right

    complex(8) :: ub1(4), v2(4), ub3(4), v4(4)
    complex(8) :: p1(4), p2(4), p3(4), p4(4)
    complex(8) :: lqcurrL(4), lqcurrR(4) ! Light quark currents
    complex(8) :: hqcurr(4) ! Heavy quark current
    complex(8) :: propfact
    
    
    ub1(1:4) = ExtParticle(4)%Pol(1:4)
    v2(1:4) = ExtParticle(3)%Pol(1:4)
    ub3(1:4) = ExtParticle(2)%Pol(1:4)
    v4(1:4) = ExtParticle(1)%Pol(1:4)
    
    p1(1:4) = ExtParticle(4)%Mom(1:4)
    p2(1:4) = ExtParticle(3)%Mom(1:4)
    p3(1:4) = ExtParticle(2)%Mom(1:4)
    p4(1:4) = ExtParticle(1)%Mom(1:4)

    lqcurrL = (-cI) * vbqq2(ub1,chir(.false.,v2))
    lqcurrR = (-cI) * vbqq2(ub1,chir(.true.,v2))

    hqcurr = (-cI) * (vbqq2(ub3,chir(.false.,v4)) * gL_Zpr(top_) + vbqq2(ub3,chir(.true.,v4)) * gR_Zpr(top_))

    propfact = -cI / (sc_(p1+p2,p1+p2)-m_Zpr**2 + cI * Ga_Zpr * M_Zpr)
    hqcurr = hqcurr * propfact
    
    LO_Res_Pol_Left  = sc_(lqcurrL,hqcurr) 
    LO_Res_Pol_Right = sc_(lqcurrR,hqcurr) 

    return

  end subroutine Tree_Zprime_tbtqbq


  subroutine Virt_Zprime_tbtqbq(Virt_Res_Pol_Left, Virt_Res_Pol_Right) ! Vertex correction and Z factor
    use ModParameters
    use ModProcess
    use ModMisc
    implicit none
    complex(8) :: p1(4), p2(4) ,p3(4), p4(4), q(4), t(4)
    real(8) :: q2, mt2, z, Cf
    complex(8) :: ub1(4), v2(4), ub3(4), v4(4)
    integer :: epv
    complex(8) :: lqcurrL(4), lqcurrR(4), hqcurr(4), scalcurr(4), axcurr(4)
    complex(8) :: ampL, ampR, ampVeL, ampVeR, ampAxL, ampAxR, propfact
    complex(8) :: massless_v, massive_v, massive_v_II, z2
    complex(8) :: qlI2, qlI3
    complex(8) :: Virt_Res_Pol_Left, Virt_Res_Pol_Right

    !!! For singularity checks: 
    complex(8) :: catani, LO_L, LO_R
    real(8) :: v34

    
    ! epsilon order
    epv = 0
    
    ! Cf
    Cf = 4d0/3d0
    
    ub1(1:4) = ExtParticle(4)%Pol(1:4)
    v2(1:4) = ExtParticle(3)%Pol(1:4)
    ub3(1:4) = ExtParticle(2)%Pol(1:4)
    v4(1:4) = ExtParticle(1)%Pol(1:4)
    
    p1(1:4) = ExtParticle(4)%Mom(1:4)
    p2(1:4) = ExtParticle(3)%Mom(1:4)
    p3(1:4) = ExtParticle(2)%Mom(1:4)
    p4(1:4) = ExtParticle(1)%Mom(1:4)

    q = p3 + p4
    t = p4 - p3
    q2 = dreal(sc_(q,q))
    mt2 = m_Top**2
    z = q2/mt2

    ! Building all currents
    propfact = -cI / (2d0 * sc_(p1,p2)-m_Zpr**2 + cI * Ga_Zpr * M_Zpr)
    
    lqcurrL = -cI * vbqq2(ub1,chir(.false.,v2))
    lqcurrR = -cI * vbqq2(ub1,chir(.true.,v2))
    
    hqcurr =  -cI * vbqq2(ub3,chir(.false.,v4)) * gL_Zpr(top_) -cI * vbqq2(ub3,chir(.true.,v4)) * gR_Zpr(top_)
    scalcurr = -cI * psp1_(ub3,v4) * t/(2d0*m_Top) * (gR_Zpr(top_)+gL_Zpr(top_))/2d0
    axcurr = -cI * (vbqq2(ub3,chir(.true.,v4))-vbqq2(ub3,chir(.false.,v4))) * (gR_Zpr(top_)-gL_Zpr(top_))/2d0
    
    
    ! Z2 correction, CDR
    z2 = -3d0 * (-epv) * (2d0+epv) - (4d0 + 3d0 * dlog(MuRen**2/mt2))*(1d0+epv)*(1d0+epv/2d0)

    ! Massless vertex correction
    ! Omit common prefactor (4pi)^ep/Gamma(1-ep)
    
    massless_v = - 2d0 * q2 * qlI3(q2,0d0,0d0,0d0,0d0,0d0,MuRen**2,epv)
    massless_v = massless_v - 3d0 * qlI2(q2,0d0,0d0,MuRen**2,epv)
!      massless_v = massless_v - 2d0 * (1d0+epv)*(1d0+epv/2d0) ! This is the CDR result, but Markus dipoles are in FDH
    massless_v = massless_v - 1d0 * (1d0+epv)*(1d0+epv/2d0) ! This is the FDH result. DO NOT TOUCH Z2 AND MASSIVE VERTEX, AS THE SCHEME CHANGE CANCELS OUT
    
    ! Massive vertex correction, term proportional to tree level, CDR
    massive_v = - 2d0 * (1d0-2d0/z) * q2 * qlI3(q2,mt2,mt2,mt2,mt2,0d0,MuRen**2,epv)
    massive_v = massive_v - 3d0 * qlI2(q2,mt2,mt2,MuRen**2,epv)
    massive_v = massive_v + 4d0 * qlI2(mt2,mt2,0d0,MuRen**2,epv)
    massive_v = massive_v - 2d0 *(1d0+epv)*(1d0+epv/2d0)
    
    ! Massive vertex, form factor for the scalar current
    massive_v_II = qlI2(q2,mt2,mt2,MuRen**2,epv)-qlI2(mt2,mt2,0d0,MuRen**2,epv)
    
    ampL = sc_(hqcurr,lqcurrL)*propfact*(massless_v+massive_v+z2) 
    ampR = sc_(hqcurr,lqcurrR)*propfact*(massless_v+massive_v+z2) 
    ampVeL = -4d0/(4d0-z)*sc_(lqcurrL,scalcurr)*propfact*massive_v_II 
    ampVeR = -4d0/(4d0-z)*sc_(lqcurrR,scalcurr)*propfact*massive_v_II 
    ampAxL = -8d0/(z-4d0)*sc_(lqcurrL,axcurr)*propfact*massive_v_II
    ampAxR = -8d0/(z-4d0)*sc_(lqcurrR,axcurr)*propfact*massive_v_II
    
    Virt_Res_Pol_Left = (ampL + ampVeL + ampAxL) * (Cf/(16d0*Pi**2))
    Virt_Res_Pol_Right = (ampR + ampVeR + ampAxR) * (Cf/(16d0*Pi**2))


    !!! Singular Part: check against Catani

    if (epv .ne. 0) then
       v34 = dsqrt(1d0 - 4d0 * m_Top**4 / (2d0 * dreal(sc_(p3,p4)))**2)
       catani = -1d0/v34 * dlog((1d0-v34)/(1d0+v34)) - 5d0 - 2d0 * dlog(MuRen**2/(2d0 * dreal(sc_(p1,p2))))
       catani = catani - 2 * cI * Pi - 2 * cI *  Pi/v34
       catani = catani * Cf / (4d0*Pi)**2

       call Tree_Zprime_tbtqbq(LO_L,LO_R)
       !    print *, 'lo CF / (4Pi)^2  ', (LO_L + LO_R)*Cf/16d0/Pi**2
       if ( epv.eq.0) then
          print *, 'res  ', Virt_Res_Pol_Left+Virt_Res_Pol_Right
       elseif ( epv.eq.-1) then
          print *, 'res  ', Virt_Res_Pol_Left+Virt_Res_Pol_Right
          print *, '1/ep ', (LO_L + LO_R)* catani
          print *, ''
       elseif (epv.eq.-2) then
          print *, 'tree ', LO_L+LO_R
          print *, 'res  ', Virt_Res_Pol_Left+Virt_Res_Pol_Right
          print *, '1/ep2', -2d0*Cf/(4d0*Pi)**2 * (LO_L+LO_R) 
          print *, ''
       endif
    endif

    !!! End of the singularity checks

    RETURN
    
  END subroutine Virt_Zprime_tbtqbq
  
  
  subroutine Tree_Zprime_tbtqbqg_i(prim1_L, prim1_R) ! Emission from the initial state
    use ModProcess
    use ModParameters
    use ModMisc
    complex(8), intent(out) :: prim1_L, prim1_R
    
    complex(8) :: ub1(4), v2(4), ub3(4), v4(4), e5(4), sptmp(4)
    complex(8) :: p1(4), p2(4), p3(4), p4(4), p5(4), p15(4), p25(4)
    complex(8) :: lqcurrL(4), lqcurrR(4) ! Light quark currents
    complex(8) :: hqcurr(4) ! Heavy quark current
    complex(8) :: propfact

    
    ub1(1:4) = ExtParticle(4)%Pol(1:4)
    v2(1:4) = ExtParticle(3)%Pol(1:4)
    ub3(1:4) = ExtParticle(2)%Pol(1:4)
    v4(1:4) = ExtParticle(1)%Pol(1:4)
    e5(1:4) = ExtParticle(5)%Pol(1:4)
    
    p1(1:4) = ExtParticle(4)%Mom(1:4)
    p2(1:4) = ExtParticle(3)%Mom(1:4)
    p3(1:4) = ExtParticle(2)%Mom(1:4)
    p4(1:4) = ExtParticle(1)%Mom(1:4)
    p5(1:4) = ExtParticle(5)%Mom(1:4)

    p15 = p1 + p5
    p25 = p2 + p5

    propfact = -cI / (sc_(p3+p4,p3+p4) - m_Zpr**2 + cI * m_Zpr * Ga_Zpr)
    hqcurr =  -cI * vbqq2(ub3,chir(.false.,v4)) * gL_Zpr(top_) -cI * vbqq2(ub3,chir(.true.,v4)) * gR_Zpr(top_)
    hqcurr = hqcurr * propfact

    ! First diagram: emission from quark line
    ! Phase factor: -cI (Vertex) * cI (Propagator) = 1
    sptmp = spb2_(ub1,e5)
    sptmp = spb2_(sptmp,p15)/sc_(p15,p15)
    lqcurrL = -cI * vbqq2(sptmp,chir(.false.,v2))
    lqcurrR = -cI * vbqq2(sptmp,chir(.true.,v2))


    ! Second diagram: emission from anti-quark line
    sptmp = spi2_(e5,v2)
    sptmp = -spi2_(p25,sptmp)/sc_(p25,p25)
    lqcurrL = lqcurrL - cI * vbqq2(ub1, chir(.false.,sptmp))
    lqcurrR = lqcurrR - cI * vbqq2(ub1, chir(.true.,sptmp))

    prim1_L = sc_(lqcurrL,hqcurr)/sqrt2 ! To match the definition of T
    prim1_R = sc_(lqcurrR,hqcurr)/sqrt2

    RETURN

  END subroutine Tree_Zprime_tbtqbqg_i

  subroutine Tree_Zprime_tbtqbqg_f(prim2_L, prim2_R) ! Emission from the final state
    use ModProcess
    use ModParameters
    use ModMisc
    complex(8), intent(out) :: prim2_L, prim2_R
    
    complex(8) :: ub1(4), v2(4), ub3(4), v4(4), e5(4), sptmp(4)
    complex(8) :: p1(4), p2(4), p3(4), p4(4), p5(4), p35(4), p45(4)
    complex(8) :: lqcurrL(4), lqcurrR(4) ! Light quark currents
    complex(8) :: hqcurr(4) ! Heavy quark current
    complex(8) :: propfact
    
    ub1(1:4) = ExtParticle(4)%Pol(1:4)
    v2(1:4) = ExtParticle(3)%Pol(1:4)
    ub3(1:4) = ExtParticle(2)%Pol(1:4)
    v4(1:4) = ExtParticle(1)%Pol(1:4)
    e5(1:4) = ExtParticle(5)%Pol(1:4)
    
    p1(1:4) = ExtParticle(4)%Mom(1:4)
    p2(1:4) = ExtParticle(3)%Mom(1:4)
    p3(1:4) = ExtParticle(2)%Mom(1:4)
    p4(1:4) = ExtParticle(1)%Mom(1:4)
    p5(1:4) = ExtParticle(5)%Mom(1:4)

    p35 = p3 + p5
    p45 = p4 + p5

    propfact = -cI / (2d0 * sc_(p1,p2) - m_Zpr**2 + cI * m_Zpr * Ga_Zpr)
    lqcurrL = -cI * vbqq2(ub1, chir(.false.,v2)) * propfact
    lqcurrR = -cI * vbqq2(ub1, chir(.true.,v2)) * propfact

    ! First diagram: emission from quark line
    ! Phase factor: -cI (Vertex) * cI (Propagator) = 1
    sptmp = spb2_(ub3,e5)
    sptmp = (spb2_(sptmp,p35)+m_Top * sptmp)/(sc_(p35,p35)-m_Top**2)
    hqcurr = -cI * vbqq2(sptmp,chir(.false.,v4)) * gL_Zpr(top_) -cI * vbqq2(sptmp,chir(.true.,v4)) * gR_Zpr(top_)

    ! Second diagram: emission from anti-quark line
    sptmp = spi2_(e5,v4)
    sptmp = (-spi2_(p45,sptmp)+m_Top * sptmp)/(sc_(p45,p45)-m_top**2)
    hqcurr = hqcurr - cI * vbqq2(ub3, chir(.false.,sptmp)) * gL_Zpr(top_) - cI * vbqq2(ub3, chir(.true.,sptmp)) * gR_Zpr(top_)

    prim2_L = sc_(lqcurrL,hqcurr)/sqrt2 ! To match the definition of T
    prim2_R = sc_(lqcurrR,hqcurr)/sqrt2 

    RETURN

  END subroutine Tree_Zprime_tbtqbqg_f


  subroutine Virt_Zprime_box(boxL, boxR)
    use ModProcess
    use ModParameters
    use ModMisc
    implicit none

    complex(8) :: p1(4), p2(4), p3(4), p4(4)
    complex(8) :: ub1(4), v2(4), ub3(4), v4(4)

    real(8) :: mt, mZp, mt2, mZp2, mur2
    real(8) :: s12, s13, s23, s14, s24, s, t, u

    complex(8) :: qlI2,qlI3,qlI4
    complex(8) :: bo1, bo2, tri1, tri2, tri3, tri4, tri5, tri6
    complex(8) :: bub1,bub2,bub3,bub4,bub5,bub6
    complex(8) :: grm3a, grm3b
    integer :: epv

    complex(8) :: spst1(-1:1), spst2(-1:1), spst3(-1:1), spst4(-1:1)
    complex(8) :: spst5(-1:1), spst6(-1:1), spst7(-1:1), spst8(-1:1)
    complex(8) :: spst9(-1:1), spst10(-1:1), spst11(-1:1)

    complex(8) :: resLL(-1:1), resLR(-1:1), boxL, boxR
    complex(8) :: LO_L, LO_R
    real(8) :: catani

    epv = -1

    mt = m_Top
    mZp = m_Zpr
    mt2 = m_Top**2
    mZp2 = m_Zpr**2
    mur2 = MuRen**2

    ub1(1:4) = ExtParticle(4)%Pol(1:4)
    v2(1:4) = ExtParticle(3)%Pol(1:4)
    ub3(1:4) = ExtParticle(2)%Pol(1:4)
    v4(1:4) = ExtParticle(1)%Pol(1:4)

    p1(1:4) = ExtParticle(4)%Mom(1:4)
    p2(1:4) = ExtParticle(3)%Mom(1:4)
    p3(1:4) = ExtParticle(2)%Mom(1:4)
    p4(1:4) = ExtParticle(1)%Mom(1:4)

    s12 = 2d0 * dreal(sc_(p1,p2))
    s13 = 2d0 * dreal(sc_(p1,p3))
    s23 = 2d0 * dreal(sc_(p2,p3))
    s14 = 2d0 * dreal(sc_(p1,p4))
    s24 = 2d0 * dreal(sc_(p2,p4))

    s = dreal(sc_(p1+p2,p1+p2))
    t = dreal(sc_(p2+p3,p2+p3))
    u = dreal(sc_(p2+p4,p2+p4))

    bo1 = qlI4(0d0,0d0,mt2,mt2,s,t,0d0,0d0,mZp2,mt2,mur2,epv)
    bo2 = qlI4(0d0,0d0,mt2,mt2,s,u,0d0,0d0,mZp2,mt2,mur2,epv)

    tri1 = qlI3(s,mt2,mt2,0d0,mZp2,mt2,mur2,epv)
    tri2 = qlI3(0d0,mt2,t,0d0,mZp2,mt2,mur2,epv)
    tri3 = qlI3(0d0,mt2,t,0d0,0d0,mt2,mur2,epv)
    tri4 = qlI3(0d0,mt2,u,0d0,0d0,mt2,mur2,epv)
    tri5 = qlI3(0d0,u,mt2,mZp2,0d0,mt2,mur2,epv)
    tri6 = qlI3(0d0,0d0,s,0d0,0d0,mZp2,mur2,epv)

    bub1 = qlI2(s,0d0,mZp2,mur2,epv)
    bub2 = qlI2(t,0d0,mt2,mur2,epv)
    bub3 = qlI2(mt2,mZp2,mt2,mur2,epv)
    bub4 = qlI2(mt2,0d0,mt2,mur2,epv)

    bub5 = qlI2(u,0d0,mt2,mur2,epv)
    bub6 = qlI2(0d0,mZp2,0d0,mur2,epv)

    spst1(-1) = spab(ub1,p2,v4)*spbb(ub3,v2)
    spst2(-1) = spab(ub1,p3,v2)*spab(ub3,p1,v4)
    spst3(-1) = spab(ub1,p3,v2)*spab(ub3,p2,v4)
    spst4(-1) = spab(ub1,p3,v2)*spbb(ub3,v4)
    spst5(-1) = spab(ub1,p3,v4)*spbb(ub3,v2)
    spst6(-1) = spab(ub3,p1,v2)*spaa(ub1,v4)
    spst7(-1) = spaa(ub1,v4)*spbb(ub3,v2)
    spst8(-1) = 1d0/2d0 * sc_(spab_4V(ub3,v4), spab_4V(ub1,v2))

    spst9(-1) = spab(ub1,p3,v2)*spba(ub3,p1,v4)
    spst10(-1) = spab(ub1,p3,v2)*spba(ub3,p2,v4)
    spst11(-1) = sc_(spab_4V(ub1,v2),spba_4V(ub3,v4))

    spst1(+1) = spba(ub1,p2,v4)*spaa(ub3,v2)
    spst2(+1) = spba(ub1,p3,v2)*spba(ub3,p1,v4)
    spst3(+1) = spba(ub1,p3,v2)*spba(ub3,p2,v4)
    spst4(+1) = spba(ub1,p3,v2)*spaa(ub3,v4)
    spst5(+1) = spba(ub1,p3,v4)*spaa(ub3,v2)
    spst6(+1) = spba(ub3,p1,v2)*spbb(ub1,v4)
    spst7(+1) = spbb(ub1,v4)*spaa(ub3,v2)
    spst8(+1) = 1d0/2d0 * sc_(spba_4V(ub3,v4), spba_4V(ub1,v2))

    spst9(+1) = spba(ub1,p3,v2)*spab(ub3,p1,v4)
    spst10(+1) = spba(ub1,p3,v2)*spab(ub3,p2,v4)
    spst11(+1) = sc_(spba_4V(ub1,v2),spab_4V(ub3,v4))

    grm3a = (s12*(-(mt**2*s12) + s13*s23))/4d0
    grm3b = (s12*(-(mt**2*s12) + s14*s24))/4d0
    grm3a = -1d0/grm3a ! To match Kirill's reduction
    grm3b = -1d0/grm3b ! To match Kirill's reduction

    resLL =  + spst8 * (  - 2d0*tri6 + tri3 + tri2 - 2d0*tri1 + 2d0*mZp**2&
     &    *bo1 + 2d0*s24*bo2 - 1d0/4d0*s23**2*mZp**2*grm3a*tri3 + 1d0/4d0*&
     &    s23**2*mZp**2*grm3a*tri2 - 1d0/4d0*s23**2*mZp**4*grm3a*bo1 + 1d0/&
     &    2d0*s12*mt**2*mZp**2*grm3a*tri2 - 1d0/2d0*s12*mt**2*mZp**2*grm3a&
     &    *tri1 + 1d0/4d0*s12*s24*mt**2*grm3b*tri5 + 1d0/4d0*s12*s24*mt**2*&
     &    grm3b*tri4 - 1d0/2d0*s12*s24*mt**2*grm3b*tri1 + 1d0/4d0*s12*s24*&
     &    mt**2*mZp**2*grm3b*bo2 - 1d0/4d0*s12*s23*mZp**2*grm3a*tri6 + 1d0/&
     &    2d0*s12*s23*mZp**2*grm3a*tri2 - 1d0/4d0*s12*s23*mZp**2*grm3a*&
     &    tri1 + 1d0/4d0*s12**2*mt**2*grm3b*tri6 - 1d0/4d0*s12**2*mt**2*&
     &    grm3b*tri1 - 1d0/4d0*s12**2*mt**2*grm3a*tri3 - 1d0/4d0*s12**2*&
     &    mt**2*grm3a*tri2 + 1d0/2d0*s12**2*mt**2*grm3a*tri1 - 1d0/2d0*&
     &    s12**2*mt**2*mZp**2*grm3a*bo1 - 1d0/4d0*s12**2*s24*mt**2*grm3b*&
     &    bo2 + 1d0/4d0*s12**2*s23*grm3a*tri6 - 1d0/4d0*s12**2*s23*grm3a*&
     &    tri3 - 1d0/4d0*s12**2*s23*grm3a*tri2 + 1d0/4d0*s12**2*s23*grm3a*&
     &    tri1 - 1d0/2d0*s12**2*s23*mZp**2*grm3a*bo1 - 1d0/4d0*s12**2*&
     &    s23**2*grm3a*bo1 )
      resLL = resLL + spst7 * (  - 1d0/8d0*s12*s23*mt**2*grm3a*tri3 - 1d0/8d0*s12*&
     &    s23*mt**2*grm3a*tri2 + 1d0/4d0*s12*s23*mt**2*grm3a*tri1 - 1d0/8d0&
     &    *s12*s23*mt**2*mZp**2*grm3a*bo1 - 1d0/8d0*s12**2*mt**2*grm3a*&
     &    tri6 + 1d0/8d0*s12**2*mt**2*grm3a*tri1 + 1d0/8d0*s12**2*s23*mt**2&
     &    *grm3a*bo1 )
      resLL = resLL + spst6 * ( mt*bo2 - mt*bo1 + 1d0/8d0*s24**2*mt*grm3b*tri5 - &
     &    1d0/8d0*s24**2*mt*grm3b*tri4 - 1d0/8d0*s24**2*mt*mZp**2*grm3b*bo2&
     &     - 1d0/8d0*s23**2*mt*grm3a*tri3 + 1d0/8d0*s23**2*mt*grm3a*tri2 - &
     &    1d0/8d0*s23**2*mt*mZp**2*grm3a*bo1 + 1d0/4d0*s12*mt**3*grm3b*tri5&
     &     - 1d0/4d0*s12*mt**3*grm3b*tri1 - 1d0/4d0*s12*mt**3*grm3a*tri3 + &
     &    1d0/4d0*s12*mt**3*grm3a*tri1 - 1d0/4d0*s12*mt**3*mZp**2*grm3a*bo1&
     &     - 1d0/8d0*s12*s24*mt*grm3b*tri6 + 1d0/4d0*s12*s24*mt*grm3b*tri5&
     &     - 1d0/8d0*s12*s24*mt*grm3b*tri1 - 1d0/8d0*s12*s24**2*mt*grm3b*&
     &    bo2 + 1d0/8d0*s12*s23*mt*grm3a*tri6 - 1d0/4d0*s12*s23*mt*grm3a*&
     &    tri3 + 1d0/8d0*s12*s23*mt*grm3a*tri1 - 1d0/4d0*s12*s23*mt*mZp**2*&
     &    grm3a*bo1 + 1d0/8d0*s12*s23**2*mt*grm3a*bo1 - 1d0/4d0*s12**2*&
     &    mt**3*grm3b*bo2 + 1d0/4d0*s12**2*mt**3*grm3a*bo1 - 1d0/4d0*s12**2&
     &    *s24*mt*grm3b*bo2 + 1d0/4d0*s12**2*s23*mt*grm3a*bo1 )
      resLL = resLL + spst5 * (  - 1d0/8d0*s12*s23*mt*grm3a*tri3 - 1d0/8d0*s12*s23*&
     &    mt*grm3a*tri2 + 1d0/4d0*s12*s23*mt*grm3a*tri1 - 1d0/8d0*s12*s23*&
     &    mt*mZp**2*grm3a*bo1 - 1d0/8d0*s12**2*mt*grm3a*tri6 + 1d0/8d0*&
     &    s12**2*mt*grm3a*tri1 + 1d0/8d0*s12**2*s23*mt*grm3a*bo1 )
      resLL = resLL + spst4 * (  - 1d0/4d0*s12*mt*grm3a*bub4 - 1d0/4d0*s12*mt*grm3a&
     &    *bub3 + 1d0/2d0*s12*mt*grm3a*bub2 + 1d0/4d0*s12*mt*mZp**2*grm3a*&
     &    tri2 + 1d0/4d0*s12*s23*mt*grm3a*tri3 + 1d0/4d0*s12*s23*mt*grm3a*&
     &    tri2 - 1d0/2d0*s12*s23*mt*grm3a*tri1 + 1d0/4d0*s12*s23*mt*mZp**2*&
     &    grm3a*bo1 + 1d0/4d0*s12**2*mt*grm3a*tri6 - 1d0/8d0*s12**2*mt*&
     &    grm3a*tri3 - 1d0/8d0*s12**2*mt*grm3a*tri2 - 1d0/4d0*s12**2*mt*&
     &    grm3a*tri1 - 1d0/4d0*s12**2*mt*mZp**2*grm3a*bo1 - 1d0/4d0*s12**2*&
     &    s23*mt*grm3a*bo1 + 1d0/16d0*s12**2*s23**2*mt*mZp**2*grm3a**2*&
     &    tri3 - 1d0/16d0*s12**2*s23**2*mt*mZp**2*grm3a**2*tri1 + 1d0/16d0*&
     &    s12**2*s23**2*mt*mZp**4*grm3a**2*bo1 - 1d0/16d0*s12**3*mt**3*&
     &    mZp**2*grm3a**2*tri2 + 1d0/16d0*s12**3*mt**3*mZp**2*grm3a**2*&
     &    tri1 + 1d0/16d0*s12**3*s23*mt*mZp**2*grm3a**2*tri6 - 1d0/16d0*&
     &    s12**3*s23*mt*mZp**2*grm3a**2*tri2 - 1d0/32d0*s12**3*s23**2*mt*&
     &    grm3a**2*tri3 - 1d0/32d0*s12**3*s23**2*mt*grm3a**2*tri2 + 1d0/16d0&
     &    *s12**3*s23**2*mt*grm3a**2*tri1 - 1d0/16d0*s12**3*s23**2*mt*&
     &    mZp**2*grm3a**2*bo1 )
      resLL = resLL + spst4 * ( 1d0/32d0*s12**4*mt**3*grm3a**2*tri3 + 1d0/32d0*&
     &    s12**4*mt**3*grm3a**2*tri2 - 1d0/16d0*s12**4*mt**3*grm3a**2*&
     &    tri1 + 1d0/16d0*s12**4*mt**3*mZp**2*grm3a**2*bo1 - 1d0/16d0*&
     &    s12**4*s23*mt*grm3a**2*tri6 + 1d0/32d0*s12**4*s23*mt*grm3a**2*&
     &    tri3 + 1d0/32d0*s12**4*s23*mt*grm3a**2*tri2 + 1d0/16d0*s12**4*s23&
     &    *mt*mZp**2*grm3a**2*bo1 + 1d0/16d0*s12**4*s23**2*mt*grm3a**2*&
     &    bo1 + 1d0/2d0/(4d0*mt**2 - s12)*s12*s23*mt*grm3a*bub4 + 1d0/2d0/(4d0&
     &    *mt**2 - s12)*s12*s23*mt*grm3a*bub3 - 1/(4d0*mt**2 - s12)*s12*&
     &    s23*mt*grm3a*bub1 + 1d0/2d0/(4d0*mt**2 - s12)*s12*s23*mt*mZp**2*&
     &    grm3a*tri1 + 1d0/4d0/(4d0*mt**2 - s12)*s12**2*mt*grm3a*bub4 + 1d0/&
     &    4d0/(4d0*mt**2 - s12)*s12**2*mt*grm3a*bub3 - 1d0/2d0/(4d0*mt**2 - &
     &    s12)*s12**2*mt*grm3a*bub1 + 1d0/4d0/(4d0*mt**2 - s12)*s12**2*mt*&
     &    mZp**2*grm3a*tri1 - 1d0/2d0/(4d0*mt**2 - s12)*s12**2*s23*mt*&
     &    grm3a*tri1 - 1d0/4d0/(4d0*mt**2 - s12)*s12**3*mt*grm3a*tri1 )
      resLL = resLL + spst3 * ( 1d0/4d0*s12*s23**(-1)*mt**2*grm3a*bub4 + 1d0/4d0*&
     &    s12*s23**(-1)*mt**2*grm3a*bub3 - 1d0/2d0*s12*s23**(-1)*mt**2*&
     &    grm3a*bub2 - 1d0/4d0*s12*s23**(-1)*mt**2*mZp**2*grm3a*tri2 - 1d0/&
     &    8d0*s12*grm3a*bub4 - 1d0/8d0*s12*grm3a*bub3 + 1d0/4d0*s12*grm3a*&
     &    bub1 - 1d0/8d0*s12*mZp**2*grm3a*tri6 - 1d0/4d0*s12*mt**2*grm3a*&
     &    tri3 - 1d0/4d0*s12*mt**2*grm3a*tri2 + 1d0/2d0*s12*mt**2*grm3a*&
     &    tri1 - 1d0/4d0*s12*mt**2*mZp**2*grm3a*bo1 + 1d0/4d0*s12*s23*grm3a&
     &    *tri6 - 1d0/4d0*s12*s23*grm3a*tri3 - 1d0/4d0*s12*s23*grm3a*tri2&
     &     - 1d0/4d0*s12*s23*grm3a*tri1 - 1d0/4d0*s12*s23*mZp**2*grm3a*bo1&
     &     - 1d0/4d0*s12*s23**2*grm3a*bo1 + 1d0/32d0*s12*s23**3*mZp**2*&
     &    grm3a**2*tri3 - 1d0/32d0*s12*s23**3*mZp**2*grm3a**2*tri2 + 1d0/&
     &    32d0*s12*s23**3*mZp**4*grm3a**2*bo1 + 1d0/4d0*s12**2*grm3a*tri6&
     &     - 1d0/8d0*s12**2*grm3a*tri3 - 1d0/8d0*s12**2*grm3a*tri2 - 1d0/4d0*&
     &    s12**2*grm3a*tri1 - 1d0/4d0*s12**2*mZp**2*grm3a*bo1 - 1d0/4d0*&
     &    s12**2*s23*grm3a*bo1 - 1d0/32d0*s12**2*s23*mt**2*mZp**2*&
     &    grm3a**2*tri3 )
      resLL = resLL + spst3 * (  - 3d0/32d0*s12**2*s23*mt**2*mZp**2*grm3a**2*tri2&
     &     + 1d0/8d0*s12**2*s23*mt**2*mZp**2*grm3a**2*tri1 - 1d0/32d0*&
     &    s12**2*s23*mt**2*mZp**4*grm3a**2*bo1 + 1d0/16d0*s12**2*s23**2*&
     &    mZp**2*grm3a**2*tri6 + 1d0/32d0*s12**2*s23**2*mZp**2*grm3a**2*&
     &    tri3 - 3d0/32d0*s12**2*s23**2*mZp**2*grm3a**2*tri2 + 1d0/32d0*&
     &    s12**2*s23**2*mZp**4*grm3a**2*bo1 - 1d0/32d0*s12**2*s23**3*&
     &    mZp**2*grm3a**2*bo1 - 1d0/16d0*s12**3*mt**2*mZp**2*grm3a**2*&
     &    tri2 + 1d0/16d0*s12**3*mt**2*mZp**2*grm3a**2*tri1 + 1d0/16d0*&
     &    s12**3*s23*mZp**2*grm3a**2*tri6 - 1d0/16d0*s12**3*s23*mZp**2*&
     &    grm3a**2*tri2 + 1d0/16d0*s12**3*s23*mt**2*grm3a**2*tri3 + 1d0/16d0&
     &    *s12**3*s23*mt**2*grm3a**2*tri2 - 1d0/8d0*s12**3*s23*mt**2*&
     &    grm3a**2*tri1 + 3d0/32d0*s12**3*s23*mt**2*mZp**2*grm3a**2*bo1&
     &     - 1d0/16d0*s12**3*s23**2*grm3a**2*tri6 + 1d0/32d0*s12**3*s23**2*&
     &    grm3a**2*tri3 + 1d0/32d0*s12**3*s23**2*grm3a**2*tri2 + 1d0/32d0*&
     &    s12**3*s23**2*mZp**2*grm3a**2*bo1 + 1d0/16d0*s12**3*s23**3*&
     &    grm3a**2*bo1 )
      resLL = resLL + spst3 * ( 1d0/32d0*s12**4*mt**2*grm3a**2*tri3 + 1d0/32d0*&
     &    s12**4*mt**2*grm3a**2*tri2 - 1d0/16d0*s12**4*mt**2*grm3a**2*&
     &    tri1 + 1d0/16d0*s12**4*mt**2*mZp**2*grm3a**2*bo1 - 1d0/16d0*&
     &    s12**4*s23*grm3a**2*tri6 + 1d0/32d0*s12**4*s23*grm3a**2*tri3 + &
     &    1d0/32d0*s12**4*s23*grm3a**2*tri2 + 1d0/16d0*s12**4*s23*mZp**2*&
     &    grm3a**2*bo1 + 1d0/16d0*s12**4*s23**2*grm3a**2*bo1 + 1d0/4d0/(4d0*&
     &    mt**2 - s12)*s12*s23*grm3a*bub4 + 1d0/4d0/(4d0*mt**2 - s12)*s12*&
     &    s23*grm3a*bub3 - 1d0/2d0/(4d0*mt**2 - s12)*s12*s23*grm3a*bub1 + &
     &    1d0/4d0/(4d0*mt**2 - s12)*s12*s23*mZp**2*grm3a*tri1 + 1/(4d0*&
     &    mt**2 - s12)*s12*s23*mt**2*grm3a*tri1 + 1d0/8d0/(4d0*mt**2 - s12&
     &    )*s12**2*grm3a*bub4 + 1d0/8d0/(4d0*mt**2 - s12)*s12**2*grm3a*&
     &    bub3 - 1d0/4d0/(4d0*mt**2 - s12)*s12**2*grm3a*bub1 + 1d0/8d0/(4d0*&
     &    mt**2 - s12)*s12**2*mZp**2*grm3a*tri1 + 1d0/2d0/(4d0*mt**2 - s12&
     &    )*s12**2*mt**2*grm3a*tri1 - 1d0/2d0/(4d0*mt**2 - s12)*s12**2*s23&
     &    *grm3a*tri1 - 1d0/4d0/(4d0*mt**2 - s12)*s12**3*grm3a*tri1 )
      resLL = resLL + spst2 * ( 2d0*bo1 - 1d0/4d0*s12*s23**(-1)*mt**2*grm3a*bub4&
     &     - 1d0/4d0*s12*s23**(-1)*mt**2*grm3a*bub3 + 1d0/2d0*s12*s23**(-1)&
     &    *mt**2*grm3a*bub2 + 1d0/4d0*s12*s23**(-1)*mt**2*mZp**2*grm3a*&
     &    tri2 - 1d0/8d0*s12*grm3a*bub4 - 1d0/8d0*s12*grm3a*bub3 + 1d0/2d0*&
     &    s12*grm3a*bub2 - 1d0/4d0*s12*grm3a*bub1 + 1d0/8d0*s12*mZp**2*&
     &    grm3a*tri6 + 1d0/4d0*s12*mZp**2*grm3a*tri2 + 1d0/4d0*s12*mt**2*&
     &    grm3a*tri3 + 1d0/4d0*s12*mt**2*grm3a*tri2 - 1d0/2d0*s12*mt**2*&
     &    grm3a*tri1 + 1d0/4d0*s12*mt**2*mZp**2*grm3a*bo1 - 1d0/4d0*s12*s23&
     &    *grm3a*tri6 + 1d0/2d0*s12*s23*grm3a*tri3 + 1d0/2d0*s12*s23*grm3a*&
     &    tri2 - 3d0/4d0*s12*s23*grm3a*tri1 + 3d0/4d0*s12*s23*mZp**2*grm3a*&
     &    bo1 - 1d0/4d0*s12*s23**2*grm3a*bo1 - 1d0/32d0*s12*s23**3*mZp**2*&
     &    grm3a**2*tri3 + 1d0/32d0*s12*s23**3*mZp**2*grm3a**2*tri2 - 1d0/&
     &    32d0*s12*s23**3*mZp**4*grm3a**2*bo1 - 1d0/4d0*s12**2*grm3a*tri1&
     &     - 1d0/2d0*s12**2*mt**2*grm3a*bo1 - 3d0/4d0*s12**2*s23*grm3a*bo1&
     &     + 1d0/32d0*s12**2*s23*mt**2*mZp**2*grm3a**2*tri3 + 3d0/32d0*&
     &    s12**2*s23*mt**2*mZp**2*grm3a**2*tri2 )
      resLL = resLL + spst2 * (  - 1d0/8d0*s12**2*s23*mt**2*mZp**2*grm3a**2*tri1&
     &     + 1d0/32d0*s12**2*s23*mt**2*mZp**4*grm3a**2*bo1 - 1d0/16d0*&
     &    s12**2*s23**2*mZp**2*grm3a**2*tri6 + 1d0/32d0*s12**2*s23**2*&
     &    mZp**2*grm3a**2*tri3 + 3d0/32d0*s12**2*s23**2*mZp**2*grm3a**2*&
     &    tri2 - 1d0/16d0*s12**2*s23**2*mZp**2*grm3a**2*tri1 + 1d0/32d0*&
     &    s12**2*s23**2*mZp**4*grm3a**2*bo1 - 1d0/32d0*s12**2*s23**3*&
     &    mZp**2*grm3a**2*bo1 - 1d0/16d0*s12**3*s23*mt**2*grm3a**2*tri3&
     &     - 1d0/16d0*s12**3*s23*mt**2*grm3a**2*tri2 + 1d0/8d0*s12**3*s23*&
     &    mt**2*grm3a**2*tri1 - 5d0/32d0*s12**3*s23*mt**2*mZp**2*grm3a**2&
     &    *bo1 + 1d0/16d0*s12**3*s23**2*grm3a**2*tri6 - 1d0/16d0*s12**3*&
     &    s23**2*grm3a**2*tri3 - 1d0/16d0*s12**3*s23**2*grm3a**2*tri2 + 1d0&
      &   /16d0*s12**3*s23**2*grm3a**2*tri1 - 5d0/32d0*s12**3*s23**2*&
     &    mZp**2*grm3a**2*bo1 + 1d0/16d0*s12**4*s23*mt**2*grm3a**2*bo1 + &
     &    1d0/16d0*s12**4*s23**2*grm3a**2*bo1 + 1d0/4d0/(4d0*mt**2 - s12)*&
     &    s12*s23*grm3a*bub4 + 1d0/4d0/(4d0*mt**2 - s12)*s12*s23*grm3a*&
     &    bub3 )
      resLL = resLL + spst2 * (  - 1d0/2d0/(4d0*mt**2 - s12)*s12*s23*grm3a*bub1 + &
     &    1d0/4d0/(4d0*mt**2 - s12)*s12*s23*mZp**2*grm3a*tri1 + 1/(4d0*&
     &    mt**2 - s12)*s12*s23*mt**2*grm3a*tri1 + 1d0/8d0/(4d0*mt**2 - s12&
     &    )*s12**2*grm3a*bub4 + 1d0/8d0/(4d0*mt**2 - s12)*s12**2*grm3a*&
     &    bub3 - 1d0/4d0/(4d0*mt**2 - s12)*s12**2*grm3a*bub1 + 1d0/8d0/(4d0*&
     &    mt**2 - s12)*s12**2*mZp**2*grm3a*tri1 + 1d0/2d0/(4d0*mt**2 - s12&
     &    )*s12**2*mt**2*grm3a*tri1 - 1d0/2d0/(4d0*mt**2 - s12)*s12**2*s23&
     &    *grm3a*tri1 - 1d0/4d0/(4d0*mt**2 - s12)*s12**3*grm3a*tri1 )
      resLL = resLL + spst1 * (  - 1d0/8d0*s24**2*mt*grm3b*tri5 + 1d0/8d0*s24**2*mt&
     &    *grm3b*tri4 + 1d0/8d0*s24**2*mt*mZp**2*grm3b*bo2 + 1d0/8d0*s23**2&
     &    *mt*grm3a*tri3 - 1d0/8d0*s23**2*mt*grm3a*tri2 + 1d0/8d0*s23**2*mt&
     &    *mZp**2*grm3a*bo1 - 1d0/4d0*s12*mt**3*grm3b*tri5 + 1d0/4d0*s12*&
     &    mt**3*grm3b*tri1 + 1d0/4d0*s12*mt**3*grm3a*tri3 - 1d0/4d0*s12*&
     &    mt**3*grm3a*tri1 + 1d0/4d0*s12*mt**3*mZp**2*grm3a*bo1 + 1d0/8d0*&
     &    s12*s24*mt*grm3b*tri6 - 1d0/4d0*s12*s24*mt*grm3b*tri5 + 1d0/8d0*&
     &    s12*s24*mt*grm3b*tri1 - 1d0/8d0*s12*s24**2*mt*grm3b*bo2 - 1d0/8d0&
     &    *s12*s23*mt*grm3a*tri6 + 1d0/8d0*s12*s23*mt*grm3a*tri3 - 1d0/8d0*&
     &    s12*s23*mt*grm3a*tri2 + 1d0/8d0*s12*s23*mt*grm3a*tri1 + 1d0/8d0*&
     &    s12*s23*mt*mZp**2*grm3a*bo1 + 1d0/8d0*s12*s23**2*mt*grm3a*bo1&
     &     - 1d0/8d0*s12**2*mt*grm3a*tri6 + 1d0/8d0*s12**2*mt*grm3a*tri1 + &
     &    1d0/8d0*s12**2*s23*mt*grm3a*bo1 )


    resLR =  + spst11 * (  - 1d0/2d0*tri5 - 1d0/2d0*tri4 - mZp**2*bo2 - 1d0/&
     &    8d0*s24**2*mZp**2*grm3b*tri5 + 1d0/8d0*s24**2*mZp**2*grm3b*tri4&
     &     + 1d0/8d0*s24**2*mZp**4*grm3b*bo2 - 1d0/4d0*s12*mt**2*mZp**2*&
     &    grm3b*tri5 + 1d0/4d0*s12*mt**2*mZp**2*grm3b*tri1 + 1d0/8d0*s12*&
     &    s24*mZp**2*grm3b*tri6 - 1d0/4d0*s12*s24*mZp**2*grm3b*tri5 + 1d0/&
     &    8d0*s12*s24*mZp**2*grm3b*tri1 + 1d0/8d0*s12**2*mt**2*grm3b*tri5&
     &     + 1d0/8d0*s12**2*mt**2*grm3b*tri4 - 1d0/4d0*s12**2*mt**2*grm3b*&
     &    tri1 + 1d0/4d0*s12**2*mt**2*mZp**2*grm3b*bo2 - 1d0/8d0*s12**2*s24&
     &    *grm3b*tri6 + 1d0/8d0*s12**2*s24*grm3b*tri5 + 1d0/8d0*s12**2*s24*&
     &    grm3b*tri4 - 1d0/8d0*s12**2*s24*grm3b*tri1 + 1d0/4d0*s12**2*s24*&
     &    mZp**2*grm3b*bo2 + 1d0/8d0*s12**2*s24**2*grm3b*bo2 )
      resLR = resLR + spst10 * (  - 1d0/2d0*s12*s24**(-1)*mt**2*grm3b*bub5 + 1d0/4d0&
     &    *s12*s24**(-1)*mt**2*grm3b*bub4 + 1d0/4d0*s12*s24**(-1)*mt**2*&
     &    grm3b*bub3 - 1d0/4d0*s12*s24**(-1)*mt**2*mZp**2*grm3b*tri5 - 1d0/&
     &    8d0*s12*grm3b*bub4 - 1d0/8d0*s12*grm3b*bub3 + 1d0/4d0*s12*grm3b*&
     &    bub1 - 1d0/8d0*s12*mZp**2*grm3b*tri6 - 1d0/4d0*s12*mt**2*grm3b*&
     &    tri5 - 1d0/4d0*s12*mt**2*grm3b*tri4 + 1d0/2d0*s12*mt**2*grm3b*&
     &    tri1 - 1d0/4d0*s12*mt**2*mZp**2*grm3b*bo2 + 1d0/4d0*s12*s24*grm3b&
     &    *tri6 - 1d0/4d0*s12*s24*grm3b*tri5 - 1d0/4d0*s12*s24*grm3b*tri4&
     &     - 1d0/4d0*s12*s24*grm3b*tri1 - 1d0/4d0*s12*s24*mZp**2*grm3b*bo2&
     &     - 1d0/4d0*s12*s24**2*grm3b*bo2 - 1d0/32d0*s12*s24**3*mZp**2*&
     &    grm3b**2*tri5 + 1d0/32d0*s12*s24**3*mZp**2*grm3b**2*tri4 + 1d0/&
     &    32d0*s12*s24**3*mZp**4*grm3b**2*bo2 + 1d0/4d0*s12**2*grm3b*tri6&
     &     - 1d0/8d0*s12**2*grm3b*tri5 - 1d0/8d0*s12**2*grm3b*tri4 - 1d0/4d0*&
     &    s12**2*grm3b*tri1 - 1d0/4d0*s12**2*mZp**2*grm3b*bo2 - 1d0/4d0*&
     &    s12**2*s24*grm3b*bo2 - 3d0/32d0*s12**2*s24*mt**2*mZp**2*&
     &    grm3b**2*tri5 )
      resLR = resLR + spst10 * (  - 1d0/32d0*s12**2*s24*mt**2*mZp**2*grm3b**2*&
     &    tri4 + 1d0/8d0*s12**2*s24*mt**2*mZp**2*grm3b**2*tri1 - 1d0/32d0*&
     &    s12**2*s24*mt**2*mZp**4*grm3b**2*bo2 + 1d0/16d0*s12**2*s24**2*&
     &    mZp**2*grm3b**2*tri6 - 3d0/32d0*s12**2*s24**2*mZp**2*grm3b**2*&
     &    tri5 + 1d0/32d0*s12**2*s24**2*mZp**2*grm3b**2*tri4 + 1d0/32d0*&
     &    s12**2*s24**2*mZp**4*grm3b**2*bo2 - 1d0/32d0*s12**2*s24**3*&
     &    mZp**2*grm3b**2*bo2 - 1d0/16d0*s12**3*mt**2*mZp**2*grm3b**2*&
     &    tri5 + 1d0/16d0*s12**3*mt**2*mZp**2*grm3b**2*tri1 + 1d0/16d0*&
     &    s12**3*s24*mZp**2*grm3b**2*tri6 - 1d0/16d0*s12**3*s24*mZp**2*&
     &    grm3b**2*tri5 + 1d0/16d0*s12**3*s24*mt**2*grm3b**2*tri5 + 1d0/16d0&
     &    *s12**3*s24*mt**2*grm3b**2*tri4 - 1d0/8d0*s12**3*s24*mt**2*&
     &    grm3b**2*tri1 + 3d0/32d0*s12**3*s24*mt**2*mZp**2*grm3b**2*bo2&
     &     - 1d0/16d0*s12**3*s24**2*grm3b**2*tri6 + 1d0/32d0*s12**3*s24**2*&
     &    grm3b**2*tri5 + 1d0/32d0*s12**3*s24**2*grm3b**2*tri4 + 1d0/32d0*&
     &    s12**3*s24**2*mZp**2*grm3b**2*bo2 + 1d0/16d0*s12**3*s24**3*&
     &    grm3b**2*bo2 )
      resLR = resLR + spst10 * ( 1d0/32d0*s12**4*mt**2*grm3b**2*tri5 + 1d0/32d0*&
     &    s12**4*mt**2*grm3b**2*tri4 - 1d0/16d0*s12**4*mt**2*grm3b**2*&
     &    tri1 + 1d0/16d0*s12**4*mt**2*mZp**2*grm3b**2*bo2 - 1d0/16d0*&
     &    s12**4*s24*grm3b**2*tri6 + 1d0/32d0*s12**4*s24*grm3b**2*tri5 + &
     &    1d0/32d0*s12**4*s24*grm3b**2*tri4 + 1d0/16d0*s12**4*s24*mZp**2*&
     &    grm3b**2*bo2 + 1d0/16d0*s12**4*s24**2*grm3b**2*bo2 + 1d0/4d0/(4d0*&
     &    mt**2 - s12)*s12*s24*grm3b*bub4 + 1d0/4d0/(4d0*mt**2 - s12)*s12*&
     &    s24*grm3b*bub3 - 1d0/2d0/(4d0*mt**2 - s12)*s12*s24*grm3b*bub1 + &
     &    1d0/4d0/(4d0*mt**2 - s12)*s12*s24*mZp**2*grm3b*tri1 + 1/(4d0*&
     &    mt**2 - s12)*s12*s24*mt**2*grm3b*tri1 + 1d0/8d0/(4d0*mt**2 - s12&
     &    )*s12**2*grm3b*bub4 + 1d0/8d0/(4d0*mt**2 - s12)*s12**2*grm3b*&
     &    bub3 - 1d0/4d0/(4d0*mt**2 - s12)*s12**2*grm3b*bub1 + 1d0/8d0/(4d0*&
     &    mt**2 - s12)*s12**2*mZp**2*grm3b*tri1 + 1d0/2d0/(4d0*mt**2 - s12&
     &    )*s12**2*mt**2*grm3b*tri1 - 1d0/2d0/(4d0*mt**2 - s12)*s12**2*s24&
     &    *grm3b*tri1 - 1d0/4d0/(4d0*mt**2 - s12)*s12**3*grm3b*tri1 )
      resLR = resLR + spst9 * ( 2d0*bo2 + 1d0/2d0*s12*s24**(-1)*mt**2*grm3b*bub5&
     &     - 1d0/4d0*s12*s24**(-1)*mt**2*grm3b*bub4 - 1d0/4d0*s12*s24**(-1)&
     &    *mt**2*grm3b*bub3 + 1d0/4d0*s12*s24**(-1)*mt**2*mZp**2*grm3b*&
     &    tri5 + 1d0/2d0*s12*grm3b*bub5 - 1d0/8d0*s12*grm3b*bub4 - 1d0/8d0*&
     &    s12*grm3b*bub3 - 1d0/4d0*s12*grm3b*bub1 + 1d0/8d0*s12*mZp**2*&
     &    grm3b*tri6 + 1d0/4d0*s12*mZp**2*grm3b*tri5 + 1d0/4d0*s12*mt**2*&
     &    grm3b*tri5 + 1d0/4d0*s12*mt**2*grm3b*tri4 - 1d0/2d0*s12*mt**2*&
     &    grm3b*tri1 + 1d0/4d0*s12*mt**2*mZp**2*grm3b*bo2 - 1d0/4d0*s12*s24&
     &    *grm3b*tri6 + 1d0/2d0*s12*s24*grm3b*tri5 + 1d0/2d0*s12*s24*grm3b*&
     &    tri4 - 3d0/4d0*s12*s24*grm3b*tri1 + 3d0/4d0*s12*s24*mZp**2*grm3b*&
     &    bo2 - 1d0/4d0*s12*s24**2*grm3b*bo2 + 1d0/32d0*s12*s24**3*mZp**2*&
     &    grm3b**2*tri5 - 1d0/32d0*s12*s24**3*mZp**2*grm3b**2*tri4 - 1d0/&
     &    32d0*s12*s24**3*mZp**4*grm3b**2*bo2 - 1d0/4d0*s12**2*grm3b*tri1&
     &     - 1d0/2d0*s12**2*mt**2*grm3b*bo2 - 3d0/4d0*s12**2*s24*grm3b*bo2&
     &     + 3d0/32d0*s12**2*s24*mt**2*mZp**2*grm3b**2*tri5 + 1d0/32d0*&
     &    s12**2*s24*mt**2*mZp**2*grm3b**2*tri4 )
      resLR = resLR + spst9 * (  - 1d0/8d0*s12**2*s24*mt**2*mZp**2*grm3b**2*tri1&
     &     + 1d0/32d0*s12**2*s24*mt**2*mZp**4*grm3b**2*bo2 - 1d0/16d0*&
     &    s12**2*s24**2*mZp**2*grm3b**2*tri6 + 3d0/32d0*s12**2*s24**2*&
     &    mZp**2*grm3b**2*tri5 + 1d0/32d0*s12**2*s24**2*mZp**2*grm3b**2*&
     &    tri4 - 1d0/16d0*s12**2*s24**2*mZp**2*grm3b**2*tri1 + 1d0/32d0*&
     &    s12**2*s24**2*mZp**4*grm3b**2*bo2 - 1d0/32d0*s12**2*s24**3*&
     &    mZp**2*grm3b**2*bo2 - 1d0/16d0*s12**3*s24*mt**2*grm3b**2*tri5&
     &     - 1d0/16d0*s12**3*s24*mt**2*grm3b**2*tri4 + 1d0/8d0*s12**3*s24*&
     &    mt**2*grm3b**2*tri1 - 5d0/32d0*s12**3*s24*mt**2*mZp**2*grm3b**2&
     &    *bo2 + 1d0/16d0*s12**3*s24**2*grm3b**2*tri6 - 1d0/16d0*s12**3*&
     &    s24**2*grm3b**2*tri5 - 1d0/16d0*s12**3*s24**2*grm3b**2*tri4 + 1d0&
      &   /16d0*s12**3*s24**2*grm3b**2*tri1 - 5d0/32d0*s12**3*s24**2*&
     &    mZp**2*grm3b**2*bo2 + 1d0/16d0*s12**4*s24*mt**2*grm3b**2*bo2 + &
     &    1d0/16d0*s12**4*s24**2*grm3b**2*bo2 + 1d0/4d0/(4d0*mt**2 - s12)*&
     &    s12*s24*grm3b*bub4 + 1d0/4d0/(4d0*mt**2 - s12)*s12*s24*grm3b*&
     &    bub3 )
      resLR = resLR + spst9 * (  - 1d0/2d0/(4d0*mt**2 - s12)*s12*s24*grm3b*bub1 + &
     &    1d0/4d0/(4d0*mt**2 - s12)*s12*s24*mZp**2*grm3b*tri1 + 1/(4d0*&
     &    mt**2 - s12)*s12*s24*mt**2*grm3b*tri1 + 1d0/8d0/(4d0*mt**2 - s12&
     &    )*s12**2*grm3b*bub4 + 1d0/8d0/(4d0*mt**2 - s12)*s12**2*grm3b*&
     &    bub3 - 1d0/4d0/(4d0*mt**2 - s12)*s12**2*grm3b*bub1 + 1d0/8d0/(4d0*&
     &    mt**2 - s12)*s12**2*mZp**2*grm3b*tri1 + 1d0/2d0/(4d0*mt**2 - s12&
     &    )*s12**2*mt**2*grm3b*tri1 - 1d0/2d0/(4d0*mt**2 - s12)*s12**2*s24&
     &    *grm3b*tri1 - 1d0/4d0/(4d0*mt**2 - s12)*s12**3*grm3b*tri1 )
      resLR = resLR + spst8 * ( 1d0/4d0*s12*s24*mt**2*grm3b*tri5 + 1d0/4d0*s12*s24*&
     &    mt**2*grm3b*tri4 - 1d0/2d0*s12*s24*mt**2*grm3b*tri1 + 1d0/4d0*s12&
     &    *s24*mt**2*mZp**2*grm3b*bo2 + 1d0/4d0*s12**2*mt**2*grm3b*tri6&
     &     - 1d0/4d0*s12**2*mt**2*grm3b*tri1 - 1d0/4d0*s12**2*s24*mt**2*&
     &    grm3b*bo2 )
      resLR = resLR + spst7 * ( 2d0*tri6 + 2d0*tri1 - 2d0*s23*bo1 - 1d0/8d0*s12*s23*&
     &    mt**2*grm3a*tri3 - 1d0/8d0*s12*s23*mt**2*grm3a*tri2 + 1d0/4d0*s12&
     &    *s23*mt**2*grm3a*tri1 - 1d0/8d0*s12*s23*mt**2*mZp**2*grm3a*bo1&
     &     - 1d0/8d0*s12**2*mt**2*grm3a*tri6 + 1d0/8d0*s12**2*mt**2*grm3a*&
     &    tri1 + 1d0/8d0*s12**2*s23*mt**2*grm3a*bo1 )
      resLR = resLR + spst6 * ( mt*bo2 - mt*bo1 - 1d0/8d0*s24**2*mt*grm3b*tri5 + &
     &    1d0/8d0*s24**2*mt*grm3b*tri4 + 1d0/8d0*s24**2*mt*mZp**2*grm3b*bo2&
     &     + 1d0/8d0*s23**2*mt*grm3a*tri3 - 1d0/8d0*s23**2*mt*grm3a*tri2 + &
     &    1d0/8d0*s23**2*mt*mZp**2*grm3a*bo1 + 1d0/4d0*s12*mt**3*grm3b*tri4&
     &     - 1d0/4d0*s12*mt**3*grm3b*tri1 - 1d0/4d0*s12*mt**3*grm3a*tri2 + &
     &    1d0/4d0*s12*mt**3*grm3a*tri1 + 1d0/4d0*s12*mt**3*mZp**2*grm3b*bo2&
     &     - 1d0/8d0*s12*s24*mt*grm3b*tri6 + 1d0/4d0*s12*s24*mt*grm3b*tri4&
     &     - 1d0/8d0*s12*s24*mt*grm3b*tri1 + 1d0/4d0*s12*s24*mt*mZp**2*&
     &    grm3b*bo2 - 1d0/8d0*s12*s24**2*mt*grm3b*bo2 + 1d0/8d0*s12*s23*mt*&
     &    grm3a*tri6 - 1d0/4d0*s12*s23*mt*grm3a*tri2 + 1d0/8d0*s12*s23*mt*&
     &    grm3a*tri1 + 1d0/8d0*s12*s23**2*mt*grm3a*bo1 - 1d0/4d0*s12**2*&
     &    mt**3*grm3b*bo2 + 1d0/4d0*s12**2*mt**3*grm3a*bo1 - 1d0/4d0*s12**2&
     &    *s24*mt*grm3b*bo2 + 1d0/4d0*s12**2*s23*mt*grm3a*bo1 )
      resLR = resLR + spst5 * (  - 1d0/8d0*s12*s23*mt*grm3a*tri3 - 1d0/8d0*s12*s23*&
     &    mt*grm3a*tri2 + 1d0/4d0*s12*s23*mt*grm3a*tri1 - 1d0/8d0*s12*s23*&
     &    mt*mZp**2*grm3a*bo1 - 1d0/8d0*s12**2*mt*grm3a*tri6 + 1d0/8d0*&
     &    s12**2*mt*grm3a*tri1 + 1d0/8d0*s12**2*s23*mt*grm3a*bo1 )
      resLR = resLR + spst4 * (  - 1d0/2d0*s12*mt*grm3b*bub5 + 1d0/4d0*s12*mt*grm3b&
     &    *bub4 + 1d0/4d0*s12*mt*grm3b*bub3 - 1d0/4d0*s12*mt*mZp**2*grm3b*&
     &    tri5 - 1d0/4d0*s12*s24*mt*grm3b*tri5 - 1d0/4d0*s12*s24*mt*grm3b*&
     &    tri4 + 1d0/2d0*s12*s24*mt*grm3b*tri1 - 1d0/4d0*s12*s24*mt*mZp**2*&
     &    grm3b*bo2 - 1d0/4d0*s12**2*mt*grm3b*tri6 + 1d0/8d0*s12**2*mt*&
     &    grm3b*tri5 + 1d0/8d0*s12**2*mt*grm3b*tri4 + 1d0/4d0*s12**2*mt*&
     &    grm3b*tri1 + 1d0/4d0*s12**2*mt*mZp**2*grm3b*bo2 + 1d0/4d0*s12**2*&
     &    s24*mt*grm3b*bo2 - 1d0/16d0*s12**2*s24**2*mt*mZp**2*grm3b**2*&
     &    tri4 + 1d0/16d0*s12**2*s24**2*mt*mZp**2*grm3b**2*tri1 - 1d0/16d0*&
     &    s12**2*s24**2*mt*mZp**4*grm3b**2*bo2 + 1d0/16d0*s12**3*mt**3*&
     &    mZp**2*grm3b**2*tri5 - 1d0/16d0*s12**3*mt**3*mZp**2*grm3b**2*&
     &    tri1 - 1d0/16d0*s12**3*s24*mt*mZp**2*grm3b**2*tri6 + 1d0/16d0*&
     &    s12**3*s24*mt*mZp**2*grm3b**2*tri5 + 1d0/32d0*s12**3*s24**2*mt*&
     &    grm3b**2*tri5 + 1d0/32d0*s12**3*s24**2*mt*grm3b**2*tri4 - 1d0/16d0&
     &    *s12**3*s24**2*mt*grm3b**2*tri1 + 1d0/16d0*s12**3*s24**2*mt*&
     &    mZp**2*grm3b**2*bo2 )
      resLR = resLR + spst4 * (  - 1d0/32d0*s12**4*mt**3*grm3b**2*tri5 - 1d0/32d0*&
     &    s12**4*mt**3*grm3b**2*tri4 + 1d0/16d0*s12**4*mt**3*grm3b**2*&
     &    tri1 - 1d0/16d0*s12**4*mt**3*mZp**2*grm3b**2*bo2 + 1d0/16d0*&
     &    s12**4*s24*mt*grm3b**2*tri6 - 1d0/32d0*s12**4*s24*mt*grm3b**2*&
     &    tri5 - 1d0/32d0*s12**4*s24*mt*grm3b**2*tri4 - 1d0/16d0*s12**4*s24&
     &    *mt*mZp**2*grm3b**2*bo2 - 1d0/16d0*s12**4*s24**2*mt*grm3b**2*&
     &    bo2 - 1d0/2d0/(4d0*mt**2 - s12)*s12*s24*mt*grm3b*bub4 - 1d0/2d0/(4d0&
     &    *mt**2 - s12)*s12*s24*mt*grm3b*bub3 + 1/(4d0*mt**2 - s12)*s12*&
     &    s24*mt*grm3b*bub1 - 1d0/2d0/(4d0*mt**2 - s12)*s12*s24*mt*mZp**2*&
     &    grm3b*tri1 - 1d0/4d0/(4d0*mt**2 - s12)*s12**2*mt*grm3b*bub4 - 1d0/&
     &    4d0/(4d0*mt**2 - s12)*s12**2*mt*grm3b*bub3 + 1d0/2d0/(4d0*mt**2 - &
     &    s12)*s12**2*mt*grm3b*bub1 - 1d0/4d0/(4d0*mt**2 - s12)*s12**2*mt*&
     &    mZp**2*grm3b*tri1 + 1d0/2d0/(4d0*mt**2 - s12)*s12**2*s24*mt*&
     &    grm3b*tri1 + 1d0/4d0/(4d0*mt**2 - s12)*s12**3*mt*grm3b*tri1 )
      resLR = resLR + spst1 * ( 1d0/8d0*s24**2*mt*grm3b*tri5 - 1d0/8d0*s24**2*mt*&
     &    grm3b*tri4 - 1d0/8d0*s24**2*mt*mZp**2*grm3b*bo2 - 1d0/8d0*s23**2*&
     &    mt*grm3a*tri3 + 1d0/8d0*s23**2*mt*grm3a*tri2 - 1d0/8d0*s23**2*mt*&
     &    mZp**2*grm3a*bo1 - 1d0/4d0*s12*mt**3*grm3b*tri4 + 1d0/4d0*s12*&
     &    mt**3*grm3b*tri1 + 1d0/4d0*s12*mt**3*grm3a*tri2 - 1d0/4d0*s12*&
     &    mt**3*grm3a*tri1 - 1d0/4d0*s12*mt**3*mZp**2*grm3b*bo2 + 1d0/8d0*&
     &    s12*s24*mt*grm3b*tri6 - 1d0/4d0*s12*s24*mt*grm3b*tri4 + 1d0/8d0*&
     &    s12*s24*mt*grm3b*tri1 - 1d0/4d0*s12*s24*mt*mZp**2*grm3b*bo2 - 1d0&
      &   /8d0*s12*s24**2*mt*grm3b*bo2 - 1d0/8d0*s12*s23*mt*grm3a*tri6 - 1d0&
      &   /8d0*s12*s23*mt*grm3a*tri3 + 1d0/8d0*s12*s23*mt*grm3a*tri2 + 1d0/&
     &    8d0*s12*s23*mt*grm3a*tri1 - 1d0/8d0*s12*s23*mt*mZp**2*grm3a*bo1&
     &     + 1d0/8d0*s12*s23**2*mt*grm3a*bo1 - 1d0/8d0*s12**2*mt*grm3a*tri6&
     &     + 1d0/8d0*s12**2*mt*grm3a*tri1 + 1d0/8d0*s12**2*s23*mt*grm3a*bo1&
     &     )


    resLL = resLL * (cI)/(4d0 * Pi)**2 ! QCDLoop normalization
    resLL = resLL * (2d0) ! Overall normalization

    resLR = resLR * (cI)/(4d0 * Pi)**2 ! QCDLoop normalization
    resLR = resLR * (2d0) ! Overall normalization

    boxL = resLL(-1) * gL_Zpr(top_) + resLR(-1) * gR_Zpr(top_)
    boxR = resLL(+1) * gR_Zpr(top_) + resLR(+1) * gL_Zpr(top_)

!   We don't use a complex mass-scheme -> 
!   Impose the 'right' s-channel propagator by hand

    boxL = boxL * (s12-mZp**2)/(s12-mZp**2+ci*Ga_Zpr*mZp)
    boxR = boxR * (s12-mZp**2)/(s12-mZp**2+ci*Ga_Zpr*mZp)

!   Singularity checks, Catani

    if ( epv .ne. 0) then
       catani = -2d0/(4d0*Pi)**2 * dlog(s13/s23)
       call Tree_Zprime_tbtqbq(LO_L,LO_R)
       print *,'tree, up quark ',(LO_L*gL_Zpr(up_)+LO_R*gR_Zpr(up_))
       print *,'res, up quark  ',boxL*gL_Zpr(up_)+boxR*gR_Zpr(up_)
       if (epv.eq.-1) then
          print *,'catani , up q  ',(LO_L*gL_Zpr(up_)+LO_R*gR_Zpr(up_))*catani
       endif
       print *, ''
    endif

    return
  end subroutine Virt_Zprime_box




  function spaa(ub,v)
    complex(8) :: spaa, ub(4), v(4)
    
    spaa = psp1_(ub,chir(.true.,v))
    
    return
  end function spaa
          
  function spbb(ub,v)
    complex(8) :: spbb, ub(4), v(4)
    
    spbb = psp1_(ub,chir(.false.,v))
    
    return
  end function spbb

  function spab(ub,p,v)
    complex(8) :: spab, ub(4), v(4), p(4)
    
    spab = psp1_(ub,spi2_(p,chir(.false.,v)))
    return

  end function spab
          
  function spba(ub,p,v)
    complex(8) :: spba, ub(4), v(4), p(4)
    
    spba = psp1_(ub,spi2_(p,chir(.true.,v)))
    return
    
  end function spba

  function spab_4V(ub,v)
    complex(8) :: spab_4V(4), ub(4), v(4)
    
    spab_4V = vbqq2(ub,chir(.false.,v))
    return
            
  end function spab_4V

  function spba_4V(ub,v)
    complex(8) :: spba_4V(4), ub(4), v(4)
    
    spba_4V = vbqq2(ub,chir(.true.,v))
    return
            
  end function spba_4V


END MODULE ModAmplitudes_Zprime

