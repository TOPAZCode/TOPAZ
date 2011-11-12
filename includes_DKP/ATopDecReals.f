
      subroutine numDecATopReals(heli,k1bb,klep,kneu,k3g,k4p,gauge,res)
      implicit none
      include 'realcommon.f'
      include 'VarAReals.f'
      integer heli(4), heli_pho, heli_glu, i, DV, gauge, heli_w,heli_b
      double complex p1tb(4),k2w(4),k1bb(4),k3g(4),k4p(4)
      double complex kneu(4), klep(4)
      double complex EpsP(4),EpsW(4), Spi_bot(4), EpsG(4)
      double complex res(4), SMEcoeff(24), SME(24,1:4), scf
      double complex k1bbp1tb,k1bbk2w,k2wp1tb
      double complex k1bbk3g, k2wk3g, k3gp1tb
      double complex EpsWk1bb, EpsWp1tb, EpsWk3g, EpsWEpsP, EpsWEpsG
      double complex EpsPk1bb, EpsPp1tb, EpsPk2w, EpsPk3g, EpsPEpsG
      double complex EpsGk1bb, EpsGk2w, EpsGp1tb
      double complex one, two, four, three
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11
      double precision GS
      Dv =4
      gs = dsqrt(4.d0*pii*alphaS)

      one = dcmplx(1.d0,0.d0)
      two = dcmplx(2.d0,0.d0)
      four = dcmplx(4.d0,0.d0)
      three = dcmplx(3.d0,0.d0)
      heli_pho = heli(3)
      heli_glu = heli(2)
      heli_w = heli(4)	
      heli_b = +1
      k2w(1:4) = kneu(1:4) + klep(1:4)
      p1tb(1:4) = k2w(1:4) + k1bb(1:4) + k3g(1:4) + k4p(1:4)
      
      call vSpi_Weyl(k1bb,heli_b,Spi_bot)
      call pol_massoff2(klep,kneu,EpsW)
C      call pol_mass(k2w,Mw,heli_w,EpsW)
      call pol_mless(k4p,heli_pho,EpsP)
      call pol_mless(k3g,heli_glu,EpsG)

      If (gauge .eq. 1) then
         EpsP(1:4) = k4p(1:4)
      endif
      If (gauge .eq. 2) then
         EpsG(1:4) = k3g(1:4)
      endif

      k1bbp1tb = scf(Dv,k1bb,p1tb)
      k1bbk2w = scf(Dv,k1bb,k2w)
      k1bbk3g = scf(Dv,k1bb,k3g)

      k2wp1tb = scf(Dv,k2w,p1tb)
      k2wk3g = scf(Dv,k2w,k3g)

      k3gp1tb = scf(Dv,p1tb,k3g)


      EpsWk1bb = scf(Dv,EpsW,k1bb)
      EpsWp1tb = scf(Dv,EpsW,p1tb)
      EpsWk3g = scf(Dv,EpsW,k3g)
      EpsWEpsG = scf(Dv,EpsW,EpsG)
      EpsWEpsP = scf(Dv,EpsW,EpsP)


      EpsPk1bb = scf(Dv,EpsP,k1bb)
      EpsPp1tb = scf(Dv,EpsP,p1tb)
      EpsPk2w = scf(Dv,EpsP,k2w)
      EpsPk3g = scf(Dv,EpsP,k3g)
      EpsPEpsG = scf(Dv,EpsP,EpsG)


      EpsGk1bb = scf(Dv,EpsG,k1bb)
      EpsGp1tb = scf(Dv,EpsG,p1tb)
      EpsGk2w = scf(Dv,EpsG,k2w)


      

      call SpATop0(Spi_bot,p1tb,mt,SME(1,1:4))

      call SpATop1(Spi_bot,EpsW,p1tb,mt,SME(2,1:4))
      call SpATop1(Spi_bot,EpsP,p1tb,mt,SME(4,1:4))
      call SpATop1(Spi_bot,k2w,p1tb,mt,SME(3,1:4))
      call SpATop1(Spi_bot,EpsG,p1tb,mt,SME(5,1:4))
      call SpATop1(Spi_bot,k3g,p1tb,mt,SME(6,1:4))

      call SpATop2(Spi_bot,EpsW,EpsP,p1tb,mt,SME(8,1:4))
      call SpATop2(Spi_bot,EpsW,EpsG,p1tb,mt,SME(7,1:4))
      call SpATop2(Spi_bot,EpsW,k3g,p1tb,mt,SME(9,1:4))

      call SpATop2(Spi_bot,EpsP,EpsG,p1tb,mt,SME(10,1:4))
      call SpATop2(Spi_bot,EpsP,k3g,p1tb,mt,SME(11,1:4))

      call SpATop2(Spi_bot,EpsG,k3g,p1tb,mt,SME(12,1:4))

      call SpATop3(Spi_bot,EpsW,k2w,EpsP,p1tb,mt,SME(13,1:4))
      call SpATop3(Spi_bot,EpsW,k2w,EpsG,p1tb,mt,SME(14,1:4))
      call SpATop3(Spi_bot,EpsW,k2w,k3g,p1tb,mt,SME(15,1:4))
      call SpATop3(Spi_bot,EpsW,EpsP,EpsG,p1tb,mt,SME(16,1:4))
      call SpATop3(Spi_bot,EpsW,EpsP,k3g,p1tb,mt,SME(17,1:4))
      call SpATop3(Spi_bot,EpsW,EpsG,k3g,p1tb,mt,SME(18,1:4))

      call SpATop3(Spi_bot,k2w,EpsP,EpsG,p1tb,mt,SME(19,1:4))
      call SpATop3(Spi_bot,k2w,EpsG,k3g,p1tb,mt,SME(21,1:4))
      call SpATop3(Spi_bot,k2w,EpsP,k3g,p1tb,mt,SME(20,1:4))

      call SpATop3(Spi_bot,EpsP,EpsG,k3g,p1tb,mt,SME(22,1:4))


      call SpATop4(Spi_bot,EpsW,EpsP,EpsG,k3g,p1tb,mt,SME(23,1:4))
      call SpATop5(Spi_bot,EpsW,k2w,EpsP,EpsG,k3g,p1tb,mt,SME(24,1:4))
      

      t1 = EL**2
      t2 = t1*EpsPk3g
      t4 = EpsWEpsG*GS*GW
      t5 = t2*t4
      t6 = mt*ne
      t7 = t6*Qqu
      t8 = 1/k3gp1tb
      t9 = two*t8
      t10 = mt**2
      t11 = Mw**2
      t12 = k1bbk2w*two
      t14 = 1/(-t10+t11+t12)
      t15 = t9*t14
      t16 = t7*t15
      t18 = t1*EpsGp1tb
      t20 = EpsWEpsP*GS*GW
      t21 = t18*t20
      t23 = t1*EpsPEpsG
      t25 = EpsWk3g*GS*GW
      t26 = t23*t25
      t28 = t1*EpsGk1bb
      t29 = t28*t20
      t30 = 1/k1bbk3g
      t31 = two*t30
      t32 = k1bbk3g*two
      t33 = k2wk3g*two
      t35 = 1/(-t10+t11+t12+t32+t33)
      t36 = t31*t35
      t40 = EpsWEpsP*four*GS
      t41 = t28*t40
      t42 = GW*mt
      t43 = t42*ne
      t45 = Qqu*t14*t35
      t46 = t43*t45
      t48 = t1*EpsGk2w
      t49 = t48*t40
      t51 = -Qqd+Qqu
      t52 = t6*t51
      t54 = k3gp1tb*two
      t56 = 1/(t10-t11+t32-k1bbp1tb*two-t54)
      t57 = t31*t56
      t60 = t9*t56
      t64 = k2wp1tb*two
      t66 = 1/(t10+t11+t33-t64-t54)
      t67 = t9*t66
      t68 = t6*Qqd*t67
      t73 = EpsPk2w*GS*GW
      t74 = t18*t73
      t75 = ne*Qqu
      t76 = t75*t15
      t80 = t48*EpsPk3g*GS*GW
      t82 = GS*GW
      t83 = t82*k2wk3g
      t84 = t23*t83
      t86 = t28*t73
      t87 = t75*t36
      t90 = EpsPk2w*four*GS
      t92 = GW*ne
      t93 = t92*t45
      t97 = ne*t51
      t98 = t97*t57
      t100 = t97*t60
      t102 = t23*t82
      t103 = ne*Qqd
      t105 = t103*two*t66
      t107 = t103*t67
      t111 = four*GS*GW
      t115 = 1/(t10+t11-t64)
      t117 = Qqd*t115*t66
      t120 = t82*t10
      t124 = t103*two*t115*t66
      t126 = t82*t11
      t129 = t74*t76-t80*t76+t84*t76-t86*t87-t28*t90*t93-t48*t90*t93-t86
     #*t98+t74*t100-t102*t105-t80*t107+t84*t107-t23*t111*k2wp1tb*ne*t117
     #+t23*t120*t124+t23*t126*t124
      t135 = t49*t93
      t139 = t26*t107
      t141 = t1*EpsWEpsG
      t142 = t141*t82
      t144 = t141*t83
      t147 = EpsWp1tb*GS*GW
      t148 = t18*t147
      t150 = t18*t25
      t152 = t48*t25
      t156 = t28*t147
      t159 = EpsWk1bb*GS*GW
      t160 = t18*t159
      t163 = EpsWp1tb*four*GS
      t165 = t92*t117
      t171 = t142*t105-t144*t107+t148*t107-t150*t107+t152*t107+t150*t100
     #-t148*t100+t156*t98+t160*t100+t48*t163*t165-t18*t163*t165-t144*t76
     #+t152*t76
      t174 = EpsWk1bb*four*GS
      t179 = t28*t159
      t181 = t28*t25
      t183 = t141*t111
      t185 = k1bbk2w*ne*t45
      t188 = EpsWk3g*four*GS
      t194 = t75*two*t14*t35
      t197 = k2wk3g*ne*t45
      t206 = t160*t76-t28*t174*t93-t48*t174*t93-t179*t87-t181*t87+t183*t
     #185-t48*t188*t93-t141*t120*t194+t183*t197+t141*t126*t194-t156*t103
     #*t31*t115-t179*t98-t181*t98
      t208 = t2*t159
      t210 = t1*EpsPk2w
      t211 = t210*t25
      t213 = t1*EpsWEpsP
      t214 = t213*t83
      t216 = t213*t82
      t222 = t213*t111
      t231 = t2*t25
      t233 = t2*t147
      t240 = -t208*t76-t211*t76+t214*t76+t216*t75*two*t35+t210*t188*t93-
     #t222*t185-t222*t197+t213*t120*t194-t213*t126*t194-t208*t100-t211*t
     #100-t231*t100+t233*t100+t214*t100+t231*t107-t233*t107+t2*t163*t165
      t241 = t210*t4
      t243 = t48*t20
      t245 = t23*t159
      t255 = t23*t147
      t260 = t241*t76-t243*t76+t245*t76-t210*EpsWEpsG*four*GS*t93+t135+t
     #241*t100-t243*t100+t245*t100+t26*t100-t255*t100-t139+t255*t107-t23
     #*t163*t165
      t261 = t2*t82
      t264 = t6*Qqu*t8*t14
      t268 = t6*Qqd*t8*t66
      t270 = t82*mt
      t274 = t18*t82
      t276 = t28*t82
      t279 = t6*Qqu*t30*t35
      t283 = t48*t270
      t299 = t1*EpsWk3g
      t300 = t299*t82
      t325 = t75*t8*t14
      t328 = t75*t30*t35
      t331 = t48*t82
      t332 = t331*t194
      t334 = t103*t30*t115
      t337 = t103*t8*t66
      t350 = t1*GS
      t352 = t350*GW*k2wk3g
      t354 = t350*GW
      t358 = t350*GW*t10
      t360 = t75*t14*t35
      t363 = t350*GW*t11
      t373 = t103*t115*t66
      t379 = -t352*t325-t354*t75*t35-t358*t360+t363*t360+t350*GW*k1bbk2w
     #*t194+t352*t194+t354*t103*t66-t352*t337-t358*t373-t363*t373+t350*G
     #W*k2wp1tb*t124
      t383 = t210*t82
      t388 = t97*t30*t56
      t391 = t97*t8*t56
      t396 = t300*t337
      t409 = t1*EpsWk1bb*t82
      t415 = t1*EpsWp1tb*t82
      t425 = -t409*t325+t409*t328+t300*t328+t409*t194+t415*t334+t409*t38
     #8+t300*t388-t415*t388-t409*t391-t300*t391+t415*t391+t396-t415*t337
     #+t415*t124
      t426 = t350*t43
      t427 = one*Qqu
      t428 = 1/two
      t429 = t8*t428
      t431 = t427*t429*t14
      t433 = t30*t428
      t435 = t427*t433*t35
      t437 = t350*t42
      t439 = one*Qqd
      t441 = t439*t433*t115
      t444 = t439*t429*t66
      t448 = t350*t92
      SMEcoeff(1) = t5*t16-t21*t16-t26*t16+t29*t7*t36+t41*t46+t49*t46+t2
     #9*t52*t57-t21*t52*t60+t5*t68-t26*t68
      SMEcoeff(2) = t129
      SMEcoeff(3) = t5*t76-t21*t76-t26*t76+t29*t87+t41*t93+t135+t29*t98-
     #t21*t100+t5*t107-t139
      SMEcoeff(4) = t171+t206
      SMEcoeff(5) = t240
      SMEcoeff(6) = t260
      SMEcoeff(7) = -t261*t264-t261*t268+t2*t270*t124
      SMEcoeff(8) = t274*t264-t276*t279-t28*t270*t194-t283*t194-t276*t6*
     #Qqd*t30*t115+t274*t268+t283*t124-t18*t270*t124
      SMEcoeff(9) = t102*t264+t102*t268-t23*t270*t124
      SMEcoeff(10) = t300*t264-t299*t270*t194+t300*t268
      SMEcoeff(11) = -t142*t264+t141*t270*t194-t142*t268
      SMEcoeff(12) = t216*t264-t216*t279-t213*t270*t194-t216*t6*t51*t30*
     #t56+t216*t6*t51*t8*t56
      SMEcoeff(13) = -t274*t325+t276*t328+t276*t194+t332+t276*t334-t274*
     #t337-t331*t124+t274*t124
      SMEcoeff(14) = t261*t325+t261*t337-t261*t124
      SMEcoeff(15) = -t102*t325-t102*t337+t102*t124
      SMEcoeff(16) = t379
      SMEcoeff(17) = t331*t325-t332+t331*t337
      SMEcoeff(18) = -t383*t325+t383*t328+t383*t194+t383*t388-t383*t391
      SMEcoeff(19) = t300*t325-t300*t194+t396
      SMEcoeff(20) = -t142*t325+t142*t194-t142*t337
      SMEcoeff(21) = t216*t325-t216*t328-t216*t194-t216*t388+t216*t391
      SMEcoeff(22) = t425
      SMEcoeff(23) = -t426*t431+t426*t435+t437*t360+t426*t441-t426*t444+
     #t437*t373
      SMEcoeff(24) = t448*t431-t448*t435-t354*t360-t448*t441+t448*t444-t
     #354*t373

      do i =1,4
         res(i) = zero
      enddo	
      do i=1,24
         res(1:4) = res(1:4) + SMEcoeff(i)*SME(i,1:4)
      enddo


      return
      end


      subroutine ATopDipole(gauge,heli,k1bb,klep,kneu,k3g,k4p,
     &alphadip,k1bb_tilde,k4p_tilde,klep_tilde,kneu_tilde,res,lo)
      implicit none
      include 'realcommon.f'
      integer heli(4), gauge, DV
      double complex p1tb(4),k2w(4),k1bb(4),k3g(4), k4p(4)
      double complex k2w_dip(4), k2w_tilde(4), kneu(4), klep(4)
      double complex k1bb_tilde(4), p1tb_tilde(4), k4p_tilde(4)
      double complex kneu_tilde(4), klep_tilde(4)
      double complex EpsP(4),EpsW(4), lo(4)
      double complex res, scf, kdum(4), kdum1(4)
      double complex k1bbk3g,k3gp1tb, k1bbp1tb
      double precision condition1, condition2, r
      double precision alphadip, ptpho, ang, rat
      double precision GS, y,z, Dipole
      logical logic
      Dv = 4
      gs = dsqrt(4.d0*pii*alphaS)
      k2w(1:4) = kneu(1:4) + klep(1:4)
      k2w_dip(1:4) = k2w(1:4) + k4p(1:4)
      p1tb(1:4) = k2w(1:4) + k1bb(1:4) + k4p(1:4) + k3g(1:4)

      lo(1:4) = zero
      k1bb_tilde(1:4) = zero
      k4p_tilde(1:4) = zero
      klep_tilde(1:4) = zero
      kneu_tilde(1:4) = zero



      k1bbp1tb = scf(Dv,k1bb,p1tb)
      k1bbk3g = scf(Dv,k1bb,k3g)
      k3gp1tb = scf(Dv,p1tb,k3g)



C This is the correct implmentation of the dipole
C Don't follow the paper of Campbell and Ellis
      r = dsqrt(dreal(scf(Dv,k2w_dip,k2w_dip)/mt**2))
      y = 2.d0*dreal(k1bbk3g)/(mt**2*(1.d0-r)**2)
C This is the correct implmentation of the dipole
C Don't follow the paper of Campbell and Ellis
      z = 1.d0 - 2.d0*dreal(k3gp1tb)/mt**2/(1.d0-r**2)




      condition1 = 1.d0-alphadip-z
      condition2 = y-alphadip*(1.d0+r)**2*z*(1.d0-z)/
     &(z+r**2*(1.d0-z))



      call DipoleBoost(k2w_dip,p1tb,kneu,kneu_tilde)
      call DipoleBoost(k2w_dip,p1tb,klep,klep_tilde)
      call DipoleBoost(k2w_dip,p1tb,k4p,k4p_tilde)

      k1bb_tilde(1:4) = p1tb(1:4) - k4p_tilde(1:4)
     & - kneu_tilde(1:4) - klep_tilde(1:4)

      call cut2(k1bb_tilde,k4p_tilde,logic,rat)
      call LOATilde(gauge,heli,k1bb_tilde
     &     ,k4p_tilde,klep_tilde,kneu_tilde,lo)

      Dipole =0.d0
      If(logic .eq. .true.) then
         If (condition1 .lt. 0.d0 .or. condition2 .lt. 0.d0) then
            Dipole = GS**2*CF*(1.d0/dreal(k1bbk3g)*
     &           (2.d0/(1.d0-z)-1.d0-z) - mt**2/dreal(k3gp1tb)**2)
         endif
      endif
      res = dcmplx(Dipole,0.d0)	

      end      

      subroutine LOATilde(gauge,heli,k1bb,k4p,klep,kneu,lo)
      implicit none
      include 'VarALO.f'
      include 'realcommon.f'
      integer heli(4), heli_pho, i, DV, gauge, heli_w, heli_b
      double complex k1bbp1tb,k1bbk2w,k2wp1tb
      double complex EpsWk1bb, EpsWp1tb, EpsWEpsP, EpsP(4)
      double complex EpsPk1bb, EpsPp1tb, EpsPk2w , EpsW(4), scf
      double complex p1tb(4),k2w(4),k1bb(4), k4p(4), klep(4), kneu(4)
      double complex lo(4), Spi_bot(4), SMEcoeff(6), SME(6,4)
      double complex one, two, four, three

      Dv =4

      one = dcmplx(1.d0,0.d0)
      two = dcmplx(2.d0,0.d0)
      four = dcmplx(4.d0,0.d0)
      three = dcmplx(3.d0,0.d0)
      
      heli_b =+1
      heli_pho = heli(3)
      heli_W = heli(4)
      k2w(1:4) = kneu(1:4) + klep(1:4)
      p1tb(1:4) = k2w(1:4) + k1bb(1:4) + k4p(1:4)


      call vSpi_Weyl(k1bb,heli_b,Spi_bot)
      call pol_massoff2(klep,kneu,EpsW)
C      call pol_mass(k2w,Mw,heli_w,EpsW)
      call pol_mless(k4p,heli_pho,EpsP)

      If (gauge .eq. 1) then
         EpsP(1:4) = k4p(1:4)
      endif


      k1bbp1tb = scf(Dv,k1bb,p1tb)
      k1bbk2w = scf(Dv,k1bb,k2w)
      k2wp1tb = scf(Dv,k2w,p1tb)
      EpsWk1bb = scf(Dv,EpsW,k1bb)
      EpsWp1tb = scf(Dv,EpsW,p1tb)
      EpsWEpsP = scf(Dv,EpsW,EpsP)
      EpsPk1bb = scf(Dv,EpsP,k1bb)
      EpsPp1tb = scf(Dv,EpsP,p1tb)
      EpsPk2w = scf(Dv,EpsP,k2w)


      call SpATop0(Spi_bot,p1tb,mt,SME(1,1:4))
      call SpATop1(Spi_bot,EpsW,p1tb,mt,SME(2,1:4))
      call SpATop1(Spi_bot,EpsP,p1tb,mt,SME(4,1:4))
      call SpATop1(Spi_bot,k2w,p1tb,mt,SME(3,1:4))

      call SpATop2(Spi_bot,EpsW,EpsP,p1tb,mt,SME(5,1:4))
      call SpATop3(Spi_bot,EpsW,k2w,EpsP,p1tb,mt,SME(6,1:4))

      t1 = EL**2
      t2 = t1*EpsWEpsP
      t6 = mt**2
      t8 = Mw**2
      t16 = -t6*Qqd+t8*Qqd+k1bbk2w*Qqd*two-k1bbk2w*Qqu*two+k1bbp1tb*Qqu*
     #two
      t19 = 1/(-t6+t8+k1bbk2w*two)
      t23 = 1/(t6-t8-k1bbp1tb*two)
      t28 = GW*ne
      t31 = t19*t23
      t32 = two*t16*t31
      t36 = t1*GW
      t49 = 1/(k1bbk2w-k1bbp1tb)
      t50 = t16*t49
      t56 = 1/two
      SMEcoeff(1) = t2*GW*mt*ne*two*t16*t19*t23
      SMEcoeff(2) = -t1*EpsPk2w*t28*t32
      SMEcoeff(3) = t2*t28*t32
      SMEcoeff(4) = t36*ne*(-EpsWp1tb*t6+EpsWp1tb*t8-EpsWk1bb*k1bbk2w*tw
     #o+EpsWp1tb*k1bbk2w*two+EpsWk1bb*k1bbp1tb*two)*t50*t31
      SMEcoeff(5) = -t36*mt*ne*one*t16*t49*t56*t19
      SMEcoeff(6) = t36*ne*one*t50*t56*t19

      do i=1,4
        lo(i) = zero
      enddo
      do i=1,6
         lo(1:4) = lo(1:4) + SMEcoeff(i)*SME(i,1:4)
      enddo
      return
      end      

       subroutine SpATop0(sp1,slash,mt,sp2)
       implicit none
       double complex sp2(4),sp1(4),slash(4),sp1n(4)
       double complex pslash(4)
       double precision mt
       integer Dv, Ds
       Dv = 4
       Ds = 4
       pslash(1:4) = -slash(1:4)
       call spi2_weyl(Dv,Ds,sp1,pslash,sp1n)
       sp2(1:4) = sp1n(1:4) + mt*sp1(1:4)

       return 
       end

       subroutine SpATop1(sp1,ap1,slash,mt,sp2)
       implicit none
       double complex sp2(4),sp1(4), ap1(4),slash(4)
       double complex pslash(4)
       double complex sp1n(4), sp2n(4)
       double precision mt
       integer Dv, Ds
       Dv = 4
       Ds = 4
       pslash(1:4) = -slash(1:4)
       call spi2_weyl(Dv,Ds,sp1,ap1,sp1n)
       call spi2_weyl(Dv,Ds,sp1n,pslash,sp2n)
       sp2(1:4) = sp2n(1:4) + mt*sp1n(1:4)
       return 
       end


       subroutine SpATop2(sp1,ap2,ap1,slash,mt,sp2) 
       implicit none
       double complex sp2(4),sp1(4), ap1(4), slash(4)
       double complex sp1n(4),ap2(4), sp2n(4)
       double complex pslash(4)
       double precision mt
       integer Dv, Ds
       Dv = 4
       Ds = 4
       pslash(1:4) = -slash(1:4)		
       call spi2_weyl(Dv,Ds,sp1,ap1,sp1n)
       call spi2_weyl(Dv,Ds,sp1n,ap2,sp2n)
       call spi2_weyl(Dv,Ds,sp2n,pslash,sp1n)
       sp2(1:4) = sp1n(1:4) + mt*sp2n(1:4)
       return 
       end


       subroutine SpATop3(sp1,ap3,ap2,ap1,slash,mt,sp2)
       implicit none
       double complex sp2(4),sp1(4), ap1(4), slash(4)
       double complex sp1n(4), sp2n(4), ap2(4), ap3(4)
       double complex pslash(4)
       double precision mt
       integer Dv, Ds
       Dv = 4
       Ds = 4
       pslash(1:4) = -slash(1:4)		

       call spi2_weyl(Dv,Ds,sp1,ap1,sp1n)
       call spi2_weyl(Dv,Ds,sp1n,ap2,sp2n)
       call spi2_weyl(Dv,Ds,sp2n,ap3,sp1n)
       call spi2_weyl(Dv,Ds,sp1n,pslash,sp2n)
       sp2(1:4) = sp2n(1:4) + mt*sp1n(1:4)
       return 
       end

       subroutine SpATop4(sp1,ap4,ap3,ap2,ap1,slash,mt,sp2)
       implicit none
       double complex sp2(4),sp1(4), ap1(4), slash(4)
       double complex sp1n(4), sp2n(4), ap2(4), ap3(4), ap4(4)
       double complex pslash(4)
       double precision mt
       integer Dv, Ds
       Dv = 4
       Ds = 4
       pslash(1:4) = -slash(1:4)		

       call spi2_weyl(Dv,Ds,sp1,ap1,sp1n)
       call spi2_weyl(Dv,Ds,sp1n,ap2,sp2n)
       call spi2_weyl(Dv,Ds,sp2n,ap3,sp1n)
       call spi2_weyl(Dv,Ds,sp1n,ap4,sp2n)
       call spi2_weyl(Dv,Ds,sp2n,pslash,sp1n)
       sp2(1:4) = sp1n(1:4) + mt*sp2n(1:4)

       return 
       end

       subroutine SpATop5(sp1,ap5,ap4,ap3,ap2,ap1,slash,mt,sp2)
       implicit none
       double complex sp2(4),sp1(4), ap1(4), ap4(4),slash(4)
       double complex sp1n(4), sp2n(4), ap2(4), ap3(4), ap5(5)
       double complex pslash(4)
       double precision mt
       integer Dv, Ds
       Dv = 4
       Ds = 4
       pslash(1:4) = -slash(1:4)		

       call spi2_weyl(Dv,Ds,sp1,ap1,sp1n)
       call spi2_weyl(Dv,Ds,sp1n,ap2,sp2n)
       call spi2_weyl(Dv,Ds,sp2n,ap3,sp1n)
       call spi2_weyl(Dv,Ds,sp1n,ap4,sp2n)
       call spi2_weyl(Dv,Ds,sp2n,ap5,sp1n)
       call spi2_weyl(Dv,Ds,sp1n,pslash,sp2n)
       sp2(1:4) = sp2n(1:4) + mt*sp1n(1:4)

       return 
       end





      subroutine vbar(p,m,i,f)
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
!     fc=dsqrt(fc2)
      
      
      if (i.eq.1) then
         f(1)=pz*fc/fc2
         f(2)=(px-(0d0,1d0)*py)*fc/fc2
         f(3)=-fc
         f(4)=dcmplx(0d0,0d0)
      endif
      
      if (i.eq.-1) then
         f(1)=(px+(0d0,1d0)*py)*fc/fc2
         f(2)=-pz*fc/fc2
         f(3)=dcmplx(0d0,0d0)
         f(4)=-fc
      endif
      return
      end 



!-------massive vector boson polarization routine, for off-shell vector bosons decaying into Neutrino and Lepton
      subroutine pol_massoff2(p_lep,p_neu,EpsW_weyl)
      implicit none
      double complex p_Lep(4), p_neu(4)
      double complex pl(4),pn(4),Wpol(4)
      double complex s_n(4), s_l(4)
      double complex  Weyl(4), EpsW_weyl(4), scf
      integer Dv,Ds,i, hl, hn
      real(8), parameter :: sqrt2 = 1.41421356237309504880168872420969d0
      Dv = 4
      Ds = 4
      do i=1,4
         EpsW_weyl(i) = (0.d0,0.d0)
         Weyl(i) = (0.d0,0.d0)
      enddo
      call ubarSpi_Weyl(p_lep,-1,s_l)
      call vSpi_Weyl(p_neu,+1,s_n)
      call vbqqWeyl(Dv,Ds,s_l,s_n,Weyl)
      EpsW_weyl(1:4) =  Weyl(1:4)*(0.d0,1.d0)
      end


