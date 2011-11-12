      MODULE ModTTBP_NLODK
      use ModParameters
      use ModMisc
      implicit none

        public  :: EvalTTBP_1LDK
        public  :: EvalTTBP_REALDK, EvalTTBP_DIPOLDK

        private :: bernuo

      contains



      SUBROUTINE EvalTTBP_DIPOLDK(Top_Atop,Momenta,MomentaTilde,dipole)
      implicit none
      include 'commondecay.f'
      complex(8) :: SpiOut(1:4)
      complex(8) :: k1b(4), kneu(4), klep(4), k4p(4), p1t(4), k3g(4)
      complex(8) :: p1t_t(4),k2w_t(4),k1b_t(4),k3g_t(4),k4p_t(4)
      complex(8) :: kneu_t(4), klep_t(4)
      real(8) :: res, dipole
      real(8) :: Momenta(1:4,1:5),  g2weak, MomentaTilde(1:4,1:5), alphadip
      real(8) :: mu, runfactor
      ! heli(1) = Bottom-Heli, heli(2) = Photon-Heli, heli(3) = Top-Heli (Top-Heli not needed)
      integer :: Top_Atop
      integer :: i, gau, Dv
      ! set parameters for Scharf's LO and NLO (Anti-)Top Decay with photon
      Dv =4
      mt = m_top
      Mw = m_w
      Qqu = Q_top
      Qqd = Q_dn
      runfactor = RunAlphaS(NLOParam,MuRen)
      alphas = alpha_S*runfactor
      GS = dsqrt(alphas*4.d0*pii)
      EL = dsqrt(4.d0*pii*alpha)
      g2weak = 4d0*dsqrt(2d0)*m_W**2*GF
      GW = 1.d0/dsqrt(2.d0)*dsqrt(g2weak)/EL
      CF = 4.d0/3.d0
      nc = 3.d0
      ! my ren-scale, result has to be independent of mu
      mu = muRen
      alphadip = alpha_DKTfi


!      if(Top_ATop .eq. 1 .and. heli(1) .eq. 1) then
!         write(*,*) "Wrong setting for Bottom-Helicity"
!         write(*,*) "Particle is top and helicity for bottom
!     &is incorrect"
!         write(*,*) "See subroutine EvalTTBP_REALDK"
!         stop
!      endif
!      if(Top_ATop .eq. -1 .and. heli(1) .eq. -1) then
!         write(*,*) "Wrong setting for Anti-Bottom-Helicity"
!         write(*,*) "Particle is anti-top and helicity for anti-bottom
!     &    is incorrect"
!         write(*,*) "See subroutine EvalTTBP_REALDK"
!         stop
!      endif


      ! Assign Momenta
      do i=1,4
       k1b(i) = dcmplx(Momenta(i,1),0.d0)
       klep(i) = dcmplx(Momenta(i,2),0.d0)
       kneu(i) = dcmplx(Momenta(i,3),0.d0)
       k3g(i) = dcmplx(Momenta(i,4),0.d0)
       k4p(i) = dcmplx(Momenta(i,5),0.d0)
       p1t(i) = k1b(i) + klep(i) + kneu(i) + k4p(i) + k3g(i)
      enddo
!      write(*,*) "Dipol -mt**2", scf(Dv,p1t,p1t)
!      write(*,*) "Dipol - k1b.k3g", scf(Dv,k1b,k3g)
!      write(*,*) "Dipol - Eg/Etop", k3g(1)/p1t(1)
      call Calc_DIPOL(Top_Atop,alphadip,k1b,klep,kneu,k3g,k4p,
     & res,k1b_t,klep_t,kneu_t,k4p_t)



      dipole = res
      do i=1,4
         MomentaTilde(i,1) = dreal(k1b_t(i))
         MomentaTilde(i,2) = dreal(klep_t(i))
         MomentaTilde(i,3) = dreal(kneu_t(i))
         MomentaTilde(i,4) = 0.d0
         MomentaTilde(i,5) = dreal(k4p_t(i))
      enddo
      RETURN
      END SUBROUTINE EvalTTBP_DIPOLDK



      subroutine Calc_DIPOL(Top_Atop,alphadip,k1b,klep,kneu,k3g,k4p,
     &   res,k1b_t,klep_t,kneu_t,k4p_t)
      implicit none
      include "commondecay.f"
      complex(8) :: p1t(4),k2w(4),k1b(4),k3g(4),k4p(4)
      complex(8) :: kneu(4), klep(4), kdum(4)
      complex(8) :: reals(4), dirac_real(4)
      complex(8) :: p1t_t(4),k2w_t(4),k1b_t(4),k3g_t(4),k4p_t(4)
      complex(8) :: kneu_t(4), klep_t(4)

      double precision alphadip, res
      integer Top_Atop, gau

      if(Top_Atop .eq. 1) then
         call TopDipole(k1b,klep,kneu,k3g,k4p,alphadip,
     &    k1b_t,k4p_t,klep_t,kneu_t,res)
      endif
      if( Top_Atop .eq. -1) then
         call ATopDipole(k1b,klep,kneu,k3g,k4p,alphadip,
     &    k1b_t,k4p_t,klep_t,kneu_t,res)
      endif
      end subroutine Calc_DIPOL





      SUBROUTINE EvalTTBP_REALDK(Top_Atop,heli,Momenta,gau,Wpar,SpiOut)
      implicit none
      include 'commondecay.f'
      complex(8) :: SpiOut(1:4)
      complex(8) :: k1b(4), kneu(4), klep(4), k4p(4), p1t(4), k3g(4)
      complex(8) :: res(4), k2w(4)
      real(8) :: Momenta(1:4,1:5),  g2weak
      real(8) :: mu, runfactor
      ! heli(1) = Bottom-Heli, heli(2) = Photon-Heli, heli(3) = Top-Heli (Top-Heli not needed)
      integer :: heli(1:3), Top_Atop
      integer :: i, gau, DV, Wpar(1:4)
      ! set parameters for Scharf's LO and NLO (Anti-)Top Decay with photon
      Dv = 4
      mt = m_top
      Mw = m_w
      Qqu = Q_top
      Qqd = Q_dn
      runfactor = RunAlphaS(NloParam,MuRen)
      alphas = alpha_S*runfactor
      GS = dsqrt(alphas*4.d0*pii)
      EL = dsqrt(4.d0*pii*alpha)
      g2weak = 4d0*dsqrt(2d0)*m_W**2*GF
      GW = 1.d0/dsqrt(2.d0)*dsqrt(g2weak)/EL
      CF = 4.d0/3.d0
      nc = 3.d0
      ! my ren-scale, result has to be independent of mu
      mu = muRen


      if(Top_ATop .eq. 1 .and. heli(1) .eq. 1) then
         write(*,*) "Wrong setting for Bottom-Helicity"
         write(*,*) "Particle is top and helicity for bottom
     &is incorrect"
         write(*,*) "See subroutine EvalTTBP_REALDK"
         stop
      endif
      if(Top_ATop .eq. -1 .and. heli(1) .eq. -1) then
         write(*,*) "Wrong setting for Anti-Bottom-Helicity"
         write(*,*) "Particle is anti-top and helicity for anti-bottom
     &    is incorrect"
         write(*,*) "See subroutine EvalTTBP_REALDK"
         stop
      endif


      ! Assign Momenta
      do i=1,4
       k1b(i) = dcmplx(Momenta(i,1),0.d0)
       klep(i) = dcmplx(Momenta(i,2),0.d0)
       kneu(i) = dcmplx(Momenta(i,3),0.d0)
       k3g(i) = dcmplx(Momenta(i,4),0.d0)
       k4p(i) = dcmplx(Momenta(i,5),0.d0)
       k2w(i) = kneu(i) + klep(i)
       p1t(i) = k1b(i) + klep(i) + kneu(i) + k4p(i) + k3g(i)
      enddo
!      write(*,*) "Real -mt**2", scf(Dv,p1t,p1t)
!      write(*,*) "Real - k1b.k3g", scf(Dv,k1b,k3g)
!      write(*,*) "Real - Eg/Etop", k3g(1)/p1t(1)

!      write(*,*) "Real p1t.p1t/mt**2", scf(Dv,p1t,p1t)/mt**2
!      write(*,*) "Real k2w.k2w/mw**2", scf(Dv,k2w,k2w)/Mw**2
!      write(*,*) "Real k3g.k3g/mt**2", scf(Dv,k3g,k3g)/mt**2
!      write(*,*) "Real k4p.k4p/mt**2", scf(Dv,k4p,k4p)/Mt**2
!      write(*,*) "Real k1b.k1b/mt**2", scf(Dv,k1b,k1b)/Mt**2
!      write(*,*) "p1t -in", p1t
      call Calc_REAL(Top_Atop,heli,Wpar,k1b,klep,kneu,k3g,k4p,gau,res)



      SpiOut(1:4) = (res(1:4))*dsqrt(CF)

      RETURN
      END SUBROUTINE EvalTTBP_REALDK



      subroutine Calc_REAL(Top_Atop,heli,Wpar,k1b,klep,kneu,k3g,k4p,gau,res)
      implicit none
      include "commondecay.f"
      double complex p1t(4),k2w(4),k1b(4),k3g(4),k4p(4)
      double complex kneu(4), klep(4), kdum(4)
      double complex reals(4), dirac_real(4), res(4), fac
      integer Top_Atop, heli(1:3), gau, Wpar(1:4)
      if(Top_Atop .eq. 1) then
         call numDecTopReals(heli,Wpar,k1b,klep,kneu,k3g,k4p,gau,reals)
         call WeyltoDiracS(reals,dirac_real)
      endif


      if( Top_Atop .eq. -1) then
         call numDecATopReals(heli,Wpar,k1b,klep,kneu,k3g,k4p,gau,reals)
         call WeyltoDiracS(reals,dirac_real)
      endif
      ! factor ne = (0.d0,1.d0) is the phase from NWA
      res = Top_Atop*ne*dirac_real

      end subroutine Calc_REAL



      SUBROUTINE EvalTTBP_1LDK(LO_NLO,Top_Atop,xe,gau,heli,Wpar,Momenta,SpiOut)
      implicit none
      include 'commondecay.f'
      complex(8) :: SpiOut(1:4)
      complex(8) :: k1b(4), kneu(4), klep(4), k3p(4), p1t(4), k2w(4)
      complex(8) :: alphadip, eta
      complex(8) :: res(4), lo(4), res_iop(-2:0)
      real(8) :: Momenta(1:4,1:4), g2weak
      real(8) :: mu, runfactor, gram_det_cut
      ! heli(1) = Bottom-Heli, heli(2) = Photon-Heli, heli(3) = Top-Heli (Top-Heli not needed)
      integer heli(1:2), Top_Atop, xe, LO_NLO, gau
      integer i, Dv, Wpar(1:4)
      ! set parameters for Scharf's LO and NLO (Anti-)Top Decay with photon
      call qlinit
      call ffini
      Dv = 4
      mt = m_top
      Mw = m_w
      Qqu = Q_top
      Qqd = Q_dn
      runfactor = RunAlphaS(NLOParam,MuRen)
      alphas = alpha_S*runfactor
      EL = dsqrt(4.d0*pii*alpha)
      g2weak = 4d0*dsqrt(2d0)*m_W**2*GF
      GW = 1.d0/dsqrt(2.d0)*dsqrt(g2weak)/EL
      CF = 4.d0/3.d0
      nc = 3.d0
C      NWAFactor_W   = 1d0/dsqrt(2d0*Ga_W*m_W)
C      NWAFactor_Top = 1d0/dsqrt(2d0*Ga_Top(0)*m_Top)

      ! left over from Scharf's debugging, eta value for calc int. Dipole in FDH or HV scheme
      eta = (1.d0,0.d0)
      ! my Dipole alpha
      alphadip = dcmplx(alpha_DKTfi,0.d0)
      ! my ren-scale, result hsa to be independent of mu
      mu = muRen

      if(Top_ATop .eq. 1 .and. heli(1) .eq. 1) then
         write(*,*) "Wrong setting for Bottom-Helicity"
         write(*,*) "Particle is top and helicity for bottom
     &is incorrect"
         write(*,*) "See subroutine EvalTTBP_1LDK"
         stop
      endif
      if(Top_ATop .eq. -1 .and. heli(1) .eq. -1) then
         write(*,*) "Wrong setting for Anti-Bottom-Helicity"
         write(*,*) "Particle is anti-top and helicity for anti-bottom
     &    is incorrect"
         write(*,*) "See subroutine EvalTTBP_1LDK"
         stop
      endif


      ! Assign Momenta & Initialize
      do i=1,4
       k1b(i) = dcmplx(Momenta(i,1),0.d0)
       klep(i) = dcmplx(Momenta(i,2),0.d0)
       kneu(i) = dcmplx(Momenta(i,3),0.d0)
       k3p(i) = dcmplx(Momenta(i,4),0.d0)
       k2w(i) =  klep(i) + kneu(i)
       p1t(i) = k1b(i) + klep(i) + kneu(i) + k3p(i)
       res(i) = zero
       lo(i) = zero
      enddo

      gram_det_cut = dabs(dreal((2.d0* scf(Dv,k1b,k2w)* scf(Dv,k1b,p1t)
     & - scf(Dv,k1b,k2w)*mt**2 + scf(Dv,k1b,p1t)*Mw**2)))/mt**4

      if (gram_det_cut .gt. 1.d-4) then
         call Calc_ALL(LO_NLO,Top_Atop,xe,gau,mu,alphadip,eta,heli,Wpar
     &                 ,p1t,k1b,klep,kneu,k3p,res,lo,res_iop)
      endif


      SpiOut(1:4) = Top_Atop*(res(1:4) + lo(1:4)*res_iop(xe))
      RETURN
      END SUBROUTINE







      ! SCHARF's virtual call: Set_kinemtics, call (A)-Top DK
      subroutine Calc_ALL(LO_NLO,Top_Atop,xe,gau,mu,alphadip,eta,heli,Wpar,
     &p1t,k1b,klep,kneu,k3p,res,lo,res_iop)
        implicit none
        include 'commondecay.f'
        double complex p1t(4),k1b(4),klep(4),kneu(4),k3p(4)
        double complex resv(4), resr(4)
        double complex resctzT(4), resctm(4)
        double precision mu
        double complex alphadip, eta, res(4), res_iop(-2:0), lo(4)
        integer i,xe, heli(1:2), LO_NLO,Top_Atop, gau, Wpar(1:4)
        call set_kinematics(Top_Atop,gau,heli,Wpar,p1t,k1b,klep,kneu,k3p)
        do i=1,4
           resv(i) = zero
           resr(i) = zero
           resctm(i) = zero
           resctzT(i) =zero
           res(i) =zero
           lo(i) =zero
        enddo
        call TopGammadecays(LO_NLO,Top_Atop,xe,mu,alphadip,eta,lo,
     &resv,resr,resctm,resctzT,res_Iop)
        if (LO_NLO .eq. 2) then
      ! factor ne = (0.d0,1.d0) is the phase from NWA
           res(1:4) = 2.d0*ne*
     &      (resv(1:4) + resr(1:4) + resctm(1:4) + resctzT(1:4))*CF
           lo(1:4) = lo(1:4)*ne*CF
        endif
        if (LO_NLO .eq. 1) then
      ! factor ne = (0.d0,1.d0) is the phase from NWA
           res(1:4) = ne*lo(1:4)
           res_iop(-2:0) = zero
C           write(*,*) "Hallo res", res(1:4)
        endif
        if (LO_NLO .ne. 1 .and. LO_NLO .ne. 2) then
           write(*,*) "The lo/virtual corrections for Photon Decay are
     &      not set correct"
           write(*,*) "Look into 'subroutine Calc_All '"
        endif
        return
      end subroutine Calc_ALL


      ! set Spinors and Polarisation vectors and load them to commonblock
      subroutine set_kinematics(Top_Atop,gau,heli,Wpar,p1t,k1b,klep,kneu,k3p)
        implicit none
        include 'commondecay.f'
        double complex k1b(4),k2w(4),k3p(4),klep(4),kneu(4), p1t(4)
        double complex e1(4),e2(4),e3(4),e4(4)
        double complex reslo
        double precision  MomPol(1:4,1:2)
        integer heli(1:2), hel_bot, hel_pho, i, Top_Atop, gau, Wpar(1:4)
        hel_bot = heli(1)
        hel_pho = heli(2)



        mom_dec(1,1:4) = k1b(1:4)
        mom_dec(2,1:4) = klep(1:4) + kneu(1:4)
        mom_dec(3,1:4) = k3p(1:4)
        mom_t(1:4) = p1t(1:4)
        do i=1,4
           MomPol(i,1) = dreal(klep(i))
           MomPol(i,2) = dreal(kneu(i))
        enddo
      ! SET kinematics for t -> b W+ gamma
        If (Top_Atop .eq. 1) then
           call ubarSpi_Weyl(k1b,hel_bot,e1)
!           call pol_massoff(klep,kneu,e2)
           call WpolvecS(Wpar(1),Wpar(2),MomPol,Wpar(3),Wpar(4),e2)
           call pol_mless(k3p,hel_pho,e3)
        endif

      ! SET kinematics for tbar -> bbar W- gamma
        If (Top_Atop .eq. -1) then
           call vSpi_Weyl(k1b,hel_bot,e1)
!           call pol_massoff2(klep,kneu,e2)
           call WpolvecS(Wpar(1),Wpar(2),MomPol,Wpar(3),Wpar(4),e2)
           call pol_mless(k3p,hel_pho,e3)
        endif

        if(gau .eq. 1) then
           e3(1:4) = k3p(1:4)
        endif

        do i=1,4
           hel_dec(1,i) = e1(i)
           hel_dec(2,i) = e2(i)
           hel_dec(3,i) = e3(i)
           hel_t(i) = e4(i)
        enddo
      end subroutine set_kinematics





      subroutine TopGammadecays(LO_NLO,Top_Atop,xe,mu,alphadip,eta,
     &res_lo,resv,resr,resctm,resctzT,Iop)
        implicit none
        include 'commondecay.f'
        double complex res_lo(4), reszw(4)
        double complex vres(4), rres(4), ctzTres(4), ctmres(4)
        double complex resv(4), resr(4)
        double complex resctzT(4), resctm(4),resIOP(4)
        double complex alphadip, eta, iop(-2:0)
        double precision mu
        integer xe,LO_NLO,Top_Atop, dummy
        dummy =0
        if (Top_Atop .eq. 1) then
           call numTopDectree(reszw)
           call WeylToDiracS(reszw,res_lo)
           call iop_me(mu,alphadip,eta,iop)
           if (LO_NLO .eq. 2) then
!CalcINTSME calcs integrals as function of xe(=eps) and mu.
! Last parameter dummy =1, then integrals are printed
              call CalcINTSME(xe,mu,dummy)
              call CalcNLO(vres,rres,ctmres,ctzTres,xe)
              call WeylToDiracS(vres,resv)
              call WeylToDiracS(rres,resr)
              call WeylToDiracS(ctmres,resctm)
              call WeylToDiracS(ctzTres,resctzT)
           endif
        endif


        if (Top_Atop .eq. -1) then
           call numATopDectree(reszw)
           call WeylToDiracS(reszw,res_lo)
           call iop_me(mu,alphadip,eta,iop)
           if (LO_NLO .eq. 2) then
              call CalcINTSMEA(xe,mu,dummy)
              call CalcANLO(vres,rres,ctmres,ctzTres,xe)
              call WeylToDiracS(vres,resv)
              call WeylToDiracS(rres,resr)
              call WeylToDiracS(ctmres,resctm)
              call WeylToDiracS(ctzTres,resctzT)
           endif
        endif
      end subroutine TopGammadecays







      subroutine Iop_me(mu,alphadip,eta,iop)
      implicit none
      include 'commondecay.f'
      double complex alphadip, eta, rli2
      double complex li2, k2w_dip(4)
      double complex rq , mtq, myln, muq, r, alphaterm
      double complex Iop(-2:0), one, two, longlog
      double precision mu
      integer i, Dv
      Dv = 4
      do i = 1,4
         k2w_dip(i) = mom_dec(2,i) +  mom_dec(3,i)
      enddo
      iop(-2:0) = (/zero,zero,zero/)

      one = dcmplx(1.d0,0.d0)
      two = dcmplx(2.d0,0.d0)
      mtq = dcmplx(mt**2,0.d0)
      muq = dcmplx(mu**2,0.d0)
      r = cdsqrt(scf(Dv,k2w_dip,k2w_dip)/mtq)
      rq = r**2
      li2 = dcmplx(1.d0,0.d0)-rq
      rli2 = dilog(li2)
      iop(-2) = dcmplx(alphaS/2.d0/pii,0.d0)
      iop(-1) = dcmplx(alphaS/2.d0/pii,0.d0)*
     &(dcmplx(5.d0/2.d0,0.d0) + cdlog(muq/mtq)
     &-two*cdlog(one-rq))
      longlog = cdlog(rq/(one-alphadip*(one-rq)))
      alphaterm = -(two*cdlog(alphadip)**2
     & - cdlog(alphadip)*(-dcmplx(7.d0/2.d0,0.d0)
     &                + dcmplx(4.d0,0.d0)*alphadip -alphadip**2/two)
     & -two*(one-alphadip)*rq/(one-rq)*longlog)
C      write(*,*) "rq, alphadip", rq,alphadip
C      write(*,*) alphaterm


      iop(0) = dcmplx(alphaS/2.d0/pii,0.d0)*
     &(dcmplx(25.d0/4.d0,0.d0)
     &+one/two*(one/(one-rq)**2 -dcmplx(8.d0,0.d0)/(one-rq)
     & + dcmplx(7.d0,0.d0))*cdlog(rq)
     & + one/two/(one-rq)
     &+ two*rli2  - dcmplx(5.d0/6.d0*pii**2,0.d0)
     &-dcmplx(5.d0,0.d0)*cdlog(one-rq) + two*cdlog(one-rq)**2 + eta/two
     &+alphaterm
     &+cdlog(muq/mtq)*(dcmplx(5.d0,0.d0)
     & -dcmplx(4.d0,0.d0)*cdlog(one-rq))
     & /two + cdlog(muq/mtq)**2/two)

      end subroutine Iop_me







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FROM THIS POINT ON IT BECOMES MAPLE-FORTRAN77-MESSY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FROM THIS POINT ON IT BECOMES MAPLE-FORTRAN77-MESSY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FROM THIS POINT ON IT BECOMES MAPLE-FORTRAN77-MESSY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FROM THIS POINT ON IT BECOMES MAPLE-FORTRAN77-MESSY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FROM THIS POINT ON IT BECOMES MAPLE-FORTRAN77-MESSY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FROM THIS POINT ON IT BECOMES MAPLE-FORTRAN77-MESSY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FROM THIS POINT ON IT BECOMES MAPLE-FORTRAN77-MESSY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FROM THIS POINT ON IT BECOMES MAPLE-FORTRAN77-MESSY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FROM THIS POINT ON IT BECOMES MAPLE-FORTRAN77-MESSY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FROM THIS POINT ON IT BECOMES MAPLE-FORTRAN77-MESSY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FROM THIS POINT ON IT BECOMES MAPLE-FORTRAN77-MESSY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FROM THIS POINT ON IT BECOMES MAPLE-FORTRAN77-MESSY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FROM THIS POINT ON IT BECOMES MAPLE-FORTRAN77-MESSY





      subroutine CalcANLO(vres,rres,ctmres,ctzTres,xe)
      implicit none
      include 'commondecay.f'
      include 'intcommonA.f'
      double complex Diags(8), DiagsRat(8), CTmass(8), CTwave(8)
      double complex vres(4), rres(4),res(4), zwres(4)
      double complex ctmres(4), ctzTres(4)
      integer j, xe

      do j=1,4
         vres(j) = (0.d0,0.d0)
         rres(j) = (0.d0,0.d0)
         ctmres(j) = (0.d0,0.d0)
         ctzTres(j) = (0.d0,0.d0)
      enddo

      call numATopDecVirt(Diags)
         zwres(1:4) = Diags(1)*AME1(1:4) + Diags(2)*AME2(1:4)
     &    	    + Diags(3)*AME3(1:4) + Diags(4)*AME4(1:4)
     &              + Diags(5)*AME5(1:4) + Diags(6)*AME6(1:4)
     &              + Diags(7)*AME7(1:4) + Diags(8)*AME8(1:4)
         vres(1:4) = zwres(1:4)
      If(xe .eq. 0) then
         call numATopDecRat(DiagsRat)
         zwres(1:4) = DiagsRat(1)*AME1(1:4) + DiagsRat(2)*AME2(1:4)
     &    	    + DiagsRat(3)*AME3(1:4) + DiagsRat(4)*AME4(1:4)
     &              + DiagsRat(5)*AME5(1:4) + DiagsRat(6)*AME6(1:4)
     &              + DiagsRat(7)*AME7(1:4) + DiagsRat(8)*AME8(1:4)
         rres(1:4) = zwres(1:4) +rres(1:4)
      endif

      If( xe .eq. -1 .or. xe .eq. 0) then
         call numATopDecCTmass(CTmass)
         ctmres(1:4) = CTmass(1)*AME1(1:4) + CTmass(2)*AME2(1:4)
     &    	     + CTmass(3)*AME3(1:4) + CTmass(4)*AME4(1:4)
     &               + CTmass(5)*AME5(1:4) + CTmass(6)*AME6(1:4)
     &               + CTmass(7)*AME7(1:4) + CTmass(8)*AME8(1:4)

         call numATopDecCTwave(CTwave)
         ctzTres(1:4) = CTwave(1)*AME1(1:4) + CTwave(2)*AME2(1:4)
     &      	      + CTwave(3)*AME3(1:4) + CTwave(4)*AME4(1:4)
     &                + CTwave(5)*AME5(1:4) + CTwave(6)*AME6(1:4)
     &                + CTwave(7)*AME7(1:4) + CTwave(8)*AME8(1:4)
      endif


      If(xe.eq. -1) then
         ctmres(1:4) = -mt*3.d0*ctmres(1:4)/4.d0/pii*alphas
	 ctzTres(1:4) = -3.d0*ctzTres(1:4)/4.d0/pii*alphas
      endif

      If(xe.eq. 0) then
         ctmres(1:4) = -mt*ctmres(1:4)*(4.d0+3.d0*log(mu**2/mt**2) )
     & *alphas/4.d0/pii
         ctzTres(1:4) = -ctzTres(1:4)*(4.d0+3.d0*log(mu**2/mt**2))
     & *alphas/4.d0/pii
      endif

      end subroutine CalcANLO

      subroutine CalcINTSMEA(xe,muRF,pri)
      implicit none
      include 'commondecay.f'
      include 'intcommonA.f'
      double precision muRq, Mwq, mtq, mzero, muFq,muRF
      double complex k1bbp1tb, k2wp1tb, k1bbk2w
      double complex qlI1,qlI2, qlI3, qlI4
      double complex p1tb(4), k2w(4), k1bb(4)
      double precision Rk1bbp1tb, Rk2wp1tb, Rk1bbk2w
      double precision p1qB, p1qC,p2qC,p1p2C
      double precision p1qD, p2qD,p3qD,p1p2D,p2p3D
      integer xe, pri, Dv, i
      mu =muRF
      Dv = 4
      do i = 1,4
         p1tb(i)=mom_t(i)
         k2w(i)=mom_dec(2,i)
         k1bb(i)=mom_dec(1,i)
      enddo
      call Calc_AME1(AME1)
      call Calc_AME2(AME2)
      call Calc_AME3(AME3)
      call Calc_AME4(AME4)
      call Calc_AME5(AME5)
      call Calc_AME6(AME6)
      call Calc_AME7(AME7)
      call Calc_AME8(AME8)
      mtq =mt**2
      muRq = mu**2
      muFq = mu**2
      Mwq = Mw**2
      mzero = 0.d0
      k1bbp1tb = scf(Dv,k1bb,p1tb)
      k1bbk2w = scf(Dv,k1bb,k2w)
      k2wp1tb = scf(Dv,k2w,p1tb)

      Rk1bbp1tb = dreal(k1bbp1tb)
      Rk1bbk2w = dreal(k1bbk2w)
      Rk2wp1tb = dreal(k2wp1tb)

C Massless Self-Energy (Diagramm 9)
      A1Sel100 = (0.d0,0.d0)
      A2Sel100 = (0.d0,0.d0)

      p1qB = mtq + mwq - 2.d0*Rk2wp1tb
      B012Sel100 = qlI2(p1qB,mzero,mzero,muRq,xe)

C Massive Self-Energy (Diagramm 10)
      A1Sel2T0 = qlI1(mtq,muRq,xe)
      A2Sel2T0 = (0.d0, 0.d0)

      p1qB = 2.d0*Rk1bbk2w + Mwq
      B012Sel2T0 = qlI2(p1qB,mtq,mzero,muRq,xe)

C Complete massless vertex (Diagramm 1)
      p1p2C = -2.d0*Rk2wp1tb + Mwq + mtq
      C0123Ver1000 = qlI3(mzero,mzero,p1p2C,mzero,mzero,mzero,muFq,xe)


C Vertex with 2 Tops (Diagramm 2)

      B013Ver2TT0 = qlI2(mtq,mtq,mzero,muRq,xe)

      p2qC = 2.d0*Rk1bbk2w + Mwq
      p1p2C = mtq
      C0123Ver2TT0 = qlI3(mzero,p2qC,p1p2C,mtq,mtq,mzero,muFq,xe)

C Vertex with 1 Top (Diagramm 3)
      B023Ver300T = qlI2(mwq,mzero,mtq,muRq,xe)

      p2qC = Mwq
      p1p2C = 2.d0*Rk1bbk2w + Mwq
      C0123Ver300T = qlI3(mzero,mwq,p1p2C,mzero,mzero,mtq,muFq,xe)

C Vertex with 1 Top (Diagramm 4)
      B013Ver4T00  = qlI2(mtq,mtq,mzero,muRq,xe)

      p1qC = Mwq
      p2qC = mtq + mwq - 2.d0*Rk2wp1tb
      p1p2C = mtq
      C0123Ver4T00 = qlI3(p1qC,p2qC,p1p2C,mtq,mzero,mzero,muFq,xe)

C Vertex with 1 Top (Diagramm 5 & 6)

      p1qB = mtq - 2.d0*Rk1bbp1tb
      B023Ver500T = qlI2(p1qB,mzero,mtq,muRq,xe)

      p2qC = mtq - 2.d0*Rk1bbp1tb
      p1p2C = mtq
      C0123Ver500T = qlI3(mzero,p2qC,p1p2C,mzero,mzero,mtq,muFq,xe)

C Box 1 with 2 Tops (Diagramm 8)

      p1p2C = mtq - 2.d0*Rk1bbp1tb
      C0234Box100TT = qlI3(mwq,mzero,p1p2C,mzero,mtq,mtq,muFq,xe)

      p1p2D = mwq + 2.d0*Rk1bbk2w
      p2p3D = mtq - 2.d0*Rk1bbp1tb
      D0Box100TT = qlI4(mzero,mwq,mzero,mtq,p1p2D,p2p3D,
     & mzero,mzero,mtq,mtq,muFq,xe)

C Box 1 with 1 Top (Diagramm 7)

      p2qC = mwq
      p1p2C = mtq - 2.d0*Rk1bbp1tb
      C0234Box2000T = qlI3(mzero,p2qC,p1p2C,mzero,mzero,mtq,muFq,xe)

      p1p2D =  mtq + mwq - 2.d0*Rk2wp1tb
      p2p3D = mtq - 2.d0*Rk1bbp1tb
      D0Box2000T = qlI4(mzero,mzero,mwq,mtq,p1p2D,p2p3D,
     & mzero,mzero,mzero,mtq,muFq,xe)

      If(pri .eq. 1 ) then

      write(*,*) "***************************** "
      write(*,*) "xe = ", xe
      write(*,*) "***************************** "

      write(*,*) "A1Sel100 = ", A1Sel100
      write(*,*) "A2Sel100 = ", A2Sel100
      write(*,*) "B012Sel100 = ", B012Sel100
      write(*,*) "A1Sel2T0 = ", A1Sel2T0
      write(*,*) "A2Sel2T0 = ", A2Sel2T0
      write(*,*) "B012Sel2T0 = ", B012Sel2T0
      write(*,*) "A1Ver1000 = ", A1Ver1000
      write(*,*) "A2Ver1000 = ", A2Ver1000
      write(*,*) "A3Ver1000 = ", A3Ver1000
      write(*,*) "B012Ver1000 = ", B012Ver1000
      write(*,*) "B013Ver1000 = ", B013Ver1000
      write(*,*) "B023Ver1000 = ", B023Ver1000
      write(*,*) "C0123Ver1000 = ", C0123Ver1000
      write(*,*) "A1Ver2TT0 = ", A1Ver2TT0
      write(*,*) "A2Ver2TT0 = ", A2Ver2TT0
      write(*,*) "A3Ver2TT0 = ", A3Ver2TT0
      write(*,*) "B012Ver2TT0 = ", B012Ver2TT0
      write(*,*) "B013Ver2TT0 = ", B013Ver2TT0
      write(*,*) "B023Ver2TT0 = ", B023Ver2TT0
      write(*,*) "C0123Ver2TT0 = ", C0123Ver2TT0
      write(*,*) "A1Ver300T = ", A1Ver300T
      write(*,*) "A2Ver300T = ", A2Ver300T
      write(*,*) "A3Ver300T = ", A3Ver300T
      write(*,*) "B012Ver300T = ", B012Ver300T
      write(*,*) "B013Ver300T = ", B013Ver300T
      write(*,*) "B023Ver300T = ", B023Ver300T
      write(*,*) "C0123Ver300T = ", C0123Ver300T
      write(*,*) "A1Ver4T00 = ", A1Ver4T00
      write(*,*) "A2Ver4T00 = ", A2Ver4T00
      write(*,*) "A3Ver4T00 = ", A3Ver4T00
      write(*,*) "B012Ver4T00 = ", B012Ver4T00
      write(*,*) "B013Ver4T00 = ", B013Ver4T00
      write(*,*) "B023Ver4T00 = ", B023Ver4T00
      write(*,*) "C0123Ver4T00 = ", C0123Ver4T00
      write(*,*) "A1Ver500T = ", A1Ver500T
      write(*,*) "A2Ver500T = ", A2Ver500T
      write(*,*) "A3Ver500T = ", A3Ver500T
      write(*,*) "B012Ver500T = ", B012Ver500T
      write(*,*) "B013Ver500T = ", B013Ver500T
      write(*,*) "B023Ver500T = ", B023Ver500T
      write(*,*) "C0123Ver500T = ", C0123Ver500T
      write(*,*) "A1Box100TT = ", A1Box100TT
      write(*,*) "A2Box100TT = ", A2Box100TT
      write(*,*) "A3Box100TT = ", A3Box100TT
      write(*,*) "A4Box100TT = ", A4Box100TT
      write(*,*) "B012Box100TT = ", B012Box100TT
      write(*,*) "B013Box100TT = ", B013Box100TT
      write(*,*) "B014Box100TT = ", B014Box100TT
      write(*,*) "B023Box100TT = ", B023Box100TT
      write(*,*) "B024Box100TT = ", B024Box100TT
      write(*,*) "B034Box100TT = ", B034Box100TT
      write(*,*) "C0123Box100TT = ", C0123Box100TT
      write(*,*) "C0124Box100TT = ", C0124Box100TT
      write(*,*) "C0134Box100TT = ", C0134Box100TT
      write(*,*) "C0234Box100TT = ", C0234Box100TT
      write(*,*) "D0Box100TT = ", D0Box100TT
      write(*,*) "A1Box2000T = ", A1Box2000T
      write(*,*) "A2Box2000T = ", A2Box2000T
      write(*,*) "A3Box2000T = ", A3Box2000T
      write(*,*) "A4Box2000T = ", A4Box2000T
      write(*,*) "B012Box2000T = ", B012Box2000T
      write(*,*) "B013Box2000T = ", B013Box2000T
      write(*,*) "B014Box2000T = ", B014Box2000T
      write(*,*) "B023Box2000T = ", B023Box2000T
      write(*,*) "B024Box2000T = ", B024Box2000T
      write(*,*) "B034Box2000T = ", B034Box2000T
      write(*,*) "C0123Box2000T = ", C0123Box2000T
      write(*,*) "C0124Box2000T = ", C0124Box2000T
      write(*,*) "C0134Box2000T = ", C0134Box2000T
      write(*,*) "C0234Box2000T = ", C0234Box2000T
      write(*,*) "D0Box2000T = ", D0Box2000T
      endif
      end subroutine CalcINTSMEA

      subroutine numATopDecRat(DiagRat)
      implicit none
      include 'VarATopEps.f'
      include 'commondecay.f'
      include 'intcommonA.f'
      double complex p1tb(4),k2w(4),k1bb(4)
      double complex EpsP(4),EpsW(4),SpB(4),SpT(4)
      double complex k1bbp1tb, k2wp1tb, k1bbk2w
      double complex EpsWk1bb, EpsWp1tb, EpsWEpsP
      double complex EpsPk1bb, EpsPp1tb, EpsPk2w
      double complex DiagRat(8)
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
      double complex s11,s12,s13,s14,s15,s16,s17,s18,s19,s20
      double complex s21,s22,s23,s24,s25,s26,s27,s28,s29,s30
      double complex s31,s32,s33,s34,s35,s36,s37,s38,s39,s40
      integer i, Dv
      Dv = 4
      do i = 1,4
         p1tb(i)=mom_t(i)
         k2w(i)=mom_dec(2,i)
         k1bb(i)=mom_dec(1,i)
         EpsP(i)=hel_dec(3,i)
         EpsW(i)=hel_dec(2,i)
         SpB(i)=hel_dec(1,i)
      enddo
      EpsWk1bb = scf(Dv,EpsW,k1bb)
      EpsWp1tb = scf(Dv,EpsW,p1tb)
      EpsWEpsP = scf(Dv,EpsW,EpsP)
      EpsPk1bb = scf(Dv,EpsP,k1bb)
      EpsPp1tb = scf(Dv,EpsP,p1tb)
      EpsPk2w = scf(Dv,EpsP,k2w)
      k1bbp1tb = scf(Dv,k1bb,p1tb)
      k1bbk2w = scf(Dv,k1bb,k2w)
      k2wp1tb = scf(Dv,k2w,p1tb)
      t1 = EL**2
      t2 = alphaS*t1
      t3 = t2*GW
      t4 = EpsPk1bb*EpsWk1bb
      t5 = mt**2
      t6 = t5**2
      t7 = t6**2
      t12 = k1bbk2w**2
      t13 = t12*k1bbk2w
      t14 = k1bbp1tb**2
      t15 = t13*t14
      t18 = EpsPk2w*EpsWk1bb
      t21 = EpsPk2w*EpsWp1tb
      t24 = t14*k1bbp1tb
      t25 = t12*t24
      t32 = t14**2
      t33 = k1bbk2w*t32
      t40 = t32*t5
      t43 = t4*k1bbk2w*t7+t4*k1bbp1tb*t7+16*t4*t15-16*t18*t15+16*t21*t15
     #-16*t4*t25+48*t18*t25-32*t21*t25-16*t4*t33-48*t18*t33+16*t21*t33-3
     #2*t4*t40
      t49 = t24*t6
      t54 = t6*t5
      t58 = t14*t54
      t63 = Mw**2
      t64 = t32*t63
      t69 = EpsPk1bb*EpsWp1tb
      t74 = t63**2
      t75 = t24*t74
      t80 = -24*t18*t40+4*t4*t13*t6+24*t4*t49+12*t18*t49+4*t4*t12*t54-8*
     #t4*t58-2*t18*t58-8*t4*t64-24*t18*t64+8*t69*t64+16*t21*t64-4*t4*t75
     #+12*t18*t75
      t86 = t74*t63
      t87 = t14*t86
      t95 = k1bbp1tb*t5
      t96 = t95*t63
      t99 = t18*t12
      t102 = t32*k1bbp1tb
      t107 = k1bbk2w*t6
      t111 = k1bbp1tb*t6
      t115 = k1bbk2w*k1bbp1tb
      t116 = t115*t86
      t121 = -8*t69*t75-24*t21*t75+2*t4*t87-2*t18*t87-6*t69*t87-4*t4*t12
     #*t96+16*t99*t96+16*t4*t102+16*t18*t102+3*t4*t107*t74+3*t4*t111*t74
     #+2*t4*t116+2*t18*t116
      t124 = t13*k1bbp1tb
      t125 = t124*t5
      t132 = t12*t14
      t133 = t132*t5
      t140 = k1bbk2w*t24
      t141 = t140*t5
      t148 = t12*k1bbp1tb
      t149 = t148*t6
      t156 = 2*t69*t116-16*t4*t125+8*t18*t125-8*t21*t125+32*t4*t133-40*t
     #18*t133+24*t21*t133+16*t4*t141+56*t18*t141-16*t21*t141-20*t4*t149+
     #8*t18*t149-4*t21*t149
      t159 = k1bbk2w*t14
      t160 = t159*t6
      t165 = t115*t54
      t170 = t124*t63
      t177 = t132*t63
      t186 = t140*t63
      t189 = -20*t18*t160+4*t21*t160-4*t4*t165+2*t18*t165+8*t4*t170+8*t1
     #8*t170-8*t21*t170+8*t4*t177-40*t18*t177+8*t69*t177+56*t21*t177-8*t
     #4*t186
      t196 = t13*t5
      t200 = t24*t5
      t201 = t200*t63
      t210 = t14*t6
      t211 = t210*t63
      t218 = k1bbk2w*t5
      t223 = 56*t18*t186-16*t69*t186-64*t21*t186-4*t4*t196*t63-4*t4*t201
     #+8*t18*t201-8*t69*t201-8*t21*t201+10*t4*t211+2*t18*t211+2*t69*t211
     #-t4*t218*t86-t4*t95*t86
      t225 = k1bbk2w*t54
      t229 = k1bbp1tb*t54
      t233 = t148*t74
      t242 = t159*t74
      t251 = t12*t5
      t255 = t14*t5
      t256 = t255*t74
      t261 = -3*t4*t225*t63-3*t4*t229*t63+8*t4*t233+8*t18*t233+4*t69*t23
     #3-4*t21*t233+4*t4*t242-20*t18*t242-4*t69*t242+28*t21*t242-4*t4*t25
     #1*t74-4*t4*t256+2*t18*t256
      t267 = t21*t12
      t270 = t4*k1bbk2w
      t271 = t255*t63
      t274 = t18*k1bbk2w
      t277 = t69*k1bbk2w
      t280 = t21*k1bbk2w
      t283 = t111*t63
      t290 = t95*t74
      t297 = 4*t69*t256-4*t69*t12*t96-24*t267*t96-4*t270*t271-24*t274*t2
     #71+20*t277*t271+32*t280*t271+10*t270*t283-2*t274*t283-6*t277*t283-
     #8*t270*t290-2*t274*t290+4*t277*t290
      t304 = 1/k1bbp1tb
      t305 = -k1bbk2w+k1bbp1tb
      t308 = 2*k1bbk2w
      t309 = 2*k1bbp1tb
      t311 = 2*mt*Mw
      t312 = -t308+t309-t5-t311-t63
      t315 = -t308+t309-t5+t311-t63
      t317 = t309-t5+t63
      t318 = 1/t317
      t321 = k1bbp1tb*t63
      t323 = 1/(2*t115-t218+t321)
      t324 = 1/pii
      t331 = t2*EpsWEpsP*GW
      t332 = Qqd-Qqu
      t333 = mt*t332
      t342 = EpsWEpsP*k1bbk2w
      t354 = mt*Qqu
      t355 = -t308+t5-t63
      t358 = t354/t355*t324
      t362 = -t308+5*t5-t63
      t364 = t355**2
      t365 = 1/t364
      t375 = 1/k1bbk2w
      t376 = -1/t355
      t378 = t375*t376*t324
      t379 = t354*t378
      t423 = 4*t4*t148+8*t18*t148-8*t21*t148-4*t4*t159-8*t18*t159+8*t69*
     #t159+8*t21*t159-2*t4*t251-4*t69*t115*t5-2*t4*t255+t4*t107+t4*t111+
     #2*t4*t12*t63+4*t69*t115*t63+2*t4*t14*t63-2*t4*t218*t63-2*t4*t96+t4
     #*k1bbk2w*t74+t4*k1bbp1tb*t74
      t427 = -1/t305
      t428 = t375*t427
      t431 = t376*t318*t324
      t439 = t427*t324
      t440 = GW*Qqd*t439
      t443 = t2*EpsPk2w
      t446 = GW*t332*t318*t324
      t453 = t365*t324
      t454 = (t308-3*t5+t63)*Qqu*t453
      t461 = GW*Qqu*t376*t324
      t491 = t12*t6
      t504 = 16*t18*t124-16*t21*t124-16*t4*t132-48*t18*t132+32*t21*t132+
     #32*t4*t140+48*t18*t140-16*t21*t140-8*t18*t196+8*t21*t196+32*t4*t20
     #0+24*t18*t200-4*t4*t491-8*t18*t491+4*t21*t491-24*t4*t210-12*t18*t2
     #10-4*t4*t225
      t511 = t24*t63
      t523 = t14*t74
      t533 = k1bbp1tb*t86
      t547 = -2*t18*t225+8*t4*t229+2*t18*t229+8*t4*t511+16*t18*t511-8*t6
     #9*t511-16*t21*t511+3*t4*t54*t63+4*t4*t523-4*t18*t523+8*t21*t523-3*
     #t4*t6*t74-2*t4*t533+2*t69*t533+t4*t5*t86-16*t18*t32-16*t4*t32-t4*t
     #7+2*t18*t290
      t551 = t218*t74
      t558 = t115*t74
      t573 = t107*t63
      t584 = t251*t63
      t589 = 4*t4*t290-4*t69*t551-2*t18*t551+4*t4*t551-8*t21*t558+4*t69*
     #t558+4*t18*t558-8*t4*t558-2*t69*t283-4*t18*t283-10*t4*t283+4*t69*t
     #573+4*t18*t573+8*t21*t271+8*t69*t271+4*t4*t271+12*t21*t584-8*t18*t
     #584
      t592 = t159*t63
      t599 = t148*t63
      t606 = t115*t6
      t613 = t159*t5
      t620 = t148*t5
      t633 = 4*t4*t584+40*t21*t592+8*t69*t592-32*t18*t592-24*t21*t599+16
     #*t18*t599-8*t4*t599-4*t21*t606+24*t4*t606+20*t18*t606-48*t4*t613-5
     #6*t18*t613+16*t21*t613+16*t4*t620+40*t18*t620-24*t21*t620+8*t274*t
     #96-12*t277*t96-20*t280*t96
      t640 = -1/t315
      t642 = -1/t312
      t648 = t2*EpsWEpsP
      t667 = k1bbk2w*t63
      t672 = t5*t63
      t676 = 4*t99-4*t267-4*t18*t115+4*t69*t115+4*t21*t115-2*t69*t218-2*
     #t4*t95+t4*t6+2*t69*t667+2*t4*t321-2*t4*t672+t4*t74
      t688 = EpsWp1tb*GW
      t701 = Qqd*t427
      t703 = t640*t642*t324
      t711 = t2*EpsWk1bb*GW
      DiagRat(1) = ne*(-t3*mt*(t43+t80+t121+t156+t189+t223+t261+t297)*Qq
     #d*t304/t305/t312/t315*t318*t323*t324/4+t331*t333/t317*t324/2+t2*GW
     #*(-2*t4-2*t18+2*t69+2*t342+5*EpsWEpsP*k1bbp1tb)*t333*t304*t318*t32
     #4/4-2*t331*t358+t331*mt*t362*Qqu*t365*t324/2+t2*GW*(-t18+4*t342)*t
     #379/2-t3*mt*t423*Qqu*t428*t304*t431/4)
      DiagRat(2) = ne*(-t2*EpsPk1bb*t440/4-t443*t446+t2*EpsPk2w*GW*t454/
     #2+t2*(EpsPk1bb-EpsPk2w)*t461/2-t443*t461)
      DiagRat(3) = ne*(-t2*GW*(t504+t547+t589+t633)*Qqd*t427*t318*t640*t
     #642*t323*t324/4+t648*t446-t331*t454/2+2*t648*t461-t2*GW*t676*Qqu*t
     #428*t431/4)
      DiagRat(4) = ne*(t2*EpsWp1tb*t440/4+t2*t688*(4*t12-8*t115+4*t14+4*
     #t218-4*t95+t6+6*t667-6*t321-3*t672+2*t74)*t701*t703/4-t2*(EpsWk1bb
     #-EpsWp1tb)*t446+t711*t454/2-t2*EpsWk1bb*t461-t711*(t308-t63)*Qqu*t
     #378/4)
      DiagRat(5) = (0.d0,0.d0)
      DiagRat(6) = ne*(3.D0/4.D0*t3*mt*Qqd*t439+2*t3*t358-t2*GW*mt*t362*
     #Qqu*t453/4)
      DiagRat(7) = ne*(-t2*t688*mt*(t308-t309+t5-t63)*Qqd*t427*t703+t711
     #*t379)/4
      DiagRat(8) = ne*(-3.D0/8.D0*t3*t701*t324-t3*t454/4+t3*Qqu*t376*t32
     #4)

      end subroutine numATopDecRat

      subroutine numATopDecCTmass(MassCT)
      implicit none
      include 'VarATopMassCT.f'
      include 'commondecay.f'
      include 'intcommonA.f'
      double complex p1tb(4),k2w(4),k1bb(4)
      double complex EpsP(4),EpsW(4),SpB(4),SpT(4)
      double complex k1bbp1tb, k2wp1tb, k1bbk2w
      double complex EpsWk1bb, EpsWp1tb, EpsWEpsP
      double complex EpsPk1bb, EpsPp1tb, EpsPk2w
      double complex MassCT(8)
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
      double complex s11,s12,s13,s14,s15,s16,s17,s18,s19,s20
      double complex s21,s22,s23,s24,s25,s26,s27,s28,s29,s30
      double complex s31,s32,s33,s34,s35,s36,s37,s38,s39,s40
      integer i, Dv
      Dv = 4
      do i = 1,4
         p1tb(i)=mom_t(i)
         k2w(i)=mom_dec(2,i)
         k1bb(i)=mom_dec(1,i)
         EpsP(i)=hel_dec(3,i)
         EpsW(i)=hel_dec(2,i)
         SpB(i)=hel_dec(1,i)
      enddo
      EpsWk1bb = scf(Dv,EpsW,k1bb)
      EpsWp1tb = scf(Dv,EpsW,p1tb)
      EpsWEpsP = scf(Dv,EpsW,EpsP)
      EpsPk1bb = scf(Dv,EpsP,k1bb)
      EpsPp1tb = scf(Dv,EpsP,p1tb)
      EpsPk2w = scf(Dv,EpsP,k2w)
      k1bbp1tb = scf(Dv,k1bb,p1tb)
      k1bbk2w = scf(Dv,k1bb,k2w)
      k2wp1tb = scf(Dv,k2w,p1tb)
      t1 = EL**2
      t3 = t1*EpsWEpsP*GW
      t4 = k1bbk2w**2
      t7 = mt**2
      t11 = t7**2
      t13 = Mw**2
      t17 = t7*t13
      t20 = t13**2
      t38 = -4*t4*Qqd+4*k1bbk2w*t7*Qqd-t11*Qqd-4*k1bbk2w*t13*Qqd+2*t17*Q
     #qd-t20*Qqd+4*t4*Qqu-8*k1bbk2w*k1bbp1tb*Qqu-12*k1bbp1tb*t7*Qqu+7*t1
     #1*Qqu-4*k1bbp1tb*t13*Qqu-6*t17*Qqu-t20*Qqu
      t43 = 2*k1bbk2w
      t45 = (t43-t7+t13)**2
      t46 = 1/t45
      t54 = mt*ne*Qqu*t46
      t63 = t1*GW
      t68 = ne*Qqu*t46
      MassCT(1) = -t3*ne*t38/(-2*k1bbp1tb+t7-t13)*t46
      MassCT(2) = 4*t1*EpsPk2w*GW*t54
      MassCT(3) = -4*t3*t54
      MassCT(4) = 4*t1*EpsWk1bb*GW*t54
      MassCT(5) = (0.d0,0.d0)
      MassCT(6) = t63*(t43+3*t7+t13)*t68
      MassCT(7) = (0.d0,0.d0)
      MassCT(8) = -2*t63*mt*t68

      end subroutine numATopDecCTmass

      subroutine numATopDecCTwave(WaveCT)
      implicit none
      include 'VarATopWaveCT.f'
      include 'commondecay.f'
      include 'intcommonA.f'
      double complex p1tb(4),k2w(4),k1bb(4)
      double complex EpsP(4),EpsW(4),SpB(4),SpT(4)

      double complex k1bbp1tb, k2wp1tb, k1bbk2w
      double complex EpsWk1bb, EpsWp1tb, EpsWEpsP
      double complex EpsPk1bb, EpsPp1tb, EpsPk2w
      double complex WaveCT(8)
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
      double complex s11,s12,s13,s14,s15,s16,s17,s18,s19,s20
      double complex s21,s22,s23,s24,s25,s26,s27,s28,s29,s30
      double complex s31,s32,s33,s34,s35,s36,s37,s38,s39,s40
      integer i, Dv
      Dv = 4
      do i = 1,4
         p1tb(i)=mom_t(i)
         k2w(i)=mom_dec(2,i)
         k1bb(i)=mom_dec(1,i)
         EpsP(i)=hel_dec(3,i)
         EpsW(i)=hel_dec(2,i)
         SpB(i)=hel_dec(1,i)
      enddo
      EpsWk1bb = scf(Dv,EpsW,k1bb)
      EpsWp1tb = scf(Dv,EpsW,p1tb)
      EpsWEpsP = scf(Dv,EpsW,EpsP)
      EpsPk1bb = scf(Dv,EpsP,k1bb)
      EpsPp1tb = scf(Dv,EpsP,p1tb)
      EpsPk2w = scf(Dv,EpsP,k2w)
      k1bbp1tb = scf(Dv,k1bb,p1tb)
      k1bbk2w = scf(Dv,k1bb,k2w)
      k2wp1tb = scf(Dv,k2w,p1tb)
      t1 = EL**2
      t2 = t1*EpsWEpsP
      t7 = mt**2
      t9 = Mw**2
      t15 = -2*k1bbk2w*Qqd+t7*Qqd-t9*Qqd+2*k1bbk2w*Qqu-2*k1bbp1tb*Qqu
      t16 = ne*t15
      t19 = 1/(-2*k1bbp1tb+t7-t9)
      t21 = 2*k1bbk2w-t7+t9
      t24 = t16*t19/t21
      t32 = t1*GW
      t45 = 1/(k1bbk2w-k1bbp1tb)
      t46 = t15*t45
      t47 = -1/t21
      WaveCT(1) = -2*t2*GW*mt*t24
      WaveCT(2) = t1*EpsPk2w*GW*t24
      WaveCT(3) = -t2*GW*t24
      WaveCT(4) = -t32*(2*EpsWk1bb*k1bbk2w-2*EpsWp1tb*k1bbk2w-2*EpsWk1bb
     #*k1bbp1tb+EpsWp1tb*t7-EpsWp1tb*t9)*ne*t46*t47*t19/2
      WaveCT(5) = (0.d0,0.d0)
      WaveCT(6) = -t32*mt*t16*t45*t47/2
      WaveCT(7) = (0.d0,0.d0)
      WaveCT(8) = t32*ne*t46*t47/4

      end subroutine numATopDecCTwave










      subroutine CalcNLO(vres,rres,ctmres,ctzTres,xe)
      implicit none
      include 'commondecay.f'
      include 'intcommon.f'
      double complex Diags(8), DiagsRat(8), CTmass(8), CTwave(8)
      double complex vres(4), rres(4),res(4), zwres(4)
      double complex ctmres(4), ctzTres(4)
      integer j, xe

      do j=1,4
         vres(j) = (0.d0,0.d0)
         rres(j) = (0.d0,0.d0)
         ctmres(j) = (0.d0,0.d0)
         ctzTres(j) = (0.d0,0.d0)
      enddo

      call numTopDecVirt(Diags)

         zwres(1:4) = Diags(1)*ME1(1:4) + Diags(2)*ME2(1:4)
     &  	    + Diags(3)*ME3(1:4) + Diags(4)*ME4(1:4)
     &              + Diags(5)*ME5(1:4) + Diags(6)*ME6(1:4)
     &              + Diags(7)*ME7(1:4) + Diags(8)*ME8(1:4)
         vres(1:4) = zwres(1:4)
      If(xe .eq. 0) then
         call numTopDecRat(DiagsRat)
         zwres(1:4) = DiagsRat(1)*ME1(1:4) + DiagsRat(2)*ME2(1:4)
     &    	    + DiagsRat(3)*ME3(1:4) + DiagsRat(4)*ME4(1:4)
     &              + DiagsRat(5)*ME5(1:4) + DiagsRat(6)*ME6(1:4)
     &              + DiagsRat(7)*ME7(1:4) + DiagsRat(8)*ME8(1:4)
         rres(1:4) = zwres(1:4)
      endif

      If( xe .eq. -1 .or. xe .eq. 0) then
         call numTopDecCTmass(CTmass)
         ctmres(1:4) = CTmass(1)*ME1(1:4) + CTmass(2)*ME2(1:4)
     &    	     + CTmass(3)*ME3(1:4) + CTmass(4)*ME4(1:4)
     &               + CTmass(5)*ME5(1:4) + CTmass(6)*ME6(1:4)
     &               + CTmass(7)*ME7(1:4) + CTmass(8)*ME8(1:4)

         call numTopDecCTwave(CTwave)
         ctzTres(1:4) = CTwave(1)*ME1(1:4) + CTwave(2)*ME2(1:4)
     &      	      + CTwave(3)*ME3(1:4) + CTwave(4)*ME4(1:4)
     &                + CTwave(5)*ME5(1:4) + CTwave(6)*ME6(1:4)
     &                + CTwave(7)*ME7(1:4) + CTwave(8)*ME8(1:4)

      endif


      If(xe.eq. -1) then
         ctmres(1:4) = -mt*3.d0*ctmres(1:4)/4.d0/pii*alphas
	 ctzTres(1:4) = -3.d0*ctzTres(1:4)/4.d0/pii*alphas
      endif

      If(xe.eq. 0) then
         ctmres(1:4) = -mt*ctmres(1:4)*(4.d0+3.d0*log(mu**2/mt**2))
     & *alphas/4.d0/pii
         ctzTres(1:4) = -ctzTres(1:4)*(4.d0+3.d0*log(mu**2/mt**2))
     & *alphas/4.d0/pii
      endif
      end subroutine CalcNLO

      subroutine CalcINTSME(xe,muRF,pri)
      implicit none
      include 'commondecay.f'
      include 'intcommon.f'
      double precision muRq, Mwq, mtq, mzero, muFq,muRF
      double complex k1bp1t, k2wp1t, k1bk2w
      double complex qlI1,qlI2, qlI3, qlI4
      double complex p1t(4), k2w(4), k1b(4)
      double precision Rk1bp1t, Rk2wp1t, Rk1bk2w
      double precision p1qB, p1qC,p2qC,p1p2C
      double precision p1qD, p2qD,p3qD,p1p2D,p2p3D
      integer xe, pri, Dv, i
      mu =muRF
      Dv = 4
      do i = 1,4
         p1t(i)=mom_t(i)
         k2w(i)=mom_dec(2,i)
         k1b(i)=mom_dec(1,i)
      enddo
      call Calc_ME1(ME1)
      call Calc_ME2(ME2)
      call Calc_ME3(ME3)
      call Calc_ME4(ME4)
      call Calc_ME5(ME5)
      call Calc_ME6(ME6)
      call Calc_ME7(ME7)
      call Calc_ME8(ME8)
      mtq =mt**2
      muRq = mu**2
      muFq = mu**2
      Mwq = Mw**2
      mzero = 0.d0
      k1bp1t = scf(Dv,k1b,p1t)
      k1bk2w = scf(Dv,k1b,k2w)
      k2wp1t = scf(Dv,k2w,p1t)

      Rk1bp1t = dreal(k1bp1t)
      Rk1bk2w = dreal(k1bk2w)
      Rk2wp1t = dreal(k2wp1t)

C Massless Self-Energy (Diagramm 9)
      A1Sel100 = (0.d0,0.d0)
      A2Sel100 = (0.d0,0.d0)

      p1qB = mtq + mwq - 2.d0*Rk2wp1t
      B012Sel100 = qlI2(p1qB,mzero,mzero,muRq,xe)

C Massive Self-Energy (Diagramm 10)
      A1Sel2T0 = qlI1(mtq,muRq,xe)
      A2Sel2T0 = (0.d0, 0.d0)

      p1qB = 2.d0*Rk1bk2w + Mwq
      B012Sel2T0 = qlI2(p1qB,mtq,mzero,muRq,xe)

C Complete massless vertex (Diagramm 1)
      p1p2C = -2.d0*Rk2wp1t + Mwq + mtq
      C0123Ver1000 = qlI3(mzero,mzero,p1p2C,mzero,mzero,mzero,muFq,xe)


C Vertex with 2 Tops (Diagramm 2)

      B013Ver2TT0 = qlI2(mtq,mtq,mzero,muRq,xe)

      p2qC = 2.d0*Rk1bk2w + Mwq
      p1p2C = mtq
      C0123Ver2TT0 = qlI3(mzero,p2qC,p1p2C,mtq,mtq,mzero,muFq,xe)

C Vertex with 1 Top (Diagramm 3)
      B023Ver300T = qlI2(mwq,mzero,mtq,muRq,xe)

      p2qC = Mwq
      p1p2C = 2.d0*Rk1bk2w + Mwq
      C0123Ver300T = qlI3(mzero,mwq,p1p2C,mzero,mzero,mtq,muFq,xe)

C Vertex with 1 Top (Diagramm 4)
      B013Ver4T00  = qlI2(mtq,mtq,mzero,muRq,xe)

      p1qC = Mwq
      p2qC = mtq + mwq - 2.d0*Rk2wp1t
      p1p2C = mtq
      C0123Ver4T00 = qlI3(p1qC,p2qC,p1p2C,mtq,mzero,mzero,muFq,xe)

C Vertex with 1 Top (Diagramm 5 & 6)

      p1qB = mtq - 2.d0*Rk1bp1t
      B023Ver500T = qlI2(p1qB,mzero,mtq,muRq,xe)

      p2qC = mtq - 2.d0*Rk1bp1t
      p1p2C = mtq
      C0123Ver500T = qlI3(mzero,p2qC,p1p2C,mzero,mzero,mtq,muFq,xe)

C Box 1 with 2 Tops (Diagramm 8)

      p1p2C = mtq - 2.d0*Rk1bp1t
      C0234Box100TT = qlI3(mwq,mzero,p1p2C,mzero,mtq,mtq,muFq,xe)

      p1p2D = mwq + 2.d0*Rk1bk2w
      p2p3D = mtq - 2.d0*Rk1bp1t
      D0Box100TT = qlI4(mzero,mwq,mzero,mtq,p1p2D,p2p3D,
     & mzero,mzero,mtq,mtq,muFq,xe)

C Box 1 with 1 Top (Diagramm 7)

      p2qC = mwq
      p1p2C = mtq - 2.d0*Rk1bp1t
      C0234Box2000T = qlI3(mzero,p2qC,p1p2C,mzero,mzero,mtq,muFq,xe)

      p1p2D =  mtq + mwq - 2.d0*Rk2wp1t
      p2p3D = mtq - 2.d0*Rk1bp1t
      D0Box2000T = qlI4(mzero,mzero,mwq,mtq,p1p2D,p2p3D,
     & mzero,mzero,mzero,mtq,muFq,xe)
      If(pri .eq. 1 ) then

      write(*,*) "***************************** "
      write(*,*) "xe = ", xe
      write(*,*) "***************************** "

      write(*,*) "A1Sel100 = ", A1Sel100
      write(*,*) "A2Sel100 = ", A2Sel100
      write(*,*) "B012Sel100 = ", B012Sel100
      write(*,*) "A1Sel2T0 = ", A1Sel2T0
      write(*,*) "A2Sel2T0 = ", A2Sel2T0
      write(*,*) "B012Sel2T0 = ", B012Sel2T0
      write(*,*) "A1Ver1000 = ", A1Ver1000
      write(*,*) "A2Ver1000 = ", A2Ver1000
      write(*,*) "A3Ver1000 = ", A3Ver1000
      write(*,*) "B012Ver1000 = ", B012Ver1000
      write(*,*) "B013Ver1000 = ", B013Ver1000
      write(*,*) "B023Ver1000 = ", B023Ver1000
      write(*,*) "C0123Ver1000 = ", C0123Ver1000
      write(*,*) "A1Ver2TT0 = ", A1Ver2TT0
      write(*,*) "A2Ver2TT0 = ", A2Ver2TT0
      write(*,*) "A3Ver2TT0 = ", A3Ver2TT0
      write(*,*) "B012Ver2TT0 = ", B012Ver2TT0
      write(*,*) "B013Ver2TT0 = ", B013Ver2TT0
      write(*,*) "B023Ver2TT0 = ", B023Ver2TT0
      write(*,*) "C0123Ver2TT0 = ", C0123Ver2TT0
      write(*,*) "A1Ver300T = ", A1Ver300T
      write(*,*) "A2Ver300T = ", A2Ver300T
      write(*,*) "A3Ver300T = ", A3Ver300T
      write(*,*) "B012Ver300T = ", B012Ver300T
      write(*,*) "B013Ver300T = ", B013Ver300T
      write(*,*) "B023Ver300T = ", B023Ver300T
      write(*,*) "C0123Ver300T = ", C0123Ver300T
      write(*,*) "A1Ver4T00 = ", A1Ver4T00
      write(*,*) "A2Ver4T00 = ", A2Ver4T00
      write(*,*) "A3Ver4T00 = ", A3Ver4T00
      write(*,*) "B012Ver4T00 = ", B012Ver4T00
      write(*,*) "B013Ver4T00 = ", B013Ver4T00
      write(*,*) "B023Ver4T00 = ", B023Ver4T00
      write(*,*) "C0123Ver4T00 = ", C0123Ver4T00
      write(*,*) "A1Ver500T = ", A1Ver500T
      write(*,*) "A2Ver500T = ", A2Ver500T
      write(*,*) "A3Ver500T = ", A3Ver500T
      write(*,*) "B012Ver500T = ", B012Ver500T
      write(*,*) "B013Ver500T = ", B013Ver500T
      write(*,*) "B023Ver500T = ", B023Ver500T
      write(*,*) "C0123Ver500T = ", C0123Ver500T
      write(*,*) "A1Box100TT = ", A1Box100TT
      write(*,*) "A2Box100TT = ", A2Box100TT
      write(*,*) "A3Box100TT = ", A3Box100TT
      write(*,*) "A4Box100TT = ", A4Box100TT
      write(*,*) "B012Box100TT = ", B012Box100TT
      write(*,*) "B013Box100TT = ", B013Box100TT
      write(*,*) "B014Box100TT = ", B014Box100TT
      write(*,*) "B023Box100TT = ", B023Box100TT
      write(*,*) "B024Box100TT = ", B024Box100TT
      write(*,*) "B034Box100TT = ", B034Box100TT
      write(*,*) "C0123Box100TT = ", C0123Box100TT
      write(*,*) "C0124Box100TT = ", C0124Box100TT
      write(*,*) "C0134Box100TT = ", C0134Box100TT
      write(*,*) "C0234Box100TT = ", C0234Box100TT
      write(*,*) "D0Box100TT = ", D0Box100TT
      write(*,*) "A1Box2000T = ", A1Box2000T
      write(*,*) "A2Box2000T = ", A2Box2000T
      write(*,*) "A3Box2000T = ", A3Box2000T
      write(*,*) "A4Box2000T = ", A4Box2000T
      write(*,*) "B012Box2000T = ", B012Box2000T
      write(*,*) "B013Box2000T = ", B013Box2000T
      write(*,*) "B014Box2000T = ", B014Box2000T
      write(*,*) "B023Box2000T = ", B023Box2000T
      write(*,*) "B024Box2000T = ", B024Box2000T
      write(*,*) "B034Box2000T = ", B034Box2000T
      write(*,*) "C0123Box2000T = ", C0123Box2000T
      write(*,*) "C0124Box2000T = ", C0124Box2000T
      write(*,*) "C0134Box2000T = ", C0134Box2000T
      write(*,*) "C0234Box2000T = ", C0234Box2000T
      write(*,*) "D0Box2000T = ", D0Box2000T
      endif
      end subroutine CalcINTSME

      subroutine numTopDecRat(DiagRat)
      implicit none
      include 'VarTopEps.f'
      include 'commondecay.f'
      include 'intcommon.f'
      double precision wrat
      double complex p1t(4),k2w(4),k1b(4), k3p(4)
      double complex EpsP(4),EpsW(4),SpB(4),SpT(4)

      double complex k1bp1t, k2wp1t, k1bk2w
      double complex EpsWk1b, EpsWp1t, EpsWEpsP
      double complex EpsPk1b, EpsPp1t, EpsPk2w
      double complex DiagRat(8)
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
      double complex s11,s12,s13,s14,s15,s16,s17,s18,s19
      integer i, Dv
      Dv = 4
      do i = 1,4
         p1t(i)=mom_t(i)
         k2w(i)=mom_dec(2,i)
         k1b(i)=mom_dec(1,i)
         EpsP(i)=hel_dec(3,i)
         EpsW(i)=hel_dec(2,i)
         SpB(i)=hel_dec(1,i)
      enddo
      EpsWk1b = scf(Dv,EpsW,k1b)
      EpsWp1t = scf(Dv,EpsW,p1t)
      EpsWEpsP = scf(Dv,EpsW,EpsP)
      EpsPk1b = scf(Dv,EpsP,k1b)
      EpsPp1t = scf(Dv,EpsP,p1t)
      EpsPk2w = scf(Dv,EpsP,k2w)
      k1bp1t = scf(Dv,k1b,p1t)
      k1bk2w = scf(Dv,k1b,k2w)
      k2wp1t = scf(Dv,k2w,p1t)
      t1 = EL**2
      t2 = alphaS*t1
      t3 = t2*EpsPk2w
      t4 = EpsWp1t*GW
      t5 = t4*mt
      t7 = 2*k1bk2w
      t8 = 2*k1bp1t
      t9 = mt**2
      t10 = Mw**2
      t14 = 1/(k1bk2w-k1bp1t)
      t17 = 2*mt*Mw
      t18 = t7-t8+t9-t17+t10
      t19 = 1/t18
      t20 = t7-t8+t9+t17+t10
      t21 = 1/t20
      t23 = 1/pii
      t24 = t19*t21*t23
      t25 = (t7-t8+t9-t10)*Qqd*t14*t24
      t28 = GW*mt
      t29 = EpsPk1b*EpsWk1b
      t30 = t10**2
      t33 = k1bk2w**2
      t34 = t33*k1bp1t
      t37 = EpsPk2w*EpsWk1b
      t40 = EpsPk2w*EpsWp1t
      t43 = k1bp1t**2
      t44 = k1bk2w*t43
      t51 = t33*t9
      t54 = t43*t9
      t59 = t9**2
      t60 = k1bk2w*t59
      t63 = k1bp1t*t59
      t68 = t43*t10
      t73 = EpsPk1b*EpsWp1t
      t78 = -t29*t9*t30+8*t29*t34+8*t37*t34-8*t40*t34-16*t29*t44-16*t37*
     #t44+8*t40*t44-4*t29*t51-12*t29*t54-8*t37*t54-4*t29*t60+6*t29*t63+2
     #*t37*t63-8*t29*t68-8*t37*t68+4*t73*t68+8*t40*t68
      t82 = k1bp1t*t30
      t89 = k1bk2w*k1bp1t
      t90 = t89*t10
      t93 = k1bk2w*t9
      t97 = k1bp1t*t9
      t98 = t97*t10
      t103 = t89*t9
      t116 = t59*t9
      t118 = t43*k1bp1t
      t123 = 2*t29*t59*t10+2*t29*t82+2*t37*t82+2*t73*t82-4*t40*t90-4*t29
     #*t93*t10-4*t37*t98-2*t73*t98+16*t29*t103+8*t37*t103-4*t40*t103+8*t
     #29*t90+8*t37*t90+4*t73*t90-t29*t116+8*t37*t118+8*t29*t118
      t127 = 1/k1bp1t
      t133 = k1bp1t*t10
      t135 = 1/(2*t89-t93+t133)
      t143 = Qqd-Qqu
      t144 = mt*t143
      t145 = -t8+t9-t10
      t161 = -1/t145
      t169 = mt*Qqu
      t170 = 1/k1bk2w
      t180 = GW*Qqd*t14*t23
      t190 = GW*t143*t161*t23
      t195 = t7-t9+t10
      t196 = 1/t195
      t198 = GW*Qqu*t196*t23
      t203 = t2*EpsWEpsP
      t206 = k1bp1t*t116
      t209 = t118*t10
      t221 = t43*t30
      t231 = t30*t10
      t232 = k1bp1t*t231
      t239 = t33*k1bk2w
      t240 = t239*k1bp1t
      t245 = t33*t43
      t252 = 2*t37*t206+8*t29*t209+16*t37*t209-8*t73*t209-16*t40*t209+3*
     #t29*t116*t10+4*t29*t221-4*t37*t221+8*t40*t221-3*t29*t59*t30-2*t29*
     #t232+2*t73*t232+t29*t9*t231+16*t37*t240-16*t40*t240-16*t29*t245-48
     #*t37*t245+32*t40*t245
      t253 = k1bk2w*t118
      t260 = t239*t9
      t265 = t118*t9
      t270 = t33*t59
      t277 = t43*t59
      t282 = k1bk2w*t116
      t298 = t93*t30
      t301 = 32*t29*t253+48*t37*t253-16*t40*t253-8*t37*t260+8*t40*t260+3
     #2*t29*t265+24*t37*t265-4*t29*t270-8*t37*t270+4*t40*t270-24*t29*t27
     #7-12*t37*t277-4*t29*t282-2*t37*t282+8*t29*t206+8*t37*k1bk2w*t98-12
     #*t73*k1bk2w*t98-20*t40*k1bk2w*t98-2*t37*t298
      t305 = t89*t30
      t314 = t63*t10
      t319 = t60*t10
      t324 = t54*t10
      t331 = t51*t10
      t340 = t34*t10
      t345 = 4*t29*t298-8*t40*t305+4*t73*t305+4*t37*t305-8*t29*t305-2*t7
     #3*t314-4*t37*t314+4*t73*t319-10*t29*t314+8*t73*t324+8*t40*t324+4*t
     #37*t319+4*t29*t331-8*t37*t331+12*t40*t331+4*t29*t324-8*t29*t340+16
     #*t37*t340
      t348 = t44*t10
      t355 = t34*t9
      t360 = t44*t9
      t367 = t89*t59
      t376 = t59**2
      t378 = t43**2
      t385 = t97*t30
      t390 = -24*t40*t340-32*t37*t348+8*t73*t348+40*t40*t348+40*t37*t355
     #-24*t40*t355-48*t29*t360-56*t37*t360+16*t40*t360+24*t29*t367+20*t3
     #7*t367-4*t40*t367+16*t29*t355-t29*t376-16*t29*t378-16*t37*t378-4*t
     #73*t298+4*t29*t385+2*t37*t385
      t419 = k1bk2w*t10
      t424 = t9*t10
      t428 = 4*t37*t33-4*t40*t33-4*t37*t89+4*t73*t89+4*t40*t89-2*t73*t93
     #-2*t29*t97+t29*t59+2*t73*t419+2*t29*t133-2*t29*t424+t29*t30
      t455 = Qqd*t14
      t463 = t2*EpsWk1b*GW
      t467 = t195**2
      t470 = (t7-3*t9+t10)*Qqu/t467*t23
      t478 = t170*t196*t23
      t497 = t2*GW
      DiagRat(1) = ne*(t3*t5*t25/2-t2*t28*(t78+t123)*Qqd*t127/t20/t18*t1
     #35*t23/4+t2*EpsWEpsP*GW*t144/t145*t23/2+t2*GW*(2*t29+2*t37-2*t73-2
     #*EpsWEpsP*k1bk2w+3*EpsWEpsP*k1bp1t)*t144*t127*t161*t23/4+t2*t29*GW
     #*t169*t170*t127*t23/4)
      DiagRat(2) = ne*(-t3*t180/4-t2*(EpsPk1b+2*EpsPk2w)*t180/4-t3*t190+
     #t2*(EpsPk1b+EpsPk2w)*t198/2)
      DiagRat(3) = ne*(3.D0/4.D0*t203*t180-t2*GW*(t252+t301+t345+t390)*Q
     #qd*t14*t161*t19*t21*t135*t23/4+t203*t190-t2*GW*t428*Qqu*t170*t14*t
     #196*t161*t23/4)
      DiagRat(4) = ne*(t2*EpsWp1t*t180/4+t2*t4*(4*t33-8*t89+4*t43+4*t93-
     #4*t97+t59+6*t419-6*t133-3*t424+2*t30)*t455*t24/4-t2*(EpsWk1b-EpsWp
     #1t)*t190+t463*t470/2-t2*EpsWk1b*t198-t463*(t7-t10)*Qqu*t478/4)
      DiagRat(5) = (0.d0,0.d0)
      DiagRat(6) = -t2*t28*ne*Qqu/t195*t23/4
      DiagRat(7) = ne*(-t2*t5*t25+t463*t169*t478)/4
      DiagRat(8) = ne*(3.D0/8.D0*t497*t455*t23+t497*t470/4-t497*Qqu*t196
     #*t23)

      end subroutine numTopDecRat

      subroutine numTopDecCTmass(MassCT)
      implicit none
      include 'VarTopMassCT.f'
      include 'commondecay.f'
      include 'intcommon.f'
      double precision wrat
      double complex p1t(4),k2w(4),k1b(4), k3p(4)
      double complex EpsP(4),EpsW(4),SpB(4),SpT(4)

      double complex k1bp1t, k2wp1t, k1bk2w
      double complex EpsWk1b, EpsWp1t, EpsWEpsP
      double complex EpsPk1b, EpsPp1t, EpsPk2w
      double complex MassCT(8)
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
      double complex s11,s12,s13,s14,s15,s16,s17,s18,s19
      integer i, Dv
      Dv = 4
      do i = 1,4
         p1t(i)=mom_t(i)
         k2w(i)=mom_dec(2,i)
         k1b(i)=mom_dec(1,i)
         EpsP(i)=hel_dec(3,i)
         EpsW(i)=hel_dec(2,i)
         SpB(i)=hel_dec(1,i)
      enddo
      EpsWk1b = scf(Dv,EpsW,k1b)
      EpsWp1t = scf(Dv,EpsW,p1t)
      EpsWEpsP = scf(Dv,EpsW,EpsP)
      EpsPk1b = scf(Dv,EpsP,k1b)
      EpsPp1t = scf(Dv,EpsP,p1t)
      EpsPk2w = scf(Dv,EpsP,k2w)
      k1bp1t = scf(Dv,k1b,p1t)
      k1bk2w = scf(Dv,k1b,k2w)
      k2wp1t = scf(Dv,k2w,p1t)
      t1 = EL**2
      t7 = mt**2
      t8 = Mw**2
      t17 = -2*k1bk2w+t7-t8
      t18 = t17**2
      t19 = 1/t18
      t24 = t1*GW
      t25 = ne*Qqu
      MassCT(1) = t1*EpsWEpsP*GW*ne*(Qqd-Qqu)/(2*k1bp1t-t7+t8)
      MassCT(2) = (0.d0,0.d0)
      MassCT(3) = (0.d0,0.d0)
      MassCT(4) = 4*t1*EpsWk1b*GW*mt*ne*Qqu*t19
      MassCT(5) = (0.d0,0.d0)
      MassCT(6) = -t24*t25/t17
      MassCT(7) = (0.d0,0.d0)
      MassCT(8) = 2*t24*mt*t25*t19

      end subroutine numTopDecCTmass

      subroutine numTopDecCTwave(WaveCT)
      implicit none
      include 'VarTopWaveCT.f'
      include 'commondecay.f'
      include 'intcommon.f'
      double precision wrat
      double complex p1t(4),k2w(4),k1b(4), k3p(4)
      double complex EpsP(4),EpsW(4),SpB(4),SpT(4)

      double complex k1bp1t, k2wp1t, k1bk2w
      double complex EpsWk1b, EpsWp1t, EpsWEpsP
      double complex EpsPk1b, EpsPp1t, EpsPk2w
      double complex WaveCT(8)
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
      double complex s11,s12,s13,s14,s15,s16,s17,s18,s19
      integer i, Dv
      Dv = 4
      do i = 1,4
         p1t(i)=mom_t(i)
         k2w(i)=mom_dec(2,i)
         k1b(i)=mom_dec(1,i)
         EpsP(i)=hel_dec(3,i)
         EpsW(i)=hel_dec(2,i)
         SpB(i)=hel_dec(1,i)
      enddo
      EpsWk1b = scf(Dv,EpsW,k1b)
      EpsWp1t = scf(Dv,EpsW,p1t)
      EpsWEpsP = scf(Dv,EpsW,EpsP)
      EpsPk1b = scf(Dv,EpsP,k1b)
      EpsPp1t = scf(Dv,EpsP,p1t)
      EpsPk2w = scf(Dv,EpsP,k2w)
      k1bp1t = scf(Dv,k1b,p1t)
      k1bk2w = scf(Dv,k1b,k2w)
      k2wp1t = scf(Dv,k2w,p1t)
      t1 = EL**2
      t6 = mt**2
      t8 = Mw**2
      t14 = -2*k1bk2w*Qqd+t6*Qqd-t8*Qqd+2*k1bk2w*Qqu-2*k1bp1t*Qqu
      t17 = 1/(k1bk2w-k1bp1t)
      t20 = 1/(-2*k1bp1t+t6-t8)
      t22 = ne*t14*t17*t20
      t29 = t1*GW
      t41 = t14*t17
      t44 = 1/(-2*k1bk2w+t6-t8)
      WaveCT(1) = (0.d0,0.d0)
      WaveCT(2) = t1*EpsPk2w*GW*t22/2
      WaveCT(3) = -t1*EpsWEpsP*GW*t22/2
      WaveCT(4) = -t29*(2*EpsWk1b*k1bk2w-2*EpsWp1t*k1bk2w-2*EpsWk1b*k1bp
     #1t+EpsWp1t*t6-EpsWp1t*t8)*ne*t41*t44*t20/2
      WaveCT(5) = (0.d0,0.d0)
      WaveCT(6) = (0.d0,0.d0)
      WaveCT(7) = (0.d0,0.d0)
      WaveCT(8) = -t29*ne*t41*t44/4

      end subroutine numTopDecCTwave















      subroutine WeyltoDiracS(sp,sp2)
      implicit none
      integer Dv, Ds
      double complex sp(4),sp2(4), sp1(4)
      sp1(1:4) = sp(1:4)
      sp2= WeylToDirac(sp1)
      return
      end subroutine WeyltoDiracS



      subroutine vbqqweylS(Dv,Ds,sp,sp2,Weyl)
      implicit none
      integer Dv, Ds
      real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0
      double complex sp(Ds),sp2(Ds), Weyl(Dv)
      Weyl = vbqq_Weyl(Dv,sp,sp2)*sqrt2/(0.d0,-1.d0)

      return
      end subroutine vbqqweylS



      subroutine spb2_weylS(Dv,Ds,sp,v,sp2)
      implicit none
      integer Dv, Ds
      double complex sp(Ds),sp2(Ds),v(Dv)
      sp2= spb2_Weyl(sp,v)
      return
      end subroutine spb2_weylS


      subroutine spi2_weylS(Dv,Ds,sp,v,sp2)
      implicit none
      integer Dv, Ds
      double complex sp(Ds),sp2(Ds),v(Dv)
      sp2= spi2_Weyl(v,sp)
      return
      end subroutine spi2_weylS


*******************************************************************************************************************************
*******************************************************************************************************************************
********************* Standard Matrix-elements for TOP-decay ******************************************************************
*******************************************************************************************************************************
*******************************************************************************************************************************

      subroutine calc_ME1(SME1)
      implicit none
      double complex SME1(4)
      integer Dv,Ds,i
      include 'commondecay.f'
      double complex sp2(4), slash(4)
      double complex sp2i(4),sp2a(4)
C Bare bottom-spinor is multiplied with (p_slash + mt)

      Dv = 4
      Ds = 4
      do i = 1,4
         sp2i(i)=hel_dec(1,i)
         slash(i) = mom_t(i)
      enddo
C      sp2i(1) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(2) = 1.d0/2.d0*(sp2(2)+sp2(4))
C      sp2i(3) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(4) = 1.d0/2.d0*(sp2(2)+sp2(4))

      call spb2_weylS(Dv,Ds,sp2i,slash,sp2a)

      do i=1,4
         SME1(i) = sp2a(i)+mt*sp2i(i)
      enddo
      return
      end subroutine calc_ME1


      subroutine calc_ME2(SME2)
      implicit none
      double complex SME2(4)
      integer Dv,Ds,i
      include 'commondecay.f'
      double complex sp2(4), slash(4), EpsW(4)
      double complex sp2i(4),sp2a(4)
C At the end bottom-spinor is multiplied with (p_slash + mt)
      Dv = 4
      Ds = 4
      do i = 1,4
         sp2i(i)=hel_dec(1,i)
         slash(i) = mom_t(i)
         EpsW(i) = hel_dec(2,i)
      enddo
C      sp2i(1) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(2) = 1.d0/2.d0*(sp2(2)+sp2(4))
C      sp2i(3) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(4) = 1.d0/2.d0*(sp2(2)+sp2(4))

      call spb2_weylS(Dv,Ds,sp2i,EpsW,sp2a)
      call spb2_weylS(Dv,Ds,sp2a,slash,sp2)


      do i=1,4
         SME2(i) = sp2(i)+mt*sp2a(i)
      enddo
      return
      end subroutine calc_ME2


      subroutine calc_ME3(SME3)
      implicit none
      double complex SME3(4)
      integer Dv,Ds,i
      include 'commondecay.f'
      double complex sp2(4), slash(4), k2w(4)
      double complex sp2i(4),sp2a(4)
      double complex sp2a2(4)
C At the end bottom-spinor is multiplied with (p_slash + mt)
      Dv = 4
      Ds = 4
      do i = 1,4
         sp2i(i) = hel_dec(1,i)
         k2w(i) = mom_dec(2,i)
         slash(i) = mom_t(i)
      enddo
C      sp2i(1) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(2) = 1.d0/2.d0*(sp2(2)+sp2(4))
C      sp2i(3) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(4) = 1.d0/2.d0*(sp2(2)+sp2(4))
      call spb2_weylS(Dv,Ds,sp2i,k2w,sp2a)
      call spb2_weylS(Dv,Ds,sp2a,slash,sp2)

      call spb2_weylS(Dv,Ds,sp2i,slash,sp2a2)


      do i=1,4
         SME3(i) = sp2(i)+mt*sp2a(i)-mt*(sp2a2(i)+mt*sp2i(i))
      enddo
      return
      end subroutine calc_ME3


C      subroutine calc_ME3(SME3)
C      implicit none
C      double complex SME3(4)
C      integer Dv,Ds,i
C      include 'commondecay.f'
C      double complex sp2(4), slash(4), k2w(4)
C      double complex sp2i(4),sp2a(4)
C      double complex sp2a2(4)
CC At the end bottom-spinor is multiplied with (p_slash + mt)
C      Dv = 4
C      Ds = 4
C      do i = 1,4
C         sp2i(i) = hel_dec(1,i)
C         k2w(i) = mom_dec(2,i)
C         slash(i) = mom_t(i)
C      enddo
CC      sp2i(1) = 1.d0/2.d0*(sp2(1)+sp2(3))
CC      sp2i(2) = 1.d0/2.d0*(sp2(2)+sp2(4))
CC      sp2i(3) = 1.d0/2.d0*(sp2(1)+sp2(3))
CC      sp2i(4) = 1.d0/2.d0*(sp2(2)+sp2(4))
C      call spb2_weylS(Dv,Ds,sp2i,k2w,sp2a)
C      call spb2_weylS(Dv,Ds,sp2a,slash,sp2)
C
C      do i=1,4
C         SME3(i) = sp2(i)+mt*sp2a(i)
C      enddo
C      return
C      end


      subroutine calc_ME4(SME4)
      implicit none
      double complex SME4(4)
      integer Dv,Ds,i
      include 'commondecay.f'
      double complex sp2(4), slash(4), EpsP(4)
      double complex sp2i(4),sp2a(4)
C At the end bottom-spinor is multiplied with (p_slash + mt)
      Dv = 4
      Ds = 4
      do i = 1,4
         sp2i(i) = hel_dec(1,i)
         EpsP(i) = hel_dec(3,i)
         slash(i) = mom_t(i)
      enddo
C      sp2i(1) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(2) = 1.d0/2.d0*(sp2(2)+sp2(4))
C      sp2i(3) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(4) = 1.d0/2.d0*(sp2(2)+sp2(4))
      call spb2_weylS(Dv,Ds,sp2i,EpsP,sp2a)
      call spb2_weylS(Dv,Ds,sp2a,slash,sp2)
      do i=1,4
         SME4(i) = sp2(i)+mt*sp2a(i)
      enddo
      return
      end subroutine calc_ME4



      subroutine calc_ME5(SME5)
      implicit none
      double complex SME5(4)
      integer Dv,Ds,i
      include 'commondecay.f'
      double complex sp2(4), slash(4), EpsW(4), k2w(4)
      double complex sp2i(4),sp2a(4),sp2b(4)
C At the end bottom-spinor is multiplied with (p_slash + mt)
      Dv = 4
      Ds = 4
      do i = 1,4
         sp2i(i) = hel_dec(1,i)
         EpsW(i) = hel_dec(2,i)
         k2w(i) = mom_dec(2,i)
         slash(i) = mom_t(i)
      enddo
C      sp2i(1) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(2) = 1.d0/2.d0*(sp2(2)+sp2(4))
C      sp2i(3) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(4) = 1.d0/2.d0*(sp2(2)+sp2(4))
      call spb2_weylS(Dv,Ds,sp2i,EpsW,sp2a)
      call spb2_weylS(Dv,Ds,sp2a,k2w,sp2b)
      call spb2_weylS(Dv,Ds,sp2b,slash,sp2)

      do i=1,4
         SME5(i) = sp2(i)+mt*sp2b(i)
      enddo
      return
      end subroutine calc_ME5




      subroutine calc_ME6(SME6)
      implicit none
      double complex SME6(4)
      integer Dv,Ds,i
      include 'commondecay.f'
      double complex sp2(4), slash(4), EpsW(4), EpsP(4)
      double complex sp2i(4),sp2a(4),sp2b(4)
C At the end bottom-spinor is multiplied with (p_slash + mt)
      Dv = 4
      Ds = 4
      do i = 1,4
         sp2i(i) = hel_dec(1,i)
         EpsW(i) = hel_dec(2,i)
         EpsP(i) = hel_dec(3,i)
         slash(i) = mom_t(i)
      enddo
C      sp2i(1) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(2) = 1.d0/2.d0*(sp2(2)+sp2(4))
C      sp2i(3) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(4) = 1.d0/2.d0*(sp2(2)+sp2(4))
      call spb2_weylS(Dv,Ds,sp2i,EpsW,sp2a)
      call spb2_weylS(Dv,Ds,sp2a,EpsP,sp2b)
      call spb2_weylS(Dv,Ds,sp2b,slash,sp2)

      do i=1,4
         SME6(i) = sp2(i)+mt*sp2b(i)
      enddo
      return
      end subroutine calc_ME6



CSMElements = {DecP0,DecP[EpsW],DecP[k2w-Mtop],DecP[EpsP],DecP[EpsW,k2w],DecP[EpsW,EpsP],DecP[k2w,EpsP],DecP[EpsW,k2w+Mtop,EpsP]}
      subroutine calc_ME7(SME7)
      implicit none
      double complex SME7(4)
      integer Dv,Ds,i
      include 'commondecay.f'
      double complex sp2(4), slash(4), k2w(4), EpsP(4)
      double complex sp2i(4),sp2a(4),sp2b(4)
C At the end bottom-spinor is multiplied with (p_slash + mt)
      Dv = 4
      Ds = 4
      do i = 1,4
         sp2i(i) = hel_dec(1,i)
         k2w(i) = mom_dec(2,i)
         EpsP(i) = hel_dec(3,i)
         slash(i) = mom_t(i)
      enddo
C      sp2i(1) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(2) = 1.d0/2.d0*(sp2(2)+sp2(4))
C      sp2i(3) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(4) = 1.d0/2.d0*(sp2(2)+sp2(4))
      call spb2_weylS(Dv,Ds,sp2i,k2w,sp2a)
      call spb2_weylS(Dv,Ds,sp2a,EpsP,sp2b)
      call spb2_weylS(Dv,Ds,sp2b,slash,sp2)

      do i=1,4
         SME7(i) = sp2(i)+mt*sp2b(i)
      enddo
      return
      end subroutine calc_ME7



      subroutine calc_ME8(SME8)
      implicit none
      double complex SME8(4)
      integer Dv,Ds,i
      include 'commondecay.f'
      double complex sp2(4), slash(4), k2w(4), EpsP(4), EpsW(4)
      double complex sp2i(4),sp2a(4),sp2b(4),sp2c(4)
      double complex sp2a2(4),sp2b2(4),sp2c2(4)
C At the end bottom-spinor is multiplied with (p_slash + mt)
      Dv = 4
      Ds = 4
      do i = 1,4
         sp2i(i) = hel_dec(1,i)
         k2w(i) = mom_dec(2,i)
         EpsP(i) = hel_dec(3,i)
         EpsW(i) = hel_dec(2,i)
         slash(i) = mom_t(i)
      enddo
C      sp2i(1) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(2) = 1.d0/2.d0*(sp2(2)+sp2(4))
C      sp2i(3) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(4) = 1.d0/2.d0*(sp2(2)+sp2(4))
      call spb2_weylS(Dv,Ds,sp2i,EpsW,sp2a)
      call spb2_weylS(Dv,Ds,sp2a,k2w,sp2b)
      call spb2_weylS(Dv,Ds,sp2b,EpsP,sp2c)
      call spb2_weylS(Dv,Ds,sp2c,slash,sp2)


      call spb2_weylS(Dv,Ds,sp2i,EpsW,sp2a2)
      call spb2_weylS(Dv,Ds,sp2a2,EpsP,sp2b2)
      call spb2_weylS(Dv,Ds,sp2b2,slash,sp2c2)

      do i=1,4
         SME8(i) = sp2(i)+mt*sp2c(i) + mt*(sp2c2(i) + mt*sp2b2(i))
      enddo
      return
      end subroutine calc_ME8

C      subroutine calc_ME8(SME8)
C      implicit none
C      double complex SME8(4)
C      integer Dv,Ds,i
C      include 'commondecay.f'
C      double complex sp2(4), slash(4), k2w(4), EpsP(4), EpsW(4)
C      double complex sp2i(4),sp2a(4),sp2b(4),sp2c(4)
C      double complex sp2a2(4),sp2b2(4),sp2c2(4)
CC At the end bottom-spinor is multiplied with (p_slash + mt)
C      Dv = 4
C      Ds = 4
C      do i = 1,4
C         sp2i(i) = hel_dec(1,i)
C         k2w(i) = mom_dec(2,i)
C         EpsP(i) = hel_dec(3,i)
C         EpsW(i) = hel_dec(2,i)
C         slash(i) = mom_t(i)
C      enddo
CC      sp2i(1) = 1.d0/2.d0*(sp2(1)+sp2(3))
CC      sp2i(2) = 1.d0/2.d0*(sp2(2)+sp2(4))
CC      sp2i(3) = 1.d0/2.d0*(sp2(1)+sp2(3))
CC      sp2i(4) = 1.d0/2.d0*(sp2(2)+sp2(4))
C      call spb2_weylS(Dv,Ds,sp2i,EpsW,sp2a)
C      call spb2_weylS(Dv,Ds,sp2a,k2w,sp2b)
C      call spb2_weylS(Dv,Ds,sp2b,EpsP,sp2c)
C      call spb2_weylS(Dv,Ds,sp2c,slash,sp2)
C
C
C
C      do i=1,4
C         SME8(i) = sp2(i)+mt*sp2c(i)
C      enddo
C      return
C      end


*******************************************************************************************************************************
*******************************************************************************************************************************
********************* Standard Matrix-elements for TBAR-decay *****************************************************************
*******************************************************************************************************************************
*******************************************************************************************************************************

C DecM0,DecM[EpsW],DecM[k2w-Mtop],DecM[EpsP],DecM[EpsW,k2w],DecM[EpsW,EpsP],DecM[k2w,EpsP],DecM[EpsW,k2w+Mtop,EpsP]
      subroutine calc_AME1(SAME1)
      implicit none
      double complex SAME1(4)
      integer Dv,Ds,i
      include 'commondecay.f'
      double complex sp2(4), slash(4)
      double complex sp2i(4),sp2a(4)
C Bare bottom-spinor is multiplied with (p_slash + mt)

      Dv = 4
      Ds = 4
      do i = 1,4
         sp2i(i)=hel_dec(1,i)
         slash(i) = -mom_t(i)
      enddo
C      sp2i(1) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(2) = 1.d0/2.d0*(sp2(2)+sp2(4))
C      sp2i(3) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(4) = 1.d0/2.d0*(sp2(2)+sp2(4))

      call spi2_weylS(Dv,Ds,sp2i,slash,sp2a)

      do i=1,4
         SAME1(i) = sp2a(i)+mt*sp2i(i)
      enddo
      return
      end subroutine calc_AME1


      subroutine calc_AME2(SAME2)
      implicit none
      double complex SAME2(4)
      integer Dv,Ds,i
      include 'commondecay.f'
      double complex sp2(4), slash(4), EpsW(4)
      double complex sp2i(4),sp2a(4), MA_weylS(4)
C At the end bottom-spinor is multiplied with (p_slash + mt)
      Dv = 4
      Ds = 4
      do i = 1,4
         sp2i(i)=hel_dec(1,i)
         slash(i) = -mom_t(i)
         EpsW(i) = hel_dec(2,i)
      enddo
C      sp2i(1) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(2) = 1.d0/2.d0*(sp2(2)+sp2(4))
C      sp2i(3) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(4) = 1.d0/2.d0*(sp2(2)+sp2(4))

      call spi2_weylS(Dv,Ds,sp2i,EpsW,sp2a)
      call spi2_weylS(Dv,Ds,sp2a,slash,sp2)


      do i=1,4
         SAME2(i) = sp2(i)+mt*sp2a(i)
      enddo
      return
      end subroutine calc_AME2


      subroutine calc_AME3(SAME3)
      implicit none
      double complex SAME3(4)
      integer Dv,Ds,i
      include 'commondecay.f'
      double complex sp2(4), slash(4), k2w(4)
      double complex sp2i(4),sp2a(4)
      double complex sp2a2(4)
C At the end bottom-spinor is multiplied with (p_slash + mt)
      Dv = 4
      Ds = 4
      do i = 1,4
         sp2i(i) = hel_dec(1,i)
         k2w(i) = mom_dec(2,i)
         slash(i) = -mom_t(i)
      enddo
C      sp2i(1) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(2) = 1.d0/2.d0*(sp2(2)+sp2(4))
C      sp2i(3) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(4) = 1.d0/2.d0*(sp2(2)+sp2(4))
      call spi2_weylS(Dv,Ds,sp2i,k2w,sp2a)
      call spi2_weylS(Dv,Ds,sp2a,slash,sp2)

      call spi2_weylS(Dv,Ds,sp2i,slash,sp2a2)


      do i=1,4
         SAME3(i) = sp2(i)+mt*sp2a(i)-mt*(sp2a2(i)+mt*sp2i(i))
      enddo
      return
      end subroutine calc_AME3

      subroutine calc_AME4(SAME4)
      implicit none
      double complex SAME4(4)
      integer Dv,Ds,i
      include 'commondecay.f'
      double complex sp2(4), slash(4), EpsP(4)
      double complex sp2i(4),sp2a(4)
C At the end bottom-spinor is multiplied with (p_slash + mt)
      Dv = 4
      Ds = 4
      do i = 1,4
         sp2i(i) = hel_dec(1,i)
         EpsP(i) = hel_dec(3,i)
         slash(i) = -mom_t(i)
      enddo
C      sp2i(1) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(2) = 1.d0/2.d0*(sp2(2)+sp2(4))
C      sp2i(3) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(4) = 1.d0/2.d0*(sp2(2)+sp2(4))
      call spi2_weylS(Dv,Ds,sp2i,EpsP,sp2a)
      call spi2_weylS(Dv,Ds,sp2a,slash,sp2)
      do i=1,4
         SAME4(i) = sp2(i)+mt*sp2a(i)
      enddo
      return
      end subroutine calc_AME4

C DecM0,DecM[EpsW],DecM[k2w-Mtop],DecM[EpsP],DecM[EpsW,k2w],DecM[EpsW,EpsP],DecM[k2w,EpsP],DecM[EpsW,k2w+Mtop,EpsP]

      subroutine calc_AME5(SAME5)
      implicit none
      double complex SAME5(4)
      integer Dv,Ds,i
      include 'commondecay.f'
      double complex sp2(4), slash(4), EpsW(4), k2w(4)
      double complex sp2i(4),sp2a(4),sp2b(4)
C At the end bottom-spinor is multiplied with (p_slash + mt)
      Dv = 4
      Ds = 4
      do i = 1,4
         sp2i(i) = hel_dec(1,i)
         EpsW(i) = hel_dec(2,i)
         k2w(i) = mom_dec(2,i)
         slash(i) = -mom_t(i)
      enddo
C      sp2i(1) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(2) = 1.d0/2.d0*(sp2(2)+sp2(4))
C      sp2i(3) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(4) = 1.d0/2.d0*(sp2(2)+sp2(4))
      call spi2_weylS(Dv,Ds,sp2i,k2w,sp2a)
      call spi2_weylS(Dv,Ds,sp2a,EpsW,sp2b)
      call spi2_weylS(Dv,Ds,sp2b,slash,sp2)

      do i=1,4
         SAME5(i) = sp2(i)+mt*sp2b(i)
      enddo
      return
      end subroutine calc_AME5




      subroutine calc_AME6(SAME6)
      implicit none
      double complex SAME6(4)
      integer Dv,Ds,i
      include 'commondecay.f'
      double complex sp2(4), slash(4), EpsW(4), EpsP(4)
      double complex sp2i(4),sp2a(4),sp2b(4)
C At the end bottom-spinor is multiplied with (p_slash + mt)
      Dv = 4
      Ds = 4
      do i = 1,4
         sp2i(i) = hel_dec(1,i)
         EpsW(i) = hel_dec(2,i)
         EpsP(i) = hel_dec(3,i)
         slash(i) = -mom_t(i)
      enddo
C      sp2i(1) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(2) = 1.d0/2.d0*(sp2(2)+sp2(4))
C      sp2i(3) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(4) = 1.d0/2.d0*(sp2(2)+sp2(4))
      call spi2_weylS(Dv,Ds,sp2i,EpsP,sp2a)
      call spi2_weylS(Dv,Ds,sp2a,EpsW,sp2b)
      call spi2_weylS(Dv,Ds,sp2b,slash,sp2)

      do i=1,4
         SAME6(i) = sp2(i)+mt*sp2b(i)
      enddo
      return
      end subroutine calc_AME6



CSAMElements = {DecM0,DecM[EpsW],DecM[k2w-Mtop],DecM[EpsP],DecM[EpsW,k2w],DecM[EpsW,EpsP],DecM[k2w,EpsP],DecM[EpsW,k2w+Mtop,EpsP]}
      subroutine calc_AME7(SAME7)
      implicit none
      double complex SAME7(4)
      integer Dv,Ds,i
      include 'commondecay.f'
      double complex sp2(4), slash(4), k2w(4), EpsP(4)
      double complex sp2i(4),sp2a(4),sp2b(4)
C At the end bottom-spinor is multiplied with (p_slash + mt)
      Dv = 4
      Ds = 4
      do i = 1,4
         sp2i(i) = hel_dec(1,i)
         k2w(i) = mom_dec(2,i)
         EpsP(i) = hel_dec(3,i)
         slash(i) = -mom_t(i)
      enddo
C      sp2i(1) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(2) = 1.d0/2.d0*(sp2(2)+sp2(4))
C      sp2i(3) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(4) = 1.d0/2.d0*(sp2(2)+sp2(4))
      call spi2_weylS(Dv,Ds,sp2i,EpsP,sp2a)
      call spi2_weylS(Dv,Ds,sp2a,k2w,sp2b)
      call spi2_weylS(Dv,Ds,sp2b,slash,sp2)

      do i=1,4
         SAME7(i) = sp2(i)+mt*sp2b(i)
      enddo
      return
      end subroutine calc_AME7


CSAMElements = {DecM0,DecM[EpsW],DecM[k2w-Mtop],DecM[EpsP],DecM[EpsW,k2w],DecM[EpsW,EpsP],DecM[k2w,EpsP],DecM[EpsW,k2w+Mtop,EpsP]}
      subroutine calc_AME8(SAME8)
      implicit none
      double complex SAME8(4)
      integer Dv,Ds,i
      include 'commondecay.f'
      double complex sp2(4), slash(4), k2w(4), EpsP(4), EpsW(4)
      double complex sp2i(4),sp2a(4),sp2b(4),sp2c(4)
      double complex sp2a2(4),sp2b2(4),sp2c2(4)
C At the end bottom-spinor is multiplied with (p_slash + mt)
      Dv = 4
      Ds = 4
      do i = 1,4
         sp2i(i) = hel_dec(1,i)
         k2w(i) = mom_dec(2,i)
         EpsP(i) = hel_dec(3,i)
         EpsW(i) = hel_dec(2,i)
         slash(i) = -mom_t(i)
      enddo
C      sp2i(1) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(2) = 1.d0/2.d0*(sp2(2)+sp2(4))
C      sp2i(3) = 1.d0/2.d0*(sp2(1)+sp2(3))
C      sp2i(4) = 1.d0/2.d0*(sp2(2)+sp2(4))
      call spi2_weylS(Dv,Ds,sp2i,EpsP,sp2a)
      call spi2_weylS(Dv,Ds,sp2a,k2w,sp2b)
      call spi2_weylS(Dv,Ds,sp2b,EpsW,sp2c)
      call spi2_weylS(Dv,Ds,sp2c,slash,sp2)


      call spi2_weylS(Dv,Ds,sp2i,EpsP,sp2a2)
      call spi2_weylS(Dv,Ds,sp2a2,EpsW,sp2b2)
      call spi2_weylS(Dv,Ds,sp2b2,slash,sp2c2)

      do i=1,4
         SAME8(i) = sp2(i)+mt*sp2c(i) + mt*(sp2c2(i) + mt*sp2b2(i))
      enddo
      return
      end subroutine calc_AME8






!-------massive vector boson polarization routine, for off-shell vector bosons decaying into Neutrino and Lepton
      subroutine pol_massoff(p_lep,p_neu,EpsW_weyl)
      implicit none
      include 'commondecay.f'
      double complex p_Lep(4), p_neu(4)
      double complex Wpol(4)
      double complex s_n(4), s_l(4)
      double complex  Weyl(4), EpsW_weyl(4)
      integer Dv,Ds,i, hl, hn
      real(8), parameter :: sqrt2 = 1.41421356237309504880168872420969d0
      Dv = 4
      Ds = 4
      do i=1,4
         EpsW_weyl(i) = (0.d0,0.d0)
         Weyl(i) = (0.d0,0.d0)
      enddo
      call ubarSpi_Weyl(p_neu,-1,s_n)
      call vSpi_Weyl(p_lep,+1,s_l)
      call vbqqweylS(Dv,Ds,s_n,s_l,Weyl)
      EpsW_weyl(1:4) =  Weyl(1:4)*(0.d0,1.d0)*EL*GW
      end subroutine pol_massoff


      subroutine pol_massoff2(p_lep,p_neu,EpsW_weyl)
      implicit none
      include 'commondecay.f'
      double complex p_Lep(4), p_neu(4)
      double complex Wpol(4)
      double complex s_n(4), s_l(4)
      double complex  Weyl(4), EpsW_weyl(4)
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
      call vbqqweylS(Dv,Ds,s_l,s_n,Weyl)
      EpsW_weyl(1:4) =  Weyl(1:4)*(0.d0,1.d0)*EL*GW
      end subroutine pol_massoff2


      double complex function scf(Dv,ap1,ap2)
        implicit none
        integer Dv
        double complex ap1(Dv),ap2(Dv)
        double complex r1

        call sc(Dv,ap1,ap2,r1)
        scf = r1
        return
      end function scf

C      subroutine sc(n,x,y,r)
C        implicit none
C        integer i,n
C        double complex x(*),y(*)
C        double complex r
C
C        r=dcmplx(0d0,0d0)
C
C        do i=1, n
C           if (i.eq.1) then
C              r = r + x(i)*y(i)
C           else
C              r = r - x(i)*y(i)
C           endif
C        enddo
C
C        return
C      end subroutine sc




      double complex function dilog(u)
      implicit none
      double complex  w, z, cdil,u
      double precision re,pi, be
      parameter(pi = 0.3141592653589793D1)
      integer i
      re = dreal(u)
      be = cdabs(u)
      if( be .lt. 1.d0) then
         if (re .gt. 0.5d0) then
            dilog = dcmplx(pi**2/6.d0,0.d0)
     &       -log(u)*log(1.d0-u)-bdil(1.d0-u)
         else
           dilog = bdil(u)
         endif
      else
         if (be .lt. 2.d0) then
            dilog = dcmplx(pi**2/6.d0,0.d0)
     &           -log(u)*log(1.d0-u)-bdil(1.d0-u)
         else
            dilog = dcmplx(-pi**2/6.d0,0.d0)-
     &        1.d0/2.d0*log(-u)**2
     &        -bdil(1.d0/u)
         endif
      endif
      return
      end function dilog









      double complex function bdil(x)
      implicit none
      double complex z, cdil,x
      double precision  r,phi,v
      integer i
      cdil = 0
      z = -log(1.d0-x)
      do i = 0,30
         cdil = cdil + bernuo(i) * z**(i+1)/fak(i+1)
      enddo
      bdil = cdil
      return
      end function bdil


      double precision function fak(j)
      implicit none
      double precision  efak
      integer i,j
      efak =1.d0
      do i=1,j
         efak= efak*dble(i)
      enddo
      fak = efak
      return
      end function fak

      double precision function bernuo(i)
      implicit none
      double precision  bern(0:30)
      integer i
      bern(0) = 1.d0
      bern(1) = -1.d0/2.d0
      bern(2) = 1.d0/6.d0
      bern(3) = 0.d0
      bern(4) = -1.d0/30.d0
      bern(5) = 0.d0
      bern(6) = 1.d0/42.d0
      bern(7) = 0.d0
      bern(8) = -1.d0/30.d0
      bern(9) = 0.d0
      bern(10) = 5.d0/66.d0
      bern(11) = 0.d0
      bern(12) = -691.d0/2730.d0
      bern(13) = 0.d0
      bern(14) = 7.d0/6.d0
      bern(15) = 0.d0
      bern(16) = -3617.d0/510.d0
      bern(17) = 0.d0
      bern(18) = 43867.d0/798.d0
      bern(19) = 0.d0
      bern(20) = -174611.d0/330.d0
      bern(21) = 0.d0
      bern(22) = 854513.d0/138.d0
      bern(23) = 0.d0
      bern(24) = -236364091.d0/2730.d0
      bern(25) = 0.d0
      bern(26) = 8553103.d0/6.d0
      bern(27) = 0.d0
      bern(28) = -23749461029.d0/870.d0
      bern(29) = 0.d0
      bern(30) = 8615841276005.d0/14322.d0

      bernuo = bern(i)
      return
      end function bernuo



      subroutine numTopDectree(res)
      implicit none
      include 'VarTopLO.f'
      include 'commondecay.f'
      double complex p1t(4),k2w(4),k1b(4)
      double complex EpsP(4),EpsW(4),SpB(4),SpT(4)
      double complex res(4)
      double complex scalar1,scalar2,scalar3,scalar4
      double complex scalar5,scalar6,scalar7,scalar8
      double complex ME1(4), ME2(4), ME3(4), ME4(4)
      double complex ME5(4), ME6(4), ME7(4), ME8(4)
      double complex k1bp1t, k2wp1t, k1bk2w, Diags(8)
      double complex EpsWk1b, EpsWp1t, EpsWEpsP
      double complex EpsPk1b, EpsPp1t, EpsPk2w
      integer i, Dv
      Dv = 4
      do i = 1,4
         p1t(i)=mom_t(i)
         k2w(i)=mom_dec(2,i)
         k1b(i)=mom_dec(1,i)
         EpsP(i)=hel_dec(3,i)
         EpsW(i)=hel_dec(2,i)
         SpB(i)=hel_dec(1,i)
	 res(i) = (0.d0,0.d0)
      enddo
      call Calc_ME1(ME1)
      call Calc_ME2(ME2)
      call Calc_ME3(ME3)
      call Calc_ME4(ME4)
      call Calc_ME5(ME5)
      call Calc_ME6(ME6)
      call Calc_ME7(ME7)
      call Calc_ME8(ME8)
      EpsWk1b = scf(Dv,EpsW,k1b)
      EpsWp1t = scf(Dv,EpsW,p1t)
      EpsWEpsP = scf(Dv,EpsW,EpsP)
      EpsPk1b = scf(Dv,EpsP,k1b)
      EpsPp1t = scf(Dv,EpsP,p1t)
      EpsPk2w = scf(Dv,EpsP,k2w)
      k1bp1t = scf(Dv,k1b,p1t)
      k1bk2w = scf(Dv,k1b,k2w)
      k2wp1t = scf(Dv,k2w,p1t)
      t1 = EL**2
      t6 = mt**2
      t8 = Mw**2
      t14 = 2*k1bk2w*Qqd-t6*Qqd+t8*Qqd-2*k1bk2w*Qqu+2*k1bp1t*Qqu
      t17 = 1/(k1bk2w-k1bp1t)
      t20 = 1/(2*k1bp1t-t6+t8)
      t22 = ne*t14*t17*t20
      t27 = t1*GW
      t39 = t14*t17
      t42 = 1/(2*k1bk2w-t6+t8)
      Diags(1) = zero
      Diags(2) = t1*EpsPk2w*GW*t22
      Diags(3) = -t1*EpsWEpsP*GW*t22
      Diags(4) = t27*(2*EpsWk1b*k1bk2w-2*EpsWp1t*k1bk2w-2*EpsWk1b*k1bp1t
     #+EpsWp1t*t6-EpsWp1t*t8)*ne*t39*t42*t20
      Diags(5) = zero
      Diags(6) = zero
      Diags(7) = zero
      Diags(8) = -t27*ne*t39*t42/2

      res(1:4) = Diags(1) * ME1(1:4)
     &         + Diags(2) * ME2(1:4)
     &         + Diags(3) * ME3(1:4)
     &         + Diags(4) * ME4(1:4)
     &         + Diags(5) * ME5(1:4)
     &         + Diags(6) * ME6(1:4)
     &         + Diags(7) * ME7(1:4)
     &         + Diags(8) * ME8(1:4)

      end subroutine numTopDectree


      subroutine numATopDectree(res)
      implicit none
      include 'VarATopLO.f'
      include 'commondecay.f'
      double complex p1tb(4),k2w(4),k1bb(4), kdum(4)
      double complex EpsP(4),EpsW(4),SpB(4),SpT(4)
      double complex res(4)
      double complex scalar1,scalar2,scalar3,scalar4
      double complex scalar5,scalar6,scalar7,scalar8
      double complex AME1(4), AME2(4), AME3(4), AME4(4)
      double complex AME5(4), AME6(4), AME7(4), AME8(4)
      double complex k1bbp1tb, k2wp1tb, k1bbk2w, Diags(8)
      double complex EpsWk1bb, EpsWp1tb, EpsWEpsP
      double complex EpsPk1bb, EpsPp1tb, EpsPk2w
      integer i, Dv
      Dv = 4
      do i = 1,4
         p1tb(i)=mom_t(i)
         k2w(i)=mom_dec(2,i)
         k1bb(i)=mom_dec(1,i)
         EpsP(i)=hel_dec(3,i)
         EpsW(i)=hel_dec(2,i)
         SpB(i)=hel_dec(1,i)
	 res(i) = (0.d0,0.d0)
      enddo
      call Calc_AME1(AME1)
      call Calc_AME2(AME2)
      call Calc_AME3(AME3)
      call Calc_AME4(AME4)
      call Calc_AME5(AME5)
      call Calc_AME6(AME6)
      call Calc_AME7(AME7)
      call Calc_AME8(AME8)

      EpsWk1bb = scf(Dv,EpsW,k1bb)
      EpsWp1tb = scf(Dv,EpsW,p1tb)
      EpsWEpsP = scf(Dv,EpsW,EpsP)
      EpsPk1bb = scf(Dv,EpsP,k1bb)
      EpsPp1tb = scf(Dv,EpsP,p1tb)
      EpsPk2w = scf(Dv,EpsP,k2w)
      k1bbp1tb = scf(Dv,k1bb,p1tb)
      k1bbk2w = scf(Dv,k1bb,k2w)
      k2wp1tb = scf(Dv,k2w,p1tb)
      t1 = EL**2
      t2 = t1*EpsWEpsP
      t7 = mt**2
      t9 = Mw**2
      t15 = 2*k1bbk2w*Qqd-t7*Qqd+t9*Qqd-2*k1bbk2w*Qqu+2*k1bbp1tb*Qqu
      t16 = ne*t15
      t19 = 1/(2*k1bbk2w-t7+t9)
      t23 = t19/(2*k1bbp1tb-t7+t9)
      t24 = t16*t23
      t34 = t1*GW
      t47 = 1/(k1bbk2w-k1bbp1tb)
      t48 = t15*t47
      Diags(1) = -4*t2*GW*mt*t24
      Diags(2) = 2*t1*EpsPk2w*GW*t24
      Diags(3) = -2*t2*GW*t24
      Diags(4) = t34*(2*EpsWk1bb*k1bbk2w-2*EpsWp1tb*k1bbk2w-2*EpsWk1bb*k
     #1bbp1tb+EpsWp1tb*t7-EpsWp1tb*t9)*ne*t48*t23
      Diags(5) = zero
      Diags(6) = -t34*mt*t16*t47*t19
      Diags(7) = zero
      Diags(8) = t34*ne*t48*t19/2

      res(1:4) = Diags(1) * AME1(1:4)
     &         + Diags(2) * AME2(1:4)
     &         + Diags(3) * AME3(1:4)
     &         + Diags(4) * AME4(1:4)
     &         + Diags(5) * AME5(1:4)
     &         + Diags(6) * AME6(1:4)
     &         + Diags(7) * AME7(1:4)
     &         + Diags(8) * AME8(1:4)
      end subroutine numATopDectree


!REALS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine numDecTopReals(heli,Wpar,k1b,klep,kneu,k3g,k4p,gauge,res)
      implicit none
      include 'commondecay.f'
      include 'VarReals.f'
      integer heli(1:3), heli_pho, heli_glu, i, DV, gauge, heli_w,heli_b
      integer Wpar(1:4)
      double precision  MomPol(1:4,1:2)
      double complex p1t(4),k2w(4),k1b(4),k3g(4),k4p(4)
      double complex kneu(4), klep(4)
      double complex EpsP(4),EpsW(4), Spi_bot(4), EpsG(4)
      double complex res(4), SMEcoeff(24), SME(24,1:4)
      double complex k1bp1t,k1bk2w,k2wp1t
      double complex k1bk3g, k2wk3g, k3gp1t
      double complex EpsWk1b, EpsWp1t, EpsWk3g, EpsWEpsP, EpsWEpsG
      double complex EpsPk1b, EpsPp1t, EpsPk2w, EpsPk3g, EpsPEpsG
      double complex EpsGk1b, EpsGk2w, EpsGp1t
      double complex one, two, four, three
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11
      Dv =4

      one = dcmplx(1.d0,0.d0)
      two = dcmplx(2.d0,0.d0)
      four = dcmplx(4.d0,0.d0)
      three = dcmplx(3.d0,0.d0)
      heli_pho = heli(3)
      heli_glu = heli(2)
      heli_b = heli(1)
      k2w(1:4) = kneu(1:4) + klep(1:4)
      p1t(1:4) = k2w(1:4) + k1b(1:4) + k3g(1:4) + k4p(1:4)

      call ubarSpi_Weyl(k1b,heli_b,Spi_bot)
!      call pol_massoff(klep,kneu,EpsW)
      do i=1,4
         MomPol(i,1) = dreal(klep(i))
         MomPol(i,2) = dreal(kneu(i))
      enddo
      call WpolvecS(Wpar(1),Wpar(2),MomPol,Wpar(3),Wpar(4),EpsW)


      call pol_mless(k4p,heli_pho,EpsP)
      call pol_mless(k3g,heli_glu,EpsG)

      If (gauge .eq. 1) then
         EpsP(1:4) = k4p(1:4)
      endif
      If (gauge .eq. 2) then
         EpsG(1:4) = k3g(1:4)
      endif

      k1bp1t = scf(Dv,k1b,p1t)
      k1bk2w = scf(Dv,k1b,k2w)
      k1bk3g = scf(Dv,k1b,k3g)

      k2wp1t = scf(Dv,k2w,p1t)
      k2wk3g = scf(Dv,k2w,k3g)

      k3gp1t = scf(Dv,p1t,k3g)


      EpsWk1b = scf(Dv,EpsW,k1b)
      EpsWp1t = scf(Dv,EpsW,p1t)
      EpsWk3g = scf(Dv,EpsW,k3g)
      EpsWEpsG = scf(Dv,EpsW,EpsG)
      EpsWEpsP = scf(Dv,EpsW,EpsP)


      EpsPk1b = scf(Dv,EpsP,k1b)
      EpsPp1t = scf(Dv,EpsP,p1t)
      EpsPk2w = scf(Dv,EpsP,k2w)
      EpsPk3g = scf(Dv,EpsP,k3g)
      EpsPEpsG = scf(Dv,EpsP,EpsG)


      EpsGk1b = scf(Dv,EpsG,k1b)
      EpsGp1t = scf(Dv,EpsG,p1t)
      EpsGk2w = scf(Dv,EpsG,k2w)




      call SpTop0(Spi_bot,p1t,mt,SME(1,1:4))

      call SpTop1(Spi_bot,EpsW,p1t,mt,SME(2,1:4))
      call SpTop1(Spi_bot,EpsP,p1t,mt,SME(4,1:4))
      call SpTop1(Spi_bot,k2w,p1t,mt,SME(3,1:4))
      call SpTop1(Spi_bot,EpsG,p1t,mt,SME(5,1:4))
      call SpTop1(Spi_bot,k3g,p1t,mt,SME(6,1:4))

      call SpTop2(Spi_bot,EpsW,EpsP,p1t,mt,SME(8,1:4))
      call SpTop2(Spi_bot,EpsW,EpsG,p1t,mt,SME(7,1:4))
      call SpTop2(Spi_bot,EpsW,k3g,p1t,mt,SME(9,1:4))

      call SpTop2(Spi_bot,EpsP,EpsG,p1t,mt,SME(10,1:4))
      call SpTop2(Spi_bot,EpsP,k3g,p1t,mt,SME(11,1:4))

      call SpTop2(Spi_bot,EpsG,k3g,p1t,mt,SME(12,1:4))

      call SpTop3(Spi_bot,EpsW,k2w,EpsP,p1t,mt,SME(13,1:4))
      call SpTop3(Spi_bot,EpsW,k2w,EpsG,p1t,mt,SME(14,1:4))
      call SpTop3(Spi_bot,EpsW,k2w,k3g,p1t,mt,SME(15,1:4))
      call SpTop3(Spi_bot,EpsW,EpsP,EpsG,p1t,mt,SME(16,1:4))
      call SpTop3(Spi_bot,EpsW,EpsP,k3g,p1t,mt,SME(17,1:4))
      call SpTop3(Spi_bot,EpsW,EpsG,k3g,p1t,mt,SME(18,1:4))

      call SpTop3(Spi_bot,k2w,EpsP,EpsG,p1t,mt,SME(19,1:4))
      call SpTop3(Spi_bot,k2w,EpsG,k3g,p1t,mt,SME(21,1:4))
      call SpTop3(Spi_bot,k2w,EpsP,k3g,p1t,mt,SME(20,1:4))

      call SpTop3(Spi_bot,EpsP,EpsG,k3g,p1t,mt,SME(22,1:4))


      call SpTop4(Spi_bot,EpsW,EpsP,EpsG,k3g,p1t,mt,SME(23,1:4))
      call SpTop5(Spi_bot,EpsW,k2w,EpsP,EpsG,k3g,p1t,mt,SME(24,1:4))


      t1 = EL**2
      t2 = t1*EpsPk3g
      t4 = EpsWEpsG*GS*GW
      t5 = t2*t4
      t6 = mt*ne
      t8 = 1/k1bk3g
      t9 = two*t8
      t10 = mt**2
      t11 = Mw**2
      t12 = k1bk2w*two
      t13 = k1bk3g*two
      t14 = k2wk3g*two
      t16 = 1/(-t10+t11+t12+t13+t14)
      t17 = t9*t16
      t18 = t6*Qqu*t17
      t20 = t1*EpsPEpsG
      t22 = EpsWk3g*GS*GW
      t23 = t20*t22
      t25 = t6*Qqd
      t26 = k2wp1t*two
      t28 = 1/(t10+t11-t26)
      t29 = t9*t28
      t30 = t25*t29
      t32 = t1*EpsGk1b
      t34 = EpsWEpsP*GS*GW
      t35 = t32*t34
      t38 = -Qqd+Qqu
      t39 = t6*t38
      t41 = k3gp1t*two
      t43 = 1/(t10-t11+t13-k1bp1t*two-t41)
      t44 = t9*t43
      t47 = t1*EpsGp1t
      t48 = t47*t34
      t49 = 1/k3gp1t
      t50 = two*t49
      t51 = t50*t43
      t55 = 1/(t10+t11+t14-t26-t41)
      t56 = t50*t55
      t59 = t1*EpsGk2w
      t61 = EpsWEpsP*four*GS
      t62 = t59*t61
      t63 = GW*mt
      t64 = t63*ne
      t66 = Qqd*t28*t55
      t67 = t64*t66
      t69 = t47*t61
      t72 = GS*GW
      t73 = t20*t72
      t74 = ne*Qqu
      t76 = t74*two*t16
      t80 = t59*EpsPk3g*GS*GW
      t81 = t74*t17
      t83 = t72*k2wk3g
      t84 = t20*t83
      t87 = four*GS*GW
      t91 = 1/(-t10+t11+t12)
      t93 = Qqu*t91*t16
      t96 = t72*t10
      t100 = t74*two*t91*t16
      t102 = t72*t11
      t106 = EpsPk2w*GS*GW
      t107 = t32*t106
      t108 = ne*Qqd
      t109 = t108*t29
      t113 = ne*t38
      t114 = t113*t44
      t116 = t47*t106
      t117 = t113*t51
      t119 = t108*t56
      t122 = EpsPk2w*four*GS
      t124 = GW*ne
      t125 = t124*t66
      t129 = t73*t76-t80*t81+t84*t81-t20*t87*k1bk2w*ne*t93+t20*t96*t100-
     #t20*t102*t100-t107*t109-t80*t109+t84*t109+t107*t114-t116*t117+t116
     #*t119+t59*t122*t125-t47*t122*t125
      t131 = t23*t81
      t138 = t62*t125
      t141 = t1*EpsWEpsG
      t142 = t141*t72
      t144 = t141*t87
      t146 = k2wp1t*ne*t66
      t149 = k2wk3g*ne*t66
      t152 = EpsWp1t*four*GS
      t158 = EpsWp1t*GS*GW
      t159 = t47*t158
      t162 = EpsWk3g*four*GS
      t166 = t47*t22
      t171 = t108*two*t28*t55
      t175 = t32*t22
      t178 = EpsWk1b*GS*GW
      t179 = t32*t178
      t181 = -t142*t76+t144*t146-t144*t149+t47*t152*t125-t59*t152*t125-t
     #159*t119+t59*t162*t125+t159*t117+t166*t119-t141*t96*t171-t141*t102
     #*t171+t175*t81+t179*t81
      t182 = t47*t178
      t187 = EpsWk1b*four*GS
      t189 = t124*t93
      t193 = t59*t22
      t195 = t141*t83
      t197 = t32*t158
      t206 = -t182*t74*t50*t91+t32*t187*t189+t59*t187*t189+t193*t81-t195
     #*t81+t197*t109+t193*t109+t179*t114-t195*t109+t175*t114-t197*t114-t
     #182*t117-t166*t117
      t208 = t2*t178
      t210 = t2*t22
      t214 = t1*EpsPk2w
      t215 = t214*t22
      t217 = t2*t158
      t219 = t1*EpsWEpsP
      t220 = t219*t83
      t227 = t219*t72
      t233 = t219*t87
      t240 = t208*t81+t210*t81+t2*t187*t189-t215*t109+t217*t109+t220*t10
     #9+t208*t114+t215*t114+t210*t114-t217*t114-t220*t114-t227*t108*two*
     #t55-t214*t162*t125+t233*t149-t233*t146+t219*t96*t171+t219*t102*t17
     #1
      t241 = t20*t178
      t245 = t214*t4
      t247 = t59*t34
      t249 = t20*t158
      t260 = -t241*t81-t131-t20*t187*t189+t245*t109-t247*t109-t249*t109-
     #t245*t114+t247*t114-t241*t114-t23*t114+t249*t114+t214*EpsWEpsG*fou
     #r*GS*t125-t138
      t261 = t2*t72
      t264 = t6*Qqu*t8*t16
      t266 = t72*mt
      t271 = t6*Qqd*t8*t28
      t274 = t47*t72
      t279 = t32*t72
      t283 = t59*t266
      t288 = t6*Qqd*t49*t55
      t299 = t1*EpsWk3g
      t300 = t299*t72
      t325 = t74*t49*t91
      t328 = t74*t8*t16
      t331 = t59*t72
      t334 = t108*t8*t28
      t337 = t108*t49*t55
      t339 = t331*t171
      t350 = t1*GS
      t351 = t350*GW
      t355 = t350*GW*k2wk3g
      t358 = t350*GW*t10
      t360 = t74*t91*t16
      t363 = t350*GW*t11
      t372 = t108*t28*t55
      t379 = -t351*t74*t16-t355*t328-t358*t360+t363*t360+t350*GW*k1bk2w*
     #t100-t355*t334+t351*t108*t55-t358*t372-t363*t372-t355*t171+t350*GW
     #*k2wp1t*t171
      t383 = t214*t72
      t386 = t113*t8*t43
      t389 = t113*t49*t43
      t394 = t300*t328
      t409 = t1*EpsWk1b*t72
      t414 = t1*EpsWp1t*t72
      t425 = -t409*t325+t409*t328+t394+t409*t100+t414*t334+t409*t386+t30
     #0*t386-t414*t386-t409*t389-t300*t389+t414*t389+t300*t337-t414*t337
     #+t414*t171
      t426 = t350*t64
      t427 = one*Qqu
      t428 = 1/two
      t429 = t49*t428
      t431 = t427*t429*t91
      t433 = t8*t428
      t435 = t427*t433*t16
      t437 = t350*t63
      t439 = one*Qqd
      t441 = t439*t433*t28
      t444 = t439*t429*t55
      t448 = t350*t124
      SMEcoeff(1) = -t5*t18+t23*t18-t5*t30-t35*t30+t23*t30+t35*t39*t44-t
     #48*t39*t51+t48*t25*t56+t62*t67-t69*t67
      SMEcoeff(2) = t129
      SMEcoeff(3) = t5*t81-t131+t5*t109+t35*t109-t23*t109-t35*t114+t48*t
     #117-t48*t119-t138+t69*t125
      SMEcoeff(4) = t181+t206
      SMEcoeff(5) = t240
      SMEcoeff(6) = t260
      SMEcoeff(7) = t261*t264+t2*t266*t100+t261*t271
      SMEcoeff(8) = -t274*t6*Qqu*t49*t91+t279*t264+t32*t266*t100+t283*t1
     #00+t279*t271-t274*t288-t283*t171+t47*t266*t171
      SMEcoeff(9) = -t73*t264-t20*t266*t100-t73*t271
      SMEcoeff(10) = -t300*t264-t300*t271-t299*t266*t171
      SMEcoeff(11) = t142*t264+t142*t271+t141*t266*t171
      SMEcoeff(12) = -t227*t271+t227*t6*t38*t8*t43-t227*t6*t38*t49*t43+t
     #227*t288-t219*t266*t171
      SMEcoeff(13) = -t274*t325+t279*t328+t279*t100+t331*t100+t279*t334-
     #t274*t337-t339+t274*t171
      SMEcoeff(14) = t261*t328+t261*t100+t261*t334
      SMEcoeff(15) = -t73*t328-t73*t100-t73*t334
      SMEcoeff(16) = t379
      SMEcoeff(17) = t331*t328+t331*t334+t339
      SMEcoeff(18) = -t383*t334+t383*t386-t383*t389+t383*t337-t383*t171
      SMEcoeff(19) = t394+t300*t334+t300*t171
      SMEcoeff(20) = -t142*t328-t142*t334-t142*t171
      SMEcoeff(21) = t227*t334-t227*t386+t227*t389-t227*t337+t227*t171
      SMEcoeff(22) = t425
      SMEcoeff(23) = -t426*t431+t426*t435+t437*t360+t426*t441-t426*t444+
     #t437*t372
      SMEcoeff(24) = -t448*t431+t448*t435+t351*t360+t448*t441-t448*t444+
     #t351*t372

      do i =1,4
         res(i) = zero
      enddo
      do i=1,24
         res(1:4) = res(1:4) + SMEcoeff(i)*SME(i,1:4)
      enddo


      return
      end subroutine numDecTopReals



      subroutine numDecATopReals(heli,Wpar,k1bb,klep,kneu,k3g,k4p,gauge,res)
      implicit none
      include 'commondecay.f'
      include 'VarAReals.f'
      integer heli(1:3), heli_pho, heli_glu, i, DV, gauge, heli_w,heli_b
      integer Wpar(1:4)
      double precision  MomPol(1:4,1:2)
      double complex p1tb(4),k2w(4),k1bb(4),k3g(4),k4p(4)
      double complex kneu(4), klep(4)
      double complex EpsP(4),EpsW(4), Spi_bot(4), EpsG(4)
      double complex res(4), SMEcoeff(24), SME(24,1:4)
      double complex k1bbp1tb,k1bbk2w,k2wp1tb
      double complex k1bbk3g, k2wk3g, k3gp1tb
      double complex EpsWk1bb, EpsWp1tb, EpsWk3g, EpsWEpsP, EpsWEpsG
      double complex EpsPk1bb, EpsPp1tb, EpsPk2w, EpsPk3g, EpsPEpsG
      double complex EpsGk1bb, EpsGk2w, EpsGp1tb
      double complex one, two, four, three
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11
      Dv =4

      one = dcmplx(1.d0,0.d0)
      two = dcmplx(2.d0,0.d0)
      four = dcmplx(4.d0,0.d0)
      three = dcmplx(3.d0,0.d0)
      heli_pho = heli(3)
      heli_glu = heli(2)
      heli_b = heli(1)
      k2w(1:4) = kneu(1:4) + klep(1:4)
      p1tb(1:4) = k2w(1:4) + k1bb(1:4) + k3g(1:4) + k4p(1:4)

      call vSpi_Weyl(k1bb,heli_b,Spi_bot)
!      call pol_massoff2(klep,kneu,EpsW)
      do i=1,4
         MomPol(i,1) = dreal(klep(i))
         MomPol(i,2) = dreal(kneu(i))
      enddo
      call WpolvecS(Wpar(1),Wpar(2),MomPol,Wpar(3),Wpar(4),EpsW)


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
      end subroutine numDecATopReals




       subroutine SpTop0(sp1,slash,mt,sp2)
       implicit none
       double complex sp2(4),sp1(4),slash(4),sp1n(4)
       double precision mt
       integer Dv, Ds
       Dv = 4
       Ds = 4
       call spb2_weylS(Dv,Ds,sp1,slash,sp1n)
       sp2(1:4) = sp1n(1:4) + mt*sp1(1:4)
       return
       end subroutine SpTop0

       subroutine SpTop1(sp1,ap1,slash,mt,sp2)
       implicit none
       double complex sp2(4),sp1(4), ap1(4),slash(4)
       double complex sp1n(4), sp2n(4)
       double precision mt
       integer Dv, Ds
       Dv = 4
       Ds = 4

       call spb2_weylS(Dv,Ds,sp1,ap1,sp1n)
       call spb2_weylS(Dv,Ds,sp1n,slash,sp2n)
       sp2(1:4) = sp2n(1:4) + mt*sp1n(1:4)
       return
       end subroutine SpTop1


       subroutine SpTop2(sp1,ap1,ap2,slash,mt,sp2)
       implicit none
       double complex sp2(4),sp1(4), ap1(4), slash(4)
       double complex sp1n(4),ap2(4), sp2n(4)
       double precision mt
       integer Dv, Ds
       Dv = 4
       Ds = 4

       call spb2_weylS(Dv,Ds,sp1,ap1,sp1n)
       call spb2_weylS(Dv,Ds,sp1n,ap2,sp2n)
       call spb2_weylS(Dv,Ds,sp2n,slash,sp1n)
       sp2(1:4) = sp1n(1:4) + mt*sp2n(1:4)
       return
       end subroutine SpTop2


       subroutine SpTop3(sp1,ap1,ap2,ap3,slash,mt,sp2)
       implicit none
       double complex sp2(4),sp1(4), ap1(4), slash(4)
       double complex sp1n(4), sp2n(4), ap2(4), ap3(4)
       double precision mt
       integer Dv, Ds
       Dv = 4
       Ds = 4

       call spb2_weylS(Dv,Ds,sp1,ap1,sp1n)
       call spb2_weylS(Dv,Ds,sp1n,ap2,sp2n)
       call spb2_weylS(Dv,Ds,sp2n,ap3,sp1n)
       call spb2_weylS(Dv,Ds,sp1n,slash,sp2n)
       sp2(1:4) = sp2n(1:4) + mt*sp1n(1:4)
       return
       end subroutine SpTop3

       subroutine SpTop4(sp1,ap1,ap2,ap3,ap4,slash,mt,sp2)
       implicit none
       double complex sp2(4),sp1(4), ap1(4), slash(4)
       double complex sp1n(4), sp2n(4), ap2(4), ap3(4), ap4(4)
       double precision mt
       integer Dv, Ds
       Dv = 4
       Ds = 4

       call spb2_weylS(Dv,Ds,sp1,ap1,sp1n)
       call spb2_weylS(Dv,Ds,sp1n,ap2,sp2n)
       call spb2_weylS(Dv,Ds,sp2n,ap3,sp1n)
       call spb2_weylS(Dv,Ds,sp1n,ap4,sp2n)
       call spb2_weylS(Dv,Ds,sp2n,slash,sp1n)
       sp2(1:4) = sp1n(1:4) + mt*sp2n(1:4)

       return
       end subroutine SpTop4

       subroutine SpTop5(sp1,ap1,ap2,ap3,ap4,ap5,slash,mt,sp2)
       implicit none
       double complex sp2(4),sp1(4), ap1(4), ap4(4),slash(4)
       double complex sp1n(4), sp2n(4), ap2(4), ap3(4), ap5(4)
       double precision mt
       integer Dv, Ds
       Dv = 4
       Ds = 4

       call spb2_weylS(Dv,Ds,sp1,ap1,sp1n)
       call spb2_weylS(Dv,Ds,sp1n,ap2,sp2n)
       call spb2_weylS(Dv,Ds,sp2n,ap3,sp1n)
       call spb2_weylS(Dv,Ds,sp1n,ap4,sp2n)
       call spb2_weylS(Dv,Ds,sp2n,ap5,sp1n)
       call spb2_weylS(Dv,Ds,sp1n,slash,sp2n)
       sp2(1:4) = sp2n(1:4) + mt*sp1n(1:4)

       return
       end subroutine SpTop5


       subroutine SpATop0(sp1,slash,mt,sp2)
       implicit none
       double complex sp2(4),sp1(4),slash(4),sp1n(4)
       double complex pslash(4)
       double precision mt
       integer Dv, Ds
       Dv = 4
       Ds = 4
       pslash(1:4) = -slash(1:4)
       call spi2_weylS(Dv,Ds,sp1,pslash,sp1n)
       sp2(1:4) = sp1n(1:4) + mt*sp1(1:4)

       return
       end subroutine SpATop0

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
       call spi2_weylS(Dv,Ds,sp1,ap1,sp1n)
       call spi2_weylS(Dv,Ds,sp1n,pslash,sp2n)
       sp2(1:4) = sp2n(1:4) + mt*sp1n(1:4)
       return
       end subroutine SpATop1


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
       call spi2_weylS(Dv,Ds,sp1,ap1,sp1n)
       call spi2_weylS(Dv,Ds,sp1n,ap2,sp2n)
       call spi2_weylS(Dv,Ds,sp2n,pslash,sp1n)
       sp2(1:4) = sp1n(1:4) + mt*sp2n(1:4)
       return
       end subroutine SpATop2


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

       call spi2_weylS(Dv,Ds,sp1,ap1,sp1n)
       call spi2_weylS(Dv,Ds,sp1n,ap2,sp2n)
       call spi2_weylS(Dv,Ds,sp2n,ap3,sp1n)
       call spi2_weylS(Dv,Ds,sp1n,pslash,sp2n)
       sp2(1:4) = sp2n(1:4) + mt*sp1n(1:4)
       return
       end subroutine SpATop3

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

       call spi2_weylS(Dv,Ds,sp1,ap1,sp1n)
       call spi2_weylS(Dv,Ds,sp1n,ap2,sp2n)
       call spi2_weylS(Dv,Ds,sp2n,ap3,sp1n)
       call spi2_weylS(Dv,Ds,sp1n,ap4,sp2n)
       call spi2_weylS(Dv,Ds,sp2n,pslash,sp1n)
       sp2(1:4) = sp1n(1:4) + mt*sp2n(1:4)

       return
       end subroutine SpATop4

       subroutine SpATop5(sp1,ap5,ap4,ap3,ap2,ap1,slash,mt,sp2)
       implicit none
       double complex sp2(4),sp1(4), ap1(4), ap4(4),slash(4)
       double complex sp1n(4), sp2n(4), ap2(4), ap3(4), ap5(4)
       double complex pslash(4)
       double precision mt
       integer Dv, Ds
       Dv = 4
       Ds = 4
       pslash(1:4) = -slash(1:4)

       call spi2_weylS(Dv,Ds,sp1,ap1,sp1n)
       call spi2_weylS(Dv,Ds,sp1n,ap2,sp2n)
       call spi2_weylS(Dv,Ds,sp2n,ap3,sp1n)
       call spi2_weylS(Dv,Ds,sp1n,ap4,sp2n)
       call spi2_weylS(Dv,Ds,sp2n,ap5,sp1n)
       call spi2_weylS(Dv,Ds,sp1n,pslash,sp2n)
       sp2(1:4) = sp2n(1:4) + mt*sp1n(1:4)

       return
       end subroutine SpATop5


!DIPOLE
!DIPOLE
      subroutine TopDipole(k1b,klep,kneu,k3g,k4p,alphadip,
     &k1b_tilde,k4p_tilde,klep_tilde,kneu_tilde,res)
      implicit none
      include 'commondecay.f'
      integer  DV
      double complex p1t(4),k2w(4),k1b(4),k3g(4), k4p(4)
      double complex k2w_dip(4), k2w_tilde(4), kneu(4), klep(4)
      double complex k1b_tilde(4), p1t_tilde(4), k4p_tilde(4)
      double complex kneu_tilde(4), klep_tilde(4)
      double complex EpsP(4),EpsW(4), lo(4)
      double complex  kdum(4), kdum1(4)
      double complex k1bk3g,k3gp1t, k1bp1t
      double precision condition1, condition2, r
      double precision alphadip, ptpho, ang, theta_y, theta_z
      double precision y,z, Dipole, y_con, z_con
      double precision condit1, condit2, rat, res
      Dv = 4
      k2w(1:4) = kneu(1:4) + klep(1:4)
      k2w_dip(1:4) = k2w(1:4) + k4p(1:4)
      p1t(1:4) = k2w(1:4) + k1b(1:4) + k4p(1:4) + k3g(1:4)


      lo(1:4) = zero
      k1B_tilde(1:4) = zero
      k4p_tilde(1:4) = zero
      klep_tilde(1:4) = zero
      kneu_tilde(1:4) = zero



      k1bp1t = scf(Dv,k1b,p1t)
      k1bk3g = scf(Dv,k1b,k3g)
      k3gp1t = scf(Dv,p1t,k3g)


      condit1 = dreal(2.d0*k3g(1)/mt)
      condit2 = dsqrt(dreal(2.d0*scf(Dv,k3g,k1b))/mt**2)

      r = dsqrt(dreal(scf(Dv,k2w_dip,k2w_dip)/mt**2))
      y = 2.d0*dreal(k1bk3g)/(mt**2*(1.d0-r)**2)
      z = 1.d0 - 2.d0*dreal(k3gp1t)/mt**2/(1.d0-r**2)


!      y_con = Theta_Y(y,r,z,alphadip)
!      z_con = Theta_Z(z,alphadip)

      condition1 = 1.d0-alphadip-z
      condition2 = y-alphadip*(1.d0+r)**2*z*(1.d0-z)/
     &(z+r**2*(1.d0-z))



      call DipoleBoost(k2w_dip,p1t,kneu,kneu_tilde)
      call DipoleBoost(k2w_dip,p1t,klep,klep_tilde)
      call DipoleBoost(k2w_dip,p1t,k4p,k4p_tilde)

      k1b_tilde(1:4) = p1t(1:4) - k4p_tilde(1:4)
     & - kneu_tilde(1:4) - klep_tilde(1:4)


      Dipole =0.d0
      If (condition1 .lt. 0.d0 .or. condition2 .lt. 0.d0) then

            Dipole = GS**2*CF*(1.d0/dreal(k1bk3g)*
     &           (2.d0/(1.d0-z)-1.d0-z) - mt**2/dreal(k3gp1t)**2)
      endif

      res = Dipole

      end subroutine TopDipole



      subroutine ATopDipole(k1bb,klep,kneu,k3g,k4p,alphadip,
     &   k1bb_tilde,k4p_tilde,klep_tilde,kneu_tilde,res)
      implicit none
      include 'commondecay.f'
      integer  DV
      double complex p1tb(4),k2w(4),k1bb(4),k3g(4), k4p(4)
      double complex k2w_dip(4), k2w_tilde(4), kneu(4), klep(4)
      double complex k1bb_tilde(4), p1tb_tilde(4), k4p_tilde(4)
      double complex kneu_tilde(4), klep_tilde(4)
      double complex EpsP(4),EpsW(4), lo(4)
      double complex  kdum(4), kdum1(4)
      double complex k1bbk3g,k3gp1tb, k1bbp1tb
      double precision condition1, condition2, r
      double precision alphadip, ptpho, ang, rat
      double precision y,z, Dipole, res
      Dv = 4
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



      r = dsqrt(dreal(scf(Dv,k2w_dip,k2w_dip)/mt**2))
      y = 2.d0*dreal(k1bbk3g)/(mt**2*(1.d0-r)**2)
      z = 1.d0 - 2.d0*dreal(k3gp1tb)/mt**2/(1.d0-r**2)




      condition1 = 1.d0-alphadip-z
      condition2 = y-alphadip*(1.d0+r)**2*z*(1.d0-z)/
     &(z+r**2*(1.d0-z))



      call DipoleBoost(k2w_dip,p1tb,kneu,kneu_tilde)
      call DipoleBoost(k2w_dip,p1tb,klep,klep_tilde)
      call DipoleBoost(k2w_dip,p1tb,k4p,k4p_tilde)

      k1bb_tilde(1:4) = p1tb(1:4) - k4p_tilde(1:4)
     & - kneu_tilde(1:4) - klep_tilde(1:4)


      Dipole =0.d0
      If (condition1 .lt. 0.d0 .or. condition2 .lt. 0.d0) then
         Dipole = GS**2*CF*(1.d0/dreal(k1bbk3g)*
     &           (2.d0/(1.d0-z)-1.d0-z) - mt**2/dreal(k3gp1tb)**2)

      endif
      res = Dipole

      end subroutine ATopDipole



      subroutine DipoleBoost(k2w,p1t,vec, vec_b)
      implicit none
      integer Dv
      double complex p1t(4),k2w(4), vec_b(4), vec(4)
      double complex sinh, cosh, k2wq, p1tq, sqt
      double complex o5, k2wp1t
      Dv =4
      k2wq = scf(Dv,k2w,k2w)
      p1tq = scf(Dv,p1t,p1t)
      k2wp1t = scf(Dv,k2w,p1t)
      sqt = sqrt(k2wp1t**2 - p1tq*k2wq)
      o5 = dcmplx(0.5d0,0.d0)

      sinh = o5/p1tq/k2wq *
     &(-(p1tq-k2wq)*k2wp1t + (p1tq+k2wq)*sqt)

      cosh = o5/p1tq/k2wq *
     &((p1tq+k2wq)*k2wp1t - (p1tq-k2wq)*sqt)

      vec_b(1:4) = vec(1:4)
     &+ sinh/sqt*(scf(Dv,k2w,vec)*p1t(1:4) - scf(Dv,p1t,vec)*k2w(1:4))
     & +(cosh-dcmplx(1.d0,0.d0))/sqt**2*
     &(k2wp1t*(scf(Dv,k2w,vec)*p1t(1:4) + scf(Dv,p1t,vec)*k2w(1:4))
     & - k2wq*scf(Dv,p1t,vec)*p1t(1:4) - p1tq*scf(Dv,k2w,vec)*k2w(1:4))
      end subroutine DipoleBoost


      subroutine WpolvecS(Wpm,WDK,MomPol,order,Wpol,EpsW)
        implicit none
        integer :: Wpm, WDK, order, Wpol
        real(8) ::  MomPol(1:4,1:2)
        complex(8) :: EpsW(1:4)
        EpsW(1:4) = WpolvecSS(Wpm,WDK,MomPol,order)
        return
      end subroutine WpolvecS

      function WpolvecSS(Wpm,WDK,MomPol,order,LamW,LamP)
        use ModParameters
        use ModMisc
        implicit none
        complex(8) :: WpolvecSS(1:4)
        complex(8) ::  Wcurr1(1:4), Wcurr2(1:4), Wmom(1:4), wcurr(1:4)
        complex(8) :: LepSpi(1:4), wcurr3(1:4), NeuSpi(1:4)
        integer :: WDK, Wpm, order, Qw
        integer, optional :: LamW, LamP
        real(8) :: MomPol(:,:)
        real(8),parameter :: sq2 =  1.414213562373095D0
        real(8),parameter :: g2_weak = 4d0*dsqrt(2d0)*m_W**2*GF
        real(8),parameter :: g_weak = dsqrt(g2_weak)
        real(8) :: sw, EL, coup, Q_f1, Q_f2
        complex(8) :: klep_dn(1:4), kneu_up(1:4), kpho(1:4)
        complex(8):: PropMom(1:4), EpsP(1:4)
        complex(8) :: SpiLep_dn(1:4), SpiNeu_up(1:4), SpiLep2(1:4)
        complex(8) :: SpiLep1(1:4), coupl_sqrt
        integer :: Dv
        ! Resolving markus coupling structure, sin(Theta_w)
        sw = dsqrt(4.d0*DblPi*alpha/g2_weak)
        EL = dsqrt(4.d0*DblPi*alpha)
        Dv =4
        Qw = Wpm
        If (order .ne. 0 ) then
           Write(*,*) "Only, LO is implemented <--> order =0"
        endif
        If(Qw .ne. -1 .and. Qw .ne. 1) then
           write(*,*) "Sorry, we can only deliver W+/W- polarisation vectors"
           stop
        endif
        if (present(LamW) .and. WDK .gt. 0) then
           write(*,*) "Incorrect call of WpolvecSS"
           write(*,*) "LamW is set to", LamW
           write(*,*) "WDK should be zero but is", WDK
           stop
        endif

        klep_dn(1:4) = dcmplx(MomPol(1:4,1),0.d0)
        kneu_up(1:4) = dcmplx(MomPol(1:4,2),0.d0)

        if(present(Lamp) ) then
           if (size(MomPol,dim=2) .ne. 3) then
              write(*,*) "Incorrect call of WpolvecSS"
              write(*,*) "LamP is set to", LamP
              write(*,*) "But photon-momentum is missing"
              stop
           endif
           kpho(1:4) = dcmplx(MomPol(1:4,3),0.d0)
        endif

        If(Q_top .eq. -4.d0/3.d0 ) then
           Qw= -Qw
        endif
        ! W+/W- vertex is I *EL /sqrt(2)/sw  times ne * sq2
        ! absorbing Factor -ne/sq2 introduced by vbqq_Weyl
        coup = -ne*EL/sq2/sw  * (ne *sq2)

        ! No W decay
        If(WDK .eq. 0) then
           WMom(1:4) = klep_dn(1:4)
           WpolvecSS(1:4) =  pol_mass(klep_dn+kneu_up,m_w,LamW)
           return
        endif


        ! ALL W+ cases
        If(Qw .eq. 1) then
           call    vSpi_Weyl(klep_dn,+1,SpiLep_dn)     ! l+ or dn_bar
           call ubarSpi_Weyl(kneu_up,-1,SpiNeu_up)  ! nu or up
           ! Anti-Lepton + Neutrino
           if(WDK .eq. 1) then
              WpolvecSS(1:4) = vbqq_Weyl(Dv,SpiNeu_up,SpiLep_dn)* coup
              return
           endif
           ! Anti-Down + Up
           if(WDK .eq. 3) then
              WpolvecSS(1:4) = vbqq_Weyl(Dv,SpiNeu_up,SpiLep_dn)* coup
     &* dsqrt(3.d0) * sq2 ! color-Factor dsqrt(3) *  Flavour  dsqrt(2)
              return
           endif

           ! Anti-Lepton + Neutrino + Photon
           ! Anti-Down + Up + Photon
           if(WDK .eq. 4 .or. WDK .eq. 2) then
              call pol_mless(kpho,LamP,EpsP)      ! photon
              ! print *, "check of gauge invariance for photon in decay 1"
              ! EpsP(1:4)=dcmplx(kpho(1:4))

              If(WDK .eq. 4) then
                 Q_f1 = -1.d0/3.d0
                 Q_f2 = 2.d0/3.d0
                 coup = sq2 * dsqrt(3.d0)! FLAVOUR * COLOUR
              endif

              If(WDK .eq. 2) then
                 Q_f1 = -1.d0
                 Q_f2 = 0.d0
                 coup =1.d0
              endif

              coupl_sqrt = (0d0,-1d0)*dsqrt(alpha4Pi)
              !       diagram 1: photon emission off anti-down
              LepSpi(1:4) = spi2_Weyl(EpsP,SpiLep_dn) * coupl_sqrt*(-Q_f1)
              PropMom(1:4) = kpho(1:4) + klep_dn(1:4)
              LepSpi(1:4) = -spi2_Weyl(PropMom(1:4),LepSpi)*(0d0,1d0)
     &        /(2d0*(kpho(1:4).dot.klep_dn(1:4)))
              WCurr1(1:4)  = vbqq_Weyl(4,SpiNeu_up,LepSpi)  * g_weak ! vbqq introduces -i/Sqrt(2)

              !       diagram 2: photon emission off up
              NeuSpi(1:4) = spb2_Weyl(SpiNeu_up,EpsP) * coupl_sqrt*(-Q_f2)
              PropMom(1:4) = kpho(1:4) + kneu_up(1:4)
              NeuSpi(1:4) = spb2_Weyl(NeuSpi,PropMom(1:4))*(0d0,1d0)
     & /(2d0*(kpho(1:4).dot.kneu_up(1:4)))
              WCurr2(1:4)  = vbqq_Weyl(4,NeuSpi,SpiLep_dn)  * g_weak ! vbqq introduces -i/Sqrt(2)

              !       diagram 3: photon emission off W+ boson
              WCurr3(1:4)  = vbqq_Weyl(4,SpiNeu_up,SpiLep_dn) * g_weak ! vbqq introduces -i/Sqrt(2)
              WMom(1:4) = klep_dn(1:4)+kneu_up(1:4)
              WCurr3(1:4) = Vertex_WPW(dreal(kpho(1:4)),dreal(WMom(1:4)),
     &EpsP,WCurr3(1:4))
     &* (0d0,-1d0)/((WMom(1:4).dot.WMom(1:4))-m_W**2) *(-coupl_sqrt*Q_Wp)
              !       connect to quark current
              WCurr(1:4) = WCurr1(1:4) + WCurr2(1:4) + Wcurr3(1:4)
              WMom(1:4) = klep_dn(1:4)+kneu_up(1:4)+kpho(1:4)

              WCurr(1:4) = (-WCurr(1:4) + (WCurr(1:4).dot.WMom(1:4))/m_W**2
     &* WMom(1:4)) *coup ! Flavour and Color

              WpolvecSS(1:4) = Wcurr(1:4)
              return
           endif

        endif

        ! All W- cases
        If(Qw .eq. -1) then
           call ubarSpi_Weyl(klep_dn,-1,SpiLep_dn)  ! l- or dn
           call    vSpi_Weyl(kneu_up,+1,SpiNeu_up)     ! nubar or up_bar

           !Lepton +  Anti-Neutrino
           if(WDK .eq. 1) then
              WpolvecSS(1:4) = vbqq_Weyl(Dv,SpiLep_dn,SpiNeu_up)* coup
              return
           endif
           !Down +  Anti-Up
           if(WDK .eq. 3) then
              WpolvecSS(1:4) =  vbqq_Weyl(Dv,SpiLep_dn,SpiNeu_up)* coup *
     &dsqrt(3.d0) * sq2 ! color-Factor dsqrt(3) *  Flavour  dsqrt(2)
              return
           endif

           ! Lepton + Anti-Neutrino + Photon
           ! Down + Anti-Up + Photon
           if(WDK .eq. 4 .or. WDK .eq. 2) then
              If(WDK .eq. 4) then
                 Q_f1 = -1.d0/3.d0
                 Q_f2 = 2.d0/3.d0
                 coup = sq2 * dsqrt(3.d0)! FLAVOUR * COLOUR
              endif

              If(WDK .eq. 2) then
                 Q_f1 = -1.d0
                 Q_f2 = 0.d0
                 coup =1.d0
              endif

              !        call ubarSpi_Weyl(klep_dn(1:4),-1,SpiLep_dn)  ! l- or dn
              !        call    vSpi_Weyl(kneu_up(1:4),+1,SpiNeu_up)     ! nubar or up_bar
              call pol_mless(kpho(1:4),LamP,EpsP)      ! photon
              !        print *, "check of gauge invariance for photon in decay 2"
              !        EpsP(1:4)=dcmplx(kpho(1:4))

              coupl_sqrt = (0d0,-1d0)*dsqrt(alpha4Pi)
              !       diagram 1: photon emission off down
              LepSpi(1:4) = spb2_Weyl(SpiLep_dn,EpsP(1:4)) * coupl_sqrt*(Q_f1)
              PropMom(1:4) = kpho(1:4) + klep_dn(1:4)
              LepSpi(1:4) = spb2_Weyl(LepSpi(1:4),PropMom(1:4))*(0d0,1d0)
     &/(2d0*(kpho(1:4).dot.klep_dn(1:4)))
              WCurr1(1:4)  = vbqq_Weyl(4,LepSpi(1:4),SpiNeu_up) * g_weak ! vbqq introduces -i/Sqrt(2)

              !       diagram 2: photon emission off anti-up
              NeuSpi(1:4) = spi2_Weyl(EpsP(1:4),SpiNeu_up) * coupl_sqrt*(Q_f2)
              PropMom(1:4) = kpho(1:4) + kneu_up(1:4)
              NeuSpi(1:4) = -spi2_Weyl(PropMom(1:4),NeuSpi(1:4))*(0d0,1d0)
     &/(2d0*(kpho(1:4).dot.kneu_up(1:4)))
              WCurr2(1:4)  = vbqq_Weyl(4,Spilep_dn,NeuSpi) * g_weak ! vbqq introduces -i/Sqrt(2)

              !       diagram 3: photon emission off W- boson!                       MINUS INTRODUCED TO ENSURE GAUGE INVARIANCE!
              WCurr3(1:4)  = vbqq_Weyl(4,SpiLep_dn,SpiNeu_up)  * g_weak ! vbqq introduces -i/Sqrt(2)
              WMom(1:4) = klep_dn(1:4)+kneu_up(1:4)
              WCurr3(1:4) = Vertex_WPW(dreal(kpho(1:4)),dreal(WMom(1:4))
     &,EpsP(1:4),WCurr3(1:4))
     &*(0d0,-1d0)/((WMom(1:4).dot.WMom(1:4))-m_W**2) *(-coupl_sqrt*Q_Wm)*(-1d0)
              !       connect to quark current
              WCurr(1:4) = (WCurr1(1:4) + WCurr2(1:4) + Wcurr3(1:4))
              WMom(1:4) = klep_dn(1:4)+kneu_up(1:4)+kpho(1:4)

              WCurr(1:4) = (-WCurr(1:4) + (WCurr(1:4).dot.WMom(1:4))/m_W**2
     &* WMom(1:4))* coup ! Flavour & Color

              WpolvecSS(1:4) = Wcurr(1:4)
              return
           endif



        endif

      end function WpolvecSS


      END MODULE ModTTBP_NLODK


