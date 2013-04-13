MODULE ModWDecay
implicit none

save

contains


  function Wpolvec(Wpm,WDK,MomPol,order,LamW,LamP)
    use ModParameters
    use ModMisc
    implicit none
    complex(8) :: Wpolvec(1:4)
    complex(8) ::  Wcurr1(1:4), Wcurr2(1:4), Wmom(1:4), wcurr(1:4), NeuSpi(1:4), LepSpi(1:4), wcurr3(1:4)
    integer :: WDK, Wpm,QW, order
    integer, optional :: LamW, LamP
    real(8) :: MomPol(:,:)
    real(8),parameter :: sq2 =  1.414213562373095D0
    real(8) :: coup, Q_f1, Q_f2, check1, check2, check3, check4, check5
    complex(8) :: klep_dn(1:4), kneu_up(1:4), kpho(1:4), PropMom(1:4), EpsP(1:4)
    complex(8) :: SpiLep_dn(1:4), SpiNeu_up(1:4), SpiLep2(1:4),SpiLep1(1:4), coupl_sqrt
    integer :: Dv
    logical nan1, nan2, nan3


    QW=Wpm
! Resolving markus coupling structure, sin(Theta_w)
    Dv =4
    If (order .ne. 0 ) then
       Write(*,*) "Only, LO is implemented <--> order =0"
    endif
    If(QW .ne. -1 .and. QW .ne. 1) then
       write(*,*) "Sorry, we can only deliver W+/W- polarisation vectors"
       stop
    endif
    if (present(LamW) .and. WDK .gt. 0) then
       write(*,*) "Incorrect call of Wpolvec"
       write(*,*) "LamW is set to", LamW
       write(*,*) "WDK should be zero but is", WDK
       stop
    endif

    klep_dn(1:4) = dcmplx(MomPol(1:4,1),0.d0)
    kneu_up(1:4) = dcmplx(MomPol(1:4,2),0.d0)

    if(present(Lamp) ) then
       if (size(MomPol,dim=2) .ne. 3) then
          write(*,*) "Incorrect call of Wpolvec"
          write(*,*) "LamP is set to", LamP
          write(*,*) "But photon-momentum is missing"
          stop
       endif
       kpho(1:4) = dcmplx(MomPol(1:4,3),0.d0)
    endif

    if(Q_top .eq. -4.d0/3.d0) QW=-QW

! W+/W- vertex is I *EL /sqrt(2)/sw  times ne * sq2
! absorbing Factor -ne/sq2 introduced by vbqq_Weyl
    coup = -ne*EL/sq2/sw  * (ne *sq2)
    ! No W decay
    If(WDK .eq. 0) then
       WMom(1:4) = klep_dn(1:4)
       Wpolvec(1:4) =  pol_mass(klep_dn+kneu_up,m_w,LamW)
       return
    endif


    ! ALL W+ cases
    If(QW .eq. 1) then
       call    vSpi_Weyl(klep_dn,+1,SpiLep_dn)     ! l+ or dn_bar
       call ubarSpi_Weyl(kneu_up,-1,SpiNeu_up)  ! nu or up
       ! Anti-Lepton + Neutrino
       if(WDK .eq. 1) then
          Wpolvec(1:4) = vbqq_Weyl(Dv,SpiNeu_up,SpiLep_dn)* coup
          return
       endif
       ! Anti-Down + Up
       if(WDK .eq. 3) then
          Wpolvec(1:4) = vbqq_Weyl(Dv,SpiNeu_up,SpiLep_dn)* coup * dsqrt(3.d0) * sq2 ! color-Factor dsqrt(3) *  Flavour  dsqrt(2)
          return
       endif

       ! Anti-Lepton + Neutrino + Photon
       ! Anti-Down + Up + Photon
       if(WDK .eq. 4 .or. WDK .eq. 2) then
!!!!!!!!!!!!!!! This ensures that we don't divide by zero
          If(dreal(klep_dn(1)) .gt. 0.d0 .and. dreal(kneu_up(1)) .gt. 0.d0) then
             call pol_mless(kpho,LamP,EpsP)      ! photon
         ! print *, "check of gauge invariance for photon in decay 1"
         ! EpsP(1:4)=dcmplx(kpho(1:4))

             If(WDK .eq. 4) then
                Q_f1 = Q_dn
                Q_f2 = Q_up
                coup = sq2 * dsqrt(3.d0)! FLAVOR * COLOR
             endif

             If(WDK .eq. 2) then
                Q_f1 = Q_el
                Q_f2 = 0.d0
                coup =1.d0
             endif

             coupl_sqrt = (0d0,-1d0)*dsqrt(alpha4Pi)
             !       diagram 1: photon emission off anti-down
             LepSpi(1:4) = spi2_Weyl(EpsP,SpiLep_dn) * coupl_sqrt*(-Q_f1)
             PropMom(1:4) = kpho(1:4) + klep_dn(1:4)
             LepSpi(1:4) = -spi2_Weyl(PropMom(1:4),LepSpi)*(0d0,1d0)/(2d0*(kpho(1:4).dot.klep_dn(1:4)))
             WCurr1(1:4)  = vbqq_Weyl(4,SpiNeu_up,LepSpi)  * g_weak ! vbqq introduces -i/Sqrt(2)

             !       diagram 2: photon emission off up
             NeuSpi(1:4) = spb2_Weyl(SpiNeu_up,EpsP) * coupl_sqrt*(-Q_f2)
             PropMom(1:4) = kpho(1:4) + kneu_up(1:4)
             NeuSpi(1:4) = spb2_Weyl(NeuSpi,PropMom(1:4))*(0d0,1d0)/(2d0*(kpho(1:4).dot.kneu_up(1:4)))
             WCurr2(1:4)  = vbqq_Weyl(4,NeuSpi,SpiLep_dn)  * g_weak ! vbqq introduces -i/Sqrt(2)

             !       diagram 3: photon emission off W+ boson
             WCurr3(1:4)  = vbqq_Weyl(4,SpiNeu_up,SpiLep_dn) * g_weak ! vbqq introduces -i/Sqrt(2)
             WMom(1:4) = klep_dn(1:4)+kneu_up(1:4)
             WCurr3(1:4) = Vertex_WPW(dreal(kpho(1:4)),dreal(WMom(1:4)),EpsP,WCurr3(1:4)) * (0d0,-1d0)/((WMom(1:4).dot.WMom(1:4))-m_W**2) *(-coupl_sqrt*Q_Wp)
             !       connect to quark current
             WCurr(1:4) = WCurr1(1:4) + WCurr2(1:4) + Wcurr3(1:4)
             WMom(1:4) = klep_dn(1:4)+kneu_up(1:4)+kpho(1:4)

             WCurr(1:4) = (-WCurr(1:4) + (WCurr(1:4).dot.WMom(1:4))/m_W**2 * WMom(1:4)) *coup ! Flavour and Color
!             check1 = dble(WCurr1.dot.WCurr1)
!             check2 = dble(WCurr2.dot.WCurr2)
!             check3 = dble(WCurr3.dot.WCurr3)
!             nan1 = IsNaN(check1)
!             nan2 = IsNaN(check2)
!             nan3 = IsNaN(check3)
!             If(nan1) then
!                write(51,*) "QW", QW
!                write(51,*) "Wcurr1", Wcurr1
!                write(51,*) "LepSpi", LepSpi
!                write(51,*) "PropMom", PropMom
!                write(51,*) "g_weak", g_weak
!                write(51,*) " vbqq_Weyl(4,SpiNeu_up,LepSpi)  * g_weak ", vbqq_Weyl(4,SpiNeu_up,LepSpi)  * g_weak
!                write(51,*) "(WCurr(1:4).dot.WMom(1:4))/m_W**2 * WMom(1:4), coup",(WCurr(1:4).dot.WMom(1:4))/m_W**2 * WMom(1:4), coup
!             endif
!             If(nan2) then
!                write(51,*) "QW", QW
!                write(51,*) "Wcurr2", Wcurr2
!                write(51,*) "LepSpi", LepSpi
!                write(51,*) "PropMom", PropMom
!                write(51,*) "g_weak", g_weak
!                write(51,*) "vbqq_Weyl(4,NeuSpi,SpiLep_dn)  * g_weak ",vbqq_Weyl(4,NeuSpi,SpiLep_dn)  * g_weak
!                write(51,*) "(WCurr(1:4).dot.WMom(1:4))/m_W**2 * WMom(1:4), coup",(WCurr(1:4).dot.WMom(1:4))/m_W**2 * WMom(1:4), coup
!             endif
!             If(nan3) then
!                write(51,*) "QW", QW
!                write(51,*) "Wcurr3", Wcurr3
!                write(51,*) "LepSpi", LepSpi
!                write(51,*) "PropMom", PropMom
!                write(51,*) "g_weak", g_weak
!                write(51,*) "Vertex_WPW(dreal(kpho(1:4)),dreal(WMom(1:4)),EpsP,WCurr3(1:4))",Vertex_WPW(dreal(kpho(1:4)),dreal(WMom(1:4)),EpsP,WCurr3(1:4))
!                write(51,*) "(WCurr(1:4).dot.WMom(1:4))/m_W**2 * WMom(1:4), coup",(WCurr(1:4).dot.WMom(1:4))/m_W**2 * WMom(1:4), coup
!             endif

          else
             Wcurr(1:4) = (0.d0,0.d0)
          endif
          Wpolvec(1:4) = Wcurr(1:4)
          return
       endif

    endif

    ! All W- cases
    If(QW .eq. -1) then
        call ubarSpi_Weyl(klep_dn,-1,SpiLep_dn)  ! l- or dn
        call    vSpi_Weyl(kneu_up,+1,SpiNeu_up)     ! nubar or up_bar

       !Lepton +  Anti-Neutrino
       if(WDK .eq. 1) then
          Wpolvec(1:4) = vbqq_Weyl(Dv,SpiLep_dn,SpiNeu_up)* coup
          return
       endif
       !Down +  Anti-Up
       if(WDK .eq. 3) then
          Wpolvec(1:4) =  vbqq_Weyl(Dv,SpiLep_dn,SpiNeu_up)* coup * dsqrt(3.d0) * sq2 ! color-Factor dsqrt(3) *  Flavour  dsqrt(2)
          return
       endif

       ! Lepton + Anti-Neutrino + Photon
       ! Down + Anti-Up + Photon
       if(WDK .eq. 4 .or. WDK .eq. 2) then
!!!!!!!!!!!!!!! This ensures that we don't divide by zero
          If(dreal(klep_dn(1)) .gt. 0.d0 .and. dreal(kneu_up(1)) .gt. 0.d0) then

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
             LepSpi(1:4) = spb2_Weyl(LepSpi(1:4),PropMom(1:4))*(0d0,1d0)/(2d0*(kpho(1:4).dot.klep_dn(1:4)))
             WCurr1(1:4)  = vbqq_Weyl(4,LepSpi(1:4),SpiNeu_up) * g_weak ! vbqq introduces -i/Sqrt(2)

             !       diagram 2: photon emission off anti-up
             NeuSpi(1:4) = spi2_Weyl(EpsP(1:4),SpiNeu_up) * coupl_sqrt*(Q_f2)
             PropMom(1:4) = kpho(1:4) + kneu_up(1:4)
             NeuSpi(1:4) = -spi2_Weyl(PropMom(1:4),NeuSpi(1:4))*(0d0,1d0)/(2d0*(kpho(1:4).dot.kneu_up(1:4)))
             WCurr2(1:4)  = vbqq_Weyl(4,Spilep_dn,NeuSpi) * g_weak ! vbqq introduces -i/Sqrt(2)

             !       diagram 3: photon emission off W- boson!                       MINUS INTRODUCED TO ENSURE GAUGE INVARIANCE!
             WCurr3(1:4)  = vbqq_Weyl(4,SpiLep_dn,SpiNeu_up)  * g_weak ! vbqq introduces -i/Sqrt(2)
             WMom(1:4) = klep_dn(1:4)+kneu_up(1:4)
             WCurr3(1:4) = Vertex_WPW(dreal(kpho(1:4)),dreal(WMom(1:4)),EpsP(1:4),WCurr3(1:4)) *(0d0,-1d0)/((WMom(1:4).dot.WMom(1:4))-m_W**2) *(-coupl_sqrt*Q_Wm)   *(-1d0)
             !       connect to quark current
             WCurr(1:4) = (WCurr1(1:4) + WCurr2(1:4) + Wcurr3(1:4))
             WMom(1:4) = klep_dn(1:4)+kneu_up(1:4)+kpho(1:4)

             WCurr(1:4) = (-WCurr(1:4) + (WCurr(1:4).dot.WMom(1:4))/m_W**2 * WMom(1:4))* coup ! Flavour & Color
!             check1 = dble(WCurr1.dot.WCurr1)
!             check2 = dble(WCurr2.dot.WCurr2)
!             check3 = dble(WCurr3.dot.WCurr3)
!             nan1 = IsNaN(check1)
!             nan2 = IsNaN(check2)
!             nan3 = IsNaN(check3)
!             If(nan1) then
!                write(51,*) "QW", QW
!                write(51,*) "Wcurr1", Wcurr1
!                write(51,*) "LepSpi", LepSpi
!                write(51,*) "PropMom", PropMom
!                write(51,*) "g_weak", g_weak
!                write(51,*) " vbqq_Weyl(4,SpiNeu_up,LepSpi)  * g_weak ", vbqq_Weyl(4,SpiNeu_up,LepSpi)  * g_weak
!                write(51,*) "(WCurr(1:4).dot.WMom(1:4))/m_W**2 * WMom(1:4), coup",(WCurr(1:4).dot.WMom(1:4))/m_W**2 * WMom(1:4), coup
!             endif
!             If(nan2) then
!                write(51,*) "QW", QW
!                write(51,*) "Wcurr2", Wcurr2
!                write(51,*) "LepSpi", LepSpi
!                write(51,*) "PropMom", PropMom
!                write(51,*) "g_weak", g_weak
!                write(51,*) "vbqq_Weyl(4,NeuSpi,SpiLep_dn)  * g_weak ",vbqq_Weyl(4,NeuSpi,SpiLep_dn)  * g_weak
!                write(51,*) "(WCurr(1:4).dot.WMom(1:4))/m_W**2 * WMom(1:4), coup",(WCurr(1:4).dot.WMom(1:4))/m_W**2 * WMom(1:4), coup
!             endif
!             If(nan3) then
!                write(51,*) "QW", QW
!                write(51,*) "Wcurr3", Wcurr3
!                write(51,*) "LepSpi", LepSpi
!                write(51,*) "PropMom", PropMom
!                write(51,*) "g_weak", g_weak
!                write(51,*) "Vertex_WPW(dreal(kpho(1:4)),dreal(WMom(1:4)),EpsP,WCurr3(1:4))",Vertex_WPW(dreal(kpho(1:4)),dreal(WMom(1:4)),EpsP,WCurr3(1:4))
!                write(51,*) "(WCurr(1:4).dot.WMom(1:4))/m_W**2 * WMom(1:4), coup",(WCurr(1:4).dot.WMom(1:4))/m_W**2 * WMom(1:4), coup
!             endif
!
          else
             Wcurr(1:4) = (0.d0,0.d0)
          endif

          Wpolvec(1:4) = Wcurr(1:4)

          return
       endif



    endif

  end function Wpolvec




END MODULE
