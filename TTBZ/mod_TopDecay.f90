MODULE ModTopDecay
implicit none
save

public :: TopDecay

contains


SUBROUTINE TopDecay(TopQuark,Topol,Mom,GluonHel,Gluon2Hel,PhotonHel,xIntDip,PartAmp)
use ModMisc
use ModProcess
use ModAmplitudes
use ModParameters
use ModHadrWDecay
use ModJPsiFrag
use ModTTBP_NLODK
use ModTTBJ_NLODK
use ModTTBJ_NLODKW
use ModDKIntDipoles
use ModWDecay
implicit none
type(Particle) :: TopQuark
integer,optional :: GluonHel,Gluon2Hel,PhotonHel,PartAmp
real(8),optional :: xIntDip(1:2)
integer :: Topol,i,Wp_DKmode,Wm_DKmode,NMom
integer,parameter :: WPlus=+1, WMinus=-1
real(8) :: Mom(:,:),PSWgt,WMom(1:4),xx,z,Cdip,C_ct,FF_xxz,FF_xx,MomT(1:4)
complex(8) :: BotSpi(1:4),Spi(1:4),BarSpi(1:4),LepSpi(1:4),GluPol(1:4),WPol(1:4),PhoPol(1:4)
complex(8) :: WCurr_new(1:4),WCurr(1:4),F0,F1,IntDip,TopProp,PropMom(1:4),WCurr1(1:4),WCurr2(1:4)
complex(8) :: SME(1:8,1:4)
real(8),parameter :: CF=4d0/3d0,CA=3d0,TR=1d0/2d0
! real(8),parameter :: g2_weak = 4d0*DblPi*alpha/sw2
!real(8),parameter :: g2_weak = 4d0*dsqrt(2d0)*m_W**2*GF
!real(8),parameter :: g_weak = dsqrt(g2_weak)
real(8),parameter :: Nc=3d0,Nflight=5d0,NFlav=2d0
real(8), parameter :: eta = 0  !1=tHV scheme, 0=FDH
complex(8) :: coupl_sqrt
real(8) :: NWAFactor_Top
real(8) :: NWAFactor_W
complex(8) :: WProp
real(8),parameter :: pisqo6 = DblPi**2/6d0
real(8) :: rsq,omrsq,wlog,rlog,Kfun,epinv,epinv2
real(8) :: s23,s24,s34,zeros(1:6)
complex(8) :: c3p2,c3p4,c4p3,c3c4
logical, save ::  init(-1:19)=.true.
type(TreeProcess),save :: ATopDKAmp(0:4),TopDKAmp(0:4)! 0:tbW, 1:tbWg, 2,3:tbWgg, 4:tbWqqb
type(Particle),save :: TopDKProd(1:7)
type(Particle),save :: ATopDKProd(1:7)
real(8), parameter :: Tt_x_Tb = 1d0/2d0*CA-CF
real(8), parameter :: Tt_x_Tg =-1d0/2d0*CA
real(8), parameter :: Tb_x_Tg =-1d0/2d0*CA
!-------------------------------------------------------KM definitions
real(8) :: ep1,ep2
real(8) :: v_factor, int_dip_factor
!-------------------------------------------------------SCHARFS inputs BEGIN
integer :: Top_Atop, LO_NLO, heli_virtuals(2), heli_reals(3), hel_pho, hel_glu, hel_bot, gau, xe, Wpar(1:4)
complex(8) :: SpiOut(1:4)
real(8) :: kb(4),kl(4),kn(4),kp(4)
!-------------------------------------------------------SCHARFS inputs END


!DEC$ IF(_CheckMomenta .EQ.1)
   NMom = size(Mom(:,:),2)
   zeros(1:4) = dble(TopQuark%Mom(1:4))
   do i=1,NMom
      zeros(1:4) = zeros(1:4) - Mom(1:4,i)
   enddo
   if( any(abs(zeros(1:4)/dble(TopQuark%Mom(1))).gt.1d-8) ) then
      print *, "ERROR: energy-momentum violation in SUBROUTINE TopDecay(): ",TopQuark%PartType,zeros(1:4)
      print *, "momentum dump:"
      print *, Mom(:,:)
   endif

   zeros(1) = dble(TopQuark%Mom(1:4).dot.TopQuark%Mom(1:4)) - m_Top**2
   do i=1,NMom
      zeros(i+1)=  Mom(1:4,i).dot.Mom(1:4,i)
   enddo
   if( any(abs(zeros(1:6)/dble(TopQuark%Mom(1))**2).gt.1d-5) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE TopDecay(): ",TopQuark%PartType,zeros(1:6)
      print *, "momentum dump:"
      print *, Mom(:,:)
   endif
!DEC$ ENDIF



if( abs(TOPDECAYS).eq.1 ) then! di-leptonic decay
    Wp_DKmode = WDK_Lep
    Wm_DKmode = WDK_Lep
elseif( TOPDECAYS.eq.2 ) then! fully hadronic decay
    Wp_DKmode = WDK_Had
    Wm_DKmode = WDK_Had
elseif( TOPDECAYS.eq.3 .or. TOPDECAYS.eq.5 ) then! Top --> hadr.
    Wp_DKmode = WDK_Had
    Wm_DKmode = WDK_Lep
elseif( TOPDECAYS.eq.4 .or. TOPDECAYS.eq.6 ) then! ATop --> hadr.
    Wp_DKmode = WDK_Lep
    Wm_DKmode = WDK_Had
endif


if( init(-1) ) then
    TopDKProd(1)%PartType = Top_
    TopDKProd(1)%ExtRef   = 1
    TopDKProd(1)%Mass = M_Top
    TopDKProd(1)%Mass2= M_Top**2

    TopDKProd(2)%PartType = Str_
    TopDKProd(2)%ExtRef   = 2
    TopDKProd(2)%Mass = 0d0
    TopDKProd(2)%Mass2= 0d0

    TopDKProd(3)%PartType = Wp_
    TopDKProd(3)%ExtRef   = 3
    TopDKProd(3)%Mass = M_W
    TopDKProd(3)%Mass2= M_W**2

    TopDKProd(4)%PartType = Glu_
    TopDKProd(4)%ExtRef   = 4
    TopDKProd(4)%Mass = 0d0
    TopDKProd(4)%Mass2= 0d0

    TopDKProd(5)%PartType = Glu_
    TopDKProd(5)%ExtRef   = 5
    TopDKProd(5)%Mass = 0d0
    TopDKProd(5)%Mass2= 0d0

    TopDKProd(6)%PartType = Str_
    TopDKProd(6)%ExtRef   = 6
    TopDKProd(6)%Mass = 0d0
    TopDKProd(6)%Mass2= 0d0

    TopDKProd(7)%PartType = AStr_
    TopDKProd(7)%ExtRef   = 7
    TopDKProd(7)%Mass = 0d0
    TopDKProd(7)%Mass2= 0d0

    ATopDKProd(1)%PartType = ATop_
    ATopDKProd(1)%ExtRef   = 1
    ATopDKProd(1)%Mass = M_Top
    ATopDKProd(1)%Mass2= M_Top**2

    ATopDKProd(2)%PartType = AStr_
    ATopDKProd(2)%ExtRef   = 2
    ATopDKProd(2)%Mass = 0d0
    ATopDKProd(2)%Mass2= 0d0

    ATopDKProd(3)%PartType = Wm_
    ATopDKProd(3)%ExtRef   = 3
    ATopDKProd(3)%Mass = M_W
    ATopDKProd(3)%Mass2= M_W**2

    ATopDKProd(4)%PartType = Glu_
    ATopDKProd(4)%ExtRef   = 4
    ATopDKProd(4)%Mass = 0d0
    ATopDKProd(4)%Mass2= 0d0

    ATopDKProd(5)%PartType = Glu_
    ATopDKProd(5)%ExtRef   = 5
    ATopDKProd(5)%Mass = 0d0
    ATopDKProd(5)%Mass2= 0d0

    ATopDKProd(6)%PartType = Str_
    ATopDKProd(6)%ExtRef   = 6
    ATopDKProd(6)%Mass = 0d0
    ATopDKProd(6)%Mass2= 0d0

    ATopDKProd(7)%PartType = AStr_
    ATopDKProd(7)%ExtRef   = 7
    ATopDKProd(7)%Mass = 0d0
    ATopDKProd(7)%Mass2= 0d0

    init(-1)=.false.
endif


NWAFactor_Top = 1d0/dsqrt(2d0*Ga_Top(0)*m_Top)
NWAFactor_W   = 1d0/dsqrt(2d0*Ga_W(0)*m_W)
WProp = (0d0,-1d0)*NWAFactor_W	

!-----------------------------------------------------------------------------
IF( Topol.eq.DK_LO ) THEN! leading order
!-----------------------------------------------------------------------------
    if( init(DK_LO) ) then
        allocate( TopDKAmp(0)%NumGlu(0:2))
        allocate( TopDKAmp(0)%PartRef(1:3))
        allocate( TopDKAmp(0)%PartType(1:3))
        allocate( TopDKAmp(0)%Quarks(1:2))
        allocate( TopDKAmp(0)%Gluons(1:0))
        allocate( ATopDKAmp(0)%NumGlu(0:2))
        allocate( ATopDKAmp(0)%PartRef(1:3))
        allocate( ATopDKAmp(0)%PartType(1:3))
        allocate( ATopDKAmp(0)%Quarks(1:2))
        allocate( ATopDKAmp(0)%Gluons(1:0))

        TopDKAmp(0)%NumPart =3; TopDKAmp(0)%NumQua =2; TopDKAmp(0)%NumGlu(0) = 0
        ATopDKAmp(0)%NumPart=3; ATopDKAmp(0)%NumQua=2; ATopDKAmp(0)%NumGlu(0)= 0
        TopDKAmp(0)%PartRef(1:3) = (/1,2,3/)
        ATopDKAmp(0)%PartRef(1:3)= (/1,2,3/)
        call LinkTreeParticles( TopDKAmp(0), TopDKProd(1:3))
        call LinkTreeParticles(ATopDKAmp(0),ATopDKProd(1:3))
        init(DK_LO)=.false.
    endif
    if( TopQuark%PartType.eq.Top_ ) then ! Top quark decay
        TopDKProd(1)%Mom(1:4) =-TopQuark%Mom(1:4)!  top
        TopDKProd(2)%Mom(1:4) = dcmplx(Mom(1:4,1))! bot
        TopDKProd(3)%Mom(1:4) = dcmplx(Mom(1:4,2)+Mom(1:4,3))!  W+
        call ubarSpi_Weyl(TopDKProd(2)%Mom(1:4),-1,TopDKProd(2)%Pol(1:4))
        TopDKProd(3)%Pol(1:4) = WPolVec(WPlus,Wp_DKmode,Mom(1:4,2:3),WDK_LO)* WProp * g_weak/dsqrt(2d0)
        call new_calc_ampl(4,4,0,0,TopDKAmp(0),TopQuark%Pol(1:4))
        TopQuark%Pol(1:4) =( spb2_Weyl(TopQuark%Pol(1:4),TopQuark%Mom(1:4)) + m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)

    elseif( TopQuark%PartType.eq.ATop_ ) then  ! Anti-Top quark decay
        ATopDKProd(1)%Mom(1:4) =-TopQuark%Mom(1:4)!  Atop
        ATopDKProd(2)%Mom(1:4) = dcmplx(Mom(1:4,1))! Abot
        ATopDKProd(3)%Mom(1:4) = dcmplx(Mom(1:4,2)+Mom(1:4,3))!  W-
        call vSpi_Weyl(ATopDKProd(2)%Mom(1:4),+1,ATopDKProd(2)%Pol(1:4))
        ATopDKProd(3)%Pol(1:4) = WPolVec(WMinus,Wm_DKmode,Mom(1:4,2:3),WDK_LO)* WProp * g_weak/dsqrt(2d0)
        call new_calc_ampl(4,4,0,0,ATopDKAmp(0),TopQuark%Pol(1:4))

        TopQuark%Pol(1:4) = ( spi2_Weyl(TopQuark%Mom(1:4),TopQuark%Pol(1:4)) - m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)
    endif

!-----------------------------------------------------------------------------
ELSEIF( Topol.eq.DK_1L_T ) THEN! 1-loop correction to top quark line + int.dipoles
!-----------------------------------------------------------------------------
!(formulae taken from Single Top Production in MCFM)
!epinv = (4*pi*musq/mtsq)^ep/Gamma(1-ep)/ep
    rsq = 1d0-2d0*(TopQuark%Mom(1:4).dot.Mom(1:4,1))/m_Top**2  ! = m_W^2 / m_Top^2 for on-shell W boson
    omrsq=1d0-rsq
    wlog=dlog(omrsq)
    rlog=dlog(rsq)
    Kfun=1d0/rsq*wlog
    epinv2 = 0d0
    epinv  = 0d0
    if( rsq.gt.1d0 ) call Error("r2>1 in virtual correction to the top decay")
!   virtual corrections
    F0 = - epinv2 -epinv*(2.5d0-2d0*wlog)  &
         - 5.5d0 - pisqo6 - 2d0*DLi2(rsq)  &  ! -5.5 corresponds to eta=0 (FHD scheme)
         + 3d0*wlog - 2d0*wlog**2 - Kfun
    F1 = 2d0*Kfun
!   integrated dipoles
    if( (TopDecays.eq.6 .and. TopQuark%PartType.eq.Top_) .or. (TopDecays.eq.5 .and. TopQuark%PartType.eq.ATop_) ) then
        xx = xIntDip(1)
        z  = xIntDip(2)
        FF_xxz = FF(xx/z)
        FF_xx  = FF(xx)
        Cdip = (epinv2-2d0*wlog*epinv+2d0*wlog**2 + 2d0-DblPi**2/6.0)*FF_xx  &
             - (epinv-2d0*wlog)*( 2d0/(1d0-z)*(FF_xxz/z - FF_xx)       &
             - (1d0+z)/z*FF_xxz -  FF_xx)                              &
             +  eta*(1d0-z)*FF_xxz/z                                    &
             + 4d0*dlog(1d0-z)/(1d0-z)*(FF_xxz/z - FF_xx)              &
             + FF_xxz/z*(-2d0*(1d0+z)*dlog(1d0-z)                       &
             - (2d0 /(1d0-z)-(1d0+z))*dlog((rsq+z*(1d0-rsq))/z**2) )      &
             - 2d0/(1d0-z)*(z/(rsq+z*(1d0-rsq))*FF_xxz/z - FF_xx)

        C_ct = (epinv+dlog(m_top**2/MuFrag**2))*( FF_xxz/z*(2d0/(1d0-z) - 1d0-z) &
             - FF_xx*(2d0/(1d0-z)- 3d0/2d0) )
        IntDip = 1d0*( Cdip + C_ct )
        F0 = F0 * FF_xx
        F1 = F1 * FF_xx
    else
        IntDip = epinv2 + epinv*(2.5d0-2d0*wlog)  &
           + 25d0/4d0 + 0.5d0*(1d0/omrsq**2-8d0/omrsq+7d0)*rlog  &
           + 0.5d0/omrsq + 2d0*DLi2(omrsq) - 5d0*pisqo6  &
           - 5d0*wlog + 2d0*wlog**2 &!    eta=0 (FHD scheme)
           - ( 2d0*dlog(alpha_DKTfi)**2 -dlog(alpha_DKTfi)*(-7d0/2d0+4d0*alpha_DKTfi-alpha_DKTfi**2/2d0) &
           - 2d0*(1d0-alpha_DKTfi)*rsq/(omrsq)*dlog(rsq/(1d0-alpha_DKTfi+rsq*alpha_DKTfi)) )
    endif


    if( TopQuark%PartType.eq.Top_ ) then ! Top quark decay
        call ubarSpi_Weyl(dcmplx(Mom(1:4,1)),-1,BotSpi(1:4))
        WCurr(1:4) = WPolVec(WPlus,Wp_DKmode,Mom(1:4,2:3),WDK_LO)* WProp * g_weak
!       connect to quark current
        TopQuark%Pol(1:4) = (F0+IntDip)*vgq_Weyl( WCurr(1:4),BotSpi(1:4) ) ! vgq introduces -i/Sqrt(2)
        TopQuark%Pol(1:4) = TopQuark%Pol(1:4) + BotSpi(1:4)*F1*(0d0,-1d0)/dsqrt(2d0)*sc_(dcmplx(Mom(1:4,1)),WCurr(1:4))/m_Top ! introduce -i/Sqrt(2) by hand
        TopQuark%Pol(1:4) =( spb2_Weyl(TopQuark%Pol(1:4),TopQuark%Mom(1:4)) + m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top * alpha_sOver2Pi*RunAlphaS(NLOParam,MuRen) * CF
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)

    elseif( TopQuark%PartType.eq.ATop_ ) then  ! Anti-Top quark decay
        call    vSpi_Weyl(dcmplx(Mom(1:4,1)),+1,BotSpi(1:4))
        WCurr(1:4) = WPolVec(WMinus,Wm_DKmode,Mom(1:4,2:3),WDK_LO)* WProp * g_weak
!       connect to quark current:
        TopQuark%Pol(1:4) = (F0+IntDip)*vbqg_Weyl( BotSpi(1:4),WCurr(1:4) )! vbqg introduces -i/Sqrt(2)
        ! INTORODUCED A MINUS SIGN FOR F1 CONTRIBUTION
        TopQuark%Pol(1:4) = TopQuark%Pol(1:4) - BotSpi(1:4)*F1*(0d0,-1d0)/dsqrt(2d0)*sc_(dcmplx(Mom(1:4,1)),WCurr(1:4))/m_Top ! introduce -i/Sqrt(2) by hand
        TopQuark%Pol(1:4) = ( spi2_Weyl(TopQuark%Mom(1:4),TopQuark%Pol(1:4)) - m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top * alpha_sOver2Pi*RunAlphaS(NLOParam,MuRen) * CF
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)
    endif

!-----------------------------------------------------------------------------
ELSEIF( Topol.eq.DK_RE_T ) THEN! real emission from top quark line!   DK_RE_T=DKJ_LO_T
!-----------------------------------------------------------------------------
    if( init(DK_RE_T) ) then
        allocate( TopDKAmp(1)%NumGlu(0:2))
        allocate( TopDKAmp(1)%PartRef(1:4))
        allocate( TopDKAmp(1)%PartType(1:4))
        allocate( TopDKAmp(1)%Quarks(1:2))
        allocate( TopDKAmp(1)%Gluons(1:1))
        allocate( ATopDKAmp(1)%NumGlu(0:2))
        allocate( ATopDKAmp(1)%PartRef(1:4))
        allocate( ATopDKAmp(1)%PartType(1:4))
        allocate( ATopDKAmp(1)%Quarks(1:2))
        allocate( ATopDKAmp(1)%Gluons(1:1))

        TopDKAmp(1)%NumPart=4;   TopDKAmp(1)%NumQua=2; TopDKAmp(1)%NumGlu(0) = 1
        ATopDKAmp(1)%NumPart=4; ATopDKAmp(1)%NumQua=2; ATopDKAmp(1)%NumGlu(0)= 1
        TopDKAmp(1)%PartRef(1:4) = (/1,2,4,3/)
        ATopDKAmp(1)%PartRef(1:4)= (/1,2,4,3/)
        call LinkTreeParticles( TopDKAmp(1), TopDKProd(1:4))
        call LinkTreeParticles(ATopDKAmp(1),ATopDKProd(1:4))
        init(DK_RE_T)=.false.
    endif

    coupl_sqrt = dsqrt(alpha_s4Pi*RunAlphaS(NLOParam,MuRen)*CF)*dsqrt(2d0)  ! dsqrt(2d0) factor is to compensate for 1/dsqrt(2d0) factors in gluon coupling
    if( TopQuark%PartType.eq.Top_ ) then ! Top quark decay   !  has a minus sign wrt. to older version
        TopDKProd(1)%Mom(1:4) =-TopQuark%Mom(1:4)!  top
        TopDKProd(2)%Mom(1:4) = dcmplx(Mom(1:4,1))! bot
        TopDKProd(3)%Mom(1:4) = dcmplx(Mom(1:4,2)+Mom(1:4,3))!  W+
        TopDKProd(4)%Mom(1:4) = dcmplx(Mom(1:4,4))! glu
        call ubarSpi_Weyl(TopDKProd(2)%Mom(1:4),-1,TopDKProd(2)%Pol(1:4))
        TopDKProd(3)%Pol(1:4) = WPolVec(WPlus,Wp_DKmode,Mom(1:4,2:3),WDK_LO)* WProp * g_weak/dsqrt(2d0)
        call pol_mless(TopDKProd(4)%Mom(1:4),GluonHel,TopDKProd(4)%Pol(1:4))
        call new_calc_ampl(4,4,0,0,TopDKAmp(1),TopQuark%Pol(1:4))

        TopQuark%Pol(1:4) =( spb2_Weyl(TopQuark%Pol(1:4),TopQuark%Mom(1:4)) + m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top * coupl_sqrt

        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)

    elseif( TopQuark%PartType.eq.ATop_ ) then  ! Anti-Top quark decay
        ATopDKProd(1)%Mom(1:4) =-TopQuark%Mom(1:4)!  top
        ATopDKProd(2)%Mom(1:4) = dcmplx(Mom(1:4,1))! bot
        ATopDKProd(3)%Mom(1:4) = dcmplx(Mom(1:4,2)+Mom(1:4,3))!  W+
        ATopDKProd(4)%Mom(1:4) = dcmplx(Mom(1:4,4))! glu
        call vSpi_Weyl(ATopDKProd(2)%Mom(1:4),+1,ATopDKProd(2)%Pol(1:4))
        ATopDKProd(3)%Pol(1:4) = WPolVec(WMinus,Wm_DKmode,Mom(1:4,2:3),WDK_LO)* WProp * g_weak/dsqrt(2d0)
        call pol_mless(ATopDKProd(4)%Mom(1:4),GluonHel,ATopDKProd(4)%Pol(1:4))
        call new_calc_ampl(4,4,0,0,ATopDKAmp(1),TopQuark%Pol(1:4))

        TopQuark%Pol(1:4) =( spi2_Weyl(TopQuark%Mom(1:4),TopQuark%Pol(1:4)) - m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top * coupl_sqrt

        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)
    endif


!-----------------------------------------------------------------------------
ELSEIF( Topol.eq.DK_1L_Q ) THEN! 1-loop correction to light quark line + int.dipoles
!-----------------------------------------------------------------------------
!--- definitions that are required
     ep2 = 0d0
     ep1 = 0d0

     v_factor = alpha_sOver2Pi*RunAlphaS(NLOParam,MuRen)*CF*(-2d0*ep2-3d0*ep1-8d0+pisq)

!---- first factor ``2'' is the number of dipoles
      int_dip_factor = 2d0*alpha_sOver2Pi*RunAlphaS(NLOParam,MuRen)*CF*(  ep2 +3d0/2d0*ep1+5d0-pisq/2d0  &
                     + 3d0/2d0*(alpha_DKWff - 1d0 -dlog(alpha_DKWff))-dlog(alpha_DKWff)**2  )

    if( TopQuark%PartType.eq.Top_ ) then ! Top quark decay
!       assemble lepton current
        call ubarSpi_Weyl(dcmplx(Mom(1:4,1)),-1,BotSpi(1:4))  ! bot
!         call    vSpi_Weyl(dcmplx(Mom(1:4,2)),+1,Spi(1:4))     ! dn_bar
!         call ubarSpi_Weyl(dcmplx(Mom(1:4,3)),-1,BarSpi(1:4))  ! up
!         call ampl_w_qbq(Mom(1:4,2:3),Wcurr)
!         Wcurr = Wcurr*WProp*g2_weak/sqrt2
!         if( TOPDECAYS.EQ.2 .or. TOPDECAYS.EQ.3 .or. TOPDECAYS.EQ.4 ) Wcurr(1:4) = Wcurr(1:4) * dsqrt(6d0)
        Wcurr = WPolVec(WPlus,Wp_DKmode,Mom(1:4,2:3),WDK_LO)  * WProp * g_weak
        Wcurr = (v_factor + int_dip_factor)*Wcurr

!       connect to quark current
        BarSpi(1:4) = BotSpi(1:4)
        TopQuark%Pol(1:4) = vgq_Weyl( WCurr(1:4),BarSpi(1:4) ) ! vgq introduces -i/Sqrt(2)
        TopQuark%Pol(1:4) =( spb2_Weyl(TopQuark%Pol(1:4),TopQuark%Mom(1:4)) + m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)

    elseif( TopQuark%PartType.eq.ATop_ ) then  ! Anti-Top quark decay
!       assemble lepton current
        call    vSpi_Weyl(dcmplx(Mom(1:4,1)),+1,BotSpi(1:4))  ! Abot
!         call ubarSpi_Weyl(dcmplx(Mom(1:4,2)),-1,BarSpi(1:4))  ! dn
!         call    vSpi_Weyl(dcmplx(Mom(1:4,3)),+1,Spi(1:4))     ! up_bar
!         call ampl_w_qbq((/Mom(1:4,3),Mom(1:4,2)/),Wcurr)
!         if( TOPDECAYS.EQ.2 .or. TOPDECAYS.EQ.3 .or. TOPDECAYS.EQ.4 ) Wcurr(1:4) = Wcurr(1:4) * dsqrt(6d0)
!         Wcurr = Wcurr*WProp*g2_weak/sqrt2
        Wcurr = WPolVec(WMinus,Wm_DKmode,Mom(1:4,2:3),WDK_LO)  * WProp * g_weak
        Wcurr = (v_factor + int_dip_factor)*Wcurr

!       connect to quark current:
        Spi(1:4) = BotSpi(1:4)
        TopQuark%Pol(1:4) = vbqg_Weyl( Spi(1:4),WCurr(1:4) )! vbqg introduces -i/Sqrt(2)
        TopQuark%Pol(1:4) = ( spi2_Weyl(TopQuark%Mom(1:4),TopQuark%Pol(1:4)) - m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)
    endif


!-----------------------------------------------------------------------------
ELSEIF( Topol.eq.DK_RE_Q ) THEN! real emission from light quark line in W decay! DK_RE_Q=DKJ_LO_Q
!-----------------------------------------------------------------------------
    if( TOPDECAYS.EQ.1 ) call Error("this topolgy is not allowed for TOPDECAYS=1")

    if( TopQuark%PartType.eq.Top_ ) then ! Top quark decay
!       assemble lepton current
        call ubarSpi_Weyl(dcmplx(Mom(1:4,1)),-1,BotSpi(1:4))  ! bot
!         call    vSpi_Weyl(dcmplx(Mom(1:4,2)),+1,Spi(1:4))     ! dn_bar
!         call ubarSpi_Weyl(dcmplx(Mom(1:4,3)),-1,BarSpi(1:4))  ! up
        call ampl_w_qbqg((/Mom(1:4,2),Mom(1:4,3),Mom(1:4,4)/),GluonHel,Wcurr(1:4))
        Wcurr(1:4) = WProp*g2_weak/sqrt2 * dsqrt(alpha_s4Pi*CF*RunAlphaS(NLOParam,MuRen))* dsqrt(NFlav*Nc) * Wcurr(1:4)

!       connect to quark current
        BarSpi(1:4) = BotSpi(1:4)
        TopQuark%Pol(1:4) = vgq_Weyl( WCurr(1:4),BarSpi(1:4) ) ! vgq introduces -i/Sqrt(2)
        TopQuark%Pol(1:4) =( spb2_Weyl(TopQuark%Pol(1:4),TopQuark%Mom(1:4)) + m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)

    elseif( TopQuark%PartType.eq.ATop_ ) then  ! Anti-Top quark decay
!       assemble lepton current
        call    vSpi_Weyl(dcmplx(Mom(1:4,1)),+1,BotSpi(1:4))  ! Abot
!         call ubarSpi_Weyl(dcmplx(Mom(1:4,2)),-1,BarSpi(1:4))  ! dn
!         call    vSpi_Weyl(dcmplx(Mom(1:4,3)),+1,Spi(1:4))     ! up_bar
        call ampl_w_qbqg((/Mom(1:4,3),Mom(1:4,2),Mom(1:4,4)/),GluonHel,Wcurr(1:4))
        Wcurr(1:4) = WProp*g2_weak/sqrt2 * dsqrt(alpha_s4Pi*CF*RunAlphaS(NLOParam,MuRen))* dsqrt(NFlav*Nc) * Wcurr(1:4)
!       connect to quark current:
        Spi(1:4) = BotSpi(1:4)
        TopQuark%Pol(1:4) = vbqg_Weyl( Spi(1:4),WCurr(1:4) )! vbqg introduces -i/Sqrt(2)
        TopQuark%Pol(1:4) = ( spi2_Weyl(TopQuark%Mom(1:4),TopQuark%Pol(1:4)) - m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)
    endif


!-----------------------------------------------------------------------------
ELSEIF( Topol.eq.DKP_LO_T ) THEN! photon emission from top quark line
!-----------------------------------------------------------------------------
    coupl_sqrt = (0d0,-1d0)*dsqrt(alpha4Pi)!   == -ie
    if( TopQuark%PartType.eq.Top_ ) then ! Top quark decay
        call ubarSpi_Weyl(dcmplx(Mom(1:4,1)),-1,BotSpi(1:4))  ! bot
        call pol_mless(dcmplx(Mom(1:4,4)),PhotonHel,PhoPol(1:4))        ! photon
        WCurr(1:4) = WPolVec(WPlus,Wp_DKmode,Mom(1:4,2:3),WDK_LO)* WProp * g_weak
!         WCurr(1:4) = WCurr_new(1:4)
!DEC$ IF(_CheckWPolVec .EQ.1)
        call    vSpi_Weyl(dcmplx(Mom(1:4,2)),+1,Spi(1:4))     ! l+ or dn_bar
        call ubarSpi_Weyl(dcmplx(Mom(1:4,3)),-1,BarSpi(1:4))  ! nu or up
        WCurr(1:4)  = vbqq_Weyl(4,BarSpi(1:4),Spi(1:4)) * WProp * g2_weak ! vbqq introduces -i/Sqrt(2)
!         print *, "check DKP_LO_T1 WPolVec routine:"
!         print *, WCurr(1:4) - WCurr_new(1:4)
!         pause
!DEC$ ENDIF
!         print *, "check of gauge invariance for photon in decay"
!         PhoPol(1:4)=dcmplx(Mom(1:4,4))

!       connect to quark current: diagram 1: photon emission off top quark
        BarSpi(1:4) = vgq_Weyl( WCurr(1:4),BotSpi(1:4) ) ! vgq introduces -i/Sqrt(2)

        PropMom(1:4) = TopQuark%Mom(1:4) - dcmplx(Mom(1:4,4))
!         TopProp = sc_(PropMom(1:4),PropMom(1:4))-m_Top**2
        TopProp = -2d0*(TopQuark%Mom(1:4).dot.Mom(1:4,4))
        BarSpi(1:4) = (0d0,1d0)*( spb2_Weyl(BarSpi(1:4),PropMom(1:4)) + m_Top*BarSpi(1:4) )/TopProp  ! top propagator

        BarSpi(1:4) = spb2_Weyl(BarSpi(1:4),PhoPol(1:4)) *coupl_sqrt*Q_top
        TopQuark%Pol(1:4) =( spb2_Weyl(BarSpi(1:4),TopQuark%Mom(1:4)) + m_Top*BarSpi(1:4) ) * NWAFactor_Top

!       adding diagram 2: photon emission off bottom quark
        BarSpi(1:4) = spb2_Weyl(BotSpi(1:4),PhoPol(1:4)) *coupl_sqrt*Q_dn

        PropMom(1:4) = dcmplx( Mom(1:4,1)+Mom(1:4,4) )
        TopProp = sc_(PropMom(1:4),PropMom(1:4))
        BarSpi(1:4) = (0d0,1d0)*( spb2_Weyl(BarSpi(1:4),PropMom(1:4)) )/TopProp  ! bot propagator

        BarSpi(1:4) = vgq_Weyl( WCurr(1:4),BarSpi(1:4) )  ! vgq introduces -i/Sqrt(2)
        TopQuark%Pol(1:4) = TopQuark%Pol(1:4) + ( spb2_Weyl(BarSpi(1:4),TopQuark%Mom(1:4)) + m_Top*BarSpi(1:4) ) * NWAFactor_Top

!       adding diagram 3: photon emission off W+ boson
        WMom(1:4) = Mom(1:4,2)+Mom(1:4,3)
        PropMom(1:4) = WMom(1:4) + Mom(1:4,4)

        WCurr1(1:4) = Vertex_WPW(Mom(1:4,4),WMom(1:4),PhoPol(1:4),WCurr(1:4)) *(-coupl_sqrt)*( Q_top-Q_dn )! = Q_Wp
        WCurr2(1:4) = -WCurr1(1:4) + (WCurr1(1:4).dot.PropMom(1:4))/m_W**2*PropMom(1:4)

        BarSpi(1:4) = vgq_Weyl( WCurr2(1:4),BotSpi(1:4) ) * (0d0,-1d0)/(2d0*(WMom(1:4).dot.Mom(1:4,4))) ! vgq introduces -i/Sqrt(2)
        TopQuark%Pol(1:4) = TopQuark%Pol(1:4) +( spb2_Weyl(BarSpi(1:4),TopQuark%Mom(1:4)) + m_Top*BarSpi(1:4) ) * NWAFactor_Top

!       convert Weyl to Dirac
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)



    elseif( TopQuark%PartType.eq.ATop_ ) then  ! Anti-Top quark decay
        call    vSpi_Weyl(dcmplx(Mom(1:4,1)),+1,BotSpi(1:4))  ! Abot
        call pol_mless(dcmplx(Mom(1:4,4)),PhotonHel,PhoPol(1:4))      ! photon
        WCurr(1:4) = WPolVec(WMinus,Wm_DKmode,Mom(1:4,2:3),WDK_LO)* WProp * g_weak
!         WCurr(1:4) = WCurr_new(1:4)

!DEC$ IF(_CheckWPolVec .EQ.1)
        call ubarSpi_Weyl(dcmplx(Mom(1:4,2)),-1,BarSpi(1:4))  ! l- or dn
        call    vSpi_Weyl(dcmplx(Mom(1:4,3)),+1,Spi(1:4))     ! nubar or up_bar
        WCurr(1:4)  = vbqq_Weyl(4,BarSpi(1:4),Spi(1:4)) * WProp * g2_weak ! vbqq introduces -i/Sqrt(2)
!         print *, "check DKP_LO_T2 WPolVec routine:"
!         print *, WCurr(1:4) - WCurr_new(1:4)
!         pause
!DEC$ ENDIF
!         print *, "check of gauge invariance for photon in decay"
!         PhoPol(1:4)=dcmplx(Mom(1:4,4))


!       connect to quark current: diagram 1: photon emission off top quark
        Spi(1:4) = vbqg_Weyl( BotSpi(1:4),WCurr(1:4) )  ! vgbq introduces -i/Sqrt(2)

        PropMom(1:4) = TopQuark%Mom(1:4) - dcmplx(Mom(1:4,4))
!         TopProp = sc_(PropMom(1:4),PropMom(1:4))-m_Top**2
        TopProp = -2d0*(TopQuark%Mom(1:4).dot.Mom(1:4,4))
        Spi(1:4) = (0d0,1d0)*(-spi2_Weyl(PropMom(1:4),Spi(1:4)) + m_Top*Spi(1:4) )/TopProp  ! top propagator

        Spi(1:4) = spi2_Weyl(PhoPol(1:4),Spi(1:4)) * coupl_sqrt*(-Q_top)
        TopQuark%Pol(1:4) =( spi2_Weyl(TopQuark%Mom(1:4),Spi(1:4)) - m_Top*Spi(1:4) ) * NWAFactor_Top

!       adding diagram 2: photon emission off bottom quark
        Spi(1:4) = spi2_Weyl(PhoPol(1:4),BotSpi(1:4)) * coupl_sqrt*(-Q_dn)

        PropMom(1:4) = dcmplx( Mom(1:4,1)+Mom(1:4,4) )
        TopProp = sc_(PropMom(1:4),PropMom(1:4))
        Spi(1:4) = (0d0,1d0)*(-spi2_Weyl(PropMom(1:4),Spi(1:4)) )/TopProp  ! bot propagator

        Spi(1:4) = vbqg_Weyl(Spi(1:4),WCurr(1:4))  ! introduces -i/Sqrt(2)
        TopQuark%Pol(1:4) = TopQuark%Pol(1:4) + ( spi2_Weyl(TopQuark%Mom(1:4),Spi(1:4)) - m_Top*Spi(1:4) ) * NWAFactor_Top


!       adding diagram 3: photon emission off W- boson
        WMom(1:4) = Mom(1:4,2)+Mom(1:4,3)
        PropMom(1:4) = WMom(1:4) + Mom(1:4,4)

        WCurr1(1:4) = Vertex_WPW(Mom(1:4,4),WMom(1:4),PhoPol(1:4),WCurr(1:4)) *(-coupl_sqrt)*( Q_dn-Q_top )! = Q_Wm
        WCurr2(1:4) = -WCurr1(1:4) + (WCurr1(1:4).dot.PropMom(1:4))/m_W**2*PropMom(1:4)
        Spi(1:4) = vbqg_Weyl(BotSpi(1:4),WCurr2(1:4)) * (0d0,+1d0)/(2d0*(WMom(1:4).dot.Mom(1:4,4)))  ! introduces -i/Sqrt(2)
        TopQuark%Pol(1:4) = TopQuark%Pol(1:4) + ( spi2_Weyl(TopQuark%Mom(1:4),Spi(1:4)) - m_Top*Spi(1:4) ) * NWAFactor_Top

! print *, "tb:",TopQuark%Pol(1:4)
! print *, "tb:",spi2_Weyl(TopQuark%Mom(1:4),TopQuark%Pol(1:4))/TopQuark%Pol(1:4)
! pause

!       convert Weyl to Dirac
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)

    endif


!         if( dabs((Mom(1:4,1).dot.Mom(1:4,4))/m_top**2) .lt.1d-2) then
!             TopQuark%Pol(1:4)= (0d0,0d0)
!         endif
!         if( cdabs((TopQuark%Mom(1:4).dot.Mom(1:4,4))/m_top**2) .lt.1d-2) then
!             TopQuark%Pol(1:4)= (0d0,0d0)
!         endif


!-----------------------------------------------------------------------------
ELSEIF( Topol.eq.DKP_LO_L ) THEN! photon emission in W decay
!-----------------------------------------------------------------------------
    coupl_sqrt = (0d0,-1d0)*dsqrt(alpha4Pi)
    if( TopQuark%PartType.eq.Top_ ) then ! Top quark decay
        call ubarSpi_Weyl(dcmplx(Mom(1:4,1)),-1,BotSpi(1:4))  ! bot
        WCurr(1:4) = WPolVec(WPlus,Wp_DKmode+1,Mom(1:4,2:4),WDK_LO,LamP=PhotonHel)* WProp * g_weak
!         WCurr(1:4) = WCurr_new(1:4)
!DEC$ IF(_CheckWPolVec .EQ.1)
        call pol_mless(dcmplx(Mom(1:4,4)),PhotonHel,PhoPol(1:4))      ! photon
        call    vSpi_Weyl(dcmplx(Mom(1:4,2)),+1,Spi(1:4))     ! l+ or dn_bar
        call ubarSpi_Weyl(dcmplx(Mom(1:4,3)),-1,BarSpi(1:4))  ! nu or up
!         print *, "check of gauge invariance for photon in decay"
!         PhoPol(1:4)=dcmplx(Mom(1:4,4))

!       diagram 1: photon emission off lepton
        LepSpi(1:4) = spi2_Weyl(PhoPol(1:4),Spi(1:4)) * coupl_sqrt*(-Q_el)
        PropMom(1:4) = Mom(1:4,4) + Mom(1:4,2)
        LepSpi(1:4) = -spi2_Weyl(PropMom(1:4),LepSpi(1:4))*(0d0,1d0)/(2d0*(Mom(1:4,4).dot.Mom(1:4,2)))
        WCurr1(1:4)  = vbqq_Weyl(4,BarSpi(1:4),LepSpi(1:4)) * WProp * g2_weak ! vbqq introduces -i/Sqrt(2)

!       diagram 2: photon emission off W+ boson
        WCurr2(1:4)  = vbqq_Weyl(4,BarSpi(1:4),Spi(1:4)) * WProp * g2_weak ! vbqq introduces -i/Sqrt(2)
        WMom(1:4) = Mom(1:4,2)+Mom(1:4,3)
        WCurr2(1:4) = Vertex_WPW(Mom(1:4,4),WMom(1:4),PhoPol(1:4),WCurr2(1:4)) * (0d0,-1d0)/((WMom(1:4).dot.WMom(1:4))-m_W**2) *(-coupl_sqrt*Q_Wp)

!       connect to quark current
        WCurr(1:4) = WCurr1(1:4) + WCurr2(1:4)
        WMom(1:4) = Mom(1:4,2)+Mom(1:4,3)+Mom(1:4,4)
        WCurr(1:4) = -WCurr(1:4) + (WCurr(1:4).dot.WMom(1:4))/m_W**2 * WMom(1:4)
!DEC$ ENDIF

        BarSpi(1:4) = vgq_Weyl( WCurr(1:4),BotSpi(1:4) ) ! introduces -i/Sqrt(2)
        TopQuark%Pol(1:4) = ( spb2_Weyl(BarSpi(1:4),TopQuark%Mom(1:4)) + m_Top*BarSpi(1:4) ) * NWAFactor_Top
!       convert Weyl to Dirac
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)


    elseif( TopQuark%PartType.eq.ATop_ ) then  ! Anti-Top quark decay
!       assemble lepton current
        call    vSpi_Weyl(dcmplx(Mom(1:4,1)),+1,BotSpi(1:4))  ! Abot
        WCurr(1:4) = WPolVec(WMinus,Wm_DKmode+1,Mom(1:4,2:4),WDK_LO,LamP=PhotonHel)* WProp * g_weak
!         WCurr(1:4) = WCurr_new(1:4)
!DEC$ IF(_CheckWPolVec .EQ.1)
        call pol_mless(dcmplx(Mom(1:4,4)),PhotonHel,PhoPol(1:4))      ! photon
        call ubarSpi_Weyl(dcmplx(Mom(1:4,2)),-1,BarSpi(1:4))  ! l- or dn
        call    vSpi_Weyl(dcmplx(Mom(1:4,3)),+1,Spi(1:4))     ! nubar or up_bar
!         print *, "check of gauge invariance for photon in decay"
!         PhoPol(1:4)=dcmplx(Mom(1:4,4))

!       diagram 1: photon emission off lepton
        LepSpi(1:4) = spb2_Weyl(BarSpi(1:4),PhoPol(1:4)) * coupl_sqrt*(Q_el)
        PropMom(1:4) = Mom(1:4,4) + Mom(1:4,2)
        LepSpi(1:4) = spb2_Weyl(LepSpi(1:4),PropMom(1:4))*(0d0,1d0)/(2d0*(Mom(1:4,4).dot.Mom(1:4,2)))
        WCurr1(1:4)  = vbqq_Weyl(4,LepSpi(1:4),Spi(1:4)) * WProp * g2_weak ! vbqq introduces -i/Sqrt(2)


!       diagram 2: photon emission off W- boson!                       MINUS INTRODUCED TO ENSURE GAUGE INVARIANCE!
        WCurr2(1:4)  = vbqq_Weyl(4,BarSpi(1:4),Spi(1:4)) * WProp * g2_weak ! vbqq introduces -i/Sqrt(2)
        WMom(1:4) = Mom(1:4,2)+Mom(1:4,3)
        WCurr2(1:4) = Vertex_WPW(Mom(1:4,4),WMom(1:4),PhoPol(1:4),WCurr2(1:4)) *(0d0,-1d0)/((WMom(1:4).dot.WMom(1:4))-m_W**2) *(-coupl_sqrt*Q_Wm)   *(-1d0)

!       connect to quark current
        WCurr(1:4) = (WCurr1(1:4) + WCurr2(1:4))
        WMom(1:4) = Mom(1:4,2)+Mom(1:4,3)+Mom(1:4,4)
        WCurr(1:4) = -WCurr(1:4) + (WCurr(1:4).dot.WMom(1:4))/m_W**2 * WMom(1:4)
!DEC$ ENDIF

        LepSpi(1:4) = vbqg_Weyl(BotSpi(1:4),WCurr(1:4))  ! introduces -i/Sqrt(2)
        TopQuark%Pol(1:4) = ( spi2_Weyl(TopQuark%Mom(1:4),LepSpi(1:4)) - m_Top*LepSpi(1:4) ) * NWAFactor_Top

!       convert Weyl to Dirac
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)
    endif





!-----------------------------------------------------------------------------
ELSEIF( Topol.eq.DKP_1L_T ) THEN! virtual corrections to photon emission from top quark line
!-----------------------------------------------------------------------------
 If(TopQuark%PartType .eq. Top_) then !  t -> b W+ gamma
! gau =0 default
! gau =1 photon gauge invariance
    gau =0
! xe = 0,-1,-2 <--> eps^0, eps^-1, eps^-2, default is 0 of course
    xe =0
    Top_Atop = 1
! LO_NLO =1,2 <--> LO,NLO
    LO_NLO = 2
    hel_bot = -1
    heli_virtuals(1) = hel_bot
    heli_virtuals(2) = PhotonHel
! Set W-decay parameters
    Wpar(1) = WPlus
    Wpar(2) = Wp_DKmode
    Wpar(3) = 0
    Wpar(4) = 0
    call EvalTTBP_1LDK(LO_NLO,Top_Atop,xe,gau,heli_virtuals,Wpar,Mom,SpiOut)
    TopQuark%Pol(1:4) = SpiOut(1:4)*NWAFactor_Top*NWAFactor_W

 elseif(TopQuark%PartType .eq. ATop_) then  !  tbar -> bbar W- gamma
! gau =0 default
! gau =1 photon gauge invariance
    gau =0
! xe = 0,-1,-2 <--> eps^0, eps^-1, eps^-2, default is 0 of course
    xe =0
    Top_Atop = -1
! LO_NLO =1,2 <--> LO,NLO
    LO_NLO = 2
    hel_bot = 1
    heli_virtuals(1) = hel_bot
    heli_virtuals(2) = PhotonHel
    Wpar(1) = WMinus
    Wpar(2) = Wm_DKmode
    Wpar(3) = 0
    Wpar(4) = 0

    call EvalTTBP_1LDK(LO_NLO,Top_Atop,xe,gau,heli_virtuals,Wpar,Mom,SpiOut)
    TopQuark%Pol(1:4) = SpiOut(1:4)*NWAFactor_Top*NWAFactor_W
 endif




!-----------------------------------------------------------------------------
ELSEIF( Topol.eq.DKP_RE_T ) THEN! virtual corrections to photon emission from top quark line
!-----------------------------------------------------------------------------
 If(TopQuark%PartType .eq. Top_) then !  t -> b W+ gamma + X
    Top_Atop = 1
    hel_bot = -1
    heli_reals(1) = hel_bot
    heli_reals(2) = GluonHel
    heli_reals(3) = PhotonHel
! gau =0 default
!
! gau =1 photon gauge invariance
!
! gau =2 gluon gauge invariance
    gau =0
! W decay mode
    Wpar(1) = WPlus
    Wpar(2) = Wp_DKmode
    Wpar(3) = 0
    Wpar(4) = 0

    call EvalTTBP_REALDK(Top_Atop,heli_reals,Mom,gau,Wpar,SpiOut)
    TopQuark%Pol(1:4) = SpiOut(1:4)*NWAFactor_Top*NWAFactor_W


 elseif(TopQuark%PartType .eq. ATop_) then  !  tbar -> bbar W- gamma + X
    Top_Atop = -1
    hel_bot = 1
! gau =0 default
!
! gau =1 photon gauge invariance
!
! gau =2 gluon gauge invariance
    gau =0
    heli_reals(1) = hel_bot
    heli_reals(2) = GluonHel
    heli_reals(3) = PhotonHel
! W decay mode
    Wpar(1) = WMinus
    Wpar(2) = Wm_DKmode
    Wpar(3) = 0
    Wpar(4) = 0

    call EvalTTBP_REALDK(Top_Atop,heli_reals,Mom,gau,Wpar,SpiOut)
    TopQuark%Pol(1:4) = SpiOut(1:4)*NWAFactor_Top*NWAFactor_W
 endif



!-----------------------------------------------------------------------------
ELSEIF( Topol.eq.DKP_1L_L ) THEN! 1-loop correction to top quark line + int.dipoles + photon emission form W
!-----------------------------------------------------------------------------
!(formulae taken from Single Top Production in MCFM)
!epinv = (4*pi*musq/mtsq)^ep/Gamma(1-ep)/ep
    rsq = 1d0-2d0*(TopQuark%Mom(1:4).dot.Mom(1:4,1))/m_Top**2  ! = m_W^2 / m_Top^2 for on-shell W boson
    omrsq=1d0-rsq
    wlog=dlog(omrsq)
    rlog=dlog(rsq)
    Kfun=1d0/rsq*wlog
    epinv2 = 0d0
    epinv  = 0d0
    if( rsq.gt.1d0 ) call Error("r2>1 in virtual correction to the top decay")
!   virtual corrections
    F0 = - epinv2 -epinv*(2.5d0-2d0*wlog)  &
         - 5.5d0 - pisqo6 - 2d0*DLi2(rsq)  &  ! -5.5 corresponds to eta=0 (FHD scheme)
         + 3d0*wlog - 2d0*wlog**2 - Kfun
    F1 = 2d0*Kfun
!   integrated dipoles
    if( (TopDecays.eq.6 .and. TopQuark%PartType.eq.Top_) .or. (TopDecays.eq.5 .and. TopQuark%PartType.eq.ATop_) ) then
        xx = xIntDip(1)
        z  = xIntDip(2)
        FF_xxz = FF(xx/z)
        FF_xx  = FF(xx)
        Cdip = (epinv2-2d0*wlog*epinv+2d0*wlog**2 + 2d0-DblPi**2/6.0)*FF_xx  &
             - (epinv-2d0*wlog)*( 2d0/(1d0-z)*(FF_xxz/z - FF_xx)       &
             - (1d0+z)/z*FF_xxz -  FF_xx)                              &
             +  eta*(1d0-z)*FF_xxz/z                                    &
             + 4d0*dlog(1d0-z)/(1d0-z)*(FF_xxz/z - FF_xx)              &
             + FF_xxz/z*(-2d0*(1d0+z)*dlog(1d0-z)                       &
             - (2d0 /(1d0-z)-(1d0+z))*dlog((rsq+z*(1d0-rsq))/z**2) )      &
             - 2d0/(1d0-z)*(z/(rsq+z*(1d0-rsq))*FF_xxz/z - FF_xx)

        C_ct = (epinv+dlog(m_top**2/MuFrag**2))*( FF_xxz/z*(2d0/(1d0-z) - 1d0-z) &
             - FF_xx*(2d0/(1d0-z)- 3d0/2d0) )
        IntDip = 1d0*( Cdip + C_ct )
        F0 = F0 * FF_xx
        F1 = F1 * FF_xx
    else
        IntDip = epinv2 + epinv*(2.5d0-2d0*wlog)  &
           + 25d0/4d0 + 0.5d0*(1d0/omrsq**2-8d0/omrsq+7d0)*rlog  &
           + 0.5d0/omrsq + 2d0*DLi2(omrsq) - 5d0*pisqo6  &
           - 5d0*wlog + 2d0*wlog**2 &
           - ( 2d0*dlog(alpha_DKTfi)**2 -dlog(alpha_DKTfi)*(-7d0/2d0+4d0*alpha_DKTfi-alpha_DKTfi**2/2d0) &
           - 2d0*(1d0-alpha_DKTfi)*rsq/(omrsq)*dlog(rsq/(1d0-alpha_DKTfi+rsq*alpha_DKTfi)) )
    endif


    if( TopQuark%PartType.eq.Top_ ) then ! Top quark decay
        call ubarSpi_Weyl(dcmplx(Mom(1:4,1)),-1,BotSpi(1:4))  ! bot
        WCurr(1:4) = WPolVec(WPlus,Wp_DKmode+1,Mom(1:4,2:4),WDK_LO,LamP=PhotonHel)* WProp * g_weak
!        WCurr(1:4) = WPolVec(WPlus,Wp_DKmode,Mom(1:4,2:3),WDK_LO)* WProp * g_weak
!         WCurr(1:4) = WCurr_new(1:4)
!DEC$ IF(_CheckWPolVec .EQ.1)
        call    vSpi_Weyl(dcmplx(Mom(1:4,2)),+1,Spi(1:4))     ! l+ or dn_bar
        call ubarSpi_Weyl(dcmplx(Mom(1:4,3)),-1,BarSpi(1:4))  ! nu or up
        WCurr(1:4)  = vbqq_Weyl(4,BarSpi(1:4),Spi(1:4)) * WProp * g2_weak ! vbqq introduces -i/Sqrt(2)
!         print *, "check DK_1L_T1 WPolVec routine:"
!         print *, WCurr(1:4) - WCurr_new(1:4)
!         pause
!DEC$ ENDIF
!       connect to quark current
        TopQuark%Pol(1:4) = (F0+IntDip)*vgq_Weyl( WCurr(1:4),BotSpi(1:4) ) ! vgq introduces -i/Sqrt(2)
        TopQuark%Pol(1:4) = TopQuark%Pol(1:4) + BotSpi(1:4)*F1*(0d0,-1d0)/dsqrt(2d0)*sc_(dcmplx(Mom(1:4,1)),WCurr(1:4))/m_Top ! introduce -i/Sqrt(2) by hand
        TopQuark%Pol(1:4) =( spb2_Weyl(TopQuark%Pol(1:4),TopQuark%Mom(1:4)) + m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top * alpha_sOver2Pi*RunAlphaS(NLOParam,MuRen) * CF
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)


    elseif( TopQuark%PartType.eq.ATop_ ) then  ! Anti-Top quark decay
        call    vSpi_Weyl(dcmplx(Mom(1:4,1)),+1,BotSpi(1:4))  ! Abot
        WCurr(1:4) = WPolVec(WMinus,Wm_DKmode+1,Mom(1:4,2:4),WDK_LO,LamP=PhotonHel)* WProp * g_weak
!        WCurr(1:4) = WPolVec(WMinus,Wm_DKmode,Mom(1:4,2:3),WDK_LO)* WProp * g_weak
!         WCurr(1:4) = WCurr_new(1:4)
!DEC$ IF(_CheckWPolVec .EQ.1)
        call ubarSpi_Weyl(dcmplx(Mom(1:4,2)),-1,BarSpi(1:4))  ! l- or dn
        call    vSpi_Weyl(dcmplx(Mom(1:4,3)),+1,Spi(1:4))     ! nubar or up_bar
        WCurr(1:4)  = vbqq_Weyl(4,BarSpi(1:4),Spi(1:4)) * WProp * g2_weak ! vbqq introduces -i/Sqrt(2)
!         print *, "check DK_1L_T2 WPolVec routine:"
!         print *, WCurr(1:4) - WCurr_new(1:4)
!         pause
!DEC$ ENDIF

!       connect to quark current:
        TopQuark%Pol(1:4) = (F0+IntDip)*vbqg_Weyl( BotSpi(1:4),WCurr(1:4) )! vbqg introduces -i/Sqrt(2)
!         call vSpi_Weyl(dcmplx(Mom(1:4,1)),-1,BotSpi(1:4))  ! Abot
        ! INTORODUCED A MINUS SIGN FOR F1 CONTRIBUTION
        TopQuark%Pol(1:4) = TopQuark%Pol(1:4) - BotSpi(1:4)*F1*(0d0,-1d0)/dsqrt(2d0)*sc_(dcmplx(Mom(1:4,1)),WCurr(1:4))/m_Top ! introduce -i/Sqrt(2) by hand
        TopQuark%Pol(1:4) = ( spi2_Weyl(TopQuark%Mom(1:4),TopQuark%Pol(1:4)) - m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top * alpha_sOver2Pi*RunAlphaS(NLOParam,MuRen) * CF
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)
    endif



!-----------------------------------------------------------------------------
ELSEIF( Topol.eq.DKP_RE_L ) THEN! real emission from top quark line
!-----------------------------------------------------------------------------
    coupl_sqrt = dsqrt(alpha_s4Pi*RunAlphaS(NLOParam,MuRen)*CF)*dsqrt(2d0)  ! dsqrt(2d0) factor is to compensate for 1/dsqrt(2d0) factors in gluon coupling
    if( TopQuark%PartType.eq.Top_ ) then ! Top quark decay
        call ubarSpi_Weyl(dcmplx(Mom(1:4,1)),-1,BotSpi(1:4))  ! bot
        call pol_mless(dcmplx(Mom(1:4,5)),GluonHel,GluPol(1:4))        ! glu
        WCurr(1:4) = WPolVec(WPlus,Wp_DKmode+1,Mom(1:4,2:4),WDK_LO,LamP=PhotonHel)* WProp * g_weak
!        WCurr(1:4) = WPolVec(WPlus,Wp_DKmode,Mom(1:4,2:3),WDK_LO)* WProp * g_weak
!         WCurr(1:4) = WCurr_new(1:4)
!DEC$ IF(_CheckWPolVec .EQ.1)
        call    vSpi_Weyl(dcmplx(Mom(1:4,2)),+1,Spi(1:4))     ! l+ or dn_bar
        call ubarSpi_Weyl(dcmplx(Mom(1:4,3)),-1,BarSpi(1:4))  ! nu or up
        WCurr(1:4)  = vbqq_Weyl(4,BarSpi(1:4),Spi(1:4)) * WProp * g2_weak ! vbqq introduces -i/Sqrt(2)
!         print *, "check DK_RE_T1 WPolVec routine:"
!         print *, WCurr(1:4) - WCurr_new(1:4)
!         pause
!DEC$ ENDIF

!       connect to quark current: diagram 1: gluon emission from top propagator
        BarSpi(1:4) = vgq_Weyl( WCurr(1:4),BotSpi(1:4) ) ! vgq introduces -i/Sqrt(2)

        PropMom(1:4) = TopQuark%Mom(1:4) - dcmplx(Mom(1:4,5))
        TopProp = sc_(PropMom(1:4),PropMom(1:4))-m_Top**2
        BarSpi(1:4) = (0d0,1d0)*( spb2_Weyl(BarSpi(1:4),PropMom(1:4)) + m_Top*BarSpi(1:4) )/TopProp  ! top propagator

        BarSpi(1:4) = vgq_Weyl( GluPol(1:4),BarSpi(1:4) ) * coupl_sqrt ! vgq introduces Sqrt(2) which is removed by coupl_sqrt
        TopQuark%Pol(1:4) =( spb2_Weyl(BarSpi(1:4),TopQuark%Mom(1:4)) + m_Top*BarSpi(1:4) ) * NWAFactor_Top


!       adding diagram 2: gluon emission from bottom propagator
        BarSpi(1:4) = vgq_Weyl( GluPol(1:4),BotSpi(1:4) ) * coupl_sqrt ! vgq introduces Sqrt(2) which is removed by coupl_sqrt

        PropMom(1:4) = dcmplx( Mom(1:4,1)+Mom(1:4,5) )
        TopProp = sc_(PropMom(1:4),PropMom(1:4))
        BarSpi(1:4) = (0d0,1d0)*( spb2_Weyl(BarSpi(1:4),PropMom(1:4)) )/TopProp  ! top propagator

        BarSpi(1:4) = vgq_Weyl( WCurr(1:4),BarSpi(1:4) )  ! vgq introduces -i/Sqrt(2)
        TopQuark%Pol(1:4) = TopQuark%Pol(1:4) + ( spb2_Weyl(BarSpi(1:4),TopQuark%Mom(1:4)) + m_Top*BarSpi(1:4) ) * NWAFactor_Top

        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)


    elseif( TopQuark%PartType.eq.ATop_ ) then  ! Anti-Top quark decay
!       assemble lepton current
        call    vSpi_Weyl(dcmplx(Mom(1:4,1)),+1,BotSpi(1:4))  ! Abot
        call pol_mless(dcmplx(Mom(1:4,5)),GluonHel,GluPol(1:4))      ! glu
        WCurr(1:4) = WPolVec(WMinus,Wm_DKmode+1,Mom(1:4,2:4),WDK_LO,LamP=PhotonHel)* WProp * g_weak
!        WCurr(1:4) = WPolVec(WMinus,Wm_DKmode,Mom(1:4,2:3),WDK_LO)* WProp * g_weak
!         WCurr(1:4) = WCurr_new(1:4)
!DEC$ IF(_CheckWPolVec.EQ.1)
        call ubarSpi_Weyl(dcmplx(Mom(1:4,2)),-1,BarSpi(1:4))  ! l- or dn
        call    vSpi_Weyl(dcmplx(Mom(1:4,3)),+1,Spi(1:4))     ! nubar or up_bar
        WCurr(1:4)  = vbqq_Weyl(4,BarSpi(1:4),Spi(1:4)) * WProp * g2_weak ! vbqq introduces -i/Sqrt(2)
!         print *, "check DK_RE_T2 WPolVec routine:"
!         print *, WCurr(1:4) - WCurr_new(1:4)
!         pause
!DEC$ ENDIF


!       connect to quark current: diagram 1: gluon emission from top propagator
        Spi(1:4) = vbqg_Weyl( BotSpi(1:4),WCurr(1:4) )  ! vgbq introduces -i/Sqrt(2)

        PropMom(1:4) = TopQuark%Mom(1:4) - dcmplx(Mom(1:4,5))
        TopProp = sc_(PropMom(1:4),PropMom(1:4))-m_Top**2
        Spi(1:4) = (0d0,1d0)*(-spi2_Weyl(PropMom(1:4),Spi(1:4)) + m_Top*Spi(1:4) )/TopProp  ! top propagator

        Spi(1:4) = vbqg_Weyl( Spi(1:4),GluPol(1:4) ) * coupl_sqrt ! vgbq introduces Sqrt(2) which is removed by coupl_sqrt
        TopQuark%Pol(1:4) =( spi2_Weyl(TopQuark%Mom(1:4),Spi(1:4)) - m_Top*Spi(1:4) ) * NWAFactor_Top

!       adding diagram 2: diagram 1: gluon emission from bottom propagator
        Spi(1:4) = vbqg_Weyl(BotSpi(1:4),GluPol(1:4)) * coupl_sqrt ! vgbq introduces Sqrt(2) which is removed by coupl_sqrt

        PropMom(1:4) = dcmplx( Mom(1:4,1)+Mom(1:4,5) )
        TopProp = sc_(PropMom(1:4),PropMom(1:4))
        Spi(1:4) = (0d0,1d0)*(-spi2_Weyl(PropMom(1:4),Spi(1:4)) )/TopProp  ! top propagator

        Spi(1:4) = vbqg_Weyl(Spi(1:4),WCurr(1:4))  ! vgq introduces -i/Sqrt(2)
        TopQuark%Pol(1:4) = TopQuark%Pol(1:4) + ( spi2_Weyl(TopQuark%Mom(1:4),Spi(1:4)) - m_Top*Spi(1:4) ) * NWAFactor_Top

        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)
    endif





!-----------------------------------------------------------------------------
ELSEIF( Topol.eq.DKJ_1L_T ) THEN! real emission from top quark line
!-----------------------------------------------------------------------------


    xe = 0
    MomT(1:4) = dble(TopQuark%Mom(1:4))
!wrong!    IntDip = 0.5d0*Tb_x_Tg * DKff_qq(xe,m_Top**2,2d0*(Mom(1:4,1).dot.Mom(1:4,4)))  &  !  cancels ff_(bb-g1,g2)        1/eps*Log checked for mu=mt       0.5=symm.fact.
!wrong!           + 0.5d0*Tb_x_Tg * DKff_qq(xe,m_Top**2,2d0*(Mom(1:4,1).dot.Mom(1:4,4)))  &  !  cancels ff_(bb-g2,g1)        1/eps*Log checked for mu=mt       0.5=symm.fact.
!wrong!           + 0.5d0*Tb_x_Tg * DKff_gg(xe,m_Top**2,2d0*(Mom(1:4,1).dot.Mom(1:4,4)))  &  !  cancels ff_(g1-g2,bb)        1/eps*Log checked for mu=mt       0.5=symm.fact.
!wrong!           + 0.5d0*Tt_x_Tb * DKfi_qq(xe,m_Top**2,-2d0*(Mom(1:4,1).dot.MomT(1:4)))  &  !  cancels fi_(bb-g1,tb)        1/eps*Log checked for mu=mt       0.5=symm.fact.
!wrong!           + 0.5d0*Tt_x_Tb * DKfi_qq(xe,m_Top**2,-2d0*(Mom(1:4,1).dot.MomT(1:4)))  &  !  cancels fi_(bb-g2,tb)        1/eps*Log checked for mu=mt       0.5=symm.fact.
!wrong!           + 0.5d0*Tt_x_Tg * DKfi_gg(xe,m_Top**2,-2d0*(Mom(1:4,4).dot.MomT(1:4))) *2d0 &  !  cancels fi_(g1-g2,tb)             1/eps*Log checked for mu=mt       0.5=symm.fact.
!wrong!           + 0.5d0*Tt_x_Tg * DKfi_gg(xe,m_Top**2,-2d0*(Mom(1:4,4).dot.MomT(1:4))) *2d0 &  !  cancels fi_(g2-g1,tb)             1/eps*Log checked for mu=mt       0.5=symm.fact.
!wrong!           - TR * Nf_light * DKff_qg(xe,m_Top**2,2d0*(Mom(1:4,1).dot.Mom(1:4,4))) !  cancels ff_(q1-q2,bb)    Nf_light=(uub+ddb+ccb+ssb)*1dipole + bbb*2dipoles*symm.fact.

    IntDip = 0.5d0*Tb_x_Tg * DKff_qq(xe,muRen**2,2d0*(Mom(1:4,1).dot.Mom(1:4,4)))  &  !  cancels ff_(bb-g1,g2)        1/eps*Log checked for mu=mt       0.5=symm.fact.
           + 0.5d0*Tb_x_Tg * DKff_qq(xe,muRen**2,2d0*(Mom(1:4,1).dot.Mom(1:4,4)))  &  !  cancels ff_(bb-g2,g1)        1/eps*Log checked for mu=mt       0.5=symm.fact.
           + 0.5d0*Tb_x_Tg * DKff_gg(xe,muRen**2,2d0*(Mom(1:4,1).dot.Mom(1:4,4)))  &  !  cancels ff_(g1-g2,bb)        1/eps*Log checked for mu=mt       0.5=symm.fact.
           + 0.5d0*Tt_x_Tb * DKfi_qq(xe,m_Top**2,muRen**2,-2d0*(Mom(1:4,1).dot.MomT(1:4)))  &  !  cancels fi_(bb-g1,tb)        1/eps*Log checked for mu=mt       0.5=symm.fact.
           + 0.5d0*Tt_x_Tb * DKfi_qq(xe,m_Top**2,muRen**2,-2d0*(Mom(1:4,1).dot.MomT(1:4)))  &  !  cancels fi_(bb-g2,tb)        1/eps*Log checked for mu=mt       0.5=symm.fact.
           + 0.5d0*Tt_x_Tg * DKfi_gg(xe,m_Top**2,muRen**2,-2d0*(Mom(1:4,4).dot.MomT(1:4))) *2d0 &  !  cancels fi_(g1-g2,tb)    1/eps*Log checked for mu=mt       0.5=symm.fact.
           + 0.5d0*Tt_x_Tg * DKfi_gg(xe,m_Top**2,muRen**2,-2d0*(Mom(1:4,4).dot.MomT(1:4))) *2d0 &  !  cancels fi_(g2-g1,tb)    1/eps*Log checked for mu=mt       0.5=symm.fact.
           - TR * Nf_light * DKff_qg(xe,muRen**2,2d0*(Mom(1:4,1).dot.Mom(1:4,4))) !  cancels ff_(q1-q2,bb)    Nf_light=(uub+ddb+ccb+ssb)*1dipole + bbb*2dipoles*symm.fact.


    IntDip =-IntDip * alpha_s*RunAlphaS(NLOParam,MuRen)/(2d0*DblPi)   ! what are these fators???
    if( xe.eq.-1 ) IntDip = IntDip + 3d0*DblPi*(0d0,1d0) * alpha_s*RunAlphaS(NLOParam,MuRen)/(2d0*DblPi)
! note: Andreas' virtual corrections come with a 3*I*Pi/eps*tree piece which I cannot reconstruct from Catani
!       therefore I add it to my integrated dipoles. anyways it doesnt contribute to the unpol.mat.el. because it's only a imag.part
! print *, "INT DIP",IntDip / (alpha_s*RunAlphaS(NLOParam,MuRen)/(2d0*DblPi))



!    this is the CATANI CHECK WITHOUT I*PI pieces
!     dummy(4) = -(CA+CF)!  eps2 poles
!    dummy(1) = CF*( -5d0/2d0-dlog(m_top**2*MuRen**2/(2d0*Mom(1:4,1).dot.MomT(1:4))**2) )
!    dummy(2) = CA*( -11d0/6d0-dlog(MuRen**2/(2d0*Mom(1:4,1).dot.Mom(1:4,4)))+0.5d0*dlog(MuRen**2*m_top**2/(2d0*Mom(1:4,1).dot.MomT(1:4))**2)-0.5d0*dlog(MuRen**2*m_top**2/(2d0*Mom(1:4,4).dot.MomT(1:4))**2) )
!    dummy(3) = Nf_light*2d0/3d0*TR
!    dummy(1:3)=dummy(1:3)*alpha_s*RunAlphaS(NLOParam,MuRen)/(4d0*DblPi)  * 2d0! factor 2 because here we normalize to alpha_s/(2pi)
!    dummy(4) = 1d0*dummy(1) + 1d0*dummy(2) + 1d0*dummy(3)
! print *, "CATANI",dummy(4) / (alpha_s*RunAlphaS(NLOParam,MuRen)/(2d0*DblPi))
! pause


!DEC$ IF(_FactCheck .EQ.1)
      if( dabs(Mom(1:4,1).dot.Mom(1:4,4))/m_Top**2 .lt. 1d-3 ) then
          TopQuark%Pol(1:4) = (1d-60,0d0)
          return
      endif
!DEC$ ENDIF

    if( TopQuark%PartType.eq.Top_ ) then
         WCurr(1:4) = WPolVec(WPlus,Wp_DKmode,Mom(1:4,2:3),WDK_LO)* WProp * g_weak/dsqrt(2d0)
         call  calc_TopDecay(MomT(1:4),Mom(1:4,1),Mom(1:4,2),Mom(1:4,3),Mom(1:4,4),GluonHel,WCurr(1:4),xe,MuRen,1,0,TopQuark%Pol(1:4) )! virt
         WCurr1(1:4) = TopQuark%Pol(1:4) * 2d0 !  Andreas' result is normalized to alpha_s/(4Pi) --> multiply by 2

         call calc_TopDecay(MomT(1:4),Mom(1:4,1),Mom(1:4,2),Mom(1:4,3),Mom(1:4,4),GluonHel,WCurr(1:4),0,MuRen,0,0,TopQuark%Pol(1:4) )!  lo
         WCurr2(1:4) = TopQuark%Pol(1:4)  * IntDip

!          print *, "SHARF",(WCurr1(1:4))
!          print *, "ID  ",(WCurr2(1:4))
!          print *, "sum  ",(WCurr2(1:4))+(WCurr1(1:4))
!          pause
!          WCurr1(1:4)=0d0!   for inv. check, virt. decaivated above, also nf deacivated

         TopQuark%Pol(1:4) = WCurr1(1:4) + WCurr2(1:4)
         TopQuark%Pol(1:4) = TopQuark%Pol(1:4) * (0d0,-1d0)! adjust phase to match DKJ_LO_T



    elseif( TopQuark%PartType.eq.ATop_ ) then  ! Anti-Top quark decay
         WCurr(1:4) = WPolVec(WMinus,Wm_DKmode,Mom(1:4,2:3),WDK_LO)* WProp * g_weak/dsqrt(2d0)
         call calc_ATopDecay(MomT(1:4),Mom(1:4,1),Mom(1:4,2),Mom(1:4,3),Mom(1:4,4),GluonHel,WCurr(1:4),xe,MuRen,1,0,TopQuark%Pol(1:4) )! virt
         WCurr1(1:4) = TopQuark%Pol(1:4) * 2d0 !  Andreas' result is normalized to alpha_s/(4Pi) --> multiply by 2

         call calc_ATopDecay(MomT(1:4),Mom(1:4,1),Mom(1:4,2),Mom(1:4,3),Mom(1:4,4),GluonHel,WCurr(1:4),0,MuRen,0,0,TopQuark%Pol(1:4) )!  lo
         WCurr2(1:4) = TopQuark%Pol(1:4) * IntDip

!          print *, "SHARF",(WCurr1(1:4))
!          print *, "ID  ",(WCurr2(1:4))
!          print *, "sum",(WCurr1(1:4))+(WCurr2(1:4))
!          print *, "ratio",(WCurr1(1:4))/(WCurr2(1:4))
!          pause
!            WCurr1(1:4)=0d0!   for inv. check, virt. decaivated above

         TopQuark%Pol(1:4) = WCurr1(1:4) + WCurr2(1:4)
         TopQuark%Pol(1:4) = TopQuark%Pol(1:4) * (0d0,1d0)! adjust phase to match DKJ_LO_T
    endif


    TopQuark%Pol(1:4) = TopQuark%Pol(1:4) * NWAFactor_Top
    TopQuark%Pol(5:16)= (0d0,0d0)



!-----------------------------------------------------------------------------
ELSEIF( Topol.eq.DKJ_RE_TT ) THEN! two gluon emission from top quark line (gg)
!-----------------------------------------------------------------------------
    if( init(DKJ_RE_TT) ) then
        allocate( TopDKAmp(2)%NumGlu(0:2))
        allocate( TopDKAmp(2)%PartRef(1:5))
        allocate( TopDKAmp(2)%PartType(1:5))
        allocate( TopDKAmp(2)%Quarks(1:2))
        allocate( TopDKAmp(2)%Gluons(1:2))
        allocate( ATopDKAmp(2)%NumGlu(0:2))
        allocate( ATopDKAmp(2)%PartRef(1:5))
        allocate( ATopDKAmp(2)%PartType(1:5))
        allocate( ATopDKAmp(2)%Quarks(1:2))
        allocate( ATopDKAmp(2)%Gluons(1:2))
        allocate( TopDKAmp(3)%NumGlu(0:2))
        allocate( TopDKAmp(3)%PartRef(1:5))
        allocate( TopDKAmp(3)%PartType(1:5))
        allocate( TopDKAmp(3)%Quarks(1:2))
        allocate( TopDKAmp(3)%Gluons(1:2))
        allocate( ATopDKAmp(3)%NumGlu(0:2))
        allocate( ATopDKAmp(3)%PartRef(1:5))
        allocate( ATopDKAmp(3)%PartType(1:5))
        allocate( ATopDKAmp(3)%Quarks(1:2))
        allocate( ATopDKAmp(3)%Gluons(1:2))

        TopDKAmp(2)%NumPart =5; TopDKAmp(2)%NumQua =2; TopDKAmp(2)%NumGlu(0) = 2
        TopDKAmp(3)%NumPart =5; TopDKAmp(3)%NumQua =2; TopDKAmp(3)%NumGlu(0) = 2
        ATopDKAmp(2)%NumPart=5; ATopDKAmp(2)%NumQua=2; ATopDKAmp(2)%NumGlu(0)= 2
        ATopDKAmp(3)%NumPart=5; ATopDKAmp(3)%NumQua=2; ATopDKAmp(3)%NumGlu(0)= 2
        TopDKAmp(2)%PartRef(1:5) = (/1,2,4,5,3/)
        TopDKAmp(3)%PartRef(1:5) = (/1,2,5,4,3/)
        ATopDKAmp(2)%PartRef(1:5)= (/1,2,4,5,3/)
        ATopDKAmp(3)%PartRef(1:5)= (/1,2,5,4,3/)
        call LinkTreeParticles( TopDKAmp(2), TopDKProd(1:5))
        call LinkTreeParticles( TopDKAmp(3), TopDKProd(1:5))
        call LinkTreeParticles(ATopDKAmp(2),ATopDKProd(1:5))
        call LinkTreeParticles(ATopDKAmp(3),ATopDKProd(1:5))
        init(DKJ_RE_TT)=.false.
    endif

!DEC$ IF(_FactCheck .EQ.1)
  if( dabs(Mom(1:4,4).dot.Mom(1:4,5)).lt.1d-3*m_Top**2 ) then
      if( ((Mom(1:4,4)+Mom(1:4,5)).dot.Mom(1:4,1)).lt.1d-3*m_Top**2 ) then
          TopQuark%Pol(1:4) = (1d-60,0d0)
          return
      endif
  endif
  if( dabs(Mom(1:4,1).dot.Mom(1:4,4)).lt.1d-3*m_Top**2 ) then
      if( ((Mom(1:4,1)+Mom(1:4,4)).dot.Mom(1:4,5)).lt.1d-3*m_Top**2 ) then
          TopQuark%Pol(1:4) = (1d-60,0d0)
          return
      endif
  endif
  if( dabs(Mom(1:4,1).dot.Mom(1:4,5)).lt.1d-3*m_Top**2 ) then
      if( ((Mom(1:4,1)+Mom(1:4,5)).dot.Mom(1:4,4)).lt.1d-3*m_Top**2 ) then
          TopQuark%Pol(1:4) = (1d-60,0d0)
          return
      endif
  endif
!DEC$ ENDIF

    coupl_sqrt = alpha_s4Pi*RunAlphaS(NLOParam,MuRen)! color factors & average have to be included outside of this routine
    if( TopQuark%PartType.eq.Top_ ) then ! Top quark decay
        TopDKProd(1)%Mom(1:4) =-TopQuark%Mom(1:4)!  top
        TopDKProd(2)%Mom(1:4) = dcmplx(Mom(1:4,1))! bot
        TopDKProd(3)%Mom(1:4) = dcmplx(Mom(1:4,2)+Mom(1:4,3))!  W+
        TopDKProd(4)%Mom(1:4) = dcmplx(Mom(1:4,4))! glu
        TopDKProd(5)%Mom(1:4) = dcmplx(Mom(1:4,5))! glu
        call ubarSpi_Weyl(TopDKProd(2)%Mom(1:4),-1,TopDKProd(2)%Pol(1:4))
        TopDKProd(3)%Pol(1:4) = WPolVec(WPlus,Wp_DKmode,Mom(1:4,2:3),WDK_LO)  * WProp * g_weak/dsqrt(2d0)
        call pol_mless(TopDKProd(4)%Mom(1:4),GluonHel,TopDKProd(4)%Pol(1:4))
        call pol_mless(TopDKProd(5)%Mom(1:4),Gluon2Hel,TopDKProd(5)%Pol(1:4))
        call new_calc_ampl(4,4,0,0,TopDKAmp(PartAmp+1),TopQuark%Pol(1:4))
        TopQuark%Pol(1:4) =( spb2_Weyl(TopQuark%Pol(1:4),TopQuark%Mom(1:4)) + m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top * coupl_sqrt
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)

    elseif( TopQuark%PartType.eq.ATop_ ) then  ! Anti-Top quark decay
        ATopDKProd(1)%Mom(1:4) =-TopQuark%Mom(1:4)!  top
        ATopDKProd(2)%Mom(1:4) = dcmplx(Mom(1:4,1))! bot
        ATopDKProd(3)%Mom(1:4) = dcmplx(Mom(1:4,2)+Mom(1:4,3))!  W+
        ATopDKProd(4)%Mom(1:4) = dcmplx(Mom(1:4,4))! glu
        ATopDKProd(5)%Mom(1:4) = dcmplx(Mom(1:4,5))! glu
        call vSpi_Weyl(ATopDKProd(2)%Mom(1:4),+1,ATopDKProd(2)%Pol(1:4))
        ATopDKProd(3)%Pol(1:4) = WPolVec(WMinus,Wm_DKmode,Mom(1:4,2:3),WDK_LO)* WProp * g_weak/dsqrt(2d0)
        call pol_mless(ATopDKProd(4)%Mom(1:4),GluonHel,ATopDKProd(4)%Pol(1:4))
        call pol_mless(ATopDKProd(5)%Mom(1:4),Gluon2Hel,ATopDKProd(5)%Pol(1:4))
        call new_calc_ampl(4,4,0,0,ATopDKAmp(PartAmp+1),TopQuark%Pol(1:4))
        TopQuark%Pol(1:4) =( spi2_Weyl(TopQuark%Mom(1:4),TopQuark%Pol(1:4)) - m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top * coupl_sqrt
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)
    endif


!-----------------------------------------------------------------------------
ELSEIF( Topol.eq.DKJ_RE2_TT ) THEN! qqbar emission from top quark line (qqb)!  careful: the nflav factor has to be multiplied outside of this routine
!-----------------------------------------------------------------------------

    coupl_sqrt = alpha_s4Pi*RunAlphaS(NLOParam,MuRen)! color factors & average have to be included outside of this routine
    if( init(DKJ_RE2_TT) ) then
        allocate( TopDKAmp(4)%NumGlu(0:2))
        allocate( TopDKAmp(4)%PartRef(1:4))
        allocate( TopDKAmp(4)%PartType(1:4))
        allocate( TopDKAmp(4)%Quarks(1:2))
        allocate( TopDKAmp(4)%Gluons(1:1))
        allocate( ATopDKAmp(4)%NumGlu(0:2))
        allocate( ATopDKAmp(4)%PartRef(1:4))
        allocate( ATopDKAmp(4)%PartType(1:4))
        allocate( ATopDKAmp(4)%Quarks(1:2))
        allocate( ATopDKAmp(4)%Gluons(1:1))

        TopDKAmp(4)%NumPart=4;   TopDKAmp(4)%NumQua=2; TopDKAmp(4)%NumGlu(0) = 1
        ATopDKAmp(4)%NumPart=4; ATopDKAmp(4)%NumQua=2; ATopDKAmp(4)%NumGlu(0)= 1
        TopDKAmp(4)%PartRef(1:4) = (/1,2,4,3/)
        ATopDKAmp(4)%PartRef(1:4)= (/1,2,4,3/)
        call LinkTreeParticles( TopDKAmp(4), TopDKProd(1:4))
        call LinkTreeParticles(ATopDKAmp(4),ATopDKProd(1:4))
        init(DKJ_RE2_TT)=.false.
    endif

! !DEC$ IF(_FactCheck .EQ.1)
!   if( dabs(Mom(1:4,4).dot.Mom(1:4,5)).lt.1d-3*m_Top**2 ) then
!       if( dabs((Mom(1:4,4)+Mom(1:4,5)).dot.Mom(1:4,1)).lt.1d-3*m_Top**2 ) then
!           TopQuark%Pol(1:4) = (1d-60,0d0)
!           return
!       endif
!   endif
! !DEC$ ENDIF

    if( TopQuark%PartType.eq.Top_ ) then ! Top quark decay   !  has a minus sign wrt. to older version
        TopDKProd(1)%Mom(1:4) =-TopQuark%Mom(1:4)!  top
        TopDKProd(2)%Mom(1:4) = dcmplx(Mom(1:4,1))! bot
        TopDKProd(3)%Mom(1:4) = dcmplx(Mom(1:4,2)+Mom(1:4,3))!  W+
        TopDKProd(4)%Mom(1:4) = dcmplx(Mom(1:4,4)+Mom(1:4,5))! glu
        call ubarSpi_Weyl(TopDKProd(2)%Mom(1:4),-1,TopDKProd(2)%Pol(1:4))! bot
        TopDKProd(3)%Pol(1:4) = WPolVec(WPlus,Wp_DKmode,Mom(1:4,2:3),WDK_LO)* WProp * g_weak/dsqrt(2d0)

        call vSpi_Weyl(dcmplx(Mom(1:4,4)),+GluonHel,Spi(1:4))
        call ubarSpi_Weyl(dcmplx(Mom(1:4,5)),-GluonHel,BarSpi(1:4))
        TopDKProd(4)%Pol(1:4) = vbqq_Weyl(4,BarSpi(1:4),Spi(1:4)) * (0d0,-1d0)/( TopDKProd(4)%Mom(1:4).dot.TopDKProd(4)%Mom(1:4) )
        call new_calc_ampl(4,4,0,0,TopDKAmp(4),TopQuark%Pol(1:4))
        TopQuark%Pol(1:4) =( spb2_Weyl(TopQuark%Pol(1:4),TopQuark%Mom(1:4)) + m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top * coupl_sqrt
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)
    elseif( TopQuark%PartType.eq.ATop_ ) then  ! Anti-Top quark decay
        ATopDKProd(1)%Mom(1:4) =-TopQuark%Mom(1:4)!  top
        ATopDKProd(2)%Mom(1:4) = dcmplx(Mom(1:4,1))! bot
        ATopDKProd(3)%Mom(1:4) = dcmplx(Mom(1:4,2)+Mom(1:4,3))!  W+
        ATopDKProd(4)%Mom(1:4) = dcmplx(Mom(1:4,4)+Mom(1:4,5))! glu
        call vSpi_Weyl(ATopDKProd(2)%Mom(1:4),+1,ATopDKProd(2)%Pol(1:4))
        ATopDKProd(3)%Pol(1:4) = WPolVec(WMinus,Wm_DKmode,Mom(1:4,2:3),WDK_LO)* WProp * g_weak/dsqrt(2d0)
        call vSpi_Weyl(dcmplx(Mom(1:4,4)),+GluonHel,Spi(1:4))
        call ubarSpi_Weyl(dcmplx(Mom(1:4,5)),-GluonHel,BarSpi(1:4))
        ATopDKProd(4)%Pol(1:4) = vbqq_Weyl(4,BarSpi(1:4),Spi(1:4)) * (0d0,-1d0)/( ATopDKProd(4)%Mom(1:4).dot.ATopDKProd(4)%Mom(1:4) )
        call new_calc_ampl(4,4,0,0,ATopDKAmp(4),TopQuark%Pol(1:4))

        TopQuark%Pol(1:4) =( spi2_Weyl(TopQuark%Mom(1:4),TopQuark%Pol(1:4)) - m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top * coupl_sqrt
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)
    endif




!-----------------------------------------------------------------------------
ELSEIF( Topol.eq.DKJ_1LQ_T ) THEN! gluon emission from top quark line and virtual correction on W
!-----------------------------------------------------------------------------
    if( init(DK_RE_T) ) then
        allocate( TopDKAmp(1)%NumGlu(0:2))
        allocate( TopDKAmp(1)%PartRef(1:4))
        allocate( TopDKAmp(1)%PartType(1:4))
        allocate( TopDKAmp(1)%Quarks(1:2))
        allocate( TopDKAmp(1)%Gluons(1:1))
        allocate( ATopDKAmp(1)%NumGlu(0:2))
        allocate( ATopDKAmp(1)%PartRef(1:4))
        allocate( ATopDKAmp(1)%PartType(1:4))
        allocate( ATopDKAmp(1)%Quarks(1:2))
        allocate( ATopDKAmp(1)%Gluons(1:1))

        TopDKAmp(1)%NumPart=4;   TopDKAmp(1)%NumQua=2; TopDKAmp(1)%NumGlu(0) = 1
        ATopDKAmp(1)%NumPart=4; ATopDKAmp(1)%NumQua=2; ATopDKAmp(1)%NumGlu(0)= 1
        TopDKAmp(1)%PartRef(1:4) = (/1,2,4,3/)
        ATopDKAmp(1)%PartRef(1:4)= (/1,2,4,3/)
        call LinkTreeParticles( TopDKAmp(1), TopDKProd(1:4))
        call LinkTreeParticles(ATopDKAmp(1),ATopDKProd(1:4))
        init(DK_RE_T)=.false.
    endif

!DEC$ IF(_FactCheck .EQ.1)
      if( dabs(Mom(1:4,1).dot.Mom(1:4,4))/m_Top**2 .lt. 1d-3 ) then
          TopQuark%Pol(1:4) = (1d-60,0d0)
          return
      endif
!DEC$ ENDIF

    coupl_sqrt = dsqrt(alpha_s4Pi*RunAlphaS(NLOParam,MuRen)*CF)*dsqrt(2d0)  ! dsqrt(2d0) factor is to compensate for 1/dsqrt(2d0) factors in gluon coupling
! is sqrt. correct because interfered with dsqrt(cf)?????????????

!--- definitions that are required
    ep2=0d0; ep1=0d0
    v_factor = alpha_sOver2Pi*RunAlphaS(NLOParam,MuRen)*CF*(-2d0*ep2-3d0*ep1-8d0+pisq)
!---- first factor ``2'' is the number of dipoles
    int_dip_factor = 2d0*alpha_sOver2Pi*RunAlphaS(NLOParam,MuRen)*CF*(ep2+3d0/2d0*ep1+5d0-pisq/2d0+3d0/2d0*(alpha_DKWff-1d0-dlog(alpha_DKWff))-dlog(alpha_DKWff)**2)

    if( TopQuark%PartType.eq.Top_ ) then ! Top quark decay
        TopDKProd(1)%Mom(1:4) =-TopQuark%Mom(1:4)!  top
        TopDKProd(2)%Mom(1:4) = dcmplx(Mom(1:4,1))! bot
        TopDKProd(3)%Mom(1:4) = dcmplx(Mom(1:4,2)+Mom(1:4,3))!  W+
        TopDKProd(4)%Mom(1:4) = dcmplx(Mom(1:4,4))! glu
        call ubarSpi_Weyl(TopDKProd(2)%Mom(1:4),-1,TopDKProd(2)%Pol(1:4))
        TopDKProd(3)%Pol(1:4) = WPolVec(WPlus,Wp_DKmode,Mom(1:4,2:3),WDK_LO)* WProp * g_weak/dsqrt(2d0)
        TopDKProd(3)%Pol(1:4) = TopDKProd(3)%Pol(1:4) * (v_factor + int_dip_factor)
        call pol_mless(TopDKProd(4)%Mom(1:4),GluonHel,TopDKProd(4)%Pol(1:4))
        call new_calc_ampl(4,4,0,0,TopDKAmp(1),TopQuark%Pol(1:4))

        TopQuark%Pol(1:4) =( spb2_Weyl(TopQuark%Pol(1:4),TopQuark%Mom(1:4)) + m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top * coupl_sqrt
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)

    elseif( TopQuark%PartType.eq.ATop_ ) then  ! Anti-Top quark decay
        ATopDKProd(1)%Mom(1:4) =-TopQuark%Mom(1:4)!  top
        ATopDKProd(2)%Mom(1:4) = dcmplx(Mom(1:4,1))! bot
        ATopDKProd(3)%Mom(1:4) = dcmplx(Mom(1:4,2)+Mom(1:4,3))!  W+
        ATopDKProd(4)%Mom(1:4) = dcmplx(Mom(1:4,4))! glu
        call vSpi_Weyl(ATopDKProd(2)%Mom(1:4),+1,ATopDKProd(2)%Pol(1:4))
        ATopDKProd(3)%Pol(1:4) = WPolVec(WMinus,Wm_DKmode,Mom(1:4,2:3),WDK_LO)* WProp * g_weak/dsqrt(2d0)
        ATopDKProd(3)%Pol(1:4) = ATopDKProd(3)%Pol(1:4) * (v_factor + int_dip_factor)

        call pol_mless(ATopDKProd(4)%Mom(1:4),GluonHel,ATopDKProd(4)%Pol(1:4))
        call new_calc_ampl(4,4,0,0,ATopDKAmp(1),TopQuark%Pol(1:4))

        TopQuark%Pol(1:4) =( spi2_Weyl(TopQuark%Mom(1:4),TopQuark%Pol(1:4)) - m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top * coupl_sqrt
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)
    endif




!-----------------------------------------------------------------------------
ELSEIF( Topol.eq.DKJ_1L_Q ) THEN! 1-loop correction to top quark line + int.dipoles and gluon emission from W
!-----------------------------------------------------------------------------
!(formulae taken from Single Top Production in MCFM)
!epinv = (4*pi*musq/mtsq)^ep/Gamma(1-ep)/ep
    rsq = 1d0-2d0*(TopQuark%Mom(1:4).dot.Mom(1:4,1))/m_Top**2  ! = m_W^2 / m_Top^2 for on-shell W boson
    omrsq=1d0-rsq
    wlog=dlog(omrsq)
    rlog=dlog(rsq)
    Kfun=1d0/rsq*wlog
    epinv2 = 0d0
    epinv  = 0d0
    if( rsq.gt.1d0 ) call Error("r2>1 in virtual correction to the top decay")
!   virtual corrections
    F0 = - epinv2 -epinv*(2.5d0-2d0*wlog)  &
         - 5.5d0 - pisqo6 - 2d0*DLi2(rsq)  &  ! -5.5 corresponds to eta=0 (FHD scheme)
         + 3d0*wlog - 2d0*wlog**2 - Kfun
    F1 = 2d0*Kfun
!   integrated dipoles
    if( (TopDecays.eq.6 .and. TopQuark%PartType.eq.Top_) .or. (TopDecays.eq.5 .and. TopQuark%PartType.eq.ATop_) ) then
        xx = xIntDip(1)
        z  = xIntDip(2)
        FF_xxz = FF(xx/z)
        FF_xx  = FF(xx)
        Cdip = (epinv2-2d0*wlog*epinv+2d0*wlog**2 + 2d0-DblPi**2/6.0)*FF_xx  &
             - (epinv-2d0*wlog)*( 2d0/(1d0-z)*(FF_xxz/z - FF_xx)       &
             - (1d0+z)/z*FF_xxz -  FF_xx)                              &
             +  eta*(1d0-z)*FF_xxz/z                                    &
             + 4d0*dlog(1d0-z)/(1d0-z)*(FF_xxz/z - FF_xx)              &
             + FF_xxz/z*(-2d0*(1d0+z)*dlog(1d0-z)                       &
             - (2d0 /(1d0-z)-(1d0+z))*dlog((rsq+z*(1d0-rsq))/z**2) )      &
             - 2d0/(1d0-z)*(z/(rsq+z*(1d0-rsq))*FF_xxz/z - FF_xx)

        C_ct = (epinv+dlog(m_top**2/MuFrag**2))*( FF_xxz/z*(2d0/(1d0-z) - 1d0-z) &
             - FF_xx*(2d0/(1d0-z)- 3d0/2d0) )
        IntDip = 1d0*( Cdip + C_ct )
        F0 = F0 * FF_xx
        F1 = F1 * FF_xx
    else
        IntDip = epinv2 + epinv*(2.5d0-2d0*wlog)  &
           + 25d0/4d0 + 0.5d0*(1d0/omrsq**2-8d0/omrsq+7d0)*rlog  &
           + 0.5d0/omrsq + 2d0*DLi2(omrsq) - 5d0*pisqo6  &
           - 5d0*wlog + 2d0*wlog**2 &!    eta=0 (FHD scheme)
           - ( 2d0*dlog(alpha_DKTfi)**2 -dlog(alpha_DKTfi)*(-7d0/2d0+4d0*alpha_DKTfi-alpha_DKTfi**2/2d0) &
           - 2d0*(1d0-alpha_DKTfi)*rsq/(omrsq)*dlog(rsq/(1d0-alpha_DKTfi+rsq*alpha_DKTfi)) )
    endif

!DEC$ IF(_FactCheck .EQ.1)
      if( dmin1(dabs(Mom(1:4,2).dot.Mom(1:4,4))/m_W**2,dabs(Mom(1:4,3).dot.Mom(1:4,4))/m_W**2) .lt. 1d-3 ) then
          TopQuark%Pol(1:4) = (1d-60,0d0)
          return
      endif
!DEC$ ENDIF

    if( TopQuark%PartType.eq.Top_ ) then ! Top quark decay
        call ubarSpi_Weyl(dcmplx(Mom(1:4,1)),-1,BotSpi(1:4))
        call ampl_w_qbqg((/Mom(1:4,2),Mom(1:4,3),Mom(1:4,4)/),GluonHel,Wcurr(1:4))
        Wcurr = Wcurr * WProp*g2_weak/dsqrt(2d0) * dsqrt(alpha_s4Pi*CF*RunAlphaS(NLOParam,MuRen))
        if( TOPDECAYS.EQ.2 .or. TOPDECAYS.EQ.3 .or. TOPDECAYS.EQ.4 ) Wcurr(1:4) = Wcurr(1:4) * dsqrt(NFlav*Nc)

!       connect to quark current
        TopQuark%Pol(1:4) = (F0+IntDip)*vgq_Weyl( WCurr(1:4),BotSpi(1:4) ) ! vgq introduces -i/Sqrt(2)
        TopQuark%Pol(1:4) = TopQuark%Pol(1:4) + BotSpi(1:4)*F1*(0d0,-1d0)/dsqrt(2d0)*sc_(dcmplx(Mom(1:4,1)),WCurr(1:4))/m_Top ! introduce -i/Sqrt(2) by hand
        TopQuark%Pol(1:4) =( spb2_Weyl(TopQuark%Pol(1:4),TopQuark%Mom(1:4)) + m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top * alpha_sOver2Pi*RunAlphaS(NLOParam,MuRen) * CF
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)

    elseif( TopQuark%PartType.eq.ATop_ ) then  ! Anti-Top quark decay
        call    vSpi_Weyl(dcmplx(Mom(1:4,1)),+1,BotSpi(1:4))
        call ampl_w_qbqg((/Mom(1:4,3),Mom(1:4,2),Mom(1:4,4)/),GluonHel,Wcurr(1:4))
        Wcurr = Wcurr * WProp*g2_weak/dsqrt(2d0) * dsqrt(alpha_s4Pi*CF*RunAlphaS(NLOParam,MuRen))
        if( TOPDECAYS.EQ.2 .or. TOPDECAYS.EQ.3 .or. TOPDECAYS.EQ.4 ) Wcurr(1:4) = Wcurr(1:4) * dsqrt(NFlav*Nc)

!       connect to quark current:
        TopQuark%Pol(1:4) = (F0+IntDip)*vbqg_Weyl( BotSpi(1:4),WCurr(1:4) )! vbqg introduces -i/Sqrt(2)
        ! INTORODUCED A MINUS SIGN FOR F1 CONTRIBUTION
        TopQuark%Pol(1:4) = TopQuark%Pol(1:4) - BotSpi(1:4)*F1*(0d0,-1d0)/dsqrt(2d0)*sc_(dcmplx(Mom(1:4,1)),WCurr(1:4))/m_Top ! introduce -i/Sqrt(2) by hand
        TopQuark%Pol(1:4) = ( spi2_Weyl(TopQuark%Mom(1:4),TopQuark%Pol(1:4)) - m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top * alpha_sOver2Pi*RunAlphaS(NLOParam,MuRen) * CF
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)
    endif





!-----------------------------------------------------------------------------
ELSEIF( Topol.eq.DKJ_1LQ_Q ) THEN! 1-loop correction on W+int.dipoles and gluon emission from W
!-----------------------------------------------------------------------------
call Wto3Jet_init
        if( init(DK_LO) ) then
        allocate( TopDKAmp(0)%NumGlu(0:2))
        allocate( TopDKAmp(0)%PartRef(1:3))
        allocate( TopDKAmp(0)%PartType(1:3))
        allocate( TopDKAmp(0)%Quarks(1:2))
        allocate( TopDKAmp(0)%Gluons(1:0))
        allocate( ATopDKAmp(0)%NumGlu(0:2))
        allocate( ATopDKAmp(0)%PartRef(1:3))
        allocate( ATopDKAmp(0)%PartType(1:3))
        allocate( ATopDKAmp(0)%Quarks(1:2))
        allocate( ATopDKAmp(0)%Gluons(1:0))

        TopDKAmp(0)%NumPart =3; TopDKAmp(0)%NumQua =2; TopDKAmp(0)%NumGlu(0) = 0
        ATopDKAmp(0)%NumPart=3; ATopDKAmp(0)%NumQua=2; ATopDKAmp(0)%NumGlu(0)= 0
        TopDKAmp(0)%PartRef(1:3) = (/1,2,3/)
        ATopDKAmp(0)%PartRef(1:3)= (/1,2,3/)
        call LinkTreeParticles( TopDKAmp(0), TopDKProd(1:3))
        call LinkTreeParticles(ATopDKAmp(0),ATopDKProd(1:3))
        init(DK_LO)=.false.
    endif

!DEC$ IF(_FactCheck .EQ.1)
      if( dmin1(dabs(Mom(1:4,2).dot.Mom(1:4,4))/m_W**2,dabs(Mom(1:4,3).dot.Mom(1:4,4))/m_W**2) .lt. 1d-3 ) then
          TopQuark%Pol(1:4) = (1d-60,0d0)
          return
      endif
!DEC$ ENDIF

! wsv (please leave that as a bookmark to search for this place)
! REMEMBER for Atop you have to change the Momneta of the quark and anti-quark to coincide with Markus Momenta-ordereing, NOT for Top
    xe = 0
    LO_NLO = 1
    if( TopQuark%PartType.eq.Top_ ) then ! Top quark decay
        TopDKProd(1)%Mom(1:4) =-TopQuark%Mom(1:4)!  top
        TopDKProd(2)%Mom(1:4) = dcmplx(Mom(1:4,1))! bot
        TopDKProd(3)%Mom(1:4) = dcmplx(Mom(1:4,2)+Mom(1:4,3)+Mom(1:4,4))!  W+
        call ubarSpi_Weyl(TopDKProd(2)%Mom(1:4),-1,TopDKProd(2)%Pol(1:4))
!         MomSV(1:4,1:3) =Mom(1:4,2:4)
! !        MomSV(1:4,1) = Mom(1:4,3)
! !        MomSV(1:4,2) = Mom(1:4,2)
        call phase_space_point_virtual((/Mom(1:4,2),Mom(1:4,3),Mom(1:4,4)/),GluonHel,LO_NLO,xe,MuRen,Wcurr(1:4))
        Wmom(1:4) = Mom(1:4,2)+Mom(1:4,3)+Mom(1:4,4)
        Wcurr(1:4) = (-WCurr(1:4) + (WCurr(1:4).dot.WMom(1:4))/m_W**2 * WMom(1:4))
        TopDKProd(3)%Pol(1:4) = Wcurr* WProp * g_weak/dsqrt(2d0) * dsqrt(2.d0)! x 2 flavor (jets )

        call new_calc_ampl(4,4,0,0,TopDKAmp(0),TopQuark%Pol(1:4))

        TopQuark%Pol(1:4) =( spb2_Weyl(TopQuark%Pol(1:4),TopQuark%Mom(1:4)) + m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)

    elseif( TopQuark%PartType.eq.ATop_ ) then  ! Anti-Top quark decay
        ATopDKProd(1)%Mom(1:4) =-TopQuark%Mom(1:4)!  Atop
        ATopDKProd(2)%Mom(1:4) = dcmplx(Mom(1:4,1))! Abot
        ATopDKProd(3)%Mom(1:4) = dcmplx(Mom(1:4,2)+Mom(1:4,3)+Mom(1:4,4))!  W-
        call vSpi_Weyl(ATopDKProd(2)%Mom(1:4),+1,ATopDKProd(2)%Pol(1:4))
!         MomSV(1:4,1:3) =Mom(1:4,2:4)
!         MomSV(1:4,1  ) = Mom(1:4,3)
!         MomSV(1:4,2)   = Mom(1:4,2)
        call phase_space_point_virtual((/Mom(1:4,3),Mom(1:4,2),Mom(1:4,4)/),GluonHel,LO_NLO,xe,MuRen,Wcurr(1:4))
        Wmom(1:4) = Mom(1:4,2)+Mom(1:4,3)+Mom(1:4,4)
!        write(*,*) "****************8"
!        write(*,*) "me-mom",  Mom(1,2),Mom(1,3),Mom(1,4)
!        write(*,*) "me-cur", Wcurr* WProp * g_weak/dsqrt(2d0) * dsqrt(2.d0)

        Wcurr(1:4) = (-WCurr(1:4) + (WCurr(1:4).dot.WMom(1:4))/m_W**2 * WMom(1:4))

        ATopDKProd(3)%Pol(1:4) = Wcurr* WProp * g_weak/dsqrt(2d0) * dsqrt(2.d0) ! x 2 jets

        call new_calc_ampl(4,4,0,0,ATopDKAmp(0),TopQuark%Pol(1:4))

        TopQuark%Pol(1:4) = ( spi2_Weyl(TopQuark%Mom(1:4),TopQuark%Pol(1:4)) - m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)
!        write(*,*) "me-hel", GluonHel, abs(wcurr(1))
!        write(*,*) "me-spi", TopQuark%Pol(1:4)
!        write(*,*) "****"

!        write(*,*) "satop spin", TopQuark%Pol(1:4)
     endif




!-----------------------------------------------------------------------------
ELSEIF( Topol.eq.DKJ_RE_TQ ) THEN! gluon emission from top quark line and from W
!-----------------------------------------------------------------------------
    coupl_sqrt = dsqrt(alpha_s4Pi*CF*RunAlphaS(NLOParam,MuRen))*dsqrt(2d0)  ! dsqrt(2d0) factor is to compensate for 1/dsqrt(2d0) factors in gluon coupling
    if( init(DK_RE_T) ) then
        allocate( TopDKAmp(1)%NumGlu(0:2))
        allocate( TopDKAmp(1)%PartRef(1:4))
        allocate( TopDKAmp(1)%PartType(1:4))
        allocate( TopDKAmp(1)%Quarks(1:2))
        allocate( TopDKAmp(1)%Gluons(1:1))
        allocate( ATopDKAmp(1)%NumGlu(0:2))
        allocate( ATopDKAmp(1)%PartRef(1:4))
        allocate( ATopDKAmp(1)%PartType(1:4))
        allocate( ATopDKAmp(1)%Quarks(1:2))
        allocate( ATopDKAmp(1)%Gluons(1:1))

        TopDKAmp(1)%NumPart=4;   TopDKAmp(1)%NumQua=2; TopDKAmp(1)%NumGlu(0) = 1
        ATopDKAmp(1)%NumPart=4; ATopDKAmp(1)%NumQua=2; ATopDKAmp(1)%NumGlu(0)= 1
        TopDKAmp(1)%PartRef(1:4) = (/1,2,4,3/)
        ATopDKAmp(1)%PartRef(1:4)= (/1,2,4,3/)
        call LinkTreeParticles( TopDKAmp(1), TopDKProd(1:4))
        call LinkTreeParticles(ATopDKAmp(1),ATopDKProd(1:4))
        init(DK_RE_T)=.false.
    endif

! ! DEC$ IF(_FactCheck .EQ.1)
!       if( dabs(Mom(1:4,1).dot.Mom(1:4,5))/m_Top**2 .lt. 1d-3 ) then
!             if( dabs(Mom(1:4,2).dot.Mom(1:4,4))/m_W**2 .lt. 1d-3 ) then
!                 TopQuark%Pol(1:4) = (1d-60,0d0)
!                 return
!             endif
!             if( dabs(Mom(1:4,3).dot.Mom(1:4,4))/m_W**2 .lt. 1d-3 ) then
!                 TopQuark%Pol(1:4) = (1d-60,0d0)
!                 return
!             endif
!       endif
!       if( dabs(Mom(1:4,2).dot.Mom(1:4,4))/m_W**2 .lt. 1d-3 ) then
!             if( dabs(Mom(1:4,1).dot.Mom(1:4,5))/m_Top**2 .lt. 1d-3 ) then
!                 TopQuark%Pol(1:4) = (1d-60,0d0)
!                 return
!             endif
!       endif
!       if( dabs(Mom(1:4,3).dot.Mom(1:4,4))/m_W**2 .lt. 1d-3 ) then
!             if( dabs(Mom(1:4,1).dot.Mom(1:4,5))/m_Top**2 .lt. 1d-3 ) then
!                 TopQuark%Pol(1:4) = (1d-60,0d0)
!                 return
!             endif
!       endif
! ! DEC$ ENDIF

    if( TopQuark%PartType.eq.Top_ ) then ! Top quark decay   !  has a minus sign wrt. to older version
        TopDKProd(1)%Mom(1:4) =-TopQuark%Mom(1:4)!  top
        TopDKProd(2)%Mom(1:4) = dcmplx(Mom(1:4,1))! bot
        TopDKProd(3)%Mom(1:4) = dcmplx(Mom(1:4,2)+Mom(1:4,3)+Mom(1:4,4))!  W+
        TopDKProd(4)%Mom(1:4) = dcmplx(Mom(1:4,5))! glu
        call ubarSpi_Weyl(TopDKProd(2)%Mom(1:4),-1,TopDKProd(2)%Pol(1:4))
        call ampl_w_qbqg((/Mom(1:4,2),Mom(1:4,3),Mom(1:4,4)/),GluonHel,Wcurr(1:4))
        TopDKProd(3)%Pol(1:4) = WProp*g2_weak/2d0 * dsqrt(alpha_s4Pi*CF*RunAlphaS(NLOParam,MuRen))* dsqrt(NFlav*Nc) * Wcurr(1:4)
        call pol_mless(TopDKProd(4)%Mom(1:4),Gluon2Hel,TopDKProd(4)%Pol(1:4))
        call new_calc_ampl(4,4,0,0,TopDKAmp(1),TopQuark%Pol(1:4))

        TopQuark%Pol(1:4) =( spb2_Weyl(TopQuark%Pol(1:4),TopQuark%Mom(1:4)) + m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top * coupl_sqrt

        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)

    elseif( TopQuark%PartType.eq.ATop_ ) then  ! Anti-Top quark decay
        ATopDKProd(1)%Mom(1:4) =-TopQuark%Mom(1:4)!  top
        ATopDKProd(2)%Mom(1:4) = dcmplx(Mom(1:4,1))! bot
        ATopDKProd(3)%Mom(1:4) = dcmplx(Mom(1:4,2)+Mom(1:4,3)+Mom(1:4,4))!  W+
        ATopDKProd(4)%Mom(1:4) = dcmplx(Mom(1:4,5))! glu
        call vSpi_Weyl(ATopDKProd(2)%Mom(1:4),+1,ATopDKProd(2)%Pol(1:4))
        call ampl_w_qbqg((/Mom(1:4,3),Mom(1:4,2),Mom(1:4,4)/),GluonHel,Wcurr(1:4))
        ATopDKProd(3)%Pol(1:4) = WProp*g2_weak/2d0 * dsqrt(alpha_s4Pi*CF*RunAlphaS(NLOParam,MuRen))* dsqrt(NFlav*Nc) * Wcurr(1:4)
        call pol_mless(ATopDKProd(4)%Mom(1:4),Gluon2Hel,ATopDKProd(4)%Pol(1:4))
        call new_calc_ampl(4,4,0,0,ATopDKAmp(1),TopQuark%Pol(1:4))

        TopQuark%Pol(1:4) =( spi2_Weyl(TopQuark%Mom(1:4),TopQuark%Pol(1:4)) - m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top * coupl_sqrt

        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)
    endif





!-----------------------------------------------------------------------------
ELSEIF( Topol.eq.DKJ_RE_QQ ) THEN! two gluon emission from W
!-----------------------------------------------------------------------------

! REMEMBER for Atop you have to change the Momneta of the quark and anti-quark to coincide with Markus Momenta-ordereing, NOT for Top
call Wto3Jet_init
    if( init(DK_LO) ) then
       allocate( TopDKAmp(0)%NumGlu(0:2))
       allocate( TopDKAmp(0)%PartRef(1:3))
       allocate( TopDKAmp(0)%PartType(1:3))
       allocate( TopDKAmp(0)%Quarks(1:2))
       allocate( TopDKAmp(0)%Gluons(1:0))
       allocate( ATopDKAmp(0)%NumGlu(0:2))
       allocate( ATopDKAmp(0)%PartRef(1:3))
       allocate( ATopDKAmp(0)%PartType(1:3))
       allocate( ATopDKAmp(0)%Quarks(1:2))
       allocate( ATopDKAmp(0)%Gluons(1:0))

       TopDKAmp(0)%NumPart =3; TopDKAmp(0)%NumQua =2; TopDKAmp(0)%NumGlu(0) = 0
       ATopDKAmp(0)%NumPart=3; ATopDKAmp(0)%NumQua=2; ATopDKAmp(0)%NumGlu(0)= 0
       TopDKAmp(0)%PartRef(1:3) = (/1,2,3/)
       ATopDKAmp(0)%PartRef(1:3)= (/1,2,3/)
       call LinkTreeParticles( TopDKAmp(0), TopDKProd(1:3))
       call LinkTreeParticles(ATopDKAmp(0),ATopDKProd(1:3))
       init(DK_LO)=.false.
    endif

! !DEC$ IF(_FactCheck .EQ.1)
!   if( (Mom(1:4,4).dot.Mom(1:4,5)).lt.1d-2*m_W**2 ) then
!       if( ((Mom(1:4,4)+Mom(1:4,5)).dot.Mom(1:4,2)).lt.1d-2*m_W**2 ) then
!           TopQuark%Pol(1:4) = (1d-60,0d0)
!           return
!       endif
!       if( ((Mom(1:4,4)+Mom(1:4,5)).dot.Mom(1:4,3)).lt.1d-2*m_W**2 ) then
!           TopQuark%Pol(1:4) = (1d-60,0d0)
!           return
!       endif
!   endif
!   if( (Mom(1:4,2).dot.Mom(1:4,4)).lt.1d-2*m_W**2 ) then
!       if( ((Mom(1:4,2)+Mom(1:4,4)).dot.Mom(1:4,5)).lt.1d-2*m_W**2 ) then
!           TopQuark%Pol(1:4) = (1d-60,0d0)
!           return
!       endif
!       if( (Mom(1:4,5).dot.Mom(1:4,3)).lt.1d-2*m_W**2 ) then
!           TopQuark%Pol(1:4) = (1d-60,0d0)
!           return
!       endif
!   endif
!   if( (Mom(1:4,2).dot.Mom(1:4,5)).lt.1d-2*m_W**2 ) then
!       if( ((Mom(1:4,2)+Mom(1:4,5)).dot.Mom(1:4,4)).lt.1d-2*m_W**2 ) then
!           TopQuark%Pol(1:4) = (1d-60,0d0)
!           return
!       endif
!       if( (Mom(1:4,4).dot.Mom(1:4,3)).lt.1d-2*m_W**2 ) then
!           TopQuark%Pol(1:4) = (1d-60,0d0)
!           return
!       endif
!   endif
!   if( (Mom(1:4,3).dot.Mom(1:4,4)).lt.1d-2*m_W**2 ) then
!       if( ((Mom(1:4,3)+Mom(1:4,4)).dot.Mom(1:4,5)).lt.1d-2*m_W**2 ) then
!           TopQuark%Pol(1:4) = (1d-60,0d0)
!           return
!       endif
!       if( (Mom(1:4,5).dot.Mom(1:4,2)).lt.1d-2*m_W**2 ) then
!           TopQuark%Pol(1:4) = (1d-60,0d0)
!           return
!       endif
!   endif
!   if( (Mom(1:4,3).dot.Mom(1:4,5)).lt.1d-2*m_W**2 ) then
!       if( ((Mom(1:4,3)+Mom(1:4,5)).dot.Mom(1:4,4)).lt.1d-2*m_W**2 ) then
!           TopQuark%Pol(1:4) = (1d-60,0d0)
!           return
!       endif
!       if( (Mom(1:4,4).dot.Mom(1:4,2)).lt.1d-2*m_W**2 ) then
!           TopQuark%Pol(1:4) = (1d-60,0d0)
!           return
!       endif
!   endif
! !DEC$ ENDIF

    if( TopQuark%PartType.eq.Top_ ) then ! Top quark decay
        TopDKProd(1)%Mom(1:4) =-TopQuark%Mom(1:4)!  top
        TopDKProd(2)%Mom(1:4) = dcmplx(Mom(1:4,1))! bot
        TopDKProd(3)%Mom(1:4) = dcmplx(Mom(1:4,2)+Mom(1:4,3)+Mom(1:4,4)+Mom(1:4,5))!  W+
        call ubarSpi_Weyl(TopDKProd(2)%Mom(1:4),-1,TopDKProd(2)%Pol(1:4))
        if(PartAmp .ne. 1 .and. PartAmp .ne. 2) then
           write(*,*) "error: search ofr this in mod_Kinematics"
        endif
!         MomSR(1:4,1:4) = Mom(1:4,2:5)
        call CalcRealWm(Mom(1:4,2:5),GluonHel,Gluon2Hel,PartAmp,Wcurr)
        Wmom(1:4) = Mom(1:4,2)+Mom(1:4,3)+Mom(1:4,4) +Mom(1:4,5)
        Wcurr(1:4) = (-WCurr(1:4) + (WCurr(1:4).dot.WMom(1:4))/m_W**2 * WMom(1:4))
        TopDKProd(3)%Pol(1:4) = Wcurr* WProp * g_weak/dsqrt(2d0) * dsqrt(2.d0) ! x 2 jets
        call new_calc_ampl(4,4,0,0,TopDKAmp(0),TopQuark%Pol(1:4))

        TopQuark%Pol(1:4) =( spb2_Weyl(TopQuark%Pol(1:4),TopQuark%Mom(1:4)) + m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)

    elseif( TopQuark%PartType.eq.ATop_ ) then  ! Anti-Top quark decay
        ATopDKProd(1)%Mom(1:4) =-TopQuark%Mom(1:4)!  Atop
        ATopDKProd(2)%Mom(1:4) = dcmplx(Mom(1:4,1))! Abot
        ATopDKProd(3)%Mom(1:4) = dcmplx(Mom(1:4,2)+Mom(1:4,3)+Mom(1:4,4)+Mom(1:4,5))!  W-
        call vSpi_Weyl(ATopDKProd(2)%Mom(1:4),+1,ATopDKProd(2)%Pol(1:4))
        if(PartAmp .ne. 1 .and. PartAmp .ne. 2) then
           write(*,*) "error: search ofr this in mod_Kinematics"
        endif

!         MomSR(1:4,1:4) = Mom(1:4,2:5)
!         MomSR(1:4,1) = Mom(1:4,3)
!         MomSR(1:4,2) = Mom(1:4,2)
        call CalcRealWm((/Mom(1:4,3),Mom(1:4,2),Mom(1:4,4),Mom(1:4,5)/),GluonHel,Gluon2Hel,PartAmp,Wcurr)
        Wmom(1:4) = Mom(1:4,2)+Mom(1:4,3)+Mom(1:4,4) +Mom(1:4,5)
        Wcurr(1:4) = (-WCurr(1:4) + (WCurr(1:4).dot.WMom(1:4))/m_W**2 * WMom(1:4))
        ATopDKProd(3)%Pol(1:4) = Wcurr* WProp * g_weak/dsqrt(2d0)* dsqrt(2.d0) ! x 2 jets

        call new_calc_ampl(4,4,0,0,ATopDKAmp(0),TopQuark%Pol(1:4))

        TopQuark%Pol(1:4) = ( spi2_Weyl(TopQuark%Mom(1:4),TopQuark%Pol(1:4)) - m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)

    endif




ELSEIF( Topol.eq.DKJ_RE2_QQ ) THEN! 4 quark emissions from  W+/-, non-identical quarks   W -> (ud+cc/ss/bb) and (cs+uu/dd/bb)
    call Wto3Jet_init
    if( init(DK_LO) ) then
       allocate( TopDKAmp(0)%NumGlu(0:2))
       allocate( TopDKAmp(0)%PartRef(1:3))
       allocate( TopDKAmp(0)%PartType(1:3))
       allocate( TopDKAmp(0)%Quarks(1:2))
       allocate( TopDKAmp(0)%Gluons(1:0))
       allocate( ATopDKAmp(0)%NumGlu(0:2))
       allocate( ATopDKAmp(0)%PartRef(1:3))
       allocate( ATopDKAmp(0)%PartType(1:3))
       allocate( ATopDKAmp(0)%Quarks(1:2))
       allocate( ATopDKAmp(0)%Gluons(1:0))

       TopDKAmp(0)%NumPart =3; TopDKAmp(0)%NumQua =2; TopDKAmp(0)%NumGlu(0) = 0
       ATopDKAmp(0)%NumPart=3; ATopDKAmp(0)%NumQua=2; ATopDKAmp(0)%NumGlu(0)= 0
       TopDKAmp(0)%PartRef(1:3) = (/1,2,3/)
       ATopDKAmp(0)%PartRef(1:3)= (/1,2,3/)
       call LinkTreeParticles( TopDKAmp(0), TopDKProd(1:3))
       call LinkTreeParticles(ATopDKAmp(0),ATopDKProd(1:3))
       init(DK_LO)=.false.
    endif

! !DEC$ IF(_FactCheck .EQ.1)
!   if( (Mom(1:4,4).dot.Mom(1:4,5)).lt.1d-3*m_W**2 ) then
!       if( ((Mom(1:4,4)+Mom(1:4,5)).dot.Mom(1:4,2)).lt.1d-3*m_W**2 ) then
!           TopQuark%Pol(1:4) = (1d-60,0d0)
!           return
!       endif
!       if( ((Mom(1:4,4)+Mom(1:4,5)).dot.Mom(1:4,3)).lt.1d-3*m_W**2 ) then
!           TopQuark%Pol(1:4) = (1d-60,0d0)
!           return
!       endif
!   endif
! !DEC$ ENDIF

    if( TopQuark%PartType.eq.Top_ ) then ! Top quark decay
        TopDKProd(1)%Mom(1:4) =-TopQuark%Mom(1:4)!  top
        TopDKProd(2)%Mom(1:4) = dcmplx(Mom(1:4,1))! bot
        TopDKProd(3)%Mom(1:4) = dcmplx(Mom(1:4,2)+Mom(1:4,3)+Mom(1:4,4)+Mom(1:4,5))!  W+
        call ubarSpi_Weyl(TopDKProd(2)%Mom(1:4),-1,TopDKProd(2)%Pol(1:4))

        call ampl_w_qbqffb((/Mom(1:4,2),Mom(1:4,3),Mom(1:4,4),Mom(1:4,5)/),GluonHel,Wcurr(1:4))!     s-channel
        Wcurr(1:4) = Wcurr(1:4) *ne * dsqrt(2.d0)  & ! * correct for use ob vbqq_weyl
                    * g_weak/dsqrt(2.d0) * (RunAlphaS(NLOParam,MuRen)*alpha_S*4.d0*Pi) !1  W-coupling & 2  QCD-vertices
        TopDKProd(3)%Pol(1:4) = Wcurr* WProp * g_weak/dsqrt(2d0) !*dsqrt(3.d0) * dsqrt(2.d0) ! x 3 qqb flavors, x 2 q'qb flavors  -> multiply this outside, now
        call new_calc_ampl(4,4,0,0,TopDKAmp(0),TopQuark%Pol(1:4))

        TopQuark%Pol(1:4) =( spb2_Weyl(TopQuark%Pol(1:4),TopQuark%Mom(1:4)) + m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)
    elseif( TopQuark%PartType.eq.ATop_ ) then  ! Anti-Top quark decay
        ATopDKProd(1)%Mom(1:4) =-TopQuark%Mom(1:4)!  Atop
        ATopDKProd(2)%Mom(1:4) = dcmplx(Mom(1:4,1))! Abot
        ATopDKProd(3)%Mom(1:4) = dcmplx(Mom(1:4,2)+Mom(1:4,3)+Mom(1:4,4)+Mom(1:4,5))!  W-
        call vSpi_Weyl(ATopDKProd(2)%Mom(1:4),+1,ATopDKProd(2)%Pol(1:4))

        call ampl_w_qbqffb((/Mom(1:4,2),Mom(1:4,3),Mom(1:4,4),Mom(1:4,5)/),GluonHel,Wcurr(1:4))!     s-channel
        Wcurr(1:4) = Wcurr(1:4) *ne * dsqrt(2.d0)  & ! * correct for use ob vbqq_weyl
                    * g_weak/dsqrt(2.d0) * (RunAlphaS(NLOParam,MuRen)*alpha_S*4.d0*Pi) !1  W-coupling & 2  QCD-vertices
        ATopDKProd(3)%Pol(1:4) = Wcurr* WProp * g_weak/dsqrt(2d0) !*dsqrt(3.d0) * dsqrt(2.d0)! 3 quark flavors, x 2 q'qb flavors -> multiply this outside, now
        call new_calc_ampl(4,4,0,0,ATopDKAmp(0),TopQuark%Pol(1:4))

        TopQuark%Pol(1:4) = ( spi2_Weyl(TopQuark%Mom(1:4),TopQuark%Pol(1:4)) - m_Top*TopQuark%Pol(1:4) ) * NWAFactor_Top
        TopQuark%Pol(1:4) = WeylToDirac(TopQuark%Pol(1:4))
        TopQuark%Pol(5:16)= (0d0,0d0)
    endif





!         dummy(2) = (0d0,0d0)
!         do hel_pho=-1,1,1!  W-pol loop: -1,0,+1
!         do hel_glu=-1,1,2!  fermion line loop: +1,-1
!               call ampl_w_qbqffb((/Mom(1:4,2),Mom(1:4,3),Mom(1:4,4),Mom(1:4,5)/),hel_glu,Wcurr(1:4))!     s-channel
!               WPol(1:4) = pol_mass(dcmplx(Mom(1:4,2)+Mom(1:4,3)+Mom(1:4,4)+Mom(1:4,5)),m_W,hel_pho)
!               dummy(1) = WPol(1:4).dot.Wcurr(1:4)
!               dummy(2) = dummy(2) + dummy(1)*dconjg(dummy(1))
!         enddo
!         enddo
!         dummy(2) = dummy(2) *(4d0*DblPi/137d0)/(2d0*0.2149577874691558d0) * (4d0*DblPi*0.13d0)**2*8d0 /3d0       /100d0**2/2d0
!
!
!        MG_MOM(0:3,1) = (Mom(1:4,2)+Mom(1:4,3)+Mom(1:4,4)+Mom(1:4,5))*100d0
!        MG_MOM(0:3,2) = Mom(1:4,2)*100d0
!        MG_MOM(0:3,3) = Mom(1:4,3)*100d0
!        MG_MOM(0:3,4) = Mom(1:4,4)*100d0
!        MG_MOM(0:3,5) = Mom(1:4,5)*100d0
!        call coupsm(0)
!        call Swm_ubdccb(MG_MOM,MadGraph_tree)
!        MadGraph_tree= MadGraph_tree
!
! print *, "me", dummy(2)
! print *, "MG", MadGraph_tree
! print *, "ratio",dummy(2)/MadGraph_tree-1d0
! pause


! !         dummy(3) = (0d0,0d0)
! !         do hel_pho=-1,1,1!  W-pol loop: -1,0,+1
! !         do hel_glu=-1,1,2!  fermion line loop: +1,-1
! !               call ampl_w_qbqffb((/Mom(1:4,2),Mom(1:4,3),Mom(1:4,4),Mom(1:4,5)/),hel_glu,Wcurr(1:4))!     s-channel
! !               call ampl_w_qbqffb((/Mom(1:4,5),Mom(1:4,3),Mom(1:4,4),Mom(1:4,2)/),hel_glu,Wcurr_new(1:4))! t-channel
! !               WPol(1:4) = pol_mass(dcmplx(Mom(1:4,2)+Mom(1:4,3)+Mom(1:4,4)+Mom(1:4,5)),m_W,hel_pho)
! !               dummy(1) = WPol(1:4).dot.Wcurr(1:4)
! !               dummy(2) = WPol(1:4).dot.Wcurr_new(1:4)
! !               dummy(3) = dummy(3) + 8d0*dummy(1)*dconjg(dummy(1)) + 8d0*dummy(2)*dconjg(dummy(2))
! !               if( hel_glu.eq.1 ) dummy(3) = dummy(3) + (8d0/3d0)*dummy(1)*dconjg(dummy(2)) + (8d0/3d0)*dummy(2)*dconjg(dummy(1))
! !         enddo
! !         enddo
! !         dummy(3) = dummy(3) *(4d0*DblPi/137d0)/(2d0*0.2149577874691558d0) * (4d0*DblPi*0.13d0)**2 /3d0 /2d0    /100d0**2/2d0
! !
! !        MG_MOM(0:3,1) = (Mom(1:4,2)+Mom(1:4,3)+Mom(1:4,4)+Mom(1:4,5))*100d0
! !        MG_MOM(0:3,2:5) = Mom(1:4,2:5)*100d0
! !        call coupsm(0); call Swm_ubduub(MG_MOM,MadGraph_tree)
! !
! !        print *, "me", dummy(3)
! !        print *, "MG", MadGraph_tree
! !        print *, "ratio",dummy(3)/MadGraph_tree-1d0
! !        pause





ELSE
      call Error("this top decay topology is not implemented",Topol)

ENDIF

return
END SUBROUTINE




END MODULE
