MODULE ModParameters
implicit none


integer, public :: Collider, Process, PDFSet, NLOParam, TopDecays, XTopDecays
integer, public :: MasterProcess, Correction, ObsSet
integer, public :: VegasIt0,VegasIt1,VegasNc0,VegasNc1,VegasMxDim,VegasSeed
integer, public :: VegasIt0_default,VegasIt1_default,VegasNc0_default,VegasNc1_default
integer, public :: NumHistograms
integer, public :: Fragm_Func_Type
real(8), public :: alpha_frag,beta_frag,delta_frag
real(8), public :: Lambda_QCD
logical, public :: unweighted
logical, public :: HelSampling
integer, public :: DKRE_switch
integer(8), public, save :: EvalCounter=0
integer(8), public, save :: PSCutCounter=0
integer(8), public, save :: SkipCounter=0
real(8),public,save :: maxWgt=0d0
integer, allocatable :: Crossing(:)
real(8), public :: MuRen, MuFac, MuFrag, AvgFactor
character, public :: HistoFile*(100)
character, public :: GridFile*(100)
integer, public :: GridIO
real(8), public :: AvgValue=0d0,MinValue=1d13,MaxValue=-1d13
real(8), public :: time_start,time_end

real(8), public, parameter :: DblPi = 3.1415926535897932384626433832795028842d0
real(8), public, parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0
real(8), public, parameter :: GeV=0.01d0


real(8), public, parameter :: alpha = 1d0/(137d0)
real(8), public, parameter :: alpha4Pi = alpha*4d0*DblPi
real(8), public, parameter :: GF = (1.16639d-5)/GeV**2
real(8), public            :: m_Top, m_SMTop
real(8), public            :: m_Bot
real(8), public, parameter :: m_Chm   = 0d0
real(8), public, parameter :: m_Str   = 0d0
real(8), public, parameter :: m_Up    = 0d0
real(8), public, parameter :: m_Dn    = 0d0
real(8), public, parameter :: m_Z     = 91.188d0*GeV
real(8), public, parameter :: m_W     = 80.419d0*GeV
real(8), public, parameter :: m_e     = 0d0
real(8), public, parameter :: m_nu    = 0d0
real(8), public, parameter :: m_HTop  = 600d0*GeV
real(8), public, parameter :: m_A0    = 50d0*GeV
real(8), public, parameter :: m_BH    = 50d0*GeV! remember: changes here require full re-compilation!
real(8), public, parameter :: g2_weak = 4d0*dsqrt(2d0)*m_W**2*GF
real(8), public, parameter :: g_weak = dsqrt(g2_weak)
real(8), public, parameter :: sw = dsqrt(4.d0*DblPi*alpha/g2_weak)
real(8), public, parameter :: sw2 = sw**2 
real(8), public, parameter :: cw = dsqrt(1d0-sw2)
real(8), public, parameter :: EL = dsqrt(4.d0*DblPi*alpha)
real(8), public            :: Ga_HTop
real(8), public            :: Ga_Htop_A0Top
real(8), public            :: Ga_Htop_BHTop
real(8), public, parameter :: m_STop  = 172d0*GeV
real(8), public            :: Ga_STop
real(8), public            :: Ga_Stop_ChiTop
real(8), public, parameter :: m_Chi   = 50d0*GeV

!!! Zprime section !!!

!real(8), public, parameter :: m_Zpr = 91.19d0*GeV ! Standard Z
!real(8), public, parameter :: Ga_Zpr = 2.4950d0*GeV ! Standard Z
real(8), public, parameter :: m_Zpr = 1500d0*GeV
real(8), public, parameter :: Ga_Zpr = 15d0*GeV
real(8), public :: gL_Zpr(6), gR_Zpr(6)

!!! End Zprime section !!!

real(8), public, parameter :: Q_up    = 2d0/3d0
real(8), public, parameter :: Q_dn    =-1d0/3d0
real(8), public, parameter :: Q_el    =-1d0
real(8), public            :: Q_top
real(8), public, parameter :: Q_Wp    =+1d0
real(8), public, parameter :: Q_Wm    =-1d0
real(8), public            :: Q_in
real(8), public :: Ga_Top(0:1),Ga_W(0:1), WidthExpansion
real(8), public :: Ga_TopExp = 1.99d0*GeV
real(8), public :: Ga_WExp   = 2.14d0*GeV
real(8), public, parameter :: fbGeV2=0.389379d12*GeV**2


real(8), public, parameter :: Nf_light=5d0
real(8), public, parameter :: Nf_heavy=1d0
real(8), public, parameter :: Nf=Nf_light+Nf_heavy
real(8), public, parameter :: sumQlight=1d0/3d0

integer, public :: nGluRadContr

real(8), public :: alpha_ii = 1d0
real(8), public :: alpha_if = 1d0
real(8), public :: alpha_fi = 1d0
real(8), public :: alpha_ff = 1d0
real(8), public :: alpha_DK = 1d0
real(8), public :: alpha_DKTfi = 1d0
real(8), public :: alpha_DKTff = 1d0
real(8), public :: alpha_DKWff = 1d0
real(8), public :: kappa_ff = 2d0/3d0  * 0d0

! Jet Algorithm
integer, public            :: AlgoType           ! 1 = kT,   -1 = anti kT
integer, public, parameter :: RecombPrescr = 0   ! 0 = 4-vector addition,   1 = Ellis-Soper prescription

! PDF Set
character, public :: PDFSetString*(100)

! particle 0 = not defined
integer, public, target :: Up_  = 1
integer, public, target :: Dn_  = 2
integer, public, target :: Chm_ = 3
integer, public, target :: Str_ = 4
integer, public, target :: Top_ = 5
integer, public, target :: Bot_ = 6
integer, public, target :: Glu_ = 10
integer, public, target :: Pho_ = 11
integer, public, target :: Z0_  = 12
integer, public, target :: Wp_  = 13
integer, public, target :: HTop_= 14
integer, public, target :: Stop_= 15
integer, public, target :: Sbot_= 16

integer, public, target :: AUp_  = -1
integer, public, target :: ADn_  = -2
integer, public, target :: AChm_ = -3
integer, public, target :: AStr_ = -4
integer, public, target :: ATop_ = -5
integer, public, target :: ABot_ = -6
integer, public, target :: Wm_   = -13
integer, public, target :: AHTop_= -14
integer, public, target :: AStop_= -15
integer, public, target :: ASBot_= -16

integer, public, parameter :: DK_LO=0
integer, public, parameter :: DK_1L_T=1
integer, public, parameter :: DK_RE_T=2
integer, public, parameter :: DK_1L_Q=3
integer, public, parameter :: DK_RE_Q=4
integer, public, parameter :: DKP_LO_T=5
integer, public, parameter :: DKP_LO_L=6
integer, public, parameter :: DKP_1L_T=8
integer, public, parameter :: DKP_RE_T=9
integer, public, parameter :: DKP_1L_L=10
integer, public, parameter :: DKP_RE_L=11

integer, public, parameter :: DKJ_LO_T=DK_RE_T
integer, public, parameter :: DKJ_LO_Q=DK_RE_Q
integer, public, parameter :: DKJ_1L_T=12
integer, public, parameter :: DKJ_1L_Q=13
integer, public, parameter :: DKJ_1LQ_T=14
integer, public, parameter :: DKJ_1LQ_Q=15
integer, public, parameter :: DKJ_RE_TT=16
integer, public, parameter :: DKJ_RE_QQ=17
integer, public, parameter :: DKJ_RE_TQ=18
integer, public, parameter :: DKJ_RE2_TT=19
integer, public, parameter :: DKJ_RE2_QQ=20

integer, public, parameter :: DKX_HTA0_LO=21
integer, public, parameter :: DKX_HTBH_LO=22
integer, public, parameter :: DKX_STChi0_LO=23

integer, public, parameter :: WDK_Lep=1
integer, public, parameter :: WDK_LepPho=2
integer, public, parameter :: WDK_Had=3
integer, public, parameter :: WDK_HadPho=4
integer, public, parameter :: WDK_LO=0

integer, public, parameter :: T_B_W=1
integer, public, parameter :: T_BG_W=2
integer, public, parameter :: T_B_WG=3
integer, public, parameter :: T_BGG_W=4
integer, public, parameter :: T_BG_WG=5
integer, public, parameter :: T_B_WGG=6
integer, public, parameter :: HT_A0_T=7
integer, public, parameter :: HT_BH_T=8
integer, public, parameter :: ST_Chi0_T=9

real(8), public :: IA0Tt(-1:+1)
real(8), public :: IBHTt(-1:+1)
real(8), public :: IChiStt(-1:+1)

real(8), public :: Collider_Energy
real(8), public :: alpha_s, alpha_s4Pi, alpha_sOver2Pi

real(8) :: ColLO_ttbgg(1:2,1:2)
real(8) :: Col1L_ttbgg(1:2,1:9)
real(8) :: Col1L_ttbggp(1:2,1:3)
real(8) :: ColLO_ttbqqb(1:1,1:1)
real(8) :: Col1L_ttbqqb(1:1,1:2)
real(8) :: ColLO_ttbggg(1:6,1:6)
real(8) :: Col1L_ttbggg(1:6,1:35)
real(8) :: ColLO_ttbqqbg(1:4,1:4)
real(8) :: Col1L_ttbqqbg(1:4,1:4)
real(8) :: ColLO_ttbgggg(1:24,1:24)

integer, public :: PrimAmp1_1234,PrimAmp1_1243
integer, public :: PrimAmp1_1324,PrimAmp1_1423
integer, public :: PrimAmp1_1342,PrimAmp1_1432
integer, public :: PrimAmp2_1234,PrimAmp2_1243
integer, public :: PrimAmp2m_1234,PrimAmp2m_1243
integer, public :: PrimAmp3_1432,PrimAmp4_1234
integer, public :: PrimAmp1_12345,PrimAmp1_12354,PrimAmp1_12435,PrimAmp1_12453,PrimAmp1_12534,PrimAmp1_12543
integer, public :: PrimAmp1_13245,PrimAmp1_13254,PrimAmp1_14235,PrimAmp1_14253,PrimAmp1_15234,PrimAmp1_15243
integer, public :: PrimAmp1_13425,PrimAmp1_13524,PrimAmp1_14325,PrimAmp1_14523,PrimAmp1_15324,PrimAmp1_15423
integer, public :: PrimAmp1_13452,PrimAmp1_13542,PrimAmp1_14352,PrimAmp1_14532,PrimAmp1_15342,PrimAmp1_15432
integer, public :: PrimAmp1_162345,PrimAmp1_123645,PrimAmp1_162534,PrimAmp1_125364,PrimAmp1_156234
integer, public :: PrimAmp1_165234,PrimAmp1_152364,PrimAmp1_162354,PrimAmp1_123564,PrimAmp1_123654

integer, public :: PrimAmp3_15432,PrimAmp3_14532,PrimAmp3_14352,PrimAmp3_14325
integer, public :: PrimAmp4_12354,PrimAmp4_12534,PrimAmp4_15234,PrimAmp4_12345
integer, public :: PrimAmp2_12345,PrimAmp2_15234,PrimAmp2_12534,PrimAmp2_12354
integer, public :: PrimAmp2_15243,PrimAmp2m_15243,PrimAmp2_12543,PrimAmp2_12453,PrimAmp2m_12435,PrimAmp2_12435,PrimAmp2m_12543,PrimAmp2m_12453
integer, public :: PrimAmp2m_12345,PrimAmp2m_15234,PrimAmp2m_12534,PrimAmp2m_12354



      real(8), public, parameter :: pi =3.141592653589793238462643383279502884197d0
      real(8), public, parameter :: pisq = pi**2

       real(8), public, parameter :: one = 1.0d0, mone = -1.0d0
       real(8), public, parameter :: half  = 0.5d0,two = 2.0d0
       real(8), public, parameter :: zero  = 0.0d0
       real(8), public, parameter :: three = 3.d0
       real(8), public, parameter :: four = 4.d0
       real(8), public, parameter :: msqrt2=-1.41421356237309504880168872420969807856967d0
       real(8), public, parameter :: sqrt3 =1.7320508075688772935274463415058723669428d0
       complex(8), parameter, public :: czero = (0d0,0d0)
       complex(8), parameter, public :: cone = 1.0d0
       complex(8), parameter, public :: ci=(0.0d0,1.0d0)
       complex(8), parameter, public :: ne=(0.0d0,1.0d0)


contains


SUBROUTINE InitKirillParameters
implicit none
include './misc/global_import'
include './misc/ucom_import'
        complex(8) one,cone,mcone,mone,sqrt2,msqrt2,csqrt2,cmsqrt2
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

!       mt=m_Top
!       mb=m_Bot
!       mc=m_Chm


END SUBROUTINE



SUBROUTINE InitParameters
use ModMisc
implicit none
real(8) :: r2, TopWidthExpansion,WWidthExpansion,WWidthChoice


m_Bot   = m_Top ! this is NOT the bottom mass! it is the mass for massive particles in closed fermion loops


IF( COLLIDER.EQ.1 ) THEN
   Collider_Energy  = 14000d0*GeV
   if( ObsSet.eq.1 ) then
        Collider_Energy  = 14000d0*GeV
   endif
   if( ObsSet.eq.3 ) then
        Collider_Energy  = 7000d0*GeV
   endif
   if( ObsSet.eq.5 ) then
        Collider_Energy  = 7000d0*GeV
   endif
   if( ObsSet.eq.6 ) then
        Collider_Energy  = 7000d0*GeV
   endif
   if( ObsSet.eq.8 ) then
        Collider_Energy  = 7000d0*GeV
   endif
   if( ObsSet.eq.11 ) then
        Collider_Energy  = 14000d0*GeV
   endif
   if( ObsSet.eq.13 ) then
        Collider_Energy  = 7000d0*GeV
   endif
   if( ObsSet.eq.15 ) then
        Collider_Energy  = 7000d0*GeV
   endif
   if( ObsSet.eq.25 ) then
        Collider_Energy  = 7000d0*GeV
   endif
ELSEIF( COLLIDER.EQ.2 ) THEN
  Collider_Energy  = 1960d0*GeV
ENDIF


IF( PDFSet.EQ.2 .AND. (NLOPARAM.EQ.1 .OR. NLOPARAM.EQ.0) ) THEN
  Lambda_QCD = 0.165d0*GeV
  alpha_s = 0.13d0  ! CTEQ6L1
!   alpha_s = 0.118d0   ! CTEQ6L
ELSEIF( PDFSet.EQ.2 .AND. (NLOPARAM.EQ.2) ) THEN
  Lambda_QCD = 0.226235d0*GeV
  alpha_s = 0.118d0
ELSEIF(  PDFSet.EQ.1 .AND. (NLOPARAM.EQ.1 .OR. NLOPARAM.EQ.0)) THEN
  Lambda_QCD = 0.167d0/100d0
  alpha_s = 0.13939d0
ELSEIF( PDFSet.EQ.1 .AND. (NLOPARAM.EQ.2) ) THEN
   Lambda_QCD = 0.226235d0/100d0
   alpha_s=0.12018D0
ELSE
  print *, "alpha_s not set"
ENDIF

alpha_s4Pi = alpha_s*4d0*DblPi
alpha_sOver2Pi = alpha_s/2d0/DblPi


!------------ top quark width -----------
!  Ga_Top(NLO) = Ga_Top(0) + Ga_Top(1)  !
!----------------------------------------
r2 = m_W**2/m_Top**2
Ga_Top(0) = GF/8d0/DblPi/dsqrt(2d0)*m_Top**3 * (1d0-r2)**2 * (1d0+2d0*r2)
! Ga_Top(1) =   GF/8d0/DblPi/dsqrt(2d0)*m_Top**3 * alpha_s/DblPi*4d0/3d0*(   &
!               5d0/4d0 - DblPi**2/3d0 + 3d0/2d0*r2 + r2**2*(-6d0+DblPi**2+3d0/2d0*dlog(r2))  &
!              +r2**3*(46d0/9d0 - 2d0/3d0*DblPi**2 - 2d0/3d0*dlog(r2)) )!+O((MW/MT)**8)    ! taken from eq.(14), hep-ph/9906273
Ga_Top(1) = Ga_Top(0) * ( - RUNALPHAS(2,MuRen)*alpha_sOver2Pi*4d0/3d0*(   &   ! this is the one-loop contribution (without LO)
                          2d0/3d0*DblPi**2 + 4d0*DLi2(r2) - 3d0/2d0 +2d0*dlog(1d0/r2-1d0)   &
                         +2d0*dlog(r2)*dlog(1d0-r2) - 4d0/3d0/(1d0-r2) + (22d0-34d0*r2)/9d0/(1d0-r2)**2*dlog(r2)  &
                         +(3d0+27d0*dlog(1d0-r2)-4d0*dlog(r2))/9d0/(1d0+2d0*r2)  &
                        ))   ! taken from eq.(26), hep-ph/0408158



WWidthChoice = 1!          0=experimental W width,    =calculated W width
IF( WWidthChoice.eq. 1 ) THEN
!   calculated W width:
    Ga_W(0) = (2d0*3d0+3d0)*GF*M_W**3/(6d0*dsqrt(2d0)*DblPi)
    Ga_W(1) = Ga_W(0) * RUNALPHAS(2,MuRen)*alpha_s/DblPi
ELSE
!   experimental W width:
    Ga_W(0) = Ga_WExp
    Ga_W(1) = 0d0
ENDIF


TopWidthExpansion = -2d0*Ga_Top(1)/Ga_Top(0) 
WWidthExpansion   = -2d0*Ga_W(1)/Ga_W(0)

IF( abs(TOPDECAYS).GE.1 .AND. NLOPARAM.EQ.2 .AND. CORRECTION.EQ.0 ) THEN
  WidthExpansion = 1d0 + TopWidthExpansion + WWidthExpansion
ELSE
  WidthExpansion = 1d0
ENDIF



!  chiral couplings for stop-Chi^0-top
   IChiStt(+1) = 1d0/50d0  *1d0
   IChiStt(-1) = 1d0/50d0  *1d0

!  stop-->Chi^0 + top partial width!
   Ga_Stop_ChiTop = SqrtLambda(m_stop**2,m_top**2,m_chi**2)/(16d0*DblPi*m_stop**3) * & 
                  ( (IChiStt(+1)**2+IChiStt(-1)**2)*(m_stop**2-m_top**2-m_chi**2) - 4d0*(IChiStt(+1)*IChiStt(-1))*m_top*m_chi )
   Ga_Stop = Ga_Stop_ChiTop! assuming no other decay channel



IF( XTOPDECAYS.EQ.2 ) THEN
!  chiral couplings for HTop-A0-top
!  IA0Tt(+1) = 1d0/50d0  *1d0
!  IA0Tt(-1) = 1d0/50d0  *1d0
   IA0Tt(+1) = M_HTop/M_A0 * 1d0/50d0  *1d0
   IA0Tt(-1) = M_Top/M_A0  * 1d0/50d0  *1d0

!  HTop-->A0 + top partial width!
   Ga_Htop_A0Top = SqrtLambda(m_Htop**2,m_top**2,m_A0**2)/(16d0*DblPi*m_Htop**3) * & 
                   1d0/2d0*( (IA0Tt(+1)**2+IA0Tt(-1)**2)*(m_Htop**2+m_top**2-m_A0**2) + 4d0*(IA0Tt(+1)*IA0Tt(-1))*m_top*m_HTop )
   Ga_HTop = Ga_Htop_A0Top! assuming no other decay channel

ELSEIF( XTOPDECAYS.EQ.1 ) THEN
!  chiral couplings for HTop-BH-top
   IBHTt(+1) = 1d0/50d0  *1d0
   IBHTt(-1) = 1d0/50d0  *1d0

!  HTop-->BH + top partial width!
   Ga_Htop_BHTop =SqrtLambda(m_Htop**2,m_top**2,m_BH**2)/(16d0*DblPi*m_Htop**3) * & 
                  1d0/2d0 *( (IBHTt(+1)**2+IBHTt(-1)**2)*( m_HTop**2+m_Top**2-2d0*m_BH**2+(m_HTop**2-m_top**2)**2/m_BH**2 ) - 12d0*(IBHTt(+1)*IBHTt(-1))*( m_Top*m_HTop ) )
   Ga_HTop = Ga_Htop_BHTop! assuming no other decay channel
ENDIF




#include "LO_ttbgg_output.f90"
#include "1L_ttbgg_output.f90"

#include "LO_ttbqqb_output.f90"
#include "1L_ttbqqb_output.f90"

#include "LO_ttbggg_output.f90"
#include "1L_ttbggg_output.f90"

#include "LO_ttbqqbg_output.f90"
#include "1L_ttbqqbg_output.f90"

#include "LO_ttbgggg_output.f90"

#include "1L_ttbggp_output.f90"


!!! Zprime section !!!

gR_Zpr(top_) = eL*(sw/cw)*2d0/3d0
gL_Zpr(top_) = -eL/(sw*cw)*( 0.5d0 - sw**2*2d0/3d0)

gL_Zpr(dn_) = -eL/(sw*cw)*(-0.5d0 + sw**2/3d0)
gR_Zpr(dn_)  = -eL*(sw/cw)/3d0

gL_Zpr(up_) = gL_Zpr(top_)
gR_Zpr(up_) = gR_Zpr(top_)
gL_Zpr(chm_) = gL_Zpr(top_)
gR_Zpr(chm_) = gR_Zpr(top_)

gL_Zpr(str_) = gL_Zpr(dn_)
gR_Zpr(str_) = gR_Zpr(dn_)
gL_Zpr(bot_) = gL_Zpr(dn_)
gR_Zpr(bot_) = gR_Zpr(dn_)

!!! End Zprime section !!!


END SUBROUTINE



FUNCTION GetMass( PartType )
implicit none
real(8) :: GetMass
integer PartType

   if    ( abs( PartType) .eq. Up_ ) then
      GetMass = m_Up
   elseif( abs( PartType) .eq. Dn_ ) then
      GetMass = m_Dn
   elseif( abs( PartType) .eq. Chm_ ) then
      GetMass = m_Chm
   elseif( abs( PartType) .eq. Str_ ) then
      GetMass = m_Str
   elseif( abs( PartType) .eq. Top_ ) then
      GetMass = m_Top
   elseif( abs( PartType) .eq. Bot_ ) then
      GetMass = m_Bot
   elseif( abs( PartType) .eq. Glu_) then
      GetMass = 0d0
   elseif( abs( PartType) .eq. Pho_) then
      GetMass = 0d0
   elseif( abs( PartType) .eq. Z0_) then
      GetMass = m_Z
   elseif( abs( PartType) .eq. Wp_) then
      GetMass = m_W
   elseif( abs( PartType) .eq. HTop_) then
      GetMass = m_HTop
   elseif( abs( PartType) .eq. STop_) then
      GetMass = m_Stop
   elseif( PartType .eq. 0 ) then
      GetMass = 0d0
   else
      print *, "GetMass particle type not available."
   endif
RETURN
END FUNCTION



FUNCTION RunAlphaS(Loop,Q)! for alphas(MZ) and Nf=5
implicit none
integer :: Loop
real(8) :: Q,w,RunAlphaS
integer, parameter :: NF=5
real(8), parameter :: beta0=11d0-2d0/3d0*NF
real(8) :: beta1=17d0*3d0-4d0/3d0*NF-5d0*NF



! !    MCFM version
    if( Loop.eq.1 ) then
      w = 1d0 - alpha_sOver2Pi*beta0*dlog(m_Z/Q)   ! one loop running
      RunAlphaS = 1d0/w
    elseif( Loop.eq.2 ) then
      w = 1d0 - alpha_sOver2Pi*beta0*dlog(m_Z/Q)
      RunAlphaS = 1d0/w * (1d0-beta1/beta0*alpha_sOver2Pi*dlog(w)/w)   ! two loop running
    else
      RunAlphaS = 1d0   ! no running
    endif


!print *, "USING DUW Alpha_S"

! ! ! !   DUW version
!     if( Loop.eq.1 ) then
!       RunAlphaS = 4d0*dblpi/beta0/dlog((Q*100d0)**2/165d-3**2)
!       RunAlphaS = RunAlphaS/alpha_s
!     elseif( Loop.eq.2 ) then
!       RunAlphaS = 4d0*dblpi/beta0/dlog((Q*100d0)**2/226d-3**2)*(1d0-(beta1*2d0)/beta0**2*dlog(dlog((Q*100d0)**2/226d-3**2))/dlog((Q*100d0)**2/226d-3**2))
!       RunAlphaS = RunAlphaS/alpha_s
!     else
!       RunAlphaS = 1d0   ! no running
!     endif


RETURN
END FUNCTION

END MODULE



