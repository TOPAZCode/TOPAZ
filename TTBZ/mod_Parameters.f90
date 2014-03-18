MODULE ModParameters
implicit none


integer, public :: Collider, Process, PDFSet, NLOParam, TopDecays, XTopDecays, ZDecays
integer, public :: MasterProcess, Correction, ObsSet
integer, public :: VegasIt0,VegasIt1,VegasNc0,VegasNc1,VegasMxDim,VegasSeed
integer, public :: VegasIt0_default,VegasIt1_default,VegasNc0_default,VegasNc1_default
integer, public :: NumHistograms
integer, public :: Fragm_Func_Type
real(8), public :: alpha_frag,beta_frag,delta_frag
real(8), public :: Lambda_QCD
logical, public :: unweighted
logical, public :: HelSampling
logical, public :: FirstLOThenVI
integer, public :: DKRE_switch
integer(8), public, save :: EvalCounter=0
integer(8), public, save :: PSCutCounter=0
integer(8), public, save :: SkipCounter=0
logical,public ,parameter :: TTBZ_SpeedUp=.false.  !  good idea but doesn't work. always set to .false.
integer, public :: TTBZ_DebugSwitch=0   !   0=disabled, 1=bosonic loops, 2=fermionic loops, 3=gamma5 renormalization

! PVegas (MPI) part
integer, public :: MPI_Rank
integer, parameter :: NUMFUNCTIONS=1
integer, parameter :: WORKERS=7
integer, parameter :: PDIM=0


integer,public, save :: pole_skipped=0
integer,public, save :: useQP=0
real(8), public      :: opp_err

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


!orig real(8), public, parameter :: alpha = 1d0/(137d0)
!orig real(8), public, parameter :: alpha4Pi = alpha*4d0*DblPi
real(8), public, parameter :: GF = (1.16639d-5)/GeV**2  !GF = (1.166379d-5)/GeV**2  ! (1.16639d-5)/GeV**2  !
real(8), public            :: m_Top, m_SMTop
real(8), public            :: m_Bot
real(8), public, parameter :: m_Chm   = 0d0
real(8), public, parameter :: m_Str   = 0d0
real(8), public, parameter :: m_Up    = 0d0
real(8), public, parameter :: m_Dn    = 0d0
real(8), public, parameter :: m_Z     = 91.1876*GeV!  91.19*GeV ! 91.18*GeV!     ! 
real(8), public, parameter :: m_W     = 80.385d0*GeV  ! 80.45*GeV  !80.385152965002916d0*GeV  !80.37576d0*GeV !
real(8), public, parameter :: m_e     = 0d0
real(8), public, parameter :: m_nu    = 0d0
real(8), public            :: m_HTop
real(8), public, parameter :: g2_weak = 4d0*dsqrt(2d0)*m_W**2*GF
real(8), public, parameter :: g_weak = dsqrt(g2_weak)
!orig real(8), public, parameter :: sw = dsqrt(4.d0*DblPi*alpha/g2_weak)
!orig real(8), public, parameter :: sw2 = sw**2
real(8), public, parameter :: sw2 = 1d0 - m_W**2/m_Z**2
real(8), public, parameter :: sw = dsqrt(sw2)
real(8), public, parameter :: alpha4Pi = g2_weak*sw2
real(8), public, parameter :: alpha = alpha4Pi/4d0/DblPi
real(8), public, parameter :: cw = dsqrt(1d0-sw2)
real(8), public, parameter :: EL = dsqrt(4.d0*DblPi*alpha)
real(8), public            :: Ga_Top(0:1)
real(8), public            :: Ga_W(0:1)
real(8), public            :: Ga_TopExp = 1.99d0*GeV
real(8), public            :: Ga_WExp   = 2.14d0*GeV
real(8), public            :: Ga_ZExp = 2.4952d0*GeV
real(8), public            :: Ga_HTop(0:1)
real(8), public            :: Ga_Htop_A0Top(0:1)
real(8), public            :: Ga_Htop_BHTop(0:1)
real(8), public            :: m_STop
real(8), public            :: m_SBot
real(8), public            :: Ga_STop(0:1)
real(8), public            :: Ga_Stop_ChiTop(0:1)
real(8), public, parameter :: m_A0    = 50d0*GeV! (scalar)
real(8), public, parameter :: m_BH    = 50d0*GeV! (vector)
real(8), public, parameter :: m_Chi   =100d0*GeV! (Majorana)
real(8), public, parameter :: Vev  = 246d0*GeV

!!!! BSM top-Z couplings
real(8), public :: DeltaF1A, DeltaF1V    !, DeltaF2A, DeltaF2V    ! later
real(8), public :: AbsDelF1A, AbsDelF1V    !, DeltaF2A, DeltaF2V    ! later
real(8), public :: RelDelF1A, RelDelF1V    !, DeltaF2A, DeltaF2V    ! later
!!!

!!! Zprime section !!!

!real(8), public, parameter :: m_Zpr = 91.19d0*GeV ! Standard Z
!real(8), public, parameter :: Ga_Zpr = 2.4950d0*GeV ! Standard Z
real(8), public :: m_Zpr
real(8), public :: Ga_Zpr
real(8), public :: gL_Zpr(6), gR_Zpr(6)
!!! End Zprime section !!!

real(8), public, parameter :: Q_up    = 2d0/3d0
real(8), public, parameter :: Q_dn    =-1d0/3d0
real(8), public, parameter :: Q_el    =-1d0
real(8), public, parameter :: Q_nu    = 0d0
real(8), public            :: Q_top
real(8), public, parameter :: Q_Wp    =+1d0
real(8), public, parameter :: Q_Wm    =-1d0
real(8), public            :: Q_in

real(8), public, parameter :: T3_up    =+1d0/2d0
real(8), public, parameter :: T3_dn    =-1d0/2d0
real(8), public, parameter :: T3_nu    =+1d0/2d0
real(8), public, parameter :: T3_el    =-1d0/2d0

real(8), public, parameter :: couplZUU_left  = -sw/cw*Q_up + 1d0/sw/cw * T3_up
real(8), public, parameter :: couplZUU_right = -sw/cw*Q_up  
real(8), public, parameter :: couplZDD_left  = -sw/cw*Q_dn + 1d0/sw/cw * T3_dn
real(8), public, parameter :: couplZDD_right = -sw/cw*Q_dn

! for check against ttb+photon
! real(8), public, parameter :: couplZUU_left  = Q_up  *1d0
! real(8), public, parameter :: couplZUU_right = Q_up  *(-1d0)
! real(8), public, parameter :: couplZDD_left  = Q_dn  *0d0
! real(8), public, parameter :: couplZDD_right = Q_dn  *0d0

real(8), public, parameter :: couplZTT_left_SM  = -sw/cw*Q_up + 1d0/sw/cw * T3_up!  treat the top separate from up quark
real(8), public, parameter :: couplZTT_right_SM = -sw/cw*Q_up
! real(8), public, parameter :: couplZTT_left_SM  = Q_up!  This is for the check against ttb+photon.
! real(8), public, parameter :: couplZTT_right_SM = Q_up
! real(8), public, parameter :: couplZTT_left_SM  = 1d0!  This is for the check against ttb+photon.
! real(8), public, parameter :: couplZTT_right_SM = 1d0
real(8), public, parameter ::  couplZTT_V_SM =  (+couplZTT_left_SM+couplZTT_right_SM)/2d0   ! our vertex is gamma^mu ( cV - cA*gamma5 ) = gamma^mu *( cL*PL + cR*PR )
real(8), public, parameter ::  couplZTT_A_SM =  (+couplZTT_left_SM-couplZTT_right_SM)/2d0   ! i.e. cV=1/2(cL+cR), cA=1/2(cL-cR)   <-->    cL=cV+cA, cR=cV-cA

real(8), public, parameter :: couplZEE_left  = -sw/cw*Q_el + 1d0/sw/cw * T3_el
real(8), public, parameter :: couplZEE_right = -sw/cw*Q_el

real(8), public, parameter :: couplZNN_left  = -sw/cw*Q_nu + 1d0/sw/cw * T3_nu
real(8), public, parameter :: couplZNN_right = -sw/cw*Q_nu

complex(8), public :: couplZTT_left,couplZTT_right,couplZTT_V,couplZTT_A,couplZTT_left_dyn,couplZTT_right_dyn,couplZTT_V_dyn,couplZTT_A_dyn!  these couplings are set dynamically (depending on whether the Z decays or not)
complex(8), public :: couplZUU_left_dyn,couplZUU_right_dyn
complex(8), public :: couplZDD_left_dyn,couplZDD_right_dyn
complex(8), public :: couplZQQ_left_dyn,couplZQQ_right_dyn

real(8), public :: WidthExpansion
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
integer, public, parameter :: DKX_STChi0_RE1=25
integer, public, parameter :: DKX_STChi0_RE2=26
integer, public, parameter :: DKX_STChi0_RE3=27
integer, public, parameter :: DKX_STChi0_1L1=28
integer, public, parameter :: DKX_STChi0_1L2=29
integer, public, parameter :: DKX_STChi0_1L3=30
integer, public, parameter :: DKX_HTBH_1L1=31
integer, public, parameter :: DKX_HTBH_1L2=32
integer, public, parameter :: DKX_HTBH_1L3=33
integer, public, parameter :: DKX_HTBH_RE1=34
integer, public, parameter :: DKX_HTBH_RE2=35
integer, public, parameter :: DKX_HTBH_RE3=36

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
integer, public, parameter :: ST_Chi0_T_G=24
integer, public, parameter :: HT_BH_T_G=10

real(8), public :: IA0Tt(-1:+1)
real(8), public :: IBHTt(-1:+1)
real(8), public :: IChiStt(-1:+1)

real(8), public :: Collider_Energy
real(8), public :: alpha_s, alpha_s4Pi, alpha_sOver2Pi

real(8) :: ColLO_ttbgg(1:2,1:2)
real(8) :: Col1L_ttbgg(1:2,1:9)
real(8) :: Col1L_ttbggp(1:2,1:3)
real(8) :: ColLO_ttbqqb(1:2,1:1)
real(8) :: Col1L_ttbqqb(1:2,1:4)
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
integer, public :: PrimAmp2_15243,PrimAmp2m_15243,PrimAmp2_12543,PrimAmp2_12453,PrimAmp2m_12435,PrimAmp2_12435,PrimAmp2m_12543,PrimAmp2m_12453,PrimAmp2_13245,PrimAmp2_13254,PrimAmp2_14235,PrimAmp2_14253,PrimAmp2m_13245,PrimAmp2m_13254,PrimAmp2m_14235,PrimAmp2m_14253
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
real(8) :: r2, TopWidthExpansion,WWidthExpansion,WWidthChoice,StopWidthExpansion,HTopWidthExpansion
real(8) :: ZWidth,BrZtoEE
real(8) :: cL,cR,omegasq,f,z
real(8) :: beta,omega,P0,P3,W0,Ppl,Pmi,Wpl,Wmi,Yp,Yw
real(8) :: term4,term7,term9

!--- Zprime section
real(8) :: myCos2thw, cot2thH, g_Zpr, f1, f2, Ga_Zpr_pref
real(8) :: ratio, gamma, rhoplus, rhominus, phi, chi, c0, d0, c1, d1
real(8) :: Ga_Zpr_TOT(0:1), Ga_Zpr_top(0:1), Ga_Zpr_up(0:1), Ga_Zpr_dn(0:1)
real(8), external :: ddilog! this is borrowed from QCDLoop
!--- End Zprime section



m_Bot  = m_Top  ! this is NOT the bottom mass! it is the mass for massive fermion in closed loops
m_SBot = m_STop ! this is NOT the sbottom mass! it is the mass for massive scalar in closed loops
! print *, "set m_SBot to zero";pause


IF( COLLIDER.EQ.1 ) THEN
   Collider_Energy  = 14000d0*GeV

ELSEIF( COLLIDER.EQ.11 ) THEN
   Collider_Energy  = 7000d0*GeV
   COLLIDER = 1

ELSEIF( COLLIDER.EQ.12 ) THEN
   Collider_Energy  = 8000d0*GeV
   COLLIDER = 1

ELSEIF( COLLIDER.EQ.13 ) THEN
   Collider_Energy  = 13000d0*GeV
   COLLIDER = 1
ELSEIF( COLLIDER.EQ.2 ) THEN
  Collider_Energy  = 1960d0*GeV
ENDIF


IF( PDFSet.EQ.2 .AND. (NLOPARAM.EQ.1 .OR. NLOPARAM.EQ.0) ) THEN
  Lambda_QCD = 0.165d0*GeV
  alpha_s = 0.13d0  ! CTEQ6L1
!   alpha_s = 0.118d0   ! CTEQ6L
! print *, "SWITCHED TO WRONG ALPHA_S FOR CHINESE CHECK"

ELSEIF( PDFSet.EQ.2 .AND. (NLOPARAM.EQ.2) ) THEN
  Lambda_QCD = 0.226235d0*GeV
  alpha_s = 0.118d0

ELSEIF(  PDFSet.EQ.1 .AND. (NLOPARAM.EQ.1 .OR. NLOPARAM.EQ.0)) THEN
  Lambda_QCD = 0.167d0/100d0
  alpha_s = 0.13939d0
! alpha_s = 0.13; print *, "change alpha_s for MSTR pdfs (comparison with Kirill)"

ELSEIF( PDFSet.EQ.1 .AND. (NLOPARAM.EQ.2) ) THEN
   Lambda_QCD = 0.226235d0/100d0
   alpha_s=0.12018D0
! alpha_s = 0.119; print *, "change alpha_s for MRST pdfs (comparison with Kirill)"
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

WWidthChoice = 1!          0=experimental W width,    1=calculated W width
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



!   top quark-Z-boson coupling


! our vertex is gamma^mu ( cV - cA*gamma5 ) = gamma^mu *( cL*PL + cR*PR )
! i.e. cV=1/2(cL+cR), cA=1/2(cL-cR)   <-->    cL=cV+cA, cR=cV-cA



! first get (possibly BSM) top-Z coupl
   couplZTT_V = couplZTT_V_SM * (1d0 + DeltaF1V)
   couplZTT_A = couplZTT_A_SM * (1d0 + DeltaF1A)
   couplZTT_left  = (couplZTT_V + couplZTT_A)/1d0     ! MARKUS:  removed 1/2
   couplZTT_right = (couplZTT_V - couplZTT_A)/1d0     ! MARKUS:  removed 1/2 and introduced a minus
   

!   the _dyn variables will be overwritten in mod_ZDecay.f90 if the Z-boson decays
   couplZTT_left_dyn  = couplZTT_left
   couplZTT_right_dyn = couplZTT_right
   
   r2 = (M_Z/(2*4.6d0*GeV))**2! mass correction for bottom quark
   ZWidth = alpha/12d0*M_Z * (  +((couplZUU_left+couplZUU_right)**2 + (couplZUU_left-couplZUU_right)**2 *1d0 )*3d0 &! up
                                +((couplZDD_left+couplZDD_right)**2 + (couplZDD_left-couplZDD_right)**2 *1d0 )*3d0 &! dn
                                +((couplZUU_left+couplZUU_right)**2 + (couplZUU_left-couplZUU_right)**2 *1d0 )*3d0 &! chm
                                +((couplZDD_left+couplZDD_right)**2 + (couplZDD_left-couplZDD_right)**2 *1d0 )*3d0 &! str
                                +((couplZDD_left+couplZDD_right)**2 + (couplZDD_left-couplZDD_right)**2 *(1d0-3d0/2d0/r2) )*3d0 &! bot
                                +((couplZEE_left+couplZEE_right)**2 + (couplZEE_left-couplZEE_right)**2 *1d0 )*3d0 &! el+mu+tau
                                +((couplZNN_left+couplZNN_right)**2 + (couplZNN_left-couplZNN_right)**2 *1d0 )*3d0 &! 3nu                      
                             )
    BrZtoEE = alpha/12d0*M_Z *( (couplZEE_left+couplZEE_right)**2 + (couplZEE_left-couplZEE_right)**2 )/ZWidth
!     print *, "Z->ee branching using calculated LO total width",BrZtoEE
    print *, "Z->ee branching using experimental  total width",alpha/12d0*M_Z *( (couplZEE_left+couplZEE_right)**2 + (couplZEE_left-couplZEE_right)**2 )/Ga_ZExp
!     print *, "sw2",sw**2


!  chiral couplings for stop-Chi^0-top
   cR = 3d0/10d0
   cL = 1d0/10d0
   IChiStt(+1) = cR * m_Top/Vev
   IChiStt(-1) = cL * m_Top/Vev

!  stop-->Chi^0 + top partial width!
   Ga_Stop_ChiTop(0) = SqrtLambda(m_stop**2,m_top**2,m_chi**2)/(16d0*DblPi*m_stop**3) * & 
                  ( (IChiStt(+1)**2+IChiStt(-1)**2)*(m_stop**2-m_top**2-m_chi**2) - 4d0*(IChiStt(+1)*IChiStt(-1))*m_top*m_chi )

IF( Process.ge.51 .and. Process.le.59 ) then
   if( cR.eq.3d0/10d0 .and. cL.eq.1d0/10d0 .and. m_top.eq.172d0*GeV .and. m_chi.eq.50d0*GeV .and. m_stop.eq.350d0*GeV ) then
      Ga_Stop_ChiTop(1) = (-0.00649d0 -0.00172d0) * GeV
      call Error("remove alphas here")
   elseif( cR.eq.3d0/10d0 .and. cL.eq.1d0/10d0 .and. m_top.eq.172d0*GeV .and. m_chi.eq.50d0*GeV .and. m_stop.eq.250d0*GeV ) then
      Ga_Stop_ChiTop(1) =  alpha_s*RunAlphaS(2,MuRen) * (-0.0062804d0 -0.00063978d0) * GeV
   elseif( cR.eq.3d0/10d0 .and. cL.eq.1d0/10d0 .and. m_top.eq.172d0*GeV .and. m_chi.eq.50d0*GeV .and. m_stop.eq.300d0*GeV ) then
      Ga_Stop_ChiTop(1) = (-0.00300d0 - 0.000595d0) * GeV
      call Error("remove alphas here")
   elseif( cR.eq.3d0/10d0 .and. cL.eq.1d0/10d0 .and. m_top.eq.172d0*GeV .and. m_chi.eq.100d0*GeV .and. m_stop.eq.500d0*GeV ) then
      Ga_Stop_ChiTop(1) =  alpha_s*RunAlphaS(2,MuRen) * (-0.179269d0 -0.069456d0) * GeV
   else
      call Error("Ga_Stop_ChiTop(1) needs to be re-calcualted for the given cR,CL,m_Stop,m_top,m_Chi")
!     ./TOPAZ Process=58 Collider=1 TopDK=4 XTopDK=3 NLOParam=2 Correction=4 ObsSet=43 VegasNc0=500000 VegasNc1=500000 MStop=3.00
   endif
ENDIF

   Ga_Stop(:) = Ga_Stop_ChiTop(:)! assuming no other decay channel
   StopWidthExpansion   = -2d0*Ga_Stop(1)/Ga_Stop(0)
IF( abs(XTOPDECAYS).EQ.3 .AND. NLOPARAM.EQ.2 .AND. CORRECTION.EQ.0 ) THEN
   WidthExpansion = WidthExpansion + StopWidthExpansion
ENDIF




IF( XTOPDECAYS.EQ.2 ) THEN
!  chiral couplings for HTop-A0-top
!  IA0Tt(+1) = 1d0/50d0  *1d0
!  IA0Tt(-1) = 1d0/50d0  *1d0
   IA0Tt(+1) = M_HTop/M_A0 * 1d0/50d0  *1d0
   IA0Tt(-1) = M_Top/M_A0  * 1d0/50d0  *1d0

!  HTop-->A0 + top partial width!
   Ga_Htop_A0Top(0) = SqrtLambda(m_Htop**2,m_top**2,m_A0**2)/(16d0*DblPi*m_Htop**3) * & 
                   1d0/2d0*( (IA0Tt(+1)**2+IA0Tt(-1)**2)*(m_Htop**2+m_top**2-m_A0**2) + 4d0*(IA0Tt(+1)*IA0Tt(-1))*m_top*m_HTop )
   Ga_Htop_A0Top(1) = 0d0
   Ga_HTop(:) = Ga_Htop_A0Top(:)! assuming no other decay channel

ELSEIF( XTOPDECAYS.EQ.1 ) THEN! HTop-BH-top


!  MY LO WIDTH for LH and RH couplings
  IBHTt(+1) = 3d0/10d0
  IBHTt(-1) = 1d0/10d0

!  HTop-->BH + top LO partial width!
   Ga_Htop_BHTop(0) = SqrtLambda(m_Htop**2,m_top**2,m_BH**2)/(16d0*DblPi*m_Htop**3) * & 
                  1d0/2d0 *( (IBHTt(+1)**2+IBHTt(-1)**2)*( m_HTop**2+m_Top**2-2d0*m_BH**2+(m_HTop**2-m_top**2)**2/m_BH**2 ) - 12d0*(IBHTt(+1)*IBHTt(-1))*( m_Top*m_HTop ) )

IF( Process.ge.41 .and. Process.le.49 ) then
   if( IBHTt(+1).eq.3d0/10d0 .and. IBHTt(-1).eq.1d0/10d0 .and. m_top.eq.172d0*GeV .and. m_BH.eq.50d0*GeV .and. m_HTop.eq.500d0*GeV ) then
      Ga_Htop_BHTop(1) =  alpha_s*RunAlphaS(2,MuRen) * (  -3.59519d0 -8.99473951d0  ) * GeV
   else
      call Error("Ga_Htop_BHTop(1) needs to be re-calcualted for the given IBHTt(+/-1),m_Htop,m_top,m_BH")
!   ./TOPAZ Process=47 Collider=1 Correction=5 TopDK=1 ObsSet=32 XTopDK=1 MHTop=5.00 MTop=1.72 NLOParam=2 MuRen=5.00 MuFac=5.00 VegasNc0=1000000 VegasNc1=1000000
   endif
ENDIF







! !  MCFM LO and NLO WIDTH for LH currents
! !    IBHTt(+1) = 0d0; print *, "setting RH coupling to zero!!"
! !    IBHTt(-1) = dsqrt(4d0*M_W**2*GF/dsqrt(2d0))
! 
!     beta = m_top/m_Htop
!     omega = m_BH/m_Htop
!     omegasq = omega**2
!     z = m_Top**2/m_Htop**2
!     P0 = 0.5d0*( 1d0 - omegasq + z )
!     P3 = 0.5d0*SqrtLambda(1d0,omegasq,z)
!     f = (1d0-z)**2 + omegasq*(1d0+z) - 2d0*omegasq**2
!     W0 = 0.5d0*( 1d0+omegasq-z )
!     Ppl = P0 + P3
!     Pmi = P0 - P3
!     Wpl = W0 + P3
!     Wmi = W0 - P3
!     Yp = 0.5d0 * dlog(Ppl/Pmi)
!     Yw = 0.5d0 * dlog(Wpl/Wmi)
! 
! 
! ! !  this is the MCFM result (actually Czarnecki)
! !    Ga_Htop_BHTop(0) = 2d0 * P3 * f
! ! print *, "remember: this is only the LH part"
! !    Ga_Htop_BHTop(0) = Ga_Htop_BHTop(0) * m_Htop**3/(8d0*DblPi*dsqrt(2d0)) * ( (IBHTt(-1)**2+IBHTt(+1)**2)/2d0/dsqrt(2d0)/m_BH**2 )! this corresponds to Gamma_0
! 
! 
! !    this is eq.(2.6) of the Campbell+Ellis paper,     seems to be in disagreement with the code below
! !    Ga_Htop_BHTop(1) = RunAlphaS(NLOParam,MuRen)*alpha_sOver2Pi * 4d0/3d0 * (Ga_Htop_BHTop(0)/(2d0*P3*f)) *  &
! !                     ( 8d0*f*P0*( DLi2(1d0-Pmi) - DLi2(1d0-Ppl) - 2d0*DLi2(1d0-Pmi/Ppl) + Yp*dlog(4d0*P3**2/Ppl**2/Wpl) + Yw*dlog(Ppl))   &
! !                      +4d0*(1d0-z)*( (1d0-z)**2 + omegasq*(1d0+z) -4d0*omegasq )*Yw   &
! !                      +(3d0 -z + 11d0*z**2 - z**3 + omegasq*(6d0-12d0*z+2d0*z**2) -omegasq**2*(21d0+5d0*z) + 12d0*omegasq**3  )*Yp   &
! !                      +8d0*f*P3*dlog(omega/4d0/P3**2) + 6d0*(1d0 -4d0*z +3d0*z**2 +omegasq*(3d0+z) -4d0*omegasq**2)*P3*dlog(beta)   &
! !                      +(5d0 -22d0*z + 5d0*z**2 + 9d0*omegasq*(1d0+z) - 6d0*omegasq**2)*P3  &
! !                     )
! 
! !    this is eq.(2.7) of the Campbell+Ellis paper, gives the correct answer for mb-->0
! !    print *, "limit",Ga_Htop_BHTop(0) + RunAlphaS(NLOParam,MuRen)*alpha_sOver2Pi * 4d0/3d0 * (Ga_Htop_BHTop(0)/(2d0*P3*f)) *  &
! !                     ( 4d0*(1d0-omegasq)**2*(1d0+2d0*omegasq)*( DLi2(1d0-omegasq)-DblPi**2/3d0 + 0.5d0*dlog(1d0-omegasq)*dlog(omegasq) )  & 
! !                      -2d0*omegasq*(1d0+omegasq)*(1d0-2d0*omegasq)*dlog(omegasq) - (1d0-omegasq)**2*(4d0*omegasq+5d0)*dlog(1d0-omegasq)  &
! !                      +0.5d0*(1d0-omegasq)*(5d0+9d0*omegasq-6d0*omegasq**2)   &
! !                     )
! 
! !     this code is stolen from MCFM, reproduces correct result for m_Htop-->m_Top, m_BH-->M_W, m_top-->m_b
!       term4=(log(Ppl)-log(beta))*dlog(4d0*P3**2*Wmi/(omegasq*Ppl**2))
!       term7=  &
!       +(3d0-z+11d0*z**2-z**3+omegasq*(6d0-12d0*z+2d0*z**2)   &
!       -omegasq**2*(21d0+5d0*z)+12d0*omegasq**3)*log(Ppl)   &
!       -(-z+11d0*z**2-z**3+omegasq*(-12d0*z+2d0*z**2)   &
!       -omegasq**2*(5d0*z))*log(beta)
!       term9=   &
!       +6d0*(1d0-4d0*z+3d0*z**2+omegasq*(3d0+z)-4d0*omegasq**2)   &
!       *(P3-0.5d0*(1d0-omegasq))*dlog(beta)   &
!       +3d0*(1d0-omegasq)*(-4d0*z+3d0*z**2+omegasq*(z))*dlog(beta) 
!  
!       Ga_Htop_BHTop(1) = RunAlphaS(NLOParam,MuRen)*alpha_sOver2Pi * 4d0/3d0 * (Ga_Htop_BHTop(0)/(2d0*P3*f)) *  &
!                         ( 8d0*f*P0*(DLi2(1d0-Pmi)-DLi2(1d0-Ppl)                &
!                             -2d0*DLi2(1d0-Pmi/Ppl)+term4    &
!                             +Yw*dlog(Ppl))   &
!                             +4d0*(1d0-z)*((1d0-z)**2+omegasq*(1d0+z)-4d0*omegasq**2)*Yw   &
!                             +term7    &
!                             +8d0*f*P3*dlog(omega/4d0/P3**2)+term9    &
!                             +(5d0-22d0*z+5d0*z**2+9d0*omegasq*(1d0+z)-6d0*omegasq**2)*P3  & 
!                         )
!       if( IBHTt(-1)*IBHTt(+1).ne.0d0 ) call Error("MCFM formula for Ga_Htop_BHTop(1) cannot be used")! only pure LH or RH currents are allowed



   Ga_HTop(:) = Ga_Htop_BHTop(:)! assuming no other decay channel

   HTopWidthExpansion   = -2d0*Ga_HTop(1)/Ga_HTop(0)
IF( abs(XTOPDECAYS).EQ.1 .AND. NLOPARAM.EQ.2 .AND. CORRECTION.EQ.0 ) THEN
   WidthExpansion = WidthExpansion + HTopWidthExpansion
ENDIF


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

!-- Sequential couplings
!  Collider=2 TopDK=0 ObsSet=67 Process=62 NLOParam=1 Correction=0  MZpr=0.9119 GaZpr=0.02 MuRen=1.73 MuFac=1.73
! SM couplings
! gR_Zpr(top_) =  eL*(sw/cw)*2d0/3d0 
! gL_Zpr(top_) = -eL/(sw*cw)*( 0.5d0  - sw**2*2d0/3d0 )
! gL_Zpr(dn_) = -eL/(sw*cw)*(-0.5d0 + sw**2/3d0)
! gR_Zpr(dn_)  = -eL*(sw/cw)/3d0
! gL_Zpr(up_) = gL_Zpr(top_)
! gR_Zpr(up_) = gR_Zpr(top_)
! gL_Zpr(chm_) = gL_Zpr(top_)
! gR_Zpr(chm_) = gR_Zpr(top_)
! gL_Zpr(str_) = gL_Zpr(dn_)
! gR_Zpr(str_) = gR_Zpr(dn_)
! gL_Zpr(bot_) = gL_Zpr(dn_)
! gR_Zpr(bot_) = gR_Zpr(dn_)


!print *, 'gR_Zpr(top_)', gR_Zpr(top_)
!print *, 'gL_Zpr(top_)', gL_Zpr(top_)
!print *, 'gR_Zpr(dn_)', gR_Zpr(dn_)
!print *, 'gL_Zpr(dn_)', gL_Zpr(dn_)

!-- Leptophobic top-color Z', see 1112.4928
myCos2thw = 0.768d0

!cot2thH = 8d0 * myCos2thw * Ga_Zpr/(alpha*m_Zpr) / &
!     (dsqrt(1d0-4d0*m_top**2/m_Zpr**2)*(2d0+4d0*m_top**2/m_Zpr**2)+4d0)
!
!Ga_Zpr_pref = alpha * M_Zpr * cot2thH / myCos2thw / 4d0
!print *, 'cot(th_H)**2', cot2thH
!g_Zpr = 1d0/2d0 * dsqrt(4d0*pi*alpha/myCos2thw * cot2thH)
!g_Zpr = dsqrt(Ga_Zpr/M_Zpr * 8d0 * Pi / (dsqrt(1d0-4d0*m_Top**2/m_Zpr**2)*(2d0+4d0*m_top**2/m_Zpr**2)+4d0))
! gL_Zpr(:) = 0d0
! gR_Zpr(:) = 0d0

!-- original experimental setup: f1 = 1, f2 = 0
!-- topcolor tilting motivated parameters: f1 > 0 and / or f2 < 0
f1 = 1d0
f2 = 0d0

gL_Zpr(up_) = -1d0
gR_Zpr(up_) = -f1 * 1d0

gL_Zpr(dn_) = -1d0
gR_Zpr(dn_) = -f2 * 1d0

gL_Zpr(top_) = 1d0
gR_Zpr(top_) = f1 * 1d0

gL_Zpr(bot_) = 1d0
gR_Zpr(bot_) = f2 * 1d0

ratio = M_Zpr**2/(4d0*M_Top**2)
gamma = dlog(2d0*dsqrt(ratio))
rhoplus = dsqrt(ratio)+dsqrt(ratio-1d0)
rhominus = dsqrt(ratio)-dsqrt(ratio-1d0)
chi = dlog(rhoplus-rhominus)
phi = dlog(rhoplus)

!-- results for the NLO corrections for a generic coupled Z' taken from
!-- Nucl.Phys. B329 (1990) 547 
!-- Overall prefactor for LO width taken from hep-ph/9911288v1,
!-- which is WRONG and must be corrected the following way:
!-- Eq. 28: -3 f1 -> +6 f1
!-- Eq. 32: +3 f1 -> -6 f1

!-- LO decay width, c0 is vector d0 is axial
 c0 = dsqrt(dabs(1d0-1d0/ratio)) * (1d0+1d0/(2d0*ratio))
 d0 = dabs(1d0-1d0/ratio)**(3d0/2d0)

!-- NLO decay, coefficients of asonpi
 c1 = 8d0/3d0 * (1d0-1d0/(4d0*ratio**2)) * ( 2d0 * ddilog(rhominus**2) &
     + ddilog(rhominus**4) + 2d0 * phi * (3d0 * phi - gamma-2d0 * chi)) &
     - 8d0/3d0 * dsqrt(1d0-1d0/ratio) * (1d0+1d0/(2d0 * ratio)) * (gamma + 2d0 * chi) &
     + 8d0 * (1d0-1d0/(6d0*ratio) - 7d0/(48d0*ratio**2))*phi + dsqrt(1d0-1d0/ratio) * (1d0+3d0/(2d0 * ratio))

 d1 = 8d0/3d0 * (1d0-3d0/(2d0 * ratio) + 1d0/(2d0 * ratio**2)) * ( 2d0 * ddilog(rhominus**2) &
      + ddilog(rhominus**4) + 2d0 * phi * (3d0 * phi-gamma-2d0 * chi)) - 8d0/3d0 * &
     (1d0-1d0/ratio)**(3d0/2d0) * (gamma+2d0 * chi) + 8d0 * (1d0-11d0/(12d0 * ratio) + 5d0/(48d0 * ratio**2) &
     + 1d0/(32d0 * ratio**3)) * phi + dsqrt(1d0-1d0/ratio) * (1d0-3d0/ratio+1d0/(4d0 * ratio**2))

 Ga_Zpr_top(0) = ( c0 * ((1d0 + f1)/2d0)**2 + d0 * ((1d0-f1)/2d0)**2 )
 Ga_Zpr_up(0) =  (((1d0+f1)/2d0)**2 + ((1d0-f1)/2d0)**2 )
 Ga_Zpr_dn(0) =  (((1d0+f2)/2d0)**2 + ((1d0-f2)/2d0)**2 )

 Ga_Zpr_top(1) = ( c1 * ((1d0 + f1)/2d0)**2 + d1 * ((1d0-f1)/2d0)**2 ) * alpha_s * RUNALPHAS(2,MuRen) / Pi
 Ga_Zpr_up(1) =  (((1d0+f1)/2d0)**2 + ((1d0-f1)/2d0)**2 ) * alpha_s * RUNALPHAS(2,MuRen) / Pi
 Ga_Zpr_dn(1) =  (((1d0+f2)/2d0)**2 + ((1d0-f2)/2d0)**2 ) * alpha_s * RUNALPHAS(2,MuRen) / Pi

 Ga_Zpr_TOT = Ga_Zpr_top + Ga_Zpr_up + 2d0 * Ga_Zpr_dn

!-- use the Ga_Zpr(0) as input for couplings
 Ga_Zpr_pref = Ga_Zpr   /Ga_Zpr_TOT(0)

!if (NLOParam.le.1) then
!   Ga_Zpr_pref = Ga_Zpr/Ga_Zpr_TOT(0)
!else
!   Ga_Zpr_pref = Ga_Zpr/(Ga_Zpr_TOT(0) + Ga_Zpr_TOT(1))
!endif

 Ga_Zpr_TOT = Ga_Zpr_TOT * Ga_Zpr_pref

if (NLOParam.le.1) then
   Ga_Zpr = Ga_Zpr_TOT(0)
else
   Ga_Zpr = Ga_Zpr_TOT(0) + Ga_Zpr_TOT(1)
endif

!-- properly normalize couplings
g_Zpr = dsqrt(4d0*pi/M_Zpr * Ga_Zpr_pref)



!  remove this rescaling if SM Z is assumed
gR_Zpr(:) = gR_Zpr(:) * g_Zpr
gL_Zpr(:) = gL_Zpr(:) * g_Zpr




!print *, 'g_Zpr', g_Zpr
!print *, 'Ga_Zpr_TOT(0)', Ga_Zpr_TOT(0)
!print *, 'Ga_Zpr_TOT(0)+Ga_Zpr_Tot(0)+Ga_Zpr_tot(1)', Ga_Zpr_TOT(0)+Ga_Zpr_Tot(1)
!print *, 'alpha_s*RUNALPHAS(2,MuRen)', alpha_s*RUNALPHAS(2,MuRen)
!print *, 'c0', c0
!print *, 'c1', c1
!print *, 'Ga_Zpr', Ga_Zpr
!stop

!!! End Zprime section !!!


END SUBROUTINE



FUNCTION GetMass( PartType )
use ModMisc
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
   elseif( abs( PartType) .eq. SBot_) then
      GetMass = m_SBot
   elseif( PartType .eq. 0 ) then
      GetMass = 0d0
   else
      call Error("GetMass particle type not available.")
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
! ! !   DUW version
!    if( Loop.eq.1 ) then
!      RunAlphaS = 4d0*dblpi/beta0/dlog((Q*100d0)**2/165d-3**2)
!      RunAlphaS = RunAlphaS/alpha_s
!    elseif( Loop.eq.2 ) then
!      RunAlphaS = 4d0*dblpi/beta0/dlog((Q*100d0)**2/226d-3**2)*(1d0-(beta1*2d0)/beta0**2*dlog(dlog((Q*100d0)**2/226d-3**2))/dlog((Q*100d0)**2/226d-3**2))
!      RunAlphaS = RunAlphaS/alpha_s
!    else
!      RunAlphaS = 1d0   ! no running
!    endif



RETURN
END FUNCTION


END MODULE



