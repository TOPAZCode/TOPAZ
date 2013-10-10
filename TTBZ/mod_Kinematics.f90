MODULE ModKinematics
implicit none
save


type :: Histogram
    integer :: NBins
    real(8) :: BinSize
    real(8) :: LowVal
    real(8) :: SetScale
    real(8),allocatable :: Value(:)
    real(8),allocatable :: Value2(:)
    integer,allocatable :: Hits(:)
    character :: Info*(50)
end type

!type :: Histogram2D
!    integer :: NBins(1:2)
!    real(8) :: BinSize(1:2)
!    real(8) :: LowVal(1:2)
!    real(8) :: SetScale(1:2)
!    real(8),allocatable :: Value(:,:)
!    real(8),allocatable :: Value2(:,:)
!    integer,allocatable :: Hits(:,:)
!    character :: Info*(50)
!end type


integer,public :: it_sav

type(Histogram),allocatable   :: Histo(:)
!type(Histogram2D),allocatable :: Histo2D(:)

real(8) :: pT_jet_cut, pT_bjet_cut, pT_lep_cut, Rsep_jet, Rsep_LepJet, pT_miss_cut, eta_sepa_cut, MInv_jets_cut, eta_lep_cut, eta_jet_cut, eta_bjet_cut, HT_cut, pT_hardestjet_cut
real(8) :: pT_pho_cut,Rsep_Pj,Rsep_Pbj,Rsep_Plep,eta_pho_cut,MTW_cut, Mttbar_cut,Rsep_jetlep

real(8),public ::MInv_LB

contains




SUBROUTINE InitPSCuts()
use ModMisc
use ModParameters
implicit none


pT_jet_cut        = 1d100
pT_bjet_cut       = 1d100
pT_lep_cut        = 1d100
Rsep_jet          = 1d100
Rsep_LepJet       = 1d100
pT_miss_cut       = 1d100
eta_sepa_cut      = 1d100
MInv_jets_cut     = 1d100
eta_lep_cut       = 1d100
eta_jet_cut       = 1d100
eta_bjet_cut      = 1d100
HT_cut            = 1d100
pT_hardestjet_cut = 1d100
pT_pho_cut        = 1d100
Rsep_Pj           = 1d100
Rsep_Pbj          = 1d100
Rsep_Plep         = 1d100
eta_pho_cut       = 1d100


IF( ObsSet.EQ.0 ) THEN! set of observables for ttb production without decays at Tevatron

ELSEIF( ObsSet.EQ.1 ) THEN! set of observables for ttb production without decays at LHC

ELSEIF( ObsSet.EQ.2 ) THEN! set of observables for ttb production as signal process at Tevatron (di-lept. decay)
    Rsep_jet    = 0.4d0         !*0d0
    pT_bjet_cut = 20d0*GeV      !*0d0
    eta_bjet_cut= 2.5d0         !*1d2
    pT_lep_cut  = 20d0*GeV      !*0d0
    pT_miss_cut = 25d0*GeV      !*0d0
    eta_lep_cut = 2.5d0         !*1d2

ELSEIF( ObsSet.EQ.3 ) THEN! set of observables for ttb production as signal process at LHC (di-lept. decay)
    Rsep_jet    = 0.4d0         !*0d0
    pT_bjet_cut = 25d0*GeV      !*0d0
    eta_bjet_cut= 2.5d0         !*1d2
    pT_lep_cut  = 25d0*GeV      !*0d0
    pT_miss_cut = 50d0*GeV      !*0d0
    eta_lep_cut = 2.5d0         !*1d2

ELSEIF( ObsSet.EQ.4 ) THEN! ! set of observables for ttb production with hadr. Atop, lept. top decay
    Rsep_jet    = 0.5d0
    pT_bjet_cut = 20d0*GeV
    eta_bjet_cut= 2.0d0
    pT_jet_cut  = 20d0*GeV
    eta_jet_cut = 2.0d0
    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 2.0d0
    pT_miss_cut = 20d0*GeV
    HT_cut      = 220d0*GeV

ELSEIF( ObsSet.EQ.5 ) THEN! set of observables for ttb production with hadr. top, lept. Atop decay at LHC

!   these are the cuts for muons
    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 2.5d0

    pT_bjet_cut = 25d0*GeV
    pT_jet_cut  = 25d0*GeV
    eta_bjet_cut= 2.5d0
    eta_jet_cut = 2.5d0

    pT_miss_cut = 25d0*GeV
!   Mwt cut is hard coded below

    Rsep_LepJet = 0.4d0
    Rsep_jet    = 0.4d0




ELSEIF( ObsSet.EQ.6 ) THEN! set of observables for ttb production with lept. top, hadr. Atop decay at LHC
    Rsep_jet    = 0.4d0         *0d0
    pT_bjet_cut = 25d0*GeV      *0d0
    eta_bjet_cut= 2.5d0         *1d2
    pT_jet_cut  = 25d0*GeV      *0d0
    eta_jet_cut = 2.5d0         *1d2
    pT_lep_cut  = 25d0*GeV      *0d0
    eta_lep_cut = 2.5d0         *1d2
    pT_miss_cut = 25d0*GeV      *0d0


ELSEIF( ObsSet.EQ.7 ) THEN! set of observables for ttb production with lept. top and J/Psi fragmentation, hadr. Atop decay at LHC
    Rsep_jet    = 0.5d0
    HT_cut      = 100d0*GeV
    pT_jet_cut  = 20d0*GeV
    pT_lep_cut  = 20d0*GeV

ELSEIF( ObsSet.EQ.8 ) THEN! set of observables for ttb spin correlations at LHC (di-lept. decay)
    pT_bjet_cut = 25d0*GeV
    pT_lep_cut  = 20d0*GeV
    pT_miss_cut = 40d0*GeV
    eta_lep_cut = 2.5d0
    Rsep_jet    = 0.4d0

ELSEIF( ObsSet.EQ.9 ) THEN! this is for the factorization check


ELSEIF( ObsSet.EQ.10 ) THEN! set of observables for ttbjet production without decays at Tevatron
    pT_jet_cut  = 20d0*GeV
    Rsep_jet    = 1d0

ELSEIF( ObsSet.EQ.11 ) THEN! set of observables for ttbjet production without decays at LHC
    pT_jet_cut  = 50d0*GeV
    Rsep_jet    = 0.4d0

ELSEIF( ObsSet.EQ.12 ) THEN! set of observables for ttbjet production as signal process at Tevatron (hadr.Atop, lept.top decay)
    pT_jet_cut  = 20d0*GeV
    pT_bjet_cut = pT_jet_cut
    eta_jet_cut = 2d0            *100d0
    eta_bjet_cut= eta_jet_cut    *100d0
    pT_lep_cut  = 20d0*GeV       *0d0
    eta_lep_cut = 1d0            *100d0
    pT_miss_cut = 20d0*GeV       *0d0
    HT_cut      = 220d0*GeV      *0d0
    Rsep_jet    = 0.5d0

ELSEIF( ObsSet.EQ.13 ) THEN! set of observables for ttbjet production as signal process at LHC (di-lept. decay)
    pT_jet_cut  = 25d0*GeV
    pT_bjet_cut = 25d0*GeV
    eta_jet_cut = 2.5d0
    eta_bjet_cut= 2.5d0
    pT_lep_cut  = 25d0*GeV
    eta_lep_cut = 2.5d0
    pT_miss_cut = 50d0*GeV
    Rsep_jet    = 0.4d0

ELSEIF( ObsSet.EQ.14 ) THEN! set of observables for ttbjet production as background process to VBF at LHC (di-lept. decay)
    pT_hardestjet_cut = 40d0*GeV
    pT_jet_cut   = 20d0*GeV
    pT_bjet_cut = pT_jet_cut
    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 2d0
    eta_sepa_cut = 3.0d0  ! 2.5d0 !  3.8d0
    MInv_jets_cut= 550d0*GeV
!     eta_jet_cut  = 2.5d0  ! REMOVED!!  this value is used in the opposite way, i.e. if(eta_jet .lt. eta_jet_cut) then cut
    eta_bjet_cut = eta_jet_cut
    Rsep_jet     = 0.5d0

ELSEIF( ObsSet.EQ.15 ) THEN! set of observables for ttbjet production as signal process at LHC (hadr.Atop, lept.top decay)
    pT_jet_cut  = 25d0*GeV
    pT_bjet_cut = pT_jet_cut
    eta_jet_cut = 2.5d0
    eta_bjet_cut= eta_jet_cut
    pT_lep_cut  = 25d0*GeV
    eta_lep_cut = 2.5d0
    pT_miss_cut = 30d0*GeV
    Rsep_jet    = 0.4d0
!   added an additional cut on MT(W)(+ET_miss) for ATLAS analysis, see kinematics_ttbjet subroutine


ELSEIF( ObsSet.EQ.19 ) THEN! for checks of ttbjet
    pT_jet_cut  = 20d0*GeV
    Rsep_jet    = 0.4d0


ELSEIF( ObsSet.EQ.20 ) THEN! set of observables for ttbgamma production without decays at Tevatron
!     Rsep_jet    = 1d0
    pT_pho_cut  = 20d0*GeV
    Rsep_Pj      = 0.4d0
!     Rsep_Plep    = 0.4d0


ELSEIF( ObsSet.EQ.21 ) THEN! set of observables for ttbgamma production without decays at LHC
!     Rsep_jet    = 1d0
    pT_pho_cut  = 20d0*GeV
    Rsep_Pj      = 0.4d0
!     Rsep_Plep    = 0.4d0


ELSEIF( ObsSet.EQ.22 ) THEN! set of observables for ttbgamma production with di-lept.decays at Tevatron
    Rsep_jet    = 0.4d0
    pT_pho_cut  = 10d0*GeV
    Rsep_Pj     = 0.4d0
    Rsep_Pbj    = 0.4d0
    Rsep_Plep   = 0.4d0

    pT_bjet_cut = 25d0*GeV
    eta_bjet_cut= 2.5d0
    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 2.5d0
    pT_miss_cut = 40d0*GeV

ELSEIF( ObsSet.EQ.23 ) THEN! set of observables for ttbgamma production with di-lept.decays at LHC
    Rsep_jet    = 0.4d0
    pT_pho_cut  = 20d0*GeV
    Rsep_Pj     = 0.4d0
    Rsep_Pbj    = 0.4d0
    Rsep_Plep   = 0.4d0

    pT_bjet_cut = 25d0*GeV
    eta_bjet_cut= 2.5d0
    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 2.5d0
    pT_miss_cut = 40d0*GeV

ELSEIF( ObsSet.EQ.24 ) THEN! set of observables for ttbgamma production with semi-lept.decays(hadr.Atop, lept.top decay) at Tevatron

    pT_pho_cut  = 10d0*GeV
    eta_pho_cut = 1.1d0
    Rsep_Plep   = 0.4d0
    Rsep_Pj     = 0.4d0
    Rsep_Pbj    = 0.4d0

    Rsep_jet    = 0.4d0
    pT_bjet_cut = 15d0*GeV
    pT_jet_cut  = 15d0*GeV
    eta_bjet_cut= 2d0
    eta_jet_cut = 2d0

    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 1.1d0

    pT_miss_cut = 20d0*GeV
    HT_cut      = 200d0*GeV



ELSEIF( ObsSet.EQ.25 ) THEN! set of observables for ttbgamma production with semi-lept.decays(hadr.Atop, lept.top decay) at LHC

!   these are the cuts for muons
    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 2.5d0

    pT_bjet_cut = 25d0*GeV
    pT_jet_cut  = 25d0*GeV
    eta_bjet_cut= 2.5d0
    eta_jet_cut = 2.5d0

    pT_pho_cut  = 15d0*GeV
    eta_pho_cut = 2.37d0
!   cracks for photon are hard coded below

    pT_miss_cut = 25d0*GeV
!   Mwt cut is hard coded below

    pT_pho_cut  = 15d0*GeV
    eta_pho_cut = 2.37d0
!   cracks for photon are hard coded below

    Rsep_LepJet = 0.4d0
    Rsep_jet    = 0.4d0
    Rsep_Pj     = 0.5d0
    Rsep_Pbj    = 0.5d0
    Rsep_Plep   = 0.4d0!  not specified




!   below is a copy of the ObsSet=28 case
!     pT_pho_cut  = 20d0*GeV
!     eta_pho_cut = 2.5d0
!     Rsep_Plep   = 0.4d0
!     Rsep_Pj     = 0.4d0
!     Rsep_Pbj    = 0.4d0
! 
!     Rsep_jet    = 0.4d0
!     pT_bjet_cut = 20d0*GeV
!     pT_jet_cut  = 20d0*GeV
!     eta_bjet_cut= 2d0
!     eta_jet_cut = 2.5d0
! 
!     pT_lep_cut  = 20d0*GeV
!     eta_lep_cut = 2.5d0
! 
!     pT_miss_cut = 20d0*GeV
!     HT_cut      = 200d0*GeV



ELSEIF( ObsSet.EQ.26 ) THEN! set of observables for ttbgamma production with semi-lept.decays(hadr.Atop, lept.top decay) at Tevatron with Qtop cuts

    pT_pho_cut  = 20d0*GeV
    eta_pho_cut = 2.5d0
    Rsep_Plep   = 0.4d0
    Rsep_Pj     = 0.4d0
    Rsep_Pbj    = 0.4d0

    Rsep_jet    = 0.4d0
    pT_bjet_cut = 20d0*GeV
    pT_jet_cut  = 20d0*GeV
    eta_bjet_cut= 2d0
    eta_jet_cut = 2.5d0

    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 2.5d0

    pT_miss_cut = 20d0*GeV
    HT_cut      = 200d0*GeV

ELSEIF( ObsSet.EQ.27 ) THEN! set of observables for ttbgamma production with semi-lept.decays(hadr.Atop, lept.top decay) at LHC with Qtop cuts

    pT_pho_cut  = 20d0*GeV
    eta_pho_cut = 2.5d0
    Rsep_Plep   = 0.4d0
    Rsep_Pj     = 0.4d0
    Rsep_Pbj    = 0.4d0

    Rsep_jet    = 0.4d0
    pT_bjet_cut = 20d0*GeV
    pT_jet_cut  = 20d0*GeV
    eta_bjet_cut= 2d0
    eta_jet_cut = 2.5d0

    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 2.5d0

    pT_miss_cut = 20d0*GeV
    HT_cut      = 200d0*GeV

ELSEIF( ObsSet.EQ.28 ) THEN! set of observables for ttbgamma production with semi-lept.decays(hadr.Atop, lept.top decay) at LHC

    pT_pho_cut  = 20d0*GeV
    eta_pho_cut = 2.5d0
    Rsep_Plep   = 0.4d0
    Rsep_Pj     = 0.4d0
    Rsep_Pbj    = 0.4d0

    Rsep_jet    = 0.4d0
    pT_bjet_cut = 20d0*GeV
    pT_jet_cut  = 20d0*GeV
    eta_bjet_cut= 2d0
    eta_jet_cut = 2.5d0

    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 2.5d0

    pT_miss_cut = 20d0*GeV
    HT_cut      = 200d0*GeV

ELSEIF( ObsSet.EQ.29 ) THEN! this is for the factorization check

    pT_pho_cut  = 20d0*GeV
    Rsep_Pj     = 0.4d0




ELSEIF( ObsSet.EQ.31 ) THEN! set of observables for HTHTbar + A0/BH production (stable)


ELSEIF( ObsSet.EQ.32 ) THEN! set of observables for HTHTbar + A0/BH production (di-lept. tops)
    Rsep_jet    = 0.4d0

    pT_bjet_cut = 30d0*GeV 
    eta_bjet_cut= 2.5d0 

    pT_lep_cut  = 20d0*GeV  
    eta_lep_cut = 2.5d0 
    pT_miss_cut = 25d0*GeV! note that this is ET and not pT


ELSEIF( ObsSet.EQ.33 ) THEN! set of observables for HTHTbar + A0/BH production (semi-hadr. tops)
    Rsep_jet    = 0.4d0

    pT_bjet_cut = 30d0*GeV
    eta_bjet_cut= 2.5d0
    pT_jet_cut = 30d0*GeV 
    eta_jet_cut= 2.5d0     
 
    pT_lep_cut  = 20d0*GeV  
    eta_lep_cut = 2.5d0      
    pT_miss_cut = 150d0*GeV! note that this is ET and not pT

    MTW_cut = 120d0*GeV



ELSEIF( ObsSet.EQ.34 ) THEN! set of observables for HTHTbar + A0/BH production (di-lept. tops) without acceptance cuts
    Rsep_jet    = 0d0

    pT_bjet_cut = 0d0*GeV 
    eta_bjet_cut= 100d0 

    pT_lep_cut  = 0d0*GeV  
    eta_lep_cut = 100d0 
    pT_miss_cut = 0d0*GeV! note that this is ET and not pT


ELSEIF( ObsSet.EQ.35 ) THEN! set of observables for HTHTbar + A0/BH production (semi-hadr. tops) without acceptance cuts
    Rsep_jet    = 0d0

    pT_bjet_cut = 0d0*GeV
    eta_bjet_cut= 100d0
    pT_jet_cut =  0d0*GeV 
    eta_jet_cut=  100d0     
 
    pT_lep_cut  = 0d0*GeV  
    eta_lep_cut = 100d0      
    pT_miss_cut = 0d0*GeV! note that this is ET and not pT

    MTW_cut = 150d0*GeV  * 0d0


ELSEIF( ObsSet.EQ.41 ) THEN! set of observables for STSTbar + Chi production (stable)



ELSEIF( ObsSet.EQ.42 ) THEN! set of observables for STSTbar + Chi production (di-lept. tops)
    Rsep_jet    = 0.4d0 
    pT_bjet_cut = 25d0*GeV 
    eta_bjet_cut= 2.5d0     
    pT_lep_cut  = 20d0*GeV       
    eta_lep_cut = 2.5d0            
    pT_miss_cut = 80d0*GeV      ! note that this is ET and not pT


ELSEIF( ObsSet.EQ.43 ) THEN! set of observables for STSTbar + Chi production (semi-hadr. tops)

! !   these are the cuts for mstop/chi = 500/100 GeV analysis at 8TeV
!  if( Collider_Energy.eq.8000d0*GeV .or. Collider_Energy.eq.14000d0*GeV ) then

    Rsep_jet    = 0.4d0
    pT_bjet_cut = 30d0*GeV
    eta_bjet_cut= 2.5d0
    pT_jet_cut = 30d0*GeV
    eta_jet_cut= 2.5d0
    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 2.5d0
    pT_miss_cut = 150d0*GeV! note that this is ET and not pT
    MTW_cut = 120d0*GeV

!  elseif( Collider_Energy.eq.7000d0*GeV ) then
! !   these are the cuts for mstop/chi = 300/100 GeV analysis at 7TeV
!     Rsep_jet    = 0.4d0
!     pT_bjet_cut = 30d0*GeV
!     eta_bjet_cut= 2.5d0
!     pT_jet_cut = 25d0*GeV
!     eta_jet_cut= 2.5d0
!     pT_lep_cut  = 20d0*GeV
!     eta_lep_cut = 2.5d0
!     pT_miss_cut = 125d0*GeV! note that this is ET and not pT
!     MTW_cut = 0d0*GeV
!  else
!     call Error("This ObsSet only supports Collider=1,11,12")
!  endif


ELSEIF( ObsSet.EQ.44 ) THEN! set of observables for STSTbar + Chi production (di-lept. tops) without acceptance cuts

    Rsep_jet    = 0d0
    pT_bjet_cut = 0d0*GeV
    eta_bjet_cut= 100d0
    pT_lep_cut  = 0d0*GeV
    eta_lep_cut = 100d0
    pT_miss_cut = 0d0*GeV! note that this is ET and not pT


ELSEIF( ObsSet.EQ.45 ) THEN! set of observables for STSTbar + Chi production (semi-hadr. tops) without acceptance cuts

    Rsep_jet    = 0d0
    pT_bjet_cut = 0d0*GeV
    eta_bjet_cut= 100d0
    pT_jet_cut = 0d0*GeV
    eta_jet_cut= 100d0
    pT_lep_cut  = 0d0*GeV
    eta_lep_cut = 100d0
    pT_miss_cut = 0d0*GeV! note that this is ET and not pT

    MTW_cut = 150d0*GeV  * 0d0


ELSEIF( ObsSet.EQ.48 ) THEN! set of observables for STOP width



ELSEIF( ObsSet.EQ.51 ) THEN! set of observables for ttb+Z (stable tops)

ELSEIF( ObsSet.EQ.52 ) THEN! set of observables for ttb+Z ( di-lept. ttbar decays and di-lept. Z decay )


    Rsep_jet    = 0.4d0 
    pT_bjet_cut = 25d0*GeV
    eta_bjet_cut= 2.5d0
    pT_lep_cut  = 25d0*GeV
    pT_miss_cut = 50d0*GeV
    eta_lep_cut = 2.5d0 



ELSEIF( ObsSet.EQ.53 ) THEN! set of observables for ttb+Z ( semi-lept. ttbar decays and di-lept. Z decay )


    Rsep_jet    = 0.4d0
    pT_bjet_cut = 20d0*GeV
    eta_bjet_cut= 2.5d0
    pT_jet_cut  = 20d0*GeV
    eta_jet_cut = 2.5d0

    pT_lep_cut  = 15d0*GeV
    pT_miss_cut = 20d0*GeV
    eta_lep_cut = 2.5d0
    Rsep_jetlep = 0.4d0


ELSEIF( ObsSet.EQ.55 ) THEN! set of observables for ttb+Z ( di-lept. ttbar decays and di-lept. Z decay ) same as 52 but no cuts


    Rsep_jet    = 0.4d0         *0d0
    pT_bjet_cut = 25d0*GeV      *0d0
    eta_bjet_cut= 2.5d0         *1d2
    pT_lep_cut  = 25d0*GeV      *0d0
    pT_miss_cut = 50d0*GeV      *0d0
    eta_lep_cut = 2.5d0         *1d2


ELSEIF( ObsSet.EQ.56 ) THEN! set of observables for ttb+Z ( semi-lept. ttbar decays and di-lept. Z decay ) same as 53 but no cuts

    Rsep_jet    = 0.4d0         *0d0
    pT_bjet_cut = 20d0*GeV      *0d0
    eta_bjet_cut= 2.5d0         *1d2
    pT_jet_cut  = 20d0*GeV      *0d0
    eta_jet_cut = 2.5d0         *1d2

    pT_lep_cut  = 15d0*GeV      *0d0 
    pT_miss_cut = 20d0*GeV      *0d0 
    eta_lep_cut = 2.5d0         *1d2
    Rsep_jetlep = 0.4d0         *0d0


ELSEIF ( ObsSet.EQ.60 ) THEN ! Zprime, stable top

   Mttbar_cut = 500d0*GeV

ELSEIF ( ObsSet.EQ.61 ) THEN ! Zprime, top decay to dileptons

   Rsep_jet = 0.5d0             !*0d0 !this removes all cuts for fact.checks

   pT_lep_cut  = 20d0*GeV       !*0d0
   eta_lep_cut = 2.5d0          !*1d6

   pT_bjet_cut = 30d0*GeV       !*0d0
   eta_bjet_cut = 2.5d0         !*1d6
   
ELSEIF ( ObsSet.EQ.62 ) THEN ! Zprime, fully hadronic top decay



ELSEIF ( ObsSet.EQ.64 ) THEN ! Zprime, semi-hadronic top decay (for factorization checks)
   
ELSEIF ( ObsSet.EQ.65 ) THEN ! Zprime, semi-hadronic top decay (for ATLAS analysis: James Ferrando)
   
!   this is for electrons and muons
    pT_lep_cut  = 25d0*GeV
    eta_lep_cut = 2.47d0

    Rsep_LepJet = 0.4d0


    pT_miss_cut = 30d0*GeV! note that this is ET and not pT
    MTW_cut = 30d0*GeV

    pT_jet_cut = 300d0*GeV
    eta_jet_cut= 2.0d0
!   more jet cuts are defined inside KinematicsZprimeTTB subroutine   




ELSEIF ( ObsSet.EQ.66 ) THEN ! Zprime, semi-hadronic top decay (for CMS analysis: Roman Kogler)

!   this is for electrons
    pT_lep_cut  = 35d0*GeV
    eta_lep_cut = 2.5d0

!   this is for muons
!    pT_lep_cut  = 45d0*GeV
!    eta_lep_cut = 2.1d0

    Rsep_LepJet = 0.5d0
!   pTrel is defined inside KinematicsZprimeTTB subroutine

    pT_miss_cut = 50d0*GeV! note that this is ET and not pT
    HT_cut = 150d0*GeV

    Rsep_jet = 0.5d0
    pT_jet_cut = 150d0*GeV
    eta_jet_cut= 2.4d0
!   more jet cuts are defined inside KinematicsZprimeTTB subroutine   

   
ENDIF

END SUBROUTINE





SUBROUTINE InitHisto()
use ModMisc
use ModParameters
implicit none
integer :: AllocStatus,NHisto

it_sav = 1


IF( ObsSet.EQ.0 ) THEN! set of observables for ttb production without decays at Tevatron
          if(Collider.ne.2)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.0) call Error("TopDecays needs to be 0!")
          NumHistograms = 6
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 50d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "eta_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

ELSEIF( ObsSet.EQ.1 ) THEN! set of observables for ttb production without decays at LHC
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.0) call Error("TopDecays needs to be 0!")
          NumHistograms = 6
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 50d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "eta_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

ELSEIF( ObsSet.EQ.2 ) THEN! set of observables for ttb production as signal process at Tevatron (di-lept. decay)

          if(Collider.ne.2)  call Error("Collider needs to be TEV!")
          if(abs(TopDecays).ne.1) call Error("TopDecays needs to be 1!")
          NumHistograms = 20
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info    = "pT_lepMinus"
          Histo(1)%NBins   = 40
          Histo(1)%BinSize = 50d0*GeV
          Histo(1)%LowVal  = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_lepMinus"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 0.2d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_lepPlus"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 50d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_lepPlus"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "HT"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 50d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "Minv_lep_bquark"
          Histo(6)%NBins  = 40
          Histo(6)%BinSize= 10d0*GeV
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "Minv_lep_bquark"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 10d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

!          Histo(7)%Info   = "Psi"
!          Histo(7)%NBins  = 60
!          Histo(7)%BinSize= 0.1d0
!          Histo(7)%LowVal = -3d0
!          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "DeltaPhi"
          Histo(8)%NBins  = 65
          Histo(8)%BinSize= 0.1d0
          Histo(8)%LowVal = -3d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "E_bj1+E_bj2"
          Histo(9)%NBins  = 200
          Histo(9)%BinSize= 5d0*GeV
          Histo(9)%LowVal = 50d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "ET_miss"
          Histo(10)%NBins  = 50
          Histo(10)%BinSize= 5d0*GeV
          Histo(10)%LowVal = 20d0*GeV
          Histo(10)%SetScale= 100d0

          Histo(11)%Info   = "E_lep1+E_lep2"
          Histo(11)%NBins  = 50
          Histo(11)%BinSize= 5d0*GeV
          Histo(11)%LowVal = 20d0*GeV
          Histo(11)%SetScale= 100d0

          Histo(12)%Info   = "pt_bj1"
          Histo(12)%NBins  = 50
          Histo(12)%BinSize= 5d0*GeV
          Histo(12)%LowVal = 20d0*GeV
          Histo(12)%SetScale= 100d0

          Histo(13)%Info   = "pt_bj2"
          Histo(13)%NBins  = 50
          Histo(13)%BinSize= 5d0*GeV
          Histo(13)%LowVal = 20d0*GeV
          Histo(13)%SetScale= 100d0

          Histo(14)%Info   = "E_bj1"
          Histo(14)%NBins  = 50
          Histo(14)%BinSize= 5d0*GeV
          Histo(14)%LowVal = 20d0*GeV
          Histo(14)%SetScale= 100d0

          Histo(15)%Info   = "E_bj2"
          Histo(15)%NBins  = 50
          Histo(15)%BinSize= 5d0*GeV
          Histo(15)%LowVal = 20d0*GeV
          Histo(15)%SetScale= 100d0

          Histo(16)%Info   = "M_lj"
          Histo(16)%NBins  = 50
          Histo(16)%BinSize= 8d0*GeV
          Histo(16)%LowVal = 0d0
          Histo(16)%SetScale= 100d0

          Histo(17)%Info   = "ET_leptons"
          Histo(17)%NBins  = 50
          Histo(17)%BinSize= 5d0*GeV
          Histo(17)%LowVal = 20d0*GeV
          Histo(17)%SetScale= 100d0

          Histo(18)%Info   = "m_T"
          Histo(18)%NBins  = 50
          Histo(18)%BinSize= 10d0*GeV
          Histo(18)%LowVal = 20d0*GeV
          Histo(18)%SetScale= 100d0

          Histo(19)%Info   = "etaFB_lepMinus"
          Histo(19)%NBins  = 2
          Histo(19)%BinSize= 5d0
          Histo(19)%LowVal =-5.0d0
          Histo(19)%SetScale= 1d0

          Histo(20)%Info   = "etaFB_lepPlus"
          Histo(20)%NBins  = 2
          Histo(20)%BinSize= 5d0
          Histo(20)%LowVal =-5.0d0
          Histo(20)%SetScale= 1d0

ELSEIF( ObsSet.EQ.3 ) THEN! set of observables for ttb production as signal process at LHC (di-lept. decay)
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(abs(TopDecays).ne.1) call Error("TopDecays needs to be |1|!")
          NumHistograms = 20
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_ATop"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_lepPlus"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 25d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_lepPlus"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "HT"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 50d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "Minv_lep_bquark"
          Histo(6)%NBins  = 40
          Histo(6)%BinSize= 10d0*GeV
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "Psi"
          Histo(7)%NBins  = 60
          Histo(7)%BinSize= 0.1d0
          Histo(7)%LowVal = -3d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "DeltaPhi"
          Histo(8)%NBins  = 65
          Histo(8)%BinSize= 0.1d0
          Histo(8)%LowVal = -3d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "E_bj1+E_bj2"
          Histo(9)%NBins  = 200
          Histo(9)%BinSize= 5d0*GeV
          Histo(9)%LowVal = 50d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "ET_miss"
          Histo(10)%NBins  = 50
          Histo(10)%BinSize= 5d0*GeV
          Histo(10)%LowVal = 20d0*GeV
          Histo(10)%SetScale= 100d0

          Histo(11)%Info   = "E_lep1+E_lep2"
          Histo(11)%NBins  = 50
          Histo(11)%BinSize= 5d0*GeV
          Histo(11)%LowVal = 20d0*GeV
          Histo(11)%SetScale= 100d0

          Histo(12)%Info   = "pt_bj1"
          Histo(12)%NBins  = 40
          Histo(12)%BinSize= 20d0*GeV
          Histo(12)%LowVal = 0d0*GeV
          Histo(12)%SetScale= 100d0

          Histo(13)%Info   = "pt_bj2"
          Histo(13)%NBins  = 50
          Histo(13)%BinSize= 5d0*GeV
          Histo(13)%LowVal = 20d0*GeV
          Histo(13)%SetScale= 100d0

          Histo(14)%Info   = "E_bj1"
          Histo(14)%NBins  = 50
          Histo(14)%BinSize= 5d0*GeV
          Histo(14)%LowVal = 20d0*GeV
          Histo(14)%SetScale= 100d0

          Histo(15)%Info   = "E_bj2"
          Histo(15)%NBins  = 50
          Histo(15)%BinSize= 5d0*GeV
          Histo(15)%LowVal = 20d0*GeV
          Histo(15)%SetScale= 100d0

          Histo(16)%Info   = "Minv_leptons"
          Histo(16)%NBins  = 50
          Histo(16)%BinSize= 10d0*GeV
          Histo(16)%LowVal = 0d0
          Histo(16)%SetScale= 100d0

          Histo(17)%Info   = "ET_leptons"
          Histo(17)%NBins  = 50
          Histo(17)%BinSize= 5d0*GeV
          Histo(17)%LowVal = 20d0*GeV
          Histo(17)%SetScale= 100d0

          Histo(18)%Info   = "m_T"
          Histo(18)%NBins  = 50
          Histo(18)%BinSize= 10d0*GeV
          Histo(18)%LowVal = 20d0*GeV
          Histo(18)%SetScale= 100d0

          Histo(19)%Info   = "etaFB_lepMinus"
          Histo(19)%NBins  = 2
          Histo(19)%BinSize= 5d0
          Histo(19)%LowVal =-5.0d0
          Histo(19)%SetScale= 1d0

          Histo(20)%Info   = "etaFB_lepPlus"
          Histo(20)%NBins  = 2
          Histo(20)%BinSize= 5d0
          Histo(20)%LowVal =-5.0d0
          Histo(20)%SetScale= 1d0

ELSEIF( ObsSet.EQ.4 ) THEN! set of observables for ttb production with semi hadronic decay
          if(Collider.ne.2)  call Error("Collider needs to be TEV!")
          if(TopDecays.ne.4) call Error("TopDecays needs to be 4!")
          NumHistograms = 12
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_ATop"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_Top"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 50d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 40
          Histo(4)%BinSize= 0.25d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Lep"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_LepP"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 20d0*GeV
          Histo(7)%LowVal =  0d0*GeV
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "eta_LepP"
          Histo(8)%NBins  = 40
          Histo(8)%BinSize= 0.25d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "pT_miss"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 20d0*GeV
          Histo(9)%LowVal =  0d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "HT(jets+lept+miss)"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 25d0*GeV
          Histo(10)%LowVal = 150d0*GeV
          Histo(10)%SetScale= 100d0

          Histo(11)%Info   = "m(lep+bjet)"
          Histo(11)%NBins  = 40
          Histo(11)%BinSize= 20d0*GeV
          Histo(11)%LowVal = 20d0*GeV
          Histo(11)%SetScale= 100d0

          Histo(12)%Info   = "pT_ttbar"
          Histo(12)%NBins  = 40
          Histo(12)%BinSize= 15d0*GeV
          Histo(12)%LowVal = 0d0
          Histo(12)%SetScale= 100d0




ELSEIF( ObsSet.EQ.5 ) THEN! set of observables for ttb production with hadr. top, lept. Atop decay
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.4) call Error("TopDecays needs to be 4!")
          NumHistograms = 11
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "y_ATop"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_Top"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 50d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "y_Top"
          Histo(4)%NBins  = 40
          Histo(4)%BinSize= 0.25d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_LepP"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 25d0*GeV
          Histo(7)%LowVal =  0d0*GeV
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "y_LepP"
          Histo(8)%NBins  = 40
          Histo(8)%BinSize= 0.25d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "pT_miss"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 25d0*GeV
          Histo(9)%LowVal =  0d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "HT(jets+lept)"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 50d0*GeV
          Histo(10)%LowVal = 150d0*GeV
          Histo(10)%SetScale= 100d0

          Histo(11)%Info   = "m(lep+bjet)"
          Histo(11)%NBins  = 90
          Histo(11)%BinSize= 5d0*GeV
          Histo(11)%LowVal = 20d0*GeV
          Histo(11)%SetScale= 100d0




ELSEIF( ObsSet.EQ.6 ) THEN! set of observables for ttb production with lept. top, hadr. Atop decay
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.4) call Error("TopDecays needs to be 4!")
          NumHistograms = 12
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif


          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_ATop"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_Top"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 50d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 40
          Histo(4)%BinSize= 0.25d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_LepP"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 25d0*GeV
          Histo(7)%LowVal =  0d0*GeV
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "eta_LepP"
          Histo(8)%NBins  = 40
          Histo(8)%BinSize= 0.25d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "pT_miss"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 25d0*GeV
          Histo(9)%LowVal =  0d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "HT(jets+lept)"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 50d0*GeV
          Histo(10)%LowVal = 150d0*GeV
          Histo(10)%SetScale= 100d0

          Histo(11)%Info   = "m(lep+bjet)"
          Histo(11)%NBins  = 40
          Histo(11)%BinSize= 20d0*GeV
          Histo(11)%LowVal = 20d0*GeV
          Histo(11)%SetScale= 100d0

          Histo(12)%Info   = "pT_ttbar"
          Histo(12)%NBins  = 40
          Histo(12)%BinSize= 15d0*GeV
          Histo(12)%LowVal = 0d0
          Histo(12)%SetScale= 100d0



ELSEIF( ObsSet.EQ.7 ) THEN! set of observables for ttb production with lept. top and J/Psi fragmentation, hadr. Atop decay at LHC
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.6) call Error("TopDecays needs to be 6!")
          NumHistograms = 6
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif
          Histo(1)%Info   = "M_LJPsi"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 5d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_JPsi"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 10d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "pT_LepPlus"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 10d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "<M_LJPsi>"
          Histo(4)%NBins  = 1
          Histo(4)%BinSize= 10000d0*GeV
          Histo(4)%LowVal = 0d0
          Histo(4)%SetScale=1d0

          Histo(5)%Info   = "<M_LJPsi^2>"
          Histo(5)%NBins  = 1
          Histo(5)%BinSize= 10000d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale=1d0

          Histo(6)%Info   = "<x>"
          Histo(6)%NBins  = 1
          Histo(6)%BinSize= 10000d0*GeV
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale=1d0


ELSEIF( ObsSet.EQ.8 ) THEN! set of observables for ttb spin correlations at LHC (di-lept. decay)

!           if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(abs(TopDecays).ne.1) call Error("TopDecays needs to be 1!")
          NumHistograms = 4
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif
          Histo(1)%Info   = "r"
          Histo(1)%NBins  = 20
          Histo(1)%BinSize= 0.05d0
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 1d0

          Histo(2)%Info    = "pT_lepMinus"
          Histo(2)%NBins   = 40
          Histo(2)%BinSize = 50d0*GeV
          Histo(2)%LowVal  = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "pT_lepPlus"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 50d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "DeltaPhi"
          Histo(4)%NBins  = 65
          Histo(4)%BinSize= 0.1d0
          Histo(4)%LowVal = -3d0
          Histo(4)%SetScale= 1d0

ELSEIF( ObsSet.EQ.9 ) THEN! this is for the factorization check
          NumHistograms = 6
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 50d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "eta_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0


ELSEIF( ObsSet.EQ.10 ) THEN! set of observables for ttbjet production without decays at Tevatron
          if(Collider.ne.2)  call Error("Collider needs to be Tevatron!")
          if(TopDecays.ne.0) call Error("TopDecays needs to be 0!")
          NumHistograms = 8
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 20d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 20d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "eta_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_jet"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 20d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "eta_jet"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 0.2d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0


ELSEIF( ObsSet.EQ.11 ) THEN! set of observables for ttbjet production without decays at LHC
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.0) call Error("TopDecays needs to be 0!")
          NumHistograms = 8
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 50d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "eta_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_jet"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 50d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "eta_jet"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 0.2d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0


ELSEIF( ObsSet.EQ.12 ) THEN! set of observables for ttbjet production as signal process at Tevatron (hadr.Atop, lept.top decay)
          if(Collider.ne.2)  call Error("Collider needs to be Tevatron!")
          if(TopDecays.ne.4) call Error("TopDecays needs to be 4!")
          NumHistograms = 44
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif
!          if( .not.allocated(Histo2D) ) then
!                allocate( Histo2D(1:3), stat=AllocStatus  )
!                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo2D")
!          endif


          Histo(1)%Info   = "pT_lepPlus"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 20d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_lepPlus"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_jet"
          Histo(3)%NBins  = 20
          Histo(3)%BinSize= 20d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_jet"
          Histo(4)%NBins  = 20
          Histo(4)%BinSize= 0.4d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "pt_5th_jet"
          Histo(5)%NBins  = 40
          Histo(5)%BinSize= 20d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "eta_5th_jet"
          Histo(6)%NBins  = 40
          Histo(6)%BinSize= 0.4d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_miss"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 20d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "HT"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 25d0*GeV
          Histo(8)%LowVal = 100d0*GeV
          Histo(8)%SetScale= 100d0

          Histo(9)%Info   = "m_lb"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 20d0*GeV
          Histo(9)%LowVal = 0d0
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "phi(l,b)"
          Histo(10)%NBins  = 15
          Histo(10)%BinSize= 0.25d0
          Histo(10)%LowVal = 0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "R(l,b)"
          Histo(11)%NBins  = 40
          Histo(11)%BinSize= 0.2d0
          Histo(11)%LowVal = 0d0
          Histo(11)%SetScale= 1d0

          Histo(12)%Info   = "eta_FBlepPlus"
          Histo(12)%NBins  = 2
          Histo(12)%BinSize= 5d0
          Histo(12)%LowVal =-5.0d0
          Histo(12)%SetScale= 1d0

          Histo(13)%Info   = "pT(ttbar)"
          Histo(13)%NBins  = 50
          Histo(13)%BinSize= 10d0*GeV
          Histo(13)%LowVal = 0d0
          Histo(13)%SetScale= 100d0

! new stuff for A_FB analysis, this is for ideally reconstructed tops
          Histo(14)%Info   = "pT(ttbar) FWD"
          Histo(14)%NBins  = 50
          Histo(14)%BinSize= 10d0*GeV
          Histo(14)%LowVal =  0d0*GeV
          Histo(14)%SetScale= 100d0

          Histo(15)%Info   = "pT(ttbar) BWD"
          Histo(15)%NBins  = 50
          Histo(15)%BinSize= 10d0*GeV
          Histo(15)%LowVal =  0d0*GeV
          Histo(15)%SetScale= 100d0

          Histo(16)%Info   = "m(ttbar) FWD"
          Histo(16)%NBins  = 50
          Histo(16)%BinSize= 50d0*GeV
          Histo(16)%LowVal =  300d0*GeV
          Histo(16)%SetScale= 100d0

          Histo(17)%Info   = "m(ttbar) BWD"
          Histo(17)%NBins  = 50
          Histo(17)%BinSize= 50d0*GeV
          Histo(17)%LowVal =  300d0*GeV
          Histo(17)%SetScale= 100d0

          Histo(18)%Info   = "y(ttbar) FWD"
          Histo(18)%NBins  = 40
          Histo(18)%BinSize= 0.5d0
          Histo(18)%LowVal =-5.0d0
          Histo(18)%SetScale= 1d0

          Histo(19)%Info   = "y(ttbar) BWD"
          Histo(19)%NBins  = 40
          Histo(19)%BinSize= 0.5d0
          Histo(19)%LowVal =-5.0d0
          Histo(19)%SetScale= 1d0

          Histo(20)%Info   = "dy(tops) FWD"
          Histo(20)%NBins  = 40
          Histo(20)%BinSize= 0.25d0
          Histo(20)%LowVal = 0.0d0
          Histo(20)%SetScale= 1d0

          Histo(21)%Info   = "dy(tops) BWD"
          Histo(21)%NBins  = 40
          Histo(21)%BinSize= 0.25d0
          Histo(21)%LowVal = 0.0d0
          Histo(21)%SetScale= 1d0

          Histo(22)%Info   = "pT(ttbar) FWD lept"
          Histo(22)%NBins  = 50
          Histo(22)%BinSize= 10d0*GeV
          Histo(22)%LowVal =  0d0*GeV
          Histo(22)%SetScale= 100d0

          Histo(23)%Info   = "pT(ttbar) BWD lept"
          Histo(23)%NBins  = 50
          Histo(23)%BinSize= 10d0*GeV
          Histo(23)%LowVal =  0d0*GeV
          Histo(23)%SetScale= 100d0

          Histo(24)%Info   = "y_FB(top)"
          Histo(24)%NBins  = 2
          Histo(24)%BinSize= 5d0
          Histo(24)%LowVal =-5.0d0
          Histo(24)%SetScale= 1d0

          Histo(25)%Info   = "pT(top)"
          Histo(25)%NBins  = 40
          Histo(25)%BinSize= 20d0*GeV
          Histo(25)%LowVal = 0d0
          Histo(25)%SetScale= 100d0

          Histo(26)%Info   = "phi(ttbar) FWD" 
          Histo(26)%NBins  = 15
          Histo(26)%BinSize= 0.25d0
          Histo(26)%LowVal = 0d0
          Histo(26)%SetScale= 1d0

          Histo(27)%Info   = "phi(ttbar) BWD" 
          Histo(27)%NBins  = 15
          Histo(27)%BinSize= 0.25d0
          Histo(27)%LowVal = 0d0
          Histo(27)%SetScale= 1d0


! new stuff for A_FB analysis, this is for realistically reconstructed tops
          Histo(28)%Info   = "pT(ttbar)"
          Histo(28)%NBins  = 50
          Histo(28)%BinSize= 10d0*GeV
          Histo(28)%LowVal = 0d0
          Histo(28)%SetScale= 100d0

          Histo(29)%Info   = "pT(ttbar) FWD"
          Histo(29)%NBins  = 50
          Histo(29)%BinSize= 10d0*GeV
          Histo(29)%LowVal =  0d0*GeV
          Histo(29)%SetScale= 100d0

          Histo(30)%Info   = "pT(ttbar) BWD"
          Histo(30)%NBins  = 50
          Histo(30)%BinSize= 10d0*GeV
          Histo(30)%LowVal =  0d0*GeV
          Histo(30)%SetScale= 100d0

          Histo(31)%Info   = "m(ttbar) FWD"
          Histo(31)%NBins  = 50
          Histo(31)%BinSize= 50d0*GeV
          Histo(31)%LowVal =  300d0*GeV
          Histo(31)%SetScale= 100d0

          Histo(32)%Info   = "m(ttbar) BWD"
          Histo(32)%NBins  = 50
          Histo(32)%BinSize= 50d0*GeV
          Histo(32)%LowVal =  300d0*GeV
          Histo(32)%SetScale= 100d0

          Histo(33)%Info   = "y(ttbar) FWD"
          Histo(33)%NBins  = 40
          Histo(33)%BinSize= 0.5d0
          Histo(33)%LowVal =-5.0d0
          Histo(33)%SetScale= 1d0

          Histo(34)%Info   = "y(ttbar) BWD"
          Histo(34)%NBins  = 40
          Histo(34)%BinSize= 0.5d0
          Histo(34)%LowVal =-5.0d0
          Histo(34)%SetScale= 1d0

          Histo(35)%Info   = "dy(tops) FWD"
          Histo(35)%NBins  = 40
          Histo(35)%BinSize= 0.25d0
          Histo(35)%LowVal = 0.0d0
          Histo(35)%SetScale= 1d0

          Histo(36)%Info   = "dy(tops) BWD"
          Histo(36)%NBins  = 40
          Histo(36)%BinSize= 0.25d0
          Histo(36)%LowVal = 0.0d0
          Histo(36)%SetScale= 1d0

          Histo(37)%Info   = "pT(ttbar) FWD lept"
          Histo(37)%NBins  = 50
          Histo(37)%BinSize= 10d0*GeV
          Histo(37)%LowVal =  0d0*GeV
          Histo(37)%SetScale= 100d0

          Histo(38)%Info   = "pT(ttbar) BWD lept"
          Histo(38)%NBins  = 50
          Histo(38)%BinSize= 10d0*GeV
          Histo(38)%LowVal =  0d0*GeV
          Histo(38)%SetScale= 100d0

          Histo(39)%Info   = "y_FB(top)"
          Histo(39)%NBins  = 2
          Histo(39)%BinSize= 5d0
          Histo(39)%LowVal =-5.0d0
          Histo(39)%SetScale= 1d0

          Histo(40)%Info   = "pT(top)"
          Histo(40)%NBins  = 40
          Histo(40)%BinSize= 20d0*GeV
          Histo(40)%LowVal = 0d0
          Histo(40)%SetScale= 100d0

          Histo(41)%Info   = "phi(ttbar) FWD" 
          Histo(41)%NBins  = 15
          Histo(41)%BinSize= 0.25d0
          Histo(41)%LowVal = 0d0
          Histo(41)%SetScale= 1d0

          Histo(42)%Info   = "phi(ttbar) BWD" 
          Histo(42)%NBins  = 15
          Histo(42)%BinSize= 0.25d0
          Histo(42)%LowVal = 0d0
          Histo(42)%SetScale= 1d0

          Histo(43)%Info   = "cos(theta_top*)" 
          Histo(43)%NBins  = 20
          Histo(43)%BinSize= 0.1d0
          Histo(43)%LowVal = -1d0
          Histo(43)%SetScale= 1d0

          Histo(44)%Info   = "beta_top*cos(theta_top*)" 
          Histo(44)%NBins  = 20
          Histo(44)%BinSize= 0.1d0
          Histo(44)%LowVal = -1d0
          Histo(44)%SetScale= 1d0

! 2D histograms for realistically reconstructed tops      
!
!          Histo2D(1)%Info   = "m_ttbar over pT_jet"  
!          Histo2D(1)%NBins(1)  = 50
!          Histo2D(1)%BinSize(1)= 50d0
!          Histo2D(1)%LowVal(1) = 300d0
!          Histo2D(1)%SetScale(1)= 100d0
!          Histo2D(1)%NBins(2)  = 50
!          Histo2D(1)%BinSize(2)= 10d0
!          Histo2D(1)%LowVal(2) = 0d0
!          Histo2D(1)%SetScale(2)= 100d0
!
!          Histo2D(2)%Info   = "m_ttbar over pT_ttbar FWD"  
!          Histo2D(2)%NBins(1)  = 50
!          Histo2D(2)%BinSize(1)= 50d0
!          Histo2D(2)%LowVal(1) = 300d0
!          Histo2D(2)%SetScale(1)= 100d0
!          Histo2D(2)%NBins(2)  = 50
!          Histo2D(2)%BinSize(2)= 20d0
!          Histo2D(2)%LowVal(2) = 0d0
!          Histo2D(2)%SetScale(2)= 100d0
!
!          Histo2D(3)%Info   = "m_ttbar over pT_ttbar BWD"  
!          Histo2D(3)%NBins(1)  = 50
!          Histo2D(3)%BinSize(1)= 50d0
!          Histo2D(3)%LowVal(1) = 300d0
!          Histo2D(3)%SetScale(1)= 100d0
!          Histo2D(3)%NBins(2)  = 50
!          Histo2D(3)%BinSize(2)= 20d0
!          Histo2D(3)%LowVal(2) = 0d0
!          Histo2D(3)%SetScale(2)= 100d0





ELSEIF( ObsSet.EQ.13 ) THEN! set of observables for ttbjet production as signal process at LHC (di-lept. decay)
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.1) call Error("TopDecays needs to be 1!")
          if(Collider_Energy.ne.7000d0*GeV) call Error("wrong collider energy for ObsSet=13!")
          NumHistograms = 16
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif


          Histo(1)%Info    = "pT_lepMinus"
          Histo(1)%NBins   = 40
          Histo(1)%BinSize = 20d0*GeV
          Histo(1)%LowVal  = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_lepMinus"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_lepPlus"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 20d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_lepPlus"
          Histo(4)%NBins  = 40
          Histo(4)%BinSize= 0.25d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "HT"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 50d0*GeV
          Histo(5)%LowVal = 100d0*GeV
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "Minv_leptons"
          Histo(6)%NBins  = 40
          Histo(6)%BinSize= 20d0*GeV
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "pT_jet"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 20d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "eta_jet"
          Histo(8)%NBins  = 40
          Histo(8)%BinSize= 0.25d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "pT_miss"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 20d0*GeV
          Histo(9)%LowVal = 0d0
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "pT_ttbar"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 15d0*GeV
          Histo(10)%LowVal = 0d0
          Histo(10)%SetScale= 100d0

          Histo(11)%Info   = "phi(l+,l-)"
          Histo(11)%NBins  = 15
          Histo(11)%BinSize= 0.25d0
          Histo(11)%LowVal = 0d0
          Histo(11)%SetScale= 1d0

          Histo(12)%Info   = "m_l+l-"
          Histo(12)%NBins  = 50
          Histo(12)%BinSize= 20d0*GeV
          Histo(12)%LowVal = 0d0
          Histo(12)%SetScale= 100d0

          Histo(13)%Info   = "pT_bjet"
          Histo(13)%NBins  = 40
          Histo(13)%BinSize= 20d0*GeV
          Histo(13)%LowVal = 0d0
          Histo(13)%SetScale= 100d0

          Histo(14)%Info   = "eta_bjet"
          Histo(14)%NBins  = 50
          Histo(14)%BinSize= 0.4d0
          Histo(14)%LowVal =-5.0d0
          Histo(14)%SetScale= 1d0

          Histo(15)%Info   = "m_bb"
          Histo(15)%NBins  = 50
          Histo(15)%BinSize= 20d0*GeV
          Histo(15)%LowVal = 0d0
          Histo(15)%SetScale= 100d0

          Histo(16)%Info   = "m_bj"
          Histo(16)%NBins  = 50
          Histo(16)%BinSize= 20d0*GeV
          Histo(16)%LowVal = 0d0
          Histo(16)%SetScale= 100d0




ELSEIF( ObsSet.EQ.14 ) THEN! set of observables for ttbjet production as background process to VBF at LHC (di-lept. decay)
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.1) call Error("TopDecays needs to be 1!")
          if(Collider_Energy.ne.14000d0*GeV) call Error("wrong collider energy for ObsSet=14!")
          NumHistograms = 11
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info    = "pT_lepMinus"
          Histo(1)%NBins   = 40
          Histo(1)%BinSize = 50d0*GeV
          Histo(1)%LowVal  = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_lepMinus"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 0.2d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_lepPlus"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 50d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_lepPlus"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "Minv_leptons"
          Histo(5)%NBins  = 40
          Histo(5)%BinSize= 50d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "pT_jet"
          Histo(6)%NBins  = 40
          Histo(6)%BinSize= 50d0*GeV
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "phi_leptons"
          Histo(7)%NBins  = 20
          Histo(7)%BinSize= 0.2d0
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "mT"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 20d0*GeV
          Histo(8)%LowVal = 0d0
          Histo(8)%SetScale= 100d0

          Histo(9)%Info   = "eta(hardest jet)"
          Histo(9)%NBins  = 50
          Histo(9)%BinSize= 0.2d0
          Histo(9)%LowVal =-5.0d0
          Histo(9)%SetScale= 1d0

          Histo(10)%Info   = "eta(2nd hardest jet)"
          Histo(10)%NBins  = 50
          Histo(10)%BinSize= 0.2d0
          Histo(10)%LowVal =-5.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "eta_Zeppi"
          Histo(11)%NBins  = 50
          Histo(11)%BinSize= 0.2d0
          Histo(11)%LowVal =-5.0d0
          Histo(11)%SetScale= 1d0



ELSEIF( ObsSet.EQ.15 ) THEN! set of observables for ttbjet production as signal process at LHC (hadr.Atop, lept.top decay)
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.4) call Error("TopDecays needs to be 4!")
          if(Collider_Energy.ne.7000d0*GeV) call Error("wrong collider energy for ObsSet=15!")
          NumHistograms = 13
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif


          Histo(1)%Info   = "pT_lepPlus"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 20d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_lepPlus"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_jet"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 20d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_jet"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.4d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "pt_5th_jet"
          Histo(5)%NBins  = 40
          Histo(5)%BinSize= 20d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "eta_5th_jet"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 0.4d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_miss"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 20d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "HT"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 50d0*GeV
          Histo(8)%LowVal = 100d0*GeV
          Histo(8)%SetScale= 100d0

          Histo(9)%Info   = "m_lb"
          Histo(9)%NBins  = 50
          Histo(9)%BinSize= 20d0*GeV
          Histo(9)%LowVal = 0d0
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "phi(l,b)"
          Histo(10)%NBins  = 15
          Histo(10)%BinSize= 0.25d0
          Histo(10)%LowVal = 0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "R(l,b)"
          Histo(11)%NBins  = 40
          Histo(11)%BinSize= 0.2d0
          Histo(11)%LowVal = 0d0
          Histo(11)%SetScale= 1d0

          Histo(12)%Info   = "eta_FBlepPlus"
          Histo(12)%NBins  = 2
          Histo(12)%BinSize= 5d0
          Histo(12)%LowVal =-5.0d0
          Histo(12)%SetScale= 1d0

          Histo(13)%Info   = "pT_ttbar"
          Histo(13)%NBins  = 40
          Histo(13)%BinSize= 15d0*GeV
          Histo(13)%LowVal = 0d0
          Histo(13)%SetScale= 100d0





ELSEIF( ObsSet.EQ.19 ) THEN! for checks of ttbjet
          NumHistograms = 10
          if(TopDecays.eq.0) call Error("TopDecays needs to be >0!")
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_lepPlus"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 20d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_lepPlus"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_jet"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 20d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_jet"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.4d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "pT_miss"
          Histo(5)%NBins  = 40
          Histo(5)%BinSize= 20d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "HT"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 20d0*GeV
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "m_lb"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 20d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "phi(l,b)"
          Histo(8)%NBins  = 40
          Histo(8)%BinSize= 0.2d0
          Histo(8)%LowVal = 0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "R(l,b)"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 0.2d0
          Histo(9)%LowVal = 0d0
          Histo(9)%SetScale= 1d0

          Histo(10)%Info   = "eta_FBlepPlus"
          Histo(10)%NBins  = 2
          Histo(10)%BinSize= 5d0
          Histo(10)%LowVal =-5.0d0
          Histo(10)%SetScale= 1d0


ELSEIF( ObsSet.EQ.20 ) THEN! set of observables for ttbgamma production without decays at Tevatron
          if(Collider.ne.2)  call Error("Collider needs to be Tevatron!")
          if(TopDecays.ne.0) call Error("TopDecays needs to be 0!")
          NumHistograms = 8
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 20d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 20d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "eta_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "pT_Photon"
          Histo(5)%NBins  = 40
          Histo(5)%BinSize= 10d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "eta_Photon"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 0.2d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "etaFB_ATop"
          Histo(7)%NBins  = 2
          Histo(7)%BinSize= 5d0
          Histo(7)%LowVal =-5.0d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "etaFB_Top"
          Histo(8)%NBins  = 2
          Histo(8)%BinSize= 5d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0


ELSEIF( ObsSet.EQ.21 ) THEN! set of observables for ttbgamma production without decays at the LHC
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.0) call Error("TopDecays needs to be 0!")
          NumHistograms = 8
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 25d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 25d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "eta_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "pT_Photon"
          Histo(5)%NBins  = 40
          Histo(5)%BinSize= 20d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "eta_Photon"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 0.2d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "etaFB_ATop"
          Histo(7)%NBins  = 2
          Histo(7)%BinSize= 5d0
          Histo(7)%LowVal =-5.0d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "etaFB_Top"
          Histo(8)%NBins  = 2
          Histo(8)%BinSize= 5d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0


ELSEIF( ObsSet.EQ.22 ) THEN! set of observables for ttbgamma production di-lept. decays at the Tevatron
          if(Collider.ne.2)  call Error("Collider needs to be TEV!")
          if(TopDecays.ne.1  ) call Error("TopDecays needs to be 1!")
          NumHistograms = 14
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 25d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 25d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "eta_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_Photon"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 20d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "eta_Photon"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 0.2d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "pT_LepP"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 20d0*GeV
          Histo(9)%LowVal = 20d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "eta_LepP"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 0.25d0
          Histo(10)%LowVal =-5.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "ET_miss"
          Histo(11)%NBins  = 40
          Histo(11)%BinSize= 20d0*GeV
          Histo(11)%LowVal = 20d0*GeV
          Histo(11)%SetScale= 100d0

          Histo(12)%Info   = "HT(jets+lept+pho+miss)"
          Histo(12)%NBins  = 40
          Histo(12)%BinSize= 20d0*GeV
          Histo(12)%LowVal = 150d0*GeV
          Histo(12)%SetScale= 100d0

          Histo(13)%Info   = "m(LepP+bjet)"
          Histo(13)%NBins  = 40
          Histo(13)%BinSize= 20d0*GeV
          Histo(13)%LowVal = 20d0*GeV
          Histo(13)%SetScale= 100d0

          Histo(14)%Info   = "phi(LepP,LepM)"
          Histo(14)%NBins  = 20
          Histo(14)%BinSize= 0.2d0
          Histo(14)%LowVal = 0d0
          Histo(14)%SetScale= 1d0


ELSEIF( ObsSet.EQ.23 ) THEN! set of observables for ttbgamma production di-lept. decays at the LHC
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.1 ) call Error("TopDecays needs to be 1!")
          NumHistograms = 14
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 25d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 25d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "eta_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_Photon"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 20d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "eta_Photon"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 0.2d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "pT_LepP"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 20d0*GeV
          Histo(9)%LowVal = 20d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "eta_LepP"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 0.25d0
          Histo(10)%LowVal =-5.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "ET_miss"
          Histo(11)%NBins  = 40
          Histo(11)%BinSize= 20d0*GeV
          Histo(11)%LowVal = 20d0*GeV
          Histo(11)%SetScale= 100d0

          Histo(12)%Info   = "HT(jets+lept+pho+miss)"
          Histo(12)%NBins  = 40
          Histo(12)%BinSize= 20d0*GeV
          Histo(12)%LowVal = 150d0*GeV
          Histo(12)%SetScale= 100d0

          Histo(13)%Info   = "m(LepP+bjet)"
          Histo(13)%NBins  = 40
          Histo(13)%BinSize= 20d0*GeV
          Histo(13)%LowVal = 20d0*GeV
          Histo(13)%SetScale= 100d0

          Histo(14)%Info   = "phi(LepP,LepM)"
          Histo(14)%NBins  = 20
          Histo(14)%BinSize= 0.2d0
          Histo(14)%LowVal = 0d0
          Histo(14)%SetScale= 1d0



ELSEIF( ObsSet.EQ.24 ) THEN! set of observables for ttbgamma production semi-lept. decays at the TEV
          if(Collider.ne.2)  call Error("Collider needs to be Tevatron!")
          if(TopDecays.ne.4  ) call Error("TopDecays needs to be 4!")
          NumHistograms = 17
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_ATop"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_Top"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 50d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 40
          Histo(4)%BinSize= 0.25d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_Lep"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Pho"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_Photon"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 20d0*GeV
          Histo(7)%LowVal =  0d0*GeV
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "eta_Photon"
          Histo(8)%NBins  = 40
          Histo(8)%BinSize= 0.25d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "pT_LepP"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 20d0*GeV
          Histo(9)%LowVal =  0d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "eta_LepP"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 0.25d0
          Histo(10)%LowVal =-5.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "ET_miss"
          Histo(11)%NBins  = 40
          Histo(11)%BinSize= 20d0*GeV
          Histo(11)%LowVal =  0d0*GeV
          Histo(11)%SetScale= 100d0

          Histo(12)%Info   = "HT(jets+lept+pho+miss)"
          Histo(12)%NBins  = 40
          Histo(12)%BinSize= 25d0*GeV
          Histo(12)%LowVal = 150d0*GeV
          Histo(12)%SetScale= 100d0

          Histo(13)%Info   = "R(pho,bjet)"
          Histo(13)%NBins  = 25
          Histo(13)%BinSize= 0.25d0
          Histo(13)%LowVal = 0d0
          Histo(13)%SetScale= 1d0

          Histo(14)%Info   = "m(lep+bjet)"
          Histo(14)%NBins  = 40
          Histo(14)%BinSize= 20d0*GeV
          Histo(14)%LowVal = 20d0*GeV
          Histo(14)%SetScale= 100d0

          Histo(15)%Info   = "mT(lep,pho;pTmiss)"
          Histo(15)%NBins  = 40
          Histo(15)%BinSize= 20d0*GeV
          Histo(15)%LowVal = 20d0*GeV
          Histo(15)%SetScale= 100d0

          Histo(16)%Info   = "phi(photon,lept)"
          Histo(16)%NBins  = 20
          Histo(16)%BinSize= 0.2d0
          Histo(16)%LowVal = 0d0
          Histo(16)%SetScale= 1d0

          Histo(17)%Info   = "ET_bjet"
          Histo(17)%NBins  = 40
          Histo(17)%BinSize= 20d0*GeV
          Histo(17)%LowVal =  0d0*GeV
          Histo(17)%SetScale= 100d0


ELSEIF( ObsSet.EQ.25 ) THEN! set of observables for ttbgamma production semi-lept. decays at the LHC
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.4) call Error("TopDecays needs to be 3 (for Qt=-4/3) or 4 (for Qt=Qup)")
          NumHistograms = 15
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_ATop"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_Top"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 50d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 40
          Histo(4)%BinSize= 0.25d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_Photon"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 25d0*GeV
          Histo(7)%LowVal = 0d0*GeV
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "eta_Photon"
          Histo(8)%NBins  = 40
          Histo(8)%BinSize= 0.25d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "pT_LepP"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 25d0*GeV
          Histo(9)%LowVal =  0d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "eta_LepP"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 0.25d0
          Histo(10)%LowVal =-5.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "pT_miss"
          Histo(11)%NBins  = 40
          Histo(11)%BinSize= 25d0*GeV
          Histo(11)%LowVal =  0d0*GeV
          Histo(11)%SetScale= 100d0

          Histo(12)%Info   = "HT(jets+lept+pho)"
          Histo(12)%NBins  = 40
          Histo(12)%BinSize= 50d0*GeV
          Histo(12)%LowVal = 150d0*GeV
          Histo(12)%SetScale= 100d0

          Histo(13)%Info   = "R(pho,bjet)"
          Histo(13)%NBins  = 25
          Histo(13)%BinSize= 0.25d0
          Histo(13)%LowVal = 0d0
          Histo(13)%SetScale= 1d0

          Histo(14)%Info   = "m(lep+bjet)"
          Histo(14)%NBins  = 90
          Histo(14)%BinSize= 5d0*GeV
          Histo(14)%LowVal = 20d0*GeV
          Histo(14)%SetScale= 100d0

          Histo(15)%Info   = "phi(photon,lept)"
          Histo(15)%NBins  = 20
          Histo(15)%BinSize= 0.2d0
          Histo(15)%LowVal = 0d0
          Histo(15)%SetScale= 1d0


ELSEIF( ObsSet.EQ.26 ) THEN! set of observables for ttbgamma production semi-lept. decays at the TEV
          if(Collider.ne.2)  call Error("Collider needs to be Tevatron!")
          if(TopDecays.ne.4  .and. TopDecays.ne.3 ) call Error("TopDecays needs to be 3 (for Qt=-4/3) or 4 (for Qt=Qup)")
          NumHistograms = 15
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_ATop"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_Top"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 50d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 40
          Histo(4)%BinSize= 0.25d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_Photon"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 20d0*GeV
          Histo(7)%LowVal = 20d0*GeV
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "eta_Photon"
          Histo(8)%NBins  = 40
          Histo(8)%BinSize= 0.25d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "pT_LepP"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 20d0*GeV
          Histo(9)%LowVal = 20d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "eta_LepP"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 0.25d0
          Histo(10)%LowVal =-5.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "pT_miss"
          Histo(11)%NBins  = 40
          Histo(11)%BinSize= 20d0*GeV
          Histo(11)%LowVal = 20d0*GeV
          Histo(11)%SetScale= 100d0

          Histo(12)%Info   = "HT(jets+lept+pho)"
          Histo(12)%NBins  = 40
          Histo(12)%BinSize= 20d0*GeV
          Histo(12)%LowVal = 150d0*GeV
          Histo(12)%SetScale= 100d0

          Histo(13)%Info   = "R(pho,bjet)"
          Histo(13)%NBins  = 25
          Histo(13)%BinSize= 0.2d0
          Histo(13)%LowVal = 0d0
          Histo(13)%SetScale= 1d0

          Histo(14)%Info   = "m(lep+bjet)"
          Histo(14)%NBins  = 40
          Histo(14)%BinSize= 20d0*GeV
          Histo(14)%LowVal = 20d0*GeV
          Histo(14)%SetScale= 100d0

          Histo(15)%Info   = "phi(photon,lept)"
          Histo(15)%NBins  = 20
          Histo(15)%BinSize= 0.2d0
          Histo(15)%LowVal = 0d0
          Histo(15)%SetScale= 1d0

ELSEIF( ObsSet.EQ.27 ) THEN! set of observables for ttbgamma production semi-lept. decays at the LHC for Q_top measurement
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.4 .and. TopDecays.ne.3) call Error("TopDecays needs to be 3 (for Qt=-4/3) or 4 (for Qt=Qup)")
          if( Q_top.ne.Q_up .and. TopDecays.ne.3 ) call Error("TopDecays needs to be 3 for Qt=-4/3")
          NumHistograms = 15
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_ATop"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_Top"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 50d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 40
          Histo(4)%BinSize= 0.25d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_Photon"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 25d0*GeV
          Histo(7)%LowVal = 0d0*GeV
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "eta_Photon"
          Histo(8)%NBins  = 40
          Histo(8)%BinSize= 0.25d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "pT_LepP"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 25d0*GeV
          Histo(9)%LowVal =  0d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "eta_LepP"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 0.25d0
          Histo(10)%LowVal =-5.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "pT_miss"
          Histo(11)%NBins  = 40
          Histo(11)%BinSize= 25d0*GeV
          Histo(11)%LowVal =  0d0*GeV
          Histo(11)%SetScale= 100d0

          Histo(12)%Info   = "HT(jets+lept+pho)"
          Histo(12)%NBins  = 40
          Histo(12)%BinSize= 50d0*GeV
          Histo(12)%LowVal = 150d0*GeV
          Histo(12)%SetScale= 100d0

          Histo(13)%Info   = "R(pho,bjet)"
          Histo(13)%NBins  = 25
          Histo(13)%BinSize= 0.25d0
          Histo(13)%LowVal = 0d0
          Histo(13)%SetScale= 1d0

          Histo(14)%Info   = "m(lep+bjet)"
          Histo(14)%NBins  = 40
          Histo(14)%BinSize= 20d0*GeV
          Histo(14)%LowVal = 20d0*GeV
          Histo(14)%SetScale= 100d0

          Histo(15)%Info   = "phi(photon,lept)"
          Histo(15)%NBins  = 20
          Histo(15)%BinSize= 0.2d0
          Histo(15)%LowVal = 0d0
          Histo(15)%SetScale= 1d0

ELSEIF( ObsSet.EQ.28 ) THEN! set of observables for ttbgamma production semi-lept. decays at the LHC for Q_top measurement
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.4 .and. TopDecays.ne.3) call Error("TopDecays needs to be 3 (for Qt=-4/3) or 4 (for Qt=Qup)")
          NumHistograms = 15
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_ATop"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_Top"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 50d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 40
          Histo(4)%BinSize= 0.25d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_Photon"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 25d0*GeV
          Histo(7)%LowVal = 0d0*GeV
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "eta_Photon"
          Histo(8)%NBins  = 40
          Histo(8)%BinSize= 0.25d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "pT_LepP"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 25d0*GeV
          Histo(9)%LowVal =  0d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "eta_LepP"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 0.25d0
          Histo(10)%LowVal =-5.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "pT_miss"
          Histo(11)%NBins  = 40
          Histo(11)%BinSize= 25d0*GeV
          Histo(11)%LowVal =  0d0*GeV
          Histo(11)%SetScale= 100d0

          Histo(12)%Info   = "HT(jets+lept+pho)"
          Histo(12)%NBins  = 40
          Histo(12)%BinSize= 50d0*GeV
          Histo(12)%LowVal = 150d0*GeV
          Histo(12)%SetScale= 100d0

          Histo(13)%Info   = "R(pho,bjet)"
          Histo(13)%NBins  = 25
          Histo(13)%BinSize= 0.25d0
          Histo(13)%LowVal = 0d0
          Histo(13)%SetScale= 1d0

          Histo(14)%Info   = "m(lep+bjet)"
          Histo(14)%NBins  = 90
          Histo(14)%BinSize= 5d0*GeV
          Histo(14)%LowVal = 20d0*GeV
          Histo(14)%SetScale= 100d0

          Histo(15)%Info   = "phi(photon,lept)"
          Histo(15)%NBins  = 20
          Histo(15)%BinSize= 0.2d0
          Histo(15)%LowVal = 0d0
          Histo(15)%SetScale= 1d0


ELSEIF( ObsSet.EQ.29 ) THEN! set of observables for ttbgamma production without decays at Tevatron
          NumHistograms = 8
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 25d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 25d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "eta_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "pT_Photon"
          Histo(5)%NBins  = 40
          Histo(5)%BinSize= 20d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "eta_Photon"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 0.2d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "etaFB_ATop"
          Histo(7)%NBins  = 2
          Histo(7)%BinSize= 5d0
          Histo(7)%LowVal =-5.0d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "etaFB_Top"
          Histo(8)%NBins  = 2
          Histo(8)%BinSize= 5d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0



ELSEIF( ObsSet.EQ.31 ) THEN! set of observables for HTHTbar + A0 production
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(XTopDecays.ne.0)  call Error("XTopDecays needs to be 0")
          NumHistograms = 1
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 25d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0



ELSEIF( ObsSet.EQ.32 .OR. ObsSet.EQ.34 ) THEN! set of observables for HTHTbar + A0/BH production
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.1  ) call Error("TopDecays needs to be 1!")
          if(XTopDecays.ne.1 .and. XTopDecays.ne.2 ) call Error("XTopDecays needs to be 1(BH) or 2(A0)!")
          NumHistograms = 7
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_LepP"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 20d0*GeV
          Histo(1)%LowVal = 20d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_LepP"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "ET_miss"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 20d0*GeV
          Histo(3)%LowVal = 20d0*GeV
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "phi(l+,l-)"
          Histo(4)%NBins  = 15
          Histo(4)%BinSize= 0.25d0
          Histo(4)%LowVal = 0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "m_l+l-"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 20d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "MT(W)"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 20d0*GeV
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "MT^eff"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 20d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0


ELSEIF( ObsSet.EQ.33 .OR. ObsSet.EQ.35 ) THEN! set of observables for HTHTbar + A0/BH production
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.4  ) call Error("TopDecays needs to be 4!")
          if(XTopDecays.ne.1 ) call Error("XTopDecays needs to be 1!")
          NumHistograms = 11
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif


          Histo(1)%Info   = "pT_LepP"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 40d0*GeV
          Histo(1)%LowVal = 20d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_LepP"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "ET_miss"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 40d0*GeV
          Histo(3)%LowVal = 20d0*GeV
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "HT"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 50d0*GeV
          Histo(4)%LowVal = 0d0
          Histo(4)%SetScale= 100d0

          Histo(5)%Info   = "pT(softest jet)"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 10d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "MT(W)"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 40d0*GeV
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "MT^eff"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 50d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "pT_Top"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 50d0*GeV
          Histo(8)%LowVal = 0d0*GeV
          Histo(8)%SetScale= 100d0

          Histo(9)%Info   = "eta_Top"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 0.25d0
          Histo(9)%LowVal =-5.0d0
          Histo(9)%SetScale= 1d0

          Histo(10)%Info   = "Log(S)"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 1.0d0
          Histo(10)%LowVal =-10.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "r_pT"
          Histo(11)%NBins  = 20
          Histo(11)%BinSize= 0.2d0
          Histo(11)%LowVal =-2.0d0
          Histo(11)%SetScale= 1d0


ELSEIF( ObsSet.EQ.41 ) THEN! set of observables for STSTbar (stable stops)
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(XTopDecays.ne.0)  call Error("XTopDecays needs to be 0")
          NumHistograms = 1
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal =  0d0*GeV
          Histo(1)%SetScale= 100d0



ELSEIF( ObsSet.EQ.42 .OR. ObsSet.EQ.44 ) THEN! set of observables for STSTbar + Chi production (di-lept. tops)

          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.1  ) call Error("TopDecays needs to be 1!")
          if(XTopDecays.ne.3  ) call Error("XTopDecays needs to be 3!")
          NumHistograms = 8
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_LepP"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 40d0*GeV
          Histo(1)%LowVal = 20d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_LepP"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "ET_miss"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 40d0*GeV
          Histo(3)%LowVal = 40d0*GeV
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "phi(l+,l-)"
          Histo(4)%NBins  = 15
          Histo(4)%BinSize= 0.25d0
          Histo(4)%LowVal = 0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "m_l+l-"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 20d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "HT"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 50d0*GeV
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "MT(W)"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 40d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "MT^eff"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 100d0*GeV
          Histo(8)%LowVal = 0d0
          Histo(8)%SetScale= 100d0


ELSEIF( ObsSet.EQ.43 .OR. ObsSet.EQ.45 ) THEN! set of observables for STSTbar + Chi production
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.4  ) call Error("TopDecays needs to be 4!")
          if(XTopDecays.ne.3 ) call Error("XTopDecays needs to be 3!")
          NumHistograms = 11
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif



          Histo(1)%Info   = "pT_LepP"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 40d0*GeV
          Histo(1)%LowVal = 20d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_LepP"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "ET_miss"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 40d0*GeV
          Histo(3)%LowVal = 20d0*GeV
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "HT"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 50d0*GeV
          Histo(4)%LowVal = 0d0
          Histo(4)%SetScale= 100d0

          Histo(5)%Info   = "pT(softest jet)"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 10d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "MT(W)"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 40d0*GeV
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "MT^eff"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 50d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "pT_Top"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 50d0*GeV
          Histo(8)%LowVal = 0d0*GeV
          Histo(8)%SetScale= 100d0

          Histo(9)%Info   = "eta_Top"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 0.25d0
          Histo(9)%LowVal =-5.0d0
          Histo(9)%SetScale= 1d0

          Histo(10)%Info   = "Log(S)"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 1.0d0
          Histo(10)%LowVal =-10.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "r_pT"
          Histo(11)%NBins  = 20
          Histo(11)%BinSize= 0.2d0
          Histo(11)%LowVal =-2.0d0
          Histo(11)%SetScale= 1d0

   

ELSEIF( ObsSet.EQ.48 ) THEN! set of observables for STOP width
          if(TopDecays.ne.0  ) call Error("TopDecays needs to be 0!")
          if(XTopDecays.ne.3 .and. XTopDecays.ne.1) call Error("XTopDecays needs to be 1 or 3!")
!           NumHistograms = 0
!           if( .not.allocated(Histo) ) then
!                 allocate( Histo(1:NumHistograms), stat=AllocStatus  )
!                 if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
!           endif



ELSEIF( ObsSet.EQ.51 ) THEN! set of observables for ttb+Z (stable tops)
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.0)  call Error("TopDecays needs to be 0")
          NumHistograms = 1
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT(top)"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal =  0d0*GeV
          Histo(1)%SetScale= 100d0





ELSEIF( ObsSet.EQ.52 .or. ObsSet.EQ.55 ) THEN! set of observables for ttb+Z ( di-lept. ttbar decays and di-lept. Z decay )
          if(abs(TopDecays).ne.1)  call Error("TopDecays needs to be 1")
          if(abs(ZDecays).ne.1)    call Error("ZDecays needs to be 1")
          NumHistograms = 23
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT(lep+)"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 20d0*GeV
          Histo(1)%LowVal =  0d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT(lep-))"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 20d0*GeV
          Histo(2)%LowVal =  0d0*GeV
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "pT(Z)"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 20d0*GeV
          Histo(3)%LowVal =  0d0*GeV
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "pT(top))"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 20d0*GeV
          Histo(4)%LowVal =  0d0*GeV
          Histo(4)%SetScale= 100d0

          Histo(5)%Info   = "pT(atop))"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 20d0*GeV
          Histo(5)%LowVal =  0d0*GeV
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "pT(j1))"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 10d0*GeV
          Histo(6)%LowVal =  0d0*GeV
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "pT(j2))"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 10d0*GeV
          Histo(7)%LowVal =  0d0*GeV
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "pT(mu+))"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 10d0*GeV
          Histo(8)%LowVal =  0d0*GeV
          Histo(8)%SetScale= 100d0

          Histo(9)%Info   = "pT(e-))"
          Histo(9)%NBins  = 50
          Histo(9)%BinSize= 10d0*GeV
          Histo(9)%LowVal =  0d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "pT(miss))"
          Histo(10)%NBins  = 50
          Histo(10)%BinSize= 10d0*GeV
          Histo(10)%LowVal =  0d0*GeV
          Histo(10)%SetScale= 100d0

          Histo(11)%Info   = "eta(lep+)"
          Histo(11)%NBins  = 50
          Histo(11)%BinSize= 0.5d0
          Histo(11)%LowVal =-5.0d0
          Histo(11)%SetScale= 1d0

          Histo(12)%Info   = "eta(lep-)"
          Histo(12)%NBins  = 50
          Histo(12)%BinSize= 0.5d0
          Histo(12)%LowVal =-5.0d0
          Histo(12)%SetScale= 1d0

          Histo(13)%Info   = "eta(Z)"
          Histo(13)%NBins  = 30
          Histo(13)%BinSize= 0.2d0
          Histo(13)%LowVal =-3.0d0
          Histo(13)%SetScale= 1d0

          Histo(14)%Info   = "eta(top)"
          Histo(14)%NBins  = 30
          Histo(14)%BinSize= 0.2d0
          Histo(14)%LowVal =-3.0d0
          Histo(14)%SetScale= 1d0

          Histo(15)%Info   = "eta(atop)"
          Histo(15)%NBins  = 30
          Histo(15)%BinSize= 0.2d0
          Histo(15)%LowVal =-3.0d0
          Histo(15)%SetScale= 1d0

          Histo(16)%Info   = "eta(j1)"
          Histo(16)%NBins  = 50
          Histo(16)%BinSize= 0.5d0
          Histo(16)%LowVal =-5.0d0
          Histo(16)%SetScale= 1d0

          Histo(17)%Info   = "eta(j2)"
          Histo(17)%NBins  = 50
          Histo(17)%BinSize= 0.5d0
          Histo(17)%LowVal =-5.0d0
          Histo(17)%SetScale= 1d0

          Histo(18)%Info   = "eta(mu+)"
          Histo(18)%NBins  = 50
          Histo(18)%BinSize= 0.5d0
          Histo(18)%LowVal =-5.0d0
          Histo(18)%SetScale= 1d0

          Histo(19)%Info   = "eta(e-)"
          Histo(19)%NBins  = 50
          Histo(19)%BinSize= 0.5d0
          Histo(19)%LowVal =-5.0d0
          Histo(19)%SetScale= 1d0

          Histo(20)%Info   = "phi(l+,l-)"
          Histo(20)%NBins  = 15    *4d0
          Histo(20)%BinSize= 0.25d0/4d0
          Histo(20)%LowVal = 0d0
          Histo(20)%SetScale= 1d0

          Histo(21)%Info   = "pseudo(Z)"
          Histo(21)%NBins  = 30
          Histo(21)%BinSize= 0.2d0
          Histo(21)%LowVal =-3.0d0
          Histo(21)%SetScale= 1d0

          Histo(22)%Info   = "pseudo(top)"
          Histo(22)%NBins  = 30
          Histo(22)%BinSize= 0.2d0
          Histo(22)%LowVal =-3.0d0
          Histo(22)%SetScale= 1d0

          Histo(23)%Info   = "pseudo(atop)"
          Histo(23)%NBins  = 30
          Histo(23)%BinSize= 0.2d0
          Histo(23)%LowVal =-3.0d0
          Histo(23)%SetScale= 1d0


ELSEIF( ObsSet.EQ.53 .or. ObsSet.EQ.56 ) THEN! set of observables for ttb+Z ( di-lept. ttbar decays and di-lept. Z decay )
          if(abs(TopDecays).ne.4)  call Error("TopDecays needs to be 4")
          if(abs(ZDecays).ne.1)    call Error("ZDecays needs to be 1")
!          NumHistograms = 3
          NumHistograms = 48
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif
          print *, 'allocating histos'

          Histo(1)%Info   = "pT(lep+)"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 20d0*GeV
          Histo(1)%LowVal =  0d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT(mu-)"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 20d0*GeV
          Histo(2)%LowVal =  0d0*GeV
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "pT(mu+)"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 20d0*GeV
          Histo(3)%LowVal =  0d0*GeV
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "pT(b1)"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 20d0*GeV
          Histo(4)%LowVal =  0d0*GeV
          Histo(4)%SetScale= 100d0

          Histo(5)%Info   = "pT(b2)"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 20d0*GeV
          Histo(5)%LowVal =  0d0*GeV
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "pT(j1)"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 20d0*GeV
          Histo(6)%LowVal =  0d0*GeV
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "pT(j2)"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 20d0*GeV
          Histo(7)%LowVal =  0d0*GeV
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "pT(miss)"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 20d0*GeV
          Histo(8)%LowVal =  0d0*GeV
          Histo(8)%SetScale= 100d0

          Histo(9)%Info   = "y(lep+)"
          Histo(9)%NBins  = 50
          Histo(9)%BinSize= 0.5d0
          Histo(9)%LowVal =-5.0d0
          Histo(9)%SetScale= 1d0

          Histo(10)%Info   = "y(mu-)"
          Histo(10)%NBins  = 50
          Histo(10)%BinSize= 0.5d0
          Histo(10)%LowVal =-5.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "y(mu+)"
          Histo(11)%NBins  = 50
          Histo(11)%BinSize= 0.5d0
          Histo(11)%LowVal =-5.0d0
          Histo(11)%SetScale= 1d0

          Histo(12)%Info   = "y(b1)"
          Histo(12)%NBins  = 50
          Histo(12)%BinSize= 0.5d0
          Histo(12)%LowVal =-5.0d0
          Histo(12)%SetScale= 1d0

          Histo(13)%Info   = "y(b2)"
          Histo(13)%NBins  = 50
          Histo(13)%BinSize= 0.5d0
          Histo(13)%LowVal =-5.0d0
          Histo(13)%SetScale= 1d0

          Histo(14)%Info   = "y(j1)"
          Histo(14)%NBins  = 50
          Histo(14)%BinSize= 0.5d0
          Histo(14)%LowVal =-5.0d0
          Histo(14)%SetScale= 1d0

          Histo(15)%Info   = "y(j2)"
          Histo(15)%NBins  = 50
          Histo(15)%BinSize= 0.5d0
          Histo(15)%LowVal =-5.0d0
          Histo(15)%SetScale= 1d0

          Histo(16)%Info   = "pT(Z)"
          Histo(16)%NBins  = 50
          Histo(16)%BinSize= 20d0*GeV
          Histo(16)%LowVal =  0d0*GeV
          Histo(16)%SetScale= 100d0

          Histo(17)%Info   = "y(Z)"
          Histo(17)%NBins  = 50
          Histo(17)%BinSize= 0.5d0
          Histo(17)%LowVal =-5.0d0
          Histo(17)%SetScale= 1d0

          Histo(18)%Info   = "pT(top)"
          Histo(18)%NBins  = 50
          Histo(18)%BinSize= 20d0*GeV
          Histo(18)%LowVal =  0d0*GeV
          Histo(18)%SetScale= 100d0

          Histo(19)%Info   = "y(top)"
          Histo(19)%NBins  = 50
          Histo(19)%BinSize= 0.5d0
          Histo(19)%LowVal =-5.0d0
          Histo(19)%SetScale= 1d0

          Histo(20)%Info   = "pT(antitop)"
          Histo(20)%NBins  = 50
          Histo(20)%BinSize= 20d0*GeV
          Histo(20)%LowVal =  0d0*GeV
          Histo(20)%SetScale= 100d0

          Histo(21)%Info   = "y(antitop)"
          Histo(21)%NBins  = 50
          Histo(21)%BinSize= 0.5d0
          Histo(21)%LowVal =-5.0d0
          Histo(21)%SetScale= 1d0

          Histo(22)%Info   = "mT_inv(lb,miss)"
          Histo(22)%NBins  = 50
          Histo(22)%BinSize= 20d0*GeV
          Histo(22)%LowVal =  0d0*GeV
          Histo(22)%SetScale= 100d0

          Histo(23)%Info   = "phi(Z,t)"
          Histo(23)%NBins  = 15    *4d0
          Histo(23)%BinSize= 0.25d0/4d0
          Histo(23)%LowVal = 0d0
          Histo(23)%SetScale= 1d0

          Histo(24)%Info   = "phi(Z,antit)"
          Histo(24)%NBins  = 15    *4d0
          Histo(24)%BinSize= 0.25d0/4d0
          Histo(24)%LowVal = 0d0
          Histo(24)%SetScale= 1d0

          Histo(25)%Info   = "phi(t,antit)"
          Histo(25)%NBins  = 15    *4d0
          Histo(25)%BinSize= 0.25d0/4d0
          Histo(25)%LowVal = 0d0
          Histo(25)%SetScale= 1d0

          Histo(26)%Info   = "phi(mu-,mu+)"
          Histo(26)%NBins  = 15    *4d0
          Histo(26)%BinSize= 0.25d0/4d0
          Histo(26)%LowVal = 0d0
          Histo(26)%SetScale= 1d0

          Histo(27)%Info   = "phi(mu-,l+)"
          Histo(27)%NBins  = 15    *4d0
          Histo(27)%BinSize= 0.25d0/4d0
          Histo(27)%LowVal = 0d0
          Histo(27)%SetScale= 1d0

          Histo(28)%Info   = "phi(mu-,b1)"
          Histo(28)%NBins  = 15    *4d0
          Histo(28)%BinSize= 0.25d0/4d0
          Histo(28)%LowVal = 0d0
          Histo(28)%SetScale= 1d0

          Histo(29)%Info   = "phi(mu-,b2)"
          Histo(29)%NBins  = 15    *4d0
          Histo(29)%BinSize= 0.25d0/4d0
          Histo(29)%LowVal = 0d0
          Histo(29)%SetScale= 1d0

          Histo(30)%Info   = "phi(mu-,j1)"
          Histo(30)%NBins  = 15    *4d0
          Histo(30)%BinSize= 0.25d0/4d0
          Histo(30)%LowVal = 0d0
          Histo(30)%SetScale= 1d0

          Histo(31)%Info   = "phi(mu-,j2)"
          Histo(31)%NBins  = 15    *4d0
          Histo(31)%BinSize= 0.25d0/4d0
          Histo(31)%LowVal = 0d0
          Histo(31)%SetScale= 1d0          
          
          Histo(32)%Info   = "phi(mu+,l+)"
          Histo(32)%NBins  = 15    *4d0
          Histo(32)%BinSize= 0.25d0/4d0
          Histo(32)%LowVal = 0d0
          Histo(32)%SetScale= 1d0

          Histo(33)%Info   = "phi(mu+,b1)"
          Histo(33)%NBins  = 15    *4d0
          Histo(33)%BinSize= 0.25d0/4d0
          Histo(33)%LowVal = 0d0
          Histo(33)%SetScale= 1d0

          Histo(34)%Info   = "phi(mu+,b2)"
          Histo(34)%NBins  = 15    *4d0
          Histo(34)%BinSize= 0.25d0/4d0
          Histo(34)%LowVal = 0d0
          Histo(34)%SetScale= 1d0

          Histo(35)%Info   = "phi(mu+,j1)"
          Histo(35)%NBins  = 15    *4d0
          Histo(35)%BinSize= 0.25d0/4d0
          Histo(35)%LowVal = 0d0
          Histo(35)%SetScale= 1d0

          Histo(36)%Info   = "phi(mu+,j2)"
          Histo(36)%NBins  = 15    *4d0
          Histo(36)%BinSize= 0.25d0/4d0
          Histo(36)%LowVal = 0d0
          Histo(36)%SetScale= 1d0

          Histo(37)%Info   = "DeltaR(Z,t)"
          Histo(37)%NBins  = 15    *4d0
          Histo(37)%BinSize= 0.25d0/4d0
          Histo(37)%LowVal = 0d0
          Histo(37)%SetScale= 1d0

          Histo(38)%Info   = "DeltaR(Z,antit)"
          Histo(38)%NBins  = 15    *4d0
          Histo(38)%BinSize= 0.25d0/4d0
          Histo(38)%LowVal = 0d0
          Histo(38)%SetScale= 1d0

          Histo(39)%Info   = "DeltaR(t,antit)"
          Histo(39)%NBins  = 15    *4d0
          Histo(39)%BinSize= 0.25d0/4d0
          Histo(39)%LowVal = 0d0
          Histo(39)%SetScale= 1d0

          Histo(40)%Info   = "DeltaR(mu-,mu+)"
          Histo(40)%NBins  = 15    *4d0
          Histo(40)%BinSize= 0.25d0/4d0
          Histo(40)%LowVal = 0d0
          Histo(40)%SetScale= 1d0

          Histo(41)%Info   = "DeltaR(mu-,l+)"
          Histo(41)%NBins  = 15    *4d0
          Histo(41)%BinSize= 0.25d0/4d0
          Histo(41)%LowVal = 0d0
          Histo(41)%SetScale= 1d0

          Histo(42)%Info   = "DeltaR(mu-,b1)"
          Histo(42)%NBins  = 15    *4d0
          Histo(42)%BinSize= 0.25d0/4d0
          Histo(42)%LowVal = 0d0
          Histo(42)%SetScale= 1d0

          Histo(43)%Info   = "DeltaR(mu-,b2)"
          Histo(43)%NBins  = 15    *4d0
          Histo(43)%BinSize= 0.25d0/4d0
          Histo(43)%LowVal = 0d0
          Histo(43)%SetScale= 1d0

          Histo(44)%Info   = "DeltaR(mu-,j1)"
          Histo(44)%NBins  = 15    *4d0
          Histo(44)%BinSize= 0.25d0/4d0
          Histo(44)%LowVal = 0d0
          Histo(44)%SetScale= 1d0

          Histo(45)%Info   = "DeltaR(mu-,j2)"
          Histo(45)%NBins  = 15    *4d0
          Histo(45)%BinSize= 0.25d0/4d0
          Histo(45)%LowVal = 0d0
          Histo(45)%SetScale= 1d0          
          
          Histo(46)%Info   = "DeltaR(mu+,l+)"
          Histo(46)%NBins  = 15    *4d0
          Histo(46)%BinSize= 0.25d0/4d0
          Histo(46)%LowVal = 0d0
          Histo(46)%SetScale= 1d0

          Histo(47)%Info   = "DeltaR(mu+,b1)"
          Histo(47)%NBins  = 15    *4d0
          Histo(47)%BinSize= 0.25d0/4d0
          Histo(47)%LowVal = 0d0
          Histo(47)%SetScale= 1d0

          Histo(48)%Info   = "DeltaR(mu+,b2)"
          Histo(48)%NBins  = 15    *4d0
          Histo(48)%BinSize= 0.25d0/4d0
          Histo(48)%LowVal = 0d0
          Histo(48)%SetScale= 1d0

          Histo(49)%Info   = "DeltaR(j1,j2)"
          Histo(49)%NBins  = 15    *4d0
          Histo(49)%BinSize= 0.25d0/4d0
          Histo(49)%LowVal = 0d0
          Histo(49)%SetScale= 1d0

          Histo(50)%Info   = "DeltaR(j1,j3)"
          Histo(50)%NBins  = 15    *4d0
          Histo(50)%BinSize= 0.25d0/4d0
          Histo(50)%LowVal = 0d0
          Histo(50)%SetScale= 1d0

          Histo(51)%Info   = "DeltaR(j1,j4)"
          Histo(51)%NBins  = 15    *4d0
          Histo(51)%BinSize= 0.25d0/4d0
          Histo(51)%LowVal = 0d0
          Histo(51)%SetScale= 1d0

          Histo(52)%Info   = "DeltaR(j2,j3)"
          Histo(52)%NBins  = 15    *4d0
          Histo(52)%BinSize= 0.25d0/4d0
          Histo(52)%LowVal = 0d0
          Histo(52)%SetScale= 1d0

          Histo(53)%Info   = "DeltaR(j2,j4)"
          Histo(53)%NBins  = 15    *4d0
          Histo(53)%BinSize= 0.25d0/4d0
          Histo(53)%LowVal = 0d0
          Histo(53)%SetScale= 1d0

          Histo(54)%Info   = "DeltaR(j3,j4)"
          Histo(54)%NBins  = 15    *4d0
          Histo(54)%BinSize= 0.25d0/4d0
          Histo(54)%LowVal = 0d0
          Histo(54)%SetScale= 1d0
!
!          Histo(49)%Info   = "DeltaR(mu+,j1)"
!          Histo(49)%NBins  = 15    *4d0
!          Histo(49)%BinSize= 0.25d0/4d0
!          Histo(49)%LowVal = 0d0
!          Histo(49)%SetScale= 1d0
!
!          Histo(50)%Info   = "DeltaR(mu+,j2)"
!          Histo(50)%NBins  = 15    *4d0
!          Histo(50)%BinSize= 0.25d0/4d0
!          Histo(50)%LowVal = 0d0
!          Histo(50)%SetScale= 1d0
!          
!          
!          print *, 'done allocating histos'


          

!          Histo(1)%Info   = "pT(lep+)"
!          Histo(1)%NBins  = 50
!          Histo(1)%BinSize= 20d0*GeV
!          Histo(1)%LowVal =  0d0*GeV
!          Histo(1)%SetScale= 100d0
!
!          Histo(2)%Info   = "pT(Z(l+,l-))"
!          Histo(2)%NBins  = 50
!          Histo(2)%BinSize= 20d0*GeV
!          Histo(2)%LowVal =  0d0*GeV
!          Histo(2)%SetScale= 100d0
!
!          Histo(3)%Info   = "phi(l+,l-)"
!          Histo(3)%NBins  = 15    *4d0
!          Histo(3)%BinSize= 0.25d0/4d0
!          Histo(3)%LowVal = 0d0
!          Histo(3)%SetScale= 1d0



ELSEIF( ObsSet.EQ.60  ) THEN! set of observables for Zprime, stable tops
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.0  ) call Error("TopDecays needs to be 0!")
          NumHistograms = 9
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif


          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 50d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "y_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.5d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "y_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.5d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "M_TTbar"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 50d0*GeV
          Histo(5)%LowVal = 350d0*GeV
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "deltaPhi_TTbar"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 0.1d0
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "CosTheta_scatter"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 0.1d0
          Histo(7)%LowVal = -1d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "CosTheta_star"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 0.1d0
          Histo(8)%LowVal = -1d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "M_TTbar+jet"
          Histo(9)%NBins  = 50
          Histo(9)%BinSize= 40d0*GeV
          Histo(9)%LowVal = 350d0*GeV
          Histo(9)%SetScale= 100d0

ELSEIF( ObsSet.EQ.61 ) THEN! set of observables for Zprime, top decaying to dileptons
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.1  ) call Error("TopDecays needs to be 1!")
          NumHistograms = 13
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_lep1"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 20d0*GeV
          Histo(1)%LowVal = 20d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_bjet1"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 20d0*GeV
          Histo(2)%LowVal = 30d0*GeV
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "pT_lep2"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 20d0*GeV
          Histo(3)%LowVal = 20d0*GeV
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "pT_bjet2"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 20d0*GeV
          Histo(4)%LowVal = 30d0*GeV
          Histo(4)%SetScale= 100d0

          Histo(5)%Info   = "M_TTbar"
          Histo(5)%NBins  = 20
          Histo(5)%BinSize= 50d0*GeV
          Histo(5)%LowVal = 1025d0*GeV
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "M_eff"
          Histo(6)%NBins  = 20
          Histo(6)%BinSize= 50d0*GeV
          Histo(6)%LowVal = 1025d0*GeV
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "y_lep1"
          Histo(7)%NBins  = 30
          Histo(7)%BinSize= 0.2d0
          Histo(7)%LowVal =-3.0d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "y_bjet1"
          Histo(8)%NBins  = 30
          Histo(8)%BinSize= 0.2d0
          Histo(8)%LowVal =-3.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "y_lep2"
          Histo(9)%NBins  = 30
          Histo(9)%BinSize= 0.2d0
          Histo(9)%LowVal =-3.0d0
          Histo(9)%SetScale= 1d0

          Histo(10)%Info   = "y_bjet2"
          Histo(10)%NBins  = 30
          Histo(10)%BinSize= 0.2d0
          Histo(10)%LowVal =-3.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "deltaPhi_LL"
          Histo(11)%NBins  = 50
          Histo(11)%BinSize= 0.08d0
          Histo(11)%LowVal = 0d0
          Histo(11)%SetScale= 1d0

          Histo(12)%Info   = "dPhiMinus"
          Histo(12)%NBins  = 50
          Histo(12)%BinSize= 0.126d0
          Histo(12)%LowVal = -DblPi
          Histo(12)%SetScale= 1d0

          Histo(13)%Info   = "dPhiPlus"
          Histo(13)%NBins  = 50
          Histo(13)%BinSize= 0.063d0
          Histo(13)%LowVal = 0
          Histo(13)%SetScale= 1d0


ELSEIF( ObsSet.EQ.62 ) THEN! set of observables for Zprime, fully hadronic top decay
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.2  ) call Error("TopDecays needs to be 2!")
          NumHistograms = 8
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif


          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 50d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "y_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "y_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "M_TTbar"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 40d0*GeV
          Histo(5)%LowVal = 350d0*GeV
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "deltaPhi_TTbar"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 0.08d0
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "CosTheta_scatter"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 0.04d0
          Histo(7)%LowVal = -1d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "CosTheta_star"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 0.04d0
          Histo(8)%LowVal = -1d0
          Histo(8)%SetScale= 1d0


ELSEIF( ObsSet.EQ.64 .OR. ObsSet.EQ.65 ) THEN ! Zprime, semi-hadronic top decay (for ATLAS analysis: James Ferrando)
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.4 .and. TopDecays.ne.3 ) call Error("TopDecays needs to be 3 or 4!")
          NumHistograms = 8
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif


          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 50d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "y_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "y_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "M_TTbar"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 40d0*GeV
          Histo(5)%LowVal = 350d0*GeV
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "deltaPhi_TTbar"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 0.08d0
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "CosTheta_scatter"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 0.04d0
          Histo(7)%LowVal = -1d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "CosTheta_star"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 0.04d0
          Histo(8)%LowVal = -1d0
          Histo(8)%SetScale= 1d0



ELSEIF( ObsSet.EQ.66 ) THEN ! Zprime, semi-hadronic top decay (for CMS analysis: Roman Kogler)
          if(Collider.ne.12)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.4  ) call Error("TopDecays needs to be 4!")
          NumHistograms = 8
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif


          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 50d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "y_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "y_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "M_TTbar"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 40d0*GeV
          Histo(5)%LowVal = 350d0*GeV
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "deltaPhi_TTbar"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 0.08d0
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "CosTheta_scatter"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 0.04d0
          Histo(7)%LowVal = -1d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "CosTheta_star"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 0.04d0
          Histo(8)%LowVal = -1d0
          Histo(8)%SetScale= 1d0



ELSE
    call Error("This ObsSet is not available ",ObsSet)
ENDIF



do NHisto=1,NumHistograms
      if( .not.allocated(Histo(NHisto)%Value) ) then
        allocate( Histo(NHisto)%Value(0:Histo(NHisto)%NBins+1), stat=AllocStatus  )
        if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
      endif
      if( .not.allocated(Histo(NHisto)%Value2) ) then
        allocate( Histo(NHisto)%Value2(0:Histo(NHisto)%NBins+1), stat=AllocStatus  )
        if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
      endif
      if( .not.allocated(Histo(NHisto)%Hits) ) then
        allocate( Histo(NHisto)%Hits(0:Histo(NHisto)%NBins+1), stat=AllocStatus  )
        if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
      endif
      Histo(NHisto)%Value(0:Histo(NHisto)%NBins+1) = 0d0
      Histo(NHisto)%Value2(0:Histo(NHisto)%NBins+1)= 0d0
      Histo(NHisto)%Hits(0:Histo(NHisto)%NBins+1)  = 0
enddo


RETURN
END SUBROUTINE





SUBROUTINE EvalPhasespace_HTopDK(Topol,HTopMom,xRndPS,MomDK,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
integer :: Topol
real(8) :: PSWgt
real(8) :: HTopMom(1:4)
real(8) :: MomDK(:,:)
real(8) :: xRndPS(:)


   if( Topol.eq.HT_A0_T ) then
      call EvalPhasespace_HTopDecay(HTopMom,xRndPS,m_A0,.false.,MomDK,PSWgt)
   elseif( Topol.eq.HT_BH_T ) then
      call EvalPhasespace_HTopDecay(HTopMom,xRndPS,m_BH,.false.,MomDK,PSWgt)
   elseif( Topol.eq.HT_BH_T_G ) then
      call EvalPhasespace_HTopDecay(HTopMom,xRndPS,m_BH,.true.,MomDK,PSWgt)
   else
      call Error("EvalPS not yet implemented")
   endif


RETURN
END SUBROUTINE



SUBROUTINE EvalPhasespace_HTopDecay(HTopMom,xRndPS,Mass,GluonRad,MomDK,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
logical GluonRad
real(8) :: PSWgt,PSWgt2,Mass
real(8) :: HTopMom(1:4),A0Mom(1:4)
real(8) :: MomDK(:,:)
real(8) :: xRndPS(:)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)
real(8),parameter :: N3=3, PiWgt3 = (2d0*Pi)**(4-N3*3) * (4d0*Pi)**(N3-1)
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,soft,coll
logical,save :: flip=.true.

    if( .not.GluonRad ) then!  no extra gluon radiation
      call genps(2,m_HTop,xRndPS(1:2),(/Mass,m_SMTop/),MomDK(1:4,1:2),PSWgt2)! Htop decay
!     boost all guys to the top frame:
      call boost(MomDK(1:4,1),HTopMom(1:4),m_HTop)
      call boost(MomDK(1:4,2),HTopMom(1:4),m_HTop)
      PSWgt = PSWgt2*PiWgt2

    else! extra gluon emission

     call genps(3,m_HTop,xRndPS(1:5),(/Mass,m_SMTop,0d0/),MomDK(1:4,1:3),PSWgt2)! Htop decay with gluon

!          Pcol1= 5 -1
!          Pcol2= 5 -1
!          SingDepth = 1e-10
!          Steps = 20
!          call gensing(3,m_HTop,(/Mass,m_SMTop,0d0/),MomDK(1:4,1:3),Pcol1,Pcol2,SingDepth,Steps); print *, "running gensing"
!          PSWgt2=1d0

!     boost all guys to the top frame:
      call boost(MomDK(1:4,1),HTopMom(1:4),m_HTop)
      call boost(MomDK(1:4,2),HTopMom(1:4),m_HTop)
      call boost(MomDK(1:4,3),HTopMom(1:4),m_HTop)
      PSWgt = PSWgt2*PiWgt3

    endif


RETURN
END SUBROUTINE





SUBROUTINE EvalPhasespace_STopDK(Topol,STopMom,xRndPS,MomDK,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
integer :: Topol
real(8) :: PSWgt
real(8) :: STopMom(1:4)
real(8) :: MomDK(:,:)
real(8) :: xRndPS(:)


   if( Topol.eq.ST_Chi0_T ) then
      call EvalPhasespace_STopDecay(STopMom,xRndPS,.false.,MomDK,PSWgt)

   elseif( Topol.eq.ST_Chi0_T_G ) then
      call EvalPhasespace_STopDecay(STopMom,xRndPS,.true.,MomDK,PSWgt)
   else
      call Error("EvalPS not yet implemented")
   endif


RETURN
END SUBROUTINE



SUBROUTINE EvalPhasespace_STopDecay(STopMom,xRndPS,GluonRad,MomDK,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
logical GluonRad
real(8) :: PSWgt,PSWgt2
real(8) :: STopMom(1:4),ChiMom(1:4)
real(8) :: MomDK(:,:)
real(8) :: xRndPS(:)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)
real(8),parameter :: N3=3, PiWgt3 = (2d0*Pi)**(4-N3*3) * (4d0*Pi)**(N3-1)
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,soft
logical,save :: flip=.true.

    if( .not.GluonRad ) then!  no extra gluon radiation
!     MomDK(1:4,i): i= 1:Chi, 2:top
      call genps(2,m_STop,xRndPS(1:2),(/m_Chi,m_Top/),MomDK(1:4,1:2),PSWgt2)! Stop decay
!     boost all guys to the top frame:
      call boost(MomDK(1:4,1),STopMom(1:4),m_STop)
      call boost(MomDK(1:4,2),STopMom(1:4),m_STop)
      PSWgt = PSWgt2*PiWgt2

    else! extra gluon emission
!     MomDK(1:4,i): i= 1:Chi, 2:top, 3:gluon
      call genps(3,m_STop,xRndPS(1:5),(/m_Chi,m_Top,0d0/),MomDK(1:4,1:3),PSWgt2)! Stop decay with gluon

!          Pcol1= 5 -1
!          Pcol2= 5 -1
!          SingDepth = 1e-10
!          Steps = 20
!          call gensing(3,m_STop,(/m_Chi,m_Top,0d0/),MomDK(1:4,1:3),Pcol1,Pcol2,SingDepth,Steps); print *, "running gensing"
!          PSWgt2=1d0

!     boost all guys to the top frame:
      call boost(MomDK(1:4,1),STopMom(1:4),m_STop)
      call boost(MomDK(1:4,2),STopMom(1:4),m_STop)
      call boost(MomDK(1:4,3),STopMom(1:4),m_STop)
      PSWgt = PSWgt2*PiWgt3

      soft = dabs(MomDK(1,3))/m_STop
      if( soft.lt.1d-6  ) PSWgt = 0d0

    endif


RETURN
END SUBROUTINE







SUBROUTINE EvalPhasespace_TopDK(Topol,TopMom,xRndPS,MomDK,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
integer :: Topol
real(8) :: PSWgt
real(8) :: TopMom(1:4)
real(8) :: MomDK(:,:)
real(8) :: xRndPS(:)



   if( Topol.eq.T_B_W ) then
      call EvalPhasespace_TopDecay(TopMom,xRndPS,.false.,MomDK,PSWgt)
   elseif( Topol.eq.T_BG_W ) then
      call EvalPhasespace_TopDecay(TopMom,xRndPS,.true.,MomDK,PSWgt)
   elseif( Topol.eq.T_B_WG ) then
      call EvalPhasespace_TopDecay2(TopMom,xRndPS,.true.,MomDK,PSWgt)
   elseif( Topol.eq.T_BGG_W ) then
      call EvalPhasespace_TopDecay3(TopMom,xRndPS,MomDK,PSWgt)
   elseif( Topol.eq.T_BG_WG) then
      call EvalPhasespace_TopDecay4(TopMom,xRndPS,MomDK,PSWgt)
   elseif( Topol.eq.T_B_WGG) then
      call EvalPhasespace_TopDecay5(TopMom,xRndPS,MomDK,PSWgt)
   endif

RETURN
END SUBROUTINE






SUBROUTINE EvalPhasespace_TopDecay(TopMom,xRndPS,GluonRad,MomDK,PSWgt)!  top quark decay phase space with/without additional massless particle
use ModProcess
use ModMisc
use ModParameters
implicit none
logical GluonRad
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: TopMom(1:4),WMom(1:4),MomChk(1:4,1:3)
real(8) :: MomDK(:,:)
real(8) :: xRndPS(:)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)
real(8),parameter :: N3=3, PiWgt3 = (2d0*Pi)**(4-N3*3) * (4d0*Pi)**(N3-1)
real(8),parameter :: N4=4, PiWgt4 = (2d0*Pi)**(4-N4*3) * (4d0*Pi)**(N4-1)
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,soft,coll
logical,save :: flip=.true.

    if( .not.GluonRad ) then!  no extra gluon radiation
!     MomDK(1:4,i): i= 1:bottom, 2:lepton, 3:neutrino
      call genps(2,m_Top,xRndPS(1:2),(/0d0,m_W/),MomDK(1:4,1:2),PSWgt2)! top decay
      WMom(1:4) = MomDK(1:4,2)
      call genps(2,m_W,xRndPS(3:4),(/0d0,0d0/),MomDK(1:4,2:3),PSWgt3)! W decay
!     boost leptons to the W frame:
      call boost(MomDK(1:4,2),WMom(1:4),m_W)
      call boost(MomDK(1:4,3),WMom(1:4),m_W)
!     boost all guys to the top frame:
      call boost(MomDK(1:4,1),TopMom(1:4),m_Top)
      call boost(MomDK(1:4,2),TopMom(1:4),m_Top)
      call boost(MomDK(1:4,3),TopMom(1:4),m_Top)
      PSWgt = PSWgt2*PiWgt2 * PSWgt3*PiWgt2

    else! extra gluon emission
!     MomDK(1:4,i): i=1:bottom, 2:lepton, 3:neutrino, 4: gluon
!      call genps(3,m_Top,xRndPS(1:5),(/0d0,0d0,m_W/),MomDK(1:4,1:3),PSWgt2)! top decay with additional gluon
!      WMom(1:4) = MomDK(1:4,3)
!      MomDK(1:4,4) = MomDK(1:4,2)
!      PSWgt = PSWgt2*PiWgt3

      call yeti3(m_Top,xRndPS(1:5),(/0d0,m_W,0d0/),MomDK(1:4,1:3),PSWgt2)
      WMom(1:4) = MomDK(1:4,2)
      MomDK(1:4,4) = MomDK(1:4,3)
      PSWgt = PSWgt2

! flip=.not.flip
! if( flip ) then!  every second call the singular event is generated
!         Pcol1= 3 -1
!         Pcol2= 4 -1
!         SingDepth = 1d-5
!         Steps = 10
!         call gensing(3,m_Top,(/0d0,0d0,m_W/),MomDK(1:4,1:3),Pcol1,Pcol2,SingDepth,Steps); print *, "gensing activated"
!         PSWgt2=1d0
!         WMom(1:4) = MomDK(1:4,3)
!         MomDK(1:4,4) = MomDK(1:4,2)
! endif

      soft = dabs(MomDK(1,4))/m_Top
      coll = dabs(MomDK(1:4,1).dot.MomDK(1:4,4))/m_Top**2
      if( soft.lt.1d-6 .or. coll.lt.1d-10 ) PSWgt = 0d0

      call genps(2,m_W,xRndPS(6:7),(/0d0,0d0/),MomDK(1:4,2:3),PSWgt3)! W decay
!     boost leptons to the W frame:
      call boost(MomDK(1:4,2),WMom(1:4),m_W)
      call boost(MomDK(1:4,3),WMom(1:4),m_W)
!     boost all guys to the top frame:
      call boost(MomDK(1:4,1),TopMom(1:4),m_Top)
      call boost(MomDK(1:4,2),TopMom(1:4),m_Top)
      call boost(MomDK(1:4,3),TopMom(1:4),m_Top)
      call boost(MomDK(1:4,4),TopMom(1:4),m_Top)
      PSWgt = PSWgt * PSWgt3*PiWgt2
    endif


RETURN
END SUBROUTINE





SUBROUTINE EvalPhasespace_TopDecay2(TopMom,xRndPS,GluonRad,MomDK,PSWgt)!  top quark decay phase space with additional massless particle from W decay
use ModProcess
use ModMisc
use ModParameters
implicit none
logical GluonRad
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: TopMom(1:4),WMom(1:4),MomChk(1:4,1:3)
real(8) :: MomDK(:,:)
real(8) :: xRndPS(:)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)
real(8),parameter :: N3=3, PiWgt3 = (2d0*Pi)**(4-N3*3) * (4d0*Pi)**(N3-1)
real(8),parameter :: N4=4, PiWgt4 = (2d0*Pi)**(4-N4*3) * (4d0*Pi)**(N4-1)
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,s_glu_qb,s_glu_q,E_glu
logical,save :: flip=.true.

!     MomDK(1:4,i): i=1:bottom, 2:lepton, 3:neutrino, 4: gluon
      call genps(2,m_Top,xRndPS(1:2),(/0d0,m_W/),MomDK(1:4,1:2),PSWgt2)! top decay
      WMom(1:4) = MomDK(1:4,2)
      PSWgt = PSWgt2*PiWgt2

!       call genps(3,m_W,xRndPS(3:7),(/0d0,0d0,0d0/),MomDK(1:4,2:4),PSWgt3)! W decay
!       PSWgt = PSWgt * PSWgt3*PiWgt3
      call yeti3(m_W,xRndPS(3:7),(/0d0,0d0,0d0/),MomDK(1:4,2:4),PSWgt3)! W decay
      PSWgt = PSWgt * PSWgt3


! flip=.not.flip
! if( flip ) then!  every second call the singular event is generated
!         Pcol1= 5 -1
!         Pcol2= 5 -1
!         SingDepth = 1e-5
!         Steps = 5
!         call gensing(3,m_W,(/0d0,0d0,0d0/),MomDK(1:4,2:4),Pcol1,Pcol2,SingDepth,Steps); print *, "gensing activated"
!         PSWgt3=1d0
! endif

!     boost leptons to the W frame:
      call boost(MomDK(1:4,2),WMom(1:4),m_W)
      call boost(MomDK(1:4,3),WMom(1:4),m_W)
      call boost(MomDK(1:4,4),WMom(1:4),m_W)
!     boost all guys to the top frame:
      call boost(MomDK(1:4,1),TopMom(1:4),m_Top)
      call boost(MomDK(1:4,2),TopMom(1:4),m_Top)
      call boost(MomDK(1:4,3),TopMom(1:4),m_Top)
      call boost(MomDK(1:4,4),TopMom(1:4),m_Top)

!     cut out very collinear/soft configurations
      s_glu_qb = dabs(MomDK(1:4,4).dot.MomDK(1:4,2))/m_W**2
      s_glu_q  = dabs(MomDK(1:4,4).dot.MomDK(1:4,3))/m_W**2
      E_glu    = dabs(MomDK(1,4))/m_W
      if( s_glu_qb.lt.1d-10 .or. s_glu_q.lt.1d-10 .or. E_glu.lt.1d-5) then
          PSWgt = 0d0
          return
      endif

RETURN
END SUBROUTINE




SUBROUTINE EvalPhasespace_TopDecay3(TopMom,xRndPS,MomDK,PSWgt)!  top quark decay phase space with two additional massless particles
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: TopMom(1:4),WMom(1:4)
real(8) :: MomDK(:,:)
real(8) :: xRndPS(:)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)
real(8),parameter :: N3=3, PiWgt3 = (2d0*Pi)**(4-N3*3) * (4d0*Pi)**(N3-1)
real(8),parameter :: N4=4, PiWgt4 = (2d0*Pi)**(4-N4*3) * (4d0*Pi)**(N4-1)
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth, soft, coll

!     MomDK(1:4,i): i=1:bottom, 2:lepton, 3:neutrino, 4: gluon, 5: gluon/photon
    call genps(4,m_Top,xRndPS(1:8),(/m_W,0d0,0d0,0d0/),MomDK(1:4,1:4),PSWgt2)

!        Pcol1= 4 -1
!        Pcol2= 5 -1
!        SingDepth = 1e-10
!        Steps = 20
!        call gensing(4,m_Top,(/m_W,0d0,0d0,0d0/),MomDK(1:4,1:4),Pcol1,Pcol2,SingDepth,Steps)
!        PSWgt2=1d0

       WMom(1:4) = MomDK(1:4,1)
       MomDK(1:4,1) = MomDK(1:4,2)
       MomDK(1:4,5) = MomDK(1:4,4)
       MomDK(1:4,4) = MomDK(1:4,3)
       call genps(2,m_W,xRndPS(9:10),(/0d0,0d0/),MomDK(1:4,2:3),PSWgt3)! W decay
       !     boost leptons to the W frame:
       call boost(MomDK(1:4,2),WMom(1:4),m_W)
       call boost(MomDK(1:4,3),WMom(1:4),m_W)
       !     boost all guys to the top frame:
       call boost(MomDK(1:4,1),TopMom(1:4),m_Top)
       call boost(MomDK(1:4,2),TopMom(1:4),m_Top)
       call boost(MomDK(1:4,3),TopMom(1:4),m_Top)
       call boost(MomDK(1:4,4),TopMom(1:4),m_Top)
       call boost(MomDK(1:4,5),TopMom(1:4),m_Top)

! cut to singular events
   soft = dmin1(dabs(MomDK(1,4)),dabs(MomDK(1,5)))/TopMom(1)
   coll = dmin1(dabs(MomDk(1:4,1).dot.MomDk(1:4,4)),dabs(MomDk(1:4,1).dot.MomDk(1:4,5)),dabs(MomDk(1:4,4).dot.MomDk(1:4,5)))/TopMom(1)**2
! print *, "soft/coll",soft,coll
   if(soft.gt.1.d-6 .and. coll.gt.1.d-10 ) then
       PSWgt = PSWgt2*PiWgt4 * PSWgt3*PiWgt2
   else
!        MomDk(1:4,1:5) = 0.d0!   why is that necessary??  it leads to NaN in sq.mat.elements
       PSWgt = 0.d0
   endif


! if( MomDK(1,4)/M_top.lt.1d-3 ) PSWgt = 0d0
! if( MomDK(1,5)/M_top.lt.1d-3 ) PSWgt = 0d0
! if( dsqrt((MomDk(1:4,1).dot.MomDk(1:4,4))/m_top**2).lt.1d-3 ) PSWgt = 0d0!  has zero effect!!
! if( dsqrt((MomDk(1:4,1).dot.MomDk(1:4,5))/m_top**2).lt.1d-3 ) PSWgt = 0d0!  has zero effect!!
! if( dsqrt((MomDk(1:4,4).dot.MomDk(1:4,5))/m_top**2).lt.1d-3 ) PSWgt = 0d0


RETURN
END SUBROUTINE





SUBROUTINE EvalPhasespace_TopDecay4(TopMom,xRndPS,MomDK,PSWgt)!  top quark decay phase space with two additional massless particles from top and W line, resp.
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: TopMom(1:4),WMom(1:4),MomChk(1:4,1:4)
real(8) :: MomDK(1:4,1:5)
real(8) :: xRndPS(:)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)
real(8),parameter :: N3=3, PiWgt3 = (2d0*Pi)**(4-N3*3) * (4d0*Pi)**(N3-1)
real(8),parameter :: N4=4, PiWgt4 = (2d0*Pi)**(4-N4*3) * (4d0*Pi)**(N4-1)
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth, soft, coll1, coll2


!     MomDK(1:4,i): i=1:bottom, 2:lepton, 3:neutrino, 4: gluon/photon 5: gluon
    call genps(3,m_Top,xRndPS(1:5),(/0d0,0d0,m_W/),MomDK(1:4,1:3),PSWgt2)! top decay with additional gluon
!        Pcol1= 4 -1
!        Pcol2= 4 -1
!        SingDepth = 1e-15
!        Steps = 20
!        call gensing(3,m_Top,(/0d0,0d0,m_W/),MomDK(1:4,1:3),Pcol1,Pcol2,SingDepth,Steps)
!        PSWgt2=1d0

    soft = dabs(MomDK(1,2))/M_top
    coll1= dabs(MomDk(1:4,2).dot.MomDk(1:4,1))/m_top**2
    if(soft.gt.1d-6 .and. coll1.gt.1d-10) then!  reject too soft/collinear gluon configurations
        WMom(1:4) = MomDK(1:4,3)
        MomDK(1:4,5) = MomDK(1:4,2)
!         call genps(3,m_W,xRndPS(6:10),(/0d0,0d0,0d0/),MomDK(1:4,2:4),PSWgt3)
!         PSWgt = PSWgt2*PiWgt3 * PSWgt3*PiWgt3
        call yeti3(m_W,xRndPS(6:10),(/0d0,0d0,0d0/),MomDK(1:4,2:4),PSWgt3)! highest sensitivity: 3-1-2
        PSWgt = PSWgt2*PiWgt3 * PSWgt3
!               Pcol1= 5 -1
!               Pcol2= 5 -1
!               SingDepth = 1e-15
!               Steps = 20
!               call gensing(3,m_W,(/0d0,0d0,0d0/),MomDK(1:4,2:4),Pcol1,Pcol2,SingDepth,Steps)
!               PSWgt3=1d0

        soft = dabs(MomDK(1,4))/m_W
        coll1= dabs(MomDk(1:4,2).dot.MomDk(1:4,4))/m_W**2
        coll2= dabs(MomDk(1:4,3).dot.MomDk(1:4,4))/m_W**2
        if(soft.gt.1d-5 .and. coll1.gt.1d-10 .and. coll2.gt.1d-10) then!  reject too soft/collinear gluon configurations
!            boost leptons to the W frame:
             call boost(MomDK(1:4,2),WMom(1:4),m_W)
             call boost(MomDK(1:4,3),WMom(1:4),m_W)
             call boost(MomDK(1:4,4),WMom(1:4),m_W)
!            boost all guys to the top frame:
             call boost(MomDK(1:4,1),TopMom(1:4),m_Top)
             call boost(MomDK(1:4,2),TopMom(1:4),m_Top)
             call boost(MomDK(1:4,3),TopMom(1:4),m_Top)
             call boost(MomDK(1:4,4),TopMom(1:4),m_Top)
             call boost(MomDK(1:4,5),TopMom(1:4),m_Top)
        else
             PSWgt = 0d0
        endif
   else
         PSWgt = 0d0
   endif


RETURN
END SUBROUTINE





SUBROUTINE EvalPhasespace_TopDecay5(TopMom,xRndPS,MomDK,PSWgt)!  top quark decay phase space with two additional massless particles from W decay
use ModProcess
use ModMisc
use ModParameters
implicit none
logical GluonRad
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: TopMom(1:4),WMom(1:4)
real(8) :: MomDK(:,:)
real(8) :: xRndPS(:)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)
real(8),parameter :: N3=3, PiWgt3 = (2d0*Pi)**(4-N3*3) * (4d0*Pi)**(N3-1)
real(8),parameter :: N4=4, PiWgt4 = (2d0*Pi)**(4-N4*3) * (4d0*Pi)**(N4-1)
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,soft,coll

!     MomDK(1:4,i): i=1:bottom, 2:lepton, 3:neutrino, 4: gluon, 5: gluon
      call genps(2,m_Top,xRndPS(1:2),(/0d0,m_W/),MomDK(1:4,1:2),PSWgt2)! top decay
      WMom(1:4) = MomDK(1:4,2)
      call genps(4,m_W,xRndPS(3:10),(/0d0,0d0,0d0,0d0/),MomDK(1:4,2:5),PSWgt3)! W decay

      soft = dmin1(dabs(MomDK(1,4)),dabs(MomDK(1,5)))/M_W
      coll = dmin1(dabs(MomDk(1:4,4).dot.MomDk(1:4,5)),dabs(MomDk(1:4,2).dot.MomDk(1:4,4)),dabs(MomDk(1:4,2).dot.MomDk(1:4,5)), &
                   dabs(MomDk(1:4,3).dot.MomDk(1:4,4)),dabs(MomDk(1:4,3).dot.MomDk(1:4,5)))/m_W**2
      if(soft.gt.1.d-6 .and. coll.gt.1.d-11) then!  reject too soft/collinear gluon configurations
!               Pcol1= 4 -1
!               Pcol2= 6 -1
!               SingDepth = 1e-15
!               Steps = 15
!               call gensing(4,m_W,(/0d0,0d0,0d0,0d0/),MomDK(1:4,2:5),Pcol1,Pcol2,SingDepth,Steps)
!               PSWgt3=1d0
!         boost leptons to the W frame:
          call boost(MomDK(1:4,2),WMom(1:4),m_W)
          call boost(MomDK(1:4,3),WMom(1:4),m_W)
          call boost(MomDK(1:4,4),WMom(1:4),m_W)
          call boost(MomDK(1:4,5),WMom(1:4),m_W)
!         boost all guys to the top frame:
          call boost(MomDK(1:4,1),TopMom(1:4),m_Top)
          call boost(MomDK(1:4,2),TopMom(1:4),m_Top)
          call boost(MomDK(1:4,3),TopMom(1:4),m_Top)
          call boost(MomDK(1:4,4),TopMom(1:4),m_Top)
          call boost(MomDK(1:4,5),TopMom(1:4),m_Top)
          PSWgt = PSWgt2*PiWgt2 * PSWgt3*PiWgt4
       else
          PSWgt = 0d0
       endif



RETURN
END SUBROUTINE







SUBROUTINE EvalPhasespace_ZDecay(mZ_inv,ZMom,xRndPS,MomDK,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: PSWgt
real(8) :: ZMom(1:4)
real(8) :: MomDK(1:4,1:2)
real(8) :: xRndPS(1:2),mZ_inv
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)


        call genps(2,mZ_inv,xRndPS(1:2),(/0d0,0d0/),MomDK(1:4,1:2),PSWgt)! top decay

!       boost leptons to the Z frame:
        call boost(MomDK(1:4,1),ZMom(1:4),mZ_inv)
        call boost(MomDK(1:4,2),ZMom(1:4),mZ_inv)
        PSWgt = PSWgt*PiWgt2


RETURN
END SUBROUTINE










SUBROUTINE EvalPhaseSpace_2tobbbWW(EHat,xRndPS,Mom,PSWgt)
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: xRndPS(1:12)
real(8) :: Mom(1:4,1:10),MomAux(1:4,3:8)
real(8),parameter :: PiWgtPr4 = (2d0*Pi)**(4-4*3) * (4d0*Pi)**(4-1)
real(8),parameter :: PiWgtPr2 = (2d0*Pi)**(4-2*3) * (4d0*Pi)**(2-1)


    call genps(4,Ehat,xRndPS(1:8),(/0d0,m_W,0d0,m_W/),MomAux(1:4,3:6),PSWgt)
    PSWgt = PSWgt*PiWgtPr4
    Mom(1:4,5)=MomAux(1:4,3)
    Mom(1:4,8)=MomAux(1:4,5)


    call genps(2,m_W,xRndPS(9:10),(/0d0,0d0/),Mom(1:4,6:7),PSWgt2)
    PSWgt = PSWgt*PSWgt2*PiWgtPr2
    call boost(Mom(1:4,6),MomAux(1:4,4),m_W)
    call boost(Mom(1:4,7),MomAux(1:4,4),m_W)

    call genps(2,m_W,xRndPS(11:12),(/0d0,0d0/),Mom(1:4,9:10),PSWgt3)
    PSWgt = PSWgt*PSWgt3*PiWgtPr2
    call boost(Mom(1:4,9),MomAux(1:4,6),m_W)
    call boost(Mom(1:4,10),MomAux(1:4,6),m_W)

    Mom(1:4,3)=Mom(1:4,5)+Mom(1:4,6)+Mom(1:4,7)
    Mom(1:4,4)=Mom(1:4,8)+Mom(1:4,9)+Mom(1:4,10)


!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0





RETURN
END SUBROUTINE












SUBROUTINE EvalPhasespace_2to2(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5
real(8) :: Mom(1:4,1:4),MomW(1:4),xRndPS(1:2)
! integer :: NPart,i
! real(8) :: vel,theta
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)
real(8),parameter :: N3=3, PiWgt3 = (2d0*Pi)**(4-N3*3) * (4d0*Pi)**(N3-1)


!  generate PS: massless + massless --> massive(anti-top) + massive(top)
   call genps(2,Ehat,xRndPS(1:2),(/m_Top,m_Top/),Mom(1:4,3:4),PSWgt)
   PSWgt = PSWgt*PiWgt2

!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0

!   call SetKinVars(4,Mom(1:4,1:4),(/0d0,0d0,m_Top,m_Top/))

!     include "kinpointDK.MCFM.f90"

return
END SUBROUTINE







SUBROUTINE EvalPhasespace_2to2HT(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5
real(8) :: Mom(1:4,1:4),MomW(1:4),xRndPS(1:2)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)


!  generate PS: massless + massless --> massive(anti-top) + massive(top)
   call genps(2,Ehat,xRndPS(1:2),(/m_HTop,m_HTop/),Mom(1:4,3:4),PSWgt)
   PSWgt = PSWgt*PiWgt2

!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0


return
END SUBROUTINE




SUBROUTINE EvalPhasespace_2to3HT(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5
real(8) :: Mom(1:4,1:5),MomW(1:4),xRndPS(1:5)
real(8),parameter :: N3=3, PiWgt3 = (2d0*Pi)**(4-N3*3) * (4d0*Pi)**(N3-1)
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,velo,parx
real(8) :: s13,s23


   call genps(3,Ehat,xRndPS(1:5),(/0d0,m_Htop,m_Htop/),Mom(1:4,3:5),PSWgt)
   PSWgt = PSWgt*PiWgt3


!     Pcol1= 2 -1
!     Pcol2= 3 -1
!     SingDepth = 1e-10
!     Steps = 15
!     PSWgt = 1d0
!     call gensing(3,EHat,(/0d0,m_HTop,m_HTop/),Mom(1:4,3:5),Pcol1,Pcol2,SingDepth,Steps); print *, "generating singular point"


!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0

   s13 = Mom(1:4,1).dot.Mom(1:4,3)
   s23 = Mom(1:4,2).dot.Mom(1:4,3)
   if( abs(s13)/EHat**2.lt.1d-9 .or. abs(s23)/EHat**2.lt.1d-9 ) PSWgt=0d0
   if( abs(Mom(1,3)/EHat).lt.1d-5  ) PSWgt=0d0


return
END SUBROUTINE








SUBROUTINE EvalPhasespace_2to2Stops(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5
real(8) :: Mom(1:4,1:4),MomW(1:4),xRndPS(1:2)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)
integer,save :: it=1
real(8) :: beta,t,u,cos13,phi


   call genps(2,Ehat,xRndPS(1:2),(/m_Stop,m_Stop/),Mom(1:4,3:4),PSWgt)
   PSWgt = PSWgt*PiWgt2



! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
! print *, "INPUT MOMENTA FOR COMPARISON WITH RADJA"
! phi = 0.231231d0
! if(it.eq.1) then
!  EHat=dsqrt(1000000.0000000d0) *GeV
!  t = -836410.161513775d0 *GeV**2
!  u = -143589.838486225d0*GeV**2
! elseif(it.eq.2) then
!  EHat=dsqrt(1440000.00000000d0) *GeV
!  t = -1324817.04595758d0*GeV**2
!  u = -95182.9540424241d0*GeV**2
! elseif(it.eq.3) then
!  EHat=dsqrt(1960000.00000000d0) *GeV
!  t = -1454974.22611929d0*GeV**2
!  u = -485025.773880714d0*GeV**2
! endif
! 
!  PSWgt=1d0
!  m_Stop = dsqrt(0.5d0*(EHat**2+t+u))
! 
!  beta=dsqrt(1d0-m_stop**2/(EHat*0.5d0)**2)
!  cos13=1d0/beta*(1d0+(t-m_stop**2)/(EHat**2*0.5d0))
! 
! 
!  Mom(1,4) = EHat*0.5d0
!  Mom(2,4) = EHat*0.5d0*beta*dsqrt(1d0-cos13**2)*dsin(phi)
!  Mom(3,4) = EHat*0.5d0*beta*dsqrt(1d0-cos13**2)*dcos(phi)
!  Mom(4,4) = EHat*0.5d0*beta*cos13
! 
!  Mom(1,3) = EHat*0.5d0
!  Mom(2,3) =-EHat*beta*0.5d0*dsqrt(1d0-cos13**2)*dsin(phi)
!  Mom(3,3) =-EHat*0.5d0*beta*dsqrt(1d0-cos13**2)*dcos(phi)
!  Mom(4,3) =-EHat*beta*0.5d0*cos13
!  it=it+1
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 




!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0


! print *, "check",((Mom(1:4,3)).dot.(Mom(1:4,3)))-m_stop**2
! print *, "check",((Mom(1:4,4)).dot.(Mom(1:4,4)))-m_stop**2
! print *, "check",Mom(1,1)+Mom(1,2)-Mom(1,3)-Mom(1,4)
! print *, "check",Mom(2,1)+Mom(2,2)-Mom(2,3)-Mom(2,4)
! print *, "check",Mom(3,1)+Mom(3,2)-Mom(3,3)-Mom(3,4)
! print *, "check",Mom(4,1)+Mom(4,2)-Mom(4,3)-Mom(4,4)
! print *, "check",ehat**2+t+u-2d0*m_stop**2
! print *, "t",(Mom(1:4,1)-Mom(1:4,3)).dot.(Mom(1:4,1)-Mom(1:4,3))
! print *, "u",(Mom(1:4,1)-Mom(1:4,4)).dot.(Mom(1:4,1)-Mom(1:4,4))
! pause


return
END SUBROUTINE






SUBROUTINE EvalPhasespace_2to3Stops(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5
real(8) :: Mom(1:4,1:5),MomW(1:4),xRndPS(1:5)
real(8),parameter :: N3=3, PiWgt3 = (2d0*Pi)**(4-N3*3) * (4d0*Pi)**(N3-1)
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,velo,parx
real(8) :: s13,s23


   call genps(3,Ehat,xRndPS(1:5),(/0d0,m_Stop,m_Stop/),Mom(1:4,3:5),PSWgt)
   PSWgt = PSWgt*PiWgt3

!     Pcol1= 3 -1
!     Pcol2= 3 -1
!     SingDepth = 1e-10
!     Steps = 15
!     PSWgt = 1d0
!     call gensing(3,EHat,(/0d0,m_sTop,m_sTop/),Mom(1:4,3:5),Pcol1,Pcol2,SingDepth,Steps); print *, "generating singular point"

!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0

   s13 = Mom(1:4,1).dot.Mom(1:4,3)
   s23 = Mom(1:4,2).dot.Mom(1:4,3)
   if( abs(s13)/EHat**2.lt.1d-9 .or. abs(s23)/EHat**2.lt.1d-9 ) PSWgt=0d0
   if( abs(Mom(1,3)/EHat).lt.1d-5  ) PSWgt=0d0


return
END SUBROUTINE






SUBROUTINE EvalPhasespace_2to3(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: xRndPS(1:5)
real(8) :: Mom(1:4,1:5),TmpMom(1:4)
! real(8) :: MomDK(1:4,1:6)
! integer :: NPart,i
! real(8) :: vel,parx,theta ! for checks
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,velo,parx
real(8),parameter :: NPr=3, PiWgtPr = (2d0*Pi)**(4-NPr*3) * (4d0*Pi)**(NPr-1)


!  generate PS: massless + massless --> massless + massive(anti-top) + massive(top)
  call genps(3,Ehat,xRndPS(1:5),(/0d0,m_Top,m_Top/),Mom(1:4,3:5),PSWgt)
  PSWgt = PSWgt*PiWgtPr

!   call yeti3(Ehat,xRndPS(1:5),(/m_Top,m_Top,0d0/),Mom(1:4,3:5),PSWgt)
!   TmpMom(1:4) = Mom(1:4,3)
!   Mom(1:4,3)  = Mom(1:4,5)
!   Mom(1:4,5)  = TmpMom(1:4)

!      Pcol1= 3 -1
!      Pcol2= 3 -1
!      SingDepth = 1e-10
!      Steps = 10
!      PSWgt = 1d0
!      call gensing(3,EHat,(/0d0,m_Top,m_Top/),Mom(1:4,3:5),Pcol1,Pcol2,SingDepth,Steps); print *, "gensing activated"

!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0

!   call SetKinVars(5,Mom(1:4,1:5),(/0d0,0d0,0d0,m_Top,m_Top/))


!      print *, "using DUW kinpoint with mtop=174"
!      include './misc/DUWkinpoint'        ! Dittm.,Uwer,Weinzierl p1+p2 --> p3+p4+p5

!     print *, "using kinpoint3."
!     include 'misc/kinpoint3.'

!    include "kinpoint."        ! 0 --> p1+..+p4
!    Mom(1:4,1) = -Mom(1:4,1)     ! p1+p2 --> p3+p4
!    Mom(1:4,2) = -Mom(1:4,2)


return
END SUBROUTINE



SUBROUTINE EvalPhasespace_2to3M(EHat,Mass,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3,Mass
real(8) :: xRndPS(1:5)
real(8) :: Mom(1:4,1:5),TmpMom(1:4)
! real(8) :: MomDK(1:4,1:6)
! integer :: NPart,i
! real(8) :: vel,parx,theta ! for checks
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,velo,parx
real(8),parameter :: NPr=3, PiWgtPr = (2d0*Pi)**(4-NPr*3) * (4d0*Pi)**(NPr-1)


!  generate PS: massless + massless --> massive(M) + massive(anti-top) + massive(top)
  call genps(3,Ehat,xRndPS(1:5),(/Mass,m_Top,m_Top/),Mom(1:4,3:5),PSWgt)
  PSWgt = PSWgt*PiWgtPr

!   call yeti3(Ehat,xRndPS(1:5),(/m_Top,m_Top,Mass/),Mom(1:4,3:5),PSWgt)
!   TmpMom(1:4) = Mom(1:4,3)
!   Mom(1:4,3)  = Mom(1:4,5)
!   Mom(1:4,5)  = TmpMom(1:4)

!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0


return
END SUBROUTINE




SUBROUTINE EvalPhasespace_2to4M(EHat,Mass,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3,Mass
real(8) :: xRndPS(1:5)
real(8) :: Mom(1:4,1:5),TmpMom(1:4)
! real(8) :: MomDK(1:4,1:6)
! integer :: NPart,i
! real(8) :: vel,parx,theta ! for checks
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,velo,parx
real(8),parameter :: NPr=4, PiWgtPr = (2d0*Pi)**(4-NPr*3) * (4d0*Pi)**(NPr-1)


!  generate PS: massless + massless --> massless + massive(M) + massive(anti-top) + massive(top)
  call genps(4,Ehat,xRndPS(1:8),(/0d0,Mass,m_Top,m_Top/),Mom(1:4,3:6),PSWgt)
  PSWgt = PSWgt*PiWgtPr

!    Pcol1= 2 -1
!    Pcol2= 3 -1
!    SingDepth = 1d-16
!    Steps = 20
!    PSWgt = 1d0
!    call gensing(4,EHat,(/0d0,Mass,m_Top,m_Top/),Mom(1:4,3:6),Pcol1,Pcol2,SingDepth,Steps); print *, "gensing"


!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0


return
END SUBROUTINE




!!! Zprime section !!!

SUBROUTINE EvalPhasespaceBWMapp(EHat,Masses,xRndPS,Mom,PSWgt)                                                                                                
use ModParameters                                                                                                                                            
use ModMisc                                                                                                                                                  
implicit none                                                                                                                                                
real(8) :: EHat,Masses(1:3),xRndPS(1:3*3-4)                                                                                                                  
real(8) :: Mom(1:4,1:5),PSWgt                                                                                                                                
real(8) :: PiWgt2,SingDepth                                                                                                                                  
integer :: N,NVar,Pcol1,Pcol2,Steps                                                                                                                          
                                                                                                                                                             
   N=3                                                                                                                                                       
   NVar=3*N-4                                                                                                                                                
   PiWgt2 = (2d0*Pi)**(-NVar) * (4d0*Pi)**(N-1)                                                                                                              
   call genpszttbg(N,Ehat,xRndPS(1:NVar),M_Zpr,Ga_Zpr,m_Top,(/0d0,m_top,m_top/),Mom(1:4,3:5),PSWgt)                                                          
   PSWgt = PSWgt*PiWgt2                                                                                                                                      
                                                                                                                                                             
   Mom(1,1) =  EHat*0.5d0                                                                                                                                    
   Mom(2,1) =  0d0                                                                                                                                           
   Mom(3,1) =  0d0                                                                                                                                           
   Mom(4,1) = +EHat*0.5d0                                                                                                                                    
                                                                                                                                                             
   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0

   call swapmom(Mom(1:4,3),Mom(1:4,5))

   if( N.eq.3 ) then
      if( dmin1( (Mom(1:4,1).dot.Mom(1:4,5))/EHat**2,(Mom(1:4,2).dot.Mom(1:4,5))/EHat**2,(Mom(1,5)/EHat)**2 ).lt.1d-10 ) PSWgt = 0d0
   endif


return
END SUBROUTINE


!!! End Zprime section !!!



SUBROUTINE EvalPhasespace_2to4(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: xRndPS(1:8)
real(8) :: Mom(1:4,1:6)!,MomTmp(1:4,3:6)
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth
real(8),parameter :: NPr=4, PiWgtPr = (2d0*Pi)**(4-NPr*3) * (4d0*Pi)**(NPr-1)


!  generate PS: massless + massless --> massless + massless + massive(anti-top) + massive(top)
    call genps(4,Ehat,xRndPS(1:8),(/0d0,0d0,m_Top,m_Top/),Mom(1:4,3:6),PSWgt)
! call genps(4,Ehat,xRndPS(1:8),(/m_Top,m_Top,0d0,0d0/),MomTmp(1:4,3:6),PSWgt)
    PSWgt = PSWgt*PiWgtPr

! Mom(1:4,3)=MomTmp(1:4,5)
! Mom(1:4,4)=MomTmp(1:4,6)
! Mom(1:4,5)=MomTmp(1:4,3)
! Mom(1:4,6)=MomTmp(1:4,4)

!    Pcol1= 1 -1
!    Pcol2= 3 -1
!    SingDepth = 1d-10
!    Steps = 10
!    PSWgt = 1d0
!    call gensing(4,EHat,(/0d0,0d0,m_Top,m_Top/),Mom(1:4,3:6),Pcol1,Pcol2,SingDepth,Steps)

!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0

!     write(*,"(4F25.16)") Mom(1:4,1)  *1d0
!     write(*,"(4F25.16)") Mom(1:4,2)  *1d0
!     write(*,"(4F25.16)") Mom(1:4,3)  *1d0
!     write(*,"(4F25.16)") Mom(1:4,4)  *1d0
!     write(*,"(4F25.16)") Mom(1:4,5)  *1d0
!     write(*,"(4F25.16)") Mom(1:4,6)  *1d0
!
! !     print *, (Mom(1:4,5).dot.Mom(1:4,6))/(Mom(1:4,1).dot.Mom(1:4,2))
!     print *, (Mom(1:4,3).dot.Mom(1:4,4))/(Mom(1:4,1).dot.Mom(1:4,2))
!     pause
!


return
END SUBROUTINE



SUBROUTINE EvalPhasespace_2to5(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: xRndPS(1:11)
real(8) :: Mom(1:4,1:7)
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth
real(8),parameter :: NPr=4, PiWgtPr = (2d0*Pi)**(4-NPr*3) * (4d0*Pi)**(NPr-1)


!  generate PS: massless + massless --> massless + massless + massive(anti-top) + massive(top)
 call genps(5,Ehat,xRndPS(1:11),(/0d0,0d0,0d0,m_Top,m_Top/),Mom(1:4,3:7),PSWgt)
 PSWgt = PSWgt*PiWgtPr

!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0


return
END SUBROUTINE







SUBROUTINE EvalPhasespace_2to6(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
use ifport
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: xRndPS(1:20)
real(8) :: Mom(1:4,1:5),MomDK(1:4,1:6)
integer :: NPart,i
real(8) :: velo,parx,theta ! for checks
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth
real(8),parameter :: NPr=6, PiWgtPr = (2d0*Pi)**(4-NPr*3) * (4d0*Pi)**(NPr-1)


!  generate PS: massless + massless --> massless + massive(anti-top) + massive(top)
  call genps(6,Ehat,xRndPS(1:14),(/0d0,0d0,0d0,0d0,0d0,0d0/),Mom(1:4,3:8),PSWgt)
  PSWgt = PSWgt*PiWgtPr

!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0


return
END SUBROUTINE




SUBROUTINE EvalPhasespace_2toN(N,EHat,xRndPS,Mom,Mass,PSWgt)
use ModProcess
use ModMisc
use ModParameters
use ifport
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: xRndPS(1:(3*N+4))
real(8) :: Mom(1:4,1:N),Mass(1:N)
integer :: NPart,i,N
real(8) :: velo,parx,theta ! for checks
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth
real(8) :: PiWgtPr


   PiWgtPr= (2d0*Pi)**(4-N*3) * (4d0*Pi)**(N-1)
!  generate PS: massless + massless --> massless + massive(anti-top) + massive(top)
  call genps(N,Ehat,xRndPS(1:(3*N+4)),Mass,Mom(1:4,3:N+2),PSWgt)
  PSWgt = PSWgt*PiWgtPr

!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0


return
END SUBROUTINE







SUBROUTINE EvalPhasespace_2to4ML(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: xRndPS(1:8)
real(8) :: Mom(1:4,1:5),MomDK(1:4,1:6)
integer :: NPart,i

!  generate PS: massless + massless --> massless +  massless +  massless +  massless
   call genps(4,Ehat,xRndPS(1:8),(/0d0,0d0,0d0,0d0/),Mom(1:4,3:6),PSWgt)
!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0


return
END SUBROUTINE



SUBROUTINE EvalPhasespace_2to4ZDK(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: xRndPS(1:8)
real(8) :: Mom(1:4,1:6),MomDK(1:4,1:6)
integer :: NPart,i

!  generate PS: massless + massless --> massless +  massless +  massive(anti-top) +  massive(Top)
   call genps(4,Ehat,xRndPS(1:8),(/0d0,0d0,m_Top,m_Top/),Mom(1:4,3:6),PSWgt)
!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0


return
END SUBROUTINE EvalPhasespace_2to4ZDK











SUBROUTINE Kinematics_TTBARJET(NPlus1PS,MomExt,MomDK,applyPSCut,NBin)
use ModMisc
use ModParameters
implicit none
integer :: NumHadr,NPlus1PS
real(8) :: MomExt(1:4,1:5+NPlus1PS),MomDK(1:4,1:6),MomJet(1:4,1:8),zeros(1:10),MomJet_ordered(1:4,1:8)
real(8) :: MomHadr(1:4,1:8),MomLept(1:4,1:4)
real(8) :: MomBoost(1:4),MomW(1:4),MomTops(1:4,1:2)
logical :: applyPSCut
integer :: NBin(:),PartList(1:8),JetList(1:8),NJet,NObsJet,n,NObsJet_Tree,nWJets
real(8) :: pT_lepM,pT_lepP,pT_miss,pT_ATop,pT_Top,HT,m_lb,R_lb,m_bb,m_bj
real(8) :: eta_ATop,eta_Top,eta_lepM,eta_lepP,eta_miss,beta,costheta
real(8) :: pT_jet(1:8),eta_jet(1:8),eta_sepa,eta_Zeppi,s34,s35,s36,s45,s46,s56,mTopHadr,mTopLept
real(8) :: R_bb,MinvJets,MinvLept,phi_Lept,pT_lept,ET_lept,ET_miss,mT,pT_x,pT_y,MTW
real(8) :: MomTTbar(1:4),pT_ttbar,m_ttbar,y_top,y_Atop,y_ttbar,dy_tops,dphi_ttbar


!DEC$ IF(_CheckMomenta .EQ.1)
   zeros(:) = 0d0
   if( NPlus1PS.eq.0 ) then
        if(TopDecays.eq.0) then
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomExt(1:4,3) - MomExt(1:4,4) - MomExt(1:4,5)
        else
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomExt(1:4,3) - MomDK(1:4,1)-MomDK(1:4,2)-MomDK(1:4,3)-MomDK(1:4,4)-MomDK(1:4,5)-MomDK(1:4,6)
        endif
   else
        if(TopDecays.eq.0) then
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomExt(1:4,3) - MomExt(1:4,4) - MomExt(1:4,5) - MomExt(1:4,6)
        else
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomExt(1:4,4) - MomExt(1:4,3) - MomDK(1:4,1)-MomDK(1:4,2)-MomDK(1:4,3)-MomDK(1:4,4)-MomDK(1:4,5)-MomDK(1:4,6)
        endif
   endif
   if( any(abs(zeros(1:4)/MomExt(1,1)).gt.1d-6) ) then
      print *, "ERROR: energy-momentum violation in SUBROUTINE Kinematics_TTBARJET(",NPlus1PS,"): ",zeros(1:4)
      print *, "momenta dump:"
      print *, MomExt(1:4,1:5+NPlus1PS)
      print *, MomDK(1:4,1:6)
   endif
   if( NPlus1PS.eq.0 ) then
        zeros(1) = (MomExt(1:4,4).dot.MomExt(1:4,4)) - m_Top**2
        zeros(2) = (MomExt(1:4,5).dot.MomExt(1:4,5)) - m_Top**2
   else
        zeros(1) = (MomExt(1:4,5).dot.MomExt(1:4,5)) - m_Top**2
        zeros(2) = (MomExt(1:4,6).dot.MomExt(1:4,6)) - m_Top**2
        zeros(10)=  MomExt(1:4,4).dot.MomExt(1:4,4)
   endif
   zeros(3) =  MomExt(1:4,3).dot.MomExt(1:4,3)
   zeros(4) =  MomDK(1:4,1).dot.MomDK(1:4,1)
   zeros(5) =  MomDK(1:4,2).dot.MomDK(1:4,2)
   zeros(6) =  MomDK(1:4,3).dot.MomDK(1:4,3)
   zeros(7) =  MomDK(1:4,4).dot.MomDK(1:4,4)
   zeros(8) =  MomDK(1:4,5).dot.MomDK(1:4,5)
   zeros(9) =  MomDK(1:4,6).dot.MomDK(1:4,6)
   if( TopDecays.eq.0 .and. any(abs(zeros(1:10)/MomExt(1,1)**2).gt.1d-6) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TTBARJET(",NPlus1PS,"): ",zeros(1:10)
      print *, "momenta dump:"
      print *, MomExt(1:4,1:5+NPlus1PS)
      print *, MomDK(1:4,1:6)
   endif
   if( TopDecays.ne.0 .and. any(abs(zeros(1:10)/MomExt(1,1)**2).gt.1d-6) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TTBARJET(",NPlus1PS,"): ",zeros(1:10)
      print *, "momenta dump:"
      print *, MomExt(1:4,1:5+NPlus1PS)
      print *, MomDK(1:4,1:6)
   endif
!DEC$ ENDIF







! NumPart = particles in the final state
! required momentum order: MomExt(:,:): 1=In_left, 2=In_right, 3,4,...=light partons, N-1=ATop, N=Top
!                          MomDK(:,1:6) : 1=ABot, 2=lep-/q, 3=ANeu/qbar, 4=Bot, 5=lep+/qbar, 6=Neu/q

applyPSCut = .false.
NBin(1:NumHistograms) = 0

MomHadr(1:4,1:8) = 0d0
MomLept(1:4,1:4) = 0d0
PartList(1:8)=(/1,2,3,4,5,6,7,8/)


! separating momenta into hadron momenta and lepton momenta to which various procedures can be applied
!-------------------------------------------------------
if( TopDecays.eq.0 ) then  ! no top decays
   MomHadr(1:4,1) = MomExt(1:4,3)   ! q/qbar/glu
   if( NPlus1PS.eq.1 ) then
      MomHadr(1:4,2) = MomExt(1:4,4)   ! q/qbar/glu
      NumHadr = 2
      MomTops(1:4,1) = MomExt(1:4,5)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,6)   ! Top
   else
      NumHadr = 1
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
   endif
!-------------------------------------------------------
elseif( TopDecays.eq.1 ) then  ! full leptonic decay
  MomHadr(1:4,1) = MomDK(1:4,1)  ! ABot
  MomHadr(1:4,2) = MomDK(1:4,4)  ! Bot
  MomHadr(1:4,3) = MomExt(1:4,3) ! q/qbar/glu
  MomLept(1:4,1) = MomDK(1:4,2)  ! lep-
  MomLept(1:4,2) = MomDK(1:4,3)  ! ANeu
  MomLept(1:4,3) = MomDK(1:4,5)  ! lep+
  MomLept(1:4,4) = MomDK(1:4,6)  ! Neu

  if( NPlus1PS.eq.1 ) then
      MomHadr(1:4,4) = MomExt(1:4,4)   ! q/qbar/glu
      NumHadr = 4
      MomTops(1:4,1) = MomExt(1:4,5)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,6)   ! Top
  else
      NumHadr = 3
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.2 ) then  ! full hadronic decay
  MomHadr(1:4,1) = MomDK(1:4,1)  ! ABot
  MomHadr(1:4,2) = MomDK(1:4,4)  ! Bot
  MomHadr(1:4,3) = MomExt(1:4,3) ! q/qbar/glu
  MomHadr(1:4,4) = MomDK(1:4,2)  ! q
  MomHadr(1:4,5) = MomDK(1:4,3)  ! qbar
  MomHadr(1:4,6) = MomDK(1:4,5)  ! qbar
  MomHadr(1:4,7) = MomDK(1:4,6)  ! q
  if( NPlus1PS.eq.1 ) then
      MomHadr(1:4,8) = MomExt(1:4,4)   ! q/qbar/glu
      NumHadr = 8
      MomTops(1:4,1) = MomExt(1:4,5)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,6)   ! Top
  else
      NumHadr = 7
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.3 ) then  ! lept. Atop, hadr. top decay
  MomHadr(1:4,1) = MomDK(1:4,1)  ! ABot
  MomHadr(1:4,2) = MomDK(1:4,4)  ! Bot
  MomHadr(1:4,3) = MomExt(1:4,3) ! q/qbar/glu
  MomLept(1:4,1) = MomDK(1:4,2)  ! lep-
  MomLept(1:4,2) = MomDK(1:4,3)  ! ANeu
  MomHadr(1:4,4) = MomDK(1:4,5)  ! qbar
  MomHadr(1:4,5) = MomDK(1:4,6)  ! q
  if( NPlus1PS.eq.1 ) then
      MomHadr(1:4,6) = MomExt(1:4,4)   ! q/qbar/glu
      NumHadr = 6
      MomTops(1:4,1) = MomExt(1:4,5)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,6)   ! Top
  else
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.4 ) then  ! hadr. Atop, lept. top decay
  MomHadr(1:4,1) = MomDK(1:4,1)  ! ABot
  MomHadr(1:4,2) = MomDK(1:4,4)  ! Bot
  MomHadr(1:4,3) = MomExt(1:4,3) ! q/qbar/glu
  MomHadr(1:4,4) = MomDK(1:4,2)  ! q
  MomHadr(1:4,5) = MomDK(1:4,3)  ! qbar
  MomLept(1:4,3) = MomDK(1:4,5)  ! lep+
  MomLept(1:4,4) = MomDK(1:4,6)  ! Neu

  if( NPlus1PS.eq.1 ) then
      MomHadr(1:4,6) = MomExt(1:4,4)   ! q/qbar/glu
      NumHadr = 6
      MomTops(1:4,1) = MomExt(1:4,5)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,6)   ! Top
  else
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  endif

else
  call Error("this TopDecay is not implemented in Kinematics_TTBARJET")
endif


!---------------------- (anti) kT jet algorithm ---------------------------------
! important: b-quarks need to be at the first two positions of MomHadr(:,:)
! recombination of MomHadr(1:4,i) and MomHadr(1:4,j) results in MomJet(1:4,i) with i<j.
! but take care: some MomJet(1:4,i) can be zero

    NJet=0
    MomJet(1:4,1:8) = MomHadr(1:4,1:8)
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr)) ! hard protojets in beam pipe are counted as jets
    call pT_order(2,MomJet(1:4,1:2))! pT-order b-jets first
    call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))! pT-order non-b jets

!-------------------------------------------------------------------------
pT_jet(1:8)  = 0d0
eta_jet(1:8) = 0d0

if( ObsSet.eq.10 .or. ObsSet.eq.11 ) then! set of observables for ttbjet production without decays at Tevatron & LHC

    call pT_order(NJet,MomJet(1:4,1:NJet))! pT ordering of jet momenta

    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))

    eta_ATop = get_ETA(MomTops(1:4,1))
    eta_Top  = get_ETA(MomTops(1:4,2))

    pT_jet(1) = get_PT(MomJet(1:4,1))
    eta_jet(1)= get_ETA(MomJet(1:4,1))

    pT_jet(2) = get_PT(MomJet(1:4,2))
    eta_jet(2)= get_ETA(MomJet(1:4,2))

! check cuts
    NObsJet_Tree = 1
    NObsJet = 0
    do n=1,NJet
        if( pT_jet(n).gt.pT_jet_cut ) NObsJet = NObsJet +1
    enddo
    if( NObsJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif

!     if( pT_Top.lt.800d0*GeV ) then!   this is for the boosted observable
!         applyPSCut = .true.
!         RETURN
!     endif


! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,eta_ATop)
    NBin(4) = WhichBin(4,eta_Top)
    NBin(5) = WhichBin(5,eta_ATop)
    NBin(6) = WhichBin(6,eta_Top)
    NBin(7) = WhichBin(7,pT_jet(1))
    NBin(8) = WhichBin(8,eta_jet(1))


!-------------------------------------------------------
elseif( ObsSet.eq.12 ) then! set of observables for ttbjet production as signal process at Tevatron (semi-lept decay)

    NObsJet_Tree = 5

!   check that there are at least NObsJet_Tree resolved jets
    if( NJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif

!   check that there are two b jets    
    if( .not.(any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif

!   check that b jets pass pT and y cut
    if( get_pT(MomJet(1:4,1)).lt.pT_bjet_cut .or. abs(get_eta(MomJet(1:4,1))).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( get_pT(MomJet(1:4,2)).lt.pT_bjet_cut .or. abs(get_eta(MomJet(1:4,2))).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif

!   calculate number of observable jets
    NObsJet = 2 != two b-jets that already passed all criteria
    do n=3,NJet
        if( get_pT(MomJet(1:4,n)).gt.pT_jet_cut .and. abs(get_eta(MomJet(1:4,n))).lt.eta_jet_cut ) then
             NObsJet = NObsJet +1
            if( n.ne.NObsJet ) MomJet(1:4,NObsJet) = MomJet(1:4,n)
        endif
    enddo
    MomJet(1:4,NObsJet+1:NJet) = 0d0!  purging the remaining array

!   check that there are at least NObsJet_Tree observable jets
    if( NObsJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif





!   evaluate kinematic variables

    pT_lepP  = get_PT(MomLept(1:4,3))
    eta_lepP = get_eta(MomLept(1:4,3))

    pT_miss  = get_PT(MomLept(1:4,4))

    phi_Lept = dabs( Get_PHI(MomLept(1:4,3)) - Get_PHI(MomJet(1:4,1)) )
    if( phi_Lept.gt.Pi ) phi_Lept=2d0*Pi-phi_Lept

    m_lb = get_MInv(MomLept(1:4,3)+MomJet(1:4,1))
    R_lb = get_R(MomLept(1:4,3),MomJet(1:4,1))

    HT = pT_lepP + pT_miss
    do n=1,NObsJet
      pT_jet(n)  = get_pT( MomJet(1:4,n))
      eta_jet(n) = get_eta(MomJet(1:4,n))
      HT = HT + pT_jet(n)
    enddo




!    check cuts
     if( pT_lepP.lt.pT_lep_cut ) then
         applyPSCut = .true.
         RETURN
     endif

     if( abs(eta_lepP).gt.eta_lep_cut ) then
         applyPSCut = .true.
         RETURN
     endif

     if( pT_miss.lt.pT_miss_cut ) then
         applyPSCut = .true.
         RETURN
     endif

     if( HT.lt.HT_cut ) then
         applyPSCut = .true.
         RETURN
     endif


! binning

    NBin(:) = 0

    NBin(1) = WhichBin(1,pT_lepP)
    NBin(2) = WhichBin(2,eta_lepP)
    NBin(3) = WhichBin(3,pT_jet(3))
    NBin(4) = WhichBin(4,eta_jet(3))

    MomJet_ordered(1:4,1:NObsJet) = MomJet(1:4,1:NObsJet)
    call pT_order(NObsJet,MomJet_ordered(1:4,1:NObsJet))! pT ordering of jet momenta for b AND non-b jets
    NBin(5) = WhichBin(5,get_pT(MomJet_ordered(1:4,5)))
    NBin(6) = WhichBin(6,get_eta(MomJet_ordered(1:4,5)))

    NBin(7) = WhichBin(7,pT_miss)
    NBin(8) = WhichBin(8,HT)
    NBin(9) = WhichBin(9,m_lb)
    NBin(10)= WhichBin(10,phi_Lept)
    NBin(11)= WhichBin(11,R_lb)
    NBin(12)= WhichBin(12,eta_lepP)




! additional histograms for A_FB analysis
! ** ideally reconstructed tops
    MomTTbar(1:4) = MomTops(1:4,1)+MomTops(1:4,2)
    pT_ttbar = get_PT(MomTTbar(1:4))
    m_ttbar = get_MInv(MomTTbar(1:4))
    y_ttbar = get_ETA(MomTTbar(1:4))
    y_top = get_ETA(MomTops(1:4,2))
    y_ATop = get_ETA(MomTops(1:4,1))
    dy_tops  = y_top-y_ATop
    pT_top = get_PT(MomTops(1:4,2))
    dphi_ttbar = dabs( Get_PHI(MomTops(1:4,1)) - Get_PHI(MomTops(1:4,2)) )
    if( dphi_ttbar.gt.Pi ) dphi_ttbar=2d0*Pi-dphi_ttbar
    beta = dsqrt( abs(1d0-m_top**2/MomTops(1,2)**2) )
    costheta = Get_CosTheta(MomTops(1:4,2))


    NBin(13)= WhichBin(13,pT_ttbar)
    if(dy_tops.ge.0d0) NBin(14)= WhichBin(14,pT_ttbar)
    if(dy_tops.lt.0d0) NBin(15)= WhichBin(15,pT_ttbar)
    if(dy_tops.ge.0d0) NBin(16)= WhichBin(16,m_ttbar)
    if(dy_tops.lt.0d0) NBin(17)= WhichBin(17,m_ttbar)
    if(dy_tops.ge.0d0) NBin(18)= WhichBin(18,y_ttbar)
    if(dy_tops.lt.0d0) NBin(19)= WhichBin(19,y_ttbar)
    if(dy_tops.ge.0d0) NBin(20)= WhichBin(20,abs(dy_tops))
    if(dy_tops.lt.0d0) NBin(21)= WhichBin(21,abs(dy_tops))
    if(eta_lepP.ge.0d0) NBin(22)= WhichBin(22,pT_ttbar)
    if(eta_lepP.lt.0d0) NBin(23)= WhichBin(23,pT_ttbar)
    NBin(24)= WhichBin(24,y_top)
    NBin(25)= WhichBin(25,pT_top)
    if(dy_tops.ge.0d0) NBin(26)= WhichBin(26,dphi_ttbar)
    if(dy_tops.lt.0d0) NBin(27)= WhichBin(27,dphi_ttbar)

    NBin(43) = WhichBin(43,costheta)
    NBin(44) = WhichBin(44,beta*costheta)



! additional histograms for A_FB analysis
! ** realistically reconstructed tops


! reconstruct top (decays leptonically)
! reconstruct anti-top (decays hadronically)
!   find two non-b jets that are closest to MW mass
    if( NObsJet.eq.6 ) then
        s34= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,4))-M_W )
        s35= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,5))-M_W )
        s36= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,6))-M_W )
        s45= dabs( get_MInv(MomJet(1:4,4)+MomJet(1:4,5))-M_W )
        s46= dabs( get_MInv(MomJet(1:4,4)+MomJet(1:4,6))-M_W )
        s56= dabs( get_MInv(MomJet(1:4,5)+MomJet(1:4,6))-M_W )
    elseif( NObsJet.eq.5 ) then
        s34= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,4))-M_W )
        s35= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,5))-M_W )
        s45= dabs( get_MInv(MomJet(1:4,4)+MomJet(1:4,5))-M_W )
        s36=1d10; s46=1d10; s56=1d10
    else
        call Error("this should not happen")
    endif
    nWJets=minloc((/s34,s35,s45,s36,s46,s56/),1)

!   construct hadr. W momentum
    if(nWJets.eq.1) then
        MomW(1:4) = MomJet(1:4,3)+MomJet(1:4,4)
    elseif(nWJets.eq.2) then
        MomW(1:4) = MomJet(1:4,3)+MomJet(1:4,5)
    elseif(nWJets.eq.3) then
        MomW(1:4) = MomJet(1:4,4)+MomJet(1:4,5)
    elseif(nWJets.eq.4) then
        MomW(1:4) = MomJet(1:4,3)+MomJet(1:4,6)
    elseif(nWJets.eq.5) then
        MomW(1:4) = MomJet(1:4,4)+MomJet(1:4,6)
    elseif(nWJets.eq.6) then
        MomW(1:4) = MomJet(1:4,5)+MomJet(1:4,6)
    else
        MomW(1:4) = 0d0
    endif

    if( get_R(MomJet(1:4,1),MomLept(1:4,3)) .lt. get_R(MomJet(1:4,2),MomLept(1:4,3))  ) then ! find smaller R-distance between lepton and bjet
        MomTops(1:4,2) = MomJet(1:4,1) + MomLept(1:4,3)+MomLept(1:4,4)
        MomTops(1:4,1) = MomJet(1:4,2) + MomW(1:4)
    else
        MomTops(1:4,2) = MomJet(1:4,2) + MomLept(1:4,3)+MomLept(1:4,4)
        MomTops(1:4,1) = MomJet(1:4,1) + MomW(1:4)
    endif

! print *, "check j",NObsJet
! print *, "check W",get_MInv(MomW(1:4))*100d0
! print *, "check T",get_MInv(MomTops(1:4,1))*100d0,get_MInv(MomTops(1:4,2))*100d0
! print *, "xxx",get_MInv(MomJet(1:4,2)+MomLept(1:4,3)+MomLept(1:4,4))*100d0
! pause

!   require a 30 GeV window around M_W
!   require a 50 GeV window around M_Top
    if( dabs(get_MInv(MomW(1:4))-M_W).lt.30d0*GeV .and. & 
        dabs(get_MInv(MomTops(1:4,1))-M_Top).lt.50d0*GeV .and. dabs(get_MInv(MomTops(1:4,2))-M_Top).lt.50d0*GeV) then

        MomTTbar(1:4) = MomTops(1:4,1)+MomTops(1:4,2)
        pT_ttbar = get_PT(MomTTbar(1:4))
        m_ttbar = get_MInv(MomTTbar(1:4))
        y_ttbar = get_ETA(MomTTbar(1:4))
        y_top = get_ETA(MomTops(1:4,2))
        y_ATop = get_ETA(MomTops(1:4,1))
        dy_tops  = y_top-y_ATop
        pT_top = get_PT(MomTops(1:4,2))
        dphi_ttbar = dabs( Get_PHI(MomTops(1:4,1)) - Get_PHI(MomTops(1:4,2)) )
        if( dphi_ttbar.gt.Pi ) dphi_ttbar=2d0*Pi-dphi_ttbar

        NBin(28)= WhichBin(28,pT_ttbar)
        if(dy_tops.ge.0d0) NBin(29)= WhichBin(29,pT_ttbar)
        if(dy_tops.lt.0d0) NBin(30)= WhichBin(30,pT_ttbar)
        if(dy_tops.ge.0d0) NBin(31)= WhichBin(31,m_ttbar)
        if(dy_tops.lt.0d0) NBin(32)= WhichBin(32,m_ttbar)
        if(dy_tops.ge.0d0) NBin(33)= WhichBin(33,y_ttbar)
        if(dy_tops.lt.0d0) NBin(34)= WhichBin(34,y_ttbar)
        if(dy_tops.ge.0d0) NBin(35)= WhichBin(35,abs(dy_tops))
        if(dy_tops.lt.0d0) NBin(36)= WhichBin(36,abs(dy_tops))
        if(eta_lepP.ge.0d0) NBin(37)= WhichBin(37,pT_ttbar)
        if(eta_lepP.lt.0d0) NBin(38)= WhichBin(38,pT_ttbar)
        NBin(39)= WhichBin(39,y_top)
        NBin(40)= WhichBin(40,pT_top)
        if(dy_tops.ge.0d0) NBin(41)= WhichBin(41,dphi_ttbar)
        if(dy_tops.lt.0d0) NBin(42)= WhichBin(42,dphi_ttbar)


    endif






!-------------------------------------------------------
elseif( ObsSet.eq.13 ) then! set of observables for ttbjet production as signal process at LHC (di-lept decay)

    NObsJet_Tree = 3

!   check that there are at least NObsJet_Tree resolved jets
    if( NJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif

!   check that there are two b jets    
    if( .not.(any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif

!   check that b jets pass pT and y cut
    if( get_pT(MomJet(1:4,1)).lt.pT_bjet_cut .or. abs(get_eta(MomJet(1:4,1))).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( get_pT(MomJet(1:4,2)).lt.pT_bjet_cut .or. abs(get_eta(MomJet(1:4,2))).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif

!   calculate number of observable jets
    NObsJet = 2 != two b-jets that already passed all criteria
    do n=3,NJet
        if( get_pT(MomJet(1:4,n)).gt.pT_jet_cut .and. abs(get_eta(MomJet(1:4,n))).lt.eta_jet_cut ) then
             NObsJet = NObsJet +1
            if( n.ne.NObsJet ) MomJet(1:4,NObsJet) = MomJet(1:4,n)
        endif
    enddo
    MomJet(1:4,NObsJet+1:NJet) = 0d0!  purging the remaining array

!   check that there are at least NObsJet_Tree observable jets
    if( NObsJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif





!   evaluate kinematic variables

    pT_lepP  = get_PT(MomLept(1:4,3))
    eta_lepP = get_eta(MomLept(1:4,3))

    pT_lepM  = get_PT(MomLept(1:4,1))
    eta_lepM = get_eta(MomLept(1:4,1))

    pT_miss  = get_PT(MomLept(1:4,2)+MomLept(1:4,4))


    HT = pT_lepP + pT_lepM + pT_miss
    do n=1,NObsJet
      pT_jet(n)  = get_pT( MomJet(1:4,n))
      eta_jet(n) = get_eta(MomJet(1:4,n))
      HT = HT + pT_jet(n)
    enddo

    MinvLept = Get_MInv(MomLept(1:4,1)+MomLept(1:4,3))
    m_lb = get_MInv(MomLept(1:4,3)+MomJet(1:4,1))
    m_bb = get_MInv(MomJet(1:4,1)+MomJet(1:4,2))
    m_bj = get_MInv(MomJet(1:4,1)+MomJet(1:4,3))

    phi_Lept = dabs( Get_PHI(MomLept(1:4,1)) - Get_PHI(MomLept(1:4,3)) )
    if( phi_Lept.gt.Pi ) phi_Lept=2d0*Pi-phi_Lept

    pT_Top = get_PT( MomJet(1:4,1)+MomJet(1:4,2)+MomLept(1:4,1)+MomLept(1:4,2)+MomLept(1:4,3)+MomLept(1:4,4) )



! check cuts

    if( pT_lepP.lt.pT_lep_cut .or. pT_lepM.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lepP).gt.eta_lep_cut .or. abs(eta_lepM).gt.eta_lep_cut) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_miss.lt.pT_miss_cut ) then
        applyPSCut = .true.
        RETURN
    endif


! binning
    NBin(1) = WhichBin(1,pT_lepM)
    NBin(2) = WhichBin(2,eta_lepM)
    NBin(3) = WhichBin(3,pT_lepP)
    NBin(4) = WhichBin(4,eta_lepP)
    NBin(5) = WhichBin(5,HT)
    NBin(6) = WhichBin(6,MinvLept)
    NBin(7) = WhichBin(7,pT_jet(3))
    NBin(8) = WhichBin(8,eta_jet(3))
    NBin(9) = WhichBin(9,pT_miss)
    NBin(10) = WhichBin(10,pT_Top)
    NBin(11) = WhichBin(11,phi_Lept)
    NBin(12) = WhichBin(12,MInvLept)
    NBin(13) = WhichBin(13,pt_jet(1))
    NBin(14) = WhichBin(14,eta_jet(1))
    NBin(15) = WhichBin(15,m_bb)
    NBin(16) = WhichBin(16,m_bj)



!-------------------------------------------------------
elseif( ObsSet.eq.14 ) then! set of observables for ttbjet production as background process to VBF at LHC (di-lept decay)

print *, "STOP! proper jet selectin needs to be implemented first! see ObsSet=12,13"; stop

    NObsJet_Tree = 3
!   request at least two b-jets and one non-b-jet
    if( .not.(NJet.ge.NObsJet_Tree .and. any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif

    call pT_order(NJet,MomJet(1:4,1:NJet))! pT ordering of jet momenta for b- and non-b-jets, Note: b-jets are no longer in position 1 & 2

    pT_lepP  = get_PT(MomLept(1:4,3))
    eta_lepP = get_eta(MomLept(1:4,3))

    pT_lepM  = get_PT(MomLept(1:4,1))
    eta_lepM = get_eta(MomLept(1:4,1))

    do n=1,NJet
      pT_jet(n)  = get_pT( MomJet(1:4,n))
      eta_jet(n) = get_eta(MomJet(1:4,n))
    enddo

    MinvJets = Get_MInv(MomJet(1:4,1)+MomJet(1:4,2))
    MinvLept = Get_MInv(MomLept(1:4,1)+MomLept(1:4,3))

!     phi_Lept = dabs( Get_PHI(MomLept(1:4,1)) - Get_PHI(MomLept(1:4,3))  )
!     if( phi_Lept.gt.Pi ) phi_Lept=2d0*DblPi-phi_Lept
    phi_Lept = dacos( (MomLept(2,1)*MomLept(2,3)+MomLept(3,1)*MomLept(3,3))/pT_lepP/pT_lepM )  ! dacos el [0,Pi]
!PRINT *, "CHECK EQUIVALENCE!"
!STOP

    pT_lept = get_PT(MomLept(1:4,1)+MomLept(1:4,3))
    pT_miss = get_PT(MomLept(1:4,2)+MomLept(1:4,4))
    ET_lept = dsqrt(pT_lept**2 + MinvLept**2)
    ET_miss = dsqrt(pT_miss**2 + MinvLept**2)
    pT_x = MomLept(2,1)+MomLept(2,2)+MomLept(2,3)+MomLept(2,4)
    pT_y = MomLept(3,1)+MomLept(3,2)+MomLept(3,3)+MomLept(3,4)
    mT = dsqrt( (ET_lept+ET_miss)**2 - (pT_x**2+pT_y**2) )
    eta_sepa = abs(eta_jet(1)-eta_jet(2))


! check cuts
    if( pT_jet(1).lt. pT_hardestjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    NObsJet = 0
    do n=1,NJet
        if( pT_jet(n).gt.pT_jet_cut ) NObsJet = NObsJet +1
    enddo

    if( NObsJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif

    if( eta_jet(1)*eta_jet(2).gt.0d0 .or. eta_sepa.lt.eta_sepa_cut ) then
        applyPSCut = .true.
        RETURN
    endif


!   a "veto jet" must lie in between the two tagged jets
    if( abs(eta_jet(3)-eta_jet(1)).lt.eta_sepa .and. abs(eta_jet(3)-eta_jet(2)).lt.eta_sepa ) then  ! jet(3) is "veto jet"
        eta_Zeppi = eta_jet(3) - 0.5d0*(eta_jet(1)+eta_jet(2))
    else  ! there is no "veto jet"
        eta_Zeppi = -100d0
    endif

!     do n=1,NJet
!         if( pT_jet(n).gt.pT_jet_cut .and. abs(eta_jet(n)).lt.eta_jet_cut ) then
!             applyPSCut = .true.
!             RETURN
!         endif
!     enddo

    if( pT_lepP.lt.pT_lep_cut .or. pT_lepM.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lepP).gt.eta_lep_cut .or. abs(eta_lepM).gt.eta_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( MinvJets.lt.Minv_jets_cut ) then
        applyPSCut = .true.
        RETURN
    endif


! binning
    NBin(1) = WhichBin(1,pT_lepM)
    NBin(2) = WhichBin(2,eta_lepM)
    NBin(3) = WhichBin(3,pT_lepP)
    NBin(4) = WhichBin(4,eta_lepP)
    NBin(5) = WhichBin(5,MinvLept)
    NBin(6) = WhichBin(6,pT_jet(1))
    NBin(7) = WhichBin(7,phi_Lept)
    NBin(8) = WhichBin(8,mT)
    NBin(9) = WhichBin(9,eta_jet(1))
    NBin(10) = WhichBin(10,eta_jet(2))
    NBin(11) = WhichBin(11,eta_Zeppi)




elseif( ObsSet.eq.15 ) then! set of observables for ttbjet production as signal process at LHC (semi-lept decay)

    NObsJet_Tree = 5

!   check that there are at least NObsJet_Tree resolved jets
    if( NJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif

!   check that there are two b jets    
    if( .not.(any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif

!   check that b jets pass pT and y cut
    if( get_pT(MomJet(1:4,1)).lt.pT_bjet_cut .or. abs(get_eta(MomJet(1:4,1))).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( get_pT(MomJet(1:4,2)).lt.pT_bjet_cut .or. abs(get_eta(MomJet(1:4,2))).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif

!   calculate number of observable jets
    NObsJet = 2 != two b-jets that already passed all criteria
    do n=3,NJet
        if( get_pT(MomJet(1:4,n)).gt.pT_jet_cut .and. abs(get_eta(MomJet(1:4,n))).lt.eta_jet_cut ) then
             NObsJet = NObsJet +1
            if( n.ne.NObsJet ) MomJet(1:4,NObsJet) = MomJet(1:4,n)
        endif
    enddo
    MomJet(1:4,NObsJet+1:NJet) = 0d0!  purging the remaining array

!   check that there are at least NObsJet_Tree observable jets
    if( NObsJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif





!   evaluate kinematic variables

    pT_lepP  = get_PT(MomLept(1:4,3))
    eta_lepP = get_eta(MomLept(1:4,3))

    pT_miss  = get_PT(MomLept(1:4,4))

    MTW = dsqrt( (MomLept(2,3)+MomLept(2,4))**2 + (MomLept(3,3)+MomLept(3,4))**2 )

    phi_Lept = dabs( Get_PHI(MomLept(1:4,3)) - Get_PHI(MomJet(1:4,1)) )
    if( phi_Lept.gt.Pi ) phi_Lept=2d0*Pi-phi_Lept

    m_lb = get_MInv(MomLept(1:4,3)+MomJet(1:4,1))
    R_lb = get_R(MomLept(1:4,3),MomJet(1:4,1))

    HT = pT_lepP + pT_miss
    do n=1,NObsJet
      pT_jet(n)  = get_pT( MomJet(1:4,n))
      eta_jet(n) = get_eta(MomJet(1:4,n))
      HT = HT + pT_jet(n)
    enddo



!    check cuts
     if( pT_lepP.lt.pT_lep_cut ) then
         applyPSCut = .true.
         RETURN
     endif

     if( abs(eta_lepP).gt.eta_lep_cut ) then
         applyPSCut = .true.
         RETURN
     endif

     if( pT_miss.lt.pT_miss_cut ) then
         applyPSCut = .true.
         RETURN
     endif

!     if( pT_miss+MTW.lt.60d0*GeV ) then!   additional cut for ATLAS muon analysis
     if( MTW.lt.30d0*GeV ) then!   additional cut for ATLAS electron analysis
         applyPSCut = .true.
         RETURN
     endif


!   find two non-b jets that are closest to MW mass
    if( NObsJet.eq.6 ) then
        s34= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,4))-M_W )
        s35= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,5))-M_W )
        s36= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,6))-M_W )
        s45= dabs( get_MInv(MomJet(1:4,4)+MomJet(1:4,5))-M_W )
        s46= dabs( get_MInv(MomJet(1:4,4)+MomJet(1:4,6))-M_W )
        s56= dabs( get_MInv(MomJet(1:4,5)+MomJet(1:4,6))-M_W )
    elseif( NObsJet.eq.5 ) then
        s34= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,4))-M_W )
        s35= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,5))-M_W )
        s45= dabs( get_MInv(MomJet(1:4,4)+MomJet(1:4,5))-M_W )
        s36=1d10; s46=1d10; s56=1d10
    endif
    nWJets=minloc((/s34,s35,s45,s36,s46,s56/),1)

!   construct hadr. W momentum
    MomTops(1:4,1) = MomJet(1:4,1)+MomJet(1:4,2)+MomLept(1:4,3)+MomLept(1:4,4)
    if(nWJets.eq.1) then
        MomW(1:4) = MomJet(1:4,3)+MomJet(1:4,4)
    elseif(nWJets.eq.2) then
        MomW(1:4) = MomJet(1:4,3)+MomJet(1:4,5)
    elseif(nWJets.eq.3) then
        MomW(1:4) = MomJet(1:4,4)+MomJet(1:4,5)
    elseif(nWJets.eq.4) then
        MomW(1:4) = MomJet(1:4,3)+MomJet(1:4,6)
    elseif(nWJets.eq.5) then
        MomW(1:4) = MomJet(1:4,4)+MomJet(1:4,6)
    elseif(nWJets.eq.6) then
        MomW(1:4) = MomJet(1:4,5)+MomJet(1:4,6)
    else
        MomW(1:4) = 0d0
    endif
    MomTops(1:4,1) = MomTops(1:4,1) + MomW(1:4)!   construct the t+bar system
    if( dmin1(s34,s35,s45,s36,s46,s56).lt.20d0*GeV ) then!   require a 20GeV window around M_W
        pT_Top = get_pT(MomTops(1:4,1))
    else
        pT_Top   = -1d0
    endif



! binning
    NBin(1) = WhichBin(1,pT_lepP)
    NBin(2) = WhichBin(2,eta_lepP)
    NBin(3) = WhichBin(3,pT_jet(3))
    NBin(4) = WhichBin(4,eta_jet(3))

    call pT_order(NObsJet,MomJet(1:4,1:NObsJet))! pT ordering of jet momenta for b AND non-b jets
    NBin(5) = WhichBin(5,get_pT(MomJet(1:4,5)))
    NBin(6) = WhichBin(6,get_eta(MomJet(1:4,5)))

    NBin(7) = WhichBin(7,pT_miss)
    NBin(8) = WhichBin(8,HT)
    NBin(9) = WhichBin(9,m_lb)
    NBin(10)= WhichBin(10,phi_Lept)
    NBin(11)= WhichBin(11,R_lb)
    NBin(12)= WhichBin(12,eta_lepP)
    NBin(13)= WhichBin(13,pT_Top)




!-------------------------------------------------------
elseif( ObsSet.eq.19 ) then! for checks of ttbjet



!DEC$ IF(_FactCheck .EQ.1)


!    this is for "checkCD"
! if( NPlus1PS.eq.0) then
!     pT_jet(3)  = get_pT(MomExt(1:4,3))  ! leading non-b jet
!     if( pT_jet(3).lt.pT_jet_cut ) then
!         applyPSCut = .true.
!         RETURN
!     endif
! else
!     if( get_R(MomExt(1:4,3),MomExt(1:4,4)).lt.0.4d0 ) then
!         pT_jet(3)  = get_pT(MomExt(1:4,3)+MomExt(1:4,4))
!         if( pT_jet(3).lt.pT_jet_cut ) then
!             applyPSCut = .true.
!             RETURN
!         endif
!     else
!         pT_jet(3)  = get_pT(MomExt(1:4,3))
!         pT_jet(4)  = get_pT(MomExt(1:4,4))
!         if( pT_jet(3).lt.pT_jet_cut .and. pT_jet(4).lt.pT_jet_cut ) then
!             applyPSCut = .true.
!             RETURN
!         endif
!     endif
! endif

RETURN!     this is needed to avoid the cuts below
!DEC$ ENDIF






    NObsJet_Tree = 3

    if( TopDecays.eq.2 ) NObsJet_Tree = NObsJet_Tree + 4
    if( TopDecays.eq.3 .or. TopDecays.eq.4 ) NObsJet_Tree = NObsJet_Tree + 2
    if( NJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif
    pT_lepP  = get_PT(MomLept(1:4,3))
    eta_lepP = get_eta(MomLept(1:4,3))

!     pT_miss  = get_PT(MomLept(1:4,4))
!     phi_Lept = dabs( Get_PHI(MomLept(1:4,3)) - Get_PHI(MomJet(1:4,1)) )
!     if( phi_Lept.gt.Pi ) phi_Lept=2d0*Pi-phi_Lept
!     m_lb = get_MInv(MomLept(1:4,3)+MomJet(1:4,1))
!     R_lb = get_R(MomLept(1:4,3),MomJet(1:4,1))


    do n=1,NJet! first two jets are always b-jets
      pT_jet(n)  = get_pT( MomJet(1:4,n))
      eta_jet(n) = get_eta(MomJet(1:4,n))
    enddo

!     HT = pT_lepP + pT_miss
!     do n=1,NJet
!         HT = HT + pT_jet(n)
!     enddo


! check cuts
    NObsJet = 0
    do n=1,NJet
        if( pT_jet(n).gt.pT_jet_cut ) NObsJet = NObsJet +1
    enddo
    if( NObsJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif
    if( TOPDECAYS.eq.1 ) then
        pT_lepP  = get_PT(MomLept(1:4,3))
        eta_lepP = get_eta(MomLept(1:4,3))
        pT_lepM  = get_PT(MomLept(1:4,1))
        eta_lepM = get_eta(MomLept(1:4,1))
        pT_miss  = get_PT (MomLept(1:4,2)+MomLept(1:4,4))
        if( pT_miss.lt.0.2d0 ) then
            applyPSCut = .true.
            RETURN
        endif
        if( pT_lepP.lt.0.2d0 .or. pT_lepM.lt.0.2d0 ) then
            applyPSCut = .true.
            RETURN
        endif
        if( abs(eta_lepP).gt.2.1d0 .or. abs(eta_lepM).gt.2.1d0 ) then
            applyPSCut = .true.
            RETURN
        endif
    endif



    if( TOPDECAYS.eq.3 ) then
        pT_lepM  = get_PT(MomLept(1:4,1))
        eta_lepM = get_eta(MomLept(1:4,1))
        pT_miss  = get_PT (MomLept(1:4,2))
        if( pT_miss.lt.0.2d0 ) then
            applyPSCut = .true.
            RETURN
        endif
        if( pT_lepM.lt.0.2d0 ) then
            applyPSCut = .true.
            RETURN
        endif
        if( abs(eta_lepM).gt.2.1d0 ) then
            applyPSCut = .true.
            RETURN
        endif
    endif



    if( TOPDECAYS.eq.4 ) then
        pT_lepP  = get_PT(MomLept(1:4,3))
        eta_lepP = get_eta(MomLept(1:4,3))
        pT_miss  = get_PT (MomLept(1:4,4))
        if( pT_miss.lt.0.2d0 ) then
            applyPSCut = .true.
            RETURN
        endif
        if( pT_lepP.lt.0.2d0 ) then
            applyPSCut = .true.
            RETURN
        endif
        if( abs(eta_lepP).gt.2.1d0 ) then
            applyPSCut = .true.
            RETURN
        endif
    endif


! binning
    NBin(1) = WhichBin(1,pT_lepP)
    NBin(2) = WhichBin(2,eta_lepP)
    NBin(3) = WhichBin(3,get_pT(MomJet(1:4,1)))
    NBin(4) = WhichBin(4,get_eta(MomJet(1:4,1)))
    NBin(5) = WhichBin(5,pT_miss)
    NBin(6) = WhichBin(6,HT)
    NBin(7) = WhichBin(7,m_lb)
    NBin(8)= WhichBin(8,phi_Lept)
    NBin(9)= WhichBin(9,R_lb)
    NBin(10)= WhichBin(10,eta_lepP)





!-------------------------------------------------------
else
  print *, "ObsSet not implemented",ObsSet
  stop
endif


return
END SUBROUTINE







SUBROUTINE Kinematics_TTBARPHOTON(NPlus1PS,Mom,MomOrder,applyPSCut,NBin)
use ModMisc
use ModParameters
implicit none
integer :: NumHadr,NPlus1PS,MomOrder(1:12)
real(8) :: Mom(1:4,1:12),zeros(1:12)
real(8) :: MomJet(1:4,1:7),MomJet_CHECK(1:4,1:7)
real(8) :: MomHadr(1:4,0:8)
real(8) :: MomBoost(1:4),MomMiss(1:4),MomObs(1:4)
logical :: applyPSCut,isolated
integer :: NBin(:),PartList(1:7),JetList(1:7),NJet,NObsJet,k,NObsJet_Tree,NJet_CHECK
real(8) :: pT_lepM,pT_lepP,ET_miss,pT_ATop,pT_Top,HT,ET_bjet
real(8) :: eta_ATop,eta_Top,eta_lepM,eta_lepP,m_lb,m_jj,mTblP,m_jjb,mT_lp
real(8) :: pT_jet(1:7),eta_jet(1:7),eta_sepa,pt_Pho,eta_Pho,Rphobjet,mT_bln(1:2),mT_blnp(1:2)
real(8) :: R_Pj(1:5),R_lj(1:5),R_PlepP,R_PlepM,pT_lept,ET_lept,mT,MInvPb1jj,mTb2lP,MInvPb2jj,mTb1lP,Phi_LP,Phi_LL
integer :: tbar,t,pho,inLeft,inRight,realp,bbar,lepM,nubar,b,lepP,nu,qdn,qbup,qbdn,qup,L,N




! momentum ordering
  tbar    = MomOrder(1)
  t       = MomOrder(2)
  pho     = MomOrder(3)
  inLeft  = MomOrder(4)
  inRight = MomOrder(5)
  realp   = MomOrder(6)
  bbar    = MomOrder(7)
  lepM    = MomOrder(8)
  nubar   = MomOrder(9)
  b       = MomOrder(10)
  lepP    = MomOrder(11)
  nu      = MomOrder(12)

  qdn    = lepM
  qbup   = nubar
  qbdn   = lepP
  qup    = nu


!DEC$ IF(_CheckMomenta .EQ.1)
   if(TopDecays.eq.0) then
      zeros(1:4) = Mom(1:4,inLeft)+Mom(1:4,inRight) - Mom(1:4,tbar) - Mom(1:4,t) - Mom(1:4,pho)
   else
      zeros(1:4) = Mom(1:4,inLeft)+Mom(1:4,inRight) - Mom(1:4,pho) - Mom(1:4,bbar)-Mom(1:4,lepM)-Mom(1:4,nubar)-Mom(1:4,b)-Mom(1:4,lepP)-Mom(1:4,nu)
   endif
   if( NPlus1PS.eq.1 ) zeros(1:4) = zeros(1:4) - Mom(1:4,realp)
   if( any(abs(zeros(1:4)/Mom(1,inLeft)).gt.1d-8) ) then
      print *, "ERROR: energy-momentum violation in SUBROUTINE Kinematics_TTBARPHOTON(): ",zeros(1:4)
      print *, "momenta dump:"
      print *, Mom(1:4,1:12)
   endif
   zeros(1) = (Mom(1:4,tbar).dot.Mom(1:4,tbar)) - m_Top**2
   zeros(2) = (Mom(1:4,t).dot.Mom(1:4,t)) - m_Top**2
   zeros(3) =  Mom(1:4,pho).dot.Mom(1:4,pho)
   zeros(4) =  Mom(1:4,bbar).dot.Mom(1:4,bbar)
   zeros(5) =  Mom(1:4,lepM).dot.Mom(1:4,lepM)
   zeros(6) =  Mom(1:4,nubar).dot.Mom(1:4,nubar)
   zeros(7) =  Mom(1:4,b).dot.Mom(1:4,b)
   zeros(8) =  Mom(1:4,lepP).dot.Mom(1:4,lepP)
   zeros(9) =  Mom(1:4,nu).dot.Mom(1:4,nu)
   if( NPlus1PS.eq.1 ) zeros(10)=  Mom(1:4,realp).dot.Mom(1:4,realp)
   if( TopDecays.eq.0 .and. any(abs(zeros(1:10)/Mom(1,inLeft)**2).gt.1d-8) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TTBARPHOTON(): ",zeros(1:10)
      print *, Mom(1:4,1:2)
   endif
   if( TopDecays.ne.0 .and. any(abs(zeros(1:10)/Mom(1,inLeft)**2).gt.1d-8) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TTBARPHOTON(): ",zeros(1:10)
      print *, "momenta dump:"
      print *, Mom(1:4,1:12)
   endif
!DEC$ ENDIF


applyPSCut = .false.
NBin(1:NumHistograms) = 0


MomHadr(1:4,1:7) = 0d0
PartList(1:7)=(/1,2,3,4,5,6,7/)

MomMiss(1:4) = 0d0
MomObs(1:4)  = 0d0

! separating momenta into hadron momenta and lepton momenta to which various procedures can be applied
!-------------------------------------------------------

IF( TOPDECAYS.EQ.0 ) THEN  ! no top decays
    if(NPlus1PS.eq.0) then
        NumHadr = 0
    else
        NumHadr = 1
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
    endif

ELSEIF( TOPDECAYS.EQ.1 ) THEN  ! di-leptonic decay
    MomHadr(1:4,1) = Mom(1:4,bbar)
    MomHadr(1:4,2) = Mom(1:4,b)     ! Bot
    if(NPlus1PS.eq.0) then
        NumHadr = 2
    else
        NumHadr = 3
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
    endif

ELSEIF( TOPDECAYS.EQ.3 ) THEN  ! lept. Atop, hadr. top decay
   MomHadr(1:4,1) = Mom(1:4,bbar)
   MomHadr(1:4,2) = Mom(1:4,b)
   MomHadr(1:4,3) = Mom(1:4,qbdn)
   MomHadr(1:4,4) = Mom(1:4,qup)
   L = LepM
   N = nubar
   if(NPlus1PS.eq.0) then
        NumHadr = 4
   else
        NumHadr = 5
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)   ! q/qbar/glu
   endif

ELSEIF( TOPDECAYS.EQ.4 ) THEN  ! hadr. Atop, lept. top decay
   MomHadr(1:4,1) = Mom(1:4,bbar)
   MomHadr(1:4,2) = Mom(1:4,b)
   MomHadr(1:4,3) = Mom(1:4,qbup)
   MomHadr(1:4,4) = Mom(1:4,qdn)
   L = LepP
   N = nu
   if(NPlus1PS.eq.0) then
        NumHadr = 4
   else
        NumHadr = 5
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
   endif

ELSE
  call Error("this decay is not yet implemented")
ENDIF



!------------------ Frixione photon isolation ----------------------------

    isolated = FrixioneIsolated(Mom(1:4,pho),Rsep_Pj,NumHadr,MomHadr(1:4,1:NumHadr))
    if( .not. isolated ) then
        applyPSCut = .true.
        RETURN
    endif


!---------------------- kT jet algorithm ---------------------------------
! important: b-quarks need to be at the first two positions of MomHadr(:,:)
! recombination of MomHadr(1:4,i) and MomHadr(1:4,j) results in MomJet(1:4,i) with i<j.
! but take care: some MomJet(1:4,i) can be zero, this is why we apply pt_order

    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr)) ! hard protojets in beam pipe are counted as jets
    call pT_order(2,MomJet(1:4,1:2))
    call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))


! call SwitchEnergyComponent(MomHadr(1:4,1:NumHadr))
! call fastjetppgenkt(MomHadr(1:4,1:NumHadr),NumHadr,Rsep_jet,-1d0,MomJet(1:4,1:NumHadr),NJet)! 4th argument:  +1d0=kT,  -1d0=anti-kT,   0d0=CA
! call SwitchEnergyComponentBack(MomJet_CHECK(1:4,1:NumHadr))


!------------------------ cuts and binning --------------------------------
if( ObsSet.eq.20 .or. ObsSet.eq.21) then! ttb+photon production without top decays at Tevatron & LHC

    pT_ATop = get_PT(Mom(1:4,tbar))
    pT_Top  = get_PT(Mom(1:4,t))

    eta_ATop = get_ETA(Mom(1:4,tbar))
    eta_Top  = get_ETA(Mom(1:4,t))

    pT_Pho  = get_PT(Mom(1:4,pho))
    eta_Pho = get_ETA(Mom(1:4,pho))

!     if(NumHadr.eq.1) then
!         if( dabs(dacos((MomHadr(2,1)*Mom(2,pho)+MomHadr(3,1)*Mom(3,pho)+MomHadr(4,1)*Mom(4,pho))/MomHadr(1,1)/Mom(1,pho))).lt.3d0/360d0*2d0*DblPi ) then !   for Chinese check
!             applyPSCut = .true.
!             RETURN
!         endif
!         R_Pj(1)  = get_R(Mom(1:4,pho),MomHadr(1:4,1))
!         if( R_Pj(1).lt.Rsep_Pj ) then
!             applyPSCut = .true.
!             RETURN
!         endif
!     endif

! check cuts
    if(pt_Pho.lt.pT_pho_cut) then
        applyPSCut = .true.
        RETURN
    endif

! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,eta_ATop)
    NBin(4) = WhichBin(4,eta_Top)
    NBin(5) = WhichBin(5,pT_Pho)
    NBin(6) = WhichBin(6,eta_Pho)
    NBin(7) = WhichBin(7,eta_ATop)
    NBin(8) = WhichBin(8,eta_Top)

!-------------------------------------------------------
elseif( ObsSet.eq.22 .or. ObsSet.eq.23 ) then! set of observables for ttb+gamma production with di-lept. decays at the Tevatron & LHC

! request at least two b-jets
    NObsJet_Tree = 2
    if( .not.(NJet.ge.NObsJet_Tree .and. any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif


    pT_ATop = get_PT(Mom(1:4,tbar))
    pT_Top  = get_PT(Mom(1:4,t))

    eta_ATop = get_ETA(Mom(1:4,tbar))
    eta_Top  = get_ETA(Mom(1:4,t))

    pT_Pho  = get_PT(Mom(1:4,pho))
    eta_Pho = get_ETA(Mom(1:4,pho))

    R_PlepM = get_R(Mom(1:4,pho),Mom(1:4,lepM))
    R_PlepP = get_R(Mom(1:4,pho),Mom(1:4,lepP))

    pT_jet(1) = get_PT(MomJet(1:4,1))
    pT_jet(2) = get_PT(MomJet(1:4,2))
    pT_jet(3) = get_PT(MomJet(1:4,3))

    pT_lepM = get_PT(Mom(1:4,lepM))
    pT_lepP = get_PT(Mom(1:4,lepP))
    eta_lepM = get_ETA(Mom(1:4,lepM))
    eta_lepP = get_ETA(Mom(1:4,lepP))

    HT = pT_Pho + pT_jet(1) + pT_jet(2) + pT_lepM + pT_lepP

    if( dabs(get_Minv(Mom(1:4,LepP)+MomJet(1:4,1))).lt.dabs(get_Minv(Mom(1:4,LepP)+MomJet(1:4,2))) ) then
        m_lb = get_Minv(Mom(1:4,LepP)+MomJet(1:4,1))
    else
        m_lb = get_Minv(Mom(1:4,LepP)+MomJet(1:4,2))
    endif

    Phi_LL = dabs( Get_PHI(Mom(1:4,LepP)) - Get_PHI(Mom(1:4,LepM)) )
    if( Phi_LL.gt.Pi ) Phi_LL=2d0*Pi-Phi_LL

! check cuts
    if(pt_Pho.lt.pT_pho_cut) then
        applyPSCut = .true.
        RETURN
    endif


    ET_miss = get_ET(Mom(1:4,nubar)+Mom(1:4,nu))
    HT = HT + ET_miss
!     if( NJet.eq.3 .and. (pT_jet(3).lt.pT_jet_cut .or. abs(eta_jet(3)).gt.eta_jet_cut) ) then
!         ET_miss = get_ET(Mom(1:4,nubar)+Mom(1:4,nu)+MomJet(1:4,3))
!         HT = HT + ET_miss
!     else
!         ET_miss = get_ET(Mom(1:4,nubar)+Mom(1:4,nu))
!         HT = HT + ET_miss
!         if(NJet.eq.3) HT = HT + pT_jet(3)
!     endif

    if( R_PlepM.lt.Rsep_Plep .or. R_PlepP.lt.Rsep_Plep ) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_jet(1).lt.pT_bjet_cut .or. pT_jet(2).lt.pT_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( abs(eta_jet(1)).gt.eta_bjet_cut .or. abs(eta_jet(2)).gt.eta_bjet_cut) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_lepM.lt.pT_lep_cut .OR. pT_lepP.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lepM).gt.eta_lep_cut .or. abs(eta_lepP).gt.eta_lep_cut) then
        applyPSCut = .true.
        RETURN
    endif

    if( ET_miss.lt.pT_miss_cut ) then
       applyPSCut = .true.
        RETURN
    endif


! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,eta_ATop)
    NBin(4) = WhichBin(4,eta_Top)
    NBin(5) = WhichBin(5,eta_ATop)
    NBin(6) = WhichBin(6,eta_Top)
    NBin(7) = WhichBin(7,pT_Pho)
    NBin(8) = WhichBin(8,eta_Pho)
    NBin(9) = WhichBin(9,pT_lepP)
    NBin(10)= WhichBin(10,eta_lepP)
    NBin(11)= WhichBin(11,ET_miss)
    NBin(12)= WhichBin(12,HT)
    NBin(13)= WhichBin(13,m_lb)
    NBin(14)= WhichBin(14,Phi_LL)


elseif( ObsSet.eq.24 .or. ObsSet.eq.25 ) then! set of observables for ttb+gamma production with semi-lept. at the Tevatron/LHC

! below is a copy of the ObsSet=28 case:
   do k=1,NJet
      pt_jet(k)  = get_PT(MomJet(1:4,k))
      eta_jet(k) = get_eta(MomJet(1:4,k))
      R_Pj(k)    = get_R(MomJet(1:4,k),Mom(1:4,pho))
      R_Lj(k)    = get_R(MomJet(1:4,k),Mom(1:4,L))
   enddo

!   request two separated b-jets
    if( .not.(any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif

!   request b-jets to be outside the Frixione cone
    if( R_Pj(1).lt.Rsep_Pbj .or. R_Pj(2).lt.Rsep_Pbj ) then
        applyPSCut = .true.
        RETURN
    endif

!   check if b-jets pass cuts
    if(  pT_jet(1).lt.pT_bjet_cut .or. abs(eta_jet(1)).gt.eta_bjet_cut .or. R_Lj(1).lt.Rsep_LepJet  ) then
        applyPSCut = .true.
        RETURN
    endif
    if( pT_jet(2).lt.pT_bjet_cut .or. abs(eta_jet(2)).gt.eta_bjet_cut .or. R_Lj(2).lt.Rsep_LepJet ) then
        applyPSCut = .true.
        RETURN
    endif
    HT = pT_jet(1) + pT_jet(2)


!   determine observable light jets
    NObsJet = 2! these are the b-jets
    do k=3,NJet
        if( R_Pj(k).gt.Rsep_Pj .and. R_Lj(k).gt.Rsep_LepJet .and. pT_jet(k).gt.pT_jet_cut .and. abs(eta_jet(k)).lt.eta_jet_cut ) then! count jets outside Frixione cone and outside beam pipe
            NObsJet = NObsJet +1
            if( k.ne.NObsJet ) MomJet(1:4,NObsJet) = MomJet(1:4,k)
            HT = HT + pT_jet(k)
        endif
    enddo

    NObsJet_Tree = 4! request two b-jets and at least two light jets
    if( NObsJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif


    pT_ATop = get_PT(Mom(1:4,tbar))
    pT_Top  = get_PT(Mom(1:4,t))

    eta_ATop = get_ETA(Mom(1:4,tbar))
    eta_Top  = get_ETA(Mom(1:4,t))

    pT_Pho  = get_PT(Mom(1:4,pho))
    eta_Pho = get_ETA(Mom(1:4,pho))

    R_PlepP  = get_R(Mom(1:4,pho),Mom(1:4,L))

    Phi_LP = dabs( Get_PHI(Mom(1:4,pho)) - Get_PHI(Mom(1:4,L)) )
    if( Phi_LP.gt.Pi ) Phi_LP=2d0*Pi-Phi_LP

    pT_lepP  = get_PT(Mom(1:4,L))
    eta_lepP = get_ETA(Mom(1:4,L))

!     MomObs(1:4) = MomObs(1:4) + Mom(1:4,pho) + Mom(1:4,L)
!     MomMiss(1:4) = Mom(1:4,inLeft) + Mom(1:4,inRight) - MomObs(1:4)
    MomMiss(1:4) = Mom(1:4,N)
    ET_miss  = get_ET(MomMiss(1:4))
    HT = HT + pT_lepP + pT_Pho + ET_miss


    m_lb = get_Minv(Mom(1:4,L)+MomJet(1:4,1))    ! these are the pT-hardest b-jets
    Rphobjet = get_R(Mom(1:4,pho),MomJet(1:4,1))


! check cuts
    if(pt_Pho.lt.pT_pho_cut) then
        applyPSCut = .true.
        RETURN
    endif

    if(abs(eta_Pho).gt.eta_pho_cut) then
        applyPSCut = .true.
        RETURN
    endif

    if(abs(eta_Pho).gt.1.37d0 .and. abs(eta_Pho).lt.1.52d0 ) then! this is the crack in the detector
        applyPSCut = .true.
        RETURN
    endif

    if( pT_lepP.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lepP).gt.eta_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( R_PlepP.lt.Rsep_Plep ) then
        applyPSCut = .true.
        RETURN
    endif

    mT_lp = get_MT(Mom(1:4,L),MomMiss(1:4))! this is the transverse W mass

    if( ET_miss.lt.pT_miss_cut .or. ET_miss+mT_lp.lt.60*GeV  ) then
        applyPSCut = .true.
        RETURN
    endif




!     print *, "old def",mt_lp
!     
!     mT_lp = dsqrt( 2d0*pT_lepP*ET_miss*(1d0 - (Mom(2,L)*MomMiss(2)+Mom(3,L)*MomMiss(3))/dsqrt(Mom(2,L)**2+Mom(3,L)**2)/dsqrt(MomMiss(2)**2+MomMiss(3)**2) )   )
!     print *, "new def",mt_lp
! 
!     Phi_LL = dabs( Get_PHI(Mom(1:4,L)) - Get_PHI(MomMiss(1:4)) )
!     if( Phi_LL.gt.Pi ) Phi_LL=2d0*Pi-Phi_LL
!     print *, dcos(Phi_LL),(Mom(2,L)*MomMiss(2)+Mom(3,L)*MomMiss(3))/dsqrt(Mom(2,L)**2+Mom(3,L)**2)/dsqrt(MomMiss(2)**2+MomMiss(3)**2)
!     pause


! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,eta_ATop)
    NBin(3) = WhichBin(3,pT_Top)
    NBin(4) = WhichBin(4,eta_Top)
    NBin(5) = WhichBin(5,eta_ATop)
    NBin(6) = WhichBin(6,eta_Top)
    NBin(7) = WhichBin(7,pT_Pho)
    NBin(8) = WhichBin(8,eta_Pho)
    NBin(9) = WhichBin(9,pT_lepP)
    NBin(10)= WhichBin(10,eta_lepP)
    NBin(11) = WhichBin(11,ET_miss)
    NBin(12) = WhichBin(12,HT)
    NBin(13) = WhichBin(13,Rphobjet)
    NBin(14) = WhichBin(14,m_lb)
    NBin(15) = WhichBin(15,Phi_LP)




elseif( ObsSet.eq.26 .or. ObsSet.eq.27 .or. ObsSet.eq.28 ) then! set of observables for ttb+gamma production with semi-lept. decays(hadr.Atop, lept.top decay) at the Tevatron/LHC
                                                               ! ObsSet 26,27 include suppression cuts for photon radiation from top decay

   do k=1,NJet
      pt_jet(k)  = get_PT(MomJet(1:4,k))
      eta_jet(k) = get_eta(MomJet(1:4,k))
      R_Pj(k)    = get_R(MomJet(1:4,k),Mom(1:4,pho))
   enddo

!   request two separated b-jets
    if( .not.(any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif

!   request b-jets to be outside the Frixione cone
    if( R_Pj(1).lt.Rsep_Pj .or. R_Pj(2).lt.Rsep_Pj ) then
        applyPSCut = .true.
        RETURN
    endif

!   check if b-jets pass cuts
    if(  pT_jet(1).lt.pT_bjet_cut .or. abs(eta_jet(1)).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( pT_jet(2).lt.pT_bjet_cut .or. abs(eta_jet(2)).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    HT = pT_jet(1) + pT_jet(2)


!   determine observable light jets
    NObsJet = 2! these are the b-jets
    do k=3,NJet
        if( R_Pj(k).gt.Rsep_Pj .and. pT_jet(k).gt.pT_jet_cut .and. abs(eta_jet(k)).lt.eta_jet_cut ) then! count jets outside Frixione cone and outside beam pipe
            NObsJet = NObsJet +1
            if( k.ne.NObsJet ) MomJet(1:4,NObsJet) = MomJet(1:4,k)
            HT = HT + pT_jet(k)
        endif
    enddo

    NObsJet_Tree = 4! request two b-jets and at least two light jets
    if( NObsJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif


    pT_ATop = get_PT(Mom(1:4,tbar))
    pT_Top  = get_PT(Mom(1:4,t))

    eta_ATop = get_ETA(Mom(1:4,tbar))
    eta_Top  = get_ETA(Mom(1:4,t))

    pT_Pho  = get_PT(Mom(1:4,pho))
    eta_Pho = get_ETA(Mom(1:4,pho))

    R_PlepP  = get_R(Mom(1:4,pho),Mom(1:4,L))

    Phi_LP = dabs( Get_PHI(Mom(1:4,pho)) - Get_PHI(Mom(1:4,L)) )
    if( Phi_LP.gt.Pi ) Phi_LP=2d0*Pi-Phi_LP

    pT_lepP  = get_PT(Mom(1:4,L))
    eta_lepP = get_ETA(Mom(1:4,L))

!     MomObs(1:4) = MomObs(1:4) + Mom(1:4,pho) + Mom(1:4,L)
!     MomMiss(1:4) = Mom(1:4,inLeft) + Mom(1:4,inRight) - MomObs(1:4)
    MomMiss(1:4) = Mom(1:4,N)
    ET_miss  = get_ET(MomMiss(1:4))
    HT = HT + pT_lepP + pT_Pho + ET_miss


    m_lb = get_Minv(Mom(1:4,L)+MomJet(1:4,1))    ! these are the pT-hardest b-jets
    Rphobjet = get_R(Mom(1:4,pho),MomJet(1:4,1))


!     if( dabs(get_Minv(Mom(1:4,pho)+MomJet(1:4,1))).lt.dabs(get_Minv(Mom(1:4,pho)+MomJet(1:4,2))) ) then
!         m_lb = get_Minv(Mom(1:4,L)+MomJet(1:4,1))
!         Rphobjet = get_R(Mom(1:4,pho),MomJet(1:4,1))
!     else
!         m_lb = get_Minv(Mom(1:4,L)+MomJet(1:4,2))
!         Rphobjet = get_R(Mom(1:4,pho),MomJet(1:4,2))
!     endif

!     Rphobjet = get_R(Mom(1:4,pho),MomJet(1:4,2))!  min.pT b-jet
!     m_lb = get_Minv(Mom(1:4,L)+MomJet(1:4,1))
!     if( get_R(Mom(1:4,pho),MomJet(1:4,1)).lt.get_R(Mom(1:4,pho),MomJet(1:4,2)) ) then
!         Rphobjet = get_R(Mom(1:4,pho),MomJet(1:4,1))
!     else
!         Rphobjet = get_R(Mom(1:4,pho),MomJet(1:4,2))
!     endif


! check cuts
    if(pt_Pho.lt.pT_pho_cut) then
        applyPSCut = .true.
        RETURN
    endif

    if(abs(eta_Pho).gt.eta_pho_cut) then
        applyPSCut = .true.
        RETURN
    endif

    if( R_PlepP.lt.Rsep_Plep ) then
        applyPSCut = .true.
        RETURN
    endif

    if( HT.lt.HT_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_lepP.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lepP).gt.eta_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( ET_miss.lt.pT_miss_cut ) then
        applyPSCut = .true.
        RETURN
    endif



IF( OBSSET.NE.28 ) THEN!   these are the cuts to suppress photon radiation from top quarks and W bosons
    m_jj = get_Minv(MomJet(1:4,3)+MomJet(1:4,4))! 3 4 is the pair
    if( dabs(get_Minv(Mom(1:4,L)+MomJet(1:4,1))).lt.dabs(get_Minv(Mom(1:4,L)+MomJet(1:4,2))) ) then
        mTblP = get_MT(MomJet(1:4,1)+Mom(1:4,pho)+Mom(1:4,L),MomMiss(1:4))
        m_jjb = get_mInv(MomJet(1:4,2)+MomJet(1:4,3)+MomJet(1:4,4))
    else
        mTblP = get_MT(MomJet(1:4,2)+Mom(1:4,pho)+Mom(1:4,L),MomMiss(1:4))
        m_jjb = get_mInv(MomJet(1:4,1)+MomJet(1:4,3)+MomJet(1:4,4))
    endif

    if(NObsJet.eq.5) then
         if( dabs(get_Minv(MomJet(1:4,3)+MomJet(1:4,5))-m_W) .lt. dabs(m_jj-m_W) ) then
              m_jj = get_Minv(MomJet(1:4,3)+MomJet(1:4,5))! 3 5 is the pair
              if( dabs(get_Minv(Mom(1:4,L)+MomJet(1:4,1))).lt.dabs(get_Minv(Mom(1:4,L)+MomJet(1:4,2))) ) then
                  m_jjb = get_mInv(MomJet(1:4,2)+MomJet(1:4,3)+MomJet(1:4,5))
              else
                  m_jjb = get_mInv(MomJet(1:4,1)+MomJet(1:4,3)+MomJet(1:4,5))
              endif
         endif
         if( dabs(get_Minv(MomJet(1:4,4)+MomJet(1:4,5))-m_W) .lt. dabs(m_jj-m_W) ) then
              m_jj = get_Minv(MomJet(1:4,4)+MomJet(1:4,5))! 4 5 is the pair
              if( dabs(get_Minv(Mom(1:4,L)+MomJet(1:4,1))).lt.dabs(get_Minv(Mom(1:4,L)+MomJet(1:4,2))) ) then
                  m_jjb = get_mInv(MomJet(1:4,2)+MomJet(1:4,4)+MomJet(1:4,5))
              else
                  m_jjb = get_mInv(MomJet(1:4,1)+MomJet(1:4,4)+MomJet(1:4,5))
              endif
         endif
    endif
    if( m_jjb.lt.160d0*GeV .or. m_jjb.gt.180d0*GeV) then
         applyPSCut = .true.
        RETURN
    endif
    if( (mTblP.lt.180d0*GeV) ) then
        applyPSCut = .true.
        RETURN
    endif


    if( m_jj.lt.70d0*GeV .or. m_jj.gt.90d0*GeV) then
            applyPSCut = .true.
            RETURN
    endif
    mT_lp = get_MT(Mom(1:4,pho)+Mom(1:4,L),MomMiss(1:4))
    if( mT_lp.lt.90d0*GeV ) then
        applyPSCut = .true.
        RETURN
    endif

ENDIF





! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,eta_ATop)
    NBin(3) = WhichBin(3,pT_Top)
    NBin(4) = WhichBin(4,eta_Top)
    NBin(5) = WhichBin(5,eta_ATop)
    NBin(6) = WhichBin(6,eta_Top)
    NBin(7) = WhichBin(7,pT_Pho)
    NBin(8) = WhichBin(8,eta_Pho)
    NBin(9) = WhichBin(9,pT_lepP)
    NBin(10)= WhichBin(10,eta_lepP)
    NBin(11) = WhichBin(11,ET_miss)
    NBin(12) = WhichBin(12,HT)
    NBin(13) = WhichBin(13,Rphobjet)
    NBin(14) = WhichBin(14,m_lb)
    NBin(15) = WhichBin(15,Phi_LP)


elseif( ObsSet.eq.29) then! ttb+photon production without top decays at Tevatron & LHC

    pT_ATop = get_PT(Mom(1:4,tbar))
    pT_Top  = get_PT(Mom(1:4,t))

    eta_ATop = get_ETA(Mom(1:4,tbar))
    eta_Top  = get_ETA(Mom(1:4,t))

    pT_Pho  = get_PT(Mom(1:4,pho))
    eta_Pho = get_ETA(Mom(1:4,pho))

! check cuts
!     if(pt_Pho.lt.pT_pho_cut) then
!         applyPSCut = .true.
!         RETURN
!     endif

!     R_Pj(1)  = get_R(Mom(1:4,pho),Mom(1:4,realp))
!     if( (Process.eq.24 .or. Process.eq.26) .and. R_Pj(1).lt.Rsep_Pj ) then
!        applyPSCut = .true.
!        RETURN
!     endif


!     if( dabs(Mom(1:4,bbar).dot.Mom(1:4,pho))/m_Top**2.lt.1d-3 .or. dabs(Mom(1:4,b).dot.Mom(1:4,pho))/m_Top**2.lt.1d-3) then
!         applyPSCut = .true.
!         RETURN
!     endif


! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,eta_ATop)
    NBin(4) = WhichBin(4,eta_Top)
    NBin(5) = WhichBin(5,pT_Pho)
    NBin(6) = WhichBin(6,eta_Pho)
    NBin(7) = WhichBin(7,eta_ATop)
    NBin(8) = WhichBin(8,eta_Top)


!-------------------------------------------------------
else
  print *, "ObsSet not implemented",ObsSet
  stop
endif


return
END SUBROUTINE





SUBROUTINE Kinematics_TTBARZ(NPlus1PS,Mom,MomOrder,applyPSCut,NBin)
use ModMisc
use ModParameters
implicit none
integer :: NumHadr,NPlus1PS,MomOrder(1:14)
real(8) :: Mom(1:4,1:14),zeros(1:14)
real(8) :: MomJet(1:4,1:7),MomJet_CHECK(1:4,1:7)
real(8) :: MomHadr(1:4,0:8)
real(8) :: MomBoost(1:4),MomMiss(1:4),MomObs(1:4)
logical :: applyPSCut
integer :: NBin(:),PartList(1:7),JetList(1:7),NJet,NObsJet,k,NObsJet_Tree,leptj(1:3),i,j
real(8) :: nZLept,s12,s13,s14,s23,s24,s34
real(8) :: pT_lep(4),ET_miss,PT_miss,pT_ATop,pT_Top,HT,ET_bjet
real(8) :: eta_ATop,eta_Top,eta_lep(1:4),pseudo_Z, pseudo_top, pseudo_tbar
real(8) :: pT_jet(1:7),eta_jet(1:7),eta_sepa,mT_bln(1:2),pT_Z
real(8) :: R_lj(1:5),R_PlepM,pT_lept,ET_lept,mT,dPhiLL
integer :: tbar,t,Zbos,inLeft,inRight,realp,bbar,lepM,nubar,b,lepP,nu,qdn,qbup,qbdn,qup,L,N,Zl,Za,ferm_Z,Aferm_Z,jlabel
! RR
real(8) :: eta_Z,Recon_M1,Recon_M2,mT_inv 
real(8) :: DphiZt,DphiZtbar,Dphittbar,Dphimumu,Dphimuml,Dphimumb1,Dphimumb2
real(8) :: Dphimumj1,Dphimumj2,Dphimupl,Dphimupb1,Dphimupb2,Dphimupj1,Dphimupj2
real(8) :: DRZt,DRZtbar,DRttbar,DRmumu,DRmuml,DRmumb1,DRmumb2
real(8) :: DRmumj1,DRmumj2,DRmupl,DRmupb1,DRmupb2,DRmupj1,DRmupj2,DRjetjet(1:6)


applyPSCut = .false.
if( Process.eq.81 ) return!  return for Z => photon


! momentum ordering
  tbar    = MomOrder(1)
  t       = MomOrder(2)
  Zbos    = MomOrder(3)
  inLeft  = MomOrder(4)
  inRight = MomOrder(5)
  realp   = MomOrder(6)
  bbar    = MomOrder(7)
  lepM    = MomOrder(8)
  nubar   = MomOrder(9)
  b       = MomOrder(10)
  lepP    = MomOrder(11)
  nu      = MomOrder(12)
  ferm_Z  = MomOrder(13)!      fermion from Z decay (lep-, q, nu)
  Aferm_Z = MomOrder(14)! anti-fermion from Z decay (lep+, qbar, nubar)

  qdn    = lepM
  qbup   = nubar
  qbdn   = lepP
  qup    = nu


!DEC$ IF(_CheckMomenta .EQ.1)
   zeros(:) = 0d0
   if(TopDecays.eq.0) then
      zeros(1:4) = Mom(1:4,inLeft)+Mom(1:4,inRight) - Mom(1:4,tbar) - Mom(1:4,t) - Mom(1:4,ZBos)
   else
      zeros(1:4) = Mom(1:4,inLeft)+Mom(1:4,inRight) - Mom(1:4,ferm_Z) - Mom(1:4,Aferm_Z) - Mom(1:4,bbar)-Mom(1:4,lepM)-Mom(1:4,nubar)-Mom(1:4,b)-Mom(1:4,lepP)-Mom(1:4,nu)
   endif
   if( NPlus1PS.eq.1 ) zeros(1:4) = zeros(1:4) - Mom(1:4,realp)
   if( any(abs(zeros(1:4)/Mom(1,inLeft)).gt.1d-8) ) then
      print *, "ERROR: energy-momentum violation in SUBROUTINE Kinematics_TTBARZ(): ",NPlus1PS,zeros(1:4)
      print *, "momenta dump:"
      do k=1,12
         print *,k, Mom(1:4,k)
      enddo
   endif
   zeros(1) = (Mom(1:4,tbar).dot.Mom(1:4,tbar)) - m_Top**2
   zeros(2) = (Mom(1:4,t).dot.Mom(1:4,t)) - m_Top**2
   if (ZDecays .le. 10) then   ! Z onshell
      zeros(3) = (Mom(1:4,ZBos).dot.Mom(1:4,ZBos)) - M_Z**2
   else ! Z off shell, so this check is meaningless
      zeros(3)=0d0
   endif
   zeros(4) =  Mom(1:4,bbar).dot.Mom(1:4,bbar)
   zeros(5) =  Mom(1:4,lepM).dot.Mom(1:4,lepM)
   zeros(6) =  Mom(1:4,nubar).dot.Mom(1:4,nubar)
   zeros(7) =  Mom(1:4,b).dot.Mom(1:4,b)
   zeros(8) =  Mom(1:4,lepP).dot.Mom(1:4,lepP)
   zeros(9) =  Mom(1:4,nu).dot.Mom(1:4,nu)
   zeros(10) = (Mom(1:4,ferm_Z).dot.Mom(1:4,ferm_Z))
   zeros(11) = (Mom(1:4,Aferm_Z).dot.Mom(1:4,Aferm_Z))


   if( NPlus1PS.eq.1 ) zeros(12)=  Mom(1:4,realp).dot.Mom(1:4,realp)
   if( TopDecays.eq.0 .and. any(abs(zeros(1:3)/Mom(1,inLeft)**2).gt.1d-8) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TTBARZ(): ",zeros(1:3)
      print *, Mom(1:4,1:3)
   endif
   if( TopDecays.ne.0 .and. any(abs(zeros(1:12)/Mom(1,inLeft)**2).gt.1d-8) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TTBARZ(): ",zeros(1:12)
      print *, "momenta dump:"
      print *, Mom(1:4,1:12)
   endif
!DEC$ ENDIF


!!! RR 2013/03/16 : havent made any changes for Z decay from here on - it looks like it shoul just be top decay though (other than some cuts on Z leptons -- much later)

NBin(1:NumHistograms) = 0
MomHadr(1:4,1:7) = 0d0
PartList(1:7)=(/1,2,3,4,5,6,7/)
MomMiss(1:4) = 0d0
MomObs(1:4)  = 0d0

! separating momenta into hadron momenta and lepton momenta to which various procedures can be applied
!-------------------------------------------------------

IF( TOPDECAYS.EQ.0 ) THEN  ! no top decays
    if(NPlus1PS.eq.0) then
        NumHadr = 0
    else
        NumHadr = 1
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
    endif

ELSEIF( TOPDECAYS.EQ.1 .OR. TOPDECAYS.eq.-1 ) THEN  ! di-leptonic decay
    MomHadr(1:4,1) = Mom(1:4,bbar)
    MomHadr(1:4,2) = Mom(1:4,b)     ! Bot
    if(NPlus1PS.eq.0) then
        NumHadr = 2
    else
        NumHadr = 3
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
    endif

ELSEIF( TOPDECAYS.EQ.3 ) THEN  ! lept. Atop, hadr. top decay
   MomHadr(1:4,1) = Mom(1:4,bbar)
   MomHadr(1:4,2) = Mom(1:4,b)
   MomHadr(1:4,3) = Mom(1:4,qbdn)
   MomHadr(1:4,4) = Mom(1:4,qup)
   L = LepM
   N = nubar
   if(NPlus1PS.eq.0) then
        NumHadr = 4
   else
        NumHadr = 5
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)   ! q/qbar/glu
   endif

ELSEIF( TOPDECAYS.EQ.4 ) THEN  ! hadr. Atop, lept. top decay
   MomHadr(1:4,1) = Mom(1:4,bbar)
   MomHadr(1:4,2) = Mom(1:4,b)
   MomHadr(1:4,3) = Mom(1:4,qbup)
   MomHadr(1:4,4) = Mom(1:4,qdn)
   L = LepP
   N = nu
   if(NPlus1PS.eq.0) then
        NumHadr = 4
   else
        NumHadr = 5
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
   endif

ELSE
  call Error("this decay is not yet implemented")
ENDIF





!---------------------- kT jet algorithm ---------------------------------
! important: b-quarks need to be at the first two positions of MomHadr(:,:)
! recombination of MomHadr(1:4,i) and MomHadr(1:4,j) results in MomJet(1:4,i) with i<j.
! but take care: some MomJet(1:4,i) can be zero, this is why we apply pt_order

    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr)) ! hard protojets in beam pipe are counted as jets
    call pT_order(2,MomJet(1:4,1:2))
    call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))


! call SwitchEnergyComponent(MomHadr(1:4,1:NumHadr))
! call fastjetppgenkt(MomHadr(1:4,1:NumHadr),NumHadr,Rsep_jet,-1d0,MomJet(1:4,1:NumHadr),NJet)! 4th argument:  +1d0=kT,  -1d0=anti-kT,   0d0=CA
! call SwitchEnergyComponentBack(MomJet_CHECK(1:4,1:NumHadr))





!------------------------ cuts and binning --------------------------------
if( ObsSet.eq.51) then! ttb+Z production without top decays at Tevatron & LHC

    pT_ATop = get_PT(Mom(1:4,tbar))
    pT_Top  = get_PT(Mom(1:4,t))


! binning
    NBin(1) = WhichBin(1,pT_Top)




elseif( ObsSet.eq.52 .or. ObsSet.eq.55 ) then! set of observables for ttb+Z ( di-lept. ttbar decays and di-lept. Z decay )


!     s12 = dabs( get_MInv(Mom(1:4,lepM)+Mom(1:4,lepP))-M_Z )
!     s13 = dabs( get_MInv(Mom(1:4,lepM)+Mom(1:4,ferm_Z))-M_Z )
!     s14 = dabs( get_MInv(Mom(1:4,lepM)+Mom(1:4,Aferm_Z))-M_Z )
!     s23 = dabs( get_MInv(Mom(1:4,lepP)+Mom(1:4,ferm_Z))-M_Z )
!     s24 = dabs( get_MInv(Mom(1:4,lepP)+Mom(1:4,Aferm_Z))-M_Z )
!     s34 = dabs( get_MInv(Mom(1:4,ferm_Z)+Mom(1:4,Aferm_Z))-M_Z )
!     nZLept=minloc((/s12,s13,s14,s23,s24,s34/),1)


    pT_jet(1) = get_PT(Mom(1:4,b))
    pT_jet(2) = get_PT(Mom(1:4,bbar))
    eta_jet(1) = get_ETA(Mom(1:4,b))
    eta_jet(2) = get_ETA(Mom(1:4,bbar))

    eta_Z = get_ETA(Mom(1:4,ferm_Z)+Mom(1:4,Aferm_Z))
    pT_top = get_PT(Mom(1:4,t))
    pT_atop = get_PT(Mom(1:4,tbar))
    eta_top = get_ETA(Mom(1:4,t))
    eta_atop = get_ETA(Mom(1:4,tbar))

    pseudo_Z=get_pseudoETA(Mom(1:4,Zbos))
    pseudo_top=get_pseudoETA(Mom(1:4,t))
    pseudo_tbar=get_pseudoETA(Mom(1:4,tbar))

    pT_Lep(1)  = get_PT(Mom(1:4,ferm_Z))
    pT_Lep(2)  = get_PT(Mom(1:4,Aferm_Z))
    pT_Lep(3)  = get_PT(Mom(1:4,LepP))
    pT_Lep(4)  = get_PT(Mom(1:4,LepM))
    eta_Lep(1)  = get_ETA(Mom(1:4,ferm_Z))
    eta_Lep(2)  = get_ETA(Mom(1:4,Aferm_Z))
    eta_Lep(3) = get_ETA(Mom(1:4,LepP))
    eta_Lep(4) = get_ETA(Mom(1:4,LepM))

    pT_miss = get_PT(Mom(1:4,nu)+Mom(1:4,nubar))

    pT_Z  = get_PT(Mom(1:4,ferm_Z)+Mom(1:4,Aferm_Z))

    DphiLL = dabs( Get_PHI(Mom(1:4,ferm_Z)) - Get_PHI(Mom(1:4,Aferm_Z))  )
    if( DphiLL.gt.Pi ) DphiLL=2d0*Pi-DphiLL


! check cuts
    if( pT_jet(1).lt.pT_bjet_cut .OR. pT_jet(2).lt.pT_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_jet(1)).gt.eta_bjet_cut .OR. abs(eta_jet(2)).gt.eta_bjet_cut) then
       applyPSCut = .true.
        RETURN
    endif

    if( pT_lep(1).lt.pT_lep_cut .OR. pT_lep(2).lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lep(1)).gt.eta_lep_cut .OR. abs(eta_lep(2)).gt.eta_lep_cut) then
       applyPSCut = .true.
        RETURN
    endif

    if( pT_miss.lt.pT_miss_cut ) then
       applyPSCut = .true.
        RETURN
    endif



! binning
    NBin(1) = WhichBin(1,pT_Lep(1))
    NBin(2) = WhichBin(2,pT_Lep(2))
    NBin(3) = WhichBin(3,pT_Z)
    NBin(4) = WhichBin(4,pT_top)
    NBin(5) = WhichBin(5,pT_atop)
    NBin(6) = WhichBin(6,pT_jet(1))
    NBin(7) = WhichBin(7,pT_jet(2))
    NBin(8) = WhichBin(8,pT_Lep(3))
    NBin(9) = WhichBin(9,pT_Lep(4))
    NBin(10) = WhichBin(10,pT_miss)
    NBin(11) = WhichBin(11,eta_Lep(1))
    NBin(12) = WhichBin(12,eta_Lep(2))
    NBin(13) = WhichBin(13,eta_Z)
    NBin(14) = WhichBin(14,eta_top)
    NBin(15) = WhichBin(15,eta_atop)
    NBin(16) = WhichBin(16,eta_jet(1))
    NBin(17) = WhichBin(17,eta_jet(2))
    NBin(18) = WhichBin(18,eta_Lep(3))
    NBin(19) = WhichBin(19,eta_Lep(4))
    NBin(20) = WhichBin(20,dphill)
    NBin(21) = WhichBin(21,pseudo_Z)
    NBin(22) = WhichBin(22,pseudo_top)
    NBin(23) = WhichBin(23,pseudo_tbar)





elseif( ObsSet.eq.53 .or. ObsSet.eq.56 ) then! set of observables for ttb+Z ( di-lept. ttbar decays and di-lept. Z decay )

!     Rsep_jetlep = 0.4d0

! request at least four jets where two are b-jets
    NObsJet_Tree = 4
    if( .not.(NJet.ge.NObsJet_Tree .and. any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif


    pT_jet(1) = get_PT(MomJet(1:4,1))
    pT_jet(2) = get_PT(MomJet(1:4,2))
    pT_jet(3) = get_PT(MomJet(1:4,3))
    pT_jet(4) = get_PT(MomJet(1:4,4))
    eta_jet(1) = get_ETA(MomJet(1:4,1))
    eta_jet(2) = get_ETA(MomJet(1:4,2))
    eta_jet(3) = get_ETA(MomJet(1:4,3))
    eta_jet(4) = get_ETA(MomJet(1:4,4))


    pT_Lep(1)  = get_PT(Mom(1:4,LepP))
    eta_Lep(1) = get_ETA(Mom(1:4,LepP))

    pT_Lep(2)  = get_PT(Mom(1:4,ferm_Z))
    eta_Lep(2) = get_ETA(Mom(1:4,ferm_Z))

    pT_Lep(3)  = get_PT(Mom(1:4,Aferm_Z))
    eta_Lep(3) = get_ETA(Mom(1:4,Aferm_Z))

    pT_miss = get_PT(Mom(1:4,nu))

    pT_Z  = get_PT(Mom(1:4,ferm_Z)+Mom(1:4,Aferm_Z))

    DphiLL = dabs( Get_PHI(Mom(1:4,ferm_Z)) - Get_PHI(Mom(1:4,Aferm_Z))  )
    if( DphiLL.gt.Pi ) DphiLL=2d0*Pi-DphiLL

! these are only needed for RR plots
    eta_Z = get_ETA(Mom(1:4,ferm_Z)+Mom(1:4,Aferm_Z))
    pT_top = get_PT(Mom(1:4,t))
    pT_atop = get_PT(Mom(1:4,tbar))
    eta_top = get_PT(Mom(1:4,t))
    eta_atop = get_PT(Mom(1:4,tbar))
    Recon_M1 = Get_MInv(Mom(1:4,lepP)+Mom(1:4,nu)+MomJet(1:4,1))
    Recon_M2 = Get_MInv(Mom(1:4,lepP)+Mom(1:4,nu)+MomJet(1:4,2))
    if ( abs(m_Top-Recon_M1) .le. abs(m_Top-Recon_M2)) then
       mT_inv = Get_MT(Mom(1:4,lepP)+MomJet(1:4,1),Mom(1:4,nu))
    else
       mT_inv = Get_MT(Mom(1:4,lepP)+MomJet(1:4,2),Mom(1:4,nu))
    endif


! delta phi
    DphiZt = dabs( Get_PHI(Mom(1:4,Zbos)) - Get_PHI(Mom(1:4,t))  )
    if( DphiZt.gt.Pi ) DphiZt=2d0*Pi-DphiZt

    DphiZtbar = dabs( Get_PHI(Mom(1:4,Zbos)) - Get_PHI(Mom(1:4,tbar))  )
    if( DphiZtbar.gt.Pi ) DphiZtbar=2d0*Pi-DphiZtbar

    Dphittbar = dabs( Get_PHI(Mom(1:4,t)) - Get_PHI(Mom(1:4,tbar))  )
    if( Dphittbar.gt.Pi ) Dphittbar=2d0*Pi-Dphittbar

    Dphimumu = dabs( Get_PHI(Mom(1:4,ferm_Z)) - Get_PHI(Mom(1:4,aferm_Z))  )
    if( Dphimumu.gt.Pi ) Dphimumu=2d0*Pi-Dphimumu

    Dphimuml = dabs( Get_PHI(Mom(1:4,ferm_Z)) - Get_PHI(Mom(1:4,lepP))  )
    if( Dphimuml.gt.Pi ) Dphimuml = 2d0*Pi-Dphimuml

    Dphimumb1 = dabs( Get_PHI(Mom(1:4,ferm_Z)) - Get_PHI(MomJet(1:4,1))  )
    if( Dphimumb1.gt.Pi ) Dphimumb1 = 2d0*Pi-Dphimumb1

    Dphimumb2 = dabs( Get_PHI(Mom(1:4,ferm_Z)) - Get_PHI(MomJet(1:4,2))  )
    if( Dphimumb2.gt.Pi ) Dphimumb2 = 2d0*Pi-Dphimumb2

    Dphimumj1 = dabs( Get_PHI(Mom(1:4,ferm_Z)) - Get_PHI(MomJet(1:4,3))  )
    if( Dphimumj1.gt.Pi ) Dphimumj1 = 2d0*Pi-Dphimumj1

    Dphimumj2 = dabs( Get_PHI(Mom(1:4,ferm_Z)) - Get_PHI(MomJet(1:4,4))  )
    if( Dphimumj2.gt.Pi ) Dphimumj2 = 2d0*Pi-Dphimumj2

    Dphimupl = dabs( Get_PHI(Mom(1:4,aferm_Z)) - Get_PHI(Mom(1:4,lepP))  )
    if( Dphimupl.gt.Pi ) Dphimupl = 2d0*Pi-Dphimupl

    Dphimupb1 = dabs( Get_PHI(Mom(1:4,aferm_Z)) - Get_PHI(MomJet(1:4,1))  )
    if( Dphimupb1.gt.Pi ) Dphimupb1 = 2d0*Pi-Dphimupb1

    Dphimupb2 = dabs( Get_PHI(Mom(1:4,aferm_Z)) - Get_PHI(MomJet(1:4,2))  )
    if( Dphimupb2.gt.Pi ) Dphimupb2 = 2d0*Pi-Dphimupb2

    Dphimupj1 = dabs( Get_PHI(Mom(1:4,aferm_Z)) - Get_PHI(MomJet(1:4,3))  )
    if( Dphimupj1.gt.Pi ) Dphimupj1 = 2d0*Pi-Dphimupj1

    Dphimupj2 = dabs( Get_PHI(Mom(1:4,aferm_Z)) - Get_PHI(MomJet(1:4,4))  )
    if( Dphimupj2.gt.Pi ) Dphimupj2 = 2d0*Pi-Dphimupj2

    ! delta R
    DRZt    = Get_R(Mom(1:4,Zbos),Mom(1:4,t))
    DRZtbar = Get_R(Mom(1:4,Zbos),Mom(1:4,tbar))
    DRttbar = Get_R(Mom(1:4,t),Mom(1:4,tbar))  
    DRmumu  = Get_R(Mom(1:4,ferm_Z),Mom(1:4,aferm_Z))
    DRmuml  = Get_R(Mom(1:4,ferm_Z),Mom(1:4,lepP))
    DRmumb1 = Get_R(Mom(1:4,ferm_Z), MomJet(1:4,1))  
    DRmumb2 = Get_R(Mom(1:4,ferm_Z), MomJet(1:4,2))  
    DRmumj1 = Get_R(Mom(1:4,ferm_Z), MomJet(1:4,3))  
    DRmumj2 = Get_R(Mom(1:4,ferm_Z), MomJet(1:4,4))  
    DRmupl  = Get_R(Mom(1:4,aferm_Z), Mom(1:4,lepP)) 
    DRmupb1 = Get_R(Mom(1:4,aferm_Z), MomJet(1:4,1)) 
    DRmupb2 = Get_R(Mom(1:4,aferm_Z), MomJet(1:4,2)) 
    DRmupj1 = Get_R(Mom(1:4,aferm_Z), MomJet(1:4,3)) 
    DRmupj2 = Get_R(Mom(1:4,aferm_Z), MomJet(1:4,4)) 

    jlabel=0
    do i=1,4
       do j=1,4
          if (i.ne.j) then
             jlabel=jlabel+1
             DRjetjet(jlabel)=Get_R(MomJet(1:4,i),MomJet(1:4,j))
          endif
       enddo
    enddo
    
    

! check cuts
    if( pT_jet(1).lt.pT_bjet_cut .OR. pT_jet(2).lt.pT_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_jet(1)).gt.eta_bjet_cut .OR. abs(eta_jet(2)).gt.eta_bjet_cut) then
       applyPSCut = .true.
        RETURN
    endif

    if( pT_jet(3).lt.pT_jet_cut .OR. pT_jet(4).lt.pT_jet_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_jet(3)).gt.eta_jet_cut .OR. abs(eta_jet(4)).gt.eta_jet_cut) then
       applyPSCut = .true.
        RETURN
    endif


    if( pT_lep(1).lt.pT_lep_cut .OR. pT_lep(2).lt.pT_lep_cut .OR. pT_lep(3).lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lep(1)).gt.eta_lep_cut .OR. abs(eta_lep(2)).gt.eta_lep_cut  .OR. abs(eta_lep(3)).gt.eta_lep_cut ) then
       applyPSCut = .true.
        RETURN
    endif

    if( pT_miss.lt.pT_miss_cut ) then
       applyPSCut = .true.
        RETURN
    endif

    leptj = (/ LepP,ferm_Z,Aferm_Z /)
    do i=1,4
    do j=1,3
      if( get_R( MomJet(1:4,i),Mom(1:4,leptj(j)) ) .lt. Rsep_jetlep ) then 
        applyPSCut = .true.
        RETURN
      endif
    enddo
    enddo

! binning
!    NBin(1) = WhichBin(1,pT_Lep(1))
!    NBin(2) = WhichBin(2,pT_Z)
!    NBin(3) = WhichBin(3,dphill)

    NBin(1) = WhichBin(1,pT_Lep(1))
    NBin(2) = WhichBin(2,pT_Lep(2))
    NBin(3) = WhichBin(3,pT_Lep(3))
    NBin(4) = WhichBin(4,pT_jet(1))
    NBin(5) = WhichBin(5,pT_jet(2))
    NBin(6) = WhichBin(6,pT_jet(3))
    NBin(7) = WhichBin(7,pT_jet(4))
    NBin(8) = WhichBin(8,pT_miss)
    NBin(9) = WhichBin(9,eta_Lep(1))
    NBin(10) = WhichBin(10,eta_Lep(2))
    NBin(11) = WhichBin(11,eta_Lep(3))
    NBin(12) = WhichBin(12,eta_jet(1))
    NBin(13) = WhichBin(13,eta_jet(2))
    NBin(14) = WhichBin(14,eta_jet(3))
    NBin(15) = WhichBin(15,eta_jet(4))
    NBin(16) = WhichBin(16,pT_Z)
    NBin(17) = WhichBin(17,eta_Z)
    NBin(18) = WhichBin(18,pT_top)
    NBin(19) = WhichBin(19,eta_top)
    NBin(20) = WhichBin(20,pT_atop)
    NBin(21) = WhichBin(21,eta_atop)
    NBin(22) = WhichBin(22,mT_inv)
!
    NBin(23) = WhichBin(23,DPhiZt)
    NBin(24) = WhichBin(24,DPhiZtbar)
    NBin(25) = WhichBin(25,DPhittbar)
    NBin(26) = WhichBin(26,DPhimumu)
    NBin(27) = WhichBin(27,DPhimuml)
    NBin(28) = WhichBin(28,DPhimumb1)
    NBin(29) = WhichBin(29,DPhimumb2)
    NBin(30) = WhichBin(30,DPhimumj1)
    NBin(31) = WhichBin(31,DPhimumj2)
    NBin(32) = WhichBin(32,DPhimupl)
    NBin(33) = WhichBin(33,DPhimupb1)
    NBin(34) = WhichBin(34,DPhimupb2)
    NBin(35) = WhichBin(35,DPhimupj1)
    NBin(36) = WhichBin(36,DPhimupj2)
    NBin(37) = WhichBin(37,DRZt)
    NBin(38) = WhichBin(38,DRZtbar)
    NBin(39) = WhichBin(39,DRttbar)
    NBin(40) = WhichBin(40,DRmumu)
    NBin(41) = WhichBin(41,DRmuml)
    NBin(42) = WhichBin(42,DRmumb1)
    NBin(43) = WhichBin(43,DRmumb2)
    NBin(44) = WhichBin(44,DRmumj1)
    NBin(45) = WhichBin(45,DRmumj2)
    NBin(46) = WhichBin(46,DRmupl)
    NBin(47) = WhichBin(47,DRmupb1)
    NBin(48) = WhichBin(48,DRmupb2)
    NBin(49) = WhichBin(49,DRjetjet(1))
    NBin(50) = WhichBin(50,DRjetjet(2))
    NBin(51) = WhichBin(51,DRjetjet(3))
    NBin(52) = WhichBin(52,DRjetjet(4))
    NBin(53) = WhichBin(53,DRjetjet(5))
    NBin(54) = WhichBin(54,DRjetjet(6))


!-------------------------------------------------------
else
  print *, "ObsSet not implemented",ObsSet
  stop
endif


return
END SUBROUTINE









SUBROUTINE Kinematics_TTBAR(NPlus1PS,MomExt,MomDK,applyPSCut,NBin,xJPsiFrag)
use ModMisc
use ModParameters
use ModSpinCorrel
implicit none
real(8), optional :: xJPsiFrag
integer :: NumHadr
real(8) :: MomExt(:,:),MomDK(:,:),MomJet(1:4,1:8)
real(8) :: MomLepP(1:4),MomLepM(1:4),MomBoost(1:4),MomZ(1:4),MomTops(1:4,1:2),MomJPsi(1:4),MomMiss(1:4),MomLepTR(1:4)
logical :: applyPSCut,NPlus1PS
integer :: NBin(:),PartList(1:7),JetList(1:7),NJet,iJet,NObsJet_Tree,NObsJet
real(8) :: pT_jet1,pT_jet2,pT_bjet1,pT_bjet2,pT_lepM,pT_lepP,pT_miss,pT_ATop,pT_Top
real(8) :: eta_ATop,eta_Top,eta_lepM,eta_lepP,eta_bjet1,eta_bjet2,R_bb
real(8) :: MInv_lepP_bjet,CosPhi_LepLep,CosPsiT_LepLep,CosPsi_LepLep,Psi_LepLep,DeltaPhi
real(8) :: MInv_Tops,CosPhi_LepPZ,InvM_Lep,CosPhi_LepPlanes,HT,E_bjet1,E_bjet2,ET_miss
real(8) :: MInv_lepP_lepM,pT_lept,ET_lept,pT_x,pT_y,mT,MinvLept,E_lept,m_lb
real(8),save :: EMin=1d12
real(8) :: MomHadr(1:4,1:7),MomLept(1:4,1:4),MomJet_aux(1:4,1:7)
real(8) :: second_pair, first_pair, pT_bjet_min, mtop_reconstr,ptaux,pt_cand(4),pij(4)
real(8) :: ptmax, xij, MomBBar(1:4),diff,cosLeBb,Ebjets,r_sc,mT_lp
integer :: i1,i,j,n,k
real(8) :: MomLeptOrd(1:4,1:2), MomJetOrd(1:4,1:8), M_eff, Minv_Lept

!-- for Zprime background
real(8) :: costheta_star, MomAux1(1:4), MomAux2(1:4), y_Top, y_Atop, m_ttbar, Mttb_cut, dphi_ttbar, costheta_scatter
real(8) :: pt_lep1, pt_lep2, eta_lep1, eta_lep2, dphi_LL, MassAux, MomTopsCMS(1:4,1:2)
real(8) :: MomBeam1(1:4), MomBeam2(1:4), MomBeam(1:4), nx(2:4), ny(2:4), nz(2:4), MomLeptTRF(1:4)
real(8) :: cosPhi, sinPhi, cosPhiBar, sinPhiBar, deltaSin, dPhiMinus, dPhiPlus

! required momentum order: MomExt(:,:): 1=In_left, 2=In_right, 3,4,...=light particles, N-1=ATop, N=Top    ! HERE WAS A BUG: tops need to be at the end
!                          MomDK(:,1:7) : 1=ABot, 2=lep-, 3=ANeu, 4=Bot, 5=lep+, 6=Neu, 7=(Glu)
! MomLept(:,1:4): 1=lep-, 2=ANeu, 3=lep+, 4=Neu


applyPSCut = .false.
NBin(1:NumHistograms) = 0


MomHadr(1:4,1:7) = 0d0
MomLept(1:4,1:4) = 0d0
PartList(1:7)=(/1,2,3,4,5,6,7/)




! separating momenta into hadron momenta and lepton momenta to which various procedures can be applied
MomHadr(1:4,1) = MomDK(1:4,1)  ! ABot
MomHadr(1:4,2) = MomDK(1:4,4)  ! Bot
!-------------------------------------------------------
if( TopDecays.eq.0 ) then  ! no top decays
   NumHadr = 0
   if( NPlus1PS ) then
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
   else
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
   endif
!-------------------------------------------------------
elseif( abs(TopDecays).eq.1 ) then  ! full leptonic decay
  MomLept(1:4,1) = MomDK(1:4,2)   ! lep-
  MomLept(1:4,2) = MomDK(1:4,3)   ! ANeu
  MomLept(1:4,3) = MomDK(1:4,5)   ! lep+
  MomLept(1:4,4) = MomDK(1:4,6)   ! Neu

  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,3) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 3
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,3) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 3
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 2
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.2 ) then  ! full hadronic decay
  MomHadr(1:4,3) = MomDK(1:4,2)   ! qbar
  MomHadr(1:4,4) = MomDK(1:4,3)   ! q
  MomHadr(1:4,5) = MomDK(1:4,5)   ! qbar
  MomHadr(1:4,6) = MomDK(1:4,6)   ! q
  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,7) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 7
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,7) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 7
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 6
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.3 ) then  ! lept. Atop, hadr. top decay
  MomLept(1:4,1) = MomDK(1:4,2)  ! lep-
  MomLept(1:4,2) = MomDK(1:4,3)  ! ANeu
  MomHadr(1:4,3) = MomDK(1:4,5)  ! qbar
  MomHadr(1:4,4) = MomDK(1:4,6)  ! q
  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,5) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,5) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 4
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.4 ) then  ! hadr. Atop, lept. top decay
  MomHadr(1:4,3) = MomDK(1:4,2)  ! qbar
  MomHadr(1:4,4) = MomDK(1:4,3)  ! q
  MomLept(1:4,3) = MomDK(1:4,5)  ! lep+
  MomLept(1:4,4) = MomDK(1:4,6)  ! Neu
  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,5) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,5) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 4
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.5 ) then  ! lept. Atop, hadr. top decay with J/Psi fragmentation

  MomHadr(1:4,1) = (1d0-xJPsiFrag) * MomDK(1:4,1)  ! b-jet remnand
  MomJPsi(1:4)   = xJPsiFrag * MomDK(1:4,1)        ! J/Psi momentum

  MomLept(1:4,1) = MomDK(1:4,2)  ! lep-
  MomLept(1:4,2) = MomDK(1:4,3)  ! ANeu
  MomHadr(1:4,3) = MomDK(1:4,5)  ! qbar
  MomHadr(1:4,4) = MomDK(1:4,6)  ! q
  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,5) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,5) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 4
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.6) then  ! hadr. Atop, lept. top decay with J/Psi fragmentation


! MomHad(1:4,2) is the b-quark momentum because J/Psi and b-remnand jet are always recombined

  MomJPsi(1:4)   = xJPsiFrag * MomDK(1:4,4)        ! J/Psi momentum
  MomHadr(1:4,3) = MomDK(1:4,2)  ! qbar
  MomHadr(1:4,4) = MomDK(1:4,3)  ! q

  MomLept(1:4,3) = MomDK(1:4,5)  ! lep+
  MomLept(1:4,4) = MomDK(1:4,6)  ! Neu
  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,5) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,5) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 4
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif

endif



!---------------------- kT jet algorithm ---------------------------------
    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)

! print *, PartList(1:NumHadr)
! print *, MomJet(1:4,1)
! print *, MomJet(1:4,2)
! print *, MomJet(1:4,3)
! print *, MomJet(1:4,4)
! print *, MomJet(1:4,5)

    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr))

! print *, "NJet",NJet
! print *, JetList(1:NJet)
! print *, MomJet(1:4,1)
! print *, MomJet(1:4,2)
! print *, MomJet(1:4,3)
! print *, MomJet(1:4,4)
! print *, MomJet(1:4,5)
! pause

!--------------------------------------------------------------------------

if( ObsSet.eq.0 .or. ObsSet.eq.1 .or. ObsSet.eq.9) then! set of observables for ttb production without decays at Tevatron & LHC
!  eval kinematic variables
    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))

    eta_ATop = get_ETA(MomTops(1:4,1))
    eta_Top  = get_ETA(MomTops(1:4,2))

!     if( pT_Top.lt.800d0*GeV .or. pT_ATop.lt.800*GeV ) then!   this is for the boosted observable
!         applyPSCut = .true.
!         RETURN
!     endif

! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,eta_ATop)
    NBin(4) = WhichBin(4,eta_Top)
    NBin(5) = WhichBin(5,eta_ATop)
    NBin(6) = WhichBin(6,eta_Top)


!-------------------------------------------------------
elseif( ObsSet.eq.2 .or. ObsSet.eq.3) then! set of observables for ttb production with di-lept. decays at TEV and LHC

! request at least two b-jets
    if( .not.(NJet.ge.2 .and. any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif
   call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))
   call pT_order(2,MomJet(1:4,1:2))


! evaluate kinematic variables
    pT_bjet1 = get_PT(MomJet(1:4,1))
    pT_bjet2 = get_PT(MomJet(1:4,2))
    eta_bjet1 = get_ETA(MomJet(1:4,1))
    eta_bjet2 = get_ETA(MomJet(1:4,2))
    E_bjet1 = MomJet(1,1)
    E_bjet2 = MomJet(1,2)

    MinvLept = Get_MInv(MomLept(1:4,1)+MomLept(1:4,3))
    pT_miss = get_PT(MomLept(1:4,2)+MomLept(1:4,4))
    pT_lept = get_PT(MomLept(1:4,1)+MomLept(1:4,3))
    pT_miss = get_PT(MomLept(1:4,2)+MomLept(1:4,4))
    ET_lept = dsqrt(pT_lept**2 + MinvLept**2)
    ET_miss = get_ET(MomLept(1:4,2)+MomLept(1:4,4))
    pT_x = MomLept(2,1)+MomLept(2,2)+MomLept(2,3)+MomLept(2,4)
    pT_y = MomLept(3,1)+MomLept(3,2)+MomLept(3,3)+MomLept(3,4)
    mT = dsqrt( (ET_lept+ET_miss)**2 - (pT_x**2+pT_y**2) )

    E_lept = MomLept(1,1)+MomLept(1,3)

    pT_lepM = get_PT(MomLept(1:4,1))
    pT_lepP = get_PT(MomLept(1:4,3))

    eta_lepM = get_ETA(MomLept(1:4,1))
    eta_lepP = get_ETA(MomLept(1:4,3))

    MInv_lepP_lepM = get_MInv( MomLept(1:4,1)+MomLept(1:4,3) )

    if( get_MInv(MomLept(1:4,3)+MomJet(1:4,1)).lt.get_MInv( MomLept(1:4,3)+MomJet(1:4,2)) ) then 
        MInv_lepP_bjet = get_MInv( MomLept(1:4,3)+MomJet(1:4,1) )
    else
        MInv_lepP_bjet = get_MInv( MomLept(1:4,3)+MomJet(1:4,2) )
    endif
!MInv_lepP_bjet = get_MInv( MomLept(1:4,3)+MomDK(1:4,4) )

    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))
    eta_ATop = get_ETA(MomTops(1:4,1))
    eta_Top  = get_ETA(MomTops(1:4,2))

    HT = pT_bjet1 + pT_bjet2 + pT_miss + pT_lepM + pT_lepP
    if( NJet.gt.2 ) HT = HT + get_PT(MomJet(1:4,3))

    MInv_Tops = get_MInv( MomTops(1:4,1)+MomTops(1:4,2) )

!     MomLepM(1:4)  = MomLept(1:4,1)
!     MomBoost(1)   =+MomTops(1,1)
!     MomBoost(2:4) =-MomTops(2:4,1)
!     call boost(MomLepM(1:4),MomBoost(1:4),m_top)
!     MomLepP(1:4)  = MomLept(1:4,3)
!     MomBoost(1)   =+MomTops(1,2)
!     MomBoost(2:4) =-MomTops(2:4,2)
!     call boost(MomLepP(1:4),MomBoost(1:4),m_top)
!     CosPhi_LepLep = (MomLepP(2)*MomLepM(2)+MomLepP(3)*MomLepM(3)+MomLepP(4)*MomLepM(4))/MomLepP(1)/MomLepM(1)

!     CosPsiT_LepLep = ( MomLept(2,1)*MomLept(2,3) + MomLept(3,1)*MomLept(3,3) )/dsqrt(MomLept(2,1)**2+MomLept(3,1)**2)/dsqrt(MomLept(2,3)**2+MomLept(3,3)**2)
    Psi_LepLep  = dacos( ( MomLept(2,1)*MomLept(2,3) + MomLept(3,1)*MomLept(3,3) + MomLept(4,1)*MomLept(4,3))/dsqrt(MomLept(2,1)**2+MomLept(3,1)**2+MomLept(4,1)**2)/dsqrt(MomLept(2,3)**2+MomLept(3,3)**2+MomLept(4,3)**2) )

   DeltaPhi = dabs( Get_PHI(MomLept(1:4,1)) - Get_PHI(MomLept(1:4,3))  )
   if( DeltaPhi.gt.Pi ) DeltaPhi=2d0*Pi-DeltaPhi


!                     MomBoost(1)   =+MomTops(1,1)
!                     MomBoost(2:4) =-MomTops(2:4,1)
!                     call boost(MomLept(1:4,1),MomBoost(1:4),m_top)
!                     MomBoost(1)   =+MomTops(1,2)
!                     MomBoost(2:4) =-MomTops(2:4,2)
!                     call boost(MomLept(1:4,3),MomBoost(1:4),m_top)
!                     Psi_LepLep = dacos( VectorProd(MomTops(2:4,1),MomLept(2:4,1))/MomLept(1,1)/dsqrt((MomTops(1,1))**2-m_top**2) )
!                     DeltaPhi   = dacos( VectorProd(MomTops(2:4,2),MomLept(2:4,3))/MomLept(1,3)/dsqrt((MomTops(1,2))**2-m_top**2) )

!                       if( get_eta(MomTops(1:4,1)).lt.0.8d0 ) then
!                           applyPSCut = .true.
!                           RETURN
!                       endif
!                     Psi_LepLep = ( MomLept(4,1)/MomLept(1,1) )
!                     DeltaPhi   = ( MomLept(4,3)/MomLept(1,3) )


   DeltaPhi = dabs( Get_PHI(MomLept(1:4,1)) - Get_PHI(MomLept(1:4,3))  )
   if( DeltaPhi.gt.Pi ) DeltaPhi=2d0*Pi-DeltaPhi
   Psi_LepLep  = dacos( ( MomLept(4,1)*MomLept(4,3) +MomLept(2,1)*MomLept(2,3) + MomLept(3,1)*MomLept(3,3) )/dsqrt(MomLept(2,1)**2+MomLept(3,1)**2+MomLept(4,1)**2)/dsqrt(MomLept(2,3)**2+MomLept(3,3)**2+MomLept(4,3)**2) )




   Ebjets = MomJet(1,1)+MomJet(1,2)

! check cuts
    if( pT_bjet1.lt.pT_bjet_cut .OR. pT_bjet2.lt.pT_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_bjet1).gt.eta_bjet_cut .OR. abs(eta_bjet2).gt.eta_bjet_cut) then
       applyPSCut = .true.
        RETURN
    endif


    if( pT_lepM.lt.pT_lep_cut .OR. pT_lepP.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lepM).gt.eta_lep_cut .OR. abs(eta_lepP).gt.eta_lep_cut) then
       applyPSCut = .true.
        RETURN
    endif

    if( pT_miss.lt.pT_miss_cut ) then
       applyPSCut = .true.
        RETURN
    endif






!  additional cuts to observe spin correlations
!     if( pT_lepM.gt.50d0*GeV .OR. pT_lepP.gt.50d0*GeV ) then
!         applyPSCut = .true.
!         RETURN
!     endif
!
!     if( MInv_lepP_lepM.gt.100d0*GeV .OR. MInv_lepP_lepM.gt.100d0*GeV ) then
!         applyPSCut = .true.
!         RETURN
!     endif

!    if( MInv_Tops.gt.400d0*GeV ) then
!      applyPSCut = .true.
!       RETURN
!    endif

!     if( pT_bjet1.gt.100d0*GeV .OR. pT_bjet2.gt.100d0*GeV ) then
!         applyPSCut = .true.
!         RETURN
!     endif





! !  additional cuts for W+dijet search
!
!     if( pT_miss.lt.25d0*GeV ) then
!        applyPSCut = .true.
!         RETURN
!     endif
!
!     if( pT_bjet1.lt.30d0*GeV .OR. pT_bjet2.lt.30d0*GeV ) then
!         applyPSCut = .true.
!         RETURN
!     endif
!
!     if( abs(eta_bjet1).gt.2.4d0 .OR. abs(eta_bjet2).gt.2.4d0 ) then
!        applyPSCut = .true.
!         RETURN
!     endif
!
!     if( Get_PT( MomJet(1:4,1)+MomJet(1:4,2) ).lt.40d0*GeV ) then
!        applyPSCut = .true.
!         RETURN
!     endif
!
!     if( dabs(Get_ETA(MomJet(1:4,1))  - Get_ETA(MomJet(1:4,2)) ).gt.2.5d0 ) then
!        applyPSCut = .true.
!         RETURN
!     endif
!
!    if( get_pt(MomJet(1:4,1)).gt.get_pt(MomJet(1:4,2)) ) then
!         DeltaPhi = dabs( Get_PHI(MomJet(1:4,1)) - Get_PHI(MomLept(1:4,2)+MomLept(1:4,4))  )
!         if( DeltaPhi.gt.Pi ) DeltaPhi=2d0*Pi-DeltaPhi
!         if( abs(DeltaPhi).lt.0.4d0 ) then
!               applyPSCut = .true.
!               RETURN
!         endif
!    else
!         DeltaPhi = dabs( Get_PHI(MomJet(1:4,2)) - Get_PHI(MomLept(1:4,2)+MomLept(1:4,4))  )
!         if( DeltaPhi.gt.Pi ) DeltaPhi=2d0*Pi-DeltaPhi
!         if( abs(DeltaPhi).lt.0.4d0 ) then
!               applyPSCut = .true.
!               RETURN
!         endif
!    endif
!
!
!
! ! now, pick the lepton
!    if( get_pt(MomLept(1:4,1)).gt.20d0*GeV .AND. get_pt(MomLept(1:4,3)).lt.10d0*GeV ) then! the lepton is MomLept(1:4,1)
!               if( get_ETA(MomLept(1:4,1)).gt.1d0 ) then
!                   applyPSCut = .true.
!                   RETURN
!               endif
!               if( get_MInv(MomLept(1:4,1)+MomLept(1:4,3)).gt.76d0*GeV .and. get_MInv(MomLept(1:4,1)+MomLept(1:4,3)).lt.106d0*GeV ) then
!                   applyPSCut = .true.
!                   RETURN
!               endif
!               if( get_R(MomLept(1:4,1),MomJet(1:4,1)).le.0.52d0 .or. get_R(MomLept(1:4,1),MomJet(1:4,2)).le.0.52d0) then
!                   applyPSCut = .true.
!                   RETURN
!               endif
!
!               MomMiss(1:4) = MomLept(1:4,2)+MomLept(1:4,4)
!               MomMiss(4) = 0d0
!               MomLepTR(1:4) = MomLept(1:4,1)
!               MomLepTR(4) = 0d0
!               mT = dsqrt( 2d0*get_PT(MomLept(1:4,1))*pT_miss*(1d0-Get_CosAlpha(MomMiss(1:4),MomLepTR(1:4))) )
!               if( mT.lt.30d0*GeV ) then
!                   applyPSCut = .true.
!                   RETURN
!               endif
!    elseif( get_pt(MomLept(1:4,3)).gt.20d0*GeV .AND. get_pt(MomLept(1:4,1)).lt.10d0*GeV ) then! the lepton is MomLept(1:4,2)
!               if( get_ETA(MomLept(1:4,3)).gt.1d0 ) then
!                   applyPSCut = .true.
!                   RETURN
!               endif
!               if( get_MInv(MomLept(1:4,1)+MomLept(1:4,3)).gt.76d0*GeV .and. get_MInv(MomLept(1:4,1)+MomLept(1:4,3)).lt.106d0*GeV ) then
!                   applyPSCut = .true.
!                   RETURN
!               endif
!               if( get_R(MomLept(1:4,3),MomJet(1:4,1)).le.0.52d0 .or. get_R(MomLept(1:4,3),MomJet(1:4,2)).le.0.52d0) then
!                   applyPSCut = .true.
!                   RETURN
!               endif
!
!               MomMiss(1:4) = MomLept(1:4,2)+MomLept(1:4,4)
!               MomMiss(4) = 0d0
!               MomLepTR(1:4) = MomLept(1:4,3)
!               MomLepTR(4) = 0d0
!               mT = dsqrt( 2d0*get_PT(MomLept(1:4,3))*pT_miss*(1d0-Get_CosAlpha(MomMiss(1:4),MomLepTR(1:4))) )
!               if( mT.lt.30d0*GeV ) then
!                   applyPSCut = .true.
!                   RETURN
!               endif
!    else
!        applyPSCut = .true.
!        RETURN
!    endif
!   MInv_lepP_lepM = get_MInv( MomJet(1:4,1)+MomJet(1:4,2) )



! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,eta_ATop)
    NBin(3) = WhichBin(3,pT_lepP)
    NBin(4) = WhichBin(4,eta_lepP)
    NBin(5) = WhichBin(5,HT)
    NBin(6) = WhichBin(6,MInv_lepP_bjet)
    NBin(7) = WhichBin(7,get_MInv(MomLept(1:4,3)+MomDK(1:4,4))) ! careful! MomDK(:,4) is not IR save
!    NBin(7) = WhichBin(7,Psi_LepLep)
    NBin(8) = WhichBin(8,DeltaPhi)
    NBin(9) = WhichBin(9,Ebjets)

    NBin(10) = WhichBin(10,ET_miss)
    NBin(11) = WhichBin(11,E_lept)
    NBin(12) = WhichBin(12,pT_bjet1)
    NBin(13) = WhichBin(13,pT_bjet2)
    NBin(14) = WhichBin(14,E_bjet1)
    NBin(15) = WhichBin(15,E_bjet2)
    NBin(16) = WhichBin(16,MInv_lepP_lepM)
    NBin(17) = WhichBin(17,ET_lept)
    NBin(18) = WhichBin(18,mT)

    NBin(19) = WhichBin(19,eta_lepM)
    NBin(20) = WhichBin(20,eta_lepP)



!-------------------------------------------------------
elseif( ObsSet.eq.4 ) then! set of observables for ttb production with hadr. Atop, lept. top decay

!   request at least two b-jets
    if( .not.(any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif
    call pT_order(2,MomJet(1:4,1:2))
    call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))

    pT_bjet1 = get_PT(MomJet(1:4,1))
    pT_bjet2 = get_PT(MomJet(1:4,2))
    eta_bjet1 = get_ETA(MomJet(1:4,1))
    eta_bjet2 = get_ETA(MomJet(1:4,2))

!   check if b-jets pass cuts
    if(  pT_bjet1.lt.pT_bjet_cut .or. abs(eta_bjet1).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( pT_bjet2.lt.pT_bjet_cut .or. abs(eta_bjet2).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    HT = pT_bjet1 + pT_bjet2

!   determine observable light jets
    NObsJet = 2! these are the b-jets
    do k=3,NJet
        if( get_PT(MomJet(1:4,k)).gt.pT_jet_cut .and. abs(get_ETA(MomJet(1:4,k))).lt.eta_jet_cut ) then! count jets outside beam pipe
            NObsJet = NObsJet +1
            if( k.ne.NObsJet ) MomJet(1:4,NObsJet) = MomJet(1:4,k)
            HT = HT + get_PT(MomJet(1:4,k))
        endif
    enddo

    NObsJet_Tree = 4! request two b-jets and at least two light jets
    if( NObsJet.lt.NObsJet_Tree ) then
!     if( NObsJet.ne.NObsJet_Tree ) then!!!    CAREFUL: this cuts out the additional hard jet: only for combination with ttbjet
        applyPSCut = .true.
        RETURN
    endif


    pT_lepP = get_PT(MomLept(1:4,3))
    eta_lepP = get_ETA(MomLept(1:4,3))

    ET_miss  = get_ET(MomLept(1:4,4))
    HT = HT + pT_lepP + ET_miss

    if( pT_bjet1.gt.pT_bjet2 ) then
        m_lb = get_Minv(MomLept(1:4,3)+MomJet(1:4,1))
    else
        m_lb = get_Minv(MomLept(1:4,3)+MomJet(1:4,2))
    endif

    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))

    eta_ATop = get_ETA(MomTops(1:4,1))
    eta_Top  = get_ETA(MomTops(1:4,2))


! check cuts
    if( HT.lt.HT_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_lepP.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lepP).gt.eta_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( ET_miss.lt.pT_miss_cut ) then
        applyPSCut = .true.
        RETURN
    endif


!   construct ttb momentum
    MomTops(1:4,1) = MomJet(1:4,1)+MomJet(1:4,2)+MomLept(1:4,3)+MomLept(1:4,4) + MomJet(1:4,3)+MomJet(1:4,4)
    if( dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,4))-M_W ).lt.20d0*GeV ) then!   require a 20GeV window around M_W
        pT_Top = get_pT(MomTops(1:4,1))
    else
        pT_Top   = -1d0
    endif



! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,eta_ATop)
    NBin(3) = WhichBin(3,pT_Top)
    NBin(4) = WhichBin(4,eta_Top)
    NBin(5) = WhichBin(5,eta_ATop)
    NBin(6) = WhichBin(6,eta_lepP)
    NBin(7) = WhichBin(7,pT_lepP)
    NBin(8) = WhichBin(8,eta_lepP)
    NBin(9) = WhichBin(9,ET_miss)
    NBin(10)= WhichBin(10,HT)
    NBin(11)= WhichBin(11,m_lb)
    NBin(12)= WhichBin(12,pT_Top)



!-------------------------------------------------------
elseif( ObsSet.eq.5 ) then! set of observables for ttb production with hadr. Atop, lept. top decay


!   request at least two b-jets
    if( .not.(any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif
    call pT_order(2,MomJet(1:4,1:2))
    call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))

    pT_bjet1 = get_PT(MomJet(1:4,1))
    pT_bjet2 = get_PT(MomJet(1:4,2))
    eta_bjet1 = get_ETA(MomJet(1:4,1))
    eta_bjet2 = get_ETA(MomJet(1:4,2))

!   check if b-jets pass cuts
    if(  pT_bjet1.lt.pT_bjet_cut .or. abs(eta_bjet1).gt.eta_bjet_cut .or. get_R(MomLept(1:4,3),MomJet(1:4,1)).lt.Rsep_LepJet ) then
        applyPSCut = .true.
        RETURN
    endif
    if( pT_bjet2.lt.pT_bjet_cut .or. abs(eta_bjet2).gt.eta_bjet_cut .or. get_R(MomLept(1:4,3),MomJet(1:4,2)).lt.Rsep_LepJet ) then
        applyPSCut = .true.
        RETURN
    endif
    HT = pT_bjet1 + pT_bjet2

!   determine observable light jets
    NObsJet = 2! these are the b-jets
    do k=3,NJet
        if( get_PT(MomJet(1:4,k)).gt.pT_jet_cut .and. abs(get_ETA(MomJet(1:4,k))).lt.eta_jet_cut .and. get_R(MomLept(1:4,3),MomJet(1:4,k)).gt.Rsep_LepJet ) then! count jets outside beam pipe
            NObsJet = NObsJet +1
            if( k.ne.NObsJet ) MomJet(1:4,NObsJet) = MomJet(1:4,k)
            HT = HT + get_PT(MomJet(1:4,k))
        endif
    enddo

    NObsJet_Tree = 4! request two b-jets and at least two light jets
    if( NObsJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif



    pT_lepP = get_PT(MomLept(1:4,3))
    eta_lepP = get_ETA(MomLept(1:4,3))

    ET_miss  = get_ET(MomLept(1:4,4))
    HT = HT + pT_lepP + ET_miss

    if( pT_bjet1.gt.pT_bjet2 ) then
        m_lb = get_Minv(MomLept(1:4,3)+MomJet(1:4,1))
    else
        m_lb = get_Minv(MomLept(1:4,3)+MomJet(1:4,2))
    endif

    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))

    eta_ATop = get_ETA(MomTops(1:4,1))
    eta_Top  = get_ETA(MomTops(1:4,2))


! check cuts
    if( pT_lepP.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lepP).gt.eta_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( ET_miss.lt.pT_miss_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    mT_lp = get_MT(MomLept(1:4,3),MomLept(1:4,4))! this is the transverse W mass

    if( ET_miss+mT_lp.lt.60*GeV  ) then
        applyPSCut = .true.
        RETURN
    endif


!   construct hadr. W momentum
    MomTops(1:4,1) = MomJet(1:4,1)+MomJet(1:4,2)+MomLept(1:4,3)+MomLept(1:4,4) + MomJet(1:4,3)+MomJet(1:4,4)
    if( dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,4))-M_W ).lt.20d0*GeV ) then!   require a 20GeV window around M_W
        pT_Top = get_pT(MomTops(1:4,1))
    else
        pT_Top   = -1d0
    endif



! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,eta_ATop)
    NBin(3) = WhichBin(3,pT_Top)
    NBin(4) = WhichBin(4,eta_Top)
    NBin(5) = WhichBin(5,eta_ATop)
    NBin(6) = WhichBin(6,eta_Top)
    NBin(7) = WhichBin(7,pT_lepP)
    NBin(8) = WhichBin(8,eta_lepP)
    NBin(9) = WhichBin(9,ET_miss)
    NBin(10)= WhichBin(10,HT)
    NBin(11)= WhichBin(11,m_lb)





!-------------------------------------------------------
elseif( ObsSet.eq.6 ) then! set of observables for ttb production with hadr. Atop, lept. top decay


!   request at least two b-jets
    if( .not.(any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif
    call pT_order(2,MomJet(1:4,1:2))
    call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))

    pT_bjet1 = get_PT(MomJet(1:4,1))
    pT_bjet2 = get_PT(MomJet(1:4,2))
    eta_bjet1 = get_ETA(MomJet(1:4,1))
    eta_bjet2 = get_ETA(MomJet(1:4,2))

!   check if b-jets pass cuts
    if(  pT_bjet1.lt.pT_bjet_cut .or. abs(eta_bjet1).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( pT_bjet2.lt.pT_bjet_cut .or. abs(eta_bjet2).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    HT = pT_bjet1 + pT_bjet2

!   determine observable light jets
    NObsJet = 2! these are the b-jets
    do k=3,NJet
        if( get_PT(MomJet(1:4,k)).gt.pT_jet_cut .and. abs(get_ETA(MomJet(1:4,k))).lt.eta_jet_cut ) then! count jets outside beam pipe
            NObsJet = NObsJet +1
            if( k.ne.NObsJet ) MomJet(1:4,NObsJet) = MomJet(1:4,k)
            HT = HT + get_PT(MomJet(1:4,k))
        endif
    enddo

    NObsJet_Tree = 4! request two b-jets and at least two light jets
    if( NObsJet.lt.NObsJet_Tree ) then
!     if( NObsJet.ne.NObsJet_Tree ) then!!!    CAREFUL: this cuts out the additional hard jet: only for combination with ttbjet
        applyPSCut = .true.
        RETURN
    endif



    pT_lepP = get_PT(MomLept(1:4,3))
    eta_lepP = get_ETA(MomLept(1:4,3))

    ET_miss  = get_ET(MomLept(1:4,4))
    HT = HT + pT_lepP + ET_miss

    if( pT_bjet1.gt.pT_bjet2 ) then
        m_lb = get_Minv(MomLept(1:4,3)+MomJet(1:4,1))
    else
        m_lb = get_Minv(MomLept(1:4,3)+MomJet(1:4,2))
    endif

    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))

    eta_ATop = get_ETA(MomTops(1:4,1))
    eta_Top  = get_ETA(MomTops(1:4,2))


! check cuts
    if( pT_lepP.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lepP).gt.eta_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( ET_miss.lt.pT_miss_cut ) then
        applyPSCut = .true.
        RETURN
    endif

!     if( HT.lt.HT_cut ) then
!         applyPSCut = .true.
!         RETURN
!     endif



!   construct hadr. ttb pT
    MomTops(1:4,1) = MomJet(1:4,1)+MomJet(1:4,2)+MomLept(1:4,3)+MomLept(1:4,4) + MomJet(1:4,3)+MomJet(1:4,4)
    if( dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,4))-M_W ).lt.20d0*GeV ) then!   require a 20GeV window around M_W
        pT_Top = get_pT(MomTops(1:4,1))
    else
        pT_Top   = -1d0
    endif



! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,eta_ATop)
    NBin(3) = WhichBin(3,pT_Top)
    NBin(4) = WhichBin(4,eta_Top)
    NBin(5) = WhichBin(5,eta_ATop)
    NBin(6) = WhichBin(6,eta_Top)
    NBin(7) = WhichBin(7,pT_lepP)
    NBin(8) = WhichBin(8,eta_lepP)
    NBin(9) = WhichBin(9,ET_miss)
    NBin(10)= WhichBin(10,HT)
    NBin(11)= WhichBin(11,m_lb)
    NBin(12)= WhichBin(12,pT_Top)





!-------------------------------------------------------
elseif( ObsSet.eq.7 ) THEN! set of observables for ttb production with lept. top and J/Psi fragmentation, hadr. Atop decay at LHC

    call pT_order(NumHadr,MomJet(1:4,1:NumHadr))

    NObsJet_Tree = 4
    NObsJet = 0
    HT = 0d0
    do n=1,NJet
        pT_jet1 = get_pT(MomJet(1:4,n))
        if( pT_jet1.gt.pT_jet_cut ) NObsJet = NObsJet +1
        HT = HT + pT_jet1
    enddo

    if( NObsJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif

    if( HT.lt.HT_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    pT_lepP = get_PT(MomLept(1:4,3))
    if( pT_lepP.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif
! print *,"removed cuts",  NBin(1:3)


    MInv_LB = get_MInv( xJPsiFrag*MomDK(1:4,4)+MomLept(1:4,3) )
    pT_lepP = get_PT(MomLept(1:4,3))
    pT_bjet1 = get_pT(MomJPsi(1:4)) ! pT(J/Psi)


! binning
    NBin(1) = WhichBin(1,MInv_LB)
    NBin(2) = WhichBin(2,pT_bjet1)
    NBin(3) = WhichBin(3,pT_lepP)




!-------------------------------------------------------
elseif( ObsSet.eq.8 ) then! set of observables for ttb spin correlations at LHC (di-lept. decay)

! request at least two b-jets
    if( .not.(NJet.ge.2 .and. any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif
   call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))


! evaluate kinematic variables
    pT_bjet1 = get_PT(MomJet(1:4,1))
    pT_bjet2 = get_PT(MomJet(1:4,2))

    pT_miss = get_PT(MomLept(1:4,2)+MomLept(1:4,4))
    pT_miss = get_PT(MomLept(1:4,2)+MomLept(1:4,4))

    pT_lepM = get_PT(MomLept(1:4,1))
    pT_lepP = get_PT(MomLept(1:4,3))

    eta_lepM = get_ETA(MomLept(1:4,1))
    eta_lepP = get_ETA(MomLept(1:4,3))

    DeltaPhi = dabs( Get_PHI(MomLept(1:4,1)) - Get_PHI(MomLept(1:4,3))  )
    if( DeltaPhi.gt.Pi ) DeltaPhi=2d0*Pi-DeltaPhi

    MInv_Tops = get_MInv( MomTops(1:4,1)+MomTops(1:4,2) )

! check cuts
!     if( pT_lepM.gt.50d0*GeV .OR. pT_lepP.gt.50d0*GeV ) then
!         applyPSCut = .true.
!         RETURN
!     endif

    if( pT_bjet1.lt.pT_bjet_cut .OR. pT_bjet2.lt.pT_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_lepM.lt.pT_lep_cut .OR. pT_lepP.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_miss.lt.pT_miss_cut ) then
       applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lepM).gt.eta_lep_cut .OR. abs(eta_lepP).gt.eta_lep_cut) then
       applyPSCut = .true.
        RETURN
    endif

!    if( MInv_Tops.gt.400d0*GeV ) then
!      applyPSCut = .true.
!       RETURN
!    endif


   MomHadr(1:4,1) = MomJet(1:4,1)
   MomHadr(1:4,2:3) = MomDK(1:4,2:3)

   MomHadr(1:4,4) = MomJet(1:4,2)
   MomHadr(1:4,5:6) = MomDK(1:4,5:6)
   if( Collider.eq.1 ) then
     r_sc = calc_rgg(MomExt(1:4,1:4),MomHadr(1:4,1:6))
   else
     r_sc = calc_rqq(MomExt(1:4,1:4),MomHadr(1:4,1:6))
   endif

! binning
    NBin(1) = WhichBin(1,r_sc)
    NBin(2) = WhichBin(2,pT_lepM)
    NBin(3) = WhichBin(3,pT_lepP)
    NBin(4) = WhichBin(4,DeltaPhi)

elseif (ObsSet.EQ.60) then ! set of observables for ttb production without decays, for Zprime background

    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))

    y_ATop = get_ETA(MomTops(1:4,1))
    y_Top  = get_ETA(MomTops(1:4,2))

    M_TTbar = get_MInv(MomTops(1:4,1)+MomTops(1:4,2))

    if (M_TTbar.lt.Mttb_cut) then
        applyPSCut = .true.
        RETURN
    endif


    Dphi_TTbar = dabs( Get_PHI(MomTops(1:4,1)) - Get_PHI(MomTops(1:4,2))  )
    if( Dphi_TTbar.gt.Pi ) Dphi_TTbar=2d0*Pi-Dphi_TTbar

!   scattering angle of the top quark in lab frame
    CosTheta_scatter = get_CosTheta( MomTops(1:4,2) )

!   see chinese paper, approximates scattering angle in ttb rest frame
    MomAux1(1:4) = MomTops(1:4,2)
    MomAux2(1:4) = MomExt(1:4,1)
    call boost(MomAux1(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)! this seems wrong!!!
    call boost(MomAux2(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)
    CosTheta_star = get_CosAlpha(MomAux1(1:4),MomAux2(1:4))



! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,y_ATop)
    NBin(4) = WhichBin(4,y_Top)
    NBin(5) = WhichBin(5,M_TTbar)
    NBin(6) = WhichBin(6,Dphi_TTbar)
    NBin(7) = WhichBin(7,CosTheta_scatter)
    NBin(8) = WhichBin(8,CosTheta_star)

    if( NPlus1PS ) NBin(9) = WhichBin(9,get_MInv(MomTops(1:4,1)+MomTops(1:4,2)+MomExt(1:4,3)))
    if( .not. NPlus1PS ) NBin(9) = WhichBin(9,get_MInv(MomTops(1:4,1)+MomTops(1:4,2)))

!-------------------------------------------------------
elseif( ObsSet.eq.61 ) then! set of observables for ttb production with di-lept. decays, for Zprime background

   !-- define ordering for leptons, without modifying existing lepton array used below for angles
   MomLeptOrd(1:4,1) = MomLept(1:4,1)
   MomLeptOrd(1:4,2) = MomLept(1:4,3)

   call pT_order(2,MomLeptOrd)! pT-order for leptons

   !-- order all jets, irrespective whether they are b-jets or not
   MomJetOrd(1:4,1:8) = MomJet(1:4,1:8)
   call pT_order(NumHadr,MomJetOrd(1:4,1:NumHadr))


!   M_LL = get_MInv(MomLept(1:4,1)+MomLept(1:4,3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--F Phase space cuts

!   check that there are two b jets    
    if( .not.(any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif

    pT_bjet2 = get_PT(MomJet(1:4,2))
    pT_lep2  = get_PT(MomLeptOrd(1:4,2))

    eta_bjet1 = get_eta(MomJet(1:4,1))
    eta_bjet2 = get_eta(MomJet(1:4,2))

    eta_lep1 = get_eta(MomLept(1:4,1))
    eta_lep2 = get_eta(MomLeptOrd(1:4,2))

    if (pT_lep2.lt.pt_lep_cut) then
       applypscut = .true.
       RETURN
    endif

    if (dabs(eta_lep1).gt.eta_lep_cut .or. dabs(eta_lep2).gt.eta_lep_cut) then
       applypscut = .true.
       RETURN
    endif

    if (pT_bjet2.lt.pt_bjet_cut) then
       applypscut = .true.
       RETURN
    endif

    if (dabs(eta_bjet1).gt.eta_bjet_cut .or. dabs(eta_bjet2).gt.eta_bjet_cut) then
       applypscut = .true.
       RETURN
    endif

!-- If we pass the cuts, compute other observables
    pT_bjet1 = get_PT(MomJet(1:4,1))
    pT_lep1 = get_PT(MomLeptOrd(1:4,1))

    M_TTbar = get_MInv(MomTops(1:4,1)+MomTops(1:4,2))

    Dphi_LL = dabs( Get_PHI(MomLept(1:4,1)) - Get_PHI(MomLept(1:4,3))  )
    if( Dphi_LL.gt.Pi ) Dphi_LL=2d0*Pi-Dphi_LL

    pT_miss = get_pT(MomLept(1:4,2) + MomLept(1:4,4))
    Minv_Lept = get_MInv(MomLept(1:4,1)+MomLept(1:4,3))
    ET_miss = dsqrt(pT_miss**2 + Minv_Lept**2)

    pT_jet1 = get_PT(MomJetOrd(1:4,1))
    pT_jet2 = get_PT(MomJetOrd(1:4,2))

    HT = pT_lep1 + pT_lep2 + pT_jet1 + pT_jet2

    M_eff = ET_miss + HT

!   construct Baumgart-Tweedie angles // for non-zero transverse momentum of TTBAR, use Collins-Soper construction

    if (dabs(M_TTbar-M_Zpr).le.(100d0*GeV)) then ! compute angles only close to the resonance

!   step 1: boost top momenta into a TTBAR CMS frame ( pT+pTbar=0 )
    MomAux1(1:4) = MomTops(1:4,1)+MomTops(1:4,2)
    MomAux1(2:4) =-MomAux1(2:4)
    MassAux = dsqrt( (MomAux1(1:4).dot.MomAux1(1:4)) +1d-12 )
    MomTopsCMS(1:4,1:2) = MomTops(1:4,1:2)
    call boost(MomTopsCMS(1:4,1),MomAux1(1:4),MassAux)
    call boost(MomTopsCMS(1:4,2),MomAux1(1:4),MassAux)

!   step 2: construct coordinate system in a TTBAR CMS frame

    MomBeam1(1:4) = (/1d0,0d0,0d0,1d0/)
    call boost(MomBeam1(1:4),MomAux1(1:4),MassAux)

    MomBeam2(1:4) = (/1d0,0d0,0d0,-1d0/)
    call boost(MomBeam2(1:4),MomAux1(1:4),MassAux)

    MomBeam(2:4) = MomBeam1(2:4)-MomBeam2(2:4)
    MomBeam(2:4) = MomBeam(2:4)/dsqrt(MomBeam(2)**2+MomBeam(3)**2+MomBeam(4)**2)     ! Collins-Soper direction

    ny(2:4) = MomTopsCMS(2:4,2).cross.MomBeam(2:4)
    ny(2:4) = ny(2:4)/dsqrt(ny(2)**2+ny(3)**2+ny(4)**2)
    nz(2:4) = MomTopsCMS(2:4,2)/dsqrt(MomTopsCMS(2,2)**2+MomTopsCMS(3,2)**2+MomTopsCMS(4,2)**2)
    nx(2:4) = ny(2:4).cross.nz(2:4)

!   step 3: boost lepton momenta first into CMS frame, then into frame with ptop=0 and call it MomLeptTRF
    MomLeptTRF(1:4) = MomLept(1:4,3)
    call boost(MomLeptTRF(1:4),MomAux1(1:4),MassAux)
    MomAux2(1:4) = MomTopsCMS(1:4,2)
    MomAux2(2:4) =-MomAux2(2:4)
    call boost(MomLeptTRF(1:4),MomAux2(1:4),m_top)

!   step 4: project MomLeptTRF onto the nx-ny plane and calculate its angle with nx

    cosPhi = VectorProd(MomLeptTRF(2:4),nx(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2  & 
    + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )

    sinPhi = VectorProd(MomLeptTRF(2:4),ny(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2  & 
    + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )


!    Phi = dacos(cosPhi)

!------- repeat the same for the anti-top

!   step 3: boost lepton momenta into frame with patop=0 and call it MomLeptTRF
    MomLeptTRF(1:4) = MomLept(1:4,1)
    call boost(MomLeptTRF(1:4),MomAux1(1:4),MassAux)

    MomAux2(1:4) = MomTopsCMS(1:4,1)
    MomAux2(2:4) =-MomAux2(2:4)
    call boost(MomLeptTRF(1:4),MomAux2(1:4),m_top)

!   step 4: project MomLeptTRF onto the nx-ny plane and calculate its angle with n
    cosPhibar = VectorProd(MomLeptTRF(2:4),nx(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2 & 
      + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )

    sinPhibar = VectorProd(MomLeptTRF(2:4),ny(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2 & 
      + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )

!    Phibar = dacos(cosPhibar)

    deltaSin = sinPhi * cosPhibar - cosPhi * sinPhibar

    if (deltaSin .gt. 0d0) then
       dPhiMinus = dacos( cosPhi*cosPhibar + sinPhi*sinPhibar)
    else
       dPhiMinus = -dacos( cosPhi*cosPhibar + sinPhi*sinPhibar)
    endif

    dPhiPlus = dabs(dacos(cosPhi*cosPhibar - sinPhi*sinPhibar))

!    dPhiMinus = Phi-Phibar          ! -pi..pi
!    dPhiPlus  = dabs(Phi+Phibar)    ! 0...2pi

!    print *, phi/DblPi,phibar/DblPi
!    print *, dPhiMinus/DblPi,dPhiPlus/DblPi
!pause

    else
       dPhiPlus=20d0
       dPhiMinus=20d0

    endif

! binning
    NBin(1) = WhichBin(1,pT_lep1)
    NBin(2) = WhichBin(2,pT_bjet1)
    NBin(3) = WhichBin(3,pT_lep2)
    NBin(4) = WhichBin(4,pT_bjet2)
    NBin(5) = WhichBin(5,M_TTbar)
    NBin(6) = WhichBin(6,M_eff)
    NBin(7) = WhichBin(7,eta_lep1)
    NBin(8) = WhichBin(8,eta_bjet1)
    NBin(9) = WhichBin(9,eta_lep2)
    NBin(10) = WhichBin(10,eta_bjet2)
    NBin(11) = WhichBin(11,dPhi_LL)
    NBin(12) = WhichBin(12,dPhiMinus)
    NBin(13) = WhichBin(13,dPhiPlus)

endif

return
END SUBROUTINE










SUBROUTINE Kinematics_TTbarETmiss(NPlus1PS,Mom,MomOrder,applyPSCut,NBin)
use ModMisc
use ModParameters
implicit none
integer :: NumHadr,MomOrder(1:13)
real(8) :: Mom(1:4,1:15),zeros(1:13)
real(8) :: MomJet(1:4,1:7)
real(8) :: MomHadr(1:4,0:8)
real(8) :: MomBoost(1:4),MomMiss(1:4),MomObs(1:4),MomAux(1:4)
logical :: applyPSCut,NPlus1PS
integer :: NBin(:),PartList(1:7),JetList(1:7),NJet,NObsJet,k,NObsJet_Tree
real(8) :: pT_lepM,pT_lepP,pT_miss,pT_ATop,pT_Top
real(8) :: eta_ATop,eta_Top,eta_lepM,eta_lepP,m_lb
real(8) :: pT_jet(1:7),eta_jet(1:7),m_X0,phi_ll,m_ll,MTW,MTeff
integer :: Htbar,Ht,X0bar,X0,tbar,t,realp,bbar,lepM,nubar,b,lepP,nu,qdn,qbup,qbdn,qup,Lep,Neu,ib_lep,ib_had
real(8) :: Stopness,rpT,H_T



! momentum ordering
  Htbar   = MomOrder(1)
  Ht      = MomOrder(2)
  X0bar   = MomOrder(3)
  X0      = MomOrder(4)
  tbar    = MomOrder(5)
  t       = MomOrder(6)
  bbar    = MomOrder(7)
  lepM    = MomOrder(8)
  nubar   = MomOrder(9)
  b       = MomOrder(10)
  lepP    = MomOrder(11)
  nu      = MomOrder(12)
  realp   = MomOrder(13)

  qdn    = lepM
  qbup   = nubar
  qbdn   = lepP
  qup    = nu

IF( ObsSet.eq.31 .or. ObsSet.eq.32 .or. ObsSet.eq.33 .or. ObsSet.eq.34 .or. ObsSet.eq.35 ) then
   IF(XTOPDECAYS.EQ.1) m_X0 = m_BH
   IF(XTOPDECAYS.EQ.2) m_X0 = m_A0
elseif( ObsSet.eq.41 .or. ObsSet.eq.42 .or. ObsSet.eq.43 .or. ObsSet.eq.44 .or. ObsSet.eq.45 ) then
   m_X0 = m_Chi
endif

!DEC$ IF(_CheckMomenta .EQ.1)
IF(XTOPDECAYS.NE.0) THEN
   zeros(1:4) = Mom(1:4,1)+Mom(1:4,2) - (Mom(1:4,X0bar)+Mom(1:4,X0)+Mom(1:4,bbar)+Mom(1:4,lepM)+Mom(1:4,nubar)+Mom(1:4,b)+Mom(1:4,lepP)+Mom(1:4,nu))
   if( NPlus1PS ) zeros(1:4) = zeros(1:4) - Mom(1:4,realp)
   if( any(abs(zeros(1:4)/m_HTop).gt.1d-8) ) then
      print *, "ERROR: energy-momentum violation 1 in SUBROUTINE Kinematics_TTbarETmiss: ",zeros(1:4),NPlus1PS
      print *, "momentum dump:"
      print *, Mom(:,:)
   endif

   if( Correction.ne.5 .or. .not.NPlus1PS ) then
          zeros(1:4) = Mom(1:4,Htbar)-Mom(1:4,X0bar)-Mom(1:4,tbar)
          zeros(5:8) = Mom(1:4,Ht)-Mom(1:4,X0)-Mom(1:4,t)
          if( any(abs(zeros(1:8)/m_HTop).gt.1d-8) ) then
              print *, "ERROR: energy-momentum violation 2 in SUBROUTINE Kinematics_TTbarETmiss: ",zeros(1:8),NPlus1PS
              print *, "momentum dump:"
              print *, Mom(:,:)
          endif

          zeros(1:4) = Mom(1:4,tbar)-Mom(1:4,bbar)-Mom(1:4,LepM)-Mom(1:4,nubar)
          zeros(5:8) = Mom(1:4,t)-Mom(1:4,b)-Mom(1:4,lepP)-Mom(1:4,nu)
          if( any(abs(zeros(1:8)/m_HTop).gt.1d-8) ) then
              print *, "ERROR: energy-momentum violation 3 in SUBROUTINE Kinematics_TTbarETmiss: ",zeros(1:8),NPlus1PS
              print *, "momentum dump:"
              print *, Mom(:,:)
          endif
    else
          zeros(1:4) = Mom(1:4,Htbar)-Mom(1:4,X0bar) + Mom(1:4,Ht)-Mom(1:4,X0) & 
                     - Mom(1:4,bbar)-Mom(1:4,LepM)-Mom(1:4,nubar) -Mom(1:4,b)-Mom(1:4,lepP)-Mom(1:4,nu) - Mom(1:4,realp)
          if( any(abs(zeros(1:4)/m_HTop).gt.1d-8) ) then
              print *, "ERROR: energy-momentum violation 3 in SUBROUTINE Kinematics_TTbarETmiss: ",zeros(1:4),NPlus1PS
              print *, "momentum dump:"
              print *, Mom(:,:)
          endif
    endif



   zeros(:) = 0d0
IF( ObsSet.eq.31 .or. ObsSet.eq.32 .or. ObsSet.eq.33 .or. ObsSet.eq.34 .or. ObsSet.eq.35 ) then
   zeros(1) = (Mom(1:4,HTbar).dot.Mom(1:4,HTbar))  - m_HTop**2
   zeros(2) = (Mom(1:4,HT).dot.Mom(1:4,HT)) - m_HTop**2
elseif( ObsSet.eq.41 .or. ObsSet.eq.42 .or. ObsSet.eq.43 .or. ObsSet.eq.44 .or. ObsSet.eq.45 ) then
   zeros(1) = (Mom(1:4,HTbar).dot.Mom(1:4,HTbar))  - m_STop**2
   zeros(2) = (Mom(1:4,HT).dot.Mom(1:4,HT)) - m_STop**2
endif
   zeros(3) = (Mom(1:4,X0bar).dot.Mom(1:4,X0bar)) - m_X0**2
   zeros(4) = (Mom(1:4,X0).dot.Mom(1:4,X0)) - m_X0**2
   zeros(5) = (Mom(1:4,tbar).dot.Mom(1:4,tbar)) - m_SMTop**2
   zeros(6) = (Mom(1:4,t).dot.Mom(1:4,t)) - m_SMTop**2
   zeros(7) = (Mom(1:4,bbar).dot.Mom(1:4,bbar))
   zeros(8) = (Mom(1:4,b).dot.Mom(1:4,b))
   zeros(9) = (Mom(1:4,lepm).dot.Mom(1:4,lepm))
   zeros(10) = (Mom(1:4,lepp).dot.Mom(1:4,lepp))
   zeros(11) = (Mom(1:4,nubar).dot.Mom(1:4,nubar))
   zeros(12) = (Mom(1:4,nu).dot.Mom(1:4,nu))
   if(NPlus1PS) zeros(13) = (Mom(1:4,realp).dot.Mom(1:4,realp))

   if( any(abs(zeros(1:13)/1d0).gt.1d-5) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TTbarETmiss: ",zeros(1:13),NPlus1PS
      print *, "momentum dump:"
      print *, Mom(:,:)
   endif
ENDIF
!DEC$ ENDIF



applyPSCut = .false.
NBin(1:NumHistograms) = 0


MomHadr(1:4,1:7) = 0d0
PartList(1:7)=(/1,2,3,4,5,6,7/)

MomMiss(1:4) = 0d0
MomObs(1:4)  = 0d0

! separating momenta into hadron momenta and lepton momenta to which various procedures can be applied
!-------------------------------------------------------

IF( TOPDECAYS.EQ.0 ) THEN  ! no decays
    if(.not.NPlus1PS) then
        NumHadr = 0
    else
        NumHadr = 1
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
    endif

ELSEIF( TOPDECAYS.EQ.1 ) THEN  ! di-leptonic decay
    MomHadr(1:4,1) = Mom(1:4,bbar)
    MomHadr(1:4,2) = Mom(1:4,b)     ! Bot
    if(.not.NPlus1PS) then
        NumHadr = 2
    else
        NumHadr = 3
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
    endif

ELSEIF( TOPDECAYS.EQ.3 ) THEN  ! lept. Atop, hadr. top decay
   MomHadr(1:4,1) = Mom(1:4,bbar)
   MomHadr(1:4,2) = Mom(1:4,b)
   MomHadr(1:4,3) = Mom(1:4,qbdn)
   MomHadr(1:4,4) = Mom(1:4,qup)
   Lep = LepM
   Neu = nubar
   if(.not.NPlus1PS) then
        NumHadr = 4
   else
        NumHadr = 5
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)   ! q/qbar/glu
   endif

ELSEIF( TOPDECAYS.EQ.4 ) THEN  ! hadr. Atop, lept. top decay
   MomHadr(1:4,1) = Mom(1:4,bbar)
   MomHadr(1:4,2) = Mom(1:4,b)
   MomHadr(1:4,3) = Mom(1:4,qbup)
   MomHadr(1:4,4) = Mom(1:4,qdn)
   Lep = LepP
   Neu = nu
   if(.not.NPlus1PS) then
        NumHadr = 4
   else
        NumHadr = 5
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
   endif

ELSE
  call Error("this decay is not yet implemented")
ENDIF



!---------------------- kT jet algorithm ---------------------------------
! important: b-quarks need to be at the first two positions of MomHadr(:,:)
! recombination of MomHadr(1:4,i) and MomHadr(1:4,j) results in MomJet(1:4,i) with i<j.
! but take care: some MomJet(1:4,i) can be zero, this is why we apply pt_order

    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr)) ! hard protojets in beam pipe are counted as jets
    call pT_order(2,MomJet(1:4,1:2))
    call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))
! call SwitchEnergyComponent(MomHadr(1:4,1:NumHadr))
! call fastjetppgenkt(MomHadr(1:4,1:NumHadr),NumHadr,Rsep_jet,-1d0,MomJet(1:4,1:NumHadr),NJet)! 4th argument:  +1d0=kT,  -1d0=anti-kT,   0d0=CA
! call SwitchEnergyComponentBack(MomJet_CHECK(1:4,1:NumHadr))

!------------------------ cuts and binning --------------------------------
if( ObsSet.eq.31 .or. ObsSet.eq.41) then! no decays


    pT_jet(1) = get_PT(Mom(1:4,HTbar))


! binning
    NBin(1)= WhichBin(1,pT_jet(1))


!-------------------------------------------------------
elseif( ObsSet.eq.32 .or. ObsSet.eq.34 .or. ObsSet.eq.42 .or. ObsSet.eq.44 )  then! set of observables for TTbar -> ttbar + ETmiss  in di-lept. top decays

! request at least two b-jets
    NObsJet_Tree = 2
    if( .not.(NJet.ge.NObsJet_Tree .and. any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif

    pT_jet(1) = get_PT(MomJet(1:4,1))
    pT_jet(2) = get_PT(MomJet(1:4,2))
    eta_jet(1)= get_ETA(MomJet(1:4,1))
    eta_jet(2)= get_ETA(MomJet(1:4,2))

    pT_lepM = get_PT(Mom(1:4,lepM))
    pT_lepP = get_PT(Mom(1:4,lepP))
    eta_lepM = get_ETA(Mom(1:4,lepM))
    eta_lepP = get_ETA(Mom(1:4,lepP))

    MomMiss(1:4) = Mom(1:4,nu)+Mom(1:4,nubar)+Mom(1:4,X0)+Mom(1:4,X0bar)
    pT_miss = get_ET( MomMiss(1:4) )! note that this is ET and not pT


    phi_ll = dabs( Get_PHI(Mom(1:4,lepM)) - Get_PHI(Mom(1:4,lepP)) )
    if( phi_ll.gt.Pi ) phi_ll=2d0*Pi-phi_ll

    m_ll = get_MInv(Mom(1:4,lepM)+Mom(1:4,lepP))


    MTW = dsqrt(  2d0*pT_lepP*get_ET( MomMiss(1:4) )*(1d0-Get_CosPhi(Mom(1:4,lepP),MomMiss(1:4))) )! let's define MTW with ET instead of pT in accordance with ATLAS

    MTeff = pt_jet(1) + pt_jet(2) + pT_lepP + pT_lepM + pT_miss
    if( NJet.eq.3 ) MTeff = MTeff + pt_jet(3)

    H_T = pt_jet(1) + pt_jet(2) + MTW ! definition in accordance with ATLAS


! check cuts
    if( pT_jet(1).lt.pT_bjet_cut .or. pT_jet(2).lt.pT_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( abs(eta_jet(1)).gt.eta_bjet_cut .or. abs(eta_jet(2)).gt.eta_bjet_cut) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_lepM.lt.pT_lep_cut .OR. pT_lepP.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( abs(eta_lepM).gt.eta_lep_cut .or. abs(eta_lepP).gt.eta_lep_cut) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_miss.lt.pT_miss_cut ) then! note that this is ET and not pT
       applyPSCut = .true.
        RETURN
    endif


! binning
    NBin(1)= WhichBin(1,pT_lepP)
    NBin(2)= WhichBin(2,eta_lepP)
    NBin(3)= WhichBin(3,pT_miss)
    NBin(4)= WhichBin(4,phi_ll)
    NBin(5)= WhichBin(5,m_ll)
    NBin(6)= WhichBin(6,H_T)
    NBin(7)= WhichBin(7,MTW)
    NBin(8)= WhichBin(8,MTeff)





!-------------------------------------------------------
elseif( ObsSet.eq.33 .or. ObsSet.eq.35 .or.ObsSet.eq.43 .or. ObsSet.eq.45 ) then! set of observables for TTbar -> ttbar + ETmiss  in semi-hadr. top decays
! request at least two b-jets and 2 jets from hadronic W decay
    NObsJet_Tree = 4
    if( .not.(NJet.ge.NObsJet_Tree .and. any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif

    pT_jet(1) = get_PT( MomJet(1:4,1))
    pT_jet(2) = get_PT( MomJet(1:4,2))
    pT_jet(3) = get_PT( MomJet(1:4,3))
    pT_jet(4) = get_PT( MomJet(1:4,4))
    if( NJet.eq.5 ) pT_jet(5) = get_PT( MomJet(1:4,5))
    eta_jet(1)= get_ETA(MomJet(1:4,1))
    eta_jet(2)= get_ETA(MomJet(1:4,2))
    eta_jet(3)= get_ETA(MomJet(1:4,3))
    eta_jet(4)= get_ETA(MomJet(1:4,4))
    if( NJet.eq.5 ) eta_jet(5)= get_ETA(MomJet(1:4,5))

    pT_lepP  = get_PT(Mom(1:4,lep))
    eta_lepP = get_ETA(Mom(1:4,lep))

    MomMiss(1:4) = Mom(1:4,nu)+Mom(1:4,X0)+Mom(1:4,X0bar)
    pT_miss = get_ET( MomMiss(1:4) )! note that this is ET and not pT

    pT_Top = get_PT(Mom(1:4,t))
    eta_Top = get_ETA(Mom(1:4,t))


!    MTW = Get_MT(Mom(1:4,lep),MomMiss(1:4))!    ==  dsqrt(  2d0*pT_lepP*get_PT( MomMiss(1:4) )*(1d0-Get_CosPhi(Mom(1:4,lep),MomMiss(1:4))) )
   MTW = dsqrt(  2d0*pT_lepP*get_ET( MomMiss(1:4) )*(1d0-Get_CosPhi(Mom(1:4,lep),MomMiss(1:4))) )! let's define MTW with ET instead of pT in accordance with ATLAS

   MTeff = pt_jet(1) + pt_jet(2) + pt_jet(3) + pt_jet(4) + pT_lepP + pT_miss
   if( NJet.eq.5 ) MTeff = MTeff + pt_jet(5)

   H_T = pt_jet(1) + pt_jet(2) + pt_jet(3) + pt_jet(4) + MTW ! definition in accoradance with ATLAS


! check cuts
    if( pT_jet(1).lt.pT_bjet_cut .or. pT_jet(2).lt.pT_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( abs(eta_jet(1)).gt.eta_bjet_cut .or. abs(eta_jet(2)).gt.eta_bjet_cut) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_jet(3).lt.pT_jet_cut .or. pT_jet(4).lt.pT_jet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( abs(eta_jet(3)).gt.eta_jet_cut .or. abs(eta_jet(4)).gt.eta_jet_cut) then
        applyPSCut = .true.
        RETURN
    endif


    if( pT_lepP.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( abs(eta_lepP).gt.eta_lep_cut) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_miss.lt.pT_miss_cut ) then! note that this is ET and not pT
       applyPSCut = .true.
        RETURN
    endif
  
   if( MTW.lt.MTW_cut ) then
      applyPSCut = .true.
       RETURN
   endif


!  topness: step 1
   if( get_MInv2(MomJet(1:4,1)+Mom(1:4,lep)) .lt. get_MInv2(MomJet(1:4,2)+Mom(1:4,lep)) ) then! select the b-jet that belongs to the leptonic side
        ib_lep=1
        ib_had=2
   else
        ib_lep=2
        ib_had=1
   endif
   Stopness = (M_W**2 - Get_MInv2(Mom(1:4,lep)+MomMiss(1:4)))**2/(5*GeV)**4  & 
            + (M_top**2 - Get_MInv2(Mom(1:4,lep)+MomMiss(1:4)+MomJet(1:4,ib_lep)))**2/(15*GeV)**4

!  topness: step 2
   if( NJet.eq.4 ) then! if there are 4 jets, we assume that all come from the decay process
      Stopness =  Stopness  +  &
                + (M_top**2 - Get_MInv2(MomJet(1:4,3)+MomJet(1:4,4)+MomJet(1:4,ib_had)))**2/(15*GeV)**4 
   elseif( NJet.eq.5 ) then! if there are 5 jets, we assume that the hardest non-bjet comes from the production process
      Stopness =  Stopness  +  &
                + (M_top**2 - Get_MInv2(MomJet(1:4,4)+MomJet(1:4,5)+MomJet(1:4,ib_had)))**2/(15*GeV)**4 
   endif

!  topness: step 3
   if( NJet.eq.4 ) then! if there are 4 jets, we assume that all come from the decay process
      Stopness =  Stopness  +  &
               + (4*M_top**2 - Get_MInv2(Mom(1:4,lep)+MomMiss(1:4)+MomJet(1:4,1)+MomJet(1:4,2)+MomJet(1:4,3)+MomJet(1:4,4)))**2/(1000*GeV)**4
   elseif( NJet.eq.5 ) then! if there are 5 jets, we assume that the hardest non-bjet comes from the production process
      Stopness =  Stopness  +  &
               + (4*M_top**2 - Get_MInv2(Mom(1:4,lep)+MomMiss(1:4)+MomJet(1:4,1)+MomJet(1:4,2)+MomJet(1:4,4)+MomJet(1:4,5)))**2/(1000*GeV)**4
   endif


    
    rpT = (pT_jet(1)-pT_lepP)/(pT_jet(1)+pT_lepP)
    



! binning
    NBin(1)= WhichBin(1,pT_lepP)
    NBin(2)= WhichBin(2,eta_lepP)
    NBin(3)= WhichBin(3,pT_miss)
    NBin(4)= WhichBin(4,H_T)
    NBin(5)= WhichBin(5,pt_jet(NJet))! softest jet
    NBin(6)= WhichBin(6,MTW)
    NBin(7)= WhichBin(7,MTeff)
    NBin(8)= WhichBin(8,pT_Top)
    NBin(9)= WhichBin(9,eta_Top)
    NBin(10)= WhichBin(10,dlog(Stopness))
    NBin(11)= WhichBin(11,rpT)




!-------------------------------------------------------
else
  print *, "ObsSet not implemented",ObsSet
  stop
endif


return
END SUBROUTINE






SUBROUTINE Kinematics_TTBARZprime(NPlus1PS,MomExt,MomDK,applyPSCut,NBin)
use ModMisc
use ModParameters
implicit none
integer :: NumHadr
real(8) :: MomExt(:,:),MomDK(:,:),MomJet(1:4,1:8),MomAux1(1:4),MomAux2(1:4)
real(8) :: MomTops(1:4,1:2),MomBoost(1:4)
logical :: applyPSCut,NPlus1PS
integer :: NBin(:),PartList(1:7),JetList(1:7),NJet,NObsJet_Tree,NObsJet
real(8) :: y_Top,y_ATop,pT_Top,pT_ATop,M_TTbar,Dphi_TTbar,pT_Lept,y_Lept,Dphi_LL,M_LL
real(8) :: pT_LeptP, pT_LeptM, eta_LeptP, eta_LeptM, pT_b, pT_bbar, eta_b, eta_bbar
real(8) :: CosTheta_scatter,CosTheta_soper,CosTheta_star,Phi,Phibar,dPhiPlus,dPhiMinus,MassAux
real(8) :: MomHadr(1:4,1:7),MomLept(1:4,1:4),zeros(1:9),cosPhi,cosPhibar
real(8) :: MomTopsCMS(1:4,1:2),nx(2:4),ny(2:4),nz(2:4),MomLeptTRF(1:4),MomBeam(2:4)
integer :: i,j,n,k
real(8) :: sinPhi, sinPhibar, deltaSin, MomBeam1(1:4), MomBeam2(1:4), Mttb_cut,MTW
real(8) :: pT_bjet1, pT_bjet2, pT_lep1, pT_lep2, eta_bjet1, eta_bjet2, eta_lep1, eta_lep2,R_LepJet, eta_jet1, eta_jet2
real(8) :: MomLeptOrd(1:4,1:2), MomJetOrd(1:4,1:8), pT_miss, ET_miss, pT_jet1, pT_jet2, M_eff, HT, Minv_Lept,MomFatJet(1:4)

!DEC$ IF(_CheckMomenta .EQ.1)
   zeros(:) = 0d0
   if( .not.NPlus1PS ) then
        if(TopDecays.eq.0) then
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomExt(1:4,3) - MomExt(1:4,4)
        else
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomDK(1:4,1)-MomDK(1:4,2)-MomDK(1:4,3)-MomDK(1:4,4)-MomDK(1:4,5)-MomDK(1:4,6)
        endif
   else
        if(TopDecays.eq.0) then
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomExt(1:4,3) - MomExt(1:4,4) - MomExt(1:4,5)
        elseif( TopDecays.ne.0 .and. Correction.eq.2 ) then 
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomExt(1:4,3) - MomDK(1:4,1)-MomDK(1:4,2)-MomDK(1:4,3)-MomDK(1:4,4)-MomDK(1:4,5)-MomDK(1:4,6)
        elseif( TopDecays.ne.0 .and. Correction.eq.5 ) then 
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomDK(1:4,1)-MomDK(1:4,2)-MomDK(1:4,3)-MomDK(1:4,4)-MomDK(1:4,5)-MomDK(1:4,6) - MomDK(1:4,7)
        endif
   endif
   if( any(abs(zeros(1:4)/MomExt(1,1)).gt.1d-6) ) then
      print *, "ERROR: energy-momentum violation in SUBROUTINE Kinematics_TTBARZprime(",NPlus1PS,"): ",zeros(1:4)
      print *, "momenta dump:"
      print *, MomExt(1:4,1:5+NPlus1PS)
      print *, MomDK(1:4,1:6)
   endif
   if( .not. NPlus1PS ) then
        zeros(1) = (MomExt(1:4,3).dot.MomExt(1:4,3)) - m_Top**2
        zeros(2) = (MomExt(1:4,4).dot.MomExt(1:4,4)) - m_Top**2
        zeros(3) = 0d0
   elseif( NPlus1PS .and. Correction.eq.2 ) then 
        zeros(1) = (MomExt(1:4,4).dot.MomExt(1:4,4)) - m_Top**2
        zeros(2) = (MomExt(1:4,5).dot.MomExt(1:4,5)) - m_Top**2
        zeros(3)=  MomExt(1:4,3).dot.MomExt(1:4,3)
   elseif( NPlus1PS .and. Correction.eq.5 ) then 
        zeros(1) = (MomExt(1:4,3).dot.MomExt(1:4,3)) - m_Top**2
        zeros(2) = (MomExt(1:4,4).dot.MomExt(1:4,4)) - m_Top**2
        zeros(3)=  MomDK(1:4,7).dot.MomDK(1:4,7)
   endif
   zeros(4) =  MomDK(1:4,1).dot.MomDK(1:4,1)
   zeros(5) =  MomDK(1:4,2).dot.MomDK(1:4,2)
   zeros(6) =  MomDK(1:4,3).dot.MomDK(1:4,3)
   zeros(7) =  MomDK(1:4,4).dot.MomDK(1:4,4)
   zeros(8) =  MomDK(1:4,5).dot.MomDK(1:4,5)
   zeros(9) =  MomDK(1:4,6).dot.MomDK(1:4,6)
   if( TopDecays.eq.0 .and. any(abs(zeros(1:9)/MomExt(1,1)**2).gt.1d-6) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TTBARZprime(",NPlus1PS,"): ",zeros(1:9)
      print *, "momenta dump:"
      print *, MomExt(1:4,1:4+NPlus1PS)
      print *, MomDK(1:4,1:6)
   endif
   if( TopDecays.ne.0 .and. any(abs(zeros(1:9)/MomExt(1,1)**2).gt.1d-6) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TTBARZprime(",NPlus1PS,"): ",zeros(1:9)
      print *, "momenta dump:"
      print *, MomExt(1:4,1:4+NPlus1PS)
      print *, MomDK(1:4,1:7)
   endif
!DEC$ ENDIF
 


! required momentum order: MomExt(:,:): 1=In_left, 2=In_right, 3,4,...=light particles, N-1=ATop, N=Top
!                          MomDK(:,1:7) : 1=ABot, 2=lep-/q, 3=ANeu/qbar, 4=Bot, 5=lep+/qbar, 6=Neu/q, 7=(Glu)
! MomLept(:,1:4): 1=lep-, 2=ANeu, 3=lep+, 4=Neu


applyPSCut = .false.
NBin(1:NumHistograms) = 0


MomHadr(1:4,1:7) = 0d0
MomLept(1:4,1:4) = 0d0
PartList(1:7)=(/1,2,3,4,5,6,7/)




! separating momenta into hadron momenta and lepton momenta to which various procedures can be applied
MomHadr(1:4,1) = MomDK(1:4,1)  ! ABot
MomHadr(1:4,2) = MomDK(1:4,4)  ! Bot
!-------------------------------------------------------
if( TopDecays.eq.0 ) then  ! no top decays
   NumHadr = 0
   if( NPlus1PS ) then
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
   else
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
   endif
!-------------------------------------------------------
elseif( TopDecays.eq.1 ) then  ! full leptonic decay
  MomLept(1:4,1) = MomDK(1:4,2)   ! lep-
  MomLept(1:4,2) = MomDK(1:4,3)   ! ANeu
  MomLept(1:4,3) = MomDK(1:4,5)   ! lep+
  MomLept(1:4,4) = MomDK(1:4,6)   ! Neu

  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,3) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 3
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,3) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 3
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 2
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.2 ) then  ! full hadronic decay
  MomHadr(1:4,3) = MomDK(1:4,2)   ! qbar
  MomHadr(1:4,4) = MomDK(1:4,3)   ! q
  MomHadr(1:4,5) = MomDK(1:4,5)   ! qbar
  MomHadr(1:4,6) = MomDK(1:4,6)   ! q
  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,7) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 7
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,7) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 7
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 6
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.3 ) then  ! lept. Atop, hadr. top decay
  MomLept(1:4,1) = MomDK(1:4,2)  ! lep-
  MomLept(1:4,2) = MomDK(1:4,3)  ! ANeu
  MomHadr(1:4,3) = MomDK(1:4,5)  ! qbar
  MomHadr(1:4,4) = MomDK(1:4,6)  ! q
  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,5) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,5) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 4
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.4 ) then  ! hadr. Atop, lept. top decay
  MomHadr(1:4,3) = MomDK(1:4,2)  ! qbar
  MomHadr(1:4,4) = MomDK(1:4,3)  ! q
  MomLept(1:4,3) = MomDK(1:4,5)  ! lep+
  MomLept(1:4,4) = MomDK(1:4,6)  ! Neu
  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,5) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,5) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 4
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif
endif





if( ObsSet.eq.60 ) then! set of observables for ttb production without decays

    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))

    y_ATop = get_ETA(MomTops(1:4,1))
    y_Top  = get_ETA(MomTops(1:4,2))

    M_TTbar = get_MInv(MomTops(1:4,1)+MomTops(1:4,2))

    if (M_TTbar.lt.Mttb_cut) then
        applyPSCut = .true.
        RETURN
    endif


    Dphi_TTbar = dabs( Get_PHI(MomTops(1:4,1)) - Get_PHI(MomTops(1:4,2))  )
    if( Dphi_TTbar.gt.Pi ) Dphi_TTbar=2d0*Pi-Dphi_TTbar

!   scattering angle of the top quark in lab frame
    CosTheta_scatter = get_CosTheta( MomTops(1:4,2) )

!   see chinese paper, approximates scattering angle in ttb rest frame
    MomAux1(1:4) = MomTops(1:4,2)
    MomAux2(1:4) = MomExt(1:4,1)
    call boost(MomAux1(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)! this seems wrong!!!
    call boost(MomAux2(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)
    CosTheta_star = get_CosAlpha(MomAux1(1:4),MomAux2(1:4))



! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,y_ATop)
    NBin(4) = WhichBin(4,y_Top)
    NBin(5) = WhichBin(5,M_TTbar)
    NBin(6) = WhichBin(6,Dphi_TTbar)
    NBin(7) = WhichBin(7,CosTheta_scatter)
    NBin(8) = WhichBin(8,CosTheta_star)

    if( NPlus1PS ) NBin(9) = WhichBin(9,get_MInv(MomTops(1:4,1)+MomTops(1:4,2)+MomExt(1:4,3)))
    if( .not. NPlus1PS ) NBin(9) = WhichBin(9,get_MInv(MomTops(1:4,1)+MomTops(1:4,2)))

!-------------------------------------------------------
elseif( ObsSet.eq.61 ) then! set of observables for ttb production with di-lept. decays

!---------------------- (anti) kT jet algorithm ---------------------------------
    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr))
    call pT_order(2,MomJet(1:4,1:2))! pT-order b-jets first
    call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))! pT-order non-b jets
!--------------------------------------------------------------------------



   !-- define ordering for leptons, without modifying existing lepton array used below for angles
   MomLeptOrd(1:4,1) = MomLept(1:4,1)
   MomLeptOrd(1:4,2) = MomLept(1:4,3)

   call pT_order(2,MomLeptOrd)! pT-order for leptons

   !-- order all jets, irrespective whether they are b-jets or not
   MomJetOrd(1:4,1:8) = MomJet(1:4,1:8)
   call pT_order(NumHadr,MomJetOrd(1:4,1:NumHadr))


!   M_LL = get_MInv(MomLept(1:4,1)+MomLept(1:4,3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--F Phase space cuts

!   check that there are two b jets    
    if( .not.(any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif

    pT_bjet2 = get_PT(MomJet(1:4,2))
    pT_lep2  = get_PT(MomLeptOrd(1:4,2))

    eta_bjet1 = get_eta(MomJet(1:4,1))
    eta_bjet2 = get_eta(MomJet(1:4,2))

    eta_lep1 = get_eta(MomLept(1:4,1))
    eta_lep2 = get_eta(MomLeptOrd(1:4,2))

    if (pT_lep2.lt.pt_lep_cut) then
       applypscut = .true.
       RETURN
    endif

    if (dabs(eta_lep1).gt.eta_lep_cut .or. dabs(eta_lep2).gt.eta_lep_cut) then
       applypscut = .true.
       RETURN
    endif

    if (pT_bjet2.lt.pt_bjet_cut) then
       applypscut = .true.
       RETURN
    endif

    if (dabs(eta_bjet1).gt.eta_bjet_cut .or. dabs(eta_bjet2).gt.eta_bjet_cut) then
       applypscut = .true.
       RETURN
    endif

!-- If we pass the cuts, compute other observables
    pT_bjet1 = get_PT(MomJet(1:4,1))
    pT_lep1 = get_PT(MomLeptOrd(1:4,1))

    M_TTbar = get_MInv(MomTops(1:4,1)+MomTops(1:4,2))

    Dphi_LL = dabs( Get_PHI(MomLept(1:4,1)) - Get_PHI(MomLept(1:4,3))  )
    if( Dphi_LL.gt.Pi ) Dphi_LL=2d0*Pi-Dphi_LL

    pT_miss = get_pT(MomLept(1:4,2) + MomLept(1:4,4))
    Minv_Lept = get_MInv(MomLept(1:4,1)+MomLept(1:4,3))
    ET_miss = dsqrt(pT_miss**2 + Minv_Lept**2)

    pT_jet1 = get_PT(MomJetOrd(1:4,1))
    pT_jet2 = get_PT(MomJetOrd(1:4,2))

    HT = pT_lep1 + pT_lep2 + pT_jet1 + pT_jet2

    M_eff = ET_miss + HT

!   construct Baumgart-Tweedie angles // for non-zero transverse momentum of TTBAR, use Collins-Soper construction

    if (dabs(M_TTbar-M_Zpr).le.(100d0*GeV)) then ! compute angles only close to the resonance

!   step 1: boost top momenta into a TTBAR CMS frame ( pT+pTbar=0 )
    MomAux1(1:4) = MomTops(1:4,1)+MomTops(1:4,2)
    MomAux1(2:4) =-MomAux1(2:4)
    MassAux = dsqrt( (MomAux1(1:4).dot.MomAux1(1:4)) +1d-12 )
    MomTopsCMS(1:4,1:2) = MomTops(1:4,1:2)
    call boost(MomTopsCMS(1:4,1),MomAux1(1:4),MassAux)
    call boost(MomTopsCMS(1:4,2),MomAux1(1:4),MassAux)

!   step 2: construct coordinate system in a TTBAR CMS frame

    MomBeam1(1:4) = (/1d0,0d0,0d0,1d0/)
    call boost(MomBeam1(1:4),MomAux1(1:4),MassAux)

    MomBeam2(1:4) = (/1d0,0d0,0d0,-1d0/)
    call boost(MomBeam2(1:4),MomAux1(1:4),MassAux)

    MomBeam(2:4) = MomBeam1(2:4)-MomBeam2(2:4)
    MomBeam(2:4) = MomBeam(2:4)/dsqrt(MomBeam(2)**2+MomBeam(3)**2+MomBeam(4)**2)     ! Collins-Soper direction

    ny(2:4) = MomTopsCMS(2:4,2).cross.MomBeam(2:4)
    ny(2:4) = ny(2:4)/dsqrt(ny(2)**2+ny(3)**2+ny(4)**2)
    nz(2:4) = MomTopsCMS(2:4,2)/dsqrt(MomTopsCMS(2,2)**2+MomTopsCMS(3,2)**2+MomTopsCMS(4,2)**2)
    nx(2:4) = ny(2:4).cross.nz(2:4)

!   step 3: boost lepton momenta first into CMS frame, then into frame with ptop=0 and call it MomLeptTRF
    MomLeptTRF(1:4) = MomLept(1:4,3)
    call boost(MomLeptTRF(1:4),MomAux1(1:4),MassAux)
    MomAux2(1:4) = MomTopsCMS(1:4,2)
    MomAux2(2:4) =-MomAux2(2:4)
    call boost(MomLeptTRF(1:4),MomAux2(1:4),m_top)

!   step 4: project MomLeptTRF onto the nx-ny plane and calculate its angle with nx

    cosPhi = VectorProd(MomLeptTRF(2:4),nx(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2  & 
    + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )

    sinPhi = VectorProd(MomLeptTRF(2:4),ny(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2  & 
    + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )


!    Phi = dacos(cosPhi)

!------- repeat the same for the anti-top

!   step 3: boost lepton momenta into frame with patop=0 and call it MomLeptTRF
    MomLeptTRF(1:4) = MomLept(1:4,1)
    call boost(MomLeptTRF(1:4),MomAux1(1:4),MassAux)

    MomAux2(1:4) = MomTopsCMS(1:4,1)
    MomAux2(2:4) =-MomAux2(2:4)
    call boost(MomLeptTRF(1:4),MomAux2(1:4),m_top)

!   step 4: project MomLeptTRF onto the nx-ny plane and calculate its angle with n
    cosPhibar = VectorProd(MomLeptTRF(2:4),nx(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2 & 
      + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )

    sinPhibar = VectorProd(MomLeptTRF(2:4),ny(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2 & 
      + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )

!    Phibar = dacos(cosPhibar)

    deltaSin = sinPhi * cosPhibar - cosPhi * sinPhibar

    if (deltaSin .gt. 0d0) then
       dPhiMinus = dacos( cosPhi*cosPhibar + sinPhi*sinPhibar)
    else
       dPhiMinus = -dacos( cosPhi*cosPhibar + sinPhi*sinPhibar)
    endif

    dPhiPlus = dabs(dacos(cosPhi*cosPhibar - sinPhi*sinPhibar))

!    dPhiMinus = Phi-Phibar          ! -pi..pi
!    dPhiPlus  = dabs(Phi+Phibar)    ! 0...2pi

!    print *, phi/DblPi,phibar/DblPi
!    print *, dPhiMinus/DblPi,dPhiPlus/DblPi
!pause

 else
    dPhiPlus = 20d0
    dPhiMinus = 20d0
 endif

! binning
    NBin(1) = WhichBin(1,pT_lep1)
    NBin(2) = WhichBin(2,pT_bjet1)
    NBin(3) = WhichBin(3,pT_lep2)
    NBin(4) = WhichBin(4,pT_bjet2)
    NBin(5) = WhichBin(5,M_TTbar)
    NBin(6) = WhichBin(6,M_eff)
    NBin(7) = WhichBin(7,eta_lep1)
    NBin(8) = WhichBin(8,eta_bjet1)
    NBin(9) = WhichBin(9,eta_lep2)
    NBin(10) = WhichBin(10,eta_bjet2)
    NBin(11) = WhichBin(11,dPhi_LL)
    NBin(12) = WhichBin(12,dPhiMinus)
    NBin(13) = WhichBin(13,dPhiPlus)

!-------------------------------------------------------
elseif( ObsSet.eq.62 ) then! set of observables for ttb production with hadr. Atop, hadr. top decay

!---------------------- (anti) kT jet algorithm ---------------------------------
    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr))
    call pT_order(2,MomJet(1:4,1:2))! pT-order b-jets first
    call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))! pT-order non-b jets
!--------------------------------------------------------------------------


    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))

    y_ATop = get_ETA(MomTops(1:4,1))
    y_Top  = get_ETA(MomTops(1:4,2))

    M_TTbar = get_MInv(MomTops(1:4,1)+MomTops(1:4,2))

    Dphi_TTbar = dabs( Get_PHI(MomTops(1:4,1)) - Get_PHI(MomTops(1:4,2))  )
    if( Dphi_TTbar.gt.Pi ) Dphi_TTbar=2d0*Pi-Dphi_TTbar

!   scattering angle of the top quark
    CosTheta_scatter = get_CosTheta( MomTops(1:4,2) )

!   see chinese paper, approximates scattering angle
    MomAux1(1:4) = MomTops(1:4,2)
    MomAux2(1:4) = MomExt(1:4,1)
    call boost(MomAux1(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)
    call boost(MomAux2(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)
    CosTheta_star = get_CosAlpha(MomAux1(1:4),MomAux2(1:4))



! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,y_ATop)
    NBin(4) = WhichBin(4,y_Top)
    NBin(5) = WhichBin(5,M_TTbar)
    NBin(6) = WhichBin(6,Dphi_TTbar)
    NBin(7) = WhichBin(7,CosTheta_scatter)
    NBin(8) = WhichBin(8,CosTheta_star)



elseif( ObsSet.eq.64 ) then ! Zprime, semi-hadronic top decay (for factorization checks)

   NBin(:) = 1

elseif( ObsSet.eq.65 ) then ! Zprime, semi-hadronic top decay (for ATLAS analysis: James Ferrando)


!---------------------- (anti) kT jet algorithm ---------------------------------

    Rsep_Jet = 1.0d0

    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr))
    call pT_order(NumHadr,MomJet(1:4,1:NumHadr))! pT-order jets
!--------------------------------------------------------------------------

    if( NJet.lt.1 ) then
       applypscut = .true.
       RETURN
    endif
    pT_jet1 = get_PT(MomJet(1:4,1))
    if( pT_jet1.lt.pT_jet_cut ) then
       applypscut = .true.
       RETURN
    endif 
    eta_jet1 = get_eta(MomJet(1:4,1))
    if( abs(eta_jet1).gt.eta_jet_cut ) then
       applypscut = .true.
       RETURN
    endif 

    MomFatJet(1:4) = MomJet(1:4,1)



!---------------------- (anti) kT jet algorithm ---------------------------------

    Rsep_Jet = 0.4d0

    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr))
    call pT_order(NumHadr,MomJet(1:4,1:NumHadr))! pT-order jets
!--------------------------------------------------------------------------

   k=0
   do i=1,NJet
       if( get_R(MomJet(1:4,i),MomFatJet(1:4)).gt.1.5d0 .and. get_PT(MomJet(1:4,i)).gt.25d0*GeV .and. abs(get_eta(MomJet(1:4,i))).lt.2.5d0 ) then
            k=k+1
       endif
   enddo
   if( k.lt.1 ) then
       applypscut = .true.
       RETURN
    endif


    pT_lep1 = get_PT(MomLept(1:4,3))
    if( pT_lep1.lt.pT_lep_cut ) then
       applypscut = .true.
       RETURN
    endif 
    eta_lep1 = get_ETA(MomLept(1:4,3))
    if( abs(eta_lep1).lt.eta_lep_cut ) then
       applypscut = .true.
       RETURN
    endif 

    pT_miss = get_pT(MomLept(1:4,4))
    if( pT_miss.lt.pT_miss_cut ) then
       applypscut = .true.
       RETURN
    endif 

    MTW = Get_MT(MomLept(1:4,3),MomLept(1:4,4))
    if( MTW.lt.MTW_cut ) then
       applypscut = .true.
       RETURN
    endif 

    do i=1,NJet
      R_LepJet = get_R(MomJet(1:4,i),MomLept(1:4,3))
      if( R_LepJet.lt.Rsep_LepJet ) then
         applypscut = .true.
         RETURN
      endif 
    enddo






!   ------------------------
!   this is just for binning
    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))

    y_ATop = get_ETA(MomTops(1:4,1))
    y_Top  = get_ETA(MomTops(1:4,2))

    M_TTbar = get_MInv(MomTops(1:4,1)+MomTops(1:4,2))

    Dphi_TTbar = dabs( Get_PHI(MomTops(1:4,1)) - Get_PHI(MomTops(1:4,2))  )
    if( Dphi_TTbar.gt.Pi ) Dphi_TTbar=2d0*Pi-Dphi_TTbar

!   scattering angle of the top quark
    CosTheta_scatter = get_CosTheta( MomTops(1:4,2) )

!   see chinese paper, approximates scattering angle
    MomAux1(1:4) = MomTops(1:4,2)
    MomAux2(1:4) = MomExt(1:4,1)
    call boost(MomAux1(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)
    call boost(MomAux2(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)
    CosTheta_star = get_CosAlpha(MomAux1(1:4),MomAux2(1:4))



! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,y_ATop)
    NBin(4) = WhichBin(4,y_Top)
    NBin(5) = WhichBin(5,M_TTbar)
    NBin(6) = WhichBin(6,Dphi_TTbar)
    NBin(7) = WhichBin(7,CosTheta_scatter)
    NBin(8) = WhichBin(8,CosTheta_star)





elseif( ObsSet.eq.66 ) then ! Zprime, semi-hadronic top decay (for CMS analysis: Roman Kogler)


!---------------------- (anti) kT jet algorithm ---------------------------------
    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr))
    call pT_order(NumHadr,MomJet(1:4,1:NumHadr))! pT-order jets
!--------------------------------------------------------------------------

   if( NJet.lt.2 ) then
       applypscut = .true.
       RETURN
    endif

    pT_jet1 = get_PT(MomJet(1:4,1))
    pT_jet2 = get_PT(MomJet(1:4,2))
    if( pT_jet1.lt.pT_jet_cut ) then
       applypscut = .true.
       RETURN
    endif 
    if( pT_jet2.lt.50d0*GeV ) then
       applypscut = .true.
       RETURN
    endif
   

    eta_jet1 = get_eta(MomJet(1:4,1))
    eta_jet2 = get_eta(MomJet(1:4,2))
    if (dabs(eta_jet1).gt.eta_jet_cut .or. dabs(eta_jet2).gt.eta_jet_cut) then
       applypscut = .true.
       RETURN
    endif

    pT_lep1 = get_PT(MomLept(1:4,3))
    if( pT_lep1.lt.pT_lep_cut ) then
       applypscut = .true.
       RETURN
    endif 
    eta_lep1 = get_ETA(MomLept(1:4,3))
    if( abs(eta_lep1).lt.eta_lep_cut ) then
       applypscut = .true.
       RETURN
    endif 

    pT_miss = get_pT(MomLept(1:4,4))
    if( pT_miss.lt.pT_miss_cut ) then
       applypscut = .true.
       RETURN
    endif 

    HT = pT_lep1 + pT_miss
    if( HT.lt.HT_cut ) then
       applypscut = .true.
       RETURN
    endif 

    do i=1,NJet
      R_LepJet = get_R(MomJet(1:4,i),MomLept(1:4,3))
      if( R_LepJet.lt.Rsep_LepJet ) then
         applypscut = .true.
         RETURN
      endif 
    enddo

!   --------------------------------
!   this is just for binning
    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))

    y_ATop = get_ETA(MomTops(1:4,1))
    y_Top  = get_ETA(MomTops(1:4,2))

    M_TTbar = get_MInv(MomTops(1:4,1)+MomTops(1:4,2))

    Dphi_TTbar = dabs( Get_PHI(MomTops(1:4,1)) - Get_PHI(MomTops(1:4,2))  )
    if( Dphi_TTbar.gt.Pi ) Dphi_TTbar=2d0*Pi-Dphi_TTbar

!   scattering angle of the top quark
    CosTheta_scatter = get_CosTheta( MomTops(1:4,2) )

!   see chinese paper, approximates scattering angle
    MomAux1(1:4) = MomTops(1:4,2)
    MomAux2(1:4) = MomExt(1:4,1)
    call boost(MomAux1(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)
    call boost(MomAux2(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)
    CosTheta_star = get_CosAlpha(MomAux1(1:4),MomAux2(1:4))



! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,y_ATop)
    NBin(4) = WhichBin(4,y_Top)
    NBin(5) = WhichBin(5,M_TTbar)
    NBin(6) = WhichBin(6,Dphi_TTbar)
    NBin(7) = WhichBin(7,CosTheta_scatter)
    NBin(8) = WhichBin(8,CosTheta_star)



endif

return
END SUBROUTINE















RECURSIVE SUBROUTINE JetAlgo_kt(Rsep_jet,PartonList,MomParton,NJet,JetList,MomJet)  ! initial call must have NJet=0 and MomJet(1:4,:) = MomPartons(1:4,:)
use ModMisc
use ModParameters
implicit none
integer :: PartonList(:), JetList(:)
real(8) :: MomParton(:,:)
integer :: NJet
real(8) :: MomJet(:,:)
real(8) :: Rsep_jet
integer :: NParton,i,j,k,dii_minp,dij_minp,ij(1:120,1:2)  ! max. partons=8, (8-1)! =5040  ! max. partons=6, (6-1)! =120
real(8) :: dii(1:6),dij(1:120),Rij,dii_min,dij_min!,eta(1:5),phi(1:5)


NParton = size(PartonList)
if( NParton.eq.0 ) then
    return
elseif( NParton.eq.1 ) then
!   print *, "HARD JET", PartonList(NParton)
   NJet = NJet +1
   JetList(NJet) = PartonList(NParton)
   return
endif

!generate dii, dij
do i=1,NParton
   dii(i) = get_PT2(MomParton(1:4,PartonList(i)))**AlgoType
enddo

k=0
do i=1,NParton-1
do j=i+1,NParton
   k = k+1
   Rij = get_R( MomParton(1:4,PartonList(i)), MomParton(1:4,PartonList(j)) )
   dij(k) = dmin1(dii(i),dii(j)) * (Rij/Rsep_jet)**2
   ij(k,1)=i
   ij(k,2)=j
enddo
enddo
!!print *, dii(1:NParton)
!!print *, dij(1:k)


! find minima
dii_min=dii(1)
dii_minp=1
do i=2,NParton
  if( dii(i).lt.dii_min ) then
    dii_min = dii(i)
    dii_minp= i
  endif
enddo

dij_min=dij(1)
dij_minp=1
do i=2,k
  if( dij(i).lt.dij_min ) then
    dij_min = dij(i)
    dij_minp= i
  endif
enddo


if( dii_min.lt.dij_min ) then ! no recombination
   NJet = NJet +1
   k=0
!   print *, "HARD JET", PartonList(dii_minp)
   JetList(NJet) = PartonList(dii_minp)
   MomJet(1:4,PartonList(dii_minp)) = MomParton(1:4,PartonList(dii_minp))
   do i=1,NParton!  remove momenta dii_minp from parton list
      if( i.eq.dii_minp ) cycle
      k=k+1
      PartonList(k) = PartonList(i)
   enddo
   call JetAlgo_kt(Rsep_jet,PartonList(1:k),MomParton(:,:),NJet,JetList(:),MomJet(:,:))
else ! recombination of dij(dij_min)
   if( RecombPrescr.eq.1 ) then
        MomJet(1:4,PartonList(ij(dij_minp,1))) = EllisSoperComb(MomJet(1:4,PartonList(ij(dij_minp,1))),MomJet(1:4,PartonList(ij(dij_minp,2))))   ! Ellis-Soper combination
   else
        MomJet(1:4,PartonList(ij(dij_minp,1))) = MomJet(1:4,PartonList(ij(dij_minp,1))) + MomJet(1:4,PartonList(ij(dij_minp,2)))                 ! four vector addition
   endif
   MomJet(1:4,PartonList(ij(dij_minp,2))) = 0d0
   k=0
!   print *, "RECOMB.",PartonList(ij(dij_minp,1)),PartonList(ij(dij_minp,2))
   do i=1,NParton!  remove momenta j from parton list
      if( i.eq.ij(dij_minp,2) ) then
         cycle
      endif
      k=k+1
      PartonList(k) = PartonList(i)
   enddo
   call JetAlgo_kt(Rsep_jet,PartonList(1:k),MomJet(:,:),NJet,JetList(:),MomJet(:,:))
endif

return
END SUBROUTINE




FUNCTION EllisSoperComb(pi,pj)
use ModMisc
implicit none
real(8), intent(in) :: pi(1:4),pj(1:4)
real(8) :: EllisSoperComb(1:4)
real(8) :: ET_i,ET_j,ET_k,eta_i,eta_j,eta_k,phi_i,phi_j,phi_k,theta_k,E_k

   ET_i = get_PT(pi(1:4))   ! ET=pT as defined in hep-ph/9305266; Ellis,Soper
   ET_j = get_PT(pj(1:4))
   eta_i = get_PseudoEta(pi(1:4))  ! check pseudo_eta = eta for m=0
   eta_j = get_PseudoEta(pj(1:4))


   phi_i = get_Phi(pi(1:4))
   phi_j = get_Phi(pj(1:4))

   ET_k = ET_i + ET_j
   eta_k = (ET_i*eta_i + ET_j*eta_j)/ET_k
   phi_k = (ET_i*phi_i + ET_j*phi_j)/ET_k
   theta_k = 2d0*datan( dexp(-eta_k) )   ! datan( (0,inf) ) = (0,pi/2)  --> theta_k = (0,pi)
   E_k = ET_k/dsin(theta_k)

   EllisSoperComb(1) = E_k
   EllisSoperComb(2) = E_k * ( dsin(theta_k) * dcos(phi_k) )     ! check E_k prefactor
   EllisSoperComb(3) = E_k * ( dsin(theta_k) * dsin(phi_k) )
   EllisSoperComb(4) = E_k * ( dcos(theta_k) )

return
END FUNCTION




FUNCTION FrixioneIsolated(MomPho,Riso,NumParton,MomParton)! Frixione photon isolation cut, arXiv:hep-ph/9801442
use ModMisc
implicit none
logical :: FrixioneIsolated
integer :: NumParton,i,j
real(8) :: MomParton(1:4,1:NumParton),MomPho(1:4),Riso
real(8) :: Chi,Esum,OneMinusCosRip,OneMinusCosRiso,Rip,Rjp

  FrixioneIsolated = .true.
  do i=1,NumParton
      Rip = dmin1(get_R(MomParton(1:4,i),MomPho(1:4)),Riso)
      OneMinusCosRip  = Rip**2/2d0*(1d0 - Rip**2/12d0 + Rip**4/360d0 - Rip**6/20160d0)
      OneMinusCosRiso = Riso**2/2d0*(1d0 - Riso**2/12d0 + Riso**4/360d0 - Riso**6/20160d0)
      Chi = get_ET(MomPho(1:4))*(OneMinusCosRip)/(OneMinusCosRiso)

      Esum = 0d0
      do j=1,NumParton
         Rjp = get_R(MomParton(1:4,j),MomPho(1:4))
         if( Rjp.gt.Rip ) cycle
         if( Rjp.gt.Riso) cycle
         Esum = Esum + get_ET(MomParton(1:4,j))
      enddo

      if( Chi.lt.Esum ) then
         FrixioneIsolated = .false.
         return
      endif
  enddo


!12345 continue
!if(  IsolateFrix(MomPho,Riso,NumParton,MomParton(1:4,1:NumParton)).eq.FrixioneIsolated   ) then
!  print *, "Andreas", .not.IsolateFrix(MomPho,Riso,NumParton,MomParton(1:4,1:NumParton))
!  print *, "Markus ",FrixioneIsolated
!  print *, Chi,Esum
!  pause
!endif

return
END FUNCTION




function IsolateFrix(kpho,delta0,numHad,MomHad)! Andreas' implementation
use ModMisc
  implicit none
  integer numHad, i, numH
  real(8) :: MomHad(4,numHad), kpho(4), delta, ET(numHad),RiPho(numHad)
  real(8) :: delta0, LHS, RHS
  logical :: IsolateFrix
  IsolateFrix = .false.

  do i=0,1000
    delta = dble(i)/1000.d0*delta0
    LHS =0.d0
    do numH=1,numHad
        ET(numH) = dsqrt(MomHad(2,numH)**2 + MomHad(3,numH)**2)
        RiPho(numH)  = get_R(MomHad(1:4,numH),kpho)
        LHS = LHS + ET(numH)*StepFunc(delta-RiPho(numH))
    enddo
    RHS = dsqrt(kpho(2)**2 + kpho(3)**2)*((1.d0-dcos(delta))/(1.d0-dcos(delta0)))
    if(LHS .gt. RHS) then
        IsolateFrix = .true.
        return
    endif
  enddo

  return
end function







FUNCTION getHelicity(yRnd)
use ModParameters
use ModProcess
implicit none
integer :: getHelicity(1:2)
real(8) :: yRnd

  IF(HelSampling) THEN
    getHelicity(1:2)=yRnd*dble(NumHelicities) + 1
  ELSE
    getHelicity(1)=1
    getHelicity(2)=NumHelicities
  ENDIF

return
END FUNCTION



FUNCTION MomCrossing(Mom)
use ModProcess
use ModParameters
implicit none
real(8) :: Mom(1:4,1:NumExtParticles),MomCrossing
integer :: NPart
!real(8) :: QuarkCrossing=-1d0, SpinAvg=1d0/4d0, QuarkColAvg=1d0/3d0, GluonColAvg=1d0/8d0

do NPart=1,NumExtParticles
   ExtParticle(NPart)%Mom(1:4) = sign(1,Crossing(NPart)) * Mom(1:4,abs(Crossing(NPart)))
   ExtParticle(NPart)%Mom(5:8) = (0d0,0d0)
enddo
MomCrossing = AvgFactor

return
END FUNCTION





SUBROUTINE HelCrossing(Hel)
use ModProcess
use ModParameters
implicit none
integer :: Hel(1:NumExtParticles),n


do n=1,NumExtParticles
  ExtParticle(n)%Helicity = Hel(n)
enddo


return
END SUBROUTINE








SUBROUTINE CheckSing(MomExt,applySingCut)
use ModMisc
use ModParameters
use ModProcess
implicit none
logical :: applySingCut
real(8) :: MomExt(:,:),s12,s13,s23,s14,s24,s34,E3,E4


!     applySingCut=.false.
!     if( dsqrt(MomExt(2,3)**2 + MomExt(3,3)**2).lt.1d0 ) then
!         applySingCut=.true.
!         return
!     endif

IF( NumExtParticles.EQ.5 .AND. CORRECTION.EQ.2 ) THEN
    s12 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,2))
    s13 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,3))
    s23 = 2d0*(MomExt(1:4,2).dot.MomExt(1:4,3))
    E3 = MomExt(1,3)
    applySingCut=.false.
    if( dabs(s13/s12).lt.1d-9 ) then
        applySingCut=.true.
        return
    endif
    if( dabs(s23/s12).lt.1d-9 ) then
        applySingCut=.true.
        return
    endif
    if( dabs(E3/(2d0*MomExt(1,1))).lt.1d-5 ) then
        applySingCut=.true.
        return
    endif
ELSEIF( NumExtParticles.EQ.6 .AND. CORRECTION.EQ.2 ) THEN
    s12 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,2))
    s13 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,3))
    s23 = 2d0*(MomExt(1:4,2).dot.MomExt(1:4,3))
    s14 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,4))
    s24 = 2d0*(MomExt(1:4,2).dot.MomExt(1:4,4))
    s34 = 2d0*(MomExt(1:4,3).dot.MomExt(1:4,4))
    E3 = MomExt(1,3)
    E4 = MomExt(1,4)
    applySingCut=.false.
    if( dabs(s13/s12).lt.1d-9 ) then
        applySingCut=.true.
        return
    endif
    if( dabs(s23/s12).lt.1d-9 ) then
        applySingCut=.true.
        return
    endif
    if( dabs(s14/s12).lt.1d-9 ) then
        applySingCut=.true.
        return
    endif
    if( dabs(s24/s12).lt.1d-9 ) then
        applySingCut=.true.
        return
    endif
    if( dabs(s34/s12).lt.1d-9 ) then
        applySingCut=.true.
        return
    endif
    if( dabs(E3/(2d0*MomExt(1,1))).lt.1d-5 ) then
        applySingCut=.true.
        return
    endif
    if( dabs(E4/(2d0*MomExt(1,1))).lt.1d-5 ) then
        applySingCut=.true.
        return
    endif


ELSEIF( CORRECTION.EQ.5 ) THEN
    s14 = MomExt(1:4,1).dot.MomExt(1:4,4)
    E4 = dabs(MomExt(1,4))
    applySingCut=.false.
    if( dabs(s14)/m_Top**2 .lt. 1d-8 ) then
!     if( dabs(s14)/m_Top**2 .lt. 1d-4 ) then
        applySingCut=.true.
        return
    endif
    if( E4/m_Top .lt. 1d-4 ) then
!     if( E4/m_Top .lt. 1d-2 ) then
        applySingCut=.true.
        return
    endif
ENDIF

return
END SUBROUTINE






FUNCTION WhichBin(NHisto,Value)
implicit none
integer :: WhichBin,NHisto
real(8) :: Value

   WhichBin = (Value-Histo(NHisto)%LowVal)/Histo(NHisto)%BinSize + 1
   if( WhichBin.lt.0 ) then
      WhichBin = 0
   elseif( WhichBin.gt.Histo(NHisto)%NBins ) then
      WhichBin = Histo(NHisto)%NBins+1
   endif

RETURN
END FUNCTION





SUBROUTINE IntoHisto(NHisto,NBin,Value)
use ModMisc
implicit none
integer :: NHisto,NBin
real(8) :: Value

    if( IsNaN(Value) ) return

     Histo(NHisto)%Value(NBin) = Histo(NHisto)%Value(NBin)  + Value
     Histo(NHisto)%Value2(NBin)= Histo(NHisto)%Value2(NBin) + Value**2
     Histo(NHisto)%Hits(NBin)  = Histo(NHisto)%Hits(NBin)+1

RETURN
END SUBROUTINE





SUBROUTINE setPolarizations()
use ModProcess
use ModMisc
use ModParameters
implicit none
integer :: NPart

   do NPart=1,NumExtParticles

      if( IsAQuark(ExtParticle(NPart)%PartType) .and. ExtParticle(NPart)%PartType.lt.0 .and. abs(ExtParticle(NPart)%Helicity).eq.1 ) then
         call vSpi(ExtParticle(NPart)%Mom(1:4),ExtParticle(NPart)%Mass,ExtParticle(NPart)%Helicity,ExtParticle(NPart)%Pol(1:4))
         ExtParticle(NPart)%Pol(5:16) = (0d0,0d0)
      elseif( IsAQuark(ExtParticle(NPart)%PartType) .and. ExtParticle(NPart)%PartType.gt.0 .and. abs(ExtParticle(NPart)%Helicity).eq.1 ) then
         call ubarSpi(ExtParticle(NPart)%Mom(1:4),ExtParticle(NPart)%Mass,ExtParticle(NPart)%Helicity,ExtParticle(NPart)%Pol(1:4))
         ExtParticle(NPart)%Pol(5:16) = (0d0,0d0)
      elseif( ExtParticle(NPart)%PartType.eq.Glu_ ) then
         call pol_mless(ExtParticle(NPart)%Mom(1:4),ExtParticle(NPart)%Helicity,ExtParticle(NPart)%Pol(1:4))
         ExtParticle(NPart)%Pol(5:16) = (0d0,0d0)
      elseif ( ExtParticle(NPart)%PartType.eq.Z0_ .and. .not. ZDecays.gt.0 ) then 
         call pol_massSR(ExtParticle(NPart)%Mom(1:4),ExtParticle(NPart)%Mass,ExtParticle(NPart)%Helicity,ExtParticle(NPart)%Pol(1:4))! on-shell Z-boson polarization
         ExtParticle(NPart)%Pol(5:16) = (0d0,0d0)
      elseif( ExtParticle(NPart)%PartType.eq.Pho_  ) then
         call pol_mless(ExtParticle(NPart)%Mom(1:4),ExtParticle(NPart)%Helicity,ExtParticle(NPart)%Pol(1:4))
         ExtParticle(NPart)%Pol(5:16) = (0d0,0d0)
      endif
   
   enddo

!       check gauge invariance
!        ExtParticle(4)%Pol(1:4) = ExtParticle(4)%Mom(1:4);       print *, "gauge invariance check"

END SUBROUTINE



SUBROUTINE SetPropagators()
use ModProcess
implicit none
integer :: NPrimAmp,Prop

   do NPrimAmp=1,NumPrimAmps
         PrimAmps(NPrimAmp)%IntPart(1)%Mom(1:5) = (0d0,0d0)
         do Prop=2,NumExtParticles
            PrimAmps(NPrimAmp)%IntPart(Prop)%Mom(1:4) = PrimAmps(NPrimAmp)%IntPart(Prop-1)%Mom(1:4)  &
                                                      + dcmplx( ExtParticle( PrimAmps(NPrimAmp)%ExtLine(Prop-1) )%Mom(1:4) )
            PrimAmps(NPrimAmp)%IntPart(Prop)%Mom( 5 ) = (0d0,0d0)
         enddo
   enddo
END SUBROUTINE






SUBROUTINE boost2Lab(x1,x2,NumPart,Mom)
implicit none
real(8) Mom(1:4,1:NumPart)
real(8) x1,x2
real(8) gamma,betagamma,MomTmp1,MomTmp4
integer :: i,NumPart
!   beta  = (x2-x1)/(x1+x2)
!   gamma = 1d0/dsqrt(1d0-beta**2)
!   betagamma = beta*gamma

  gamma     = (x1+x2)/2d0/dsqrt(x1*x2)
  betagamma = (x2-x1)/2d0/dsqrt(x1*x2)

!   Mom(1,1)= x1*ColliderEnergy/2d0
!   Mom(2,1)= 0d0
!   Mom(3,1)= 0d0
!   Mom(4,1)= x1*ColliderEnergy/2d0
!
!   Mom(1,2)= x2*ColliderEnergy/2d0
!   Mom(2,2)= 0d0
!   Mom(3,2)= 0d0
!   Mom(4,2)=-x2*ColliderEnergy/2d0

  do i=1,NumPart
      MomTmp1=Mom(1,i)
      MomTmp4=Mom(4,i)
      Mom(1,i)= gamma*MomTmp1 - betagamma*MomTmp4
      Mom(4,i)= gamma*MomTmp4 - betagamma*MomTmp1
  enddo

RETURN
END SUBROUTINE






SUBROUTINE PDFMapping(MapType,yRnd,eta1,eta2,Ehat,sHatJacobi)
use ModParameters
use ModMisc
implicit none
integer :: MapType
real(8) :: yRnd(1:2),eta1,eta2,EHat,sHatJacobi,tau,nPotMap,z,sbar,fmax

  if( MapType.eq.1 ) then!  no mapping
      eta1 = yRnd(1)
      eta2 = yRnd(2)
      sHatJacobi = 1d0
  elseif( MapType.eq.2 ) then!  exponential mapping
      tau = (2d0*m_Top/Collider_Energy)**2
      eta1 = tau**yRnd(1)
      eta2 = tau**( (1d0-yRnd(1))*yRnd(2) )
      sHatJacobi = dlog(tau)**2*(1d0-yRnd(1))*eta1*eta2
  elseif( MapType.eq.3 ) then!  linear mapping
      tau = (2d0*m_Top/Collider_Energy)**2
      eta1 = (1d0-tau)*yRnd(1) + tau
      eta2 = ((1d0-tau)*yRnd(1))/((1d0-tau)*yRnd(1)+tau)*yRnd(2) + tau/((1d0-tau)*yRnd(1)+tau)
      sHatJacobi = (1d0-tau)*((1d0-tau)*yRnd(1))/((1d0-tau)*yRnd(1)+tau)
  elseif( MapType.eq.4 ) then!  MCFM mapping
      tau = dexp(dlog(((2d0*m_Top/Collider_Energy)**2))*yRnd(1))
      eta1 = dsqrt(tau)*dexp(0.5d0*dlog(tau)*(1d0-2d0*yRnd(2)))
      eta2 = dsqrt(tau)/dexp(0.5d0*dlog(tau)*(1d0-2d0*yRnd(2)))
      sHatJacobi = dlog(((2d0*m_Top/Collider_Energy)**2))*tau*dlog(tau)
  elseif( MapType.eq.5 ) then!  nPotMap mapping
      nPotMap = 0.5d0
      tau = (2d0*m_Top/Collider_Energy)**2
      yRnd(1) = yRnd(1)**nPotMap
      yRnd(2) = yRnd(2)**nPotMap
      eta1 = (1d0-tau) * yRnd(1) + tau
      eta2 = ((1d0-tau)*yRnd(1))/((1d0-tau)*yRnd(1)+tau)*yRnd(2) + tau/((1d0-tau)*yRnd(1)+tau)
      sHatJacobi=(1d0-tau)*((1d0-tau)*yRnd(1))/((1d0-tau)*yRnd(1)+tau)*nPotMap**2*((yRnd(1)*yRnd(2))**(1d0/nPotMap))**(nPotMap-1d0)
  elseif( MapType.eq.10 ) then!  Breit-Wigner mapping
!       fmax = 1d0/m_Grav/Ga_Grav * ( datan((Collider_Energy**2-m_Grav**2)/m_Grav/Ga_Grav) - datan(-m_Grav/Ga_Grav) )
!       sbar = m_Grav*Ga_Grav * dtan(fmax*yRnd(1)*m_Grav*Ga_Grav - atan(m_Grav/Ga_Grav) ) + m_Grav**2
!       z = sbar/Collider_Energy**2
!       eta1 = z + (1d0-z)*yRnd(2)
!       eta2 = z/eta1
!       sHatJacobi = fmax/Collider_Energy**2 * (1d0-z)/eta1  * ( (Collider_Energy**2*eta1*eta2 - m_Grav**2)**2 + m_Grav**2*Ga_Grav**2 )

!!! Zprime section !!!
  elseif ( MapType.eq.62) then ! BW for Zprime
     fmax = 1d0/M_Zpr/Ga_Zpr * ( datan((Collider_Energy**2-M_Zpr**2)/M_Zpr/Ga_Zpr) - datan(-M_Zpr/Ga_Zpr) )
     sbar = M_Zpr*Ga_Zpr * dtan(fmax*yRnd(1)*M_Zpr*Ga_Zpr - atan(M_Zpr/Ga_Zpr) ) + M_Zpr**2
     z = sbar/Collider_Energy**2
     eta1 = z + (1d0-z)*yRnd(2)
     eta2 = z/eta1
     sHatJacobi = fmax/Collider_Energy**2 * (1d0-z)/eta1  * ( (Collider_Energy**2*eta1*eta2 - M_Zpr**2)**2 + M_Zpr**2*Ga_Zpr**2 )
!!! End Zprime section !!!

  else

      call Error("PDF mapping not available")
  endif
  EHat = Collider_Energy*dsqrt(eta1*eta2)

RETURN
END SUBROUTINE



SUBROUTINE setPDFs(x1,x2,MuFac,pdf)
use ModParameters
implicit none
real(8) :: x1,x2,PDFScale,MuFac
real(8) :: upv(1:2),dnv(1:2),usea(1:2),dsea(1:2),str(1:2),sbar(1:2),chm(1:2),cbar(1:2),bot(1:2),bbar(1:2),glu(1:2),phot
integer,parameter :: swPDF_u=1, swPDF_d=1, swPDF_c=1, swPDF_s=1, swPDF_b=1, swPDF_g=1
real(8) :: pdf(-6:6,1:2)

        PDFScale=MuFac*100d0


!   MRSW PDFS
IF( PDFSET.EQ.1 .AND. NLOPARAM.LE.1) THEN
        if( x1.lt.1d0 ) then ! this is needed for integrated dipole routines, where eta/z appears
            call mrstlo(x1,PDFScale,1,upv(1),dnv(1),usea(1),dsea(1),str(1),chm(1),bot(1),glu(1))
! RR added
            sbar(1)=str(1)
            cbar(1)=chm(1)
            bbar(1)=bot(1)
!             call mrs96(x1,PDFScale,2,upv(1),dnv(1),usea(1),dsea(1),str(1),chm(1),bot(1),glu(1))
!             call GetAllPDFs("mstw2008lo",0,x1,PDFScale,upv(1),dnv(1),usea(1),dsea(1),str(1),sbar(1),chm(1),cbar(1),bot(1),bbar(1),glu(1),phot)
        else
            upv(1) = 0d0
            dnv(1) = 0d0
            usea(1)= 0d0
            dsea(1)= 0d0
            str(1) = 0d0
            sbar(1)= 0d0
            chm(1) = 0d0
            cbar(1)= 0d0
            bot(1) = 0d0
            bbar(1)= 0d0
            glu(1) = 0d0
        endif
        if( x2.lt.1d0 ) then
            call mrstlo(x2,PDFScale,1,upv(2),dnv(2),usea(2),dsea(2),str(2),chm(2),bot(2),glu(2))
! RR added
            sbar(2)=str(2)
            cbar(2)=chm(2)
            bbar(2)=bot(2)

!             call mrs96(x2,PDFScale,2,upv(2),dnv(2),usea(2),dsea(2),str(2),chm(2),bot(2),glu(2))
!             call GetAllPDFs("mstw2008lo",0,x2,PDFScale,upv(2),dnv(2),usea(2),dsea(2),str(2),sbar(2),chm(2),cbar(2),bot(2),bbar(2),glu(2),phot)
        else
            upv(2) = 0d0
            dnv(2) = 0d0
            usea(2)= 0d0
            dsea(2)= 0d0
            str(2) = 0d0
            sbar(2)= 0d0
            chm(2) = 0d0
            cbar(2)= 0d0
            bot(2) = 0d0
            bbar(2)= 0d0
            glu(2) = 0d0
        endif
ELSEIF( PDFSET.EQ.1 .AND. NLOPARAM.EQ.2) THEN
        if( x1.lt.1d0 ) then ! this is needed for integrated dipole routines, where eta/z appears
!             call mrst2004(x1,PDFScale,1,upv(1),dnv(1),usea(1),dsea(1),str(1),chm(1),bot(1),glu(1))
            call mrst2001(x1,PDFScale,1,upv(1),dnv(1),usea(1),dsea(1),str(1),chm(1),bot(1),glu(1))
!             call GetAllPDFs("mstw2008nlo",0,x1,PDFScale,upv(1),dnv(1),usea(1),dsea(1),str(1),sbar(1),chm(1),cbar(1),bot(1),bbar(1),glu(1),phot)
        else
            upv(1) = 0d0
            dnv(1) = 0d0
            usea(1)= 0d0
            dsea(1)= 0d0
            str(1) = 0d0
            sbar(1)= 0d0
            chm(1) = 0d0
            cbar(1)= 0d0
            bot(1) = 0d0
            bbar(1)= 0d0
            glu(1) = 0d0
        endif
        if( x2.lt.1d0 ) then
!             call mrst2004(x2,PDFScale,1,upv(2),dnv(2),usea(2),dsea(2),str(2),chm(2),bot(2),glu(2))
            call mrst2001(x2,PDFScale,1,upv(2),dnv(2),usea(2),dsea(2),str(2),chm(2),bot(2),glu(2))
!             call GetAllPDFs("mstw2008nlo",0,x2,PDFScale,upv(2),dnv(2),usea(2),dsea(2),str(2),sbar(2),chm(2),cbar(2),bot(2),bbar(2),glu(2),phot)
        else
            upv(2) = 0d0
            dnv(2) = 0d0
            usea(2)= 0d0
            dsea(2)= 0d0
            str(2) = 0d0
            sbar(2)= 0d0
            chm(2) = 0d0
            cbar(2)= 0d0
            bot(2) = 0d0
            bbar(2)= 0d0
            glu(2) = 0d0
        endif


!   CTEQ PDFS
ELSEIF( PDFSET.EQ.2 ) THEN
        if( x1.lt.1d0 ) then ! this is needed for integrated dipole routines
            call cteq6(x1,PDFScale,99,upv(1),dnv(1),usea(1),dsea(1),str(1),chm(1),bot(1),glu(1))
            sbar(1) = str(1)
            cbar(1) = chm(1)
            bbar(1) = bot(1)
        else
            upv(1) = 0d0
            dnv(1) = 0d0
            usea(1)= 0d0
            dsea(1)= 0d0
            str(1) = 0d0
            sbar(1)= 0d0
            chm(1) = 0d0
            cbar(1)= 0d0
            bot(1) = 0d0
            bbar(1)= 0d0
            glu(1) = 0d0
        endif
        if( x2.lt.1d0 )then
            call cteq6(x2,PDFScale,99,upv(2),dnv(2),usea(2),dsea(2),str(2),chm(2),bot(2),glu(2))
            sbar(2) = str(2)
            cbar(2) = chm(2)
            bbar(2) = bot(2)
        else
            upv(2) = 0d0
            dnv(2) = 0d0
            usea(2)= 0d0
            dsea(2)= 0d0
            str(2) = 0d0
            sbar(2)= 0d0
            chm(2) = 0d0
            cbar(2)= 0d0
            bot(2) = 0d0
            bbar(2)= 0d0
            glu(2) = 0d0
        endif
ENDIF




IF( COLLIDER.EQ.1 ) THEN
!       PROTON CONTENT
        pdf(Up_,1)   = (upv(1) + usea(1))  * swPDF_u / x1
        pdf(AUp_,1)  = usea(1)             * swPDF_u / x1
        pdf(Dn_,1)   = (dnv(1) + dsea(1))  * swPDF_d / x1
        pdf(ADn_,1)  = dsea(1)             * swPDF_d / x1
        pdf(Chm_,1)  = chm(1)              * swPDF_c / x1
        pdf(AChm_,1) = cbar(1)             * swPDF_c / x1
        pdf(Str_,1)  = str(1)              * swPDF_s / x1
        pdf(AStr_,1) = sbar(1)             * swPDF_s / x1
        pdf(Bot_,1)  = bot(1)              * swPDF_b / x1
        pdf(ABot_,1) = bbar(1)             * swPDF_b / x1
        pdf(0,1)     = glu(1)              * swPDF_g / x1

!       PROTON CONTENT
        pdf(Up_,2)   = (upv(2) + usea(2))  * swPDF_u / x2
        pdf(AUp_,2)  = usea(2)             * swPDF_u / x2
        pdf(Dn_,2)   = (dnv(2) + dsea(2))  * swPDF_d / x2
        pdf(ADn_,2)  = dsea(2)             * swPDF_d / x2
        pdf(Chm_,2)  = chm(2)              * swPDF_c / x2
        pdf(AChm_,2) = cbar(2)             * swPDF_c / x2
        pdf(Str_,2)  = str(2)              * swPDF_s / x2
        pdf(AStr_,2) = sbar(2)             * swPDF_s / x2
        pdf(Bot_,2)  = bot(2)              * swPDF_b / x2
        pdf(ABot_,2) = bbar(2)             * swPDF_b / x2
        pdf(0,2)     = glu(2)              * swPDF_g / x2

ELSEIF( COLLIDER.EQ.2 ) THEN
!       PROTON CONTENT
        pdf(Up_,1)   = (upv(1) + usea(1))  * swPDF_u / x1
        pdf(AUp_,1)  = usea(1)             * swPDF_u / x1
        pdf(Dn_,1)   = (dnv(1) + dsea(1))  * swPDF_d / x1
        pdf(ADn_,1)  = dsea(1)             * swPDF_d / x1
        pdf(Chm_,1)  = chm(1)              * swPDF_c / x1
        pdf(AChm_,1) = cbar(1)             * swPDF_c / x1
        pdf(Str_,1)  = str(1)              * swPDF_s / x1
        pdf(AStr_,1) = sbar(1)             * swPDF_s / x1
        pdf(Bot_,1)  = bot(1)              * swPDF_b / x1
        pdf(ABot_,1) = bbar(1)             * swPDF_b / x1
        pdf(0,1)     = glu(1)              * swPDF_g / x1

!       ANTI-PROTON CONTENT
        pdf(Up_,2)   = usea(2)             * swPDF_u / x2
        pdf(AUp_,2)  = (upv(2) + usea(2))  * swPDF_u / x2
        pdf(Dn_,2)   = dsea(2)             * swPDF_d / x2
        pdf(ADn_,2)  = (dnv(2) + dsea(2))  * swPDF_d / x2
        pdf(Chm_,2)  = chm(2)              * swPDF_c / x2
        pdf(AChm_,2) = cbar(2)             * swPDF_c / x2
        pdf(Str_,2)  = str(2)              * swPDF_s / x2
        pdf(AStr_,2) = sbar(2)             * swPDF_s / x2
        pdf(Bot_,2)  = bot(2)              * swPDF_b / x2
        pdf(ABot_,2) = bbar(2)             * swPDF_b / x2
        pdf(0,2)     = glu(2)              * swPDF_g / x2
ENDIF


RETURN
END SUBROUTINE














END MODULE
