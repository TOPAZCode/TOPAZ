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


integer,public :: it_sav

type(Histogram),allocatable :: Histo(:)
real(8) :: pT_jet_cut, pT_bjet_cut, pT_lep_cut, Rsep_jet, Rsep_LepJet, pT_miss_cut, eta_sepa_cut, MInv_jets_cut, eta_lep_cut, eta_jet_cut, eta_bjet_cut, HT_cut, pT_hardestjet_cut
real(8) :: pT_pho_cut,Rsep_Pj,Rsep_Pbj,Rsep_Plep,eta_pho_cut


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
    Rsep_jet    = 0.4d0
    pT_bjet_cut = 20d0*GeV
    eta_bjet_cut= 2.5d0
    pT_lep_cut  = 20d0*GeV
    pT_miss_cut = 25d0*GeV
    eta_lep_cut = 2.5d0

ELSEIF( ObsSet.EQ.3 ) THEN! set of observables for ttb production as signal process at LHC (di-lept. decay)
    Rsep_jet    = 0.4d0
    pT_bjet_cut = 25d0*GeV
    eta_bjet_cut= 2.5d0
    pT_lep_cut  = 25d0*GeV
    pT_miss_cut = 50d0*GeV
    eta_lep_cut = 2.5d0

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
    Rsep_jet    = 0.4d0
    pT_bjet_cut = 20d0*GeV
    eta_bjet_cut= 2.0d0
    pT_jet_cut  = 20d0*GeV
    eta_jet_cut = 2.5d0
    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 2.5d0
    pT_miss_cut = 20d0*GeV


ELSEIF( ObsSet.EQ.6 ) THEN! set of observables for ttb production with lept. top, hadr. Atop decay at LHC
    Rsep_jet    = 0.5d0
    pT_bjet_cut = 20d0*GeV
    eta_bjet_cut= 2.0d0
    pT_jet_cut  = 20d0*GeV
    eta_jet_cut = 2.0d0
    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 2.5d0
    pT_miss_cut = 20d0*GeV


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
    pT_jet_cut  = 25d0*GeV
    Rsep_jet    = 0.4d0

ELSEIF( ObsSet.EQ.12 ) THEN! set of observables for ttbjet production as signal process at Tevatron (hadr.Atop, lept.top decay)
    pT_jet_cut  = 20d0*GeV
    pT_bjet_cut = pT_jet_cut
    eta_jet_cut = 2d0
    eta_bjet_cut= eta_jet_cut
    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 2d0
    pT_miss_cut = 20d0*GeV
    HT_cut      = 220d0*GeV
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
    pT_miss_cut = 25d0*GeV
    Rsep_jet    = 0.4d0
    Rsep_LepJet = 0.4d0



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




ELSEIF( ObsSet.EQ.31 ) THEN! set of observables for TTbar + A0/BH production (stable)


ELSEIF( ObsSet.EQ.32 ) THEN! set of observables for TTbar + A0/BH production (di-lept. tops)
    Rsep_jet    = 0.4d0

    pT_bjet_cut = 30d0*GeV 
    eta_bjet_cut= 2.5d0 

    pT_lep_cut  = 20d0*GeV  
    eta_lep_cut = 2.5d0 
    pT_miss_cut = 25d0*GeV


ELSEIF( ObsSet.EQ.33 ) THEN! set of observables for TTbar + A0/BH production (semi-hadr. tops)
    Rsep_jet    = 0.4d0

    pT_bjet_cut = 30d0*GeV
    eta_bjet_cut= 2.5d0
    pT_jet_cut = 25d0*GeV 
    eta_jet_cut= 2.5d0     
 
    pT_lep_cut  = 20d0*GeV  
    eta_lep_cut = 2.5d0      
    pT_miss_cut = 25d0*GeV  


ELSEIF( ObsSet.EQ.41 ) THEN! set of observables for STSTbar + Chi production (stable)



ELSEIF( ObsSet.EQ.42 ) THEN! set of observables for STSTbar + Chi production (di-lept. tops)
    Rsep_jet    = 0.4d0

    pT_bjet_cut = 25d0*GeV
    eta_bjet_cut= 2.5d0
    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 2.5d0
    pT_miss_cut = 80d0*GeV


ELSEIF( ObsSet.EQ.43 ) THEN! set of observables for STSTbar + Chi production (semi-hadr. tops)
    Rsep_jet    = 0.4d0

    pT_bjet_cut = 30d0*GeV
    eta_bjet_cut= 2.5d0
    pT_jet_cut = 25d0*GeV
    eta_jet_cut= 2.5d0

    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 2.5d0 
    pT_miss_cut = 25d0*GeV 



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
          if(TopDecays.ne.3) call Error("TopDecays needs to be 3!")
          NumHistograms = 11
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif
          call Error("I think the binning routine is not set up properly for TopDK=3")


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
          NumHistograms = 15
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

          Histo(13)%Info   = "pT_ttbar"
          Histo(13)%NBins  = 40
          Histo(13)%BinSize= 15d0*GeV
          Histo(13)%LowVal = 0d0
          Histo(13)%SetScale= 100d0

          Histo(14)%Info   = "reconst. mTop(hadr)"
          Histo(14)%NBins  = 50
          Histo(14)%BinSize= 1d0*GeV
          Histo(14)%LowVal = 140d0*GeV
          Histo(14)%SetScale= 100d0

          Histo(15)%Info   = "reconst. mTop(lept)"
          Histo(15)%NBins  = 50
          Histo(15)%BinSize= 1d0*GeV
          Histo(15)%LowVal = 140d0*GeV
          Histo(15)%SetScale= 100d0


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
          NumHistograms = 15
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

          Histo(14)%Info   = "reconst. mTop(hadr)"
          Histo(14)%NBins  = 50
          Histo(14)%BinSize= 1d0*GeV
          Histo(14)%LowVal = 140d0*GeV
          Histo(14)%SetScale= 100d0

          Histo(15)%Info   = "reconst. mTop(lept)"
          Histo(15)%NBins  = 50
          Histo(15)%BinSize= 1d0*GeV
          Histo(15)%LowVal = 140d0*GeV
          Histo(15)%SetScale= 100d0




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
          if(TopDecays.ne.4  .and. TopDecays.ne.3 ) call Error("TopDecays needs to be 3 (for Qt=-4/3) or 4 (for Qt=Qup)")
          NumHistograms = 19
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

          Histo(15)%Info   = "m(jet+jet+pho)"
          Histo(15)%NBins  = 40
          Histo(15)%BinSize= 20d0*GeV
          Histo(15)%LowVal = 20d0*GeV
          Histo(15)%SetScale= 100d0

          Histo(16)%Info   = "mT(lep,pho;pTmiss)"
          Histo(16)%NBins  = 40
          Histo(16)%BinSize= 20d0*GeV
          Histo(16)%LowVal = 20d0*GeV
          Histo(16)%SetScale= 100d0

          Histo(17)%Info   = "m(bjet+jet+jet+pho)"
          Histo(17)%NBins  = 40
          Histo(17)%BinSize= 20d0*GeV
          Histo(17)%LowVal = 20d0*GeV
          Histo(17)%SetScale= 100d0

          Histo(18)%Info   = "mT(bjet+lep+pho;pTmiss)"
          Histo(18)%NBins  = 40
          Histo(18)%BinSize= 20d0*GeV
          Histo(18)%LowVal = 20d0*GeV
          Histo(18)%SetScale= 100d0

          Histo(19)%Info   = "phi(photon,lept)"
          Histo(19)%NBins  = 20
          Histo(19)%BinSize= 0.2d0
          Histo(19)%LowVal = 0d0
          Histo(19)%SetScale= 1d0


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



ELSEIF( ObsSet.EQ.31 ) THEN! set of observables for TTbar + A0 production
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          NumHistograms = 1
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

ELSEIF( ObsSet.EQ.32 ) THEN! set of observables for TTbar + A0/BH production
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.1  ) call Error("TopDecays needs to be 1!")
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


ELSEIF( ObsSet.EQ.33 ) THEN! set of observables for TTbar + A0/BH production
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.4  ) call Error("TopDecays needs to be 4!")
          NumHistograms = 9
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_LepP"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 40d0*GeV
          Histo(1)%LowVal = 00d0*GeV
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
          Histo(8)%BinSize= 40d0*GeV
          Histo(8)%LowVal = 0d0*GeV
          Histo(8)%SetScale= 100d0

          Histo(9)%Info   = "eta_Top"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 0.25d0
          Histo(9)%LowVal =-5.0d0
          Histo(9)%SetScale= 1d0


ELSEIF( ObsSet.EQ.41 ) THEN! set of observables for STSTbar + Chi production (stable tops)
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          NumHistograms = 1
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 20d0*GeV
          Histo(1)%LowVal = 20d0*GeV
          Histo(1)%SetScale= 100d0



ELSEIF( ObsSet.EQ.42 ) THEN! set of observables for STSTbar + Chi production (di-lept. tops)

          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.1  ) call Error("TopDecays needs to be 1!")
          NumHistograms = 5
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


ELSEIF( ObsSet.EQ.43 ) THEN! set of observables for TTbar + A0/BH production
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.4  ) call Error("TopDecays needs to be 4!")
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
          Histo(8)%BinSize= 20d0*GeV
          Histo(8)%LowVal = 0d0*GeV
          Histo(8)%SetScale= 100d0


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
      call genps(2,m_HTop,xRndPS(1:2),(/Mass,m_Top/),MomDK(1:4,1:2),PSWgt2)! Htop decay
!     boost all guys to the top frame:
      call boost(MomDK(1:4,1),HTopMom(1:4),m_HTop)
      call boost(MomDK(1:4,2),HTopMom(1:4),m_HTop)
      PSWgt = PSWgt2*PiWgt2

    else! extra gluon emission

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
real(8) :: SingDepth,soft,coll
logical,save :: flip=.true.

    if( .not.GluonRad ) then!  no extra gluon radiation
!     MomDK(1:4,i): i= 1:bottom, 2:lepton, 3:neutrino
      call genps(2,m_STop,xRndPS(1:2),(/m_Chi,m_Top/),MomDK(1:4,1:2),PSWgt2)! Stop decay
!     boost all guys to the top frame:
      call boost(MomDK(1:4,1),STopMom(1:4),m_STop)
      call boost(MomDK(1:4,2),STopMom(1:4),m_STop)
      PSWgt = PSWgt2*PiWgt2

    else! extra gluon emission

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
!         Pcol1= 4 -1
!         Pcol2= 4 -1
!         SingDepth = 1e-10
!         Steps = 20
!         call gensing(3,m_Top,(/0d0,0d0,m_W/),MomDK(1:4,1:3),Pcol1,Pcol2,SingDepth,Steps)
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
!         Pcol1= 3 -1
!         Pcol2= 5 -1
!         SingDepth = 1e-10
!         Steps = 20
!         call gensing(3,m_W,(/0d0,0d0,0d0/),MomDK(1:4,2:4),Pcol1,Pcol2,SingDepth,Steps)
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





SUBROUTINE EvalPhasespace_2to2Stops(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5
real(8) :: Mom(1:4,1:4),MomW(1:4),xRndPS(1:2)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)


!  generate PS: massless + massless --> massive(anti-top) + massive(top)
   call genps(2,Ehat,xRndPS(1:2),(/m_Stop,m_Stop/),Mom(1:4,3:4),PSWgt)
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
!      Steps = 15
!      PSWgt = 1d0
!      call gensing(3,EHat,(/0d0,m_Top,m_Top/),Mom(1:4,3:5),Pcol1,Pcol2,SingDepth,Steps)

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











SUBROUTINE Kinematics_TTBARJET(NPlus1PS,MomExt,MomDK,applyPSCut,NBin)
use ModMisc
use ModParameters
implicit none
integer :: NumHadr,NPlus1PS
real(8) :: MomExt(1:4,1:5+NPlus1PS),MomDK(1:4,1:6),MomJet(1:4,1:8),zeros(1:10)
real(8) :: MomHadr(1:4,1:8),MomLept(1:4,1:4)
real(8) :: MomBoost(1:4),MomW(1:4),MomTops(1:4,1:2)
logical :: applyPSCut
integer :: NBin(:),PartList(1:8),JetList(1:8),NJet,NObsJet,n,NObsJet_Tree,nWJets
real(8) :: pT_lepM,pT_lepP,pT_miss,pT_ATop,pT_Top,HT,m_lb,R_lb,m_bb,m_bj
real(8) :: eta_ATop,eta_Top,eta_lepM,eta_lepP,eta_miss
real(8) :: pT_jet(1:8),eta_jet(1:8),eta_sepa,eta_Zeppi,s34,s35,s36,s45,s46,s56,mTopHadr,mTopLept
real(8) :: R_bb,MinvJets,MinvLept,phi_Lept,pT_lept,ET_lept,ET_miss,mT,pT_x,pT_y



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

elseif( TopDecays.eq.5 ) then  ! hadr.  top, lept. top decay with J/Psi fragmentation
  print *, "do something"
  stop
elseif( TopDecays.eq.6 ) then  ! hadr. Atop, lept. top decay with J/Psi fragmentation
  print *, "do something"
  stop
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

!       reconstruct mTop
        if( get_MInv(MomJet(1:4,1)+MomLept(1:4,3)) .lt. get_MInv(MomJet(1:4,2)+MomLept(1:4,3)) ) then! pair bjet and lepton wrt. inv.mass
            mTopLept = get_MInv( MomJet(1:4,1)+MomLept(1:4,3)+MomLept(1:4,4) )
            mTopHadr = get_MInv( MomJet(1:4,2)+MomW(1:4) )
        else
            mTopLept = get_MInv( MomJet(1:4,2)+MomLept(1:4,3)+MomLept(1:4,4) )
            mTopHadr = get_MInv( MomJet(1:4,1)+MomW(1:4) )
        endif
    else
        pT_Top   = -1d0
        mTopHadr = -1d0
        mTopLept = -1d0
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
    NBin(14)= WhichBin(14,mTopHadr)
    NBin(15)= WhichBin(15,mTopLept)


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

print *, "STOP! proper jet selectin needs to be implemented first! see ObsSet=12,13"; stop


    NObsJet_Tree = 5
!   request at least two b-jets and three non-b-jet
    if( .not.(NJet.ge.NObsJet_Tree .and. any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif

    pT_lepP  = get_PT(MomLept(1:4,3))
    eta_lepP = get_eta(MomLept(1:4,3))

    pT_miss  = get_PT(MomLept(1:4,4))

    phi_Lept = dabs( Get_PHI(MomLept(1:4,3)) - Get_PHI(MomJet(1:4,1)) )
    if( phi_Lept.gt.Pi ) phi_Lept=2d0*Pi-phi_Lept

    m_lb = get_MInv(MomLept(1:4,3)+MomJet(1:4,1))
    R_lb = get_R(MomLept(1:4,3),MomJet(1:4,1))


    do n=1,NJet! first two jets are always b-jets
      pT_jet(n)  = get_pT( MomJet(1:4,n))
      eta_jet(n) = get_eta(MomJet(1:4,n))
    enddo

    HT = pT_lepP + pT_miss

! check cuts
    if( pT_jet(1).lt.pT_bjet_cut .or. pT_jet(2).lt.pT_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( abs(eta_jet(1)).gt.eta_bjet_cut .or. abs(eta_jet(2)).gt.eta_bjet_cut) then
        applyPSCut = .true.
        RETURN
    endif
    if( get_R(MomJet(1:4,1),MomLept(1:4,3)).lt.Rsep_LepJet .or. get_R(MomJet(1:4,2),MomLept(1:4,3)).lt.Rsep_LepJet ) then
        applyPSCut = .true.
        RETURN
    endif

    HT = HT + pT_jet(1) + pT_jet(2)

    NObsJet = 0
    do n=3,NJet
        if( pT_jet(n).gt.pT_jet_cut ) NObsJet = NObsJet +1
    enddo
    if( NObsJet.lt.NObsJet_Tree-2 ) then
        applyPSCut = .true.
        RETURN
    endif
    do n=3,NJet
        if( pT_jet(n).gt.pT_jet_cut .and. ( abs(eta_jet(n)).gt.eta_jet_cut .or. get_R(MomJet(1:4,n),MomLept(1:4,3)).lt.Rsep_LepJet ) ) then
            applyPSCut = .true.
            RETURN
        endif
        if( pT_jet(n).gt.pT_jet_cut .and. abs(eta_jet(n)).lt.eta_jet_cut ) HT = HT + pT_jet(n)
    enddo

!   construct the t+bar system
    MomTops(1:4,1) = MomJet(1:4,1)+MomJet(1:4,2)+MomLept(1:4,3)+MomLept(1:4,4)
!   find two non-b jets that are closest to MW mass
    if( NObsJet.eq.4 ) then
        s34= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,4))-M_W )
        s35= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,5))-M_W )
        s36= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,6))-M_W )
        s45= dabs( get_MInv(MomJet(1:4,4)+MomJet(1:4,5))-M_W )
        s46= dabs( get_MInv(MomJet(1:4,4)+MomJet(1:4,6))-M_W )
        s56= dabs( get_MInv(MomJet(1:4,5)+MomJet(1:4,6))-M_W )
    elseif( NObsJet.eq.3 ) then
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
    if( dmin1(s34,s35,s45,s36,s46,s56).lt.20d0*GeV ) then!   require a 20GeV window for M_W
        pT_Top = get_pt(MomTops(1:4,1))

!       reconstruct mTop
        if( get_MInv(MomJet(1:4,1)+MomLept(1:4,3)) .lt. get_MInv(MomJet(1:4,2)+MomLept(1:4,3)) ) then! pair bjet and lepton wrt. inv.mass
            mTopLept = get_MInv( MomJet(1:4,1)+MomLept(1:4,3)+MomLept(1:4,4) )
            mTopHadr = get_MInv( MomJet(1:4,2)+MomW(1:4) )
        else
            mTopLept = get_MInv( MomJet(1:4,2)+MomLept(1:4,3)+MomLept(1:4,4) )
            mTopHadr = get_MInv( MomJet(1:4,1)+MomW(1:4) )
        endif
    else
        pT_Top   = -1d0
        mTopHadr = -1d0
        mTopLept = -1d0
    endif


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


! binning
    call pT_order(NJet,MomJet(1:4,1:NJet))! pT ordering of jet momenta for b- AND non-b-jets
    pT_jet(5)  = get_pT(MomJet(1:4,5))
    eta_jet(5) = get_eta(MomJet(1:4,5))

    NBin(1) = WhichBin(1,pT_lepP)
    NBin(2) = WhichBin(2,eta_lepP)
    NBin(3) = WhichBin(3,get_pT(MomJet(1:4,1)))
    NBin(4) = WhichBin(4,get_eta(MomJet(1:4,1)))
    NBin(5) = WhichBin(5,get_pT(MomJet(1:4,5)))
    NBin(6) = WhichBin(6,get_eta(MomJet(1:4,5)))
    NBin(7) = WhichBin(7,pT_miss)
    NBin(8) = WhichBin(8,HT)
    NBin(9) = WhichBin(9,m_lb)
    NBin(10)= WhichBin(10,phi_Lept)
    NBin(11)= WhichBin(11,R_lb)
    NBin(12)= WhichBin(12,eta_lepP)
    NBin(13)= WhichBin(13,pT_Top)
    NBin(14)= WhichBin(14,mTopHadr)
    NBin(15)= WhichBin(15,mTopLept)



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
real(8) :: eta_ATop,eta_Top,eta_lepM,eta_lepP,m_lb,mT_lpmiss,mT_blpmiss,m_jj,mTblP,m_jjb,mT_lp
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
            HT = HT + get_pT(MomJet(1:4,k))
        endif
    enddo

    NObsJet_Tree = 3! request two b-jets and at least one light jet
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

    if( get_pT(MomJet(1:4,1)).gt.get_pT(MomJet(1:4,2)) ) then
        ET_bjet  = get_ET(MomJet(1:4,1))
        m_lb = get_Minv(Mom(1:4,L)+MomJet(1:4,1))
        Rphobjet = get_R(Mom(1:4,pho),MomJet(1:4,1))
    else
        ET_bjet  = get_ET(MomJet(1:4,2))
        m_lb = get_Minv(Mom(1:4,L)+MomJet(1:4,1))
        Rphobjet = get_R(Mom(1:4,pho),MomJet(1:4,1))
    endif
    pT_lepP  = get_PT(Mom(1:4,L))
    eta_lepP = get_ETA(Mom(1:4,L))

!    MomObs(1:4) = MomObs(1:4) + Mom(1:4,pho) + Mom(1:4,L)
!    MomMiss(1:4) = Mom(1:4,inLeft) + Mom(1:4,inRight) - MomObs(1:4)
    MomMiss(1:4) = Mom(1:4,N)
    ET_miss  = get_ET(MomMiss(1:4))
    HT = HT + pT_lepP + ET_miss + pT_Pho



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


! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,eta_ATop)
    NBin(3) = WhichBin(3,pT_Top)
    NBin(4) = WhichBin(4,eta_Top)
    NBin(5) = WhichBin(5,eta_lepP)
    NBin(6) = WhichBin(6,eta_Pho)
    NBin(7) = WhichBin(7,pT_Pho)
    NBin(8) = WhichBin(8,eta_Pho)
    NBin(9) = WhichBin(9,pT_lepP)
    NBin(10)= WhichBin(10,eta_lepP)
    NBin(11) = WhichBin(11,ET_miss)
    NBin(12) = WhichBin(12,HT)
    NBin(13) = WhichBin(13,Rphobjet)
    NBin(14) = WhichBin(14,m_lb)
    NBin(15) = WhichBin(15,mT_lpmiss)
    NBin(16) = WhichBin(16,Phi_LP)
    NBin(17) = WhichBin(17,ET_bjet)




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
real(8) :: ptmax, xij, MomBBar(1:4),diff,cosLeBb,Ebjets,r_sc
integer :: i1,i,j,n,k

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
    NBin(1) = WhichBin(1,pT_lepM)
    NBin(2) = WhichBin(2,eta_lepM)
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

    if( HT.lt.HT_cut ) then
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
integer :: Htbar,Ht,X0bar,X0,tbar,t,realp,bbar,lepM,nubar,b,lepP,nu,qdn,qbup,qbdn,qup,Lep,Neu




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

IF( ObsSet.eq.31 .or. ObsSet.eq.32 .or. ObsSet.eq.33 ) then
   IF(XTOPDECAYS.EQ.1) m_X0 = m_BH
   IF(XTOPDECAYS.EQ.2) m_X0 = m_A0
elseif( ObsSet.eq.41 .or. ObsSet.eq.42 .or. ObsSet.eq.43 ) then
   m_X0 = m_Chi
endif

!DEC$ IF(_CheckMomenta .EQ.1)
IF(XTOPDECAYS.NE.0) THEN
   zeros(1:4) = Mom(1:4,1)+Mom(1:4,2) - (Mom(1:4,X0bar)+Mom(1:4,X0)+Mom(1:4,bbar)+Mom(1:4,lepM)+Mom(1:4,nubar)+Mom(1:4,b)+Mom(1:4,lepP)+Mom(1:4,nu))
   if( NPlus1PS ) zeros(1:4) = zeros(1:4) + Mom(1:4,realp)
   if( any(abs(zeros(1:4)/m_HTop).gt.1d-8) ) then
      print *, "ERROR: energy-momentum violation 1 in SUBROUTINE Kinematics_HTTBARDM: ",zeros(1:4)
      print *, "momentum dump:"
      print *, Mom(:,:)
   endif

   zeros(1:4) = Mom(1:4,Htbar)-Mom(1:4,X0bar)-Mom(1:4,tbar)
   zeros(5:8) = Mom(1:4,Ht)-Mom(1:4,X0)-Mom(1:4,t)
   if( any(abs(zeros(1:8)/m_HTop).gt.1d-8) ) then
      print *, "ERROR: energy-momentum violation 2 in SUBROUTINE Kinematics_HTTBARDM: ",zeros(1:8)
      print *, "momentum dump:"
      print *, Mom(:,:)
      pause
   endif

   zeros(1:4) = Mom(1:4,tbar)-Mom(1:4,bbar)-Mom(1:4,LepM)-Mom(1:4,nubar)
   zeros(5:8) = Mom(1:4,t)-Mom(1:4,b)-Mom(1:4,lepP)-Mom(1:4,nu)
   if( any(abs(zeros(1:8)/m_HTop).gt.1d-8) ) then
      print *, "ERROR: energy-momentum violation 3 in SUBROUTINE Kinematics_HTTBARDM: ",zeros(1:8)
      print *, "momentum dump:"
      print *, Mom(:,:)
      pause
   endif

   zeros(:) = 0d0
IF( ObsSet.eq.31 .or. ObsSet.eq.32 .or. ObsSet.eq.33 ) then
   zeros(1) = (Mom(1:4,HTbar).dot.Mom(1:4,HTbar))  - m_HTop**2
   zeros(2) = (Mom(1:4,HT).dot.Mom(1:4,HT)) - m_HTop**2
elseif( ObsSet.eq.41 .or. ObsSet.eq.42 .or. ObsSet.eq.43 ) then
   zeros(1) = (Mom(1:4,HTbar).dot.Mom(1:4,HTbar))  - m_STop**2
   zeros(2) = (Mom(1:4,HT).dot.Mom(1:4,HT)) - m_STop**2
endif
   zeros(3) = (Mom(1:4,X0bar).dot.Mom(1:4,X0bar)) - m_X0**2
   zeros(4) = (Mom(1:4,X0).dot.Mom(1:4,X0)) - m_X0**2
   zeros(5) = (Mom(1:4,tbar).dot.Mom(1:4,tbar)) - m_Top**2
   zeros(6) = (Mom(1:4,t).dot.Mom(1:4,t)) - m_Top**2
   zeros(7) = (Mom(1:4,bbar).dot.Mom(1:4,bbar))
   zeros(8) = (Mom(1:4,b).dot.Mom(1:4,b))
   zeros(9) = (Mom(1:4,lepm).dot.Mom(1:4,lepm))
   zeros(10) = (Mom(1:4,lepp).dot.Mom(1:4,lepp))
   zeros(11) = (Mom(1:4,nubar).dot.Mom(1:4,nubar))
   zeros(12) = (Mom(1:4,nu).dot.Mom(1:4,nu))
   if(NPlus1PS) zeros(13) = (Mom(1:4,realp).dot.Mom(1:4,realp))

   if( any(abs(zeros(1:13)/1d0).gt.1d-5) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_HTTBARDM: ",zeros(1:13)
      print *, "momentum dump:"
      print *, Mom(:,:)
      pause
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




!-------------------------------------------------------
elseif( ObsSet.eq.32 .or. ObsSet.eq.42) then! set of observables for TTbar -> ttbar + ETmiss  in di-lept. top decays

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

    pT_miss = get_PT( Mom(1:4,nu)+Mom(1:4,nubar)+Mom(1:4,X0)+Mom(1:4,X0bar) )

    phi_ll = dabs( Get_PHI(Mom(1:4,lepM)) - Get_PHI(Mom(1:4,lepP)) )
    if( phi_ll.gt.Pi ) phi_ll=2d0*Pi-phi_ll

    m_ll = get_MInv(Mom(1:4,lepM)+Mom(1:4,lepP))


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

    if( pT_miss.lt.pT_miss_cut ) then
       applyPSCut = .true.
        RETURN
    endif


! binning
    NBin(1)= WhichBin(1,pT_lepP)
    NBin(2)= WhichBin(2,eta_lepP)
    NBin(3)= WhichBin(3,pT_miss)
    NBin(4)= WhichBin(4,phi_ll)
    NBin(5)= WhichBin(5,m_ll)





!-------------------------------------------------------
elseif( ObsSet.eq.33 .or. ObsSet.eq.43) then! set of observables for TTbar -> ttbar + ETmiss  in semi-hadr. top decays

! request at least two b-jets
    NObsJet_Tree = 4
    if( .not.(NJet.ge.NObsJet_Tree .and. any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif

    pT_jet(1) = get_PT( MomJet(1:4,1))
    pT_jet(2) = get_PT( MomJet(1:4,2))
    pT_jet(3) = get_PT( MomJet(1:4,3))
    pT_jet(4) = get_PT( MomJet(1:4,4))
    eta_jet(1)= get_ETA(MomJet(1:4,1))
    eta_jet(2)= get_ETA(MomJet(1:4,2))
    eta_jet(3)= get_ETA(MomJet(1:4,3))
    eta_jet(4)= get_ETA(MomJet(1:4,4))

    pT_lepP  = get_PT(Mom(1:4,lep))
    eta_lepP = get_ETA(Mom(1:4,lep))

    pT_miss = get_PT( Mom(1:4,nu)+Mom(1:4,X0)+Mom(1:4,X0bar) )

    phi_ll = dabs( Get_PHI(Mom(1:4,lepM)) - Get_PHI(Mom(1:4,lepP)) )
    if( phi_ll.gt.Pi ) phi_ll=2d0*Pi-phi_ll

    m_ll = get_MInv(Mom(1:4,lepM)+Mom(1:4,lepP))

    pT_Top = get_PT(Mom(1:4,t))
    eta_Top = get_ETA(Mom(1:4,t))

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

    if( pT_miss.lt.pT_miss_cut ) then
       applyPSCut = .true.
        RETURN
    endif

   MomAux(2) = Mom(2,lep)+Mom(2,neu)
   MomAux(3) = Mom(3,lep)+Mom(3,neu)
!   MTW = dsqrt( (pT_lepP+get_PT(Mom(1:4,neu)))**2  - (MomAux(2)**2+MomAux(3)**2) )
   MTW = Get_MT(Mom(1:4,lep),Mom(1:4,nu)+Mom(1:4,X0)+Mom(1:4,X0bar))
   MTeff = pt_jet(1)+pt_jet(2)+pt_jet(3)+pt_jet(4)+pT_lepP + pT_miss


!    if( MTW.lt.220*GeV ) then
!       applyPSCut = .true.
!        RETURN
!    endif


!    if( MTeff.lt.1000*GeV ) then
!       applyPSCut = .true.
!        RETURN
!    endif


! binning
    NBin(1)= WhichBin(1,pT_lepP)
    NBin(2)= WhichBin(2,eta_lepP)
    NBin(3)= WhichBin(3,pT_miss)
    NBin(4)= WhichBin(4,phi_ll)
    NBin(5)= WhichBin(5,m_ll)
    NBin(6)= WhichBin(6,MTW)
    NBin(7)= WhichBin(7,MTeff)
    NBin(8)= WhichBin(8,pT_Top)
    NBin(9)= WhichBin(9,eta_Top)




!-------------------------------------------------------
else
  print *, "ObsSet not implemented",ObsSet
  stop
endif


return
END SUBROUTINE













RECURSIVE SUBROUTINE JetAlgo_kt(Rsep_jet,PartonList,MomParton,NJet,JetList,MomJet)  ! initial call must have NJet=0
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
complex(8) e(1:NumExtParticles,1:4), tmp(1:4)
integer :: NPart

   do NPart=1,NumExtParticles
      if(abs(ExtParticle(NPart)%Helicity).ne.1) cycle ! this prevents overriding of top decay spinors,  (ExtParticle(NPart)%Helicity=0)
      if( ExtParticle(NPart)%PartType .eq. Glu_ ) then
         call pol_mless(ExtParticle(NPart)%Mom(1:4),ExtParticle(NPart)%Helicity,ExtParticle(NPart)%Pol(1:4))
         ExtParticle(NPart)%Pol(5:16) = (0d0,0d0)
      elseif( IsAQuark(ExtParticle(NPart)%PartType) .and. ExtParticle(NPart)%PartType.gt.0 ) then
         call ubarSpi(ExtParticle(NPart)%Mom(1:4),ExtParticle(NPart)%Mass,ExtParticle(NPart)%Helicity,ExtParticle(NPart)%Pol(1:4))
         ExtParticle(NPart)%Pol(5:16) = (0d0,0d0)
      elseif( IsAQuark(ExtParticle(NPart)%PartType) .and. ExtParticle(NPart)%PartType.lt.0 ) then
         call vSpi(ExtParticle(NPart)%Mom(1:4),ExtParticle(NPart)%Mass,ExtParticle(NPart)%Helicity,ExtParticle(NPart)%Pol(1:4))
         ExtParticle(NPart)%Pol(5:16) = (0d0,0d0)
      endif
   enddo

!       check gauge invariance
        ExtParticle(3)%Pol(1:4) = ExtParticle(3)%Mom(1:4);       print *, "gauge invariance check"
!         ExtParticle(4)%Pol(1:4) = ExtParticle(4)%Mom(1:4);       print *, "gauge invariance check"

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
!             call mrstlo(x1,PDFScale,1,upv(1),dnv(1),usea(1),dsea(1),str(1),chm(1),bot(1),glu(1))
!             call mrs96(x1,PDFScale,2,upv(1),dnv(1),usea(1),dsea(1),str(1),chm(1),bot(1),glu(1))
            call GetAllPDFs("mstw2008lo",0,x1,PDFScale,upv(1),dnv(1),usea(1),dsea(1),str(1),sbar(1),chm(1),cbar(1),bot(1),bbar(1),glu(1),phot)
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
!             call mrstlo(x2,PDFScale,1,upv(2),dnv(2),usea(2),dsea(2),str(2),chm(2),bot(2),glu(2))
!             call mrs96(x2,PDFScale,2,upv(2),dnv(2),usea(2),dsea(2),str(2),chm(2),bot(2),glu(2))
            call GetAllPDFs("mstw2008lo",0,x2,PDFScale,upv(2),dnv(2),usea(2),dsea(2),str(2),sbar(2),chm(2),cbar(2),bot(2),bbar(2),glu(2),phot)
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
            call GetAllPDFs("mstw2008nlo",0,x1,PDFScale,upv(1),dnv(1),usea(1),dsea(1),str(1),sbar(1),chm(1),cbar(1),bot(1),bbar(1),glu(1),phot)
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
            call GetAllPDFs("mstw2008nlo",0,x2,PDFScale,upv(2),dnv(2),usea(2),dsea(2),str(2),sbar(2),chm(2),cbar(2),bot(2),bbar(2),glu(2),phot)
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
