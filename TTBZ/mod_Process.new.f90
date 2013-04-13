MODULE ModProcess
implicit none


!! important convention for 1-loop amplitudes:
!! external quarks need to be Str/AStr
!! massive quarks (i.e. tops) in closed fermion loops are Bot
!! massless quarks in closed fermion loops are Chm


type :: Particle
   integer  :: PartType
   integer  :: ExtRef
   integer  :: Helicity
   real(8)  :: Mass
   real(8)  :: Mass2
   complex(8) :: Mom(1:8)
   complex(8) :: Pol(1:16)
end type

type :: Propagator
   integer  :: PartType
   integer  :: ExtRef
   real(8)  :: Mass
   real(8)  :: Mass2
   complex(8) :: Mom(1:8)
end type

type :: PtrToParticle
   integer,pointer  :: PartType
   integer,pointer  :: ExtRef
   real(8),pointer  :: Mass
   real(8),pointer  :: Mass2
   integer ,pointer :: Helicity
   complex(8),pointer :: Mom(:)
   complex(8),pointer :: Pol(:)
end type

type :: TreeProcess
   integer :: NumPart
   integer :: NumQua
   integer,allocatable :: NumGlu(:)
   integer,allocatable :: PartRef(:)
   integer,allocatable :: PartType(:)
   type(PtrToParticle),allocatable :: Quarks(:)
   type(PtrToParticle),allocatable :: Gluons(:)
   type(PtrToParticle) :: Boson
end type


type :: AMatch
  integer :: NumMatch
  integer, allocatable :: MatchHiCuts(:)
  integer, allocatable :: FirstHiProp(:)
  integer, allocatable :: MissHiProp(:,:)
end type

type :: UCutMatch
  type(AMatch) :: Subt(2:5)
end type


type :: UnitarityCut
   integer              :: CutType              ! 5,...,1: pent.,...,sing. cut
   integer              :: NumCuts              ! number of cuts
   integer, allocatable :: CutProp(:,:)         ! CutProp(cut number, 1..CutType): cutted prop. in highest-level diagr.
   complex(8), allocatable :: Coeff(:,:)        ! coefficients of these cuts
   type(TreeProcess), allocatable :: TreeProcess(:,:)  ! process for a cut at a vertex, first and last number are the prop.ID, the rest corresp. to the ExtParticle ID
   real(8), allocatable :: KMom(:,:,:)
   complex(8), allocatable :: NMom(:,:,:)
   type(UCutMatch),allocatable :: Match(:)      ! for matching a cut with higher cuts
end type


type :: PrimitiveAmplitude
   integer :: NPoint                            ! highest level of loop diagram
   integer :: AmpType                             ! J=0/1/2: subset of closed scalar loop/all other/fermion loop
   integer :: FermLoopPart                      ! fermion type in closed loop
   integer, allocatable :: ExtLine(:)        ! sequence of ext. particles
   type(Propagator), allocatable :: IntPart(:)   ! internal particle, i.e. propagator particle
   integer :: FermionLines                      ! number of quark lines
   integer :: FermLine1In                       ! vertex number of beginning of 1st quark line
   integer :: FermLine1Out                      ! vertex number of termination of 1st quark line
   integer :: FermLine2In                       ! vertex number of beginning of 2nd quark line
   integer :: FermLine2Out                      ! vertex number of termination of 2nd quark line
   type(UnitarityCut) :: UCuts(1:5)             ! unitarity cuts for this prim.ampl.
   complex(8) :: Result(-2:1)
!    integer :: NumInsertions1                    ! number of insertions of colorless particles in quark line 1 (corresp. to NCut)
!    integer :: NumInsertions2                    ! number of insertions of colorless particles in quark line 2
!    type(TreeProcess) :: TreeProc                ! tree amplitude
end type

type :: BornAmplitude
   integer, allocatable :: ExtLine(:)        ! sequence of ext. particles
   type(TreeProcess) :: TreeProc             ! tree amplitude
   complex(8) :: Result
end type


type(Particle),allocatable, target :: ExtParticle(:)
integer,allocatable :: Helicities(:,:)
type(PrimitiveAmplitude),allocatable, target :: PrimAmps(:)
type(BornAmplitude),allocatable, target :: BornAmps(:)

integer, public :: NumExtParticles,NumHelicities,NumPrimAmps,NumBornAmps

integer h1,h2,h3,h4,h5,h6,ih

contains




SUBROUTINE InitProcess()
use ModParameters
use ModMisc
implicit none
real(8) :: QuarkCrossing=-1d0, SpinAvg=1d0/4d0, QuarkColAvg=1d0/3d0, GluonColAvg=1d0/8d0, SymmFact=1d0/2d0
include "vegas_common.f"


NDim = 0
IF( abs(TOPDECAYS).GE.1 ) THEN
  NDim = NDim + 4+4
ENDIF

IF( TOPDECAYS.EQ.5 .OR. TOPDECAYS.EQ.6 ) THEN
  NDim = NDim + 1
  IF(CORRECTION.EQ.4) THEN
      NDim = NDim + 1
  ENDIF
ENDIF

IF(HelSampling) THEN
  NDim = NDim + 1
ENDIF


IF( PROCESS.EQ.0 ) THEN !   3_Glu  + 4_Glu  --> 5_Glu  + 1_Glu  + 2_Glu + 6_Glu
  IF( CORRECTION.EQ.0 ) THEN
      NumExtParticles = 6
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/5,6,-1,-2,3,4/)
      MasterProcess=0
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 8
      VegasNc0_default = 0
      VegasNc1_default = 2
  ELSEIF( CORRECTION.EQ.1 ) THEN
      NumExtParticles = 7
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/6,7,-1,-2,3,4,5/)
      MasterProcess=0
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 11
      VegasNc0_default = 0
      VegasNc1_default = 10
  ELSE
      call Error("Correction to this process is not available")
  ENDIF



ELSEIF( PROCESS.EQ.1 ) THEN !   3_Glu  + 4_Glu  --> 1_ATop + 2_Top
  IF( CORRECTION.EQ.0 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=1
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSEIF( CORRECTION.EQ.1 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=1
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 50000
      VegasNc1_default = 50000
  ELSEIF( CORRECTION.EQ.4 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=1
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 1000000
      VegasNc1_default = 1000000
  ELSEIF( CORRECTION.EQ.5 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=1
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 3+3  ! additional gluons in the top decay
!                                      NDIM = NDIM + 8  ! ADDITIONAL GLUONS IN THE TOP DECAY WITH SEPARATE YRND!!!!!!!
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF



ELSEIF( PROCESS.EQ.2 ) THEN !   3_Str  + 4_AStr --> 1_ATop + 2_Top
  IF( CORRECTION.EQ.0 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSEIF( CORRECTION.EQ.1 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 50000
      VegasNc1_default = 50000
  ELSEIF( CORRECTION.EQ.4 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 1000000
      VegasNc1_default = 1000000
  ELSEIF( CORRECTION.EQ.5 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 3+3  ! additional gluons in the top decay
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF



ELSEIF( PROCESS.EQ.3 ) THEN !   3_Str  + 5_Glu  --> 4_Str  + 1_ATop + 2_Top
  IF( CORRECTION.EQ.0 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,3,-2/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg*GluonColAvg
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000
      VegasNc1_default = 10000
  ELSEIF( CORRECTION.EQ.1 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,3,-2/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg*GluonColAvg
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 0
      VegasNc1_default = 10000000
  ELSEIF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,3,-2/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg * GluonColAvg
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg * GluonColAvg
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSEIF( CORRECTION.EQ.4 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,3,-2/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg*GluonColAvg
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 1000000
      VegasNc1_default = 1000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF



ELSEIF( PROCESS.EQ.4 ) THEN !   4_AStr + 5_Glu  --> 3_AStr + 1_ATop + 2_Top
  IF( CORRECTION.EQ.0 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,3,-2,-1/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg*GluonColAvg
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000
      VegasNc1_default = 10000
  ELSEIF( CORRECTION.EQ.1 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,3,-2,-1/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg*GluonColAvg
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 0
      VegasNc1_default = 10000000
  ELSEIF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,3,-1,-2/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg * GluonColAvg
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg * GluonColAvg
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSEIF( CORRECTION.EQ.4 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,3,-2,-1/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg*GluonColAvg
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 1000000
      VegasNc1_default = 1000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF



ELSEIF( PROCESS.EQ.5 ) THEN !   3_Glu  + 4_Glu  --> 5_Glu  + 1_ATop + 2_Top
  IF( CORRECTION.EQ.0 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=3
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000
      VegasNc1_default = 10000
  ELSEIF( CORRECTION.EQ.1 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=3
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 0
      VegasNc1_default = 10
  ELSEIF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=3
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 20000000
      VegasNc1_default = 20000000
  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=1
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSEIF( CORRECTION.EQ.4 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=3
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 100000
      VegasNc1_default = 100000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF



ELSEIF( PROCESS.EQ.6 ) THEN !   3_Str  + 4_AStr --> 5_Glu  + 1_ATop + 2_Top
  IF( CORRECTION.EQ.0 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000
      VegasNc1_default = 10000

!       NumExtParticles = 4!!!   this is for run of ttbar + jet from decay
!       allocate(Crossing(1:NumExtParticles))
!       allocate(ExtParticle(1:NumExtParticles))
!       Crossing(:) = (/3,4,-1,-2/)
!       MasterProcess=2
!       AvgFactor = SpinAvg * QuarkColAvg**2
!       NDim = NDim + 2    ! t tbar PS integration
!       NDim = NDim + 2    ! shat integration
!       ndim =ndim + 3 ! jet
!       VegasNc0_default = 10000000
!       VegasNc1_default = 10000000


  ELSEIF( CORRECTION.EQ.1 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 0
      VegasNc1_default = 10

  ELSEIF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 20000000
      VegasNc1_default = 20000000
  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSEIF( CORRECTION.EQ.4 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 100000
      VegasNc1_default = 100000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF



ELSEIF( PROCESS.EQ.9 ) THEN !   3_Glu  + 4_Glu  --> 5_Glu  + 6_Glu  + 1_ATop + 2_Top
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 6
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/5,6,-1,-2,3,4/)
      MasterProcess=5
      AvgFactor = SpinAvg * GluonColAvg**2 * SymmFact
      NDim = NDim + 8    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 2000000
      VegasNc1_default = 2000000

  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=3
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 5    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 500000
      VegasNc1_default = 500000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF

ELSEIF( PROCESS.EQ.10 ) THEN !   3_Str  + 4_AStr --> 5_Glu  + 6_Glu  + 1_ATop + 2_Top
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 6
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/5,6,-1,-2,3,4/)
      MasterProcess=6
      AvgFactor = SpinAvg * QuarkColAvg**2  * SymmFact
      NDim = NDim + 8    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 2000000
      VegasNc1_default = 2000000

  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 2000000
      VegasNc1_default = 2000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF



ELSEIF( PROCESS.EQ.11 ) THEN !   5_Glu  + 6_Glu --> 3_Str  + 4_AStr + 1_ATop + 2_Top
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 6
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/5,6,3,4,-1,-2/)
      MasterProcess=6
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 8    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 2000000
      VegasNc1_default = 2000000

  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,3,-1,-2/)
      MasterProcess=3
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 5    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 500000
      VegasNc1_default = 500000

  ELSE
      call Error("Correction to this process is not available")
  ENDIF


ELSEIF( PROCESS.EQ.12 ) THEN !   3_Str + 5_Glu --> 4_Str  + 1_ATop + 2_Top + 6_Glu
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 6
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/5,6,-1,3,-2,4/)
      MasterProcess=6
      AvgFactor = SpinAvg * QuarkColAvg * GluonColAvg
      NDim = NDim + 8    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 2000000
      VegasNc1_default = 2000000

      VegasNc1_default = 2000000

  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,3,-2/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg * GluonColAvg
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 2000000
      VegasNc1_default = 2000000

  ELSE
      call Error("Correction to this process is not available")
  ENDIF



ELSEIF( PROCESS.EQ.13 ) THEN !   4_AStr + 5_Glu --> 3_Str  + 1_ATop + 2_Top + 6_Glu
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 6
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/5,6,3,-1,-2,4/)
      MasterProcess=6
      AvgFactor = SpinAvg * QuarkColAvg * GluonColAvg
      NDim = NDim + 8    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 2000000
      VegasNc1_default = 2000000

  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,3,-1,-2/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg*GluonColAvg
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 2000000
      VegasNc1_default = 2000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF


ELSEIF( PROCESS.EQ.14 ) THEN !   3_Str + 4_AStr -->  5_Chm + 6_AChm + 1_ATop + 2_Top //  5_Str + 6_AStr + 1_ATop + 2_Top
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 6
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/5,6,-1,-2,3,4/)
      MasterProcess=7
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 8    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000

  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,3,-1,-2/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg*QuarkColAvg
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 2000000
      VegasNc1_default = 2000000

  ELSE
      call Error("Correction to this process is not available")
  ENDIF


ELSEIF( PROCESS.EQ.15 ) THEN !   3_Str + 6_AChm -->  4_Str + 5_AChm + 1_ATop + 2_Top
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 6
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/5,6,-1,3,4,-2/)
      MasterProcess=7
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 8    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000000
      VegasNc1_default = 50

  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,3,-2/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg*QuarkColAvg
      NDim = NDim + 5    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 2000000
      VegasNc1_default = 2000000

  ELSE
      call Error("Correction to this process is not available")
  ENDIF


ELSEIF( PROCESS.EQ.16 ) THEN !   3_Str + 6_Chm -->  4_Str + 5_Chm + 1_ATop + 2_Top
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 6
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/5,6,-1,3,-2,4/)
      MasterProcess=7
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 8    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000000
      VegasNc1_default = 50

  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,3,-2/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg*QuarkColAvg
      NDim = NDim + 5    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 2000000
      VegasNc1_default = 2000000

  ELSE
      call Error("Correction to this process is not available")
  ENDIF


ELSEIF( PROCESS.EQ.17 ) THEN !   4_AStr + 5_AChm -->  3_AStr + 5_AChm + 1_ATop + 2_Top
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 6
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/5,6,3,-1,4,-2/)
      MasterProcess=7
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 8    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000000
      VegasNc1_default = 50

  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,3,-1,-2/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg*QuarkColAvg
      NDim = NDim + 5    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 2000000
      VegasNc1_default = 2000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF


ELSEIF( PROCESS.EQ.18 ) THEN !   3_Chm + 5_Chm -->  4_Chm + 6_Chm + 1_ATop + 2_Top
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 6
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/5,6,-1,3,-2,4/)
      MasterProcess=7
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 8    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000000
      VegasNc1_default = 50

  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,3,-2/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg*QuarkColAvg
      NDim = NDim + 5    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 2000000
      VegasNc1_default = 2000000

  ELSE
      call Error("Correction to this process is not available")
  ENDIF



ELSEIF( PROCESS.EQ.19 ) THEN !   4_AChm + 6_AChm -->  3_AChm + 5_AChm + 1_ATop + 2_Top
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 6
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/5,6,3,-1,4,-2/)
      MasterProcess=7
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 8    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000000
      VegasNc1_default = 50

  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-3,-1,-2/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg*QuarkColAvg
      NDim = NDim + 5    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 2000000
      VegasNc1_default = 2000000

  ELSE
      call Error("Correction to this process is not available")
  ENDIF



ELSEIF( PROCESS.EQ.20 ) THEN !   3_Glu  + 4_Glu  --> 1_ATop + 2_Top + 5_Pho(in production)
  IF( CORRECTION.EQ.0 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=8
      NDim = NDim + 5    ! PS integration
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 100000
      VegasNc1_default = 100000
  ELSEIF( CORRECTION.EQ.1 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=8
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 5    ! t tbar photon PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 100000
      VegasNc1_default = 100000
  ELSEIF( CORRECTION.EQ.4 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=8
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 5    ! t tbar photon PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 100000
      VegasNc1_default = 100000
  ELSEIF( CORRECTION.EQ.5 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=8
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 5    ! t tbar photon PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 3  ! additional gluons in the top decay
      VegasNc0_default = 1000000
      VegasNc1_default = 1000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF


ELSEIF( PROCESS.EQ.21 ) THEN !   3_Glu  + 4_Glu  --> 1_ATop + 2_Top + 5_Pho(in decay)
  IF( CORRECTION.LE.1 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=1
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 3    ! photon phase space
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 100000
      VegasNc1_default = 100000
  ELSEIF( CORRECTION.EQ.4 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=1
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 5    ! t tbar photon PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 1000000
      VegasNc1_default = 1000000
  ELSEIF( CORRECTION.EQ.5 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=1
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 5    ! t tbar photon PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 3    ! additional gluons in the top decay
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF


ELSEIF( PROCESS.EQ.22 ) THEN !   3_Str  + 4_AStr --> 1_ATop + 2_Top + 5_Pho(in production)
  IF( CORRECTION.EQ.0 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=9
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 5    ! t tbar photon PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSEIF( CORRECTION.EQ.1 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=9
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 5    ! t tbar photon PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 100000
      VegasNc1_default = 100000
  ELSEIF( CORRECTION.EQ.4 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=9
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 5    ! t tbar photon PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 1000000
      VegasNc1_default = 1000000
  ELSEIF( CORRECTION.EQ.5 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=9
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 5    ! t tbar photon PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 3  ! additional gluons in the top decay
      VegasNc0_default = 1000000
      VegasNc1_default = 1000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF



ELSEIF( PROCESS.EQ.23 ) THEN !   3_Str  + 4_AStr --> 1_ATop + 2_Top + 5_Pho(in decay)
  IF( CORRECTION.LE.1 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 3    ! photon phase space
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 100000
      VegasNc1_default = 100000
  ELSEIF( CORRECTION.EQ.4 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 5    ! t tbar photon PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 1000000
      VegasNc1_default = 1000000
  ELSEIF( CORRECTION.EQ.5 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 5    ! t tbar photon PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 3    ! additional gluons in the top decay
      VegasNc0_default = 1000000
      VegasNc1_default = 1000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF



ELSEIF( PROCESS.EQ.24 ) THEN !   3_Str  + 5_Glu  --> 4_Str  + 1_ATop + 2_Top + 6_Pho(in production)
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 6
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/5,6,-1,3,-2,4/)
      MasterProcess=11
      AvgFactor = SpinAvg * QuarkColAvg * GluonColAvg
      NDim = NDim + 8    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=9
      AvgFactor = SpinAvg * QuarkColAvg * GluonColAvg
      NDim = NDim + 5    ! t tbar photon PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF



ELSEIF( PROCESS.EQ.25 ) THEN  !   3_Str  + 5_Glu  --> 4_Str  + 1_ATop + 2_Top + 6_Pho(in decay)
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,3,-2/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg * GluonColAvg
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 3    ! photon phase space
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg * GluonColAvg
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 3    ! photon phase space
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF



ELSEIF( PROCESS.EQ.26 ) THEN !   4_AStr + 5_Glu  --> 3_AStr + 1_ATop + 2_Top + 6_Pho(in production)
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 6
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/5,6,3,-1,-2,4/)
      MasterProcess=11
      AvgFactor = SpinAvg * QuarkColAvg * GluonColAvg
      NDim = NDim + 8    ! t tbar glu photon PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=9
      AvgFactor = SpinAvg * QuarkColAvg * GluonColAvg
      NDim = NDim + 5    ! t tbar photon PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF



ELSEIF( PROCESS.EQ.27 ) THEN !   4_AStr + 5_Glu  --> 3_AStr + 1_ATop + 2_Top + 6_Pho(in decay)
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,3,-1,-2/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg * GluonColAvg
      NDim = NDim + 8    ! t tbar glu photon PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg * GluonColAvg
      NDim = NDim + 5    ! t tbar photon PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF



ELSEIF( PROCESS.EQ.28 ) THEN !   3_Glu  + 4_Glu  --> 5_Glu  + 1_ATop + 2_Top + 6_Pho(in production)
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 6
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/5,6,-1,-2,3,4/)
      MasterProcess=10
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 8    ! t tbar glu photon PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=8
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 5    ! t tbar photon PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF


ELSEIF( PROCESS.EQ.29 ) THEN !   3_Glu  + 4_Glu  --> 5_Glu  + 1_ATop + 2_Top + 6_Pho(in decay)
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=3
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 8    ! t tbar glu photon PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=1
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 5    ! t tbar PS photon integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF



ELSEIF( PROCESS.EQ.30 ) THEN !   3_Str  + 4_AStr --> 5_Glu  + 1_ATop + 2_Top + 6_Pho(in production)
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 6
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/5,6,-1,-2,3,4/)
      MasterProcess=11
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 8    ! t tbar glu photon PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=9
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 5    ! t tbar photon PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF


ELSEIF( PROCESS.EQ.31 ) THEN !   3_Str  + 4_AStr --> 5_Glu  + 1_ATop + 2_Top + 6_Pho(in decay)
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 8    ! t tbar glu photon PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 1000000
      VegasNc1_default = 1000000
  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 5    ! t tbar photon PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 1000000
      VegasNc1_default = 1000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF


ELSEIF( PROCESS.EQ.33 ) THEN !   3_Glu  + 4_Glu  --> 1_ATop + 2_Top + 5_Glu(in decay)
  IF( CORRECTION.LE.1 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=1
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 3    ! gluon phase space
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 100000
      VegasNc1_default = 100000
  ELSEIF( CORRECTION.EQ.4 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=1
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 3    ! gluon phase space
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 1000000
      VegasNc1_default = 1000000
  ELSEIF( CORRECTION.EQ.5 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=1
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 5    ! t tbar gluon PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 3    ! additional gluon in the top decay
      VegasNc0_default = 1000000
      VegasNc1_default = 1000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF



ELSEIF( PROCESS.EQ.34 ) THEN !   3_Str  + 4_AStr --> 1_ATop + 2_Top + 5_Glu(in decay)
  IF( CORRECTION.LE.1 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 3    ! gluon phase space
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 100000
      VegasNc1_default = 100000
  ELSEIF( CORRECTION.EQ.4 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 5    ! t tbar gluon PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 1000000
      VegasNc1_default = 1000000
  ELSEIF( CORRECTION.EQ.5 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 5    ! t tbar gluon PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 3    ! additional gluon in the top decay
      VegasNc0_default = 1000000
      VegasNc1_default = 1000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF




ELSEIF( PROCESS.EQ.35 ) THEN  !   3_Str  + 5_Glu  --> 4_Str  + 1_ATop + 2_Top + 6_Glu(in decay)
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,3,-2/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg * GluonColAvg
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 3    ! gluon phase space
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg * GluonColAvg
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 3    ! gluon phase space
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF


ELSEIF( PROCESS.EQ.36 ) THEN !   4_AStr + 5_Glu  --> 3_AStr + 1_ATop + 2_Top + 6_Glu(in decay)
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,3,-1,-2/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg * GluonColAvg
      NDim = NDim + 8    ! t tbar glu glu PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg * GluonColAvg
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF




ELSEIF( PROCESS.EQ.37 ) THEN !   3_Glu  + 4_Glu  --> 5_Glu  + 1_ATop + 2_Top + 6_Glu(in decay)
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=3
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 8    ! t tbar glu glu PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=1
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 5    ! t tbar PS glu integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF




ELSEIF( PROCESS.EQ.38 ) THEN !   3_Str  + 4_AStr --> 5_Glu  + 1_ATop + 2_Top + Glu(in decay)
  IF( CORRECTION.EQ.2 ) THEN
      NumExtParticles = 5
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/4,5,-1,-2,3/)
      MasterProcess=4
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 8    ! t tbar glu glu PS integration
      NDim = NDim + 2    ! shat integration
      VegasNc0_default = 1000000
      VegasNc1_default = 1000000
  ELSEIF( CORRECTION.EQ.3 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 5    ! t tbar glu PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! x integration
      VegasNc0_default = 1000000
      VegasNc1_default = 1000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF




ELSEIF( PROCESS.EQ.41 ) THEN !   3_Glu  + 4_Glu  --> 1_AHeavyTop + 2_HeavyTop

! temporarily reset m_Top for InitMasterprocess and InitProcess
 m_SMTop = m_Top
 m_Top = m_HTop
! will be restored in StartVegas

  IF( CORRECTION.EQ.0 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=1
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 4    ! T -> A0+t decays
      VegasNc0_default = 100000
      VegasNc1_default = 100000
  ELSEIF( CORRECTION.EQ.1 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=1
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 4    ! T -> A0+t decays
      VegasNc0_default = 50000
      VegasNc1_default = 50000
  ELSEIF( CORRECTION.EQ.4 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=1
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 4    ! T -> A0+t decays
      VegasNc0_default = 1000000
      VegasNc1_default = 1000000
  ELSEIF( CORRECTION.EQ.5 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=1
      AvgFactor = SpinAvg * GluonColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 4    ! T -> A0+t decays
      NDim = NDim + 3+3  ! additional gluons in the top decay
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF



ELSEIF( PROCESS.EQ.42 ) THEN !   3_Str  + 4_AStr --> 1_AHeavyTop + 2_HeavyTop

! temporarily reset m_Top for InitMasterprocess and InitProcess
 m_SMTop = m_Top
 m_Top = m_HTop
! will be restored in StartVegas

  IF( CORRECTION.EQ.0 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 4    ! T -> A0+t decays
      VegasNc0_default = 100000
      VegasNc1_default = 100000
  ELSEIF( CORRECTION.EQ.1 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 4    ! T -> A0+t decays
      VegasNc0_default = 50000
      VegasNc1_default = 50000
  ELSEIF( CORRECTION.EQ.4 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 4    ! T -> A0+t decays
      VegasNc0_default = 1000000
      VegasNc1_default = 1000000
  ELSEIF( CORRECTION.EQ.5 ) THEN
      NumExtParticles = 4
      allocate(Crossing(1:NumExtParticles))
      allocate(ExtParticle(1:NumExtParticles))
      Crossing(:) = (/3,4,-1,-2/)
      MasterProcess=2
      AvgFactor = SpinAvg * QuarkColAvg**2
      NDim = NDim + 2    ! t tbar PS integration
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 4    ! T -> A0+t decays
      NDim = NDim + 3+3  ! additional gluons in the top decay
      VegasNc0_default = 10000000
      VegasNc1_default = 10000000
  ELSE
      call Error("Correction to this process is not available")
  ENDIF






ELSE
    call Error("Process not available")
ENDIF

  call InitMasterProcess()

END SUBROUTINE







SUBROUTINE InitMasterProcess()
use ModParameters
use ModMisc
implicit none
integer :: NPart,sig_tb,sig_t,NAmp

! print *, ""
! print *, "Initializing Process   ", Process
! print *, "             Master    ", MasterProcess
! print *, "             Correction", Correction


IF( MASTERPROCESS.EQ.0 ) THEN

    ExtParticle(1)%PartType = ATop_
    ExtParticle(2)%PartType = Top_
    ExtParticle(3)%PartType = Glu_
    ExtParticle(4)%PartType = Glu_
    ExtParticle(5)%PartType = Glu_
    ExtParticle(6)%PartType = Glu_
    ExtParticle(7)%PartType = Glu_
    NumPrimAmps = 1
    NumBornAmps = 1
    allocate(PrimAmps(1:NumPrimAmps))
    allocate(BornAmps(1:NumPrimAmps))
    do NAmp=1,NumPrimAmps
        allocate(BornAmps(NAmp)%ExtLine(1:NumExtParticles))
        allocate(PrimAmps(NAmp)%ExtLine(1:NumExtParticles))
        allocate(PrimAmps(NAmp)%IntPart(1:NumExtParticles))
    enddo
    NumHelicities = 1
    allocate(Helicities(1:NumHelicities,1:NumExtParticles))
    Helicities(1,1:NumExtParticles) = (/+1,+1,-1,-1,-1,-1,+1/)


ELSEIF( MASTERPROCESS.EQ.1 ) THEN

    ExtParticle(1)%PartType = ATop_
    ExtParticle(2)%PartType = Top_
    ExtParticle(3)%PartType = Glu_
    ExtParticle(4)%PartType = Glu_
    IF( Correction.EQ.0 .OR. Correction.GE.4 ) THEN
      NumPrimAmps = 2
      NumBornAmps = 2
    ELSEIF( Correction.EQ.1 ) THEN
      NumPrimAmps = 10
      NumBornAmps = 2
    ELSEIF( Correction.EQ.3 ) THEN
      NumPrimAmps = 0
      NumBornAmps = 0
    ENDIF
    allocate(PrimAmps(1:NumPrimAmps))
    allocate(BornAmps(1:NumPrimAmps))
    do NAmp=1,NumPrimAmps
        allocate(BornAmps(NAmp)%ExtLine(1:NumExtParticles))
        allocate(PrimAmps(NAmp)%ExtLine(1:NumExtParticles))
        allocate(PrimAmps(NAmp)%IntPart(1:NumExtParticles))
    enddo

    IF( TOPDECAYS.GE.1 ) THEN
              NumHelicities = 4
              allocate(Helicities(1:NumHelicities,1:NumExtParticles))
              Helicities(1,1:4) = (/0,0,+1,+1/)
              Helicities(2,1:4) = (/0,0,+1,-1/)
              Helicities(3,1:4) = (/0,0,-1,+1/)
              Helicities(4,1:4) = (/0,0,-1,-1/)
    ELSE
              NumHelicities = 8
              allocate(Helicities(1:NumHelicities,1:NumExtParticles))
              sig_tb=+1; sig_t =+1;
              Helicities(1,1:NumExtParticles) = (/sig_tb,sig_t,+1,+1/)
              Helicities(2,1:NumExtParticles) = (/sig_tb,sig_t,+1,-1/)
              sig_tb=+1; sig_t =-1;
              Helicities(3,1:NumExtParticles) = (/sig_tb,sig_t,+1,+1/)
              Helicities(4,1:NumExtParticles) = (/sig_tb,sig_t,+1,-1/)
              sig_tb=-1; sig_t =+1;
              Helicities(5,1:NumExtParticles) = (/sig_tb,sig_t,+1,+1/)
              Helicities(6,1:NumExtParticles) = (/sig_tb,sig_t,+1,-1/)
              sig_tb=-1; sig_t =-1;
              Helicities(7,1:NumExtParticles) = (/sig_tb,sig_t,+1,+1/)
              Helicities(8,1:NumExtParticles) = (/sig_tb,sig_t,+1,-1/)
    !  additional helicities when parity inversion is not applied:    changes affect also EvalCS_ttb_NLODK_noSC
    !         sig_tb=+1; sig_t =+1;
    !         Helicities(9 ,1:NumExtParticles) = (/sig_tb,sig_t,-1,-1/)
    !         Helicities(10,1:NumExtParticles) = (/sig_tb,sig_t,-1,+1/)
    !         sig_tb=+1; sig_t =-1;
    !         Helicities(11,1:NumExtParticles) = (/sig_tb,sig_t,-1,-1/)
    !         Helicities(12,1:NumExtParticles) = (/sig_tb,sig_t,-1,+1/)
    !         sig_tb=-1; sig_t =+1;
    !         Helicities(13,1:NumExtParticles) = (/sig_tb,sig_t,-1,-1/)
    !         Helicities(14,1:NumExtParticles) = (/sig_tb,sig_t,-1,+1/)
    !         sig_tb=-1; sig_t =-1;
    !         Helicities(15,1:NumExtParticles) = (/sig_tb,sig_t,-1,-1/)
    !         Helicities(16,1:NumExtParticles) = (/sig_tb,sig_t,-1,+1/)
    ENDIF



ELSEIF( MASTERPROCESS.EQ.2 ) THEN

    ExtParticle(1)%PartType = ATop_
    ExtParticle(2)%PartType = Top_
    ExtParticle(3)%PartType = AStr_
    ExtParticle(4)%PartType = Str_
    IF( Correction.EQ.0 .OR. Correction.GE.4 ) THEN
      NumPrimAmps = 1
      NumBornAmps = 1
    ELSEIF( Correction.EQ.1 ) THEN
      NumPrimAmps = 6
      NumBornAmps = 1
    ELSEIF( Correction.EQ.3 ) THEN
      NumPrimAmps = 0
      NumBornAmps = 0
    ENDIF
    allocate(PrimAmps(1:NumPrimAmps))
    allocate(BornAmps(1:NumPrimAmps))
    do NAmp=1,NumPrimAmps
        allocate(BornAmps(NAmp)%ExtLine(1:NumExtParticles))
        allocate(PrimAmps(NAmp)%ExtLine(1:NumExtParticles))
        allocate(PrimAmps(NAmp)%IntPart(1:NumExtParticles))
    enddo

    IF( TOPDECAYS.GE.1 ) THEN
              NumHelicities = 4
              allocate(Helicities(1:NumHelicities,1:NumExtParticles))
              Helicities(1,1:4) = (/0,0,+1,+1/)
              Helicities(2,1:4) = (/0,0,+1,-1/)
              Helicities(3,1:4) = (/0,0,-1,+1/)
              Helicities(4,1:4) = (/0,0,-1,-1/)
    ELSE
      NumHelicities = 4
      allocate(Helicities(1:NumHelicities,1:NumExtParticles))
      sig_tb=+1; sig_t =+1;
  !    Helicities(1,1:NumExtParticles) = (/sig_tb,sig_t,+1,+1/)  ! the x,x,+1,+1 helicities lead to vanishing tree contribution
      Helicities(1,1:NumExtParticles) = (/sig_tb,sig_t,+1,-1/)
      sig_tb=+1; sig_t =-1;
  !    Helicities(3,1:NumExtParticles) = (/sig_tb,sig_t,+1,+1/)
      Helicities(2,1:NumExtParticles) = (/sig_tb,sig_t,+1,-1/)
      sig_tb=-1; sig_t =+1;
  !    Helicities(5,1:NumExtParticles) = (/sig_tb,sig_t,+1,+1/)
      Helicities(3,1:NumExtParticles) = (/sig_tb,sig_t,+1,-1/)
      sig_tb=-1; sig_t =-1;
  !    Helicities(7,1:NumExtParticles) = (/sig_tb,sig_t,+1,+1/)
      Helicities(4,1:NumExtParticles) = (/sig_tb,sig_t,+1,-1/)
  !   additional helicities when parity inversion is not applied: changes affect also EvalCS_ttb_NLODK_noSC
  !     sig_tb=+1; sig_t =+1;
  !     Helicities(9 ,1:NumExtParticles) = (/sig_tb,sig_t,-1,-1/)
  !     Helicities(10,1:NumExtParticles) = (/sig_tb,sig_t,-1,+1/)
  !     sig_tb=+1; sig_t =-1;
  !     Helicities(11,1:NumExtParticles) = (/sig_tb,sig_t,-1,-1/)
  !     Helicities(12,1:NumExtParticles) = (/sig_tb,sig_t,-1,+1/)
  !     sig_tb=-1; sig_t =+1;
  !     Helicities(13,1:NumExtParticles) = (/sig_tb,sig_t,-1,-1/)
  !     Helicities(14,1:NumExtParticles) = (/sig_tb,sig_t,-1,+1/)
  !     sig_tb=-1; sig_t =-1;
  !     Helicities(15,1:NumExtParticles) = (/sig_tb,sig_t,-1,-1/)
  !     Helicities(16,1:NumExtParticles) = (/sig_tb,sig_t,-1,+1/)
    ENDIF



ELSEIF( MASTERPROCESS.EQ.3 ) THEN

    ExtParticle(1)%PartType = ATop_
    ExtParticle(2)%PartType = Top_
    ExtParticle(3)%PartType = Glu_
    ExtParticle(4)%PartType = Glu_
    ExtParticle(5)%PartType = Glu_

    IF( Correction.EQ.0 .OR.  Correction.EQ.2 .OR.  Correction.EQ.4 ) THEN
      NumPrimAmps = 6
      NumBornAmps = 6
    ELSEIF( Correction.EQ.1 ) THEN
      NumPrimAmps = 48
      NumBornAmps = 6
    ENDIF
    allocate(PrimAmps(1:NumPrimAmps))
    allocate(BornAmps(1:NumPrimAmps))
    do NAmp=1,NumPrimAmps
        allocate(BornAmps(NAmp)%ExtLine(1:NumExtParticles))
        allocate(PrimAmps(NAmp)%ExtLine(1:NumExtParticles))
        allocate(PrimAmps(NAmp)%IntPart(1:NumExtParticles))
    enddo

    IF( TOPDECAYS.GE.1 ) THEN
        NumHelicities = 8
        allocate(Helicities(1:NumHelicities,1:NumExtParticles))
        Helicities(1,1:5) = (/0,0,+1,+1,+1/)
        Helicities(2,1:5) = (/0,0,+1,+1,-1/)
        Helicities(3,1:5) = (/0,0,+1,-1,+1/)
        Helicities(4,1:5) = (/0,0,+1,-1,-1/)
        Helicities(5,1:5) = (/0,0,-1,+1,+1/)
        Helicities(6,1:5) = (/0,0,-1,+1,-1/)
        Helicities(7,1:5) = (/0,0,-1,-1,+1/)
        Helicities(8,1:5) = (/0,0,-1,-1,-1/)
    ELSE
        NumHelicities = 32
        allocate(Helicities(1:NumHelicities,1:NumExtParticles))
        sig_tb=+1; sig_t =+1;
        Helicities( 1,1:5) = (/sig_tb,sig_t,+1,+1,+1/)
        Helicities( 2,1:5) = (/sig_tb,sig_t,+1,+1,-1/)
        Helicities( 3,1:5) = (/sig_tb,sig_t,+1,-1,-1/)
        Helicities( 4,1:5) = (/sig_tb,sig_t,+1,-1,+1/)
        sig_tb=+1; sig_t =-1;
        Helicities( 5,1:5) = (/sig_tb,sig_t,+1,+1,+1/)
        Helicities( 6,1:5) = (/sig_tb,sig_t,+1,+1,-1/)
        Helicities( 7,1:5) = (/sig_tb,sig_t,+1,-1,-1/)
        Helicities( 8,1:5) = (/sig_tb,sig_t,+1,-1,+1/)
        sig_tb=-1; sig_t =+1;
        Helicities( 9,1:5) = (/sig_tb,sig_t,+1,+1,+1/)
        Helicities(10,1:5) = (/sig_tb,sig_t,+1,+1,-1/)
        Helicities(11,1:5) = (/sig_tb,sig_t,+1,-1,-1/)
        Helicities(12,1:5) = (/sig_tb,sig_t,+1,-1,+1/)
        sig_tb=-1; sig_t =-1;
        Helicities(13,1:5) = (/sig_tb,sig_t,+1,+1,+1/)
        Helicities(14,1:5) = (/sig_tb,sig_t,+1,+1,-1/)
        Helicities(15,1:5) = (/sig_tb,sig_t,+1,-1,-1/)
        Helicities(16,1:5) = (/sig_tb,sig_t,+1,-1,+1/)
    !   additional helicities when parity inversion is not applied:
        sig_tb=-1; sig_t =-1;
        Helicities(17,1:5) = (/sig_tb,sig_t,-1,-1,-1/)
        Helicities(18,1:5) = (/sig_tb,sig_t,-1,-1,+1/)
        Helicities(19,1:5) = (/sig_tb,sig_t,-1,+1,+1/)
        Helicities(20,1:5) = (/sig_tb,sig_t,-1,+1,-1/)
        sig_tb=-1; sig_t =+1;
        Helicities(21,1:5) = (/sig_tb,sig_t,-1,-1,-1/)
        Helicities(22,1:5) = (/sig_tb,sig_t,-1,-1,+1/)
        Helicities(23,1:5) = (/sig_tb,sig_t,-1,+1,+1/)
        Helicities(24,1:5) = (/sig_tb,sig_t,-1,+1,-1/)
        sig_tb=+1; sig_t =-1;
        Helicities(25,1:5) = (/sig_tb,sig_t,-1,-1,-1/)
        Helicities(26,1:5) = (/sig_tb,sig_t,-1,-1,+1/)
        Helicities(27,1:5) = (/sig_tb,sig_t,-1,+1,+1/)
        Helicities(28,1:5) = (/sig_tb,sig_t,-1,+1,-1/)
        sig_tb=+1; sig_t =+1;
        Helicities(29,1:5) = (/sig_tb,sig_t,-1,-1,-1/)
        Helicities(30,1:5) = (/sig_tb,sig_t,-1,-1,+1/)
        Helicities(31,1:5) = (/sig_tb,sig_t,-1,+1,+1/)
        Helicities(32,1:5) = (/sig_tb,sig_t,-1,+1,-1/)
    ENDIF



ELSEIF( MASTERPROCESS.EQ.4 ) THEN

    ExtParticle(1)%PartType = ATop_
    ExtParticle(2)%PartType = Top_
    ExtParticle(3)%PartType = AStr_
    ExtParticle(4)%PartType = Str_
    ExtParticle(5)%PartType = Glu_

    IF( Correction.EQ.0 .OR.  Correction.EQ.2 .OR.  Correction.EQ.4 ) THEN
      NumPrimAmps = 4
      NumBornAmps = 4
    ELSEIF( Correction.EQ.1 ) THEN
      NumPrimAmps = 24
      NumBornAmps = 4
    ENDIF
    allocate(PrimAmps(1:NumPrimAmps))
    allocate(BornAmps(1:NumPrimAmps))
    do NAmp=1,NumPrimAmps
        allocate(BornAmps(NAmp)%ExtLine(1:NumExtParticles))
        allocate(PrimAmps(NAmp)%ExtLine(1:NumExtParticles))
        allocate(PrimAmps(NAmp)%IntPart(1:NumExtParticles))
    enddo

    IF( TOPDECAYS.GE.1 ) THEN
        NumHelicities = 8
        allocate(Helicities(1:NumHelicities,1:NumExtParticles))
        Helicities(1,1:5) = (/0,0,+1,+1,+1/)
        Helicities(2,1:5) = (/0,0,+1,+1,-1/)
        Helicities(3,1:5) = (/0,0,+1,-1,+1/)
        Helicities(4,1:5) = (/0,0,+1,-1,-1/)
        Helicities(5,1:5) = (/0,0,-1,+1,+1/)
        Helicities(6,1:5) = (/0,0,-1,+1,-1/)
        Helicities(7,1:5) = (/0,0,-1,-1,+1/)
        Helicities(8,1:5) = (/0,0,-1,-1,-1/)
    !    print *, "remove vanishing helicities!"
    !    stop

    ELSE
!     NumHelicities = 16
    NumHelicities = 32
    allocate(Helicities(1:NumHelicities,1:NumExtParticles))
!     sig_tb=+1; sig_t =+1;         ! works for qqb initial state but not for qg and qbg because ++/-- helicities don't vanish
! !    Helicities( 1,1:5) = (/sig_tb,sig_t,+1,+1,+1/)!
! !    Helicities( 2,1:5) = (/sig_tb,sig_t,+1,+1,-1/)!
!     Helicities( 1,1:5) = (/sig_tb,sig_t,+1,-1,-1/)
!     Helicities( 2,1:5) = (/sig_tb,sig_t,+1,-1,+1/)
!     sig_tb=+1; sig_t =-1;
! !    Helicities( 5,1:5) = (/sig_tb,sig_t,+1,+1,+1/)!
! !    Helicities( 6,1:5) = (/sig_tb,sig_t,+1,+1,-1/)!
!     Helicities( 3,1:5) = (/sig_tb,sig_t,+1,-1,-1/)
!     Helicities( 4,1:5) = (/sig_tb,sig_t,+1,-1,+1/)
!     sig_tb=-1; sig_t =+1;
! !    Helicities( 9,1:5) = (/sig_tb,sig_t,+1,+1,+1/)!
! !    Helicities(10,1:5) = (/sig_tb,sig_t,+1,+1,-1/)!
!     Helicities(5,1:5) = (/sig_tb,sig_t,+1,-1,-1/)
!     Helicities(6,1:5) = (/sig_tb,sig_t,+1,-1,+1/)
!     sig_tb=-1; sig_t =-1;
! !    Helicities(13,1:5) = (/sig_tb,sig_t,+1,+1,+1/)!
! !    Helicities(14,1:5) = (/sig_tb,sig_t,+1,+1,-1/)!
!     Helicities(7,1:5) = (/sig_tb,sig_t,+1,-1,-1/)
!     Helicities(8,1:5) = (/sig_tb,sig_t,+1,-1,+1/)
! ! !   additional helicities when parity inversion is not applied:
!       sig_tb=-1; sig_t =-1;
! !      Helicities(17,1:5) = (/sig_tb,sig_t,-1,-1,-1/)!
! !      Helicities(18,1:5) = (/sig_tb,sig_t,-1,-1,+1/)!
!       Helicities(9,1:5)  = (/sig_tb,sig_t,-1,+1,+1/)
!       Helicities(10,1:5) = (/sig_tb,sig_t,-1,+1,-1/)
!       sig_tb=-1; sig_t =+1;
! !      Helicities(21,1:5) = (/sig_tb,sig_t,-1,-1,-1/)!
! !      Helicities(22,1:5) = (/sig_tb,sig_t,-1,-1,+1/)!
!       Helicities(11,1:5) = (/sig_tb,sig_t,-1,+1,+1/)
!       Helicities(12,1:5) = (/sig_tb,sig_t,-1,+1,-1/)
!       sig_tb=+1; sig_t =-1;
! !      Helicities(25,1:5) = (/sig_tb,sig_t,-1,-1,-1/)!
! !      Helicities(26,1:5) = (/sig_tb,sig_t,-1,-1,+1/)!
!       Helicities(13,1:5) = (/sig_tb,sig_t,-1,+1,+1/)
!       Helicities(14,1:5) = (/sig_tb,sig_t,-1,+1,-1/)
!       sig_tb=+1; sig_t =+1;
! !      Helicities(29,1:5) = (/sig_tb,sig_t,-1,-1,-1/)!
! !      Helicities(30,1:5) = (/sig_tb,sig_t,-1,-1,+1/)!
!       Helicities(15,1:5) = (/sig_tb,sig_t,-1,+1,+1/)
!       Helicities(16,1:5) = (/sig_tb,sig_t,-1,+1,-1/)



    sig_tb=+1; sig_t =+1;
    Helicities( 1,1:5) = (/sig_tb,sig_t,+1,-1,-1/)
    Helicities( 2,1:5) = (/sig_tb,sig_t,+1,-1,+1/)
    sig_tb=+1; sig_t =-1;
    Helicities( 3,1:5) = (/sig_tb,sig_t,+1,-1,-1/)
    Helicities( 4,1:5) = (/sig_tb,sig_t,+1,-1,+1/)
    sig_tb=-1; sig_t =+1;
    Helicities( 5,1:5) = (/sig_tb,sig_t,+1,-1,-1/)
    Helicities( 6,1:5) = (/sig_tb,sig_t,+1,-1,+1/)
    sig_tb=-1; sig_t =-1;
    Helicities( 7,1:5) = (/sig_tb,sig_t,+1,-1,-1/)
    Helicities( 8,1:5) = (/sig_tb,sig_t,+1,-1,+1/)

    sig_tb=-1; sig_t =-1;
    Helicities( 9,1:5) = (/sig_tb,sig_t,-1,+1,+1/)
    Helicities(10,1:5) = (/sig_tb,sig_t,-1,+1,-1/)
    sig_tb=-1; sig_t =+1;
    Helicities(11,1:5) = (/sig_tb,sig_t,-1,+1,+1/)
    Helicities(12,1:5) = (/sig_tb,sig_t,-1,+1,-1/)
    sig_tb=+1; sig_t =-1;
    Helicities(13,1:5) = (/sig_tb,sig_t,-1,+1,+1/)
    Helicities(14,1:5) = (/sig_tb,sig_t,-1,+1,-1/)
    sig_tb=+1; sig_t =+1;
    Helicities(15,1:5) = (/sig_tb,sig_t,-1,+1,+1/)
    Helicities(16,1:5) = (/sig_tb,sig_t,-1,+1,-1/)

!   these helicities don't contribute for qqb:
    if( PROCESS.EQ.6 ) NumHelicities = 16

    sig_tb=+1; sig_t =+1;
    Helicities(17,1:5) = (/sig_tb,sig_t,+1,+1,+1/)!
    Helicities(18,1:5) = (/sig_tb,sig_t,+1,+1,-1/)!
    sig_tb=+1; sig_t =-1;
    Helicities(19,1:5) = (/sig_tb,sig_t,+1,+1,+1/)!
    Helicities(20,1:5) = (/sig_tb,sig_t,+1,+1,-1/)!
    sig_tb=-1; sig_t =+1;
    Helicities(21,1:5) = (/sig_tb,sig_t,+1,+1,+1/)!
    Helicities(22,1:5) = (/sig_tb,sig_t,+1,+1,-1/)!
    sig_tb=-1; sig_t =-1;
    Helicities(23,1:5) = (/sig_tb,sig_t,+1,+1,+1/)!
    Helicities(24,1:5) = (/sig_tb,sig_t,+1,+1,-1/)!

    sig_tb=-1; sig_t =-1;
    Helicities(25,1:5) = (/sig_tb,sig_t,-1,-1,-1/)!
    Helicities(26,1:5) = (/sig_tb,sig_t,-1,-1,+1/)!
    sig_tb=-1; sig_t =+1;
    Helicities(27,1:5) = (/sig_tb,sig_t,-1,-1,-1/)!
    Helicities(28,1:5) = (/sig_tb,sig_t,-1,-1,+1/)!
    sig_tb=+1; sig_t =-1;
    Helicities(29,1:5) = (/sig_tb,sig_t,-1,-1,-1/)!
    Helicities(30,1:5) = (/sig_tb,sig_t,-1,-1,+1/)!
    sig_tb=+1; sig_t =+1;
    Helicities(31,1:5) = (/sig_tb,sig_t,-1,-1,-1/)!
    Helicities(32,1:5) = (/sig_tb,sig_t,-1,-1,+1/)!

    ENDIF


ELSEIF( MASTERPROCESS.EQ.5 ) THEN

    ExtParticle(1)%PartType = ATop_
    ExtParticle(2)%PartType = Top_
    ExtParticle(3)%PartType = Glu_
    ExtParticle(4)%PartType = Glu_
    ExtParticle(5)%PartType = Glu_
    ExtParticle(6)%PartType = Glu_

    IF( Correction.EQ.2 ) THEN
      NumPrimAmps = 24
      NumBornAmps = 24
    ENDIF
    allocate(PrimAmps(1:NumPrimAmps))
    allocate(BornAmps(1:NumPrimAmps))
    do NAmp=1,NumPrimAmps
        allocate(BornAmps(NAmp)%ExtLine(1:NumExtParticles))
        allocate(PrimAmps(NAmp)%ExtLine(1:NumExtParticles))
        allocate(PrimAmps(NAmp)%IntPart(1:NumExtParticles))
    enddo

    IF( TOPDECAYS.GE.1 ) THEN
    NumHelicities = 16
    allocate(Helicities(1:NumHelicities,1:NumExtParticles))
      ih=1
      do h3=-1,1,2
      do h4=-1,1,2
      do h5=-1,1,2
      do h6=-1,1,2
          if( ih.ge.17 ) cycle
          Helicities(ih,1:6) = (/0,0,h3,h4,h5,h6/)
          ih=ih+1
      enddo
      enddo
      enddo
      enddo
    ELSEIF( TOPDECAYS.EQ.0 ) then
        NumHelicities = 32  ! uses the parity flip
        allocate(Helicities(1:NumHelicities,1:NumExtParticles))
          ih=1
          do h1=-1,1,2
          do h2=-1,1,2
          do h3=-1,1,2
          do h4=-1,1,2
          do h5=-1,1,2
          do h6=-1,1,2
              if( ih.ge.33 ) cycle
              Helicities(ih,1:6) = (/h1,h2,h3,h4,h5,h6/)
              ih=ih+1
          enddo
          enddo
          enddo
          enddo
          enddo
          enddo
    ENDIF


ELSEIF( MASTERPROCESS.EQ.6 ) THEN

    ExtParticle(1)%PartType = ATop_
    ExtParticle(2)%PartType = Top_
    ExtParticle(3)%PartType = AStr_
    ExtParticle(4)%PartType = Str_
    ExtParticle(5)%PartType = Glu_
    ExtParticle(6)%PartType = Glu_

    IF( Correction.EQ.2 ) THEN
       NumPrimAmps = 12
       NumBornAmps = 12
    ENDIF
    allocate(PrimAmps(1:NumPrimAmps))
    allocate(BornAmps(1:NumPrimAmps))
    do NAmp=1,NumPrimAmps
        allocate(BornAmps(NAmp)%ExtLine(1:NumExtParticles))
        allocate(PrimAmps(NAmp)%ExtLine(1:NumExtParticles))
        allocate(PrimAmps(NAmp)%IntPart(1:NumExtParticles))
    enddo

    IF( TOPDECAYS.GE.1 ) THEN
    NumHelicities = 16
    allocate(Helicities(1:NumHelicities,1:NumExtParticles))
      ih=1
      do h3=-1,1,2
      do h4=-1,1,2
      do h5=-1,1,2
      do h6=-1,1,2
          if( ih.ge.17 ) cycle
          Helicities(ih,1:6) = (/0,0,h3,h4,h5,h6/)
          ih=ih+1
      enddo
      enddo
      enddo
      enddo
    ELSEIF( TOPDECAYS.EQ.0 ) then
    NumHelicities = 64
    allocate(Helicities(1:NumHelicities,1:NumExtParticles))
      ih=1
      do h1=-1,1,2
      do h2=-1,1,2
      do h3=-1,1,2
      do h4=-1,1,2
      do h5=-1,1,2
      do h6=-1,1,2
          if( ih.ge.65 ) cycle
          Helicities(ih,1:6) = (/h1,h2,h3,h4,h5,h6/)
          ih=ih+1
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
    ENDIF

ELSEIF( MASTERPROCESS.EQ.7 ) THEN

    ExtParticle(1)%PartType = ATop_
    ExtParticle(2)%PartType = Top_
    ExtParticle(3)%PartType = AStr_
    ExtParticle(4)%PartType = Str_
    ExtParticle(5)%PartType = AChm_
    ExtParticle(6)%PartType = Chm_


ELSEIF( MASTERPROCESS.EQ.8 ) THEN

    ExtParticle(1)%PartType = ATop_
    ExtParticle(2)%PartType = Top_
    ExtParticle(3)%PartType = Glu_
    ExtParticle(4)%PartType = Glu_
    ExtParticle(5)%PartType = Glu_  ! this is the photon!

    IF( Correction.EQ.0 .OR. Correction.EQ.4 .OR.Correction.EQ.5 ) THEN
      NumPrimAmps = 2
      NumBornAmps = 2
    ELSEIF( Correction.EQ.1 ) THEN
      NumPrimAmps = 28
      NumBornAmps = 2
    ENDIF
    allocate(PrimAmps(1:NumPrimAmps))
    allocate(BornAmps(1:NumPrimAmps))
    do NAmp=1,NumPrimAmps
        allocate(BornAmps(NAmp)%ExtLine(1:NumExtParticles))
        allocate(PrimAmps(NAmp)%ExtLine(1:NumExtParticles))
        allocate(PrimAmps(NAmp)%IntPart(1:NumExtParticles))
    enddo

    IF( TOPDECAYS.GE.1 ) THEN
        NumHelicities = 8
        allocate(Helicities(1:NumHelicities,1:NumExtParticles))
        Helicities(1,1:5) = (/0,0,+1,+1,+1/)
        Helicities(2,1:5) = (/0,0,+1,+1,-1/)
        Helicities(3,1:5) = (/0,0,+1,-1,+1/)
        Helicities(4,1:5) = (/0,0,+1,-1,-1/)
        Helicities(5,1:5) = (/0,0,-1,+1,+1/)
        Helicities(6,1:5) = (/0,0,-1,+1,-1/)
        Helicities(7,1:5) = (/0,0,-1,-1,+1/)
        Helicities(8,1:5) = (/0,0,-1,-1,-1/)
    ELSE
        NumHelicities = 32
        allocate(Helicities(1:NumHelicities,1:NumExtParticles))
        sig_tb=+1; sig_t =+1;
        Helicities( 1,1:5) = (/sig_tb,sig_t,+1,+1,+1/)
        Helicities( 2,1:5) = (/sig_tb,sig_t,+1,+1,-1/)
        Helicities( 3,1:5) = (/sig_tb,sig_t,+1,-1,-1/)
        Helicities( 4,1:5) = (/sig_tb,sig_t,+1,-1,+1/)
        sig_tb=+1; sig_t =-1;
        Helicities( 5,1:5) = (/sig_tb,sig_t,+1,+1,+1/)
        Helicities( 6,1:5) = (/sig_tb,sig_t,+1,+1,-1/)
        Helicities( 7,1:5) = (/sig_tb,sig_t,+1,-1,-1/)
        Helicities( 8,1:5) = (/sig_tb,sig_t,+1,-1,+1/)
        sig_tb=-1; sig_t =+1;
        Helicities( 9,1:5) = (/sig_tb,sig_t,+1,+1,+1/)
        Helicities(10,1:5) = (/sig_tb,sig_t,+1,+1,-1/)
        Helicities(11,1:5) = (/sig_tb,sig_t,+1,-1,-1/)
        Helicities(12,1:5) = (/sig_tb,sig_t,+1,-1,+1/)
        sig_tb=-1; sig_t =-1;
        Helicities(13,1:5) = (/sig_tb,sig_t,+1,+1,+1/)
        Helicities(14,1:5) = (/sig_tb,sig_t,+1,+1,-1/)
        Helicities(15,1:5) = (/sig_tb,sig_t,+1,-1,-1/)
        Helicities(16,1:5) = (/sig_tb,sig_t,+1,-1,+1/)
    !   additional helicities when parity inversion is not applied:
        sig_tb=-1; sig_t =-1;
        Helicities(17,1:5) = (/sig_tb,sig_t,-1,-1,-1/)
        Helicities(18,1:5) = (/sig_tb,sig_t,-1,-1,+1/)
        Helicities(19,1:5) = (/sig_tb,sig_t,-1,+1,+1/)
        Helicities(20,1:5) = (/sig_tb,sig_t,-1,+1,-1/)
        sig_tb=-1; sig_t =+1;
        Helicities(21,1:5) = (/sig_tb,sig_t,-1,-1,-1/)
        Helicities(22,1:5) = (/sig_tb,sig_t,-1,-1,+1/)
        Helicities(23,1:5) = (/sig_tb,sig_t,-1,+1,+1/)
        Helicities(24,1:5) = (/sig_tb,sig_t,-1,+1,-1/)
        sig_tb=+1; sig_t =-1;
        Helicities(25,1:5) = (/sig_tb,sig_t,-1,-1,-1/)
        Helicities(26,1:5) = (/sig_tb,sig_t,-1,-1,+1/)
        Helicities(27,1:5) = (/sig_tb,sig_t,-1,+1,+1/)
        Helicities(28,1:5) = (/sig_tb,sig_t,-1,+1,-1/)
        sig_tb=+1; sig_t =+1;
        Helicities(29,1:5) = (/sig_tb,sig_t,-1,-1,-1/)
        Helicities(30,1:5) = (/sig_tb,sig_t,-1,-1,+1/)
        Helicities(31,1:5) = (/sig_tb,sig_t,-1,+1,+1/)
        Helicities(32,1:5) = (/sig_tb,sig_t,-1,+1,-1/)
    ENDIF


ELSEIF( MASTERPROCESS.EQ.9 ) THEN

    ExtParticle(1)%PartType = ATop_
    ExtParticle(2)%PartType = Top_
    ExtParticle(3)%PartType = AStr_
    ExtParticle(4)%PartType = Str_
    ExtParticle(5)%PartType = Glu_  ! this is the photon!

    IF( Correction.EQ.0 .OR. Correction.EQ.4 .OR.Correction.EQ.5) THEN
      NumPrimAmps = 2
      NumBornAmps = 2
    ELSEIF( Correction.EQ.1 ) THEN
      NumPrimAmps = 14
      NumBornAmps = 2
    ENDIF
    allocate(PrimAmps(1:NumPrimAmps))
    allocate(BornAmps(1:NumPrimAmps))
    do NAmp=1,NumPrimAmps
        allocate(BornAmps(NAmp)%ExtLine(1:NumExtParticles))
        allocate(PrimAmps(NAmp)%ExtLine(1:NumExtParticles))
        allocate(PrimAmps(NAmp)%IntPart(1:NumExtParticles))
    enddo

    IF( TOPDECAYS.GE.1 ) THEN
        NumHelicities = 8
        allocate(Helicities(1:NumHelicities,1:NumExtParticles))
        Helicities(1,1:5) = (/0,0,+1,+1,+1/)
        Helicities(2,1:5) = (/0,0,+1,+1,-1/)
        Helicities(3,1:5) = (/0,0,+1,-1,+1/)
        Helicities(4,1:5) = (/0,0,+1,-1,-1/)
        Helicities(5,1:5) = (/0,0,-1,+1,+1/)
        Helicities(6,1:5) = (/0,0,-1,+1,-1/)
        Helicities(7,1:5) = (/0,0,-1,-1,+1/)
        Helicities(8,1:5) = (/0,0,-1,-1,-1/)
    ELSE
        NumHelicities = 32
        allocate(Helicities(1:NumHelicities,1:NumExtParticles))
        sig_tb=+1; sig_t =+1;
        Helicities( 1,1:5) = (/sig_tb,sig_t,+1,+1,+1/)
        Helicities( 2,1:5) = (/sig_tb,sig_t,+1,+1,-1/)
        Helicities( 3,1:5) = (/sig_tb,sig_t,+1,-1,-1/)
        Helicities( 4,1:5) = (/sig_tb,sig_t,+1,-1,+1/)
        sig_tb=+1; sig_t =-1;
        Helicities( 5,1:5) = (/sig_tb,sig_t,+1,+1,+1/)
        Helicities( 6,1:5) = (/sig_tb,sig_t,+1,+1,-1/)
        Helicities( 7,1:5) = (/sig_tb,sig_t,+1,-1,-1/)
        Helicities( 8,1:5) = (/sig_tb,sig_t,+1,-1,+1/)
        sig_tb=-1; sig_t =+1;
        Helicities( 9,1:5) = (/sig_tb,sig_t,+1,+1,+1/)
        Helicities(10,1:5) = (/sig_tb,sig_t,+1,+1,-1/)
        Helicities(11,1:5) = (/sig_tb,sig_t,+1,-1,-1/)
        Helicities(12,1:5) = (/sig_tb,sig_t,+1,-1,+1/)
        sig_tb=-1; sig_t =-1;
        Helicities(13,1:5) = (/sig_tb,sig_t,+1,+1,+1/)
        Helicities(14,1:5) = (/sig_tb,sig_t,+1,+1,-1/)
        Helicities(15,1:5) = (/sig_tb,sig_t,+1,-1,-1/)
        Helicities(16,1:5) = (/sig_tb,sig_t,+1,-1,+1/)
    !   additional helicities when parity inversion is not applied:
        sig_tb=-1; sig_t =-1;
        Helicities(17,1:5) = (/sig_tb,sig_t,-1,-1,-1/)
        Helicities(18,1:5) = (/sig_tb,sig_t,-1,-1,+1/)
        Helicities(19,1:5) = (/sig_tb,sig_t,-1,+1,+1/)
        Helicities(20,1:5) = (/sig_tb,sig_t,-1,+1,-1/)
        sig_tb=-1; sig_t =+1;
        Helicities(21,1:5) = (/sig_tb,sig_t,-1,-1,-1/)
        Helicities(22,1:5) = (/sig_tb,sig_t,-1,-1,+1/)
        Helicities(23,1:5) = (/sig_tb,sig_t,-1,+1,+1/)
        Helicities(24,1:5) = (/sig_tb,sig_t,-1,+1,-1/)
        sig_tb=+1; sig_t =-1;
        Helicities(25,1:5) = (/sig_tb,sig_t,-1,-1,-1/)
        Helicities(26,1:5) = (/sig_tb,sig_t,-1,-1,+1/)
        Helicities(27,1:5) = (/sig_tb,sig_t,-1,+1,+1/)
        Helicities(28,1:5) = (/sig_tb,sig_t,-1,+1,-1/)
        sig_tb=+1; sig_t =+1;
        Helicities(29,1:5) = (/sig_tb,sig_t,-1,-1,-1/)
        Helicities(30,1:5) = (/sig_tb,sig_t,-1,-1,+1/)
        Helicities(31,1:5) = (/sig_tb,sig_t,-1,+1,+1/)
        Helicities(32,1:5) = (/sig_tb,sig_t,-1,+1,-1/)
    ENDIF




ELSEIF( MASTERPROCESS.EQ.10 ) THEN

    ExtParticle(1)%PartType = ATop_
    ExtParticle(2)%PartType = Top_
    ExtParticle(3)%PartType = Glu_
    ExtParticle(4)%PartType = Glu_
    ExtParticle(5)%PartType = Glu_
    ExtParticle(6)%PartType = Glu_  ! this is the photon!

    IF( Correction.EQ.2 ) THEN
      NumPrimAmps = 6
      NumBornAmps = 6
    ENDIF
    allocate(PrimAmps(1:NumPrimAmps))
    allocate(BornAmps(1:NumPrimAmps))
    do NAmp=1,NumPrimAmps
        allocate(BornAmps(NAmp)%ExtLine(1:NumExtParticles))
        allocate(PrimAmps(NAmp)%ExtLine(1:NumExtParticles))
        allocate(PrimAmps(NAmp)%IntPart(1:NumExtParticles))
    enddo



    IF( TOPDECAYS.GE.1 ) THEN
    NumHelicities = 16
    allocate(Helicities(1:NumHelicities,1:NumExtParticles))
      ih=1
      do h3=-1,1,2
      do h4=-1,1,2
      do h5=-1,1,2
      do h6=-1,1,2
          if( ih.ge.17 ) cycle
          Helicities(ih,1:6) = (/0,0,h3,h4,h5,h6/)
          ih=ih+1
      enddo
      enddo
      enddo
      enddo
    ELSEIF( TOPDECAYS.EQ.0 ) then
    NumHelicities = 64
    allocate(Helicities(1:NumHelicities,1:NumExtParticles))
      ih=1
      do h1=-1,1,2
      do h2=-1,1,2
      do h3=-1,1,2
      do h4=-1,1,2
      do h5=-1,1,2
      do h6=-1,1,2
          if( ih.ge.65 ) cycle
          Helicities(ih,1:6) = (/h1,h2,h3,h4,h5,h6/)
          ih=ih+1
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
    ENDIF






ELSEIF( MASTERPROCESS.EQ.11 ) THEN

    ExtParticle(1)%PartType = ATop_
    ExtParticle(2)%PartType = Top_
    ExtParticle(3)%PartType = AChm_
    ExtParticle(4)%PartType = Chm_
    ExtParticle(5)%PartType = Glu_
    ExtParticle(6)%PartType = Glu_  ! this is the photon!

    IF( Correction.EQ.2 ) THEN
      NumPrimAmps = 10
      NumBornAmps = 10
    ENDIF
    allocate(PrimAmps(1:NumPrimAmps))
    allocate(BornAmps(1:NumPrimAmps))
    do NAmp=1,NumPrimAmps
        allocate(BornAmps(NAmp)%ExtLine(1:NumExtParticles))
        allocate(PrimAmps(NAmp)%ExtLine(1:NumExtParticles))
        allocate(PrimAmps(NAmp)%IntPart(1:NumExtParticles))
    enddo

    IF( TOPDECAYS.GE.1 ) THEN
    NumHelicities = 16
    allocate(Helicities(1:NumHelicities,1:NumExtParticles))
      ih=1
      do h3=-1,1,2
      do h4=-1,1,2
      do h5=-1,1,2
      do h6=-1,1,2
          if( ih.ge.17 ) cycle
          Helicities(ih,1:6) = (/0,0,h3,h4,h5,h6/)
          ih=ih+1
      enddo
      enddo
      enddo
      enddo
    ELSEIF( TOPDECAYS.EQ.0 ) then
    NumHelicities = 64
    allocate(Helicities(1:NumHelicities,1:NumExtParticles))
      ih=1
      do h1=-1,1,2
      do h2=-1,1,2
      do h3=-1,1,2
      do h4=-1,1,2
      do h5=-1,1,2
      do h6=-1,1,2
          if( ih.ge.65 ) cycle
          Helicities(ih,1:6) = (/h1,h2,h3,h4,h5,h6/)
          ih=ih+1
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
    ENDIF





ELSE
    call Error("MasterProcess not implemented in InitMasterProcess")

ENDIF


   do NPart=1,NumExtParticles
         ExtParticle(NPart)%ExtRef = NPart
         ExtParticle(NPart)%Mass = GetMass( ExtParticle(NPart)%PartType )
         ExtParticle(NPart)%Mass2 = (ExtParticle(NPart)%Mass)**2
         ExtParticle(NPart)%Mom(1:8) = (0d0,0d0)
   enddo
END SUBROUTINE








SUBROUTINE InitAmps()
use ModMisc
use ModParameters
implicit none
integer :: Vertex,Propa,PropaMinus1,ExtPartType,NPrimAmp,k
integer :: AllocStatus,counter,counterQ,counterG,QuarkPos(1:6),NPart
logical :: ColorLessParticles
type(PrimitiveAmplitude),pointer :: ThePrimAmp
type(BornAmplitude),pointer :: TheBornAmp
type(TreeProcess),pointer :: TheTree


! convention for labeling particles with 4 quarks: (1_tb, 2_t, 3_qb, 4_q)
! AmpType = 1/a:  1,2,3,4 (both fermion lines in loop)
! AmpType = 3/b:  1,4,3,2 (top line in loop)
! AmpType = 4/c:  1,2,3,4 (light fermion line in loop)
! AmpType = 2/d:  1,2,3,4 (closed fermion loop)
!  the ordering of the first PrimAmps should match the BornAmps because of the MassCT contribution



IF( MASTERPROCESS.EQ.0 ) THEN
   BornAmps(1)%ExtLine = (/1,2,3,4,5,6,7/)
   PrimAmps(1)%ExtLine = (/1,2,3,4,5,6,7/)
   PrimAmps(1)%AmpType = 1



ELSEIF( MASTERPROCESS.EQ.1 ) THEN
   IF ( Correction.EQ.1 ) THEN
   BornAmps(1)%ExtLine = (/1,2,3,4/)
   BornAmps(2)%ExtLine = (/1,2,4,3/)

   PrimAmps(1)%ExtLine = (/1,2,3,4/)
   PrimAmp1_1234 = 2
   PrimAmps(1)%AmpType = 1

   PrimAmps(2)%ExtLine = (/1,2,4,3/)
   PrimAmp1_1243 = 1
   PrimAmps(2)%AmpType = 1

   PrimAmps(3)%ExtLine = (/1,3,2,4/)
   PrimAmp1_1324 = 3
   PrimAmps(3)%AmpType = 1

   PrimAmps(4)%ExtLine = (/1,4,2,3/)
   PrimAmp1_1423 = 4
   PrimAmps(4)%AmpType = 1

   PrimAmps(5)%ExtLine = (/1,3,4,2/)
   PrimAmp1_1342 = 6
   PrimAmps(5)%AmpType = 1

   PrimAmps(6)%ExtLine = (/1,4,3,2/)
   PrimAmp1_1432 = 5
   PrimAmps(6)%AmpType = 1

   PrimAmps(7)%ExtLine = (/1,2,3,4/)
   PrimAmp2_1234 = 1
   PrimAmps(7)%AmpType = 2
   PrimAmps(7)%FermLoopPart = Chm_

   PrimAmps(8)%ExtLine = (/1,2,4,3/)
   PrimAmp2_1243 = 2
   PrimAmps(8)%AmpType = 2
   PrimAmps(8)%FermLoopPart = Chm_

   PrimAmps(9)%ExtLine = (/1,2,3,4/)
   PrimAmp2m_1234 = 1
   PrimAmps(9)%AmpType = 2
   PrimAmps(9)%FermLoopPart = Bot_

   PrimAmps(10)%ExtLine = (/1,2,4,3/)
   PrimAmp2m_1243 = 2
   PrimAmps(10)%AmpType = 2
   PrimAmps(10)%FermLoopPart = Bot_

   ELSEIF( Correction.EQ.0 .OR. Correction.GE.4 ) THEN
   BornAmps(1)%ExtLine = (/1,2,3,4/)
   BornAmps(2)%ExtLine = (/1,2,4,3/)

   PrimAmps(1)%ExtLine = (/1,2,3,4/)
   PrimAmps(2)%ExtLine = (/1,2,4,3/)
   ENDIF


ELSEIF( MasterProcess.EQ.2) THEN
   IF( Correction.EQ.1 ) THEN
   BornAmps(1)%ExtLine = (/1,2,3,4/)
   BornAmps(2)%ExtLine = (/1,2,3,4/)
   BornAmps(3)%ExtLine = (/1,4,3,2/)
   BornAmps(4)%ExtLine = (/1,2,3,4/)
   BornAmps(5)%ExtLine = (/1,2,3,4/)
   BornAmps(6)%ExtLine = (/1,2,3,4/)

   PrimAmps(1)%ExtLine = (/1,2,3,4/)
   PrimAmp1_1234 = 1
   PrimAmps(1)%AmpType = 1

   PrimAmps(2)%ExtLine = (/1,2,4,3/)
   PrimAmp1_1243 = 2
   PrimAmps(2)%AmpType = 1

!    PrimAmps(3)%ExtLine = (/1,3,4,2/)
   PrimAmps(3)%ExtLine = (/1,4,3,2/)
   PrimAmp3_1432 = 3
   PrimAmps(3)%AmpType = 3

!    PrimAmps(4)%ExtLine = (/1,2,4,3/)
   PrimAmps(4)%ExtLine = (/1,2,3,4/)
   PrimAmp4_1234 = 4
   PrimAmps(4)%AmpType = 4

   PrimAmps(5)%ExtLine = (/1,2,3,4/)
   PrimAmp2_1234 = 5
   PrimAmps(5)%AmpType = 2
   PrimAmps(5)%FermLoopPart = Chm_

   PrimAmps(6)%ExtLine = (/1,2,3,4/)
   PrimAmp2m_1234 = 6
   PrimAmps(6)%AmpType = 2
   PrimAmps(6)%FermLoopPart = Bot_

   ELSEIF ( Correction.EQ.0 .OR. Correction.GE.4 ) THEN
   BornAmps(1)%ExtLine = (/1,2,3,4/)
   PrimAmps(1)%ExtLine = (/1,2,3,4/)
   ENDIF



ELSEIF( MasterProcess.EQ.3 ) THEN
   IF( Correction.EQ.1 ) THEN
   BornAmps(1)%ExtLine = (/1,2,3,4,5/)
   BornAmps(2)%ExtLine = (/1,2,3,5,4/)
   BornAmps(3)%ExtLine = (/1,2,4,3,5/)
   BornAmps(4)%ExtLine = (/1,2,4,5,3/)
   BornAmps(5)%ExtLine = (/1,2,5,3,4/)
   BornAmps(6)%ExtLine = (/1,2,5,4,3/)

   PrimAmps(1)%ExtLine = (/1,2,3,4,5/)
   PrimAmp1_12345 = 1
   PrimAmps(1)%AmpType = 1

   PrimAmps(2)%ExtLine = (/1,2,3,5,4/)
   PrimAmp1_12354 = 2
   PrimAmps(2)%AmpType = 1

   PrimAmps(3)%ExtLine = (/1,2,4,3,5/)
   PrimAmp1_12435 = 3
   PrimAmps(3)%AmpType = 1

   PrimAmps(4)%ExtLine = (/1,2,4,5,3/)
   PrimAmp1_12453 = 4
   PrimAmps(4)%AmpType = 1

   PrimAmps(5)%ExtLine = (/1,2,5,3,4/)
   PrimAmp1_12534 = 5
   PrimAmps(5)%AmpType = 1

   PrimAmps(6)%ExtLine = (/1,2,5,4,3/)
   PrimAmp1_12543 = 6
   PrimAmps(6)%AmpType = 1

!-------

   PrimAmps(7)%ExtLine = (/1,3,2,4,5/)
   PrimAmp1_13245 = 7
   PrimAmps(7)%AmpType = 1

   PrimAmps(8)%ExtLine = (/1,3,2,5,4/)
   PrimAmp1_13254 = 8
   PrimAmps(8)%AmpType = 1

   PrimAmps(9)%ExtLine = (/1,4,2,3,5/)
   PrimAmp1_14235 = 9
   PrimAmps(9)%AmpType = 1

   PrimAmps(10)%ExtLine = (/1,4,2,5,3/)
   PrimAmp1_14253 = 10
   PrimAmps(10)%AmpType = 1

   PrimAmps(11)%ExtLine = (/1,5,2,3,4/)
   PrimAmp1_15234 = 11
   PrimAmps(11)%AmpType = 1

   PrimAmps(12)%ExtLine = (/1,5,2,4,3/)
   PrimAmp1_15243 = 12
   PrimAmps(12)%AmpType = 1

!-------

   PrimAmps(13)%ExtLine = (/1,3,4,2,5/)
   PrimAmp1_13425 = 13
   PrimAmps(13)%AmpType = 1

   PrimAmps(14)%ExtLine = (/1,3,5,2,4/)
   PrimAmp1_13524 = 14
   PrimAmps(14)%AmpType = 1

   PrimAmps(15)%ExtLine = (/1,4,3,2,5/)
   PrimAmp1_14325 = 15
   PrimAmps(15)%AmpType = 1

   PrimAmps(16)%ExtLine = (/1,4,5,2,3/)
   PrimAmp1_14523 = 16
   PrimAmps(16)%AmpType = 1

   PrimAmps(17)%ExtLine = (/1,5,3,2,4/)
   PrimAmp1_15324 = 17
   PrimAmps(17)%AmpType = 1

   PrimAmps(18)%ExtLine = (/1,5,4,2,3/)
   PrimAmp1_15423 = 18
   PrimAmps(18)%AmpType = 1

!-------

   PrimAmps(19)%ExtLine = (/1,3,4,5,2/)
   PrimAmp1_13452 = 19
   PrimAmps(19)%AmpType = 1

   PrimAmps(20)%ExtLine = (/1,3,5,4,2/)
   PrimAmp1_13542 = 20
   PrimAmps(20)%AmpType = 1

   PrimAmps(21)%ExtLine = (/1,4,3,5,2/)
   PrimAmp1_14352 = 21
   PrimAmps(21)%AmpType = 1

   PrimAmps(22)%ExtLine = (/1,4,5,3,2/)
   PrimAmp1_14532 = 22
   PrimAmps(22)%AmpType = 1

   PrimAmps(23)%ExtLine = (/1,5,3,4,2/)
   PrimAmp1_15342 = 23
   PrimAmps(23)%AmpType = 1

   PrimAmps(24)%ExtLine = (/1,5,4,3,2/)
   PrimAmp1_15432 = 24
   PrimAmps(24)%AmpType = 1


!------- fermion loops (massless)


   PrimAmps(25)%ExtLine = (/1,2,3,4,5/)
!    PrimAmp1_12345 = 1
   PrimAmps(25)%AmpType = 2
   PrimAmps(25)%FermLoopPart = Chm_

   PrimAmps(26)%ExtLine = (/1,2,3,5,4/)
!    PrimAmp1_12354 = 2
   PrimAmps(26)%AmpType = 2
   PrimAmps(26)%FermLoopPart = Chm_

   PrimAmps(27)%ExtLine = (/1,2,4,3,5/)
!    PrimAmp1_12435 = 3
   PrimAmps(27)%AmpType = 2
   PrimAmps(27)%FermLoopPart = Chm_

   PrimAmps(28)%ExtLine = (/1,2,4,5,3/)
!    PrimAmp1_12453 = 4
   PrimAmps(28)%AmpType = 2
   PrimAmps(28)%FermLoopPart = Chm_

   PrimAmps(29)%ExtLine = (/1,2,5,3,4/)
!    PrimAmp1_12534 = 5
   PrimAmps(29)%AmpType = 2
   PrimAmps(29)%FermLoopPart = Chm_

   PrimAmps(30)%ExtLine = (/1,2,5,4,3/)
!    PrimAmp1_12543 = 6
   PrimAmps(30)%AmpType = 2
   PrimAmps(30)%FermLoopPart = Chm_

!-------

   PrimAmps(31)%ExtLine = (/1,3,2,4,5/)
!    PrimAmp1_13245 = 7
   PrimAmps(31)%AmpType = 2
   PrimAmps(31)%FermLoopPart = Chm_

   PrimAmps(32)%ExtLine = (/1,3,2,5,4/)
!    PrimAmp1_13254 = 8
   PrimAmps(32)%AmpType = 2
   PrimAmps(32)%FermLoopPart = Chm_

   PrimAmps(33)%ExtLine = (/1,4,2,3,5/)
!    PrimAmp1_14235 = 9
   PrimAmps(33)%AmpType = 2
   PrimAmps(33)%FermLoopPart = Chm_

   PrimAmps(34)%ExtLine = (/1,4,2,5,3/)
!    PrimAmp1_14253 = 10
   PrimAmps(34)%AmpType = 2
   PrimAmps(34)%FermLoopPart = Chm_

   PrimAmps(35)%ExtLine = (/1,5,2,3,4/)
!    PrimAmp1_15234 = 11
   PrimAmps(35)%AmpType = 2
   PrimAmps(35)%FermLoopPart = Chm_

   PrimAmps(36)%ExtLine = (/1,5,2,4,3/)
!    PrimAmp1_15243 = 12
   PrimAmps(36)%AmpType = 2
   PrimAmps(36)%FermLoopPart = Chm_


!------- fermion loops (massive)


   PrimAmps(37)%ExtLine = (/1,2,3,4,5/)
!    PrimAmp1_12345 = 1
   PrimAmps(37)%AmpType = 2
   PrimAmps(37)%FermLoopPart = Bot_

   PrimAmps(38)%ExtLine = (/1,2,3,5,4/)
!    PrimAmp1_12354 = 2
   PrimAmps(38)%AmpType = 2
   PrimAmps(38)%FermLoopPart = Bot_

   PrimAmps(39)%ExtLine = (/1,2,4,3,5/)
!    PrimAmp1_12435 = 3
   PrimAmps(39)%AmpType = 2
   PrimAmps(39)%FermLoopPart = Bot_

   PrimAmps(40)%ExtLine = (/1,2,4,5,3/)
!    PrimAmp1_12453 = 4
   PrimAmps(40)%AmpType = 2
   PrimAmps(40)%FermLoopPart = Bot_

   PrimAmps(41)%ExtLine = (/1,2,5,3,4/)
!    PrimAmp1_12534 = 5
   PrimAmps(41)%AmpType = 2
   PrimAmps(41)%FermLoopPart = Bot_

   PrimAmps(42)%ExtLine = (/1,2,5,4,3/)
!    PrimAmp1_12543 = 6
   PrimAmps(42)%AmpType = 2
   PrimAmps(42)%FermLoopPart = Bot_

!-------

   PrimAmps(43)%ExtLine = (/1,3,2,4,5/)
!    PrimAmp1_13245 = 7
   PrimAmps(43)%AmpType = 2
   PrimAmps(43)%FermLoopPart = Bot_

   PrimAmps(44)%ExtLine = (/1,3,2,5,4/)
!    PrimAmp1_13254 = 8
   PrimAmps(44)%AmpType = 2
   PrimAmps(44)%FermLoopPart = Bot_

   PrimAmps(45)%ExtLine = (/1,4,2,3,5/)
!    PrimAmp1_14235 = 9
   PrimAmps(45)%AmpType = 2
   PrimAmps(45)%FermLoopPart = Bot_

   PrimAmps(46)%ExtLine = (/1,4,2,5,3/)
!    PrimAmp1_14253 = 10
   PrimAmps(46)%AmpType = 2
   PrimAmps(46)%FermLoopPart = Bot_

   PrimAmps(47)%ExtLine = (/1,5,2,3,4/)
!    PrimAmp1_15234 = 11
   PrimAmps(47)%AmpType = 2
   PrimAmps(47)%FermLoopPart = Bot_

   PrimAmps(48)%ExtLine = (/1,5,2,4,3/)
!    PrimAmp1_15243 = 12
   PrimAmps(48)%AmpType = 2
   PrimAmps(48)%FermLoopPart = Bot_


   ELSEIF( Correction.EQ.0 .OR. Correction.EQ.2 .OR. Correction.EQ.4 ) THEN
   BornAmps(1)%ExtLine = (/1,2,3,4,5/)
   BornAmps(2)%ExtLine = (/1,2,3,5,4/)
   BornAmps(3)%ExtLine = (/1,2,4,3,5/)
   BornAmps(4)%ExtLine = (/1,2,4,5,3/)
   BornAmps(5)%ExtLine = (/1,2,5,3,4/)
   BornAmps(6)%ExtLine = (/1,2,5,4,3/)

   PrimAmps(1)%ExtLine = (/1,2,3,4,5/)
   PrimAmps(2)%ExtLine = (/1,2,3,5,4/)
   PrimAmps(3)%ExtLine = (/1,2,4,3,5/)
   PrimAmps(4)%ExtLine = (/1,2,4,5,3/)
   PrimAmps(5)%ExtLine = (/1,2,5,3,4/)
   PrimAmps(6)%ExtLine = (/1,2,5,4,3/)
   ENDIF



ELSEIF( MasterProcess.EQ.4 ) THEN
   IF( Correction.EQ.0 .OR. Correction.EQ.2 .OR. Correction.EQ.4 ) THEN
   PrimAmps(1)%ExtLine = (/1,2,3,4,5/)
   PrimAmps(2)%ExtLine = (/1,5,2,3,4/)
   PrimAmps(3)%ExtLine = (/1,2,5,3,4/)
   PrimAmps(4)%ExtLine = (/1,2,3,5,4/)

   ELSEIF( Correction.EQ.1 ) THEN
!    BornAmps(1)%ExtLine = (/1,2,3,4,5/)  ! will be overwritten anyways
!    BornAmps(2)%ExtLine = (/1,5,2,3,4/)
!    BornAmps(3)%ExtLine = (/1,2,5,3,4/)
!    BornAmps(4)%ExtLine = (/1,2,3,5,4/)

   PrimAmps(1)%ExtLine = (/1,2,3,4,5/)
   PrimAmp1_12345 = 1
   PrimAmps(1)%AmpType = 1

   PrimAmps(2)%ExtLine = (/1,5,2,3,4/)
   PrimAmp1_15234 = 2
   PrimAmps(2)%AmpType = 1

   PrimAmps(3)%ExtLine = (/1,2,5,3,4/)
   PrimAmp1_12534 = 3
   PrimAmps(3)%AmpType = 1

   PrimAmps(4)%ExtLine = (/1,2,3,5,4/)
   PrimAmp1_12354 = 4
   PrimAmps(4)%AmpType = 1

   PrimAmps(5)%ExtLine = (/1,2,5,4,3/)
   PrimAmp1_12543 = 5
   PrimAmps(5)%AmpType = 1

   PrimAmps(6)%ExtLine = (/1,2,4,3,5/)
   PrimAmp1_12435 = 6
   PrimAmps(6)%AmpType = 1

   PrimAmps(7)%ExtLine = (/1,2,4,5,3/)
   PrimAmp1_12453 = 7
   PrimAmps(7)%AmpType = 1

   PrimAmps(8)%ExtLine = (/1,5,2,4,3/)
   PrimAmp1_15243 = 8
   PrimAmps(8)%AmpType = 1



   PrimAmps(9)%ExtLine = (/1,5,4,3,2/)
   PrimAmp3_15432 = 9
   PrimAmps(9)%AmpType = 3

   PrimAmps(10)%ExtLine = (/1,4,5,3,2/)
   PrimAmp3_14532 = 10
   PrimAmps(10)%AmpType = 3

   PrimAmps(11)%ExtLine = (/1,4,3,5,2/)
   PrimAmp3_14352 = 11
   PrimAmps(11)%AmpType = 3

   PrimAmps(12)%ExtLine = (/1,4,3,2,5/)
   PrimAmp3_14325 = 12
   PrimAmps(12)%AmpType = 3



   PrimAmps(13)%ExtLine = (/1,2,3,5,4/)
   PrimAmp4_12354 = 13
   PrimAmps(13)%AmpType = 4

   PrimAmps(14)%ExtLine = (/1,2,5,3,4/)
   PrimAmp4_12534 = 14
   PrimAmps(14)%AmpType = 4

   PrimAmps(15)%ExtLine = (/1,5,2,3,4/)
   PrimAmp4_15234 = 15
   PrimAmps(15)%AmpType = 4

   PrimAmps(16)%ExtLine = (/1,2,3,4,5/)
   PrimAmp4_12345 = 16
   PrimAmps(16)%AmpType = 4


!------- fermion loops (massless)

   PrimAmps(17)%ExtLine = (/1,2,3,4,5/)
   PrimAmp2_12345 = 17
   PrimAmps(17)%AmpType = 2
   PrimAmps(17)%FermLoopPart = Chm_

   PrimAmps(18)%ExtLine = (/1,5,2,3,4/)
   PrimAmp2_15234 = 18
   PrimAmps(18)%AmpType = 2
   PrimAmps(18)%FermLoopPart = Chm_

   PrimAmps(19)%ExtLine = (/1,2,5,3,4/)
   PrimAmp2_12534 = 19
   PrimAmps(19)%AmpType = 2
   PrimAmps(19)%FermLoopPart = Chm_

   PrimAmps(20)%ExtLine = (/1,2,3,5,4/)
   PrimAmp2_12354 = 20
   PrimAmps(20)%AmpType = 2
   PrimAmps(20)%FermLoopPart = Chm_

!------- fermion loops (massive)

   PrimAmps(21)%ExtLine = (/1,2,3,4,5/)
   PrimAmp2m_12345 = 21
   PrimAmps(21)%AmpType = 2
   PrimAmps(21)%FermLoopPart = Bot_

   PrimAmps(22)%ExtLine = (/1,5,2,3,4/)
   PrimAmp2m_15234 = 22
   PrimAmps(22)%AmpType = 2
   PrimAmps(22)%FermLoopPart = Bot_

   PrimAmps(23)%ExtLine = (/1,2,5,3,4/)
   PrimAmp2m_12534 = 23
   PrimAmps(23)%AmpType = 2
   PrimAmps(23)%FermLoopPart = Bot_

   PrimAmps(24)%ExtLine = (/1,2,3,5,4/)
   PrimAmp2m_12354 = 24
   PrimAmps(24)%AmpType = 2
   PrimAmps(24)%FermLoopPart = Bot_
   ENDIF

ELSEIF( MasterProcess.EQ.5 ) THEN
   IF( Correction.EQ.2 ) THEN
   BornAmps( 1)%ExtLine = (/1,2,3, 4, 5, 6/)
   BornAmps( 2)%ExtLine = (/1,2,3, 4, 6, 5/)
   BornAmps( 3)%ExtLine = (/1,2,3, 5, 4, 6/)
   BornAmps( 4)%ExtLine = (/1,2,3, 5, 6, 4/)
   BornAmps( 5)%ExtLine = (/1,2,3, 6, 4, 5/)
   BornAmps( 6)%ExtLine = (/1,2,3, 6, 5, 4/)
   BornAmps( 7)%ExtLine = (/1,2,4, 3, 5, 6/)
   BornAmps( 8)%ExtLine = (/1,2,4, 3, 6, 5/)
   BornAmps( 9)%ExtLine = (/1,2,4, 5, 3, 6/)
   BornAmps(10)%ExtLine = (/1,2,4, 5, 6, 3/)
   BornAmps(11)%ExtLine = (/1,2,4, 6, 3, 5/)
   BornAmps(12)%ExtLine = (/1,2,4, 6, 5, 3/)
   BornAmps(13)%ExtLine = (/1,2,5, 3, 4, 6/)
   BornAmps(14)%ExtLine = (/1,2,5, 3, 6, 4/)
   BornAmps(15)%ExtLine = (/1,2,5, 4, 3, 6/)
   BornAmps(16)%ExtLine = (/1,2,5, 4, 6, 3/)
   BornAmps(17)%ExtLine = (/1,2,5, 6, 3, 4/)
   BornAmps(18)%ExtLine = (/1,2,5, 6, 4, 3/)
   BornAmps(19)%ExtLine = (/1,2,6, 3, 4, 5/)
   BornAmps(20)%ExtLine = (/1,2,6, 3, 5, 4/)
   BornAmps(21)%ExtLine = (/1,2,6, 4, 3, 5/)
   BornAmps(22)%ExtLine = (/1,2,6, 4, 5, 3/)
   BornAmps(23)%ExtLine = (/1,2,6, 5, 3, 4/)
   BornAmps(24)%ExtLine = (/1,2,6, 5, 4, 3/)


   PrimAmps( 1)%ExtLine = (/1,2,3, 4, 5, 6/)
   PrimAmps( 2)%ExtLine = (/1,2,3, 4, 6, 5/)
   PrimAmps( 3)%ExtLine = (/1,2,3, 5, 4, 6/)
   PrimAmps( 4)%ExtLine = (/1,2,3, 5, 6, 4/)
   PrimAmps( 5)%ExtLine = (/1,2,3, 6, 4, 5/)
   PrimAmps( 6)%ExtLine = (/1,2,3, 6, 5, 4/)
   PrimAmps( 7)%ExtLine = (/1,2,4, 3, 5, 6/)
   PrimAmps( 8)%ExtLine = (/1,2,4, 3, 6, 5/)
   PrimAmps( 9)%ExtLine = (/1,2,4, 5, 3, 6/)
   PrimAmps(10)%ExtLine = (/1,2,4, 5, 6, 3/)
   PrimAmps(11)%ExtLine = (/1,2,4, 6, 3, 5/)
   PrimAmps(12)%ExtLine = (/1,2,4, 6, 5, 3/)
   PrimAmps(13)%ExtLine = (/1,2,5, 3, 4, 6/)
   PrimAmps(14)%ExtLine = (/1,2,5, 3, 6, 4/)
   PrimAmps(15)%ExtLine = (/1,2,5, 4, 3, 6/)
   PrimAmps(16)%ExtLine = (/1,2,5, 4, 6, 3/)
   PrimAmps(17)%ExtLine = (/1,2,5, 6, 3, 4/)
   PrimAmps(18)%ExtLine = (/1,2,5, 6, 4, 3/)
   PrimAmps(19)%ExtLine = (/1,2,6, 3, 4, 5/)
   PrimAmps(20)%ExtLine = (/1,2,6, 3, 5, 4/)
   PrimAmps(21)%ExtLine = (/1,2,6, 4, 3, 5/)
   PrimAmps(22)%ExtLine = (/1,2,6, 4, 5, 3/)
   PrimAmps(23)%ExtLine = (/1,2,6, 5, 3, 4/)
   PrimAmps(24)%ExtLine = (/1,2,6, 5, 4, 3/)
   ENDIF

ELSEIF( MasterProcess.EQ.6 ) THEN
   IF( Correction.EQ.2 ) THEN
   BornAmps( 1)%ExtLine = (/1,2,3,4,5,6/)
   BornAmps( 2)%ExtLine = (/1,2,3,4,6,5/)
   BornAmps( 3)%ExtLine = (/1,2,5,6,3,4/)
   BornAmps( 4)%ExtLine = (/1,2,6,5,3,4/)
   BornAmps( 5)%ExtLine = (/1,2,5,3,4,6/)
   BornAmps( 6)%ExtLine = (/1,2,6,3,4,5/)
   BornAmps( 7)%ExtLine = (/1,5,6,2,3,4/)
   BornAmps( 8)%ExtLine = (/1,6,5,2,3,4/)
   BornAmps( 9)%ExtLine = (/1,5,2,3,6,4/)
   BornAmps(10)%ExtLine = (/1,6,2,3,5,4/)
   BornAmps(11)%ExtLine = (/1,2,3,5,6,4/)
   BornAmps(12)%ExtLine = (/1,2,3,6,5,4/)

   PrimAmps( 1)%ExtLine = (/1,2,3,4,5,6/)
   PrimAmps( 2)%ExtLine = (/1,2,3,4,6,5/)
   PrimAmps( 3)%ExtLine = (/1,2,5,6,3,4/)
   PrimAmps( 4)%ExtLine = (/1,2,6,5,3,4/)
   PrimAmps( 5)%ExtLine = (/1,2,5,3,4,6/)
   PrimAmps( 6)%ExtLine = (/1,2,6,3,4,5/)
   PrimAmps( 7)%ExtLine = (/1,5,6,2,3,4/)
   PrimAmps( 8)%ExtLine = (/1,6,5,2,3,4/)
   PrimAmps( 9)%ExtLine = (/1,5,2,3,6,4/)
   PrimAmps(10)%ExtLine = (/1,6,2,3,5,4/)
   PrimAmps(11)%ExtLine = (/1,2,3,5,6,4/)
   PrimAmps(12)%ExtLine = (/1,2,3,6,5,4/)
   ENDIF

ELSEIF( MasterProcess.EQ.7 ) THEN
! there's nothing to do here

ELSEIF( MasterProcess.EQ.8 ) THEN! tb t g g pho

   IF( Correction.EQ.0 .OR. Correction.EQ.4 .OR.Correction.EQ.5 ) THEN
      BornAmps(1)%ExtLine = (/1,5,2,3,4/)
      BornAmps(2)%ExtLine = (/1,5,2,4,3/)

      PrimAmps(1)%ExtLine = (/1,5,2,3,4/)
      PrimAmps(2)%ExtLine = (/1,5,2,4,3/)
   ELSEIF( Correction.EQ.1 ) THEN
      BornAmps(1)%ExtLine = (/1,5,2,3,4/)
      BornAmps(2)%ExtLine = (/1,5,2,4,3/)

      PrimAmps(1)%ExtLine = (/1,5,2,3,4/)
      PrimAmp1_15234 = 1
      PrimAmps(1)%AmpType = 1

      PrimAmps(2)%ExtLine = (/1,5,2,4,3/)
      PrimAmp1_15243 = 2
      PrimAmps(2)%AmpType = 1

      PrimAmps(3)%ExtLine = (/1,3,5,4,2/)
      PrimAmp1_13542 = 3
      PrimAmps(3)%AmpType = 1

      PrimAmps(4)%ExtLine = (/1,3,4,5,2/)
      PrimAmp1_13452 = 4
      PrimAmps(4)%AmpType = 1

      PrimAmps(5)%ExtLine = (/1,5,3,4,2/)
      PrimAmp1_15342 = 5
      PrimAmps(5)%AmpType = 1

      PrimAmps(6)%ExtLine = (/1,5,4,3,2/)
      PrimAmp1_15432 = 6
      PrimAmps(6)%AmpType = 1

      PrimAmps(7)%ExtLine = (/1,4,5,3,2/)
      PrimAmp1_14532 = 7
      PrimAmps(7)%AmpType = 1

      PrimAmps(8)%ExtLine = (/1,4,3,5,2/)
      PrimAmp1_14352 = 8
      PrimAmps(8)%AmpType = 1

      PrimAmps(9)%ExtLine = (/1,5,3,2,4/)
      PrimAmp1_15324 = 9
      PrimAmps(9)%AmpType = 1

      PrimAmps(10)%ExtLine = (/1,3,5,2,4/)
      PrimAmp1_13524 = 10
      PrimAmps(10)%AmpType = 1

      PrimAmps(11)%ExtLine = (/1,5,4,2,3/)
      PrimAmp1_15423 = 11
      PrimAmps(11)%AmpType = 1

      PrimAmps(12)%ExtLine = (/1,4,5,2,3/)
      PrimAmp1_14523 = 12
      PrimAmps(12)%AmpType = 1




      PrimAmps(13)%ExtLine = (/1,5,2,3,4/)
      PrimAmp2_15234 = 13
      PrimAmps(13)%AmpType = 2
      PrimAmps(13)%FermLoopPart = Chm_

      PrimAmps(14)%ExtLine = (/1,2,3,4,5/)
      PrimAmp2_12345 = 14
      PrimAmps(14)%AmpType = 2
      PrimAmps(14)%FermLoopPart = Chm_

      PrimAmps(15)%ExtLine = (/1,2,3,5,4/)
      PrimAmp2_12354 = 15
      PrimAmps(15)%AmpType = 2
      PrimAmps(15)%FermLoopPart = Chm_

      PrimAmps(16)%ExtLine = (/1,2,5,3,4/)
      PrimAmp2_12534 = 16
      PrimAmps(16)%AmpType = 2
      PrimAmps(16)%FermLoopPart = Chm_

      PrimAmps(17)%ExtLine = (/1,5,2,4,3/)
      PrimAmp2_15243 = 17
      PrimAmps(17)%AmpType = 2
      PrimAmps(17)%FermLoopPart = Chm_

      PrimAmps(18)%ExtLine = (/1,2,5,4,3/)
      PrimAmp2_12543 = 18
      PrimAmps(18)%AmpType = 2
      PrimAmps(18)%FermLoopPart = Chm_

      PrimAmps(19)%ExtLine = (/1,2,4,5,3/)
      PrimAmp2_12453 = 19
      PrimAmps(19)%AmpType = 2
      PrimAmps(19)%FermLoopPart = Chm_

      PrimAmps(20)%ExtLine = (/1,2,4,3,5/)
      PrimAmp2_12435 = 20
      PrimAmps(20)%AmpType = 2
      PrimAmps(20)%FermLoopPart = Chm_




      PrimAmps(21)%ExtLine = (/1,5,2,3,4/)
      PrimAmp2m_15234 = 21
      PrimAmps(21)%AmpType = 2
      PrimAmps(21)%FermLoopPart = Bot_

      PrimAmps(22)%ExtLine = (/1,2,3,4,5/)
      PrimAmp2m_12345 = 22
      PrimAmps(22)%AmpType = 2
      PrimAmps(22)%FermLoopPart = Bot_

      PrimAmps(23)%ExtLine = (/1,2,3,5,4/)
      PrimAmp2m_12354 = 23
      PrimAmps(23)%AmpType = 2
      PrimAmps(23)%FermLoopPart = Bot_

      PrimAmps(24)%ExtLine = (/1,2,5,3,4/)
      PrimAmp2m_12534 = 24
      PrimAmps(24)%AmpType = 2
      PrimAmps(24)%FermLoopPart = Bot_

      PrimAmps(25)%ExtLine = (/1,5,2,4,3/)
      PrimAmp2m_15243 = 25
      PrimAmps(25)%AmpType = 2
      PrimAmps(25)%FermLoopPart = Bot_

      PrimAmps(26)%ExtLine = (/1,2,5,4,3/)
      PrimAmp2m_12543 = 26
      PrimAmps(26)%AmpType = 2
      PrimAmps(26)%FermLoopPart = Bot_

      PrimAmps(27)%ExtLine = (/1,2,4,5,3/)
      PrimAmp2m_12453 = 27
      PrimAmps(27)%AmpType = 2
      PrimAmps(27)%FermLoopPart = Bot_

      PrimAmps(28)%ExtLine = (/1,2,4,3,5/)
      PrimAmp2m_12435 = 28
      PrimAmps(28)%AmpType = 2
      PrimAmps(28)%FermLoopPart = Bot_


   ENDIF


ELSEIF( MASTERPROCESS.EQ.9 ) THEN! tb t qb q pho

   IF( Correction.EQ.0  .OR. Correction.EQ.4 .OR.Correction.EQ.5 ) THEN
      BornAmps(1)%ExtLine = (/1,5,2,3,4/)
      BornAmps(2)%ExtLine = (/1,2,3,5,4/)

      PrimAmps(1)%ExtLine = (/1,5,2,3,4/)
      PrimAmps(2)%ExtLine = (/1,2,3,5,4/)
   ELSEIF( Correction.EQ.1 ) THEN
      BornAmps(1)%ExtLine = (/1,5,2,3,4/)
      BornAmps(2)%ExtLine = (/1,2,3,5,4/)

      PrimAmps(1)%ExtLine = (/1,5,2,3,4/)
      PrimAmp1_15234 = 1
      PrimAmps(1)%AmpType = 1

      PrimAmps(2)%ExtLine = (/1,2,3,5,4/)
      PrimAmp1_12354 = 2
      PrimAmps(2)%AmpType = 1

      PrimAmps(3)%ExtLine = (/1,5,2,4,3/)
      PrimAmp1_15243 = 3
      PrimAmps(3)%AmpType = 1

      PrimAmps(4)%ExtLine = (/1,2,4,5,3/)
      PrimAmp1_12453 = 4
      PrimAmps(4)%AmpType = 1

      PrimAmps(5)%ExtLine = (/1,5,4,3,2/)
      PrimAmp3_15432 = 5
      PrimAmps(5)%AmpType = 3

      PrimAmps(6)%ExtLine = (/1,4,5,3,2/)
      PrimAmp3_14532 = 6
      PrimAmps(6)%AmpType = 3

      PrimAmps(7)%ExtLine = (/1,4,3,5,2/)
      PrimAmp3_14352 = 7
      PrimAmps(7)%AmpType = 3

      PrimAmps(8)%ExtLine = (/1,5,2,3,4/)
      PrimAmp4_15234 = 8
      PrimAmps(8)%AmpType = 4

      PrimAmps(9)%ExtLine = (/1,2,5,3,4/)
      PrimAmp4_12534 = 9
      PrimAmps(9)%AmpType = 4

      PrimAmps(10)%ExtLine = (/1,2,3,4,5/)
      PrimAmp4_12345 = 10
      PrimAmps(10)%AmpType = 4



      PrimAmps(11)%ExtLine = (/1,5,2,3,4/)
      PrimAmp2_15234 = 11
      PrimAmps(11)%AmpType = 2
      PrimAmps(11)%FermLoopPart = Chm_

      PrimAmps(12)%ExtLine = (/1,2,3,5,4/)
      PrimAmp2_12354 = 12
      PrimAmps(12)%AmpType = 2
      PrimAmps(12)%FermLoopPart = Chm_

!       PrimAmps(13)%ExtLine = (/1,2,5,3,4/)
!       PrimAmp2_12534 = 13
!       PrimAmps(13)%AmpType = 2
!       PrimAmps(13)%FermLoopPart = Chm_
!
!       PrimAmps(14)%ExtLine = (/1,2,3,4,5/)
!       PrimAmp2_12345 = 14
!       PrimAmps(14)%AmpType = 2
!       PrimAmps(14)%FermLoopPart = Chm_


      PrimAmps(13)%ExtLine = (/1,5,2,3,4/)
      PrimAmp2m_15234 = 13
      PrimAmps(13)%AmpType = 2
      PrimAmps(13)%FermLoopPart = Bot_

      PrimAmps(14)%ExtLine = (/1,2,3,5,4/)
      PrimAmp2m_12354 = 14
      PrimAmps(14)%AmpType = 2
      PrimAmps(14)%FermLoopPart = Bot_

!       PrimAmps(17)%ExtLine = (/1,2,5,3,4/)
!       PrimAmp2m_12534 = 17
!       PrimAmps(17)%AmpType = 2
!       PrimAmps(17)%FermLoopPart = Bot_
!
!       PrimAmps(18)%ExtLine = (/1,2,3,4,5/)
!       PrimAmp2m_12345 = 18
!       PrimAmps(18)%AmpType = 2
!       PrimAmps(18)%FermLoopPart = Bot_

   ENDIF


ELSEIF( MasterProcess.EQ.10 ) THEN

   IF( Correction.EQ.2 ) THEN
      PrimAmps(1)%ExtLine = (/1,6,2,3,4,5/)
      BornAmps(1)%ExtLine = (/1,6,2,3,4,5/)

      PrimAmps(2)%ExtLine = (/1,6,2,3,5,4/)
      BornAmps(2)%ExtLine = (/1,6,2,3,5,4/)

      PrimAmps(3)%ExtLine = (/1,6,2,4,3,5/)
      BornAmps(3)%ExtLine = (/1,6,2,4,3,5/)

      PrimAmps(4)%ExtLine = (/1,6,2,4,5,3/)
      BornAmps(4)%ExtLine = (/1,6,2,4,5,3/)

      PrimAmps(5)%ExtLine = (/1,6,2,5,3,4/)
      BornAmps(5)%ExtLine = (/1,6,2,5,3,4/)

      PrimAmps(6)%ExtLine = (/1,6,2,5,4,3/)
      BornAmps(6)%ExtLine = (/1,6,2,5,4,3/)
   ENDIF



ELSEIF( MasterProcess.EQ.11 ) THEN

   IF( Correction.EQ.2 ) THEN
      PrimAmps( 1)%ExtLine = (/1,6,2,3,4,5/)
      BornAmps( 1)%ExtLine = (/1,6,2,3,4,5/)
      PrimAmp1_162345 = 1

      PrimAmps( 2)%ExtLine = (/1,2,3,6,4,5/)
      BornAmps( 2)%ExtLine = (/1,2,3,6,4,5/)
      PrimAmp1_123645 = 2

      PrimAmps( 3)%ExtLine = (/1,6,2,5,3,4/)
      BornAmps( 3)%ExtLine = (/1,6,2,5,3,4/)
      PrimAmp1_162534 = 3

      PrimAmps( 4)%ExtLine = (/1,2,5,3,6,4/)
      BornAmps( 4)%ExtLine = (/1,2,5,3,6,4/)
      PrimAmp1_125364 = 4

      PrimAmps( 5)%ExtLine = (/1,5,6,2,3,4/)
      BornAmps( 5)%ExtLine = (/1,5,6,2,3,4/)
      PrimAmp1_156234 = 5

      PrimAmps( 6)%ExtLine = (/1,6,5,2,3,4/)
      BornAmps( 6)%ExtLine = (/1,6,5,2,3,4/)
      PrimAmp1_165234 = 6

      PrimAmps( 7)%ExtLine = (/1,5,2,3,6,4/)
      BornAmps( 7)%ExtLine = (/1,5,2,3,6,4/)
      PrimAmp1_152364 = 7

      PrimAmps( 8)%ExtLine = (/1,6,2,3,5,4/)
      BornAmps( 8)%ExtLine = (/1,6,2,3,5,4/)
      PrimAmp1_162354 = 8

      PrimAmps( 9)%ExtLine = (/1,2,3,5,6,4/)
      BornAmps( 9)%ExtLine = (/1,2,3,5,6,4/)
      PrimAmp1_123564 = 9

      PrimAmps(10)%ExtLine = (/1,2,3,6,5,4/)
      BornAmps(10)%ExtLine = (/1,2,3,6,5,4/)
      PrimAmp1_123654 = 10
   ENDIF




ELSE
    call Error("MasterProcess not implemented in InitAmps")

ENDIF





!  loop over all born amplitudes
!    do NPrimAmp=1,NumBornAmps
   do NPrimAmp=1,NumPrimAmps
!!!!!! overwriting bornamps initializations  !!!!
            BornAmps(NPrimAmp)%ExtLine = PrimAmps(NPrimAmp)%ExtLine
!!!!!!
            TheBornAmp => BornAmps(NPrimAmp)
            TheTree => TheBornAmp%TreeProc
            TheTree%NumPart = NumExtParticles
            allocate( TheTree%PartType(1:NumExtParticles), stat=AllocStatus )
            if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartType for Born")
            allocate( TheTree%PartRef(1:NumExtParticles), stat=AllocStatus )
            if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartRef for Born")
            TheTree%PartRef(1:NumExtParticles) = BornAmps(NPrimAmp)%ExtLine(1:NumExtParticles)
!           set number of quarks and gluons
            TheTree%NumQua = 0
            counter = 0

            do NPart=1,TheTree%NumPart
                  TheTree%PartType(NPart) = ExtParticle( TheBornAmp%ExtLine(NPart) )%PartType
                  if( IsAQuark(TheTree%PartType(NPart)) ) then
                     TheTree%NumQua = TheTree%NumQua + 1
                     counter = counter + 1
                     QuarkPos(counter) = NPart
                  endif
            enddo

            if( IsAQuark(TheTree%PartType(1)) ) then
               allocate( TheTree%NumGlu(0:TheTree%NumQua), stat=AllocStatus )
               TheTree%NumGlu(0:TheTree%NumQua) = 0
            elseif( TheTree%PartType(1).eq.Glu_ ) then
               allocate( TheTree%NumGlu(0:TheTree%NumQua+1), stat=AllocStatus )
               TheTree%NumGlu(0:TheTree%NumQua) = 0
            else
               call Error("TheTree%NumGlu")
            endif
            if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%NumGlu")
            do NPart=1,TheTree%NumPart
                  if( TheTree%PartType(NPart) .eq. Glu_ ) then
                     TheTree%NumGlu(0) = TheTree%NumGlu(0) + 1
                  endif
            enddo

!           set number of gluons between quark lines
            if( IsAQuark(TheTree%PartType(1)) ) then
            if( TheTree%NumQua .eq. 2 ) then
                  TheTree%NumGlu(1) = QuarkPos(2) - QuarkPos(1) - 1
                  TheTree%NumGlu(2) = TheTree%NumPart - QuarkPos(2)
            endif
            if( TheTree%NumQua .eq. 4 ) then
                  TheTree%NumGlu(1) = QuarkPos(2) - QuarkPos(1) - 1
                  TheTree%NumGlu(2) = QuarkPos(3) - QuarkPos(2) - 1
                  TheTree%NumGlu(3) = QuarkPos(4) - QuarkPos(3) - 1
                  TheTree%NumGlu(4) = TheTree%NumPart - QuarkPos(4)
            endif
            if( TheTree%NumQua .eq. 6 ) then
                  TheTree%NumGlu(1) = QuarkPos(2) - QuarkPos(1) - 1
                  TheTree%NumGlu(2) = QuarkPos(3) - QuarkPos(2) - 1
                  TheTree%NumGlu(3) = QuarkPos(4) - QuarkPos(3) - 1
                  TheTree%NumGlu(4) = QuarkPos(5) - QuarkPos(4) - 1
                  TheTree%NumGlu(5) = QuarkPos(6) - QuarkPos(5) - 1
                  TheTree%NumGlu(6) = TheTree%NumPart - QuarkPos(6)
            endif
            elseif( TheTree%PartType(1).eq.Glu_ ) then
            if( TheTree%NumQua .eq. 2 ) then
                  TheTree%NumGlu(1) = QuarkPos(1) - 2
                  TheTree%NumGlu(2) = QuarkPos(2) - QuarkPos(1) - 1
                  TheTree%NumGlu(3) = TheTree%NumPart - QuarkPos(2)
            endif
            if( TheTree%NumQua .eq. 4 ) then
                  TheTree%NumGlu(1) = QuarkPos(1) - 2
                  TheTree%NumGlu(2) = QuarkPos(2) - QuarkPos(1) - 1
                  TheTree%NumGlu(3) = QuarkPos(3) - QuarkPos(2) - 1
                  TheTree%NumGlu(4) = QuarkPos(4) - QuarkPos(3) - 1
                  TheTree%NumGlu(5) = TheTree%NumPart - QuarkPos(4)
            endif
            if( TheTree%NumQua .eq. 6 ) then
                  TheTree%NumGlu(1) = QuarkPos(1) - 2
                  TheTree%NumGlu(2) = QuarkPos(2) - QuarkPos(1) - 1
                  TheTree%NumGlu(3) = QuarkPos(3) - QuarkPos(2) - 1
                  TheTree%NumGlu(4) = QuarkPos(4) - QuarkPos(3) - 1
                  TheTree%NumGlu(5) = QuarkPos(5) - QuarkPos(4) - 1
                  TheTree%NumGlu(6) = QuarkPos(6) - QuarkPos(5) - 1
                  TheTree%NumGlu(7) = TheTree%NumPart - QuarkPos(6)
            endif
            endif

!          allocate memory for pointer to quarks
           allocate( TheTree%Quarks(1:TheTree%NumQua), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%Quarks")
!          allocate memory for pointer to gluons
           allocate( TheTree%Gluons(1:TheTree%NumGlu(0)), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%Gluons")

           counterQ = 0
           counterG = 0

           do NPart=1,TheTree%NumPart
               if( IsAQuark(TheTree%PartType(NPart)) ) then
                     counterQ = counterQ + 1
                     TheTree%Quarks(counterQ)%PartType => ExtParticle( TheBornAmp%ExtLine(NPart) )%PartType
                     TheTree%Quarks(counterQ)%ExtRef => ExtParticle( TheBornAmp%ExtLine(NPart) )%ExtRef
                     TheTree%Quarks(counterQ)%Mass => ExtParticle( TheBornAmp%ExtLine(NPart) )%Mass
                     TheTree%Quarks(counterQ)%Mass2 => ExtParticle( TheBornAmp%ExtLine(NPart) )%Mass2
                     TheTree%Quarks(counterQ)%Helicity => ExtParticle( TheBornAmp%ExtLine(NPart) )%Helicity
                     TheTree%Quarks(counterQ)%Mom => ExtParticle( TheBornAmp%ExtLine(NPart) )%Mom
                     TheTree%Quarks(counterQ)%Pol => ExtParticle( TheBornAmp%ExtLine(NPart) )%Pol
               endif
               if( TheTree%PartType(NPart) .eq. Glu_ ) then
                     counterG = counterG + 1
                     TheTree%Gluons(counterG)%PartType => ExtParticle( TheBornAmp%ExtLine(NPart) )%PartType
                     TheTree%Gluons(counterG)%ExtRef => ExtParticle( TheBornAmp%ExtLine(NPart) )%ExtRef
                     TheTree%Gluons(counterG)%Mass => ExtParticle( TheBornAmp%ExtLine(NPart) )%Mass
                     TheTree%Gluons(counterG)%Mass2 => ExtParticle( TheBornAmp%ExtLine(NPart) )%Mass2
                     TheTree%Gluons(counterG)%Helicity => ExtParticle( TheBornAmp%ExtLine(NPart) )%Helicity
                     TheTree%Gluons(counterG)%Mom => ExtParticle( TheBornAmp%ExtLine(NPart) )%Mom
                     TheTree%Gluons(counterG)%Pol => ExtParticle( TheBornAmp%ExtLine(NPart) )%Pol
               endif
           enddo
   enddo


IF( Correction.EQ.1 ) THEN
!  loop over all primitive amplitudes
   do NPrimAmp=1,NumPrimAmps

         ThePrimAmp => PrimAmps(NPrimAmp)
         ColorLessParticles = .false.
         ThePrimAmp%IntPart(1:NumExtParticles)%PartType = 99

         ExtPartType = ExtParticle( ThePrimAmp%ExtLine(1) )%PartType
!  set internal lines: associate each int.line with the ext.line at the next vertex
!  negativ  IntPart()%PartType <--> fermion flow along ascending  propagators
!  positive IntPart()%PartType <--> fermion flow along descending propagators


!  determine number of quark lines
         ThePrimAmp%FermLine1In = 0
         ThePrimAmp%FermLine1Out= 0
         ThePrimAmp%FermLine2In = 0
         ThePrimAmp%FermLine2Out= 0
         do Vertex=1,NumExtParticles
            Propa       = Vertex + 1
            PropaMinus1 = Propa - 1
            ExtPartType = ExtParticle( ThePrimAmp%ExtLine(Vertex) )%PartType
            if( Vertex .eq. NumExtParticles ) then
               Propa       = 1
               PropaMinus1 = NumExtParticles
            endif

            if ( Vertex.eq.1 ) then
               ThePrimAmp%IntPart(Propa)%PartType = ExtPartType
               if( IsAQuark(ExtPartType) ) then
                  ThePrimAmp%FermionLines = 1
                  ThePrimAmp%FermLine1In  = 1
               else
                  ThePrimAmp%FermionLines = 0
               endif

            elseif( IsAQuark(ExtPartType) ) then
               if( ThePrimAmp%AmpType.eq.1 ) then
                   if( ThePrimAmp%FermLine1Out.eq.0 ) then
                      ThePrimAmp%FermLine1Out = Vertex
                      ThePrimAmp%IntPart(Propa)%PartType = Glu_
                   elseif( ThePrimAmp%FermLine2In.eq.0 ) then
                      ThePrimAmp%FermionLines = 2
                      ThePrimAmp%FermLine2In = Vertex
                      ThePrimAmp%IntPart(Propa)%PartType = ExtPartType
                   elseif( ThePrimAmp%FermLine2Out.eq.0 ) then
                      ThePrimAmp%FermLine2Out = Vertex
                      ThePrimAmp%IntPart(Propa)%PartType = Glu_
                   endif

               elseif( ThePrimamp%AmpType.eq.2 ) then
                   if( ThePrimAmp%FermLine1Out.eq.0 ) then
                            ThePrimAmp%FermLine1Out = Vertex
                            ThePrimAmp%IntPart(1)%PartType = -abs(ThePrimAmp%FermLoopPart)
                            ThePrimAmp%IntPart(Propa)%PartType = -abs(ThePrimAmp%FermLoopPart)
                            do k=ThePrimAmp%FermLine1In+1,ThePrimAmp%FermLine1Out
                                ThePrimAmp%IntPart(k)%PartType = 0
                            enddo
                   elseif( ThePrimAmp%FermLine2In.eq.0 ) then
                            ThePrimAmp%FermionLines = 2
                            ThePrimAmp%FermLine2In = Vertex
                   elseif( ThePrimAmp%FermLine2Out.eq.0 ) then
                            ThePrimAmp%FermLine2Out = Vertex
                            do k=ThePrimAmp%FermLine2In+1,ThePrimAmp%FermLine2Out
                                ThePrimAmp%IntPart(k)%PartType = 0
                            enddo
                            do k=ThePrimAmp%FermLine2Out+1,NumExtParticles
                                ThePrimAmp%IntPart(k)%PartType = -abs(ThePrimAmp%FermLoopPart)
                            enddo
                   endif

               elseif( ThePrimamp%AmpType.eq.3 ) then
                       if( ThePrimAmp%FermLine2In.eq.0 ) then
                          ThePrimAmp%FermionLines = 2
                          ThePrimAmp%FermLine2In = Vertex
                       elseif( ThePrimAmp%FermLine2Out.eq.0 ) then
                          ThePrimAmp%FermLine2Out = Vertex
                          ThePrimAmp%IntPart(Propa)%PartType = ThePrimAmp%IntPart(ThePrimAmp%FermLine2In)%PartType
                          do k=ThePrimAmp%FermLine2In+1,ThePrimAmp%FermLine2Out
                              ThePrimAmp%IntPart(k)%PartType = 0
                          enddo
                       elseif( ThePrimAmp%FermLine1Out.eq.0 ) then
                           ThePrimAmp%FermLine1Out = Vertex
                           ThePrimAmp%IntPart(Propa)%PartType = Glu_
                       endif

               elseif( ThePrimamp%AmpType.eq.4 ) then
                       if( ThePrimAmp%FermLine1Out.eq.0 ) then
                           ThePrimAmp%FermLine1Out = Vertex
                           do k=ThePrimAmp%FermLine1In+1,ThePrimAmp%FermLine1Out
                              ThePrimAmp%IntPart(k)%PartType = 0
                           enddo
                       elseif( ThePrimAmp%FermLine2In.eq.0 ) then
                           ThePrimAmp%FermionLines = 2
                           ThePrimAmp%FermLine2In = Vertex
                            ThePrimAmp%IntPart(Propa)%PartType = Glu_
                            do k=ThePrimAmp%FermLine1Out+1,ThePrimAmp%FermLine2In
                                ThePrimAmp%IntPart(k)%PartType = -ExtPartType
                            enddo
                            ThePrimAmp%IntPart(1)%PartType = -ExtPartType
                       elseif( ThePrimAmp%FermLine2Out.eq.0 ) then
                            ThePrimAmp%FermLine2Out = Vertex
                            do k=ThePrimAmp%FermLine2Out+1,NumExtParticles
                                ThePrimAmp%IntPart(k)%PartType = ExtPartType
                            enddo
                       endif
                       print *, "remember: check again this code, int/ext particles and trees"
               endif

            elseif( (ExtPartType .eq. Glu_) ) then
               ThePrimAmp%IntPart(Propa)%PartType = ThePrimAmp%IntPart(PropaMinus1)%PartType

            elseif( (ExtPartType .eq. Pho_) .or. (ExtPartType .eq. Z0_) ) then
               ThePrimAmp%IntPart(Propa)%PartType = ThePrimAmp%IntPart(PropaMinus1)%PartType
               ColorLessParticles = .true.

            elseif( ExtPartType .eq. Wp_ ) then
               ThePrimAmp%IntPart(Propa)%PartType = abs( ThePrimAmp%IntPart(PropaMinus1)%PartType -1)
               ColorLessParticles = .true.

            elseif( ExtPartType .eq. Wm_ ) then
               ThePrimAmp%IntPart(Propa)%PartType = abs( ThePrimAmp%IntPart(PropaMinus1)%PartType +1)
               ColorLessParticles = .true.
            endif
         enddo

         do Propa=1,NumExtParticles
            if( ThePrimAmp%IntPart(Propa)%PartType .eq. 99 ) call Error("internal particle type is 99")
            ThePrimAmp%IntPart(Propa)%Mass  = GetMass( ThePrimAmp%IntPart(Propa)%PartType )
            ThePrimAmp%IntPart(Propa)%Mass2 = (ThePrimAmp%IntPart(Propa)%Mass)**2
            ThePrimAmp%IntPart(Propa)%ExtRef = -1
         enddo

!        set number of possible insertions of colorless particles into quark lines
!          if ( ColorLessParticles ) then
!             ThePrimAmp%NumInsertions1 = ThePrimAmp%FermLine1Out - ThePrimAmp%FermLine1In - 1
!             if( ThePrimAmp%FermionLines .eq. 2 ) then
!                ThePrimAmp%NumInsertions2 = ThePrimAmp%FermLine2Out - ThePrimAmp%FermLine2In - 1
!             endif
!          endif

         call InitUCuts(ThePrimAmp)
   enddo
ENDIF

RETURN
END SUBROUTINE






SUBROUTINE SetKirill(ThePrimAmp)
use ModMisc
use ModParameters
implicit none
integer :: Vertex,Propa,PropaMinus1,ExtPartType,NPrimAmp,k
integer :: AllocStatus,counter,counterQ,counterG,QuarkPos(1:6),NPart
logical :: ColorLessParticles
type(PrimitiveAmplitude) :: ThePrimAmp
type(BornAmplitude),pointer :: TheBornAmp
type(TreeProcess),pointer :: TheTree
integer :: h1,h2,h3,h4,h5,h6
complex(8) e(1:NumExtParticles,1:4)
integer :: i
include 'misc/global_import'

!     conversion to kirills conv for last prim. ampl.
   NPoint = NumExtParticles

!      if(NumExtParticles.ge.4) then
!         h1=ExtParticle(1)%Helicity
!         h2=ExtParticle(2)%Helicity
!         h3=ExtParticle(3)%Helicity
!         h4=ExtParticle(4)%Helicity
!      endif
!      if(NumExtParticles.eq.5) then
! !         h5=ExtParticle(5)%Helicity
! !         print *, "set h5!"
!      endif
!      if(NumExtParticles.eq.6) then
! !         h6=ExtParticle(6)%Helicity
!           print *, "set h6!"
!      endif

      do i=1,NumExtParticles
         mom(i,1:4) = ExtParticle( ThePrimAmp%ExtLine(i) )%Mom(1:4)
         hel(i,1:4) = ExtParticle( ThePrimAmp%ExtLine(i) )%Pol(1:4)
      enddo


!       the list of ``propagator momenta''
      do i=1,4
         momline(1,i)=dcmplx(0d0,0d0)
         momline(2,i)=mom(1,i)
         momline(3,i)=mom(1,i)+mom(2,i)
         momline(4,i)=mom(1,i)+mom(2,i)+mom(3,i)
         momline(5,i)=mom(1,i)+mom(2,i)+mom(3,i)+mom(4,i)
         momline(6,i)=mom(1,i)+mom(2,i)+mom(3,i)+mom(4,i)+mom(5,i)
      enddo


   do Vertex=1,NumExtParticles
            Propa       = Vertex + 1
            PropaMinus1 = Propa - 1
            ExtPartType = ExtParticle( ThePrimAmp%ExtLine(Vertex) )%PartType
            if( Vertex .eq. NumExtParticles ) then
               Propa       = 1
               PropaMinus1 = NumExtParticles
            endif

      if(ExtPartType.eq.Top_ .or. ExtPartType.eq.ATop_ ) then
         Lab_ex(Vertex)='top'
      elseif(ExtPartType.eq.Bot_ .or. ExtPartType.eq.ABot_ ) then
         Lab_ex(Vertex)='bot'
      elseif(ExtPartType.eq.Chm_ .or. ExtPartType.eq.AChm_ ) then
         Lab_ex(Vertex)='chm'
      elseif(ExtPartType.eq.Str_ .or. ExtPartType.eq.AStr_ ) then
         Lab_ex(Vertex)='str'
      elseif(ExtPartType.eq.Glu_ ) then
         Lab_ex(Vertex)='glu'
      else
         print *, "error in kirills conv, ExtPartType=", ExtPartType
      endif

      if( ThePrimAmp%IntPart(Propa)%PartType.eq.Top_ .or. ThePrimAmp%IntPart(Propa)%PartType.eq.ATop_ ) then
         Lab_in(Propa)='top'
      elseif( ThePrimAmp%IntPart(Propa)%PartType.eq.Bot_ .or. ThePrimAmp%IntPart(Propa)%PartType.eq.ABot_ ) then
         Lab_in(Propa)='bot'
      elseif( ThePrimAmp%IntPart(Propa)%PartType.eq.Chm_ .or. ThePrimAmp%IntPart(Propa)%PartType.eq.AChm_ ) then
         Lab_in(Propa)='chm'
      elseif( ThePrimAmp%IntPart(Propa)%PartType.eq.Str_ .or. ThePrimAmp%IntPart(Propa)%PartType.eq.AStr_ ) then
         Lab_in(Propa)='str'
      elseif( ThePrimAmp%IntPart(Propa)%PartType.eq.Glu_ ) then
         Lab_in(Propa)='glu'
      else
         Lab_in(Propa)='notset'
      endif
  enddo

       N5=ThePrimAmp%UCuts(5)%NumCuts
       N4=ThePrimAmp%UCuts(4)%NumCuts
       N3=ThePrimAmp%UCuts(3)%NumCuts
       N2=ThePrimAmp%UCuts(2)%NumCuts
       N1=ThePrimAmp%UCuts(1)%NumCuts

       do k=1,N5
         Lc5(k,1)=ThePrimAmp%UCuts(5)%CutProp(k,1)
         Lc5(k,2)=ThePrimAmp%UCuts(5)%CutProp(k,2)
         Lc5(k,3)=ThePrimAmp%UCuts(5)%CutProp(k,3)
         Lc5(k,4)=ThePrimAmp%UCuts(5)%CutProp(k,4)
         Lc5(k,5)=ThePrimAmp%UCuts(5)%CutProp(k,5)
       enddo
       do k=1,N4
         Lc4(k,1)=ThePrimAmp%UCuts(4)%CutProp(k,1)
         Lc4(k,2)=ThePrimAmp%UCuts(4)%CutProp(k,2)
         Lc4(k,3)=ThePrimAmp%UCuts(4)%CutProp(k,3)
         Lc4(k,4)=ThePrimAmp%UCuts(4)%CutProp(k,4)
       enddo
       do k=1,N3
         Lc3(k,1)=ThePrimAmp%UCuts(3)%CutProp(k,1)
         Lc3(k,2)=ThePrimAmp%UCuts(3)%CutProp(k,2)
         Lc3(k,3)=ThePrimAmp%UCuts(3)%CutProp(k,3)
       enddo
       do k=1,N2
         Lc2(k,1)=ThePrimAmp%UCuts(2)%CutProp(k,1)
         Lc2(k,2)=ThePrimAmp%UCuts(2)%CutProp(k,2)
       enddo
       do k=1,N1
         Lc1(k,1)=ThePrimAmp%UCuts(1)%CutProp(k,1)
       enddo

return
END SUBROUTINE




SUBROUTINE InitUCuts(ThePrimAmp)
use ModMisc
use ModParameters
implicit none
integer :: AllocStatus,NCut
integer :: i1,i2,i3,i4,i5
type(PrimitiveAmplitude),target :: ThePrimAmp
integer :: NumVertPart,NPart,NPoint,NTree
logical :: MasslessExtLeg,MasslessIntParticles
integer :: QuarkPos(1:6),counter,counterQ,counterG
type(TreeProcess),pointer :: TheTree



   if( ThePrimAmp%AmpType.eq.1 ) then
      ThePrimAmp%NPoint = NumExtParticles
   elseif( ThePrimAmp%AmpType.eq.2 ) then
      ThePrimAmp%NPoint = NumExtParticles-(ThePrimAmp%FermLine1Out-ThePrimAmp%FermLine1In)-(ThePrimAmp%FermLine2Out-ThePrimAmp%FermLine2In)
   elseif( ThePrimAmp%AmpType.eq.3 ) then
      ThePrimAmp%NPoint = NumExtParticles-(ThePrimAmp%FermLine2Out-ThePrimAmp%FermLine2In)
   elseif( ThePrimAmp%AmpType.eq.4 ) then
      ThePrimAmp%NPoint = NumExtParticles-(ThePrimAmp%FermLine1Out-ThePrimAmp%FermLine1In)
   endif

!  set number of cuts
   ThePrimAmp%UCuts(5)%CutType = 5
   if ( NumExtParticles .ge. 5) then
      ThePrimAmp%UCuts(5)%NumCuts = Binomial(ThePrimAmp%NPoint,5)
   else
      ThePrimAmp%UCuts(5)%NumCuts = 0
   endif

   ThePrimAmp%UCuts(4)%CutType = 4
   ThePrimAmp%UCuts(4)%NumCuts = Binomial(ThePrimAmp%NPoint,4)

   ThePrimAmp%UCuts(3)%CutType = 3
   ThePrimAmp%UCuts(3)%NumCuts = Binomial(ThePrimAmp%NPoint,3)

   ThePrimAmp%UCuts(2)%CutType = 2
   ThePrimAmp%UCuts(2)%NumCuts = Binomial(ThePrimAmp%NPoint,2)

   ThePrimAmp%UCuts(1)%CutType = 1
   ThePrimAmp%UCuts(1)%NumCuts = Binomial(ThePrimAmp%NPoint,1)


!  allocate memory for CutProp
   allocate( ThePrimAmp%UCuts(5)%CutProp(1:ThePrimAmp%UCuts(5)%NumCuts,1:5), stat=AllocStatus)
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in ThePrimAmp%UCuts(5)")
   allocate( ThePrimAmp%UCuts(5)%Coeff(1:ThePrimAmp%UCuts(5)%NumCuts,0:0),     stat=AllocStatus)
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in ThePrimAmp%UCuts(5)")
   allocate( ThePrimAmp%UCuts(5)%KMom(1:ThePrimAmp%UCuts(5)%NumCuts,1:4,1:4) )
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in KMom 5")

   allocate( ThePrimAmp%UCuts(4)%CutProp(1:ThePrimAmp%UCuts(4)%NumCuts,1:4), stat=AllocStatus)
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in ThePrimAmp%UCuts(4)")
   allocate( ThePrimAmp%UCuts(4)%Coeff(1:ThePrimAmp%UCuts(4)%NumCuts,0:4),   stat=AllocStatus)
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in ThePrimAmp%UCuts(4)")
   allocate( ThePrimAmp%UCuts(4)%KMom(1:ThePrimAmp%UCuts(4)%NumCuts,1:3,1:4) )
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in KMom 4")
   allocate( ThePrimAmp%UCuts(4)%NMom(1:ThePrimAmp%UCuts(4)%NumCuts,1:1,1:4) )
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in NMom 4")

   allocate( ThePrimAmp%UCuts(3)%CutProp(1:ThePrimAmp%UCuts(3)%NumCuts,1:3), stat=AllocStatus)
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in ThePrimAmp%UCuts(3)")
   allocate( ThePrimAmp%UCuts(3)%Coeff(1:ThePrimAmp%UCuts(3)%NumCuts,0:9),   stat=AllocStatus)
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in ThePrimAmp%UCuts(3)")
   allocate( ThePrimAmp%UCuts(3)%KMom(1:ThePrimAmp%UCuts(3)%NumCuts,1:2,1:4) )
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in KMom 3")
   allocate( ThePrimAmp%UCuts(3)%NMom(1:ThePrimAmp%UCuts(3)%NumCuts,1:2,1:4) )
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in NMom 3")

   allocate( ThePrimAmp%UCuts(2)%CutProp(1:ThePrimAmp%UCuts(2)%NumCuts,1:2), stat=AllocStatus)
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in ThePrimAmp%UCuts(2)")
   allocate( ThePrimAmp%UCuts(2)%Coeff(1:ThePrimAmp%UCuts(2)%NumCuts,0:9),   stat=AllocStatus)
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in ThePrimAmp%UCuts(2)")
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in ThePrimAmp%UCuts(3)")
   allocate( ThePrimAmp%UCuts(2)%KMom(1:ThePrimAmp%UCuts(2)%NumCuts,1:1,1:4) )
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in KMom 2")
   allocate( ThePrimAmp%UCuts(2)%NMom(1:ThePrimAmp%UCuts(2)%NumCuts,1:3,1:4) )
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in NMom 2")

   allocate( ThePrimAmp%UCuts(1)%CutProp(1:ThePrimAmp%UCuts(1)%NumCuts,1:1), stat=AllocStatus)
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in ThePrimAmp%UCuts(1)")
   allocate( ThePrimAmp%UCuts(1)%Coeff(1:ThePrimAmp%UCuts(1)%NumCuts,0:0),     stat=AllocStatus)
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in ThePrimAmp%UCuts(1)")



!  init pentcuts
   allocate(ThePrimAmp%UCuts(5)%TreeProcess(1:ThePrimAmp%UCuts(5)%NumCuts,1:5), stat=AllocStatus)
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in ThePrimAmp%UCuts(5)%TreeProcess(ThePrimAmp%UCuts(5)%NumCuts,5)")
   NCut = 1
   do i1 = 1,    NumExtParticles-4
   do i2 = i1+1, NumExtParticles-3
   do i3 = i2+1, NumExtParticles-2
   do i4 = i3+1, NumExtParticles-1
   do i5 = i4+1, NumExtParticles

!          skip cuts that are marked for discarding
           if( ThePrimAmp%IntPart(i1)%PartType .eq. 0 ) cycle
           if( ThePrimAmp%IntPart(i2)%PartType .eq. 0 ) cycle
           if( ThePrimAmp%IntPart(i3)%PartType .eq. 0 ) cycle
           if( ThePrimAmp%IntPart(i4)%PartType .eq. 0 ) cycle
           if( ThePrimAmp%IntPart(i5)%PartType .eq. 0 ) cycle

!          set tree processes
!          vertex 1
           TheTree => ThePrimAmp%UCuts(5)%TreeProcess(NCut,1)
           NumVertPart = i2-i1
           TheTree%NumPart = NumVertPart+2
           allocate( TheTree%PartRef(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartRef 5-1")
           allocate( TheTree%PartType(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartType 5-1")
!          set particle reference (wrt.prim.amp) and flavor
           TheTree%PartRef(1) = i1
           TheTree%PartType(1) = (ThePrimAmp%IntPart(i1)%PartType)
           TheTree%PartRef(NumVertPart+2) = i2
           TheTree%PartType(NumVertPart+2) = ChargeConj(ThePrimAmp%IntPart(i2)%PartType)
           do NPart=0,NumVertPart-1
               TheTree%PartRef(NPart+2) = i1+NPart
               TheTree%PartType(NPart+2) = ExtParticle( ThePrimAmp%ExtLine(i1+NPart) )%PartType
           enddo

!          vertex 2
           TheTree => ThePrimAmp%UCuts(5)%TreeProcess(NCut,2)
           NumVertPart = i3-i2
           TheTree%NumPart = NumVertPart+2
           allocate( TheTree%PartRef(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartRef 5-2")
           allocate( TheTree%PartType(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartType 5-2")
!          set particle reference (wrt.prim.amp) and flavor
           TheTree%PartRef(1) = i2
           TheTree%PartType(1) = (ThePrimAmp%IntPart(i2)%PartType)
           TheTree%PartRef(NumVertPart+2) = i3
           TheTree%PartType(NumVertPart+2) = ChargeConj(ThePrimAmp%IntPart(i3)%PartType)
           do NPart=0,NumVertPart-1
               TheTree%PartRef(NPart+2) = i2+NPart
               TheTree%PartType(NPart+2) = ExtParticle( ThePrimAmp%ExtLine(i2+NPart) )%PartType
           enddo


!          vertex 3
           TheTree => ThePrimAmp%UCuts(5)%TreeProcess(NCut,3)
           NumVertPart = i4-i3
           TheTree%NumPart = NumVertPart+2
           allocate( TheTree%PartRef(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartRef 5-3")
           allocate( TheTree%PartType(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartType 5-3")
!          set particle reference (wrt.prim.amp) and flavor
           TheTree%PartRef(1) = i3
           TheTree%PartType(1) = (ThePrimAmp%IntPart(i3)%PartType)
           TheTree%PartRef(NumVertPart+2) = i4
           TheTree%PartType(NumVertPart+2) = ChargeConj(ThePrimAmp%IntPart(i4)%PartType)
           do NPart=0,NumVertPart-1
               TheTree%PartRef(NPart+2) = i3+NPart
               TheTree%PartType(NPart+2) = ExtParticle( ThePrimAmp%ExtLine(i3+NPart) )%PartType
           enddo


!          vertex 4
           TheTree => ThePrimAmp%UCuts(5)%TreeProcess(NCut,4)
           NumVertPart = i5-i4
           TheTree%NumPart = NumVertPart+2
           allocate( TheTree%PartRef(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartRef 5-4")
           allocate( TheTree%PartType(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartType 5-4")
!          set particle reference (wrt.prim.amp) and flavor
           TheTree%PartRef(1) = i4
           TheTree%PartType(1) = (ThePrimAmp%IntPart(i4)%PartType)
           TheTree%PartRef(NumVertPart+2) = i5
           TheTree%PartType(NumVertPart+2) = ChargeConj(ThePrimAmp%IntPart(i5)%PartType)
           do NPart=0,NumVertPart-1
               TheTree%PartRef(NPart+2) = i4+NPart
               TheTree%PartType(NPart+2) = ExtParticle( ThePrimAmp%ExtLine(i4+NPart) )%PartType
           enddo


!          vertex 5
           TheTree => ThePrimAmp%UCuts(5)%TreeProcess(NCut,5)
           NumVertPart = i1-(i5-NumExtParticles)
           TheTree%NumPart = NumVertPart+2
           allocate( TheTree%PartRef(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartRef 5-5")
           allocate( TheTree%PartType(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartType 5-5")
!          set particle reference (wrt.prim.amp) and flavor
           TheTree%PartRef(1) = i5
           TheTree%PartType(1) = (ThePrimAmp%IntPart(i5)%PartType)
           TheTree%PartRef(NumVertPart+2) = i1
           TheTree%PartType(NumVertPart+2) = ChargeConj(ThePrimAmp%IntPart(i1)%PartType)
           do NPart=0,NumVertPart-1
             if( i5+NPart.le.NumExtParticles ) then
                  TheTree%PartRef(NPart+2) = i5+NPart
                  TheTree%PartType(NPart+2) = ExtParticle( ThePrimAmp%ExtLine(i5+NPart) )%PartType
               else
                  TheTree%PartRef(NPart+2) = i5+NPart-NumExtParticles
                  TheTree%PartType(NPart+2) = ExtParticle( ThePrimAmp%ExtLine(i5+NPart-NumExtParticles) )%PartType
               endif
           enddo

           ThePrimAmp%UCuts(5)%CutProp(NCut,1) = i1
           ThePrimAmp%UCuts(5)%CutProp(NCut,2) = i2
           ThePrimAmp%UCuts(5)%CutProp(NCut,3) = i3
           ThePrimAmp%UCuts(5)%CutProp(NCut,4) = i4
           ThePrimAmp%UCuts(5)%CutProp(NCut,5) = i5

           NCut = NCut + 1
   enddo
   enddo
   enddo
   enddo
   enddo
   if ( NCut-1 .ne. ThePrimAmp%UCuts(5)%NumCuts ) call Error("Something went wrong while setting pent-cuts.")




!  init quadcuts
   allocate(ThePrimAmp%UCuts(4)%TreeProcess(1:ThePrimAmp%UCuts(4)%NumCuts,1:4), stat=AllocStatus)
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in ThePrimAmp%UCuts(4)%TreeProcess(ThePrimAmp%UCuts(4)%NumCuts,4).")
   NCut = 1
   do i2 = 1, NumExtParticles-3
   do i3 = i2+1, NumExtParticles-2
   do i4 = i3+1, NumExtParticles-1
   do i5 = i4+1, NumExtParticles


!          skip cuts that are marked for discarding
           if( ThePrimAmp%IntPart(i2)%PartType .eq. 0 ) cycle
           if( ThePrimAmp%IntPart(i3)%PartType .eq. 0 ) cycle
           if( ThePrimAmp%IntPart(i4)%PartType .eq. 0 ) cycle
           if( ThePrimAmp%IntPart(i5)%PartType .eq. 0 ) cycle

!          set tree processes
!          vertex 1
           TheTree => ThePrimAmp%UCuts(4)%TreeProcess(NCut,1)
           NumVertPart = i3-i2
           TheTree%NumPart = NumVertPart+2
           allocate( TheTree%PartRef(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartRef 4-1")
           allocate( TheTree%PartType(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartType 4-1")
!          set particle reference (wrt.prim.amp) and flavor
           TheTree%PartRef(1) = i2
           TheTree%PartType(1) = (ThePrimAmp%IntPart(i2)%PartType)
           TheTree%PartRef(NumVertPart+2) = i3
           TheTree%PartType(NumVertPart+2) = ChargeConj(ThePrimAmp%IntPart(i3)%PartType)
           do NPart=0,NumVertPart-1
               TheTree%PartRef(NPart+2) = i2+NPart
               TheTree%PartType(NPart+2) = ExtParticle( ThePrimAmp%ExtLine(i2+NPart) )%PartType
           enddo

!          vertex 2
           TheTree => ThePrimAmp%UCuts(4)%TreeProcess(NCut,2)
           NumVertPart = i4-i3
           TheTree%NumPart = NumVertPart+2
           allocate( TheTree%PartRef(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartRef 4-2")
           allocate( TheTree%PartType(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartType 4-2")
!          set particle reference (wrt.prim.amp) and flavor
           TheTree%PartRef(1) = i3
           TheTree%PartType(1) = (ThePrimAmp%IntPart(i3)%PartType)
           TheTree%PartRef(NumVertPart+2) = i4
           TheTree%PartType(NumVertPart+2) = ChargeConj(ThePrimAmp%IntPart(i4)%PartType)
           do NPart=0,NumVertPart-1
               TheTree%PartRef(NPart+2) = i3+NPart
               TheTree%PartType(NPart+2) = ExtParticle( ThePrimAmp%ExtLine(i3+NPart) )%PartType
           enddo

!          vertex 3
           TheTree => ThePrimAmp%UCuts(4)%TreeProcess(NCut,3)
           NumVertPart = i5-i4
           TheTree%NumPart = NumVertPart+2
           allocate( TheTree%PartRef(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartRef 4-3")
           allocate( TheTree%PartType(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartType 4-3")
!          set particle reference (wrt.prim.amp) and flavor
           TheTree%PartRef(1) = i4
           TheTree%PartType(1) = (ThePrimAmp%IntPart(i4)%PartType)
           TheTree%PartRef(NumVertPart+2) = i5
           TheTree%PartType(NumVertPart+2) = ChargeConj(ThePrimAmp%IntPart(i5)%PartType)
           do NPart=0,NumVertPart-1
               TheTree%PartRef(NPart+2) = i4+NPart
               TheTree%PartType(NPart+2) = ExtParticle( ThePrimAmp%ExtLine(i4+NPart) )%PartType
           enddo

!          vertex 4
           TheTree => ThePrimAmp%UCuts(4)%TreeProcess(NCut,4)
           NumVertPart = i2-(i5-NumExtParticles)
           TheTree%NumPart = NumVertPart+2
           allocate( TheTree%PartRef(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartRef 4-4")
           allocate( TheTree%PartType(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartType 4-4")
!          set particle reference (wrt.prim.amp) and flavor
           TheTree%PartRef(1) = i5
           TheTree%PartType(1) = (ThePrimAmp%IntPart(i5)%PartType)
           TheTree%PartRef(NumVertPart+2) = i2
           TheTree%PartType(NumVertPart+2) = ChargeConj(ThePrimAmp%IntPart(i2)%PartType)
           do NPart=0,NumVertPart-1
               if( i5+NPart.le.NumExtParticles ) then
                  TheTree%PartRef(NPart+2) = i5+NPart
                  TheTree%PartType(NPart+2) = ExtParticle( ThePrimAmp%ExtLine(i5+NPart) )%PartType
               else
                  TheTree%PartRef(NPart+2) = i5+NPart-NumExtParticles
                  TheTree%PartType(NPart+2) = ExtParticle( ThePrimAmp%ExtLine(i5+NPart-NumExtParticles) )%PartType
               endif
           enddo

           ThePrimAmp%UCuts(4)%CutProp(NCut,1) = i2
           ThePrimAmp%UCuts(4)%CutProp(NCut,2) = i3
           ThePrimAmp%UCuts(4)%CutProp(NCut,3) = i4
           ThePrimAmp%UCuts(4)%CutProp(NCut,4) = i5

           NCut = NCut + 1
   enddo
   enddo
   enddo
   enddo

   if ( NCut-1 .ne. ThePrimAmp%UCuts(4)%NumCuts ) call Error("Something went wrong while setting quad-cuts.")





!  init tripcuts
   allocate(ThePrimAmp%UCuts(3)%TreeProcess(1:ThePrimAmp%UCuts(3)%NumCuts,1:3), stat=AllocStatus)
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in ThePrimAmp%UCuts(3)%TreeProcess(ThePrimAmp%UCuts(3)%NumCuts,1:3).")
   NCut = 1
   do i3 = 1,    NumExtParticles-2
   do i4 = i3+1, NumExtParticles-1
   do i5 = i4+1, NumExtParticles

!          skip cuts that are marked for discarding
           if( ThePrimAmp%IntPart(i3)%PartType .eq. 0 ) cycle
           if( ThePrimAmp%IntPart(i4)%PartType .eq. 0 ) cycle
           if( ThePrimAmp%IntPart(i5)%PartType .eq. 0 ) cycle

!          set tree processes
!          vertex 1
           TheTree => ThePrimAmp%UCuts(3)%TreeProcess(NCut,1)
           NumVertPart = i4-i3
           TheTree%NumPart = NumVertPart+2
           allocate( TheTree%PartRef(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartRef 3-1")
           allocate( TheTree%PartType(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartType 3-1")
!          set particle reference (wrt.prim.amp) and flavor
           TheTree%PartRef(1) = i3
           TheTree%PartType(1) = (ThePrimAmp%IntPart(i3)%PartType)
           TheTree%PartRef(NumVertPart+2) = i4
           TheTree%PartType(NumVertPart+2) = ChargeConj(ThePrimAmp%IntPart(i4)%PartType)
           do NPart=0,NumVertPart-1
               TheTree%PartRef(NPart+2) = i3+NPart
               TheTree%PartType(NPart+2) = ExtParticle( ThePrimAmp%ExtLine(i3+NPart) )%PartType
           enddo

!          vertex 2
           TheTree => ThePrimAmp%UCuts(3)%TreeProcess(NCut,2)
           NumVertPart = i5-i4
           TheTree%NumPart = NumVertPart+2
           allocate( TheTree%PartRef(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartRef 3-2")
           allocate( TheTree%PartType(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartType 3-2")
!          set particle reference (wrt.prim.amp) and flavor
           TheTree%PartRef(1) = i4
           TheTree%PartType(1) = (ThePrimAmp%IntPart(i4)%PartType)
           TheTree%PartRef(NumVertPart+2) = i5
           TheTree%PartType(NumVertPart+2) = ChargeConj(ThePrimAmp%IntPart(i5)%PartType)
           do NPart=0,NumVertPart-1
               TheTree%PartRef(NPart+2) = i4+NPart
               TheTree%PartType(NPart+2) = ExtParticle( ThePrimAmp%ExtLine(i4+NPart) )%PartType
           enddo

!          vertex 3
           TheTree => ThePrimAmp%UCuts(3)%TreeProcess(NCut,3)
           NumVertPart = i3-(i5-NumExtParticles)
           TheTree%NumPart = NumVertPart+2
           allocate( TheTree%PartRef(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartRef 3-3")
           allocate( TheTree%PartType(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartType 3-3")
!          set particle reference (wrt.prim.amp) and flavor
           TheTree%PartRef(1) = i5
           TheTree%PartType(1) = (ThePrimAmp%IntPart(i5)%PartType)
           TheTree%PartRef(NumVertPart+2) = i3
           TheTree%PartType(NumVertPart+2) = ChargeConj(ThePrimAmp%IntPart(i3)%PartType)
           do NPart=0,NumVertPart-1
               if( i5+NPart.le.NumExtParticles ) then
                  TheTree%PartRef(NPart+2) = i5+NPart
                  TheTree%PartType(NPart+2) = ExtParticle( ThePrimAmp%ExtLine(i5+NPart) )%PartType
               else
                  TheTree%PartRef(NPart+2) = i5+NPart-NumExtParticles
                  TheTree%PartType(NPart+2) = ExtParticle( ThePrimAmp%ExtLine(i5+NPart-NumExtParticles) )%PartType
               endif
           enddo

           ThePrimAmp%UCuts(3)%CutProp(NCut,1) = i3
           ThePrimAmp%UCuts(3)%CutProp(NCut,2) = i4
           ThePrimAmp%UCuts(3)%CutProp(NCut,3) = i5

           NCut = NCut + 1
   enddo
   enddo
   enddo
   if ( NCut-1 .ne. ThePrimAmp%UCuts(3)%NumCuts ) call Error("Something went wrong while setting trip-cuts.")





!  set doub-cuts
   allocate(ThePrimAmp%UCuts(2)%TreeProcess(1:ThePrimAmp%UCuts(2)%NumCuts,1:2), stat=AllocStatus)
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in ThePrimAmp%UCuts(2)%TreeProcess(ThePrimAmp%UCuts(2)%NumCuts,1:2).")
   NCut = 1
   do i4 = 1,    NumExtParticles-1
   do i5 = i4+1, NumExtParticles

!          skip cuts that are marked for discarding
           if( ThePrimAmp%IntPart(i4)%PartType .eq. 0 ) cycle
           if( ThePrimAmp%IntPart(i5)%PartType .eq. 0 ) cycle

!          set tree processes
!          vertex 1
           TheTree => ThePrimAmp%UCuts(2)%TreeProcess(NCut,1)
           NumVertPart = i5-i4
           TheTree%NumPart = NumVertPart+2
           allocate( TheTree%PartRef(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartRef 2-1")
           allocate( TheTree%PartType(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartType 2-1")
!          set particle reference (wrt.prim.amp) and flavor
           TheTree%PartRef(1) = i4
           TheTree%PartType(1) = (ThePrimAmp%IntPart(i4)%PartType)
           TheTree%PartRef(NumVertPart+2) = i5
           TheTree%PartType(NumVertPart+2) = ChargeConj(ThePrimAmp%IntPart(i5)%PartType)
           do NPart=0,NumVertPart-1
               TheTree%PartRef(NPart+2) = i4+NPart
               TheTree%PartType(NPart+2) = ExtParticle( ThePrimAmp%ExtLine(i4+NPart) )%PartType
           enddo

!          check for massless external leg at vertex 1
           MasslessExtLeg = .false.
           if( NumVertPart.eq.1 .and. ExtParticle( TheTree%PartRef(2) )%Mass .le. 1d-10 ) then
               MasslessExtLeg = .true.
           endif


!          vertex 2
           TheTree => ThePrimAmp%UCuts(2)%TreeProcess(NCut,2)
           NumVertPart = i4-(i5-NumExtParticles)
           TheTree%NumPart = NumVertPart+2
           allocate( TheTree%PartRef(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartRef 2-2")
           allocate( TheTree%PartType(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartType 2-2")
!          set particle reference (wrt.prim.amp) and flavor
           TheTree%PartRef(1) = i5
           TheTree%PartType(1) = (ThePrimAmp%IntPart(i5)%PartType)
           TheTree%PartRef(NumVertPart+2) = i4
           TheTree%PartType(NumVertPart+2) = ChargeConj(ThePrimAmp%IntPart(i4)%PartType)
           do NPart=0,NumVertPart-1
               if( i5+NPart.le.NumExtParticles ) then
                  TheTree%PartRef(NPart+2) = i5+NPart
                  TheTree%PartType(NPart+2) = ExtParticle( ThePrimAmp%ExtLine(i5+NPart) )%PartType
               else
                  TheTree%PartRef(NPart+2) = i5+NPart-NumExtParticles
                  TheTree%PartType(NPart+2) = ExtParticle( ThePrimAmp%ExtLine(i5+NPart-NumExtParticles) )%PartType
               endif
           enddo

!          check for massless external leg at vertex 2
           if( NumVertPart.eq.1 .and. ExtParticle( TheTree%PartRef(2) )%Mass .le. 1d-10 ) then
               MasslessExtLeg = .true.
           endif

!          check for massless internal particles
           MasslessIntParticles = .false.
           if( ThePrimAmp%IntPart(i4)%Mass .le. 1d-10 .and. ThePrimAmp%IntPart(i5)%Mass .le. 1d-10 ) then
               MasslessIntParticles = .true.
           endif


!          check for massless bubbles
           if( MasslessExtLeg .and. MasslessIntParticles ) then
               deallocate( ThePrimAmp%UCuts(2)%TreeProcess(NCut,1)%PartRef )
               deallocate( ThePrimAmp%UCuts(2)%TreeProcess(NCut,1)%PartType )
               deallocate( ThePrimAmp%UCuts(2)%TreeProcess(NCut,2)%PartRef )
               deallocate( ThePrimAmp%UCuts(2)%TreeProcess(NCut,2)%PartType )
               ThePrimAmp%UCuts(2)%NumCuts = ThePrimAmp%UCuts(2)%NumCuts - 1
               cycle
           else
               ThePrimAmp%UCuts(2)%CutProp(NCut,1) = i4
               ThePrimAmp%UCuts(2)%CutProp(NCut,2) = i5
               NCut = NCut + 1
           endif


   enddo
   enddo
   if ( NCut-1 .ne. ThePrimAmp%UCuts(2)%NumCuts ) call Error("Something went wrong while setting doub-cuts.")




!  set sing-cuts
   allocate(ThePrimAmp%UCuts(1)%TreeProcess(1:ThePrimAmp%UCuts(1)%NumCuts,1:1), stat=AllocStatus)
   if( AllocStatus .ne. 0 ) call Error("Memory allocation in ThePrimAmp%UCuts(1)%TreeProcess(ThePrimAmp%UCuts(1)%NumCuts,1).")

   NCut = 1
   do i5 = 1, NumExtParticles

!          skip cuts that are marked for discarding
           if( ThePrimAmp%IntPart(i5)%PartType .eq. 0 ) cycle

!          set tree processes
!          vertex 1
           TheTree => ThePrimAmp%UCuts(1)%TreeProcess(NCut,1)
           NumVertPart = NumExtParticles
           TheTree%NumPart = NumVertPart+2
           allocate( TheTree%PartRef(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartRef 1-1")
           allocate( TheTree%PartType(1:NumVertPart+2), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%PartType 1-1")
!          set particle reference (wrt.prim.amp) and flavor
           TheTree%PartRef(1) = i5
           TheTree%PartType(1) = (ThePrimAmp%IntPart(i5)%PartType)
           TheTree%PartRef(NumVertPart+2) = i5
           TheTree%PartType(NumVertPart+2) = ChargeConj(ThePrimAmp%IntPart(i5)%PartType)
           do NPart=0,NumVertPart-1
               if( i5+NPart.le.NumExtParticles ) then
                  TheTree%PartRef(NPart+2) = i5+NPart
                  TheTree%PartType(NPart+2) = ExtParticle( ThePrimAmp%ExtLine(i5+NPart) )%PartType
               else
                  TheTree%PartRef(NPart+2) = i5+NPart-NumExtParticles
                  TheTree%PartType(NPart+2) = ExtParticle( ThePrimAmp%ExtLine(i5+NPart-NumExtParticles) )%PartType
               endif
           enddo

!          check for massless internal particles
           MasslessIntParticles = .false.
           if( ThePrimAmp%IntPart(i5)%Mass .le. 1d-10 ) then
               MasslessIntParticles = .true.
           endif

!          check for massless tadpole
           if( MasslessIntParticles ) then
               deallocate( ThePrimAmp%UCuts(1)%TreeProcess(NCut,1)%PartRef )
               deallocate( ThePrimAmp%UCuts(1)%TreeProcess(NCut,1)%PartType )
               ThePrimAmp%UCuts(1)%NumCuts = ThePrimAmp%UCuts(1)%NumCuts - 1
           else
!              set sing cut
               ThePrimAmp%UCuts(1)%CutProp(NCut,1) = i5
               NCut = NCut + 1
           endif
   enddo
   if ( NCut-1 .ne. ThePrimAmp%UCuts(1)%NumCuts ) call Error("Something went wrong while setting sing-cuts.")



! set cut matchings
  allocate( ThePrimAmp%UCuts(4)%Match(1:ThePrimAmp%UCuts(4)%NumCuts)  , stat=AllocStatus )
  if( AllocStatus .ne. 0 ) call Error("Memory allocation for quad cut matching")

  allocate( ThePrimAmp%UCuts(3)%Match(1:ThePrimAmp%UCuts(3)%NumCuts)  , stat=AllocStatus )
  if( AllocStatus .ne. 0 ) call Error("Memory allocation for trip cut matching")

  allocate( ThePrimAmp%UCuts(2)%Match(1:ThePrimAmp%UCuts(2)%NumCuts)  , stat=AllocStatus )
  if( AllocStatus .ne. 0 ) call Error("Memory allocation for doub cut matching")

  allocate( ThePrimAmp%UCuts(1)%Match(1:ThePrimAmp%UCuts(1)%NumCuts)  , stat=AllocStatus )
  if( AllocStatus .ne. 0 ) call Error("Memory allocation for sing cut matching")

  do NCut=1,ThePrimAmp%UCuts(4)%NumCuts
     call MatchUCuts(ThePrimAmp,4,5,NCut)
  enddo

  do NCut=1,ThePrimAmp%UCuts(3)%NumCuts
     call MatchUCuts(ThePrimAmp,3,5,NCut)
     call MatchUCuts(ThePrimAmp,3,4,NCut)
  enddo

  do NCut=1,ThePrimAmp%UCuts(2)%NumCuts
     call MatchUCuts(ThePrimAmp,2,5,NCut)
     call MatchUCuts(ThePrimAmp,2,4,NCut)
     call MatchUCuts(ThePrimAmp,2,3,NCut)
  enddo

  do NCut=1,ThePrimAmp%UCuts(1)%NumCuts
     call MatchUCuts(ThePrimAmp,1,5,NCut)
     call MatchUCuts(ThePrimAmp,1,4,NCut)
     call MatchUCuts(ThePrimAmp,1,3,NCut)
     call MatchUCuts(ThePrimAmp,1,2,NCut)
  enddo



! set tree particles
   do NPoint=1,5
      do NCut=1,ThePrimAmp%UCuts(NPoint)%NumCuts
         do NTree=1,NPoint

            TheTree => ThePrimAmp%UCuts(NPoint)%TreeProcess(NCut,NTree)
!           set number of quarks and gluons
            TheTree%NumQua = 0
            counter = 0
            do NPart=1,TheTree%NumPart
                  if( IsAQuark(TheTree%PartType(NPart)) ) then
                     TheTree%NumQua = TheTree%NumQua + 1
                     counter = counter + 1
                     QuarkPos(counter) = NPart
                  endif
            enddo

            if( IsAQuark(TheTree%PartType(1)) ) then
               allocate( TheTree%NumGlu(0:TheTree%NumQua), stat=AllocStatus )
               TheTree%NumGlu(0:TheTree%NumQua) = 0
            elseif( TheTree%PartType(1).eq.Glu_ ) then
               allocate( TheTree%NumGlu(0:TheTree%NumQua+1), stat=AllocStatus )
               TheTree%NumGlu(0:TheTree%NumQua) = 0
            else
               call Error("TheTree%NumGlu")
            endif
            if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%NumGlu")
            do NPart=1,TheTree%NumPart
                  if( TheTree%PartType(NPart) .eq. Glu_ ) then
                     TheTree%NumGlu(0) = TheTree%NumGlu(0) + 1
                  endif
            enddo

!           set number of gluons between quark lines
            if( IsAQuark(TheTree%PartType(1)) ) then
            if( TheTree%NumQua .eq. 2 ) then
                  TheTree%NumGlu(1) = QuarkPos(2) - QuarkPos(1) - 1
                  TheTree%NumGlu(2) = TheTree%NumPart - QuarkPos(2)
            endif
            if( TheTree%NumQua .eq. 4 ) then
                  TheTree%NumGlu(1) = QuarkPos(2) - QuarkPos(1) - 1
                  TheTree%NumGlu(2) = QuarkPos(3) - QuarkPos(2) - 1
                  TheTree%NumGlu(3) = QuarkPos(4) - QuarkPos(3) - 1
                  TheTree%NumGlu(4) = TheTree%NumPart - QuarkPos(4)
            endif
            if( TheTree%NumQua .eq. 6 ) then
                  TheTree%NumGlu(1) = QuarkPos(2) - QuarkPos(1) - 1
                  TheTree%NumGlu(2) = QuarkPos(3) - QuarkPos(2) - 1
                  TheTree%NumGlu(3) = QuarkPos(4) - QuarkPos(3) - 1
                  TheTree%NumGlu(4) = QuarkPos(5) - QuarkPos(4) - 1
                  TheTree%NumGlu(5) = QuarkPos(6) - QuarkPos(5) - 1
                  TheTree%NumGlu(6) = TheTree%NumPart - QuarkPos(6)
            endif
            elseif( TheTree%PartType(1).eq.Glu_ ) then
            if( TheTree%NumQua .eq. 2 ) then
                  TheTree%NumGlu(1) = QuarkPos(1) - 2
                  TheTree%NumGlu(2) = QuarkPos(2) - QuarkPos(1) - 1
                  TheTree%NumGlu(3) = TheTree%NumPart - QuarkPos(2)
            endif
            if( TheTree%NumQua .eq. 4 ) then
                  TheTree%NumGlu(1) = QuarkPos(1) - 2
                  TheTree%NumGlu(2) = QuarkPos(2) - QuarkPos(1) - 1
                  TheTree%NumGlu(3) = QuarkPos(3) - QuarkPos(2) - 1
                  TheTree%NumGlu(4) = QuarkPos(4) - QuarkPos(3) - 1
                  TheTree%NumGlu(5) = TheTree%NumPart - QuarkPos(4)
            endif
            if( TheTree%NumQua .eq. 6 ) then
                  TheTree%NumGlu(1) = QuarkPos(1) - 2
                  TheTree%NumGlu(2) = QuarkPos(2) - QuarkPos(1) - 1
                  TheTree%NumGlu(3) = QuarkPos(3) - QuarkPos(2) - 1
                  TheTree%NumGlu(4) = QuarkPos(4) - QuarkPos(3) - 1
                  TheTree%NumGlu(5) = QuarkPos(5) - QuarkPos(4) - 1
                  TheTree%NumGlu(6) = QuarkPos(6) - QuarkPos(5) - 1
                  TheTree%NumGlu(7) = TheTree%NumPart - QuarkPos(6)
            endif
            endif


!          allocate memory for pointer to quarks
           allocate( TheTree%Quarks(1:TheTree%NumQua), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%Quarks")
!          allocate memory for pointer to gluons
           allocate( TheTree%Gluons(1:TheTree%NumGlu(0)), stat=AllocStatus )
           if( AllocStatus .ne. 0 ) call Error("Memory allocation in TheTree%Gluons")

           counterQ = 0
           counterG = 0

           do NPart=1,TheTree%NumPart
               if( IsAQuark(TheTree%PartType(NPart)) ) then
                     counterQ = counterQ + 1
                     if( NPart.eq.1 .or. NPart.eq.TheTree%NumPart) then    ! first and last particles are in the loop
                        TheTree%Quarks(counterQ)%PartType => TheTree%PartType(NPart)
                        TheTree%Quarks(counterQ)%ExtRef => ThePrimAmp%IntPart( TheTree%PartRef(NPart) )%ExtRef
                        TheTree%Quarks(counterQ)%Mass => ThePrimAmp%IntPart( TheTree%PartRef(NPart) )%Mass
                        TheTree%Quarks(counterQ)%Mass2 => ThePrimAmp%IntPart( TheTree%PartRef(NPart) )%Mass2
                        TheTree%Quarks(counterQ)%Helicity => Null()
                        TheTree%Quarks(counterQ)%Mom => Null()
                        TheTree%Quarks(counterQ)%Pol => Null()
                     else
                        TheTree%Quarks(counterQ)%PartType => TheTree%PartType(NPart)
                        TheTree%Quarks(counterQ)%ExtRef => ExtParticle( ThePrimAmp%ExtLine(TheTree%PartRef(NPart)) )%ExtRef
                        TheTree%Quarks(counterQ)%Mass => ExtParticle( ThePrimAmp%ExtLine(TheTree%PartRef(NPart)) )%Mass
                        TheTree%Quarks(counterQ)%Mass2 => ExtParticle( ThePrimAmp%ExtLine(TheTree%PartRef(NPart)) )%Mass2
                        TheTree%Quarks(counterQ)%Helicity => ExtParticle( ThePrimAmp%ExtLine(TheTree%PartRef(NPart)) )%Helicity
                        TheTree%Quarks(counterQ)%Mom => ExtParticle( ThePrimAmp%ExtLine(TheTree%PartRef(NPart)) )%Mom
                        TheTree%Quarks(counterQ)%Pol => ExtParticle( ThePrimAmp%ExtLine(TheTree%PartRef(NPart)) )%Pol
                     endif
               endif
               if( TheTree%PartType(NPart) .eq. Glu_ ) then
                     counterG = counterG + 1
                     if( NPart.eq.1 .or. NPart.eq.TheTree%NumPart) then    ! first and last particles are in the loop
                        TheTree%Gluons(counterG)%PartType => TheTree%PartType(NPart)
                        TheTree%Gluons(counterG)%ExtRef => ThePrimAmp%IntPart( TheTree%PartRef(NPart) )%ExtRef
                        TheTree%Gluons(counterG)%Mass => ThePrimAmp%IntPart( TheTree%PartRef(NPart) )%Mass
                        TheTree%Gluons(counterG)%Mass2 => ThePrimAmp%IntPart( TheTree%PartRef(NPart) )%Mass2
                        TheTree%Gluons(counterG)%Helicity => Null()
                        TheTree%Gluons(counterG)%Mom => Null()
                        TheTree%Gluons(counterG)%Pol => Null()
                     else
                        TheTree%Gluons(counterG)%PartType => TheTree%PartType(NPart)
                        TheTree%Gluons(counterG)%ExtRef => ExtParticle( ThePrimAmp%ExtLine(TheTree%PartRef(NPart)) )%ExtRef
                        TheTree%Gluons(counterG)%Mass => ExtParticle( ThePrimAmp%ExtLine(TheTree%PartRef(NPart)) )%Mass
                        TheTree%Gluons(counterG)%Mass2 => ExtParticle( ThePrimAmp%ExtLine(TheTree%PartRef(NPart)) )%Mass2
                        TheTree%Gluons(counterG)%Helicity => ExtParticle( ThePrimAmp%ExtLine(TheTree%PartRef(NPart)) )%Helicity
                        TheTree%Gluons(counterG)%Mom => ExtParticle( ThePrimAmp%ExtLine(TheTree%PartRef(NPart)) )%Mom
                        TheTree%Gluons(counterG)%Pol => ExtParticle( ThePrimAmp%ExtLine(TheTree%PartRef(NPart)) )%Pol
                     endif
               endif
           enddo

         enddo
      enddo
   enddo


END SUBROUTINE



SUBROUTINE MatchUCuts(PrimAmp,Cut,HiCut,CutNum)
use ModMisc
implicit none
type(PrimitiveAmplitude),target :: PrimAmp
integer :: Cut,CutNum,HiCut,NumMatch,MatchHiCuts(1:20),FirstHiProp(1:20),MissHiProp(1:20,1:4)
integer :: Prop(1:5),HiProp(1:5),INTER(1:5),COMPL(1:5),FirstInterPos
integer :: i,numcheck,AllocStatus
type(UCutMatch),pointer :: TheMatch

NumMatch=0
do i=1,PrimAmp%UCuts(HiCut)%NumCuts
    Prop(1:Cut) = PrimAmp%UCuts(Cut)%CutProp(CutNum,1:Cut)
    HiProp(1:HiCut) = PrimAmp%UCuts(HiCut)%CutProp(i,1:HiCut)
    numcheck= MatchSets(HiProp(1:HiCut),Prop(1:Cut),INTER,COMPL,FirstInterPos)
    if( numcheck.eq.Cut ) then
        NumMatch=NumMatch+1
        MatchHiCuts(NumMatch)=i
        FirstHiProp(NumMatch)=FirstInterPos
        MissHiProp(NumMatch,1:HiCut-Cut)=COMPL(1:HiCut-Cut)
    endif
enddo

  TheMatch => PrimAmp%UCuts(Cut)%Match(CutNum)

  TheMatch%Subt(HiCut)%NumMatch = NumMatch

  allocate(TheMatch%Subt(HiCut)%MatchHiCuts(1:NumMatch), stat=AllocStatus )
  if( AllocStatus .ne. 0 ) call Error("Memory allocation in MatchUCuts 1")
  TheMatch%Subt(HiCut)%MatchHiCuts(1:NumMatch) = MatchHiCuts(1:NumMatch)

  allocate(TheMatch%Subt(HiCut)%FirstHiProp(1:NumMatch), stat=AllocStatus )
  if( AllocStatus .ne. 0 ) call Error("Memory allocation in MatchUCuts 2")
  TheMatch%Subt(HiCut)%FirstHiProp(1:NumMatch) = FirstHiProp(1:NumMatch)

  allocate(TheMatch%Subt(HiCut)%MissHiProp(1:NumMatch,1:HiCut-Cut), stat=AllocStatus )
  if( AllocStatus .ne. 0 ) call Error("Memory allocation in MatchUCuts 3")
  TheMatch%Subt(HiCut)%MissHiProp(1:NumMatch,1:HiCut-Cut) = MissHiProp(1:NumMatch,1:HiCut-Cut)

RETURN
END SUBROUTINE





SUBROUTINE InfoPrimAmps(Filename)
use ModMisc
implicit none
type(PrimitiveAmplitude),pointer :: ThePrimAmp
integer :: NPrimAmp,NCut,NPart
character :: Filename*(*)


   open(unit=14,file=Filename,form='formatted',access= 'sequential',status='replace')


   write(14,"(A31,I3)") "Number of primitive amplitudes:",NumPrimAmps
   write(14,*) ""

   do NPrimAmp=1,NumBornAmps
      write(14,"(A31,10I3)") "LO amplitude: ",BornAmps(NPrimAmp)%TreeProc%PartType(1:NumExtParticles)
      write(14,"(A32,A100)") " ", WritePartType( BornAmps(NPrimAmp)%TreeProc%PartType(1:NumExtParticles) )
      write(14,*) ""
   enddo

   do NPrimAmp=1,NumPrimAmps
      ThePrimAmp => PrimAmps(NPrimAmp)
      write(14,"(A)") "----------------------------------------------------------------------------"
      write(14,"(A31,I3)") "Primitive amplitude: ",NPrimAmp
      if(ThePrimAmp%AmpType .eq.1) then
         write(14,"(A31,I3,A30)") "Amplitude type: ",ThePrimAmp%AmpType," (class (a) diagram)"
      elseif(ThePrimAmp%AmpType .eq.2) then
         write(14,"(A31,I3,A30)") "Amplitude type: ",ThePrimAmp%AmpType," (fermion loop diagram)"
      elseif(ThePrimAmp%AmpType .eq.3) then
         write(14,"(A31,I3,A30)") "Amplitude type: ",ThePrimAmp%AmpType," (class (b) diagram)"
      elseif(ThePrimAmp%AmpType .eq.4) then
         write(14,"(A31,I3,A30)") "Amplitude type: ",ThePrimAmp%AmpType," (class (c) diagram)"
      else
         write(14,"(A31,I3)") "Amplitude type: ",ThePrimAmp%AmpType
      endif

      write(14,"(A31,I3)") "N-Point: ",ThePrimAmp%NPoint
      write(14,"(A31,10I3)") "External ordering: ",ThePrimAmp%ExtLine(1:NumExtParticles)
      write(14,"(A32,A100)") " ", WritePartType( ExtParticle(ThePrimAmp%ExtLine(1:NumExtParticles))%PartType )
      write(14,"(A31,10I3)") "Internal ordering: ",ThePrimAmp%IntPart(1:NumExtParticles)%PartType
      write(14,"(A32,A100)") " ", WritePartType( ThePrimAmp%IntPart(1:NumExtParticles)%PartType )
      write(14,*) ""
      write(14,"(A)") "--------------------------------------"

      write(14,"(A31,I3)") "Number of pent-cuts: ",ThePrimAmp%UCuts(5)%NumCuts
      do NCut=1,ThePrimAmp%UCuts(5)%NumCuts
         write(14,"(A31,I3,A1,5I3)") "cut (pointer)",NCut,":",ThePrimAmp%UCuts(5)%CutProp(NCut,1:5)
         NPart = ThePrimAmp%UCuts(5)%TreeProcess(NCut,1)%NumPart
         write(14,"(A40,10I3)") "Tree 1:",ThePrimAmp%UCuts(5)%TreeProcess(NCut,1)%PartRef(1:NPart)
         write(14,"(A41,A100)") "       ",WritePartType( ThePrimAmp%UCuts(5)%TreeProcess(NCut,1)%PartType(1:NPart) )
         NPart = ThePrimAmp%UCuts(5)%TreeProcess(NCut,2)%NumPart
         write(14,"(A40,10I3)") "Tree 2:",ThePrimAmp%UCuts(5)%TreeProcess(NCut,2)%PartRef(1:NPart)
         write(14,"(A41,A100)") "       ",WritePartType( ThePrimAmp%UCuts(5)%TreeProcess(NCut,2)%PartType(1:NPart) )
         NPart = ThePrimAmp%UCuts(5)%TreeProcess(NCut,3)%NumPart
         write(14,"(A40,10I3)") "Tree 3:",ThePrimAmp%UCuts(5)%TreeProcess(NCut,3)%PartRef(1:NPart)
         write(14,"(A41,A100)") "       ",WritePartType( ThePrimAmp%UCuts(5)%TreeProcess(NCut,3)%PartType(1:NPart) )
         NPart = ThePrimAmp%UCuts(5)%TreeProcess(NCut,4)%NumPart
         write(14,"(A40,10I3)") "Tree 4:",ThePrimAmp%UCuts(5)%TreeProcess(NCut,4)%PartRef(1:NPart)
         write(14,"(A41,A100)") "       ",WritePartType( ThePrimAmp%UCuts(5)%TreeProcess(NCut,4)%PartType(1:NPart) )
         NPart = ThePrimAmp%UCuts(5)%TreeProcess(NCut,5)%NumPart
         write(14,"(A40,10I3)") "Tree 5:",ThePrimAmp%UCuts(5)%TreeProcess(NCut,5)%PartRef(1:NPart)
         write(14,"(A41,A100)") "       ",WritePartType( ThePrimAmp%UCuts(5)%TreeProcess(NCut,5)%PartType(1:NPart) )
      enddo

      write(14,"(A31,I3)") "Number of quad-cuts: ",ThePrimAmp%UCuts(4)%NumCuts
      do NCut=1,ThePrimAmp%UCuts(4)%NumCuts
         write(14,"(A31,I3,A1,4I3)") "cut (pointer)",NCut,":",ThePrimAmp%UCuts(4)%CutProp(NCut,1:4)
         NPart = ThePrimAmp%UCuts(4)%TreeProcess(NCut,1)%NumPart
         write(14,"(A40,10I3)") "Tree 1:",ThePrimAmp%UCuts(4)%TreeProcess(NCut,1)%PartRef(1:NPart)
         write(14,"(A41,A100)") "       ",WritePartType( ThePrimAmp%UCuts(4)%TreeProcess(NCut,1)%PartType(1:NPart) )
         NPart = ThePrimAmp%UCuts(4)%TreeProcess(NCut,2)%NumPart
         write(14,"(A40,10I3)") "Tree 2:",ThePrimAmp%UCuts(4)%TreeProcess(NCut,2)%PartRef(1:NPart)
         write(14,"(A41,A100)") "       ",WritePartType( ThePrimAmp%UCuts(4)%TreeProcess(NCut,2)%PartType(1:NPart) )
         NPart = ThePrimAmp%UCuts(4)%TreeProcess(NCut,3)%NumPart
         write(14,"(A40,10I3)") "Tree 3:",ThePrimAmp%UCuts(4)%TreeProcess(NCut,3)%PartRef(1:NPart)
         write(14,"(A41,A100)") "       ",WritePartType( ThePrimAmp%UCuts(4)%TreeProcess(NCut,3)%PartType(1:NPart) )
         NPart = ThePrimAmp%UCuts(4)%TreeProcess(NCut,4)%NumPart
         write(14,"(A40,10I3)") "Tree 4:",ThePrimAmp%UCuts(4)%TreeProcess(NCut,4)%PartRef(1:NPart)
         write(14,"(A41,A100)") "       ",WritePartType( ThePrimAmp%UCuts(4)%TreeProcess(NCut,4)%PartType(1:NPart) )
      enddo

      write(14,"(A31,I3)") "Number of trip-cuts: ",ThePrimAmp%UCuts(3)%NumCuts
      do NCut=1,ThePrimAmp%UCuts(3)%NumCuts
         write(14,"(A31,I3,A1,4I3)") "cut (pointer)",NCut,":",ThePrimAmp%UCuts(3)%CutProp(NCut,1:3)
         NPart = ThePrimAmp%UCuts(3)%TreeProcess(NCut,1)%NumPart
         write(14,"(A40,10I3)") "Tree 1:",ThePrimAmp%UCuts(3)%TreeProcess(NCut,1)%PartRef(1:NPart)
         write(14,"(A41,A100)") "       ",WritePartType( ThePrimAmp%UCuts(3)%TreeProcess(NCut,1)%PartType(1:NPart) )
         NPart = ThePrimAmp%UCuts(3)%TreeProcess(NCut,2)%NumPart
         write(14,"(A40,10I3)") "Tree 2:",ThePrimAmp%UCuts(3)%TreeProcess(NCut,2)%PartRef(1:NPart)
         write(14,"(A41,A100)") "       ",WritePartType( ThePrimAmp%UCuts(3)%TreeProcess(NCut,2)%PartType(1:NPart) )
         NPart = ThePrimAmp%UCuts(3)%TreeProcess(NCut,3)%NumPart
         write(14,"(A40,10I3)") "Tree 3:",ThePrimAmp%UCuts(3)%TreeProcess(NCut,3)%PartRef(1:NPart)
         write(14,"(A41,A100)") "       ",WritePartType( ThePrimAmp%UCuts(3)%TreeProcess(NCut,3)%PartType(1:NPart) )
      enddo

      write(14,"(A31,I3)") "Number of doub-cuts: ",ThePrimAmp%UCuts(2)%NumCuts
      do NCut=1,ThePrimAmp%UCuts(2)%NumCuts
         write(14,"(A31,I3,A1,4I3)") "cut (pointer)",NCut,":",ThePrimAmp%UCuts(2)%CutProp(NCut,1:2)
         NPart = ThePrimAmp%UCuts(2)%TreeProcess(NCut,1)%NumPart
         write(14,"(A40,10I3)") "Tree 1:",ThePrimAmp%UCuts(2)%TreeProcess(NCut,1)%PartRef(1:NPart)
         write(14,"(A41,A100)") "       ",WritePartType( ThePrimAmp%UCuts(2)%TreeProcess(NCut,1)%PartType(1:NPart) )
         NPart = ThePrimAmp%UCuts(2)%TreeProcess(NCut,2)%NumPart
         write(14,"(A40,10I3)") "Tree 2:",ThePrimAmp%UCuts(2)%TreeProcess(NCut,2)%PartRef(1:NPart)
         write(14,"(A41,A100)") "       ",WritePartType( ThePrimAmp%UCuts(2)%TreeProcess(NCut,2)%PartType(1:NPart) )
      enddo

      write(14,"(A31,I3)") "Number of sing-cuts: ",ThePrimAmp%UCuts(1)%NumCuts
      do NCut=1,ThePrimAmp%UCuts(1)%NumCuts
         write(14,"(A31,I3,A1,4I3)") "cut (pointer)",NCut,":",ThePrimAmp%UCuts(1)%CutProp(NCut,1:1)
         NPart = ThePrimAmp%UCuts(1)%TreeProcess(NCut,1)%NumPart
         write(14,"(A40,10I3)") "Tree 1:",ThePrimAmp%UCuts(1)%TreeProcess(NCut,1)%PartRef(1:NPart)
         write(14,"(A41,A100)") "       ",WritePartType( ThePrimAmp%UCuts(1)%TreeProcess(NCut,1)%PartType(1:NPart) )
      enddo
   enddo
   write(14,*) ""
   close(14)

END SUBROUTINE







END MODULE
