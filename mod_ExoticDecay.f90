MODULE ModExoticDecay
implicit none
save

public :: HTopA0Decay ! heavy top --> A0(scalar)+top
public :: HTopBHDecay ! heavy top --> BS(vector)+top
public :: StopDecay   ! stop --> Chi0 + top

contains






SUBROUTINE HTopA0Decay(HeavyTop,Topol,Mom)
use ModMisc
use ModProcess
use ModParameters
use ModTopDecay
implicit none
type(Particle) :: HeavyTop,TopQuark
integer :: Topol,NMom,i
real(8) :: NWAFactor_HTop,zeros(1:7)
real(8) :: Mom(:,:)! order: A0,t,b,lep,neu
integer,parameter :: A0=1,top=2,bot=3,lep=4,neu=5



!DEC$ IF(_CheckMomenta .EQ.1)
   NMom = size(Mom(:,:),2)
   zeros(1:4) = dble(HeavyTop%Mom(1:4)) - Mom(1:4,A0) - Mom(1:4,top)

   if( any(abs(zeros(1:4)/dble(HeavyTop%Mom(1))).gt.1d-8) ) then
      print *, "ERROR: energy-momentum violation in SUBROUTINE HTopA0Decay(): ",HeavyTop%PartType,zeros(1:4)
      print *, "momentum dump:"
      print *, Mom(:,:)
      pause
   endif

   zeros(:) = 0d0
   zeros(1) = dble(HeavyTop%Mom(1:4).dot.HeavyTop%Mom(1:4)) - m_HTop**2
   zeros(2) = (Mom(1:4,A0).dot.Mom(1:4,A0)) - m_A0**2
   zeros(3) = (Mom(1:4,top).dot.Mom(1:4,top)) - m_Top**2
   do i=3,NMom
      zeros(i+1)=  (Mom(1:4,i).dot.Mom(1:4,i))
   enddo

   if( any(abs(zeros(1:NMom+1)/dble(HeavyTop%Mom(1))**2).gt.1d-5) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE HTopA0Decay(): ",HeavyTop%PartType,zeros(1:NMom+1)
      print *, "momentum dump:"
      print *, Mom(:,:)
      pause
   endif
!DEC$ ENDIF



NWAFactor_HTop = 1d0/dsqrt(2d0*Ga_HTop(0)*M_HTop)
!-----------------------------------------------------------------------------
IF( Topol.eq.DKX_HTA0_LO ) THEN! leading order
!-----------------------------------------------------------------------------
    if( HeavyTop%PartType.gt.0 ) then! HTop quark decay, remember: PartType=Top and ATop for HTop and AHTop!!
          TopQuark%PartType = Top_
          TopQuark%Mass = m_Top
          TopQuark%Mass2 = m_Top**2
          TopQuark%Mom(1:4) = dcmplx(Mom(1:4,bot)+Mom(1:4,lep)+Mom(1:4,neu))
          call TopDecay(TopQuark,DK_LO,Mom(1:4,bot:neu))
          HeavyTop%Pol(1:4) = ( IA0Tt(+1)*Chir(.true.,TopQuark%Pol(1:4)) + IA0Tt(-1)*Chir(.false.,TopQuark%Pol(1:4)) ) * (0d0,1d0)
  
          HeavyTop%Pol(1:4) =( spb2_(HeavyTop%Pol(1:4),HeavyTop%Mom(1:4)) + M_HTop*HeavyTop%Pol(1:4) ) * NWAFactor_HTop
          HeavyTop%Pol(5:16)= (0d0,0d0)

    else! anti HTop quark decay
          TopQuark%PartType = ATop_
          TopQuark%Mass = m_Top
          TopQuark%Mass2 = m_Top**2
          TopQuark%Mom(1:4) = dcmplx(Mom(1:4,bot)+Mom(1:4,lep)+Mom(1:4,neu))
          call TopDecay(TopQuark,DK_LO,Mom(1:4,bot:neu))
          HeavyTop%Pol(1:4) = ( IA0Tt(+1)*Chir(.false.,TopQuark%Pol(1:4)) + IA0Tt(-1)*Chir(.true.,TopQuark%Pol(1:4)) ) * (0d0,1d0)

          HeavyTop%Pol(1:4) = ( spi2_(HeavyTop%Mom(1:4),HeavyTop%Pol(1:4)) - M_HTop*HeavyTop%Pol(1:4) ) * NWAFactor_HTop
          HeavyTop%Pol(5:16)= (0d0,0d0)
    endif
ENDIF


END SUBROUTINE









SUBROUTINE HTopBHDecay(HeavyTop,Topol,BHHel,Mom,MomGlu,HelGlu)
use ModMisc
use ModProcess
use ModParameters
use ModTopDecay
implicit none
type(Particle) :: HeavyTop,TopQuark
integer :: Topol,NMom,BHHel,i
real(8),optional :: MomGlu(1:4)
integer,optional :: HelGlu
real(8) :: NWAFactor_HTop,zeros(1:7)
real(8) :: Mom(:,:)! order: BH,t,b,lep,neu
complex(8) :: BarSpi(1:4,-1,1), Spi(1:4,-1,+1),BHPol(1:4)
integer,parameter :: BH=1,top=2,bot=3,lep=4,neu=5



!DEC$ IF(_CheckMomenta .EQ.1)
   NMom = size(Mom(:,:),2)
   zeros(1:4) = dble(HeavyTop%Mom(1:4)) - Mom(1:4,BH) - Mom(1:4,top)

   if( any(abs(zeros(1:4)/dble(HeavyTop%Mom(1))).gt.1d-8) ) then
      print *, "ERROR: energy-momentum violation in SUBROUTINE HTopBHDecay(): ",HeavyTop%PartType,zeros(1:4)
      print *, "momentum dump:"
      print *, Mom(:,:)
      pause
   endif

   zeros(:) = 0d0
   zeros(1) = dble(HeavyTop%Mom(1:4).dot.HeavyTop%Mom(1:4)) - m_HTop**2
   zeros(2) = (Mom(1:4,BH).dot.Mom(1:4,BH)) - m_BH**2
   zeros(3) = (Mom(1:4,top).dot.Mom(1:4,top)) - m_Top**2
   do i=3,NMom
      zeros(i+1)=  (Mom(1:4,i).dot.Mom(1:4,i))
   enddo

   if( any(abs(zeros(1:NMom+1)/dble(HeavyTop%Mom(1))**2).gt.1d-5) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE HTopBHDecay(): ",HeavyTop%PartType,zeros(1:NMom+1)
      print *, "momentum dump:"
      print *, Mom(:,:)
      pause
   endif
!DEC$ ENDIF



NWAFactor_HTop = 1d0/dsqrt(2d0*Ga_HTop(0)*M_HTop)
!-----------------------------------------------------------------------------
IF( Topol.eq.DKX_HTBH_LO ) THEN! leading order
!-----------------------------------------------------------------------------

    BHPol(1:4) = pol_mass(dcmplx(Mom(1:4,BH)),m_BH,BHHel)
    if( HeavyTop%PartType.gt.0 ) then! HTop quark decay, remember: PartType=Top and ATop for HTop and AHTop!!
          TopQuark%PartType = Top_
          TopQuark%Mass = m_Top
          TopQuark%Mass2= m_Top**2
          TopQuark%Mom(1:4) = dcmplx(Mom(1:4,bot)+Mom(1:4,lep)+Mom(1:4,neu))
          call TopDecay(TopQuark,DK_LO,Mom(1:4,bot:neu))
          HeavyTop%Pol(1:4) = spb2_(TopQuark%Pol(1:4),BHPol(1:4))
          HeavyTop%Pol(1:4) = ( IBHTt(+1)*Chir(.true.,HeavyTop%Pol(1:4)) + IBHTt(-1)*Chir(.false.,HeavyTop%Pol(1:4)) ) * (0d0,1d0)
          HeavyTop%Pol(1:4) =( spb2_(HeavyTop%Pol(1:4),HeavyTop%Mom(1:4)) + M_HTop*HeavyTop%Pol(1:4) ) * NWAFactor_HTop
          HeavyTop%Pol(5:16)= (0d0,0d0)
    else! anti HTop quark decay
          TopQuark%PartType = ATop_
          TopQuark%Mass = m_Top
          TopQuark%Mass2= m_Top**2
          TopQuark%Mom(1:4) = dcmplx(Mom(1:4,bot)+Mom(1:4,lep)+Mom(1:4,neu))
          call TopDecay(TopQuark,DK_LO,Mom(1:4,bot:neu))
          HeavyTop%Pol(1:4) = ( IBHTt(-1)*Chir(.false.,TopQuark%Pol(1:4)) + IBHTt(+1)*Chir(.true.,TopQuark%Pol(1:4)) ) * (0d0,1d0)
          HeavyTop%Pol(1:4) = spi2_(BHPol(1:4),HeavyTop%Pol(1:4))
          HeavyTop%Pol(1:4) = ( spi2_(HeavyTop%Mom(1:4),HeavyTop%Pol(1:4)) - M_HTop*HeavyTop%Pol(1:4) ) * NWAFactor_HTop
          HeavyTop%Pol(5:16)= (0d0,0d0)
    endif
ENDIF


END SUBROUTINE










SUBROUTINE StopDecay(StopQuark,Topol,ChiHel,Mom,MomGlu,HelGlu,HelTop)
use ModMisc
use ModProcess
use ModParameters
use ModTopDecay
implicit none
type(Particle) :: StopQuark
type(Particle) :: TopQuark
integer :: Topol,ChiHel,i,epMI
real(8) :: Mom(:,:)! chi top bot lep neu
real(8),optional :: MomGlu(1:4)
integer,optional :: HelGlu
integer,optional :: HelTop! this is for stable tops
complex(8) :: Spi(1:4),BarSpi(1:4),Diagram1,Diagram2,PolGlu(1:4)
complex(8) :: SpiPlus(1:4),SpiMinus(1:4),MI(1:4),FF1,FF2,ID(-1:0)
real(8) :: NWAFactor_STop,zeros(1:7),MomIn(1:4,1:4),Norm
real(8) :: beta,omega,P3,P0,Ppl,Pmi,Yp,Yw,z,W0,Wpl,Wmi,v13,Catani
real(8),parameter :: CF=4d0/3d0
complex(8) :: qlI2,qlI3
if( XTopDecays.eq.0 ) then
  StopQuark%Pol(1)    = (1d0,0d0)
  StopQuark%Pol(2:16) = (0d0,0d0)
  return
endif



!DEC$ IF(_CheckMomenta .EQ.1)
   zeros(:) = 0d0
   zeros(1:4) = dble(STopQuark%Mom(1:4)) - Mom(1:4,1) - Mom(1:4,3) - Mom(1:4,4) - Mom(1:4,5)
   if( present(MomGlu) ) zeros(1:4) = zeros(1:4) - MomGlu(1:4)
   if( any(abs(zeros(1:4)/dble(STopQuark%Mom(1))).gt.1d-8) ) then
      print *, "ERROR: energy-momentum violation in SUBROUTINE STopDecay(): ",STopQuark%PartType,zeros(1:4)
      print *, "momentum dump:"
      print *, Mom(:,:)
   endif

   zeros(1) = dble(STopQuark%Mom(1:4).dot.STopQuark%Mom(1:4)) - m_STop**2
   zeros(2) = (Mom(1:4,1).dot.Mom(1:4,1)) - m_Chi**2
   zeros(3) = (Mom(1:4,2).dot.Mom(1:4,2)) - m_Top**2
   zeros(4) =  Mom(1:4,3).dot.Mom(1:4,3)
   zeros(5) =  Mom(1:4,4).dot.Mom(1:4,4)
   zeros(6) =  Mom(1:4,5).dot.Mom(1:4,5)
   if( present(MomGlu) ) zeros(7) = MomGlu(1:4).dot.MomGlu(1:4)
   if( any(abs(zeros(1:7)/dble(STopQuark%Mom(1))**2).gt.1d-5) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE STopDecay(): ",STopQuark%PartType,zeros(1:7)
      print *, "momentum dump:"
      print *, Mom(:,:)
   endif
!DEC$ ENDIF






   NWAFactor_STop = 1d0/dsqrt(2d0*Ga_STop(0)*M_STop)


IF( Topol.eq.DKX_STChi0_LO ) THEN! leading order

   TopQuark%PartType= sign(1,StopQuark%PartType) * Top_
   TopQuark%Mass = m_Top
   TopQuark%Mass2= m_Top**2
   TopQuark%Mom(1:4)=dcmplx(Mom(1:4,2))
   if( present(HelTop) ) then! stable top
      TopQuark%Helicity = HelTop
      if( TopQuark%PartType.gt.0 ) then 
          call ubarSpi(TopQuark%Mom(1:4),TopQuark%Mass,TopQuark%Helicity,TopQuark%Pol(1:4))
      else
          call vSpi(TopQuark%Mom(1:4),TopQuark%Mass,TopQuark%Helicity,TopQuark%Pol(1:4))
      endif
   else
      call TopDecay(TopQuark,DK_LO,Mom(1:4,3:5))
   endif

   if( StopQuark%PartType.eq.Stop_ ) then
      call vMajoSpi(dcmplx(Mom(1:4,1)),m_Chi,ChiHel,Spi(1:4))
      BarSpi(1:4) = TopQuark%Pol(1:4)
      Spi(1:4) = (0d0,1d0)*( IChiStt(+1)*Chir(.true.,Spi(1:4)) + IChiStt(-1)*Chir(.false.,Spi(1:4)) )

! StopQuark%Pol(1) = 1d0
! RETURN

   elseif( StopQuark%PartType.eq.AStop_ ) then
      call ubarMajoSpi(dcmplx(Mom(1:4,1)),m_Chi,ChiHel,BarSpi(1:4))
      Spi(1:4) = TopQuark%Pol(1:4)
      Spi(1:4) = (0d0,1d0)*( IChiStt(-1)*Chir(.true.,Spi(1:4)) + IChiStt(+1)*Chir(.false.,Spi(1:4)) )

! StopQuark%Pol(1) = 1d0
! RETURN

   endif
   StopQuark%Pol(1) = psp1_(BarSpi(1:4),Spi(1:4)) * NWAFactor_STop


! print *, "check", StopQuark%PartType
! print *,  barspi(1:4)
! print *, spi(1:4)
! print *, "check"
! call uMajoSpi_Weyl(dcmplx(-Mom(1:4,4)),+1,Spi(1:4))
! call ubarMajoSpi_Weyl(dcmplx(Mom(1:4,5)),+1,BarSpi(1:4))
! print *, vbqq_Weyl(4,BarSpi,Spi)
! call vMajoSpi_Weyl(dcmplx(Mom(1:4,5)),+1,Spi(1:4))
! call vbarMajoSpi_Weyl(dcmplx(-Mom(1:4,4)),+1,BarSpi(1:4))
! print *, vbqq_Weyl(4,BarSpi,Spi)
!pause






ELSEIF( Topol.eq.DKX_STChi0_RE1 ) THEN! gluon emission off (anti)stop

   TopQuark%PartType= sign(1,StopQuark%PartType) * Top_
   TopQuark%Mass = m_Top
   TopQuark%Mass2= m_Top**2
   TopQuark%Mom(1:4)=dcmplx(Mom(1:4,2))
   call TopDecay(TopQuark,DK_LO,Mom(1:4,3:5))
   if( present(HelTop) ) then! stable top
      TopQuark%Helicity = HelTop
      if( TopQuark%PartType.gt.0 ) then 
          call ubarSpi(TopQuark%Mom(1:4),TopQuark%Mass,TopQuark%Helicity,TopQuark%Pol(1:4))
      else
          call vSpi(TopQuark%Mom(1:4),TopQuark%Mass,TopQuark%Helicity,TopQuark%Pol(1:4))
      endif
   endif

   call pol_mless(dcmplx(MomGlu(1:4)),HelGlu,PolGlu(1:4))
!    PolGlu(1:4)= MomGlu(1:4); print *, "Stop decay gauge invariance check"


   if( StopQuark%PartType.eq.Stop_ ) then
      call vMajoSpi(dcmplx(Mom(1:4,1)),m_Chi,ChiHel,Spi(1:4))
      BarSpi(1:4) = TopQuark%Pol(1:4)
      Spi(1:4) = (0d0,1d0)*( IChiStt(+1)*Chir(.true.,Spi(1:4)) + IChiStt(-1)*Chir(.false.,Spi(1:4)) )
      Diagram1 = psp1_(BarSpi(1:4),Spi(1:4)) * cgs(PolGlu(1:4),dcmplx(MomGlu(1:4)),StopQuark%Mom(1:4)-MomGlu(1:4)) * sqrt2! sqrt2 removes 1/sqrt2 from cgs
      Diagram1 = Diagram1 * (0d0,1d0)/( ((StopQuark%Mom(1:4)-MomGlu(1:4)).dot.(StopQuark%Mom(1:4)-MomGlu(1:4))) -m_stop**2 )

      BarSpi(1:4) = vgq(PolGlu(1:4),TopQuark%Pol(1:4)) * sqrt2! sqrt2 removes 1/sqrt2 from vgq
      BarSpi(1:4) = ( spb2_(BarSpi(1:4),TopQuark%Mom(1:4)+MomGlu(1:4)) + m_top*BarSpi(1:4) ) * (0d0,1d0)/( ((TopQuark%Mom(1:4)+MomGlu(1:4)).dot.(TopQuark%Mom(1:4)+MomGlu(1:4))) -m_top**2 )
      Diagram2 = psp1_(BarSpi(1:4),Spi(1:4))


   elseif( StopQuark%PartType.eq.AStop_ ) then
      call ubarMajoSpi(dcmplx(Mom(1:4,1)),m_Chi,ChiHel,BarSpi(1:4))
      Spi(1:4) = TopQuark%Pol(1:4)
      Spi(1:4) = (0d0,1d0)*( IChiStt(+1)*Chir(.true.,Spi(1:4)) + IChiStt(-1)*Chir(.false.,Spi(1:4)) )
      Diagram1 = psp1_(BarSpi(1:4),Spi(1:4)) * cgbs(PolGlu(1:4),dcmplx(MomGlu(1:4)),StopQuark%Mom(1:4)-MomGlu(1:4)) * sqrt2! sqrt2 removes 1/sqrt2 from cgbs
      Diagram1 = Diagram1 * (0d0,1d0)/( ((StopQuark%Mom(1:4)-MomGlu(1:4)).dot.(StopQuark%Mom(1:4)-MomGlu(1:4))) -m_stop**2 )

      Spi(1:4) = vgbq(PolGlu(1:4),TopQuark%Pol(1:4)) * sqrt2! sqrt2 removes 1/sqrt2 from vgbq
      Spi(1:4) = ( -spi2_(TopQuark%Mom(1:4)+MomGlu(1:4),Spi(1:4)) + m_top*Spi(1:4) ) * (0d0,1d0)/( ((TopQuark%Mom(1:4)+MomGlu(1:4)).dot.(TopQuark%Mom(1:4)+MomGlu(1:4))) -m_top**2 )
      Spi(1:4) = (0d0,1d0)*( IChiStt(+1)*Chir(.true.,Spi(1:4)) + IChiStt(-1)*Chir(.false.,Spi(1:4)) )
      Diagram2 = psp1_(BarSpi(1:4),Spi(1:4))

   endif

   StopQuark%Pol(1) = ( Diagram1 + Diagram2 ) * NWAFactor_STop * dsqrt( alpha_s4Pi*RunAlphaS(NLOParam,MuRen) * CF )


ELSEIF( Topol.eq.DKX_STChi0_RE2 ) THEN! gluon emission off (anti)top

   TopQuark%PartType= sign(1,StopQuark%PartType) * Top_
   TopQuark%Mass = m_Top
   TopQuark%Mass2= m_Top**2
   TopQuark%Mom(1:4)=dcmplx(Mom(1:4,2))
   MomIn(1:4,1:3) = Mom(1:4,3:5)
   MomIn(1:4,4)   = MomGlu(1:4)
   call TopDecay(TopQuark,DK_RE_T,MomIn(1:4,1:4),GluonHel=HelGlu)

   if( StopQuark%PartType.eq.Stop_ ) then
      call vMajoSpi(dcmplx(Mom(1:4,1)),m_Chi,ChiHel,Spi(1:4))
      BarSpi(1:4) = TopQuark%Pol(1:4)
      Spi(1:4) = (0d0,1d0)*( IChiStt(+1)*Chir(.true.,Spi(1:4)) + IChiStt(-1)*Chir(.false.,Spi(1:4)) )
   elseif( StopQuark%PartType.eq.AStop_ ) then
      call ubarMajoSpi(dcmplx(Mom(1:4,1)),m_Chi,ChiHel,BarSpi(1:4))
      Spi(1:4) = TopQuark%Pol(1:4)      
      Spi(1:4) = (0d0,1d0)*( IChiStt(+1)*Chir(.true.,Spi(1:4)) + IChiStt(-1)*Chir(.false.,Spi(1:4)) )
   endif
   StopQuark%Pol(1) = psp1_(BarSpi(1:4),Spi(1:4)) * NWAFactor_STop



ELSEIF( Topol.eq.DKX_STChi0_RE3 ) THEN! gluon emission off W+/-


   TopQuark%PartType= sign(1,StopQuark%PartType) * Top_
   TopQuark%Mass = m_Top
   TopQuark%Mass2= m_Top**2
   TopQuark%Mom(1:4)=dcmplx(Mom(1:4,2))
   MomIn(1:4,1:3) = Mom(1:4,3:5)
   MomIn(1:4,4)   = MomGlu(1:4)
   call TopDecay(TopQuark,DK_RE_Q,MomIn(1:4,1:4),GluonHel=HelGlu)

   if( StopQuark%PartType.eq.Stop_ ) then
      call vMajoSpi(dcmplx(Mom(1:4,1)),m_Chi,ChiHel,Spi(1:4))
      BarSpi(1:4) = TopQuark%Pol(1:4)
      Spi(1:4) = (0d0,1d0)*( IChiStt(+1)*Chir(.true.,Spi(1:4)) + IChiStt(-1)*Chir(.false.,Spi(1:4)) )
   elseif( StopQuark%PartType.eq.AStop_ ) then
      call ubarMajoSpi(dcmplx(Mom(1:4,1)),m_Chi,ChiHel,BarSpi(1:4))
      Spi(1:4) = TopQuark%Pol(1:4)      
      Spi(1:4) = (0d0,1d0)*( IChiStt(+1)*Chir(.true.,Spi(1:4)) + IChiStt(-1)*Chir(.false.,Spi(1:4)) )
   endif
   StopQuark%Pol(1) = psp1_(BarSpi(1:4),Spi(1:4)) * NWAFactor_STop






ELSEIF( Topol.eq.DKX_STChi0_1L1 ) THEN! virtual correction and integr.dipole on (anti)stop

   TopQuark%PartType= sign(1,StopQuark%PartType) * Top_
   TopQuark%Mass = m_Top
   TopQuark%Mass2= m_Top**2
   TopQuark%Mom(1:4)=dcmplx(Mom(1:4,2))
   call TopDecay(TopQuark,DK_LO,Mom(1:4,3:5))
   if( present(HelTop) ) then! stable top
      TopQuark%Helicity = HelTop
      if( TopQuark%PartType.gt.0 ) then 
          call ubarSpi(TopQuark%Mom(1:4),TopQuark%Mass,TopQuark%Helicity,TopQuark%Pol(1:4))
      else
          call vSpi(TopQuark%Mom(1:4),TopQuark%Mass,TopQuark%Helicity,TopQuark%Pol(1:4))
      endif
   endif




   if( StopQuark%PartType.eq.Stop_ ) then
      call vMajoSpi(dcmplx(Mom(1:4,1)),m_Chi,ChiHel,Spi(1:4))
      BarSpi(1:4) = TopQuark%Pol(1:4)
      SpiPlus(1:4)  = (0d0,1d0)*( IChiStt(+1)*Chir(.true.,Spi(1:4)) + IChiStt(-1)*Chir(.false.,Spi(1:4)) )!  =Tree(+1)
      SpiMinus(1:4) = (0d0,1d0)*( IChiStt(-1)*Chir(.true.,Spi(1:4)) + IChiStt(+1)*Chir(.false.,Spi(1:4)) )!  =Tree(-1)


      Norm = CF *  alpha_sOver2Pi*RunAlphaS(NLOParam,MuRen)! = CF * gs^2/(16pi^2) * 2*Re[..]
      epMI = 0

      MI(1) = 1 * qlI2(M_Chi**2,M_Stop**2,M_top**2,MuRen**2,epMI) !MI(2,M_Chi**2,M_Stop**2,M_top**2)
      MI(2) = 1 * qlI2(M_Stop**2,0d0,M_Stop**2,MuRen**2,epMI) !MI(2,M_Stop**2,0,M_Stop**2)
      MI(3) = 1 * qlI2(M_top**2,0d0,M_top**2,MuRen**2,epMI) !MI(2,M_top**2,0,M_top**2)
      MI(4) = 1 * qlI3(M_Chi**2,M_top**2,M_Stop**2,M_Stop**2,M_top**2,0d0,MuRen**2,epMI) !MI(3,M_Chi**2,M_top**2,M_Stop**2,M_Stop**2,M_top**2,0)

      ! virtual spinor
      ff2= ((m_Chi*m_Top*(m_Chi**2 + m_Stop**2 - m_Top**2)*MI(1))/ &
             (m_Chi**4 + (m_Stop**2 - m_Top**2)**2 - 2*m_Chi**2*(m_Stop**2 + m_Top**2)) &
             - (2*m_Chi*m_Stop**2*m_Top*MI(2))/ &
             (m_Chi**4 + (m_Stop**2 - m_Top**2)**2 - 2*m_Chi**2*(m_Stop**2 + m_Top**2)) &
             + (m_Chi*m_Top*(-m_Chi**2 + m_Stop**2 + m_Top**2)*MI(3))/ &
             (m_Chi**4 + (m_Stop**2 - m_Top**2)**2 - 2*m_Chi**2*(m_Stop**2 + m_Top**2)))

      ff1= ((m_Chi**2*(-m_Chi**2 + m_Stop**2 + m_Top**2)*MI(1))/ &
             (m_Chi**4 + (m_Stop**2 - m_Top**2)**2 - 2*m_Chi**2*(m_Stop**2 + m_Top**2)) &
             + ((-m_Stop**4 + (m_Chi**2 - m_Top**2)**2)*MI(2))/ &
             (2.0d0*(m_Chi**4 + (m_Stop**2 - m_Top**2)**2 - &
                 2*m_Chi**2*(m_Stop**2 + m_Top**2))) + &
            (((m_Chi**2 - m_Stop**2)**2 - (m_Chi**2 + m_Stop**2)*m_Top**2)*MI(3))/ &
             (m_Chi**4 + (m_Stop**2 - m_Top**2)**2 - 2*m_Chi**2*(m_Stop**2 + m_Top**2)) &
             + (-m_Chi**2 + m_Stop**2 + m_Top**2)*MI(4))

      FF1 = FF1*dsqrt(2d0)**2
      FF2 = FF2*dsqrt(2d0)**2
      if(epMI.eq.-1) FF1 = FF1 - 3d0
      if(epMI.eq.0)  FF1 = FF1 - ( 4d0 +3d0*dlog(MuRen**2/m_Top**2) )
      Spi(1:4) = ( FF1*SpiPlus(1:4) + FF2*SpiMinus(1:4) ) * Norm



      ! integrated dipole
      beta = m_top/m_STop
      omega = m_Chi/m_Stop
      z = m_Top**2/m_Stop**2
      P0 = 0.5d0*( 1d0 - omega**2 + z )
      P3 = 0.5d0*SqrtLambda(1d0,omega**2,z)
      W0 = 0.5d0*( 1d0+omega**2-z )
      Ppl = P0 + P3
      Pmi = P0 - P3
      Wpl = W0 + P3
      Wmi = W0 - P3
      Yp = 0.5d0 * dlog(Ppl/Pmi)
      Yw = 0.5d0 * dlog(Wpl/Wmi)

      ID(-1) = 1d0*( 2d0*(1d0-P0/P3*Yp) )
      ID(0)  = 1d0*( -4d0*dlog(4d0*P3**2/omega/beta) * (1d0-P0/P3*Yp) + 4d0 + 2d0/P3*( (1d0-omega**2)*Yp + (1d0-beta**2)*Yw )   & 
                     +P0/P3*( 2d0*Yp-6d0*Yp**2 + 4d0*Yw*dlog(beta) - 6d0*DLi2(1d0-Pmi/Ppl) - 2d0*DLi2(1d0-Ppl) + 2d0*DLi2(1d0-Pmi) ) )
      Spi(1:4) = Spi(1:4) + ID(epMI)*SpiPlus(1:4)*Norm

! spi=spiplus! this has to give the LO result
! print *, "spi=spiplus"

! StopQuark%Pol(1) = 1d0
! RETURN

!print *, "check stop"
!print *, "integrals",MI(1)
!print *, "integrals",MI(2)
!print *, "integrals",MI(3)
!print *, "integrals",MI(4)
!print *, "virt",FF1,FF2
!print *, "intD",ID(epMI)
!print *, "sum ",FF1 + ID(epMI)
!print *, "ratio",FF1/ID(epMI)
! v13 = dsqrt( 1d0-4d0*m_Stop**2*m_Top**2/(m_Stop**2+m_Top**2-m_Chi**2)**2 )
! Catani = -2d0 - 1/v13 * dlog((1d0-v13)/(1d0+v13))
! print *, "Catani",Catani
!pause





   elseif( StopQuark%PartType.eq.AStop_ ) then
      call ubarMajoSpi(dcmplx(Mom(1:4,1)),m_Chi,ChiHel,BarSpi(1:4))
      Spi(1:4) = TopQuark%Pol(1:4)      
      SpiPlus(1:4)  = (0d0,1d0)*( IChiStt(+1)*Chir(.false.,Spi(1:4)) + IChiStt(-1)*Chir(.true.,Spi(1:4)) )!  =Tree(1)
      SpiMinus(1:4) = (0d0,1d0)*( IChiStt(-1)*Chir(.false.,Spi(1:4)) + IChiStt(+1)*Chir(.true.,Spi(1:4)) )!  =Tree(-1)

      Norm = CF   *  alpha_sOver2Pi*RunAlphaS(NLOParam,MuRen)! = CF * gs^2/(16pi^2) * 2*Re[..]
      epMI = 0

      MI(1) = 1  *  qlI2(M_Chi**2,M_Stop**2,M_top**2,MuRen**2,epMI) !SI(2,M_Chi**2,M_Stop**2,M_top**2)
      MI(2) = 1  *  qlI2(M_Stop**2,0d0,M_Stop**2,MuRen**2,epMI) !SI(2,M_Stop**2,0,M_Stop**2)
      MI(3) = 1  *  qlI2(M_top**2,0d0,M_top**2,MuRen**2,epMI) !SI(2,M_top**2,0,M_top**2)
      MI(4) = 1  *  qlI3(M_Chi**2,M_top**2,M_Stop**2,M_Stop**2,M_top**2,0d0,MuRen**2,epMI) !SI(3,M_Chi**2,M_top**2,M_Stop**2,M_Stop**2,M_top**2,0)
      

      ! virtual spinor
      ff2 = ((m_Chi*m_Top*(m_Chi**2 + m_Stop**2 - m_Top**2)*MI(1))/ &
             (m_Chi**4 + (m_Stop**2 - m_Top**2)**2 - 2*m_Chi**2*(m_Stop**2 + m_Top**2)) &
             - (2*m_Chi*m_Stop**2*m_Top*MI(2))/ &
             (m_Chi**4 + (m_Stop**2 - m_Top**2)**2 - 2*m_Chi**2*(m_Stop**2 + m_Top**2)) &
             + (m_Chi*m_Top*(-m_Chi**2 + m_Stop**2 + m_Top**2)*MI(3))/ &
             (m_Chi**4 + (m_Stop**2 - m_Top**2)**2 - 2*m_Chi**2*(m_Stop**2 + m_Top**2)))

      ff1 = ((m_Chi**2*(-m_Chi**2 + m_Stop**2 + m_Top**2)*MI(1))/ &
             (m_Chi**4 + (m_Stop**2 - m_Top**2)**2 - 2*m_Chi**2*(m_Stop**2 + m_Top**2)) &
             + ((-m_Stop**4 + (m_Chi**2 - m_Top**2)**2)*MI(2))/ &
             (2.0d0*(m_Chi**4 + (m_Stop**2 - m_Top**2)**2 - &
                 2*m_Chi**2*(m_Stop**2 + m_Top**2))) + &
            (((m_Chi**2 - m_Stop**2)**2 - (m_Chi**2 + m_Stop**2)*m_Top**2)*MI(3))/ &
             (m_Chi**4 + (m_Stop**2 - m_Top**2)**2 - 2*m_Chi**2*(m_Stop**2 + m_Top**2)) &
             + (-m_Chi**2 + m_Stop**2 + m_Top**2)*MI(4))

      FF1 = FF1*dsqrt(2d0)**2
      FF2 = FF2*dsqrt(2d0)**2
      if(epMI.eq.-1) FF1 = FF1 - 3d0
      if(epMI.eq.0)  FF1 = FF1 - ( 4d0 +3d0*dlog(MuRen**2/m_Top**2) )
      Spi(1:4) = ( FF1*SpiPlus(1:4) + FF2*SpiMinus(1:4) ) * Norm



      ! integrated dipole
      beta = m_top/m_STop
      omega = m_Chi/m_Stop
      z = m_Top**2/m_Stop**2
      P0 = 0.5d0*( 1d0 - omega**2 + z )
      P3 = 0.5d0*SqrtLambda(1d0,omega**2,z)
      W0 = 0.5d0*( 1d0+omega**2-z )
      Ppl = P0 + P3
      Pmi = P0 - P3
      Wpl = W0 + P3
      Wmi = W0 - P3
      Yp = 0.5d0 * dlog(Ppl/Pmi)
      Yw = 0.5d0 * dlog(Wpl/Wmi)

      ID(-1) = 1d0*( 2d0*(1d0-P0/P3*Yp) )
      ID(0)  = 1d0*( -4d0*dlog(4d0*P3**2/omega/beta) * (1d0-P0/P3*Yp) + 4d0 + 2d0/P3*( (1d0-omega**2)*Yp + (1d0-beta**2)*Yw )   & 
                     +P0/P3*( 2d0*Yp-6d0*Yp**2 + 4d0*Yw*dlog(beta) - 6d0*DLi2(1d0-Pmi/Ppl) - 2d0*DLi2(1d0-Ppl) + 2d0*DLi2(1d0-Pmi) ) )
      Spi(1:4) = Spi(1:4) + ID(epMI)*SpiPlus(1:4)*Norm

! spi=spiminus! this has to give the LO result
! print *, "spi=spiplus"



!print *, "check astop"
!print *, "integrals",MI(1)
!print *, "integrals",MI(2)
!print *, "integrals",MI(3)
!print *, "integrals",MI(4)
!print *, "virt",FF1,FF2
!print *, "intD",ID(epMI)
!print *, "sum ",FF1 + ID(epMI)
!print *, "ratio",FF1/ID(epMI)
! v13 = dsqrt( 1d0-4d0*m_Stop**2*m_Top**2/(m_Stop**2+m_Top**2-m_Chi**2)**2 )
! Catani = -2d0 - 1/v13 * dlog((1d0-v13)/(1d0+v13))
! print *, "Catani",Catani
!pause


   endif
   StopQuark%Pol(1) = psp1_(BarSpi(1:4),Spi(1:4)) * NWAFactor_STop






ELSEIF( Topol.eq.DKX_STChi0_1L2 ) THEN! virtual correction and integr.dipole on (anti)top


   TopQuark%PartType= sign(1,StopQuark%PartType) * Top_
   TopQuark%Mass = m_Top
   TopQuark%Mass2= m_Top**2
   TopQuark%Mom(1:4)=dcmplx(Mom(1:4,2))
   call TopDecay(TopQuark,DK_1L_T,Mom(1:4,3:5))
   if( present(HelTop) ) then! stable top
      TopQuark%Helicity = HelTop
      if( TopQuark%PartType.gt.0 ) then 
          call ubarSpi(TopQuark%Mom(1:4),TopQuark%Mass,TopQuark%Helicity,TopQuark%Pol(1:4))
      else
          call vSpi(TopQuark%Mom(1:4),TopQuark%Mass,TopQuark%Helicity,TopQuark%Pol(1:4))
      endif
   endif

   if( StopQuark%PartType.eq.Stop_ ) then
      call vMajoSpi(dcmplx(Mom(1:4,1)),m_Chi,ChiHel,Spi(1:4))
      BarSpi(1:4) = TopQuark%Pol(1:4)
      Spi(1:4) = (0d0,1d0)*( IChiStt(+1)*Chir(.true.,Spi(1:4)) + IChiStt(-1)*Chir(.false.,Spi(1:4)) )
   elseif( StopQuark%PartType.eq.AStop_ ) then
      call ubarMajoSpi(dcmplx(Mom(1:4,1)),m_Chi,ChiHel,BarSpi(1:4))
      Spi(1:4) = TopQuark%Pol(1:4)
      Spi(1:4) = (0d0,1d0)*( IChiStt(+1)*Chir(.true.,Spi(1:4)) + IChiStt(-1)*Chir(.false.,Spi(1:4)) )
   endif
   StopQuark%Pol(1) = psp1_(BarSpi(1:4),Spi(1:4)) * NWAFactor_STop



ELSEIF( Topol.eq.DKX_STChi0_1L3 ) THEN! virtual correction and integr.dipole on W+/W-


   TopQuark%PartType= sign(1,StopQuark%PartType) * Top_
   TopQuark%Mass = m_Top
   TopQuark%Mass2= m_Top**2
   TopQuark%Mom(1:4)=dcmplx(Mom(1:4,2))
   call TopDecay(TopQuark,DK_1L_Q,Mom(1:4,3:5))
   if( present(HelTop) ) then! stable top
      TopQuark%Helicity = HelTop
      if( TopQuark%PartType.gt.0 ) then 
          call ubarSpi(TopQuark%Mom(1:4),TopQuark%Mass,TopQuark%Helicity,TopQuark%Pol(1:4))
      else
          call vSpi(TopQuark%Mom(1:4),TopQuark%Mass,TopQuark%Helicity,TopQuark%Pol(1:4))
      endif
   endif

   if( StopQuark%PartType.eq.Stop_ ) then
      call vMajoSpi(dcmplx(Mom(1:4,1)),m_Chi,ChiHel,Spi(1:4))
      BarSpi(1:4) = TopQuark%Pol(1:4)
      Spi(1:4) = (0d0,1d0)*( IChiStt(+1)*Chir(.true.,Spi(1:4)) + IChiStt(-1)*Chir(.false.,Spi(1:4)) )
   elseif( StopQuark%PartType.eq.AStop_ ) then
      call ubarMajoSpi(dcmplx(Mom(1:4,1)),m_Chi,ChiHel,BarSpi(1:4))
      Spi(1:4) = TopQuark%Pol(1:4)
      Spi(1:4) = (0d0,1d0)*( IChiStt(+1)*Chir(.true.,Spi(1:4)) + IChiStt(-1)*Chir(.false.,Spi(1:4)) )
   endif
   StopQuark%Pol(1) = psp1_(BarSpi(1:4),Spi(1:4)) * NWAFactor_STop


ELSE
    call Error("Unknown topology in SUBROUTINE StopDecay")
ENDIF


END SUBROUTINE





END MODULE
