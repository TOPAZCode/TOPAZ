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



NWAFactor_HTop = 1d0/dsqrt(2d0*Ga_HTop*M_HTop)
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









SUBROUTINE HTopBHDecay(HeavyTop,Topol,BHHel,Mom)
use ModMisc
use ModProcess
use ModParameters
use ModTopDecay
implicit none
type(Particle) :: HeavyTop,TopQuark
integer :: Topol,NMom,BHHel,i
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



NWAFactor_HTop = 1d0/dsqrt(2d0*Ga_HTop*M_HTop)
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










SUBROUTINE StopDecay(StopQuark,Topol,ChiHel,Mom,MomGlu,HelGlu)
use ModMisc
use ModProcess
use ModParameters
use ModTopDecay
implicit none
type(Particle) :: StopQuark
type(Particle) :: TopQuark
integer :: Topol,ChiHel
real(8) :: Mom(:,:)! chi top bot lep neu
real(8),optional :: MomGlu(1:4)
integer,optional :: HelGlu
complex(8) :: Spi(1:4),BarSpi(1:4),Diagram1,Diagram2,PolGlu(1:4)
real(8) :: NWAFactor_STop

if( XTopDecays.eq.0 ) then
  StopQuark%Pol(1)    = (1d0,0d0)
  StopQuark%Pol(2:16) = (0d0,0d0)
  return
endif

   NWAFactor_STop = 1d0/dsqrt(2d0*Ga_STop*M_STop)


IF( Topol.eq.DKX_STChi0_LO ) THEN

   TopQuark%PartType= sign(1,StopQuark%PartType) * Top_
   TopQuark%Mass = m_Top
   TopQuark%Mass2= m_Top**2
   TopQuark%Mom(1:4)=dcmplx(Mom(1:4,2))
   call TopDecay(TopQuark,DK_LO,Mom(1:4,3:5))

   if( StopQuark%PartType.eq.Stop_ ) then
      call vMajoSpi(dcmplx(Mom(1:4,1)),m_Chi,ChiHel,Spi(1:4))
      BarSpi(1:4) = TopQuark%Pol(1:4)
      Spi(1:4) = (0d0,1d0)*( IChiStt(+1)*Chir(.true.,Spi(1:4)) + IChiStt(-1)*Chir(.false.,Spi(1:4)) )
   elseif( StopQuark%PartType.eq.AStop_ ) then
      call ubarMajoSpi(dcmplx(Mom(1:4,1)),m_Chi,ChiHel,BarSpi(1:4))
      Spi(1:4) = TopQuark%Pol(1:4)      
      Spi(1:4) = (0d0,1d0)*( IChiStt(+1)*Chir(.false.,Spi(1:4)) + IChiStt(-1)*Chir(.true.,Spi(1:4)) )
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







ELSEIF( Topol.eq.DKX_STChi0_RE1 ) THEN

   TopQuark%PartType= sign(1,StopQuark%PartType) * Top_
   TopQuark%Mass = m_Top
   TopQuark%Mass2= m_Top**2
   TopQuark%Mom(1:4)=dcmplx(Mom(1:4,2))
   call TopDecay(TopQuark,DK_LO,Mom(1:4,3:5))

   call pol_mless(dcmplx(MomGlu(1:4)),HelGlu,PolGlu(1:4))

   if( StopQuark%PartType.eq.Stop_ ) then
      call vMajoSpi(dcmplx(Mom(1:4,1)),m_Chi,ChiHel,Spi(1:4))
      BarSpi(1:4) = TopQuark%Pol(1:4)
      Spi(1:4) = (0d0,1d0)*( IChiStt(+1)*Chir(.true.,Spi(1:4)) + IChiStt(-1)*Chir(.false.,Spi(1:4)) )
      Diagram1 = psp1_(BarSpi(1:4),Spi(1:4)) * cgs(PolGlu(1:4),dcmplx(MomGlu(1:4)),StopQuark%Mom(1:4)-MomGlu(1:4))
      Diagram1 = Diagram1 * (0d0,1d0)/( ((StopQuark%Mom(1:4)-MomGlu(1:4)).dot.(StopQuark%Mom(1:4)-MomGlu(1:4))) -m_stop**2 )

      BarSpi(1:4) = vgq(PolGlu(1:4),TopQuark%Pol(1:4))
      BarSpi(1:4) = ( spb2_(TopQuark%Mom(1:4)+MomGlu(1:4),Spi(1:4)) + m_top*Spi(1:4) ) * (0d0,1d0)/( ((TopQuark%Mom(1:4)+MomGlu(1:4)).dot.(TopQuark%Mom(1:4)+MomGlu(1:4))) -m_top**2 )
      Diagram2 = psp1_(BarSpi(1:4),Spi(1:4))

   elseif( StopQuark%PartType.eq.AStop_ ) then
      call ubarMajoSpi(dcmplx(Mom(1:4,1)),m_Chi,ChiHel,BarSpi(1:4))
      Spi(1:4) = TopQuark%Pol(1:4)
      Spi(1:4) = (0d0,1d0)*( IChiStt(+1)*Chir(.false.,Spi(1:4)) + IChiStt(-1)*Chir(.true.,Spi(1:4)) )
      Diagram1 = psp1_(BarSpi(1:4),Spi(1:4)) * cgbs(PolGlu(1:4),dcmplx(MomGlu(1:4)),StopQuark%Mom(1:4)-MomGlu(1:4))
      Diagram1 = Diagram1 * (0d0,1d0)/( ((StopQuark%Mom(1:4)-MomGlu(1:4)).dot.(StopQuark%Mom(1:4)-MomGlu(1:4))) -m_stop**2 )

      Spi(1:4) = vgbq(PolGlu(1:4),TopQuark%Pol(1:4))
      Spi(1:4) = ( -spi2_(TopQuark%Mom(1:4)+MomGlu(1:4),Spi(1:4)) + m_top*Spi(1:4) ) * (0d0,1d0)/( ((TopQuark%Mom(1:4)+MomGlu(1:4)).dot.(TopQuark%Mom(1:4)+MomGlu(1:4))) -m_top**2 )
      Spi(1:4) = (0d0,1d0)*( IChiStt(+1)*Chir(.false.,Spi(1:4)) + IChiStt(-1)*Chir(.true.,Spi(1:4)) )
      Diagram2 = psp1_(BarSpi(1:4),Spi(1:4))
   endif

   StopQuark%Pol(1) = ( Diagram1 + Diagram2 ) * NWAFactor_STop


ELSEIF( Topol.eq.DKX_STChi0_RE2 ) THEN



ELSEIF( Topol.eq.DKX_STChi0_RE3 ) THEN



ELSE
    call Error("Unknown topology in SUBROUTINE StopDecay")
ENDIF


END SUBROUTINE





END MODULE
