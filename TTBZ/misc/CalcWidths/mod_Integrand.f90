MODULE ModIntegrand
implicit none




contains








FUNCTION Integrand_01(yRnd,VgsWgt)!   T --> W B
use ModParameters
use ModKinematics
use ModMisc
implicit none
real(8) :: Integrand_01
real(8) :: yRnd(1:VegasMxDim),VgsWgt
include 'vegas_common.f'
real(8) :: Mom(1:4,1:3),PSWgt,MG_Result,Norm

  Integrand_01 = 0d0
  Norm = 1d0/2d0/m_Top * yRnd(1)*2d0
  Mom(1:4,1) = (/m_Top,0d0,0d0,0d0/)
  call EvalPhasespace(2,m_Top,(/0d0,m_W/),yRnd(1:2),Mom(1:4,2:3),PSWgt)
  call ST_BWP(Mom(1:4,1:3),MG_Result)

  Integrand_01 = PSWgt * MG_Result * Norm

RETURN
END FUNCTION






FUNCTION Integrand_02(yRnd,VgsWgt)!   W --> E NU   or    W --> Q QB
use ModParameters
use ModKinematics
use ModMisc
implicit none
real(8) :: Integrand_02
real(8) :: yRnd(1:VegasMxDim),VgsWgt
include 'vegas_common.f'
real(8) :: Mom(1:4,1:3),PSWgt,MG_Result,Norm

  Integrand_02 = 0d0
  Norm = 1d0/2d0/m_W * yRnd(1)*2d0
  Mom(1:4,1) = (/m_W,0d0,0d0,0d0/)
  call EvalPhasespace(2,m_W,(/0d0,0d0/),yRnd(1:2),Mom(1:4,2:3),PSWgt)
!   call SWP_EPVE(Mom(1:4,1:3),MG_Result)
  call SWP_UDB(Mom(1:4,1:3),MG_Result)

  Integrand_02 = PSWgt * MG_Result * Norm

RETURN
END FUNCTION







FUNCTION Integrand_03(yRnd,VgsWgt)!   T --> W B G
use ModParameters
use ModKinematics
use ModMisc
implicit none
real(8) :: Integrand_03
real(8) :: yRnd(1:VegasMxDim),VgsWgt
include 'vegas_common.f'
real(8) :: Mom(1:4,1:4),PSWgt,MG_Result,Norm,alpha_sCorr
real(8) :: sbg

  Integrand_03 = 0d0

!   alpha_sCorr = 0.1258113d0/0.13d0!  1-loop running
  alpha_sCorr = 0.1095170d0/0.13d0!  2-loop running
  Norm = 1d0/2d0/m_Top * alpha_sCorr
  Mom(1:4,1) = (/m_Top,0d0,0d0,0d0/)
  call EvalPhasespace(3,m_Top,(/m_W,0d0,0d0/),yRnd(1:5),Mom(1:4,2:4),PSWgt)
  call swapMom(Mom(1:4,2),Mom(1:4,3))

  sbg = Mom(1:4,2).dot.Mom(1:4,4)
  if( sbg/m_Top**2.lt.1d-3 ) return

  call ST_BWPG(Mom(1:4,1:4),MG_Result)

  Integrand_03 = PSWgt * MG_Result * Norm

RETURN
END FUNCTION






FUNCTION Integrand_04(yRnd,VgsWgt)!   W --> Q QB G
use ModParameters
use ModKinematics
use ModMisc
implicit none
real(8) :: Integrand_04
real(8) :: yRnd(1:VegasMxDim),VgsWgt
include 'vegas_common.f'
real(8) :: Mom(1:4,1:4),PSWgt,MG_Result,Norm,alpha_sCorr
real(8) :: sqg

  Integrand_04 = 0d0


!   alpha_sCorr = 0.1258113d0/0.13d0!  1-loop running
  alpha_sCorr = 0.1095170d0/0.13d0!  2-loop running
  Norm = 1d0/2d0/m_W * alpha_sCorr
  Mom(1:4,1) = (/m_W,0d0,0d0,0d0/)
  call EvalPhasespace(3,m_W,(/0d0,0d0,0d0/),yRnd(1:5),Mom(1:4,2:4),PSWgt)

  sqg = dmin1(Mom(1:4,2).dot.Mom(1:4,4),Mom(1:4,3).dot.Mom(1:4,4))
  if( sqg/m_W**2.lt.1d-3 ) return

  call SWP_UDBG(Mom(1:4,1:3),MG_Result)

  Integrand_04 = PSWgt * MG_Result * Norm

RETURN
END FUNCTION








END MODULE
