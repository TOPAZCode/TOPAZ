MODULE ModDipolesDKJTTB
use ModTopdecay
implicit none
public :: DipolesDKJ
integer,private,parameter :: NumMaxHisto=45

contains



SUBROUTINE DipolesDKJ(Contr,nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
use ModProcess
implicit none
real(8) :: MomExt(1:4,1:12)
integer :: Contr,nDipole,NBin(1:NumMaxHisto),TopCorr
complex(8) :: Dipole(0:1,0:1)
logical :: applyPSCut
integer, parameter :: plus_ = 0, minus_= 1
integer, parameter :: NumDKPrimAmp = 2
type(Particle) :: ATop(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)! indices: Contribution, GluHel1, GluHel2, ColFacDK
type(Particle) ::  Top(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)


   Dipole(:,:) = (0d0,0d0)
   if(Contr.eq. 1) call DipolesDKJ_Contr1(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
   if(Contr.eq. 2) call DipolesDKJ_Contr2(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)! ANDREAS
   if(Contr.eq. 3) call DipolesDKJ_Contr3(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
   if(Contr.eq. 4) call DipolesDKJ_Contr4(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
   if(Contr.eq. 5) call DipolesDKJ_Contr5(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
   if(Contr.eq. 6) call DipolesDKJ_Contr6(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
   if(Contr.eq. 7) call DipolesDKJ_Contr7(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
   if(Contr.eq. 8) call DipolesDKJ_Contr8(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
   if(Contr.eq. 9) call DipolesDKJ_Contr9(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)! ANDREAS
   if(Contr.eq.10) call DipolesDKJ_Contr10(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
   if(Contr.eq.11) call DipolesDKJ_Contr11(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
   if(Contr.eq.12) call DipolesDKJ_Contr12(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
   if(Contr.eq.13) call DipolesDKJ_Contr13(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
   if(Contr.eq.14) call DipolesDKJ_Contr14(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
   if(Contr.eq.15) call DipolesDKJ_Contr15(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)! ANDREAS!!!! factor two wrong!!! + swapped u-d momenta
   if(Contr.eq.16) call DipolesDKJ_Contr16(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)! ANDREAS!!!  swapped u-d momenta
   if(Contr.eq.17) call DipolesDKJ_Contr17(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)! ANDREAS!!!  swapped u-d momenta
   if(Contr.eq.18) call DipolesDKJ_Contr18(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)! ANDREAS!!!! factor two wrong!!!
   if(Contr.eq.19) call DipolesDKJ_Contr19(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)! ANDREAS
   if(Contr.eq.20) call DipolesDKJ_Contr20(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)! ANDREAS


END SUBROUTINE









SUBROUTINE DipolesDKJ_Contr1(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
use ModProcess
use ModMisc
use ModParameters
use ModKinematics
implicit none
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,5:11)
real(8) :: sij,ColorCorr,Ria,fAux,PAux(1:4),PAux2(1:4)
real(8) :: z,y,RunFactor,ymax
logical :: applyPSCut
integer :: NBin(1:NumMaxHisto),TopCorr
integer :: i,j,k,nDipole
complex(8) :: Dipole(0:1,0:1),Splitting
complex(8) :: xp,xm,GluPol(1:4)
integer, parameter :: plus_ = 0, minus_= 1
real(8), parameter :: CA=3d0, CF=4d0/3d0
real(8), parameter :: Tt_x_Tb = 1d0/2d0*CA-CF
real(8), parameter :: Tt_x_Tg =-1d0/2d0*CA
real(8), parameter :: Tb_x_Tg =-1d0/2d0*CA
integer, parameter :: NumDKPrimAmp = 2
type(Particle) :: ATop(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)! indices: Contribution, GluHel1, GluHel2, ColFacDK
type(Particle) ::  Top(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)
#define InvertAlphaContr1 0


    applyPSCut = .false.
    TopCorr = ATop_
    call TopDecay(Top(1,1,1),DK_LO,MomExt(1:4,10:12))! evaluate this spinor once and for all dipoles before acceptance cuts
    RunFactor = RunAlphaS(NLOParam,MuRen)


IF( nDipole.eq.1 ) then
      i=5; j=8; k=9;
      call DKDip_FFmapping(MomExt(1:4,5:9),i-4,j-4,k-4,y,z,MomExtTd(1:4,5:8))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExtTd(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
!DEC$ IF( InvertAlphaContr1==0)
      if( y.gt.alpha_DKTff ) applyPSCut=.true.!       regular alpha cut
!DEC$ ELSE
      if( y.lt.alpha_DKTff ) applyPSCut=.true.!  inverted alpha cut for checks
!DEC$ ENDIF
      if( applyPSCut ) return
      call TopDecay(ATop(plus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**minus_))
      ColorCorr = Tb_x_Tg
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))
      Splitting = ( 2d0/(1d0-z*(1d0-y)) -1d0-z )
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0


ELSEIF( nDipole.eq.2 ) then
      i=5; j=9; k=8;
      call DKDip_FFmapping(MomExt(1:4,5:9),i-4,j-4,k-4,y,z,MomExtTd(1:4,5:8))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExtTd(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
!DEC$ IF( InvertAlphaContr1==0)
      if( y.gt.alpha_DKTff ) applyPSCut=.true.!       regular alpha cut
!DEC$ ELSE
      if( y.lt.alpha_DKTff ) applyPSCut=.true.!  inverted alpha cut for checks
!DEC$ ENDIF
      if( applyPSCut ) return
      call TopDecay(ATop(plus_,plus_,1), DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**minus_))
      ColorCorr = Tb_x_Tg
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))
      Splitting = ( 2d0/(1d0-z*(1d0-y)) -1d0-z )
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0


ELSEIF( nDipole.eq.3 ) then
      i=8; j=9; k=5;
      call DKDip_FFmapping(MomExt(1:4,5:9),i-4,j-4,k-4,y,z,MomExtTd(1:4,5:8))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExtTd(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
!DEC$ IF( InvertAlphaContr1==0)
      if( y.gt.alpha_DKTff ) applyPSCut=.true.!       regular alpha cut
!DEC$ ELSE
      if( y.lt.alpha_DKTff ) applyPSCut=.true.!  inverted alpha cut for checks
!DEC$ ENDIF
      if( applyPSCut ) return
      call TopDecay(ATop(plus_,plus_,1), DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**minus_))
      ColorCorr = Tb_x_Tg
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))

      PAux(1:4) = z*MomExt(1:4,i) - (1d0-z)*MomExt(1:4,j)
      call pol_mless(dcmplx(MomExtTd(1:4,8)),-1,GluPol(1:4))
      xm = dcmplx(PAux(1:4)).dot.GluPol(1:4)
      call pol_mless(dcmplx(MomExtTd(1:4,8)),+1,GluPol(1:4))
      xp = dcmplx(PAux(1:4)).dot.GluPol(1:4)

      Splitting = 2d0*(1d0/(1d0-z*(1d0-y)) + 1d0/(1d0-(1d0-z)*(1d0-y))-2d0 ) + 4d0/sij*dconjg(xp)*xp
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 2d0*(1d0/(1d0-z*(1d0-y)) + 1d0/(1d0-(1d0-z)*(1d0-y))-2d0 ) + 4d0/sij*dconjg(xm)*xm
      Dipole(minus_,minus_) = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 4d0/sij*dconjg(xp)*xm
      Dipole(plus_,minus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 4d0/sij*dconjg(xm)*xp
      Dipole(minus_,plus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij


ELSEIF( nDipole.eq.4 ) then
      i=5; j=8; k=3;
      call DKDip_FImapping(MomExt(1:4,5:9),i-4,j-4,y,z,Ria,ymax,MomExtTd(1:4,5:8))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExtTd(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
!DEC$ IF( InvertAlphaContr1==0)
      if( (1d0-alpha_DKTfi-z).ge.0d0 .and. (y-alpha_DKTfi*ymax).ge.0d0 ) applyPSCut=.true.!       regular alpha cut
!DEC$ ELSE
      if( (1d0-alpha_DKTfi-z).lt.0d0 .or. (y-alpha_DKTfi*ymax).lt.0d0 ) applyPSCut=.true.!       inverted alpha cut
!DEC$ ENDIF
      if( applyPSCut ) return
      call TopDecay(ATop(plus_,plus_,1), DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**minus_))
      ColorCorr = Tt_x_Tb
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))
      Splitting = (1d0+z**2)/(1d0-z) - (MomExt(1:4,i).dot.MomExt(1:4,j))*m_Top**2/(MomExt(1:4,k).dot.MomExt(1:4,j))**2
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0


ELSEIF( nDipole.eq.5 ) then
      i=5; j=9; k=3;
      call DKDip_FImapping(MomExt(1:4,5:9),i-4,j-4,y,z,Ria,ymax,MomExtTd(1:4,5:8))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExtTd(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
!DEC$ IF( InvertAlphaContr1==0)
      if( (1d0-alpha_DKTfi-z).ge.0d0 .and. (y-alpha_DKTfi*ymax).ge.0d0 ) applyPSCut=.true.!       regular alpha cut
!DEC$ ELSE
      if( (1d0-alpha_DKTfi-z).lt.0d0 .or. (y-alpha_DKTfi*ymax).lt.0d0 ) applyPSCut=.true.!       inverted alpha cut
!DEC$ ENDIF
      if( applyPSCut ) return
      call TopDecay(ATop(plus_,plus_,1), DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**minus_))
      ColorCorr = Tt_x_Tb
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))
      Splitting = (1d0+z**2)/(1d0-z) - (MomExt(1:4,i).dot.MomExt(1:4,j))*m_Top**2/(MomExt(1:4,k).dot.MomExt(1:4,j))**2
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0


ELSEIF( nDipole.eq.6 ) then
      i=8; j=9; k=3;
      call DKDip_FImapping(MomExt(1:4,5:9),i-4,j-4,y,z,Ria,ymax,MomExtTd(1:4,5:8))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExtTd(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
!DEC$ IF( InvertAlphaContr1==0)
      if( (1d0-alpha_DKTfi-z).ge.0d0 .and. (y-alpha_DKTfi*ymax).ge.0d0 ) applyPSCut=.true.!       regular alpha cut
!DEC$ ELSE
      if( (1d0-alpha_DKTfi-z).lt.0d0 .or. (y-alpha_DKTfi*ymax).lt.0d0 ) applyPSCut=.true.!       inverted alpha cut
!DEC$ ENDIF
      if( applyPSCut ) return
      call TopDecay(ATop(plus_,plus_,1), DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**minus_))
      ColorCorr = Tt_x_Tg
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))


      PAux(1:4) = MomExt(1:4,5) + MomExt(1:4,6) + MomExt(1:4,7)
      fAux = 4d0/Ria**2 * ((MomExt(1:4,k).dot.PAux(1:4))**2- m_Top**2*(m_Top**2-Ria) )
      PAux2(1:4)= ( (MomExt(1:4,k).dot.MomExt(1:4,i))*MomExt(1:4,i) - (MomExt(1:4,k).dot.MomExt(1:4,j))*MomExt(1:4,j) )*2d0/Ria
      PAux(1:4) = 1d0/(MomExt(1:4,k).dot.MomExtTd(1:4,8))*( (MomExt(1:4,k).dot.MomExtTd(1:4,8))*PAux2(1:4)  - (PAux2(1:4).dot.MomExtTd(1:4,8))*MomExt(1:4,k)  ) * fAux
      call pol_mless(dcmplx(MomExtTd(1:4,8)),-1,GluPol(1:4))
      xm = dcmplx(PAux(1:4)).dot.GluPol(1:4)
      call pol_mless(dcmplx(MomExtTd(1:4,8)),+1,GluPol(1:4))
      xp = dcmplx(PAux(1:4)).dot.GluPol(1:4)

      Splitting = 2d0*( z/(1d0-z)   - m_Top**2/2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))/(MomExt(1:4,k).dot.MomExt(1:4,j))**2  ) + 2d0/sij*dconjg(xp)*xp
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 2d0*( z/(1d0-z)   - m_Top**2/2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))/(MomExt(1:4,k).dot.MomExt(1:4,j))**2  ) + 2d0/sij*dconjg(xm)*xm
      Dipole(minus_,minus_) = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 2d0/sij*dconjg(xp)*xm
      Dipole(plus_,minus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 2d0/sij*dconjg(xm)*xp
      Dipole(minus_,plus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij


ELSEIF( nDipole.eq.7 ) then
      i=9; j=8; k=3;
      call DKDip_FImapping(MomExt(1:4,5:9),i-4,j-4,y,z,Ria,ymax,MomExtTd(1:4,5:8))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExtTd(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
!DEC$ IF( InvertAlphaContr1==0)
      if( (1d0-alpha_DKTfi-z).ge.0d0 .and. (y-alpha_DKTfi*ymax).ge.0d0 ) applyPSCut=.true.!       regular alpha cut
!DEC$ ELSE
      if( (1d0-alpha_DKTfi-z).lt.0d0 .or. (y-alpha_DKTfi*ymax).lt.0d0 ) applyPSCut=.true.!       inverted alpha cut
!DEC$ ENDIF
      if( applyPSCut ) return
      call TopDecay(ATop(plus_,plus_,1), DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**minus_))
      ColorCorr = Tt_x_Tg
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))


      PAux(1:4) = MomExt(1:4,5) + MomExt(1:4,6) + MomExt(1:4,7)
      fAux = 4d0/Ria**2 * ((MomExt(1:4,k).dot.PAux(1:4))**2- m_Top**2*(m_Top**2-Ria) )
      PAux2(1:4)= ( (MomExt(1:4,k).dot.MomExt(1:4,i))*MomExt(1:4,i) - (MomExt(1:4,k).dot.MomExt(1:4,j))*MomExt(1:4,j) )*2d0/Ria
      PAux(1:4) = 1d0/(MomExt(1:4,k).dot.MomExtTd(1:4,8))*( (MomExt(1:4,k).dot.MomExtTd(1:4,8))*PAux2(1:4)  - (PAux2(1:4).dot.MomExtTd(1:4,8))*MomExt(1:4,k)  ) * fAux
      call pol_mless(dcmplx(MomExtTd(1:4,8)),-1,GluPol(1:4))
      xm = dcmplx(PAux(1:4)).dot.GluPol(1:4)
      call pol_mless(dcmplx(MomExtTd(1:4,8)),+1,GluPol(1:4))
      xp = dcmplx(PAux(1:4)).dot.GluPol(1:4)

      Splitting = 2d0*( z/(1d0-z)   - m_Top**2/2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))/(MomExt(1:4,k).dot.MomExt(1:4,j))**2  ) + 2d0/sij*dconjg(xp)*xp
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 2d0*( z/(1d0-z)   - m_Top**2/2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))/(MomExt(1:4,k).dot.MomExt(1:4,j))**2  ) + 2d0/sij*dconjg(xm)*xm
      Dipole(minus_,minus_) = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 2d0/sij*dconjg(xp)*xm
      Dipole(plus_,minus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 2d0/sij*dconjg(xm)*xp
      Dipole(minus_,plus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

ENDIF! dipoles


END SUBROUTINE





SUBROUTINE DipolesDKJ_Contr2(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
use ModProcess
use ModMisc
use ModParameters
use ModKinematics
use ModTTBJ_NLODKW
implicit none
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,5:9)
real(8) :: sij,ColorCorr,Ria,fAux,PAux(1:4),PAux2(1:4)
real(8) :: z,y,RunFactor,ymax, muS
logical :: applyPSCut
integer :: NBin(1:NumMaxHisto),TopCorr
integer :: i,j,k,nDipole
complex(8) :: Dipole(0:1,0:1),Splitting, dip
complex(8) :: xp,xm,GluPol(1:4)
integer, parameter :: plus_ = 0, minus_= 1
real(8), parameter :: CA=3d0, CF=4d0/3d0
real(8), parameter :: Tt_x_Tb = 1d0/2d0*CA-CF
real(8), parameter :: Tt_x_Tg =-1d0/2d0*CA
real(8), parameter :: Tb_x_Tg =-1d0/2d0*CA
integer, parameter :: NumDKPrimAmp = 2
integer :: xe, LO_NLO
type(Particle) :: ATop(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)! indices: Contribution, GluHel1, GluHel2, ColFacDK
type(Particle) ::  Top(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)


    applyPSCut = .false.
    TopCorr = ATop_
    call TopDecay(Top(1,1,1),DK_LO,MomExt(1:4,10:12))! evaluate this spinor once and for all dipoles before acceptance cuts
    RunFactor = RunAlphaS(NLOParam,MuRen)
    Dipole(0:1,0:1) = czero
! Markus order for tb -> bb W- -> bb q qb g g
!    5-9 == bb q qb g1 g2
! my order
!    1-4 == qb q g1 g2
call Wto3Jet_init

IF( nDipole.le. 8 ) then

   call WmDipolesQG(ndipole,MomExt(1:4,6:9),dip,MomExtTd(1:4,6:8)) ! order: qb , q, g1, g2 -> qb, q,  g
   MomExtTd(1:4,5) = MomExt(1:4,5) ! Set bottom-mom for Tilde-kienmatics in TopDecay
   if(dip .eq. czero) return
   call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
   if( applyPSCut ) return
   call TopDecay(ATop(plus_,plus_,1), DKJ_LO_Q,MomExtTd(1:4,5:8),GluonHel=int((-1)**plus_))
   call TopDecay(ATop(minus_,plus_,1),DKJ_LO_Q,MomExtTd(1:4,5:8),GluonHel=int((-1)**minus_))
   Dipole(plus_,plus_)   = -dip
   Dipole(minus_,minus_) = -dip
   Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0
elseif(nDipole .gt. 8 .and. nDipole .le. 10) then
   call WmDipolesGG(ndipole,MomExt(1:4,6:9),dip,1,1,MomExtTd(1:4,6:8))! order: qb , q, g1, g2 -> qb, q,  g

   MomExtTd(1:4,5) = MomExt(1:4,5) ! Set bottom-mom for Tilde-kienmatics in TopDecay

   if(dip .eq. czero) return
   call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
   if( applyPSCut ) return
   call TopDecay(ATop(plus_,plus_,1), DKJ_LO_Q,MomExtTd(1:4,5:8),GluonHel=int((-1)**plus_))
   call TopDecay(ATop(minus_,plus_,1),DKJ_LO_Q,MomExtTd(1:4,5:8),GluonHel=int((-1)**minus_))
   Dipole(plus_,plus_)   = -dip

   call WmDipolesGG(ndipole,MomExt(1:4,6:9),dip,1,-1,MomExtTd(1:4,6:8))

   Dipole(plus_,minus_) = -dip
   call WmDipolesGG(ndipole,MomExt(1:4,6:9),dip,-1,1,MomExtTd(1:4,6:8))

   Dipole(minus_,plus_) = -dip
   call WmDipolesGG(ndipole,MomExt(1:4,6:9),dip,-1,-1,MomExtTd(1:4,6:8))

   Dipole(minus_,minus_) = -dip
else
   write(*,*) "ERROR IN DipolesDKJ_Contr2 (SCHARF's), nDipole = ", nDipole
endif
! If(MomExt(1,9)/M_W .lt. 1.d-3) then
!    write(*,*) "ndip", nDipole
!    write(*,*) "Dipole", Dipole
! endif

END SUBROUTINE






SUBROUTINE DipolesDKJ_Contr3(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
use ModProcess
use ModMisc
use ModParameters
use ModKinematics
use ModHadrWDecay
implicit none
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,5:12)
real(8) :: sij,ColorCorr,Ria,fAux,PAux(1:4),PAux2(1:4)
real(8) :: z,y,RunFactor,ymax
logical :: applyPSCut
integer :: NBin(1:NumMaxHisto),TopCorr
integer :: i,j,k,nDipole
complex(8) :: Dipole(0:1,0:1)
real(8) :: Splitting
complex(8) :: xp,xm,GluPol(1:4)
integer, parameter :: plus_ = 0, minus_= 1
real(8), parameter :: CA=3d0, CF=4d0/3d0
integer, parameter :: NumDKPrimAmp = 2
type(Particle) :: ATop(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)! indices: Contribution, GluHel1, GluHel2, ColFacDK
type(Particle) ::  Top(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)



    applyPSCut = .false.
    RunFactor = RunAlphaS(NLOParam,MuRen)

IF( nDipole.eq.1 ) then
      call TopDecay(Top(1,1,1),DK_LO,MomExt(1:4,10:12))

      TopCorr = ATop_
      i=5; j=9; k=3;
      call DKDip_FImapping(MomExt(1:4,5:9),i-4,j-4,y,z,Ria,ymax,MomExtTd(1:4,5:8))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExtTd(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
      if( (1d0-alpha_DKTfi-z).ge.0d0 .and. (y-alpha_DKTfi*ymax).ge.0d0 ) applyPSCut=.true.
!DEC$ IF(_FactCheck .EQ.1)
      if( dabs(MomExtTd(1:4,6).dot.MomExtTd(1:4,8))/m_W**2.lt.1d-3 .or. dabs(MomExtTd(1:4,7).dot.MomExtTd(1:4,8))/m_W**2.lt.1d-3) applyPSCut=.true.
!DEC$ ENDIF
      if( applyPSCut ) return
      call TopDecay(ATop(plus_,plus_,1), DKJ_LO_Q,MomExtTd(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,plus_,1),DKJ_LO_Q,MomExtTd(1:4,5:8),GluonHel=int((-1)**minus_))
      ColorCorr =-CF
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))
      Splitting = (1d0+z**2)/(1d0-z) - (MomExt(1:4,i).dot.MomExt(1:4,j))*m_Top**2/(MomExt(1:4,k).dot.MomExt(1:4,j))**2
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0
ELSEIF( nDipole.eq.2 ) then
      call TopDecay(Top(1,1,1),DK_LO,MomExt(1:4,10:12))

      TopCorr = ATop_
      call wdec_trans(1,MomExt(1:4,5:8),MomExtTd(1:4,5:7),alpha_DKWff,Splitting)
      MomExtTd(1:4,8) = MomExt(1:4,9)
      if( Splitting.eq.0d0 ) then
          applyPSCut=.true.
          return
      endif
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,9),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExtTd(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
!DEC$ IF(_FactCheck .EQ.1)
      if( dabs(MomExtTd(1:4,5).dot.MomExtTd(1:4,8))/m_Top**2 .lt. 1d-3 ) applyPSCut=.true.
!DEC$ ENDIF
      if( applyPSCut ) return
      call TopDecay(ATop(plus_,plus_,1), DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**minus_))

      ColorCorr = -CF
      Dipole(plus_,plus_) = 4d0*DblPi*(alpha_s*RunFactor)* ColorCorr *Splitting
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0
ELSEIF( nDipole.eq.3 ) then
      call TopDecay(Top(1,1,1),DK_LO,MomExt(1:4,10:12))

      TopCorr = ATop_
      call wdec_trans(2,MomExt(1:4,5:8),MomExtTd(1:4,5:7),alpha_DKWff,Splitting)
      MomExtTd(1:4,8) = MomExt(1:4,9)
      if( Splitting.eq.0d0 ) then
          applyPSCut=.true.
          return
      endif
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,9),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExtTd(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
!DEC$ IF(_FactCheck .EQ.1)
      if( dabs(MomExtTd(1:4,5).dot.MomExtTd(1:4,8))/m_Top**2 .lt. 1d-3 ) applyPSCut=.true.
!DEC$ ENDIF
      if( applyPSCut ) return
      call TopDecay(ATop(plus_,plus_,1), DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**minus_))

      ColorCorr = -CF
      Dipole(plus_,plus_) = 4d0*DblPi*(alpha_s*RunFactor)* ColorCorr *Splitting
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0
ENDIF! dipoles



END SUBROUTINE









SUBROUTINE DipolesDKJ_Contr4(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
use ModProcess
use ModMisc
use ModParameters
use ModKinematics
implicit none
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,5:12)
real(8) :: sij,ColorCorr,Ria,fAux,PAux(1:4),PAux2(1:4)
real(8) :: z,y,RunFactor,ymax
logical :: applyPSCut
integer :: NBin(1:NumMaxHisto),TopCorr
integer :: i,j,k,nDipole
complex(8) :: Dipole(0:1,0:1),Splitting
complex(8) :: xp,xm,GluPol(1:4)
integer, parameter :: plus_ = 0, minus_= 1
real(8), parameter :: CA=3d0, CF=4d0/3d0
integer, parameter :: NumDKPrimAmp = 2
type(Particle) :: ATop(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)! indices: Contribution, GluHel1, GluHel2, ColFacDK
type(Particle) ::  Top(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)

    applyPSCut = .false.
    RunFactor = RunAlphaS(NLOParam,MuRen)

IF( nDipole.eq.1 ) then
      call TopDecay(Top(plus_,plus_,1), DKJ_LO_T,MomExt(1:4,9:12),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,plus_,1),DKJ_LO_T,MomExt(1:4,9:12),GluonHel=int((-1)**minus_))

      TopCorr = Top_
      i=5; j=8; k=3;
      call DKDip_FImapping(MomExt(1:4,5:8),i-4,j-4,y,z,Ria,ymax,MomExtTd(1:4,5:7))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,12),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExtTd(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,9),MomExt(1:4,10),MomExt(1:4,11)/),applyPSCut,NBin)
      if( (1d0-alpha_DKTfi-z).ge.0d0 .and. (y-alpha_DKTfi*ymax).ge.0d0 ) applyPSCut=.true.
!DEC$ IF(_FactCheck .EQ.1)
      if( dabs(MomExt(1:4,9).dot.MomExt(1:4,12))/m_Top**2 .lt. 1d-3 ) applyPSCut=.true.
!DEC$ ENDIF
      if( applyPSCut ) return
      call TopDecay(ATop(1,1,1),DK_LO,MomExtTd(1:4,5:7))
      ColorCorr =-CF
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))
      Splitting = (1d0+z**2)/(1d0-z) - (MomExt(1:4,i).dot.MomExt(1:4,j))*m_Top**2/(MomExt(1:4,k).dot.MomExt(1:4,j))**2
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0


ELSEIF( nDipole.eq.2 ) then
      call TopDecay(ATop(plus_,plus_,1), DKJ_LO_T,MomExt(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,plus_,1),DKJ_LO_T,MomExt(1:4,5:8),GluonHel=int((-1)**minus_))

      TopCorr = ATop_
      i=9; j=12; k=4;
      call DKDip_FImapping(MomExt(1:4,9:12),i-8,j-8,y,z,Ria,ymax,MomExtTd(1:4,9:11))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExtTd(1:4,9),MomExtTd(1:4,10),MomExtTd(1:4,11)/),applyPSCut,NBin)
      if( (1d0-alpha_DKTfi-z).ge.0d0 .and. (y-alpha_DKTfi*ymax).ge.0d0 ) applyPSCut=.true.
!DEC$ IF(_FactCheck .EQ.1)
      if( dabs(MomExt(1:4,5).dot.MomExt(1:4,8))/m_Top**2 .lt. 1d-3 ) applyPSCut=.true.
!DEC$ ENDIF
      if( applyPSCut ) return
      call TopDecay(Top(1,1,1),DK_LO,MomExtTd(1:4,9:11))
      ColorCorr =-CF
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))
      Splitting = (1d0+z**2)/(1d0-z) - (MomExt(1:4,i).dot.MomExt(1:4,j))*m_Top**2/(MomExt(1:4,k).dot.MomExt(1:4,j))**2
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0

ENDIF! dipoles


END SUBROUTINE









SUBROUTINE DipolesDKJ_Contr5(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
use ModProcess
use ModMisc
use ModParameters
use ModKinematics
use ModHadrWDecay
implicit none
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,5:12)
real(8) :: sij,ColorCorr,Ria,fAux,PAux(1:4),PAux2(1:4)
real(8) :: z,y,RunFactor,ymax
logical :: applyPSCut
integer :: NBin(1:NumMaxHisto),TopCorr
integer :: i,j,k,nDipole
complex(8) :: Dipole(0:1,0:1)
real(8) :: Splitting
complex(8) :: xp,xm,GluPol(1:4)
integer, parameter :: plus_ = 0, minus_= 1
real(8), parameter :: CA=3d0, CF=4d0/3d0
integer, parameter :: NumDKPrimAmp = 2
type(Particle) :: ATop(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)! indices: Contribution, GluHel1, GluHel2, ColFacDK
type(Particle) ::  Top(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)



    applyPSCut = .false.
    RunFactor = RunAlphaS(NLOParam,MuRen)

IF( nDipole.eq.1 ) then
      call TopDecay(Top(plus_,0,1),  DKJ_LO_Q,MomExt(1:4,9:12),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,0,1), DKJ_LO_Q,MomExt(1:4,9:12),GluonHel=int((-1)**minus_))

      TopCorr = Top_
      i=5; j=8; k=3;
      call DKDip_FImapping(MomExt(1:4,5:8),i-4,j-4,y,z,Ria,ymax,MomExtTd(1:4,5:7))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,12),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExtTd(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,9),MomExt(1:4,10),MomExt(1:4,11)/),applyPSCut,NBin)
      if( (1d0-alpha_DKTfi-z).ge.0d0 .and. (y-alpha_DKTfi*ymax).ge.0d0 ) applyPSCut=.true.
      if( applyPSCut ) return
      call TopDecay(ATop(1,1,1),DK_LO,MomExtTd(1:4,5:7))
      ColorCorr =-CF
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))
      Splitting = (1d0+z**2)/(1d0-z) - (MomExt(1:4,i).dot.MomExt(1:4,j))*m_Top**2/(MomExt(1:4,k).dot.MomExt(1:4,j))**2
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0


ELSEIF( nDipole.eq.2 ) then
      call TopDecay(ATop(plus_,0,1),  DKJ_LO_T,MomExt(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,0,1), DKJ_LO_T,MomExt(1:4,5:8),GluonHel=int((-1)**minus_))

      TopCorr = ATop_
      call wdec_trans(1,MomExt(1:4,9:12),MomExtTd(1:4,9:11),alpha_DKWff,Splitting)
      if( Splitting.eq.0d0 ) then
          applyPSCut=.true.
          return
      endif
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExtTd(1:4,9),MomExtTd(1:4,10),MomExtTd(1:4,11)/),applyPSCut,NBin)
      if( applyPSCut ) return
      call TopDecay(Top(1,1,1),DK_LO,MomExtTd(1:4,9:11))

      ColorCorr = -CF
      Dipole(plus_,plus_) = 4d0*DblPi*(alpha_s*RunFactor)* ColorCorr *Splitting
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0


ELSEIF( nDipole.eq.3 ) then
      call TopDecay(ATop(plus_,0,1),  DKJ_LO_T,MomExt(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,0,1), DKJ_LO_T,MomExt(1:4,5:8),GluonHel=int((-1)**minus_))

      TopCorr = ATop_
      call wdec_trans(2,MomExt(1:4,9:12),MomExtTd(1:4,9:11),alpha_DKWff,Splitting)
      if( Splitting.eq.0d0 ) then
          applyPSCut=.true.
          return
      endif
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExtTd(1:4,9),MomExtTd(1:4,10),MomExtTd(1:4,11)/),applyPSCut,NBin)
      if( applyPSCut ) return
      call TopDecay(Top(1,1,1),DK_LO,MomExtTd(1:4,9:11))

      ColorCorr = -CF
      Dipole(plus_,plus_) = 4d0*DblPi*(alpha_s*RunFactor)* ColorCorr *Splitting
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0
ENDIF! dipoles



END SUBROUTINE






SUBROUTINE DipolesDKJ_Contr6(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
use ModProcess
use ModMisc
use ModParameters
use ModKinematics
use ModHadrWDecay
implicit none
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,5:12)
real(8) :: sij,ColorCorr,Ria,fAux,PAux(1:4),PAux2(1:4)
real(8) :: z,y,RunFactor,ymax
logical :: applyPSCut
integer :: NBin(1:NumMaxHisto),TopCorr
integer :: i,j,k,nDipole
complex(8) :: Dipole(0:1,0:1)
real(8) :: Splitting
complex(8) :: xp,xm,GluPol(1:4)
integer, parameter :: plus_ = 0, minus_= 1
real(8), parameter :: CA=3d0, CF=4d0/3d0
integer, parameter :: NumDKPrimAmp = 2
type(Particle) :: ATop(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)! indices: Contribution, GluHel1, GluHel2, ColFacDK
type(Particle) ::  Top(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)



    applyPSCut = .false.
    RunFactor = RunAlphaS(NLOParam,MuRen)

IF( nDipole.eq.1 ) then
      call TopDecay(ATop(plus_,0,1),  DKJ_LO_Q,MomExt(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,0,1), DKJ_LO_Q,MomExt(1:4,5:8),GluonHel=int((-1)**minus_))

      TopCorr = ATop_
      i=9; j=12; k=4;
      call DKDip_FImapping(MomExt(1:4,9:12),i-8,j-8,y,z,Ria,ymax,MomExtTd(1:4,9:11))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExtTd(1:4,9),MomExtTd(1:4,10),MomExtTd(1:4,11)/),applyPSCut,NBin)
      if( (1d0-alpha_DKTfi-z).ge.0d0 .and. (y-alpha_DKTfi*ymax).ge.0d0 ) applyPSCut=.true.
!DEC$ IF(_FactCheck .EQ.1)
      if( dabs(MomExt(1:4,6).dot.MomExt(1:4,8))/m_W**2.lt.1d-3 .or. dabs(MomExt(1:4,7).dot.MomExt(1:4,8))/m_W**2.lt. 1d-3 ) applyPSCut=.true.
!DEC$ ENDIF
      if( applyPSCut ) return
      call TopDecay(Top(1,1,1),DK_LO,MomExtTd(1:4,9:11))
      ColorCorr =-CF
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))
      Splitting = (1d0+z**2)/(1d0-z) - (MomExt(1:4,i).dot.MomExt(1:4,j))*m_Top**2/(MomExt(1:4,k).dot.MomExt(1:4,j))**2
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0


ELSEIF( nDipole.eq.2 ) then
      call TopDecay(Top(plus_,0,1),  DKJ_LO_T,MomExt(1:4,9:12),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,0,1), DKJ_LO_T,MomExt(1:4,9:12),GluonHel=int((-1)**minus_))

      TopCorr = Top_
      call wdec_trans(1,MomExt(1:4,5:8),MomExtTd(1:4,5:7),alpha_DKWff,Splitting)
      if( Splitting.eq.0d0 ) then
          applyPSCut=.true.
          return
      endif
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,12),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExtTd(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,9),MomExt(1:4,10),MomExt(1:4,11)/),applyPSCut,NBin)
!DEC$ IF(_FactCheck .EQ.1)
      if( dabs(MomExt(1:4,9).dot.MomExt(1:4,12))/m_Top**2 .lt. 1d-3 ) applyPSCut=.true.
!DEC$ ENDIF
      if( applyPSCut ) return
      call TopDecay(ATop(1,1,1),DK_LO,MomExtTd(1:4,5:7))

      ColorCorr = -CF
      Dipole(plus_,plus_) = 4d0*DblPi*(alpha_s*RunFactor)* ColorCorr *Splitting
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0


ELSEIF( nDipole.eq.3 ) then
      call TopDecay(Top(plus_,0,1),  DKJ_LO_T,MomExt(1:4,9:12),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,0,1), DKJ_LO_T,MomExt(1:4,9:12),GluonHel=int((-1)**minus_))

      TopCorr = Top_
      call wdec_trans(2,MomExt(1:4,5:8),MomExtTd(1:4,5:7),alpha_DKWff,Splitting)
      if( Splitting.eq.0d0 ) then
          applyPSCut=.true.
          return
      endif
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,12),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExtTd(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,9),MomExt(1:4,10),MomExt(1:4,11)/),applyPSCut,NBin)
!DEC$ IF(_FactCheck .EQ.1)
      if( dabs(MomExt(1:4,9).dot.MomExt(1:4,12))/m_Top**2 .lt. 1d-3 ) applyPSCut=.true.
!DEC$ ENDIF
      if( applyPSCut ) return
      call TopDecay(ATop(1,1,1),DK_LO,MomExtTd(1:4,5:7))

      ColorCorr = -CF
      Dipole(plus_,plus_) = 4d0*DblPi*(alpha_s*RunFactor)* ColorCorr *Splitting
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0
ENDIF! dipoles




END SUBROUTINE






SUBROUTINE DipolesDKJ_Contr7(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
use ModProcess
use ModMisc
use ModParameters
use ModKinematics
use ModHadrWDecay
implicit none
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,5:12)
real(8) :: sij,ColorCorr,Ria,fAux,PAux(1:4),PAux2(1:4)
real(8) :: z,y,RunFactor,ymax
logical :: applyPSCut
integer :: NBin(1:NumMaxHisto),TopCorr
integer :: i,j,k,nDipole
complex(8) :: Dipole(0:1,0:1)
real(8) :: Splitting
complex(8) :: xp,xm,GluPol(1:4)
integer, parameter :: plus_ = 0, minus_= 1
real(8), parameter :: CA=3d0, CF=4d0/3d0
integer, parameter :: NumDKPrimAmp = 2
type(Particle) :: ATop(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)! indices: Contribution, GluHel1, GluHel2, ColFacDK
type(Particle) ::  Top(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)



    applyPSCut = .false.
    RunFactor = RunAlphaS(NLOParam,MuRen)

IF( nDipole.eq.1 ) then
      call TopDecay(Top(plus_,0,1),  DKJ_LO_Q,MomExt(1:4,9:12),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,0,1), DKJ_LO_Q,MomExt(1:4,9:12),GluonHel=int((-1)**minus_))

      TopCorr = Top_
      call wdec_trans(1,MomExt(1:4,5:8),MomExtTd(1:4,5:7),alpha_DKWff,Splitting)
      if( Splitting.eq.0d0 ) then
          applyPSCut=.true.
          return
      endif
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,12),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExtTd(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,9),MomExt(1:4,10),MomExt(1:4,11)/),applyPSCut,NBin)
      if( applyPSCut ) return
      call TopDecay(ATop(1,1,1),DK_LO,MomExtTd(1:4,5:7))

      ColorCorr = -CF
      Dipole(plus_,plus_) = 4d0*DblPi*(alpha_s*RunFactor)* ColorCorr *Splitting
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0

ELSEIF( nDipole.eq.2 ) then
      call TopDecay(Top(plus_,0,1),  DKJ_LO_Q,MomExt(1:4,9:12),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,0,1), DKJ_LO_Q,MomExt(1:4,9:12),GluonHel=int((-1)**minus_))

      TopCorr = Top_
      call wdec_trans(2,MomExt(1:4,5:8),MomExtTd(1:4,5:7),alpha_DKWff,Splitting)
      if( Splitting.eq.0d0 ) then
          applyPSCut=.true.
          return
      endif
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,12),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExtTd(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,9),MomExt(1:4,10),MomExt(1:4,11)/),applyPSCut,NBin)
      if( applyPSCut ) return
      call TopDecay(ATop(1,1,1),DK_LO,MomExtTd(1:4,5:7))

      ColorCorr = -CF
      Dipole(plus_,plus_) = 4d0*DblPi*(alpha_s*RunFactor)* ColorCorr *Splitting
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0

ELSEIF( nDipole.eq.3 ) then
      call TopDecay(ATop(plus_,0,1),  DKJ_LO_Q,MomExt(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,0,1), DKJ_LO_Q,MomExt(1:4,5:8),GluonHel=int((-1)**minus_))

      TopCorr = ATop_
      call wdec_trans(1,MomExt(1:4,9:12),MomExtTd(1:4,9:11),alpha_DKWff,Splitting)
      if( Splitting.eq.0d0 ) then
          applyPSCut=.true.
          return
      endif
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExtTd(1:4,9),MomExtTd(1:4,10),MomExtTd(1:4,11)/),applyPSCut,NBin)
      if( applyPSCut ) return
      call TopDecay(Top(1,1,1),DK_LO,MomExtTd(1:4,9:11))

      ColorCorr = -CF
      Dipole(plus_,plus_) = 4d0*DblPi*(alpha_s*RunFactor)* ColorCorr *Splitting
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0

ELSEIF( nDipole.eq.4 ) then
      call TopDecay(ATop(plus_,0,1),  DKJ_LO_Q,MomExt(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,0,1), DKJ_LO_Q,MomExt(1:4,5:8),GluonHel=int((-1)**minus_))

      TopCorr = ATop_
      call wdec_trans(2,MomExt(1:4,9:12),MomExtTd(1:4,9:11),alpha_DKWff,Splitting)
      if( Splitting.eq.0d0 ) then
          applyPSCut=.true.
          return
      endif
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExtTd(1:4,9),MomExtTd(1:4,10),MomExtTd(1:4,11)/),applyPSCut,NBin)
      if( applyPSCut ) return
      call TopDecay(Top(1,1,1),DK_LO,MomExtTd(1:4,9:11))

      ColorCorr = -CF
      Dipole(plus_,plus_) = 4d0*DblPi*(alpha_s*RunFactor)* ColorCorr *Splitting
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0


ENDIF! dipoles



END SUBROUTINE











SUBROUTINE DipolesDKJ_Contr8(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
use ModProcess
use ModMisc
use ModParameters
use modKinematics
implicit none
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,5:11)
real(8) :: sij,ColorCorr,Ria,fAux,PAux(1:4),PAux2(1:4)
real(8) :: z,y,RunFactor,ymax
integer :: i,j,k,nDipole
logical :: applyPSCut
integer :: NBin(1:NumMaxHisto),TopCorr
complex(8) :: Dipole(0:1,0:1),Splitting
complex(8) :: xp,xm,GluPol(1:4)
integer, parameter :: plus_ = 0, minus_= 1
real(8), parameter :: CA=3d0, CF=4d0/3d0
real(8), parameter :: Tt_x_Tb = 1d0/2d0*CA-CF
real(8), parameter :: Tt_x_Tg =-1d0/2d0*CA
real(8), parameter :: Tb_x_Tg =-1d0/2d0*CA
integer, parameter :: NumDKPrimAmp = 2
type(Particle) :: ATop(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)! indices: Contribution, GluHel1, GluHel2, ColFacDK
type(Particle) ::  Top(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)


    applyPSCut = .false.
    TopCorr = Top_
    call TopDecay(ATop(1,1,1),DK_LO,MomExt(1:4,5:7))! evaluate this spinor once and for all dipoles before acceptance cuts
    RunFactor = RunAlphaS(NLOParam,MuRen)


! if( ndipole.ne.1  .and. ndipole.ne.3 .and. ndipole.ne.4 .and. ndipole.ne.7 )  return


IF( nDipole.eq.1 ) then
      i=8; j=11; k=12;
      call DKDip_FFmapping(MomExt(1:4,8:12),i-7,j-7,k-7,y,z,MomExtTd(1:4,8:11))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExtTd(1:4,8),MomExtTd(1:4,9),MomExtTd(1:4,10)/),applyPSCut,NBin)
      if( y.gt.alpha_DKTff ) applyPSCut=.true.
      if( applyPSCut ) return
      call TopDecay(Top(plus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**minus_))
      ColorCorr = Tb_x_Tg
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))
      Splitting = ( 2d0/(1d0-z*(1d0-y)) -1d0-z )
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0


ELSEIF( nDipole.eq.2 ) then
      i=8; j=12; k=11;
      call DKDip_FFmapping(MomExt(1:4,8:12),i-7,j-7,k-7,y,z,MomExtTd(1:4,8:11))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExtTd(1:4,8),MomExtTd(1:4,9),MomExtTd(1:4,10)/),applyPSCut,NBin)
      if( y.gt.alpha_DKTff ) applyPSCut=.true.
      if( applyPSCut ) return
      call TopDecay(Top(plus_,plus_,1), DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**minus_))
      ColorCorr = Tb_x_Tg
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))
      Splitting = ( 2d0/(1d0-z*(1d0-y)) -1d0-z )
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0


ELSEIF( nDipole.eq.3 ) then
      i=11; j=12; k=8;
      call DKDip_FFmapping(MomExt(1:4,8:12),i-7,j-7,k-7,y,z,MomExtTd(1:4,8:11))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExtTd(1:4,8),MomExtTd(1:4,9),MomExtTd(1:4,10)/),applyPSCut,NBin)
      if( y.gt.alpha_DKTff ) applyPSCut=.true.
      if( applyPSCut ) return
      call TopDecay(Top(plus_,plus_,1), DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**minus_))
      ColorCorr = Tb_x_Tg
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))

      PAux(1:4) = z*MomExt(1:4,i) - (1d0-z)*MomExt(1:4,j)
      call pol_mless(dcmplx(MomExtTd(1:4,11)),-1,GluPol(1:4))
      xm = dcmplx(PAux(1:4)).dot.GluPol(1:4)
      call pol_mless(dcmplx(MomExtTd(1:4,11)),+1,GluPol(1:4))
      xp = dcmplx(PAux(1:4)).dot.GluPol(1:4)

      Splitting = 2d0*(1d0/(1d0-z*(1d0-y)) + 1d0/(1d0-(1d0-z)*(1d0-y))-2d0 ) + 4d0/sij*dconjg(xp)*xp
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 2d0*(1d0/(1d0-z*(1d0-y)) + 1d0/(1d0-(1d0-z)*(1d0-y))-2d0 ) + 4d0/sij*dconjg(xm)*xm
      Dipole(minus_,minus_) = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 4d0/sij*dconjg(xp)*xm
      Dipole(plus_,minus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 4d0/sij*dconjg(xm)*xp
      Dipole(minus_,plus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij


ELSEIF( nDipole.eq.4 ) then
      i=8; j=11; k=4;
      call DKDip_FImapping(MomExt(1:4,8:12),i-7,j-7,y,z,Ria,ymax,MomExtTd(1:4,8:11))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExtTd(1:4,8),MomExtTd(1:4,9),MomExtTd(1:4,10)/),applyPSCut,NBin)
      if( (1d0-alpha_DKTfi-z).ge.0d0 .and. (y-alpha_DKTfi*ymax).ge.0d0 ) applyPSCut=.true.
      if( applyPSCut ) return
      call TopDecay(Top(plus_,plus_,1), DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**minus_))
      ColorCorr = Tt_x_Tb
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))
      Splitting = (1d0+z**2)/(1d0-z) - (MomExt(1:4,i).dot.MomExt(1:4,j))*m_Top**2/(MomExt(1:4,k).dot.MomExt(1:4,j))**2
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0


ELSEIF( nDipole.eq.5 ) then
      i=8; j=12; k=4;
      call DKDip_FImapping(MomExt(1:4,8:12),i-7,j-7,y,z,Ria,ymax,MomExtTd(1:4,8:11))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExtTd(1:4,8),MomExtTd(1:4,9),MomExtTd(1:4,10)/),applyPSCut,NBin)
      if( (1d0-alpha_DKTfi-z).ge.0d0 .and. (y-alpha_DKTfi*ymax).ge.0d0 ) applyPSCut=.true.
      if( applyPSCut ) return
      call TopDecay(Top(plus_,plus_,1), DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**minus_))
      ColorCorr = Tt_x_Tb
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))
      Splitting = (1d0+z**2)/(1d0-z) - (MomExt(1:4,i).dot.MomExt(1:4,j))*m_Top**2/(MomExt(1:4,k).dot.MomExt(1:4,j))**2
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0


ELSEIF( nDipole.eq.6 ) then
      i=11; j=12; k=4;
      call DKDip_FImapping(MomExt(1:4,8:12),i-7,j-7,y,z,Ria,ymax,MomExtTd(1:4,8:11))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExtTd(1:4,8),MomExtTd(1:4,9),MomExtTd(1:4,10)/),applyPSCut,NBin)
      if( (1d0-alpha_DKTfi-z).ge.0d0 .and. (y-alpha_DKTfi*ymax).ge.0d0 ) applyPSCut=.true.
      if( applyPSCut ) return
      call TopDecay(Top(plus_,plus_,1), DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**minus_))
      ColorCorr = Tt_x_Tg
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))
      PAux(1:4) = MomExt(1:4,8) + MomExt(1:4,9) + MomExt(1:4,10)
      fAux = 4d0/Ria**2 * ((MomExt(1:4,k).dot.PAux(1:4))**2- m_Top**2*(m_Top**2-Ria) )
      PAux2(1:4)= ( (MomExt(1:4,k).dot.MomExt(1:4,i))*MomExt(1:4,i) - (MomExt(1:4,k).dot.MomExt(1:4,j))*MomExt(1:4,j) )*2d0/Ria
      PAux(1:4) = 1d0/(MomExt(1:4,k).dot.MomExtTd(1:4,11))*( (MomExt(1:4,k).dot.MomExtTd(1:4,11))*PAux2(1:4)  - (PAux2(1:4).dot.MomExtTd(1:4,11))*MomExt(1:4,k)  ) * fAux
      call pol_mless(dcmplx(MomExtTd(1:4,11)),-1,GluPol(1:4))
      xm = dcmplx(PAux(1:4)).dot.GluPol(1:4)
      call pol_mless(dcmplx(MomExtTd(1:4,11)),+1,GluPol(1:4))
      xp = dcmplx(PAux(1:4)).dot.GluPol(1:4)

      Splitting = 2d0*( z/(1d0-z)   - m_Top**2/2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))/(MomExt(1:4,k).dot.MomExt(1:4,j))**2  ) + 2d0/sij*dconjg(xp)*xp
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 2d0*( z/(1d0-z)   - m_Top**2/2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))/(MomExt(1:4,k).dot.MomExt(1:4,j))**2  ) + 2d0/sij*dconjg(xm)*xm
      Dipole(minus_,minus_) = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 2d0/sij*dconjg(xp)*xm
      Dipole(plus_,minus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 2d0/sij*dconjg(xm)*xp
      Dipole(minus_,plus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij


ELSEIF( nDipole.eq.7 ) then
      i=12; j=11; k=4;
      call DKDip_FImapping(MomExt(1:4,8:12),i-7,j-7,y,z,Ria,ymax,MomExtTd(1:4,8:11))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExtTd(1:4,8),MomExtTd(1:4,9),MomExtTd(1:4,10)/),applyPSCut,NBin)
      if( (1d0-alpha_DKTfi-z).ge.0d0 .and. (y-alpha_DKTfi*ymax).ge.0d0 ) applyPSCut=.true.
      if( applyPSCut ) return
      call TopDecay(Top(plus_,plus_,1), DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**minus_))
      ColorCorr = Tt_x_Tg
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))


      PAux(1:4) = MomExt(1:4,8) + MomExt(1:4,9) + MomExt(1:4,10)
      fAux = 4d0/Ria**2 * ((MomExt(1:4,k).dot.PAux(1:4))**2- m_Top**2*(m_Top**2-Ria) )
      PAux2(1:4)= ( (MomExt(1:4,k).dot.MomExt(1:4,i))*MomExt(1:4,i) - (MomExt(1:4,k).dot.MomExt(1:4,j))*MomExt(1:4,j) )*2d0/Ria
      PAux(1:4) = 1d0/(MomExt(1:4,k).dot.MomExtTd(1:4,11))*( (MomExt(1:4,k).dot.MomExtTd(1:4,11))*PAux2(1:4)  - (PAux2(1:4).dot.MomExtTd(1:4,11))*MomExt(1:4,k)  ) * fAux
      call pol_mless(dcmplx(MomExtTd(1:4,11)),-1,GluPol(1:4))
      xm = dcmplx(PAux(1:4)).dot.GluPol(1:4)
      call pol_mless(dcmplx(MomExtTd(1:4,11)),+1,GluPol(1:4))
      xp = dcmplx(PAux(1:4)).dot.GluPol(1:4)

      Splitting = 2d0*( z/(1d0-z)   - m_Top**2/2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))/(MomExt(1:4,k).dot.MomExt(1:4,j))**2  ) + 2d0/sij*dconjg(xp)*xp
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 2d0*( z/(1d0-z)   - m_Top**2/2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))/(MomExt(1:4,k).dot.MomExt(1:4,j))**2  ) + 2d0/sij*dconjg(xm)*xm
      Dipole(minus_,minus_) = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 2d0/sij*dconjg(xp)*xm
      Dipole(plus_,minus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 2d0/sij*dconjg(xm)*xp
      Dipole(minus_,plus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

ENDIF! dipoles


END SUBROUTINE




SUBROUTINE DipolesDKJ_Contr9(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
use ModProcess
use ModMisc
use ModParameters
use ModKinematics
use ModTTBJ_NLODKW
implicit none
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,8:11)
real(8) :: sij,ColorCorr,Ria,fAux,PAux(1:4),PAux2(1:4), MomSR(1:4,1:4)
real(8) :: z,y,RunFactor,ymax, muS
logical :: applyPSCut
integer :: NBin(1:NumMaxHisto),TopCorr
integer :: i,j,k,nDipole
complex(8) :: Dipole(0:1,0:1),Splitting, dip
complex(8) :: xp,xm,GluPol(1:4)
integer, parameter :: plus_ = 0, minus_= 1
real(8), parameter :: CA=3d0, CF=4d0/3d0
real(8), parameter :: Tt_x_Tb = 1d0/2d0*CA-CF
real(8), parameter :: Tt_x_Tg =-1d0/2d0*CA
real(8), parameter :: Tb_x_Tg =-1d0/2d0*CA
integer, parameter :: NumDKPrimAmp = 2
integer :: xe, LO_NLO
type(Particle) :: ATop(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)! indices: Contribution, GluHel1, GluHel2, ColFacDK
type(Particle) ::  Top(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)

call Wto3Jet_init
    applyPSCut = .false.
    TopCorr = Top_
    call TopDecay(ATop(1,1,1),DK_LO,MomExt(1:4,5:7))! evaluate this spinor once and for all dipoles before acceptance cuts
    Dipole(0:1,0:1) = czero
! Markus order for t -> b W+ -> b qb q g g
!    9-12 == bb  qb q g1 g2
! my order
!    1-4 == qb q g1 g2


IF( nDipole.le. 8 ) then
   call WmDipolesQG(ndipole,MomExt(1:4,9:12),dip,MomExtTd(1:4,9:11)) ! order qb, q,g,g -> qb, q, g
   MomExtTd(1:4,8) = MomExt(1:4,8) ! Set b-momentum for TOP-decays with tilde kinematics
   if(dip .eq. czero) return
!   write(*,*) "t",(MomExt(1:4,8)+MomExtTd(1:4,9)+MomExtTd(1:4,10)+MomExtTd(1:4,11)).dot.(MomExt(1:4,8)+MomExtTd(1:4,9)+MomExtTd(1:4,10)+MomExtTd(1:4,11))/m_top**2
!   write(*,*) "W",(MomExtTd(1:4,9)+MomExtTd(1:4,10)+MomExtTd(1:4,11)).dot.(MomExtTd(1:4,9)+MomExtTd(1:4,10)+MomExtTd(1:4,11))/m_W**2
   call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExtTd(1:4,9),MomExtTd(1:4,10)/),applyPSCut,NBin)

   if( applyPSCut ) return

   call TopDecay(Top(plus_,plus_,1), DKJ_LO_Q,MomExtTd(1:4,8:11),GluonHel=int((-1)**plus_))
   call TopDecay(Top(minus_,plus_,1),DKJ_LO_Q,MomExtTd(1:4,8:11),GluonHel=int((-1)**minus_))


   Dipole(plus_,plus_)   = -dip
   Dipole(minus_,minus_) = -dip
   Dipole(plus_,minus_) = czero
   Dipole(minus_,plus_) = czero
elseif(nDipole .gt. 8 .and. nDipole .le. 10) then
   call WmDipolesGG(ndipole,MomExt(1:4,9:12),dip,1,1,MomExtTd(1:4,9:11))

   MomExtTd(1:4,8) = MomExt(1:4,8) ! Set b-momentum for TOP-decays with tilde kinematics
   if(dip .eq. czero) return

   call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExtTd(1:4,8),MomExtTd(1:4,9),MomExtTd(1:4,10)/),applyPSCut,NBin)

   if( applyPSCut ) return


   call TopDecay(Top(plus_,plus_,1),DKJ_LO_Q,MomExtTd(1:4,8:11),GluonHel=int((-1)**plus_))
   call TopDecay(Top(minus_,plus_,1),DKJ_LO_Q,MomExtTd(1:4,8:11),GluonHel=int((-1)**minus_))
   Dipole(plus_,plus_)   = -dip

   call WmDipolesGG(ndipole,MomExt(1:4,9:12),dip,1,-1,MomExtTd(1:4,9:11))
   Dipole(plus_,minus_) = -dip
   call WmDipolesGG(ndipole,MomExt(1:4,9:12),dip,-1,1,MomExtTd(1:4,9:11))
   Dipole(minus_,plus_) = -dip
   call WmDipolesGG(ndipole,MomExt(1:4,9:12),dip,-1,-1,MomExtTd(1:4,9:11))
   Dipole(minus_,minus_) = -dip
else
   write(*,*) "ERROR IN DipolesDKJ_Contr9 (SCHARF's), nDipole = ", nDipole
endif
!write(*,*) "ndipole", ndipole
!write(*,*) "++", Dipole(plus_,plus_)
!write(*,*) "+-", Dipole(plus_,minus_)
!write(*,*) "-+", Dipole(minus_,plus_)
! write(*,*) "--",Dipole(minus_,s_)
END SUBROUTINE





SUBROUTINE DipolesDKJ_Contr10(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
use ModProcess
use ModMisc
use ModParameters
use ModKinematics
use ModHadrWDecay
implicit none
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,5:12)
real(8) :: sij,ColorCorr,Ria,fAux,PAux(1:4),PAux2(1:4)
real(8) :: z,y,RunFactor,ymax
logical :: applyPSCut
integer :: NBin(1:NumMaxHisto),TopCorr
integer :: i,j,k,nDipole
complex(8) :: Dipole(0:1,0:1)
real(8) :: Splitting
complex(8) :: xp,xm,GluPol(1:4)
integer, parameter :: plus_ = 0, minus_= 1
real(8), parameter :: CA=3d0, CF=4d0/3d0
integer, parameter :: NumDKPrimAmp = 2
type(Particle) :: ATop(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)! indices: Contribution, GluHel1, GluHel2, ColFacDK
type(Particle) ::  Top(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)



    applyPSCut = .false.
    RunFactor = RunAlphaS(NLOParam,MuRen)

IF( nDipole.eq.1 ) then
      call TopDecay(ATop(1,1,1),DK_LO,MomExt(1:4,5:7))

      TopCorr = Top_
      i=8; j=12; k=4;
      call DKDip_FImapping(MomExt(1:4,8:12),i-7,j-7,y,z,Ria,ymax,MomExtTd(1:4,8:11))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExtTd(1:4,8),MomExtTd(1:4,9),MomExtTd(1:4,10)/),applyPSCut,NBin)
      if( (1d0-alpha_DKTfi-z).ge.0d0 .and. (y-alpha_DKTfi*ymax).ge.0d0 ) applyPSCut=.true.
      if( applyPSCut ) return
      call TopDecay(Top(plus_,plus_,1), DKJ_LO_Q,MomExtTd(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,plus_,1),DKJ_LO_Q,MomExtTd(1:4,8:11),GluonHel=int((-1)**minus_))
      ColorCorr =-CF
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))
      Splitting = (1d0+z**2)/(1d0-z) - (MomExt(1:4,i).dot.MomExt(1:4,j))*m_Top**2/(MomExt(1:4,k).dot.MomExt(1:4,j))**2
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0


ELSEIF( nDipole.eq.2 ) then
      call TopDecay(ATop(1,1,1),DK_LO,MomExt(1:4,5:7))

      TopCorr = Top_
      call wdec_trans(1,MomExt(1:4,8:11),MomExtTd(1:4,8:10),alpha_DKWff,Splitting)
      if( Splitting.eq.0d0 ) return
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,12),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExtTd(1:4,8),MomExtTd(1:4,9),MomExtTd(1:4,10)/),applyPSCut,NBin)
      if( applyPSCut ) return
      MomExtTd(1:4,11) = MomExt(1:4,12)
      call TopDecay(Top(plus_,plus_,1), DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**minus_))

      ColorCorr = -CF
      Dipole(plus_,plus_) = 4d0*DblPi*(alpha_s*RunFactor)* ColorCorr *Splitting
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0

ELSEIF( nDipole.eq.3 ) then
      call TopDecay(ATop(1,1,1),DK_LO,MomExt(1:4,5:7))

      TopCorr = Top_
      call wdec_trans(2,MomExt(1:4,8:11),MomExtTd(1:4,8:10),alpha_DKWff,Splitting)
      if( Splitting.eq.0d0 ) return
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,12),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExtTd(1:4,8),MomExtTd(1:4,9),MomExtTd(1:4,10)/),applyPSCut,NBin)
      if( applyPSCut ) return
      MomExtTd(1:4,11) = MomExt(1:4,12)
      call TopDecay(Top(plus_,plus_,1), DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**minus_))

      ColorCorr = -CF
      Dipole(plus_,plus_) = 4d0*DblPi*(alpha_s*RunFactor)* ColorCorr *Splitting
      Dipole(minus_,minus_) = Dipole(plus_,plus_)
      Dipole(plus_,minus_) = 0d0; Dipole(minus_,plus_) = 0d0
ENDIF! dipoles



END SUBROUTINE






SUBROUTINE DipolesDKJ_Contr11(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
use ModProcess
use ModMisc
use ModParameters
use ModKinematics
implicit none
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,5:11)
real(8) :: sij,ColorCorr,Ria,fAux,PAux(1:4),PAux2(1:4)
real(8) :: z,y,RunFactor,ymax
logical :: applyPSCut
integer :: NBin(1:NumMaxHisto),TopCorr
integer :: i,j,k,nDipole
complex(8) :: Dipole(0:1,0:1),Splitting
complex(8) :: xp,xm,GluPol(1:4)
integer, parameter :: plus_ = 0, minus_= 1
real(8), parameter :: CA=3d0, CF=4d0/3d0
integer, parameter :: NumDKPrimAmp = 2
type(Particle) :: ATop(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)! indices: Contribution, GluHel1, GluHel2, ColFacDK
type(Particle) ::  Top(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)
#define InvertAlphaContr11 0

    applyPSCut = .false.
    TopCorr = ATop_
    call TopDecay(Top(1,1,1),DK_LO,MomExt(1:4,10:12))! evaluate this spinor once and for all dipoles before acceptance cuts
    RunFactor = RunAlphaS(NLOParam,MuRen)


IF( nDipole.eq.1 ) then
      i=8; j=9; k=5;
      call DKDip_FFmapping(MomExt(1:4,5:9),i-4,j-4,k-4,y,z,MomExtTd(1:4,5:8))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExtTd(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
!DEC$ IF( InvertAlphaContr11==0)
      if( y.gt.alpha_DKTff ) applyPSCut=.true.!  regular alpha cut
!DEC$ ELSE
      if( y.lt.alpha_DKTff ) applyPSCut=.true.!  inverted alpha cut for checks
!DEC$ ENDIF
      if( applyPSCut ) return
      call TopDecay(ATop(plus_,plus_,1), DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**minus_))
      ColorCorr = -0.5d0
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))

      PAux(1:4) = z*MomExt(1:4,i) - (1d0-z)*MomExt(1:4,j)
      call pol_mless(dcmplx(MomExtTd(1:4,8)),-1,GluPol(1:4))
      xm = dcmplx(PAux(1:4)).dot.GluPol(1:4)
      call pol_mless(dcmplx(MomExtTd(1:4,8)),+1,GluPol(1:4))
      xp = dcmplx(PAux(1:4)).dot.GluPol(1:4)

      Splitting = 1d0 - 4d0/sij*dconjg(xp)*xp
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 1d0 - 4d0/sij*dconjg(xm)*xm
      Dipole(minus_,minus_) = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting =-4d0/sij*dconjg(xp)*xm
      Dipole(plus_,minus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting =-4d0/sij*dconjg(xm)*xp
      Dipole(minus_,plus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

ENDIF! dipoles


END SUBROUTINE



SUBROUTINE DipolesDKJ_Contr12(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
use ModProcess
use ModMisc
use ModParameters
use ModKinematics
implicit none
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,5:11)
real(8) :: sij,ColorCorr,Ria,fAux,PAux(1:4),PAux2(1:4)
real(8) :: z,y,RunFactor,ymax
logical :: applyPSCut
integer :: NBin(1:NumMaxHisto),TopCorr
integer :: i,j,k,nDipole
complex(8) :: Dipole(0:1,0:1),Splitting
complex(8) :: xp,xm,GluPol(1:4)
integer, parameter :: plus_ = 0, minus_= 1
real(8), parameter :: CA=3d0, CF=4d0/3d0
integer, parameter :: NumDKPrimAmp = 2
type(Particle) :: ATop(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)! indices: Contribution, GluHel1, GluHel2, ColFacDK
type(Particle) ::  Top(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)
#define InvertAlphaContr12 0

    applyPSCut = .false.
    TopCorr = ATop_
    call TopDecay(Top(1,1,1),DK_LO,MomExt(1:4,10:12))! evaluate this spinor once and for all dipoles before acceptance cuts
    RunFactor = RunAlphaS(NLOParam,MuRen)


IF( nDipole.eq.1 ) then
      i=8; j=9; k=5;
      call DKDip_FFmapping(MomExt(1:4,5:9),i-4,j-4,k-4,y,z,MomExtTd(1:4,5:8))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExtTd(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
!DEC$ IF( InvertAlphaContr12==0)
      if( y.gt.alpha_DKTff ) applyPSCut=.true.!  regular alpha cut
!DEC$ ELSE
      if( y.lt.alpha_DKTff ) applyPSCut=.true.!  inverted alpha cut for checks
!DEC$ ENDIF
      if( applyPSCut ) return
      call TopDecay(ATop(plus_,plus_,1), DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**minus_))
      ColorCorr = -0.5d0
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))

      PAux(1:4) = z*MomExt(1:4,i) - (1d0-z)*MomExt(1:4,j)
      call pol_mless(dcmplx(MomExtTd(1:4,8)),-1,GluPol(1:4))
      xm = dcmplx(PAux(1:4)).dot.GluPol(1:4)
      call pol_mless(dcmplx(MomExtTd(1:4,8)),+1,GluPol(1:4))
      xp = dcmplx(PAux(1:4)).dot.GluPol(1:4)

      Splitting = 1d0 - 4d0/sij*dconjg(xp)*xp
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 1d0 - 4d0/sij*dconjg(xm)*xm
      Dipole(minus_,minus_) = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting =-4d0/sij*dconjg(xp)*xm
      Dipole(plus_,minus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting =-4d0/sij*dconjg(xm)*xp
      Dipole(minus_,plus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij



ELSEIF( nDipole.eq.2 ) then
      i=5; j=9; k=8;
      call DKDip_FFmapping(MomExt(1:4,5:9),i-4,j-4,k-4,y,z,MomExtTd(1:4,5:8))
      call swapMom(MomExtTd(1:4,5),MomExtTd(1:4,8))! because (5,9) is recombined into 5 which is the "reduced: gluon that has to be in position 8
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExtTd(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
!DEC$ IF( InvertAlphaContr12==0)
      if( y.gt.alpha_DKTff ) applyPSCut=.true.!  regular alpha cut
!DEC$ ELSE
      if( y.lt.alpha_DKTff ) applyPSCut=.true.!  inverted alpha cut for checks
!DEC$ ENDIF
      if( applyPSCut ) return
      call TopDecay(ATop(plus_,plus_,1), DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,5:8),GluonHel=int((-1)**minus_))
      ColorCorr = -0.5d0
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))

      PAux(1:4) = z*MomExt(1:4,i) - (1d0-z)*MomExt(1:4,j)
      call pol_mless(dcmplx(MomExtTd(1:4,8)),-1,GluPol(1:4))
      xm = dcmplx(PAux(1:4)).dot.GluPol(1:4)
      call pol_mless(dcmplx(MomExtTd(1:4,8)),+1,GluPol(1:4))
      xp = dcmplx(PAux(1:4)).dot.GluPol(1:4)

      Splitting = 1d0 - 4d0/sij*dconjg(xp)*xp
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 1d0 - 4d0/sij*dconjg(xm)*xm
      Dipole(minus_,minus_) = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting =-4d0/sij*dconjg(xp)*xm
      Dipole(plus_,minus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting =-4d0/sij*dconjg(xm)*xp
      Dipole(minus_,plus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij
ENDIF! dipoles


END SUBROUTINE







SUBROUTINE DipolesDKJ_Contr13(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
use ModProcess
use ModMisc
use ModParameters
use modKinematics
implicit none
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,5:11)
real(8) :: sij,ColorCorr,Ria,fAux,PAux(1:4),PAux2(1:4)
real(8) :: z,y,RunFactor,ymax
integer :: i,j,k,nDipole
logical :: applyPSCut
integer :: NBin(1:NumMaxHisto),TopCorr
complex(8) :: Dipole(0:1,0:1),Splitting
complex(8) :: xp,xm,GluPol(1:4)
integer, parameter :: plus_ = 0, minus_= 1
real(8), parameter :: CA=3d0, CF=4d0/3d0
integer, parameter :: NumDKPrimAmp = 2
type(Particle) :: ATop(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)! indices: Contribution, GluHel1, GluHel2, ColFacDK
type(Particle) ::  Top(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)


    applyPSCut = .false.
    TopCorr = Top_
    call TopDecay(ATop(1,1,1),DK_LO,MomExt(1:4,5:7))! evaluate this spinor once and for all dipoles before acceptance cuts
    RunFactor = RunAlphaS(NLOParam,MuRen)



IF( nDipole.eq.1 ) then
      i=11; j=12; k=8;
      call DKDip_FFmapping(MomExt(1:4,8:12),i-7,j-7,k-7,y,z,MomExtTd(1:4,8:11))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExtTd(1:4,8),MomExtTd(1:4,9),MomExtTd(1:4,10)/),applyPSCut,NBin)
      if( y.gt.alpha_DKTff ) applyPSCut=.true.
      if( applyPSCut ) return
      call TopDecay(Top(plus_,plus_,1), DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**minus_))
      ColorCorr = -0.5d0
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))

      PAux(1:4) = z*MomExt(1:4,i) - (1d0-z)*MomExt(1:4,j)
      call pol_mless(dcmplx(MomExtTd(1:4,11)),-1,GluPol(1:4))
      xm = dcmplx(PAux(1:4)).dot.GluPol(1:4)
      call pol_mless(dcmplx(MomExtTd(1:4,11)),+1,GluPol(1:4))
      xp = dcmplx(PAux(1:4)).dot.GluPol(1:4)

      Splitting = 1d0 - 4d0/sij*dconjg(xp)*xp
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 1d0 - 4d0/sij*dconjg(xm)*xm
      Dipole(minus_,minus_) = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting =-4d0/sij*dconjg(xp)*xm
      Dipole(plus_,minus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting =-4d0/sij*dconjg(xm)*xp
      Dipole(minus_,plus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

ENDIF! dipoles
END SUBROUTINE







SUBROUTINE DipolesDKJ_Contr14(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
use ModProcess
use ModMisc
use ModParameters
use modKinematics
implicit none
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,5:11)
real(8) :: sij,ColorCorr,Ria,fAux,PAux(1:4),PAux2(1:4)
real(8) :: z,y,RunFactor,ymax
integer :: i,j,k,nDipole
logical :: applyPSCut
integer :: NBin(1:NumMaxHisto),TopCorr
complex(8) :: Dipole(0:1,0:1),Splitting
complex(8) :: xp,xm,GluPol(1:4)
integer, parameter :: plus_ = 0, minus_= 1
real(8), parameter :: CA=3d0, CF=4d0/3d0
integer, parameter :: NumDKPrimAmp = 2
type(Particle) :: ATop(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)! indices: Contribution, GluHel1, GluHel2, ColFacDK
type(Particle) ::  Top(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)


    applyPSCut = .false.
    TopCorr = Top_
    call TopDecay(ATop(1,1,1),DK_LO,MomExt(1:4,5:7))! evaluate this spinor once and for all dipoles before acceptance cuts
    RunFactor = RunAlphaS(NLOParam,MuRen)



IF( nDipole.eq.1 ) then
      i=11; j=12; k=8;
      call DKDip_FFmapping(MomExt(1:4,8:12),i-7,j-7,k-7,y,z,MomExtTd(1:4,8:11))
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExtTd(1:4,8),MomExtTd(1:4,9),MomExtTd(1:4,10)/),applyPSCut,NBin)
      if( y.gt.alpha_DKTff ) applyPSCut=.true.
      if( applyPSCut ) return
      call TopDecay(Top(plus_,plus_,1), DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**minus_))
      ColorCorr = -0.5d0
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))

      PAux(1:4) = z*MomExt(1:4,i) - (1d0-z)*MomExt(1:4,j)
      call pol_mless(dcmplx(MomExtTd(1:4,11)),-1,GluPol(1:4))
      xm = dcmplx(PAux(1:4)).dot.GluPol(1:4)
      call pol_mless(dcmplx(MomExtTd(1:4,11)),+1,GluPol(1:4))
      xp = dcmplx(PAux(1:4)).dot.GluPol(1:4)

      Splitting = 1d0 - 4d0/sij*dconjg(xp)*xp
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 1d0 - 4d0/sij*dconjg(xm)*xm
      Dipole(minus_,minus_) = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting =-4d0/sij*dconjg(xp)*xm
      Dipole(plus_,minus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting =-4d0/sij*dconjg(xm)*xp
      Dipole(minus_,plus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij



ELSEIF( nDipole.eq.2 ) then
      i=8; j=11; k=12;
      call DKDip_FFmapping(MomExt(1:4,8:12),i-7,j-7,k-7,y,z,MomExtTd(1:4,8:11))
      call swapMom(MomExtTd(1:4,8),MomExtTd(1:4,11))! because (8,11) is recombined into 8 which is the "reduced: gluon that has to be in position 11
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExtTd(1:4,8),MomExtTd(1:4,9),MomExtTd(1:4,10)/),applyPSCut,NBin)
      if( y.gt.alpha_DKTff ) applyPSCut=.true.
      if( applyPSCut ) return
      call TopDecay(Top(plus_,plus_,1), DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,plus_,1),DKJ_LO_T,MomExtTd(1:4,8:11),GluonHel=int((-1)**minus_))
      ColorCorr = -0.5d0
      sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))

      PAux(1:4) = z*MomExt(1:4,i) - (1d0-z)*MomExt(1:4,j)
      call pol_mless(dcmplx(MomExtTd(1:4,11)),-1,GluPol(1:4))
      xm = dcmplx(PAux(1:4)).dot.GluPol(1:4)
      call pol_mless(dcmplx(MomExtTd(1:4,11)),+1,GluPol(1:4))
      xp = dcmplx(PAux(1:4)).dot.GluPol(1:4)

      Splitting = 1d0 - 4d0/sij*dconjg(xp)*xp
      Dipole(plus_,plus_)   = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting = 1d0 - 4d0/sij*dconjg(xm)*xm
      Dipole(minus_,minus_) = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting =-4d0/sij*dconjg(xp)*xm
      Dipole(plus_,minus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij

      Splitting =-4d0/sij*dconjg(xm)*xp
      Dipole(minus_,plus_)  = 8d0*DblPi*(alpha_s*RunFactor)*ColorCorr*Splitting/sij


ENDIF! dipoles
END SUBROUTINE







SUBROUTINE DipolesDKJ_Contr15(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)! originally 13
use ModProcess
use ModMisc
use ModParameters
use modKinematics
use ModTTBJ_NLODKW
implicit none
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,5:11)
real(8) :: sij,ColorCorr,Ria,fAux,PAux(1:4),PAux2(1:4)
real(8) :: z,y,RunFactor,ymax
integer :: i,j,k,nDipole
logical :: applyPSCut
integer :: NBin(1:NumMaxHisto),TopCorr
complex(8) :: Dipole(0:1,0:1),Splitting, dip
complex(8) :: xp,xm,GluPol(1:4)
integer, parameter :: plus_ = 0, minus_= 1
real(8), parameter :: CA=3d0, CF=4d0/3d0
integer, parameter :: NumDKPrimAmp = 2
type(Particle) :: ATop(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)! indices: Contribution, GluHel1, GluHel2, ColFacDK
type(Particle) ::  Top(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)

call Wto3Jet_init
! DIPOLE for W-decay with W- -> ub d q qb, q = b,s,c
    applyPSCut = .false.
    TopCorr = ATop_
    call TopDecay(Top(1,1,1),DK_LO,MomExt(1:4,10:12))! evaluate this spinor once and for all dipoles before acceptance cuts
    RunFactor = RunAlphaS(NLOParam,MuRen)
    Dipole(0:1,0:1) = czero

   call WmDipolesQQB(ndipole,MomExt(1:4,6:9),dip,1,1,MomExtTd(1:4,6:8))

   MomExtTd(1:4,5) = MomExt(1:4,5) ! Set bottom-mom for Tilde-kienmatics in TopDecay

   call swapmom(MomExtTd(1:4,6),MomExtTd(1:4,7))! MARKUS:  introduced this to assure correct limit

   if(dip .eq. czero) return
   call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
   if( applyPSCut ) return

   call TopDecay(ATop(plus_,plus_,1), DKJ_LO_Q,MomExtTd(1:4,5:8),GluonHel=int((-1)**plus_))
   call TopDecay(ATop(minus_,plus_,1),DKJ_LO_Q,MomExtTd(1:4,5:8),GluonHel=int((-1)**minus_))
   Dipole(plus_,plus_)   = -dip

   call WmDipolesQQB(ndipole,MomExt(1:4,6:9),dip,1,-1,MomExtTd(1:4,6:8))
   Dipole(plus_,minus_) = -dip

   call WmDipolesQQB(ndipole,MomExt(1:4,6:9),dip,-1,1,MomExtTd(1:4,6:8))
   Dipole(minus_,plus_) = -dip

   call WmDipolesQQB(ndipole,MomExt(1:4,6:9),dip,-1,-1,MomExtTd(1:4,6:8))
   Dipole(minus_,minus_) = -dip


!MARKUS:: introduced fudge factor of 0.5.... has to be clarified!!
dipole(:,:)= 0.5d0*dipole(:,:)



If(ndipole .gt. 2) then
   write(*,*) "ERROR IN DipolesDKJ_Contr13 (SCHARF's), nDipole = ", nDipole
endif

END SUBROUTINE






SUBROUTINE DipolesDKJ_Contr16(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)! originally 15
use ModProcess
use ModMisc
use ModParameters
use modKinematics
use ModTTBJ_NLODKW
implicit none
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,5:11)
real(8) :: sij,ColorCorr,Ria,fAux,PAux(1:4),PAux2(1:4)
real(8) :: z,y,RunFactor,ymax
integer :: i,j,k,nDipole
logical :: applyPSCut
integer :: NBin(1:NumMaxHisto),TopCorr
complex(8) :: Dipole(0:1,0:1),Splitting, dip
complex(8) :: xp,xm,GluPol(1:4)
integer, parameter :: plus_ = 0, minus_= 1
real(8), parameter :: CA=3d0, CF=4d0/3d0
integer, parameter :: NumDKPrimAmp = 2
type(Particle) :: ATop(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)! indices: Contribution, GluHel1, GluHel2, ColFacDK
type(Particle) ::  Top(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)

call Wto3Jet_init
! DIPOLE for W-decay with W- -> ub d q qb, q = b,s,c
    applyPSCut = .false.
    TopCorr = ATop_
    call TopDecay(Top(1,1,1),DK_LO,MomExt(1:4,10:12))! evaluate this spinor once and for all dipoles before acceptance cuts
    RunFactor = RunAlphaS(NLOParam,MuRen)
    Dipole(0:1,0:1) = czero
    call WmDipolesUUB(ndipole,MomExt(1:4,6:9),dip,1,1,MomExtTd(1:4,6:8))

   MomExtTd(1:4,5) = MomExt(1:4,5) ! Set bottom-mom for Tilde-kienmatics in TopDecay

   call swapmom(MomExtTd(1:4,6),MomExtTd(1:4,7))! MARKUS:  introduced this to assure correct limit

   if(dip .eq. czero) return
   call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
   if( applyPSCut ) return

   call TopDecay(ATop(plus_,plus_,1), DKJ_LO_Q,MomExtTd(1:4,5:8),GluonHel=int((-1)**plus_))
   call TopDecay(ATop(minus_,plus_,1),DKJ_LO_Q,MomExtTd(1:4,5:8),GluonHel=int((-1)**minus_))
   Dipole(plus_,plus_)   = -dip

   call WmDipolesUUB(ndipole,MomExt(1:4,6:9),dip,1,-1,MomExtTd(1:4,6:8))

   Dipole(plus_,minus_) = -dip
   call WmDipolesUUB(ndipole,MomExt(1:4,6:9),dip,-1,1,MomExtTd(1:4,6:8))

   Dipole(minus_,plus_) = -dip
   call WmDipolesUUB(ndipole,MomExt(1:4,6:9),dip,-1,-1,MomExtTd(1:4,6:8))

   Dipole(minus_,minus_) = -dip
If(ndipole .gt. 4) then
   write(*,*) "ERROR IN DipolesDKJ_Contr15 (SCHARF's), nDipole = ", nDipole
endif

END SUBROUTINE





SUBROUTINE DipolesDKJ_Contr17(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)
use ModProcess
use ModMisc
use ModParameters
use modKinematics
use ModTTBJ_NLODKW
implicit none
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,5:11)
real(8) :: sij,ColorCorr,Ria,fAux,PAux(1:4),PAux2(1:4)
real(8) :: z,y,RunFactor,ymax
integer :: i,j,k,nDipole
logical :: applyPSCut
integer :: NBin(1:NumMaxHisto),TopCorr
complex(8) :: Dipole(0:1,0:1),Splitting, dip
complex(8) :: xp,xm,GluPol(1:4)
integer, parameter :: plus_ = 0, minus_= 1
real(8), parameter :: CA=3d0, CF=4d0/3d0
integer, parameter :: NumDKPrimAmp = 2
type(Particle) :: ATop(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)! indices: Contribution, GluHel1, GluHel2, ColFacDK
type(Particle) ::  Top(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)

call Wto3Jet_init
! DIPOLE for W-decay with W- -> ub d q qb, q = b,s,c
    applyPSCut = .false.
    TopCorr = ATop_
    call TopDecay(Top(1,1,1),DK_LO,MomExt(1:4,10:12))! evaluate this spinor once and for all dipoles before acceptance cuts
    RunFactor = RunAlphaS(NLOParam,MuRen)
    Dipole(0:1,0:1) = czero

   call WmDipolesDDB(ndipole,MomExt(1:4,6:9),dip,1,1,MomExtTd(1:4,6:8))

   MomExtTd(1:4,5) = MomExt(1:4,5) ! Set bottom-mom for Tilde-kienmatics in TopDecay

   call swapmom(MomExtTd(1:4,6),MomExtTd(1:4,7))! MARKUS:  introduced this to assure correct limit

   if(dip .eq. czero) return
   call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExtTd(1:4,6),MomExtTd(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
   if( applyPSCut ) return

   call TopDecay(ATop(plus_,plus_,1), DKJ_LO_Q,MomExtTd(1:4,5:8),GluonHel=int((-1)**plus_))
   call TopDecay(ATop(minus_,plus_,1),DKJ_LO_Q,MomExtTd(1:4,5:8),GluonHel=int((-1)**minus_))
   Dipole(plus_,plus_)   = -dip

   call WmDipolesDDB(ndipole,MomExt(1:4,6:9),dip,1,-1,MomExtTd(1:4,6:8))

   Dipole(plus_,minus_) = -dip
   call WmDipolesDDB(ndipole,MomExt(1:4,6:9),dip,-1,1,MomExtTd(1:4,6:8))

   Dipole(minus_,plus_) = -dip
   call WmDipolesDDB(ndipole,MomExt(1:4,6:9),dip,-1,-1,MomExtTd(1:4,6:8))

   Dipole(minus_,minus_) = -dip
If(ndipole .gt. 4) then
   write(*,*) "ERROR IN DipolesDKJ_Contr15 (SCHARF's), nDipole = ", nDipole
endif

END SUBROUTINE






SUBROUTINE DipolesDKJ_Contr18(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)! originally 14
use ModProcess
use ModMisc
use ModParameters
use ModKinematics
use ModTTBJ_NLODKW
implicit none
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,5:11)
real(8) :: sij,ColorCorr,Ria,fAux,PAux(1:4),PAux2(1:4)
real(8) :: z,y,RunFactor,ymax, muS
logical :: applyPSCut
integer :: NBin(1:NumMaxHisto),TopCorr
integer :: i,j,k,nDipole
complex(8) :: Dipole(0:1,0:1),Splitting, dip
complex(8) :: xp,xm,GluPol(1:4)
integer, parameter :: plus_ = 0, minus_= 1
real(8), parameter :: CA=3d0, CF=4d0/3d0
real(8), parameter :: Tt_x_Tb = 1d0/2d0*CA-CF
real(8), parameter :: Tt_x_Tg =-1d0/2d0*CA
real(8), parameter :: Tb_x_Tg =-1d0/2d0*CA
integer, parameter :: NumDKPrimAmp = 2
integer :: xe, LO_NLO
type(Particle) :: ATop(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)! indices: Contribution, GluHel1, GluHel2, ColFacDK
type(Particle) ::  Top(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)
call Wto3Jet_init
    applyPSCut = .false.
    TopCorr = Top_
    call TopDecay(ATop(1,1,1),DK_LO,MomExt(1:4,5:7))! evaluate this spinor once and for all dipoles before acceptance cuts
    RunFactor = RunAlphaS(NLOParam,MuRen)
    Dipole(0:1,0:1) = czero
    call WmDipolesQQB(ndipole,MomExt(1:4,9:12),dip,1,1,MomExtTd(1:4,9:11))

    MomExtTd(1:4,8) = MomExt(1:4,8) ! Set b-momentum for TOP-decays with tilde kinematics
    if(dip .eq. czero) return

   call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExtTd(1:4,9),MomExtTd(1:4,10)/),applyPSCut,NBin)

   if( applyPSCut ) return

   call TopDecay(Top(plus_,plus_,1), DKJ_LO_Q,MomExtTd(1:4,8:11),GluonHel=int((-1)**plus_))
   call TopDecay(Top(minus_,plus_,1),DKJ_LO_Q,MomExtTd(1:4,8:11),GluonHel=int((-1)**minus_))
   Dipole(plus_,plus_)   = -dip

   call WmDipolesQQB(ndipole,MomExt(1:4,9:12),dip,1,-1,MomExtTd(1:4,9:11))
   Dipole(plus_,minus_) = -dip
   call WmDipolesQQB(ndipole,MomExt(1:4,9:12),dip,-1,1,MomExtTd(1:4,9:11))
   Dipole(minus_,plus_) = -dip
   call WmDipolesQQB(ndipole,MomExt(1:4,9:12),dip,-1,-1,MomExtTd(1:4,9:11))
   Dipole(minus_,minus_) = -dip



!MARKUS:: introduced fudge factor of 0.5.... has to be clarified!!
dipole(:,:)= 0.5d0*dipole(:,:)


If(ndipole .gt. 2) then
   write(*,*) "ERROR IN DipolesDKJ_Contr9 (SCHARF's), nDipole = ", nDipole
endif

END SUBROUTINE






SUBROUTINE DipolesDKJ_Contr19(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)! originally 16
use ModProcess
use ModMisc
use ModParameters
use ModKinematics
use ModTTBJ_NLODKW
implicit none
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,5:11)
real(8) :: sij,ColorCorr,Ria,fAux,PAux(1:4),PAux2(1:4)
real(8) :: z,y,RunFactor,ymax, muS
logical :: applyPSCut
integer :: NBin(1:NumMaxHisto),TopCorr
integer :: i,j,k,nDipole
complex(8) :: Dipole(0:1,0:1),Splitting, dip
complex(8) :: xp,xm,GluPol(1:4)
integer, parameter :: plus_ = 0, minus_= 1
real(8), parameter :: CA=3d0, CF=4d0/3d0
real(8), parameter :: Tt_x_Tb = 1d0/2d0*CA-CF
real(8), parameter :: Tt_x_Tg =-1d0/2d0*CA
real(8), parameter :: Tb_x_Tg =-1d0/2d0*CA
integer, parameter :: NumDKPrimAmp = 2
integer :: xe, LO_NLO
type(Particle) :: ATop(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)! indices: Contribution, GluHel1, GluHel2, ColFacDK
type(Particle) ::  Top(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)
call Wto3Jet_init
    applyPSCut = .false.
    TopCorr = Top_
    call TopDecay(ATop(1,1,1),DK_LO,MomExt(1:4,5:7))! evaluate this spinor once and for all dipoles before acceptance cuts
    RunFactor = RunAlphaS(NLOParam,MuRen)
    Dipole(0:1,0:1) = czero
    call WmDipolesUUB(ndipole,MomExt(1:4,9:12),dip,1,1,MomExtTd(1:4,9:11))

    MomExtTd(1:4,8) = MomExt(1:4,8) ! Set b-momentum for TOP-decays with tilde kinematics
    if(dip .eq. czero) return

   call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExtTd(1:4,9),MomExtTd(1:4,10)/),applyPSCut,NBin)

   if( applyPSCut ) return

   call TopDecay(Top(plus_,plus_,1), DKJ_LO_Q,MomExtTd(1:4,8:11),GluonHel=int((-1)**plus_))
   call TopDecay(Top(minus_,plus_,1),DKJ_LO_Q,MomExtTd(1:4,8:11),GluonHel=int((-1)**minus_))
   Dipole(plus_,plus_)   = -dip

   call WmDipolesUUB(ndipole,MomExt(1:4,9:12),dip,1,-1,MomExtTd(1:4,9:11))
   Dipole(plus_,minus_) = -dip
   call WmDipolesUUB(ndipole,MomExt(1:4,9:12),dip,-1,1,MomExtTd(1:4,9:11))
   Dipole(minus_,plus_) = -dip
   call WmDipolesUUB(ndipole,MomExt(1:4,9:12),dip,-1,-1,MomExtTd(1:4,9:11))
   Dipole(minus_,minus_) = -dip

If(ndipole .gt. 4) then
   write(*,*) "ERROR IN DipolesDKJ_Contr9 (SCHARF's), nDipole = ", nDipole
endif

END SUBROUTINE






SUBROUTINE DipolesDKJ_Contr20(nDipole,ATop,Top,TopCorr,MomExt,Dipole,applyPSCut,NBin)! originally 18
use ModProcess
use ModMisc
use ModParameters
use ModKinematics
use ModTTBJ_NLODKW
implicit none
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,5:11)
real(8) :: sij,ColorCorr,Ria,fAux,PAux(1:4),PAux2(1:4)
real(8) :: z,y,RunFactor,ymax, muS
logical :: applyPSCut
integer :: NBin(1:NumMaxHisto),TopCorr
integer :: i,j,k,nDipole
complex(8) :: Dipole(0:1,0:1),Splitting, dip
complex(8) :: xp,xm,GluPol(1:4)
integer, parameter :: plus_ = 0, minus_= 1
real(8), parameter :: CA=3d0, CF=4d0/3d0
real(8), parameter :: Tt_x_Tb = 1d0/2d0*CA-CF
real(8), parameter :: Tt_x_Tg =-1d0/2d0*CA
real(8), parameter :: Tb_x_Tg =-1d0/2d0*CA
integer, parameter :: NumDKPrimAmp = 2
integer :: xe, LO_NLO
type(Particle) :: ATop(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)! indices: Contribution, GluHel1, GluHel2, ColFacDK
type(Particle) ::  Top(plus_:minus_,plus_:minus_,1:NumDKPrimAmp)
call Wto3Jet_init
    applyPSCut = .false.
    TopCorr = Top_
    call TopDecay(ATop(1,1,1),DK_LO,MomExt(1:4,5:7))! evaluate this spinor once and for all dipoles before acceptance cuts
    RunFactor = RunAlphaS(NLOParam,MuRen)
    Dipole(0:1,0:1) = czero
    call WmDipolesDDB(ndipole,MomExt(1:4,9:12),dip,1,1,MomExtTd(1:4,9:11))

    MomExtTd(1:4,8) = MomExt(1:4,8) ! Set b-momentum for TOP-decays with tilde kinematics
    if(dip .eq. czero) return

   call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExtTd(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExtTd(1:4,9),MomExtTd(1:4,10)/),applyPSCut,NBin)

   if( applyPSCut ) return

   call TopDecay(Top(plus_,plus_,1), DKJ_LO_Q,MomExtTd(1:4,8:11),GluonHel=int((-1)**plus_))
   call TopDecay(Top(minus_,plus_,1),DKJ_LO_Q,MomExtTd(1:4,8:11),GluonHel=int((-1)**minus_))
   Dipole(plus_,plus_)   = -dip

   call WmDipolesDDB(ndipole,MomExt(1:4,9:12),dip,1,-1,MomExtTd(1:4,9:11))
   Dipole(plus_,minus_) = -dip
   call WmDipolesDDB(ndipole,MomExt(1:4,9:12),dip,-1,1,MomExtTd(1:4,9:11))
   Dipole(minus_,plus_) = -dip
   call WmDipolesDDB(ndipole,MomExt(1:4,9:12),dip,-1,-1,MomExtTd(1:4,9:11))
   Dipole(minus_,minus_) = -dip

If(ndipole .gt. 4) then
   write(*,*) "ERROR IN DipolesDKJ_Contr9 (SCHARF's), nDipole = ", nDipole
endif

END SUBROUTINE











SUBROUTINE DKDip_FFmapping(MomExt,i,j,k,y,z,MomExtTd)!  i has to be smaller than j
use ModMisc
implicit none
real(8) :: MomExt(:,:)!   order: Bot,Lep,Lep,Glu,Glu
real(8) :: MomExtTd(:,:)! order: Bot,Lep,Lep,Glu
real(8) :: y,z
integer :: i,j,k,NReal


          NReal = size(MomExt(:,:),2)
          MomExtTd(1:4,1:NReal-1)=MomExt(1:4,1:NReal-1)

          y = (MomExt(1:4,i).dot.MomExt(1:4,j))/((MomExt(1:4,i).dot.MomExt(1:4,j))+(MomExt(1:4,i).dot.MomExt(1:4,k))+(MomExt(1:4,j).dot.MomExt(1:4,k)))
          z = (MomExt(1:4,i).dot.MomExt(1:4,k))/((MomExt(1:4,i).dot.MomExt(1:4,k))+(MomExt(1:4,k).dot.MomExt(1:4,j)))
          if( k.eq.NReal) then! this means i<j<k
              MomExtTd(1:4,NReal-1) = 1d0/(1d0-y)*MomExt(1:4,k)
          else
              MomExtTd(1:4,k) = 1d0/(1d0-y)*MomExt(1:4,k)
          endif
          MomExtTd(1:4,i) = MomExt(1:4,i)+MomExt(1:4,j)-y/(1d0-y)*MomExt(1:4,k)


! print *, "FF Mapping",i,j,k
! print *, dsqrt((MomExtTd(1:4,5)+MomExtTd(1:4,6)+MomExtTd(1:4,7)+MomExtTd(1:4,8)).dot.(MomExtTd(1:4,5)+MomExtTd(1:4,6)+MomExtTd(1:4,7)+MomExtTd(1:4,8)))
! print *, dsqrt((MomExtTd(1:4,6)+MomExtTd(1:4,7)).dot.(MomExtTd(1:4,6)+MomExtTd(1:4,7)))
! print *, (MomExt(1:4,5)+MomExt(1:4,6)+MomExt(1:4,7)+MomExt(1:4,8)+MomExt(1:4,9))-(MomExtTd(1:4,5)+MomExtTd(1:4,6)+MomExtTd(1:4,7)+MomExtTd(1:4,8))
! print *, MomExtTd(1:4,5).dot.MomExtTd(1:4,5)
! print *, MomExtTd(1:4,6).dot.MomExtTd(1:4,6)
! print *, MomExtTd(1:4,7).dot.MomExtTd(1:4,7)
! print *, MomExtTd(1:4,8).dot.MomExtTd(1:4,8)
! print *, "----"
! pause

END SUBROUTINE





SUBROUTINE DKDip_FImapping(MomExt,i,j,y,z,Ria,ymax,MomExtTd)
use ModMisc
use ModParameters
implicit none
real(8) :: MomExt(:,:)!   order: Bot,Lep,Lep,Glu,Glu
real(8) :: MomExtTd(:,:)! order: Bot,Lep,Lep,Glu
real(8) :: MomTop(1:4),MomK(1:4)
real(8) :: z,y,Ria,ymax
integer :: i,j,n,NReal


         NReal = size(MomExt(:,:),2)
         MomTop(1:4)=0d0; MomK(1:4) = 0d0
         do n=1,NReal
              MomTop(1:4) = MomTop(1:4) + MomExt(1:4,n)
              if( n.ne.i .and. n.ne.j ) MomK(1:4) = MomK(1:4) + MomExt(1:4,n)
         enddo
         z = 1d0-(    (MomTop(1:4).dot.MomExt(1:4,j))/((MomExt(1:4,i).dot.MomTop(1:4))+(MomTop(1:4).dot.MomExt(1:4,j))-(MomExt(1:4,i).dot.MomExt(1:4,j)))     )
         Ria = 2d0*(MomExt(1:4,j).dot.MomTop(1:4))/(1d0-z)! = mt^2 * (1-r^2)
         y = (MomExt(1:4,i).dot.MomExt(1:4,j))*2d0/m_Top**2 / (1d0-dsqrt(1d0-Ria/m_Top**2))**2
         ymax = (1d0+dsqrt(1d0-Ria/m_Top**2))**2/(z+(1d0-z)*(1d0-Ria/m_Top**2))*z*(1d0-z)

         MomExtTd(1:4,min(i,j)) = MomTop(1:4)
         do n=1,NReal
            if( n.ne.i .and. n.ne.j ) then
                if( n.ne.NReal ) then
                      call DipoleTrafo(MomTop(1:4),MomK(1:4),MomExt(1:4,n),MomExtTd(1:4,n))
                      MomExtTd(1:4,min(i,j)) = MomExtTd(1:4,min(i,j)) - MomExtTd(1:4,n)
                else
                      call DipoleTrafo(MomTop(1:4),MomK(1:4),MomExt(1:4,n),MomExtTd(1:4,NReal-1))
                      MomExtTd(1:4,min(i,j)) = MomExtTd(1:4,min(i,j)) - MomExtTd(1:4,NReal-1)
                endif
            endif
         enddo

! print *, "FI Mapping",i,j
! print *, dsqrt((MomExtTd(1:4,5)+MomExtTd(1:4,6)+MomExtTd(1:4,7)+MomExtTd(1:4,8)).dot.(MomExtTd(1:4,5)+MomExtTd(1:4,6)+MomExtTd(1:4,7)+MomExtTd(1:4,8)))
! print *, dsqrt((MomExtTd(1:4,6)+MomExtTd(1:4,7)).dot.(MomExtTd(1:4,6)+MomExtTd(1:4,7)))
! print *, (MomExt(1:4,5)+MomExt(1:4,6)+MomExt(1:4,7)+MomExt(1:4,8)+MomExt(1:4,9))-(MomExtTd(1:4,5)+MomExtTd(1:4,6)+MomExtTd(1:4,7)+MomExtTd(1:4,8))
! print *, MomExtTd(1:4,5).dot.MomExtTd(1:4,5)
! print *, MomExtTd(1:4,6).dot.MomExtTd(1:4,6)
! print *, MomExtTd(1:4,7).dot.MomExtTd(1:4,7)
! print *, MomExtTd(1:4,8).dot.MomExtTd(1:4,8)
! print *, "----"

END SUBROUTINE











END MODULE
