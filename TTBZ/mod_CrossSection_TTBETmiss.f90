MODULE ModCrossSection_TTBETmiss
use ModExoticDecay
implicit none

integer,private,parameter :: NumMaxHisto=45

integer,private,parameter :: RemoveHTopClosedLoop = 1!   0=remove, 1=include   closed heavy HT loop


 contains




FUNCTION EvalCS_1L_HtHtbgg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModIntDipoles_GGHTHTBG
implicit none
real(8) ::  EvalCS_1L_HtHtbgg,yRnd(1:VegasMxDim),VgsWgt,HOp(1:3)
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1),FermionLoopPartAmp(7:8,-2:1)
integer :: iHel,jHel,kHel,iPrimAmp,jPrimAmp
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5,ISFac
real(8) :: MomExt(1:4,1:15),MomP(1:4,1:4)
logical :: applyPSCut
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,AccPoles
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2),xE
integer :: NBin(1:NumMaxHisto),NHisto,nHel(1:2),ParityFlip,BH1Hel,BH2Hel,BHMaxHel,iQuark,iCut
include 'misc/global_import'
include 'vegas_common.f'

IF( XTopDecays.ge.1 ) THEN
  ParityFlip=1
ELSE
  ParityFlip=2
ENDIF

  EvalCS_1L_HtHtbgg = 0d0
  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_HTop ) then
      EvalCS_1L_HtHtbgg = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

!   MuRen = 0.5d0*Ehat
!   MuFac = MuRen

   call EvalPhaseSpace_2to2HT(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt)! HTbar, HT
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
   ISFac = MomCrossing(MomExt)
   BHMaxHel = -1
   IF(XTOPDECAYS.EQ.1) THEN
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  BH top
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,4),yRnd(7:8),MomExt(1:4,10:11),PSWgt3)
      m_Top = m_SMTop
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd( 9:12),MomExt(1:4,7:9),PSWgt4)! bot lep neu
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(13:16),MomExt(1:4,12:14),PSWgt5)
      m_Top = m_HTop
      PSWgt = PSWgt * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
      BHMaxHel = +1     
      nHel(1:2) = getHelicity(yrnd(17))

   ELSEIF(XTOPDECAYS.EQ.2) THEN
      call EvalPhasespace_HTopDK(HT_A0_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  A0 top
      call EvalPhasespace_HTopDK(HT_A0_T,MomExt(1:4,4),yRnd(7:8),MomExt(1:4,10:11),PSWgt3)
      m_Top = m_SMTop
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd( 9:12),MomExt(1:4,7:9),PSWgt4)! bot lep neu
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(13:16),MomExt(1:4,12:14),PSWgt5)
      m_Top = m_HTop
      PSWgt = PSWgt * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
      call HTopA0Decay(ExtParticle(1),DKX_HTA0_LO,MomExt(1:4,5:9))
      call HTopA0Decay(ExtParticle(2),DKX_HTA0_LO,MomExt(1:4,10:14))
      nHel(1:2) = getHelicity(yrnd(17))
  
   ELSEIF(XTOPDECAYS.EQ.0) THEN
      nHel(1:2) = getHelicity(yrnd(5))
   ENDIF

   call Kinematics_TTbarETmiss(.false.,MomExt,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_1L_HtHtbgg = 0d0
      return
   endif

   call InitCurrCache()
   call SetPropagators()

   call SetPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   RunFactor = RunAlphaS(NLOParam,MuRen)
   PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))

   LO_Res_Unpol             = (0d0,0d0)
   NLO_Res_Unpol(-2:1)      = (0d0,0d0)
   NLO_Res_Unpol_Ferm(-2:1) = (0d0,0d0)
!------------ LO ----------------
IF( Correction.EQ.0 ) THEN
   do iHel=nHel(1),nHel(2)
   do BH1Hel=-1,BHMaxHel
   do BH2Hel=-1,BHMaxHel
      IF( XTOPDECAYS.EQ.1 ) THEN
         call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExt(1:4,5:9))
         call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExt(1:4,10:14))
      ENDIF

      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo
      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * ParityFlip*BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
   enddo!helicity loop
   enddo!helicity loop
   enddo!helicity loop

!------------ 1 LOOP --------------
ELSEIF( Correction.EQ.1 ) THEN
   do iHel=nHel(1),nHel(2)
   do BH1Hel=-1,BHMaxHel
   do BH2Hel=-1,BHMaxHel
      IF( XTOPDECAYS.EQ.1 ) THEN
         call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExt(1:4,5:9))
         call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExt(1:4,10:14))
      ENDIF

      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      do iPrimAmp=1,6
          call EvalTree(BornAmps(iPrimAmp))
      enddo
      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * ParityFlip*BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
!------------ bosonic loops --------------
       do iPrimAmp=1,6
          call SetKirill(PrimAmps(iPrimAmp))
          call PentCut(PrimAmps(iPrimAmp))
          call QuadCut(PrimAmps(iPrimAmp))
          call TripCut(PrimAmps(iPrimAmp))
          call DoubCut(PrimAmps(iPrimAmp))
          call SingCut(PrimAmps(iPrimAmp))
          call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
          call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen**2)
          PrimAmps(iPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(iPrimAmp)%Result(-2:1)
          call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,2,rdiv(2),rdiv(1))
!           call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
! pause

!DEC$ IF (_QuadPrecImpr==1)
          AccPoles = CheckPoles(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv(1:2))
          if( AccPoles.gt.1d-4 ) then
              coeff4_128(:,:) = qcmplx( coeff4(:,:) )
              coeff5_128(:,:) = qcmplx( coeff5(:,:) )
              call TripCut_128(PrimAmps(iPrimAmp))
              call DoubCut_128(PrimAmps(iPrimAmp))
              call SingCut_128(PrimAmps(iPrimAmp))
              call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
              call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen*2)
              PrimAmps(iPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(iPrimAmp)%Result(-2:1)
              AccPoles = CheckPoles(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv(1:2))
              if( AccPoles.gt.5d-2 ) then
                  print *, "SKIP",AccPoles
                  call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
                  EvalCS_1L_HtHtbgg = 0d0
                  SkipCounter = SkipCounter + 1
                  return
              endif
          endif
!DEC$ ENDIF
      enddo

      NLO_Res_Pol(-2:1) = (0d0,0d0)
      do jPrimAmp=1,6
      do iPrimAmp=1,NumBornAmps
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + ParityFlip*Col1L_ttbgg(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(PrimAmps(jPrimAmp)%Result(-2:1)) )
      enddo
      enddo
      NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1) + NLO_Res_Pol(-2:1)


! ------------ fermionic loops --------------
      do iPrimAmp=7,12
          call SetKirill(PrimAmps(iPrimAmp))
          call PentCut(PrimAmps(iPrimAmp))
          call QuadCut(PrimAmps(iPrimAmp))
          call TripCut(PrimAmps(iPrimAmp))
          call DoubCut(PrimAmps(iPrimAmp))
          call SingCut(PrimAmps(iPrimAmp))
          call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
          PrimAmps(iPrimAmp)%Result(-2:1) = -(0d0,1d0)*PrimAmps(iPrimAmp)%Result(-2:1) !minus is from closed fermion loop
      enddo
      FermionLoopPartAmp(7,-2:1) = Nf_light*PrimAmps(7)%Result(-2:1) + PrimAmps(9)%Result(-2:1)  + PrimAmps(11)%Result(-2:1)*RemoveHTopClosedLoop
      FermionLoopPartAmp(8,-2:1) = Nf_light*PrimAmps(8)%Result(-2:1) + PrimAmps(10)%Result(-2:1) + PrimAmps(12)%Result(-2:1)*RemoveHTopClosedLoop


      NLO_Res_Pol(-2:1) = (0d0,0d0)
      do jPrimAmp=7,8
      do iPrimAmp=1,NumBornAmps
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + ParityFlip*Col1L_ttbgg(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result*dconjg(FermionLoopPartAmp(jPrimAmp,-2:1)) )
      enddo
      enddo
      NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) + NLO_Res_Pol(-2:1)
   enddo! helicity loop
   enddo! helicity loop
   enddo! helicity loop
ENDIF


IF( Correction.EQ.0 ) THEN
!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * WidthExpansion
   EvalCS_1L_HtHtbgg = LO_Res_Unpol * PreFac

ELSEIF( Correction.EQ.1 ) THEN
!  overall normalization: (4*Pi)^eps/Gamma(1-eps)
!  CT contributions                           ! beta        !top WFRC
   NLO_Res_UnPol(-1) = NLO_Res_UnPol(-1) + (-11d0/3d0*3d0 - 3d0*4d0/3d0 )*LO_Res_Unpol
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-3d0*4d0/3d0)*2d0*dlog(MuRen/m_top)*LO_Res_Unpol  ! finite log(mu) contrib. from  top WFRC

   NLO_Res_UnPol_Ferm(-1) = NLO_Res_UnPol_Ferm(-1) - (-2d0/3d0*Nf_light)*LO_Res_Unpol

   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-5d0/2d0*8d0/3d0 )*LO_Res_Unpol   ! finite contribution from top WFRC's
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + LO_Res_Unpol  ! shift alpha_s^DR --> alpha_s^MSbar

!  factor out (Mu2/mTop**2)^eps
!    NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + NLO_Res_UnPol(-1)*2d0*dlog(m_top/MuRen) + NLO_Res_UnPol(-2)*dlog(m_top/MuRen)**2
!    NLO_Res_UnPol(-1) = NLO_Res_UnPol(-1) + NLO_Res_UnPol(-2)*2d0*dlog(m_top/MuRen)
!    NLO_Res_UnPol_Ferm(0) = NLO_Res_UnPol_Ferm(0) + NLO_Res_UnPol_Ferm(-1)*2d0*dlog(m_top/MuRen)

!  normalization
   LO_Res_Unpol = LO_Res_Unpol                         * ISFac * (alpha_s4Pi*RunFactor)**2
   NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1)           * ISFac * (alpha_s4Pi*RunFactor)**2 * alpha_sOver2Pi*RunFactor
   NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) * ISFac * (alpha_s4Pi*RunFactor)**2 * alpha_sOver2Pi*RunFactor


   EvalCS_1L_HtHtbgg = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac


ELSEIF( Correction.EQ.3 ) THEN


   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt

   if( XTopDecays.eq.0 ) then
      xE = yRnd(5)
   else
      xE = yRnd(17)
   endif
! xE=0.4d0

   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
   call EvalIntDipoles_GGHTHTBG((/MomExt(1:4,3),MomExt(1:4,4),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,5:14),xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
   EvalCS_1L_HtHtbgg = HOp(1)    * pdf(0,1)  * pdf(0,2)   &
                     + HOp(2)/xE * pdf_z(0,1)* pdf(0,2)   &
                     + HOp(3)/xE * pdf(0,1)  * pdf_z(0,2)

! checked for XTopDK=1 TopDK=1 and two mu-scales and for XTopDK=TopDK=0
! print *, "1L check",NLO_Res_UnPol(-2)                            * PreFac/(alpha_sOver2Pi*RunFactor)/(LO_Res_Unpol*PreFac)
! print *, "1L check",( NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1) ) * PreFac/(alpha_sOver2Pi*RunFactor)/(LO_Res_Unpol*PreFac)
! print *, "res1=", HOp(1)/(alpha_sOver2Pi*RunFactor)/(LO_Res_Unpol*PreFac)
! pause

ENDIF


   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_HtHtbgg)
   enddo
   EvalCounter = EvalCounter + 1

   EvalCS_1L_HtHtbgg = EvalCS_1L_HtHtbgg/VgsWgt

return
END FUNCTION










FUNCTION EvalCS_1L_HtHtbqqb(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModIntDipoles_QQBHTHTBG
use ModIntDipoles_QGHTHTBQ
implicit none
real(8) ::  EvalCS_1L_HtHtbqqb,yRnd(1:VegasMxDim),VgsWgt,xE,HOp(1:3)
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1)
complex(8) :: BosonicPartAmp(1:2,-2:1),FermionPartAmp(1:2,-2:1)
integer :: iHel,iPrimAmp,jPrimAmp
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5,ISFac
real(8) :: MomExt(1:4,1:15),MomP(1:4,1:4)
logical :: applyPSCut
real(8),parameter :: Nc=3d0
real(8) :: tau,eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac_a,PDFFac_b,PDFFac,AccPoles
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2),mydummy(1:10)
integer :: NHisto,NBin(1:NumMaxHisto),npdf,nHel(1:2),NRndHel,ParityFlip,BH1Hel,BH2Hel,BHMaxHel
include 'misc/global_import'
include "vegas_common.f"


IF( XTopDecays.ge.1 ) THEN
  ParityFlip=1
ELSE
  ParityFlip=2
ENDIF

  EvalCS_1L_HtHtbqqb=0d0
  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_HTop ) then
      EvalCS_1L_HtHtbqqb = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

!   MuRen = 0.5d0*Ehat
!   MuFac = MuRen


   call EvalPhaseSpace_2to2HT(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt)! HTbar, HT
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
   ISFac = MomCrossing(MomExt)

   BHMaxHel = -1
   IF(XTOPDECAYS.EQ.1) THEN
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  BH top
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,4),yRnd(7:8),MomExt(1:4,10:11),PSWgt3)

      m_Top = m_SMTop
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd( 9:12),MomExt(1:4,7:9),PSWgt4)! bot lep neu
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(13:16),MomExt(1:4,12:14),PSWgt5)
      m_Top = m_HTop

      PSWgt = PSWgt * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
      BHMaxHel = +1
      nHel(1:2) = getHelicity(yrnd(17))
   ELSEIF(XTOPDECAYS.EQ.2) THEN
      call EvalPhasespace_HTopDK(HT_A0_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  A0 top
      call EvalPhasespace_HTopDK(HT_A0_T,MomExt(1:4,4),yRnd(7:8),MomExt(1:4,10:11),PSWgt3)

      m_Top = m_SMTop
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd( 9:12),MomExt(1:4,7:9),PSWgt4)! bot lep neu
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(13:16),MomExt(1:4,12:14),PSWgt5)
      m_Top = m_HTop
      PSWgt = PSWgt * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
      nHel(1:2) = getHelicity(yrnd(17))
   ELSEIF(XTOPDECAYS.EQ.0) THEN
      nHel(1:2) = getHelicity(yrnd(5))
   ENDIF

   call Kinematics_TTbarETmiss(.false.,MomExt,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_1L_HtHtbqqb = 0d0
      return
   endif

   call setPDFs(eta1,eta2,MuFac,pdf)
   PDFFac_a = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
            + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
            + pdf(Bot_,1)*pdf(ABot_,2)
   PDFFac_b = pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
            + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
            + pdf(Bot_,2)*pdf(ABot_,1)

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   RunFactor = RunAlphaS(NLOParam,MuRen)
   PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))

   LO_Res_Unpol             = (0d0,0d0)
   NLO_Res_Unpol(-2:1)      = (0d0,0d0)
   NLO_Res_Unpol_Ferm(-2:1) = (0d0,0d0)
!------------ LO --------------
IF( CORRECTION.EQ.0 ) THEN
  do npdf=1,2
    if(npdf.eq.1) then
        PDFFac = PDFFac_a
    elseif(npdf.eq.2) then
        PDFFac = PDFFac_b
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
    endif
    ISFac = MomCrossing(MomExt)
    call InitCurrCache()
    call SetPropagators()

    do iHel=nHel(1),nHel(2)
    do BH1Hel=-1,BHMaxHel
    do BH2Hel=-1,BHMaxHel
        if( Helicities(iHel,3).eq.Helicities(iHel,4) ) cycle ! remove helicities where IS quarks have ++ / -- helicities because tree vanishes
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
        IF(XTOPDECAYS.eq.1) THEN
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExt(1:4,5:9))
          call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExt(1:4,10:14))
        ELSEIF(XTOPDECAYS.EQ.2) THEN
          call HTopA0Decay(ExtParticle(1),DKX_HTA0_LO,MomExt(1:4,5:9))
          call HTopA0Decay(ExtParticle(2),DKX_HTA0_LO,MomExt(1:4,10:14))
        ENDIF
        call SetPolarizations()
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
        enddo
        LO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * ParityFlip*BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
        enddo
        enddo
        LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol*PDFFac
    enddo!helicity loop
    enddo!helicity loop
    enddo!helicity loop
  enddo! npdf loop

!------------ 1 LOOP --------------
ELSEIF( CORRECTION.EQ.1 ) THEN
  do npdf=1,2
!   do npdf=2,2; print *, "npdf=2"
    if(npdf.eq.1) then
        PDFFac = PDFFac_a
    elseif(npdf.eq.2) then
        PDFFac = PDFFac_b
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
    endif
    ISFac = MomCrossing(MomExt)
    call InitCurrCache()
    call SetPropagators()
    do iHel=nHel(1),nHel(2)
    do BH1Hel=-1,BHMaxHel
    do BH2Hel=-1,BHMaxHel
      if( Helicities(iHel,3).eq.Helicities(iHel,4) ) cycle ! remove helicities where IS quarks have ++ / -- helicities because tree vanishes
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      IF(XTOPDECAYS.eq.1) THEN
        call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExt(1:4,5:9))
        call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExt(1:4,10:14))
      ELSEIF(XTOPDECAYS.EQ.2) THEN
        call HTopA0Decay(ExtParticle(1),DKX_HTA0_LO,MomExt(1:4,5:9))
        call HTopA0Decay(ExtParticle(2),DKX_HTA0_LO,MomExt(1:4,10:14))
      ENDIF
      call SetPolarizations()
      do iPrimAmp=1,7
          call EvalTree(BornAmps(iPrimAmp))
      enddo

      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * ParityFlip*BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol*PDFFac


! ------------ bosonic loops --------------
      do iPrimAmp=1,4
          call SetKirill(PrimAmps(iPrimAmp))
          call PentCut(PrimAmps(iPrimAmp))
          call QuadCut(PrimAmps(iPrimAmp))
          call TripCut(PrimAmps(iPrimAmp))
          call DoubCut(PrimAmps(iPrimAmp))
          call SingCut(PrimAmps(iPrimAmp))
          call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)

          call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen**2)
          PrimAmps(iPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(iPrimAmp)%Result(-2:1)
!           call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,2,rdiv(2),rdiv(1))
!           call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
!           pause
      enddo
      BosonicPartAmp(1,-2:1) =  &
                             + Nc * PrimAmps(PrimAmp1_1234)%Result(-2:1) &
                             - 2d0/Nc*( PrimAmps(PrimAmp1_1234)%Result(-2:1) + PrimAmps(PrimAmp1_1243)%Result(-2:1) ) &
                             - 1d0/Nc * ( PrimAmps(PrimAmp3_1432)%Result(-2:1) + PrimAmps(PrimAmp4_1234)%Result(-2:1) )

      BosonicPartAmp(2,-2:1) = PrimAmps(PrimAmp1_1243)%Result(-2:1)  &
                             + 1d0/Nc**2 * ( PrimAmps(PrimAmp1_1234)%Result(-2:1) + PrimAmps(PrimAmp1_1243)%Result(-2:1) )  &
                             + 1d0/Nc**2 * ( PrimAmps(PrimAmp3_1432)%Result(-2:1) + PrimAmps(PrimAmp4_1234)%Result(-2:1) )


      NLO_Res_Pol(-2:1) = (0d0,0d0)
      do jPrimAmp=1,2
      do iPrimAmp=1,NumBornAmps
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + Col1L_ttbqqb(iPrimAmp,jPrimAmp) * ParityFlip*dreal( BornAmps(iPrimAmp)%Result * dconjg(BosonicPartAmp(jPrimAmp,-2:1)) )
      enddo
      enddo
      NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1) + NLO_Res_Pol(-2:1)*PDFFac

! ------------ fermionic loops --------------
      do iPrimAmp=5,7
          call SetKirill(PrimAmps(iPrimAmp))
          call PentCut(PrimAmps(iPrimAmp))
          call QuadCut(PrimAmps(iPrimAmp))
          call TripCut(PrimAmps(iPrimAmp))
          call DoubCut(PrimAmps(iPrimAmp))
          call SingCut(PrimAmps(iPrimAmp))
          call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
          PrimAmps(iPrimAmp)%Result(-2:1) = -(0d0,1d0)*PrimAmps(iPrimAmp)%Result(-2:1) !minus if from closed fermion loop
!           call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,2,rdiv(2),rdiv(1))
!           call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
      enddo

      FermionPartAmp(1,-2:1) =          ( Nf_light*PrimAmps(5)%Result(-2:1) + PrimAmps(6)%Result(-2:1) + PrimAmps(7)%Result(-2:1)*RemoveHTopClosedLoop )
      FermionPartAmp(2,-2:1) = -1d0/Nc *( Nf_light*PrimAmps(5)%Result(-2:1) + PrimAmps(6)%Result(-2:1) + PrimAmps(7)%Result(-2:1)*RemoveHTopClosedLoop )

      NLO_Res_Pol(-2:1) = (0d0,0d0)
      do jPrimAmp=1,2
      do iPrimAmp=1,NumBornAmps
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + Col1L_ttbqqb(iPrimAmp,jPrimAmp) * ParityFlip*dreal( BornAmps(iPrimAmp)%Result * dconjg(FermionPartAmp(jPrimAmp,-2:1)) )
      enddo
      enddo
      NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) + NLO_Res_Pol(-2:1)*PDFFac

   enddo!helicity loop
   enddo!helicity loop
   enddo!helicity loop
  enddo! npdf loop
  call swapMom(MomExt(1:4,1),MomExt(1:4,2))   ! swap back to original order, for ID below
! print *, "mom swp disabled"
ENDIF



IF( CORRECTION.EQ.0 ) THEN
!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * WidthExpansion
   EvalCS_1L_HtHtbqqb = LO_Res_Unpol * PreFac

ELSEIF( CORRECTION.EQ.1 ) THEN
!  overall normalization: (4*Pi)^eps/Gamma(1-eps)
!  CT contributions                           ! beta        !Htop WFRC
   NLO_Res_UnPol(-1) = NLO_Res_UnPol(-1) + (-11d0/3d0*3d0 - 3d0*4d0/3d0)*LO_Res_Unpol
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-3d0*4d0/3d0)*2d0*dlog(MuRen/m_Htop)*LO_Res_Unpol  ! finite log(mu2) contrib. from  Htop WFRC
   NLO_Res_UnPol_Ferm(-1) = NLO_Res_UnPol_Ferm(-1) + (2d0/3d0*Nf_light+2d0/3d0*Nf_heavy)*LO_Res_Unpol
   NLO_Res_UnPol_Ferm(-1) = NLO_Res_UnPol_Ferm(-1) + (2d0/3d0*RemoveHTopClosedLoop     )*LO_Res_Unpol! this is from the closed HTop loop
   NLO_Res_UnPol_Ferm( 0) = NLO_Res_UnPol_Ferm( 0) + (2d0/3d0*Nf_heavy)*2d0*dlog(MuRen/m_SMtop)*LO_Res_Unpol  ! finite log(mu2) contrib. from heavy flavor in alpha_s ren.
   NLO_Res_UnPol_Ferm( 0) = NLO_Res_UnPol_Ferm( 0) + (2d0/3d0*RemoveHTopClosedLoop)*2d0*dlog(MuRen/m_Htop)*LO_Res_Unpol   ! finite log(mu2) contrib. from HTop in alpha_s ren.
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-5d0/2d0*8d0/3d0 )*LO_Res_Unpol   ! finite contribution from top WFRC's
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + LO_Res_Unpol ! shift alpha_s^DR --> alpha_s^MSbar

!  factor out (Mu2/mTop**2)^eps
!    NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + NLO_Res_UnPol(-1)*2d0*dlog(m_top/MuRen) + NLO_Res_UnPol(-2)*dlog(m_top/MuRen)**2
!    NLO_Res_UnPol(-1) = NLO_Res_UnPol(-1) + NLO_Res_UnPol(-2)*2d0*dlog(m_top/MuRen)
!    NLO_Res_UnPol_Ferm(0) = NLO_Res_UnPol_Ferm(0) + NLO_Res_UnPol_Ferm(-1)*2d0*dlog(m_top/MuRen)

!  normalization
   LO_Res_Unpol = LO_Res_Unpol                         * ISFac * (alpha_s4Pi*RunFactor)**2
   NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1)           * ISFac * (alpha_s4Pi*RunFactor)**2 * alpha_sOver2Pi*RunFactor
   NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) * ISFac * (alpha_s4Pi*RunFactor)**2 * alpha_sOver2Pi*RunFactor

   EvalCS_1L_HtHtbqqb = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac

ELSEIF( CORRECTION.EQ.3 ) THEN


   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   ISFac = MomCrossing(MomExt)
   if( XTopDecays.eq.0 ) then
      xE = yRnd(5)
   else
      xE = yRnd(17)
   endif


IF( PROCESS.EQ.46 ) THEN
! xE=0.3d0
   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
   call EvalIntDipoles_QQBHTHTBG((/MomExt(1:4,3),MomExt(1:4,4),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,5:14),xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
   EvalCS_1L_HtHtbqqb = HOp(1)    * (pdf(Up_,1)*pdf(AUp_,2)+pdf(Dn_,1)*pdf(ADn_,2)+pdf(Chm_,1)*pdf(AChm_,2)+pdf(Str_,1)*pdf(AStr_,2)+pdf(Bot_,1)*pdf(ABot_,2) ) &
                      + HOp(2)/xE * (pdf_z(Up_,1)*pdf(AUp_,2)+pdf_z(Dn_,1)*pdf(ADn_,2)+pdf_z(Chm_,1)*pdf(AChm_,2)+pdf_z(Str_,1)*pdf(AStr_,2)+pdf_z(Bot_,1)*pdf(ABot_,2) ) &
                      + HOp(3)/xE * (pdf(Up_,1)*pdf_z(AUp_,2)+pdf(Dn_,1)*pdf_z(ADn_,2)+pdf(Chm_,1)*pdf_z(AChm_,2)+pdf(Str_,1)*pdf_z(AStr_,2)+pdf(Bot_,1)*pdf_z(ABot_,2) )


    call EvalIntDipoles_QQBHTHTBG((/MomExt(1:4,3),MomExt(1:4,4),-MomExt(1:4,2),-MomExt(1:4,1)/),MomExt(1:4,5:14),xE,HOp(1:3))
    HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
    EvalCS_1L_HtHtbqqb = EvalCS_1L_HtHtbqqb  &
                     + HOp(1)    * (pdf(Up_,2)*pdf(AUp_,1)+pdf(Dn_,2)*pdf(ADn_,1)+pdf(Chm_,2)*pdf(AChm_,1)+pdf(Str_,2)*pdf(AStr_,1)+pdf(Bot_,2)*pdf(ABot_,1) ) &
                     + HOp(2)/xE * (pdf_z(Up_,2)*pdf(AUp_,1)+pdf_z(Dn_,2)*pdf(ADn_,1)+pdf_z(Chm_,2)*pdf(AChm_,1)+pdf_z(Str_,2)*pdf(AStr_,1)+pdf_z(Bot_,2)*pdf(ABot_,1) ) &
                     + HOp(3)/xE * (pdf(Up_,2)*pdf_z(AUp_,1)+pdf(Dn_,2)*pdf_z(ADn_,1)+pdf(Chm_,2)*pdf_z(AChm_,1)+pdf(Str_,2)*pdf_z(AStr_,1)+pdf(Bot_,2)*pdf_z(ABot_,1) )


! checked for    XTopDK=1 TopDK=1 and two mu-scales and for XTopDK=TopDK=0
! print *, "1L check",NLO_Res_UnPol(-2)                            * PreFac/(alpha_sOver2Pi*RunFactor)/(LO_Res_Unpol*PreFac)
! print *, "1L check",( NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1) ) * PreFac/(alpha_sOver2Pi*RunFactor)/(LO_Res_Unpol*PreFac)
! print *, "res1=", (EvalCS_1L_HtHtbqqb)/(alpha_sOver2Pi*RunFactor)/(LO_Res_Unpol*PreFac)
! pause


ELSEIF( PROCESS.EQ.43 ) THEN
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

! checked for    XTopDK=1 TopDK=1 and two mu-scales and for XTopDK=TopDK=0
      call EvalIntDipoles_QGHTHTBQ((/MomExt(1:4,3),MomExt(1:4,4),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,5:14),xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_1L_HtHtbqqb = HOp(1)    * (pdf(Up_,1)*pdf(0,2)+pdf(Dn_,1)*pdf(0,2)+pdf(Chm_,1)*pdf(0,2)+pdf(Str_,1)*pdf(0,2)+pdf(Bot_,1)*pdf(0,2) ) &
                         + HOp(2)/xE * (pdf_z(Up_,1)*pdf(0,2)+pdf_z(Dn_,1)*pdf(0,2)+pdf_z(Chm_,1)*pdf(0,2)+pdf_z(Str_,1)*pdf(0,2)+pdf_z(Bot_,1)*pdf(0,2) ) &
                         + HOp(3)/xE * (pdf(Up_,1)*pdf_z(0,2)+pdf(Dn_,1)*pdf_z(0,2)+pdf(Chm_,1)*pdf_z(0,2)+pdf(Str_,1)*pdf_z(0,2)+pdf(Bot_,1)*pdf_z(0,2) )
 

      call EvalIntDipoles_QGHTHTBQ((/MomExt(1:4,3),MomExt(1:4,4),-MomExt(1:4,2),-MomExt(1:4,1)/),MomExt(1:4,5:14),xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_1L_HtHtbqqb = EvalCS_1L_HtHtbqqb  &
                      + HOp(1)    * (pdf(Up_,2)*pdf(0,1)+pdf(Dn_,2)*pdf(0,1)+pdf(Chm_,2)*pdf(0,1)+pdf(Str_,2)*pdf(0,1)+pdf(Bot_,2)*pdf(0,1) ) &
                      + HOp(2)/xE * (pdf_z(Up_,2)*pdf(0,1)+pdf_z(Dn_,2)*pdf(0,1)+pdf_z(Chm_,2)*pdf(0,1)+pdf_z(Str_,2)*pdf(0,1)+pdf_z(Bot_,2)*pdf(0,1) ) &
                      + HOp(3)/xE * (pdf(Up_,2)*pdf_z(0,1)+pdf(Dn_,2)*pdf_z(0,1)+pdf(Chm_,2)*pdf_z(0,1)+pdf(Str_,2)*pdf_z(0,1)+pdf(Bot_,2)*pdf_z(0,1) )


ELSEIF( PROCESS.EQ.44 ) THEN
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

! checked for    XTopDK=1 TopDK=1 and two mu-scales and for XTopDK=TopDK=0
      call EvalIntDipoles_QGHTHTBQ((/MomExt(1:4,3),MomExt(1:4,4),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,5:14),xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_1L_HtHtbqqb = HOp(1)    * (pdf(AUp_,1)*pdf(0,2)+pdf(ADn_,1)*pdf(0,2)+pdf(AChm_,1)*pdf(0,2)+pdf(AStr_,1)*pdf(0,2)+pdf(ABot_,1)*pdf(0,2) ) &
                         + HOp(2)/xE * (pdf_z(AUp_,1)*pdf(0,2)+pdf_z(ADn_,1)*pdf(0,2)+pdf_z(AChm_,1)*pdf(0,2)+pdf_z(AStr_,1)*pdf(0,2)+pdf_z(ABot_,1)*pdf(0,2) ) &
                         + HOp(3)/xE * (pdf(AUp_,1)*pdf_z(0,2)+pdf(ADn_,1)*pdf_z(0,2)+pdf(AChm_,1)*pdf_z(0,2)+pdf(AStr_,1)*pdf_z(0,2)+pdf(ABot_,1)*pdf_z(0,2) )
 

      call EvalIntDipoles_QGHTHTBQ((/MomExt(1:4,3),MomExt(1:4,4),-MomExt(1:4,2),-MomExt(1:4,1)/),MomExt(1:4,5:14),xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_1L_HtHtbqqb = EvalCS_1L_HtHtbqqb &
                        + HOp(1)    * (pdf(AUp_,2)*pdf(0,1)+pdf(ADn_,2)*pdf(0,1)+pdf(AChm_,2)*pdf(0,1)+pdf(AStr_,2)*pdf(0,1)+pdf(ABot_,2)*pdf(0,1) ) &
                        + HOp(2)/xE * (pdf_z(AUp_,2)*pdf(0,1)+pdf_z(ADn_,2)*pdf(0,1)+pdf_z(AChm_,2)*pdf(0,1)+pdf_z(AStr_,2)*pdf(0,1)+pdf_z(ABot_,2)*pdf(0,1) ) &
                        + HOp(3)/xE * (pdf(AUp_,2)*pdf_z(0,1)+pdf(ADn_,2)*pdf_z(0,1)+pdf(AChm_,2)*pdf_z(0,1)+pdf(AStr_,2)*pdf_z(0,1)+pdf(ABot_,2)*pdf_z(0,1) )
ENDIF
ENDIF


   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_HtHtbqqb)
   enddo
   EvalCounter = EvalCounter + 1


   EvalCS_1L_HtHtbqqb = EvalCS_1L_HtHtbqqb/VgsWgt

return
END FUNCTION






FUNCTION EvalCS_Real_HtHtbggg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModDipoles_GGHtHtbG
implicit none
real(8) :: EvalCS_Real_HtHtbggg,DipoleResult,yRnd(1:VegasMxDim),VgsWgt
complex(8) :: LO_Res_Pol,LO_Res_Unpol
integer :: iHel,jHel,kHel,iPrimAmp,jPrimAmp
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5,ISFac
real(8) :: MomExt(1:4,1:15),MomBoost(1:4)
logical :: applyPSCut
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac
real(8) :: pdf(-6:6,1:2),s12,s13
integer :: NBin(1:NumMaxHisto),NHisto,nHel(1:2),BHMaxHel,BH1Hel,BH2Hel
!include 'misc/global_import'
include 'vegas_common.f'


   EvalCS_Real_HtHtbggg = 0d0
   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if( EHat.le.2d0*m_HTop ) then
      EvalCS_Real_HtHtbggg = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to3HT(EHat,yRnd(3:7),MomExt(1:4,1:5),PSWgt)! Glu HTbar, HT
   if( PSWgt.eq.0d0 ) then
      EvalCS_Real_HtHtbggg = 0d0
      SkipCounter = SkipCounter + 1
      return
   endif
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))
   ISFac = MomCrossing(MomExt)

   BHMaxHel = -1
   IF(XTOPDECAYS.EQ.1) THEN
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,4),yRnd(8:9),MomExt(1:4,6:7),PSWgt2)!  BH top
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,5),yRnd(10:11),MomExt(1:4,11:12),PSWgt3)

      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,7), yRnd(12:15),MomExt(1:4,8:10),PSWgt4)! bot lep neu
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,12),yRnd(16:19),MomExt(1:4,13:15),PSWgt5)

      PSWgt = PSWgt * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
      BHMaxHel = +1
   ELSEIF(XTOPDECAYS.EQ.2) THEN
      call Error("XTOPDECAYS.EQ.2 is not yet supported")
   ENDIF

   call SetPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   RunFactor = RunAlphaS(2,MuRen)
   call Kinematics_TTbarETmiss(.true.,MomExt,(/4,5,6,11,7,12,8,9,10,13,14,15,3/),applyPSCut,NBin)

! applyPScut = .true.! this is for inverted check

if(  applyPSCut ) then
   EvalCS_Real_HtHtbggg = 0d0
else
   call InitCurrCache()
   call SetPropagators()
   PreFac = PreFac
   LO_Res_Unpol = (0d0,0d0)
   do iHel=1,NumHelicities
   do BH1Hel=-1,BHMaxHel
   do BH2Hel=-1,BHMaxHel
      IF( XTOPDECAYS.EQ.1 ) THEN
         call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExt(1:4,6:10))
         call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExt(1:4,11:15))
      ENDIF
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      m_Top=m_HTop
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo
      m_Top=m_SMTop
      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbggg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
   enddo!helicity loop
   enddo!helicity loop
   enddo!helicity loop

!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**3
   EvalCS_Real_HtHtbggg = LO_Res_Unpol * PreFac

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_Real_HtHtbggg)
   enddo
   EvalCounter = EvalCounter + 1

! !     gg->ttbg
!       MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!       MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!       MG_MOM(0:3,3) = MomExt(1:4,5)*100d0
!       MG_MOM(0:3,4) = MomExt(1:4,4)*100d0
!       MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!       call coupsm(0)
!       call SGG_TTBG(MG_MOM,MadGraph_tree)
!       print *, ""
!       print *, "Unpolarized LO result:"
!       print *, "My tree:          ", LO_Res_Unpol
!       print *, "MadGraph hel.amp:", MadGraph_tree*(100d0)**2
!       print *, "ratio: ", MadGraph_tree*(100d0)**2/dble((LO_Res_Unpol))
!       pause
endif! applyPSCut


!    cancellation checked for XTopDK=0 and (XTopDK=1,TopDK=4,1)
     call EvalDipoles_GGHtHtbG((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_HTop**2,m_HTop**2,0d0/),yRnd(8:19),PreFac,DipoleResult)

!   s12 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,2))
!   s13 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,3))
!   write(* ,"(1PE23.16,3X,1PE23.16,3X,1PE23.16,3X,1PE23.16)") s13/s12,EvalCS_Real_HtHtbggg,DipoleResult, (EvalCS_Real_HtHtbggg/(-DipoleResult)-1d0)
!   pause


   if( IsNan(EvalCS_Real_HtHtbggg) .or. IsNan(DipoleResult) ) then
        print *, "NAN:",EvalCS_Real_HtHtbggg,DipoleResult
        print *, yRnd(:)
        print *, PSWgt , VgsWgt , PDFFac, sHatJacobi
        print *, eta1,eta2,MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,1)
        print *, MomExt(1:4,2)
        print *, MomExt(1:4,3)
        print *, MomExt(1:4,4)
        print *, MomExt(1:4,5)
!         stop
   endif

        EvalCS_Real_HtHtbggg = (EvalCS_Real_HtHtbggg + DipoleResult)/VgsWgt

return
END FUNCTION









FUNCTION EvalCS_Real_HtHtbqqbg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModDipoles_QQBHtHtbG
use ModDipoles_QGHtHtbQ
implicit none
real(8) ::  EvalCS_Real_HtHtbqqbg,EvalCS_Dips_ttbqqbgp,DipoleResult,yRnd(1:VegasMxDim),VgsWgt
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol
integer :: iHel,jHel,kHel,iPrimAmp,jPrimAmp
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5,ISFac,PreFacDip
real(8) :: MomExt(1:4,1:15)
logical :: applyPSCut
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac_a,PDFFac_b,PDFFac
real(8) :: pdf(-6:6,1:2),s12,s13
integer :: NBin(1:NumMaxHisto),NHisto,nHel(1:2),npdf,BHMaxHel,BH1Hel,BH2Hel
real(8),parameter :: Nc=3d0
include 'vegas_common.f'



   EvalCS_Real_HtHtbqqbg = 0d0
   EvalCS_Dips_ttbqqbgp = 0d0
   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if( EHat.le.2d0*m_HTop ) then
      EvalCS_Real_HtHtbqqbg = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to3HT(EHat,yRnd(3:7),MomExt(1:4,1:5),PSWgt)! Glu HTbar, HT
   if( PSWgt.eq.0d0 ) then
      EvalCS_Real_HtHtbqqbg = 0d0
      SkipCounter = SkipCounter + 1
      return
   endif
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))
   ISFac = MomCrossing(MomExt)
   BHMaxHel = -1
   IF(XTOPDECAYS.EQ.1) THEN
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,4),yRnd(8:9),MomExt(1:4,6:7),PSWgt2)!  BH top
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,5),yRnd(10:11),MomExt(1:4,11:12),PSWgt3)
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,7), yRnd(12:15),MomExt(1:4,8:10),PSWgt4)! bot lep neu
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,12),yRnd(16:19),MomExt(1:4,13:15),PSWgt5)
      PSWgt = PSWgt * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
      BHMaxHel = +1
   ELSEIF(XTOPDECAYS.EQ.2) THEN
      call Error("XTOPDECAYS.EQ.2 is not yet supported")
   ENDIF


   call SetPDFs(eta1,eta2,MuFac,pdf)
   IF( PROCESS.EQ.46 ) THEN
          PDFFac_a = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
                   + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
                   + pdf(Bot_,1)*pdf(ABot_,2)
          PDFFac_b = pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
                   + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
                   + pdf(Bot_,2)*pdf(ABot_,1)
    ELSEIF( PROCESS.EQ.43 ) THEN
          PDFFac_a = pdf(Up_,1) *pdf(0,2) + pdf(Dn_,1) *pdf(0,2)  &
                   + pdf(Chm_,1)*pdf(0,2) + pdf(Str_,1)*pdf(0,2)  &
                   + pdf(Bot_,1)*pdf(0,2)
          PDFFac_b = pdf(Up_,2) *pdf(0,1) + pdf(Dn_,2) *pdf(0,1)  &
                   + pdf(Chm_,2)*pdf(0,1) + pdf(Str_,2)*pdf(0,1)  &
                   + pdf(Bot_,2)*pdf(0,1)
    ELSEIF( PROCESS.EQ.44 ) THEN
          PDFFac_a = pdf(AUp_,2) *pdf(0,1) + pdf(ADn_,2) *pdf(0,1)  &
                   + pdf(AChm_,2)*pdf(0,1) + pdf(AStr_,2)*pdf(0,1)  &
                   + pdf(ABot_,2)*pdf(0,1)
          PDFFac_b = pdf(AUp_,1) *pdf(0,2) + pdf(ADn_,1) *pdf(0,2)  &
                   + pdf(AChm_,1)*pdf(0,2) + pdf(AStr_,1)*pdf(0,2)  &
                   + pdf(ABot_,1)*pdf(0,2)
    ENDIF


DO NPDF=1,2   !1; print *, "fixed npdf"
    call Kinematics_TTbarETmiss(.true.,MomExt,(/4,5,6,11,7,12,8,9,10,13,14,15,3/),applyPSCut,NBin)

! applyPScut = .true.! this is for inverted check

    if(npdf.eq.1) then
        PDFFac = PDFFac_a
    elseif(npdf.eq.2) then
        PDFFac = PDFFac_b
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
    endif
    PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
    RunFactor = RunAlphaS(2,MuRen)
    ISFac = MomCrossing(MomExt)
    call InitCurrCache()
    call SetPropagators()

if(  applyPSCut ) then
    EvalCS_Real_HtHtbqqbg = 0d0
else

!------------ LO ----------------
   LO_Res_Unpol = (0d0,0d0)
    do iHel=1,NumHelicities
    do BH1Hel=-1,BHMaxHel
    do BH2Hel=-1,BHMaxHel
        IF(XTOPDECAYS.eq.1) THEN
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExt(1:4,6:10))
          call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExt(1:4,11:15))
        ELSEIF(XTOPDECAYS.EQ.2) THEN
          call HTopA0Decay(ExtParticle(1),DKX_HTA0_LO,MomExt(1:4,6:10))
          call HTopA0Decay(ExtParticle(2),DKX_HTA0_LO,MomExt(1:4,11:15))
        ENDIF
      call HelCrossing(Helicities(iHel,1:5))
      call SetPolarizations()
      
      m_Top=m_HTop
       do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo
      m_Top=m_SMTop

      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqbg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
   enddo!helicity loop
   enddo!helicity loop
   enddo!helicity loop

!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**3 * PreFac *  PDFFac 
   EvalCS_Real_HtHtbqqbg = EvalCS_Real_HtHtbqqbg + LO_Res_Unpol

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol))
   enddo
   EvalCounter = EvalCounter + 1
endif! applyPSCut



  PreFacDip = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
  IF( PROCESS.EQ.46 ) THEN
      call EvalDipoles_QQBHtHtbG((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_HTop**2,m_HTop**2,0d0/),yRnd(8:19),PreFacDip,DipoleResult)
  ELSEIF( PROCESS.EQ.43 ) THEN
      call EvalDipoles_QGHtHtbQ((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_HTop**2,m_HTop**2,0d0/),yRnd(8:19),PreFacDip,DipoleResult)
  ELSEIF( PROCESS.EQ.44 ) THEN
      call EvalDipoles_QGHtHtbQ((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_HTop**2,m_HTop**2,0d0/),yRnd(8:19),PreFacDip,DipoleResult)
  ENDIF
  EvalCS_Dips_ttbqqbgp = EvalCS_Dips_ttbqqbgp + DipoleResult

ENDDO! npdf loop
! call swapMom(MomExt(1:4,1),MomExt(1:4,2))   ! swap back, not necessary here


!    cancellations checked for XTopDK=0 and (XTopDK=1,TopDK=4,1)
!  s12 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,2))
!  s13 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,3))
!  write(* ,"(1PE23.16,3X,1PE23.16,3X,1PE23.16,3X,1PE23.16)") s13/s12,EvalCS_Real_HtHtbqqbg,EvalCS_Dips_ttbqqbgp, (EvalCS_Real_HtHtbqqbg/(-EvalCS_Dips_ttbqqbgp)- 1d0)
!  pause

   if( IsNan(EvalCS_Real_HtHtbqqbg) .or. IsNan(EvalCS_Dips_ttbqqbgp) ) then
        print *, "NAN:",EvalCS_Real_HtHtbqqbg,EvalCS_Dips_ttbqqbgp
        print *, yRnd(:)
        print *, PSWgt , VgsWgt , PDFFac, sHatJacobi
        print *, eta1,eta2,MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,1)
        print *, MomExt(1:4,2)
        print *, MomExt(1:4,3)
        print *, MomExt(1:4,4)
        print *, MomExt(1:4,5)
   endif

   EvalCS_Real_HtHtbqqbg = (EvalCS_Real_HtHtbqqbg + EvalCS_Dips_ttbqqbgp)/VgsWgt

END FUNCTION











FUNCTION EvalCS_1L_ststbgg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModIntDipoles_GGSTSTBG
use ModTopdecay
implicit none
real(8) ::  EvalCS_1L_ststbgg,yRnd(1:VegasMxDim),VgsWgt,HOp(1:3)
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1),FermionLoopPartAmp(7:8,-2:1),CataniRes(-2:-1)
complex(8) :: LO_Decay_UnPol,Spi(1:4),BarSpi(1:4),M_T_BWENU,Msq_T_BWENU
integer :: iHel,jHel,kHel,lHel,iPrimAmp,jPrimAmp
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5,ISFac
real(8) :: MomExt(1:4,1:15),MomP(1:4,1:4),MomBoost(1:4)
logical :: applyPSCut
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,AccPoles
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2),xE,sigmaTot,beta
integer :: NBin(1:NumMaxHisto),NHisto,nHel(1:2)
type(Particle) :: TopQuark(1:2)
!include 'misc/global_import'
include 'vegas_common.f'

! yrnd(1)=0.8d0; yrnd(2)=0.65d0

   EvalCS_1L_ststbgg = 0d0
   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if( EHat.le.2d0*m_STop ) then
      EvalCS_1L_ststbgg = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to2Stops(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt)! AStop, Stop
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
! !!!!!!!!!!!!!!!!!!!!! this is just for checking gauge invariance with gluon off the beam pipe
! print *, "Test Boost activated"
! MomBoost(1:4) = 1.8d0*MomExt(1:4,3)
! call boost(MomExt(1:4,1),MomBoost(1:4),1.8d0*m_sTop)
! call boost(MomExt(1:4,2),MomBoost(1:4),1.8d0*m_sTop)
! call boost(MomExt(1:4,3),MomBoost(1:4),1.8d0*m_sTop)
! call boost(MomExt(1:4,4),MomBoost(1:4),1.8d0*m_sTop)
! !!!!!!!!!!!!!!!!!!!!
   ISFac = MomCrossing(MomExt)

   IF(XTOPDECAYS.EQ.3) THEN
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  Chi top
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,4),yRnd(7:8),MomExt(1:4,10:11),PSWgt3)
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd( 9:12),MomExt(1:4,7:9),PSWgt4)! bot lep neu
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(13:16),MomExt(1:4,12:14),PSWgt5)
      PSWgt = PSWgt * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
      nHel(1:2) = getHelicity(yrnd(17))
   ELSE
      nHel(1:2) = getHelicity(yrnd(5))
   ENDIF

   call Kinematics_TTbarETmiss(.false.,MomExt,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_1L_ststbgg = 0d0
      return
   endif
   call InitCurrCache()
   call SetPropagators()

   call SetPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   RunFactor = RunAlphaS(NLOParam,MuRen)
   PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))


!!!!! for of total cross section (fix y1,y2!)
!   beta = dsqrt(1d0-4d0*m_stop**2/EHat**2)
!   SigmaTot = (alpha_s*RunFactor)**2*DblPi/EHat**2*( beta*(5d0/48d0+31d0*m_stop**2/24d0/EHat**2) + &
!                                                     dlog((1d0-beta)/(1d0+beta))*(2d0*m_stop**2/3d0/EHat**2+m_stop**4/6d0/Ehat**4)  )
!   SigmaTot = SigmaTot * fbGeV2*SHatJacobi*PDFFac 
!   PreFac = PreFac/SigmaTot
!!!!!!
   LO_Res_Unpol             = (0d0,0d0)
   NLO_Res_Unpol(-2:1)      = (0d0,0d0)
   NLO_Res_Unpol_Ferm(-2:1) = (0d0,0d0)
!------------ LO ----------------
IF( Correction.EQ.0 ) THEN


!    do iHel=nHel(1),nHel(2)
!       call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExt(1:4,5:9))
!       call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,10:14))
!       call HelCrossing(Helicities(iHel,1:4))
!       call SetPolarizations()
!       do iPrimAmp=1,NumBornAmps
!           call EvalTree(BornAmps(iPrimAmp))
!       enddo
! 
!       LO_Res_Pol = (0d0,0d0)
!       do jPrimAmp=1,NumBornAmps
!       do iPrimAmp=1,NumBornAmps
!           LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
!       enddo
!       enddo
!       LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
!    enddo!helicity loop





! ! !  this LO code block is equivalent to the above but it is 4 times faster
   LO_Res_Unpol = (0d0,0d0)
   if( XTopDecays.eq.0 ) then
      jHel=1
   else
      jHel=4
   endif

   do iHel=nHel(1),nHel(2),jHel  ! ,4 because we jump over stop polarizations which are summed over outside
      ExtParticle(1)%Pol(1)=(1d0,0d0)
      ExtParticle(2)%Pol(1)=(1d0,0d0)
      call HelCrossing(Helicities(iHel,1:4))
      call SetPolarizations()
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo

      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
   enddo!helicity loop


!  spin correlated top quarks
   if( XTopDecays.ne.0) then
      LO_Decay_UnPol = (0d0,0d0)
      do jHel=-1,1,2
      do kHel=-1,1,2
          call STopDecay(ExtParticle(1),DKX_STChi0_LO,jHel,MomExt(1:4,5:9))
          call STopDecay(ExtParticle(2),DKX_STChi0_LO,kHel,MomExt(1:4,10:14))
          LO_Decay_UnPol = LO_Decay_UnPol + ExtParticle(1)%Pol(1)*dconjg(ExtParticle(1)%Pol(1)) * ExtParticle(2)%Pol(1)*dconjg(ExtParticle(2)%Pol(1))
      enddo
      enddo

      LO_Res_UnPol = LO_Res_UnPol * LO_Decay_UnPol
   endif

! ! ! !  spin uncorrelated top quarks
!    LO_Decay_UnPol = (0d0,0d0)
!    do iHel=-1,1,2
!    do jHel=-1,1,2
!    do kHel=-1,1,2
!    do lHel=-1,1,2
!       call STopDecay(ExtParticle(1),DKX_STChi0_LO,jHel,MomExt(1:4,5:9),HelTop=iHel)
!       call STopDecay(ExtParticle(2),DKX_STChi0_LO,kHel,MomExt(1:4,10:14),HelTop=lHel)
!       LO_Decay_UnPol = LO_Decay_UnPol + ExtParticle(1)%Pol(1)*dconjg(ExtParticle(1)%Pol(1)) * ExtParticle(2)%Pol(1)*dconjg(ExtParticle(2)%Pol(1))
!    enddo
!    enddo
!    enddo
!    enddo
! 
! ! spherical decay of ATop
!       TopQuark(1)%Mom(1:4) = MomExt(1:4,6)
!       TopQuark(1)%PartType = ATop_
!       TopQuark(1)%Mass = m_Top
!       TopQuark(1)%Mass2 = m_Top**2
! 
!       TopQuark(2)%Mom(1:4) = MomExt(1:4,11)
!       TopQuark(2)%PartType = Top_
!       TopQuark(2)%Mass = m_Top
!       TopQuark(2)%Mass2 = m_Top**2
! 
!       call TopDecay(TopQuark(1),DK_LO,MomExt(1:4,7:9))
!       Spi(1:4) = TopQuark(1)%Pol(1:4)
!       call vBarSpi(TopQuark(1)%Mom(1:4),m_Top,-1,BarSpi(1:4))
!       M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
!       Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0
! 
!       call vBarSpi(TopQuark(1)%Mom(1:4),m_Top,+1,BarSpi(1:4))
!       M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
!       Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
! 
!       LO_Decay_UnPol = LO_Decay_UnPol * Msq_T_BWENU
! 
! ! spherical decay of Top
!       call TopDecay(TopQuark(2),DK_LO,MomExt(1:4,12:14))
!       BarSpi(1:4) = TopQuark(2)%Pol(1:4)
!       call uSpi(TopQuark(2)%Mom(1:4),m_Top,-1,Spi(1:4))
!       M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
!       Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0
! 
!       call uSpi(TopQuark(2)%Mom(1:4),m_Top,+1,Spi(1:4))
!       M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
!       Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
! 
!       LO_Decay_UnPol = LO_Decay_UnPol * Msq_T_BWENU
! 
! 
!       LO_Res_UnPol = LO_Res_UnPol * LO_Decay_UnPol


!------------ 1 LOOP --------------
ELSEIF( Correction.EQ.1 ) THEN
   do iHel=nHel(1),nHel(2)
!    do iHel=1,1; print *, "fixed heli"
      call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExt(1:4,5:9))
      call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,10:14))
      call HelCrossing(Helicities(iHel,1:4))
      call SetPolarizations()
      do iPrimAmp=1,6
          call EvalTree(BornAmps(iPrimAmp))
      enddo

      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol



!------------ bosonic loops --------------
       do iPrimAmp=1,6
! print *, "hel",ihel,"primamp no",iPrimAmp
          call SetKirill(PrimAmps(iPrimAmp))
          call PentCut(PrimAmps(iPrimAmp))
          call QuadCut(PrimAmps(iPrimAmp))
          call TripCut(PrimAmps(iPrimAmp))
          call DoubCut(PrimAmps(iPrimAmp))
          call SingCut(PrimAmps(iPrimAmp))
          call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
          call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen**2)
          PrimAmps(iPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(iPrimAmp)%Result(-2:1)
          call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,2,rdiv(2),rdiv(1))
!           call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv)
! print *, "gauge inv. fails for XTopDK=3", PrimAmps(iPrimAmp)%Result(-2:-1)
! print *, "------"
! print *, "LO", BornAmps(iPrimAmp)%Result
! print *, "NLO",PrimAmps(iPrimAmp)%Result(-2)
! print *, "NLO",PrimAmps(iPrimAmp)%Result(-1)
! print *, "NLO",PrimAmps(iPrimAmp)%Result(0)+PrimAmps(iPrimAmp)%Result(1)
! print *, "ratio",PrimAmps(iPrimAmp)%Result(-2)/BornAmps(iPrimAmp)%Result
! print *, "ratio",PrimAmps(iPrimAmp)%Result(-1)/BornAmps(iPrimAmp)%Result
! pause
!DEC$ IF (_QuadPrecImpr==1)
          AccPoles = CheckPoles(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv(1:2))
! print *, "DP AccPoles",AccPoles
          if( AccPoles.gt.1d-4 ) then
!               call PentCut_128(PrimAmps(iPrimAmp))
              call QuadCut_128(PrimAmps(iPrimAmp))
              call TripCut_128(PrimAmps(iPrimAmp))
              call DoubCut_128(PrimAmps(iPrimAmp))
              call SingCut_128(PrimAmps(iPrimAmp))
              call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
              call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen*2)
              PrimAmps(iPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(iPrimAmp)%Result(-2:1)
              AccPoles = CheckPoles(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv(1:2))
! print *, "QP AccPoles",AccPoles;pause
              if( AccPoles.gt.5d-2 ) then
                  print *, "SKIP",AccPoles
                  call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
                  EvalCS_1L_ststbgg = 0d0
                  SkipCounter = SkipCounter + 1
                  return
              endif
          endif
!DEC$ ENDIF
      enddo


      NLO_Res_Pol(-2:1) = (0d0,0d0)
      do jPrimAmp=1,6
      do iPrimAmp=1,NumBornAmps
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + Col1L_ttbgg(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(PrimAmps(jPrimAmp)%Result(-2:1)) )
      enddo
      enddo
      NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1) + NLO_Res_Pol(-2:1)


! ------------ fermionic loops --------------
      do iPrimAmp=7,12; !12; print *, "only primamp",iPrimAmp
          call SetKirill(PrimAmps(iPrimAmp))

          call PentCut(PrimAmps(iPrimAmp))
          call QuadCut(PrimAmps(iPrimAmp))
          call TripCut(PrimAmps(iPrimAmp))
          call DoubCut(PrimAmps(iPrimAmp))
          call SingCut(PrimAmps(iPrimAmp))
          call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)

          PrimAmps(iPrimAmp)%Result(-2:1) = -(0d0,1d0)*PrimAmps(iPrimAmp)%Result(-2:1) !minus is from closed fermion loop
          call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,2,rdiv(2),rdiv(1))
!           call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv)
! print *, "------",iPrimAmp
! print *, "LO", BornAmps(1)%Result
! print *, "NLO(-2)",PrimAmps(iPrimAmp)%Result(-2)
! print *, "NLO(-1)",PrimAmps(iPrimAmp)%Result(-1)
! print *, "NLO(0,1)",PrimAmps(iPrimAmp)%Result(0),PrimAmps(iPrimAmp)%Result(1)
! pause
      enddo
      FermionLoopPartAmp(7,-2:1) = Nf_light*PrimAmps(7)%Result(-2:1) + PrimAmps(9)%Result(-2:1)  - PrimAmps(11)%Result(-2:1) ! minus sign to remove minus from closed fermion loop
      FermionLoopPartAmp(8,-2:1) = Nf_light*PrimAmps(8)%Result(-2:1) + PrimAmps(10)%Result(-2:1) - PrimAmps(12)%Result(-2:1) 

      NLO_Res_Pol(-2:1) = (0d0,0d0)
      do jPrimAmp=7,8
      do iPrimAmp=1,NumBornAmps
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + Col1L_ttbgg(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result*dconjg(FermionLoopPartAmp(jPrimAmp,-2:1)) )
      enddo
      enddo
      NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) + NLO_Res_Pol(-2:1)

! print *, "hel",ihel
! print *, "LO",BornAmps(1)%Result
! print *, "NLO",PrimAmps(7)%Result(-2:1)
! print *, "ratio",PrimAmps(7)%Result(-1)/BornAmps(1)%Result
! pause

   enddo! helicity loop

ENDIF


IF( Correction.EQ.0 ) THEN
!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * WidthExpansion
   EvalCS_1L_ststbgg = LO_Res_Unpol * PreFac

ELSEIF( Correction.EQ.1 ) THEN
!  overall normalization: (4*Pi)^eps/Gamma(1-eps)
!  CT contributions                           ! beta        !stop WFRC
   NLO_Res_UnPol(-1) = NLO_Res_UnPol(-1) + (-11d0/3d0*3d0   +  0d0  )*LO_Res_Unpol
   NLO_Res_UnPol_Ferm(-1) = NLO_Res_UnPol_Ferm(-1) + (+2d0/3d0*Nf_light+2d0/3d0+1d0/6d0)*LO_Res_Unpol! coupling renorm.
   NLO_Res_UnPol_Ferm(-1) = NLO_Res_UnPol_Ferm(-1) + (-2d0/6d0         -2d0/3d0)*LO_Res_Unpol! gluon self energy

   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-1d0/6d0 )*2d0*dlog(MuRen/m_stop)*LO_Res_Unpol  ! finite log(mu) contrib. from  stop WFRC
!    NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-7d0/2d0*8d0/3d0 )*LO_Res_Unpol   ! finite contribution from stop WFRC's   ! set to zero in MSBAR
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + 1d0*LO_Res_Unpol  ! shift alpha_s^DR --> alpha_s^MSbar


! call Catani_1L_ststbgg(MomExt(1:4,1:4),CataniRes(-2:-1))!  alpha_s renormalization needs to be switched off
! print *, "LO",LO_Res_Unpol
! print *, "bare NLO(-2)",(NLO_Res_UnPol(-2)/LO_Res_Unpol)
! print *, "bare NLO(-1)",(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/LO_Res_Unpol
! print *, "bare NLO( 0)",(NLO_Res_UnPol(0)+NLO_Res_UnPol_Ferm(0))/LO_Res_Unpol
! print *, "bare NLO( 1)",(NLO_Res_UnPol(1)+NLO_Res_UnPol_Ferm(1))/LO_Res_Unpol
! print *, "Catani(-2)",dreal(CataniRes(-2))/LO_Res_Unpol
! print *, "Catani(-1)",real(CataniRes(-1))/LO_Res_Unpol
! print *, "ratio",real(CataniRes(-1)/LO_Res_Unpol)/((NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/LO_Res_Unpol)
! print *, "diff ",real(CataniRes(-1)/LO_Res_Unpol)-((NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/LO_Res_Unpol)
! pause


!  normalization
   LO_Res_Unpol = LO_Res_Unpol                         * ISFac * (alpha_s4Pi*RunFactor)**2
   NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1)           * ISFac * (alpha_s4Pi*RunFactor)**2 * alpha_sOver2Pi*RunFactor
   NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) * ISFac * (alpha_s4Pi*RunFactor)**2 * alpha_sOver2Pi*RunFactor

   EvalCS_1L_ststbgg = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac

!print *, "1L check",NLO_Res_UnPol(-2)/(alpha_sOver2Pi*RunFactor)/LO_Res_Unpol
!print *, "1L check",(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/(alpha_sOver2Pi*RunFactor)/LO_Res_Unpol


ELSEIF( Correction.EQ.3 ) THEN

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   IF( XTOPDECAYS.GE.1 ) THEN
       xE = yRnd(17+HelSampling)
   ELSEIF( XTOPDECAYS.EQ.0 ) THEN
       xE = yRnd(5+HelSampling)
   ENDIF
!    xe = 0.35d0; print *, "fixed xe"

   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
   call EvalIntDipoles_GGSTSTBG((/MomExt(1:4,3),MomExt(1:4,4),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,5:14),xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**3
   EvalCS_1L_ststbgg = HOp(1)    * pdf(0,1)  * pdf(0,2)   &
                     + HOp(2)/xE * pdf_z(0,1)* pdf(0,2)   &
                     + HOp(3)/xE * pdf(0,1)  * pdf_z(0,2)

!    print *, "LO cross",LO_Res_Unpol
!    print *, "HOp check",HOp(1)/(alpha_sOver2Pi*RunFactor)/LO_Res_Unpol
!    pause

    EvalCS_1L_ststbgg = EvalCS_1L_ststbgg * PreFac
ENDIF



   if( IsNan(EvalCS_1L_ststbgg) ) then
        print *, "NAN:",EvalCS_1L_ststbgg
        print *, yRnd(:)
        print *, NLO_Res_UnPol(0),NLO_Res_UnPol(1),NLO_Res_UnPol_Ferm(0),NLO_Res_UnPol_Ferm(1)
        print *, PSWgt , VgsWgt , PDFFac, sHatJacobi
        print *, eta1,eta2,((1d0-eta1)*xE+eta1),((1d0-eta2)*xE+eta2),MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,1:11)
        print *, "SKIP EVENT!!!!!"
        EvalCS_1L_ststbgg = 0d0
        return
   endif

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_ststbgg)
   enddo
   EvalCounter = EvalCounter + 1

   EvalCS_1L_ststbgg = EvalCS_1L_ststbgg/VgsWgt

return
END FUNCTION




SUBROUTINE Catani_1L_ststbgg(Mom,res)
use ModParameters
use ModMisc
use ModAmplitudes
use ModIntDipoles_GGSTSTBG
implicit none
real(8) :: Mom(1:4,4)
real(8) :: T1T2(1:2),T1T3(1:2),T1T4(1:2),TreeCol(1:2),beta0
complex(8) :: res(-2:-1),s12,s13,s14,s34,rdiv2,rdiv1
real(8) :: Tree_ij(0:6),Tree,Tree_12,Tree_13,Tree_14,Tree_23,Tree_24,Tree_34
complex(8) :: TreeMom(1:4,1:4)


  res(:) = (0d0,0d0)


!  tree momenta for g g -> t tb
  TreeMom(1:4,1) = dcmplx( mom(1:4,1) )
  TreeMom(1:4,2) = dcmplx( mom(1:4,2) )
  TreeMom(1:4,3) = dcmplx( mom(1:4,4) )
  TreeMom(1:4,4) = dcmplx( mom(1:4,3) )
  Tree_ij= Tree_GG_TTb_ij(TreeMom,(/0d0,0d0,m_STop**2,m_STop**2/))  /(alpha_s4Pi**2)*(4d0*8d0*8d0)
  Tree_12  =  Tree_ij(1) * 0.5d0
  Tree_13  =  Tree_ij(2) * 0.5d0
  Tree_14  =  Tree_ij(3) * 0.5d0
  Tree = Tree_ij(0)
!   Tree_23  =  Tree_ij(4) * 0.5d0
!   Tree_24  =  Tree_ij(5) * 0.5d0
   Tree_34  =  Tree_ij(6) * 0.5d0

! print *, "checker",Tree_12
! print *, "checker",Tree_34
! print *, "checker",(Tree_12+Tree_13+Tree_14)/Tree

  s12=2d0*(Mom(1:4,1).dot.Mom(1:4,2)) + (0D0,1D-15)
  s13=2d0*(Mom(1:4,1).dot.Mom(1:4,4)) + (0D0,1D-15)
  s14=2d0*(Mom(1:4,1).dot.Mom(1:4,3)) + (0D0,1D-15)
  s34=2d0*(Mom(1:4,3).dot.Mom(1:4,4)) + (0D0,1D-15)

  rdiv1=0d0; rdiv2=0d0
  call J_sing(s12,0d0,0d0,MuRen**2,rdiv2,rdiv1)
  res(-2:-1) = res(-2:-1) + 2*(/-rdiv2,-rdiv1/)*Tree_12

  rdiv1=0d0; rdiv2=0d0
  call J_sing(s34,m_stop,m_stop,MuRen**2,rdiv2,rdiv1)
  res(-2:-1) = res(-2:-1) + 2*(/-rdiv2,-rdiv1/)*Tree_34

  rdiv1=0d0; rdiv2=0d0
  call J_sing(s13,m_stop,0d0,MuRen**2,rdiv2,rdiv1)
  res(-2:-1) = res(-2:-1) + 4*(/-rdiv2,-rdiv1/)*Tree_13

  rdiv1=0d0; rdiv2=0d0
  call J_sing(s14,m_stop,0d0,MuRen**2,rdiv2,rdiv1)
  res(-2:-1) = res(-2:-1) + 4*(/-rdiv2,-rdiv1/)*Tree_14



 rdiv1=-11d0/6d0*3d0 + 2d0/3d0/2d0*5d0
 res(-1) = res(-1) + rdiv1*Tree

 rdiv1=-11d0/6d0*3d0 + 2d0/3d0/2d0*5d0
 res(-1) = res(-1) + rdiv1*Tree

 rdiv1=-4d0/3d0 
 res(-1) = res(-1) + rdiv1*Tree

 rdiv1=-4d0/3d0 
 res(-1) = res(-1) + rdiv1*Tree


 beta0 = 11d0/3d0*3d0 -2d0/3d0*6d0 - 1d0/6d0
  res(-1) = res(-1) + beta0*Tree


RETURN
END SUBROUTINE











FUNCTION EvalCS_1L_ststbqqb(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModIntDipoles_QQBSTSTBG
use ModIntDipoles_QGSTSTBQ
implicit none
real(8) ::  EvalCS_1L_ststbqqb,yRnd(1:VegasMxDim),VgsWgt,HOp(1:3)
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1),FermionLoopPartAmp(1:2,-2:1)
complex(8) :: LO_Decay_UnPol,Spi(1:4),BarSpi(1:4),M_T_BWENU,Msq_T_BWENU
integer :: iHel,jHel,kHel,lHel,iPrimAmp,jPrimAmp
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5,ISFac
complex(8) :: BosonicPartAmp(1:2,-2:1),FermionPartAmp(1:2,-2:1),CataniRes(-2:-1),dZStop(-1:0)
real(8) :: MomExt(1:4,1:15),MomP(1:4,1:4),MomBoost(1:4)
logical :: applyPSCut
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac_a,PDFFac_b,PDFFac,AccPoles
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2),xE,sigmaTot,beta
integer :: NBin(1:NumMaxHisto),NHisto,nHel(1:2),npdf
real(8),parameter :: Nc=3d0
type(Particle) :: TopQuark(1:2)
!include 'misc/global_import'
include 'vegas_common.f'

   EvalCS_1L_ststbqqb = 0d0
   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if( EHat.le.2d0*m_STop ) then
      EvalCS_1L_ststbqqb = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to2Stops(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt)! AStop, Stop
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
! print *, "switch off boost for comparison with Radja"

   IF(XTOPDECAYS.EQ.3) THEN
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  Chi top
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,4),yRnd(7:8),MomExt(1:4,10:11),PSWgt3)
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd( 9:12),MomExt(1:4,7:9),PSWgt4)! bot lep neu
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(13:16),MomExt(1:4,12:14),PSWgt5)
      PSWgt = PSWgt * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
      nHel(1:2) = getHelicity(yrnd(17))
   ELSE
      nHel(1:2) = getHelicity(yrnd(5))
   ENDIF
 
   call Kinematics_TTbarETmiss(.false.,MomExt,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_1L_ststbqqb = 0d0
      return
   endif

   call SetPDFs(eta1,eta2,MuFac,pdf)
   PDFFac_a = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
            + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
            + pdf(Bot_,1)*pdf(ABot_,2)
   PDFFac_b = pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
            + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
            + pdf(Bot_,2)*pdf(ABot_,1)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   RunFactor = RunAlphaS(NLOParam,MuRen)
   PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))

!!!!! for of total cross section (fix y1,y2!)
!   beta = dsqrt(1d0-4d0*m_stop**2/EHat**2)
!   SigmaTot = (alpha_s*RunFactor)**2*DblPi/EHat**2*( 2d0/27d0 * beta**3  )
!   SigmaTot = SigmaTot * fbGeV2*SHatJacobi
!   PreFac = PreFac/SigmaTot
!   pdffac_a=0d0
!   pdffac_b=1d0
!!!!!!

  LO_Res_Unpol             = (0d0,0d0)
  NLO_Res_Unpol(-2:1)      = (0d0,0d0)
  NLO_Res_Unpol_Ferm(-2:1) = (0d0,0d0)

 do npdf=1,2;! print *, "no npdf loop"
    if(npdf.eq.1) then
        PDFFac = PDFFac_a
    elseif(npdf.eq.2) then
        PDFFac = PDFFac_b
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
    endif
! PDFFac=1d0; print *, "pdfs set to one"
    ISFac = MomCrossing(MomExt)
    call InitCurrCache()
    call SetPropagators()

!------------ LO ----------------
IF( Correction.EQ.0 ) THEN


!    do iHel=nHel(1),nHel(2)
! 
!       if( Helicities(iHel,3).eq.Helicities(iHel,4) ) cycle ! remove helicities where IS quarks have ++ / -- helicities because tree vanishes
! 
!       call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExt(1:4,5:9))
!       call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,10:14))
! 
!       call HelCrossing(Helicities(iHel,1:4))
!       call SetPolarizations()
!       do iPrimAmp=1,NumBornAmps
!           call EvalTree(BornAmps(iPrimAmp))
!       enddo
! 
!       LO_Res_Pol = (0d0,0d0)
!       do jPrimAmp=1,NumBornAmps
!       do iPrimAmp=1,NumBornAmps
!           LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
!       enddo
!       enddo
!       LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol*PDFFac
!    enddo!helicity loop





! ! !  this LO code block is equivalent to the above but it is 4 times faster
!  spin correlated top quarks
   if( XTopDecays.ne.0 ) then
        LO_Decay_UnPol = (0d0,0d0)
        do jHel=-1,1,2
        do kHel=-1,1,2
            call STopDecay(ExtParticle(1),DKX_STChi0_LO,jHel,MomExt(1:4,5:9))
            call STopDecay(ExtParticle(2),DKX_STChi0_LO,kHel,MomExt(1:4,10:14))
            LO_Decay_UnPol = LO_Decay_UnPol + ExtParticle(1)%Pol(1)*dconjg(ExtParticle(1)%Pol(1)) * ExtParticle(2)%Pol(1)*dconjg(ExtParticle(2)%Pol(1))
        enddo
        enddo
   else
        LO_Decay_UnPol = (1d0,0d0)
   endif


! ! ! !  spin uncorrelated top quarks
!    LO_Decay_UnPol = (0d0,0d0)
!    do iHel=-1,1,2
!    do jHel=-1,1,2
!    do kHel=-1,1,2
!    do lHel=-1,1,2
!       call STopDecay(ExtParticle(1),DKX_STChi0_LO,jHel,MomExt(1:4,5:9),HelTop=iHel)
!       call STopDecay(ExtParticle(2),DKX_STChi0_LO,kHel,MomExt(1:4,10:14),HelTop=lHel)
!       LO_Decay_UnPol = LO_Decay_UnPol + ExtParticle(1)%Pol(1)*dconjg(ExtParticle(1)%Pol(1)) * ExtParticle(2)%Pol(1)*dconjg(ExtParticle(2)%Pol(1))
!    enddo
!    enddo
!    enddo
!    enddo
! 
! ! spherical decay of ATop
!       TopQuark(1)%Mom(1:4) = MomExt(1:4,6)
!       TopQuark(1)%PartType = ATop_
!       TopQuark(1)%Mass = m_Top
!       TopQuark(1)%Mass2 = m_Top**2
! 
!       TopQuark(2)%Mom(1:4) = MomExt(1:4,11)
!       TopQuark(2)%PartType = Top_
!       TopQuark(2)%Mass = m_Top
!       TopQuark(2)%Mass2 = m_Top**2
! 
!       call TopDecay(TopQuark(1),DK_LO,MomExt(1:4,7:9))
!       Spi(1:4) = TopQuark(1)%Pol(1:4)
!       call vBarSpi(TopQuark(1)%Mom(1:4),m_Top,-1,BarSpi(1:4))
!       M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
!       Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0
! 
!       call vBarSpi(TopQuark(1)%Mom(1:4),m_Top,+1,BarSpi(1:4))
!       M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
!       Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
! 
!       LO_Decay_UnPol = LO_Decay_UnPol * Msq_T_BWENU
! 
! ! spherical decay of Top
!       call TopDecay(TopQuark(2),DK_LO,MomExt(1:4,12:14))
!       BarSpi(1:4) = TopQuark(2)%Pol(1:4)
!       call uSpi(TopQuark(2)%Mom(1:4),m_Top,-1,Spi(1:4))
!       M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
!       Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0
! 
!       call uSpi(TopQuark(2)%Mom(1:4),m_Top,+1,Spi(1:4))
!       M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
!       Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
! 
!       LO_Decay_UnPol = LO_Decay_UnPol * Msq_T_BWENU


! !  production process
   if( XTopDecays.eq.0 ) then
      jHel=1
   else
      jHel=4
   endif

   do iHel=nHel(1),nHel(2),jHel! ,4 because we jump over stop polarizations which are summed over outside
      if( Helicities(iHel,3).eq.Helicities(iHel,4) ) cycle ! remove helicities where IS quarks have ++ / -- helicities because tree vanishes

      ExtParticle(1)%Pol(1)=(1d0,0d0)
      ExtParticle(2)%Pol(1)=(1d0,0d0)
      call HelCrossing(Helicities(iHel,1:4))
      call SetPolarizations()
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo

      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol * LO_Decay_UnPol * PDFFac


   enddo!helicity loop



!------------ 1 LOOP --------------
ELSEIF( Correction.EQ.1 ) THEN
    do iHel=nHel(1),nHel(2)  !; print *, "eval hel2 only"
      if( Helicities(iHel,3).eq.Helicities(iHel,4) ) cycle ! remove helicities where IS quarks have ++ / -- helicities because tree vanishes

      call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExt(1:4,5:9))
      call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,10:14))

      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      do iPrimAmp=1,4
          call EvalTree(BornAmps(iPrimAmp))
      enddo
      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol*PDFFac


! ------------ bosonic loops --------------
      do iPrimAmp=1,4; !print *, "Evaluate bosonic loop primamp",iPrimAmp
          call SetKirill(PrimAmps(iPrimAmp))
          call PentCut(PrimAmps(iPrimAmp))
          call QuadCut(PrimAmps(iPrimAmp))
          call TripCut(PrimAmps(iPrimAmp))
          call DoubCut(PrimAmps(iPrimAmp))
          call SingCut(PrimAmps(iPrimAmp))
          call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)

!         call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen**2)
          PrimAmps(iPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(iPrimAmp)%Result(-2:1)
          call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,2,rdiv(2),rdiv(1))!    this fails for iPrimAmp=3 because ssss vertex is not accounted for in OneLoopDiv
!           call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv)
! print *, "LO",BornAmps(iPrimAmp)%Result
! print *, "NLO",PrimAmps(iPrimAmp)%Result(-2:1)
! print *, "ratio",PrimAmps(iPrimAmp)%Result(-2)/BornAmps(iPrimAmp)%Result
! print *, "ratio",PrimAmps(iPrimAmp)%Result(-1)/BornAmps(iPrimAmp)%Result
! print *, "ratio",PrimAmps(iPrimAmp)%Result(0)/BornAmps(1)%Result
! print *, "ratio",PrimAmps(iPrimAmp)%Result(1)/BornAmps(1)%Result
! pause

!DEC$ IF (_QuadPrecImpr==1)
          AccPoles = CheckPoles(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv(1:2))

 if( iPrimAmp.eq.3 ) AccPoles = 0d0!   remove checking PrimAmp no.3 because ssss vertex is not accounted for in OneLoopDiv

          if( AccPoles.gt.1d-4 ) then
              call PentCut_128(PrimAmps(iPrimAmp))
              call QuadCut_128(PrimAmps(iPrimAmp))
              call TripCut_128(PrimAmps(iPrimAmp))
              call DoubCut_128(PrimAmps(iPrimAmp))
              call SingCut_128(PrimAmps(iPrimAmp))
              call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
!               call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen*2)
              PrimAmps(iPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(iPrimAmp)%Result(-2:1)
              AccPoles = CheckPoles(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv(1:2))
              if( AccPoles.gt.1d-2 ) then
                  print *, "SKIP",AccPoles
                  call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
                  EvalCS_1L_ststbqqb = 0d0
                  SkipCounter = SkipCounter + 1
                  return
              endif
          endif
!DEC$ ENDIF
      enddo


      BosonicPartAmp(1,-2:1) =  &
                             + Nc * PrimAmps(PrimAmp1_1234)%Result(-2:1) &
                             - 2d0/Nc*( PrimAmps(PrimAmp1_1234)%Result(-2:1) + PrimAmps(PrimAmp1_1243)%Result(-2:1) ) &
                             - 1d0/Nc * ( PrimAmps(PrimAmp3_1432)%Result(-2:1) + PrimAmps(PrimAmp4_1234)%Result(-2:1) )

      BosonicPartAmp(2,-2:1) = PrimAmps(PrimAmp1_1243)%Result(-2:1)  &
                             + 1d0/Nc**2 * ( PrimAmps(PrimAmp1_1234)%Result(-2:1) + PrimAmps(PrimAmp1_1243)%Result(-2:1) )  &
                             + 1d0/Nc**2 * ( PrimAmps(PrimAmp3_1432)%Result(-2:1) + PrimAmps(PrimAmp4_1234)%Result(-2:1) )


      NLO_Res_Pol(-2:1) = (0d0,0d0)
      do jPrimAmp=1,2
      do iPrimAmp=1,NumBornAmps
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + Col1L_ttbqqb(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(BosonicPartAmp(jPrimAmp,-2:1)) )
      enddo
      enddo
      NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1) + NLO_Res_Pol(-2:1)*PDFFac

! ------------ fermionic loops --------------
      do iPrimAmp=5,7; !print *, "Evaluate fermion loop primamp",iPrimAmp
          call SetKirill(PrimAmps(iPrimAmp))
          call PentCut(PrimAmps(iPrimAmp))
          call QuadCut(PrimAmps(iPrimAmp))
          call TripCut(PrimAmps(iPrimAmp))
          call DoubCut(PrimAmps(iPrimAmp))
          call SingCut(PrimAmps(iPrimAmp))
          call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
          PrimAmps(iPrimAmp)%Result(-2:1) = -(0d0,1d0)*PrimAmps(iPrimAmp)%Result(-2:1) !minus from closed fermion loop
          call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,2,rdiv(2),rdiv(1))
!           call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(1),rdiv)
! print *, "LO",BornAmps(1)%Result
! print *, "NLO",PrimAmps(iPrimAmp)%Result(-2:1)
! print *, "ratio",PrimAmps(iPrimAmp)%Result(-2)/BornAmps(1)%Result
! print *, "ratio",PrimAmps(iPrimAmp)%Result(-1),BornAmps(1)%Result
! print *, "ratio",PrimAmps(iPrimAmp)%Result(0),BornAmps(1)%Result
! print *, "ratio",PrimAmps(iPrimAmp)%Result(1),BornAmps(1)%Result
! pause
      enddo
                         ! minus to remove minus from closed fermion loop
      FermionLoopPartAmp(1,-2:1) =           ( Nf_light*PrimAmps(5)%Result(-2:1) + PrimAmps(6)%Result(-2:1) - PrimAmps(7)%Result(-2:1) )
      FermionLoopPartAmp(2,-2:1) = -1d0/Nc * ( Nf_light*PrimAmps(5)%Result(-2:1) + PrimAmps(6)%Result(-2:1) - PrimAmps(7)%Result(-2:1) )

! print *, "For comparison with RADJA set mt=0 in closed loops"
! FermionLoopPartAmp(1,-2:1) =           ( (Nf_light)*PrimAmps(5)%Result(-2:1) - PrimAmps(7)%Result(-2:1) )
! FermionLoopPartAmp(2,-2:1) = -1d0/Nc * ( (Nf_light)*PrimAmps(5)%Result(-2:1) - PrimAmps(7)%Result(-2:1) )


      NLO_Res_Pol(-2:1) = (0d0,0d0)
      do jPrimAmp=1,2
      do iPrimAmp=1,NumBornAmps
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + Col1L_ttbqqb(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(FermionLoopPartAmp(jPrimAmp,-2:1)) )
      enddo
      enddo
      NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) + NLO_Res_Pol(-2:1)*PDFFac

   enddo!helicity loop
ENDIF! Correction loop
  enddo! npdf loop
  call swapMom(MomExt(1:4,1),MomExt(1:4,2))   ! swap back to original order, for ID below
! print *, "npdf swap off"

! Radja's results with nf=5:
! E1 = -56.3557171738789
! E2 = -75.9457406360220
! E3 = -76.5210321061188
! 
! ns = 0.575291470097(R)
! ns = 0.575291470094(M)

! nf = 19.5900234621(R)
! nf = 19.5900234621448(M)

! rest = -76.5210321061188(R)


IF( Correction.EQ.0 ) THEN
!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * WidthExpansion
   EvalCS_1L_ststbqqb = LO_Res_Unpol * PreFac

ELSEIF( Correction.EQ.1 ) THEN

! NLO_Res_UnPol(0) = NLO_Res_UnPol(0) +  1d0*LO_Res_Unpol! shift alpha_s from FDH to CDR
! NLO_Res_UnPol(0) = NLO_Res_UnPol(0) -  4d0/3d0*LO_Res_Unpol! shift anom.dim. from FDH to CDR
! NLO_Res_UnPol(0) = NLO_Res_UnPol(0) -  NLO_Res_UnPol(-2)*dblPi**2/12d0!  convert from my normalization to radja's 

! print *, ""
! print *, "LO",LO_Res_Unpol* ISFac/5.3333333333333333d-2
! print *, "LO",LO_Res_Unpol* ISFac/2.700617283950616d-002
! print *, "LO",LO_Res_Unpol* ISFac/8.163265306122448d-002
! print *, ""
! print *, "bare NLO(-2)",2d0*(NLO_Res_UnPol(-2)/LO_Res_Unpol)
! print *, "bare NLO(-1)",2d0*(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/LO_Res_Unpol-44.2491195362494d0
! print *, "bare NLO(-1)",2d0*(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/LO_Res_Unpol-47.0012561334728d0
! print *, "bare NLO(-1)",2d0*(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/LO_Res_Unpol-51.5866549563775d0
! print *, ""
! print *, "bare NLO(0)",2d0*(NLO_Res_UnPol(0)+NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol(1)+NLO_Res_UnPol_Ferm(1))/LO_Res_Unpol+56.3557171738789d0
! print *, "bare NLO(0)",2d0*(NLO_Res_UnPol(0)+NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol(1)+NLO_Res_UnPol_Ferm(1))/LO_Res_Unpol+67.6225159647733d0
! print *, "bare NLO(0)",2d0*(NLO_Res_UnPol(0)+NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol(1)+NLO_Res_UnPol_Ferm(1))/LO_Res_Unpol+84.9030724385563d0
! print *, ""



!  overall normalization: (4*Pi)^eps/Gamma(1-eps)

   call deltaZ_Stop(dZStop(-1:0))

!  CT contributions
   NLO_Res_UnPol(-1) = NLO_Res_UnPol(-1) + (-11d0/3d0*3d0   +  2d0*dZStop(-1) )*LO_Res_Unpol
                                                                                          ! the last term is from the scalar loop
   NLO_Res_UnPol_Ferm(-1) = NLO_Res_UnPol_Ferm(-1) + (2d0/3d0*Nf_light+  2d0/3d0*Nf_heavy + 1d0/6d0)*LO_Res_Unpol
   NLO_Res_UnPol_Ferm( 0) = NLO_Res_UnPol_Ferm( 0) + (2d0/3d0*Nf_heavy)*2d0*dlog(MuRen/m_top)*LO_Res_Unpol  ! finite log(mu2) contrib. from heavy flavor in alpha_s ren.
   NLO_Res_UnPol_Ferm( 0) = NLO_Res_UnPol_Ferm( 0) +          (1d0/6d0)*2d0*dlog(MuRen/m_Stop)*LO_Res_Unpol  ! finite log(mu2) contrib. from heavy flavor in alpha_s ren.
!    NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + dZStop(0)*LO_Res_Unpol ! finite contribution from dZStop !  set to zero in MSBAR


   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + LO_Res_Unpol ! shift alpha_s^DR --> alpha_s^MSbar

!  print *, "ren. NLO(-2)",2d0*NLO_Res_UnPol(-2)/LO_Res_Unpol
!  print *, "ren. NLO(-1)",2d0*(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/LO_Res_Unpol-29.2491195362494d0
!  print *, "ren. NLO(-1)",2d0*(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/LO_Res_Unpol-32.0012561334728d0
!  print *, "ren. NLO(-1)",2d0*(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/LO_Res_Unpol-36.5866549563775d0
!  print *, "ren. NLO(0)",2d0*(NLO_Res_UnPol(0)+NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol(1)+NLO_Res_UnPol_Ferm(1))/LO_Res_Unpol+56.3557171738789d0
!  print *, "ren. NLO(0)",2d0*(NLO_Res_UnPol(0)+NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol(1)+NLO_Res_UnPol_Ferm(1))/LO_Res_Unpol+67.6225159647733d0
!  print *, "ren. NLO(0)",2d0*(NLO_Res_UnPol(0)+NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol(1)+NLO_Res_UnPol_Ferm(1))/LO_Res_Unpol+84.9030724385563d0
!  pause


! call Catani_1L_ststbqqb(MomExt(1:4,1:4),CataniRes(-2:-1))!  alpha_s renormalization needs to be switched off
! print *, "LO",LO_Res_Unpol
! print *, "bare NLO(-2)",(NLO_Res_UnPol(-2)/LO_Res_Unpol)
! print *, "bare NLO(-1)",(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/LO_Res_Unpol
! print *, "Catani(-2)",dreal(CataniRes(-2))
! print *, "Catani(-1)",real(CataniRes(-1))
! print *, "ratio",real(CataniRes(-1))/((NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/LO_Res_Unpol)
! print *, "diff ",real(CataniRes(-1))-((NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/LO_Res_Unpol)
! pause


!  normalization
   LO_Res_Unpol = LO_Res_Unpol                         * ISFac * (alpha_s4Pi*RunFactor)**2
   NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1)           * ISFac * (alpha_s4Pi*RunFactor)**2 * alpha_sOver2Pi*RunFactor
   NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) * ISFac * (alpha_s4Pi*RunFactor)**2 * alpha_sOver2Pi*RunFactor

! print *, "virt 1/eps2",NLO_Res_UnPol(-2)/(alpha_sOver2Pi*RunFactor)/LO_Res_Unpol
! print *, "virt 1/eps",(dble(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1)))/(alpha_sOver2Pi*RunFactor)/LO_Res_Unpol
! pause

   EvalCS_1L_ststbqqb = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac



ELSEIF( Correction.EQ.3 ) THEN

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   IF( XTOPDECAYS.GE.1 ) THEN
       xE = yRnd(17+HelSampling)
   ELSEIF( XTOPDECAYS.EQ.0 ) THEN
       xE = yRnd(5+HelSampling)
   ENDIF
! xe = 0.35d0; print *, "fixed xe"

   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
   IF( PROCESS.eq.56 ) THEN
      call EvalIntDipoles_QQBSTSTBG((/MomExt(1:4,3),MomExt(1:4,4),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,5:14),xE,HOp(1:3))
! print *, "IntDip",HOp(1)*RunFactor**3/(alpha_sOver2Pi*RunFactor)/LO_Res_Unpol
! pause
! IF(PROCESS.eq.99999999999999 ) THEN
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ststbqqb = HOp(1) * (pdf(Up_,1)*pdf(AUp_,2)+pdf(Dn_,1)*pdf(ADn_,2)+pdf(Chm_,1)*pdf(AChm_,2)+pdf(Str_,1)*pdf(AStr_,2)+pdf(Bot_,1)*pdf(ABot_,2) ) &
                         + HOp(2)/xE * (pdf_z(Up_,1)*pdf(AUp_,2)+pdf_z(Dn_,1)*pdf(ADn_,2)+pdf_z(Chm_,1)*pdf(AChm_,2)+pdf_z(Str_,1)*pdf(AStr_,2)+pdf_z(Bot_,1)*pdf(ABot_,2) ) &
                         + HOp(3)/xE * (pdf(Up_,1)*pdf_z(AUp_,2)+pdf(Dn_,1)*pdf_z(ADn_,2)+pdf(Chm_,1)*pdf_z(AChm_,2)+pdf(Str_,1)*pdf_z(AStr_,2)+pdf(Bot_,1)*pdf_z(ABot_,2) )

      call EvalIntDipoles_QQBSTSTBG((/MomExt(1:4,3),MomExt(1:4,4),-MomExt(1:4,2),-MomExt(1:4,1)/),MomExt(1:4,5:14),xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ststbqqb = EvalCS_1L_ststbqqb &
                         + HOp(1)    * (pdf(Up_,2)*pdf(AUp_,1)+pdf(Dn_,2)*pdf(ADn_,1)+pdf(Chm_,2)*pdf(AChm_,1)+pdf(Str_,2)*pdf(AStr_,1)+pdf(Bot_,2)*pdf(ABot_,1) ) &
                         + HOp(2)/xE * (pdf_z(Up_,2)*pdf(AUp_,1)+pdf_z(Dn_,2)*pdf(ADn_,1)+pdf_z(Chm_,2)*pdf(AChm_,1)+pdf_z(Str_,2)*pdf(AStr_,1)+pdf_z(Bot_,2)*pdf(ABot_,1) ) &
                         + HOp(3)/xE * (pdf(Up_,2)*pdf_z(AUp_,1)+pdf(Dn_,2)*pdf_z(ADn_,1)+pdf(Chm_,2)*pdf_z(AChm_,1)+pdf(Str_,2)*pdf_z(AStr_,1)+pdf(Bot_,2)*pdf_z(ABot_,1) )


   ELSEIF( PROCESS.eq.53 ) THEN
      call EvalIntDipoles_QGSTSTBQ((/MomExt(1:4,3),MomExt(1:4,4),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,5:14),xE,HOp(1:3))

      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ststbqqb = HOp(1) * (pdf(Up_,1)*pdf(0,2)+pdf(Dn_,1)*pdf(0,2)+pdf(Chm_,1)*pdf(0,2)+pdf(Str_,1)*pdf(0,2)+pdf(Bot_,1)*pdf(0,2) ) &
                         + HOp(2)/xE * (pdf_z(Up_,1)*pdf(0,2)+pdf_z(Dn_,1)*pdf(0,2)+pdf_z(Chm_,1)*pdf(0,2)+pdf_z(Str_,1)*pdf(0,2)+pdf_z(Bot_,1)*pdf(0,2) ) &
                         + HOp(3)/xE * (pdf(Up_,1)*pdf_z(0,2)+pdf(Dn_,1)*pdf_z(0,2)+pdf(Chm_,1)*pdf_z(0,2)+pdf(Str_,1)*pdf_z(0,2)+pdf(Bot_,1)*pdf_z(0,2) )

      call EvalIntDipoles_QGSTSTBQ((/MomExt(1:4,3),MomExt(1:4,4),-MomExt(1:4,2),-MomExt(1:4,1)/),MomExt(1:4,5:14),xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac

      EvalCS_1L_ststbqqb = EvalCS_1L_ststbqqb  &
                         + HOp(1)    * (pdf(Up_,2)*pdf(0,1)+pdf(Dn_,2)*pdf(0,1)+pdf(Chm_,2)*pdf(0,1)+pdf(Str_,2)*pdf(0,1)+pdf(Bot_,2)*pdf(0,1) ) &
                         + HOp(2)/xE * (pdf_z(Up_,2)*pdf(0,1)+pdf_z(Dn_,2)*pdf(0,1)+pdf_z(Chm_,2)*pdf(0,1)+pdf_z(Str_,2)*pdf(0,1)+pdf_z(Bot_,2)*pdf(0,1) ) &
                         + HOp(3)/xE * (pdf(Up_,2)*pdf_z(0,1)+pdf(Dn_,2)*pdf_z(0,1)+pdf(Chm_,2)*pdf_z(0,1)+pdf(Str_,2)*pdf_z(0,1)+pdf(Bot_,2)*pdf_z(0,1) )



   ELSEIF( PROCESS.eq.54 ) THEN
      call EvalIntDipoles_QGSTSTBQ((/MomExt(1:4,3),MomExt(1:4,4),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,5:14),xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ststbqqb = HOp(1) * (pdf(AUp_,1)*pdf(0,2)+pdf(ADn_,1)*pdf(0,2)+pdf(AChm_,1)*pdf(0,2)+pdf(AStr_,1)*pdf(0,2)+pdf(ABot_,1)*pdf(0,2) ) &
                         + HOp(2)/xE * (pdf_z(AUp_,1)*pdf(0,2)+pdf_z(ADn_,1)*pdf(0,2)+pdf_z(AChm_,1)*pdf(0,2)+pdf_z(AStr_,1)*pdf(0,2)+pdf_z(ABot_,1)*pdf(0,2) ) &
                         + HOp(3)/xE * (pdf(AUp_,1)*pdf_z(0,2)+pdf(ADn_,1)*pdf_z(0,2)+pdf(AChm_,1)*pdf_z(0,2)+pdf(AStr_,1)*pdf_z(0,2)+pdf(ABot_,1)*pdf_z(0,2) )

      call EvalIntDipoles_QGSTSTBQ((/MomExt(1:4,3),MomExt(1:4,4),-MomExt(1:4,2),-MomExt(1:4,1)/),MomExt(1:4,5:14),xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ststbqqb = EvalCS_1L_ststbqqb &
                         + HOp(1) * (pdf(AUp_,1)*pdf(0,2)+pdf(ADn_,1)*pdf(0,2)+pdf(AChm_,1)*pdf(0,2)+pdf(AStr_,1)*pdf(0,2)+pdf(ABot_,1)*pdf(0,2) ) &
                         + HOp(2)/xE * (pdf_z(AUp_,1)*pdf(0,2)+pdf_z(ADn_,1)*pdf(0,2)+pdf_z(AChm_,1)*pdf(0,2)+pdf_z(AStr_,1)*pdf(0,2)+pdf_z(ABot_,1)*pdf(0,2) ) &
                         + HOp(3)/xE * (pdf(AUp_,1)*pdf_z(0,2)+pdf(ADn_,1)*pdf_z(0,2)+pdf(AChm_,1)*pdf_z(0,2)+pdf(AStr_,1)*pdf_z(0,2)+pdf(ABot_,1)*pdf_z(0,2) )
   ENDIF

!    print *, "LO cross",LO_Res_Unpol
!    print *, "HOp check",HOp(1)/(alpha_sOver2Pi*RunFactor)
!    pause

ENDIF


   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_ststbqqb)
   enddo
   EvalCounter = EvalCounter + 1


   EvalCS_1L_ststbqqb = EvalCS_1L_ststbqqb/VgsWgt

return
END FUNCTION





SUBROUTINE Catani_1L_ststbqqb(Mom,res)
use ModParameters
use ModMisc
use ModAmplitudes
implicit none
real(8) :: Mom(1:4,4)
real(8) :: T1T2,T1T3,T1T4,TreeCol,beta0
complex(8) :: res(-2:-1),s12,s13,s14,s34,rdiv2,rdiv1

  res(:) = (0d0,0d0)

  TreeCol=8d0
  T1T2 = 4d0/3d0
  T1T3 =-28d0/3d0
  T1T4 =-8d0/3d0

! print *,"check",(t1t2+t1t3+t1t4)/TreeCol,-4.0/3

  s12=2d0*(Mom(1:4,1).dot.Mom(1:4,2)) + (0D0,1D-15)
  s13=2d0*(Mom(1:4,1).dot.Mom(1:4,4)) + (0D0,1D-15)
  s14=2d0*(Mom(1:4,1).dot.Mom(1:4,3)) + (0D0,1D-15)
  s34=2d0*(Mom(1:4,3).dot.Mom(1:4,4)) + (0D0,1D-15)

  rdiv1=0d0; rdiv2=0d0
  call J_sing(s12,0d0,0d0,MuRen**2,rdiv2,rdiv1)
  res(-2:-1) = res(-2:-1) + 2*T1T2*(/-rdiv2,-rdiv1/)/TreeCol

  rdiv1=0d0; rdiv2=0d0
  call J_sing(s34,m_stop,m_stop,MuRen**2,rdiv2,rdiv1)
  res(-2:-1) = res(-2:-1) + 2*T1T2*(/-rdiv2,-rdiv1/)/TreeCol

  rdiv1=0d0; rdiv2=0d0
  call J_sing(s13,m_stop,0d0,MuRen**2,rdiv2,rdiv1)
  res(-2:-1) = res(-2:-1) + 4*T1T3*(/-rdiv2,-rdiv1/)/TreeCol

  rdiv1=0d0; rdiv2=0d0
  call J_sing(s14,m_stop,0d0,MuRen**2,rdiv2,rdiv1)
  res(-2:-1) = res(-2:-1) + 4*T1T4*(/-rdiv2,-rdiv1/)/TreeCol

 rdiv1=-4d0/3d0 *3d0/2d0
 res(-1) = res(-1) + rdiv1

 rdiv1=-4d0/3d0 *3d0/2d0
 res(-1) = res(-1) + rdiv1

 rdiv1=-4d0/3d0 
 res(-1) = res(-1) + rdiv1

 rdiv1=-4d0/3d0 
 res(-1) = res(-1) + rdiv1

                        ! Nf+NF=6
 beta0 = 11d0/3d0*3d0 -2d0/3d0*6d0 - 1d0/6d0
 res(-1) = res(-1) + beta0


RETURN
END SUBROUTINE




SUBROUTINE deltaZ_Stop(res)
use ModParameters
use ModMisc
use ModAmplitudes
implicit none
real(8) :: Mom(1:4,4)
complex(8) :: res(-1:0)
complex(8) :: qlI2

  res(-1) = -2d0 * qlI2(m_stop**2,0d0,m_stop**2,MuRen**2,-1)  &
            -4d0 * 0  &
            +2d0

  res( 0) = -2d0 * qlI2(m_stop**2,0d0,m_stop**2,MuRen**2,0)  &
            -4d0 * (-1d0)  &
            -1d0

  res(:) = res(:) * 4d0/3d0

RETURN
END SUBROUTINE






FUNCTION EvalCS_1L_ststbgggg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
implicit none
real(8) ::  EvalCS_1L_ststbgggg,yRnd(1:VegasMxDim),VgsWgt,HOp(1:3)
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1),FermionLoopPartAmp(7:8,-2:1)
integer :: iHel,jHel,kHel,iPrimAmp,jPrimAmp
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5,ISFac
real(8) :: MomExt(1:4,1:15),MomP(1:4,1:4),MomBoost(1:4)
logical :: applyPSCut
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,AccPoles
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2),xE,sigmaTot,beta
integer :: NBin(1:NumMaxHisto),NHisto,nHel(1:2)
!include 'misc/global_import'
include 'vegas_common.f'


   

   EvalCS_1L_ststbgggg = 0d0
   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if( EHat.le.2d0*m_STop ) then
      EvalCS_1L_ststbgggg = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)

!    call EvalPhaseSpace_2to2Stops(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt)! AStop, Stop
!    call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
   call EvalPhasespace_2toN(4,EHat,yRnd(3:10),MomExt(1:4,1:6),(/m_Stop,m_Stop,m_Stop,m_Stop/),PSWgt)! AStop, Stop
   call boost2Lab(eta1,eta2,6,MomExt(1:4,1:6))

! write(*,"(4F15.6)") momext(1:4,1)
! write(*,"(4F15.6)") momext(1:4,2)
! write(*,"(4F15.6)") momext(1:4,3)
! write(*,"(4F15.6)") momext(1:4,4)
! write(*,"(4F15.6)") momext(1:4,5)
! write(*,"(4F15.6)") momext(1:4,6)
! print *, momext(1:4,1).dot.momext(1:4,1)
! print *, momext(1:4,2).dot.momext(1:4,2)
! print *, momext(1:4,3).dot.momext(1:4,3)
! print *, momext(1:4,4).dot.momext(1:4,4)
! print *, momext(1:4,5).dot.momext(1:4,5)
! print *, momext(1:4,6).dot.momext(1:4,6)
! pause
!    call Kinematics_TTbarETmiss(.false.,MomExt,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
!    if( applyPSCut ) then
!       EvalCS_1L_ststbgggg = 0d0
!       return
!    endif

!    call SetPDFs(eta1,eta2,MuFac,pdf)
!    PDFFac_a = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
!             + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
!             + pdf(Bot_,1)*pdf(ABot_,2)
!    PDFFac_b = pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
!             + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
!             + pdf(Bot_,2)*pdf(ABot_,1)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   RunFactor = RunAlphaS(NLOParam,MuRen)
   nHel(1:2) = getHelicity(yrnd(17))
   PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))

!!!!! for of total cross section (fix y1,y2!)
!   beta = dsqrt(1d0-4d0*m_stop**2/EHat**2)
!   SigmaTot = (alpha_s*RunFactor)**2*DblPi/EHat**2*( 2d0/27d0 * beta**3  )
!   SigmaTot = SigmaTot * fbGeV2*SHatJacobi
!   PreFac = PreFac/SigmaTot
!   pdffac_a=0d0
!   pdffac_b=1d0
!!!!!!

  LO_Res_Unpol             = (0d0,0d0)
  NLO_Res_Unpol(-2:1)      = (0d0,0d0)
  NLO_Res_Unpol_Ferm(-2:1) = (0d0,0d0)
!   do npdf=1,2
!     if(npdf.eq.1) then
!         PDFFac = PDFFac_a
!     elseif(npdf.eq.2) then
!         PDFFac = PDFFac_b
!         call swapMom(MomExt(1:4,1),MomExt(1:4,2))
!     endif
    ISFac = MomCrossing(MomExt)
    call InitCurrCache()
    call SetPropagators()

!------------ LO ----------------
IF( Correction.EQ.0 ) THEN
   do iHel=nHel(1),nHel(2)
      call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,1),MomExt(1:4,5:9))
      call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,2),MomExt(1:4,10:14))
      call STopDecay(ExtParticle(5),DKX_STChi0_LO,Helicities(iHel,5),MomExt(1:4,5:9))
      call STopDecay(ExtParticle(6),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,10:14))

      call HelCrossing(Helicities(iHel,1:4))
      call SetPolarizations()

      do iPrimAmp=1,NumBornAmps
print *,"EVAL PRIMAMP",iPrimAmp
          call EvalTree(BornAmps(iPrimAmp))
print *, BornAmps(iPrimAmp)%Result
pause
      enddo
      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
! print *, "removed for cross check"
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol*PDFFac
pause

   enddo!helicity loop

!------------ 1 LOOP --------------
ELSEIF( Correction.EQ.1 ) THEN
    do iHel=nHel(1),nHel(2)
      call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,1),MomExt(1:4,5:9))
      call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,2),MomExt(1:4,10:14))

      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      do iPrimAmp=1,1
          call EvalTree(BornAmps(iPrimAmp))
      enddo
      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol*PDFFac

! ------------ bosonic loops --------------
      do iPrimAmp=1,1

print *, "hel",ihel,"primamp no",iPrimAmp
print *, "ordering ",PrimAmps(iPrimAmp)%ExtLine(1:5)
!           call SetKirill(PrimAmps(iPrimAmp))
!           call PentCut(PrimAmps(iPrimAmp))
!           call QuadCut(PrimAmps(iPrimAmp))
!           call TripCut(PrimAmps(iPrimAmp))
!           call DoubCut(PrimAmps(iPrimAmp))
!           call SingCut(PrimAmps(iPrimAmp))
!           call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)

!         call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen**2)
          PrimAmps(iPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(iPrimAmp)%Result(-2:1)
          call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,2,rdiv(2),rdiv(1))
          call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv)

print *, "----------------"
print *, "LO",BornAmps(iPrimAmp)%Result
print *, "NLO",PrimAmps(iPrimAmp)%Result(-2:1)
print *, "ratio",PrimAmps(iPrimAmp)%Result(-2)/BornAmps(iPrimAmp)%Result
print *, "ratio",PrimAmps(iPrimAmp)%Result(-1)/BornAmps(iPrimAmp)%Result

! print *, "check",PrimAmps(iPrimAmp)%UCuts(2)%Coeff(1:4,:)
! print *, "check",PrimAmps(iPrimAmp)%UCuts(1)%Coeff(1:1,:)

pause
      enddo



      enddo
ENDIF








!    EvalCS_1L_ststbgggg = 0d0
!    call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
!    if( EHat.le.4d0*m_STop ) then!   NOTE this is 4*mtop for check
!       EvalCS_1L_ststbgggg = 0d0
!       return
!    endif
!    FluxFac = 1d0/(2d0*EHat**2)
! 
!    call EvalPhasespace_2toN(4,EHat,yRnd(3:10),MomExt(1:4,1:6),(/m_Stop,m_Stop,m_Stop,m_Stop,0d0,0d0/),PSWgt)! AStop, Stop
!    call boost2Lab(eta1,eta2,6,MomExt(1:4,1:6))
! !    call EvalPhasespace_2toN(3,EHat,yRnd(3:7),MomExt(1:4,1:5),(/m_Stop,m_Stop,0d0,0d0,0d0/),PSWgt)! AStop, Stop
! !    call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))
! 
! ! !!!!!!!!!!!!!!!!!!!!! this is just for checking gauge invariance with gluon off the beam pipe
! print *, "Test Boost activated"
! MomBoost(1:4) = 1.8d0*MomExt(1:4,3)
! call boost(MomExt(1:4,1),MomBoost(1:4),1.8d0*m_sTop)
! call boost(MomExt(1:4,2),MomBoost(1:4),1.8d0*m_sTop)
! call boost(MomExt(1:4,3),MomBoost(1:4),1.8d0*m_sTop)
! call boost(MomExt(1:4,4),MomBoost(1:4),1.8d0*m_sTop)
! call boost(MomExt(1:4,5),MomBoost(1:4),1.8d0*m_sTop)
! call boost(MomExt(1:4,6),MomBoost(1:4),1.8d0*m_sTop)
! ! !!!!!!!!!!!!!!!!!!!!
!    ISFac = MomCrossing(MomExt)
! !    print *, MomExt(1:4,1).dot.MomExt(1:4,1)
! !    print *, MomExt(1:4,2).dot.MomExt(1:4,2)
! !    print *, MomExt(1:4,3).dot.MomExt(1:4,3)
! !    print *, MomExt(1:4,4).dot.MomExt(1:4,4)
! !    print *, MomExt(1:4,5).dot.MomExt(1:4,5)
! !    print *, MomExt(1:4,6).dot.MomExt(1:4,6)
! !    print *, MomExt(1,1)+MomExt(1,2)-MomExt(1,3)-MomExt(1,4)-MomExt(1,5)-MomExt(1,6)
! !    print *, MomExt(2,1)+MomExt(2,2)-MomExt(2,3)-MomExt(2,4)-MomExt(2,5)-MomExt(2,6)
! !    print *, MomExt(3,1)+MomExt(3,2)-MomExt(3,3)-MomExt(3,4)-MomExt(3,5)-MomExt(3,6)
! !    print *, MomExt(4,1)+MomExt(4,2)-MomExt(4,3)-MomExt(4,4)-MomExt(4,5)-MomExt(4,6)
! ! pause
! !    IF(XTOPDECAYS.EQ.3) THEN
! !       call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  Chi top
! !       call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,4),yRnd(7:8),MomExt(1:4,10:11),PSWgt3)
! !       call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd( 9:12),MomExt(1:4,7:9),PSWgt4)! bot lep neu
! !       call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(13:16),MomExt(1:4,12:14),PSWgt5)
! !       PSWgt = PSWgt * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
! !    ENDIF
!  
! !    call Kinematics_TTbarETmiss(.false.,MomExt,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
! !    if( applyPSCut ) then
! !       EvalCS_1L_ststbgggg = 0d0
! !       return
! !    endif
! 
!    call InitCurrCache()
!    call SetPropagators()
! 
!    call SetPDFs(eta1,eta2,MuFac,pdf)
!    PDFFac = pdf(0,1) * pdf(0,2)
!    PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
!    RunFactor = RunAlphaS(NLOParam,MuRen)
!    nHel(1:2) = getHelicity(yrnd(17))
!    PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))
! 
! !!!!! for of total cross section (fix y1,y2!)
! !   beta = dsqrt(1d0-4d0*m_stop**2/EHat**2)
! !   SigmaTot = (alpha_s*RunFactor)**2*DblPi/EHat**2*( beta*(5d0/48d0+31d0*m_stop**2/24d0/EHat**2) + &
! !                                                     dlog((1d0-beta)/(1d0+beta))*(2d0*m_stop**2/3d0/EHat**2+m_stop**4/6d0/Ehat**4)  )
! !   SigmaTot = SigmaTot * fbGeV2*SHatJacobi*PDFFac 
! !   PreFac = PreFac/SigmaTot
! !!!!!!
! 
!    LO_Res_Unpol = (0d0,0d0)
! !------------ LO ----------------
! IF( Correction.EQ.0 ) THEN
!    do iHel=nHel(1),nHel(2)
! !       call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExt(1:4,5:9))
! !       call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,10:14))
!       ExtParticle(1)%Pol(1:4) = (/1d0,0d0,0d0,0d0/)
!       ExtParticle(2)%Pol(1:4) = (/1d0,0d0,0d0,0d0/)
!       ExtParticle(3)%Pol(1:4) = (/1d0,0d0,0d0,0d0/)
!       ExtParticle(4)%Pol(1:4) = (/1d0,0d0,0d0,0d0/)
!       call HelCrossing(Helicities(iHel,1:6))
! !       call HelCrossing(Helicities(iHel,1:5))
!       call SetPolarizations()
!       do iPrimAmp=1,1!  NumBornAmps
!           call EvalTree(BornAmps(iPrimAmp))
!       enddo
!       print *, "Helicity",iHel
!       print *, "Result",BornAmps(1)%Result
! !       print *, "Result",BornAmps(2)%Result
!       pause
! 
!       LO_Res_Pol = (0d0,0d0)
!       do jPrimAmp=1,NumBornAmps
!       do iPrimAmp=1,NumBornAmps
!           LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
!       enddo
!       enddo
!       LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
!    enddo!helicity loop
! ENDIF
! 
! 
! IF( Correction.EQ.0 ) THEN
! !  normalization
!    LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * WidthExpansion
!    EvalCS_1L_ststbgggg = LO_Res_Unpol * PreFac
! ENDIF
! 
!    do NHisto=1,NumHistograms
!       call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_ststbgggg)
!    enddo
! 
!    EvalCS_1L_ststbgggg = EvalCS_1L_ststbgggg/VgsWgt




return
END FUNCTION







FUNCTION EvalCS_Real_ststbggg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModDipoles_GGSTSTBG
implicit none
real(8) :: EvalCS_Real_ststbggg,DipoleResult,yRnd(1:VegasMxDim),VgsWgt
complex(8) :: LO_Res_Pol,LO_Res_Unpol
integer :: iHel,jHel,kHel,iPrimAmp,jPrimAmp
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5,ISFac
real(8) :: MomExt(1:4,1:15),MomBoost(1:4)
logical :: applyPSCut
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac
real(8) :: pdf(-6:6,1:2),s12,s13
integer :: NBin(1:NumMaxHisto),NHisto,nHel(1:2)
!include 'misc/global_import'
include 'vegas_common.f'

! yrnd(1)= 0.1d0
! yrnd(2)= 0.22d0
! print *, "fixed yrnd"

   EvalCS_Real_ststbggg = 0d0
   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if( EHat.le.2d0*m_STop ) then
      EvalCS_Real_ststbggg = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to3Stops(EHat,yRnd(3:7),MomExt(1:4,1:5),PSWgt)! Glu AStop, Stop
   if( PSWgt.eq.0d0 )  then
      EvalCS_Real_ststbggg = 0d0
      SkipCounter = SkipCounter + 1
      return
   endif
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))
! !!!!!!!!!!!!!!!!!!!!! this is just for checking gauge invariance with gluon off the beam pipe
! print *, "Test Boost activated"
! MomBoost(1:4) = 1.8d0*MomExt(1:4,3)
! call boost(MomExt(1:4,1),MomBoost(1:4),1.8d0*m_sTop)
! call boost(MomExt(1:4,2),MomBoost(1:4),1.8d0*m_sTop)
! call boost(MomExt(1:4,3),MomBoost(1:4),1.8d0*m_sTop)
! call boost(MomExt(1:4,4),MomBoost(1:4),1.8d0*m_sTop)
! !!!!!!!!!!!!!!!!!!!!
   ISFac = MomCrossing(MomExt)

   IF(XTOPDECAYS.EQ.3) THEN
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,4),yRnd(8:9),MomExt(1:4,6:7),PSWgt2)!  Chi top
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,5),yRnd(10:11),MomExt(1:4,11:12),PSWgt3)
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,7),yRnd(12:15),MomExt(1:4,8:10),PSWgt4)! bot lep neu
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,12),yRnd(16:19),MomExt(1:4,13:15),PSWgt5)
      PSWgt = PSWgt * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
   ENDIF
 
   call SetPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   RunFactor = RunAlphaS(2,MuRen)

   call Kinematics_TTbarETmiss(.true.,MomExt,(/4,5,6,11,7,12,8,9,10,13,14,15,3/),applyPSCut,NBin)

! applyPScut = .true.! this is for inverted check

if(  applyPSCut ) then
   EvalCS_Real_ststbggg = 0d0
else
   call InitCurrCache()
   call SetPropagators()
   PreFac = PreFac
   LO_Res_Unpol = (0d0,0d0)
   do iHel=1,NumHelicities
      call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,6:10))
      call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,7),MomExt(1:4,11:15))
      call HelCrossing(Helicities(iHel,1:5))
      call SetPolarizations()
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo

      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbggg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
   enddo!helicity loop


!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**3
   EvalCS_Real_ststbggg = LO_Res_Unpol * PreFac

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_Real_ststbggg)
   enddo
   EvalCounter = EvalCounter + 1

! !     gg->ttbg
!       MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!       MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!       MG_MOM(0:3,3) = MomExt(1:4,5)*100d0
!       MG_MOM(0:3,4) = MomExt(1:4,4)*100d0
!       MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!       call coupsm(0)
!       call SGG_TTBG(MG_MOM,MadGraph_tree)
!       print *, ""
!       print *, "Unpolarized LO result:"
!       print *, "My tree:          ", LO_Res_Unpol
!       print *, "MadGraph hel.amp:", MadGraph_tree*(100d0)**2
!       print *, "ratio: ", MadGraph_tree*(100d0)**2/dble((LO_Res_Unpol))
!       pause
endif! applyPSCut


! IF( TopDecays.GE.1 ) THEN
     call EvalDipoles_GGSTSTBG((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_STop**2,m_STop**2,0d0/),yRnd(8:19),PreFac,DipoleResult)
! ELSE
!      call EvalDipoles_GGSTSTBG_noDK((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:19),PreFac,DipoleResult)
! ENDIF


!   s12 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,2))
!   s13 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,3))
!   write(* ,"(1PE23.16,3X,1PE23.16,3X,1PE23.16,3X,1PE23.16)") s13/s12,EvalCS_Real_ststbggg,DipoleResult, (EvalCS_Real_ststbggg/(-DipoleResult)-1d0)
!   pause


   if( IsNan(EvalCS_Real_ststbggg) .or. IsNan(DipoleResult) ) then
        print *, "NAN:",EvalCS_Real_ststbggg,DipoleResult
        print *, yRnd(:)
        print *, PSWgt , VgsWgt , PDFFac, sHatJacobi
        print *, eta1,eta2,MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,1)
        print *, MomExt(1:4,2)
        print *, MomExt(1:4,3)
        print *, MomExt(1:4,4)
        print *, MomExt(1:4,5)
        stop
   endif

        EvalCS_Real_ststbggg = (EvalCS_Real_ststbggg + DipoleResult)/VgsWgt



return
END FUNCTION









FUNCTION EvalCS_Real_ststbqqbg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModDipoles_QQBSTSTBG
use ModDipoles_QGSTSTBQ
implicit none
real(8) ::  EvalCS_Real_ststbqqbg,EvalCS_Dips_ttbqqbgp,DipoleResult,yRnd(1:VegasMxDim),VgsWgt
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol
integer :: iHel,jHel,kHel,iPrimAmp,jPrimAmp
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5,ISFac,PreFacDip
real(8) :: MomExt(1:4,1:15)
logical :: applyPSCut
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac_a,PDFFac_b,PDFFac
real(8) :: pdf(-6:6,1:2),s12,s13
integer :: NBin(1:NumMaxHisto),NHisto,nHel(1:2),npdf
real(8),parameter :: Nc=3d0
include 'vegas_common.f'



   EvalCS_Real_ststbqqbg = 0d0
   EvalCS_Dips_ttbqqbgp = 0d0
   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if( EHat.le.2d0*m_STop ) then
      EvalCS_Real_ststbqqbg = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to3Stops(EHat,yRnd(3:7),MomExt(1:4,1:5),PSWgt)! Glu AStop, Stop
   if( PSWgt.eq.0d0 )  then
      EvalCS_Real_ststbqqbg = 0d0
      SkipCounter = SkipCounter + 1
      return
   endif
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))

   IF(XTOPDECAYS.EQ.3) THEN
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,4),yRnd(8:9),MomExt(1:4,6:7),PSWgt2)!  Chi top
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,5),yRnd(10:11),MomExt(1:4,11:12),PSWgt3)
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,7),yRnd(12:15),MomExt(1:4,8:10),PSWgt4)! bot lep neu
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,12),yRnd(16:19),MomExt(1:4,13:15),PSWgt5)
      PSWgt = PSWgt * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
   ENDIF

   call SetPDFs(eta1,eta2,MuFac,pdf)
   IF( PROCESS.EQ.56 ) THEN
          PDFFac_a = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
                   + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
                   + pdf(Bot_,1)*pdf(ABot_,2)
          PDFFac_b = pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
                   + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
                   + pdf(Bot_,2)*pdf(ABot_,1)
    ELSEIF( PROCESS.EQ.53 ) THEN
          PDFFac_a = pdf(Up_,1) *pdf(0,2) + pdf(Dn_,1) *pdf(0,2)  &
                   + pdf(Chm_,1)*pdf(0,2) + pdf(Str_,1)*pdf(0,2)  &
                   + pdf(Bot_,1)*pdf(0,2)
          PDFFac_b = pdf(Up_,2) *pdf(0,1) + pdf(Dn_,2) *pdf(0,1)  &
                   + pdf(Chm_,2)*pdf(0,1) + pdf(Str_,2)*pdf(0,1)  &
                   + pdf(Bot_,2)*pdf(0,1)
    ELSEIF( PROCESS.EQ.54 ) THEN
          PDFFac_a = pdf(AUp_,2) *pdf(0,1) + pdf(ADn_,2) *pdf(0,1)  &
                   + pdf(AChm_,2)*pdf(0,1) + pdf(AStr_,2)*pdf(0,1)  &
                   + pdf(ABot_,2)*pdf(0,1)
          PDFFac_b = pdf(AUp_,1) *pdf(0,2) + pdf(ADn_,1) *pdf(0,2)  &
                   + pdf(AChm_,1)*pdf(0,2) + pdf(AStr_,1)*pdf(0,2)  &
                   + pdf(ABot_,1)*pdf(0,2)
    ENDIF

DO NPDF=1,2
    call Kinematics_TTbarETmiss(.true.,MomExt,(/4,5,6,11,7,12,8,9,10,13,14,15,3/),applyPSCut,NBin)

! applyPScut = .true.! this is for inverted check

    if(npdf.eq.1) then
        PDFFac = PDFFac_a
    elseif(npdf.eq.2) then
        PDFFac = PDFFac_b
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
    endif
    PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
    RunFactor = RunAlphaS(2,MuRen)
    ISFac = MomCrossing(MomExt)
    call InitCurrCache()
    call SetPropagators()


if(  applyPSCut ) then
    EvalCS_Real_ststbqqbg = 0d0
else


!------------ LO ----------------
   LO_Res_Unpol = (0d0,0d0)
    do iHel=1,NumHelicities
      call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,6:10))
      call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,7),MomExt(1:4,11:15))
      call HelCrossing(Helicities(iHel,1:5))
      call SetPolarizations()
      
       do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
!           print *, "amp",BornAmps(iPrimAmp)%ExtLine
!           print *, "res",BornAmps(iPrimAmp)%Result
!           pause
      enddo

      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqbg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
   enddo!helicity loop

!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**3 * PreFac *  PDFFac 
   EvalCS_Real_ststbqqbg = EvalCS_Real_ststbqqbg + LO_Res_Unpol

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol))
   enddo
   EvalCounter = EvalCounter + 1
endif! applyPSCut



  PreFacDip = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
  IF( PROCESS.EQ.56 ) THEN
      call EvalDipoles_QQBSTSTBG((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_STop**2,m_STop**2,0d0/),yRnd(8:19),PreFacDip,DipoleResult)
  ELSEIF( PROCESS.EQ.53 ) THEN
      call EvalDipoles_QGSTSTBQ((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_STop**2,m_STop**2,0d0/),yRnd(8:19),PreFacDip,DipoleResult)
  ELSEIF( PROCESS.EQ.54 ) THEN
      call EvalDipoles_QGSTSTBQ((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_STop**2,m_STop**2,0d0/),yRnd(8:19),PreFacDip,DipoleResult)
  ENDIF
  EvalCS_Dips_ttbqqbgp = EvalCS_Dips_ttbqqbgp + DipoleResult

ENDDO! npdf loop
! call swapMom(MomExt(1:4,1),MomExt(1:4,2))   ! swap back, not necessary here

!  s12 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,2))
!  s13 = 2d0*(MomExt(1:4,2).dot.MomExt(1:4,3))
!  write(* ,"(1PE23.16,3X,1PE23.16,3X,1PE23.16,3X,1PE23.16)") s13/s12,EvalCS_Real_ststbqqbg,EvalCS_Dips_ttbqqbgp, (EvalCS_Real_ststbqqbg/(-EvalCS_Dips_ttbqqbgp)- 1d0)
!  pause

   if( IsNan(EvalCS_Real_ststbqqbg) .or. IsNan(EvalCS_Dips_ttbqqbgp) ) then
        print *, "NAN:",EvalCS_Real_ststbqqbg,EvalCS_Dips_ttbqqbgp
        print *, yRnd(:)
        print *, PSWgt , VgsWgt , PDFFac, sHatJacobi
        print *, eta1,eta2,MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,1)
        print *, MomExt(1:4,2)
        print *, MomExt(1:4,3)
        print *, MomExt(1:4,4)
        print *, MomExt(1:4,5)
        stop
   endif

   EvalCS_Real_ststbqqbg = (EvalCS_Real_ststbqqbg + EvalCS_Dips_ttbqqbgp)/VgsWgt

END FUNCTION







! should be named DKReal
FUNCTION EvalCS_DKJ_Real_ststbgg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes
use ModParameters
use ModExoticDecay
use ModHadrWDecay
use ModMisc
implicit none
real(8) ::  EvalCS_DKJ_Real_ststbgg,EvalCS_DKJ_Real_ststbgg_1,EvalCS_DKJ_Real_ststbgg_2
real(8) :: yRnd(1:VegasMxDim),VgsWgt
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol
integer :: iHel,jHel,kHel,iPrimAmp,jPrimAmp
real(8) :: EHat,RunFactor,PSWgt,PSWgt0,PSWgt2,PSWgt3,PSWgt4,PSWgt5,ISFac
real(8) :: MomExt(1:4,1:15),MomExtTd(1:4,1:15),MomExtTdIn(1:4,1:4),MomBoost(1:4)
logical :: applyPSCut
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,TheDipole,DipoleResult
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2), mydummy(1:10)
integer :: NBin(1:NumMaxHisto),NHisto,nHel(1:2)
integer :: nJetRad,nJetRad1,nJetRad2,nJetRad3,nJetRad4,GluHel,ndip
real(8) :: pbDpg,ptDpg,ptDpb,omz,rsq,z,y
real(8),parameter :: CF=4d0/3d0
include 'vegas_common.f'



   EvalCS_DKJ_Real_ststbgg   = 0d0
   EvalCS_DKJ_Real_ststbgg_1 = 0d0
   EvalCS_DKJ_Real_ststbgg_2 = 0d0

! yRnd( 1)=  0.3385585941088194d0
! yRnd( 2)=  0.2799513116385563d0
! yRnd( 3)=  0.012473622342792d0
! yRnd( 4)=  0.2879364093709448d0
! yRnd( 5)=  0.1334328211068331d0
! yRnd( 6)=  0.7829718273519412d0
! yRnd( 7)=  0.3479862101366653d0
! yRnd( 8)=  0.1332233664734401d0
! yRnd( 9)=  0.2332185946559626d0
! yRnd(10)=  0.7471774192280964d0
! yRnd(11)=  0.438860277849650d0
! yRnd(12)=  0.02669250668338154d0
! yRnd(13)=  0.1692973097643342d0
! yRnd(14)=  0.6054243425398250d0
! yRnd(15)=  0.1081832926292833d0
! yRnd(16)=  0.1284279972912875d0
! yRnd(17)=  0.4608715537287632d0
! yRnd(18)=  0.3695206294159972d0
! yRnd(19)=  0.4436338681465174d0
! print *, "fixing yrnds"


   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if( EHat.le.2d0*m_STop ) then
      EvalCS_DKJ_Real_ststbgg = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to2Stops(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt0)! AStop, Stop
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))

   ISFac = MomCrossing(MomExt)
   call SetPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
   RunFactor = RunAlphaS(NLOParam,MuRen)


   nJetRad1=1
   nJetRad3=1
   if( TOPDECAYS.eq.1 ) then
      nJetRad2=2;  nJetRad4=2
   elseif( TOPDECAYS.eq.2 ) then
      nJetRad2=3;  nJetRad4=3
   elseif( TOPDECAYS.eq.3 ) then
      nJetRad2=2;  nJetRad4=3
   elseif( TOPDECAYS.eq.4 ) then
      nJetRad2=3;  nJetRad4=2
   endif

! ------------------- for checks --------------------------
!      nJetRad1=3
!      nJetRad2=3
!      nJetRad3=99
!      nJetRad4=2
! ---------------------------------------------



!----------------------------------
!  gluon  emission off anti-stop   |
!----------------------------------
do nJetRad=nJetRad1,nJetRad2!   nJetRad=1: gluon radiation off stop line,nJetRad=2: gluon radiation off top line,nJetRad=3: gluon radiation off W line

   EvalCS_DKJ_Real_ststbgg = 0d0
   IF(XTOPDECAYS.EQ.3) THEN
   if( nJetRad.eq.1 ) then
      call EvalPhasespace_StopDK(ST_Chi0_T_G,MomExt(1:4,3),yRnd(5:9),MomExt(1:4,5:7),PSWgt2)!  Chi(5) top(6) glu(7->15)
      MomExt(1:4,15) = MomExt(1:4,7)
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd(10:13),MomExt(1:4,7:9),PSWgt3)!  b(7) el(8) nu(9)
   elseif( nJetRad.eq.2 ) then
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  Chi(5) top(6)
      call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,6),yRnd(7:13),MomExt(1:4,7:10),PSWgt3)!  b(7) el(8) nu(9) glu(10->15)
      MomExt(1:4,15) = MomExt(1:4,10)
   elseif( nJetRad.eq.3 ) then
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  Chi(5) top(6)
      call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,6),yRnd(7:13),MomExt(1:4,7:10),PSWgt3)!  b(7) el(8) nu(9) glu(10->15)
      MomExt(1:4,15) = MomExt(1:4,10)
   endif
   call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,4),yRnd(14:15),MomExt(1:4,10:11),PSWgt4)! Chi(10) Top(11)
   call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(16:19),MomExt(1:4,12:14),PSWgt5)! b(12) el(13) nu(14)
   ENDIF
   call Kinematics_TTbarETmiss(.true.,MomExt,(/3,4,5,10,6,11,7,8,9,12,13,14,15/),applyPSCut,NBin)

   PSWgt = PSWgt0 * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   if( PreFac.eq.0d0 ) then
        SkipCounter = SkipCounter + 1
        cycle
   endif
   if( applyPSCut ) then
      EvalCS_DKJ_Real_ststbgg = 0d0
      PSCutCounter = PSCutCounter + 1
      goto 2000
   endif


   LO_Res_Unpol = (0d0,0d0)
   do iHel=1,NumHelicities
   do GluHel=1,-1,-2
      call HelCrossing(Helicities(iHel,1:6))
      call SetPolarizations()
      if( nJetRad.eq.1 ) then
          call STopDecay(ExtParticle(1),DKX_STChi0_RE1,Helicities(iHel,5),MomExt(1:4,5:9),MomExt(1:4,15),GluHel)
      elseif( nJetRad.eq.2 ) then
          call STopDecay(ExtParticle(1),DKX_STChi0_RE2,Helicities(iHel,5),MomExt(1:4,5:9),MomExt(1:4,15),GluHel)
      else
          call STopDecay(ExtParticle(1),DKX_STChi0_RE3,Helicities(iHel,5),MomExt(1:4,5:9),MomExt(1:4,15),GluHel)
      endif
      call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,10:14))

      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo

      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol

   enddo!helicity loop
   enddo!helicity loop


   EvalCS_DKJ_Real_ststbgg = LO_Res_UnPol * PreFac * ISFac * (alpha_s4Pi*RunFactor)**2
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_DKJ_Real_ststbgg)
   enddo
   EvalCounter = EvalCounter + 1

   EvalCS_DKJ_Real_ststbgg_1 = EvalCS_DKJ_Real_ststbgg_1 + EvalCS_DKJ_Real_ststbgg

mydummy(1) = EvalCS_DKJ_Real_ststbgg


2000 continue!! dipoles for gluon emission off anti-stop



   if( nJetRad.eq.1 ) then

        call WTransform3(MomExt(1:4,5:6),MomExt(1:4,15),MomExtTd(1:4,5:6),pbDpg,ptDpg,ptDpb)
        call EvalPhasespace_TopDK(T_B_W,MomExtTd(1:4,6),yRnd(10:13),MomExtTd(1:4,7:9),PSWgt3)
        MomExtTd(1:4,1:4) = MomExt(1:4,1:4)
        MomExtTd(1:4,10:14) = MomExt(1:4,10:14)
        
        TheDipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg/ptDpg*( m_Stop**2+m_Top**2-(MomExt(1:4,5).dot.MomExt(1:4,5)) ) - (m_STop/ptDpg)**2 - (m_Top/pbDpg)**2 )

        call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
        if( applyPSCut .or. TheDipole.eq.0d0 ) then
            cycle
        endif

        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
            call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExtTd(1:4,5:9))
            call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExtTd(1:4,10:14))
            call HelCrossing(Helicities(iHel,1:4))
            call SetPolarizations()
            do iPrimAmp=1,NumBornAmps
                call EvalTree(BornAmps(iPrimAmp))
            enddo

            LO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
            do iPrimAmp=1,NumBornAmps
                LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
            enddo
            enddo
            LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop
        DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

        do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),DipoleResult)
        enddo
        EvalCS_DKJ_Real_ststbgg_1 = EvalCS_DKJ_Real_ststbgg_1 + DipoleResult

! mydummy(2) = DipoleResult
! print *, "E(glu)/E(stop)",MomExt(1,15)/MomExt(1,5)
! print *, "real    result",mydummy(1)
! print *, "dipole  result",mydummy(2)
! print *, "ratios        ",mydummy(1)/mydummy(2), mydummy(1)/mydummy(2)+1d0
! pause
! return

   elseif( nJetRad.eq.2 ) then

        MomExtTd(1:4,1:6) = MomExt(1:4,1:6)
        MomExtTd(1:4,10:14) = MomExt(1:4,10:14)
        MomExtTdIn(1:4,1:3) = MomExt(1:4,7:9)
        MomExtTdIn(1:4,4)  = MomExt(1:4,15)
        call WTransform(MomExtTdIn(1:4,1:4),MomExtTd(1:4,7:9),pbDpg,ptDpg,ptDpb)
        
        omz=ptDpg/(ptDpb+ptDpg-pbDpg)  !  for some reason this is not (1-z) as defined in the paper...
        rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
        z=1d0-omz
        y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2

        TheDipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
        TheDipole = TheDipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )


        call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
        if( applyPSCut .or. TheDipole.eq.0d0) then
            cycle
        endif


        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
            call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExtTd(1:4,5:9))
            call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExtTd(1:4,10:14))
            call HelCrossing(Helicities(iHel,1:4))
            call SetPolarizations()
            do iPrimAmp=1,NumBornAmps
                call EvalTree(BornAmps(iPrimAmp))
            enddo

            LO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
            do iPrimAmp=1,NumBornAmps
                LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
            enddo
            enddo
            LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop
        DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

        do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),DipoleResult)
        enddo
        EvalCS_DKJ_Real_ststbgg_1 = EvalCS_DKJ_Real_ststbgg_1 + DipoleResult


! mydummy(2) = DipoleResult
! print *, "sij/mtop",(MomExt(1:4,15).dot.MomExt(1:4,7))/m_top
! print *, "E(glu)/E(stop)",MomExt(1,15)/MomExt(1,7)
! print *, "real    result",mydummy(1)
! print *, "dipole  result",mydummy(2)
! print *, "ratios        ",mydummy(1)/mydummy(2), mydummy(1)/mydummy(2)+1d0
! pause
! return

   elseif( nJetRad.eq.3 ) then
          MomExtTd(1:4,1:6) = MomExt(1:4,1:6)
          MomExtTd(1:4,10:14) = MomExt(1:4,10:14)
          MomExtTdIn(1:4,1:3) = MomExt(1:4,7:9)
          MomExtTdIn(1:4,4) = MomExt(1:4,15)
          do ndip=1,2
              call wdec_trans(ndip,MomExtTdIn(1:4,1:4),MomExtTd(1:4,7:9),alpha_DKWff,TheDipole)
              TheDipole = - alpha_s4Pi*RunFactor * CF * TheDipole

              call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
              if( applyPSCut .or. TheDipole.eq.0d0 ) then
                  cycle
              endif

              LO_Res_Unpol = (0d0,0d0)
              do iHel=1,NumHelicities! helicity summation
                  call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExtTd(1:4,5:9))
                  call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExtTd(1:4,10:14))
                  call HelCrossing(Helicities(iHel,1:4))
                  call SetPolarizations()
                  do iPrimAmp=1,NumBornAmps
                      call EvalTree(BornAmps(iPrimAmp))
                  enddo

                  LO_Res_Pol = (0d0,0d0)
                  do jPrimAmp=1,NumBornAmps
                  do iPrimAmp=1,NumBornAmps
                      LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
                  enddo
                  enddo
                  LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
              enddo!helicity loop
              DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

              do NHisto=1,NumHistograms
                  call intoHisto(NHisto,NBin(NHisto),DipoleResult)
              enddo
              EvalCS_DKJ_Real_ststbgg_1 = EvalCS_DKJ_Real_ststbgg_1 + DipoleResult
mydummy(2) = mydummy(2) + DipoleResult
          enddo   !dipole loop

! print *, "sij/shat",(MomExt(1:4,15).dot.MomExt(1:4,7))/EHat**2
! print *, "E(glu)/E(stop)",MomExt(1,15)/MomExt(1,7)
! print *, "real    result",mydummy(1)
! print *, "dipole  result",mydummy(2)
! print *, "ratios        ",mydummy(1)/mydummy(2), mydummy(1)/mydummy(2)+1d0
! pause
! return

  endif! nJetRad for dipoles

enddo! nJetRad loop







!----------------------------------
!  gluon  emission off stop   |
!----------------------------------
do nJetRad=nJetRad3,nJetRad4!   nJetRad=1: gluon radiation off stop line,nJetRad=2: gluon radiation off top line,nJetRad=3: gluon radiation off W line

   EvalCS_DKJ_Real_ststbgg = 0d0
   IF(XTOPDECAYS.EQ.3) THEN
   if( nJetRad.eq.1 ) then
      call EvalPhasespace_StopDK(ST_Chi0_T_G,MomExt(1:4,4),yRnd(11:15),MomExt(1:4,10:12),PSWgt2)!  Chi(10) top(11) glu(12->15)
      MomExt(1:4,15) = MomExt(1:4,12)
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(16:19),MomExt(1:4,12:14),PSWgt3)!  b(12) el(13) nu(14)
   elseif( nJetRad.eq.2 ) then
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,4),yRnd(11:12),MomExt(1:4,10:11),PSWgt2)!  Chi(10) top(11)
      call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,11),yRnd(13:19),MomExt(1:4,12:15),PSWgt3)!  b(12) el(13) nu(14) glu(15)
   elseif( nJetRad.eq.3 ) then
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,4),yRnd(11:12),MomExt(1:4,10:11),PSWgt2)!  Chi(10) top(11)
      call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,11),yRnd(13:19),MomExt(1:4,12:15),PSWgt3)!  b(12) el(13) nu(14) glu(15)
   endif   
   call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt4)! Chi(5) Top(6)
   call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd(7:10),MomExt(1:4,7:9),PSWgt5)! b(7) el(8) nu(9)
   ENDIF
   call Kinematics_TTbarETmiss(.true.,MomExt,(/3,4,5,10,6,11,7,8,9,12,13,14,15/),applyPSCut,NBin)
   PSWgt = PSWgt0 * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   if( PreFac.eq.0d0 ) then
      SkipCounter = SkipCounter + 1
      cycle
   endif
   if( applyPSCut ) then
      EvalCS_DKJ_Real_ststbgg = 0d0
      PSCutCounter = PSCutCounter + 1
      goto 2001
   endif


   LO_Res_Unpol = (0d0,0d0)
   do iHel=1,NumHelicities
   do GluHel=1,-1,-2

      call HelCrossing(Helicities(iHel,1:6))
      call SetPolarizations()
      call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExt(1:4,5:9))
      if( nJetRad.eq.1 ) then
          call STopDecay(ExtParticle(2),DKX_STChi0_RE1,Helicities(iHel,6),MomExt(1:4,10:14),MomExt(1:4,15),GluHel)
      elseif( nJetRad.eq.2 ) then
          call STopDecay(ExtParticle(2),DKX_STChi0_RE2,Helicities(iHel,6),MomExt(1:4,10:14),MomExt(1:4,15),GluHel)
      else
          call STopDecay(ExtParticle(2),DKX_STChi0_RE3,Helicities(iHel,6),MomExt(1:4,10:14),MomExt(1:4,15),GluHel)
      endif

      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo

      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol

   enddo!helicity loop
   enddo!helicity loop

   EvalCS_DKJ_Real_ststbgg = LO_Res_UnPol * PreFac * ISFac * (alpha_s4Pi*RunFactor)**2
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_DKJ_Real_ststbgg)
   enddo
   EvalCounter = EvalCounter + 1

   EvalCS_DKJ_Real_ststbgg_2 = EvalCS_DKJ_Real_ststbgg_2 + EvalCS_DKJ_Real_ststbgg

mydummy(3) = EvalCS_DKJ_Real_ststbgg


2001 continue!! dipoles for gluon emission off stop




   if( nJetRad.eq.1 ) then

        call WTransform3(MomExt(1:4,10:11),MomExt(1:4,15),MomExtTd(1:4,10:11),pbDpg,ptDpg,ptDpb)
        call EvalPhasespace_TopDK(T_B_W,MomExtTd(1:4,11),yRnd(16:19),MomExtTd(1:4,12:14),PSWgt3)
        MomExtTd(1:4,1:9) = MomExt(1:4,1:9)

        TheDipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg/ptDpg*( m_Stop**2+m_Top**2-(MomExt(1:4,10).dot.MomExt(1:4,10)) ) - (m_STop/ptDpg)**2 - (m_Top/pbDpg)**2 )

        call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
        if( applyPSCut .or. TheDipole.eq.0d0 ) then
            cycle
        endif


        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
            call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExtTd(1:4,5:9))
            call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExtTd(1:4,10:14))
            call HelCrossing(Helicities(iHel,1:4))
            call SetPolarizations()
            do iPrimAmp=1,NumBornAmps
                call EvalTree(BornAmps(iPrimAmp))
            enddo

            LO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
            do iPrimAmp=1,NumBornAmps
                LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
            enddo
            enddo
            LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop
        DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

        do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),DipoleResult)
        enddo
        EvalCS_DKJ_Real_ststbgg_2 = EvalCS_DKJ_Real_ststbgg_2 + DipoleResult


! mydummy(4) = DipoleResult
! print *, "E(glu)/E(stop)",MomExt(1,15)/MomExt(1,10)
! print *, "real    result",mydummy(3)
! print *, "dipole  result",mydummy(4)
! print *, "ratios        ",mydummy(3)/mydummy(4), mydummy(3)/mydummy(4)+1d0
! pause
! return


   elseif( nJetRad.eq.2 ) then

        MomExtTd(1:4,1:11) = MomExt(1:4,1:11)
        MomExtTdIn(1:4,1:4) = MomExt(1:4,12:15)
        call WTransform(MomExtTdIn(1:4,1:4),MomExtTd(1:4,12:14),pbDpg,ptDpg,ptDpb)
        
        omz=ptDpg/(ptDpb+ptDpg-pbDpg)  !  for some reason this is not (1-z) as defined in the paper...
        rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
        z=1d0-omz
        y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2

        TheDipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
        TheDipole = TheDipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )


        call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
        if( applyPSCut .or. TheDipole.eq.0d0) then
            cycle
        endif


        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
            call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExtTd(1:4,5:9))
            call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExtTd(1:4,10:14))
            call HelCrossing(Helicities(iHel,1:4))
            call SetPolarizations()
            do iPrimAmp=1,NumBornAmps
                call EvalTree(BornAmps(iPrimAmp))
            enddo

            LO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
            do iPrimAmp=1,NumBornAmps
                LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
            enddo
            enddo
            LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop
        DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

        do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),DipoleResult)
        enddo
        EvalCS_DKJ_Real_ststbgg_2 = EvalCS_DKJ_Real_ststbgg_2 + DipoleResult


! mydummy(4) = DipoleResult
! print *, "sij/shat",(MomExt(1:4,15).dot.MomExt(1:4,12))/EHat**2
! print *, "E(glu)/E(top)",MomExt(1,15)/MomExt(1,12)
! print *, "real    result",mydummy(3)
! print *, "dipole  result",mydummy(4)
! print *, "ratios        ",mydummy(3)/mydummy(4), mydummy(3)/mydummy(4)+1d0
! pause
! return

   elseif( nJetRad.eq.3 ) then
          MomExtTd(1:4,1:11) = MomExt(1:4,1:11)
          MomExtTdIn(1:4,1:4) = MomExt(1:4,12:15)
          do ndip=1,2
              call wdec_trans(ndip,MomExtTdIn(1:4,1:4),MomExtTd(1:4,12:14),alpha_DKWff,TheDipole)
              TheDipole = - alpha_s4Pi*RunFactor * CF * TheDipole

              call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
              if( applyPSCut .or. TheDipole.eq.0d0 ) then
                  cycle
              endif

              LO_Res_Unpol = (0d0,0d0)
              do iHel=1,NumHelicities! helicity summation
                  call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExtTd(1:4,5:9))
                  call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExtTd(1:4,10:14))
                  call HelCrossing(Helicities(iHel,1:4))
                  call SetPolarizations()
                  do iPrimAmp=1,NumBornAmps
                      call EvalTree(BornAmps(iPrimAmp))
                  enddo

                  LO_Res_Pol = (0d0,0d0)
                  do jPrimAmp=1,NumBornAmps
                  do iPrimAmp=1,NumBornAmps
                      LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
                  enddo
                  enddo
                  LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
              enddo!helicity loop
              DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

              do NHisto=1,NumHistograms
                  call intoHisto(NHisto,NBin(NHisto),DipoleResult)
              enddo
              EvalCS_DKJ_Real_ststbgg_2 = EvalCS_DKJ_Real_ststbgg_2 + DipoleResult
          enddo   !dipole loop
  endif! nJetRad for dipoles

enddo! nJetRad loop





   EvalCS_DKJ_Real_ststbgg = (EvalCS_DKJ_Real_ststbgg_1+EvalCS_DKJ_Real_ststbgg_2)/VgsWgt

   if( IsNan(EvalCS_DKJ_Real_ststbgg) ) then
        print *, "NAN:",EvalCS_DKJ_Real_ststbgg
        print *, yRnd(:)
        print *, PSWgt , VgsWgt , PDFFac, sHatJacobi
        print *, eta1,eta2,MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,1:11)
        print *, "SKIP EVENT!!!!!"
        EvalCS_DKJ_Real_ststbgg = 0d0
        return
   endif


return
END FUNCTION









FUNCTION EvalCS_DKJ_Real_ststbqqb(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes
use ModParameters
use ModExoticDecay
use ModHadrWDecay
use ModMisc
implicit none
real(8) ::  EvalCS_DKJ_Real_ststbqqb,EvalCS_DKJ_Real_ststbqqb_1,EvalCS_DKJ_Real_ststbqqb_2
real(8) :: yRnd(1:VegasMxDim),VgsWgt
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol
integer :: iHel,jHel,kHel,iPrimAmp,jPrimAmp
real(8) :: EHat,RunFactor,PSWgt,PSWgt0,PSWgt2,PSWgt3,PSWgt4,PSWgt5,ISFac
real(8) :: MomExt(1:4,1:15),MomExtTd(1:4,1:15),MomExtTdIn(1:4,1:4),MomBoost(1:4)
logical :: applyPSCut
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,TheDipole,DipoleResult,PDFFac_a,PDFFac_b
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2)
integer :: NBin(1:NumMaxHisto),NHisto,nHel(1:2)
integer :: nJetRad,nJetRad1,nJetRad2,nJetRad3,nJetRad4,GluHel,ndip
real(8) :: pbDpg,ptDpg,ptDpb,omz,rsq,z,y
real(8),parameter :: CF=4d0/3d0
include 'vegas_common.f'



   EvalCS_DKJ_Real_ststbqqb   = 0d0
   EvalCS_DKJ_Real_ststbqqb_1 = 0d0
   EvalCS_DKJ_Real_ststbqqb_2 = 0d0


   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if( EHat.le.2d0*m_STop ) then
      EvalCS_DKJ_Real_ststbqqb = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to2Stops(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt0)! AStop, Stop
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))

   ISFac = MomCrossing(MomExt)
   call SetPDFs(eta1,eta2,MuFac,pdf)
   PDFFac_a = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
            + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
            + pdf(Bot_,1)*pdf(ABot_,2)
   PDFFac_b = pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
            + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
            + pdf(Bot_,2)*pdf(ABot_,1)
   PDFFac = PDFFac_a + PDFFac_b
   RunFactor = RunAlphaS(NLOParam,MuRen)


   nJetRad1=1
   nJetRad3=1
   if( TOPDECAYS.eq.1 ) then
      nJetRad2=2;  nJetRad4=2
   elseif( TOPDECAYS.eq.2 ) then
      nJetRad2=3;  nJetRad4=3
   elseif( TOPDECAYS.eq.3 ) then
      nJetRad2=2;  nJetRad4=3
   elseif( TOPDECAYS.eq.4 ) then
      nJetRad2=3;  nJetRad4=2
   endif


! ------------------- for checks --------------------------
!    nJetRad1=1
!    nJetRad2=-1
!    nJetRad3=2
!    nJetRad4=2
! ---------------------------------------------






!-----------------------------------
!  gluon  emission off anti-stop   |
!----------------------------------
do nJetRad=nJetRad1,nJetRad2!   nJetRad=1: gluon radiation off stop line,nJetRad=2: gluon radiation off top line,nJetRad=3: gluon radiation off W line

   EvalCS_DKJ_Real_ststbqqb = 0d0
   IF(XTOPDECAYS.EQ.3) THEN
   if( nJetRad.eq.1 ) then
      call EvalPhasespace_StopDK(ST_Chi0_T_G,MomExt(1:4,3),yRnd(5:9),MomExt(1:4,5:7),PSWgt2)!  Chi(5) top(6) glu(7->15)
      MomExt(1:4,15) = MomExt(1:4,7)
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd(10:13),MomExt(1:4,7:9),PSWgt3)!  b(7) el(8) nu(9)
   elseif( nJetRad.eq.2 ) then
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  Chi(5) top(6)
      call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,6),yRnd(7:13),MomExt(1:4,7:10),PSWgt3)!  b(7) el(8) nu(9) glu(10->15)
      MomExt(1:4,15) = MomExt(1:4,10)
   elseif( nJetRad.eq.3 ) then
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  Chi(5) top(6)
      call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,6),yRnd(7:13),MomExt(1:4,7:10),PSWgt3)!  b(7) el(8) nu(9) glu(10->15)
      MomExt(1:4,15) = MomExt(1:4,10)
   endif   
   call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,4),yRnd(14:15),MomExt(1:4,10:11),PSWgt4)! Chi(10) Top(11)
   call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(16:19),MomExt(1:4,12:14),PSWgt5)! b(12) el(13) nu(14)
   ENDIF
   call Kinematics_TTbarETmiss(.true.,MomExt,(/3,4,5,10,6,11,7,8,9,12,13,14,15/),applyPSCut,NBin)
   PSWgt = PSWgt0 * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   if( PreFac.eq.0d0 ) then
        SkipCounter = SkipCounter + 1
        cycle
   endif
   if( applyPSCut ) then
      EvalCS_DKJ_Real_ststbqqb = 0d0
      PSCutCounter = PSCutCounter + 1
      goto 2000
   endif


   LO_Res_Unpol = (0d0,0d0)
   do iHel=1,NumHelicities
   do GluHel=1,-1,-2

      call HelCrossing(Helicities(iHel,1:6))
      call SetPolarizations()
      if( nJetRad.eq.1 ) then
          call STopDecay(ExtParticle(1),DKX_STChi0_RE1,Helicities(iHel,5),MomExt(1:4,5:9),MomExt(1:4,15),GluHel)
      elseif( nJetRad.eq.2 ) then
          call STopDecay(ExtParticle(1),DKX_STChi0_RE2,Helicities(iHel,5),MomExt(1:4,5:9),MomExt(1:4,15),GluHel)
      else
          call STopDecay(ExtParticle(1),DKX_STChi0_RE3,Helicities(iHel,5),MomExt(1:4,5:9),MomExt(1:4,15),GluHel)
      endif
      call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,10:14))

      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo

      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol

   enddo!helicity loop
   enddo!helicity loop

   EvalCS_DKJ_Real_ststbqqb = LO_Res_UnPol * PreFac * ISFac * (alpha_s4Pi*RunFactor)**2
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_DKJ_Real_ststbqqb)
   enddo
   EvalCounter = EvalCounter + 1

   EvalCS_DKJ_Real_ststbqqb_1 = EvalCS_DKJ_Real_ststbqqb_1 + EvalCS_DKJ_Real_ststbqqb




2000 continue!! dipoles for gluon emission off anti-stop



! print *, "check helicity arrays"


   if( nJetRad.eq.1 ) then
       
        call WTransform3(MomExt(1:4,5:6),MomExt(1:4,15),MomExtTd(1:4,5:6),pbDpg,ptDpg,ptDpb)
        call EvalPhasespace_TopDK(T_B_W,MomExtTd(1:4,6),yRnd(10:13),MomExtTd(1:4,7:9),PSWgt3)
        MomExtTd(1:4,1:4) = MomExt(1:4,1:4)
        MomExtTd(1:4,10:14) = MomExt(1:4,10:14)

        TheDipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg/ptDpg*( m_Stop**2+m_Top**2-(MomExt(1:4,5).dot.MomExt(1:4,5)) ) - (m_STop/ptDpg)**2 - (m_Top/pbDpg)**2 )

        call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
        if( applyPSCut .or. TheDipole.eq.0d0 ) then
            cycle
        endif



        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
            call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExtTd(1:4,5:9))
            call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExtTd(1:4,10:14))
            call HelCrossing(Helicities(iHel,1:4))
            call SetPolarizations()
            do iPrimAmp=1,NumBornAmps
                call EvalTree(BornAmps(iPrimAmp))
            enddo

            LO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
            do iPrimAmp=1,NumBornAmps
                LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
            enddo
            enddo
            LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop
        DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

        do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),DipoleResult)
        enddo
        EvalCS_DKJ_Real_ststbqqb_1 = EvalCS_DKJ_Real_ststbqqb_1 + DipoleResult


   elseif( nJetRad.eq.2 ) then

        MomExtTd(1:4,1:6) = MomExt(1:4,1:6)
        MomExtTd(1:4,10:14) = MomExt(1:4,10:14)
        MomExtTdIn(1:4,1:3) = MomExt(1:4,7:9)
        MomExtTdIn(1:4,4)  = MomExt(1:4,15)
        call WTransform(MomExtTdIn(1:4,1:4),MomExtTd(1:4,7:9),pbDpg,ptDpg,ptDpb)
        
        omz=ptDpg/(ptDpb+ptDpg-pbDpg)  !  for some reason this is not (1-z) as defined in the paper...
        rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
        z=1d0-omz
        y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2

        TheDipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
        TheDipole = TheDipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )


        call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
        if( applyPSCut .or. TheDipole.eq.0d0) then
            cycle
        endif


        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
            call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExtTd(1:4,5:9))
            call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExtTd(1:4,10:14))
            call HelCrossing(Helicities(iHel,1:4))
            call SetPolarizations()
            do iPrimAmp=1,NumBornAmps
                call EvalTree(BornAmps(iPrimAmp))
            enddo

            LO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
            do iPrimAmp=1,NumBornAmps
                LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
            enddo
            enddo
            LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop
        DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

        do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),DipoleResult)
        enddo
        EvalCS_DKJ_Real_ststbqqb_1 = EvalCS_DKJ_Real_ststbqqb_1 + DipoleResult


   elseif( nJetRad.eq.3 ) then
          MomExtTd(1:4,1:6) = MomExt(1:4,1:6)
          MomExtTd(1:4,10:14) = MomExt(1:4,10:14)
          MomExtTdIn(1:4,1:3) = MomExt(1:4,7:9)
          MomExtTdIn(1:4,4) = MomExt(1:4,15)
          do ndip=1,2
              call wdec_trans(ndip,MomExtTdIn(1:4,1:4),MomExtTd(1:4,7:9),alpha_DKWff,TheDipole)
              TheDipole = - alpha_s4Pi*RunFactor * CF * TheDipole

              call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
              if( applyPSCut .or. TheDipole.eq.0d0 ) then
                  cycle
              endif

              LO_Res_Unpol = (0d0,0d0)
              do iHel=1,NumHelicities! helicity summation
                  call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExtTd(1:4,5:9))
                  call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExtTd(1:4,10:14))
                  call HelCrossing(Helicities(iHel,1:4))
                  call SetPolarizations()
                  do iPrimAmp=1,NumBornAmps
                      call EvalTree(BornAmps(iPrimAmp))
                  enddo

                  LO_Res_Pol = (0d0,0d0)
                  do jPrimAmp=1,NumBornAmps
                  do iPrimAmp=1,NumBornAmps
                      LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
                  enddo
                  enddo
                  LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
              enddo!helicity loop
              DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

              do NHisto=1,NumHistograms
                  call intoHisto(NHisto,NBin(NHisto),DipoleResult)
              enddo
              EvalCS_DKJ_Real_ststbqqb_1 = EvalCS_DKJ_Real_ststbqqb_1 + DipoleResult
          enddo   !dipole loop
  endif! nJetRad for dipoles

enddo! nJetRad loop







!----------------------------------
!  gluon  emission off stop   |
!----------------------------------
do nJetRad=nJetRad3,nJetRad4!   nJetRad=1: gluon radiation off stop line,nJetRad=2: gluon radiation off top line,nJetRad=3: gluon radiation off W line

   EvalCS_DKJ_Real_ststbqqb = 0d0
   IF(XTOPDECAYS.EQ.3) THEN
   if( nJetRad.eq.1 ) then
      call EvalPhasespace_StopDK(ST_Chi0_T_G,MomExt(1:4,4),yRnd(11:15),MomExt(1:4,10:12),PSWgt2)!  Chi(10) top(11) glu(12->15)
      MomExt(1:4,15) = MomExt(1:4,12)
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(16:19),MomExt(1:4,12:14),PSWgt3)!  b(12) el(13) nu(14)
   elseif( nJetRad.eq.2 ) then
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,4),yRnd(11:12),MomExt(1:4,10:11),PSWgt2)!  Chi(10) top(11)
      call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,11),yRnd(13:19),MomExt(1:4,12:15),PSWgt3)!  b(12) el(13) nu(14) glu(15)
   elseif( nJetRad.eq.3 ) then
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,4),yRnd(11:12),MomExt(1:4,10:11),PSWgt2)!  Chi(10) top(11)
      call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,11),yRnd(13:19),MomExt(1:4,12:15),PSWgt3)!  b(12) el(13) nu(14) glu(15)
   endif   
   call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt4)! Chi(5) Top(6)
   call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd(7:10),MomExt(1:4,7:9),PSWgt5)! b(7) el(8) nu(9)
   ENDIF
   call Kinematics_TTbarETmiss(.true.,MomExt,(/3,4,5,10,6,11,7,8,9,12,13,14,15/),applyPSCut,NBin)
   PSWgt = PSWgt0 * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   if( PreFac.eq.0d0 ) then
        SkipCounter = SkipCounter + 1
        cycle
   endif
   if( applyPSCut ) then
      EvalCS_DKJ_Real_ststbqqb = 0d0
      PSCutCounter = PSCutCounter + 1
      goto 2001
   endif


   LO_Res_Unpol = (0d0,0d0)
   do iHel=1,NumHelicities
   do GluHel=1,-1,-2

      call HelCrossing(Helicities(iHel,1:6))
      call SetPolarizations()
      call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExt(1:4,5:9))
      if( nJetRad.eq.1 ) then
          call STopDecay(ExtParticle(2),DKX_STChi0_RE1,Helicities(iHel,6),MomExt(1:4,10:14),MomExt(1:4,15),GluHel)
      elseif( nJetRad.eq.2 ) then
          call STopDecay(ExtParticle(2),DKX_STChi0_RE2,Helicities(iHel,6),MomExt(1:4,10:14),MomExt(1:4,15),GluHel)
      else
          call STopDecay(ExtParticle(2),DKX_STChi0_RE3,Helicities(iHel,6),MomExt(1:4,10:14),MomExt(1:4,15),GluHel)
      endif

      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo

      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol

   enddo!helicity loop
   enddo!helicity loop

   EvalCS_DKJ_Real_ststbqqb = LO_Res_UnPol * PreFac * ISFac * (alpha_s4Pi*RunFactor)**2
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_DKJ_Real_ststbqqb)
   enddo
   EvalCounter = EvalCounter + 1

   EvalCS_DKJ_Real_ststbqqb_2 = EvalCS_DKJ_Real_ststbqqb_2 + EvalCS_DKJ_Real_ststbqqb




2001 continue!! dipoles for gluon emission off stop




   if( nJetRad.eq.1 ) then
        call WTransform3(MomExt(1:4,10:11),MomExt(1:4,15),MomExtTd(1:4,10:11),pbDpg,ptDpg,ptDpb)
        call EvalPhasespace_TopDK(T_B_W,MomExtTd(1:4,11),yRnd(16:19),MomExtTd(1:4,12:14),PSWgt3)
        MomExtTd(1:4,1:9) = MomExt(1:4,1:9)

        TheDipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg/ptDpg*( m_Stop**2+m_Top**2-(MomExt(1:4,10).dot.MomExt(1:4,10)) ) - (m_STop/ptDpg)**2 - (m_Top/pbDpg)**2 )

        call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
        if( applyPSCut .or. TheDipole.eq.0d0 ) then
            cycle
        endif


        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
            call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExtTd(1:4,5:9))
            call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExtTd(1:4,10:14))
            call HelCrossing(Helicities(iHel,1:4))
            call SetPolarizations()
            do iPrimAmp=1,NumBornAmps
                call EvalTree(BornAmps(iPrimAmp))
            enddo

            LO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
            do iPrimAmp=1,NumBornAmps
                LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
            enddo
            enddo
            LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop
        DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

        do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),DipoleResult)
        enddo
        EvalCS_DKJ_Real_ststbqqb_2 = EvalCS_DKJ_Real_ststbqqb_2 + DipoleResult


   elseif( nJetRad.eq.2 ) then

        MomExtTd(1:4,1:11) = MomExt(1:4,1:11)
        MomExtTdIn(1:4,1:4) = MomExt(1:4,12:15)
        call WTransform(MomExtTdIn(1:4,1:4),MomExtTd(1:4,12:14),pbDpg,ptDpg,ptDpb)
        
        omz=ptDpg/(ptDpb+ptDpg-pbDpg)  !  for some reason this is not (1-z) as defined in the paper...
        rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
        z=1d0-omz
        y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2

        TheDipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
        TheDipole = TheDipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )


        call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
        if( applyPSCut .or. TheDipole.eq.0d0) then
            cycle
        endif


        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
            call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExtTd(1:4,5:9))
            call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExtTd(1:4,10:14))
            call HelCrossing(Helicities(iHel,1:4))
            call SetPolarizations()
            do iPrimAmp=1,NumBornAmps
                call EvalTree(BornAmps(iPrimAmp))
            enddo

            LO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
            do iPrimAmp=1,NumBornAmps
                LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
            enddo
            enddo
            LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop
        DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

        do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),DipoleResult)
        enddo
        EvalCS_DKJ_Real_ststbqqb_2 = EvalCS_DKJ_Real_ststbqqb_2 + DipoleResult


   elseif( nJetRad.eq.3 ) then
          MomExtTd(1:4,1:11) = MomExt(1:4,1:11)
          MomExtTdIn(1:4,1:4) = MomExt(1:4,12:15)
          do ndip=1,2
              call wdec_trans(ndip,MomExtTdIn(1:4,1:4),MomExtTd(1:4,12:14),alpha_DKWff,TheDipole)
              TheDipole = - alpha_s4Pi*RunFactor * CF * TheDipole

              call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
              if( applyPSCut .or. TheDipole.eq.0d0 ) then
                  cycle
              endif

              LO_Res_Unpol = (0d0,0d0)
              do iHel=1,NumHelicities! helicity summation
                  call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExtTd(1:4,5:9))
                  call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExtTd(1:4,10:14))
                  call HelCrossing(Helicities(iHel,1:4))
                  call SetPolarizations()
                  do iPrimAmp=1,NumBornAmps
                      call EvalTree(BornAmps(iPrimAmp))
                  enddo

                  LO_Res_Pol = (0d0,0d0)
                  do jPrimAmp=1,NumBornAmps
                  do iPrimAmp=1,NumBornAmps
                      LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
                  enddo
                  enddo
                  LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
              enddo!helicity loop
              DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

              do NHisto=1,NumHistograms
                  call intoHisto(NHisto,NBin(NHisto),DipoleResult)
              enddo
              EvalCS_DKJ_Real_ststbqqb_2 = EvalCS_DKJ_Real_ststbqqb_2 + DipoleResult
          enddo   !dipole loop
  endif! nJetRad for dipoles

enddo! nJetRad loop





   EvalCS_DKJ_Real_ststbqqb = (EvalCS_DKJ_Real_ststbqqb_1+EvalCS_DKJ_Real_ststbqqb_2)/VgsWgt



   if( IsNan(EvalCS_DKJ_Real_ststbqqb) ) then
        print *, "NAN:",EvalCS_DKJ_Real_ststbqqb
        print *, yRnd(:)
        print *, PSWgt , VgsWgt , PDFFac, sHatJacobi
        print *, eta1,eta2,MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,1:11)
        print *, "SKIP EVENT!!!!!"
        EvalCS_DKJ_Real_ststbqqb = 0d0
        return
   endif


return
END FUNCTION






! should be named DK1L
FUNCTION EvalCS_DKJ_1L_ststbgg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
implicit none
real(8) ::  EvalCS_DKJ_1L_ststbgg,yRnd(1:VegasMxDim),VgsWgt
complex(8) :: NLO_Res_Pol,NLO_Res_Unpol,TreeResult(1:6),VirtResult(1:6),LO_Res_Unpol
integer :: iHel,jHel,iPrimAmp,jPrimAmp,nHel(1:2)
integer :: NBin(1:NumMaxHisto),NHisto
real(8) :: SpinAvg,ColorAvg,EHat,PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5,ISFac,RunFactor
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac
real(8) :: MomExt(1:4,1:15)
real(8) :: pdf(-6:6,1:2)
logical :: applyPSCut
include 'vegas_common.f'


   EvalCS_DKJ_1L_ststbgg = 0d0
   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if( EHat.le.2d0*m_STop ) then
      EvalCS_DKJ_1L_ststbgg = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)
   call EvalPhaseSpace_2to2Stops(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt)! AStop, Stop
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
   ISFac = MomCrossing(MomExt)

   IF(XTOPDECAYS.EQ.3) THEN
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  Chi top
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,4),yRnd(7:8),MomExt(1:4,10:11),PSWgt3)
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd( 9:12),MomExt(1:4,7:9),PSWgt4)! bot lep neu
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(13:16),MomExt(1:4,12:14),PSWgt5)
      PSWgt = PSWgt * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
   ENDIF
 

   call Kinematics_TTbarETmiss(.false.,MomExt,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_DKJ_1L_ststbgg = 0d0
      return
   endif

   call SetPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   RunFactor = RunAlphaS(NLOParam,MuRen)
   nHel(1:2) = getHelicity(yrnd(17))
   PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))


!----------------------------------------
! one loop correction to anti-stop decay |
!----------------------------------------
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
        call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExt(1:4,5:9))
        call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,10:14))
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call STopDecay(ExtParticle(1),DKX_STChi0_1L1,Helicities(iHel,5),MomExt(1:4,5:9))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol

    enddo! iHel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_ststbgg = NLO_Res_Unpol




!----------------------------------------
! one loop correction to anti-top decay |
!----------------------------------------
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
        call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExt(1:4,5:9))
        call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,10:14))
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call STopDecay(ExtParticle(1),DKX_STChi0_1L2,Helicities(iHel,5),MomExt(1:4,5:9))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol

    enddo! iHel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_ststbgg = EvalCS_DKJ_1L_ststbgg + NLO_Res_Unpol



!----------------------------------------
! one loop correction to W- decay       |
!----------------------------------------
  IF( TOPDECAYS.EQ.4 .OR. TOPDECAYS.EQ.2 ) THEN
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
        call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExt(1:4,5:9))
        call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,10:14))
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call STopDecay(ExtParticle(1),DKX_STChi0_1L3,Helicities(iHel,5),MomExt(1:4,5:9))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol

    enddo! iHel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_ststbgg = EvalCS_DKJ_1L_ststbgg + NLO_Res_Unpol
  ENDIF






!----------------------------------------
! one loop correction on stop decay     |
! ----------------------------------------
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
        call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExt(1:4,5:9))
        call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,10:14))
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call STopDecay(ExtParticle(2),DKX_STChi0_1L1,Helicities(iHel,6),MomExt(1:4,10:14))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol

    enddo! iHel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_ststbgg = EvalCS_DKJ_1L_ststbgg + NLO_Res_Unpol


! !!!!!!!!!!!!! COMPARISON WITH RADJA !!!!!!!!!!!!!!!
! 
! 
! LO_Res_UnPol = (0d0,0d0)
! NLO_Res_UnPol = (0d0,0d0)
! do ihel=-1,1,2
! do jhel=-1,1,2
! 
!  print *, "hel",ihel,jhel
!  call STopDecay(ExtParticle(2),DKX_STChi0_LO,iHel,MomExt(1:4,10:14),HelTop=jhel)
!  TreeResult(1)=ExtParticle(2)%Pol(1)    *dsqrt(2d0*Ga_STop*M_STop)
!  call STopDecay(ExtParticle(2),DKX_STChi0_1L1,iHel,MomExt(1:4,10:14),HelTop=jhel)
!  VirtResult(1)=ExtParticle(2)%Pol(1)    *dsqrt(2d0*Ga_STop*M_STop)/(alpha_sOver2Pi*RunAlphaS(NLOParam,MuRen))
!  print *, "tree",treeresult(1)
!  print *, "virt",virtresult(1)
!  
!  LO_Res_UnPol  = LO_Res_UnPol + TreeResult(1)*dconjg(TreeResult(1))
!  NLO_Res_UnPol = NLO_Res_UnPol + dreal(  TreeResult(1)*dconjg(VirtResult(1))   )
! 
! enddo
! enddo
! 
! print *, "--unpol--"
! write(*,"(A,3F10.4)") "mstop,mchi,mtop",m_stop*100d0,m_chi*100d0,m_top*100d0
! write(*,"(A,3F10.4)") "muren",MuRen*100d0
! ga_stop_chitop=( (IChiStt(+1)**2+IChiStt(-1)**2)*(m_stop**2-m_top**2-m_chi**2) - 4d0*(IChiStt(+1)*IChiStt(-1))*m_top*m_chi )
! 
! print *, "tree",LO_Res_UnPol
! !print *, "my LO width",Ga_Stop_ChiTop
! !print *, "LO ratio",LO_Res_UnPol/(Ga_Stop_ChiTop)
! 
! nlo_res_unpol=nlo_res_unpol*2d0
! print *, "my nlo( 0)",NLO_Res_UnPol/LO_Res_UnPol
! print *, "radja(-1)",6.35546544566007d0!5.19664786378765d0 !3.45067452120216d0
! print *, "ratio to radja(-1)",6.35546544566007d0/(NLO_Res_UnPol/LO_Res_UnPol)
! print *, "radja( 0)",11.4104415737049d0!7.94514801248857!2.47361046117698d0
! print *, "ratio to radja( 0)",11.4104415737049d0/(NLO_Res_UnPol/LO_Res_UnPol)
! pause
! 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!----------------------------------------
! one loop correction on top decay      |
!----------------------------------------
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
        call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExt(1:4,5:9))
        call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,10:14))
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call STopDecay(ExtParticle(2),DKX_STChi0_1L2,Helicities(iHel,6),MomExt(1:4,10:14))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol

    enddo! iHel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_ststbgg = EvalCS_DKJ_1L_ststbgg + NLO_Res_Unpol




!----------------------------------------
! one loop correction on W+ decay       |
!----------------------------------------
  IF( TOPDECAYS.EQ.3 .OR. TOPDECAYS.EQ.2 ) THEN
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
        call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExt(1:4,5:9))
        call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,10:14))
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call STopDecay(ExtParticle(2),DKX_STChi0_1L3,Helicities(iHel,6),MomExt(1:4,10:14))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol

    enddo! iHel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_ststbgg = EvalCS_DKJ_1L_ststbgg + NLO_Res_Unpol
  ENDIF




   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_DKJ_1L_ststbgg)
   enddo


   EvalCS_DKJ_1L_ststbgg = EvalCS_DKJ_1L_ststbgg/VgsWgt
   if( IsNan(EvalCS_DKJ_1L_ststbgg) ) then
        print *, "NAN:",EvalCS_DKJ_1L_ststbgg
        call printYrnd(yrnd(:))
        EvalCS_DKJ_1L_ststbgg = 0d0
        return
   endif

return
END FUNCTION











FUNCTION EvalCS_DKJ_1L_ststbqqb(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
implicit none
real(8) ::  EvalCS_DKJ_1L_ststbqqb
real(8) :: yRnd(1:VegasMxDim),VgsWgt
complex(8) :: NLO_Res_Pol,NLO_Res_Unpol,TreeResult(1:6),VirtResult(1:6),LO_Res_Unpol
integer :: iHel,jHel,iPrimAmp,jPrimAmp,nHel(1:2)
integer :: NBin(1:NumMaxHisto),NHisto
real(8) :: SpinAvg,ColorAvg,EHat,PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5,ISFac,RunFactor
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,PDFFac_a,PDFFac_b
real(8) :: MomExt(1:4,1:15)
real(8) :: pdf(-6:6,1:2)
logical :: applyPSCut
include 'vegas_common.f'




EvalCS_DKJ_1L_ststbqqb= 0d0


   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if( EHat.le.2d0*m_STop ) then
      EvalCS_DKJ_1L_ststbqqb = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)
   call EvalPhaseSpace_2to2Stops(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt)! AStop, Stop
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
   ISFac = MomCrossing(MomExt)

   IF(XTOPDECAYS.EQ.3) THEN
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  Chi top
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,4),yRnd(7:8),MomExt(1:4,10:11),PSWgt3)
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd( 9:12),MomExt(1:4,7:9),PSWgt4)! bot lep neu
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(13:16),MomExt(1:4,12:14),PSWgt5)
      PSWgt = PSWgt * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
   ENDIF
 
   call Kinematics_TTbarETmiss(.false.,MomExt,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_DKJ_1L_ststbqqb = 0d0
      return
   endif

   call SetPDFs(eta1,eta2,MuFac,pdf)
   PDFFac_a = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
            + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
            + pdf(Bot_,1)*pdf(ABot_,2)
   PDFFac_b = pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
            + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
            + pdf(Bot_,2)*pdf(ABot_,1)
   PDFFac = PDFFac_a + PDFFac_b
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   RunFactor = RunAlphaS(NLOParam,MuRen)
   nHel(1:2) = getHelicity(yrnd(17))
   PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))


!----------------------------------------
! one loop correction to anti-stop decay |
!----------------------------------------
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
        call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExt(1:4,5:9))
        call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,10:14))
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call STopDecay(ExtParticle(1),DKX_STChi0_1L1,Helicities(iHel,5),MomExt(1:4,5:9))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol

    enddo! iHel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_ststbqqb = NLO_Res_Unpol




!----------------------------------------
! one loop correction to anti-top decay |
! ----------------------------------------
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
        call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExt(1:4,5:9))
        call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,10:14))
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call STopDecay(ExtParticle(1),DKX_STChi0_1L2,Helicities(iHel,5),MomExt(1:4,5:9))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol

    enddo! iHel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_ststbqqb = EvalCS_DKJ_1L_ststbqqb + NLO_Res_Unpol



!----------------------------------------
! one loop correction to W- decay       |
!----------------------------------------
  IF( TOPDECAYS.EQ.2 .OR. TOPDECAYS.EQ.4 ) THEN
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
        call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExt(1:4,5:9))
        call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,10:14))
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call STopDecay(ExtParticle(1),DKX_STChi0_1L3,Helicities(iHel,5),MomExt(1:4,5:9))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol

    enddo! iHel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_ststbqqb = EvalCS_DKJ_1L_ststbqqb + NLO_Res_Unpol
  ENDIF






!----------------------------------------
! one loop correction on stop decay     |
!----------------------------------------
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
        call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExt(1:4,5:9))
        call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,10:14))
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call STopDecay(ExtParticle(2),DKX_STChi0_1L1,Helicities(iHel,6),MomExt(1:4,10:14))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol

    enddo! iHel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_ststbqqb = EvalCS_DKJ_1L_ststbqqb + NLO_Res_Unpol





!----------------------------------------
! one loop correction on top decay      |
!----------------------------------------
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
        call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExt(1:4,5:9))
        call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,10:14))
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call STopDecay(ExtParticle(2),DKX_STChi0_1L2,Helicities(iHel,6),MomExt(1:4,10:14))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol

    enddo! iHel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_ststbqqb = EvalCS_DKJ_1L_ststbqqb + NLO_Res_Unpol




!----------------------------------------
! one loop correction on W+ decay       |
!----------------------------------------
  IF( TOPDECAYS.EQ.2 .OR. TOPDECAYS.EQ.3 ) THEN
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
        call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExt(1:4,5:9))
        call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,10:14))
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call STopDecay(ExtParticle(2),DKX_STChi0_1L3,Helicities(iHel,6),MomExt(1:4,10:14))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol

    enddo! iHel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_ststbqqb = EvalCS_DKJ_1L_ststbqqb + NLO_Res_Unpol
  ENDIF




   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_DKJ_1L_ststbqqb)
   enddo


   EvalCS_DKJ_1L_ststbqqb = EvalCS_DKJ_1L_ststbqqb/VgsWgt
   if( IsNan(EvalCS_DKJ_1L_ststbqqb) ) then
        print *, "NAN:",EvalCS_DKJ_1L_ststbqqb
        call printYrnd(yrnd(:))
        EvalCS_DKJ_1L_ststbqqb = 0d0
        return
   endif



END FUNCTION











! should be named DKReal
FUNCTION EvalCS_DKJ_Real_HtHtbgg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes
use ModParameters
use ModExoticDecay
use ModHadrWDecay
use ModMisc
implicit none
real(8) ::  EvalCS_DKJ_Real_HtHtbgg,EvalCS_DKJ_Real_HtHtbgg_1,EvalCS_DKJ_Real_HtHtbgg_2
real(8) :: yRnd(1:VegasMxDim),VgsWgt
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol
integer :: iHel,jHel,kHel,iPrimAmp,jPrimAmp,BH1Hel,BH2Hel
real(8) :: EHat,RunFactor,PSWgt,PSWgt0,PSWgt2,PSWgt3,PSWgt4,PSWgt5,ISFac
real(8) :: MomExt(1:4,1:15),MomExtTd(1:4,1:15),MomExtTdIn(1:4,1:4),MomBoost(1:4)
logical :: applyPSCut
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,TheDipole,DipoleResult
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2), mydummy(1:10)
integer :: NBin(1:NumMaxHisto),NHisto,nHel(1:2)
integer :: nJetRad,nJetRad1,nJetRad2,nJetRad3,nJetRad4,GluHel,ndip
real(8) :: pbDpg,ptDpg,ptDpb,omz,rsq,z,y
real(8),parameter :: CF=4d0/3d0
include 'vegas_common.f'


   EvalCS_DKJ_Real_HtHtbgg   = 0d0
   EvalCS_DKJ_Real_HtHtbgg_1 = 0d0
   EvalCS_DKJ_Real_HtHtbgg_2 = 0d0

! yRnd( 1)=  0.3385585941088194d0
! yRnd( 2)=  0.2799513116385563d0
! yRnd( 3)=  0.012473622342792d0
! yRnd( 4)=  0.2879364093709448d0
! yRnd( 5)=  0.1334328211068331d0
! yRnd( 6)=  0.7829718273519412d0
! yRnd( 7)=  0.3479862101366653d0
! yRnd( 8)=  0.1332233664734401d0
! yRnd( 9)=  0.2332185946559626d0
! yRnd(10)=  0.7471774192280964d0
! yRnd(11)=  0.438860277849650d0
! yRnd(12)=  0.02669250668338154d0
! yRnd(13)=  0.1692973097643342d0
! yRnd(14)=  0.6054243425398250d0
! yRnd(15)=  0.1081832926292833d0
! yRnd(16)=  0.1284279972912875d0
! yRnd(17)=  0.4608715537287632d0
! yRnd(18)=  0.3695206294159972d0
! yRnd(19)=  0.4436338681465174d0
! print *, "fixing yrnds"


   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if( EHat.le.2d0*m_HTop ) then
      EvalCS_DKJ_Real_HtHtbgg = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)
   call EvalPhaseSpace_2to2HT(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt0)! HTbar, HT
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
   ISFac = MomCrossing(MomExt)

   call SetPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
   RunFactor = RunAlphaS(NLOParam,MuRen)


   nJetRad1=1
   nJetRad3=1
   if( TOPDECAYS.eq.1 ) then
      nJetRad2=2;  nJetRad4=2
   elseif( TOPDECAYS.eq.2 ) then
      nJetRad2=3;  nJetRad4=3
   elseif( TOPDECAYS.eq.3 ) then
      nJetRad2=2;  nJetRad4=3
   elseif( TOPDECAYS.eq.4 ) then
      nJetRad2=3;  nJetRad4=2
   endif


! ------------------- for checks --------------------------
!    print *, "fixing nJetRad"
!    nJetRad1=2
!    nJetRad2=2
!    nJetRad3=2
!    nJetRad4=2
! ---------------------------------------------


!----------------------------------
!  gluon  emission off HTbar   |
!----------------------------------
do nJetRad=nJetRad1,nJetRad2!   nJetRad=1: gluon radiation off stop line,nJetRad=2: gluon radiation off top line,nJetRad=3: gluon radiation off W line
   EvalCS_DKJ_Real_HtHtbgg = 0d0
   IF(XTOPDECAYS.EQ.1) THEN
   if( nJetRad.eq.1 ) then
      call EvalPhasespace_HTopDK(HT_BH_T_G,MomExt(1:4,3),yRnd(5:9),MomExt(1:4,5:7),PSWgt2)!  BH(5) top(6) glu(7->15)
      MomExt(1:4,15) = MomExt(1:4,7)
      m_Top = m_SMTop
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd(10:13),MomExt(1:4,7:9),PSWgt3)!  b(7) el(8) nu(9)
      m_Top = m_HTop
   elseif( nJetRad.eq.2 ) then
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  BH(5) top(6)
      m_Top = m_SMTop
      call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,6),yRnd(7:13),MomExt(1:4,7:10),PSWgt3)!  b(7) el(8) nu(9) glu(10->15)
      m_Top = m_HTop
      MomExt(1:4,15) = MomExt(1:4,10)
   elseif( nJetRad.eq.3 ) then
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  BH(5) top(6)
      m_Top = m_SMTop
      call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,6),yRnd(7:13),MomExt(1:4,7:10),PSWgt3)!  b(7) el(8) nu(9) glu(10->15)
      m_Top = m_HTop
      MomExt(1:4,15) = MomExt(1:4,10)
   endif
   call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,4),yRnd(14:15),MomExt(1:4,10:11),PSWgt4)!  BH(10) top(11)
   m_Top = m_SMTop
   call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(16:19),MomExt(1:4,12:14),PSWgt5)! b(12) el(13) nu(14)
      m_Top = m_HTop
   ENDIF
   call Kinematics_TTbarETmiss(.true.,MomExt,(/3,4,5,10,6,11,7,8,9,12,13,14,15/),applyPSCut,NBin)

   PSWgt = PSWgt0 * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   if( PreFac.eq.0d0 ) then
        SkipCounter = SkipCounter + 1
        cycle
   endif
   if( applyPSCut ) then
      EvalCS_DKJ_Real_HtHtbgg = 0d0
      PSCutCounter = PSCutCounter + 1
      goto 2000
   endif


   LO_Res_Unpol = (0d0,0d0)
   do iHel=1,NumHelicities
   do GluHel=1,-1,-2
   do BH1Hel=-1,1
   do BH2Hel=-1,1
      call HelCrossing(Helicities(iHel,1:6))
      call SetPolarizations()
      if( nJetRad.eq.1 ) then
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_RE1,BH1Hel,MomExt(1:4,5:9),MomExt(1:4,15),GluHel)
      elseif( nJetRad.eq.2 ) then
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_RE2,BH1Hel,MomExt(1:4,5:9),MomExt(1:4,15),GluHel)
      else
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_RE3,BH1Hel,MomExt(1:4,5:9),MomExt(1:4,15),GluHel)
      endif
      call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExt(1:4,10:14))

      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo

      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol

   enddo!helicity loop
   enddo!helicity loop
   enddo!helicity loop
   enddo!helicity loop


   EvalCS_DKJ_Real_HtHtbgg = LO_Res_UnPol * PreFac * ISFac * (alpha_s4Pi*RunFactor)**2
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_DKJ_Real_HtHtbgg)
   enddo
   EvalCounter = EvalCounter + 1

   EvalCS_DKJ_Real_HtHtbgg_1 = EvalCS_DKJ_Real_HtHtbgg_1 + EvalCS_DKJ_Real_HtHtbgg

mydummy(1) = EvalCS_DKJ_Real_HtHtbgg


2000 continue!! dipoles for gluon emission off anti-stop


   if( nJetRad.eq.1 ) then

        call WTransform3(MomExt(1:4,5:6),MomExt(1:4,15),MomExtTd(1:4,5:6),pbDpg,ptDpg,ptDpb)
        m_Top=m_SMTop
        call EvalPhasespace_TopDK(T_B_W,MomExtTd(1:4,6),yRnd(10:13),MomExtTd(1:4,7:9),PSWgt3)
        m_Top=m_HTop
        MomExtTd(1:4,1:4) = MomExt(1:4,1:4)
        MomExtTd(1:4,10:14) = MomExt(1:4,10:14)
        
        TheDipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg/ptDpg*( m_Htop**2+m_SMTop**2-(MomExt(1:4,5).dot.MomExt(1:4,5)) ) - (m_HTop/ptDpg)**2 - (m_SMTop/pbDpg)**2 )

        call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
        if( applyPSCut .or. TheDipole.eq.0d0 ) then
            cycle
        endif

        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
        do BH1Hel=-1,1
        do BH2Hel=-1,1
            call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExtTd(1:4,5:9))
            call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExtTd(1:4,10:14))
            call HelCrossing(Helicities(iHel,1:4))
            call SetPolarizations()
            do iPrimAmp=1,NumBornAmps
                call EvalTree(BornAmps(iPrimAmp))
            enddo

            LO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
            do iPrimAmp=1,NumBornAmps
                LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
            enddo
            enddo
            LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop
        enddo!helicity loop
        enddo!helicity loop
        DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

        do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),DipoleResult)
        enddo
        EvalCS_DKJ_Real_HtHtbgg_1 = EvalCS_DKJ_Real_HtHtbgg_1 + DipoleResult

! mydummy(2) = DipoleResult
! print *, "E(glu)/E(stop)",MomExt(1,15)/MomExt(1,5)
! print *, "real    result",mydummy(1)
! print *, "dipole  result",mydummy(2)
! print *, "ratios        ",mydummy(1)/mydummy(2), mydummy(1)/mydummy(2)+1d0
! pause
! return

   elseif( nJetRad.eq.2 ) then

        MomExtTd(1:4,1:6) = MomExt(1:4,1:6)
        MomExtTd(1:4,10:14) = MomExt(1:4,10:14)
        MomExtTdIn(1:4,1:3) = MomExt(1:4,7:9)
        MomExtTdIn(1:4,4)  = MomExt(1:4,15)
        call WTransform(MomExtTdIn(1:4,1:4),MomExtTd(1:4,7:9),pbDpg,ptDpg,ptDpb)
       
        omz=ptDpg/(ptDpb+ptDpg-pbDpg)  !  for some reason this is not (1-z) as defined in the paper...
        rsq = 1d0 - 2d0/m_SMtop**2*(ptDpb+ptDpg-pbDpg)
        z=1d0-omz
        y=pbDpg*2d0/m_SMtop**2/(1d0-dsqrt(rsq))**2

        TheDipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_SMTop/ptDpg)**2 )
        TheDipole = TheDipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )


        call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
        if( applyPSCut .or. TheDipole.eq.0d0) then
            cycle
        endif

        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
        do BH1Hel=-1,1
        do BH2Hel=-1,1
            call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExtTd(1:4,5:9))
            call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExtTd(1:4,10:14))
            call HelCrossing(Helicities(iHel,1:4))
            call SetPolarizations()
            do iPrimAmp=1,NumBornAmps
                call EvalTree(BornAmps(iPrimAmp))
            enddo

            LO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
            do iPrimAmp=1,NumBornAmps
                LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
            enddo
            enddo
            LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop
        enddo!helicity loop
        enddo!helicity loop
        DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

        do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),DipoleResult)
        enddo
        EvalCS_DKJ_Real_HtHtbgg_1 = EvalCS_DKJ_Real_HtHtbgg_1 + DipoleResult


! mydummy(2) = DipoleResult
! print *, "sij/shat",(MomExt(1:4,15).dot.MomExt(1:4,7))/EHat**2
! print *, "E(glu)/E(stop)",MomExt(1,15)/MomExt(1,7)
! print *, "real    result",mydummy(1)
! print *, "dipole  result",mydummy(2)
! print *, "ratios        ",mydummy(1)/mydummy(2), mydummy(1)/mydummy(2)+1d0
! pause
! return

   elseif( nJetRad.eq.3 ) then
          MomExtTd(1:4,1:6) = MomExt(1:4,1:6)
          MomExtTd(1:4,10:14) = MomExt(1:4,10:14)
          MomExtTdIn(1:4,1:3) = MomExt(1:4,7:9)
          MomExtTdIn(1:4,4) = MomExt(1:4,15)
          do ndip=1,2
              call wdec_trans(ndip,MomExtTdIn(1:4,1:4),MomExtTd(1:4,7:9),alpha_DKWff,TheDipole)
              TheDipole = - alpha_s4Pi*RunFactor * CF * TheDipole

              call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
              if( applyPSCut .or. TheDipole.eq.0d0 ) then
                  cycle
              endif

              LO_Res_Unpol = (0d0,0d0)
              do iHel=1,NumHelicities! helicity summation
              do BH1Hel=-1,1
              do BH2Hel=-1,1
                  call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExtTd(1:4,5:9))
                  call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExtTd(1:4,10:14))
                  call HelCrossing(Helicities(iHel,1:4))
                  call SetPolarizations()
                  do iPrimAmp=1,NumBornAmps
                      call EvalTree(BornAmps(iPrimAmp))
                  enddo

                  LO_Res_Pol = (0d0,0d0)
                  do jPrimAmp=1,NumBornAmps
                  do iPrimAmp=1,NumBornAmps
                      LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
                  enddo
                  enddo
                  LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
              enddo!helicity loop
              enddo!helicity loop
              enddo!helicity loop
              DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

              do NHisto=1,NumHistograms
                  call intoHisto(NHisto,NBin(NHisto),DipoleResult)
              enddo
              EvalCS_DKJ_Real_HtHtbgg_1 = EvalCS_DKJ_Real_HtHtbgg_1 + DipoleResult
mydummy(2) = mydummy(2) + DipoleResult
          enddo   !dipole loop

! print *, "sij/shat",(MomExt(1:4,15).dot.MomExt(1:4,7))/EHat**2
! print *, "E(glu)/E(stop)",MomExt(1,15)/MomExt(1,7)
! print *, "real    result",mydummy(1)
! print *, "dipole  result",mydummy(2)
! print *, "ratios        ",mydummy(1)/mydummy(2), mydummy(1)/mydummy(2)+1d0
! pause
! return

  endif! nJetRad for dipoles

enddo! nJetRad loop






!----------------------------------
!  gluon  emission off HT        |
!----------------------------------
do nJetRad=nJetRad3,nJetRad4!   nJetRad=1: gluon radiation off stop line,nJetRad=2: gluon radiation off top line,nJetRad=3: gluon radiation off W line

   EvalCS_DKJ_Real_HtHtbgg = 0d0
   IF(XTOPDECAYS.EQ.1) THEN
   if( nJetRad.eq.1 ) then
      call EvalPhasespace_HTopDK(HT_BH_T_G,MomExt(1:4,4),yRnd(11:15),MomExt(1:4,10:12),PSWgt2)!  BH(10) top(11) glu(12->15)
      MomExt(1:4,15) = MomExt(1:4,12)
      m_Top = m_SMTop
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(16:19),MomExt(1:4,12:14),PSWgt3)!  b(12) el(13) nu(14)
      m_Top = m_HTop
   elseif( nJetRad.eq.2 ) then
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,4),yRnd(11:12),MomExt(1:4,10:11),PSWgt2)!  BH(10) top(11)
      m_Top = m_SMTop
      call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,11),yRnd(13:19),MomExt(1:4,12:15),PSWgt3)!  b(12) el(13) nu(14) glu(15)
      m_Top = m_HTop
   elseif( nJetRad.eq.3 ) then
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,4),yRnd(11:12),MomExt(1:4,10:11),PSWgt2)!  BH(10) top(11)
      m_Top = m_SMTop
      call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,11),yRnd(13:19),MomExt(1:4,12:15),PSWgt3)!  b(12) el(13) nu(14) glu(15)
      m_Top = m_HTop
   endif   
   call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt4)! BH(5) top(6)
   m_Top = m_SMTop
   call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd(7:10),MomExt(1:4,7:9),PSWgt5)! b(7) el(8) nu(9)
   m_Top = m_HTop
   ENDIF
   call Kinematics_TTbarETmiss(.true.,MomExt,(/3,4,5,10,6,11,7,8,9,12,13,14,15/),applyPSCut,NBin)
   PSWgt = PSWgt0 * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   if( PreFac.eq.0d0 ) then
        SkipCounter = SkipCounter + 1
        cycle
   endif
   if( applyPSCut ) then
      EvalCS_DKJ_Real_HtHtbgg = 0d0
      PSCutCounter = PSCutCounter + 1
      goto 2001
   endif


   LO_Res_Unpol = (0d0,0d0)
   do iHel=1,NumHelicities
   do GluHel=1,-1,-2
   do BH1Hel=-1,1
   do BH2Hel=-1,1
      call HelCrossing(Helicities(iHel,1:6))
      call SetPolarizations()
      call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExt(1:4,5:9))
      if( nJetRad.eq.1 ) then
          call HTopBHDecay(ExtParticle(2),DKX_HTBH_RE1,BH2Hel,MomExt(1:4,10:14),MomExt(1:4,15),GluHel)
      elseif( nJetRad.eq.2 ) then
          call HTopBHDecay(ExtParticle(2),DKX_HTBH_RE2,BH2Hel,MomExt(1:4,10:14),MomExt(1:4,15),GluHel)
      else
          call HTopBHDecay(ExtParticle(2),DKX_HTBH_RE3,BH2Hel,MomExt(1:4,10:14),MomExt(1:4,15),GluHel)
      endif
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo

      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol

   enddo!helicity loop
   enddo!helicity loop
   enddo!helicity loop
   enddo!helicity loop

   EvalCS_DKJ_Real_HtHtbgg = LO_Res_UnPol * PreFac * ISFac * (alpha_s4Pi*RunFactor)**2
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_DKJ_Real_HtHtbgg)
   enddo
   EvalCounter = EvalCounter + 1

   EvalCS_DKJ_Real_HtHtbgg_2 = EvalCS_DKJ_Real_HtHtbgg_2 + EvalCS_DKJ_Real_HtHtbgg

mydummy(3) = EvalCS_DKJ_Real_HtHtbgg


2001 continue!! dipoles for gluon emission off stop





   if( nJetRad.eq.1 ) then

        call WTransform3(MomExt(1:4,10:11),MomExt(1:4,15),MomExtTd(1:4,10:11),pbDpg,ptDpg,ptDpb)
        m_Top=m_SMTop
        call EvalPhasespace_TopDK(T_B_W,MomExtTd(1:4,11),yRnd(16:19),MomExtTd(1:4,12:14),PSWgt3)
        m_Top=m_HTop
        MomExtTd(1:4,1:9) = MomExt(1:4,1:9)

        TheDipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg/ptDpg*( m_Htop**2+m_SMTop**2-(MomExt(1:4,10).dot.MomExt(1:4,10)) ) - (m_HTop/ptDpg)**2 - (m_SMTop/pbDpg)**2 )

        call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
        if( applyPSCut .or. TheDipole.eq.0d0 ) then
            cycle
        endif


        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
        do BH1Hel=-1,1
        do BH2Hel=-1,1
            call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExtTd(1:4,5:9))
            call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExtTd(1:4,10:14))
            call HelCrossing(Helicities(iHel,1:4))
            call SetPolarizations()
            do iPrimAmp=1,NumBornAmps
                call EvalTree(BornAmps(iPrimAmp))
            enddo

            LO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
            do iPrimAmp=1,NumBornAmps
                LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
            enddo
            enddo
            LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop
        enddo!helicity loop
        enddo!helicity loop
        DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

        do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),DipoleResult)
        enddo
        EvalCS_DKJ_Real_HtHtbgg_2 = EvalCS_DKJ_Real_HtHtbgg_2 + DipoleResult


! mydummy(4) = DipoleResult
! print *, "E(glu)/E(stop)",MomExt(1,15)/MomExt(1,10)
! print *, "real    result",mydummy(3)
! print *, "dipole  result",mydummy(4)
! print *, "ratios        ",mydummy(3)/mydummy(4), mydummy(3)/mydummy(4)+1d0
! pause
! return


   elseif( nJetRad.eq.2 ) then

        MomExtTd(1:4,1:11) = MomExt(1:4,1:11)
        MomExtTdIn(1:4,1:4) = MomExt(1:4,12:15)
        call WTransform(MomExtTdIn(1:4,1:4),MomExtTd(1:4,12:14),pbDpg,ptDpg,ptDpb)
        
        omz=ptDpg/(ptDpb+ptDpg-pbDpg)  !  for some reason this is not (1-z) as defined in the paper...
        rsq = 1d0 - 2d0/m_SMtop**2*(ptDpb+ptDpg-pbDpg)
        z=1d0-omz
        y=pbDpg*2d0/m_SMtop**2/(1d0-dsqrt(rsq))**2

        TheDipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_SMTop/ptDpg)**2 )
        TheDipole = TheDipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )


        call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
        if( applyPSCut .or. TheDipole.eq.0d0) then
            cycle
        endif


        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
        do BH1Hel=-1,1
        do BH2Hel=-1,1
            call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExtTd(1:4,5:9))
            call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExtTd(1:4,10:14))
            call HelCrossing(Helicities(iHel,1:4))
            call SetPolarizations()
            do iPrimAmp=1,NumBornAmps
                call EvalTree(BornAmps(iPrimAmp))
            enddo

            LO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
            do iPrimAmp=1,NumBornAmps
                LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
            enddo
            enddo
            LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop
        enddo!helicity loop
        enddo!helicity loop
        DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

        do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),DipoleResult)
        enddo
        EvalCS_DKJ_Real_HtHtbgg_2 = EvalCS_DKJ_Real_HtHtbgg_2 + DipoleResult


! mydummy(4) = DipoleResult
! print *, "sij/shat",(MomExt(1:4,15).dot.MomExt(1:4,12))/EHat**2
! print *, "E(glu)/E(top)",MomExt(1,15)/MomExt(1,12)
! print *, "real    result",mydummy(3)
! print *, "dipole  result",mydummy(4)
! print *, "ratios        ",mydummy(3)/mydummy(4), mydummy(3)/mydummy(4)+1d0
! pause
! return

   elseif( nJetRad.eq.3 ) then
          MomExtTd(1:4,1:11) = MomExt(1:4,1:11)
          MomExtTdIn(1:4,1:4) = MomExt(1:4,12:15)
          do ndip=1,2
              call wdec_trans(ndip,MomExtTdIn(1:4,1:4),MomExtTd(1:4,12:14),alpha_DKWff,TheDipole)
              TheDipole = - alpha_s4Pi*RunFactor * CF * TheDipole

              call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
              if( applyPSCut .or. TheDipole.eq.0d0 ) then
                  cycle
              endif

              LO_Res_Unpol = (0d0,0d0)
              do iHel=1,NumHelicities! helicity summation
              do BH1Hel=-1,1
              do BH2Hel=-1,1
                  call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExtTd(1:4,5:9))
                  call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExtTd(1:4,10:14))
                  call HelCrossing(Helicities(iHel,1:4))
                  call SetPolarizations()
                  do iPrimAmp=1,NumBornAmps
                      call EvalTree(BornAmps(iPrimAmp))
                  enddo

                  LO_Res_Pol = (0d0,0d0)
                  do jPrimAmp=1,NumBornAmps
                  do iPrimAmp=1,NumBornAmps
                      LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
                  enddo
                  enddo
                  LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
              enddo!helicity loop
              enddo!helicity loop
              enddo!helicity loop
              DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

              do NHisto=1,NumHistograms
                  call intoHisto(NHisto,NBin(NHisto),DipoleResult)
              enddo
              EvalCS_DKJ_Real_HtHtbgg_2 = EvalCS_DKJ_Real_HtHtbgg_2 + DipoleResult
          enddo   !dipole loop
  endif! nJetRad for dipoles

enddo! nJetRad loop





   EvalCS_DKJ_Real_HtHtbgg = (EvalCS_DKJ_Real_HtHtbgg_1+EvalCS_DKJ_Real_HtHtbgg_2)/VgsWgt



   if( IsNan(EvalCS_DKJ_Real_HtHtbgg) ) then
        print *, "NAN:",EvalCS_DKJ_Real_HtHtbgg
        print *, yRnd(:)
        print *, PSWgt , VgsWgt , PDFFac, sHatJacobi
        print *, eta1,eta2,MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,1:11)
        print *, "SKIP EVENT!!!!!"
        EvalCS_DKJ_Real_HtHtbgg = 0d0
        return
   endif


return
END FUNCTION






! should be named DKReal
FUNCTION EvalCS_DKJ_Real_HtHtbqqb(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes
use ModParameters
use ModExoticDecay
use ModHadrWDecay
use ModMisc
implicit none
real(8) ::  EvalCS_DKJ_Real_HtHtbqqb,EvalCS_DKJ_Real_HtHtbqqb_1,EvalCS_DKJ_Real_HtHtbqqb_2
real(8) :: yRnd(1:VegasMxDim),VgsWgt
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol
integer :: iHel,jHel,kHel,iPrimAmp,jPrimAmp,BH1Hel,BH2Hel
real(8) :: EHat,RunFactor,PSWgt,PSWgt0,PSWgt2,PSWgt3,PSWgt4,PSWgt5,ISFac
real(8) :: MomExt(1:4,1:15),MomExtTd(1:4,1:15),MomExtTdIn(1:4,1:4),MomBoost(1:4)
logical :: applyPSCut
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,TheDipole,DipoleResult,PDFFac_a,PDFFac_b
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2), mydummy(1:10)
integer :: NBin(1:NumMaxHisto),NHisto,nHel(1:2)
integer :: nJetRad,nJetRad1,nJetRad2,nJetRad3,nJetRad4,GluHel,ndip
real(8) :: pbDpg,ptDpg,ptDpb,omz,rsq,z,y
real(8),parameter :: CF=4d0/3d0
include 'vegas_common.f'


   EvalCS_DKJ_Real_HtHtbqqb   = 0d0
   EvalCS_DKJ_Real_HtHtbqqb_1 = 0d0
   EvalCS_DKJ_Real_HtHtbqqb_2 = 0d0

! yRnd( 1)=  0.3385585941088194d0
! yRnd( 2)=  0.2799513116385563d0
! yRnd( 3)=  0.012473622342792d0
! yRnd( 4)=  0.2879364093709448d0
! yRnd( 5)=  0.1334328211068331d0
! yRnd( 6)=  0.7829718273519412d0
! yRnd( 7)=  0.3479862101366653d0
! yRnd( 8)=  0.1332233664734401d0
! yRnd( 9)=  0.2332185946559626d0
! yRnd(10)=  0.7471774192280964d0
! yRnd(11)=  0.438860277849650d0
! yRnd(12)=  0.02669250668338154d0
! yRnd(13)=  0.1692973097643342d0
! yRnd(14)=  0.6054243425398250d0
! yRnd(15)=  0.1081832926292833d0
! yRnd(16)=  0.1284279972912875d0
! yRnd(17)=  0.4608715537287632d0
! yRnd(18)=  0.3695206294159972d0
! yRnd(19)=  0.4436338681465174d0
! print *, "fixing yrnds"


   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if( EHat.le.2d0*m_HTop ) then
      EvalCS_DKJ_Real_HtHtbqqb = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)
   call EvalPhaseSpace_2to2HT(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt0)! HTbar, HT
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
   ISFac = MomCrossing(MomExt)

   call SetPDFs(eta1,eta2,MuFac,pdf)
   PDFFac_a = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
            + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
            + pdf(Bot_,1)*pdf(ABot_,2)
   PDFFac_b = pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
            + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
            + pdf(Bot_,2)*pdf(ABot_,1)
   PDFFac = PDFFac_a + PDFFac_b

   RunFactor = RunAlphaS(NLOParam,MuRen)


   nJetRad1=1
   nJetRad3=1
   if( TOPDECAYS.eq.1 ) then
      nJetRad2=2;  nJetRad4=2
   elseif( TOPDECAYS.eq.2 ) then
      nJetRad2=3;  nJetRad4=3
   elseif( TOPDECAYS.eq.3 ) then
      nJetRad2=2;  nJetRad4=3
   elseif( TOPDECAYS.eq.4 ) then
      nJetRad2=3;  nJetRad4=2
   endif


! ------------------- for checks --------------------------
!    print *, "fixing nJetRad"
!    nJetRad1=1
!    nJetRad2=-1
!    nJetRad3=2
!    nJetRad4=2
! ---------------------------------------------



!----------------------------------
!  gluon  emission off HTbar   |
!----------------------------------
do nJetRad=nJetRad1,nJetRad2!   nJetRad=1: gluon radiation off stop line,nJetRad=2: gluon radiation off top line,nJetRad=3: gluon radiation off W line
   EvalCS_DKJ_Real_HtHtbqqb = 0d0
   IF(XTOPDECAYS.EQ.1) THEN
   if( nJetRad.eq.1 ) then
      call EvalPhasespace_HTopDK(HT_BH_T_G,MomExt(1:4,3),yRnd(5:9),MomExt(1:4,5:7),PSWgt2)!  BH(5) top(6) glu(7->15)
      MomExt(1:4,15) = MomExt(1:4,7)
      m_Top = m_SMTop
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd(10:13),MomExt(1:4,7:9),PSWgt3)!  b(7) el(8) nu(9)
      m_Top = m_HTop
   elseif( nJetRad.eq.2 ) then
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  BH(5) top(6)
      m_Top = m_SMTop
      call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,6),yRnd(7:13),MomExt(1:4,7:10),PSWgt3)!  b(7) el(8) nu(9) glu(10->15)
      m_Top = m_HTop
      MomExt(1:4,15) = MomExt(1:4,10)
   elseif( nJetRad.eq.3 ) then
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  BH(5) top(6)
      m_Top = m_SMTop
      call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,6),yRnd(7:13),MomExt(1:4,7:10),PSWgt3)!  b(7) el(8) nu(9) glu(10->15)
      m_Top = m_HTop
      MomExt(1:4,15) = MomExt(1:4,10)
   endif
   call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,4),yRnd(14:15),MomExt(1:4,10:11),PSWgt4)!  BH(10) top(11)
   m_Top = m_SMTop
   call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(16:19),MomExt(1:4,12:14),PSWgt5)! b(12) el(13) nu(14)
   m_Top = m_HTop
   ENDIF
   call Kinematics_TTbarETmiss(.true.,MomExt,(/3,4,5,10,6,11,7,8,9,12,13,14,15/),applyPSCut,NBin)

   PSWgt = PSWgt0 * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   if( PreFac.eq.0d0 ) then
        SkipCounter = SkipCounter + 1
        cycle
   endif
   if( applyPSCut ) then
      EvalCS_DKJ_Real_HtHtbqqb = 0d0
      PSCutCounter = PSCutCounter + 1
      goto 2000
   endif


   LO_Res_Unpol = (0d0,0d0)
   do iHel=1,NumHelicities
   do GluHel=1,-1,-2
   do BH1Hel=-1,1
   do BH2Hel=-1,1
      call HelCrossing(Helicities(iHel,1:6))
      call SetPolarizations()
      if( nJetRad.eq.1 ) then
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_RE1,BH1Hel,MomExt(1:4,5:9),MomExt(1:4,15),GluHel)
      elseif( nJetRad.eq.2 ) then
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_RE2,BH1Hel,MomExt(1:4,5:9),MomExt(1:4,15),GluHel)
      else
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_RE3,BH1Hel,MomExt(1:4,5:9),MomExt(1:4,15),GluHel)
      endif
      call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExt(1:4,10:14))

      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo

      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol

   enddo!helicity loop
   enddo!helicity loop
   enddo!helicity loop
   enddo!helicity loop


   EvalCS_DKJ_Real_HtHtbqqb = LO_Res_UnPol * PreFac * ISFac * (alpha_s4Pi*RunFactor)**2
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_DKJ_Real_HtHtbqqb)
   enddo
   EvalCounter = EvalCounter + 1

   EvalCS_DKJ_Real_HtHtbqqb_1 = EvalCS_DKJ_Real_HtHtbqqb_1 + EvalCS_DKJ_Real_HtHtbqqb

mydummy(1) = EvalCS_DKJ_Real_HtHtbqqb


2000 continue!! dipoles for gluon emission off anti-stop



   if( nJetRad.eq.1 ) then

        call WTransform3(MomExt(1:4,5:6),MomExt(1:4,15),MomExtTd(1:4,5:6),pbDpg,ptDpg,ptDpb)
        m_Top=m_SMTop
        call EvalPhasespace_TopDK(T_B_W,MomExtTd(1:4,6),yRnd(10:13),MomExtTd(1:4,7:9),PSWgt3)
        m_Top=m_HTop
        MomExtTd(1:4,1:4) = MomExt(1:4,1:4)
        MomExtTd(1:4,10:14) = MomExt(1:4,10:14)
        
        TheDipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg/ptDpg*( m_Htop**2+m_SMTop**2-(MomExt(1:4,5).dot.MomExt(1:4,5)) ) - (m_HTop/ptDpg)**2 - (m_SMTop/pbDpg)**2 )

        call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
        if( applyPSCut .or. TheDipole.eq.0d0 ) then
            cycle
        endif

        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
        do BH1Hel=-1,1
        do BH2Hel=-1,1
            call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExtTd(1:4,5:9))
            call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExtTd(1:4,10:14))
            call HelCrossing(Helicities(iHel,1:4))
            call SetPolarizations()
            do iPrimAmp=1,NumBornAmps
                call EvalTree(BornAmps(iPrimAmp))
            enddo

            LO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
            do iPrimAmp=1,NumBornAmps
                LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
            enddo
            enddo
            LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop
        enddo!helicity loop
        enddo!helicity loop
        DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

        do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),DipoleResult)
        enddo
        EvalCS_DKJ_Real_HtHtbqqb_1 = EvalCS_DKJ_Real_HtHtbqqb_1 + DipoleResult

! mydummy(2) = DipoleResult
! print *, "E(glu)/E(stop)",MomExt(1,15)/MomExt(1,5)
! print *, "real    result",mydummy(1)
! print *, "dipole  result",mydummy(2)
! print *, "ratios        ",mydummy(1)/mydummy(2), mydummy(1)/mydummy(2)+1d0
! pause
! return

   elseif( nJetRad.eq.2 ) then

        MomExtTd(1:4,1:6) = MomExt(1:4,1:6)
        MomExtTd(1:4,10:14) = MomExt(1:4,10:14)
        MomExtTdIn(1:4,1:3) = MomExt(1:4,7:9)
        MomExtTdIn(1:4,4)  = MomExt(1:4,15)
        call WTransform(MomExtTdIn(1:4,1:4),MomExtTd(1:4,7:9),pbDpg,ptDpg,ptDpb)
        
        omz=ptDpg/(ptDpb+ptDpg-pbDpg)  !  for some reason this is not (1-z) as defined in the paper...
        rsq = 1d0 - 2d0/m_SMtop**2*(ptDpb+ptDpg-pbDpg)
        z=1d0-omz
        y=pbDpg*2d0/m_SMtop**2/(1d0-dsqrt(rsq))**2

        TheDipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_SMTop/ptDpg)**2 )
        TheDipole = TheDipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )


        call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
        if( applyPSCut .or. TheDipole.eq.0d0) then
            cycle
        endif


        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
        do BH1Hel=-1,1
        do BH2Hel=-1,1
            call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExtTd(1:4,5:9))
            call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExtTd(1:4,10:14))
            call HelCrossing(Helicities(iHel,1:4))
            call SetPolarizations()
            do iPrimAmp=1,NumBornAmps
                call EvalTree(BornAmps(iPrimAmp))
            enddo

            LO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
            do iPrimAmp=1,NumBornAmps
                LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
            enddo
            enddo
            LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop
        enddo!helicity loop
        enddo!helicity loop
        DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

        do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),DipoleResult)
        enddo
        EvalCS_DKJ_Real_HtHtbqqb_1 = EvalCS_DKJ_Real_HtHtbqqb_1 + DipoleResult


! mydummy(2) = DipoleResult
! print *, "sij/shat",(MomExt(1:4,15).dot.MomExt(1:4,7))/EHat**2
! print *, "E(glu)/E(stop)",MomExt(1,15)/MomExt(1,7)
! print *, "real    result",mydummy(1)
! print *, "dipole  result",mydummy(2)
! print *, "ratios        ",mydummy(1)/mydummy(2), mydummy(1)/mydummy(2)+1d0
! pause
! return

   elseif( nJetRad.eq.3 ) then
          MomExtTd(1:4,1:6) = MomExt(1:4,1:6)
          MomExtTd(1:4,10:14) = MomExt(1:4,10:14)
          MomExtTdIn(1:4,1:3) = MomExt(1:4,7:9)
          MomExtTdIn(1:4,4) = MomExt(1:4,15)
          do ndip=1,2
              call wdec_trans(ndip,MomExtTdIn(1:4,1:4),MomExtTd(1:4,7:9),alpha_DKWff,TheDipole)
              TheDipole = - alpha_s4Pi*RunFactor * CF * TheDipole

              call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
              if( applyPSCut .or. TheDipole.eq.0d0 ) then
                  cycle
              endif

              LO_Res_Unpol = (0d0,0d0)
              do iHel=1,NumHelicities! helicity summation
              do BH1Hel=-1,1
              do BH2Hel=-1,1
                  call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExtTd(1:4,5:9))
                  call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExtTd(1:4,10:14))
                  call HelCrossing(Helicities(iHel,1:4))
                  call SetPolarizations()
                  do iPrimAmp=1,NumBornAmps
                      call EvalTree(BornAmps(iPrimAmp))
                  enddo

                  LO_Res_Pol = (0d0,0d0)
                  do jPrimAmp=1,NumBornAmps
                  do iPrimAmp=1,NumBornAmps
                      LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
                  enddo
                  enddo
                  LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
              enddo!helicity loop
              enddo!helicity loop
              enddo!helicity loop
              DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

              do NHisto=1,NumHistograms
                  call intoHisto(NHisto,NBin(NHisto),DipoleResult)
              enddo
              EvalCS_DKJ_Real_HtHtbqqb_1 = EvalCS_DKJ_Real_HtHtbqqb_1 + DipoleResult
mydummy(2) = mydummy(2) + DipoleResult
          enddo   !dipole loop

! print *, "sij/shat",(MomExt(1:4,15).dot.MomExt(1:4,7))/EHat**2
! print *, "E(glu)/E(stop)",MomExt(1,15)/MomExt(1,7)
! print *, "real    result",mydummy(1)
! print *, "dipole  result",mydummy(2)
! print *, "ratios        ",mydummy(1)/mydummy(2), mydummy(1)/mydummy(2)+1d0
! pause
! return

  endif! nJetRad for dipoles

enddo! nJetRad loop






!----------------------------------
!  gluon  emission off HT        |
!----------------------------------
do nJetRad=nJetRad3,nJetRad4!   nJetRad=1: gluon radiation off stop line,nJetRad=2: gluon radiation off top line,nJetRad=3: gluon radiation off W line

   EvalCS_DKJ_Real_HtHtbqqb = 0d0
   IF(XTOPDECAYS.EQ.1) THEN
   if( nJetRad.eq.1 ) then
      call EvalPhasespace_HTopDK(HT_BH_T_G,MomExt(1:4,4),yRnd(11:15),MomExt(1:4,10:12),PSWgt2)!  BH(10) top(11) glu(12->15)
      MomExt(1:4,15) = MomExt(1:4,12)
      m_Top = m_SMTop
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(16:19),MomExt(1:4,12:14),PSWgt3)!  b(12) el(13) nu(14)
      m_Top = m_HTop
   elseif( nJetRad.eq.2 ) then
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,4),yRnd(11:12),MomExt(1:4,10:11),PSWgt2)!  BH(10) top(11)
      m_Top = m_SMTop
      call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,11),yRnd(13:19),MomExt(1:4,12:15),PSWgt3)!  b(12) el(13) nu(14) glu(15)
      m_Top = m_HTop
   elseif( nJetRad.eq.3 ) then
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,4),yRnd(11:12),MomExt(1:4,10:11),PSWgt2)!  BH(10) top(11)
      m_Top = m_SMTop
      call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,11),yRnd(13:19),MomExt(1:4,12:15),PSWgt3)!  b(12) el(13) nu(14) glu(15)
      m_Top = m_HTop
   endif   
   call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt4)! BH(5) top(6)
   m_Top = m_SMTop
   call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd(7:10),MomExt(1:4,7:9),PSWgt5)! b(7) el(8) nu(9)
   m_Top = m_HTop
   ENDIF
   call Kinematics_TTbarETmiss(.true.,MomExt,(/3,4,5,10,6,11,7,8,9,12,13,14,15/),applyPSCut,NBin)
   PSWgt = PSWgt0 * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   if( PreFac.eq.0d0 ) then
        SkipCounter = SkipCounter + 1
        cycle
   endif
   if( applyPSCut ) then
      EvalCS_DKJ_Real_HtHtbqqb = 0d0
      PSCutCounter = PSCutCounter + 1
      goto 2001
   endif


   LO_Res_Unpol = (0d0,0d0)
   do iHel=1,NumHelicities
   do GluHel=1,-1,-2
   do BH1Hel=-1,1
   do BH2Hel=-1,1
      call HelCrossing(Helicities(iHel,1:6))
      call SetPolarizations()
      call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExt(1:4,5:9))
      if( nJetRad.eq.1 ) then
          call HTopBHDecay(ExtParticle(2),DKX_HTBH_RE1,BH2Hel,MomExt(1:4,10:14),MomExt(1:4,15),GluHel)
      elseif( nJetRad.eq.2 ) then
          call HTopBHDecay(ExtParticle(2),DKX_HTBH_RE2,BH2Hel,MomExt(1:4,10:14),MomExt(1:4,15),GluHel)
      else
          call HTopBHDecay(ExtParticle(2),DKX_HTBH_RE3,BH2Hel,MomExt(1:4,10:14),MomExt(1:4,15),GluHel)
      endif


      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo

      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol

   enddo!helicity loop
   enddo!helicity loop
   enddo!helicity loop
   enddo!helicity loop

   EvalCS_DKJ_Real_HtHtbqqb = LO_Res_UnPol * PreFac * ISFac * (alpha_s4Pi*RunFactor)**2
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_DKJ_Real_HtHtbqqb)
   enddo
   EvalCounter = EvalCounter + 1

   EvalCS_DKJ_Real_HtHtbqqb_2 = EvalCS_DKJ_Real_HtHtbqqb_2 + EvalCS_DKJ_Real_HtHtbqqb

mydummy(3) = EvalCS_DKJ_Real_HtHtbqqb


2001 continue!! dipoles for gluon emission off stop




   if( nJetRad.eq.1 ) then

        call WTransform3(MomExt(1:4,10:11),MomExt(1:4,15),MomExtTd(1:4,10:11),pbDpg,ptDpg,ptDpb)
        m_Top=m_SMTop
        call EvalPhasespace_TopDK(T_B_W,MomExtTd(1:4,11),yRnd(16:19),MomExtTd(1:4,12:14),PSWgt3)
        m_Top=m_HTop
        MomExtTd(1:4,1:9) = MomExt(1:4,1:9)

        TheDipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg/ptDpg*( m_Htop**2+m_SMTop**2-(MomExt(1:4,10).dot.MomExt(1:4,10)) ) - (m_HTop/ptDpg)**2 - (m_SMTop/pbDpg)**2 )

        call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
        if( applyPSCut .or. TheDipole.eq.0d0 ) then
            cycle
        endif


        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
        do BH1Hel=-1,1
        do BH2Hel=-1,1
            call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExtTd(1:4,5:9))
            call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExtTd(1:4,10:14))
            call HelCrossing(Helicities(iHel,1:4))
            call SetPolarizations()
            do iPrimAmp=1,NumBornAmps
                call EvalTree(BornAmps(iPrimAmp))
            enddo

            LO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
            do iPrimAmp=1,NumBornAmps
                LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
            enddo
            enddo
            LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop
        enddo!helicity loop
        enddo!helicity loop
        DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

        do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),DipoleResult)
        enddo
        EvalCS_DKJ_Real_HtHtbqqb_2 = EvalCS_DKJ_Real_HtHtbqqb_2 + DipoleResult


! mydummy(4) = DipoleResult
! print *, "E(glu)/E(stop)",MomExt(1,15)/MomExt(1,10)
! print *, "real    result",mydummy(3)
! print *, "dipole  result",mydummy(4)
! print *, "ratios        ",mydummy(3)/mydummy(4), mydummy(3)/mydummy(4)+1d0
! pause
! return


   elseif( nJetRad.eq.2 ) then

        MomExtTd(1:4,1:11) = MomExt(1:4,1:11)
        MomExtTdIn(1:4,1:4) = MomExt(1:4,12:15)
        call WTransform(MomExtTdIn(1:4,1:4),MomExtTd(1:4,12:14),pbDpg,ptDpg,ptDpb)
        
        omz=ptDpg/(ptDpb+ptDpg-pbDpg)  !  for some reason this is not (1-z) as defined in the paper...
        rsq = 1d0 - 2d0/m_SMtop**2*(ptDpb+ptDpg-pbDpg)
        z=1d0-omz
        y=pbDpg*2d0/m_SMtop**2/(1d0-dsqrt(rsq))**2

        TheDipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_SMTop/ptDpg)**2 )
        TheDipole = TheDipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )


        call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
        if( applyPSCut .or. TheDipole.eq.0d0) then
            cycle
        endif


        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
        do BH1Hel=-1,1
        do BH2Hel=-1,1
            call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExtTd(1:4,5:9))
            call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExtTd(1:4,10:14))
            call HelCrossing(Helicities(iHel,1:4))
            call SetPolarizations()
            do iPrimAmp=1,NumBornAmps
                call EvalTree(BornAmps(iPrimAmp))
            enddo

            LO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
            do iPrimAmp=1,NumBornAmps
                LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
            enddo
            enddo
            LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop
        enddo!helicity loop
        enddo!helicity loop
        DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

        do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),DipoleResult)
        enddo
        EvalCS_DKJ_Real_HtHtbqqb_2 = EvalCS_DKJ_Real_HtHtbqqb_2 + DipoleResult


! mydummy(4) = DipoleResult
! print *, "sij/shat",(MomExt(1:4,15).dot.MomExt(1:4,12))/EHat**2
! print *, "E(glu)/E(top)",MomExt(1,15)/MomExt(1,12)
! print *, "real    result",mydummy(3)
! print *, "dipole  result",mydummy(4)
! print *, "ratios        ",mydummy(3)/mydummy(4), mydummy(3)/mydummy(4)+1d0
! pause
! return

   elseif( nJetRad.eq.3 ) then
          MomExtTd(1:4,1:11) = MomExt(1:4,1:11)
          MomExtTdIn(1:4,1:4) = MomExt(1:4,12:15)
          do ndip=1,2
              call wdec_trans(ndip,MomExtTdIn(1:4,1:4),MomExtTd(1:4,12:14),alpha_DKWff,TheDipole)
              TheDipole = - alpha_s4Pi*RunFactor * CF * TheDipole

              call Kinematics_TTbarETmiss(.false.,MomExtTd,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
              if( applyPSCut .or. TheDipole.eq.0d0 ) then
                  cycle
              endif

              LO_Res_Unpol = (0d0,0d0)
              do iHel=1,NumHelicities! helicity summation
              do BH1Hel=-1,1
              do BH2Hel=-1,1
                  call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExtTd(1:4,5:9))
                  call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExtTd(1:4,10:14))
                  call HelCrossing(Helicities(iHel,1:4))
                  call SetPolarizations()
                  do iPrimAmp=1,NumBornAmps
                      call EvalTree(BornAmps(iPrimAmp))
                  enddo

                  LO_Res_Pol = (0d0,0d0)
                  do jPrimAmp=1,NumBornAmps
                  do iPrimAmp=1,NumBornAmps
                      LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
                  enddo
                  enddo
                  LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
              enddo!helicity loop
              enddo!helicity loop
              enddo!helicity loop
              DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**2 * ISFac

              do NHisto=1,NumHistograms
                  call intoHisto(NHisto,NBin(NHisto),DipoleResult)
              enddo
              EvalCS_DKJ_Real_HtHtbqqb_2 = EvalCS_DKJ_Real_HtHtbqqb_2 + DipoleResult
          enddo   !dipole loop
  endif! nJetRad for dipoles

enddo! nJetRad loop





   EvalCS_DKJ_Real_HtHtbqqb = (EvalCS_DKJ_Real_HtHtbqqb_1+EvalCS_DKJ_Real_HtHtbqqb_2)/VgsWgt



   if( IsNan(EvalCS_DKJ_Real_HtHtbqqb) ) then
        print *, "NAN:",EvalCS_DKJ_Real_HtHtbqqb
        print *, yRnd(:)
        print *, PSWgt , VgsWgt , PDFFac, sHatJacobi
        print *, eta1,eta2,MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,1:11)
        print *, "SKIP EVENT!!!!!"
        EvalCS_DKJ_Real_HtHtbqqb = 0d0
        return
   endif


return
END FUNCTION







! should be named DK1L
FUNCTION EvalCS_DKJ_1L_HtHtbgg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
implicit none
real(8) ::  EvalCS_DKJ_1L_HtHtbgg,yRnd(1:VegasMxDim),VgsWgt
complex(8) :: NLO_Res_Pol,NLO_Res_Unpol,TreeResult(1:6),VirtResult(1:6),LO_Res_Unpol
integer :: iHel,jHel,iPrimAmp,jPrimAmp,nHel(1:2)
integer :: NBin(1:NumMaxHisto),NHisto,BHMaxHel,BH1Hel,BH2Hel
real(8) :: SpinAvg,ColorAvg,EHat,PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5,ISFac,RunFactor
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac
real(8) :: MomExt(1:4,1:15)
real(8) :: pdf(-6:6,1:2)
logical :: applyPSCut
include 'vegas_common.f'


   EvalCS_DKJ_1L_HtHtbgg = 0d0
   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if( EHat.le.2d0*m_HTop ) then
      EvalCS_DKJ_1L_HtHtbgg = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)
   call EvalPhaseSpace_2to2HT(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt)! HTbar, HT
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
   ISFac = MomCrossing(MomExt)

   BHMaxHel = -1
   IF(XTOPDECAYS.EQ.1) THEN
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  BH top
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,4),yRnd(7:8),MomExt(1:4,10:11),PSWgt3)

      m_Top = m_SMTop
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd( 9:12),MomExt(1:4,7:9),PSWgt4)! bot lep neu
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(13:16),MomExt(1:4,12:14),PSWgt5)
      m_Top = m_HTop

      PSWgt = PSWgt * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
      BHMaxHel = +1     
   ELSEIF(XTOPDECAYS.EQ.2) THEN
      call Error("XTOPDECAYS.EQ.2 is not supported")
   ENDIF



   call Kinematics_TTbarETmiss(.false.,MomExt,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_DKJ_1L_HtHtbgg = 0d0
      return
   endif

   call SetPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   RunFactor = RunAlphaS(NLOParam,MuRen)
   nHel(1:2) = getHelicity(yrnd(17))
   PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))


!----------------------------------------
! one loop correction to HTbar decay |
!----------------------------------------
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
    do BH1Hel=-1,BHMaxHel
    do BH2Hel=-1,BHMaxHel
        IF( XTOPDECAYS.EQ.1 ) THEN
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExt(1:4,5:9))
          call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExt(1:4,10:14))
        ENDIF
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call HTopBHDecay(ExtParticle(1),DKX_HTBH_1L1,BH1Hel,MomExt(1:4,5:9))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol

    enddo! iHel loop
    enddo! Hel loop
    enddo! Hel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_HtHtbgg = NLO_Res_Unpol




!----------------------------------------
! one loop correction to anti-top decay |
!----------------------------------------
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
    do BH1Hel=-1,BHMaxHel
    do BH2Hel=-1,BHMaxHel
        IF( XTOPDECAYS.EQ.1 ) THEN
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExt(1:4,5:9))
          call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExt(1:4,10:14))
        ENDIF
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call HTopBHDecay(ExtParticle(1),DKX_HTBH_1L2,BH1Hel,MomExt(1:4,5:9))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol

    enddo! iHel loop
    enddo! Hel loop
    enddo! Hel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_HtHtbgg = EvalCS_DKJ_1L_HtHtbgg + NLO_Res_Unpol



!----------------------------------------
! one loop correction to W- decay       |
!----------------------------------------
  IF( TOPDECAYS.EQ.4 .OR. TOPDECAYS.EQ.2 ) THEN
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
    do BH1Hel=-1,BHMaxHel
    do BH2Hel=-1,BHMaxHel
        IF( XTOPDECAYS.EQ.1 ) THEN
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExt(1:4,5:9))
          call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExt(1:4,10:14))
        ENDIF
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call HTopBHDecay(ExtParticle(1),DKX_HTBH_1L3,BH1Hel,MomExt(1:4,5:9))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol

    enddo! iHel loop
    enddo! Hel loop
    enddo! Hel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_HtHtbgg = EvalCS_DKJ_1L_HtHtbgg + NLO_Res_Unpol
  ENDIF






!----------------------------------------
! one loop correction on HT decay     |
! ----------------------------------------
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
    do BH1Hel=-1,BHMaxHel
    do BH2Hel=-1,BHMaxHel
        IF( XTOPDECAYS.EQ.1 ) THEN
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExt(1:4,5:9))
          call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExt(1:4,10:14))
        ENDIF
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call HTopBHDecay(ExtParticle(2),DKX_HTBH_1L1,BH1Hel,MomExt(1:4,10:14))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
    enddo! iHel loop
    enddo! Hel loop
    enddo! Hel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_HtHtbgg = EvalCS_DKJ_1L_HtHtbgg + NLO_Res_Unpol





!----------------------------------------
! one loop correction on top decay      |
!----------------------------------------
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
    do BH1Hel=-1,BHMaxHel
    do BH2Hel=-1,BHMaxHel
        IF( XTOPDECAYS.EQ.1 ) THEN
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExt(1:4,5:9))
          call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExt(1:4,10:14))
        ENDIF
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call HTopBHDecay(ExtParticle(2),DKX_HTBH_1L2,BH1Hel,MomExt(1:4,10:14))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
    enddo! iHel loop
    enddo! Hel loop
    enddo! Hel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_HtHtbgg = EvalCS_DKJ_1L_HtHtbgg + NLO_Res_Unpol




!----------------------------------------
! one loop correction on W+ decay       |
!----------------------------------------
  IF( TOPDECAYS.EQ.3 .OR. TOPDECAYS.EQ.2 ) THEN
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
    do BH1Hel=-1,BHMaxHel
    do BH2Hel=-1,BHMaxHel
        IF( XTOPDECAYS.EQ.1 ) THEN
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExt(1:4,5:9))
          call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExt(1:4,10:14))
        ENDIF
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call HTopBHDecay(ExtParticle(2),DKX_HTBH_1L3,BH1Hel,MomExt(1:4,10:14))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
    enddo! iHel loop
    enddo! Hel loop
    enddo! Hel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_HtHtbgg = EvalCS_DKJ_1L_HtHtbgg + NLO_Res_Unpol
  ENDIF




   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_DKJ_1L_HtHtbgg)
   enddo


   EvalCS_DKJ_1L_HtHtbgg = EvalCS_DKJ_1L_HtHtbgg/VgsWgt
   if( IsNan(EvalCS_DKJ_1L_HtHtbgg) ) then
        print *, "NAN:",EvalCS_DKJ_1L_HtHtbgg
        call printYrnd(yrnd(:))
        EvalCS_DKJ_1L_HtHtbgg = 0d0
        return
   endif

return
END FUNCTION











FUNCTION EvalCS_DKJ_1L_HtHtbqqb(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
implicit none
real(8) ::  EvalCS_DKJ_1L_HtHtbqqb
real(8) :: yRnd(1:VegasMxDim),VgsWgt
complex(8) :: NLO_Res_Pol,NLO_Res_Unpol,TreeResult(1:6),VirtResult(1:6),LO_Res_Unpol
integer :: iHel,jHel,iPrimAmp,jPrimAmp,nHel(1:2)
integer :: NBin(1:NumMaxHisto),NHisto,BHMaxHel,BH1Hel,BH2Hel
real(8) :: SpinAvg,ColorAvg,EHat,PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5,ISFac,RunFactor
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,PDFFac_a,PDFFac_b
real(8) :: MomExt(1:4,1:15)
real(8) :: pdf(-6:6,1:2)
logical :: applyPSCut
include 'vegas_common.f'


   EvalCS_DKJ_1L_HtHtbqqb= 0d0


   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if( EHat.le.2d0*m_HTop ) then
      EvalCS_DKJ_1L_HtHtbqqb = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)
   call EvalPhaseSpace_2to2HT(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt)! HTbar, HT
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
   ISFac = MomCrossing(MomExt)

   BHMaxHel = -1
   IF(XTOPDECAYS.EQ.1) THEN
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  BH top
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,4),yRnd(7:8),MomExt(1:4,10:11),PSWgt3)

      m_Top = m_SMTop
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd( 9:12),MomExt(1:4,7:9),PSWgt4)! bot lep neu
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(13:16),MomExt(1:4,12:14),PSWgt5)
      m_Top = m_HTop

      PSWgt = PSWgt * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
      BHMaxHel = +1     
   ELSEIF(XTOPDECAYS.EQ.2) THEN
      call Error("XTOPDECAYS.EQ.2 is not supported")
   ENDIF

   call Kinematics_TTbarETmiss(.false.,MomExt,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_DKJ_1L_HtHtbqqb = 0d0
      return
   endif

   call SetPDFs(eta1,eta2,MuFac,pdf)
   PDFFac_a = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
            + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
            + pdf(Bot_,1)*pdf(ABot_,2)
   PDFFac_b = pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
            + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
            + pdf(Bot_,2)*pdf(ABot_,1)
   PDFFac = PDFFac_a + PDFFac_b
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   RunFactor = RunAlphaS(NLOParam,MuRen)
   nHel(1:2) = getHelicity(yrnd(17))
   PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))


! -------------------------------------
! one loop correction to HTbar decay  |
!-------------------------------------
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
    do BH1Hel=-1,BHMaxHel
    do BH2Hel=-1,BHMaxHel
        IF( XTOPDECAYS.EQ.1 ) THEN
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExt(1:4,5:9))
          call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExt(1:4,10:14))
        ENDIF
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call HTopBHDecay(ExtParticle(1),DKX_HTBH_1L1,BH1Hel,MomExt(1:4,5:9))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol

    enddo! iHel loop
    enddo! Hel loop
    enddo! Hel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_HtHtbqqb = NLO_Res_Unpol




!----------------------------------------
! one loop correction to anti-top decay |
!----------------------------------------
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
    do BH1Hel=-1,BHMaxHel
    do BH2Hel=-1,BHMaxHel
        IF( XTOPDECAYS.EQ.1 ) THEN
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExt(1:4,5:9))
          call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExt(1:4,10:14))
        ENDIF
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call HTopBHDecay(ExtParticle(1),DKX_HTBH_1L2,BH1Hel,MomExt(1:4,5:9))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol

    enddo! iHel loop
    enddo! Hel loop
    enddo! Hel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_HtHtbqqb = EvalCS_DKJ_1L_HtHtbqqb + NLO_Res_Unpol



!----------------------------------------
! one loop correction to W- decay       |
!----------------------------------------
  IF( TOPDECAYS.EQ.4 .OR. TOPDECAYS.EQ.2 ) THEN
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
    do BH1Hel=-1,BHMaxHel
    do BH2Hel=-1,BHMaxHel
        IF( XTOPDECAYS.EQ.1 ) THEN
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExt(1:4,5:9))
          call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExt(1:4,10:14))
        ENDIF
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call HTopBHDecay(ExtParticle(1),DKX_HTBH_1L3,BH1Hel,MomExt(1:4,5:9))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol

    enddo! iHel loop
    enddo! Hel loop
    enddo! Hel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_HtHtbqqb = EvalCS_DKJ_1L_HtHtbqqb + NLO_Res_Unpol
  ENDIF






!----------------------------------------
! one loop correction on HT decay     |
! ----------------------------------------
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
    do BH1Hel=-1,BHMaxHel
    do BH2Hel=-1,BHMaxHel
        IF( XTOPDECAYS.EQ.1 ) THEN
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExt(1:4,5:9))
          call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExt(1:4,10:14))
        ENDIF
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call HTopBHDecay(ExtParticle(2),DKX_HTBH_1L1,BH1Hel,MomExt(1:4,10:14))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
    enddo! iHel loop
    enddo! Hel loop
    enddo! Hel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_HtHtbqqb = EvalCS_DKJ_1L_HtHtbqqb + NLO_Res_Unpol





!----------------------------------------
! one loop correction on top decay      |
!----------------------------------------
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
    do BH1Hel=-1,BHMaxHel
    do BH2Hel=-1,BHMaxHel
        IF( XTOPDECAYS.EQ.1 ) THEN
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExt(1:4,5:9))
          call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExt(1:4,10:14))
        ENDIF
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call HTopBHDecay(ExtParticle(2),DKX_HTBH_1L2,BH1Hel,MomExt(1:4,10:14))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
    enddo! iHel loop
    enddo! Hel loop
    enddo! Hel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_HtHtbqqb = EvalCS_DKJ_1L_HtHtbqqb + NLO_Res_Unpol




!----------------------------------------
! one loop correction on W+ decay       |
!----------------------------------------
  IF( TOPDECAYS.EQ.3 .OR. TOPDECAYS.EQ.2 ) THEN
    NLO_Res_Unpol = (0d0,0d0)
    do iHel=nHel(1),nHel(2)
    do BH1Hel=-1,BHMaxHel
    do BH2Hel=-1,BHMaxHel
        IF( XTOPDECAYS.EQ.1 ) THEN
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExt(1:4,5:9))
          call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExt(1:4,10:14))
        ENDIF
        call HelCrossing(Helicities(iHel,1:4))
        call SetPolarizations()
        
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call HTopBHDecay(ExtParticle(2),DKX_HTBH_1L3,BH1Hel,MomExt(1:4,10:14))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
    enddo! iHel loop
    enddo! Hel loop
    enddo! Hel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
    EvalCS_DKJ_1L_HtHtbqqb = EvalCS_DKJ_1L_HtHtbqqb + NLO_Res_Unpol
  ENDIF




   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_DKJ_1L_HtHtbqqb)
   enddo


   EvalCS_DKJ_1L_HtHtbqqb = EvalCS_DKJ_1L_HtHtbqqb/VgsWgt
   if( IsNan(EvalCS_DKJ_1L_HtHtbqqb) ) then
        print *, "NAN:",EvalCS_DKJ_1L_HtHtbqqb
        call printYrnd(yrnd(:))
        EvalCS_DKJ_1L_HtHtbqqb = 0d0
        return
   endif



END FUNCTION








FUNCTION EvalCS_StopWidth(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes
use ModMyRecurrence
use ModParameters
implicit none
real(8) ::  EvalCS_StopWidth,yRnd(1:VegasMxDim),VgsWgt,HOp(1:3)
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1)
integer :: iHel,jHel,kHel,iPrimAmp,jPrimAmp
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5,ISFac
real(8) :: MomExt(1:4,1:15),MomP(1:4,1:4),MomBoost(1:4),MomExtTd(1:4,1:9),pbDpg,ptDpg,ptDpb,TheDipole
logical :: applyPSCut
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,AccPoles
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2),xE,sigmaTot,beta
integer :: NBin(1:NumMaxHisto),NHisto,nHel(1:2)
real(8), parameter :: CF=4d0/3d0
include 'vegas_common.f'


   EvalCS_StopWidth = 0d0
  
! note: top decay is included but the results are independent of TopDK=1,2,3,4

   IF(XTOPDECAYS.EQ.3) THEN
      MomExt(1:4,3) = (/m_Stop,0d0,0d0,0d0/)
      ExtParticle(1)%Mom(1:4) = dcmplx( MomExt(1:4,3) )
      IF( CORRECTION.EQ.0 .OR.CORRECTION.EQ.4  ) THEN
          call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,3),yRnd(1:2),MomExt(1:4,5:6),PSWgt2)!  Chi top
          call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd(3:6),MomExt(1:4,7:9),PSWgt4)! bot lep neu
      ELSEIF( CORRECTION.EQ.5 ) THEN
          call EvalPhasespace_StopDK(ST_Chi0_T_G,MomExt(1:4,3),yRnd(1:5),MomExt(1:4,5:7),PSWgt2)!  Chi(5) top(6) glu(7->15)
          MomExt(1:4,15) = MomExt(1:4,7)
          call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd(6:9),MomExt(1:4,7:9),PSWgt4)!  b(7) el(8) nu(9)
      ENDIF
      PSWgt = PSWgt2
      if( PSWgt.eq.0d0 ) then
            SkipCounter = SkipCounter + 1
            return
      endif
   ENDIF

 

   PreFac = PSWgt * VgsWgt /(2d0*m_Stop) * (2d0*Ga_STop(0)*M_STop)! this removed the factors from STopDecay
   RunFactor = RunAlphaS(NLOParam,MuRen)

   IF( CORRECTION.EQ.0 ) THEN
      LO_Res_Unpol = (0d0,0d0)
      do iHel=-1,+1,2
      do jHel=-1,+1,2
          call STopDecay(ExtParticle(1),DKX_STChi0_LO,iHel,MomExt(1:4,5:9),HelTop=jHel)
          LO_Res_Pol = ExtParticle(1)%Pol(1)
          LO_Res_UnPol = LO_Res_UnPol +  LO_Res_Pol*dconjg( LO_Res_Pol )
      enddo!helicity loop
      enddo!helicity loop
!     normalization
      EvalCS_StopWidth = LO_Res_Unpol * PreFac




   ELSEIF( CORRECTION.EQ.4 ) THEN
      LO_Res_Unpol = (0d0,0d0)
      do iHel=-1,+1,2
      do jHel=-1,+1,2
          call STopDecay(ExtParticle(1),DKX_STChi0_LO,iHel,MomExt(1:4,5:9),HelTop=jHel)
          LO_Res_Pol = ExtParticle(1)%Pol(1)
          call STopDecay(ExtParticle(1),DKX_STChi0_1L1,iHel,MomExt(1:4,5:9),HelTop=jHel)
          NLO_Res_Pol(0) = ExtParticle(1)%Pol(1)

          LO_Res_UnPol = LO_Res_UnPol +  dreal( LO_Res_Pol*dconjg(NLO_Res_Pol(0)) )
      enddo!helicity loop
      enddo!helicity loop
!     normalization
      EvalCS_StopWidth = LO_Res_Unpol * PreFac



   ELSEIF( CORRECTION.EQ.5 ) THEN
      LO_Res_Unpol = (0d0,0d0)
      do iHel=-1,+1,2
      do jHel=-1,+1,2
      do kHel=-1,+1,2
          call STopDecay(ExtParticle(1),DKX_STChi0_RE1,iHel,MomExt(1:4,5:9),MomExt(1:4,15),HelGlu=kHel,HelTop=jHel)
          LO_Res_Pol = ExtParticle(1)%Pol(1)
          LO_Res_UnPol = LO_Res_UnPol +  LO_Res_Pol*dconjg( LO_Res_Pol )
      enddo!helicity loop
      enddo!helicity loop
      enddo!helicity loop
!     normalization
      EvalCS_StopWidth = LO_Res_Unpol * PreFac
! print *, LO_Res_Unpol * PreFac




      call WTransform3(MomExt(1:4,5:6),MomExt(1:4,15),MomExtTd(1:4,5:6),pbDpg,ptDpg,ptDpb)
      call EvalPhasespace_TopDK(T_B_W,MomExtTd(1:4,6),yRnd(6:9),MomExtTd(1:4,7:9),PSWgt4)    
      LO_Res_Unpol = (0d0,0d0)
      do iHel=-1,+1,2
      do jHel=-1,+1,2
          call STopDecay(ExtParticle(1),DKX_STChi0_LO,iHel,MomExtTd(1:4,5:9),HelTop=jHel)
          LO_Res_Pol = ExtParticle(1)%Pol(1)
          LO_Res_UnPol = LO_Res_UnPol +  LO_Res_Pol*dconjg( LO_Res_Pol )
      enddo!helicity loop
      enddo!helicity loop

      TheDipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg/ptDpg*( m_Stop**2+m_Top**2-(MomExt(1:4,5).dot.MomExt(1:4,5)) ) - (m_STop/ptDpg)**2 - (m_Top/pbDpg)**2 )
!       TheDipole = - 4d0*Pi * CF * ( 1d0/pbDpg/ptDpg*( m_Stop**2+m_Top**2-(MomExt(1:4,5).dot.MomExt(1:4,5)) ) - (m_STop/ptDpg)**2 - (m_Top/pbDpg)**2 )! alpha_s removed
! print *, LO_Res_Unpol * TheDipole * PreFac
! pause
!     normalization
      EvalCS_StopWidth = EvalCS_StopWidth + LO_Res_Unpol * TheDipole * PreFac



   ENDIF








   if( IsNan(EvalCS_StopWidth) ) then
        print *, "NAN:",EvalCS_StopWidth
        print *, yRnd(:)
        print *, NLO_Res_UnPol(0),NLO_Res_UnPol(1)
        print *, PSWgt , VgsWgt , PDFFac, sHatJacobi
        print *, eta1,eta2,((1d0-eta1)*xE+eta1),((1d0-eta2)*xE+eta2),MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,1:11)
        print *, "SKIP EVENT!!!!!"
        EvalCS_StopWidth = 0d0
        return
   endif

   EvalCounter = EvalCounter + 1
   EvalCS_StopWidth = EvalCS_StopWidth/VgsWgt

return
END FUNCTION













FUNCTION EvalCS_HTopWidth(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes
use ModMyRecurrence
use ModParameters
implicit none
real(8) ::  EvalCS_HTopWidth,yRnd(1:VegasMxDim),VgsWgt,HOp(1:3)
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),Amp,Amp2,Spi(1:4),BarSpi(1:4)
integer :: iHel,jHel,kHel,lHel,iPrimAmp,jPrimAmp
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5,ISFac
real(8) :: MomExt(1:4,1:15),MomP(1:4,1:4),MomBoost(1:4),MomExtTd(1:4,1:9),pbDpg,ptDpg,ptDpb,TheDipole
logical :: applyPSCut
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,AccPoles
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2),xE,sigmaTot,beta
integer :: NBin(1:NumMaxHisto),NHisto,nHel(1:2)
real(8), parameter :: CF=4d0/3d0
include 'vegas_common.f'


   EvalCS_HTopWidth = 0d0

   IF(XTOPDECAYS.EQ.1) THEN
      MomExt(1:4,3) = (/m_HTop,0d0,0d0,0d0/)
      ExtParticle(1)%Mom(1:4) = dcmplx( MomExt(1:4,3) )
      IF( CORRECTION.EQ.0 .OR.CORRECTION.EQ.4  ) THEN
          call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,3),yRnd(1:2),MomExt(1:4,5:6),PSWgt2)!  BH top
          call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd(3:6),MomExt(1:4,7:9),PSWgt4)! bot lep neu
      ELSEIF( CORRECTION.EQ.5 ) THEN
          call EvalPhasespace_HTopDK(HT_BH_T_G,MomExt(1:4,3),yRnd(1:5),MomExt(1:4,5:7),PSWgt2)!  BH(5) top(6) glu(7->15)
          MomExt(1:4,15) = MomExt(1:4,7)
          call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd(6:9),MomExt(1:4,7:9),PSWgt4)!  b(7) el(8) nu(9)
      ENDIF
      PSWgt = PSWgt2
   ENDIF
 

   PreFac = PSWgt * VgsWgt /(2d0*m_HTop) * 2d0*Ga_HTop(0)*M_HTop
   RunFactor = RunAlphaS(NLOParam,MuRen)


   IF( CORRECTION.EQ.0 ) THEN
      LO_Res_Unpol = (0d0,0d0)
      do iHel=-1,+1,1
      do jHel=-1,+1,2
      do kHel=-1,+1,2
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,iHel,MomExt(1:4,5:6),HelTop=jHel)
          Spi(1:4) = ExtParticle(1)%Pol(1:4)
          IF( ExtParticle(1)%PartType.lt.0d0 ) THEN
                call vBarSpi(ExtParticle(1)%Mom(1:4),m_HTop,kHel,BarSpi(1:4))
          ELSE
                call uSpi(ExtParticle(1)%Mom(1:4),m_HTop,kHel,BarSpi(1:4))
          ENDIF
          Amp = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_HTop)
          LO_Res_UnPol = LO_Res_UnPol +  Amp*dconjg(Amp) * 0.5d0
      enddo!helicity loop
      enddo!helicity loop
      enddo!helicity loop

!  normalization
   EvalCS_HTopWidth = LO_Res_Unpol * PreFac


   ELSEIF( CORRECTION.EQ.4 ) THEN


      LO_Res_Unpol = (0d0,0d0)
      do iHel=-1,+1,1
      do jHel=-1,+1,2
      do kHel=-1,+1,2
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,iHel,MomExt(1:4,5:6),HelTop=jHel)
          Spi(1:4) = ExtParticle(1)%Pol(1:4)
          IF( ExtParticle(1)%PartType.lt.0d0 ) THEN
                call vBarSpi(ExtParticle(1)%Mom(1:4),m_HTop,kHel,BarSpi(1:4))
          ELSE
                call uSpi(ExtParticle(1)%Mom(1:4),m_HTop,kHel,BarSpi(1:4))
          ENDIF
          Amp = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_HTop)

          call HTopBHDecay(ExtParticle(1),DKX_HTBH_1L1,iHel,MomExt(1:4,5:6),HelTop=jHel)
          Spi(1:4) = ExtParticle(1)%Pol(1:4)
          Amp2= psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_HTop)

          LO_Res_UnPol = LO_Res_UnPol +  dreal( Amp*dconjg(Amp2) ) * 0.5d0
      enddo!helicity loop
      enddo!helicity loop
      enddo!helicity loop

!  normalization
   EvalCS_HTopWidth = LO_Res_Unpol * PreFac




   ELSEIF( CORRECTION.EQ.5 ) THEN

      LO_Res_Unpol = (0d0,0d0)
      do iHel=-1,+1,1
      do jHel=-1,+1,2
      do kHel=-1,+1,2
      do lHel=-1,+1,2
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_RE1,iHel,MomExt(1:4,5:6),MomGlu=MomExt(1:4,15),GluHel=lHel,HelTop=jHel)
          Spi(1:4) = ExtParticle(1)%Pol(1:4)
          IF( ExtParticle(1)%PartType.lt.0d0 ) THEN
                call vBarSpi(ExtParticle(1)%Mom(1:4),m_HTop,kHel,BarSpi(1:4))
          ELSE
                call uSpi(ExtParticle(1)%Mom(1:4),m_HTop,kHel,BarSpi(1:4))
          ENDIF
          Amp = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_HTop)
          LO_Res_UnPol = LO_Res_UnPol +  Amp*dconjg(Amp) * 0.5d0
      enddo!helicity loop
      enddo!helicity loop
      enddo!helicity loop
      enddo!helicity loop
!     normalization
      EvalCS_HTopWidth = LO_Res_Unpol * PreFac


      call WTransform3(MomExt(1:4,5:6),MomExt(1:4,15),MomExtTd(1:4,5:6),pbDpg,ptDpg,ptDpb)
      call EvalPhasespace_TopDK(T_B_W,MomExtTd(1:4,6),yRnd(6:9),MomExtTd(1:4,7:9),PSWgt4) 
      if( dabs(ptDpg/m_HTop**2) .lt. 1d-6 )  then
          EvalCS_HTopWidth = 0d0
          return
      endif
      LO_Res_Unpol = (0d0,0d0)
      do iHel=-1,+1,1
      do jHel=-1,+1,2
      do kHel=-1,+1,2
          call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,iHel,MomExtTd(1:4,5:6),HelTop=jHel)
          Spi(1:4) = ExtParticle(1)%Pol(1:4)
          IF( ExtParticle(1)%PartType.lt.0d0 ) THEN
                call vBarSpi(ExtParticle(1)%Mom(1:4),m_HTop,kHel,BarSpi(1:4))
          ELSE
                call uSpi(ExtParticle(1)%Mom(1:4),m_HTop,kHel,BarSpi(1:4))
          ENDIF
          Amp = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_HTop)
          LO_Res_UnPol = LO_Res_UnPol +  Amp*dconjg(Amp) * 0.5d0
      enddo!helicity loop
      enddo!helicity loop
      enddo!helicity loop
      TheDipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg/ptDpg*( m_HTop**2+m_SMTop**2-(MomExt(1:4,5).dot.MomExt(1:4,5)) ) - (m_HTop/ptDpg)**2 - (m_SMTop/pbDpg)**2 )

! print *, "Sing",ptDpg/m_HTop**2, pbDpg/m_HTop**2
! print *, "Real",EvalCS_HTopWidth
! print *, "Dipo",LO_Res_Unpol*TheDipole*PreFac
! print *, "rati",(LO_Res_Unpol*TheDipole*PreFac)/EvalCS_HTopWidth + 1d0
! pause

!     normalization
      EvalCS_HTopWidth = EvalCS_HTopWidth + LO_Res_Unpol * TheDipole * PreFac



   ENDIF





   if( IsNan(EvalCS_HTopWidth) ) then
        print *, "NAN:",EvalCS_HTopWidth
        print *, yRnd(:)
        print *, NLO_Res_UnPol(0),NLO_Res_UnPol(1)
        print *, PSWgt , VgsWgt , PDFFac, sHatJacobi
        print *, eta1,eta2,((1d0-eta1)*xE+eta1),((1d0-eta2)*xE+eta2),MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,1:11)
        print *, "SKIP EVENT!!!!!"
        EvalCS_HTopWidth = 0d0
        return
   endif

   EvalCounter = EvalCounter + 1
   EvalCS_HTopWidth = EvalCS_HTopWidth/VgsWgt

return
END FUNCTION





END MODULE ModCrossSection_TTBETmiss


