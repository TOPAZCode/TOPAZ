MODULE ModCrossSection_TTBETmiss
use ModExoticDecay
implicit none

integer,private,parameter :: NumMaxHisto=20




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
use ModIntDipoles_GGTTBG
implicit none
real(8) ::  EvalCS_1L_HtHtbgg,yRnd(1:VegasMxDim),VgsWgt,HOp(1:3)
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1),FermionLoopPartAmp(7:8,-2:1)
integer :: iHel,jHel,kHel,iPrimAmp,jPrimAmp
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5,ISFac
real(8) :: MomExt(1:4,1:15),MomP(1:4,1:4)
logical :: applyPSCut
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,AccPoles
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2),xE
integer :: NBin(1:NumMaxHisto),NHisto,nHel(1:2),ParityFlip,BH1Hel,BH2Hel,BHMaxHel
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

  MuRen = 0.5d0*Ehat
  MuFac = MuRen

! ! run into threshold:
!    countr=countr+0.25d0m
!    EHAT = 2D0*M_TOP * (1d0+ 10d0**(-countr))
!    EHAT = 2D0*M_TOP * (1.0000000001d0)

   call EvalPhaseSpace_2to2HT(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt)! HTbar, HT
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
   ISFac = MomCrossing(MomExt)
   BHMaxHel = -1
   IF(XTOPDECAYS.EQ.1) THEN
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  BH top
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,4),yRnd(7:8),MomExt(1:4,10:11),PSWgt3)

      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd( 9:12),MomExt(1:4,7:9),PSWgt4)! bot lep neu
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(13:16),MomExt(1:4,12:14),PSWgt5)
      PSWgt = PSWgt * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
      BHMaxHel = +1     

   ELSEIF(XTOPDECAYS.EQ.2) THEN
      call EvalPhasespace_HTopDK(HT_A0_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  A0 top
      call EvalPhasespace_HTopDK(HT_A0_T,MomExt(1:4,4),yRnd(7:8),MomExt(1:4,10:11),PSWgt3)

      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd( 9:12),MomExt(1:4,7:9),PSWgt4)! bot lep neu
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(13:16),MomExt(1:4,12:14),PSWgt5)
      PSWgt = PSWgt * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
     
      call HTopA0Decay(ExtParticle(1),DKX_HTA0_LO,MomExt(1:4,5:9))
      call HTopA0Decay(ExtParticle(2),DKX_HTA0_LO,MomExt(1:4,10:14))
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
   nHel(1:2) = getHelicity(yrnd(17))
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
!	  call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))

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
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + Col1L_ttbgg(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(PrimAmps(jPrimAmp)%Result(-2:1)) )
      enddo
      enddo
      NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1) + NLO_Res_Pol(-2:1)


! ------------ fermionic loops --------------
      do iPrimAmp=7,10
          call SetKirill(PrimAmps(iPrimAmp))

          call PentCut(PrimAmps(iPrimAmp))
          call QuadCut(PrimAmps(iPrimAmp))
          call TripCut(PrimAmps(iPrimAmp))
          call DoubCut(PrimAmps(iPrimAmp))
          call SingCut(PrimAmps(iPrimAmp))
          call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)

          PrimAmps(iPrimAmp)%Result(-2:1) = -(0d0,1d0)*PrimAmps(iPrimAmp)%Result(-2:1) !minus is from closed fermion loop
!           call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,rdiv(2),rdiv(1))
!           call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv)
      enddo
      FermionLoopPartAmp(7,-2:1) = Nf_light*PrimAmps(7)%Result(-2:1) + PrimAmps(9)%Result(-2:1)
      FermionLoopPartAmp(8,-2:1) = Nf_light*PrimAmps(8)%Result(-2:1) + PrimAmps(10)%Result(-2:1)

      NLO_Res_Pol(-2:1) = (0d0,0d0)
      do jPrimAmp=7,8
      do iPrimAmp=1,NumBornAmps
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + Col1L_ttbgg(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result*dconjg(FermionLoopPartAmp(jPrimAmp,-2:1)) )
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


!   check Coulomb singularity
!   print *, "EHat,beta:",EHat,dsqrt(1d0-4d0*m_top**2/EHat**2)
!   CoulombFac = LO_Res_Unpol*(alpha_s*RunFactor)*DblPi/2d0/dsqrt(1d0-4d0*m_top**2/EHat**2) * 4d0/3d0/2d0
!   print *, "ratio",dble(NLO_Res_UnPol(0))/CoulombFac
!   print *, dsqrt(1d0-4d0*m_top**2/EHat**2), dble(NLO_Res_UnPol(0))/CoulombFac
!    STOP

   EvalCS_1L_HtHtbgg = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac


ELSEIF( Correction.EQ.3 ) THEN

   MomP(1:4,1) = MomExt(1:4,3)
   MomP(1:4,2) = MomExt(1:4,4)
   MomP(1:4,3) =-MomExt(1:4,1)
   MomP(1:4,4) =-MomExt(1:4,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt

   xE = yRnd(13)
   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
!   call EvalIntDipoles_GGTTBG(MomP(1:4,1:4),(/Mom(1:4,),/),xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
   EvalCS_1L_HtHtbgg = HOp(1)    * pdf(0,1)  * pdf(0,2)   &
                   + HOp(2)/xE * pdf_z(0,1)* pdf(0,2)   &
                   + HOp(3)/xE * pdf(0,1)  * pdf_z(0,2)
ENDIF


   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_HtHtbgg)
   enddo

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
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2)
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

  MuRen = 0.5d0*Ehat
  MuFac = MuRen


   call EvalPhaseSpace_2to2HT(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt)! HTbar, HT
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
   ISFac = MomCrossing(MomExt)

   BHMaxHel = -1
   IF(XTOPDECAYS.EQ.1) THEN
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  BH top
      call EvalPhasespace_HTopDK(HT_BH_T,MomExt(1:4,4),yRnd(7:8),MomExt(1:4,10:11),PSWgt3)

      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd( 9:12),MomExt(1:4,7:9),PSWgt4)! bot lep neu
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(13:16),MomExt(1:4,12:14),PSWgt5)
      PSWgt = PSWgt * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
      BHMaxHel = +1
   ELSEIF(XTOPDECAYS.EQ.2) THEN
      call EvalPhasespace_HTopDK(HT_A0_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  A0 top
      call EvalPhasespace_HTopDK(HT_A0_T,MomExt(1:4,4),yRnd(7:8),MomExt(1:4,10:11),PSWgt3)

      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd( 9:12),MomExt(1:4,7:9),PSWgt4)! bot lep neu
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(13:16),MomExt(1:4,12:14),PSWgt5)
      PSWgt = PSWgt * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
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
   nHel(1:2) = getHelicity(yrnd(17))
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
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      IF(XTOPDECAYS.eq.1) THEN
        call HTopBHDecay(ExtParticle(1),DKX_HTBH_LO,BH1Hel,MomExt(1:4,5:9))
        call HTopBHDecay(ExtParticle(2),DKX_HTBH_LO,BH2Hel,MomExt(1:4,10:14))
      ELSEIF(XTOPDECAYS.EQ.2) THEN
        call HTopA0Decay(ExtParticle(1),DKX_HTA0_LO,MomExt(1:4,5:9))
        call HTopA0Decay(ExtParticle(2),DKX_HTA0_LO,MomExt(1:4,10:14))
      ENDIF
      call SetPolarizations()
      do iPrimAmp=1,6
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
!           call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,rdiv(2),rdiv(1))
!           call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv)
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
      do iPrimAmp=5,6
          call SetKirill(PrimAmps(iPrimAmp))
          call PentCut(PrimAmps(iPrimAmp))
          call QuadCut(PrimAmps(iPrimAmp))
          call TripCut(PrimAmps(iPrimAmp))
          call DoubCut(PrimAmps(iPrimAmp))
          call SingCut(PrimAmps(iPrimAmp))
          call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
          PrimAmps(iPrimAmp)%Result(-2:1) = -(0d0,1d0)*PrimAmps(iPrimAmp)%Result(-2:1) !minus if from closed fermion loop
!           call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,rdiv(2),rdiv(1))
!           call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv)
      enddo

      FermionPartAmp(1,-2:1) = ( Nf_light*PrimAmps(PrimAmp2_1234)%Result(-2:1) +  PrimAmps(PrimAmp2m_1234)%Result(-2:1) )
      FermionPartAmp(2,-2:1) = -1d0/Nc * ( Nf_light*PrimAmps(PrimAmp2_1234)%Result(-2:1) + PrimAmps(PrimAmp2_1234)%Result(-2:1) )

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
ENDIF




IF( CORRECTION.EQ.0 ) THEN
!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * WidthExpansion
   EvalCS_1L_HtHtbqqb = LO_Res_Unpol * PreFac

ELSEIF( CORRECTION.EQ.1 ) THEN
!  overall normalization: (4*Pi)^eps/Gamma(1-eps)
!  CT contributions                           ! beta           !top WFRC
   NLO_Res_UnPol(-1) = NLO_Res_UnPol(-1) + (-11d0/3d0*3d0 - 3d0*4d0/3d0)*LO_Res_Unpol
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-3d0*4d0/3d0)*2d0*dlog(MuRen/m_top)*LO_Res_Unpol  ! finite log(mu2) contrib. from  top WFRC

   NLO_Res_UnPol_Ferm(-1) = NLO_Res_UnPol_Ferm(-1) + (2d0/3d0*Nf_light+2d0/3d0*Nf_heavy)*LO_Res_Unpol
   NLO_Res_UnPol_Ferm( 0) = NLO_Res_UnPol_Ferm( 0) + (2d0/3d0*Nf_heavy)*2d0*dlog(MuRen/m_top)*LO_Res_Unpol  ! finite log(mu2) contrib. from heavy flavor in alpha_s ren.

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
   MomP(1:4,1) = MomExt(1:4,3)
   MomP(1:4,2) = MomExt(1:4,4)
   MomP(1:4,3) =-MomExt(1:4,1)
   MomP(1:4,4) =-MomExt(1:4,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   ISFac = MomCrossing(MomExt)

IF( PROCESS.EQ.6 ) THEN
      xE = yRnd(13)
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
!      call EvalIntDipoles_QQBTTBG(MomP(1:4,1:4),MomDK(1:4,1:6),xE,HOp(1:3))

   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
   EvalCS_1L_HtHtbqqb = HOp(1)    * (pdf(Up_,1)*pdf(AUp_,2)+pdf(Dn_,1)*pdf(ADn_,2)+pdf(Chm_,1)*pdf(AChm_,2)+pdf(Str_,1)*pdf(AStr_,2)+pdf(Bot_,1)*pdf(ABot_,2) ) &
                    + HOp(2)/xE * (pdf_z(Up_,1)*pdf(AUp_,2)+pdf_z(Dn_,1)*pdf(ADn_,2)+pdf_z(Chm_,1)*pdf(AChm_,2)+pdf_z(Str_,1)*pdf(AStr_,2)+pdf_z(Bot_,1)*pdf(ABot_,2) ) &
                    + HOp(3)/xE * (pdf(Up_,1)*pdf_z(AUp_,2)+pdf(Dn_,1)*pdf_z(ADn_,2)+pdf(Chm_,1)*pdf_z(AChm_,2)+pdf(Str_,1)*pdf_z(AStr_,2)+pdf(Bot_,1)*pdf_z(ABot_,2) )

    call swapMom(MomP(1:4,3),MomP(1:4,4))
!      call EvalIntDipoles_QQBTTBG(MomP(1:4,1:4),MomDK(1:4,1:6),xE,HOp(1:3))
    HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
    EvalCS_1L_HtHtbqqb = EvalCS_1L_HtHtbqqb  &
                     + HOp(1)    * (pdf(Up_,2)*pdf(AUp_,1)+pdf(Dn_,2)*pdf(ADn_,1)+pdf(Chm_,2)*pdf(AChm_,1)+pdf(Str_,2)*pdf(AStr_,1)+pdf(Bot_,2)*pdf(ABot_,1) ) &
                     + HOp(2)/xE * (pdf_z(Up_,2)*pdf(AUp_,1)+pdf_z(Dn_,2)*pdf(ADn_,1)+pdf_z(Chm_,2)*pdf(AChm_,1)+pdf_z(Str_,2)*pdf(AStr_,1)+pdf_z(Bot_,2)*pdf(ABot_,1) ) &
                     + HOp(3)/xE * (pdf(Up_,2)*pdf_z(AUp_,1)+pdf(Dn_,2)*pdf_z(ADn_,1)+pdf(Chm_,2)*pdf_z(AChm_,1)+pdf(Str_,2)*pdf_z(AStr_,1)+pdf(Bot_,2)*pdf_z(ABot_,1) )


ELSEIF( PROCESS.EQ.3 ) THEN
      xE = yRnd(13)
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
!      call EvalIntDipoles_QGTTBQ(MomP(1:4,1:4),MomDK(1:4,1:6),xE,HOp(1:3))
    HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
    EvalCS_1L_HtHtbqqb = HOp(1)    * (pdf(Up_,1)*pdf(0,2)+pdf(Dn_,1)*pdf(0,2)+pdf(Chm_,1)*pdf(0,2)+pdf(Str_,1)*pdf(0,2)+pdf(Bot_,1)*pdf(0,2) ) &
                     + HOp(2)/xE * (pdf_z(Up_,1)*pdf(0,2)+pdf_z(Dn_,1)*pdf(0,2)+pdf_z(Chm_,1)*pdf(0,2)+pdf_z(Str_,1)*pdf(0,2)+pdf_z(Bot_,1)*pdf(0,2) ) &
                     + HOp(3)/xE * (pdf(Up_,1)*pdf_z(0,2)+pdf(Dn_,1)*pdf_z(0,2)+pdf(Chm_,1)*pdf_z(0,2)+pdf(Str_,1)*pdf_z(0,2)+pdf(Bot_,1)*pdf_z(0,2) )


    call swapMom(MomP(1:4,3),MomP(1:4,4))
!      call EvalIntDipoles_QGTTBQ(MomP(1:4,1:4),MomDK(1:4,1:6),xE,HOp(1:3))
    HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
    EvalCS_1L_HtHtbqqb = EvalCS_1L_HtHtbqqb  &
                    + HOp(1)    * (pdf(Up_,2)*pdf(0,1)+pdf(Dn_,2)*pdf(0,1)+pdf(Chm_,2)*pdf(0,1)+pdf(Str_,2)*pdf(0,1)+pdf(Bot_,2)*pdf(0,1) ) &
                    + HOp(2)/xE * (pdf_z(Up_,2)*pdf(0,1)+pdf_z(Dn_,2)*pdf(0,1)+pdf_z(Chm_,2)*pdf(0,1)+pdf_z(Str_,2)*pdf(0,1)+pdf_z(Bot_,2)*pdf(0,1) ) &
                    + HOp(3)/xE * (pdf(Up_,2)*pdf_z(0,1)+pdf(Dn_,2)*pdf_z(0,1)+pdf(Chm_,2)*pdf_z(0,1)+pdf(Str_,2)*pdf_z(0,1)+pdf(Bot_,2)*pdf_z(0,1) )


ELSEIF( PROCESS.EQ.4 ) THEN
      xE = yRnd(13)
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
!      call EvalIntDipoles_QGTTBQ(MomP(1:4,1:4),MomDK(1:4,1:6),xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
   EvalCS_1L_HtHtbqqb =HOp(1)    * (pdf(AUp_,1)*pdf(0,2)+pdf(ADn_,1)*pdf(0,2)+pdf(AChm_,1)*pdf(0,2)+pdf(AStr_,1)*pdf(0,2)+pdf(ABot_,1)*pdf(0,2) ) &
                    +HOp(2)/xE * (pdf_z(AUp_,1)*pdf(0,2)+pdf_z(ADn_,1)*pdf(0,2)+pdf_z(AChm_,1)*pdf(0,2)+pdf_z(AStr_,1)*pdf(0,2)+pdf_z(ABot_,1)*pdf(0,2) ) &
                    +HOp(3)/xE * (pdf(AUp_,1)*pdf_z(0,2)+pdf(ADn_,1)*pdf_z(0,2)+pdf(AChm_,1)*pdf_z(0,2)+pdf(AStr_,1)*pdf_z(0,2)+pdf(ABot_,1)*pdf_z(0,2) )


   call swapMom(MomP(1:4,3),MomP(1:4,4))
!   call EvalIntDipoles_QGTTBQ(MomP(1:4,1:4),MomDK(1:4,1:6),xE,HOp(1:3))
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


   EvalCS_1L_HtHtbqqb = EvalCS_1L_HtHtbqqb/VgsWgt

return
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
implicit none
real(8) ::  EvalCS_1L_ststbgg,yRnd(1:VegasMxDim),VgsWgt,HOp(1:3)
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
   nHel(1:2) = getHelicity(yrnd(17))
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
   do iHel=nHel(1),nHel(2)
      call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExt(1:4,5:9))
      call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,10:14))
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


!------------ 1 LOOP --------------
ELSEIF( Correction.EQ.1 ) THEN
!    do iHel=nHel(1),nHel(2)
   do iHel=1,1; print *, "fixed heli"
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
!           call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,2,rdiv(2),rdiv(1))
!           call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv)

! print *, "------"
! print *, "LO", BornAmps(iPrimAmp)%Result
! print *, "NLO",PrimAmps(iPrimAmp)%Result(-2)
! print *, "NLO",PrimAmps(iPrimAmp)%Result(-1)
! print *, "NLO",PrimAmps(iPrimAmp)%Result(0)+PrimAmps(iPrimAmp)%Result(1)
! print *, "ratio",PrimAmps(iPrimAmp)%Result(-2)/BornAmps(iPrimAmp)%Result
! print *, "ratio",PrimAmps(iPrimAmp)%Result(-1)/BornAmps(iPrimAmp)%Result
! pause

! !DEC$ IF (_QuadPrecImpr==1)
!           AccPoles = CheckPoles(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv(1:2))
!           if( AccPoles.gt.1d-4 ) then
!               coeff4_128(:,:) = qcmplx( coeff4(:,:) )
!               coeff5_128(:,:) = qcmplx( coeff5(:,:) )
!               call TripCut_128(PrimAmps(iPrimAmp))
!               call DoubCut_128(PrimAmps(iPrimAmp))
!               call SingCut_128(PrimAmps(iPrimAmp))
!               call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
!               call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen*2)
!               PrimAmps(iPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(iPrimAmp)%Result(-2:1)
!               AccPoles = CheckPoles(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv(1:2))
!               if( AccPoles.gt.5d-2 ) then
!                   print *, "SKIP",AccPoles
!                   call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
!                   EvalCS_1L_ttbgg = 0d0
!                   SkipCounter = SkipCounter + 1
!                   return
!               endif
!           endif
! !DEC$ ENDIF
      enddo

      NLO_Res_Pol(-2:1) = (0d0,0d0)
      do jPrimAmp=1,6
      do iPrimAmp=1,NumBornAmps
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + Col1L_ttbgg(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(PrimAmps(jPrimAmp)%Result(-2:1)) )
      enddo
      enddo
      NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1) + NLO_Res_Pol(-2:1)


! ------------ fermionic loops --------------
      do iPrimAmp=7,10; print *, "only primamp",iPrimAmp
          call SetKirill(PrimAmps(iPrimAmp))

          call PentCut(PrimAmps(iPrimAmp))
          call QuadCut(PrimAmps(iPrimAmp))
          call TripCut(PrimAmps(iPrimAmp))
          call DoubCut(PrimAmps(iPrimAmp))
          call SingCut(PrimAmps(iPrimAmp))
          call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)

          PrimAmps(iPrimAmp)%Result(-2:1) = -(0d0,1d0)*PrimAmps(iPrimAmp)%Result(-2:1) !minus is from closed fermion loop
!           call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,2,rdiv(2),rdiv(1))
!           call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv)
! print *, "------"
! print *, "LO", BornAmps(iPrimAmp)%Result
! print *, "NLO",PrimAmps(iPrimAmp)%Result(-2)
! print *, "NLO",PrimAmps(iPrimAmp)%Result(-1)
! print *, "NLO",PrimAmps(iPrimAmp)%Result(0)+PrimAmps(iPrimAmp)%Result(1)
! pause
      enddo
      FermionLoopPartAmp(7,-2:1) = Nf_light*PrimAmps(7)%Result(-2:1) + PrimAmps(9)%Result(-2:1)
      FermionLoopPartAmp(8,-2:1) = Nf_light*PrimAmps(8)%Result(-2:1) + PrimAmps(10)%Result(-2:1)

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
!  CT contributions                           ! beta        !top WFRC
   NLO_Res_UnPol(-1) = NLO_Res_UnPol(-1) + (-11d0/3d0*3d0 - 3d0*4d0/3d0 )*LO_Res_Unpol
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-3d0*4d0/3d0)*2d0*dlog(MuRen/m_top)*LO_Res_Unpol  ! finite log(mu) contrib. from  top WFRC

   NLO_Res_UnPol_Ferm(-1) = NLO_Res_UnPol_Ferm(-1) - (-2d0/3d0*Nf_light)*LO_Res_Unpol

   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-5d0/2d0*8d0/3d0 )*LO_Res_Unpol   ! finite contribution from top WFRC's
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + LO_Res_Unpol  ! shift alpha_s^DR --> alpha_s^MSbar


!  normalization
   LO_Res_Unpol = LO_Res_Unpol                         * ISFac * (alpha_s4Pi*RunFactor)**2
   NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1)           * ISFac * (alpha_s4Pi*RunFactor)**2 * alpha_sOver2Pi*RunFactor
   NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) * ISFac * (alpha_s4Pi*RunFactor)**2 * alpha_sOver2Pi*RunFactor

   EvalCS_1L_ststbgg = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac




ELSEIF( Correction.EQ.3 ) THEN

ENDIF


   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_ststbgg)
   enddo

   EvalCS_1L_ststbgg = EvalCS_1L_ststbgg/VgsWgt

return
END FUNCTION











FUNCTION EvalCS_1L_ststbqqb(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
implicit none
real(8) ::  EvalCS_1L_ststbqqb,yRnd(1:VegasMxDim),VgsWgt,HOp(1:3)
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1),FermionLoopPartAmp(1:2,-2:1)
integer :: iHel,jHel,kHel,iPrimAmp,jPrimAmp
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5,ISFac
real(8) :: MomExt(1:4,1:15),MomP(1:4,1:4),MomBoost(1:4)
logical :: applyPSCut
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac_a,PDFFac_b,PDFFac,AccPoles
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2),xE,sigmaTot,beta
integer :: NBin(1:NumMaxHisto),NHisto,nHel(1:2),npdf
real(8),parameter :: Nc=3d0
!include 'misc/global_import'
include 'vegas_common.f'

! yrnd(1)=0.8d0; yrnd(2)=0.65d0

   EvalCS_1L_ststbqqb = 0d0
   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if( EHat.le.2d0*m_STop ) then
      EvalCS_1L_ststbqqb = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to2Stops(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt)! AStop, Stop
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))

   IF(XTOPDECAYS.EQ.3) THEN
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,3),yRnd(5:6),MomExt(1:4,5:6),PSWgt2)!  Chi top
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,4),yRnd(7:8),MomExt(1:4,10:11),PSWgt3)
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd( 9:12),MomExt(1:4,7:9),PSWgt4)! bot lep neu
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(13:16),MomExt(1:4,12:14),PSWgt5)
      PSWgt = PSWgt * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
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

!------------ LO ----------------
IF( Correction.EQ.0 ) THEN
   do iHel=nHel(1),nHel(2)
      call STopDecay(ExtParticle(1),DKX_STChi0_LO,Helicities(iHel,5),MomExt(1:4,5:9))
      call STopDecay(ExtParticle(2),DKX_STChi0_LO,Helicities(iHel,6),MomExt(1:4,10:14))

      call HelCrossing(Helicities(iHel,1:4))
      call SetPolarizations()
      do iPrimAmp=1,NumBornAmps
! print *,"EVAL PRIMAMP",iPrimAmp
          call EvalTree(BornAmps(iPrimAmp))
! print *, BornAmps(iPrimAmp)%Result
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
      do iPrimAmp=1,2!         4

print *, "hel",ihel,"primamp no",iPrimAmp
print *, "ordering ",PrimAmps(iPrimAmp)%ExtLine(1:4)
          call SetKirill(PrimAmps(iPrimAmp))
          call PentCut(PrimAmps(iPrimAmp))
          call QuadCut(PrimAmps(iPrimAmp))
          call TripCut(PrimAmps(iPrimAmp))
          call DoubCut(PrimAmps(iPrimAmp))
          call SingCut(PrimAmps(iPrimAmp))
          call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)

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


!       BosonicPartAmp(1,-2:1) =  &
!                              + Nc * PrimAmps(PrimAmp1_1234)%Result(-2:1) &
!                              - 2d0/Nc*( PrimAmps(PrimAmp1_1234)%Result(-2:1) + PrimAmps(PrimAmp1_1243)%Result(-2:1) ) &
!                              - 1d0/Nc * ( PrimAmps(PrimAmp3_1432)%Result(-2:1) + PrimAmps(PrimAmp4_1234)%Result(-2:1) )
! 
!       BosonicPartAmp(2,-2:1) = PrimAmps(PrimAmp1_1243)%Result(-2:1)  &
!                              + 1d0/Nc**2 * ( PrimAmps(PrimAmp1_1234)%Result(-2:1) + PrimAmps(PrimAmp1_1243)%Result(-2:1) )  &
!                              + 1d0/Nc**2 * ( PrimAmps(PrimAmp3_1432)%Result(-2:1) + PrimAmps(PrimAmp4_1234)%Result(-2:1) )
! 
! 
!       NLO_Res_Pol(-2:1) = (0d0,0d0)
!       do jPrimAmp=1,2
!       do iPrimAmp=1,NumBornAmps
!           NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + Col1L_ttbqqb(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(BosonicPartAmp(jPrimAmp,-2:1)) )
!       enddo
!       enddo
!       NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1) + NLO_Res_Pol(-2:1)*PDFFac

! ------------ fermionic loops --------------
!       do iPrimAmp=1,2; print *, "Evaluate primamp only",iPrimAmp
!           call SetKirill(PrimAmps(iPrimAmp))
!           call PentCut(PrimAmps(iPrimAmp))
!           call QuadCut(PrimAmps(iPrimAmp))
!           call TripCut(PrimAmps(iPrimAmp))
!           call DoubCut(PrimAmps(iPrimAmp))
!           call SingCut(PrimAmps(iPrimAmp))
!           call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
!           PrimAmps(iPrimAmp)%Result(-2:1) = -(0d0,1d0)*PrimAmps(iPrimAmp)%Result(-2:1) !minus from closed fermion loop
! !           call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,rdiv(2),rdiv(1))
! !           call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv)
!       enddo

!       FermionLoopPartAmp(1,-2:1) =           ( Nf_light*PrimAmps(5)%Result(-2:1) + PrimAmps(6)%Result(-2:1) ) 
!       FermionLoopPartAmp(2,-2:1) = -1d0/Nc * ( Nf_light*PrimAmps(5)%Result(-2:1) + PrimAmps(6)%Result(-2:1) )

      NLO_Res_Pol(-2:1) = (0d0,0d0)
      do jPrimAmp=1,2
      do iPrimAmp=1,NumBornAmps
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + Col1L_ttbqqb(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(FermionLoopPartAmp(jPrimAmp,-2:1)) )
      enddo
      enddo
      NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) + NLO_Res_Pol(-2:1)*PDFFac

! print *, "hel",ihel
! print *, "LO",BornAmps(1)%Result
! print *, "NLO",PrimAmps(2)%Result(-2:1)
! print *, "ratio",PrimAmps(2)%Result(-1)/BornAmps(1)%Result
! pause

   enddo!helicity loop
ENDIF! Correction loop
  enddo! npdf loop
  call swapMom(MomExt(1:4,1),MomExt(1:4,2))   ! swap back to original order, for ID below



IF( Correction.EQ.0 ) THEN
!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * WidthExpansion
   EvalCS_1L_ststbqqb = LO_Res_Unpol * PreFac

ELSEIF( Correction.EQ.1 ) THEN

ELSEIF( Correction.EQ.3 ) THEN

ENDIF


   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_ststbqqb)
   enddo

   EvalCS_1L_ststbqqb = EvalCS_1L_ststbqqb/VgsWgt

return
END FUNCTION












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
   call EvalPhasespace_2toN(3,EHat,yRnd(3:7),MomExt(1:4,1:5),(/m_Stop,m_Stop,0d0,0d0,0d0/),PSWgt)! AStop, Stop
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))


 
   call Kinematics_TTbarETmiss(.false.,MomExt,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_1L_ststbgggg = 0d0
      return
   endif

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

      call HelCrossing(Helicities(iHel,1:4))
      call SetPolarizations()
      do iPrimAmp=1,NumBornAmps
! print *,"EVAL PRIMAMP",iPrimAmp
          call EvalTree(BornAmps(iPrimAmp))
! print *, BornAmps(iPrimAmp)%Result
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
use ModDipoles_DKP_GGSTSTBG
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


   EvalCS_Real_ststbggg = 0d0
   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if( EHat.le.2d0*m_STop ) then
      EvalCS_Real_ststbggg = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to3Stops(EHat,yRnd(3:7),MomExt(1:4,1:5),PSWgt)! Glu AStop, Stop
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
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,3),yRnd(8:9),MomExt(1:4,6:7),PSWgt2)!  Chi top
      call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,4),yRnd(10:11),MomExt(1:4,11:12),PSWgt3)
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRnd(12:15),MomExt(1:4,8:10),PSWgt4)! bot lep neu
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRnd(16:19),MomExt(1:4,13:15),PSWgt5)
      PSWgt = PSWgt * PSWgt2*PSWgt3 * PSWgt4*PSWgt5
   ENDIF
 
   call SetPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   RunFactor = RunAlphaS(2,MuRen)

   call Kinematics_TTbarETmiss(.true.,MomExt,(/4,5,6,11,7,12,8,9,10,13,14,15,3/),applyPSCut,NBin)
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
     call EvalDipoles_GGSTSTBG((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:19),PreFac,DipoleResult)
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









END MODULE ModCrossSection_TTBETmiss


