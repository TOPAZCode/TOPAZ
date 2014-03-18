MODULE ModCrossSection_TTB
use ModTopDecay
implicit none

integer,private,parameter :: NumMaxHisto=45

contains


FUNCTION EvalCS_1L_ttbgg_MPI(yRnd,VgsWgt,res)
implicit none
integer :: EvalCS_1L_ttbgg_MPI
real(8) ::  yRnd(*),res(*),VgsWgt

res(1) = EvalCS_1L_ttbgg(yRnd,VgsWgt)
EvalCS_1L_ttbgg_MPI=0
RETURN
END FUNCTION


FUNCTION EvalCS_1L_ttbgg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModIntDipoles_GGTTBG
use ModIntDipoles_GGTTBG_noDK
use ModJPsiFrag
implicit none
real(8) ::  EvalCS_1L_ttbgg,yRnd(1:VegasMxDim),VgsWgt,IOp(-2:0),HOp(1:3)
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1),FermionLoopPartAmp(7:8,-2:1)
complex(8) ::  Msq_T_BWENU,M_T_BWENU,Spi(1:4),BarSpi(1:4)
integer :: iHel,jHel,kHel,iPrimAmp,jPrimAmp,WHel
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,ISFac,xFrag
real(8) :: MomExt(1:4,1:4),MomDK(1:4,1:6),MomP(1:4,1:4),MomJPsi(1:4)
logical :: applyPSCut
real(8) :: MG_MOM(0:3,1:4),MG_MOM2(0:3,1:8)
real(8) :: MadGraph_tree,GG_TTB,Msq_T_BW,Msq_W_ENU
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,AccPoles
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2),xE,r_sc
integer :: NBin(1:NumMaxHisto),NHisto,ParityFlip,nHel(1:2),NRndHel
real(8) :: ThresholdCutOff = 1.0d0
! real(8) :: countr=0d0,CoulombFac
include 'misc/global_import'
include 'vegas_common.f'




IF( CORRECTION.EQ.1 ) ThresholdCutOff = 1.01d0
IF( TopDecays.ge.1 ) THEN
  ParityFlip=1
ELSE
  ParityFlip=2
ENDIF

  EvalCS_1L_ttbgg = 0d0

  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top * ThresholdCutOff ) then
      EvalCS_1L_ttbgg = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

! ! run into threshold:
!    countr=countr+0.25d0
!    EHAT = 2D0*M_TOP * (1d0+ 10d0**(-countr))
!    EHAT = 2D0*M_TOP * (1.0000000001d0)

   call EvalPhaseSpace_2to2(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt)
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
   ISFac = MomCrossing(MomExt)

!    muren= dsqrt( m_top**2 + get_pt(MomExt(1:4,3))**2 )
!    mufac= dsqrt( m_top**2 + get_pt(MomExt(1:4,3))**2 )

   NRndHel=5
IF( TOPDECAYS.NE.0 ) THEN
   call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:8),.false.,MomDK(1:4,1:3),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(9:12),.false.,MomDK(1:4,4:6),PSWgt3)
   PSWgt = PSWgt * PSWgt2*PSWgt3

   NRndHel=13
   IF( TOPDECAYS.GE.1 ) THEN
! top decay with spin correlations
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
   ELSE
! spherical decay of ATop
!       Msq_T_BW = dsqrt(2d0)*GF*m_top**4*(1d0-m_w**2/m_top**2)*(1d0+2d0*m_w**2/m_top**2)  ! sq.mat.el. for T -> b W
!       Msq_W_ENU = dsqrt(2d0)*GF*m_W**4/3d0*4d0  ! sq.mat.el. for W -> l nu
!       PSWgt = PSWgt * Msq_T_BW/(2d0*Ga_Top(0)*m_Top)
!       PSWgt = PSWgt * Msq_W_ENU/(2d0*Ga_W*m_W)
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      Spi(1:4) = ExtParticle(1)%Pol(1:4)
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,-1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0

      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,+1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0

      PSWgt = PSWgt * Msq_T_BWENU

! spherical decay of Top
!       PSWgt = PSWgt * Msq_T_BW/(2d0*Ga_Top(0)*m_Top)
!       PSWgt = PSWgt * Msq_W_ENU/(2d0*Ga_W*m_W)
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
      BarSpi(1:4) = ExtParticle(2)%Pol(1:4)
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,-1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0

      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,+1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
      PSWgt = PSWgt * Msq_T_BWENU
   ENDIF

!    MuFac = get_pt(MomDK(1:4,2))
!    MuRen = get_pt(MomDK(1:4,2))

   IF( TOPDECAYS.EQ.5 ) THEN
          if(Correction.EQ.3 ) then
              xFrag= yRnd(14)
          else
              xFrag= yRnd(13)
          endif
          PSWgt = PSWgt * FF(xFrag)
   ELSEIF( TOPDECAYS.EQ.6 ) THEN
          if(Correction.EQ.3 ) then
              xFrag= yRnd(14)
          else
              xFrag= yRnd(13)
          endif
          PSWgt = PSWgt * FF(xFrag)
   ENDIF
ENDIF

   call Kinematics_TTBAR(.false.,MomExt,MomDK,applyPSCut,NBin,xJPsiFrag=xFrag)
   if( applyPSCut ) then
      EvalCS_1L_ttbgg = 0d0
      return
   endif

   call InitCurrCache()
   call SetPropagators()

   call SetPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   RunFactor = RunAlphaS(NLOParam,MuRen)
   nHel(1:2) = getHelicity(yrnd(NRndHel))
   PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))

   LO_Res_Unpol             = (0d0,0d0)
   NLO_Res_Unpol(-2:1)      = (0d0,0d0)
   NLO_Res_Unpol_Ferm(-2:1) = (0d0,0d0)

!------------ LO --------------
IF( Correction.EQ.0 ) THEN
   do iHel=nHel(1),nHel(2)
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

!------------ 1 LOOP --------------
ELSEIF( Correction.EQ.1 ) THEN
   do iHel=nHel(1),nHel(2)
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
!              coeff4_128(:,:) = qcmplx( coeff4(:,:) )
              coeff5_128(:,:) = qcmplx( coeff5(:,:) )
!print *, "before QUAD PREC",AccPoles
              call QuadCut_128(PrimAmps(iPrimAmp))
              call TripCut_128(PrimAmps(iPrimAmp))
              call DoubCut_128(PrimAmps(iPrimAmp))
              call SingCut_128(PrimAmps(iPrimAmp))
!PrimAmps(iPrimAmp)%UCuts(5)%Coeff(:,:)= dcmplx(PrimAmps(iPrimAmp)%UCuts(5)%Coeff_128(:,:) ) 
!PrimAmps(iPrimAmp)%UCuts(4)%Coeff(:,:)= dcmplx(PrimAmps(iPrimAmp)%UCuts(4)%Coeff_128(:,:) ) 
!PrimAmps(iPrimAmp)%UCuts(3)%Coeff(:,:)= dcmplx(PrimAmps(iPrimAmp)%UCuts(3)%Coeff_128(:,:) ) 
!PrimAmps(iPrimAmp)%UCuts(2)%Coeff(:,:)= dcmplx(PrimAmps(iPrimAmp)%UCuts(2)%Coeff_128(:,:) ) 
!PrimAmps(iPrimAmp)%UCuts(1)%Coeff(:,:)= dcmplx(PrimAmps(iPrimAmp)%UCuts(1)%Coeff_128(:,:) ) 
              call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
              call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen*2)
              PrimAmps(iPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(iPrimAmp)%Result(-2:1)
              AccPoles = CheckPoles(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv(1:2))
!print *, "after QUAD PREC",AccPoles
!pause
              if( AccPoles.gt.5d-2 ) then
                  print *, "SKIP",AccPoles
                  call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
                  EvalCS_1L_ttbgg = 0d0
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
      do iPrimAmp=7,10
          call SetKirill(PrimAmps(iPrimAmp))
          call PentCut(PrimAmps(iPrimAmp))
          call QuadCut(PrimAmps(iPrimAmp))
          call TripCut(PrimAmps(iPrimAmp))
          call DoubCut(PrimAmps(iPrimAmp))
          call SingCut(PrimAmps(iPrimAmp))
          call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)

          PrimAmps(iPrimAmp)%Result(-2:1) = -(0d0,1d0)*PrimAmps(iPrimAmp)%Result(-2:1) !minus is from closed fermion loop
          call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,2,rdiv(2),rdiv(1))
!           call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp-6),rdiv,(/EHat/))
      enddo
      FermionLoopPartAmp(7,-2:1) = Nf_light*PrimAmps(7)%Result(-2:1) + PrimAmps(9)%Result(-2:1)
      FermionLoopPartAmp(8,-2:1) = Nf_light*PrimAmps(8)%Result(-2:1) + PrimAmps(10)%Result(-2:1)

      NLO_Res_Pol(-2:1) = (0d0,0d0)
      do jPrimAmp=7,8
      do iPrimAmp=1,NumBornAmps
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + ParityFlip*Col1L_ttbgg(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result*dconjg(FermionLoopPartAmp(jPrimAmp,-2:1)) )
      enddo
      enddo
      NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) + NLO_Res_Pol(-2:1)


   enddo! helicity loop
ENDIF


IF( Correction.EQ.0 ) THEN
!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * WidthExpansion
   EvalCS_1L_ttbgg = LO_Res_Unpol * PreFac

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
!   22 continue

!   check Coulomb singularity
!   print *, "EHat,beta:",EHat,dsqrt(1d0-4d0*m_top**2/EHat**2)
!   CoulombFac = LO_Res_Unpol*(alpha_s*RunFactor)*DblPi/2d0/dsqrt(1d0-4d0*m_top**2/EHat**2) * 4d0/3d0/2d0
!   print *, "ratio",dble(NLO_Res_UnPol(0))/CoulombFac
!   print *, dsqrt(1d0-4d0*m_top**2/EHat**2), dble(NLO_Res_UnPol(0))/CoulombFac
!    STOP

   EvalCS_1L_ttbgg = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac


!    xE = 0.123d0
!    MomP(1:4,1) = MomExt(1:4,3)
!    MomP(1:4,2) = MomExt(1:4,4)
!    MomP(1:4,3) =-MomExt(1:4,1)
!    MomP(1:4,4) =-MomExt(1:4,2)
!    PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
!    call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
!    call EvalIntDipoles_GGTTBG(MomP(1:4,1:4),MomDK(1:4,1:6),xE,HOp(1:3))
!    HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
!    print *, ( NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2) ) * PreFac
!    print *, ( NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1) ) * PreFac
!    print *, "1L check",NLO_Res_UnPol(-2)/(alpha_sOver2Pi*RunFactor)/LO_Res_Unpol
!    print *, "res1=", HOp(1)
!    print *, "res2=", HOp(2)
!    print *, "res3=", HOp(3)
!    stop

ELSEIF( Correction.EQ.3 ) THEN

   MomP(1:4,1) = MomExt(1:4,3)
   MomP(1:4,2) = MomExt(1:4,4)
   MomP(1:4,3) =-MomExt(1:4,1)
   MomP(1:4,4) =-MomExt(1:4,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   IF( TOPDECAYS.GE.1 ) THEN
      xE = yRnd(13)
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
      call EvalIntDipoles_GGTTBG(MomP(1:4,1:4),MomDK(1:4,1:6),xE,HOp(1:3))
   ELSEIF( TOPDECAYS.LE.0 ) THEN
      xE = yRnd(5)
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
      call EvalIntDipoles_GGTTBG_noDK(MomP(1:4,1:4),xE,HOp(1:3))
   ENDIF
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
   EvalCS_1L_ttbgg = HOp(1)    * pdf(0,1)  * pdf(0,2)   &
                   + HOp(2)/xE * pdf_z(0,1)* pdf(0,2)   &
                   + HOp(3)/xE * pdf(0,1)  * pdf_z(0,2)
ENDIF


   if( IsNan(EvalCS_1L_ttbgg) ) then
        print *, "NAN:",EvalCS_1L_ttbgg
        print *, yRnd(:)
        print *, NLO_Res_UnPol(0),NLO_Res_UnPol(1),NLO_Res_UnPol_Ferm(0),NLO_Res_UnPol_Ferm(1),IOp(0)
        print *, PSWgt , VgsWgt , PDFFac, sHatJacobi
        print *, eta1,eta2,((1d0-eta1)*xE+eta1),((1d0-eta2)*xE+eta2),MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,1)
        print *, MomExt(1:4,2)
        print *, MomExt(1:4,3)
        print *, MomExt(1:4,4)
        print *, MomDK(1:4,1)
        print *, MomDK(1:4,2)
        print *, MomDK(1:4,3)
        print *, MomDK(1:4,4)
        print *, MomDK(1:4,5)
        print *, MomDK(1:4,6)
        print *, "SKIP EVENT!!!!!"
        EvalCS_1L_ttbgg = 0d0
        return
   endif



IF( ObsSet.EQ.7 ) THEN
   do NHisto=1,3
      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_ttbgg)
   enddo
   call intoHisto(4,1,MInv_LB*EvalCS_1L_ttbgg)
   call intoHisto(5,1,MInv_LB**2*EvalCS_1L_ttbgg)
   call intoHisto(6,1,xFrag*EvalCS_1L_ttbgg)
ELSE
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_ttbgg)
   enddo

! if(EvalCS_1L_ttbgg.lt.0d0) then
!     print *, EvalCS_1L_ttbgg
!     print *, pdf(0,1) , pdf(0,2),fluxFac , sHatJacobi , PSWgt , VgsWgt , PDFFac,LO_Res_Unpol, ISFac
!     pause
! endif


ENDIF

   EvalCS_1L_ttbgg = EvalCS_1L_ttbgg/VgsWgt
   EvalCounter = EvalCounter + 1

   if( IsNan(EvalCS_1L_ttbgg) ) then
        print *, "NAN:",EvalCS_1L_ttbgg
        print *, yRnd(:)
        print *, NLO_Res_UnPol(0),NLO_Res_UnPol(1),NLO_Res_UnPol_Ferm(0),NLO_Res_UnPol_Ferm(1),IOp(0)
        print *, PSWgt , VgsWgt , PDFFac, sHatJacobi
        print *, eta1,eta2,((1d0-eta1)*xE+eta1),((1d0-eta2)*xE+eta2),MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,1)
        print *, MomExt(1:4,2)
        print *, MomExt(1:4,3)
        print *, MomExt(1:4,4)
        print *, MomDK(1:4,1)
        print *, MomDK(1:4,2)
        print *, MomDK(1:4,3)
        print *, MomDK(1:4,4)
        print *, MomDK(1:4,5)
        print *, MomDK(1:4,6)
   endif

! !      MADGRAPH CHECK: gg->ttb
!        MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!        MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!        MG_MOM(0:3,3) = MomExt(1:4,3)*100d0
!        MG_MOM(0:3,4) = MomExt(1:4,4)*100d0
!        call coupsm(0)
!        call SGG_TTB(MG_MOM,MadGraph_tree)
!        MadGraph_tree= MadGraph_tree
!        print *, ""
!        print *, "My tree:         ", LO_Res_Unpol
!        print *, "MadGraph hel.amp:", MadGraph_tree
!        print *, "MG/ME ratio: ", MadGraph_tree/dble(LO_Res_Unpol)
!        PAUSE
return
END FUNCTION





FUNCTION EvalCS_1L_ttbqqb_MPI(yRnd,VgsWgt,res)
implicit none
integer :: EvalCS_1L_ttbqqb_MPI
real(8) ::  yRnd(*),res(*),VgsWgt

res(1) = EvalCS_1L_ttbqqb(yRnd,VgsWgt)
EvalCS_1L_ttbqqb_MPI=0
RETURN
END FUNCTION



FUNCTION EvalCS_1L_ttbqqb(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModIntDipoles_QQBTTBG
use ModIntDipoles_QQBTTBG_noDK
use ModIntDipoles_QGTTBQ
use ModIntDipoles_QGTTBQ_noDK
use ModJPsiFrag
implicit none
real(8) ::  EvalCS_1L_ttbqqb,yRnd(1:VegasMxDim),VgsWgt,xE,IOp(-2:0),HOp(1:3)
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1)
complex(8) :: BosonicPartAmp(1:2,-2:1),FermionPartAmp(1:2,-2:1)
complex(8) ::  Msq_T_BWENU,M_T_BWENU,Spi(1:4),BarSpi(1:4)
integer :: iHel,iPrimAmp,jPrimAmp
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,ISFac,xFrag,SpinDecorr
real(8) :: MomExt(1:4,1:NumExtParticles),MomDK(1:4,1:6),MomP(1:4,1:4),MomJPsi(1:4)
logical :: applyPSCut
real(8) :: MG_MOM(0:3,1:NumExtParticles),MomTmp(1:4)
real(8) :: MadGraph_tree,UUB_TTB,Msq_T_BW,Msq_W_ENU
real(8),parameter :: Nc=3d0
real(8) :: tau,eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac_a,PDFFac_b,PDFFac,AccPoles
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2),r_sc
integer :: NHisto,NBin(1:NumMaxHisto),ParityFlip,npdf,nHel(1:2),NRndHel
real(8) :: ThresholdCutOff = 1.0d0
include 'misc/global_import'
include "vegas_common.f"



IF( CORRECTION.EQ.1 ) ThresholdCutOff = 1.01d0
IF( TopDecays.ge.1 ) THEN
  ParityFlip=1
ELSE
  ParityFlip=2
ENDIF

   EvalCS_1L_ttbqqb = 0d0

  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top * ThresholdCutOff ) then
      EvalCS_1L_ttbqqb = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)


   call EvalPhaseSpace_2to2(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt)
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
   NRndHel=5
   SpinDecorr = 1d0
IF( TOPDECAYS.NE.0 ) THEN
   call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:8),.false.,MomDK(1:4,1:3),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(9:12),.false.,MomDK(1:4,4:6),PSWgt3)
   PSWgt = PSWgt * PSWgt2*PSWgt3
   NRndHel=13

   IF( TOPDECAYS.EQ.5 ) THEN
          if(Correction.EQ.3 ) then
              xFrag= yRnd(14)
          else
              xFrag= yRnd(13)
          endif
          PSWgt = PSWgt * FF(xFrag)
   ELSEIF( TOPDECAYS.EQ.6 ) THEN
          if(Correction.EQ.3 ) then
              xFrag= yRnd(14)
          else
              xFrag= yRnd(13)
          endif
          PSWgt = PSWgt * FF(xFrag)
   ENDIF
ENDIF


   call Kinematics_TTBAR(.false.,MomExt,MomDK,applyPSCut,NBin,xJPsiFrag=xFrag)
   if( applyPSCut ) then
      EvalCS_1L_ttbqqb = 0d0
      return
   endif

!    muren= dsqrt( m_top**2 + get_pt(MomExt(1:4,3))**2 )
!    mufac= dsqrt( m_top**2 + get_pt(MomExt(1:4,3))**2 )

   call setPDFs(eta1,eta2,MuFac,pdf)
   PDFFac_a = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
            + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
            + pdf(Bot_,1)*pdf(ABot_,2)
   PDFFac_b = pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
            + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
            + pdf(Bot_,2)*pdf(ABot_,1)

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   RunFactor = RunAlphaS(NLOParam,MuRen)
   nHel(1:2) = getHelicity(yrnd(NRndHel))
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

IF( TOPDECAYS.GE.1 ) THEN
!     top decays with spin correlations
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
ELSEIF( TOPDECAYS.LE.-1 ) THEN
! spherical decay of ATop
!      Msq_T_BW = dsqrt(2d0)*GF*m_top**4*(1d0-m_w**2/m_top**2)*(1d0+2d0*m_w**2/m_top**2)  ! sq.mat.el. for T -> b W
!      Msq_W_ENU = dsqrt(2d0)*GF*m_W**4/3d0*4d0  ! sq.mat.el. for W -> l nu
!      PSWgt = PSWgt * Msq_T_BW/(2d0*Ga_Top(0)*m_Top)
!      PSWgt = PSWgt * Msq_W_ENU/(2d0*Ga_W*m_W)
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      Spi(1:4) = ExtParticle(1)%Pol(1:4)
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,-1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0

      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,+1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
      SpinDecorr =  Msq_T_BWENU

! spherical decay of Top
!      PSWgt = PSWgt * Msq_T_BW/(2d0*Ga_Top(0)*m_Top)
!      PSWgt = PSWgt * Msq_W_ENU/(2d0*Ga_W*m_W)
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
      BarSpi(1:4) = ExtParticle(2)%Pol(1:4)
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,-1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0

      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,+1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
      SpinDecorr = SpinDecorr * Msq_T_BWENU
ENDIF
    call InitCurrCache()
    call SetPropagators()

    do iHel=nHel(1),nHel(2)
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
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
        LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol*PDFFac*SpinDecorr
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
! PDFFac=1d0; print *, "set pdfs to one"
    ISFac = MomCrossing(MomExt)

IF( TOPDECAYS.GE.1 ) THEN
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
ELSEIF( TOPDECAYS.LE.-1 ) THEN
! spherical decay of ATop
!      Msq_T_BW = dsqrt(2d0)*GF*m_top**4*(1d0-m_w**2/m_top**2)*(1d0+2d0*m_w**2/m_top**2)  ! sq.mat.el. for T -> b W
!      Msq_W_ENU = dsqrt(2d0)*GF*m_W**4/3d0*4d0  ! sq.mat.el. for W -> l nu
!      PSWgt = PSWgt * Msq_T_BW/(2d0*Ga_Top(0)*m_Top)
!      PSWgt = PSWgt * Msq_W_ENU/(2d0*Ga_W*m_W)
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      Spi(1:4) = ExtParticle(1)%Pol(1:4)
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,-1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0

      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,+1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
      SpinDecorr =  Msq_T_BWENU

! spherical decay of Top
!      PSWgt = PSWgt * Msq_T_BW/(2d0*Ga_Top(0)*m_Top)
!      PSWgt = PSWgt * Msq_W_ENU/(2d0*Ga_W*m_W)
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
      BarSpi(1:4) = ExtParticle(2)%Pol(1:4)
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,-1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0

      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,+1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
      SpinDecorr = SpinDecorr * Msq_T_BWENU
ENDIF
    call InitCurrCache()
    call SetPropagators()
    do iHel=nHel(1),nHel(2)
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
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
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol*PDFFac*SpinDecorr

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
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + ParityFlip*Col1L_ttbqqb(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(BosonicPartAmp(jPrimAmp,-2:1)) )
      enddo
      enddo
      NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1) + NLO_Res_Pol(-2:1)*PDFFac*SpinDecorr

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

      FermionPartAmp(1,-2:1) = ( Nf_light*PrimAmps(PrimAmp2_1234)%Result(-2:1) + PrimAmps(PrimAmp2m_1234)%Result(-2:1) )
      FermionPartAmp(2,-2:1) = -1d0/Nc * ( Nf_light*PrimAmps(PrimAmp2_1234)%Result(-2:1) + PrimAmps(PrimAmp2m_1234)%Result(-2:1) )

      NLO_Res_Pol(-2:1) = (0d0,0d0)
      do jPrimAmp=1,2
      do iPrimAmp=1,NumBornAmps
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + ParityFlip*Col1L_ttbqqb(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(FermionPartAmp(jPrimAmp,-2:1)) )
      enddo
      enddo
      NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) + NLO_Res_Pol(-2:1)*PDFFac*SpinDecorr
   enddo!helicity loop
  enddo! npdf loop
  call swapMom(MomExt(1:4,1),MomExt(1:4,2))   ! swap back to original order, for ID below
! print *, "npdf swapping off"
ENDIF




IF( CORRECTION.EQ.0 ) THEN
!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * WidthExpansion
   EvalCS_1L_ttbqqb = LO_Res_Unpol * PreFac

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

!    MCFM comparison
!    print *, "res(-2): ",dble(NLO_Res_UnPol(-2))
!    print *, "res(-1): ",dble(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))
!    print *, "res(0)+res(1): ", dble(NLO_Res_UnPol(0))+dble(NLO_Res_UnPol(1))+dble(NLO_Res_UnPol_Ferm(0))+dble(NLO_Res_UnPol_Ferm(1))
!    print *, "res(0+1)_ferm: ", dble( NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1))

! print *, "virt 1/eps2",NLO_Res_UnPol(-2)/(alpha_sOver2Pi*RunFactor)/LO_Res_Unpol
! print *, "virt 1/eps",(dble(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1)))/(alpha_sOver2Pi*RunFactor)/LO_Res_Unpol


   EvalCS_1L_ttbqqb = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac

ELSEIF( CORRECTION.EQ.3 ) THEN
   MomP(1:4,1) = MomExt(1:4,3)
   MomP(1:4,2) = MomExt(1:4,4)
   MomP(1:4,3) =-MomExt(1:4,1)
   MomP(1:4,4) =-MomExt(1:4,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   ISFac = MomCrossing(MomExt)

IF( PROCESS.EQ.6 ) THEN
   IF( TOPDECAYS.GE.1 ) THEN
      xE = yRnd(13)
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
      call EvalIntDipoles_QQBTTBG(MomP(1:4,1:4),MomDK(1:4,1:6),xE,HOp(1:3))
   ELSEIF( TopDecays.EQ.0 ) THEN
      xE = yRnd(5)
! xE=0.46d0; print *, "fixed xE"
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
      call EvalIntDipoles_QQBTTBG_noDK(MomP(1:4,1:4),xE,HOp(1:3))
! print *, "IntDip",HOp(1)/(alpha_sOver2Pi*RunFactor)/LO_Res_Unpol
! pause
   ELSE
      xE = yRnd(5)
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
      call EvalIntDipoles_QQBTTBG_noDK(MomP(1:4,1:4),xE,HOp(1:3))
! spherical decay of ATop
!      Msq_T_BW = dsqrt(2d0)*GF*m_top**4*(1d0-m_w**2/m_top**2)*(1d0+2d0*m_w**2/m_top**2)  ! sq.mat.el. for T -> b W
!      Msq_W_ENU = dsqrt(2d0)*GF*m_W**4/3d0*4d0  ! sq.mat.el. for W -> l nu
!      PSWgt = PSWgt * Msq_T_BW/(2d0*Ga_Top(0)*m_Top)
!      PSWgt = PSWgt * Msq_W_ENU/(2d0*Ga_W*m_W)
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      Spi(1:4) = ExtParticle(1)%Pol(1:4)
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,-1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0

      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,+1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
      SpinDecorr =  Msq_T_BWENU

! spherical decay of Top
!      PSWgt = PSWgt * Msq_T_BW/(2d0*Ga_Top(0)*m_Top)
!      PSWgt = PSWgt * Msq_W_ENU/(2d0*Ga_W*m_W)
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
      BarSpi(1:4) = ExtParticle(2)%Pol(1:4)
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,-1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0

      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,+1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
      SpinDecorr = SpinDecorr * Msq_T_BWENU
   ENDIF
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac * SpinDecorr
   EvalCS_1L_ttbqqb = HOp(1)    * (pdf(Up_,1)*pdf(AUp_,2)+pdf(Dn_,1)*pdf(ADn_,2)+pdf(Chm_,1)*pdf(AChm_,2)+pdf(Str_,1)*pdf(AStr_,2)+pdf(Bot_,1)*pdf(ABot_,2) ) &
                    + HOp(2)/xE * (pdf_z(Up_,1)*pdf(AUp_,2)+pdf_z(Dn_,1)*pdf(ADn_,2)+pdf_z(Chm_,1)*pdf(AChm_,2)+pdf_z(Str_,1)*pdf(AStr_,2)+pdf_z(Bot_,1)*pdf(ABot_,2) ) &
                    + HOp(3)/xE * (pdf(Up_,1)*pdf_z(AUp_,2)+pdf(Dn_,1)*pdf_z(ADn_,2)+pdf(Chm_,1)*pdf_z(AChm_,2)+pdf(Str_,1)*pdf_z(AStr_,2)+pdf(Bot_,1)*pdf_z(ABot_,2) )

    call swapMom(MomP(1:4,3),MomP(1:4,4))
    IF( TOPDECAYS.GE.1 ) THEN
      call EvalIntDipoles_QQBTTBG(MomP(1:4,1:4),MomDK(1:4,1:6),xE,HOp(1:3))
    ELSE
      call EvalIntDipoles_QQBTTBG_noDK(MomP(1:4,1:4),xE,HOp(1:3))
    ENDIF
    HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac * SpinDecorr
    EvalCS_1L_ttbqqb = EvalCS_1L_ttbqqb  &
                     + HOp(1)    * (pdf(Up_,2)*pdf(AUp_,1)+pdf(Dn_,2)*pdf(ADn_,1)+pdf(Chm_,2)*pdf(AChm_,1)+pdf(Str_,2)*pdf(AStr_,1)+pdf(Bot_,2)*pdf(ABot_,1) ) &
                     + HOp(2)/xE * (pdf_z(Up_,2)*pdf(AUp_,1)+pdf_z(Dn_,2)*pdf(ADn_,1)+pdf_z(Chm_,2)*pdf(AChm_,1)+pdf_z(Str_,2)*pdf(AStr_,1)+pdf_z(Bot_,2)*pdf(ABot_,1) ) &
                     + HOp(3)/xE * (pdf(Up_,2)*pdf_z(AUp_,1)+pdf(Dn_,2)*pdf_z(ADn_,1)+pdf(Chm_,2)*pdf_z(AChm_,1)+pdf(Str_,2)*pdf_z(AStr_,1)+pdf(Bot_,2)*pdf_z(ABot_,1) )


ELSEIF( PROCESS.EQ.3 ) THEN
   IF( TopDecays.GE.1 ) THEN
      xE = yRnd(13)
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
      call EvalIntDipoles_QGTTBQ(MomP(1:4,1:4),MomDK(1:4,1:6),xE,HOp(1:3))
   ELSEIF( TopDecays.EQ.0 ) THEN
      xE = yRnd(5)
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
      call EvalIntDipoles_QGTTBQ_noDK(MomP(1:4,1:4),xE,HOp(1:3))
   ELSE
      xE = yRnd(5)
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
      call EvalIntDipoles_QGTTBQ_noDK(MomP(1:4,1:4),xE,HOp(1:3))

! spherical decay of ATop
!      Msq_T_BW = dsqrt(2d0)*GF*m_top**4*(1d0-m_w**2/m_top**2)*(1d0+2d0*m_w**2/m_top**2)  ! sq.mat.el. for T -> b W
!      Msq_W_ENU = dsqrt(2d0)*GF*m_W**4/3d0*4d0  ! sq.mat.el. for W -> l nu
!      PSWgt = PSWgt * Msq_T_BW/(2d0*Ga_Top(0)*m_Top)
!      PSWgt = PSWgt * Msq_W_ENU/(2d0*Ga_W*m_W)
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      Spi(1:4) = ExtParticle(1)%Pol(1:4)
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,-1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0

      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,+1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
      SpinDecorr =  Msq_T_BWENU

! spherical decay of Top
!      PSWgt = PSWgt * Msq_T_BW/(2d0*Ga_Top(0)*m_Top)
!      PSWgt = PSWgt * Msq_W_ENU/(2d0*Ga_W*m_W)
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
      BarSpi(1:4) = ExtParticle(2)%Pol(1:4)
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,-1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0

      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,+1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
      SpinDecorr = SpinDecorr * Msq_T_BWENU
   ENDIF
    HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac * SpinDecorr
    EvalCS_1L_ttbqqb = HOp(1)    * (pdf(Up_,1)*pdf(0,2)+pdf(Dn_,1)*pdf(0,2)+pdf(Chm_,1)*pdf(0,2)+pdf(Str_,1)*pdf(0,2)+pdf(Bot_,1)*pdf(0,2) ) &
                     + HOp(2)/xE * (pdf_z(Up_,1)*pdf(0,2)+pdf_z(Dn_,1)*pdf(0,2)+pdf_z(Chm_,1)*pdf(0,2)+pdf_z(Str_,1)*pdf(0,2)+pdf_z(Bot_,1)*pdf(0,2) ) &
                     + HOp(3)/xE * (pdf(Up_,1)*pdf_z(0,2)+pdf(Dn_,1)*pdf_z(0,2)+pdf(Chm_,1)*pdf_z(0,2)+pdf(Str_,1)*pdf_z(0,2)+pdf(Bot_,1)*pdf_z(0,2) )


    call swapMom(MomP(1:4,3),MomP(1:4,4))
    IF( TopDecays.GE.1 ) THEN
      call EvalIntDipoles_QGTTBQ(MomP(1:4,1:4),MomDK(1:4,1:6),xE,HOp(1:3))
    ELSE
      call EvalIntDipoles_QGTTBQ_noDK(MomP(1:4,1:4),xE,HOp(1:3))
    ENDIF
    HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac * SpinDecorr
    EvalCS_1L_ttbqqb = EvalCS_1L_ttbqqb  &
                    + HOp(1)    * (pdf(Up_,2)*pdf(0,1)+pdf(Dn_,2)*pdf(0,1)+pdf(Chm_,2)*pdf(0,1)+pdf(Str_,2)*pdf(0,1)+pdf(Bot_,2)*pdf(0,1) ) &
                    + HOp(2)/xE * (pdf_z(Up_,2)*pdf(0,1)+pdf_z(Dn_,2)*pdf(0,1)+pdf_z(Chm_,2)*pdf(0,1)+pdf_z(Str_,2)*pdf(0,1)+pdf_z(Bot_,2)*pdf(0,1) ) &
                    + HOp(3)/xE * (pdf(Up_,2)*pdf_z(0,1)+pdf(Dn_,2)*pdf_z(0,1)+pdf(Chm_,2)*pdf_z(0,1)+pdf(Str_,2)*pdf_z(0,1)+pdf(Bot_,2)*pdf_z(0,1) )


ELSEIF( PROCESS.EQ.4 ) THEN
   IF( TopDecays.GE.1 ) THEN
      xE = yRnd(13)
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
      call EvalIntDipoles_QGTTBQ(MomP(1:4,1:4),MomDK(1:4,1:6),xE,HOp(1:3))
   ELSEIF( TopDecays.EQ.0 ) THEN
      xE = yRnd(5)
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
      call EvalIntDipoles_QGTTBQ_noDK(MomP(1:4,1:4),xE,HOp(1:3))
   ELSE
      xE = yRnd(5)
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
      call EvalIntDipoles_QGTTBQ_noDK(MomP(1:4,1:4),xE,HOp(1:3))
! spherical decay of ATop
!      Msq_T_BW = dsqrt(2d0)*GF*m_top**4*(1d0-m_w**2/m_top**2)*(1d0+2d0*m_w**2/m_top**2)  ! sq.mat.el. for T -> b W
!      Msq_W_ENU = dsqrt(2d0)*GF*m_W**4/3d0*4d0  ! sq.mat.el. for W -> l nu
!      PSWgt = PSWgt * Msq_T_BW/(2d0*Ga_Top(0)*m_Top)
!      PSWgt = PSWgt * Msq_W_ENU/(2d0*Ga_W*m_W)
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      Spi(1:4) = ExtParticle(1)%Pol(1:4)
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,-1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0

      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,+1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
      SpinDecorr =  Msq_T_BWENU

! spherical decay of Top
!      PSWgt = PSWgt * Msq_T_BW/(2d0*Ga_Top(0)*m_Top)
!      PSWgt = PSWgt * Msq_W_ENU/(2d0*Ga_W*m_W)
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
      BarSpi(1:4) = ExtParticle(2)%Pol(1:4)
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,-1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0

      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,+1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
      SpinDecorr = SpinDecorr * Msq_T_BWENU
   ENDIF
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac * SpinDecorr
   EvalCS_1L_ttbqqb =HOp(1)    * (pdf(AUp_,1)*pdf(0,2)+pdf(ADn_,1)*pdf(0,2)+pdf(AChm_,1)*pdf(0,2)+pdf(AStr_,1)*pdf(0,2)+pdf(ABot_,1)*pdf(0,2) ) &
                    +HOp(2)/xE * (pdf_z(AUp_,1)*pdf(0,2)+pdf_z(ADn_,1)*pdf(0,2)+pdf_z(AChm_,1)*pdf(0,2)+pdf_z(AStr_,1)*pdf(0,2)+pdf_z(ABot_,1)*pdf(0,2) ) &
                    +HOp(3)/xE * (pdf(AUp_,1)*pdf_z(0,2)+pdf(ADn_,1)*pdf_z(0,2)+pdf(AChm_,1)*pdf_z(0,2)+pdf(AStr_,1)*pdf_z(0,2)+pdf(ABot_,1)*pdf_z(0,2) )


    call swapMom(MomP(1:4,3),MomP(1:4,4))
    IF( TopDecays.GE.1 ) THEN
      call EvalIntDipoles_QGTTBQ(MomP(1:4,1:4),MomDK(1:4,1:6),xE,HOp(1:3))
    ELSE
      call EvalIntDipoles_QGTTBQ_noDK(MomP(1:4,1:4),xE,HOp(1:3))
   ENDIF
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac * SpinDecorr
   EvalCS_1L_ttbqqb = EvalCS_1L_ttbqqb &
                    + HOp(1)    * (pdf(AUp_,2)*pdf(0,1)+pdf(ADn_,2)*pdf(0,1)+pdf(AChm_,2)*pdf(0,1)+pdf(AStr_,2)*pdf(0,1)+pdf(ABot_,2)*pdf(0,1) ) &
                    + HOp(2)/xE * (pdf_z(AUp_,2)*pdf(0,1)+pdf_z(ADn_,2)*pdf(0,1)+pdf_z(AChm_,2)*pdf(0,1)+pdf_z(AStr_,2)*pdf(0,1)+pdf_z(ABot_,2)*pdf(0,1) ) &
                    + HOp(3)/xE * (pdf(AUp_,2)*pdf_z(0,1)+pdf(ADn_,2)*pdf_z(0,1)+pdf(AChm_,2)*pdf_z(0,1)+pdf(AStr_,2)*pdf_z(0,1)+pdf(ABot_,2)*pdf_z(0,1) )
ENDIF
ENDIF

   if( IsNan(EvalCS_1L_ttbqqb) ) then
        print *, "NAN:",EvalCS_1L_ttbqqb
        print *, yRnd(:)
        print *, NLO_Res_UnPol(0),NLO_Res_UnPol(1),NLO_Res_UnPol_Ferm(0),NLO_Res_UnPol_Ferm(1),IOp(0)
        print *, PSWgt , VgsWgt , PDFFac_a,PDFFac_b, sHatJacobi
        print *, eta1,eta2,((1d0-eta1)*xE+eta1),((1d0-eta2)*xE+eta2),MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,1)
        print *, MomExt(1:4,2)
        print *, MomExt(1:4,3)
        print *, MomExt(1:4,4)
        print *, MomDK(1:4,1)
        print *, MomDK(1:4,2)
        print *, MomDK(1:4,3)
        print *, MomDK(1:4,4)
        print *, MomDK(1:4,5)
        print *, MomDK(1:4,6)
        print *, "SKIP EVENT!!!!!"
        EvalCS_1L_ttbqqb = 0d0
        return
   endif


IF( ObsSet.EQ.7 ) THEN
   do NHisto=1,3
      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_ttbqqb)
   enddo
   call intoHisto(4,1,MInv_LB*EvalCS_1L_ttbqqb)
   call intoHisto(5,1,MInv_LB**2*EvalCS_1L_ttbqqb)
   call intoHisto(6,1,xFrag*EvalCS_1L_ttbqqb)
ELSE
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_ttbqqb)
   enddo
ENDIF


   EvalCS_1L_ttbqqb = EvalCS_1L_ttbqqb/VgsWgt
   EvalCounter = EvalCounter + 1


! !      MADGRAPH CHECK: uub->ttb
!        MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!        MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!        MG_MOM(0:3,3) = MomExt(1:4,3)*100d0
!        MG_MOM(0:3,4) = MomExt(1:4,4)*100d0
!        call coupsm(0)
!        call SUUB_TTB(MG_MOM,MadGraph_tree)
!        MadGraph_tree = MadGraph_tree *100d0**2
!        print *, ""
!        print *, "My tree:         ", LO_Res_Unpol
!        print *, "MadGraph hel.amp:", MadGraph_tree
!        print *, "MG/ME ratio: ", MadGraph_tree/dble(LO_Res_Unpol)
!        pause
return
END FUNCTION





FUNCTION EvalCS_Real_ttbggg_MPI(yRnd,VgsWgt,res)
implicit none
integer :: EvalCS_Real_ttbggg_MPI
real(8) ::  yRnd(*),res(*),VgsWgt

res(1) = EvalCS_Real_ttbggg(yRnd,VgsWgt)
EvalCS_Real_ttbggg_MPI=0
RETURN
END FUNCTION




FUNCTION EvalCS_Real_ttbggg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModDipoles_GGTTBG
use ModDipoles_GGTTBG_noDK
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModJPsiFrag
implicit none
real(8) ::  EvalCS_Real_ttbggg,yRnd(1:VegasMxDim),VgsWgt
complex(8) :: LO_Res_Pol,LO_Res_Unpol,Spi(1:4),BarSpi(1:4),M_T_BWENU,Msq_T_BWENU
integer :: iHel,iPrimAmp,jPrimAmp
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,DipoleResult,ISFac
real(8) :: MomExt(1:4,1:5),MomDK(1:4,1:6)
real(8) :: MG_MOM(0:3,1:5),xFrag,r_sc
real(8) :: MadGraph_tree,GG_TTBG,GG_TTB,GG_UUB,Msq_T_BW,Msq_W_ENU
logical :: applyPSCut,applySingCut
real(8) :: tau,eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac
real(8) :: pdf(-6:6,1:2)
integer :: NBin(1:NumMaxHisto),NHisto,ParityFlip
! real(8) :: MomBoost(1:4),MomLep1(1:4),MomLep2(1:4),MomPart(1:4,7),MomJet(1:4,1:7)
real(8) :: s13,s23,s12,s55
real(8),parameter :: ThresholdCutOff = 1.00d0
include "vegas_common.f"


! IF( TopDecays.ge.1 ) THEN
  ParityFlip=1
! ELSE
!   ParityFlip=2
! ENDIF

!  for gensing
!    print *, "fixed yRnd"
!   yRnd(1) =0.669596162750198
!   yRnd(2) =0.221727079861671
!   yRnd(3) =6.927929705440973d-002
!   yRnd(4) =8.997565581927815d-002
!   yRnd(5) =6.236811171396189d-004
!   yRnd(6) =0.393968204685472
!   yRnd(7) =0.116716410553417
!   yRnd(8) =0.441485913675971
!   yRnd(9) =0.489933605068333
!   yRnd(10) =0.216611683236720
!   yRnd(11) =0.316609297327981
!   yRnd(12) =7.358870961404818d-002
!   yRnd(13) =2.194301389248251d-002
!   yRnd(14) =0.133462533416908
!   yRnd(15) =0.284648654882167
!   yRnd(16) =0.6152712171269913

  EvalCS_Real_ttbggg = 0d0

  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top * ThresholdCutOff ) then
      EvalCS_Real_ttbggg = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to3(EHat,yRnd(3:7),MomExt(1:4,1:5),PSWgt)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))
   ISFac = MomCrossing(MomExt)

IF( TOPDECAYS.NE.0 ) THEN
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomDK(1:4,1:3),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:15),.false.,MomDK(1:4,4:6),PSWgt3)
   PSWgt = PSWgt * PSWgt2*PSWgt3

   IF( TOPDECAYS.GE.1 ) THEN
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
   ELSE
! spherical decay of ATop
!       Msq_T_BW = dsqrt(2d0)*GF*m_top**4*(1d0-m_w**2/m_top**2)*(1d0+2d0*m_w**2/m_top**2)  ! sq.mat.el. for T -> b W
!       Msq_W_ENU = dsqrt(2d0)*GF*m_W**4/3d0*4d0  ! sq.mat.el. for W -> l nu
!       PSWgt = PSWgt * Msq_T_BW/(2d0*Ga_Top(0)*m_Top)
!       PSWgt = PSWgt * Msq_W_ENU/(2d0*Ga_W*m_W)
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      Spi(1:4) = ExtParticle(1)%Pol(1:4)
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,-1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0

      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,+1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
      PSWgt = PSWgt * Msq_T_BWENU

! spherical decay of Top
!       PSWgt = PSWgt * Msq_T_BW/(2d0*Ga_Top(0)*m_Top)
!       PSWgt = PSWgt * Msq_W_ENU/(2d0*Ga_W*m_W)
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
      BarSpi(1:4) = ExtParticle(2)%Pol(1:4)
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,-1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0

      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,+1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
      PSWgt = PSWgt * Msq_T_BWENU
   ENDIF

   IF( TOPDECAYS.EQ.5 .OR. TOPDECAYS.EQ.6) THEN
          xFrag= yRnd(16)
          PSWgt = PSWgt * FF(xFrag)
   ENDIF

ENDIF

   call CheckSing(MomExt,applySingCut)
   if( applySingCut ) then
       EvalCS_Real_ttbggg = 0d0
       return
   endif
   call Kinematics_TTBAR(.true.,MomExt,MomDK,applyPSCut,NBin,xJPsiFrag=xFrag)

   call InitCurrCache()
   call setPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   RunFactor = RunAlphaS(2,MuRen)
   if( applyPSCut ) then
        EvalCS_Real_ttbggg = 0d0
   else
        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
          call HelCrossing(Helicities(iHel,1:NumExtParticles))
          call SetPolarizations()
          do iPrimAmp=1,NumPrimAmps
              call EvalTree(BornAmps(iPrimAmp))
          enddo
          LO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,NumBornAmps
          do iPrimAmp=1,NumBornAmps
              LO_Res_Pol = LO_Res_Pol + ColLO_ttbggg(iPrimAmp,jPrimAmp) * ParityFlip*BornAmps(iPrimAmp)%Result * dconjg(BornAmps(jPrimAmp)%Result)
          enddo
          enddo
          LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
      enddo!helicity loop

      LO_Res_UnPol = LO_Res_UnPol * (alpha_s4Pi*RunFactor)**3 * ISFac
      EvalCS_Real_ttbggg = LO_Res_UnPol * PreFac

IF( ObsSet.EQ.7 ) THEN
   do NHisto=1,3
      call intoHisto(NHisto,NBin(NHisto),EvalCS_Real_ttbggg)
   enddo
   call intoHisto(4,1,MInv_LB*EvalCS_Real_ttbggg)
   call intoHisto(5,1,MInv_LB**2*EvalCS_Real_ttbggg)
   call intoHisto(6,1,xFrag*EvalCS_Real_ttbggg)

ELSE
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_Real_ttbggg)
   enddo
ENDIF
endif!applyPSCut

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

IF( TopDecays.GE.1 ) THEN
     call EvalDipoles_GGTTBG((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:16),PreFac,DipoleResult)
ELSE
     call EvalDipoles_GGTTBG_noDK((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:15),PreFac,DipoleResult)
ENDIF

!   s12 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,2))
!   s13 = 2d0*(MomExt(1:4,2).dot.MomExt(1:4,3))
!   write(* ,"(1PE23.16,3X,1PE23.16,3X,1PE23.16,3X,1PE23.16)") s13/s12,EvalCS_Real_ttbggg,DipoleResult, (EvalCS_Real_ttbggg/(-DipoleResult)-1d0)
!   pause


   if( IsNan(EvalCS_Real_ttbggg) .or. IsNan(DipoleResult) ) then
        print *, "NAN:",EvalCS_Real_ttbggg,DipoleResult
        print *, yRnd(:)
        print *, PSWgt , VgsWgt , PDFFac, sHatJacobi
        print *, eta1,eta2,MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,1)
        print *, MomExt(1:4,2)
        print *, MomExt(1:4,3)
        print *, MomExt(1:4,4)
        print *, MomExt(1:4,5)
        print *, MomDK(1:4,1)
        print *, MomDK(1:4,2)
        print *, MomDK(1:4,3)
        print *, MomDK(1:4,4)
        print *, MomDK(1:4,5)
        print *, MomDK(1:4,6)
        stop
   endif

        EvalCS_Real_ttbggg = (EvalCS_Real_ttbggg + DipoleResult)/VgsWgt
return
END FUNCTION






FUNCTION EvalCS_Real_ttbqqbg_MPI(yRnd,VgsWgt,res)
implicit none
integer :: EvalCS_Real_ttbqqbg_MPI
real(8) ::  yRnd(*),res(*),VgsWgt

res(1) = EvalCS_Real_ttbqqbg(yRnd,VgsWgt)
EvalCS_Real_ttbqqbg_MPI=0
RETURN
END FUNCTION






FUNCTION EvalCS_Real_ttbqqbg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModDipoles_QQBTTBG
use ModDipoles_QQBTTBG_noDK
use ModDipoles_QGTTBQ
use ModDipoles_QGTTBQ_noDK
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModJPsiFrag
implicit none
real(8) ::  EvalCS_Real_ttbqqbg,yRnd(1:VegasMxDim),VgsWgt,DipoleResult,EvalCS_Dips_ttbqqbg
complex(8) :: LO_Res_Pol,LO_Res_Unpol,Spi(1:4),BarSpi(1:4),M_T_BWENU,Msq_T_BWENU
integer :: iHel,iPrimAmp,jPrimAmp,NPDF
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,s12,s13,xFrag
real(8) :: PDFFac,PDFFac_a,PDFFac_b,pdf(-6:6,1:2),tau,eta1,eta2,FluxFac,sHatJacobi,ISFac,PreFac
real(8) :: MomExt(1:4,1:5),MomDK(1:4,1:6),AccPoles,r_sc
real(8) :: MG_MOM(0:3,1:5)
integer :: NBin(1:NumMaxHisto),NHisto,ParityFlip
real(8) :: MadGraph_tree,SpinDecorr
logical :: applyPSCut,applySingCut
real(8),parameter :: ThresholdCutOff = 1.00d0
include "vegas_common.f"


  ParityFlip=1
  EvalCS_Real_ttbqqbg = 0d0
  EvalCS_Dips_ttbqqbg = 0d0

  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top * ThresholdCutOff ) then
      EvalCS_Real_ttbqqbg = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to3(EHat,yRnd(3:7),MomExt(1:4,1:5),PSWgt)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))

   call CheckSing(MomExt,applySingCut)
   if( applySingCut ) then
      EvalCS_Real_ttbqqbg = 0d0
      return
   endif
   SpinDecorr = 1d0
IF( TOPDECAYS.NE.0 ) THEN
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomDK(1:4,1:3),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:15),.false.,MomDK(1:4,4:6),PSWgt3)
   PSWgt = PSWgt * PSWgt2*PSWgt3
   IF( TOPDECAYS.EQ.5 .OR. TOPDECAYS.EQ.6) THEN
       xFrag= yRnd(16)
       PSWgt = PSWgt * FF(xFrag)
   ENDIF
ENDIF

   call Kinematics_TTBAR(.true.,MomExt,MomDK,applyPSCut,NBin,xJPsiFrag=xFrag)

   call setPDFs(eta1,eta2,MuFac,pdf)
IF( PROCESS.EQ.6 ) THEN
   PDFFac_a = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
            + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
            + pdf(Bot_,1)*pdf(ABot_,2)
   PDFFac_b = pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
            + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
            + pdf(Bot_,2)*pdf(ABot_,1)
ELSEIF( PROCESS.EQ.3 ) THEN
   PDFFac_a = pdf(Up_,1) *pdf(0,2) + pdf(Dn_,1) *pdf(0,2)  &
            + pdf(Chm_,1)*pdf(0,2) + pdf(Str_,1)*pdf(0,2)  &
            + pdf(Bot_,1)*pdf(0,2)
   PDFFac_b = pdf(Up_,2) *pdf(0,1) + pdf(Dn_,2) *pdf(0,1)  &
            + pdf(Chm_,2)*pdf(0,1) + pdf(Str_,2)*pdf(0,1)  &
            + pdf(Bot_,2)*pdf(0,1)
ELSEIF( PROCESS.EQ.4 ) THEN
   PDFFac_a = pdf(AUp_,2) *pdf(0,1) + pdf(ADn_,2) *pdf(0,1)  &
            + pdf(AChm_,2)*pdf(0,1) + pdf(AStr_,2)*pdf(0,1)  &
            + pdf(ABot_,2)*pdf(0,1)
   PDFFac_b = pdf(AUp_,1) *pdf(0,2) + pdf(ADn_,1) *pdf(0,2)  &
            + pdf(AChm_,1)*pdf(0,2) + pdf(AStr_,1)*pdf(0,2)  &
            + pdf(ABot_,1)*pdf(0,2)
ENDIF
   RunFactor = RunAlphaS(NLOParam,MuRen)

DO NPDF=1,2
   if(npdf.eq.1) then
        PDFFac = PDFFac_a
   elseif(npdf.eq.2) then
        PDFFac = PDFFac_b
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
   endif
   ISFac = MomCrossing(MomExt)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt

   IF( TOPDECAYS.GE.1 ) THEN
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
   ELSEIF( TOPDECAYS.LE.-1 ) THEN
! spherical decay of ATop
!       Msq_T_BW = dsqrt(2d0)*GF*m_top**4*(1d0-m_w**2/m_top**2)*(1d0+2d0*m_w**2/m_top**2)  ! sq.mat.el. for T -> b W
!       Msq_W_ENU = dsqrt(2d0)*GF*m_W**4/3d0*4d0  ! sq.mat.el. for W -> l nu
!       PSWgt = PSWgt * Msq_T_BW/(2d0*Ga_Top(0)*m_Top)
!       PSWgt = PSWgt * Msq_W_ENU/(2d0*Ga_W*m_W)
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      Spi(1:4) = ExtParticle(1)%Pol(1:4)
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,-1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0

      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,+1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
      SpinDecorr =  Msq_T_BWENU

! spherical decay of Top
!       PSWgt = PSWgt * Msq_T_BW/(2d0*Ga_Top(0)*m_Top)
!       PSWgt = PSWgt * Msq_W_ENU/(2d0*Ga_W*m_W)
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
      BarSpi(1:4) = ExtParticle(2)%Pol(1:4)
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,-1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0

      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,+1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
      SpinDecorr =  SpinDecorr * Msq_T_BWENU
   ENDIF

   if( applyPSCut ) then
        EvalCS_Real_ttbqqbg = 0d0
   else
        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities
          call HelCrossing(Helicities(iHel,1:NumExtParticles))
          call SetPolarizations()
          do iPrimAmp=1,NumPrimAmps
              call EvalTree(BornAmps(iPrimAmp))
          enddo
          LO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,NumBornAmps
          do iPrimAmp=1,NumBornAmps
              LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqbg(iPrimAmp,jPrimAmp) * ParityFlip*BornAmps(iPrimAmp)%Result * dconjg(BornAmps(jPrimAmp)%Result)
          enddo
          enddo
          LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
      enddo!helicity loop
      LO_Res_UnPol = LO_Res_UnPol * (alpha_s4Pi*RunFactor)**3 * ISFac * PreFac *PDFFac * SpinDecorr


      IF( ObsSet.EQ.7 ) THEN
        do NHisto=1,3
            call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol))
        enddo
        call intoHisto(4,1,MInv_LB*dble(LO_Res_Unpol))
        call intoHisto(5,1,MInv_LB**2*dble(LO_Res_Unpol))
        call intoHisto(6,1,xFrag*dble(LO_Res_Unpol))
      ELSE
        do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol))
        enddo
      ENDIF
      EvalCS_Real_ttbqqbg = EvalCS_Real_ttbqqbg + dble(LO_Res_Unpol)
   endif!applyPSCut


!       MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!       MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!       MG_MOM(0:3,3) = MomExt(1:4,5)*100d0
!       MG_MOM(0:3,4) = MomExt(1:4,4)*100d0
!       MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!       call coupsm(0)
!       call SUUB_TTBG(MG_MOM,MadGraph_tree)
!       call SUG_TTBU(MG_MOM,MadGraph_tree)
!       call SUBG_TTBUB(MG_MOM,MadGraph_tree)
!       print *, ""
!       print *, "Unpolarized LO result:"
!       print *, "My tree:          ", LO_Res_Unpol
!       print *, "MadGraph hel.amp:", MadGraph_tree*(100d0)**2
!       print *, "ratio x: ", MadGraph_tree*(100d0)**2/dble((LO_Res_Unpol))
!       pause


PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac * SpinDecorr
IF( PROCESS.EQ.6 ) THEN
  IF( TOPDECAYS.GE.1 ) THEN
    call EvalDipoles_QQBTTBG((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:16),PreFac,DipoleResult)
  ELSE
    call EvalDipoles_QQBTTBG_noDK((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:15),PreFac,DipoleResult)
  ENDIF
ELSEIF( PROCESS.EQ.3 ) THEN
  IF( TOPDECAYS.GE.1 ) THEN
    call EvalDipoles_QGTTBQ((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:16),PreFac,DipoleResult)
  ELSE
    call EvalDipoles_QGTTBQ_noDK((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:15),PreFac,DipoleResult)
  ENDIF
ELSEIF( PROCESS.EQ.4 ) THEN
  IF( TOPDECAYS.GE.1 ) THEN
    call EvalDipoles_QGTTBQ((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:16),PreFac,DipoleResult)
  ELSE
    call EvalDipoles_QGTTBQ_noDK((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:15),PreFac,DipoleResult)
  ENDIF
ENDIF
    EvalCS_Dips_ttbqqbg = EvalCS_Dips_ttbqqbg + DipoleResult
ENDDO! loop over a<-->b pdfs


!     s12 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,2))
!     s13 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,3))
!     write(* ,"(1PE23.16,3X,1PE23.16,3X,1PE23.16,3X,1PE23.16)") s13/s12,EvalCS_Real_ttbqqbg,EvalCS_Dips_ttbqqbg, (EvalCS_Real_ttbqqbg/(-EvalCS_Dips_ttbqqbg)-1d0)
!     pause

    EvalCS_Real_ttbqqbg = (EvalCS_Real_ttbqqbg + EvalCS_Dips_ttbqqbg) /VgsWgt

   if( IsNan(EvalCS_Real_ttbqqbg) .or. IsNan(EvalCS_Dips_ttbqqbg) ) then
        print *, "NAN:",EvalCS_Real_ttbqqbg,EvalCS_Dips_ttbqqbg
        print *, yRnd(:)
        print *, PSWgt , VgsWgt , PDFFac, sHatJacobi
        print *, eta1,eta2,MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,1)
        print *, MomExt(1:4,2)
        print *, MomExt(1:4,3)
        print *, MomExt(1:4,4)
        print *, MomExt(1:4,5)
        print *, MomDK(1:4,1)
        print *, MomDK(1:4,2)
        print *, MomDK(1:4,3)
        print *, MomDK(1:4,4)
        print *, MomDK(1:4,5)
        print *, MomDK(1:4,6)
        stop
   endif


return
END FUNCTION







FUNCTION EvalCS_NLODK_ttb_noSC(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes
use ModMyRecurrence
use ModParameters
implicit none
real(8) ::  EvalCS_NLODK_ttb_noSC,yRnd(1:VegasMxDim),VgsWgt,xE,IOp(-2:0),HOp(1:3)
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1)
complex(8) :: BosonicPartAmp(1:2,-2:1),FermionPartAmp(1:2,-2:1)
complex(8) :: Spi(1:4),BarSpi(1:4),M_T_BWENU,Msq_T_BWENU,MsqLO_T_BW,MsqVI_T_BW,Msq_W_ENU,MsqRE_T_BW
integer :: iHel,iPrimAmp,jPrimAmp
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,ISFac,xFrag
real(8) :: MomExt(1:4,1:NumExtParticles),MomDK(1:4,1:7),MomDKTd(1:4,1:3),MomDKx(1:4,1:7),MomP(1:4)
logical :: applyPSCut,applySingCut
real(8) :: MG_MOM(0:3,1:NumExtParticles),ptb,ptg,pbg,pbDpg,ptDpg,ptDpb,omz,z,y,Dipole
real(8) :: MadGraph_tree,UUB_TTB,DKFac_Top,DKFac_ATop
real(8),parameter :: Nc=3d0
real(8) :: tau,eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac_a,PDFFac_b,PDFFac
real(8) :: pdf(-6:6,1:2),ColCorrLO(1:NumBornAmps,1:NumBornAmps),r_sc
integer :: NHisto,NBin(1:NumMaxHisto),ParityFlip,npdf
real(8) :: ThresholdCutOff = 1.0d0
real(8),parameter :: pisqo6 = DblPi**2/6d0
real(8) :: rsq,omrsq,wlog,rlog,Kfun,epinv,epinv2,F0,F1,IntDip
real(8),parameter :: CF=4d0/3d0
include "vegas_common.f"


! print *, "fixed yrnd"
! yrnd(1)= 0.120417488562138
! yrnd(2)=0.467848163315956
! yrnd(3)=0.314100383461500
! yrnd(4)=  4.800610991567660E-002
! yrnd(5)=  6.554446535443172E-002
! yrnd(6)=  0.105829211932528
! yrnd(7)=  0.492114467077011
! yrnd(8)=  2.928262810655061E-002
! yrnd(9)=  0.153130586796035
! yrnd(10)=  0.426867230994100
! yrnd(11)=  0.381774851764448
! yrnd(12)=  6.696560190383624E-003
! yrnd(13)=  4.668414455218438E-002
! yrnd(14)=  0.356730263846335
! yrnd(15)=  0.351952110348247
! yrnd(16)=  8.698269659047150E-002
! yrnd(17)=  0.171564949989116
! yrnd(18)=  0.357551317828499



  ParityFlip=2
  EvalCS_NLODK_ttb_noSC = 0d0
IF( Process.EQ.1 ) THEN
    ColCorrLO(1:NumBornAmps,1:NumBornAmps) = ColLO_ttbgg(1:NumBornAmps,1:NumBornAmps)
ELSEIF( Process.EQ.2 ) THEN
    ColCorrLO(1:NumBornAmps,1:NumBornAmps) = ColLO_ttbqqb(1:NumBornAmps,1:NumBornAmps)
ENDIF

  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top * ThresholdCutOff ) then
      EvalCS_NLODK_ttb_noSC = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to2(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt)
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))

   call setPDFs(eta1,eta2,MuFac,pdf)
IF( PROCESS.EQ.1 ) THEN
   PDFFac_a = pdf(0,1) * pdf(0,2)
   PDFFac_b = 1d300
ELSEIF( PROCESS.EQ.2 ) THEN
   PDFFac_a = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
            + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
            + pdf(Bot_,1)*pdf(ABot_,2)
   PDFFac_b = pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
            + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
            + pdf(Bot_,2)*pdf(ABot_,1)
ENDIF

   PreFac = fbGeV2 * FluxFac * sHatJacobi * VgsWgt
   RunFactor = RunAlphaS(NLOParam,MuRen)

   LO_Res_Unpol             = (0d0,0d0)
   NLO_Res_Unpol(-2:1)      = (0d0,0d0)
   NLO_Res_Unpol_Ferm(-2:1) = (0d0,0d0)

  do npdf=1,2
    if(npdf.eq.1) then
        PDFFac = PDFFac_a
    elseif(npdf.eq.2) then
        if( PROCESS.EQ.1 ) cycle
        PDFFac = PDFFac_b
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
    endif
    ISFac = MomCrossing(MomExt)
    call InitCurrCache()
    call SetPropagators()

    do iHel=1,NumHelicities
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
        call SetPolarizations()
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
        enddo

        LO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            LO_Res_Pol = LO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
        enddo
        enddo
        LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol*PDFFac
    enddo!helicity loop
  enddo! npdf loop
  call swapMom(MomExt(1:4,1),MomExt(1:4,2))   ! swap back to original order

!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2




IF( CORRECTION.EQ.4 ) THEN


   call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:8),.false.,MomDK(1:4,1:3),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(9:12),.false.,MomDK(1:4,4:6),PSWgt3)

   call Kinematics_TTBAR(.false.,MomExt,MomDK,applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_NLODK_ttb_noSC = 0d0
      return
   endif


if( DKRE_switch.eq.0 .or. DKRE_switch.eq.1 ) then
!----------------------------------------
! one loop correction to Anti-top decay |
!----------------------------------------
      call TopDecay(ExtParticle(1),DK_1L_T,MomDK(1:4,1:3))
      Spi(1:4) = ExtParticle(1)%Pol(1:4)
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,-1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      Spi(1:4) = ExtParticle(1)%Pol(1:4)
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,-1,BarSpi(1:4))
      Msq_T_BWENU =  dreal(dconjg(M_T_BWENU)*psp1_(Spi(1:4),BarSpi(1:4)))/(2d0*m_Top)

      call TopDecay(ExtParticle(1),DK_1L_T,MomDK(1:4,1:3))
      Spi(1:4) = ExtParticle(1)%Pol(1:4)
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,+1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      Spi(1:4) = ExtParticle(1)%Pol(1:4)
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,+1,BarSpi(1:4))
      Msq_T_BWENU =  Msq_T_BWENU + dreal(dconjg(M_T_BWENU)*psp1_(Spi(1:4),BarSpi(1:4)))/(2d0*m_Top)
      DKFac_ATop = Msq_T_BWENU/2d0


      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
      BarSpi(1:4) = ExtParticle(2)%Pol(1:4)
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,-1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)

      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,+1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)
      DKFac_Top = Msq_T_BWENU/2d0


      EvalCS_NLODK_ttb_noSC = dble(LO_Res_Unpol)*PreFac*PSWgt*PSWgt2*PSWgt3 * DKFac_ATop*DKFac_Top

      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),EvalCS_NLODK_ttb_noSC)
      enddo

endif



if( DKRE_switch.eq.0 .or. DKRE_switch.eq.2 ) then
!-----------------------------------
! one loop correction to top decay |
!-----------------------------------
!       MsqLO_T_BW = dsqrt(2d0)*GF*m_top**4*(1d0-m_w**2/m_top**2)*(1d0+2d0*m_w**2/m_top**2)  ! LO sq.mat.el. for T -> b W
!       MsqVI_T_BW = MsqLO_T_BW * ( F0+IntDip + 0.5d0*F1*(1d0-rsq)/(1d0+2d0*rsq) ) * alpha_sOver2Pi*RunAlphaS(2,MuRen) * CF   ! note: F1 comes with same sign!
!       Msq_W_ENU = dsqrt(2d0)*GF*m_W**4/3d0*4d0  ! sq.mat.el. for W -> l nu
!       DKFac_Top  = MsqVI_T_BW/(2d0*Ga_Top(0)*m_Top) * Msq_W_ENU/(2d0*Ga_W*m_W)          ! virt. corr. to Top
!       DKFac_ATop = MsqLO_T_BW/(2d0*Ga_Top(0)*m_Top) * Msq_W_ENU/(2d0*Ga_W*m_W)          ! LO corr. to ATop


      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      Spi(1:4) = ExtParticle(1)%Pol(1:4)
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,-1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)

      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,+1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)
      DKFac_ATop = Msq_T_BWENU/2d0




      call TopDecay(ExtParticle(2),DK_1L_T,MomDK(1:4,4:6))
      BarSpi(1:4) = ExtParticle(2)%Pol(1:4)
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,-1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)

      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
      BarSpi(1:4) = ExtParticle(2)%Pol(1:4)
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,-1,Spi(1:4))
      Msq_T_BWENU =  dreal(dconjg(M_T_BWENU)*psp1_(Spi(1:4),BarSpi(1:4)))/(2d0*m_Top)


      call TopDecay(ExtParticle(2),DK_1L_T,MomDK(1:4,4:6))
      BarSpi(1:4) = ExtParticle(2)%Pol(1:4)
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,+1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)

      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
      BarSpi(1:4) = ExtParticle(2)%Pol(1:4)
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,+1,Spi(1:4))
      Msq_T_BWENU =  Msq_T_BWENU + dreal(dconjg(M_T_BWENU)*psp1_(Spi(1:4),BarSpi(1:4)))/(2d0*m_Top)
      DKFac_Top = Msq_T_BWENU/2d0

      EvalCS_NLODK_ttb_noSC = EvalCS_NLODK_ttb_noSC + dble(LO_Res_Unpol)*PreFac*PSWgt*PSWgt2*PSWgt3 * DKFac_ATop*DKFac_Top

      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol)*PreFac*PSWgt*PSWgt2*PSWgt3 * DKFac_ATop*DKFac_Top)
      enddo
endif


ELSEIF( CORRECTION.EQ.5 ) THEN
if( DKRE_switch.eq.0 .or. DKRE_switch.eq.1 ) then
!----------------------------------------
! real gluon emission for Anti-top decay |
!----------------------------------------
   call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:11),.true.,MomDK(1:4,1:4),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(12:15),.false.,MomDK(1:4,5:7),PSWgt3)

   call CheckSing(MomDK(1:4,1:4),applySingCut)
   if( applySingCut) then
      EvalCS_NLODK_ttb_noSC = 0d0
      goto 25
   endif

   MomDKx(1:4,1:3) = MomDK(1:4,1:3)
   MomDKx(1:4,4:6) = MomDK(1:4,5:7)
   MomDKx(1:4,7)   = MomDK(1:4,4)
   call Kinematics_TTBAR(.true.,MomExt,MomDKx,applyPSCut,NBin)

   if( applyPSCut ) then
      goto 24
   endif


!     MomP(1:4) = MomExt(1:4,3)-MomDK(1:4,1)
!     rsq = (MomP(1:4).dot.MomP(1:4))/m_Top**2  ! this does not work
!       rsq = m_W**2/m_Top**2
!       ptb = MomExt(1:4,3).dot.MomDK(1:4,1)
!       pbg = MomDK(1:4,1).dot.MomDK(1:4,4)
!       ptg = MomExt(1:4,3).dot.MomDK(1:4,4)
!       Msq_W_ENU = dsqrt(2d0)*GF*m_W**4/3d0*4d0  ! sq.mat.el. for W -> l nu
!       MsqLO_T_BW = dsqrt(2d0)*GF*m_top**4*(1d0-m_w**2/m_top**2)*(1d0+2d0*m_w**2/m_top**2)  ! LO sq.mat.el. for T -> b W
!       MsqRE_T_BW = 2d0*(alpha_s4Pi*RunFactor)*CF*GF/dsqrt(2d0)*m_top**2 *(           &
!                  + m_top**2*(1d0-rsq)*(1d0+2d0*rsq)*(2d0*ptb/pbg/ptg-m_top**2/ptg**2)  &
!                  + 8d0*rsq + 2d0*(ptg-pbg)**2/pbg/ptg*(1d0+2d0*rsq)        )
!       DKFac_ATop = MsqRE_T_BW/(2d0*Ga_Top(0)*m_Top) * Msq_W_ENU/(2d0*Ga_W*m_W)          ! real corr. to ATop
!       DKFac_Top  = MsqLO_T_BW/(2d0*Ga_Top(0)*m_Top) * Msq_W_ENU/(2d0*Ga_W*m_W)          ! LO corr. to Top



      call TopDecay(ExtParticle(1),DK_RE_T,MomDK(1:4,1:4),GluonHel=+1)
      Spi(1:4) = ExtParticle(1)%Pol(1:4)
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,-1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,+1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)

      call TopDecay(ExtParticle(1),DK_RE_T,MomDK(1:4,1:4),GluonHel=-1)
      Spi(1:4) = ExtParticle(1)%Pol(1:4)
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,-1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,+1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)
      DKFac_ATop =  Msq_T_BWENU/2d0

      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,5:7))
      BarSpi(1:4) = ExtParticle(2)%Pol(1:4)
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,-1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)

      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,+1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)
      DKFac_Top = Msq_T_BWENU/2d0


      EvalCS_NLODK_ttb_noSC = dble(LO_Res_Unpol)*PreFac*PSWgt*PSWgt2*PSWgt3 *DKFac_ATop*DKFac_Top

      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),EvalCS_NLODK_ttb_noSC)
      enddo

!       print *, "Ma",dble(LO_Res_Unpol)*PreFac*PSWgt*PSWgt2*PSWgt3 *DKFac_ATop*DKFac_Top


! ENDIF

24 continue

!-------------------------------------
! dipole subtraction for Atop-decay  |
!-------------------------------------
      call WTransform(MomDK(1:4,1:4),MomDKTd(1:4,1:3),pbDpg,ptDpg,ptDpb)
      omz=ptDpg/(ptDpb+ptDpg-pbDpg)  ! for some reason this is not (1-z) as defined in the paper...
      rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
      z=1d0-omz
      y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2
      Dipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
      Dipole = Dipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )

      MomDKx(1:4,1:3) = MomDKTd(1:4,1:3)
      MomDKx(1:4,4:6) = MomDK(1:4,5:7)
      call Kinematics_TTBAR(.false.,MomExt,MomDKx,applyPSCut,NBin)
      if( applyPSCut ) then
          goto 25
      endif


      call TopDecay(ExtParticle(1),DK_LO,MomDKTd(1:4,1:4))
      Spi(1:4) = ExtParticle(1)%Pol(1:4)
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,-1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0

      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,+1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
      DKFac_ATop = Msq_T_BWENU * Dipole

      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,5:7))
      BarSpi(1:4) = ExtParticle(2)%Pol(1:4)
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,-1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0

      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,+1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
      DKFac_Top = Msq_T_BWENU


      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol*PreFac*PSWgt*PSWgt2*PSWgt3 * DKFac_ATop*DKFac_Top))
      enddo
! ENDIF
      EvalCS_NLODK_ttb_noSC = EvalCS_NLODK_ttb_noSC + dble(LO_Res_Unpol*PreFac*PSWgt*PSWgt2*PSWgt3 * DKFac_ATop*DKFac_Top)


!       print *, "Di",dble(LO_Res_Unpol*PreFac*PSWgt*PSWgt2*PSWgt3 * DKFac_ATop*DKFac_Top)
!       pause


25 continue
endif



if( DKRE_switch.eq.0 .or. DKRE_switch.eq.2 ) then
!----------------------------------------
! real gluon emission for top decay     |
!----------------------------------------
   call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:8),.false.,MomDK(1:4,1:3),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(12:18),.true.,MomDK(1:4,4:7),PSWgt3)


   call CheckSing(MomDK(1:4,4:7),applySingCut)
   if( applySingCut) then
      EvalCS_NLODK_ttb_noSC = EvalCS_NLODK_ttb_noSC/VgsWgt
      return
   endif
   call Kinematics_TTBAR(.true.,MomExt,MomDK,applyPSCut,NBin)
   if( applyPSCut ) then
      goto 26
   endif

! !       MomP(1:4) = MomExt(1:4,4)-MomDK(1:4,4)
! !       rsq = (MomP(1:4).dot.MomP(1:4))/m_Top**2
!       rsq = m_W**2/m_Top**2
!       ptb = MomExt(1:4,4).dot.MomDK(1:4,4)
!       pbg = MomDK(1:4,4).dot.MomDK(1:4,7)
!       ptg = MomExt(1:4,4).dot.MomDK(1:4,7)
!
!       MsqRE_T_BW = 2d0*(alpha_s4Pi*RunFactor)*CF*GF/dsqrt(2d0)*m_top**2 *(           &
!                  + m_top**2*(1d0-rsq)*(1d0+2d0*rsq)*(2d0*ptb/pbg/ptg-m_top**2/ptg**2)  &
!                  + 8d0*rsq + 2d0*(ptg-pbg)**2/pbg/ptg*(1d0+2d0*rsq)        )
!
!       MsqLO_T_BW = dsqrt(2d0)*GF*m_top**4*(1d0-m_w**2/m_top**2)*(1d0+2d0*m_w**2/m_top**2)  ! LO sq.mat.el. for T -> b W
!       Msq_W_ENU  = dsqrt(2d0)*GF*m_W**4/3d0*4d0  ! sq.mat.el. for W -> l nu
!
!       DKFac_Top  = MsqRE_T_BW/(2d0*Ga_Top(0)*m_Top) * Msq_W_ENU/(2d0*Ga_W*m_W)          ! real corr. to Top
!       DKFac_ATop = MsqLO_T_BW/(2d0*Ga_Top(0)*m_Top) * Msq_W_ENU/(2d0*Ga_W*m_W)          ! LO corr. to ATop

      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      Spi(1:4) = ExtParticle(1)%Pol(1:4)
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,-1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)

      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,+1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)
      DKFac_ATop =  Msq_T_BWENU/2d0



      call TopDecay(ExtParticle(2),DK_RE_T,MomDK(1:4,4:7),GluonHel=-1)
      BarSpi(1:4) = ExtParticle(2)%Pol(1:4)
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,-1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,+1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)

      call TopDecay(ExtParticle(2),DK_RE_T,MomDK(1:4,4:7),GluonHel=+1)
      BarSpi(1:4) = ExtParticle(2)%Pol(1:4)
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,-1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,+1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)
      DKFac_Top = Msq_T_BWENU/2d0


      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol*PreFac*PSWgt*PSWgt2*PSWgt3 * DKFac_Top*DKFac_ATop))
      enddo
! ENDIF

!       print *, "Ma",dble(LO_Res_Unpol*PreFac*PSWgt*PSWgt2*PSWgt3 * DKFac_Top*DKFac_ATop)

      EvalCS_NLODK_ttb_noSC = EvalCS_NLODK_ttb_noSC + dble(LO_Res_Unpol*PreFac*PSWgt*PSWgt2*PSWgt3 * DKFac_Top*DKFac_ATop)


26 continue
!-------------------------------------
! dipole subtraction for top-decay   |
!-------------------------------------
  call WTransform(MomDK(1:4,4:7),MomDKTd(1:4,1:3),pbDpg,ptDpg,ptDpb)
  omz=ptDpg/(ptDpb+ptDpg-pbDpg)
  rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
  z=1d0-omz
  y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2
  Dipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
  Dipole = Dipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )


   MomDKx(1:4,1:3) = MomDK(1:4,1:3)
   MomDKx(1:4,4:6) = MomDKTd(1:4,1:3)
   call Kinematics_TTBAR(.false.,MomExt,MomDKx,applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_NLODK_ttb_noSC = EvalCS_NLODK_ttb_noSC/VgsWgt
      return
   endif

!       MsqLO_T_BW = dsqrt(2d0)*GF*m_top**4*(1d0-m_w**2/m_top**2)*(1d0+2d0*m_w**2/m_top**2)  ! LO sq.mat.el. for T -> b W
!       Msq_W_ENU = dsqrt(2d0)*GF*m_W**4/3d0*4d0  ! sq.mat.el. for W -> l nu
!
!       DKFac_Top  = Dipole* MsqLO_T_BW/(2d0*Ga_Top(0)*m_Top) * Msq_W_ENU/(2d0*Ga_W*m_W)  ! dipole for to Top
!       DKFac_ATop = MsqLO_T_BW/(2d0*Ga_Top(0)*m_Top) * Msq_W_ENU/(2d0*Ga_W*m_W)          ! LO corr. to ATop

      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:4))
      Spi(1:4) = ExtParticle(1)%Pol(1:4)
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,-1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0

      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,+1,BarSpi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
      DKFac_ATop = Msq_T_BWENU

      call TopDecay(ExtParticle(2),DK_LO,MomDKTd(1:4,1:3))
      BarSpi(1:4) = ExtParticle(2)%Pol(1:4)
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,-1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = M_T_BWENU*dconjg(M_T_BWENU)/2d0

      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,+1,Spi(1:4))
      M_T_BWENU = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BWENU = Msq_T_BWENU + M_T_BWENU*dconjg(M_T_BWENU)/2d0
      DKFac_Top = Msq_T_BWENU * Dipole


      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol*PreFac*PSWgt*PSWgt2*PSWgt3 * DKFac_Top*DKFac_ATop))
      enddo
! ENDIF

!       print *, "Di",dble(LO_Res_Unpol*PreFac*PSWgt*PSWgt2*PSWgt3 * DKFac_ATop*DKFac_Top)
!       pause

      EvalCS_NLODK_ttb_noSC = EvalCS_NLODK_ttb_noSC + dble(LO_Res_Unpol)*PreFac*PSWgt*PSWgt2*PSWgt3 * DKFac_Top*DKFac_ATop


ENDIF
endif


   EvalCS_NLODK_ttb_noSC = EvalCS_NLODK_ttb_noSC/VgsWgt


return
END FUNCTION







FUNCTION EvalCS_NLODK_ttb(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModHadrWDecay
use ModJPsiFrag
implicit none
real(8) ::  EvalCS_NLODK_ttb,yRnd(1:VegasMxDim),VgsWgt
complex(8) :: LO_Res_Pol,LO_Res_Unpol,Dip_Res_Unpol,NLO_Res_Pol,NLO_Res_UnPol
complex(8) :: TreeResult(1:NumBornAmps),DKResult(1:NumBornAmps)
integer :: iHel,jHel,kHel,GluHel,iPrimAmp,jPrimAmp,ndip
real(8) :: EHat,PSWgt1,PSWgt2,PSWgt3,ISFac,dip_res_w,xFrag
real(8) :: MomExt(1:4,1:NumExtParticles),MomDKTd(1:4,1:3),MomDK(1:4,1:7),MomDKx(1:4,1:7)
logical :: applyPSCut,applySingCut
real(8) :: MG_MOM(0:3,1:NumExtParticles)
integer, parameter :: ParityFlip=1
real(8) :: tau,eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,RunFactor,r_sc
real(8) :: pdf(-6:6,1:2),ColCorrLO(1:NumBornAmps,1:NumBornAmps)
integer :: NBin(1:NumMaxHisto),NHisto
real(8) :: pbDpg,ptDpg,ptDpb,z,omz,Dipole,rsq,y
real(8), parameter :: CF=4d0/3d0
real(8) :: MomBoost(1:4),MomLep1(1:4),MomLep2(1:4)
real(8) :: ThresholdCutOff = 1.0d0
include "vegas_common.f"
logical, parameter :: include_WCorr=.true.

! print *, "fixed yrnd"
! yrnd(1)= 0.120417488562138
! yrnd(2)=0.467848163315956
! yrnd(3)=0.314100383461500
! yrnd(4)=  4.800610991567660E-002
! yrnd(5)=  6.554446535443172E-002
! yrnd(6)=  0.105829211932528
! yrnd(7)=  0.492114467077011
! yrnd(8)=  2.928262810655061E-002
! yrnd(9)=  0.153130586796035
! yrnd(10)=  0.426867230994100
! yrnd(11)=  0.381774851764448
! yrnd(12)=  6.696560190383624E-003
! yrnd(13)=  4.668414455218438E-002
! yrnd(14)=  0.356730263846335
! yrnd(15)=  0.351952110348247
! yrnd(16)=  8.698269659047150E-002
! yrnd(17)=  0.171564949989116
! yrnd(18)=  0.357551317828499



EvalCS_NLODK_ttb = 0d0
IF( Process.EQ.1 ) THEN
    ColCorrLO(1:NumBornAmps,1:NumBornAmps) = ColLO_ttbgg(1:NumBornAmps,1:NumBornAmps)
ELSEIF( Process.EQ.2 ) THEN
    ColCorrLO(1:NumBornAmps,1:NumBornAmps) = ColLO_ttbqqb(1:NumBornAmps,1:NumBornAmps)
ENDIF

  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top * ThresholdCutOff ) then
      EvalCS_NLODK_ttb = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to2(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt1)
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
   ISFac = MomCrossing(MomExt)



IF( CORRECTION.EQ.4 ) THEN
if( DKRE_switch.eq.0 .or. DKRE_switch.eq.1 ) then
!----------------------------------------
! one loop correction to Anti-top decay |
!----------------------------------------
   call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:8),.false.,MomDK(1:4,1:3),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(9:12),.false.,MomDK(1:4,4:6),PSWgt3)
   xFrag = yRnd(13)
   if( TopDecays.eq.6 ) then ! evaluate LO fragmentation function
      PSWgt2 = PSWgt2 * FF(xFrag)
   endif
   call Kinematics_TTBAR(.false.,MomExt,MomDK,applyPSCut,NBin,xJPsiFrag=xFrag)
   if( applyPSCut ) then
      EvalCS_NLODK_ttb = 0d0
      goto 12
   endif

   !call InitCurrCache()
   call setPDFs(eta1,eta2,MuFac,pdf)
IF( PROCESS.EQ.1 ) THEN
   PDFFac = pdf(0,1) * pdf(0,2)
ELSEIF( PROCESS.EQ.2 ) THEN
   PDFFac =       ( pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
                  + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
                  + pdf(Bot_,1)*pdf(ABot_,2) ) +                         &
                  ( pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
                  + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
                  + pdf(Bot_,2)*pdf(ABot_,1) )
ENDIF
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt * PDFFac
   RunFactor = RunAlphaS(2,MuRen)

   LO_Res_Unpol = (0d0,0d0)
   NLO_Res_UnPol= (0d0,0d0)
   do iHel=1,NumHelicities ! loop over initial state chiralities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()

      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
          TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
      enddo

      call TopDecay(ExtParticle(1),DK_1L_T,MomDK(1:4,1:3),xIntDip=yRnd(13:14))   ! xIntDip is dummy in case of TopDecays.ne.5
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
          DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
      enddo
      NLO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          NLO_Res_Pol = NLO_Res_Pol + ParityFlip*ColCorrLO(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(DKResult(jPrimAmp)) )
      enddo
      enddo

      if( include_WCorr .and. (TopDecays.eq.2 .or. TopDecays.eq.4 .or. TopDecays.eq.6) ) then  !  virt.corr. to hadronic W decay
          call TopDecay(ExtParticle(1),DK_1L_Q,MomDK(1:4,1:3))
          do iPrimAmp=1,NumBornAmps
              call EvalTree(BornAmps(iPrimAmp))
              DKResult(iPrimAmp) =  BornAmps(iPrimAmp)%Result
          enddo
          do jPrimAmp=1,NumBornAmps
          do iPrimAmp=1,NumBornAmps
              NLO_Res_Pol = NLO_Res_Pol + ParityFlip*ColCorrLO(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(DKResult(jPrimAmp)) )
          enddo
          enddo
      endif

      NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
   enddo!helicity loop
!  normalization
   NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
   EvalCS_NLODK_ttb = dble(NLO_Res_Unpol)



IF( ObsSet.EQ.7 ) THEN
   do NHisto=1,3
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo
   call intoHisto(4,1,MInv_LB*dble(NLO_Res_Unpol))
   call intoHisto(5,1,MInv_LB**2*dble(NLO_Res_Unpol))
   call intoHisto(6,1,xFrag*dble(NLO_Res_Unpol))
ELSE
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo
ENDIF


endif


12 continue


!-------------------------------------
! one loop correction to top-decay   |
!-------------------------------------
if( DKRE_switch.eq.0 .or. DKRE_switch.eq.2 ) then

   call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:8),.false.,MomDK(1:4,1:3),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(9:12),.false.,MomDK(1:4,4:6),PSWgt3)
   xFrag = yRnd(13)
   if( TopDecays.eq.5 ) then ! evaluate LO fragmentation function
      PSWgt3 = PSWgt3 * FF(xFrag)
   endif
   call Kinematics_TTBAR(.false.,MomExt,MomDK,applyPSCut,NBin,xJPsiFrag=xFrag)
   if( applyPSCut ) then
      EvalCS_NLODK_ttb = EvalCS_NLODK_ttb/VgsWgt
      return
   endif

   !call InitCurrCache()
   call setPDFs(eta1,eta2,MuFac,pdf)
IF( PROCESS.EQ.1 ) THEN
   PDFFac = pdf(0,1) * pdf(0,2)
ELSEIF( PROCESS.EQ.2 ) THEN
   PDFFac =       ( pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
                  + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
                  + pdf(Bot_,1)*pdf(ABot_,2) ) +                         &
                  ( pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
                  + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
                  + pdf(Bot_,2)*pdf(ABot_,1) )
ENDIF
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt * PDFFac
   RunFactor = RunAlphaS(2,MuRen)

   LO_Res_Unpol = (0d0,0d0)
   NLO_Res_UnPol= (0d0,0d0)
   do iHel=1,NumHelicities ! loop over initial state chiralities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
          TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
      enddo

      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_1L_T,MomDK(1:4,4:6),xIntDip=yRnd(13:14))   ! xIntDip is dummy in case of TopDecays.ne.6
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
          DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
      enddo
      NLO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          NLO_Res_Pol = NLO_Res_Pol + ParityFlip*ColCorrLO(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(DKResult(jPrimAmp)) )
      enddo
      enddo

      if( include_WCorr .and. (TopDecays.eq.2 .or. TopDecays.eq.3 .or. TopDecays.eq.5) ) then  !  virt.corr. to hadronic W decay
          call TopDecay(ExtParticle(2),DK_1L_Q,MomDK(1:4,4:6))
          do iPrimAmp=1,NumBornAmps
              call EvalTree(BornAmps(iPrimAmp))
              DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
          enddo
          do jPrimAmp=1,NumBornAmps
          do iPrimAmp=1,NumBornAmps
              NLO_Res_Pol = NLO_Res_Pol + ParityFlip*ColCorrLO(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(DKResult(jPrimAmp)) )
          enddo
          enddo
      endif

      NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
   enddo!helicity loop
!  normalization
   NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
   EvalCS_NLODK_ttb = EvalCS_NLODK_ttb + NLO_Res_Unpol


IF( ObsSet.EQ.7 ) THEN

   do NHisto=1,3
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo
   call intoHisto(4,1,MInv_LB*dble(NLO_Res_Unpol))
   call intoHisto(5,1,MInv_LB**2*dble(NLO_Res_Unpol))
   call intoHisto(6,1,xFrag*dble(NLO_Res_Unpol))

ELSE
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo
ENDIF

endif




ELSEIF( CORRECTION.EQ.5 ) THEN
if( DKRE_switch.eq.0 .or. DKRE_switch.eq.1 ) then
!----------------------------------------
! real gluon emission for Anti-top decay |
!----------------------------------------
   call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:11),.true.,MomDK(1:4,1:4),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(12:15),.false.,MomDK(1:4,5:7),PSWgt3)
   if( TopDecays.eq.5 .or. TopDecays.eq.6 ) then ! evaluate LO fragmentation function
      xFrag = yRnd(19)
      PSWgt2 = PSWgt2 * FF(xFrag)
   endif



   !call InitCurrCache()
   call setPDFs(eta1,eta2,MuFac,pdf)
IF( PROCESS.EQ.1 ) THEN
   PDFFac = pdf(0,1) * pdf(0,2)
ELSEIF( PROCESS.EQ.2 ) THEN
   PDFFac =       ( pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
                  + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
                  + pdf(Bot_,1)*pdf(ABot_,2) ) +                         &
                  ( pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
                  + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
                  + pdf(Bot_,2)*pdf(ABot_,1) )
ENDIF
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt * PDFFac
   RunFactor = RunAlphaS(2,MuRen)

   call CheckSing(MomDK(1:4,1:4),applySingCut)
   if( applySingCut) then
      EvalCS_NLODK_ttb = 0d0
      goto 13
   endif

   MomDKx(1:4,1:3) = MomDK(1:4,1:3)
   MomDKx(1:4,4:6) = MomDK(1:4,5:7)
   MomDKx(1:4,7)   = MomDK(1:4,4)
   call Kinematics_TTBAR(.true.,MomExt,MomDKx,applyPSCut,NBin,xJPsiFrag=xFrag)
   if( applyPSCut ) then
      goto 14
   endif

   LO_Res_Unpol = (0d0,0d0)
   do iHel=1,NumHelicities ! loop over initial state chiralities
   do GluHel=1,-1,-2 ! loop over additional gluon chiralities from decay
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_RE_T,MomDK(1:4,1:4),GluonHel=GluHel)
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,5:7))

      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo
      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol

   enddo!helicity loop
   enddo!helicity loop

!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
   EvalCS_NLODK_ttb = dble(LO_Res_Unpol)

IF( ObsSet.EQ.7 ) THEN
   do NHisto=1,3
      call intoHisto(NHisto,NBin(NHisto),EvalCS_NLODK_ttb)
   enddo
   call intoHisto(4,1,MInv_LB*EvalCS_NLODK_ttb)
   call intoHisto(5,1,MInv_LB**2*EvalCS_NLODK_ttb)
   call intoHisto(6,1,xFrag*EvalCS_NLODK_ttb)
ELSE

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_NLODK_ttb)
   enddo
ENDIF


14 continue
!-------------------------------------
! dipole subtraction for Atop-decay |
!-------------------------------------
  call WTransform(MomDK(1:4,1:4),MomDKTd(1:4,1:3),pbDpg,ptDpg,ptDpb)
  omz=ptDpg/(ptDpb+ptDpg-pbDpg)  !  for some reason this is not (1-z) as defined in the paper...
  rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
  z=1d0-omz
  y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2

  Dipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
  Dipole = Dipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )

   MomDKx(1:4,1:3) = MomDKTd(1:4,1:3)
   MomDKx(1:4,4:6) = MomDK(1:4,5:7)
   call Kinematics_TTBAR(.false.,MomExt,MomDKx,applyPSCut,NBin,xJPsiFrag=xFrag)
   if( applyPSCut ) then
      goto 13
   endif

   Dip_Res_Unpol= (0d0,0d0)
   do iHel=1,NumHelicities ! loop over initial state chiralities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()

      call TopDecay(ExtParticle(1),DK_LO,MomDKTd(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,5:7))

      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo
      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      Dip_Res_UnPol = Dip_Res_UnPol + LO_Res_Pol
   enddo!helicity loop

!  normalization
   Dip_Res_Unpol = Dip_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * Dipole * PreFac
   EvalCS_NLODK_ttb = EvalCS_NLODK_ttb + Dip_Res_Unpol


IF( ObsSet.EQ.7 ) THEN

   do NHisto=1,3
      call intoHisto(NHisto,NBin(NHisto),dble(Dip_Res_Unpol))
   enddo
   call intoHisto(4,1,MInv_LB*dble(Dip_Res_Unpol))
   call intoHisto(5,1,MInv_LB**2*dble(Dip_Res_Unpol))
   call intoHisto(6,1,xFrag*dble(Dip_Res_Unpol))
ELSE
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(Dip_Res_Unpol))
   enddo
ENDIF





13 continue


!             print *, MomDK(1,4)/EHat,EvalCS_NLODK_ttb/dble(Dip_Res_Unpol)
!             print *, (MomDK(1:4,1).dot.MomDK(1:4,4))/EHat**2,EvalCS_NLODK_ttb/dble(Dip_Res_Unpol)
!             EvalCS_NLODK_ttb=EvalCS_NLODK_ttb/Vgswgt
!             pause
!             return


   if( include_WCorr .and. (TopDecays.eq.2 .or. TopDecays.eq.4 .or. TopDecays.eq.6) ) then
          call EvalPhasespace_TopDecay2(MomExt(1:4,3),yRnd(5:11),.true.,MomDK(1:4,1:4),PSWgt2)
          call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(12:15),.false.,MomDK(1:4,5:7),PSWgt3)
          if( TopDecays.eq.6 ) then ! evaluate LO fragmentation function
              xFrag = yRnd(19)
              PSWgt2 = PSWgt2 * FF(xFrag)
          endif
          if( PSWgt2 .eq. 0d0 ) then  ! this rejects "too singular" events, similar to CheckSing
              EvalCS_NLODK_ttb = 0d0
              goto 19
          endif

          PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt * PDFFac
          RunFactor = RunAlphaS(2,MuRen)

          MomDKx(1:4,1:3) = MomDK(1:4,1:3)
          MomDKx(1:4,4:6) = MomDK(1:4,5:7)
          MomDKx(1:4,7)   = MomDK(1:4,4)
          call Kinematics_TTBAR(.true.,MomExt,MomDKx,applyPSCut,NBin,xJPsiFrag=xFrag)
          if( applyPSCut ) then
              goto 16
          endif

          LO_Res_Unpol = (0d0,0d0)
          do iHel=1,NumHelicities ! loop over initial state chiralities
          do GluHel=1,-1,-2 ! loop over additional gluon chiralities from decay
              call HelCrossing(Helicities(iHel,1:NumExtParticles))
              call SetPolarizations()
                  call TopDecay(ExtParticle(1),DK_RE_Q,MomDK(1:4,1:4),GluonHel=GluHel)
                  call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,5:7))
                  do iPrimAmp=1,NumBornAmps
                      call EvalTree(BornAmps(iPrimAmp))
                  enddo
                  LO_Res_Pol = (0d0,0d0)
                  do jPrimAmp=1,NumBornAmps
                  do iPrimAmp=1,NumBornAmps
                      LO_Res_Pol = LO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
                  enddo
                  enddo
                  LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
          enddo!helicity loop
          enddo!helicity loop

        !  normalization
          LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
          EvalCS_NLODK_ttb = EvalCS_NLODK_ttb + dble(LO_Res_Unpol)


IF( ObsSet.EQ.7 ) THEN

   do NHisto=1,3
      call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol))
   enddo
   call intoHisto(4,1,MInv_LB*dble(LO_Res_Unpol))
   call intoHisto(5,1,MInv_LB**2*dble(LO_Res_Unpol))
   call intoHisto(6,1,xFrag*dble(LO_Res_Unpol))
ELSE
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol))
   enddo
ENDIF


16 continue
!             print *, (MomDK(1:4,3).dot.MomDK(1:4,4))/EHat**2,   EvalCS_NLODK_ttb
          do ndip=1,2        !---there are two dipoles
              call wdec_trans(ndip,MomDK(1:4,1:4),MomDKTd(1:4,1:3),alpha_DKWff,dip_res_w)
              if( dip_res_w.eq.0d0 ) goto 18
              Dipole = - alpha_s4Pi*RunFactor * CF * dip_res_w
              MomDKx(1:4,1:3) = MomDKTd(1:4,1:3)
              MomDKx(1:4,4:6) = MomDK(1:4,5:7)

              call Kinematics_TTBAR(.false.,MomExt,MomDKx,applyPSCut,NBin,xJPsiFrag=xFrag)    ! this needs to be changed; for extra partons
              if( applyPSCut ) then
                  goto 18
              endif

              Dip_Res_Unpol= (0d0,0d0)
              do iHel=1,NumHelicities ! loop over initial state chiralities
                  call HelCrossing(Helicities(iHel,1:NumExtParticles))
                  call SetPolarizations()
                  call TopDecay(ExtParticle(1),DK_LO,MomDKTd(1:4,1:3))
                  call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,5:7))

                  do iPrimAmp=1,NumBornAmps
                      call EvalTree(BornAmps(iPrimAmp))
                  enddo
                  LO_Res_Pol = (0d0,0d0)
                  do jPrimAmp=1,NumBornAmps
                  do iPrimAmp=1,NumBornAmps
                      LO_Res_Pol = LO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
                  enddo
                  enddo
                  Dip_Res_UnPol = Dip_Res_UnPol + LO_Res_Pol
              enddo!helicity loop

!  normalization
            Dip_Res_Unpol = Dip_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * Dipole * PreFac
            EvalCS_NLODK_ttb = EvalCS_NLODK_ttb + Dip_Res_Unpol

IF( ObsSet.EQ.7 ) THEN
   do NHisto=1,3
      call intoHisto(NHisto,NBin(NHisto),dble(Dip_Res_Unpol))
   enddo
   call intoHisto(4,1,MInv_LB*dble(Dip_Res_Unpol))
   call intoHisto(5,1,MInv_LB**2*dble(Dip_Res_Unpol))
   call intoHisto(6,1,xFrag*dble(Dip_Res_Unpol))
ELSE
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(Dip_Res_Unpol))
   enddo
ENDIF

        18 continue
        enddo   !dipole loop

!             print *, MomDK(1,4)/EHat,  EvalCS_NLODK_ttb,  dble(Dip_Res_Unpol),   EvalCS_NLODK_ttb/dble(Dip_Res_Unpol)
!             print *, (MomDK(1:4,3).dot.MomDK(1:4,4))/EHat**2,   dble(Dip_Res_Unpol)
!             pause

   endif

endif

19 continue



if( DKRE_switch.eq.0 .or. DKRE_switch.eq.2 ) then
!-------------------------------------
! real gluon emission for top-decay  |
!-------------------------------------
   call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:8),.false.,MomDK(1:4,1:3),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(12:18),.true.,MomDK(1:4,4:7),PSWgt3)

   !call InitCurrCache()
   call setPDFs(eta1,eta2,MuFac,pdf)
IF( PROCESS.EQ.1 ) THEN
   PDFFac = pdf(0,1) * pdf(0,2)
ELSEIF( PROCESS.EQ.2 ) THEN
   PDFFac =       ( pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
                  + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
                  + pdf(Bot_,1)*pdf(ABot_,2) ) +                         &
                  ( pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
                  + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
                  + pdf(Bot_,2)*pdf(ABot_,1) )
ENDIF
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt * PDFFac
   RunFactor = RunAlphaS(2,MuRen)

   if( TopDecays.eq.5 ) then ! evaluate LO fragmentation function
      xFrag = yRnd(19)
      PreFac = PreFac * FF(xFrag)
   elseif( TopDecays.eq.6 ) then
      z = 2d0*(MomExt(1:4,4).dot.MomDK(1:4,4))/m_top**2/(one-M_W**2/m_Top**2)
      xFrag = yRnd(19)/z
      PreFac  = PreFac * FF(xFrag)/z
   endif

   call CheckSing(MomDK(1:4,4:7),applySingCut)
   if( applySingCut ) then
      goto 17
   endif

   call Kinematics_TTBAR(.true.,MomExt,MomDK,applyPSCut,NBin,xJPsiFrag=xFrag)
   if( applyPSCut ) then
      goto 15
   endif

   LO_Res_Unpol = (0d0,0d0)
   do iHel=1,NumHelicities ! loop over initial state chiralities
   do GluHel=1,-1,-2 ! loop over additional gluon chiralities from decay
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_RE_T,MomDK(1:4,4:7),GluonHel=GluHel)
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo
      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
   enddo!helicity loop
   enddo!helicity loop

!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
   EvalCS_NLODK_ttb = EvalCS_NLODK_ttb + LO_Res_Unpol

IF( ObsSet.EQ.7 ) THEN
   do NHisto=1,3
      call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol))
   enddo
   call intoHisto(4,1,MInv_LB*dble(LO_Res_Unpol))
   call intoHisto(5,1,MInv_LB**2*dble(LO_Res_Unpol))
   call intoHisto(6,1,xFrag*dble(LO_Res_Unpol))
ELSE
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol))
   enddo
ENDIF

15 continue

!-------------------------------------
! dipole subtraction for top-decay   |
!-------------------------------------
  call WTransform(MomDK(1:4,4:7),MomDKTd(1:4,1:3),pbDpg,ptDpg,ptDpb)
  omz=ptDpg/(ptDpb+ptDpg-pbDpg)
  rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
  z=1d0-omz
  y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2
  Dipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
  Dipole = Dipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt * PDFFac
   if( TopDecays.eq.5 ) then ! evaluate LO fragmentation function
      xFrag = yRnd(19)
      PreFac = PreFac * FF(xFrag)
   elseif( TopDecays.eq.6 ) then
      xFrag = yRnd(19)/z
      PreFac  = PreFac * FF(xFrag)/z
   endif

   MomDKx(1:4,1:3) = MomDK(1:4,1:3)
   MomDKx(1:4,4:6) = MomDKTd(1:4,1:3)
   call Kinematics_TTBAR(.false.,MomExt,MomDKx,applyPSCut,NBin,xJPsiFrag=yRnd(19))  ! MomDKx(4)= MomDKTd(1:4,1) includes already z factor, so p_JPsi=yRnd(19)*ptilde_b
   if( applyPSCut ) then
      goto 17
   endif

   Dip_Res_Unpol= (0d0,0d0)
   do iHel=1,NumHelicities ! loop over initial state chiralities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()

      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_LO,MomDKTd(1:4,1:3))

      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo
      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      Dip_Res_UnPol = Dip_Res_UnPol + LO_Res_Pol
   enddo!helicity loop

!  normalization
   Dip_Res_Unpol = Dip_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * Dipole * PreFac
   EvalCS_NLODK_ttb = EvalCS_NLODK_ttb + dble(Dip_Res_Unpol)

!             print *, MomDK(1,7)/EHat,pbDpg/EHat**2,EvalCS_NLODK_ttb/dble(Dip_Res_Unpol)
!             print *, dble(LO_Res_Unpol),dble(Dip_Res_Unpol),FF(xFrag)/z
!             pause


IF( ObsSet.EQ.7 ) THEN
   do NHisto=1,3
      call intoHisto(NHisto,NBin(NHisto),dble(Dip_Res_Unpol))
   enddo
   call intoHisto(4,1,MInv_LB*dble(Dip_Res_Unpol))
   call intoHisto(5,1,MInv_LB**2*dble(Dip_Res_Unpol))
   call intoHisto(6,1,xFrag*dble(Dip_Res_Unpol))
ELSE
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(Dip_Res_Unpol))
   enddo
ENDIF



17 continue



      if( include_WCorr .and. (TopDecays.eq.2 .or. TopDecays.eq.3 .or. TopDecays.eq.5) ) then
          call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:8),.false.,MomDK(1:4,1:3),PSWgt2)
          call EvalPhasespace_TopDecay2(MomExt(1:4,4),yRnd(12:18),.true.,MomDK(1:4,4:7),PSWgt3)
          if( TopDecays.eq.5 ) then ! evaluate LO fragmentation function
              xFrag = yRnd(19)
              PSWgt3 = PSWgt3 * FF(xFrag)
          endif
          if( PSWgt3 .eq. 0d0 ) then
              EvalCS_NLODK_ttb = EvalCS_NLODK_ttb/VgsWgt
              return
          endif

          PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt * PDFFac
          RunFactor = RunAlphaS(2,MuRen)

          call Kinematics_TTBAR(.true.,MomExt,MomDK,applyPSCut,NBin,xJPsiFrag=xFrag)
          if( applyPSCut ) then
              goto 21
          endif

          LO_Res_Unpol = (0d0,0d0)
          do iHel=1,NumHelicities ! loop over initial state chiralities
          do GluHel=1,-1,-2 ! loop over additional gluon chiralities from decay
              call HelCrossing(Helicities(iHel,1:NumExtParticles))
              call SetPolarizations()
              call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
              call TopDecay(ExtParticle(2),DK_RE_Q,MomDK(1:4,4:7),GluonHel=GluHel)

              do iPrimAmp=1,NumBornAmps
                  call EvalTree(BornAmps(iPrimAmp))
              enddo
              LO_Res_Pol = (0d0,0d0)
              do jPrimAmp=1,NumBornAmps
              do iPrimAmp=1,NumBornAmps
                  LO_Res_Pol = LO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
              enddo
              enddo
              LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
          enddo!helicity loop
          enddo!helicity loop
        !  normalization
          LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
          EvalCS_NLODK_ttb = EvalCS_NLODK_ttb + LO_Res_Unpol

          IF( ObsSet.EQ.7 ) THEN
                  do NHisto=1,3
                      call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol))
                  enddo
                  call intoHisto(4,1,MInv_LB*dble(LO_Res_Unpol))
                  call intoHisto(5,1,MInv_LB**2*dble(LO_Res_Unpol))
                  call intoHisto(6,1,xFrag*dble(LO_Res_Unpol))
           ELSE
                  do NHisto=1,NumHistograms
                      call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol))
                  enddo
           ENDIF

21 continue

          do ndip=1,2
                call wdec_trans(ndip,MomDK(1:4,4:7),MomDKTd(1:4,1:3),alpha_DKWff,dip_res_w)
                if( dip_res_w.eq.0d0 ) goto 20
                Dipole = - alpha_s4Pi*RunFactor * CF * dip_res_w

                MomDKx(1:4,1:3) = MomDK(1:4,1:3)
                MomDKx(1:4,4:6) = MomDKTd(1:4,1:3)
                call Kinematics_TTBAR(.false.,MomExt,MomDKx,applyPSCut,NBin,xJPsiFrag=xFrag)
                if( applyPSCut ) then
                    goto 20
                    return
                endif

                Dip_Res_Unpol= (0d0,0d0)
                do iHel=1,NumHelicities ! loop over initial state chiralities
                    call HelCrossing(Helicities(iHel,1:NumExtParticles))
                    call SetPolarizations()

                    call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
                    call TopDecay(ExtParticle(2),DK_LO,MomDKTd(1:4,1:3))

                    do iPrimAmp=1,NumBornAmps
                        call EvalTree(BornAmps(iPrimAmp))
                    enddo
                    LO_Res_Pol = (0d0,0d0)
                    do jPrimAmp=1,NumBornAmps
                    do iPrimAmp=1,NumBornAmps
                        LO_Res_Pol = LO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
                    enddo
                    enddo
                    Dip_Res_UnPol = Dip_Res_UnPol + LO_Res_Pol
                enddo!helicity loop

              !  normalization
                Dip_Res_Unpol = Dip_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * Dipole * PreFac
                EvalCS_NLODK_ttb = EvalCS_NLODK_ttb + Dip_Res_Unpol

                IF( ObsSet.EQ.7 ) THEN
                      do NHisto=1,3
                          call intoHisto(NHisto,NBin(NHisto),dble(Dip_Res_Unpol))
                      enddo
                      call intoHisto(4,1,MInv_LB*dble(Dip_Res_Unpol))
                      call intoHisto(5,1,MInv_LB**2*dble(Dip_Res_Unpol))
                      call intoHisto(6,1,xFrag*dble(Dip_Res_Unpol))
                ELSE
                      do NHisto=1,NumHistograms
                          call intoHisto(NHisto,NBin(NHisto),dble(Dip_Res_Unpol))
                      enddo
                ENDIF

                20 continue
          enddo ! dipole loop

!             print *, (MomDK(1:4,7).dot.MomDK(1:4,6))/EHat**2,EvalCS_NLODK_ttb/LO_Res_Unpol
!             pause
   endif


endif
ENDIF


   EvalCS_NLODK_ttb = EvalCS_NLODK_ttb/VgsWgt
return
END FUNCTION







FUNCTION EvalCS_1L_gggggg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
implicit none
real(8) ::  EvalCS_1L_gggggg,yRnd(1:VegasMxDim),VgsWgt
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1)
integer :: iHel,iPrimAmp,jPrimAmp
real(8) :: SpinAvg,ColorAvg,EHat,Mu,PSWgt,CrossingFactor,AccPoles
real(8) :: MomExt(1:4,1:NumExtParticles)
logical :: applyPSCut

   include 'misc/global_import'

   EHat = 10.d0
   Mu   = EHat

!   call EvalPhaseSpace_2to2(EHat,yRnd(3:4),MomExt(1:4,1:NumExtParticles),PSWgt)
!   call EvalPhaseSpace_2to3(EHat,yRnd(3:7),MomExt(1:4,1:NumExtParticles),PSWgt)
!   call EvalPhaseSpace_2to4(EHat,yRnd(1:8),MomExt(1:4,1:NumExtParticles),PSWgt)
   call EvalPhaseSpace_2to5(EHat,yRnd(1:11),MomExt(1:4,1:NumExtParticles),PSWgt)
!    call CheckPSCuts(MomExt,applyPSCut)
!    if( applyPSCut ) return
!    call InitCurrCache()
   CrossingFactor = MomCrossing(MomExt)
   call SetPropagators()

!    LO_Res_Unpol             = (0d0,0d0)
!    NLO_Res_Unpol(-2:1)      = (0d0,0d0)
!    NLO_Res_Unpol_Ferm(-2:1) = (0d0,0d0)

!------------ Helicity summation --------------
   do iHel=1,NumHelicities
!      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()

!------------ LO --------------
      do iPrimAmp=1,1
          call EvalTree(BornAmps(iPrimAmp))
      enddo
      print *, BornAmps(iPrimAmp)%Result
!       LO_Res_Pol = (0d0,0d0)
!       do jPrimAmp=1,NumBornAmps
!       do iPrimAmp=1,NumBornAmps
!           LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * 2d0*BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
!       enddo
!       enddo
!       LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol

!------------ NLO bosonic --------------
      do iPrimAmp=1,1
          call SetKirill(PrimAmps(iPrimAmp))

          BornAmps(iPrimAmp)%Result = (0d0,1d0)*BornAmps(iPrimAmp)%Result
          call PentCut(PrimAmps(iPrimAmp))
          call QuadCut(PrimAmps(iPrimAmp))
          call TripCut(PrimAmps(iPrimAmp))
          call DoubCut(PrimAmps(iPrimAmp))
          call SingCut(PrimAmps(iPrimAmp))
          call EvalMasterIntegrals(PrimAmps(iPrimAmp),Mu**2)

!         UV renormalization (only mass CT in the moment)
!           call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),Mu)
!           if(iPrimAmp.le.6) then
!               PrimAmps(iPrimAmp)%Result(-1) = PrimAmps(iPrimAmp)%Result(-1) -0d0*3d0/2d0*(11d0/3d0- 0d0*2d0/3d0*Nf/3d0)*BornAmps(iPrimAmp)%Result
!           endif
!          call OneLoopDiv(PrimAmps(iPrimAmp),Mu,rdiv(2),rdiv(1))
!          call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,EHat)
      enddo

      print *, PrimAmps(1)%Result(0)+PrimAmps(1)%Result(1)

!       NLO_Res_Pol(-2:1) = (0d0,0d0)
!       do jPrimAmp=1,6
!       do iPrimAmp=1,NumBornAmps
!           NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + 4d0*Col1L_ttbgg(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(PrimAmps(jPrimAmp)%Result(-2:1)) )
!       enddo
!       enddo
!       NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1) + NLO_Res_Pol(-2:1)

    enddo
    EvalCS_1L_gggggg = PrimAmps(1)%Result(0)+PrimAmps(1)%Result(1)

return
END FUNCTION








FUNCTION EvalCS_LO_bbbWWgg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes
use ModMisc
use ModParameters
implicit none
real(8) ::  EvalCS_LO_bbbWWgg,yRnd(1:VegasMxDim),VgsWgt
complex(8) :: LO_Res_Pol,LO_Res_Unpol
integer :: iHel,jHel,kHel,iPrimAmp,jPrimAmp
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,ISFac
real(8) :: MomExt(1:4,1:10)
logical :: applyPSCut
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac
real(8) :: pdf(-6:6,1:2)
integer :: NBin(1:NumMaxHisto),NHisto,nHel(1:2),NRndHel
include 'vegas_common.f'

!   note:   gauge-invariance fails

   TopDecays = 1! reset 101-->1 to ensure normal running

   EvalCS_LO_bbbWWgg = 0d0
   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if( EHat.le.2d0*m_Top ) then
      EvalCS_LO_bbbWWgg = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)


   call EvalPhaseSpace_2tobbbWW(EHat,yRnd(3:14),MomExt(1:4,1:10),PSWgt)
   call boost2Lab(eta1,eta2,8,MomExt(1:4,1:8))
   ISFac = 0.5d0**2
   ExtParticle(1)%Mom(1:4) = dcmplx(MomExt(1:4,3))
   ExtParticle(2)%Mom(1:4) = dcmplx(MomExt(1:4,4))
   ExtParticle(3)%Mom(1:4) =-dcmplx(MomExt(1:4,1))
   ExtParticle(4)%Mom(1:4) =-dcmplx(MomExt(1:4,2))

   call Kinematics_TTBAR(.false.,MomExt(1:4,1:4),MomExt(1:4,5:10),applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_LO_bbbWWgg = 0d0
      return
   endif

   call SetPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   RunFactor = RunAlphaS(NLOParam,MuRen)
   nHel(1:2) = getHelicity(yrnd(15))
   PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))

   call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,5:7))
   call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,8:10))

   ExtParticle(1)%Pol(1:4) = ExtParticle(1)%Pol(1:4) * dsqrt(2d0*Ga_Top(0)*m_Top) /cdsqrt( (MomExt(1:4,3).dot.MomExt(1:4,3))-m_Top**2+(0d0,1d0)*Ga_Top(0)*0.1d-4 )
   ExtParticle(2)%Pol(1:4) = ExtParticle(2)%Pol(1:4) * dsqrt(2d0*Ga_Top(0)*m_Top) /cdsqrt( (MomExt(1:4,4).dot.MomExt(1:4,4))-m_Top**2+(0d0,1d0)*Ga_Top(0)*0.1d-4 )


   LO_Res_Unpol = (0d0,0d0)
   do iHel=nHel(1),nHel(2)
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
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


!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * WidthExpansion
   EvalCS_LO_bbbWWgg = LO_Res_Unpol * PreFac

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_LO_bbbWWgg)
   enddo
   EvalCS_LO_bbbWWgg = EvalCS_LO_bbbWWgg/VgsWgt


return
END FUNCTION











END MODULE ModCrossSection_TTB


