MODULE ModCrossSection_TTBJ
use ModTopDecay
implicit none

integer,private,parameter :: NumMaxHisto=45


contains



FUNCTION EvalCS_1L_ttbggg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_new
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModIntDipoles_GGTTBGG
use ModIntDipoles_GGTTBQQB
implicit none
real(8) ::  EvalCS_1L_ttbggg,yRnd(1:VegasMxDim),VgsWgt,HOp(1:3),HOp2(1:3)
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1)
integer :: iHel,iPrimAmp,jPrimAmp,nHel(1:2)
integer :: NBin(1:NumMaxHisto),NHisto,ParityFlip,NRndHel
real(8) :: SpinAvg,ColorAvg,EHat,PSWgt,PSWgt2,PSWgt3,ISFac,RunFactor
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac
real(8) :: MomExt(1:4,1:5),MomDK(1:4,1:7),MomP(1:4,1:5)
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2),xE,AccPoles
real(8) :: MG_MOM(0:3,1:NumExtParticles)
real(8) :: MadGraph_tree,GG_TTBG,GG_TTB,GG_UUB
complex(8) :: FermionLoopPartAmp(25:35,-2:1)
logical :: applyPSCut
real(8) :: ThresholdCutOff = 1d0
include 'misc/global_import'
include 'vegas_common.f'



ParityFlip=1
IF( CORRECTION.EQ.1 ) ThresholdCutOff = 1.01d0


  EvalCS_1L_ttbggg = 0d0
  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)

  if( EHat.le.2d0*m_Top * ThresholdCutOff ) then
      EvalCS_1L_ttbggg = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)
! !DEC$ IF(_WriteTmpHisto .EQ.1)
!    if( it.ne.it_sav ) then
!       it_sav=it
! !       call WriteInstantHisto()
!       call WriteHisto(14,0d0,0d0,0d0)
!    endif
! !DEC$ ENDIF

   call EvalPhaseSpace_2to3(EHat,yRnd(3:7),MomExt(1:4,1:5),PSWgt)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))

   ISFac = MomCrossing(MomExt)
IF( TOPDECAYS.GE.1 ) THEN
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomDK(1:4,1:3),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:15),.false.,MomDK(1:4,4:6),PSWgt3)
   PSWgt = PSWgt * PSWgt2*PSWgt3
ENDIF

   call Kinematics_TTBARJET(0,MomExt,MomDK,applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_1L_ttbggg = 0d0
      PSCutCounter = PSCutCounter + 1
      return
   endif

   NRndHel=8
IF( TOPDECAYS.GE.1 ) THEN
   call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
   call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
   NRndHel=16
ENDIF


   call InitCurrCache()
   ISFac = MomCrossing(MomExt)
   call SetPropagators()
   call setPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   RunFactor = RunAlphaS(NLOParam,MuRen)

   LO_Res_Unpol             = (0d0,0d0)
   NLO_Res_Unpol(-2:1)      = (0d0,0d0)
   NLO_Res_Unpol_Ferm(-2:1) = (0d0,0d0)

IF( Correction.EQ.0 ) THEN
    nHel(1:2) = getHelicity(yrnd(NRndHel))
    PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))
    do iHel=nHel(1),nHel(2)
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()

      do iPrimAmp=1,6
          call EvalTree(BornAmps(iPrimAmp))
      enddo

      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ParityFlip*ColLO_ttbggg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result * dconjg(BornAmps(jPrimAmp)%Result)
!           LO_Res_Pol = LO_Res_Pol + ColLO_ttbggg(iPrimAmp,jPrimAmp) *                     &
!                  2d0*( dreal(BornAmps(iPrimAmp)%Result)*dreal(BornAmps(jPrimAmp)%Result) &
!                      + dimag(BornAmps(iPrimAmp)%Result)*dimag(BornAmps(jPrimAmp)%Result) )
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
    enddo

ELSEIF( Correction.EQ.1 ) THEN
    nHel(1:2) = getHelicity(yrnd(NRndHel))
    PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))
    do iHel=nHel(1),nHel(2)
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()

      do iPrimAmp=1,24
          call EvalTree(BornAmps(iPrimAmp))
      enddo
      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ParityFlip*ColLO_ttbggg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result * dconjg(BornAmps(jPrimAmp)%Result)
!           LO_Res_Pol = LO_Res_Pol + ColLO_ttbggg(iPrimAmp,jPrimAmp) *                    &
!                  2d0*( dreal(BornAmps(iPrimAmp)%Result)*dreal(BornAmps(jPrimAmp)%Result) &
!                      + dimag(BornAmps(iPrimAmp)%Result)*dimag(BornAmps(jPrimAmp)%Result) )
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol

!------------ 1 LOOP bosonic --------------
        do iPrimAmp=1,24
          call SetKirill(PrimAmps(iPrimAmp))
          call PentCut(PrimAmps(iPrimAmp))
          call QuadCut(PrimAmps(iPrimAmp))
          call TripCut(PrimAmps(iPrimAmp))
          call DoubCut(PrimAmps(iPrimAmp))
          call SingCut(PrimAmps(iPrimAmp))
          call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
! print *, "before",iPrimAmp,PrimAmps(iPrimAmp)%Result(-2:1)
!           call SetKirill(PrimAmps(iPrimAmp))
!           call PentCut_new(PrimAmps(:),iPrimAmp)
!           call QuadCut_new(PrimAmps(:),iPrimAmp)
!           call TripCut_new(PrimAmps(:),iPrimAmp)
!           call DoubCut_new(PrimAmps(:),iPrimAmp)
!           call SingCut_new(PrimAmps(:),iPrimAmp)
!           call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
! print *, "after ",iPrimAmp,PrimAmps(iPrimAmp)%Result(-2:1)
! pause

          call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen**2)
          PrimAmps(iPrimAmp)%Result(-2:1) = -PrimAmps(iPrimAmp)%Result(-2:1)    ! A_R -> A_L factor
          PrimAmps(iPrimAmp)%Result(-2:1) = -PrimAmps(iPrimAmp)%Result(-2:1)    ! correction to compensate in color factor
          PrimAmps(iPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(iPrimAmp)%Result(-2:1)
          call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,3,rdiv(2),rdiv(1))
!         call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
!         call WriteLatexOutput(iPrimAmp,PrimAmps(iPrimAmp),BornAmps(iPrimAmp))
!DEC$ IF (_QuadPrecImpr==1)
          AccPoles = CheckPoles(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv(1:2))
          if( AccPoles.gt.1d-4 ) then
!              print *, ""
!              print *, "QP re-run",AccPoles
!             coeff4_128(:,:) = qcmplx( coeff4(:,:) )
!             coeff5_128(:,:) = qcmplx( coeff5(:,:) )
              call PentCut_128(PrimAmps(iPrimAmp))
              call QuadCut_128(PrimAmps(iPrimAmp))
              call TripCut_128(PrimAmps(iPrimAmp))
              call DoubCut_128(PrimAmps(iPrimAmp))
              call SingCut_128(PrimAmps(iPrimAmp))
              call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
              call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen*2)
              PrimAmps(iPrimAmp)%Result(-2:1) = -PrimAmps(iPrimAmp)%Result(-2:1)    ! A_R -> A_L factor
              PrimAmps(iPrimAmp)%Result(-2:1) = -PrimAmps(iPrimAmp)%Result(-2:1)    ! correction to compensate in color factor
              PrimAmps(iPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(iPrimAmp)%Result(-2:1)
!              call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,rdiv(2),rdiv(1))
!              call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
              AccPoles = CheckPoles(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv(1:2))
              if( AccPoles.gt.5d-2 ) then
                  print *, "SKIP",AccPoles
                  call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
                  EvalCS_1L_ttbggg = 0d0
                  SkipCounter = SkipCounter + 1
                  return
              endif
          endif
!DEC$ ENDIF
      enddo

      NLO_Res_Pol(-2:1) = (0d0,0d0)
      do jPrimAmp=1,24
      do iPrimAmp=1,NumBornAmps
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) +  ParityFlip*Col1L_ttbggg(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(PrimAmps(jPrimAmp)%Result(-2:1)) )  *(-1d0)   ! sign correction to compensate in color factor
      enddo
      enddo
      NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1) + NLO_Res_Pol(-2:1)

! ------------ 1 LOOP fermionic loops  --------------
       do iPrimAmp=25,36
           call SetKirill(PrimAmps(iPrimAmp))
           call PentCut(PrimAmps(iPrimAmp))
           call QuadCut(PrimAmps(iPrimAmp))
           call TripCut(PrimAmps(iPrimAmp))
           call DoubCut(PrimAmps(iPrimAmp))
           call SingCut(PrimAmps(iPrimAmp))
           call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
           call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen**2)
           PrimAmps(iPrimAmp)%Result(-2:1) = -(0d0,1d0)*PrimAmps(iPrimAmp)%Result(-2:1)  ! minus is for closed fermion loop
!            call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),-rdiv,(/EHat/))
!            call WriteLatexOutput(iPrimAmp,PrimAmps(iPrimAmp),BornAmps(iPrimAmp))
       enddo
       FermionLoopPartAmp(25,-2:1) = PrimAmps(25)%Result(-2:1)
       FermionLoopPartAmp(26,-2:1) = PrimAmps(26)%Result(-2:1)
       FermionLoopPartAmp(27,-2:1) = PrimAmps(27)%Result(-2:1)
       FermionLoopPartAmp(28,-2:1) = PrimAmps(28)%Result(-2:1)
       FermionLoopPartAmp(29,-2:1) = PrimAmps(29)%Result(-2:1)
       FermionLoopPartAmp(30,-2:1) = PrimAmps(30)%Result(-2:1)

       FermionLoopPartAmp(31,-2:1) = PrimAmps(35)%Result(-2:1)  &
                                   + PrimAmps(36)%Result(-2:1)

       FermionLoopPartAmp(32,-2:1) = PrimAmps(33)%Result(-2:1)  &
                                   + PrimAmps(34)%Result(-2:1)

       FermionLoopPartAmp(33,-2:1) = PrimAmps(31)%Result(-2:1)  &
                                   + PrimAmps(32)%Result(-2:1)

       FermionLoopPartAmp(34,-2:1) =-PrimAmps(25)%Result(-2:1)  &
                                   - PrimAmps(28)%Result(-2:1)  &
                                   - PrimAmps(29)%Result(-2:1)  &
                                   - PrimAmps(31)%Result(-2:1)  &
                                   - PrimAmps(34)%Result(-2:1)  &
                                   - PrimAmps(35)%Result(-2:1)

       FermionLoopPartAmp(35,-2:1) =-PrimAmps(26)%Result(-2:1)  &
                                   - PrimAmps(27)%Result(-2:1)  &
                                   - PrimAmps(30)%Result(-2:1)  &
                                   - PrimAmps(32)%Result(-2:1)  &
                                   - PrimAmps(33)%Result(-2:1)  &
                                   - PrimAmps(36)%Result(-2:1)

       NLO_Res_Pol(-2:1) = (0d0,0d0)
       do jPrimAmp=25,35
       do iPrimAmp=1,NumBornAmps
           NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + ParityFlip*Col1L_ttbggg(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(FermionLoopPartAmp(jPrimAmp,-2:1)) )
       enddo
       enddo
       NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) + Nf_light*NLO_Res_Pol(-2:1)


       do iPrimAmp=37,48
           call SetKirill(PrimAmps(iPrimAmp))
           call PentCut(PrimAmps(iPrimAmp))
           call QuadCut(PrimAmps(iPrimAmp))
           call TripCut(PrimAmps(iPrimAmp))
           call DoubCut(PrimAmps(iPrimAmp))
           call SingCut(PrimAmps(iPrimAmp))
           call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
           call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen**2)
           PrimAmps(iPrimAmp)%Result(-2:1) = -(0d0,1d0)*PrimAmps(iPrimAmp)%Result(-2:1)  ! minus is for closed fermion loop
!            call WriteLatexOutput(iPrimAmp,PrimAmps(iPrimAmp),BornAmps(iPrimAmp))
!            call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,rdiv(2),rdiv(1))
!            call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),-rdiv,(/EHat/))
       enddo
       FermionLoopPartAmp(25,-2:1) = PrimAmps(37)%Result(-2:1)
       FermionLoopPartAmp(26,-2:1) = PrimAmps(38)%Result(-2:1)
       FermionLoopPartAmp(27,-2:1) = PrimAmps(39)%Result(-2:1)
       FermionLoopPartAmp(28,-2:1) = PrimAmps(40)%Result(-2:1)
       FermionLoopPartAmp(29,-2:1) = PrimAmps(41)%Result(-2:1)
       FermionLoopPartAmp(30,-2:1) = PrimAmps(42)%Result(-2:1)

       FermionLoopPartAmp(31,-2:1) = PrimAmps(47)%Result(-2:1)  &
                                   + PrimAmps(48)%Result(-2:1)

       FermionLoopPartAmp(32,-2:1) = PrimAmps(45)%Result(-2:1)  &
                                   + PrimAmps(46)%Result(-2:1)

       FermionLoopPartAmp(33,-2:1) = PrimAmps(43)%Result(-2:1)  &
                                   + PrimAmps(44)%Result(-2:1)

       FermionLoopPartAmp(34,-2:1) =-PrimAmps(37)%Result(-2:1)  &
                                   - PrimAmps(40)%Result(-2:1)  &
                                   - PrimAmps(41)%Result(-2:1)  &
                                   - PrimAmps(43)%Result(-2:1)  &
                                   - PrimAmps(46)%Result(-2:1)  &
                                   - PrimAmps(47)%Result(-2:1)

       FermionLoopPartAmp(35,-2:1) =-PrimAmps(38)%Result(-2:1)  &
                                   - PrimAmps(39)%Result(-2:1)  &
                                   - PrimAmps(42)%Result(-2:1)  &
                                   - PrimAmps(44)%Result(-2:1)  &
                                   - PrimAmps(45)%Result(-2:1)  &
                                   - PrimAmps(48)%Result(-2:1)

       NLO_Res_Pol(-2:1) = (0d0,0d0)
       do jPrimAmp=25,35
       do iPrimAmp=1,NumBornAmps
           NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + ParityFlip*Col1L_ttbggg(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(FermionLoopPartAmp(jPrimAmp,-2:1)) )
       enddo
       enddo
       NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) + Nf_heavy*NLO_Res_Pol(-2:1)
   enddo! helicity loop
ENDIF

IF( CORRECTION.EQ.0 ) THEN
!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**3 * WidthExpansion
   EvalCS_1L_ttbggg = LO_Res_Unpol * PreFac
ELSEIF( CORRECTION.EQ.1 ) THEN
!  overall normalization: (4*Pi)^eps/Gamma(1-eps)
!  CT contributions                                   ! beta fct            top WFRC
   NLO_Res_UnPol(-1) = NLO_Res_UnPol(-1) + (3d0/2d0*(-11d0 +2d0/3d0*Nf_light) -4d0   )*LO_Res_Unpol
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-3d0*4d0/3d0)*2d0*dlog(MuRen/m_top)*LO_Res_Unpol  ! finite log(mu) contrib. from  top WFRC

   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-5d0/2d0*8d0/3d0 )*LO_Res_Unpol   ! finite contribution from top WFRC's
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + 3d0/2d0*LO_Res_Unpol  ! shift alpha_s^DR --> alpha_s^MSbar

!    NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) - 3d0/2d0*LO_Res_Unpol     ! scheme shift to t'HV for DUW comparison, commented out because already in int.dips
!    print *, "this is for comparison with DUW, replace by alpha_s shift only"

!  factor out (Mu2/mTop**2)^eps
!    NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + NLO_Res_UnPol(-1)*2d0*dlog(m_top/MuRen) + NLO_Res_UnPol(-2)*dlog(m_top/MuRen)**2
!    NLO_Res_UnPol(-1) = NLO_Res_UnPol(-1) + NLO_Res_UnPol(-2)*2d0*dlog(m_top/MuRen)
!    NLO_Res_UnPol_Ferm(0) = NLO_Res_UnPol_Ferm(0) + NLO_Res_UnPol_Ferm(-1)*2d0*dlog(m_top/MuRen)

!  normalization
   LO_Res_Unpol = LO_Res_Unpol                         * ISFac * (alpha_s4Pi*RunFactor)**3
   NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1)           * ISFac * (alpha_s4Pi*RunFactor)**3 * alpha_sOver2Pi*RunFactor
   NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) * ISFac * (alpha_s4Pi*RunFactor)**3 * alpha_sOver2Pi*RunFactor

   EvalCS_1L_ttbggg = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac

!   HOp(1:3)=0d0
!   xE = 0.123d0
!   MomP(1:4,1) = MomExt(1:4,4)
!   MomP(1:4,2) = MomExt(1:4,5)
!   MomP(1:4,3) =-MomExt(1:4,1)
!   MomP(1:4,4) =-MomExt(1:4,2)
!   MomP(1:4,5) = MomExt(1:4,3)
!   call EvalIntDipoles_GGTTBGG(MomP(1:4,1:5),xE,HOp(1:3))
!   HOp(1:3) = HOp(1:3)*RunFactor**4 * 2d0
!    print *, "HOp(1)=", HOp(1)
!    print *, "HOp(2)=", HOp(2)
!    print *, "HOp(3)=", HOp(3)
!    print *, "NLO(-2)=",( NLO_Res_UnPol(-2) + NLO_Res_UnPol_Ferm(-2) )
!    print *, "NLO(-1)=",( NLO_Res_UnPol(-1) + NLO_Res_UnPol_Ferm(-1) )


!    HOp2(1:3)=0d0
!    xE = 0.123d0
!    MomP(1:4,1) = MomExt(1:4,4)
!    MomP(1:4,2) = MomExt(1:4,5)
!    MomP(1:4,3) =-MomExt(1:4,1)
!    MomP(1:4,4) =-MomExt(1:4,2)
!    MomP(1:4,5) = MomExt(1:4,3)
!    call EvalIntDipoles_GGTTBQQB(MomP(1:4,1:5),xE,HOp2(1:3))
!    HOp2(1:3) = HOp2(1:3)*RunFactor**4 * 2d0


!     print *, "HOp(1)=", HOp2(1)
!     print *, "HOp(2)=", HOp2(2)
!     print *, "HOp(3)=", HOp2(3)
!     print *, "NLO(-2)=",( NLO_Res_UnPol(-2) + NLO_Res_UnPol_Ferm(-2) )
!     print *, "NLO(-1)=",( NLO_Res_UnPol(-1) + NLO_Res_UnPol_Ferm(-1) )
!      NLO_Res_UnPol(0) = NLO_Res_UnPol(0) + HOp(1) + HOp2(1)


ELSEIF( CORRECTION.EQ.3 ) THEN

   MomP(1:4,1) = MomExt(1:4,4)
   MomP(1:4,2) = MomExt(1:4,5)
   MomP(1:4,3) =-MomExt(1:4,1)
   MomP(1:4,4) =-MomExt(1:4,2)
   MomP(1:4,5) = MomExt(1:4,3)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   IF( TOPDECAYS.GE.1 ) THEN
       xE = yRnd(16)
   ELSEIF( TOPDECAYS.EQ.0 ) THEN
       xE = yRnd(8)
   ENDIF
   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
   IF( PROCESS.EQ.9 ) THEN
          call EvalIntDipoles_GGTTBGG(MomP(1:4,1:5),MomDK(1:4,1:6),xE,HOp(1:3))
   ELSEIF( PROCESS.EQ.11 ) THEN
          call EvalIntDipoles_GGTTBQQB(MomP(1:4,1:5),MomDK(1:4,1:6),xE,HOp(1:3))
   ENDIF

   HOp(1:3) = HOp(1:3)*RunFactor**4 * PreFac
   EvalCS_1L_ttbggg = HOp(1)    * pdf(0,1)  * pdf(0,2)   &
                    + HOp(2)/xE * pdf_z(0,1)* pdf(0,2)   &
                    + HOp(3)/xE * pdf(0,1)  * pdf_z(0,2)
ENDIF


   if( IsNan(EvalCS_1L_ttbggg) .or. dabs(EvalCS_1L_ttbggg).eq.1d0/0d0  ) then
        print *, "NAN:",EvalCS_1L_ttbggg
        call printYrnd(yrnd(:))
        print *, PSWgt , VgsWgt , PDFFac, sHatJacobi
        print *, eta1,eta2,MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,1)
        print *, MomExt(1:4,2)
        print *, MomExt(1:4,3)
        print *, MomExt(1:4,4)
        print *, MomExt(1:4,5)
!         stop
        EvalCS_1L_ttbggg = 0d0
        return
   endif

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_ttbggg)
   enddo
   EvalCounter = EvalCounter + 1
   EvalCS_1L_ttbggg = EvalCS_1L_ttbggg/VgsWgt


!     gg->ttbg
!       ColorAvg=1d0/64d0
!       SpinAvg =0.25d0
!       MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!       MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!       MG_MOM(0:3,3) = MomExt(1:4,5)*100d0
!       MG_MOM(0:3,4) = MomExt(1:4,4)*100d0
!       MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!       call coupsm(0)
!       call SGG_TTBG(MG_MOM,MadGraph_tree)
!       print *, "MG  ratio: ", MadGraph_tree/dble(LO_Res_Unpol)
!       pause

!       print *, "./ttbarjets Collider=1 TopDK=0 Correction=1  Process=5 ObsSet=11 NLOParam=2 MTop=1.74 PDFSet=1"
!       print *, "remember: set DUW alpha_s,mt=1.74, t'H-V shift, NLOParam=2"
!       print *, ""
!       print *, "Unpolarized LO result:"
!       print *, "My tree:          ", LO_Res_Unpol/(4d0*DblPi*alpha_s*RunFactor)**3
!       print *, "DUW:             ",0.6566843362709775d-3 *(100d0)**2
!       print *, "ratio: ", 0.6566843362709775d-3 *(100d0)**2/(LO_Res_Unpol/(4d0*DblPi*alpha_s*RunFactor)**3)
!       print *, ""
!       print *, "Unpolarized NLO result:",NLO_Res_UnPol(-2:1)
!       print *, "c_{-2)_bos:",NLO_Res_UnPol(-2)/LO_Res_UnPol
!       print *, "c_{-2)_fer:",NLO_Res_UnPol_Ferm(-2)/LO_Res_UnPol
!       print *, "c_{-2)_DUW:",(-0.1540118421074573d0,0d0)
!       print *, "ratio to DUW:",(NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2))/LO_Res_UnPol / (-0.1540118421074573d0)
!       print *, ""
!       print *, "c_{-1)_bos:",NLO_Res_UnPol(-1)/LO_Res_UnPol
!       print *, "c_{-1)_fer:",NLO_Res_UnPol_Ferm(-1)/LO_Res_UnPol
!       print *, "c_{-1)_DUW:",(0.0731096894943437d0,0d0)
!       print *, "ratio to DUW:",(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/LO_Res_UnPol /(0.0731096894943437d0)
!       print *, ""
!       NLO_Res_UnPol(0) = NLO_Res_UnPol(0) - NLO_Res_UnPol(-2)*DblPi**2/6d0  ! different gamma function
!       print *, "c_{0)_bos:",(NLO_Res_UnPol(0)+NLO_Res_UnPol(1))/LO_Res_UnPol
!       print *, "c_{0)_fer:",(NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1))/LO_Res_UnPol
!       print *, "c_{0)_DUW:",(0.5295183452413002d0,0d0)
!       print *, "ratio to DUW:",(NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1)+NLO_Res_UnPol(0)+NLO_Res_UnPol(1))/LO_Res_UnPol/(0.5295183452413002d0)
!       stop


return
END FUNCTION









FUNCTION EvalCS_1L_ttbqqbg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use modIntDipoles_QQBTTBGG
use ModIntDipoles_QGTTBQG
use ModIntDipoles_QBGTTBQBG
use ModSixFermionIntDip
implicit none
real(8) ::  EvalCS_1L_ttbqqbg,yRnd(1:VegasMxDim),VgsWgt,xE,IOp(-2:0),HOp(1:3),HOp2(1:3)
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1)
complex(8) :: BosonicPartAmp(1:4,-2:1),FermionPartAmp(1:4,-2:1),PrimAmpsHel(1:16,1:16)
integer :: iHel,iPrimAmp,jPrimAmp,NRndHel,nHel(1:2)
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,ISFac
real(8) :: MomExt(1:4,1:NumExtParticles),MomDK(1:4,1:7),MomP(1:4,1:5),MomTmp(1:4)
logical :: applyPSCut
real(8) :: MG_MOM(0:3,1:5)
real(8) :: MadGraph_tree,UUB_TTB
real(8),parameter :: Nc=3d0
real(8) :: tau,eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,PDFFac_a,PDFFac_b,AccPoles
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2)
real(8) :: ThresholdCutOff = 1.0d0
integer :: NBin(1:NumMaxHisto),NHisto,Bin_PT,Bin_M,ParityFlip,npdf
include 'misc/global_import'
include 'vegas_common.f'



  EvalCS_1L_ttbqqbg = 0d0
  ParityFlip=1
  IF( CORRECTION.EQ.1 ) ThresholdCutOff = 1.01d0

  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top * ThresholdCutOff ) then
      EvalCS_1L_ttbqqbg = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to3(EHat,yRnd(3:7),MomExt(1:4,1:5),PSWgt)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))

IF( TOPDECAYS.GE.1 ) THEN
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomDK(1:4,1:3),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:15),.false.,MomDK(1:4,4:6),PSWgt3)
   PSWgt = PSWgt * PSWgt2*PSWgt3
ENDIF
!    call Kinematics(3,MomExt,MomDK,applyPSCut,NBin)
   call Kinematics_TTBARJET(0,MomExt,MomDK,applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_1L_ttbqqbg = 0d0
      PSCutCounter = PSCutCounter + 1
      return
   endif

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
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   RunFactor = RunAlphaS(NLOParam,MuRen)

   LO_Res_Unpol             = (0d0,0d0)
   NLO_Res_Unpol(-2:1)      = (0d0,0d0)
   NLO_Res_Unpol_Ferm(-2:1) = (0d0,0d0)

IF( TOPDECAYS.GE.1 ) THEN
   NRndHel=16
ELSE
   NRndHel=8
ENDIF
   nHel(1:2) = getHelicity(yrnd(NRndHel))
   PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))

DO NPDF=1,2
   if(npdf.eq.1) then
        PDFFac = PDFFac_a
   elseif(npdf.eq.2) then
        PDFFac = PDFFac_b
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
   endif

   ISFac = MomCrossing(MomExt)
IF( TOPDECAYS.GE.1 ) THEN
   call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
   call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
ENDIF
   call InitCurrCache()
   call SetPropagators()

   do iHel=nHel(1),nHel(2)
IF( CORRECTION.EQ.0 ) THEN
!------------ LO --------------
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo
      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqbg(iPrimAmp,jPrimAmp) * ParityFlip*BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol *PDFFac

ELSEIF( CORRECTION.EQ.1 ) THEN
!------------ 1 LOOP bosonic --------------
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      do iPrimAmp=1,24
          call EvalTree(BornAmps(iPrimAmp))
      enddo

      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqbg(iPrimAmp,jPrimAmp) * ParityFlip*BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol *PDFFac
      do iPrimAmp=1,16
           call SetKirill(PrimAmps(iPrimAmp))
           call PentCut(PrimAmps(iPrimAmp))
           call QuadCut(PrimAmps(iPrimAmp))
           call TripCut(PrimAmps(iPrimAmp))
           call DoubCut(PrimAmps(iPrimAmp))
           call SingCut(PrimAmps(iPrimAmp))
           call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
           call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen**2)
           PrimAmps(iPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(iPrimAmp)%Result(-2:1)
           call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,3,rdiv(2),rdiv(1))
!            call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
!            call WriteLatexOutput(iPrimAmp,PrimAmps(iPrimAmp),BornAmps(iPrimAmp))
!DEC$ IF (_QuadPrecImpr==1)
          AccPoles = CheckPoles(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv(1:2))
          if( AccPoles.gt.1d-4 ) then
!              print *, ""
!               print *, "QP re-run",AccPoles
!             coeff4_128(:,:) = qcmplx( coeff4(:,:) )
!             coeff5_128(:,:) = qcmplx( coeff5(:,:) )
              call PentCut_128(PrimAmps(iPrimAmp))
              call QuadCut_128(PrimAmps(iPrimAmp))
              call TripCut_128(PrimAmps(iPrimAmp))
              call DoubCut_128(PrimAmps(iPrimAmp))
              call SingCut_128(PrimAmps(iPrimAmp))
              call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
              call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen*2)
              PrimAmps(iPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(iPrimAmp)%Result(-2:1)
!              call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,rdiv(2),rdiv(1))
!              call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
              AccPoles = CheckPoles(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv(1:2))
          if( AccPoles.gt.5d-2 ) then
            print *, "SKIP",AccPoles
            call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
            EvalCS_1L_ttbqqbg = 0d0
            SkipCounter = SkipCounter + 1
            return
          endif
          endif
!DEC$ ENDIF
      enddo



!           iPrimAmp=7
!           ExtParticle(3)%Helicity=+1
!           ExtParticle(4)%Helicity=-1
!           call SetPolarizations()
!           call SetKirill(PrimAmps(iPrimAmp))
!           call PentCut(PrimAmps(iPrimAmp))
!           call QuadCut(PrimAmps(iPrimAmp))
!           call TripCut(PrimAmps(iPrimAmp))
!           call DoubCut(PrimAmps(iPrimAmp))
!           call SingCut(PrimAmps(iPrimAmp))
!           call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
!           call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen**2)
!           PrimAmps(iPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(iPrimAmp)%Result(-2:1)
!           print *, PrimAmps(iPrimAmp)%Result(-2)
!
!           iPrimAmp=4
!           ExtParticle(3)%Helicity=+1
!           ExtParticle(4)%Helicity=-1
!           MomTmp(1:4) = ExtParticle(3)%Mom(1:4)
!           ExtParticle(3)%Mom(1:4) = ExtParticle(4)%Mom(1:4)
!           ExtParticle(4)%Mom(1:4) = MomTmp(1:4)
!           call SetPropagators()
!           call SetPolarizations()
!
!           call SetKirill(PrimAmps(iPrimAmp))
!           call PentCut(PrimAmps(iPrimAmp))
!           call QuadCut(PrimAmps(iPrimAmp))
!           call TripCut(PrimAmps(iPrimAmp))
!           call DoubCut(PrimAmps(iPrimAmp))
!           call SingCut(PrimAmps(iPrimAmp))
!           call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
!           call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen**2)
!           PrimAmps(iPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(iPrimAmp)%Result(-2:1)
!           print *, PrimAmps(iPrimAmp)%Result(-2)
!           stop

!         print *, ExtParticle(1)%Helicity
!         print *, ExtParticle(2)%Helicity
!         print *, ExtParticle(3)%Helicity
!         print *, ExtParticle(4)%Helicity
!         print *, ExtParticle(5)%Helicity
!
!          write(*,"(A,PE20.10,PE20.10,PE20.10,PE20.10)")  "p1",dreal(ExtParticle(1)%Mom(1:4))
!          write(*,"(A,PE20.10,PE20.10,PE20.10,PE20.10)")  "p2",dreal(ExtParticle(2)%Mom(1:4))
!          write(*,"(A,PE20.10,PE20.10,PE20.10,PE20.10)")  "p3",dreal(ExtParticle(3)%Mom(1:4))
!          write(*,"(A,PE20.10,PE20.10,PE20.10,PE20.10)")  "p4",dreal(ExtParticle(4)%Mom(1:4))
!          write(*,"(A,PE20.10,PE20.10,PE20.10,PE20.10)")  "p5",dreal(ExtParticle(5)%Mom(1:4))
!
!          write(*,"(A,PE20.10,PE20.10,PE20.10,PE20.10,PE20.10,PE20.10,PE20.10,PE20.10)")  "e1",ExtParticle(1)%Pol(1:4)
!          write(*,"(A,PE20.10,PE20.10,PE20.10,PE20.10,PE20.10,PE20.10,PE20.10,PE20.10)")  "e2",ExtParticle(2)%Pol(1:4)
!          write(*,"(A,PE20.10,PE20.10,PE20.10,PE20.10,PE20.10,PE20.10,PE20.10,PE20.10)")  "e3",ExtParticle(3)%Pol(1:4)
!          write(*,"(A,PE20.10,PE20.10,PE20.10,PE20.10,PE20.10,PE20.10,PE20.10,PE20.10)")  "e4",ExtParticle(4)%Pol(1:4)
!          write(*,"(A,PE20.10,PE20.10,PE20.10,PE20.10,PE20.10,PE20.10,PE20.10,PE20.10)")  "e5",ExtParticle(5)%Pol(1:4)




!   set:   EHat=6d0,   MuRen=m_Top,   MuFac=m_Top, m_Top=175
!           YRnd(1:10)=0.1d0
!       LO_Res_Pol = BornAmps(PrimAmp1_12345)%Result   ! B7,1
!       BosonicPartAmp(1,-2:1) =    &    ! leading color
!                              + PrimAmps(PrimAmp1_12345)%Result(-2:1) * Nc

!       BosonicPartAmp(1,-2:1) =  BosonicPartAmp(1,-2:1) *0d0 +   &    ! sub-leading color
!                              ! (a) contribution
!                              - PrimAmps(PrimAmp1_12345)%Result(-2:1)/Nc   &
!                              + PrimAmps(PrimAmp1_15234)%Result(-2:1)/Nc   &
!                              + PrimAmps(PrimAmp1_15243)%Result(-2:1)/Nc   &
!                              + PrimAmps(PrimAmp1_12534)%Result(-2:1)/Nc   &
!                              + PrimAmps(PrimAmp1_12543)%Result(-2:1)/Nc   &
!                              + PrimAmps(PrimAmp1_12354)%Result(-2:1)/Nc   &
!                              - PrimAmps(PrimAmp1_12453)%Result(-2:1)/Nc   &
!                              - PrimAmps(PrimAmp1_12435)%Result(-2:1)/Nc   &
!                              ! (b) contribution
!                              + PrimAmps(PrimAmp3_15432)%Result(-2:1)/Nc   &
!                              ! (c) contribution
!                              - PrimAmps(PrimAmp4_12345)%Result(-2:1)/Nc
!
!       BosonicPartAmp(1,-2:1) = BosonicPartAmp(1,-2:1) /Nc
! B71, leading color   /Nc
!    -2.0000000000E+00    8.6257443404E-13
!     1.0461309737E+01   -3.1415926493E+00
!    -4.0680737845E+00    1.1348659129E+01

! B71, sub-leading color  /Nc
!     1.0000000000E+00   -7.0430729253E-13
!    -7.9474715420E+00    6.5044923245E+00
!     2.7941688916E+00   -1.7873699180E+01

!          write(*,"(A)") "B_7;1"
!          write(*,"(PE20.10,PE20.10)") BosonicPartAmp(1,-2)/LO_Res_Pol
!          write(*,"(PE20.10,PE20.10)") BosonicPartAmp(1,-1)/LO_Res_Pol
!          write(*,"(PE20.10,PE20.10)") (BosonicPartAmp(1,0)+BosonicPartAmp(1,1))/LO_Res_Pol
!          pause





!       LO_Res_Pol = BornAmps(PrimAmp1_15234)%Result   ! B7,2
!       BosonicPartAmp(2,-2:1) =  &   ! leading color
!                              ! (a) contribution
!                              + PrimAmps(PrimAmp1_12543)%Result(-2:1) &
!                              + PrimAmps(PrimAmp1_12453)%Result(-2:1) &
!                              + PrimAmps(PrimAmp1_12435)%Result(-2:1) &
!                              ! (c) contribution
!                              + PrimAmps(PrimAmp4_12345)%Result(-2:1) &
!                              + PrimAmps(PrimAmp4_12354)%Result(-2:1) &
!                              + PrimAmps(PrimAmp4_12534)%Result(-2:1) &
!                              + PrimAmps(PrimAmp4_15234)%Result(-2:1)
!
!       BosonicPartAmp(2,-2:1) = BosonicPartAmp(2,-2:1)  &   ! sub-leading color
!                              ! (a) contribution
!                              - PrimAmps(PrimAmp1_15234)%Result(-2:1) /Nc**2 &
!                              - PrimAmps(PrimAmp1_15243)%Result(-2:1) /Nc**2 &
!                              ! (b) contribution
!                              - PrimAmps(PrimAmp3_15432)%Result(-2:1) /Nc**2 &
!                              - PrimAmps(PrimAmp3_14532)%Result(-2:1) /Nc**2 &
!                              - PrimAmps(PrimAmp3_14352)%Result(-2:1) /Nc**2 &
!                              ! (c) contribution
!                              - PrimAmps(PrimAmp4_15234)%Result(-2:1)/Nc**2

!       BosonicPartAmp(2,-2:1) =       BosonicPartAmp(2,-2:1) *nc**2
! B_7;2  leading color
!    -2.0000000000E+00    5.7145785068E-13
!     9.7888518787E+00   -3.0980211637E+00
!    -4.2401138788E+00    1.2405890447E+01
! B_7;2  subleading color /Nc**2
!     1.0000000000E+00    1.4369248109E-14
!    -6.5010255505E+00    6.2862901119E+00
!     3.0346660424E+00   -1.5605472257E+01


!          write(*,"(A)") "B_7;2"
!          write(*,"(PE20.10,PE20.10)") BosonicPartAmp(2,-2)/LO_Res_Pol
!          write(*,"(PE20.10,PE20.10)") BosonicPartAmp(2,-1)/LO_Res_Pol
!          write(*,"(PE20.10,PE20.10)") (BosonicPartAmp(2,0)+BosonicPartAmp(2,1))/LO_Res_Pol
!          pause




!       LO_Res_Pol = BornAmps(PrimAmp1_12534)%Result   ! B7,3
!       BosonicPartAmp(3,-2:1) =  &    ! leading color
!                              ! (a) contribution
!                              + PrimAmps(PrimAmp1_12534)%Result(-2:1) * Nc
!
!       BosonicPartAmp(3,-2:1) = BosonicPartAmp(3,-2:1) +  &    ! sub-leading color
!                              ! (a) contribution
!                              - PrimAmps(PrimAmp1_12534)%Result(-2:1)/Nc  &
!                              + PrimAmps(PrimAmp1_15234)%Result(-2:1)/Nc  &
!                              + PrimAmps(PrimAmp1_15243)%Result(-2:1)/Nc  &
!                              + PrimAmps(PrimAmp1_12354)%Result(-2:1)/Nc  &
!                              - PrimAmps(PrimAmp1_12453)%Result(-2:1)/Nc  &
!                              + PrimAmps(PrimAmp1_12345)%Result(-2:1)/Nc  &
!                              + PrimAmps(PrimAmp1_12435)%Result(-2:1)/Nc  &
!                              - PrimAmps(PrimAmp1_12543)%Result(-2:1)/Nc  &
!                              ! (b) contribution
!                              + PrimAmps(PrimAmp3_14352)%Result(-2:1)/Nc  &
!                              ! (c) contribution
!                              - PrimAmps(PrimAmp4_12534)%Result(-2:1)/Nc

!       BosonicPartAmp(3,-2:1) =       BosonicPartAmp(3,-2:1) *nc
! B_7;3 leading color
!   -2.0000000000E+00   -1.0802418871E-13
!    1.2079989629E+01   -3.1415926001E+00
!   -4.7467636899E+00    1.7712577764E+01
! B_7;3 subleading color
!     1.0000000000E+00    2.6139445995E-12
!    -8.8206114163E+00    5.7630506569E+00
!     1.1110410001E+01   -2.7766413908E+01

!          write(*,"(A)") "B_7;3"
!          write(*,"(PE20.10,PE20.10)") BosonicPartAmp(3,-2)/LO_Res_Pol
!          write(*,"(PE20.10,PE20.10)") BosonicPartAmp(3,-1)/LO_Res_Pol
!          write(*,"(PE20.10,PE20.10)") (BosonicPartAmp(3,0)+BosonicPartAmp(3,1))/LO_Res_Pol
!          pause




!       LO_Res_Pol = BornAmps(PrimAmp1_12354)%Result   ! B7,4
!       BosonicPartAmp(4,-2:1) =  &   ! leading color
!                              ! (a) contribution
!                              - PrimAmps(PrimAmp1_15243)%Result(-2:1) &
!                              - PrimAmps(PrimAmp1_12543)%Result(-2:1) &
!                              - PrimAmps(PrimAmp1_12435)%Result(-2:1) &
!                              ! (b) contribution
!                              - PrimAmps(PrimAmp3_15432)%Result(-2:1) &
!                              - PrimAmps(PrimAmp3_14352)%Result(-2:1) &
!                              - PrimAmps(PrimAmp3_14325)%Result(-2:1) &
!                              - PrimAmps(PrimAmp3_14532)%Result(-2:1)
!
!       BosonicPartAmp(4,-2:1) = BosonicPartAmp(4,-2:1)   &   ! sub-leading color
!                              ! (a) contribution
!                              + PrimAmps(PrimAmp1_12453)%Result(-2:1) /Nc**2 &
!                              - PrimAmps(PrimAmp1_12354)%Result(-2:1) /Nc**2      &
!                              ! (b) contribution
!                              + PrimAmps(PrimAmp3_14532)%Result(-2:1) /Nc**2      &
! !                              (c) contribution
!                              + PrimAmps(PrimAmp4_12345)%Result(-2:1) /Nc**2   &
!                              + PrimAmps(PrimAmp4_12534)%Result(-2:1) /Nc**2   &
!                              + PrimAmps(PrimAmp4_15234)%Result(-2:1) /Nc**2

!       BosonicPartAmp(4,-2:1) = BosonicPartAmp(4,-2:1) * Nc**2
! B_7;4  leading color
!    -2.0000000000E+00   -1.4012219994E-12
!     4.3603698103E+01    3.5794091571E+01
!    -1.5327304453E+02   -1.4419326844E+02
! B_7;4 subleading color
!     1.0000000000E+00    8.2651394543E-14
!    -6.5010276769E+00    6.2862893559E+00
!     6.9060737787E+00   -1.3972184221E+01

!          write(*,"(A)") "B_7;4"
!          write(*,"(PE20.10,PE20.10)") BosonicPartAmp(4,-2)/LO_Res_Pol
!          write(*,"(PE20.10,PE20.10)") BosonicPartAmp(4,-1)/LO_Res_Pol
!          write(*,"(PE20.10,PE20.10)") (BosonicPartAmp(4,0)+BosonicPartAmp(4,1))/LO_Res_Pol
!          pause
!!! check end

!       print *, "uncrossed box"
!       print *,  PrimAmps(PrimAmp1_15234)%Result(-2)
!       print *,  PrimAmps(PrimAmp1_15234)%Result(-1)
!       print *,  PrimAmps(PrimAmp1_15234)%Result(0)
!       print *,  PrimAmps(PrimAmp1_15234)%Result(1)
!       print *, "crossed box"
!       print *,  PrimAmps(PrimAmp1_15243)%Result(-2)
!       print *,  PrimAmps(PrimAmp1_15243)%Result(-1)
!       print *,  PrimAmps(PrimAmp1_15243)%Result(0)
!       print *,  PrimAmps(PrimAmp1_15243)%Result(1)
!       print *, "hand crossed box"
!
!       MomTmp(1:4)=extparticle(3)%Mom(1:4)
!       extparticle(3)%Mom(1:4)=extparticle(4)%Mom(1:4)
!       extparticle(4)%Mom(1:4)=MomTmp(1:4)
!
! !       MomTmp(1)=extparticle(3)%Helicity
! !       extparticle(3)%Helicity=extparticle(4)%Helicity
! !       extparticle(4)%Helicity=MomTmp(1)
!
!       call SetPropagators()
!       call SetPolarizations()
!            call SetKirill(PrimAmps(PrimAmp1_15234))
!            call PentCut(PrimAmps(PrimAmp1_15234))
!            call QuadCut(PrimAmps(PrimAmp1_15234))
!            call TripCut(PrimAmps(PrimAmp1_15234))
!            call DoubCut(PrimAmps(PrimAmp1_15234))
!            call SingCut(PrimAmps(PrimAmp1_15234))
!            call EvalMasterIntegrals(PrimAmps(PrimAmp1_15234),MuRen**2)
!            call RenormalizeUV(PrimAmps(PrimAmp1_15234),BornAmps(PrimAmp1_15234),MuRen**2)
!            PrimAmps(PrimAmp1_15234)%Result(-2:1) = (0d0,1d0) * PrimAmps(PrimAmp1_15234)%Result(-2:1)
!
!       print *, - PrimAmps(PrimAmp1_15234)%Result(-2)
!       print *, - PrimAmps(PrimAmp1_15234)%Result(-1)
!       print *, - PrimAmps(PrimAmp1_15234)%Result(0)
!       print *, - PrimAmps(PrimAmp1_15234)%Result(1)
!         pause
!!! check crossed box

      BosonicPartAmp(1,-2:1) =    &    ! leading color
                             + PrimAmps(PrimAmp1_12345)%Result(-2:1) * Nc
      BosonicPartAmp(1,-2:1) =  BosonicPartAmp(1,-2:1)  +   &    ! sub-leading color
                             ! (a) contribution
                             - PrimAmps(PrimAmp1_12345)%Result(-2:1)/Nc   &
                             + PrimAmps(PrimAmp1_15234)%Result(-2:1)/Nc   &
                             + PrimAmps(PrimAmp1_15243)%Result(-2:1)/Nc   &
                             + PrimAmps(PrimAmp1_12534)%Result(-2:1)/Nc   &
                             + PrimAmps(PrimAmp1_12543)%Result(-2:1)/Nc   &
                             + PrimAmps(PrimAmp1_12354)%Result(-2:1)/Nc   &
                             - PrimAmps(PrimAmp1_12453)%Result(-2:1)/Nc   &
                             - PrimAmps(PrimAmp1_12435)%Result(-2:1)/Nc   &
                             ! (b) contribution
                             + PrimAmps(PrimAmp3_15432)%Result(-2:1)/Nc   &
                             ! (c) contribution
                             - PrimAmps(PrimAmp4_12345)%Result(-2:1)/Nc


      BosonicPartAmp(2,-2:1) =  &   ! leading color
                             ! (a) contribution
                             + PrimAmps(PrimAmp1_12543)%Result(-2:1) &
                             + PrimAmps(PrimAmp1_12453)%Result(-2:1) &
                             + PrimAmps(PrimAmp1_12435)%Result(-2:1) &
                             ! (c) contribution
                             + PrimAmps(PrimAmp4_12345)%Result(-2:1) &
                             + PrimAmps(PrimAmp4_12354)%Result(-2:1) &
                             + PrimAmps(PrimAmp4_12534)%Result(-2:1) &
                             + PrimAmps(PrimAmp4_15234)%Result(-2:1)
      BosonicPartAmp(2,-2:1) = BosonicPartAmp(2,-2:1)  &   ! sub-leading color
                             ! (a) contribution
                             - PrimAmps(PrimAmp1_15234)%Result(-2:1) /Nc**2 &
                             - PrimAmps(PrimAmp1_15243)%Result(-2:1) /Nc**2 &
                             ! (b) contribution
                             - PrimAmps(PrimAmp3_15432)%Result(-2:1) /Nc**2 &
                             - PrimAmps(PrimAmp3_14532)%Result(-2:1) /Nc**2 &
                             - PrimAmps(PrimAmp3_14352)%Result(-2:1) /Nc**2 &
                             ! (c) contribution
                             - PrimAmps(PrimAmp4_15234)%Result(-2:1)/Nc**2



      BosonicPartAmp(3,-2:1) =  &    ! leading color
                             ! (a) contribution
                             + PrimAmps(PrimAmp1_12534)%Result(-2:1) * Nc
      BosonicPartAmp(3,-2:1) = BosonicPartAmp(3,-2:1) +  &    ! sub-leading color
                             ! (a) contribution
                             - PrimAmps(PrimAmp1_12534)%Result(-2:1)/Nc  &
                             + PrimAmps(PrimAmp1_15234)%Result(-2:1)/Nc  &
                             + PrimAmps(PrimAmp1_15243)%Result(-2:1)/Nc  &
                             + PrimAmps(PrimAmp1_12354)%Result(-2:1)/Nc  &
                             - PrimAmps(PrimAmp1_12453)%Result(-2:1)/Nc  &
                             + PrimAmps(PrimAmp1_12345)%Result(-2:1)/Nc  &
                             + PrimAmps(PrimAmp1_12435)%Result(-2:1)/Nc  &
                             - PrimAmps(PrimAmp1_12543)%Result(-2:1)/Nc  &
                             ! (b) contribution
                             + PrimAmps(PrimAmp3_14352)%Result(-2:1)/Nc  &
                             ! (c) contribution
                             - PrimAmps(PrimAmp4_12534)%Result(-2:1)/Nc



      BosonicPartAmp(4,-2:1) =  &   ! leading color
                             ! (a) contribution
                             - PrimAmps(PrimAmp1_15243)%Result(-2:1) &
                             - PrimAmps(PrimAmp1_12543)%Result(-2:1) &
                             - PrimAmps(PrimAmp1_12435)%Result(-2:1) &
                             ! (b) contribution
                             - PrimAmps(PrimAmp3_15432)%Result(-2:1) &
                             - PrimAmps(PrimAmp3_14352)%Result(-2:1) &
                             - PrimAmps(PrimAmp3_14325)%Result(-2:1) &
                             - PrimAmps(PrimAmp3_14532)%Result(-2:1)
      BosonicPartAmp(4,-2:1) = BosonicPartAmp(4,-2:1)   &   ! sub-leading color
                             ! (a) contribution
                             + PrimAmps(PrimAmp1_12453)%Result(-2:1) /Nc**2 &
                             - PrimAmps(PrimAmp1_12354)%Result(-2:1) /Nc**2      &
                             ! (b) contribution
                             + PrimAmps(PrimAmp3_14532)%Result(-2:1) /Nc**2      &
!                              (c) contribution
                             + PrimAmps(PrimAmp4_12345)%Result(-2:1) /Nc**2   &
                             + PrimAmps(PrimAmp4_12534)%Result(-2:1) /Nc**2   &
                             + PrimAmps(PrimAmp4_15234)%Result(-2:1) /Nc**2


      NLO_Res_Pol(-2:1) = (0d0,0d0)
      do jPrimAmp=1,4
      do iPrimAmp=1,NumBornAmps
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + ParityFlip*Col1L_ttbqqbg(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(BosonicPartAmp(jPrimAmp,-2:1)) )
      enddo
      enddo
      NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1) + NLO_Res_Pol(-2:1)*PDFFac


! ------------ 1-loop fermionic --------------
      do iPrimAmp=17,24
          call SetKirill(PrimAmps(iPrimAmp))

          call PentCut(PrimAmps(iPrimAmp))
          call QuadCut(PrimAmps(iPrimAmp))
          call TripCut(PrimAmps(iPrimAmp))
          call DoubCut(PrimAmps(iPrimAmp))
          call SingCut(PrimAmps(iPrimAmp))
          call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)

          PrimAmps(iPrimAmp)%Result(-2:1) = -(0d0,1d0)*PrimAmps(iPrimAmp)%Result(-2:1) !minus if from closed fermion loop
!           call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,3,rdiv(2),rdiv(1))
!           call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv); pause
!           call WriteLatexOutput(iPrimAmp,PrimAmps(iPrimAmp),BornAmps(iPrimAmp))
      enddo

      FermionPartAmp(1,-2:1) =  (Nf_light*PrimAmps(PrimAmp2_12345)%Result(-2:1) + Nf_heavy*PrimAmps(PrimAmp2m_12345)%Result(-2:1))
      FermionPartAmp(2,-2:1) =  (Nf_light*PrimAmps(PrimAmp2_15234)%Result(-2:1) + Nf_heavy*PrimAmps(PrimAmp2m_15234)%Result(-2:1))/Nc
      FermionPartAmp(3,-2:1) =  (Nf_light*PrimAmps(PrimAmp2_12534)%Result(-2:1) + Nf_heavy*PrimAmps(PrimAmp2m_12534)%Result(-2:1))
      FermionPartAmp(4,-2:1) =  (Nf_light*PrimAmps(PrimAmp2_12354)%Result(-2:1) + Nf_heavy*PrimAmps(PrimAmp2m_12354)%Result(-2:1))/Nc

      NLO_Res_Pol(-2:1) = (0d0,0d0)
      do jPrimAmp=1,4
      do iPrimAmp=1,4
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + ParityFlip*Col1L_ttbqqbg(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(FermionPartAmp(jPrimAmp,-2:1)) )
      enddo
      enddo
      NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) + NLO_Res_Pol(-2:1)*PDFFac
ENDIF

ENDDO! helicity loop
ENDDO! loop over a<-->b pdfs

call swapMom(MomExt(1:4,1),MomExt(1:4,2))   ! swap back to original order, for ID below



IF( CORRECTION.EQ.0 ) THEN
!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**3 * WidthExpansion
   EvalCS_1L_ttbqqbg = LO_Res_Unpol * PreFac

ELSEIF( CORRECTION.EQ.1 ) THEN

!  overall normalization: (4*Pi)^eps/Gamma(1-eps)
!  CT contributions                           ! beta           !top WFRC
   NLO_Res_UnPol(-1) = NLO_Res_UnPol(-1) + (-3d0/2d0*11d0 - 4d0  )*LO_Res_Unpol
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-3d0*4d0/3d0)*2d0*dlog(MuRen/m_top)*LO_Res_Unpol  ! finite log(mu2) contrib. from  top WFRC

   NLO_Res_UnPol_Ferm(-1) = NLO_Res_UnPol_Ferm(-1) + (Nf_light+Nf_heavy-Nf_heavy/3d0)*LO_Res_Unpol
   NLO_Res_UnPol_Ferm( 0) = NLO_Res_UnPol_Ferm( 0) + (2d0/3d0*Nf_heavy)*2d0*dlog(MuRen/m_top)*LO_Res_Unpol  ! finite log(mu2) contrib. from heavy flavor in alpha_s ren.

   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-5d0/2d0*8d0/3d0 )*LO_Res_Unpol   ! finite contribution from top WFRC's
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + 3d0/2d0*LO_Res_Unpol  ! shift alpha_s^DR --> alpha_s^MSbar

!    print *, "this is for comparison with DUW, replace by alpha_s shift only:"
!    NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-11d0/6d0)*LO_Res_Unpol ! scheme shift to t'HV for DUW comparison

!  factor out (Mu2/mTop**2)^eps
!    NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + NLO_Res_UnPol(-1)*2d0*dlog(m_top/MuRen) + NLO_Res_UnPol(-2)*dlog(m_top/MuRen)**2
!    NLO_Res_UnPol(-1) = NLO_Res_UnPol(-1) + NLO_Res_UnPol(-2)*2d0*dlog(m_top/MuRen)
!    NLO_Res_UnPol_Ferm(0) = NLO_Res_UnPol_Ferm(0) + NLO_Res_UnPol_Ferm(-1)*2d0*dlog(m_top/MuRen)

!  normalization
   LO_Res_Unpol = LO_Res_Unpol                         * ISFac * (alpha_s4Pi*RunFactor)**3
   NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1)           * ISFac * (alpha_s4Pi*RunFactor)**3 * alpha_sOver2Pi*RunFactor
   NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) * ISFac * (alpha_s4Pi*RunFactor)**3 * alpha_sOver2Pi*RunFactor

   EvalCS_1L_ttbqqbg = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac

   IF( PROCESS.EQ.3 ) THEN
!         xE = 0.123d0
!         MomP(1:4,1) = MomExt(1:4,4)
!         MomP(1:4,2) = MomExt(1:4,5)
!         MomP(1:4,3) =-MomExt(1:4,1)
!         MomP(1:4,4) = MomExt(1:4,3)
!         MomP(1:4,5) =-MomExt(1:4,2)
!        call EvalIntDipoles_QGTTBQG(MomP(1:4,1:5),xE,HOp(1:3))
!         HOp(1:3) = HOp(1:3)*RunFactor**4
!         call sixquark_intdip(MomP(1:4,1:5),xE,+1,-1,HOp2(1:3))
!         HOp2(1:3) = HOp2(1:3)*RunFactor**4
    ELSEIF( PROCESS.EQ.4 ) THEN
!         xE = 0.123d0
!         MomP(1:4,1) = MomExt(1:4,4)
!         MomP(1:4,2) = MomExt(1:4,5)
!         MomP(1:4,3) = MomExt(1:4,3)   ! check this
!         MomP(1:4,4) =-MomExt(1:4,1)
!         MomP(1:4,5) =-MomExt(1:4,2)
!        call EvalIntDipoles_QBGTTBQBG(MomP(1:4,1:5),xE,HOp(1:3))
!         HOp(1:3) = HOp(1:3)*RunFactor**4
!         call sixquark_intdip(MomP(1:4,1:5),xE,+1,-1,HOp2(1:3))
!         HOp2(1:3) = HOp2(1:3)*RunFactor**4
    ELSEIF( PROCESS.EQ.6 ) THEN
!        xE = 0.123d0
!        MomP(1:4,1) = MomExt(1:4,4)
!        MomP(1:4,2) = MomExt(1:4,5)
!        MomP(1:4,3) =-MomExt(1:4,1)
!        MomP(1:4,4) =-MomExt(1:4,2)
!        MomP(1:4,5) = MomExt(1:4,3)
!        call EvalIntDipoles_QQBTTBGG(MomP(1:4,1:5),xE,HOp(1:3))
!        HOp(1:3) = HOp(1:3)*RunFactor**4
!        call sixquark_intdip(MomP(1:4,1:5),xE,+1,-1,HOp2(1:3))
!        HOp2(1:3) = HOp2(1:3)*RunFactor**4
    ENDIF

!   print *, "HOp(1)=", (HOp(1)+HOp2(1)*0d0)/LO_Res_Unpol/(alpha_sOver2Pi*RunFactor)
!   print *, "NLO(-2)=",( NLO_Res_UnPol(-2) + NLO_Res_UnPol_Ferm(-2) )/LO_Res_Unpol/(alpha_sOver2Pi*RunFactor)
!   print *, "NLO(-1)=",( NLO_Res_UnPol(-1) + NLO_Res_UnPol_Ferm(-1) )/LO_Res_Unpol/(alpha_sOver2Pi*RunFactor)
!   pause

!      NLO_Res_UnPol(0) = NLO_Res_UnPol(0) + HOp(1) + HOp2(1)



ELSEIF( CORRECTION.EQ.3 ) THEN


IF( PROCESS.EQ.10 ) THEN

   MomP(1:4,1) = MomExt(1:4,4)
   MomP(1:4,2) = MomExt(1:4,5)
   MomP(1:4,3) =-MomExt(1:4,1)
   MomP(1:4,4) =-MomExt(1:4,2)
   MomP(1:4,5) = MomExt(1:4,3)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   IF( TOPDECAYS.GE.1 ) THEN
      xE = yRnd(16)
   ELSEIF( TOPDECAYS.EQ.0 ) THEN
      xE = yRnd(8)
   ENDIF

   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
   call EvalIntDipoles_QQBTTBGG(MomP(1:4,1:5),MomDK(1:4,1:6),xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**4 * PreFac
   EvalCS_1L_ttbqqbg=HOp(1)    * (pdf(Up_,1)*pdf(AUp_,2)+pdf(Dn_,1)*pdf(ADn_,2)+pdf(Chm_,1)*pdf(AChm_,2)+pdf(Str_,1)*pdf(AStr_,2)+pdf(Bot_,1)*pdf(ABot_,2) ) &
                    +HOp(2)/xE * (pdf_z(Up_,1)*pdf(AUp_,2)+pdf_z(Dn_,1)*pdf(ADn_,2)+pdf_z(Chm_,1)*pdf(AChm_,2)+pdf_z(Str_,1)*pdf(AStr_,2)+pdf_z(Bot_,1)*pdf(ABot_,2) ) &
                    +HOp(3)/xE * (pdf(Up_,1)*pdf_z(AUp_,2)+pdf(Dn_,1)*pdf_z(ADn_,2)+pdf(Chm_,1)*pdf_z(AChm_,2)+pdf(Str_,1)*pdf_z(AStr_,2)+pdf(Bot_,1)*pdf_z(ABot_,2) )


   call swapMom(MomP(1:4,3),MomP(1:4,4))
   call EvalIntDipoles_QQBTTBGG(MomP(1:4,1:5),MomDK(1:4,1:6),xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**4 * PreFac
   EvalCS_1L_ttbqqbg= EvalCS_1L_ttbqqbg +      &
                    +HOp(1)    * (pdf(Up_,2)*pdf(AUp_,1)+pdf(Dn_,2)*pdf(ADn_,1)+pdf(Chm_,2)*pdf(AChm_,1)+pdf(Str_,2)*pdf(AStr_,1)+pdf(Bot_,2)*pdf(ABot_,1) ) &
                    +HOp(2)/xE * (pdf_z(Up_,2)*pdf(AUp_,1)+pdf_z(Dn_,2)*pdf(ADn_,1)+pdf_z(Chm_,2)*pdf(AChm_,1)+pdf_z(Str_,2)*pdf(AStr_,1)+pdf_z(Bot_,2)*pdf(ABot_,1) ) &
                    +HOp(3)/xE * (pdf(Up_,2)*pdf_z(AUp_,1)+pdf(Dn_,2)*pdf_z(ADn_,1)+pdf(Chm_,2)*pdf_z(AChm_,1)+pdf(Str_,2)*pdf_z(AStr_,1)+pdf(Bot_,2)*pdf_z(ABot_,1) )

ELSEIF( PROCESS.EQ.12 ) THEN

   MomP(1:4,1) = MomExt(1:4,4)
   MomP(1:4,2) = MomExt(1:4,5)
   MomP(1:4,3) =-MomExt(1:4,1)
   MomP(1:4,4) = MomExt(1:4,3)
   MomP(1:4,5) =-MomExt(1:4,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   IF( TOPDECAYS.GE.1 ) THEN
      xE = yRnd(16)
   ELSEIF( TOPDECAYS.EQ.0 ) THEN
      xE = yRnd(8)
   ENDIF

   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
   call EvalIntDipoles_QGTTBQG(MomP(1:4,1:5),MomDK(1:4,1:6),xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**4 * PreFac
   EvalCS_1L_ttbqqbg=HOp(1)    * (pdf(Up_,1)*pdf(0,2)+pdf(Dn_,1)*pdf(0,2)+pdf(Chm_,1)*pdf(0,2)+pdf(Str_,1)*pdf(0,2)+pdf(Bot_,1)*pdf(0,2) ) &
                    +HOp(2)/xE * (pdf_z(Up_,1)*pdf(0,2)+pdf_z(Dn_,1)*pdf(0,2)+pdf_z(Chm_,1)*pdf(0,2)+pdf_z(Str_,1)*pdf(0,2)+pdf_z(Bot_,1)*pdf(0,2) ) &
                    +HOp(3)/xE * (pdf(Up_,1)*pdf_z(0,2)+pdf(Dn_,1)*pdf_z(0,2)+pdf(Chm_,1)*pdf_z(0,2)+pdf(Str_,1)*pdf_z(0,2)+pdf(Bot_,1)*pdf_z(0,2) )


   call swapMom(MomP(1:4,3),MomP(1:4,5))
   call EvalIntDipoles_QGTTBQG(MomP(1:4,1:5),MomDK(1:4,1:6),xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**4 * PreFac
   EvalCS_1L_ttbqqbg=EvalCS_1L_ttbqqbg +    &
                    +HOp(1)    * (pdf(Up_,2)*pdf(0,1)+pdf(Dn_,2)*pdf(0,1)+pdf(Chm_,2)*pdf(0,1)+pdf(Str_,2)*pdf(0,1)+pdf(Bot_,2)*pdf(0,1) ) &
                    +HOp(2)/xE * (pdf_z(Up_,2)*pdf(0,1)+pdf_z(Dn_,2)*pdf(0,1)+pdf_z(Chm_,2)*pdf(0,1)+pdf_z(Str_,2)*pdf(0,1)+pdf_z(Bot_,2)*pdf(0,1) ) &
                    +HOp(3)/xE * (pdf(Up_,2)*pdf_z(0,1)+pdf(Dn_,2)*pdf_z(0,1)+pdf(Chm_,2)*pdf_z(0,1)+pdf(Str_,2)*pdf_z(0,1)+pdf(Bot_,2)*pdf_z(0,1) )


ELSEIF( PROCESS.EQ.13 ) THEN

   MomP(1:4,1) = MomExt(1:4,4)
   MomP(1:4,2) = MomExt(1:4,5)
   MomP(1:4,3) = MomExt(1:4,3)
   MomP(1:4,4) =-MomExt(1:4,1)
   MomP(1:4,5) =-MomExt(1:4,2)

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   IF( TOPDECAYS.GE.1 ) THEN
      xE = yRnd(16)
   ELSEIF( TOPDECAYS.EQ.0 ) THEN
      xE = yRnd(8)
   ENDIF

   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
   call EvalIntDipoles_QBGTTBQBG(MomP(1:4,1:5),MomDK(1:4,1:6),xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**4 * PreFac
   EvalCS_1L_ttbqqbg=HOp(1)    * (pdf(AUp_,1)*pdf(0,2)+pdf(ADn_,1)*pdf(0,2)+pdf(AChm_,1)*pdf(0,2)+pdf(AStr_,1)*pdf(0,2)+pdf(ABot_,1)*pdf(0,2) ) &
                    +HOp(2)/xE * (pdf_z(AUp_,1)*pdf(0,2)+pdf_z(ADn_,1)*pdf(0,2)+pdf_z(AChm_,1)*pdf(0,2)+pdf_z(AStr_,1)*pdf(0,2)+pdf_z(ABot_,1)*pdf(0,2) ) &
                    +HOp(3)/xE * (pdf(AUp_,1)*pdf_z(0,2)+pdf(ADn_,1)*pdf_z(0,2)+pdf(AChm_,1)*pdf_z(0,2)+pdf(AStr_,1)*pdf_z(0,2)+pdf(ABot_,1)*pdf_z(0,2) )

   call swapMom(MomP(1:4,4),MomP(1:4,5))
   call EvalIntDipoles_QBGTTBQBG(MomP(1:4,1:5),MomDK(1:4,1:6),xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**4 * PreFac
   EvalCS_1L_ttbqqbg=EvalCS_1L_ttbqqbg +          &
                    +HOp(1)    * (pdf(AUp_,2)*pdf(0,1)+pdf(ADn_,2)*pdf(0,1)+pdf(AChm_,2)*pdf(0,1)+pdf(AStr_,2)*pdf(0,1)+pdf(ABot_,2)*pdf(0,1) ) &
                    +HOp(2)/xE * (pdf_z(AUp_,2)*pdf(0,1)+pdf_z(ADn_,2)*pdf(0,1)+pdf_z(AChm_,2)*pdf(0,1)+pdf_z(AStr_,2)*pdf(0,1)+pdf_z(ABot_,2)*pdf(0,1) ) &
                    +HOp(3)/xE * (pdf(AUp_,2)*pdf_z(0,1)+pdf(ADn_,2)*pdf_z(0,1)+pdf(AChm_,2)*pdf_z(0,1)+pdf(AStr_,2)*pdf_z(0,1)+pdf(ABot_,2)*pdf_z(0,1) )




!
! !----------------------------------
!
!
!
!    MomP(1:4,1) = MomExt(1:4,4)
!    MomP(1:4,2) = MomExt(1:4,5)
!    MomP(1:4,3) =-MomExt(1:4,1)
!    MomP(1:4,4) = MomExt(1:4,3)
!    MomP(1:4,5) =-MomExt(1:4,2)
!    PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
!    IF( TOPDECAYS.GE.1 ) THEN
!       xE = yRnd(16)
!    ELSEIF( TOPDECAYS.EQ.0 ) THEN
!       xE = yRnd(8)
!    ENDIF
!
!    call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
!    call EvalIntDipoles_QGTTBQG(MomP(1:4,1:5),MomDK(1:4,1:6),xE,HOp(1:3))
!    HOp(1:3) = HOp(1:3)*RunFactor**4 * PreFac
! !    EvalCS_1L_ttbqqbg=HOp(1)    * (pdf(AUp_,1)*pdf(0,2)+pdf(ADn_,1)*pdf(0,2)+pdf(AChm_,1)*pdf(0,2)+pdf(AStr_,1)*pdf(0,2)+pdf(ABot_,1)*pdf(0,2) ) &
! !                     +HOp(2)/xE * (pdf_z(AUp_,1)*pdf(0,2)+pdf_z(ADn_,1)*pdf(0,2)+pdf_z(AChm_,1)*pdf(0,2)+pdf_z(AStr_,1)*pdf(0,2)+pdf_z(ABot_,1)*pdf(0,2) ) &
! !                     +HOp(3)/xE * (pdf(AUp_,1)*pdf_z(0,2)+pdf(ADn_,1)*pdf_z(0,2)+pdf(AChm_,1)*pdf_z(0,2)+pdf(AStr_,1)*pdf_z(0,2)+pdf(ABot_,1)*pdf_z(0,2) )
!
!    EvalCS_1L_ttbqqbg=HOp(1)    * (pdf(ABot_,1)*pdf(0,2) ) &
!                     +HOp(2)/xE * (pdf_z(ABot_,1)*pdf(0,2) ) &
!                     +HOp(3)/xE * (pdf(ABot_,1)*pdf_z(0,2) )
!
! !   fails for ABtom, works for Bot


ELSEIF( PROCESS.EQ.14 ) THEN

   MomP(1:4,1) = MomExt(1:4,4)
   MomP(1:4,2) = MomExt(1:4,5)
   MomP(1:4,3) =-MomExt(1:4,1)
   MomP(1:4,4) =-MomExt(1:4,2)
   MomP(1:4,5) = MomExt(1:4,3)

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   IF( TOPDECAYS.GE.1 ) THEN
      xE = yRnd(16)
   ELSEIF( TOPDECAYS.EQ.0 ) THEN
      xE = yRnd(8)
   ENDIF

   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

   call sixquark_intdip(MomP(1:4,1:5),MomDK(1:4,1:6),xE,+1,-1,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**4 * PreFac
   EvalCS_1L_ttbqqbg = HOp(1) *( pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
                               + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
                               + pdf(Bot_,1)*pdf(ABot_,2)  ) &
!                             &
                    + HOp(2)/xE *( pdf_z(Up_,1) *pdf(AUp_,2)  + pdf_z(Dn_,1) *pdf(ADn_,2)   &
                                 + pdf_z(Chm_,1)*pdf(AChm_,2) + pdf_z(Str_,1)*pdf(AStr_,2)  &
                                 + pdf_z(Bot_,1)*pdf(ABot_,2) ) &
!                             &
                    + HOp(3)/xE *( pdf(Up_,1) *pdf_z(AUp_,2)  + pdf(Dn_,1) *pdf_z(ADn_,2)   &
                                 + pdf(Chm_,1)*pdf_z(AChm_,2) + pdf(Str_,1)*pdf_z(AStr_,2)  &
                                 + pdf(Bot_,1)*pdf_z(ABot_,2) )


   call sixquark_intdip(MomP(1:4,1:5),MomDK(1:4,1:6),xE,-1,+1,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**4 * PreFac
   EvalCS_1L_ttbqqbg = EvalCS_1L_ttbqqbg +  &
                       HOp(1) *( pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
                               + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
                               + pdf(Bot_,2)*pdf(ABot_,1)  ) &
!                             &
                    + HOp(2)/xE *( pdf_z(Up_,2) *pdf(AUp_,1)  + pdf_z(Dn_,2) *pdf(ADn_,1)   &
                                 + pdf_z(Chm_,2)*pdf(AChm_,1) + pdf_z(Str_,2)*pdf(AStr_,1)  &
                                 + pdf_z(Bot_,2)*pdf(ABot_,1) ) &
!                             &
                    + HOp(3)/xE *( pdf(Up_,2) *pdf_z(AUp_,1)  + pdf(Dn_,2) *pdf_z(ADn_,1)   &
                                 + pdf(Chm_,2)*pdf_z(AChm_,1) + pdf(Str_,2)*pdf_z(AStr_,1)  &
                                 + pdf(Bot_,2)*pdf_z(ABot_,1)  )



ELSEIF( PROCESS.EQ.15 ) THEN

   MomP(1:4,1) = MomExt(1:4,4)
   MomP(1:4,2) = MomExt(1:4,5)
   MomP(1:4,3) =-MomExt(1:4,1)
   MomP(1:4,4) =-MomExt(1:4,2)
   MomP(1:4,5) = MomExt(1:4,3)

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   IF( TOPDECAYS.GE.1 ) THEN
      xE = yRnd(16)
   ELSEIF( TOPDECAYS.EQ.0 ) THEN
      xE = yRnd(8)
   ENDIF

   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

   call sixquark_intdip(MomP(1:4,1:5),MomDK(1:4,1:6),xE,+1,-2,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**4 * PreFac
   EvalCS_1L_ttbqqbg = HOp(1) *( pdf(Up_,1) * ( pdf(ADn_,2) +pdf(AChm_,2) +pdf(AStr_,2) +pdf(ABot_,2) )    &
                              +  pdf(Dn_,1) * ( pdf(AUp_,2) +pdf(AChm_,2) +pdf(AStr_,2) +pdf(ABot_,2) )    &
                              +  pdf(Chm_,1)* ( pdf(AUp_,2) +pdf(ADn_,2)  +pdf(AStr_,2) +pdf(ABot_,2) )    &
                              +  pdf(Str_,1)* ( pdf(AUp_,2) +pdf(ADn_,2)  +pdf(AChm_,2) +pdf(ABot_,2) )    &
                              +  pdf(Bot_,1)* ( pdf(AUp_,2) +pdf(ADn_,2)  +pdf(AChm_,2) +pdf(AStr_,2) )  ) &
!                             &
                    + HOp(2)/xE *( pdf_z(Up_,1) * ( pdf(ADn_,2) +pdf(AChm_,2) +pdf(AStr_,2) +pdf(ABot_,2) )    &
                                +  pdf_z(Dn_,1) * ( pdf(AUp_,2) +pdf(AChm_,2) +pdf(AStr_,2) +pdf(ABot_,2) )    &
                                +  pdf_z(Chm_,1)* ( pdf(AUp_,2) +pdf(ADn_,2)  +pdf(AStr_,2) +pdf(ABot_,2) )    &
                                +  pdf_z(Str_,1)* ( pdf(AUp_,2) +pdf(ADn_,2)  +pdf(AChm_,2) +pdf(ABot_,2) )    &
                                +  pdf_z(Bot_,1)* ( pdf(AUp_,2) +pdf(ADn_,2)  +pdf(AChm_,2) +pdf(AStr_,2) ) ) &
!                                 &
                    + HOp(3)/xE *( pdf(Up_,1) * ( pdf_z(ADn_,2) +pdf_z(AChm_,2) +pdf_z(AStr_,2) +pdf_z(ABot_,2) )    &
                                +  pdf(Dn_,1) * ( pdf_z(AUp_,2) +pdf_z(AChm_,2) +pdf_z(AStr_,2) +pdf_z(ABot_,2) )    &
                                +  pdf(Chm_,1)* ( pdf_z(AUp_,2) +pdf_z(ADn_,2)  +pdf_z(AStr_,2) +pdf_z(ABot_,2) )    &
                                +  pdf(Str_,1)* ( pdf_z(AUp_,2) +pdf_z(ADn_,2)  +pdf_z(AChm_,2) +pdf_z(ABot_,2) )    &
                                +  pdf(Bot_,1)* ( pdf_z(AUp_,2) +pdf_z(ADn_,2)  +pdf_z(AChm_,2) +pdf_z(AStr_,2) )  )


   call sixquark_intdip(MomP(1:4,1:5),MomDK(1:4,1:6),xE,-2,+1,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**4 * PreFac
   EvalCS_1L_ttbqqbg = EvalCS_1L_ttbqqbg +  &
                       HOp(1) *( pdf(Up_,2) * ( pdf(ADn_,1) +pdf(AChm_,1) +pdf(AStr_,1) +pdf(ABot_,1) )    &
                              +  pdf(Dn_,2) * ( pdf(AUp_,1) +pdf(AChm_,1) +pdf(AStr_,1) +pdf(ABot_,1) )    &
                              +  pdf(Chm_,2)* ( pdf(AUp_,1) +pdf(ADn_,1)  +pdf(AStr_,1) +pdf(ABot_,1) )    &
                              +  pdf(Str_,2)* ( pdf(AUp_,1) +pdf(ADn_,1)  +pdf(AChm_,1) +pdf(ABot_,1) )    &
                              +  pdf(Bot_,2)* ( pdf(AUp_,1) +pdf(ADn_,1)  +pdf(AChm_,1) +pdf(AStr_,1) )  ) &
!                             &
                    + HOp(2)/xE *( pdf_z(Up_,2) * ( pdf(ADn_,1) +pdf(AChm_,1) +pdf(AStr_,1) +pdf(ABot_,1) )    &
                                +  pdf_z(Dn_,2) * ( pdf(AUp_,1) +pdf(AChm_,1) +pdf(AStr_,1) +pdf(ABot_,1) )    &
                                +  pdf_z(Chm_,2)* ( pdf(AUp_,1) +pdf(ADn_,1)  +pdf(AStr_,1) +pdf(ABot_,1) )    &
                                +  pdf_z(Str_,2)* ( pdf(AUp_,1) +pdf(ADn_,1)  +pdf(AChm_,1) +pdf(ABot_,1) )    &
                                +  pdf_z(Bot_,2)* ( pdf(AUp_,1) +pdf(ADn_,1)  +pdf(AChm_,1) +pdf(AStr_,1) ) ) &
!                                 &
                    + HOp(3)/xE *( pdf(Up_,2) * ( pdf_z(ADn_,1) +pdf_z(AChm_,1) +pdf_z(AStr_,1) +pdf_z(ABot_,1) )    &
                                +  pdf(Dn_,2) * ( pdf_z(AUp_,1) +pdf_z(AChm_,1) +pdf_z(AStr_,1) +pdf_z(ABot_,1) )    &
                                +  pdf(Chm_,2)* ( pdf_z(AUp_,1) +pdf_z(ADn_,1)  +pdf_z(AStr_,1) +pdf_z(ABot_,1) )    &
                                +  pdf(Str_,2)* ( pdf_z(AUp_,1) +pdf_z(ADn_,1)  +pdf_z(AChm_,1) +pdf_z(ABot_,1) )    &
                                +  pdf(Bot_,2)* ( pdf_z(AUp_,1) +pdf_z(ADn_,1)  +pdf_z(AChm_,1) +pdf_z(AStr_,1) )  )


ELSEIF( PROCESS.EQ.16 ) THEN

   MomP(1:4,1) = MomExt(1:4,4)
   MomP(1:4,2) = MomExt(1:4,5)
   MomP(1:4,3) =-MomExt(1:4,1)
   MomP(1:4,4) =-MomExt(1:4,2)
   MomP(1:4,5) = MomExt(1:4,3)

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   IF( TOPDECAYS.GE.1 ) THEN
      xE = yRnd(16)
   ELSEIF( TOPDECAYS.EQ.0 ) THEN
      xE = yRnd(8)
   ENDIF

   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

   call sixquark_intdip(MomP(1:4,1:5),MomDK(1:4,1:6),xE,+1,+2,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**4 * PreFac
   EvalCS_1L_ttbqqbg = HOp(1) *( pdf(Up_,1) * ( pdf(Dn_,2) +pdf(Chm_,2) +pdf(Str_,2) +pdf(Bot_,2) )    &
                              +  pdf(Dn_,1) * ( pdf(Up_,2) +pdf(Chm_,2) +pdf(Str_,2) +pdf(Bot_,2) )    &
                              +  pdf(Chm_,1)* ( pdf(Up_,2) +pdf(Dn_,2)  +pdf(Str_,2) +pdf(Bot_,2) )    &
                              +  pdf(Str_,1)* ( pdf(Up_,2) +pdf(Dn_,2)  +pdf(Chm_,2) +pdf(Bot_,2) )    &
                              +  pdf(Bot_,1)* ( pdf(Up_,2) +pdf(Dn_,2)  +pdf(Chm_,2) +pdf(Str_,2) )  ) &
!                             &
                    + HOp(2)/xE *( pdf_z(Up_,1) * ( pdf(Dn_,2) +pdf(Chm_,2) +pdf(Str_,2) +pdf(Bot_,2) )    &
                                +  pdf_z(Dn_,1) * ( pdf(Up_,2) +pdf(Chm_,2) +pdf(Str_,2) +pdf(Bot_,2) )    &
                                +  pdf_z(Chm_,1)* ( pdf(Up_,2) +pdf(Dn_,2)  +pdf(Str_,2) +pdf(Bot_,2) )    &
                                +  pdf_z(Str_,1)* ( pdf(Up_,2) +pdf(Dn_,2)  +pdf(Chm_,2) +pdf(Bot_,2) )    &
                                +  pdf_z(Bot_,1)* ( pdf(Up_,2) +pdf(Dn_,2)  +pdf(Chm_,2) +pdf(Str_,2) ) ) &
!                                 &
                    + HOp(3)/xE *( pdf(Up_,1) * ( pdf_z(Dn_,2) +pdf_z(Chm_,2) +pdf_z(Str_,2) +pdf_z(Bot_,2) )    &
                                +  pdf(Dn_,1) * ( pdf_z(Up_,2) +pdf_z(Chm_,2) +pdf_z(Str_,2) +pdf_z(Bot_,2) )    &
                                +  pdf(Chm_,1)* ( pdf_z(Up_,2) +pdf_z(Dn_,2)  +pdf_z(Str_,2) +pdf_z(Bot_,2) )    &
                                +  pdf(Str_,1)* ( pdf_z(Up_,2) +pdf_z(Dn_,2)  +pdf_z(Chm_,2) +pdf_z(Bot_,2) )    &
                                +  pdf(Bot_,1)* ( pdf_z(Up_,2) +pdf_z(Dn_,2)  +pdf_z(Chm_,2) +pdf_z(Str_,2) )  )



!    call sixquark_intdip(MomP(1:4,1:5),MomDK(1:4,1:6),xE,+2,+1,HOp(1:3))
!    HOp(1:3) = HOp(1:3)*RunFactor**4 * PreFac
!    EvalCS_1L_ttbqqbg = EvalCS_1L_ttbqqbg +  &
!                        HOp(1) *( pdf(Up_,2) * ( pdf(Dn_,1) +pdf(Chm_,1) +pdf(Str_,1) +pdf(Bot_,1) )    &
!                               +  pdf(Dn_,2) * ( pdf(Up_,1) +pdf(Chm_,1) +pdf(Str_,1) +pdf(Bot_,1) )    &
!                               +  pdf(Chm_,2)* ( pdf(Up_,1) +pdf(Dn_,1)  +pdf(Str_,1) +pdf(Bot_,1) )    &
!                               +  pdf(Str_,2)* ( pdf(Up_,1) +pdf(Dn_,1)  +pdf(Chm_,1) +pdf(Bot_,1) )    &
!                               +  pdf(Bot_,2)* ( pdf(Up_,1) +pdf(Dn_,1)  +pdf(Chm_,1) +pdf(Str_,1) )  ) &
! !                             &
!                     + HOp(2)/xE *( pdf_z(Up_,2) * ( pdf(Dn_,1) +pdf(Chm_,1) +pdf(Str_,1) +pdf(Bot_,1) )    &
!                                 +  pdf_z(Dn_,2) * ( pdf(Up_,1) +pdf(Chm_,1) +pdf(Str_,1) +pdf(Bot_,1) )    &
!                                 +  pdf_z(Chm_,2)* ( pdf(Up_,1) +pdf(Dn_,1)  +pdf(Str_,1) +pdf(Bot_,1) )    &
!                                 +  pdf_z(Str_,2)* ( pdf(Up_,1) +pdf(Dn_,1)  +pdf(Chm_,1) +pdf(Bot_,1) )    &
!                                 +  pdf_z(Bot_,2)* ( pdf(Up_,1) +pdf(Dn_,1)  +pdf(Chm_,1) +pdf(Str_,1) ) ) &
! !                                 &
!                     + HOp(3)/xE *( pdf(Up_,2) * ( pdf_z(Dn_,1) +pdf_z(Chm_,1) +pdf_z(Str_,1) +pdf_z(Bot_,1) )    &
!                                 +  pdf(Dn_,2) * ( pdf_z(Up_,1) +pdf_z(Chm_,1) +pdf_z(Str_,1) +pdf_z(Bot_,1) )    &
!                                 +  pdf(Chm_,2)* ( pdf_z(Up_,1) +pdf_z(Dn_,1)  +pdf_z(Str_,1) +pdf_z(Bot_,1) )    &
!                                 +  pdf(Str_,2)* ( pdf_z(Up_,1) +pdf_z(Dn_,1)  +pdf_z(Chm_,1) +pdf_z(Bot_,1) )    &
!                                 +  pdf(Bot_,2)* ( pdf_z(Up_,1) +pdf_z(Dn_,1)  +pdf_z(Chm_,1) +pdf_z(Str_,1) )  )

ELSEIF( PROCESS.EQ.17 ) THEN

   MomP(1:4,1) = MomExt(1:4,4)
   MomP(1:4,2) = MomExt(1:4,5)
   MomP(1:4,3) =-MomExt(1:4,1)
   MomP(1:4,4) =-MomExt(1:4,2)
   MomP(1:4,5) = MomExt(1:4,3)

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   IF( TOPDECAYS.GE.1 ) THEN
      xE = yRnd(16)
   ELSEIF( TOPDECAYS.EQ.0 ) THEN
      xE = yRnd(8)
   ENDIF

   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

   call sixquark_intdip(MomP(1:4,1:5),MomDK(1:4,1:6),xE,-1,-2,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**4 * PreFac
   EvalCS_1L_ttbqqbg = HOp(1) *( pdf(AUp_,1) * ( pdf(ADn_,2) +pdf(AChm_,2) +pdf(AStr_,2) +pdf(ABot_,2) )    &
                              +  pdf(ADn_,1) * ( pdf(AUp_,2) +pdf(AChm_,2) +pdf(AStr_,2) +pdf(ABot_,2) )    &
                              +  pdf(AChm_,1)* ( pdf(AUp_,2) +pdf(ADn_,2)  +pdf(AStr_,2) +pdf(ABot_,2) )    &
                              +  pdf(AStr_,1)* ( pdf(AUp_,2) +pdf(ADn_,2)  +pdf(AChm_,2) +pdf(ABot_,2) )    &
                              +  pdf(ABot_,1)* ( pdf(AUp_,2) +pdf(ADn_,2)  +pdf(AChm_,2) +pdf(AStr_,2) )  ) &
!                             &
                    + HOp(2)/xE *( pdf_z(AUp_,1) * ( pdf(ADn_,2) +pdf(AChm_,2) +pdf(AStr_,2) +pdf(ABot_,2) )    &
                                +  pdf_z(ADn_,1) * ( pdf(AUp_,2) +pdf(AChm_,2) +pdf(AStr_,2) +pdf(ABot_,2) )    &
                                +  pdf_z(AChm_,1)* ( pdf(AUp_,2) +pdf(ADn_,2)  +pdf(AStr_,2) +pdf(ABot_,2) )    &
                                +  pdf_z(AStr_,1)* ( pdf(AUp_,2) +pdf(ADn_,2)  +pdf(AChm_,2) +pdf(ABot_,2) )    &
                                +  pdf_z(ABot_,1)* ( pdf(AUp_,2) +pdf(ADn_,2)  +pdf(AChm_,2) +pdf(AStr_,2) ) ) &
!                                 &
                    + HOp(3)/xE *( pdf(AUp_,1) * ( pdf_z(ADn_,2) +pdf_z(AChm_,2) +pdf_z(AStr_,2) +pdf_z(ABot_,2) )    &
                                +  pdf(ADn_,1) * ( pdf_z(AUp_,2) +pdf_z(AChm_,2) +pdf_z(AStr_,2) +pdf_z(ABot_,2) )    &
                                +  pdf(AChm_,1)* ( pdf_z(AUp_,2) +pdf_z(ADn_,2)  +pdf_z(AStr_,2) +pdf_z(ABot_,2) )    &
                                +  pdf(AStr_,1)* ( pdf_z(AUp_,2) +pdf_z(ADn_,2)  +pdf_z(AChm_,2) +pdf_z(ABot_,2) )    &
                                +  pdf(ABot_,1)* ( pdf_z(AUp_,2) +pdf_z(ADn_,2)  +pdf_z(AChm_,2) +pdf_z(AStr_,2) )  )



!    call sixquark_intdip(MomP(1:4,1:5),MomDK(1:4,1:6),xE,-2,-1,HOp(1:3))
!    HOp(1:3) = HOp(1:3)*RunFactor**4 * PreFac
!    EvalCS_1L_ttbqqbg = EvalCS_1L_ttbqqbg  + &
!                        HOp(1) *( pdf(AUp_,2) * ( pdf(ADn_,1) +pdf(AChm_,1) +pdf(AStr_,1) +pdf(ABot_,1) )    &
!                               +  pdf(ADn_,2) * ( pdf(AUp_,1) +pdf(AChm_,1) +pdf(AStr_,1) +pdf(ABot_,1) )    &
!                               +  pdf(AChm_,2)* ( pdf(AUp_,1) +pdf(ADn_,1)  +pdf(AStr_,1) +pdf(ABot_,1) )    &
!                               +  pdf(AStr_,2)* ( pdf(AUp_,1) +pdf(ADn_,1)  +pdf(AChm_,1) +pdf(ABot_,1) )    &
!                               +  pdf(ABot_,2)* ( pdf(AUp_,1) +pdf(ADn_,1)  +pdf(AChm_,1) +pdf(AStr_,1) )  ) &
! !                             &
!                     + HOp(2)/xE *( pdf_z(AUp_,2) * ( pdf(ADn_,1) +pdf(AChm_,1) +pdf(AStr_,1) +pdf(ABot_,1) )    &
!                                 +  pdf_z(ADn_,2) * ( pdf(AUp_,1) +pdf(AChm_,1) +pdf(AStr_,1) +pdf(ABot_,1) )    &
!                                 +  pdf_z(AChm_,2)* ( pdf(AUp_,1) +pdf(ADn_,1)  +pdf(AStr_,1) +pdf(ABot_,1) )    &
!                                 +  pdf_z(AStr_,2)* ( pdf(AUp_,1) +pdf(ADn_,1)  +pdf(AChm_,1) +pdf(ABot_,1) )    &
!                                 +  pdf_z(ABot_,2)* ( pdf(AUp_,1) +pdf(ADn_,1)  +pdf(AChm_,1) +pdf(AStr_,1) ) ) &
! !                                 &
!                     + HOp(3)/xE *( pdf(AUp_,2) * ( pdf_z(ADn_,1) +pdf_z(AChm_,1) +pdf_z(AStr_,1) +pdf_z(ABot_,1) )    &
!                                 +  pdf(ADn_,2) * ( pdf_z(AUp_,1) +pdf_z(AChm_,1) +pdf_z(AStr_,1) +pdf_z(ABot_,1) )    &
!                                 +  pdf(AChm_,2)* ( pdf_z(AUp_,1) +pdf_z(ADn_,1)  +pdf_z(AStr_,1) +pdf_z(ABot_,1) )    &
!                                 +  pdf(AStr_,2)* ( pdf_z(AUp_,1) +pdf_z(ADn_,1)  +pdf_z(AChm_,1) +pdf_z(ABot_,1) )    &
!                                 +  pdf(ABot_,2)* ( pdf_z(AUp_,1) +pdf_z(ADn_,1)  +pdf_z(AChm_,1) +pdf_z(AStr_,1) )  )


ELSEIF( PROCESS.EQ.18 ) THEN

   MomP(1:4,1) = MomExt(1:4,4)
   MomP(1:4,2) = MomExt(1:4,5)
   MomP(1:4,3) =-MomExt(1:4,1)
   MomP(1:4,4) =-MomExt(1:4,2)
   MomP(1:4,5) = MomExt(1:4,3)

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   IF( TOPDECAYS.GE.1 ) THEN
      xE = yRnd(16)
   ELSEIF( TOPDECAYS.EQ.0 ) THEN
      xE = yRnd(8)
   ENDIF

   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

   call sixquark_intdip(MomP(1:4,1:5),MomDK(1:4,1:6),xE,+1,+1,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**4 * PreFac

   EvalCS_1L_ttbqqbg = HOp(1) *( pdf(Up_,1) *pdf(Up_,2)  + pdf(Dn_,1) *pdf(Dn_,2)   &
                              +  pdf(Chm_,1)*pdf(Chm_,2) + pdf(Str_,1)*pdf(Str_,2)  &
                              +  pdf(Bot_,1)*pdf(Bot_,2)  ) &
!                             &
                    + HOp(2)/xE *( pdf_z(Up_,1) *pdf(Up_,2)  + pdf_z(Dn_,1) *pdf(Dn_,2)   &
                                +  pdf_z(Chm_,1)*pdf(Chm_,2) + pdf_z(Str_,1)*pdf(Str_,2)  &
                                +  pdf_z(Bot_,1)*pdf(Bot_,2) ) &
!                               &
                    + HOp(3)/xE *( pdf(Up_,1) *pdf_z(Up_,2)  + pdf(Dn_,1) *pdf_z(Dn_,2)   &
                                +  pdf(Chm_,1)*pdf_z(Chm_,2) + pdf(Str_,1)*pdf_z(Str_,2)  &
                                +  pdf(Bot_,1)*pdf_z(Bot_,2)  )



ELSEIF( PROCESS.EQ.19 ) THEN

   MomP(1:4,1) = MomExt(1:4,4)
   MomP(1:4,2) = MomExt(1:4,5)
   MomP(1:4,3) =-MomExt(1:4,1)
   MomP(1:4,4) =-MomExt(1:4,2)
   MomP(1:4,5) = MomExt(1:4,3)

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   IF( TOPDECAYS.GE.1 ) THEN
      xE = yRnd(16)
   ELSEIF( TOPDECAYS.EQ.0 ) THEN
      xE = yRnd(8)
   ENDIF

   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

   call sixquark_intdip(MomP(1:4,1:5),MomDK(1:4,1:6),xE,-1,-1,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**4 * PreFac

   EvalCS_1L_ttbqqbg = HOp(1) *( pdf(AUp_,1) *pdf(AUp_,2)  + pdf(ADn_,1) *pdf(ADn_,2)   &
                              +  pdf(AChm_,1)*pdf(AChm_,2) + pdf(AStr_,1)*pdf(AStr_,2)  &
                              +  pdf(ABot_,1)*pdf(ABot_,2)  ) &
!                             &
                    + HOp(2)/xE *( pdf_z(AUp_,1) *pdf(AUp_,2)  + pdf_z(ADn_,1) *pdf(ADn_,2)   &
                                +  pdf_z(AChm_,1)*pdf(AChm_,2) + pdf_z(AStr_,1)*pdf(AStr_,2)  &
                                +  pdf_z(ABot_,1)*pdf(ABot_,2) ) &
!                               &
                    + HOp(3)/xE *( pdf(AUp_,1) *pdf_z(AUp_,2)  + pdf(ADn_,1) *pdf_z(ADn_,2)   &
                                +  pdf(AChm_,1)*pdf_z(AChm_,2) + pdf(AStr_,1)*pdf_z(AStr_,2)  &
                                +  pdf(ABot_,1)*pdf_z(ABot_,2)  )


ENDIF
ENDIF


   if( IsNan(EvalCS_1L_ttbqqbg)  .or. dabs(EvalCS_1L_ttbqqbg).eq.1d0/0d0  ) then
        print *, "NAN:",EvalCS_1L_ttbqqbg
        call printYrnd(yrnd(:))
        print *, PSWgt , VgsWgt , PDFFac, sHatJacobi
        print *, eta1,eta2,MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,1)
        print *, MomExt(1:4,2)
        print *, MomExt(1:4,3)
        print *, MomExt(1:4,4)
        print *, MomExt(1:4,5)
!         stop
        EvalCS_1L_ttbqqbg = 0d0
        return
   endif


   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_ttbqqbg)
   enddo

   EvalCounter = EvalCounter + 1
   EvalCS_1L_ttbqqbg = EvalCS_1L_ttbqqbg/VgsWgt


!        MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!        MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!        MG_MOM(0:3,3) = MomExt(1:4,5)*100d0
!        MG_MOM(0:3,4) = MomExt(1:4,4)*100d0
!        MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!        call coupsm(0)
!       call SUUB_TTBG(MG_MOM,MadGraph_tree)
! !       call SUG_TTBU(MG_MOM,MadGraph_tree)
! !       call SUBG_TTBUB(MG_MOM,MadGraph_tree)
!       print *, ""
!       print *, "Unpolarized LO result:"
!       print *, "My tree:          ", LO_Res_Unpol  /10000d0
!       print *, "MadGraph hel.amp:", MadGraph_tree*(100d0)**2
!       print *, "ratio: ", MadGraph_tree*(100d0)**2/dble((LO_Res_Unpol))
!       pause

!       IF( PROCESS.EQ.6 ) THEN
!        print *, ""
!        print *, "Unpolarized LO result:",alpha_s*RunFactor
!        print *, "My tree:          ", LO_Res_Unpol/(4d0*DblPi*alpha_s*RunFactor)**3
!        print *, "DUW:             ",0.57903680015509d-4 *(100d0)**2
!        print *, "ratio: ", 0.57903680015509d-4*(100d0)**2/(LO_Res_Unpol/(4d0*DblPi*alpha_s*RunFactor)**3)
!        print *, ""
!
!        print *, "c_{-2)_bos:",NLO_Res_UnPol(-2)/LO_Res_UnPol
!        print *, "c_{-2)_fer:",NLO_Res_UnPol_Ferm(-2)/LO_Res_UnPol
!        print *, "c_{-2)_DUW:",(-0.096970419104695d0,0d0)
!        print *, "ratio to DUW:",(NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2))/LO_Res_UnPol / (-0.096970419104695d0)
!        print *, ""
!        print *, "c_{-1)_bos:",NLO_Res_UnPol(-1)/LO_Res_UnPol
!        print *, "c_{-1)_fer:",NLO_Res_UnPol_Ferm(-1)/LO_Res_UnPol
!        print *, "c_{-1)_DUW:",(-0.0126983208241662,0d0)
!        print *, "ratio to DUW:",(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/LO_Res_UnPol /(-0.0126983208241662d0)
! !        print *, "diff  to DUW:",((NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/LO_Res_UnPol-(-0.0126983208241662d0))/(alpha_sOver2Pi*RunFactor)
!        print *, ""
!        NLO_Res_UnPol(0) = NLO_Res_UnPol(0) - NLO_Res_UnPol(-2)*DblPi**2/6d0  ! different gamma function
!        print *, "c_{0)_bos:",(NLO_Res_UnPol(0)+NLO_Res_UnPol(1))/LO_Res_UnPol
!        print *, "c_{0)_fer:",(NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1))/LO_Res_UnPol
!        print *, "c_{0)_DUW:",(0.243567243908d0,0d0)
!        print *, "ratio to DUW:",(NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1)+NLO_Res_UnPol(0)+NLO_Res_UnPol(1))/LO_Res_UnPol/(0.243567243908d0)
! !       print *, "diff to DUW:",((NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1)+NLO_Res_UnPol(0)+NLO_Res_UnPol(1))/LO_Res_UnPol-(0.243567243908d0))/(alpha_sOver2Pi*RunFactor)
!        pause
!
!       ELSEIF( PROCESS.EQ.3 ) THEN
!       print *, ""
!       print *, "Unpolarized LO result:",alpha_s*RunFactor
!       print *, "My tree:          ", LO_Res_Unpol/(4d0*DblPi*alpha_s*RunFactor)**3
!       print *, "DUW:             ",0.1607845322071585d-4 *(100d0)**2
!       print *, "ratio: ", 0.1607845322071585d-4*(100d0)**2/(LO_Res_Unpol/(4d0*DblPi*alpha_s*RunFactor)**3)
!       print *, ""
!
!       print *, "c_{-2)_bos:",NLO_Res_UnPol(-2)/LO_Res_UnPol
!       print *, "c_{-2)_fer:",NLO_Res_UnPol_Ferm(-2)/LO_Res_UnPol
!       print *, "c_{-2)_DUW:",(-0.096970419104695d0,0d0)
!       print *, "ratio to DUW:",(NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2))/LO_Res_UnPol / (-0.096970419104695d0)
!       print *, ""
!       print *, "c_{-1)_bos:",NLO_Res_UnPol(-1)/LO_Res_UnPol
!       print *, "c_{-1)_fer:",NLO_Res_UnPol_Ferm(-1)/LO_Res_UnPol
!       print *, "c_{-1)_DUW:",(-0.0056430956994203d0,0d0)
!       print *, "ratio to DUW:",(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/LO_Res_UnPol /(-0.0056430956994203d0)
! !        print *, "diff  to DUW:",((NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/LO_Res_UnPol-(-0.0056430956994203d0))/(alpha_sOver2Pi*RunFactor)
!       print *, ""
!       NLO_Res_UnPol(0) = NLO_Res_UnPol(0) - NLO_Res_UnPol(-2)*DblPi**2/6d0  ! different gamma function
!       print *, "c_{0)_bos:",(NLO_Res_UnPol(0)+NLO_Res_UnPol(1))/LO_Res_UnPol
!       print *, "c_{0)_fer:",(NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1))/LO_Res_UnPol
!       print *, "c_{0)_DUW:",(0.4003849386477017d0,0d0)
!       print *, "ratio to DUW:",(NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1)+NLO_Res_UnPol(0)+NLO_Res_UnPol(1))/LO_Res_UnPol/0.4003849386477017d0
! !       print *, "diff to DUW:",((NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1)+NLO_Res_UnPol(0)+NLO_Res_UnPol(1))/LO_Res_UnPol-(0.4003849386477017d0))/(alpha_sOver2Pi*RunFactor)
!
!       pause
!
!       ELSEIF( PROCESS.EQ.4 ) THEN
!        print *, ""
!        print *, "Unpolarized LO result:",alpha_s*RunFactor
!        print *, "My tree:          ", LO_Res_Unpol/(4d0*DblPi*alpha_s*RunFactor)**3
!        print *, "DUW:             ",0.2603527972645622d-3 *(100d0)**2
!        print *, "ratio: ", 0.2603527972645622d-3*(100d0)**2/(LO_Res_Unpol/(4d0*DblPi*alpha_s*RunFactor)**3)
!        print *, ""
!
!        print *, "c_{-2)_bos:",NLO_Res_UnPol(-2)/LO_Res_UnPol
!        print *, "c_{-2)_fer:",NLO_Res_UnPol_Ferm(-2)/LO_Res_UnPol
!        print *, "c_{-2)_DUW:",(-0.096970419104695d0,0d0)
!        print *, "ratio to DUW:",(NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2))/LO_Res_UnPol / (-0.096970419104695d0)
!        print *, ""
!        print *, "c_{-1)_bos:",NLO_Res_UnPol(-1)/LO_Res_UnPol
!        print *, "c_{-1)_fer:",NLO_Res_UnPol_Ferm(-1)/LO_Res_UnPol
!        print *, "c_{-1)_DUW:",(0.083336273912803d0,0d0)
!        print *, "ratio to DUW:",(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/LO_Res_UnPol /(0.083336273912803d0)
! !        print *, "diff  to DUW:",((NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/LO_Res_UnPol-(0.083336273912803d0))/(alpha_sOver2Pi*RunFactor)
!        print *, ""
!        NLO_Res_UnPol(0) = NLO_Res_UnPol(0) - NLO_Res_UnPol(-2)*DblPi**2/6d0  ! different gamma function
!        print *, "c_{0)_bos:",(NLO_Res_UnPol(0)+NLO_Res_UnPol(1))/LO_Res_UnPol
!        print *, "c_{0)_fer:",(NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1))/LO_Res_UnPol
!        print *, "c_{0)_DUW:",(0.5384721403213878d0,0d0)
!        print *, "ratio to DUW:",(NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1)+NLO_Res_UnPol(0)+NLO_Res_UnPol(1))/LO_Res_UnPol/0.5384721403213878d0
! !       print *, "diff to DUW:",((NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1)+NLO_Res_UnPol(0)+NLO_Res_UnPol(1))/LO_Res_UnPol-(0.5384721403213878d0))/(alpha_sOver2Pi*RunFactor)
!
!        pause
!       ENDIF



return
END FUNCTION







FUNCTION EvalCS_Real_ttbgggg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModDipoles_GGTTBGG
use ifport
implicit none
real(8) ::  EvalCS_Real_ttbgggg,yRnd(1:VegasMxDim),VgsWgt,DipoleResult
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol
integer :: iHel,iPrimAmp,jPrimAmp,ParityFlip,NHisto,NBin(1:NumMaxHisto)
real(8) :: SpinAvg,ColorAvg,EHat,PSWgt,PSWgt2,PSWgt3,ISFac,RunFactor,PreFac
real(8) :: eta1,eta2,sHatJacobi,FluxFac,PDFFac
real(8) :: MomExt(1:4,1:6),MomDK(1:4,1:6),s12,sij,pdf(-6:6,1:2)
real(8) :: MG_MOM(0:3,1:6)
real(8) :: MadGraph_tree,GG_TTBGG,GG_TTBG,GG_TTB,GG_UUB
logical :: applyPSCut,applySingCut
real(8) :: weight,pjetout(1:4,1:6)
real(8),parameter :: ThresholdCutOff = 1.0d0
integer :: Njet
include "vegas_common.f"

IF( TopDecays.ge.1 ) THEN
  ParityFlip=1
ELSEIF( TopDecays.eq.0 ) THEN
  ParityFlip=2
ENDIF

!     print *, "fixed yRnd(:)"
!     print *,  yrnd(1:18)

!   yrnd(1) = 0.343119567652755
!   yrnd(2) = 0.338689352776245
!   yrnd(3) = 0.305442483306603
!   yrnd(4) = 0.220847354606189
!   yrnd(5) = 0.210255306777663
!   yrnd(6) = 0.185989333403292
!   yrnd(7) = 0.108586014764656
!   yrnd(8) = 0.367722897263115
!   yrnd(9) = 5.478499925452522E-002
!   yrnd(10) = 5.648925204597872E-003
!   yrnd(11) = 0.271482470804585
!   yrnd(12) = 0.494176612931386
!   yrnd(13) = 5.375115203380175E-002
!   yrnd(14) = 0.467031235092800
!   yrnd(15) = 0.391219816352809
!   yrnd(16) = 0.281488866210677
!   yrnd(17) = 1.933753002404123E-002
!   yrnd(18) = 3.966021772458232E-002


!   yrnd(1:18) = (/   9.491522648504568E-002,  0.325382507613596 ,      8.343304449358904E-002,  &
!   0.830449263260874    ,   0.891891869707152   ,    0.998966039346925, &
!   0.421468289226620   ,    0.171773186890813   ,    0.265246821430602, &
!   0.461662659062379    ,   0.826638621796715   ,    0.447198954078335, &
!   0.598441700876272    ,   0.627273513217870  ,     0.487814010579232, &
!   5.323968406189610E-002 , 0.479819143094906    ,   0.557825871948250 &
! /)



  EvalCS_Real_ttbgggg = 0d0
  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top * ThresholdCutOff ) then
      EvalCS_Real_ttbgggg = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)
! !DEC$ IF(_WriteTmpHisto .EQ.1)
!    if( it.ne.it_sav ) then
!       it_sav=it
! !       call WriteInstantHisto()
!       call WriteHisto(14,0d0,0d0,0d0)
!    endif
! !DEC$ ENDIF
   call EvalPhaseSpace_2to4(EHat,yRnd(3:10),MomExt(1:4,1:6),PSWgt)
   call boost2Lab(eta1,eta2,6,MomExt(1:4,1:6))
   ISFac = MomCrossing(MomExt)
   PSWgt2 = 1d0
   PSWgt3 = 1d0
IF( TopDecays.GE.1 ) THEN
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(11:14),.false.,MomDK(1:4,1:3),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,6),yRnd(15:18),.false.,MomDK(1:4,4:6),PSWgt3)
   PSWgt = PSWgt * PSWgt2*PSWgt3
   call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
   call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
ENDIF
   call CheckSing(MomExt,applySingCut)
   if( applySingCut ) then
       EvalCS_Real_ttbgggg = 0d0
       SkipCounter = SkipCounter + 1
       return
   endif
!    call Kinematics(4,MomExt,MomDK,applyPSCut,NBin)
   call Kinematics_TTBARJET(1,MomExt,MomDK,applyPSCut,NBin)

   call setPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   RunFactor = RunAlphaS(NLOParam,MuRen)
   if( applyPSCut ) then
       EvalCS_Real_ttbgggg = 0d0
   else
        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities
          call HelCrossing(Helicities(iHel,1:NumExtParticles))
          call SetPolarizations()

          do iPrimAmp=1,NumBornAmps
              call EvalTree(BornAmps(iPrimAmp))
          enddo

          LO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,NumBornAmps
          do iPrimAmp=1,NumBornAmps
              LO_Res_Pol = LO_Res_Pol + ParityFlip*ColLO_ttbgggg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result * dconjg(BornAmps(jPrimAmp)%Result)
          enddo
          enddo
          LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop

        LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**4
        EvalCS_Real_ttbgggg = LO_Res_Unpol * PreFac
        do NHisto=1,NumHistograms
               call intoHisto(NHisto,NBin(NHisto),EvalCS_Real_ttbgggg)
        enddo
        EvalCounter = EvalCounter + 1

endif!applyPSCut

    PreFac = PreFac * ISFac * (alpha_s4Pi*RunFactor)**4 / PSWgt2/PSWgt3
    call EvalDipoles_GGTTBGG((/MomExt(1:4,5),MomExt(1:4,6),-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4)/),yRnd(11:18),PreFac,DipoleResult)


!     call EvalDipoles_GGTTBGG((/MomExt(1:4,5),MomExt(1:4,6),-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4)/),DipoleResult)
!     DipoleResult = - PreFac * DipoleResult * ISFac * (alpha_s4Pi*RunFactor)**4

! !     gg->ttbg
!       MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!       MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!       MG_MOM(0:3,3) = MomExt(1:4,6)*100d0
!       MG_MOM(0:3,4) = MomExt(1:4,5)*100d0
!       MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!       MG_MOM(0:3,6) = MomExt(1:4,4)*100d0
!       call coupsm(0)
!       call SGG_TTBGG(MG_MOM,MadGraph_tree)
!       MadGraph_tree = MadGraph_tree*(100d0)**4
!        print *, alpha_s*RunFactor,m_Top
!        print *, "Unpolarized LO result:"
!        print *, "My tree:", LO_Res_Unpol
!        print *, "MadGraph:", MadGraph_tree
!        print *, "MG  ratio: ", MadGraph_tree/dble(LO_Res_Unpol)
!        pause
!         LO_Res_Unpol = MadGraph_tree
!         EvalCS_Real_ttbgggg = LO_Res_Unpol * PreFac

!      s12 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,2))
!      sij = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,3))

!      print *, "E4: ",MomExt(1,4)/dsqrt(s12)
!      print *, "si12",2d0*(MomExt(1:4,1).dot.MomExt(1:4,4))/s12,2d0*(MomExt(1:4,2).dot.MomExt(1:4,4))/s12
!      print *, "R: ",get_R(MomExt(1:4,4),MomExt(1:4,3))

!      write(* ,"(1PE23.16,3X,1PE23.16,3X,1PE23.16,3X,1PE23.16)") MomExt(1,4)/dsqrt(s12),EvalCS_Real_ttbgggg,DipoleResult, dabs(1d0+EvalCS_Real_ttbgggg/(DipoleResult))
!        write(* ,"(1PE23.16,3X,1PE23.16,3X,1PE23.16,3X,1PE23.16)") sij/s12,EvalCS_Real_ttbgggg,DipoleResult, dabs(1d0+EvalCS_Real_ttbgggg/(DipoleResult))
!        pause


    if( IsNan(EvalCS_Real_ttbgggg) ) then
          print *, "NAN:",EvalCS_Real_ttbgggg,DipoleResult
          call printYrnd(yrnd(:))
          print *, PSWgt , VgsWgt , PDFFac, sHatJacobi
          print *, eta1,eta2,MuFac,EHat
          print *, "Mom"
          print *, MomExt(1:4,1)
          print *, MomExt(1:4,2)
          print *, MomExt(1:4,3)
          print *, MomExt(1:4,4)
          print *, MomExt(1:4,5)
          print *, MomExt(1:4,6)
          EvalCS_Real_ttbgggg = 0d0
          return
    endif

    EvalCS_Real_ttbgggg = (EvalCS_Real_ttbgggg + DipoleResult)/VgsWgt


return
END FUNCTION









FUNCTION EvalCS_Real_ttbqqbgg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes
use ModParameters
use ModMyRecurrence
use ModDipoles_GGTTBQQB
use ModDipoles_QQBTTBGG
use ModDipoles_QGTTBQG
use ModDipoles_QBGTTBQBG
implicit none
real(8) ::  EvalCS_Real_ttbqqbgg,EvalCS_Dips_ttbqqbgg,yRnd(1:VegasMxDim),VgsWgt,DipoleResult
complex(8) :: A(1:12),Ac(1:12)
complex(8) :: LO_Res_Pol,LO_Res_Unpol
integer :: iHel,iPrimAmp,jPrimAmp,ParityFlip,NPDF
real(8) :: SpinAvg,ColorAvg,EHat,PSWgt,ISFac,RunFactor,PreFac,PSWgt2,PSWgt3
real(8) :: eta1,eta2,sHatJacobi,FluxFac,PDFFac_a,PDFFac_b,PDFFac,pdf(-6:6,1:2)
real(8) :: MomExt(1:4,1:6),MomDK(1:4,1:6),s12,sij
real(8) :: MG_MOM(0:3,1:6)
real(8) :: MadGraph_tree
logical :: applyPSCut,applySingCut
integer :: NBin(1:NumMaxHisto),NHisto
real(8),parameter :: ThresholdCutOff=1d0
include "vegas_common.f"

  ParityFlip=1
! print *, "yRnd fixed"


  EvalCS_Real_ttbqqbgg = 0d0
  EvalCS_Dips_ttbqqbgg = 0d0

  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top * ThresholdCutOff ) then
      EvalCS_Real_ttbqqbgg = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)
! !DEC$ IF(_WriteTmpHisto .EQ.1)
!    if( it.ne.it_sav ) then
!       it_sav=it
! !       call WriteInstantHisto()
!       call WriteHisto(14,0d0,0d0,0d0)
!    endif
! !DEC$ ENDIF

   call EvalPhaseSpace_2to4(EHat,yRnd(3:10),MomExt(1:4,1:6),PSWgt)
   call boost2Lab(eta1,eta2,6,MomExt(1:4,1:6))

   PSWgt2 = 1d0
   PSWgt3 = 1d0
IF( TopDecays.GE.1 ) THEN
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(11:14),.false.,MomDK(1:4,1:3),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,6),yRnd(15:18),.false.,MomDK(1:4,4:6),PSWgt3)
   PSWgt = PSWgt * PSWgt2*PSWgt3
ENDIF
   call CheckSing(MomExt,applySingCut)
   if( applySingCut ) then
      EvalCS_Real_ttbqqbgg = 0d0
      return
   endif

!    call Kinematics(4,MomExt,MomDK,applyPSCut,NBin)
   call Kinematics_TTBARJET(1,MomExt,MomDK,applyPSCut,NBin)

   call setPDFs(eta1,eta2,MuFac,pdf)
   IF( PROCESS.EQ.10 ) THEN
      PDFFac_a = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
               + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
               + pdf(Bot_,1)*pdf(ABot_,2)
      PDFFac_b = pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
               + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
               + pdf(Bot_,2)*pdf(ABot_,1)
   ELSEIF( PROCESS.EQ.11 ) THEN
      PDFFac = pdf(0,1) * pdf(0,2)
      PDFFac = PDFFac * nf_light   ! nf_light accounts for all light quarks in the final state
      PDFFac_a = PDFFac
      PDFFac_b = 1d300   ! unused
   ELSEIF( PROCESS.EQ.12 ) THEN
      PDFFac_a = pdf(Up_,1) *pdf(0,2) + pdf(Dn_,1) *pdf(0,2)   &
                    + pdf(Chm_,1)*pdf(0,2) + pdf(Str_,1)*pdf(0,2)   &
                    + pdf(Bot_,1)*pdf(0,2)
      PDFFac_b =      pdf(Up_,2) *pdf(0,1) + pdf(Dn_,2) *pdf(0,1)   &
                    + pdf(Chm_,2)*pdf(0,1) + pdf(Str_,2)*pdf(0,1)   &
                    + pdf(Bot_,2)*pdf(0,1)
   ELSEIF( PROCESS.EQ.13 ) THEN
      PDFFac_a = pdf(AUp_,1) *pdf(0,2) + pdf(ADn_,1) *pdf(0,2)   &
               + pdf(AChm_,1)*pdf(0,2) + pdf(AStr_,1)*pdf(0,2)   &
               + pdf(ABot_,1)*pdf(0,2)
      PDFFac_b = pdf(AUp_,2) *pdf(0,1) + pdf(ADn_,2) *pdf(0,1)   &
               + pdf(AChm_,2)*pdf(0,1) + pdf(AStr_,2)*pdf(0,1)   &
               + pdf(ABot_,2)*pdf(0,1)
   ENDIF
   RunFactor = RunAlphaS(NLOParam,MuRen)

DO NPDF=1,2
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   if(npdf.eq.1) then
        PDFFac = PDFFac_a
   elseif(npdf.eq.2) then
        if( PROCESS.EQ.11 ) cycle
        PDFFac = PDFFac_b
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
   endif
   ISFac = MomCrossing(MomExt)
IF( TopDecays.GE.1 ) THEN
   call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
   call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
ENDIF

if( .not. applyPSCut ) then
    LO_Res_Unpol = (0d0,0d0)
    do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()

      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
          A(iPrimAmp) = BornAmps(iPrimAmp)%Result
          Ac(iPrimAmp)= dconjg( A(iPrimAmp) )
      enddo
      A(7:12)  = - A(7:12)
      Ac(7:12) = -Ac(7:12)

      LO_Res_Pol = (0d0,0d0)
      LO_Res_Pol = LO_Res_Pol + A(1)*Ac(1) * ( 64.D0 )
      LO_Res_Pol = LO_Res_Pol + A(1)*Ac(2) * (  - 8.D0 )
      LO_Res_Pol = LO_Res_Pol + A(1)*Ac(3) * ( 8.D0 )
      LO_Res_Pol = LO_Res_Pol + A(1)*Ac(4) * ( 8.D0 )
      LO_Res_Pol = LO_Res_Pol + A(1)*Ac(7) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(1)*Ac(8) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(1)*Ac(9) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(1)*Ac(10) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(1)*Ac(11) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(1)*Ac(12) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(2)*Ac(1) * (  - 8.D0 )
      LO_Res_Pol = LO_Res_Pol + A(2)*Ac(2) * ( 64.D0 )
      LO_Res_Pol = LO_Res_Pol + A(2)*Ac(3) * ( 8.D0 )
      LO_Res_Pol = LO_Res_Pol + A(2)*Ac(4) * ( 8.D0 )
      LO_Res_Pol = LO_Res_Pol + A(2)*Ac(7) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(2)*Ac(8) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(2)*Ac(9) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(2)*Ac(10) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(2)*Ac(11) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(2)*Ac(12) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(3)*Ac(1) * ( 8.D0 )
      LO_Res_Pol = LO_Res_Pol + A(3)*Ac(2) * ( 8.D0 )
      LO_Res_Pol = LO_Res_Pol + A(3)*Ac(3) * ( 64.D0 )
      LO_Res_Pol = LO_Res_Pol + A(3)*Ac(4) * (  - 8.D0 )
      LO_Res_Pol = LO_Res_Pol + A(3)*Ac(7) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(3)*Ac(8) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(3)*Ac(9) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(3)*Ac(10) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(3)*Ac(11) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(3)*Ac(12) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(4)*Ac(1) * ( 8.D0 )
      LO_Res_Pol = LO_Res_Pol + A(4)*Ac(2) * ( 8.D0 )
      LO_Res_Pol = LO_Res_Pol + A(4)*Ac(3) * (  - 8.D0 )
      LO_Res_Pol = LO_Res_Pol + A(4)*Ac(4) * ( 64.D0 )
      LO_Res_Pol = LO_Res_Pol + A(4)*Ac(7) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(4)*Ac(8) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(4)*Ac(9) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(4)*Ac(10) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(4)*Ac(11) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(4)*Ac(12) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(5)*Ac(5) * ( 64.D0 )
      LO_Res_Pol = LO_Res_Pol + A(5)*Ac(6) * ( 8.D0 )
      LO_Res_Pol = LO_Res_Pol + A(5)*Ac(7) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(5)*Ac(8) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(5)*Ac(9) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(5)*Ac(10) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(5)*Ac(11) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(5)*Ac(12) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(6)*Ac(5) * ( 8.D0 )
      LO_Res_Pol = LO_Res_Pol + A(6)*Ac(6) * ( 64.D0 )
      LO_Res_Pol = LO_Res_Pol + A(6)*Ac(7) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(6)*Ac(8) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(6)*Ac(9) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(6)*Ac(10) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(6)*Ac(11) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(6)*Ac(12) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(7)*Ac(1) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(7)*Ac(2) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(7)*Ac(3) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(7)*Ac(4) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(7)*Ac(5) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(7)*Ac(6) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(7)*Ac(7) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(7)*Ac(8) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(7)*Ac(11) * ( 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(7)*Ac(12) * ( 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(8)*Ac(1) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(8)*Ac(2) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(8)*Ac(3) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(8)*Ac(4) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(8)*Ac(5) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(8)*Ac(6) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(8)*Ac(7) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(8)*Ac(8) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(8)*Ac(11) * ( 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(8)*Ac(12) * ( 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(9)*Ac(1) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(9)*Ac(2) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(9)*Ac(3) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(9)*Ac(4) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(9)*Ac(5) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(9)*Ac(6) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(9)*Ac(9) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(9)*Ac(10) * ( 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(10)*Ac(1) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(10)*Ac(2) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(10)*Ac(3) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(10)*Ac(4) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(10)*Ac(5) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(10)*Ac(6) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(10)*Ac(9) * ( 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(10)*Ac(10) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(11)*Ac(1) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(11)*Ac(2) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(11)*Ac(3) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(11)*Ac(4) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(11)*Ac(5) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(11)*Ac(6) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(11)*Ac(7) * ( 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(11)*Ac(8) * ( 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(11)*Ac(11) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(11)*Ac(12) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(12)*Ac(1) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(12)*Ac(2) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(12)*Ac(3) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(12)*Ac(4) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(12)*Ac(5) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(12)*Ac(6) * ( 64.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(12)*Ac(7) * ( 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(12)*Ac(8) * ( 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(12)*Ac(11) * (  - 8.D0/9.D0 )
      LO_Res_Pol = LO_Res_Pol + A(12)*Ac(12) * ( 64.D0/9.D0 )
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol * PDFFac
    enddo! helicity loop

    LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**4 * PreFac
    do NHisto=1,NumHistograms
        call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol))
    enddo
    EvalCounter = EvalCounter + 1
    EvalCS_Real_ttbqqbgg = EvalCS_Real_ttbqqbgg + LO_Res_Unpol
endif!applyPSCut

    PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * ISFac * (alpha_s4Pi*RunFactor)**4   * PDFFac  / PSWgt2/PSWgt3
    IF( PROCESS.EQ.10 ) THEN
        call EvalDipoles_QQBTTBGG((/MomExt(1:4,5),MomExt(1:4,6),-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4)/),yRnd(11:18),PreFac,DipoleResult)
    ELSEIF( PROCESS.EQ.11 ) THEN
        call EvalDipoles_GGTTBQQB((/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,3),MomExt(1:4,4),-MomExt(1:4,1),-MomExt(1:4,2)/),yRnd(11:18),PreFac,DipoleResult)
    ELSEIF( PROCESS.EQ.12 ) THEN
        call EvalDipoles_QGTTBQG((/MomExt(1:4,5),MomExt(1:4,6),-MomExt(1:4,1),MomExt(1:4,3),-MomExt(1:4,2),MomExt(1:4,4)/),yRnd(11:18),PreFac,DipoleResult)
    ELSEIF( PROCESS.EQ.13 ) THEN
        call EvalDipoles_QBGTTBQBG((/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,3),-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,4)/),yRnd(11:18),PreFac,DipoleResult)
    ENDIF
    EvalCS_Dips_ttbqqbgg = EvalCS_Dips_ttbqqbgg + DipoleResult
ENDDO! loop over a<-->b pdfs

!       MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!       MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!       MG_MOM(0:3,3) = MomExt(1:4,6)*100d0
!       MG_MOM(0:3,4) = MomExt(1:4,5)*100d0
!       MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!       MG_MOM(0:3,6) = MomExt(1:4,4)*100d0
!       call coupsm(0)
!       call SUUB_TTBGG(MG_MOM,MadGraph_tree)
!       MadGraph_tree =MadGraph_tree*1d8

!       MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!       MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!       MG_MOM(0:3,3) = MomExt(1:4,6)*100d0
!       MG_MOM(0:3,4) = MomExt(1:4,5)*100d0
!       MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!       MG_MOM(0:3,6) = MomExt(1:4,4)*100d0
!       call coupsm(0)
!       call SUBG_TTBUBG(MG_MOM,MadGraph_tree)
!       MadGraph_tree =MadGraph_tree*1d8

!       MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!       MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!       MG_MOM(0:3,3) = MomExt(1:4,6)*100d0
!       MG_MOM(0:3,4) = MomExt(1:4,5)*100d0
!       MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!       MG_MOM(0:3,6) = MomExt(1:4,4)*100d0
!       call coupsm(0)
!       call SUG_TTBUG(MG_MOM,MadGraph_tree)
!       MadGraph_tree =MadGraph_tree*1d8

!       MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!       MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!       MG_MOM(0:3,3) = MomExt(1:4,6)*100d0
!       MG_MOM(0:3,4) = MomExt(1:4,5)*100d0
!       MG_MOM(0:3,5) = MomExt(1:4,4)*100d0
!       MG_MOM(0:3,6) = MomExt(1:4,3)*100d0
!       call coupsm(0)
!       call SGG_TTBUUB(MG_MOM,MadGraph_tree)
!       MadGraph_tree = MadGraph_tree*(100d0)**4



!        print *, ""
!        print *, alpha_s*RunFactor,m_Top
!        print *, "Unpolarized LO result:"
!        print *, "My tree:", LO_Res_Unpol
!        print *, "MadGraph:", MadGraph_tree
!        print *, "MG  ratio: ", MadGraph_tree/dble(LO_Res_Unpol)
!        pause

!      s12 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,2))
!      sij = 2d0*(MomExt(1:4,2).dot.MomExt(1:4,4))
!      print *, "E4: ",MomExt(1,4)/dsqrt(s12)
!      print *, "si12",2d0*(MomExt(1:4,1).dot.MomExt(1:4,4))/s12,2d0*(MomExt(1:4,2).dot.MomExt(1:4,4))/s12
!      print *, "R: ",get_R(MomExt(1:4,4),MomExt(1:4,3))

!      write(* ,"(1PE23.16,3X,1PE23.16,3X,1PE23.16,3X,1PE23.16)") MomExt(1,4)/dsqrt(s12),EvalCS_Real_ttbqqbgg,EvalCS_Dips_ttbqqbgg,dabs(1d0+EvalCS_Real_ttbqqbgg/(EvalCS_Dips_ttbqqbgg))
!      write(* ,"(1PE23.16,3X,1PE23.16,3X,1PE23.16,3X,1PE23.16)") sij/s12,EvalCS_Real_ttbqqbgg,EvalCS_Dips_ttbqqbgg, dabs(1d0+EvalCS_Real_ttbqqbgg/(EvalCS_Dips_ttbqqbgg))
!      pause

    EvalCS_Real_ttbqqbgg = (EvalCS_Real_ttbqqbgg + EvalCS_Dips_ttbqqbgg) /VgsWgt

   if( IsNan(EvalCS_Real_ttbqqbgg) ) then
        print *, "NAN:",EvalCS_Real_ttbqqbgg,EvalCS_Dips_ttbqqbgg
        call printYrnd(yrnd(:))
        print *, PSWgt , VgsWgt , PDFFac, sHatJacobi
        print *, eta1,eta2,MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,1)
        print *, MomExt(1:4,2)
        print *, MomExt(1:4,3)
        print *, MomExt(1:4,4)
        print *, MomExt(1:4,5)
        print *, MomExt(1:4,6)
        EvalCS_Real_ttbqqbgg = 0d0
        return
   endif

return
END FUNCTION





FUNCTION EvalCS_Real_ttbqqbqqb(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes
use ModParameters
use ModMyRecurrence
use ModSixFermionProcs
use ifport
! use ModDipoles_QQpBTTBQQpB
! use ModDipoles_QQpTTBQQp
! use ModDipoles_QBQpBTTBQBQpB
implicit none
real(8) ::  EvalCS_Real_ttbqqbqqb,yRnd(1:VegasMxDim),VgsWgt,DipoleResult,res
real(8) :: SpinAvg,ColorAvg,EHat,PSWgt,ISFac,RunFactor,PreFac
real(8) :: LO_Res_Unpol
real(8) :: MomExt(1:4,1:6)
real(8) :: s12,sij,eta1,eta2,FluxFac,sHatJacobi
real(8) :: pdf(-6:6,1:2),PDFFac
real(8) :: MG_MOM(0:3,1:6)
real(8) :: MadGraph_tree
logical :: applyPSCut,applySingCut
! include 'global_import'
real(8) :: weight,PJetOut(1:4,1:6)
real(8),parameter :: ThresholdCutOff = 1.0d0
integer :: NJet,NBin(1:NumMaxHisto),NHisto,ParityFlip
include 'vegas_common.f'


!   print *, "fixed yRnd"
!   yrnd(1) = 0.292000321574510d0
!   yrnd(2) = 0.454017554621221d0
!   yrnd(3) = 0.121794757024290d0
!   yrnd(4) = 0.235149903332419d0
!   yrnd(5) = 3.478175263608888d-002
!   yrnd(6) = 0.308248961255071d0
!   yrnd(7) = 0.173040518850573d0
!   yrnd(8) = 0.466578881240720d0
!   yrnd(9) = 0.138669799612216d0
!   yrnd(10) = 0.164425307961379d0
!   yrnd(11) = 0.308292174622553d0
!   yrnd(12) = 0.184182640949349d0
!   yrnd(13) = 5.764643571229954d-002
!   yrnd(14) = 0.240291813966023d0
!   yrnd(15) = 0.174240712623224d0
!   yrnd(16) = 0.381112783160579d0
!   yrnd(17) = 0.368773503167915d0
!   yrnd(18) = 0.118900980157266d0


  EvalCS_Real_ttbqqbqqb = 0d0

  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top * ThresholdCutOff ) then
      EvalCS_Real_ttbqqbqqb = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)
! !DEC$ IF(_WriteTmpHisto .EQ.1)
!    if( it.ne.it_sav ) then
!       it_sav=it
! !       call WriteInstantHisto()
!       call WriteHisto(14,0d0,0d0,0d0)
!    endif
! !DEC$ ENDIF


   call EvalPhaseSpace_2to4(EHat,yRnd(3:10),MomExt(1:4,1:6),PSWgt)
   call boost2Lab(eta1,eta2,6,MomExt(1:4,1:6))
   ISFac = MomCrossing(MomExt)
!    PSWgt2 = 1d0
!    PSWgt3 = 1d0
! IF( TOPDECAYS.EQ.1 ) THEN     ! this is done inside of the sixfermion routines
!    call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(11:14),.false.,MomDK(1:4,1:3),PSWgt2)
!    call EvalPhasespace_TopDecay(MomExt(1:4,6),yRnd(15:18),.false.,MomDK(1:4,4:6),PSWgt3)
! !    PSWgt = PSWgt * PSWgt2*PSWgt3  is multiplied inside of the sixfermion routines
!    call TopDecay(ExtParticle(1),.false.,0,MomDK(1:4,1:3))
!    call TopDecay(ExtParticle(2),.false.,0,MomDK(1:4,4:6))
! ENDIF
   call CheckSing(MomExt,applySingCut)
   if( applySingCut ) then
      EvalCS_Real_ttbqqbqqb = 0d0
      return
   endif

!    call Kinematics(4,MomExt,MomDK,applyPSCut,NBin)
    call setPDFs(eta1,eta2,MuFac,pdf)



    LO_Res_Unpol = 0d0
    IF( PROCESS.EQ.14 ) THEN
       PDFFac =  pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
               + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
               + pdf(Bot_,1)*pdf(ABot_,2)
       RunFactor = RunAlphaS(NLOParam,MuRen)
       PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * ISFac * PDFFac * (alpha_s4Pi*RunFactor)**4
       call sixquark((/MomExt(1:4,5),MomExt(1:4,6),-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4)/),yRnd(11:18),PreFac,+1,-1,.true.,Res,DipoleResult)
       LO_Res_Unpol = Res


       PDFFac = ( pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
                + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
                + pdf(Bot_,2)*pdf(ABot_,1) )
       PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * ISFac * PDFFac * (alpha_s4Pi*RunFactor)**4
       call sixquark((/MomExt(1:4,5),MomExt(1:4,6),-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4)/),yRnd(11:18),PreFac,-1,+1,.true.,Res,DipoleResult)
       LO_Res_Unpol = LO_Res_Unpol + Res

       EvalCS_Real_ttbqqbqqb = LO_Res_Unpol / VgsWgt
!        print *, EvalCS_Real_ttbqqbqqb

!         MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!         MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!         MG_MOM(0:3,3) = MomExt(1:4,6)*100d0
!         MG_MOM(0:3,4) = MomExt(1:4,5)*100d0
!         MG_MOM(0:3,5) = MomExt(1:4,4)*100d0
!         MG_MOM(0:3,6) = MomExt(1:4,3)*100d0
!         call coupsm(0)
!         call SUUB_TTBUUB(MG_MOM,PreFac)
!         MadGraph_tree = PreFac
!         MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!         MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!         MG_MOM(0:3,3) = MomExt(1:4,6)*100d0
!         MG_MOM(0:3,4) = MomExt(1:4,5)*100d0
!         MG_MOM(0:3,5) = MomExt(1:4,4)*100d0
!         MG_MOM(0:3,6) = MomExt(1:4,3)*100d0
!         call coupsm(0)
!         call SUUB_TTBDDB(MG_MOM,PreFac)
!         MadGraph_tree = MadGraph_tree + PreFac*4d0
!         MadGraph_tree = MadGraph_tree*(100d0)**4

    ELSEIF( PROCESS.EQ.15 ) THEN
       PDFFac =  pdf(Up_,1) * ( pdf(ADn_,2) +pdf(AChm_,2) +pdf(AStr_,2) +pdf(ABot_,2) )    &
              +  pdf(Dn_,1) * ( pdf(AUp_,2) +pdf(AChm_,2) +pdf(AStr_,2) +pdf(ABot_,2) )    &
              +  pdf(Chm_,1)* ( pdf(AUp_,2) +pdf(ADn_,2)  +pdf(AStr_,2) +pdf(ABot_,2) )    &
              +  pdf(Str_,1)* ( pdf(AUp_,2) +pdf(ADn_,2)  +pdf(AChm_,2) +pdf(ABot_,2) )    &
              +  pdf(Bot_,1)* ( pdf(AUp_,2) +pdf(ADn_,2)  +pdf(AChm_,2) +pdf(AStr_,2) )
       RunFactor = RunAlphaS(NLOParam,MuRen)
       PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac * ISFac * (alpha_s4Pi*RunFactor)**4
       call sixquark((/MomExt(1:4,5),MomExt(1:4,6),-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,4),MomExt(1:4,3)/),yRnd(11:18),PreFac,+1,-2,.true.,res,DipoleResult)
       LO_Res_Unpol = Res

       PDFFac =  pdf(Up_,2) * ( pdf(ADn_,1) +pdf(AChm_,1) +pdf(AStr_,1) +pdf(ABot_,1) )    &
              +  pdf(Dn_,2) * ( pdf(AUp_,1) +pdf(AChm_,1) +pdf(AStr_,1) +pdf(ABot_,1) )    &
              +  pdf(Chm_,2)* ( pdf(AUp_,1) +pdf(ADn_,1)  +pdf(AStr_,1) +pdf(ABot_,1) )    &
              +  pdf(Str_,2)* ( pdf(AUp_,1) +pdf(ADn_,1)  +pdf(AChm_,1) +pdf(ABot_,1) )    &
              +  pdf(Bot_,2)* ( pdf(AUp_,1) +pdf(ADn_,1)  +pdf(AChm_,1) +pdf(AStr_,1) )
       PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac * ISFac * (alpha_s4Pi*RunFactor)**4
       call sixquark((/MomExt(1:4,5),MomExt(1:4,6),-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,4),MomExt(1:4,3)/),yRnd(11:18),PreFac,-2,+1,.true.,res,DipoleResult)
       LO_Res_Unpol = LO_Res_Unpol + Res
       EvalCS_Real_ttbqqbqqb = LO_Res_Unpol / VgsWgt
!         LO_Res_Unpol = LO_Res_Unpol/PreFac * ISFac * (alpha_s4Pi*RunFactor)**4
!         MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!         MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!         MG_MOM(0:3,3) = MomExt(1:4,6)*100d0
!         MG_MOM(0:3,4) = MomExt(1:4,5)*100d0
!         MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!         MG_MOM(0:3,6) = MomExt(1:4,4)*100d0
!         call coupsm(0)
!         call SUDB_TTBUDB(MG_MOM,MadGraph_tree)
!         MadGraph_tree = MadGraph_tree*(100d0)**4

    ELSEIF( PROCESS.EQ.16 ) THEN
       PDFFac =  pdf(Up_,1) * ( pdf(Dn_,2) +pdf(Chm_,2) +pdf(Str_,2) +pdf(Bot_,2) )    &
              +  pdf(Dn_,1) * ( pdf(Up_,2) +pdf(Chm_,2) +pdf(Str_,2) +pdf(Bot_,2) )    &
              +  pdf(Chm_,1)* ( pdf(Up_,2) +pdf(Dn_,2)  +pdf(Str_,2) +pdf(Bot_,2) )    &
              +  pdf(Str_,1)* ( pdf(Up_,2) +pdf(Dn_,2)  +pdf(Chm_,2) +pdf(Bot_,2) )    &
              +  pdf(Bot_,1)* ( pdf(Up_,2) +pdf(Dn_,2)  +pdf(Chm_,2) +pdf(Str_,2) )
       RunFactor = RunAlphaS(NLOParam,MuRen)
       PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac * ISFac * (alpha_s4Pi*RunFactor)**4
       call SixQuark((/MomExt(1:4,5),MomExt(1:4,6),-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4)/),yRnd(11:18),PreFac,+1,+2,.true.,Res,DipoleResult)
       LO_Res_Unpol = Res

!        PDFFac =  pdf(Up_,2) * ( pdf(Dn_,1) +pdf(Chm_,1) +pdf(Str_,1) +pdf(Bot_,1) )    &
!               +  pdf(Dn_,2) * ( pdf(Up_,1) +pdf(Chm_,1) +pdf(Str_,1) +pdf(Bot_,1) )    &
!               +  pdf(Chm_,2)* ( pdf(Up_,1) +pdf(Dn_,1)  +pdf(Str_,1) +pdf(Bot_,1) )    &
!               +  pdf(Str_,2)* ( pdf(Up_,1) +pdf(Dn_,1)  +pdf(Chm_,1) +pdf(Bot_,1) )    &
!               +  pdf(Bot_,2)* ( pdf(Up_,1) +pdf(Dn_,1)  +pdf(Chm_,1) +pdf(Str_,1) )
!        PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac * ISFac * (alpha_s4Pi*RunFactor)**4
!        call SixQuark((/MomExt(1:4,5),MomExt(1:4,6),-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4)/),yRnd(11:18),PreFac,+2,+1,.true.,Res,DipoleResult)
!        LO_Res_Unpol = LO_Res_Unpol + Res

       EvalCS_Real_ttbqqbqqb = LO_Res_Unpol / VgsWgt
!         MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!         MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!         MG_MOM(0:3,3) = MomExt(1:4,6)*100d0
!         MG_MOM(0:3,4) = MomExt(1:4,5)*100d0
!         MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!         MG_MOM(0:3,6) = MomExt(1:4,4)*100d0
!         call coupsm(0)
!         call SUD_TTBUD(MG_MOM,MadGraph_tree)
!         MadGraph_tree = MadGraph_tree*(100d0)**4
    ELSEIF( PROCESS.EQ.17 ) THEN
       PDFFac =  pdf(AUp_,1) * ( pdf(ADn_,2) +pdf(AChm_,2) +pdf(AStr_,2) +pdf(ABot_,2) )    &
              +  pdf(ADn_,1) * ( pdf(AUp_,2) +pdf(AChm_,2) +pdf(AStr_,2) +pdf(ABot_,2) )    &
              +  pdf(AChm_,1)* ( pdf(AUp_,2) +pdf(ADn_,2)  +pdf(AStr_,2) +pdf(ABot_,2) )    &
              +  pdf(AStr_,1)* ( pdf(AUp_,2) +pdf(ADn_,2)  +pdf(AChm_,2) +pdf(ABot_,2) )    &
              +  pdf(ABot_,1)* ( pdf(AUp_,2) +pdf(ADn_,2)  +pdf(AChm_,2) +pdf(AStr_,2) )
       RunFactor = RunAlphaS(NLOParam,MuRen)
       PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac * ISFac * (alpha_s4Pi*RunFactor)**4
       call sixquark((/MomExt(1:4,5),MomExt(1:4,6),-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,4),MomExt(1:4,3)/),yRnd(11:18),PreFac,-1,-2,.true.,Res,DipoleResult)
       LO_Res_Unpol = Res

!        PDFFac =  pdf(AUp_,2) * ( pdf(ADn_,1) +pdf(AChm_,1) +pdf(AStr_,1) +pdf(ABot_,1) )    &
!               +  pdf(ADn_,2) * ( pdf(AUp_,1) +pdf(AChm_,1) +pdf(AStr_,1) +pdf(ABot_,1) )    &
!               +  pdf(AChm_,2)* ( pdf(AUp_,1) +pdf(ADn_,1)  +pdf(AStr_,1) +pdf(ABot_,1) )    &
!               +  pdf(AStr_,2)* ( pdf(AUp_,1) +pdf(ADn_,1)  +pdf(AChm_,1) +pdf(ABot_,1) )    &
!               +  pdf(ABot_,2)* ( pdf(AUp_,1) +pdf(ADn_,1)  +pdf(AChm_,1) +pdf(AStr_,1) )
!        PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac * ISFac * (alpha_s4Pi*RunFactor)**4
!        call sixquark((/MomExt(1:4,5),MomExt(1:4,6),-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,4),MomExt(1:4,3)/),yRnd(11:18),PreFac,-2,-1,.true.,Res,DipoleResult)
!        LO_Res_Unpol = LO_Res_Unpol + Res

       EvalCS_Real_ttbqqbqqb = LO_Res_Unpol / VgsWgt
!         MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!         MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!         MG_MOM(0:3,3) = MomExt(1:4,6)*100d0
!         MG_MOM(0:3,4) = MomExt(1:4,5)*100d0
!         MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!         MG_MOM(0:3,6) = MomExt(1:4,4)*100d0
!         call coupsm(0)
!         call SUBDB_TTBUBDB(MG_MOM,MadGraph_tree)
!         MadGraph_tree = MadGraph_tree*(100d0)**4
    ELSEIF( PROCESS.EQ.18 ) THEN
       PDFFac = ( pdf(Up_,1) *pdf(Up_,2)  + pdf(Dn_,1) *pdf(Dn_,2)   &
                + pdf(Chm_,1)*pdf(Chm_,2) + pdf(Str_,1)*pdf(Str_,2)  &
                + pdf(Bot_,1)*pdf(Bot_,2) )
       RunFactor = RunAlphaS(NLOParam,MuRen)
       PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac * ISFac * (alpha_s4Pi*RunFactor)**4 * 0.5d0
       call sixquark((/MomExt(1:4,5),MomExt(1:4,6),-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,4),MomExt(1:4,3)/),yRnd(11:18),PreFac,1,1,.true.,res,DipoleResult)
       EvalCS_Real_ttbqqbqqb = Res / VgsWgt
!         MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!         MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!         MG_MOM(0:3,3) = MomExt(1:4,6)*100d0
!         MG_MOM(0:3,4) = MomExt(1:4,5)*100d0
!         MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!         MG_MOM(0:3,6) = MomExt(1:4,4)*100d0
!         call coupsm(0)
!         call SUU_TTBUU(MG_MOM,MadGraph_tree)
!         MadGraph_tree = MadGraph_tree*(100d0)**4


    ELSEIF( PROCESS.EQ.19 ) THEN
       PDFFac = ( pdf(AUp_,1) *pdf(AUp_,2)  + pdf(ADn_,1) *pdf(ADn_,2)   &
                + pdf(AChm_,1)*pdf(AChm_,2) + pdf(AStr_,1)*pdf(AStr_,2)  &
                + pdf(ABot_,1)*pdf(ABot_,2) )
       RunFactor = RunAlphaS(NLOParam,MuRen)
       PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac * ISFac * (alpha_s4Pi*RunFactor)**4 * 0.5d0
       call sixquark((/MomExt(1:4,5),MomExt(1:4,6),-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,4),MomExt(1:4,3)/),yRnd(11:18),PreFac,-1,-1,.true.,res,DipoleResult)
!        res = res*0.5d0
       EvalCS_Real_ttbqqbqqb = Res / VgsWgt
!         MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!         MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!         MG_MOM(0:3,3) = MomExt(1:4,6)*100d0
!         MG_MOM(0:3,4) = MomExt(1:4,5)*100d0
!         MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!         MG_MOM(0:3,6) = MomExt(1:4,4)*100d0
!         call coupsm(0)
!         call SUBUB_TTBUBUB(MG_MOM,MadGraph_tree)
!         MadGraph_tree = MadGraph_tree*(100d0)**4
    ENDIF

!    call jetktalg((/MomExt(1:4,5),MomExt(1:4,6),-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4)/),6,pjetout,Njet,weight)
!     IF( WEIGHT.EQ.0D0) RETURN  ! THIS MUST BE REMOVED AFTER CHECKING PHASE! EVERYWHERE!!!
!     DipoleResult= DipoleResult * PreFac * ISFac * (alpha_s4Pi*RunFactor)**4


!        print *, ""
!        print *, alpha_s*RunFactor,m_Top
!        print *, "Unpolarized LO result:"
!        print *, "My tree:", LO_Res_Unpol
!        print *, "MadGraph:", MadGraph_tree
!        print *, "MG  ratio: ", MadGraph_tree/dble(LO_Res_Unpol)
!        pause
!

!       s12 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,2))
!       sij = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,3))
!      print *, "E4: ",MomExt(1,4)/dsqrt(s12)
!      print *, "si12",2d0*(MomExt(1:4,1).dot.MomExt(1:4,4))/s12,2d0*(MomExt(1:4,2).dot.MomExt(1:4,4))/s12
!      print *, "R: ",get_R(MomExt(1:4,4),MomExt(1:4,3))
!       write(* ,"(1PE23.16,3X,1PE23.16,3X,1PE23.16,3X,1PE23.16)") sij/s12,EvalCS_Real_ttbqqbqqb/(DipoleResult/VgsWgt),DipoleResult
!       pause


return
END FUNCTION








FUNCTION EvalCS_DKJ_1L_ttbgg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModIntDipoles_DKJ_GGTTBG
implicit none
real(8) ::  EvalCS_DKJ_1L_ttbgg,EvalCS_DKJ_1L_ttbgg_1,EvalCS_DKJ_1L_ttbgg_2
real(8) ::  yRnd(1:VegasMxDim),VgsWgt,IOp(-2:0),HOp(1:3)
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1),FermionLoopPartAmp(7:8,-2:1)
integer :: iHel,jHel,kHel,iPrimAmp,jPrimAmp,nJetRad,GluHel
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,ISFac
real(8) :: MomExt(1:4,1:12),MomP(1:4,1:4)
logical :: applyPSCut
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,AccPoles
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2),xE
integer :: NBin(1:NumMaxHisto),NHisto,nHel(1:2)
integer :: nJetRad1=1,nJetRad3=1
integer :: nJetRad2,nJetRad4
include 'misc/global_import'
include 'vegas_common.f'

  EvalCS_DKJ_1L_ttbgg = 0d0
  EvalCS_DKJ_1L_ttbgg_1 = 0d0
  EvalCS_DKJ_1L_ttbgg_2 = 0d0

  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top ) then
      EvalCS_DKJ_1L_ttbgg = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to2(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt)
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
   ISFac = MomCrossing(MomExt)

   call SetPropagators()

   call SetPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
   RunFactor = RunAlphaS(NLOParam,MuRen)
   nHel(1:2) = getHelicity(yRnd(16))


   if( TOPDECAYS.eq.1 ) then
      nJetRad2=1;  nJetRad4=1
   elseif( TOPDECAYS.eq.2 ) then
      nJetRad2=2;  nJetRad4=2
   elseif( TOPDECAYS.eq.3 ) then
      nJetRad2=1;  nJetRad4=2
   elseif( TOPDECAYS.eq.4 ) then
      nJetRad2=2;  nJetRad4=1
   endif
!----------------------------------
!    jet emission off anti-top    |
!----------------------------------
do nJetRad=nJetRad1,nJetRad2!   nJetRad=1: jet radiation off top/bot/W, nJetRad=2: jet radiation off W/lep
   LO_Res_Unpol             = (0d0,0d0)
   NLO_Res_Unpol(-2:1)      = (0d0,0d0)
   NLO_Res_Unpol_Ferm(-2:1) = (0d0,0d0)
   if( nJetRad.eq.1 ) then
      call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:11),.true.,MomExt(1:4,5:8),PSWgt2)
   else
      call EvalPhasespace_TopDecay2(MomExt(1:4,3),yRnd(5:11),.true.,MomExt(1:4,5:8),PSWgt2)
   endif
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(12:15),.false.,MomExt(1:4,9:11),PSWgt3)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt*PSWgt2*PSWgt3 * VgsWgt * PDFFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))

   call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,9),MomExt(1:4,10),MomExt(1:4,11)/),applyPSCut,NBin)
   if( applyPSCut ) then
      cycle
   endif
!------------ LO --------------
IF( Correction.EQ.0 ) THEN
   do iHel=nHel(1),nHel(2)
   do GluHel=1,-1,-2 ! loop over additional gluon polarization
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      if( nJetRad.eq.1 ) then
          call TopDecay(ExtParticle(1),DKJ_LO_T,MomExt(1:4,5:8),GluonHel=GluHel)
      else
          call TopDecay(ExtParticle(1),DKJ_LO_Q,MomExt(1:4,5:8),GluonHel=GluHel)
      endif
      call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))

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

!------------ 1 LOOP --------------
ELSEIF( Correction.EQ.1 ) THEN
   do iHel=nHel(1),nHel(2)
   do GluHel=1,-1,-2 ! loop over additional gluon polarization
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      if( nJetRad.eq.1 ) then
         call TopDecay(ExtParticle(1),DKJ_LO_T,MomExt(1:4,5:8),GluonHel=GluHel)
      else
         call TopDecay(ExtParticle(1),DKJ_LO_Q,MomExt(1:4,5:8),GluonHel=GluHel)
      endif
      call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))

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
!    call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))

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
              if( AccPoles.gt.1d-3 ) then
                  print *, "SKIP",AccPoles
!                  call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
                  EvalCS_DKJ_1L_ttbgg = 0d0
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
ENDIF


IF( Correction.EQ.0 ) THEN
!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * WidthExpansion
   EvalCS_DKJ_1L_ttbgg = LO_Res_Unpol * PreFac

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

   EvalCS_DKJ_1L_ttbgg = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac


ELSEIF( Correction.EQ.3 ) THEN
! print *, nJetRad
! print *, "1-loop eps2:",(NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2) )* PreFac,  (NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2))/(alpha_sOver2Pi*RunFactor*LO_Res_Unpol)
! print *, "1-loop eps1:",(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1) )* PreFac,  (NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/(alpha_sOver2Pi*RunFactor*LO_Res_Unpol)
! print *, "tree virt",LO_Res_Unpol * alpha_sOver2Pi*RunFactor* PreFac

   if( PROCESS.NE.37 ) call Error("PROCESS has to be 37")
   MomP(1:4,1) = MomExt(1:4,3)
   MomP(1:4,2) = MomExt(1:4,4)
   MomP(1:4,3) =-MomExt(1:4,1)
   MomP(1:4,4) =-MomExt(1:4,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt*PSWgt2*PSWgt3 * VgsWgt
   xE = yRnd(16+HelSampling); ! xE=0.45d0
   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
   call EvalIntDipoles_DKJ_GGTTBG(MomP(1:4,1:4),MomExt(1:4,5:11),ATop_,nJetRad,xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
! print *, "HOp(1:3)",HOp(1)*pdf(0,1)*pdf(0,2),HOp(2)/xE*pdf_z(0,1)*pdf(0,2),HOp(3)/xE*pdf(0,1)*pdf_z(0,2)
! pause
   EvalCS_DKJ_1L_ttbgg = HOp(1)    * pdf(0,1)  * pdf(0,2)   &
                       + HOp(2)/xE * pdf_z(0,1)* pdf(0,2)   &
                       + HOp(3)/xE * pdf(0,1)  * pdf_z(0,2)
ENDIF

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_DKJ_1L_ttbgg)
   enddo

EvalCS_DKJ_1L_ttbgg_1 = EvalCS_DKJ_1L_ttbgg_1 + EvalCS_DKJ_1L_ttbgg

enddo! nJetRad loop


! EvalCS_DKJ_1L_ttbgg = (EvalCS_DKJ_1L_ttbgg_1)/VgsWgt
! return
! -----------------------------------------------------------------------------


!----------------------------------
!    jet emission off top         |
!----------------------------------
do nJetRad=nJetRad3,nJetRad4!   nJetRad=1: photon radiation off top/bot/W,nJetRad=2: photon radiation off W/lep
   LO_Res_Unpol             = (0d0,0d0)
   NLO_Res_Unpol(-2:1)      = (0d0,0d0)
   NLO_Res_Unpol_Ferm(-2:1) = (0d0,0d0)
   call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:8),.false.,MomExt(1:4,5:7),PSWgt2)
   if( nJetRad.eq.1 ) then
      call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(9:15),.true.,MomExt(1:4,8:11),PSWgt3)
   else
      call EvalPhasespace_TopDecay2(MomExt(1:4,4),yRnd(9:15),.true.,MomExt(1:4,8:11),PSWgt3)
   endif
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt*PSWgt2*PSWgt3 * VgsWgt * PDFFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))


   call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,10)/),applyPSCut,NBin)
   if( applyPSCut ) then
      cycle
   endif

!------------ LO --------------
IF( Correction.EQ.0 ) THEN
   do iHel=nHel(1),nHel(2)
   do GluHel=1,-1,-2 ! loop over additional gluon polarization
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,5:7))
      if( nJetRad.eq.1 ) then
         call TopDecay(ExtParticle(2),DKJ_LO_T,MomExt(1:4,8:11),GluonHel=GluHel)
      else
         call TopDecay(ExtParticle(2),DKJ_LO_Q,MomExt(1:4,8:11),GluonHel=GluHel)
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

!------------ 1 LOOP --------------
ELSEIF( Correction.EQ.1 ) THEN
   do iHel=nHel(1),nHel(2)
   do GluHel=1,-1,-2 ! loop over additional gluon polarization
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,5:7))
      if( nJetRad.eq.1 ) then
         call TopDecay(ExtParticle(2),DKJ_LO_T,MomExt(1:4,8:11),GluonHel=GluHel)
      else
         call TopDecay(ExtParticle(2),DKJ_LO_Q,MomExt(1:4,8:11),GluonHel=GluHel)
      endif

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
!    call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))

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
              if( AccPoles.gt.1d-3 ) then
                  print *, "SKIP",AccPoles
!                  call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
                  EvalCS_DKJ_1L_ttbgg = 0d0
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
ENDIF


IF( Correction.EQ.0 ) THEN
!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * WidthExpansion
   EvalCS_DKJ_1L_ttbgg = LO_Res_Unpol * PreFac

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

   EvalCS_DKJ_1L_ttbgg = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac


ELSEIF( Correction.EQ.3 ) THEN
! print *, nJetRad
! print *, "1-loop eps2:",(NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2) )* PreFac,  (NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2))/(alpha_sOver2Pi*RunFactor*LO_Res_Unpol)
! print *, "1-loop eps1:",(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1) )* PreFac,  (NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/(alpha_sOver2Pi*RunFactor*LO_Res_Unpol)
! print *, "tree virt",LO_Res_Unpol * alpha_sOver2Pi*RunFactor* PreFac


   MomP(1:4,1) = MomExt(1:4,3)
   MomP(1:4,2) = MomExt(1:4,4)
   MomP(1:4,3) =-MomExt(1:4,1)
   MomP(1:4,4) =-MomExt(1:4,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt*PSWgt2*PSWgt3 * VgsWgt

      xE = yRnd(16+HelSampling); !xE=0.055d0
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

   call EvalIntDipoles_DKJ_GGTTBG(MomP(1:4,1:4),MomExt(1:4,5:11),Top_,nJetRad,xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
! print *, "HOp(1:3)",HOp(1)*pdf(0,1)*pdf(0,2),HOp(2)/xE*pdf_z(0,1)*pdf(0,2),HOp(3)/xE*pdf(0,1)*pdf_z(0,2)
! pause
   EvalCS_DKJ_1L_ttbgg = HOp(1)    * pdf(0,1)  * pdf(0,2)   &
                       + HOp(2)/xE * pdf_z(0,1)* pdf(0,2)   &
                       + HOp(3)/xE * pdf(0,1)  * pdf_z(0,2)
ENDIF

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_DKJ_1L_ttbgg)
   enddo

   EvalCS_DKJ_1L_ttbgg_2 = EvalCS_DKJ_1L_ttbgg_2 + EvalCS_DKJ_1L_ttbgg

! print *,EvalCS_DKJ_1L_ttbgg
! pause

enddo! nJetRad loop



   EvalCS_DKJ_1L_ttbgg = (EvalCS_DKJ_1L_ttbgg_1+EvalCS_DKJ_1L_ttbgg_2)/VgsWgt
return
END FUNCTION






FUNCTION EvalCS_DKJ_1L_ttbqqb(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModIntDipoles_DKJ_QQBTTBG
use ModIntDipoles_DKJ_QGTTBQ
implicit none
real(8) :: EvalCS_DKJ_1L_ttbqqb,EvalCS_DKJ_1L_ttbqqb_1,EvalCS_DKJ_1L_ttbqqb_2
real(8) :: yRnd(1:VegasMxDim),VgsWgt,xE,HOp(1:3)
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1)
complex(8) :: BosonicPartAmp(1:2,-2:1),FermionPartAmp(1:2,-2:1)
integer :: iHel,iPrimAmp,jPrimAmp,nJetRad,GluHel
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,ISFac
real(8) :: MomExt(1:4,1:12),MomP(1:4,1:4)
logical :: applyPSCut
real(8),parameter :: Nc=3d0
real(8) :: tau,eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac_a,PDFFac_b,PDFFac,AccPoles
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2)
integer :: NHisto,NBin(1:NumMaxHisto),npdf,nHel(1:2),NRndHel
integer,parameter :: nJetRad1=1,nJetRad3=1
integer :: nJetRad2,nJetRad4
include 'misc/global_import'
include "vegas_common.f"


  EvalCS_DKJ_1L_ttbqqb = 0d0
  EvalCS_DKJ_1L_ttbqqb_1 = 0d0
  EvalCS_DKJ_1L_ttbqqb_2 = 0d0

  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top ) then
      EvalCS_DKJ_1L_ttbqqb = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to2(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt)
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))

   call setPDFs(eta1,eta2,MuFac,pdf)
   PDFFac_a = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
            + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
            + pdf(Bot_,1)*pdf(ABot_,2)
   PDFFac_b = pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
            + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
            + pdf(Bot_,2)*pdf(ABot_,1)
   RunFactor = RunAlphaS(NLOParam,MuRen)
   nHel(1:2) = getHelicity(yrnd(16))


   if( TOPDECAYS.eq.1 ) then
      nJetRad2=1;  nJetRad4=1
   elseif( TOPDECAYS.eq.2 ) then
      nJetRad2=2;  nJetRad4=2
   elseif( TOPDECAYS.eq.3 ) then
      nJetRad2=1;  nJetRad4=2
   elseif( TOPDECAYS.eq.4 ) then
      nJetRad2=2;  nJetRad4=1
   endif
!----------------------------------
! jet emission off anti-top       |
!----------------------------------
do nJetRad=nJetRad1,nJetRad2!   nJetRad=1: photon radiation off top/bot/W, nJetRad=2: photon radiation off W/lep
   LO_Res_Unpol             = (0d0,0d0)
   NLO_Res_Unpol(-2:1)      = (0d0,0d0)
   NLO_Res_Unpol_Ferm(-2:1) = (0d0,0d0)
   if( nJetRad.eq.1 ) then
      call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:11),.true.,MomExt(1:4,5:8),PSWgt2)
   else
      call EvalPhasespace_TopDecay2(MomExt(1:4,3),yRnd(5:11),.true.,MomExt(1:4,5:8),PSWgt2)
   endif
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(12:15),.false.,MomExt(1:4,9:11),PSWgt3)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt*PSWgt2*PSWgt3 * VgsWgt * dble(NumHelicities/(nHel(2)-nHel(1)+1))

   call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,9),MomExt(1:4,10),MomExt(1:4,11)/),applyPSCut,NBin)
   if( applyPSCut ) then
      cycle
   endif
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
!     call InitCurrCache()
    call SetPropagators()

    do iHel=nHel(1),nHel(2)
    do GluHel=1,-1,-2 ! loop over additional gluon polarization
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
        call SetPolarizations()
        if( nJetRad.eq.1 ) then
          call TopDecay(ExtParticle(1),DKJ_LO_T,MomExt(1:4,5:8),GluonHel=GluHel)
        else
          call TopDecay(ExtParticle(1),DKJ_LO_Q,MomExt(1:4,5:8),GluonHel=GluHel)
        endif
        call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))

        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
        enddo
        LO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
        enddo
        enddo
        LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol*PDFFac
    enddo!helicity loop
    enddo!helicity loop
  enddo! npdf loop
  call swapMom(MomExt(1:4,1),MomExt(1:4,2))   ! swap back

!------------ 1 LOOP --------------
ELSEIF( CORRECTION.EQ.1 ) THEN
  do npdf=1,2
    if(npdf.eq.1) then
        PDFFac = PDFFac_a
!         PDFFac = 1d0
    elseif(npdf.eq.2) then
        PDFFac = PDFFac_b
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
    endif
    ISFac = MomCrossing(MomExt)
!     call InitCurrCache()
    call SetPropagators()

    do iHel=nHel(1),nHel(2)
    do GluHel=1,-1,-2 ! loop over additional gluon polarization
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      if( nJetRad.eq.1 ) then
         call TopDecay(ExtParticle(1),DKJ_LO_T,MomExt(1:4,5:8),GluonHel=GluHel)
      else
         call TopDecay(ExtParticle(1),DKJ_LO_Q,MomExt(1:4,5:8),GluonHel=GluHel)
      endif
      call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))
      do iPrimAmp=1,6
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
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + Col1L_ttbqqb(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(BosonicPartAmp(jPrimAmp,-2:1)) )
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
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + Col1L_ttbqqb(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(FermionPartAmp(jPrimAmp,-2:1)) )
      enddo
      enddo
      NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) + NLO_Res_Pol(-2:1)*PDFFac
   enddo!helicity loop
   enddo!helicity loop
  enddo! npdf loop
  call swapMom(MomExt(1:4,1),MomExt(1:4,2))   ! swap back
ENDIF


IF( CORRECTION.EQ.0 ) THEN
!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * WidthExpansion
   EvalCS_DKJ_1L_ttbqqb = LO_Res_Unpol * PreFac

ELSEIF( CORRECTION.EQ.1 ) THEN
!  overall normalization: (4*Pi)^eps/Gamma(1-eps)
!  CT contributions                           ! beta           !top WFRC
   NLO_Res_UnPol(-1) = NLO_Res_UnPol(-1) + (-11d0/3d0*3d0 - 3d0*4d0/3d0)*LO_Res_Unpol
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-3d0*4d0/3d0)*2d0*dlog(MuRen/m_top)*LO_Res_Unpol  ! finite log(mu2) contrib. from  top WFRC

   NLO_Res_UnPol_Ferm(-1) = NLO_Res_UnPol_Ferm(-1) + (+2d0/3d0*Nf_light+2d0/3d0*Nf_heavy)*LO_Res_Unpol
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

   EvalCS_DKJ_1L_ttbqqb = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac

ELSEIF( CORRECTION.EQ.3 ) THEN

! print *, nJetRad
! print *, "1-loop eps2:",(NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2) )* PreFac,  (NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2))/(alpha_sOver2Pi*RunFactor*LO_Res_Unpol)
! print *, "1-loop eps1:",(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1) )* PreFac,  (NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/(alpha_sOver2Pi*RunFactor*LO_Res_Unpol)
! print *, "tree virt",LO_Res_Unpol * alpha_sOver2Pi*RunFactor* PreFac

   MomP(1:4,1) = MomExt(1:4,3)
   MomP(1:4,2) = MomExt(1:4,4)
   MomP(1:4,3) =-MomExt(1:4,1)
   MomP(1:4,4) =-MomExt(1:4,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt*PSWgt2*PSWgt3 * VgsWgt
   xE = yRnd(16+HelSampling); ! xe=0.12d0


IF( PROCESS.EQ.38 ) THEN
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

      call EvalIntDipoles_DKJ_QQBTTBG(MomP(1:4,1:4),MomExt(1:4,5:11),ATop_,nJetRad,xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_DKJ_1L_ttbqqb = HOp(1)    * (pdf(Up_,1)*pdf(AUp_,2)+pdf(Dn_,1)*pdf(ADn_,2)+pdf(Chm_,1)*pdf(AChm_,2)+pdf(Str_,1)*pdf(AStr_,2)+pdf(Bot_,1)*pdf(ABot_,2) ) &
                           + HOp(2)/xE * (pdf_z(Up_,1)*pdf(AUp_,2)+pdf_z(Dn_,1)*pdf(ADn_,2)+pdf_z(Chm_,1)*pdf(AChm_,2)+pdf_z(Str_,1)*pdf(AStr_,2)+pdf_z(Bot_,1)*pdf(ABot_,2) ) &
                           + HOp(3)/xE * (pdf(Up_,1)*pdf_z(AUp_,2)+pdf(Dn_,1)*pdf_z(ADn_,2)+pdf(Chm_,1)*pdf_z(AChm_,2)+pdf(Str_,1)*pdf_z(AStr_,2)+pdf(Bot_,1)*pdf_z(ABot_,2) )
! print *, "HOp1(1:3)",HOp(1) * (pdf(Up_,1)*pdf(AUp_,2)+pdf(Dn_,1)*pdf(ADn_,2)+pdf(Chm_,1)*pdf(AChm_,2)+pdf(Str_,1)*pdf(AStr_,2)+pdf(Bot_,1)*pdf(ABot_,2) ) &
!                     , HOp(2)/xE * (pdf_z(Up_,1)*pdf(AUp_,2)+pdf_z(Dn_,1)*pdf(ADn_,2)+pdf_z(Chm_,1)*pdf(AChm_,2)+pdf_z(Str_,1)*pdf(AStr_,2)+pdf_z(Bot_,1)*pdf(ABot_,2) ) &
!                     , HOp(3)/xE * (pdf(Up_,1)*pdf_z(AUp_,2)+pdf(Dn_,1)*pdf_z(ADn_,2)+pdf(Chm_,1)*pdf_z(AChm_,2)+pdf(Str_,1)*pdf_z(AStr_,2)+pdf(Bot_,1)*pdf_z(ABot_,2) )

      call swapMom(MomP(1:4,3),MomP(1:4,4))
      call EvalIntDipoles_DKJ_QQBTTBG(MomP(1:4,1:4),MomExt(1:4,5:11),ATop_,nJetRad,xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_DKJ_1L_ttbqqb = EvalCS_DKJ_1L_ttbqqb &
                           + HOp(1)    * (pdf(Up_,2)*pdf(AUp_,1)+pdf(Dn_,2)*pdf(ADn_,1)+pdf(Chm_,2)*pdf(AChm_,1)+pdf(Str_,2)*pdf(AStr_,1)+pdf(Bot_,2)*pdf(ABot_,1) ) &
                           + HOp(2)/xE * (pdf_z(Up_,2)*pdf(AUp_,1)+pdf_z(Dn_,2)*pdf(ADn_,1)+pdf_z(Chm_,2)*pdf(AChm_,1)+pdf_z(Str_,2)*pdf(AStr_,1)+pdf_z(Bot_,2)*pdf(ABot_,1) ) &
                           + HOp(3)/xE * (pdf(Up_,2)*pdf_z(AUp_,1)+pdf(Dn_,2)*pdf_z(ADn_,1)+pdf(Chm_,2)*pdf_z(AChm_,1)+pdf(Str_,2)*pdf_z(AStr_,1)+pdf(Bot_,2)*pdf_z(ABot_,1) )
      call swapMom(MomP(1:4,3),MomP(1:4,4))! swap back
! print *, "HOp2(1:3)",HOp(1)    * (pdf(Up_,2)*pdf(AUp_,1)+pdf(Dn_,2)*pdf(ADn_,1)+pdf(Chm_,2)*pdf(AChm_,1)+pdf(Str_,2)*pdf(AStr_,1)+pdf(Bot_,2)*pdf(ABot_,1) ) &
!                    , HOp(2)/xE * (pdf_z(Up_,2)*pdf(AUp_,1)+pdf_z(Dn_,2)*pdf(ADn_,1)+pdf_z(Chm_,2)*pdf(AChm_,1)+pdf_z(Str_,2)*pdf(AStr_,1)+pdf_z(Bot_,2)*pdf(ABot_,1) ) &
!                    , HOp(3)/xE * (pdf(Up_,2)*pdf_z(AUp_,1)+pdf(Dn_,2)*pdf_z(ADn_,1)+pdf(Chm_,2)*pdf_z(AChm_,1)+pdf(Str_,2)*pdf_z(AStr_,1)+pdf(Bot_,2)*pdf_z(ABot_,1) )
! pause


ELSEIF( PROCESS.EQ.35 ) THEN
   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

   call EvalIntDipoles_DKJ_QGTTBQ(MomP(1:4,1:4),MomExt(1:4,5:11),ATop_,nJetRad,xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
   EvalCS_DKJ_1L_ttbqqb = HOp(1) * (pdf(Up_,1)*pdf(0,2)+pdf(Dn_,1)*pdf(0,2)+pdf(Chm_,1)*pdf(0,2)+pdf(Str_,1)*pdf(0,2)+pdf(Bot_,1)*pdf(0,2) ) &
                       + HOp(2)/xE * (pdf_z(Up_,1)*pdf(0,2)+pdf_z(Dn_,1)*pdf(0,2)+pdf_z(Chm_,1)*pdf(0,2)+pdf_z(Str_,1)*pdf(0,2)+pdf_z(Bot_,1)*pdf(0,2) ) &
                       + HOp(3)/xE * (pdf(Up_,1)*pdf_z(0,2)+pdf(Dn_,1)*pdf_z(0,2)+pdf(Chm_,1)*pdf_z(0,2)+pdf(Str_,1)*pdf_z(0,2)+pdf(Bot_,1)*pdf_z(0,2) )

   call swapMom(MomP(1:4,3),MomP(1:4,4))
   call EvalIntDipoles_DKJ_QGTTBQ(MomP(1:4,1:4),MomExt(1:4,5:11),ATop_,nJetRad,xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
   EvalCS_DKJ_1L_ttbqqb = EvalCS_DKJ_1L_ttbqqb &
                       + HOp(1)    * (pdf(Up_,2)*pdf(0,1)+pdf(Dn_,2)*pdf(0,1)+pdf(Chm_,2)*pdf(0,1)+pdf(Str_,2)*pdf(0,1)+pdf(Bot_,2)*pdf(0,1) ) &
                       + HOp(2)/xE * (pdf_z(Up_,2)*pdf(0,1)+pdf_z(Dn_,2)*pdf(0,1)+pdf_z(Chm_,2)*pdf(0,1)+pdf_z(Str_,2)*pdf(0,1)+pdf_z(Bot_,2)*pdf(0,1) ) &
                       + HOp(3)/xE * (pdf(Up_,2)*pdf_z(0,1)+pdf(Dn_,2)*pdf_z(0,1)+pdf(Chm_,2)*pdf_z(0,1)+pdf(Str_,2)*pdf_z(0,1)+pdf(Bot_,2)*pdf_z(0,1) )
   call swapMom(MomP(1:4,3),MomP(1:4,4))! swap back




ELSEIF( PROCESS.EQ.36 ) THEN
   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

   call EvalIntDipoles_DKJ_QGTTBQ(MomP(1:4,1:4),MomExt(1:4,5:11),ATop_,nJetRad,xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
   EvalCS_DKJ_1L_ttbqqb = HOp(1) * (pdf(AUp_,1)*pdf(0,2)+pdf(ADn_,1)*pdf(0,2)+pdf(AChm_,1)*pdf(0,2)+pdf(AStr_,1)*pdf(0,2)+pdf(ABot_,1)*pdf(0,2) ) &
                       + HOp(2)/xE * (pdf_z(AUp_,1)*pdf(0,2)+pdf_z(ADn_,1)*pdf(0,2)+pdf_z(AChm_,1)*pdf(0,2)+pdf_z(AStr_,1)*pdf(0,2)+pdf_z(ABot_,1)*pdf(0,2) ) &
                       + HOp(3)/xE * (pdf(AUp_,1)*pdf_z(0,2)+pdf(ADn_,1)*pdf_z(0,2)+pdf(AChm_,1)*pdf_z(0,2)+pdf(AStr_,1)*pdf_z(0,2)+pdf(ABot_,1)*pdf_z(0,2) )

   call swapMom(MomP(1:4,3),MomP(1:4,4))
   call EvalIntDipoles_DKJ_QGTTBQ(MomP(1:4,1:4),MomExt(1:4,5:11),ATop_,nJetRad,xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
   EvalCS_DKJ_1L_ttbqqb = EvalCS_DKJ_1L_ttbqqb &
                       + HOp(1)    * (pdf(AUp_,2)*pdf(0,1)+pdf(ADn_,2)*pdf(0,1)+pdf(AChm_,2)*pdf(0,1)+pdf(AStr_,2)*pdf(0,1)+pdf(ABot_,2)*pdf(0,1) ) &
                       + HOp(2)/xE * (pdf_z(AUp_,2)*pdf(0,1)+pdf_z(ADn_,2)*pdf(0,1)+pdf_z(AChm_,2)*pdf(0,1)+pdf_z(AStr_,2)*pdf(0,1)+pdf_z(ABot_,2)*pdf(0,1) ) &
                       + HOp(3)/xE * (pdf(AUp_,2)*pdf_z(0,1)+pdf(ADn_,2)*pdf_z(0,1)+pdf(AChm_,2)*pdf_z(0,1)+pdf(AStr_,2)*pdf_z(0,1)+pdf(ABot_,2)*pdf_z(0,1) )
   call swapMom(MomP(1:4,3),MomP(1:4,4))! swap back

ENDIF
ENDIF

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_DKJ_1L_ttbqqb)
   enddo

EvalCS_DKJ_1L_ttbqqb_1 = EvalCS_DKJ_1L_ttbqqb_1 + EvalCS_DKJ_1L_ttbqqb
enddo! nJetRad loop


!    EvalCS_DKJ_1L_ttbqqb = (EvalCS_DKJ_1L_ttbqqb_1)/VgsWgt
!    return
!-----------------------------------------------------------------------------









!----------------------------------
!         jet emission off top    |
!----------------------------------
do nJetRad=nJetRad3,nJetRad4!   nJetRad=1: photon radiation off top/bot/W, nJetRad=2: photon radiation off W/lep
   LO_Res_Unpol             = (0d0,0d0)
   NLO_Res_Unpol(-2:1)      = (0d0,0d0)
   NLO_Res_Unpol_Ferm(-2:1) = (0d0,0d0)
   call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:8),.false.,MomExt(1:4,5:7),PSWgt2)
   if( nJetRad.eq.1 ) then
      call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(9:15),.true.,MomExt(1:4,8:11),PSWgt3)
   else
      call EvalPhasespace_TopDecay2(MomExt(1:4,4),yRnd(9:15),.true.,MomExt(1:4,8:11),PSWgt3)
   endif
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt*PSWgt2*PSWgt3 * VgsWgt * dble(NumHelicities/(nHel(2)-nHel(1)+1))

   call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,10)/),applyPSCut,NBin)
   if( applyPSCut ) then
      cycle
   endif

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
    call SetPropagators()

    do iHel=nHel(1),nHel(2)
    do GluHel=1,-1,-2 ! loop over additional gluon polarization
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
        call SetPolarizations()
        call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,5:7))
        if( nJetRad.eq.1 ) then
          call TopDecay(ExtParticle(2),DKJ_LO_T,MomExt(1:4,8:11),GluonHel=GluHel)
        else
          call TopDecay(ExtParticle(2),DKJ_LO_Q,MomExt(1:4,8:11),GluonHel=GluHel)
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
        LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol*PDFFac
    enddo!helicity loop
    enddo!helicity loop
  enddo! npdf loop
  call swapMom(MomExt(1:4,1),MomExt(1:4,2))!  swap back

!------------ 1 LOOP --------------
ELSEIF( CORRECTION.EQ.1 ) THEN
  do npdf=1,2
    if(npdf.eq.1) then
        PDFFac = PDFFac_a
!         PDFFac = 1d0
    elseif(npdf.eq.2) then
        PDFFac = PDFFac_b
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
    endif
    ISFac = MomCrossing(MomExt)
    call SetPropagators()

    do iHel=nHel(1),nHel(2)
    do GluHel=1,-1,-2 ! loop over additional gluon polarization
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,5:7))
      if( nJetRad.eq.1 ) then
          call TopDecay(ExtParticle(2),DKJ_LO_T,MomExt(1:4,8:11),GluonHel=GluHel)
      else
          call TopDecay(ExtParticle(2),DKJ_LO_Q,MomExt(1:4,8:11),GluonHel=GluHel)
      endif

      do iPrimAmp=1,6
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
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + Col1L_ttbqqb(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(BosonicPartAmp(jPrimAmp,-2:1)) )
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
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + Col1L_ttbqqb(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(FermionPartAmp(jPrimAmp,-2:1)) )
      enddo
      enddo
      NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) + NLO_Res_Pol(-2:1)*PDFFac
   enddo!helicity loop
   enddo!helicity loop
  enddo! npdf loop
  call swapMom(MomExt(1:4,1),MomExt(1:4,2))   ! swap back
ENDIF


IF( CORRECTION.EQ.0 ) THEN
!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * WidthExpansion
   EvalCS_DKJ_1L_ttbqqb = LO_Res_Unpol * PreFac

ELSEIF( CORRECTION.EQ.1 ) THEN
!  overall normalization: (4*Pi)^eps/Gamma(1-eps)
!  CT contributions                           ! beta           !top WFRC
   NLO_Res_UnPol(-1) = NLO_Res_UnPol(-1) + (-11d0/3d0*3d0 - 3d0*4d0/3d0)*LO_Res_Unpol
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-3d0*4d0/3d0)*2d0*dlog(MuRen/m_top)*LO_Res_Unpol  ! finite log(mu2) contrib. from  top WFRC

   NLO_Res_UnPol_Ferm(-1) = NLO_Res_UnPol_Ferm(-1) + (+2d0/3d0*Nf_light+2d0/3d0*Nf_heavy)*LO_Res_Unpol
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

   EvalCS_DKJ_1L_ttbqqb = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac


ELSEIF( CORRECTION.EQ.3 ) THEN
! print *, nJetRad
! print *, "1-loop eps2:",(NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2) )* PreFac,  (NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2))/(alpha_sOver2Pi*RunFactor*LO_Res_Unpol)
! print *, "1-loop eps1:",(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1) )* PreFac,  (NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/(alpha_sOver2Pi*RunFactor*LO_Res_Unpol)
! print *, "tree virt",LO_Res_Unpol * alpha_sOver2Pi*RunFactor* PreFac


   MomP(1:4,1) = MomExt(1:4,3)
   MomP(1:4,2) = MomExt(1:4,4)
   MomP(1:4,3) =-MomExt(1:4,1)
   MomP(1:4,4) =-MomExt(1:4,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt*PSWgt2*PSWgt3 * VgsWgt
   xE = yRnd(16+HelSampling); ! xE=0.123d0



IF( PROCESS.EQ.38 ) THEN
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

      call EvalIntDipoles_DKJ_QQBTTBG(MomP(1:4,1:4),MomExt(1:4,5:11),Top_,nJetRad,xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_DKJ_1L_ttbqqb = HOp(1)    * (pdf(Up_,1)*pdf(AUp_,2)+pdf(Dn_,1)*pdf(ADn_,2)+pdf(Chm_,1)*pdf(AChm_,2)+pdf(Str_,1)*pdf(AStr_,2)+pdf(Bot_,1)*pdf(ABot_,2) ) &
                           + HOp(2)/xE * (pdf_z(Up_,1)*pdf(AUp_,2)+pdf_z(Dn_,1)*pdf(ADn_,2)+pdf_z(Chm_,1)*pdf(AChm_,2)+pdf_z(Str_,1)*pdf(AStr_,2)+pdf_z(Bot_,1)*pdf(ABot_,2) ) &
                           + HOp(3)/xE * (pdf(Up_,1)*pdf_z(AUp_,2)+pdf(Dn_,1)*pdf_z(ADn_,2)+pdf(Chm_,1)*pdf_z(AChm_,2)+pdf(Str_,1)*pdf_z(AStr_,2)+pdf(Bot_,1)*pdf_z(ABot_,2) )
! print *, "HOp1(1:3)",HOp(1) * (pdf(Up_,1)*pdf(AUp_,2)+pdf(Dn_,1)*pdf(ADn_,2)+pdf(Chm_,1)*pdf(AChm_,2)+pdf(Str_,1)*pdf(AStr_,2)+pdf(Bot_,1)*pdf(ABot_,2) ) &
!                     , HOp(2)/xE * (pdf_z(Up_,1)*pdf(AUp_,2)+pdf_z(Dn_,1)*pdf(ADn_,2)+pdf_z(Chm_,1)*pdf(AChm_,2)+pdf_z(Str_,1)*pdf(AStr_,2)+pdf_z(Bot_,1)*pdf(ABot_,2) ) &
!                     , HOp(3)/xE * (pdf(Up_,1)*pdf_z(AUp_,2)+pdf(Dn_,1)*pdf_z(ADn_,2)+pdf(Chm_,1)*pdf_z(AChm_,2)+pdf(Str_,1)*pdf_z(AStr_,2)+pdf(Bot_,1)*pdf_z(ABot_,2) )


      call swapMom(MomP(1:4,3),MomP(1:4,4))
      call EvalIntDipoles_DKJ_QQBTTBG(MomP(1:4,1:4),MomExt(1:4,5:11),Top_,nJetRad,xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_DKJ_1L_ttbqqb = EvalCS_DKJ_1L_ttbqqb &
                           + HOp(1)    * (pdf(Up_,2)*pdf(AUp_,1)+pdf(Dn_,2)*pdf(ADn_,1)+pdf(Chm_,2)*pdf(AChm_,1)+pdf(Str_,2)*pdf(AStr_,1)+pdf(Bot_,2)*pdf(ABot_,1) ) &
                           + HOp(2)/xE * (pdf_z(Up_,2)*pdf(AUp_,1)+pdf_z(Dn_,2)*pdf(ADn_,1)+pdf_z(Chm_,2)*pdf(AChm_,1)+pdf_z(Str_,2)*pdf(AStr_,1)+pdf_z(Bot_,2)*pdf(ABot_,1) ) &
                           + HOp(3)/xE * (pdf(Up_,2)*pdf_z(AUp_,1)+pdf(Dn_,2)*pdf_z(ADn_,1)+pdf(Chm_,2)*pdf_z(AChm_,1)+pdf(Str_,2)*pdf_z(AStr_,1)+pdf(Bot_,2)*pdf_z(ABot_,1) )
      call swapMom(MomP(1:4,3),MomP(1:4,4))! swap back
! print *, "HOp2(1:3)",HOp(1)    * (pdf(Up_,2)*pdf(AUp_,1)+pdf(Dn_,2)*pdf(ADn_,1)+pdf(Chm_,2)*pdf(AChm_,1)+pdf(Str_,2)*pdf(AStr_,1)+pdf(Bot_,2)*pdf(ABot_,1) ) &
!                    , HOp(2)/xE * (pdf_z(Up_,2)*pdf(AUp_,1)+pdf_z(Dn_,2)*pdf(ADn_,1)+pdf_z(Chm_,2)*pdf(AChm_,1)+pdf_z(Str_,2)*pdf(AStr_,1)+pdf_z(Bot_,2)*pdf(ABot_,1) ) &
!                    , HOp(3)/xE * (pdf(Up_,2)*pdf_z(AUp_,1)+pdf(Dn_,2)*pdf_z(ADn_,1)+pdf(Chm_,2)*pdf_z(AChm_,1)+pdf(Str_,2)*pdf_z(AStr_,1)+pdf(Bot_,2)*pdf_z(ABot_,1) )
! pause


ELSEIF( PROCESS.EQ.35 ) THEN
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

      call EvalIntDipoles_DKJ_QGTTBQ(MomP(1:4,1:4),MomExt(1:4,5:11),Top_,nJetRad,xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_DKJ_1L_ttbqqb = HOp(1) * (pdf(Up_,1)*pdf(0,2)+pdf(Dn_,1)*pdf(0,2)+pdf(Chm_,1)*pdf(0,2)+pdf(Str_,1)*pdf(0,2)+pdf(Bot_,1)*pdf(0,2) ) &
                          + HOp(2)/xE * (pdf_z(Up_,1)*pdf(0,2)+pdf_z(Dn_,1)*pdf(0,2)+pdf_z(Chm_,1)*pdf(0,2)+pdf_z(Str_,1)*pdf(0,2)+pdf_z(Bot_,1)*pdf(0,2) ) &
                          + HOp(3)/xE * (pdf(Up_,1)*pdf_z(0,2)+pdf(Dn_,1)*pdf_z(0,2)+pdf(Chm_,1)*pdf_z(0,2)+pdf(Str_,1)*pdf_z(0,2)+pdf(Bot_,1)*pdf_z(0,2) )

      call swapMom(MomP(1:4,3),MomP(1:4,4))
      call EvalIntDipoles_DKJ_QGTTBQ(MomP(1:4,1:4),MomExt(1:4,5:11),Top_,nJetRad,xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_DKJ_1L_ttbqqb = EvalCS_DKJ_1L_ttbqqb &
                          + HOp(1)    * (pdf(Up_,2)*pdf(0,1)+pdf(Dn_,2)*pdf(0,1)+pdf(Chm_,2)*pdf(0,1)+pdf(Str_,2)*pdf(0,1)+pdf(Bot_,2)*pdf(0,1) ) &
                          + HOp(2)/xE * (pdf_z(Up_,2)*pdf(0,1)+pdf_z(Dn_,2)*pdf(0,1)+pdf_z(Chm_,2)*pdf(0,1)+pdf_z(Str_,2)*pdf(0,1)+pdf_z(Bot_,2)*pdf(0,1) ) &
                          + HOp(3)/xE * (pdf(Up_,2)*pdf_z(0,1)+pdf(Dn_,2)*pdf_z(0,1)+pdf(Chm_,2)*pdf_z(0,1)+pdf(Str_,2)*pdf_z(0,1)+pdf(Bot_,2)*pdf_z(0,1) )
      call swapMom(MomP(1:4,3),MomP(1:4,4))! swap back

ELSEIF( PROCESS.EQ.36 ) THEN
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

      call EvalIntDipoles_DKJ_QGTTBQ(MomP(1:4,1:4),MomExt(1:4,5:11),Top_,nJetRad,xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_DKJ_1L_ttbqqb = HOp(1) * (pdf(AUp_,1)*pdf(0,2)+pdf(ADn_,1)*pdf(0,2)+pdf(AChm_,1)*pdf(0,2)+pdf(AStr_,1)*pdf(0,2)+pdf(ABot_,1)*pdf(0,2) ) &
                          + HOp(2)/xE * (pdf_z(AUp_,1)*pdf(0,2)+pdf_z(ADn_,1)*pdf(0,2)+pdf_z(AChm_,1)*pdf(0,2)+pdf_z(AStr_,1)*pdf(0,2)+pdf_z(ABot_,1)*pdf(0,2) ) &
                          + HOp(3)/xE * (pdf(AUp_,1)*pdf_z(0,2)+pdf(ADn_,1)*pdf_z(0,2)+pdf(AChm_,1)*pdf_z(0,2)+pdf(AStr_,1)*pdf_z(0,2)+pdf(ABot_,1)*pdf_z(0,2) )

      call swapMom(MomP(1:4,3),MomP(1:4,4))
      call EvalIntDipoles_DKJ_QGTTBQ(MomP(1:4,1:4),MomExt(1:4,5:11),Top_,nJetRad,xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_DKJ_1L_ttbqqb = EvalCS_DKJ_1L_ttbqqb &
                          + HOp(1)    * (pdf(AUp_,2)*pdf(0,1)+pdf(ADn_,2)*pdf(0,1)+pdf(AChm_,2)*pdf(0,1)+pdf(AStr_,2)*pdf(0,1)+pdf(ABot_,2)*pdf(0,1) ) &
                          + HOp(2)/xE * (pdf_z(AUp_,2)*pdf(0,1)+pdf_z(ADn_,2)*pdf(0,1)+pdf_z(AChm_,2)*pdf(0,1)+pdf_z(AStr_,2)*pdf(0,1)+pdf_z(ABot_,2)*pdf(0,1) ) &
                          + HOp(3)/xE * (pdf(AUp_,2)*pdf_z(0,1)+pdf(ADn_,2)*pdf_z(0,1)+pdf(AChm_,2)*pdf_z(0,1)+pdf(AStr_,2)*pdf_z(0,1)+pdf(ABot_,2)*pdf_z(0,1) )
      call swapMom(MomP(1:4,3),MomP(1:4,4))! swap back

ENDIF
ENDIF

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_DKJ_1L_ttbqqb)
   enddo

   EvalCS_DKJ_1L_ttbqqb_2 = EvalCS_DKJ_1L_ttbqqb_2 + EvalCS_DKJ_1L_ttbqqb
enddo! nJetRad loop

   EvalCS_DKJ_1L_ttbqqb = (EvalCS_DKJ_1L_ttbqqb_1+EvalCS_DKJ_1L_ttbqqb_2)/VgsWgt

return
END FUNCTION





FUNCTION EvalCS_1LDK_ttbggg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
implicit none
real(8) ::  EvalCS_1LDK_ttbggg,yRnd(1:VegasMxDim),VgsWgt
complex(8) :: NLO_Res_Pol,NLO_Res_Unpol,TreeResult(1:6),VirtResult(1:6)
integer :: iHel,iPrimAmp,jPrimAmp,nHel(1:2)
integer :: NBin(1:NumMaxHisto),NHisto
real(8) :: SpinAvg,ColorAvg,EHat,PSWgt,PSWgt2,PSWgt3,ISFac,RunFactor
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac
real(8) :: MomExt(1:4,1:5),MomDK(1:4,1:7)
real(8) :: pdf(-6:6,1:2)
logical :: applyPSCut
include 'vegas_common.f'

   EvalCS_1LDK_ttbggg = 0d0

   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if( EHat.le.2d0*m_Top ) then
      EvalCS_1LDK_ttbggg = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)
   call EvalPhaseSpace_2to3(EHat,yRnd(3:7),MomExt(1:4,1:5),PSWgt)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))
   ISFac = MomCrossing(MomExt)

   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomDK(1:4,1:3),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:15),.false.,MomDK(1:4,4:6),PSWgt3)
   PSWgt = PSWgt * PSWgt2*PSWgt3
   call Kinematics_TTBARJET(0,MomExt,MomDK,applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_1LDK_ttbggg = 0d0
      return
   endif

   ISFac = MomCrossing(MomExt)
   call setPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   RunFactor = RunAlphaS(NLOParam,MuRen)


!----------------------------------------
! one loop correction on anti-top decay |
!----------------------------------------
    NLO_Res_Unpol = (0d0,0d0)
    nHel(1:2) = getHelicity(yrnd(16))
    PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))
    do iHel=nHel(1),nHel(2)
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
        call SetPolarizations()
        call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
        call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))

        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call TopDecay(ExtParticle(1),DK_1L_T,MomDK(1:4,1:3))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbggg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol


        if( TOPDECAYS.eq.2 .or. TOPDECAYS.eq.4 ) then
              call TopDecay(ExtParticle(1),DK_1L_Q,MomDK(1:4,1:3))
              do iPrimAmp=1,NumBornAmps
                  call EvalTree(BornAmps(iPrimAmp))
                  VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
              enddo

              NLO_Res_Pol = (0d0,0d0)
              do jPrimAmp=1,NumBornAmps
              do iPrimAmp=1,NumBornAmps
                  NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbggg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
              enddo
              enddo
              NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
        endif

    enddo! iHel loop
    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**3 * PreFac
    EvalCS_1LDK_ttbggg = NLO_Res_Unpol





!-------------------------------------
! one loop correction to top-decay   |
!-------------------------------------
    NLO_Res_Unpol = (0d0,0d0)
    nHel(1:2) = getHelicity(yrnd(16))
    PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))
    do iHel=nHel(1),nHel(2)
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
        call SetPolarizations()
        call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
        call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))

        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call TopDecay(ExtParticle(2),DK_1L_T,MomDK(1:4,4:6))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbggg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol

        if( TOPDECAYS.eq.2 .or. TOPDECAYS.eq.3 ) then
              call TopDecay(ExtParticle(2),DK_1L_Q,MomDK(1:4,4:6))
              do iPrimAmp=1,NumBornAmps
                  call EvalTree(BornAmps(iPrimAmp))
                  VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
              enddo

              NLO_Res_Pol = (0d0,0d0)
              do jPrimAmp=1,NumBornAmps
              do iPrimAmp=1,NumBornAmps
                  NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbggg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
              enddo
              enddo
              NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
        endif
    enddo

   NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**3 * PreFac
   EvalCS_1LDK_ttbggg = EvalCS_1LDK_ttbggg + NLO_Res_Unpol


   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_1LDK_ttbggg)
   enddo


   EvalCS_1LDK_ttbggg = EvalCS_1LDK_ttbggg/VgsWgt
   if( IsNan(EvalCS_1LDK_ttbggg) ) then
        print *, "NAN:",EvalCS_1LDK_ttbggg
        call printYrnd(yrnd(:))
        EvalCS_1LDK_ttbggg = 0d0
        return
   endif

return
END FUNCTION









FUNCTION EvalCS_1LDK_ttbqqbg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes
use ModParameters
use ModMisc
implicit none
real(8) ::  EvalCS_1LDK_ttbqqbg,yRnd(1:VegasMxDim),VgsWgt
complex(8) :: NLO_Res_Pol,NLO_Res_Unpol,TreeResult(1:6),VirtResult(1:6)
integer :: iHel,iPrimAmp,jPrimAmp,nHel(1:2)
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,ISFac
real(8) :: MomExt(1:4,1:NumExtParticles),MomDK(1:4,1:7)
logical :: applyPSCut
real(8),parameter :: Nc=3d0
real(8) :: tau,eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,PDFFac_a,PDFFac_b
real(8) :: pdf(-6:6,1:2)
integer :: NBin(1:NumMaxHisto),NHisto,npdf
include 'vegas_common.f'


  EvalCS_1LDK_ttbqqbg = 0d0

  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top ) then
      EvalCS_1LDK_ttbqqbg = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to3(EHat,yRnd(3:7),MomExt(1:4,1:5),PSWgt)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))

   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomDK(1:4,1:3),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:15),.false.,MomDK(1:4,4:6),PSWgt3)
   PSWgt = PSWgt * PSWgt2*PSWgt3

   call Kinematics_TTBARJET(0,MomExt,MomDK,applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_1LDK_ttbqqbg = 0d0
      PSCutCounter = PSCutCounter + 1
      return
   endif

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
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   RunFactor = RunAlphaS(NLOParam,MuRen)


!----------------------------------------
! one loop correction on anti-top decay |
!----------------------------------------
   NLO_Res_Unpol = (0d0,0d0)
   nHel(1:2) = getHelicity(yrnd(16))
   PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))
DO NPDF=1,2
   if(npdf.eq.1) then
        PDFFac = PDFFac_a
   elseif(npdf.eq.2) then
        PDFFac = PDFFac_b
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
   endif

   ISFac = MomCrossing(MomExt)
   do iHel=nHel(1),nHel(2)
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
        call SetPolarizations()
        call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
        call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))

        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call TopDecay(ExtParticle(1),DK_1L_T,MomDK(1:4,1:3))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbqqbg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol *PDFFac

        if( TOPDECAYS.eq.2 .or. TOPDECAYS.eq.4 ) then
              call TopDecay(ExtParticle(1),DK_1L_Q,MomDK(1:4,1:3))
              do iPrimAmp=1,NumBornAmps
                  call EvalTree(BornAmps(iPrimAmp))
                  VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
              enddo

              NLO_Res_Pol = (0d0,0d0)
              do jPrimAmp=1,NumBornAmps
              do iPrimAmp=1,NumBornAmps
                  NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbqqbg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
              enddo
              enddo
              NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol *PDFFac
        endif

    enddo! helicity loop
ENDDO! loop over a<-->b pdfs
call swapMom(MomExt(1:4,1),MomExt(1:4,2))   ! swap back to original order

    NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**3 * PreFac
    EvalCS_1LDK_ttbqqbg = NLO_Res_Unpol






!----------------------------------------
! one loop correction on top decay      |
!----------------------------------------
   NLO_Res_Unpol = (0d0,0d0)
   nHel(1:2) = getHelicity(yrnd(16))
   PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))
DO NPDF=1,2
   if(npdf.eq.1) then
        PDFFac = PDFFac_a
   elseif(npdf.eq.2) then
        PDFFac = PDFFac_b
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
   endif

   ISFac = MomCrossing(MomExt)
   do iHel=nHel(1),nHel(2)
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
        call SetPolarizations()
        call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
        call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))

        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        call TopDecay(ExtParticle(2),DK_1L_T,MomDK(1:4,4:6))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo

        NLO_Res_Pol = (0d0,0d0)
        do jPrimAmp=1,NumBornAmps
        do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbqqbg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
        enddo
        enddo
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol *PDFFac

        if( TOPDECAYS.eq.2 .or. TOPDECAYS.eq.3 ) then
              call TopDecay(ExtParticle(2),DK_1L_Q,MomDK(1:4,4:6))
              do iPrimAmp=1,NumBornAmps
                  call EvalTree(BornAmps(iPrimAmp))
                  VirtResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
              enddo

              NLO_Res_Pol = (0d0,0d0)
              do jPrimAmp=1,NumBornAmps
              do iPrimAmp=1,NumBornAmps
                  NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbqqbg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(VirtResult(jPrimAmp)) )
              enddo
              enddo
              NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol *PDFFac
        endif

    enddo! helicity loop
ENDDO! loop over a<-->b pdfs
! call swapMom(MomExt(1:4,1),MomExt(1:4,2))   ! swap back to original order

   NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**3 * PreFac
   EvalCS_1LDK_ttbqqbg = EvalCS_1LDK_ttbqqbg + NLO_Res_Unpol


   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_1LDK_ttbqqbg)
   enddo




   EvalCS_1LDK_ttbqqbg = EvalCS_1LDK_ttbqqbg/VgsWgt
   if( IsNan(EvalCS_1LDK_ttbqqbg) ) then
        print *, "NAN:",EvalCS_1LDK_ttbqqbg
        call printYrnd(yrnd(:))
        EvalCS_1LDK_ttbqqbg = 0d0
        return
   endif

return
END FUNCTION






! something very funny happens here: sometimes there are NaN events when compiled with O2
! everything works fine with O1, compiler is ifort V.10.1.  l_fc_p_10.1.018 (64bit)
! it is probably only the dipole module
FUNCTION EvalCS_DKJ_Real_ttbggg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModDipoles_DKJ_GGTTBG
use ModAmplitudes
use ModParameters
use ModMisc
use ModHadrWDecay
implicit none
real(8) :: EvalCS_DKJ_Real_ttbggg,EvalCS_DKJ_Real_ttbggg_1,EvalCS_DKJ_Real_ttbggg_2
real(8) :: yRnd(1:VegasMxDim),VgsWgt
complex(8) :: LO_Res_Pol,LO_Res_Unpol
integer :: iHel,iPrimAmp,jPrimAmp,nJetRad,GluHel,nDip
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,DipoleResult,ISFac
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,6:11),dummy(1:5)
real(8) :: pbDpg,ptDpg,ptDpb,omz,rsq,z,y,TheDipole
logical :: applyPSCut,applySingCut
real(8) :: tau,eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac
real(8) :: pdf(-6:6,1:2)
integer :: NBin(1:NumMaxHisto),NHisto
real(8) :: CF=4d0/3d0
integer,parameter :: nJetRad1=1,nJetRad3=1
integer :: nJetRad2,nJetRad4
include "vegas_common.f"


dummy(:)=0d0
! print *, "yrnd fixed"
!  yrnd(1)=  0.440979372449676d0
!  yrnd(2)= 0.322643898112068d0
!  yrnd(3)= 0.203127593827959d0
!  yrnd(4)=   4.031276169247592d-002
!  yrnd(5)=   0.396959786255360d0
!  yrnd(6)=   0.266077416607215d0
!  yrnd(7)=   0.311592539917488d0
!  yrnd(8)=    9.021701574801334d-002
!  yrnd(9)=     0.119413255070994d0
!  yrnd(10)=   0.460079616615586d0
!  yrnd(11)=   0.250086936051998d0
!  yrnd(12)=   0.414713387105015d0
!  yrnd(13)=   0.435818393219178d0
!  yrnd(14)=   0.237879456597324d0
!  yrnd(15)=   0.465404917935564d0
!  yrnd(16)=   2.032386908322755d-002
!  yrnd(17)=    0.532935153242637d0
!  yrnd(18)=     0.772143932186135d0


  EvalCS_DKJ_Real_ttbggg = 0d0
  EvalCS_DKJ_Real_ttbggg_1 = 0d0
  EvalCS_DKJ_Real_ttbggg_2 = 0d0

  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top ) then
      EvalCS_DKJ_Real_ttbggg = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to3(EHat,yRnd(3:7),MomExt(1:4,1:5),PSWgt)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))
   ISFac = MomCrossing(MomExt)


   call CheckSing(MomExt,applySingCut)
   if( applySingCut ) then
       EvalCS_DKJ_Real_ttbggg = 0d0
       SkipCounter = SkipCounter + 1
       return
   endif

   call setPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
   RunFactor = RunAlphaS(2,MuRen)



   if( TOPDECAYS.eq.1 ) then
      nJetRad2=1;  nJetRad4=1
   elseif( TOPDECAYS.eq.2 ) then
      nJetRad2=2;  nJetRad4=2
   elseif( TOPDECAYS.eq.3 ) then
      nJetRad2=1;  nJetRad4=2
   elseif( TOPDECAYS.eq.4 ) then
      nJetRad2=2;  nJetRad4=1
   endif
!----------------------------------
!  gluon  emission off anti-top   |
!----------------------------------
do nJetRad=nJetRad1,nJetRad2!   nJetRad=1: gluon radiation off top line,nJetRad=2: gluon radiation off W line
   EvalCS_DKJ_Real_ttbggg = 0d0
   if( nJetRad.eq.1 ) then
      call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:14),.true.,MomExt(1:4,6:9),PSWgt2)
   else
      call EvalPhasespace_TopDecay2(MomExt(1:4,4),yRnd(8:14),.true.,MomExt(1:4,6:9),PSWgt2)
   endif
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(15:18),.false.,MomExt(1:4,10:12),PSWgt3)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt*PSWgt2*PSWgt3 * VgsWgt * PDFFac
   if( PreFac.eq.0d0 ) cycle

   call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,9),MomExt(1:4,4),MomExt(1:4,5)/),(/MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_DKJ_Real_ttbggg = 0d0
       PSCutCounter = PSCutCounter + 1
      goto 50
   endif

        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
        do GluHel=1,-1,-2 ! loop over additional gluon polarization
          call HelCrossing(Helicities(iHel,1:NumExtParticles))
          call SetPolarizations()
          if( nJetRad.eq.1 ) then
            call TopDecay(ExtParticle(1),DKJ_LO_T,MomExt(1:4,6:9),GluonHel=GluHel)
          else
            call TopDecay(ExtParticle(1),DKJ_LO_Q,MomExt(1:4,6:9),GluonHel=GluHel)
          endif
          call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,10:12))

          do iPrimAmp=1,NumPrimAmps
              call EvalTree(BornAmps(iPrimAmp))
          enddo
          LO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,NumBornAmps
          do iPrimAmp=1,NumBornAmps
              LO_Res_Pol = LO_Res_Pol + ColLO_ttbggg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result * dconjg(BornAmps(jPrimAmp)%Result)
          enddo
          enddo
          LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
      enddo!helicity loop
      enddo!helicity loop
      EvalCS_DKJ_Real_ttbggg = LO_Res_UnPol * PreFac * (alpha_s4Pi*RunFactor)**3 * ISFac; dummy(1)=dummy(1)+EvalCS_DKJ_Real_ttbggg

      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),EvalCS_DKJ_Real_ttbggg)
      enddo
      EvalCounter = EvalCounter + 1
      EvalCS_DKJ_Real_ttbggg_1 = EvalCS_DKJ_Real_ttbggg_1 + EvalCS_DKJ_Real_ttbggg


50 continue!!   Dipoles for gluon  emission off anti-top

      call EvalDipoles_DKJ_GGTTBG((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:18),PreFac,ATop_,nJetRad,DipoleResult)
      EvalCS_DKJ_Real_ttbggg_1 = EvalCS_DKJ_Real_ttbggg_1 + DipoleResult;  dummy(2)=dummy(2)+DipoleResult

   if( nJetRad.eq.1 ) then
        call WTransform(MomExt(1:4,6:9),MomExtTd(1:4,6:8),pbDpg,ptDpg,ptDpb)
        omz=ptDpg/(ptDpb+ptDpg-pbDpg)  !  for some reason this is not (1-z) as defined in the paper...
        rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
        z=1d0-omz
        y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2

        TheDipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
        TheDipole = TheDipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )

        call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4),MomExt(1:4,5)/),(/MomExtTd(1:4,6),MomExtTd(1:4,7),MomExtTd(1:4,8),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
        if( applyPSCut .or. TheDipole.eq.0d0 ) then
            cycle
        endif

        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
            call HelCrossing(Helicities(iHel,1:NumExtParticles))
            call SetPolarizations()
            call TopDecay(ExtParticle(1),DK_LO,MomExtTd(1:4,6:8))
            call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,10:12))

            do iPrimAmp=1,NumPrimAmps
                call EvalTree(BornAmps(iPrimAmp))
            enddo
            LO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
            do iPrimAmp=1,NumBornAmps
                LO_Res_Pol = LO_Res_Pol + ColLO_ttbggg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result * dconjg(BornAmps(jPrimAmp)%Result)
            enddo
            enddo
            LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop
        DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**3 * ISFac

        do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),DipoleResult)
        enddo
        EvalCS_DKJ_Real_ttbggg_1 = EvalCS_DKJ_Real_ttbggg_1 + DipoleResult; !dummy(3)=DipoleResult


   elseif( nJetRad.eq.2 ) then
          do ndip=1,2
              call wdec_trans(ndip,MomExt(1:4,6:9),MomExtTd(1:4,6:8),alpha_DKWff,TheDipole)
              TheDipole = - alpha_s4Pi*RunFactor * CF * TheDipole

              call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4),MomExt(1:4,5)/),(/MomExtTd(1:4,6),MomExtTd(1:4,7),MomExtTd(1:4,8),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
              if( applyPSCut .or. TheDipole.eq.0d0 ) then
                  cycle
              endif

              LO_Res_Unpol = (0d0,0d0)
              do iHel=1,NumHelicities! helicity summation
                  call HelCrossing(Helicities(iHel,1:NumExtParticles))
                  call SetPolarizations()
                  call TopDecay(ExtParticle(1),DK_LO,MomExtTd(1:4,6:8))
                  call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,10:12))

                  do iPrimAmp=1,NumPrimAmps
                      call EvalTree(BornAmps(iPrimAmp))
                  enddo
                  LO_Res_Pol = (0d0,0d0)
                  do jPrimAmp=1,NumBornAmps
                  do iPrimAmp=1,NumBornAmps
                      LO_Res_Pol = LO_Res_Pol + ColLO_ttbggg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result * dconjg(BornAmps(jPrimAmp)%Result)
                  enddo
                  enddo
                  LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
              enddo!helicity loop
              DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**3 * ISFac

              do NHisto=1,NumHistograms
                  call intoHisto(NHisto,NBin(NHisto),DipoleResult)
              enddo
              EvalCS_DKJ_Real_ttbggg_1 = EvalCS_DKJ_Real_ttbggg_1 + DipoleResult; !dummy(3+ndip)=DipoleResult
        enddo   !dipole loop
  endif! nJetRad for dipoles

enddo! nJetRad loop

! print *, ""
! print *, "Real",dummy(1)
! print *, "DipolesPR",dummy(2)
! print *, "DipolesDK1",dummy(3)
! print *, "DipolesDK2",dummy(4)+dummy(5)
! 
! print *, "ratio", dummy(1)/(dummy(2)+dummy(3)+dummy(4)+dummy(5)), dummy(1)/(dummy(2)+dummy(3)+dummy(4)+dummy(5)) +1d0
! print *, " Sing",(MomExt(1:4,1).dot.MomExt(1:4,3))/EHat**2
! print *, " Sing",(MomExt(1:4,5).dot.MomExt(1:4,9))/EHat**2
! print *, " Sing",(MomExt(1:4,8).dot.MomExt(1:4,9))/EHat**2

! if( (MomExt(1:4,1).dot.MomExt(1:4,3))/EHat**2.lt.1d-5  .and. abs(dummy(1)+dummy(2)).gt.1d0 ) then
!     print *, ""
!     print *, yrnd(1:18)
!     print *, "Real",dummy(1)
!     print *, "DipolesPR",dummy(2)
!     print *, "ratio", dummy(1)/(dummy(2)+dummy(3)+dummy(4)+dummy(5)), dummy(1)/(dummy(2)+dummy(3)+dummy(4)+dummy(5)) +1d0
!     print *, " Sing",(MomExt(1:4,1).dot.MomExt(1:4,3))/EHat**2
!     pause
! endif
! 
! 
! if( (MomExt(1:4,2).dot.MomExt(1:4,3))/EHat**2.lt.1d-5  .and. abs(dummy(1)+dummy(2)).gt.1d0 ) then
!     print *, ""
!     print *, yrnd(1:18)
!     print *, "Real",dummy(1)
!     print *, "DipolesPR",dummy(2)
!     print *, "ratio", dummy(1)/(dummy(2)+dummy(3)+dummy(4)+dummy(5)), dummy(1)/(dummy(2)+dummy(3)+dummy(4)+dummy(5)) +1d0
!     print *, " Sing",(MomExt(1:4,1).dot.MomExt(1:4,3))/EHat**2
!     pause
! endif
! 
! 
! if( (MomExt(1,3))/EHat.lt.1d-4 .and. abs(dummy(1)+dummy(2)).gt.1d0 ) then
!     print *, ""
!     print *, yrnd(1:18)
!     print *, "Real",dummy(1)
!     print *, "DipolesPR",dummy(2)
!     print *, "ratio", dummy(1)/(dummy(2)+dummy(3)+dummy(4)+dummy(5)), dummy(1)/(dummy(2)+dummy(3)+dummy(4)+dummy(5)) +1d0
!     print *, " Sing",(MomExt(1:4,1).dot.MomExt(1:4,3))/EHat**2
!     pause
! endif







!----------------------------------
!    gluon emission off top       |
!----------------------------------
do nJetRad=nJetRad3,nJetRad4!   nJetRad=1: gluon radiation off top line,nJetRad=2: gluon radiation off W line
! dummy(:) = 0d0
   EvalCS_DKJ_Real_ttbggg = 0d0
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomExt(1:4,6:8),PSWgt2)
   if( nJetRad.eq.1 ) then
      call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:18),.true.,MomExt(1:4,9:12),PSWgt3)
   else
      call EvalPhasespace_TopDecay2(MomExt(1:4,5),yRnd(12:18),.true.,MomExt(1:4,9:12),PSWgt3)
   endif
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt*PSWgt2*PSWgt3 * VgsWgt * PDFFac
   if( PreFac.eq.0d0 ) cycle
   call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,12),MomExt(1:4,4),MomExt(1:4,5)/),(/MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,10),MomExt(1:4,11)/),applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_DKJ_Real_ttbggg = 0d0
      PSCutCounter = PSCutCounter + 1
      goto 51
   endif

        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
        do GluHel=1,-1,-2 ! loop over additional gluon polarization
          call HelCrossing(Helicities(iHel,1:NumExtParticles))
          call SetPolarizations()

          call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
          if( nJetRad.eq.1 ) then
            call TopDecay(ExtParticle(2),DKJ_LO_T,MomExt(1:4,9:12),GluonHel=GluHel)
          else
            call TopDecay(ExtParticle(2),DKJ_LO_Q,MomExt(1:4,9:12),GluonHel=GluHel)
          endif

          do iPrimAmp=1,NumPrimAmps
              call EvalTree(BornAmps(iPrimAmp))
          enddo
          LO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,NumBornAmps
          do iPrimAmp=1,NumBornAmps
              LO_Res_Pol = LO_Res_Pol + ColLO_ttbggg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result * dconjg(BornAmps(jPrimAmp)%Result)
          enddo
          enddo
          LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
      enddo!helicity loop
      enddo!helicity loop
      EvalCS_DKJ_Real_ttbggg = LO_Res_UnPol * PreFac * (alpha_s4Pi*RunFactor)**3 * ISFac;    !dummy(1)=EvalCS_DKJ_Real_ttbggg !+ dummy(1)

      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),EvalCS_DKJ_Real_ttbggg)
      enddo
      EvalCounter = EvalCounter + 1
      EvalCS_DKJ_Real_ttbggg_2 = EvalCS_DKJ_Real_ttbggg_2 + EvalCS_DKJ_Real_ttbggg

51 continue

    call EvalDipoles_DKJ_GGTTBG((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:18),PreFac,Top_,nJetRad,DipoleResult)
    EvalCS_DKJ_Real_ttbggg_2 = EvalCS_DKJ_Real_ttbggg_2 + DipoleResult;    !dummy(2)=DipoleResult !+ dummy(2)

   if( nJetRad.eq.1 ) then
        call WTransform(MomExt(1:4,9:12),MomExtTd(1:4,9:11),pbDpg,ptDpg,ptDpb)
        omz=ptDpg/(ptDpb+ptDpg-pbDpg)  !  for some reason this is not (1-z) as defined in the paper...
        rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
        z=1d0-omz
        y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2

        TheDipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
        TheDipole = TheDipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )
        call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4),MomExt(1:4,5)/),(/MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExtTd(1:4,9),MomExtTd(1:4,10),MomExtTd(1:4,11)/),applyPSCut,NBin)
        if( applyPSCut .or. TheDipole.eq.0d0 ) then
            cycle
        endif
        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
            call HelCrossing(Helicities(iHel,1:NumExtParticles))
            call SetPolarizations()
            call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
            call TopDecay(ExtParticle(2),DK_LO,MomExtTd(1:4,9:11))

            do iPrimAmp=1,NumPrimAmps
                call EvalTree(BornAmps(iPrimAmp))
            enddo
            LO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
            do iPrimAmp=1,NumBornAmps
                LO_Res_Pol = LO_Res_Pol + ColLO_ttbggg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result * dconjg(BornAmps(jPrimAmp)%Result)
            enddo
            enddo
            LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop
        DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**3 * ISFac
        do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),DipoleResult)
        enddo
        EvalCS_DKJ_Real_ttbggg_2 = EvalCS_DKJ_Real_ttbggg_2 + DipoleResult;    !dummy(3)=DipoleResult

   elseif( nJetRad.eq.2 ) then

          do ndip=1,2
              call wdec_trans(ndip,MomExt(1:4,9:12),MomExtTd(1:4,9:11),alpha_DKWff,TheDipole)
              TheDipole = - alpha_s4Pi*RunFactor * CF * TheDipole

        call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4),MomExt(1:4,5)/),(/MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExtTd(1:4,9),MomExtTd(1:4,10),MomExtTd(1:4,11)/),applyPSCut,NBin)
              if( applyPSCut .or. TheDipole.eq.0d0 ) then
                  cycle
              endif

              LO_Res_Unpol = (0d0,0d0)
              do iHel=1,NumHelicities! helicity summation
                  call HelCrossing(Helicities(iHel,1:NumExtParticles))
                  call SetPolarizations()
                  call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
                  call TopDecay(ExtParticle(2),DK_LO,MomExtTd(1:4,9:11))

                  do iPrimAmp=1,NumPrimAmps
                      call EvalTree(BornAmps(iPrimAmp))
                  enddo
                  LO_Res_Pol = (0d0,0d0)
                  do jPrimAmp=1,NumBornAmps
                  do iPrimAmp=1,NumBornAmps
                      LO_Res_Pol = LO_Res_Pol + ColLO_ttbggg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result * dconjg(BornAmps(jPrimAmp)%Result)
                  enddo
                  enddo
                  LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
              enddo!helicity loop
              DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**3 * ISFac

              do NHisto=1,NumHistograms
                  call intoHisto(NHisto,NBin(NHisto),DipoleResult)
              enddo
              EvalCS_DKJ_Real_ttbggg_2 = EvalCS_DKJ_Real_ttbggg_2 + DipoleResult;   !dummy(3+ndip)=DipoleResult
        enddo   !dipole loop
  endif! nJetRad for dipoles

! print *, "nJetRad",nJetRad
! print *, "Real",dummy(1)
! print *, "DipolesPR",dummy(2)
! print *, "DipolesDK1",dummy(3)
! print *, "DipolesDK2",dummy(4)+dummy(5)
!
! print *, "ratio", dummy(1)/(dummy(2)+dummy(3)+dummy(4)+dummy(5)), dummy(1)/(dummy(2)+dummy(3)+dummy(4)+dummy(5)) +1d0
! print *, " Sing",(MomExt(1:4,12).dot.MomExt(1:4,11))/EHat**2

enddo! nJetRad loop


! print *, "check",EvalCS_DKJ_Real_ttbggg_1+EvalCS_DKJ_Real_ttbggg_2
! print *, "Real",dummy(1)
! print *, "DipolesPR",dummy(2)
! print *, "DipolesDK1",dummy(3)
! print *, "DipolesDK2",dummy(4)+dummy(5)
!
! print *, "ratio", dummy(1)/(dummy(2)+dummy(3)+dummy(4)+dummy(5)), dummy(1)/(dummy(2)+dummy(3)+dummy(4)+dummy(5)) +1d0
! print *, " Sing",(MomExt(1:4,2).dot.MomExt(1:4,3))/EHat**2
! print *, " Sing",(MomExt(1:4,9).dot.MomExt(1:4,12))/EHat**2
! print *, " Sing",(MomExt(1:4,11).dot.MomExt(1:4,12))/EHat**2
! if ((MomExt(1:4,12).dot.MomExt(1:4,11))/EHat**2.lt.1d-5) pause



EvalCS_DKJ_Real_ttbggg = (EvalCS_DKJ_Real_ttbggg_1+EvalCS_DKJ_Real_ttbggg_2)/VgsWgt
if( isNaN(EvalCS_DKJ_Real_ttbggg) ) then
    print *, "NaN event"
    print *,  EvalCS_DKJ_Real_ttbggg
    print *, EvalCS_DKJ_Real_ttbggg_1,EvalCS_DKJ_Real_ttbggg_2,VgsWgt
    call printYrnd(yrnd(:))
    print *, ""
    EvalCS_DKJ_Real_ttbggg = 0d0
endif


RETURN
END FUNCTION








FUNCTION EvalCS_DKJ_Real_ttbqqbg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModMisc
use ModDipoles_DKJ_QQBTTBG
use ModDipoles_DKJ_QGTTBQ
use ModAmplitudes
use ModParameters
use ModHadrWDecay
implicit none
real(8) ::  EvalCS_DKJ_Real_ttbqqbg,EvalCS_DKJ_Real_ttbqqbg_1,EvalCS_DKJ_Real_ttbqqbg_2
real(8) ::  yRnd(1:VegasMxDim),VgsWgt,DipoleResult
complex(8) :: LO_Res_Pol,LO_Res_Unpol
integer :: iHel,iPrimAmp,jPrimAmp,nJetRad,GluHel,ndip
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3
real(8) :: PDFFac_a,PDFFac_b,PDFFac,pdf(-6:6,1:2),tau,eta1,eta2,FluxFac,sHatJacobi,ISFac,PreFac
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,6:11),dummy(1:5)
real(8) :: pbDpg,ptDpg,ptDpb,omz,rsq,z,y,TheDipole
integer :: NBin(1:NumMaxHisto),NHisto,npdf
logical :: applyPSCut,applySingCut
real(8) :: CF=4d0/3d0
integer,parameter :: nJetRad1=1,nJetRad3=1
integer :: nJetRad2,nJetRad4
include "vegas_common.f"

  EvalCS_DKJ_Real_ttbqqbg = 0d0
  EvalCS_DKJ_Real_ttbqqbg_1 = 0d0
  EvalCS_DKJ_Real_ttbqqbg_2 = 0d0

  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top ) then
      EvalCS_DKJ_Real_ttbqqbg = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to3(EHat,yRnd(3:7),MomExt(1:4,1:5),PSWgt)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))

   call CheckSing(MomExt,applySingCut)
   if( applySingCut ) then
      EvalCS_DKJ_Real_ttbqqbg = 0d0
      SkipCounter = SkipCounter + 1
      return
   endif


   RunFactor = RunAlphaS(2,MuRen)
   call setPDFs(eta1,eta2,MuFac,pdf)
IF( PROCESS.EQ.38 ) THEN
   PDFFac_a = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
            + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
            + pdf(Bot_,1)*pdf(ABot_,2)
   PDFFac_b = pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
            + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
            + pdf(Bot_,2)*pdf(ABot_,1)
ELSEIF( PROCESS.EQ.35 ) THEN
   PDFFac_a = pdf(Up_,1) *pdf(0,2) + pdf(Dn_,1) *pdf(0,2)  &
            + pdf(Chm_,1)*pdf(0,2) + pdf(Str_,1)*pdf(0,2)  &
            + pdf(Bot_,1)*pdf(0,2)
   PDFFac_b = pdf(Up_,2) *pdf(0,1) + pdf(Dn_,2) *pdf(0,1)  &
            + pdf(Chm_,2)*pdf(0,1) + pdf(Str_,2)*pdf(0,1)  &
            + pdf(Bot_,2)*pdf(0,1)
ELSEIF( PROCESS.EQ.36 ) THEN
   PDFFac_a = pdf(AUp_,1) *pdf(0,2) + pdf(ADn_,1) *pdf(0,2)  &
            + pdf(AChm_,1)*pdf(0,2) + pdf(AStr_,1)*pdf(0,2)  &
            + pdf(ABot_,1)*pdf(0,2)
   PDFFac_b = pdf(AUp_,2) *pdf(0,1) + pdf(ADn_,2) *pdf(0,1)  &
            + pdf(AChm_,2)*pdf(0,1) + pdf(AStr_,2)*pdf(0,1)  &
            + pdf(ABot_,2)*pdf(0,1)
ENDIF



   if( TOPDECAYS.eq.1 ) then
      nJetRad2=1;  nJetRad4=1
   elseif( TOPDECAYS.eq.2 ) then
      nJetRad2=2;  nJetRad4=2
   elseif( TOPDECAYS.eq.3 ) then
      nJetRad2=1;  nJetRad4=2
   elseif( TOPDECAYS.eq.4 ) then
      nJetRad2=2;  nJetRad4=1
   endif
!----------------------------------
! photon emission off anti-top    |
!----------------------------------
EvalCS_DKJ_Real_ttbqqbg_1 = 0d0
do nJetRad=nJetRad1,nJetRad2!   nJetRad=1: photon radiation off top/bot/W, nJetRad=2: photon radiation off W/lep
! dummy(:) = 0d0
   if( nJetRad.eq.1 ) then
      call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:14),.true.,MomExt(1:4,6:9),PSWgt2)
   else
      call EvalPhasespace_TopDecay2(MomExt(1:4,4),yRnd(8:14),.true.,MomExt(1:4,6:9),PSWgt2)
   endif
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(15:18),.false.,MomExt(1:4,10:12),PSWgt3)

   do npdf=1,2
        if(npdf.eq.1) then
            PDFFac = PDFFac_a
        elseif(npdf.eq.2) then
            PDFFac = PDFFac_b
            call swapMom(MomExt(1:4,1),MomExt(1:4,2))
        endif
        ISFac = MomCrossing(MomExt)
        PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt*PSWgt2*PSWgt3 * VgsWgt * PDFFac
        if( PreFac.eq.0d0 ) cycle

        call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,9),MomExt(1:4,4),MomExt(1:4,5)/),(/MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
        if( applyPSCut ) then
            EvalCS_DKJ_Real_ttbqqbg = 0d0
              PSCutCounter = PSCutCounter + 1
            goto 60
        endif

        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities
        do GluHel=1,-1,-2 ! loop over additional gluon polarization
          call HelCrossing(Helicities(iHel,1:NumExtParticles))
          call SetPolarizations()

          if( nJetRad.eq.1 ) then
            call TopDecay(ExtParticle(1),DKJ_LO_T,MomExt(1:4,6:9),GluonHel=GluHel)
          else
            call TopDecay(ExtParticle(1),DKJ_LO_Q,MomExt(1:4,6:9),GluonHel=GluHel)
          endif
          call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,10:12))

          do iPrimAmp=1,NumPrimAmps
              call EvalTree(BornAmps(iPrimAmp))
          enddo
          LO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,NumBornAmps
          do iPrimAmp=1,NumBornAmps
              LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqbg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result * dconjg(BornAmps(jPrimAmp)%Result)
          enddo
          enddo
          LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
      enddo!helicity loop
      enddo!helicity loop
      EvalCS_DKJ_Real_ttbqqbg = LO_Res_UnPol * PreFac * (alpha_s4Pi*RunFactor)**3 * ISFac

      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),EvalCS_DKJ_Real_ttbqqbg)
      enddo
      EvalCounter = EvalCounter + 1
      EvalCS_DKJ_Real_ttbqqbg_1 = EvalCS_DKJ_Real_ttbqqbg_1 + EvalCS_DKJ_Real_ttbqqbg
! dummy(1) = EvalCS_DKJ_Real_ttbqqbg + dummy(1)

60 continue

      IF( PROCESS.EQ.38 ) THEN
          call EvalDipoles_DKJ_QQBTTBG((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:18),PreFac,ATop_,nJetRad,DipoleResult)
      ELSEIF( PROCESS.EQ.35 ) THEN
          call EvalDipoles_DKJ_QGTTBQ((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:18),PreFac,ATop_,nJetRad,DipoleResult)
      ELSEIF( PROCESS.EQ.36 ) THEN
          call EvalDipoles_DKJ_QGTTBQ((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:18),PreFac,ATop_,nJetRad,DipoleResult)
      ENDIF

      EvalCS_DKJ_Real_ttbqqbg_1 = EvalCS_DKJ_Real_ttbqqbg_1 + DipoleResult
! dummy(2) = DipoleResult + dummy(2)






   if( nJetRad.eq.1 ) then

        call WTransform(MomExt(1:4,6:9),MomExtTd(1:4,6:8),pbDpg,ptDpg,ptDpb)
        omz=ptDpg/(ptDpb+ptDpg-pbDpg)  !  for some reason this is not (1-z) as defined in the paper...
        rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
        z=1d0-omz
        y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2

        TheDipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
        TheDipole = TheDipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )

        call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4),MomExt(1:4,5)/),(/MomExtTd(1:4,6),MomExtTd(1:4,7),MomExtTd(1:4,8),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
        if( applyPSCut .or. TheDipole.eq.0d0 ) then
            cycle
        endif

        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
            call HelCrossing(Helicities(iHel,1:NumExtParticles))
            call SetPolarizations()
            call TopDecay(ExtParticle(1),DK_LO,MomExtTd(1:4,6:8))
            call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,10:12))

            do iPrimAmp=1,NumPrimAmps
                call EvalTree(BornAmps(iPrimAmp))
            enddo
            LO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
            do iPrimAmp=1,NumBornAmps
                LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqbg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result * dconjg(BornAmps(jPrimAmp)%Result)
            enddo
            enddo
            LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop
        DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**3 * ISFac

        do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),DipoleResult)
        enddo
        EvalCS_DKJ_Real_ttbqqbg_1 = EvalCS_DKJ_Real_ttbqqbg_1 + DipoleResult
! dummy(3) = DipoleResult + dummy(3)

   elseif( nJetRad.eq.2 ) then

          do ndip=1,2
              call wdec_trans(ndip,MomExt(1:4,6:9),MomExtTd(1:4,6:8),alpha_DKWff,TheDipole)
              TheDipole = - alpha_s4Pi*RunFactor * CF * TheDipole

              call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4),MomExt(1:4,5)/),(/MomExtTd(1:4,6),MomExtTd(1:4,7),MomExtTd(1:4,8),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
              if( applyPSCut .or. TheDipole.eq.0d0 ) then
                  cycle
              endif

              LO_Res_Unpol = (0d0,0d0)
              do iHel=1,NumHelicities! helicity summation
                  call HelCrossing(Helicities(iHel,1:NumExtParticles))
                  call SetPolarizations()
                  call TopDecay(ExtParticle(1),DK_LO,MomExtTd(1:4,6:8))
                  call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,10:12))

                  do iPrimAmp=1,NumPrimAmps
                      call EvalTree(BornAmps(iPrimAmp))
                  enddo
                  LO_Res_Pol = (0d0,0d0)
                  do jPrimAmp=1,NumBornAmps
                  do iPrimAmp=1,NumBornAmps
                      LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqbg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result * dconjg(BornAmps(jPrimAmp)%Result)
                  enddo
                  enddo
                  LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
              enddo!helicity loop
              DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**3 * ISFac

              do NHisto=1,NumHistograms
                  call intoHisto(NHisto,NBin(NHisto),DipoleResult)
              enddo
              EvalCS_DKJ_Real_ttbqqbg_1 = EvalCS_DKJ_Real_ttbqqbg_1 + DipoleResult
! dummy(4) = DipoleResult + dummy(4)
          enddo!dipole loop
  endif! nJetRad for dipoles

enddo! nPdf loop
call swapMom(MomExt(1:4,1),MomExt(1:4,2))! swap back

! print *, "nJetRad",nJetRad
! print *, "Real",dummy(1)
! print *, "DipolesPR",dummy(2)
! print *, "DipolesDK1",dummy(3)
! print *, "DipolesDK2",dummy(4)+dummy(5)
!
! print *, "ratio", dummy(1)/(dummy(2)+dummy(3)+dummy(4)+dummy(5)), dummy(1)/(dummy(2)+dummy(3)+dummy(4)+dummy(5)) +1d0
! print *, " Sing",(MomExt(1:4,8).dot.MomExt(1:4,9))/EHat**2
! pause

enddo! nJetRad loop









!----------------------------------
! photon emission off top         |
!----------------------------------
  EvalCS_DKJ_Real_ttbqqbg_2 = 0d0
do nJetRad=nJetRad3,nJetRad4!   nJetRad=1: photon radiation off top/bot/W, nJetRad=2: photon radiation off W/lep
! dummy(:) = 0d0
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomExt(1:4,6:8),PSWgt2)
   if( nJetRad.eq.1 ) then
      call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:18),.true.,MomExt(1:4,9:12),PSWgt3)
   else
      call EvalPhasespace_TopDecay2(MomExt(1:4,5),yRnd(12:18),.true.,MomExt(1:4,9:12),PSWgt3)
   endif
   do npdf=1,2
        if(npdf.eq.1) then
            PDFFac = PDFFac_a
        elseif(npdf.eq.2) then
            PDFFac = PDFFac_b
            call swapMom(MomExt(1:4,1),MomExt(1:4,2))
        endif
        ISFac = MomCrossing(MomExt)
        PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt*PSWgt2*PSWgt3 * VgsWgt * PDFFac
        if( PreFac.eq.0d0 ) then
            SkipCounter = SkipCounter + 1
            cycle
        endif

        call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,12),MomExt(1:4,4),MomExt(1:4,5)/),(/MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,10),MomExt(1:4,11)/),applyPSCut,NBin)
        if( applyPSCut ) then
            EvalCS_DKJ_Real_ttbqqbg = 0d0
            PSCutCounter = PSCutCounter + 1
            goto 61
        endif

        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities
        do GluHel=1,-1,-2 ! loop over additional gluon polarization
          call HelCrossing(Helicities(iHel,1:NumExtParticles))
          call SetPolarizations()
          call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
          if( nJetRad.eq.1 ) then
            call TopDecay(ExtParticle(2),DKJ_LO_T,MomExt(1:4,9:12),GluonHel=GluHel)
          else
            call TopDecay(ExtParticle(2),DKJ_LO_Q,MomExt(1:4,9:12),GluonHel=GluHel)
          endif

          do iPrimAmp=1,NumPrimAmps
              call EvalTree(BornAmps(iPrimAmp))
          enddo
          LO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,NumBornAmps
          do iPrimAmp=1,NumBornAmps
              LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqbg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result * dconjg(BornAmps(jPrimAmp)%Result)
          enddo
          enddo
          LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
      enddo!helicity loop
      enddo!helicity loop
      EvalCS_DKJ_Real_ttbqqbg = LO_Res_UnPol * PreFac * (alpha_s4Pi*RunFactor)**3 * ISFac

      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),EvalCS_DKJ_Real_ttbqqbg)
      enddo
      EvalCS_DKJ_Real_ttbqqbg_2 = EvalCS_DKJ_Real_ttbqqbg_2 + EvalCS_DKJ_Real_ttbqqbg
      EvalCounter = EvalCounter + 1
! dummy(1) = EvalCS_DKJ_Real_ttbqqbg + dummy(1)

61 continue

      IF( PROCESS.EQ.38 ) THEN
          call EvalDipoles_DKJ_QQBTTBG((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:18),PreFac,Top_,nJetRad,DipoleResult)
      ELSEIF( PROCESS.EQ.35 ) THEN
          call EvalDipoles_DKJ_QGTTBQ((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:18),PreFac,Top_,nJetRad,DipoleResult)
      ELSEIF( PROCESS.EQ.36 ) THEN
          call EvalDipoles_DKJ_QGTTBQ((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:18),PreFac,Top_,nJetRad,DipoleResult)
      ENDIF
      EvalCS_DKJ_Real_ttbqqbg_2 = EvalCS_DKJ_Real_ttbqqbg_2 + DipoleResult
! dummy(2) = DipoleResult + dummy(2)



   if( nJetRad.eq.1 ) then

        call WTransform(MomExt(1:4,9:12),MomExtTd(1:4,9:11),pbDpg,ptDpg,ptDpb)
        omz=ptDpg/(ptDpb+ptDpg-pbDpg)  !  for some reason this is not (1-z) as defined in the paper...
        rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
        z=1d0-omz
        y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2

        TheDipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
        TheDipole = TheDipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )

        call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4),MomExt(1:4,5)/),(/MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExtTd(1:4,9),MomExtTd(1:4,10),MomExtTd(1:4,11)/),applyPSCut,NBin)
        if( applyPSCut .or. TheDipole.eq.0d0 ) then
            cycle
        endif

        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
            call HelCrossing(Helicities(iHel,1:NumExtParticles))
            call SetPolarizations()
            call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
            call TopDecay(ExtParticle(2),DK_LO,MomExtTd(1:4,9:11))

            do iPrimAmp=1,NumPrimAmps
                call EvalTree(BornAmps(iPrimAmp))
            enddo
            LO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
            do iPrimAmp=1,NumBornAmps
                LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqbg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result * dconjg(BornAmps(jPrimAmp)%Result)
            enddo
            enddo
            LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop
        DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**3 * ISFac

        do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),DipoleResult)
        enddo
        EvalCS_DKJ_Real_ttbqqbg_2 = EvalCS_DKJ_Real_ttbqqbg_2 + DipoleResult
! dummy(3) = DipoleResult + dummy(3)

   elseif( nJetRad.eq.2 ) then

          do ndip=1,2
              call wdec_trans(ndip,MomExt(1:4,9:12),MomExtTd(1:4,9:11),alpha_DKWff,TheDipole)
              TheDipole = - alpha_s4Pi*RunFactor * CF * TheDipole

              call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4),MomExt(1:4,5)/),(/MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExtTd(1:4,9),MomExtTd(1:4,10),MomExtTd(1:4,11)/),applyPSCut,NBin)
              if( applyPSCut .or. TheDipole.eq.0d0 ) then
                  cycle
              endif

              LO_Res_Unpol = (0d0,0d0)
              do iHel=1,NumHelicities! helicity summation
                  call HelCrossing(Helicities(iHel,1:NumExtParticles))
                  call SetPolarizations()
                  call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
                  call TopDecay(ExtParticle(2),DK_LO,MomExtTd(1:4,9:11))

                  do iPrimAmp=1,NumPrimAmps
                      call EvalTree(BornAmps(iPrimAmp))
                  enddo
                  LO_Res_Pol = (0d0,0d0)
                  do jPrimAmp=1,NumBornAmps
                  do iPrimAmp=1,NumBornAmps
                      LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqbg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result * dconjg(BornAmps(jPrimAmp)%Result)
                  enddo
                  enddo
                  LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
              enddo!helicity loop
              DipoleResult = LO_Res_UnPol * PreFac * TheDipole * (alpha_s4Pi*RunFactor)**3 * ISFac

              do NHisto=1,NumHistograms
                  call intoHisto(NHisto,NBin(NHisto),DipoleResult)
              enddo
              EvalCS_DKJ_Real_ttbqqbg_2 = EvalCS_DKJ_Real_ttbqqbg_2 + DipoleResult
! dummy(4) = DipoleResult + dummy(4)
        enddo   !dipole loop
  endif! nJetRad for dipoles

enddo! nPdf loop
call swapMom(MomExt(1:4,1),MomExt(1:4,2))! swap back

! print *, "nJetRad",nJetRad
! print *, "Real",dummy(1)
! print *, "DipolesPR",dummy(2)
! print *, "DipolesDK1",dummy(3)
! print *, "DipolesDK2",dummy(4)+dummy(5)
!
! print *, "ratio", dummy(1)/(dummy(2)+dummy(3)+dummy(4)+dummy(5)), dummy(1)/(dummy(2)+dummy(3)+dummy(4)+dummy(5)) +1d0
! print *, " Sing",(MomExt(1:4,11).dot.MomExt(1:4,12))/EHat**2
! if((MomExt(1:4,11).dot.MomExt(1:4,12))/EHat**2.lt.1d-6) pause

enddo! nJetRad loop






EvalCS_DKJ_Real_ttbqqbg = (EvalCS_DKJ_Real_ttbqqbg_1+EvalCS_DKJ_Real_ttbqqbg_2)/VgsWgt

if( isNaN(EvalCS_DKJ_Real_ttbqqbg) ) then
    print *, "NaN event"
    print *, EvalCS_DKJ_Real_ttbqqbg
    print *, EvalCS_DKJ_Real_ttbqqbg_1,EvalCS_DKJ_Real_ttbqqbg_2,VgsWgt
    call printYrnd(yrnd(:))
    print *, ""
    EvalCS_DKJ_Real_ttbqqbg = 0d0
endif
RETURN
END FUNCTION








FUNCTION EvalCS_1LDKJ_ttb(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModDipolesDKJTTB
implicit none
real(8) ::  EvalCS_1LDKJ_ttb,yRnd(1:VegasMxDim),VgsWgt
complex(8) :: LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol,NLO_Res_UnPol
complex(8) :: TreeResult(1:NumBornAmps), LoopResult(1:NumBornAmps)
integer :: mu,iHel,GluHel,iPrimAmp,jPrimAmp, Contr
real(8) :: EHat,PSWgt1,PSWgt2,PSWgt3,ISFac
real(8) :: MomExt(1:4,1:12)
logical :: applyPSCut,applySingCut
real(8) :: tau,eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,RunFactor
real(8) :: pdf(-6:6,1:2),ColCorrLO(1:NumBornAmps,1:NumBornAmps)
integer :: NBin(1:NumMaxHisto),NHisto,nPhoRad
type(Particle) :: ATop(0:1,0:1)
type(Particle) ::  Top(0:1,0:1)
integer :: ATopCorr,TopCorr
integer, parameter :: tree_ = 0
integer, parameter :: virt_ = 1
integer, parameter :: plus_ = 0
integer, parameter :: minus_= 1
integer, parameter :: NumContr=16
include "vegas_common.f"



EvalCS_1LDKJ_ttb = 0d0


  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top) then
      EvalCS_1LDKJ_ttb = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to2(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt1)
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))

   ISFac = MomCrossing(MomExt)
   RunFactor = RunAlphaS(NLOParam,MuRen)
   call setPDFs(eta1,eta2,MuFac,pdf)
IF( PROCESS.EQ.33 ) THEN
   ColCorrLO(1:NumBornAmps,1:NumBornAmps) = ColLO_ttbgg(1:NumBornAmps,1:NumBornAmps)
   PDFFac = pdf(0,1) * pdf(0,2)
ELSEIF( PROCESS.EQ.34 ) THEN
   ColCorrLO(1:NumBornAmps,1:NumBornAmps) = ColLO_ttbqqb(1:NumBornAmps,1:NumBornAmps)
   PDFFac =       ( pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
                  + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
                  + pdf(Bot_,1)*pdf(ABot_,2) ) +                         &
                  ( pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
                  + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
                  + pdf(Bot_,2)*pdf(ABot_,1) )
ENDIF


    ATop(0:1,0:1)%PartType = ATop_
    ATop(0:1,0:1)%Mass = M_Top
    ATop(0:1,0:1)%Mass2 = M_Top**2
    Top(0:1,0:1)%PartType = Top_
    Top(0:1,0:1)%Mass = M_Top
    Top(0:1,0:1)%Mass2 = M_Top**2
    do mu=1,4
        ATop(0:1,0:1)%Mom(mu) = dcmplx(MomExt(mu,3))
        Top(0:1,0:1)%Mom(mu)  = dcmplx(MomExt(mu,4))
    enddo



 Contr=0
DO WHILE ( Contr.LT.NumContr )
    Contr = Contr + 1

                          IF( CONTR.ne.nGluRadContr ) CYCLE; !print *, "nGluRadContr=",nGluRadContr

! note: it turns out that contr. 6=10 and 2=14 because the virt.corr. on W is prop. to tree level, i.e. a form factor
!------------------------------------------------------------------------
!     virtual correction on anti-top:  jet radiation from anti-top      |
!------------------------------------------------------------------------
IF( Contr.eq.1 ) THEN
      call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,3),yRnd(5:11),MomExt(1:4,5:8),PSWgt2)
      call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,4),yRnd(12:15),MomExt(1:4,9:11),PSWgt3)
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,9),MomExt(1:4,10),MomExt(1:4,11)/),applyPSCut,NBin)
      if( applyPSCut ) then
          cycle
      endif
      call TopDecay(ATop(tree_,plus_), DKJ_LO_T,MomExt(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(tree_,minus_),DKJ_LO_T,MomExt(1:4,5:8),GluonHel=int((-1)**minus_))
      call TopDecay(ATop(virt_,plus_), DKJ_1L_T,MomExt(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(virt_,minus_),DKJ_1L_T,MomExt(1:4,5:8),GluonHel=int((-1)**minus_))

      call TopDecay(Top(tree_,plus_),DK_LO,MomExt(1:4,9:11));  Top(tree_,minus_)%Pol(1:4) = Top(tree_,plus_)%Pol(1:4)
      ATopCorr = virt_
      TopCorr  = tree_

!------------------------------------------------------------------
!     virtual correction on W-:  jet radiation from anti-top      |
!------------------------------------------------------------------
ELSEIF( Contr.eq.2 ) THEN
      if( TopDecays.ne.2 .and. TopDecays.ne.4 ) cycle
            call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,3),yRnd(5:11),MomExt(1:4,5:8),PSWgt2)
            call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,4),yRnd(12:15),MomExt(1:4,9:11),PSWgt3)
            call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,9),MomExt(1:4,10),MomExt(1:4,11)/),applyPSCut,NBin)
            if( applyPSCut ) then
      !           Contr=2
                cycle
            endif

      call TopDecay(ATop(tree_,plus_), DKJ_LO_T,MomExt(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(tree_,minus_),DKJ_LO_T,MomExt(1:4,5:8),GluonHel=int((-1)**minus_))
      call TopDecay(ATop(virt_,plus_), DKJ_1LQ_T,MomExt(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(virt_,minus_),DKJ_1LQ_T,MomExt(1:4,5:8),GluonHel=int((-1)**minus_))

      call TopDecay(Top(tree_,plus_),DK_LO,MomExt(1:4,9:11));  Top(tree_,minus_)%Pol(1:4) = Top(tree_,plus_)%Pol(1:4)
      ATopCorr = virt_
      TopCorr  = tree_



!------------------------------------------------------------
!     virtual correction on W-:  jet radiation from W-      |
!------------------------------------------------------------
ELSEIF( Contr.eq.3 ) THEN
      if( TopDecays.ne.2 .and. TopDecays.ne.4 ) cycle
      call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,3),yRnd(5:11),MomExt(1:4,5:8),PSWgt2)
      call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,4),yRnd(12:15),MomExt(1:4,9:11),PSWgt3)
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,9),MomExt(1:4,10),MomExt(1:4,11)/),applyPSCut,NBin)
      if( applyPSCut )  then
!           Contr=4
          cycle
      endif
      call TopDecay(ATop(tree_,plus_), DKJ_LO_Q,MomExt(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(tree_,minus_),DKJ_LO_Q,MomExt(1:4,5:8),GluonHel=int((-1)**minus_))
      call TopDecay(ATop(virt_,plus_), DKJ_1LQ_Q,MomExt(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(virt_,minus_),DKJ_1LQ_Q,MomExt(1:4,5:8),GluonHel=int((-1)**minus_))

      call TopDecay(Top(tree_,plus_),DK_LO,MomExt(1:4,9:11));  Top(tree_,minus_)%Pol(1:4) = Top(tree_,plus_)%Pol(1:4)
      ATopCorr = virt_
      TopCorr  = tree_



!------------------------------------------------------------------
!     virtual correction on anti-top:  jet radiation from W-      |
!------------------------------------------------------------------
ELSEIF( Contr.eq.4 ) THEN
      if( TopDecays.ne.2 .and. TopDecays.ne.4 ) cycle
            call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,3),yRnd(5:11),MomExt(1:4,5:8),PSWgt2)
            call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,4),yRnd(12:15),MomExt(1:4,9:11),PSWgt3)
            call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,9),MomExt(1:4,10),MomExt(1:4,11)/),applyPSCut,NBin)
            if( applyPSCut )  then
      !           Contr=4
                cycle
            endif

      call TopDecay(ATop(tree_,plus_), DKJ_LO_Q,MomExt(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(tree_,minus_),DKJ_LO_Q,MomExt(1:4,5:8),GluonHel=int((-1)**minus_))
      call TopDecay(ATop(virt_,plus_), DKJ_1L_Q,MomExt(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(virt_,minus_),DKJ_1L_Q,MomExt(1:4,5:8),GluonHel=int((-1)**minus_))

      call TopDecay(Top(tree_,plus_),DK_LO,MomExt(1:4,9:11));  Top(tree_,minus_)%Pol(1:4) = Top(tree_,plus_)%Pol(1:4)
      ATopCorr = virt_
      TopCorr  = tree_



!-------------------------------------------------------------------
!     virtual correction on anti-top:  jet radiation from top      |
!-------------------------------------------------------------------
ELSEIF( Contr.eq.5 ) THEN
      call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,3),yRnd(5:8),MomExt(1:4,5:7),PSWgt2)
      call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,4),yRnd(9:15),MomExt(1:4,8:11),PSWgt3)
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,10)/),applyPSCut,NBin)
      if( applyPSCut )  then
!           Contr=6
          cycle
      endif
      call TopDecay(ATop(tree_,plus_),DK_LO,MomExt(1:4,5:7));   ATop(tree_,minus_)%Pol(1:4) = ATop(tree_,plus_)%Pol(1:4)
      call TopDecay(ATop(virt_,plus_),DK_1L_T,MomExt(1:4,5:7)); ATop(virt_,minus_)%Pol(1:4) = ATop(virt_,plus_)%Pol(1:4)

!DEC$ IF(_FactCheck .EQ.1)
      if( dabs(MomExt(1:4,8).dot.MomExt(1:4,11))/m_Top**2 .lt. 1d-3 ) then
          cycle
      endif
!DEC$ ENDIF

      call TopDecay(Top(tree_,plus_), DKJ_LO_T,MomExt(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(tree_,minus_),DKJ_LO_T,MomExt(1:4,8:11),GluonHel=int((-1)**minus_))
      ATopCorr = virt_
      TopCorr  = tree_



!-------------------------------------------------------------
!     virtual correction on W-:  jet radiation from top      |
!-------------------------------------------------------------
ELSEIF( Contr.eq.6 ) THEN
      if( TopDecays.ne.2 .and. TopDecays.ne.4 ) cycle
            call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,3),yRnd(5:8),MomExt(1:4,5:7),PSWgt2)
            call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,4),yRnd(9:15),MomExt(1:4,8:11),PSWgt3)
            call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,10)/),applyPSCut,NBin)
            if( applyPSCut )  then
      !           Contr=6
                cycle
            endif
      call TopDecay(ATop(tree_,plus_),DK_LO,MomExt(1:4,5:7));   ATop(tree_,minus_)%Pol(1:4) = ATop(tree_,plus_)%Pol(1:4)
      call TopDecay(ATop(virt_,plus_),DK_1L_Q,MomExt(1:4,5:7)); ATop(virt_,minus_)%Pol(1:4) = ATop(virt_,plus_)%Pol(1:4)

      call TopDecay(Top(tree_,plus_), DKJ_LO_T,MomExt(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(tree_,minus_),DKJ_LO_T,MomExt(1:4,8:11),GluonHel=int((-1)**minus_))
      ATopCorr = virt_
      TopCorr  = tree_



!------------------------------------------------------------------
!     virtual correction on anti-top:  jet radiation from W+      |
!------------------------------------------------------------------
ELSEIF( Contr.eq.7 ) THEN
      if( TopDecays.ne.2 .and. TopDecays.ne.3 ) cycle
      call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,3),yRnd(5:8),MomExt(1:4,5:7),PSWgt2)
      call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,4),yRnd(9:15),MomExt(1:4,8:11),PSWgt3)
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,10)/),applyPSCut,NBin)
      if( applyPSCut )  then
!           Contr=8
          cycle
      endif
      call TopDecay(ATop(tree_,plus_),DK_LO,MomExt(1:4,5:7));   ATop(tree_,minus_)%Pol(1:4) = ATop(tree_,plus_)%Pol(1:4)
      call TopDecay(ATop(virt_,plus_),DK_1L_T,MomExt(1:4,5:7)); ATop(virt_,minus_)%Pol(1:4) = ATop(virt_,plus_)%Pol(1:4)

      call TopDecay(Top(tree_,plus_), DKJ_LO_Q,MomExt(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(tree_,minus_),DKJ_LO_Q,MomExt(1:4,8:11),GluonHel=int((-1)**minus_))
      ATopCorr = virt_
      TopCorr  = tree_



!------------------------------------------------------------
!     virtual correction on W-:  jet radiation from W+      |
!------------------------------------------------------------
ELSEIF( Contr.eq.8 ) THEN
      if( TopDecays.ne.2 ) cycle
              call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,3),yRnd(5:8),MomExt(1:4,5:7),PSWgt2)
              call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,4),yRnd(9:15),MomExt(1:4,8:11),PSWgt3)
              call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,10)/),applyPSCut,NBin)
              if( applyPSCut )  then
!                   Contr=8
                  cycle
              endif

      call TopDecay(ATop(tree_,plus_),DK_LO,MomExt(1:4,5:7));   ATop(tree_,minus_)%Pol(1:4) = ATop(tree_,plus_)%Pol(1:4)
      call TopDecay(ATop(virt_,plus_),DK_1L_Q,MomExt(1:4,5:7)); ATop(virt_,minus_)%Pol(1:4) = ATop(virt_,plus_)%Pol(1:4)
      call TopDecay(Top(tree_,plus_), DKJ_LO_Q,MomExt(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(tree_,minus_),DKJ_LO_Q,MomExt(1:4,8:11),GluonHel=int((-1)**minus_))
      ATopCorr = virt_
      TopCorr  = tree_



!--------------------------------------------------------------
!     virtual correction on top:  jet radiation from top      |
!--------------------------------------------------------------
ELSEIF( Contr.eq.9 ) THEN
      call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,3),yRnd(5:8),MomExt(1:4,5:7),PSWgt2)
      call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,4),yRnd(9:15),MomExt(1:4,8:11),PSWgt3)
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,10)/),applyPSCut,NBin)
      if( applyPSCut )  then
!           Contr=10
          cycle
      endif

      call TopDecay(ATop(tree_,plus_),DK_LO,MomExt(1:4,5:7)); ATop(tree_,minus_)%Pol(1:4) = ATop(tree_,plus_)%Pol(1:4)
      call TopDecay(Top(tree_,plus_), DKJ_LO_T,MomExt(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(tree_,minus_),DKJ_LO_T,MomExt(1:4,8:11),GluonHel=int((-1)**minus_))
      call TopDecay(Top(virt_,plus_), DKJ_1L_T,MomExt(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(virt_,minus_),DKJ_1L_T,MomExt(1:4,8:11),GluonHel=int((-1)**minus_))
      ATopCorr = tree_
      TopCorr  = virt_

!-------------------------------------------------------------
!     virtual correction on W+:  jet radiation from top      |
!-------------------------------------------------------------
ELSEIF( Contr.eq.10 ) THEN
      if( TopDecays.ne.2 .and. TopDecays.ne.3 ) cycle
            call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,3),yRnd(5:8),MomExt(1:4,5:7),PSWgt2)
            call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,4),yRnd(9:15),MomExt(1:4,8:11),PSWgt3)
            call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,10)/),applyPSCut,NBin)
            if( applyPSCut )  then
!                 Contr=10
                cycle
            endif

      call TopDecay(ATop(tree_,plus_),DK_LO,MomExt(1:4,5:7)); ATop(tree_,minus_)%Pol(1:4) = ATop(tree_,plus_)%Pol(1:4)

      call TopDecay(Top(tree_,plus_), DKJ_LO_T,MomExt(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(tree_,minus_),DKJ_LO_T,MomExt(1:4,8:11),GluonHel=int((-1)**minus_))
      call TopDecay(Top(virt_,plus_), DKJ_1LQ_T,MomExt(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(virt_,minus_),DKJ_1LQ_T,MomExt(1:4,8:11),GluonHel=int((-1)**minus_))
      ATopCorr = tree_
      TopCorr  = virt_



!------------------------------------------------------------
!     virtual correction on W+:  jet radiation from W+      |
!------------------------------------------------------------
ELSEIF( Contr.eq.11 ) THEN
      if( TopDecays.ne.2 .and. TopDecays.ne.3 ) cycle
      call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,3),yRnd(5:8),MomExt(1:4,5:7),PSWgt2)
      call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,4),yRnd(9:15),MomExt(1:4,8:11),PSWgt3)
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,10)/),applyPSCut,NBin)
      if( applyPSCut )  then
!           Contr=12
          cycle
      endif

      call TopDecay(ATop(tree_,plus_),DK_LO,MomExt(1:4,5:7)); ATop(tree_,minus_)%Pol(1:4) = ATop(tree_,plus_)%Pol(1:4)

      call TopDecay(Top(tree_,plus_), DKJ_LO_Q,MomExt(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(tree_,minus_),DKJ_LO_Q,MomExt(1:4,8:11),GluonHel=int((-1)**minus_))
      call TopDecay(Top(virt_,plus_), DKJ_1LQ_Q,MomExt(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(virt_,minus_),DKJ_1LQ_Q,MomExt(1:4,8:11),GluonHel=int((-1)**minus_))
      ATopCorr = tree_
      TopCorr  = virt_



!-------------------------------------------------------------
!     virtual correction on top:  jet radiation from W+      |
!-------------------------------------------------------------
ELSEIF( Contr.eq.12 ) THEN
      if( TopDecays.ne.2 .and. TopDecays.ne.3 ) cycle
            call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,3),yRnd(5:8),MomExt(1:4,5:7),PSWgt2)
            call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,4),yRnd(9:15),MomExt(1:4,8:11),PSWgt3)
            call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,11),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,10)/),applyPSCut,NBin)
            if( applyPSCut )  then
!                 Contr=12
                cycle
            endif

      call TopDecay(ATop(tree_,plus_),DK_LO,MomExt(1:4,5:7)); ATop(tree_,minus_)%Pol(1:4) = ATop(tree_,plus_)%Pol(1:4)

      call TopDecay(Top(tree_,plus_), DKJ_LO_Q,MomExt(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(tree_,minus_),DKJ_LO_Q,MomExt(1:4,8:11),GluonHel=int((-1)**minus_))
      call TopDecay(Top(virt_,plus_), DKJ_1L_Q,MomExt(1:4,8:11),GluonHel=int((-1)**plus_))
      call TopDecay(Top(virt_,minus_),DKJ_1L_Q,MomExt(1:4,8:11),GluonHel=int((-1)**minus_))
      ATopCorr = tree_
      TopCorr  = virt_



!-------------------------------------------------------------------
!     virtual correction on top:  jet radiation from anti-top      |
!-------------------------------------------------------------------
ELSEIF( Contr.eq.13 ) THEN
      call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,3),yRnd(5:11),MomExt(1:4,5:8),PSWgt2)
      call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,4),yRnd(12:15),MomExt(1:4,9:11),PSWgt3)
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,9),MomExt(1:4,10),MomExt(1:4,11)/),applyPSCut,NBin)
      if( applyPSCut )  then
!           Contr=14
          cycle
      endif

      call TopDecay(ATop(tree_,plus_), DKJ_LO_T,MomExt(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(tree_,minus_),DKJ_LO_T,MomExt(1:4,5:8),GluonHel=int((-1)**minus_))

!DEC$ IF(_FactCheck .EQ.1)
      if( dabs(MomExt(1:4,5).dot.MomExt(1:4,8))/m_Top**2 .lt. 1d-3 ) then
          cycle
      endif
!DEC$ ENDIF

      call TopDecay(Top(tree_,plus_),DK_LO,MomExt(1:4,9:11));   Top(tree_,minus_)%Pol(1:4) = Top(tree_,plus_)%Pol(1:4)
      call TopDecay(Top(virt_,plus_),DK_1L_T,MomExt(1:4,9:11)); Top(virt_,minus_)%Pol(1:4) = Top(virt_,plus_)%Pol(1:4)
      ATopCorr = tree_
      TopCorr  = virt_



!------------------------------------------------------------------
!     virtual correction on W+:  jet radiation from anti-top      |
!------------------------------------------------------------------
ELSEIF( Contr.eq.14 ) THEN
      if( TopDecays.ne.2 .and. TopDecays.ne.3 ) cycle
            call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,3),yRnd(5:11),MomExt(1:4,5:8),PSWgt2)
            call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,4),yRnd(12:15),MomExt(1:4,9:11),PSWgt3)
            call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,9),MomExt(1:4,10),MomExt(1:4,11)/),applyPSCut,NBin)
            if( applyPSCut )  then
!                 Contr=14
                cycle
            endif

      call TopDecay(ATop(tree_,plus_), DKJ_LO_T,MomExt(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(tree_,minus_),DKJ_LO_T,MomExt(1:4,5:8),GluonHel=int((-1)**minus_))

      call TopDecay(Top(tree_,plus_),DK_LO,MomExt(1:4,9:11));   Top(tree_,minus_)%Pol(1:4) = Top(tree_,plus_)%Pol(1:4)
      call TopDecay(Top(virt_,plus_),DK_1L_Q,MomExt(1:4,9:11)); Top(virt_,minus_)%Pol(1:4) = Top(virt_,plus_)%Pol(1:4)
      ATopCorr = tree_
      TopCorr  = virt_



!-------------------------------------------------------------
!     virtual correction on top:  jet radiation from W-      |
!-------------------------------------------------------------
ELSEIF( Contr.eq.15 ) THEN
      if( TopDecays.ne.2 .and. TopDecays.ne.4 ) cycle
      call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,3),yRnd(5:11),MomExt(1:4,5:8),PSWgt2)
      call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,4),yRnd(12:15),MomExt(1:4,9:11),PSWgt3)
      call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,9),MomExt(1:4,10),MomExt(1:4,11)/),applyPSCut,NBin)
      if( applyPSCut )  then
!           Contr=16
          cycle
      endif

      call TopDecay(ATop(tree_,plus_), DKJ_LO_Q,MomExt(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(tree_,minus_),DKJ_LO_Q,MomExt(1:4,5:8),GluonHel=int((-1)**minus_))

      call TopDecay(Top(tree_,plus_),DK_LO,MomExt(1:4,9:11));    Top(tree_,minus_)%Pol(1:4) = Top(tree_,plus_)%Pol(1:4)
      call TopDecay(Top(virt_,plus_),DK_1L_T,MomExt(1:4,9:11)); Top(virt_,minus_)%Pol(1:4) = Top(virt_,plus_)%Pol(1:4)
      ATopCorr = tree_
      TopCorr  = virt_


!------------------------------------------------------------
!     virtual correction on W+:  jet radiation from W-      |
!------------------------------------------------------------
ELSEIF( Contr.eq.16 ) THEN
      if( TopDecays.ne.2 ) cycle
            call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,3),yRnd(5:11),MomExt(1:4,5:8),PSWgt2)
            call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,4),yRnd(12:15),MomExt(1:4,9:11),PSWgt3)
            call Kinematics_TTBARJET(0,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,9),MomExt(1:4,10),MomExt(1:4,11)/),applyPSCut,NBin)
            if( applyPSCut )  then
!                 Contr=16
                cycle
            endif

      call TopDecay(ATop(tree_,plus_), DKJ_LO_Q,MomExt(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(tree_,minus_),DKJ_LO_Q,MomExt(1:4,5:8),GluonHel=int((-1)**minus_))

      call TopDecay(Top(tree_,plus_),DK_LO,MomExt(1:4,9:11));   Top(tree_,minus_)%Pol(1:4) = Top(tree_,plus_)%Pol(1:4)
      call TopDecay(Top(virt_,plus_),DK_1L_Q,MomExt(1:4,9:11)); Top(virt_,minus_)%Pol(1:4) = Top(virt_,plus_)%Pol(1:4)
      ATopCorr = tree_
      TopCorr  = virt_


ENDIF! Contr





      PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt * PDFFac
      NLO_Res_Unpol = (0d0,0d0)
      do iHel=1,NumHelicities ! loop over polarizations
      do GluHel=0,1           ! loop over additional gluon polarization
            call HelCrossing(Helicities(iHel,1:NumExtParticles))
            call SetPolarizations()

            ExtParticle(1)%Pol(1:4) = ATop(tree_,GluHel)%Pol(1:4)
            ExtParticle(2)%Pol(1:4) =  Top(tree_,GluHel)%Pol(1:4)
            do iPrimAmp=1,NumBornAmps
               call EvalTree(BornAmps(iPrimAmp))
               TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
            enddo

            ExtParticle(1)%Pol(1:4) = ATop(ATopCorr,GluHel)%Pol(1:4)
            ExtParticle(2)%Pol(1:4) =  Top( TopCorr,GluHel)%Pol(1:4)
            do iPrimAmp=1,NumBornAmps
               call EvalTree(BornAmps(iPrimAmp))
               LoopResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
            enddo

            NLO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
               do iPrimAmp=1,NumBornAmps
                  NLO_Res_Pol = NLO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * dreal(TreeResult(jPrimAmp)*dconjg(LoopResult(iPrimAmp)))
               enddo
            enddo
            NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
      enddo!helicity loop
      enddo!helicity loop

      NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
      EvalCS_1LDKJ_ttb = EvalCS_1LDKJ_ttb + dble(NLO_Res_Unpol)

! print *, "Virt",Contr,NLO_Res_Unpol
! pause
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),dreal(NLO_Res_UnPol))
      enddo

ENDDO! Contr



   EvalCS_1LDKJ_ttb = EvalCS_1LDKJ_ttb/VgsWgt
RETURN
END FUNCTION









FUNCTION EvalCS_REDKJ_ttb(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModDipolesDKJTTB
implicit none
real(8) ::  EvalCS_REDKJ_ttb,yRnd(1:VegasMxDim),VgsWgt
complex(8) :: LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol,NLO_Res_UnPol
complex(8) :: TreeResult(1:NumBornAmps,0:1)! the last index enumerates either the DKPrimAmp (in real emission) OR helicity correlation (in dioples)
integer ::  mu,iHel,GluHel1,GluHel2,iPrimAmp,jPrimAmp,iPrimAmpDK,jPrimAmpDK,Contr
real(8) :: EHat,PSWgt,PSWgt1,PSWgt2,PSWgt3,ISFac
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,5:11)
complex(8) :: Dipole(0:1,0:1)
logical :: applyPSCut,applySingCut
real(8) :: tau,eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,RunFactor
real(8) :: pdf(-6:6,1:2),ColCorrLO(1:NumBornAmps,1:NumBornAmps)
integer :: NBin(1:NumMaxHisto),NHisto,nDipole
real(8) :: dummy(1:5)=0d0
integer, parameter :: NumContr = 20
integer, parameter :: plus_ = 0, minus_= 1
real(8), parameter :: SymmFac = 0.5d0
type(Particle) :: ATop(plus_:minus_,plus_:minus_,1:2)! indices: Contribution, GluHel1, GluHel2, ColFacDK
type(Particle) ::  Top(plus_:minus_,plus_:minus_,1:2)
real(8) :: ColCorrDK(1:2,1:2),HelCorr(1:2,1:2)
integer :: TopCorr!                                                  |
integer, parameter :: NumDKPrimAmp(1:NumContr)=(/2,2, 1,1,1,1,1,2,2, 1,1,2,2,2,1,2,2,1,2,2/)
integer, parameter :: NumDipoles(1:NumContr)  =(/7,10,3,2,3,3,4,7,10,3,1,2,2,2,2,4,4,2,4,4/)
#DEFINE TheDipolesREDKJ 3789
include "vegas_common.f"



! print *,"fixed yRnd"
! yrnd(1)=   5.158616905837767d-003
! yrnd(2)=   0.226384259038744d0
! yrnd(3)=   0.639774218319495d0
! yrnd(4)=   0.363389155775710d0
      ! yrnd(5)=   0.899851546503038d0
      ! yrnd(6)=   0.116913734572266d0
      ! yrnd(7)=   0.284107823003173d0
      ! yrnd(8)=   0.748669132523022d0
      ! yrnd(9)=   0.494439351378572d0
      ! yrnd(10)=   0.182568801269892d0
      ! yrnd(11)=   3.504963000886823d-002
      ! yrnd(12)=   0.819452873633415d0
      ! yrnd(13)=   0.898317497701006d0
      ! yrnd(14)=   6.124850891079013d-002
! yrnd(15)=   0.898218286552350d0
! yrnd(16)=   6.977205951980933d-002
! yrnd(17)=   0.553658300509745d0
! yrnd(18)=   0.567500606778878d0






  EvalCS_REDKJ_ttb = 0d0

  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top ) then
      EvalCS_REDKJ_ttb = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to2(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt1)
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))

   ISFac = MomCrossing(MomExt)
   RunFactor = RunAlphaS(NLOParam,MuRen)
   call setPDFs(eta1,eta2,MuFac,pdf)
IF( PROCESS.EQ.33 ) THEN
   ColCorrLO(1:NumBornAmps,1:NumBornAmps) = ColLO_ttbgg(1:NumBornAmps,1:NumBornAmps)
   PDFFac = pdf(0,1) * pdf(0,2)
ELSEIF( PROCESS.EQ.34 ) THEN
   ColCorrLO(1:NumBornAmps,1:NumBornAmps) = ColLO_ttbqqb(1:NumBornAmps,1:NumBornAmps)
   PDFFac =       ( pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
                  + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
                  + pdf(Bot_,1)*pdf(ABot_,2) ) +                         &
                  ( pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
                  + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
                  + pdf(Bot_,2)*pdf(ABot_,1) )
ENDIF



    ATop(0:1,0:1,1:2)%PartType = ATop_
    ATop(0:1,0:1,1:2)%Mass = M_Top
    ATop(0:1,0:1,1:2)%Mass2 = M_Top**2
    Top(0:1,0:1,1:2)%PartType = Top_
    Top(0:1,0:1,1:2)%Mass = M_Top
    Top(0:1,0:1,1:2)%Mass2 = M_Top**2
    do mu=1,4
        ATop(0:1,0:1,1:2)%Mom(mu) = dcmplx(MomExt(mu,3))
        Top(0:1,0:1,1:2)%Mom(mu)  = dcmplx(MomExt(mu,4))
    enddo

dummy(:) = 0d0

!------------------------------------------------------------
!                   real emission spinors                   |
!------------------------------------------------------------
 Contr=0
DO WHILE ( Contr.LT.NumContr )
 Contr = Contr + 1

                          IF( CONTR.ne.nGluRadContr ) CYCLE;        !PRINT *, "CONTRIBUTION ",CONTR
IF( Contr.eq.1 ) THEN
      TopCorr = ATop_
      call EvalPhasespace_TopDK(T_BGG_W,MomExt(1:4,3),yRnd(5:14),MomExt(1:4,5:9),PSWgt2); PSWgt2=PSWgt2*SymmFac
      call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,4),yRnd(15:18),MomExt(1:4,10:12),PSWgt3)
      PSWgt = PSWgt1*PSWgt2*PSWgt3
      if( PSWgt.eq.0d0 ) cycle
      call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)

      if( applyPSCut ) goto TheDipolesREDKJ
      do iPrimAmpDK=1,NumDKPrimAmp(Contr)
            call TopDecay(ATop(plus_,plus_,iPrimAmpDK),  DKJ_RE_TT,MomExt(1:4,5:9),PartAmp=iPrimAmpDK,GluonHel=int((-1)**plus_),Gluon2Hel=int((-1)**plus_))
            call TopDecay(ATop(plus_,minus_,iPrimAmpDK), DKJ_RE_TT,MomExt(1:4,5:9),PartAmp=iPrimAmpDK,GluonHel=int((-1)**plus_),Gluon2Hel=int((-1)**minus_))
            call TopDecay(ATop(minus_,plus_,iPrimAmpDK), DKJ_RE_TT,MomExt(1:4,5:9),PartAmp=iPrimAmpDK,GluonHel=int((-1)**minus_),Gluon2Hel=int((-1)**plus_))
            call TopDecay(ATop(minus_,minus_,iPrimAmpDK),DKJ_RE_TT,MomExt(1:4,5:9),PartAmp=iPrimAmpDK,GluonHel=int((-1)**minus_),Gluon2Hel=int((-1)**minus_))
      enddo
      call TopDecay(Top(0,0,1),DK_LO,MomExt(1:4,10:12))
      ColCorrDK(1:2,1:2) = ColLO_ttbgg(1:2,1:2) * 1d0/3d0


ELSEIF( Contr.eq.2 ) THEN
      if( TopDecays.ne.2 .and. TopDecays.ne.4 ) cycle
      TopCorr = ATop_
      call EvalPhasespace_TopDK(T_B_WGG,MomExt(1:4,3),yRnd(5:14),MomExt(1:4,5:9),PSWgt2); !PSWgt2=PSWgt2*SymmFac removed because Andreas has it inside
      call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,4),yRnd(15:18),MomExt(1:4,10:12),PSWgt3)
      PSWgt = PSWgt1*PSWgt2*PSWgt3
      if( PSWgt.eq.0d0 ) cycle
      call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
      if( applyPSCut ) goto TheDipolesREDKJ
      do iPrimAmpDK=1,NumDKPrimAmp(Contr)
            call TopDecay(ATop(plus_,plus_,iPrimAmpDK),  DKJ_RE_QQ,MomExt(1:4,5:9),PartAmp=iPrimAmpDK,GluonHel=int((-1)**plus_),Gluon2Hel=int((-1)**plus_))
            call TopDecay(ATop(plus_,minus_,iPrimAmpDK), DKJ_RE_QQ,MomExt(1:4,5:9),PartAmp=iPrimAmpDK,GluonHel=int((-1)**plus_),Gluon2Hel=int((-1)**minus_))
            call TopDecay(ATop(minus_,plus_,iPrimAmpDK), DKJ_RE_QQ,MomExt(1:4,5:9),PartAmp=iPrimAmpDK,GluonHel=int((-1)**minus_),Gluon2Hel=int((-1)**plus_))
            call TopDecay(ATop(minus_,minus_,iPrimAmpDK),DKJ_RE_QQ,MomExt(1:4,5:9),PartAmp=iPrimAmpDK,GluonHel=int((-1)**minus_),Gluon2Hel=int((-1)**minus_))
      enddo
      call TopDecay(Top(0,0,1),DK_LO,MomExt(1:4,10:12))
      ColCorrDK(1,1) = 16d0/3d0;    ColCorrDK(2,2) = 16d0/3d0;
      ColCorrDK(1,2) =- 2d0/3d0;    ColCorrDK(2,1) =- 2d0/3d0;

!        check factor 2 with andreas!! symm. fact??


ELSEIF( Contr.eq.3 ) THEN!    note: it seems that alpha_DKWff should not be bigger than 0.1, otherwise it's instable
      if( TopDecays.ne.2 .and. TopDecays.ne.4 ) cycle
      TopCorr = ATop_
      call EvalPhasespace_TopDK(T_BG_WG,MomExt(1:4,3),yRnd(5:14),MomExt(1:4,5:9),PSWgt2)
      call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,4),yRnd(15:18),MomExt(1:4,10:12),PSWgt3)
      PSWgt = PSWgt1*PSWgt2*PSWgt3
      if( PSWgt.eq.0d0 ) cycle
      call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
!DEC$ IF(_FactCheck .EQ.1)
  if( (MomExt(1:4,5).dot.MomExt(1:4,9))/m_Top**2.lt.1d-3 .and. ( (MomExt(1:4,6).dot.MomExt(1:4,8))/m_W**2.lt.1d-3 .or. (MomExt(1:4,7).dot.MomExt(1:4,8))/m_W**2.lt.1d-3) ) applyPSCut=.true.
!DEC$ ENDIF
      if( applyPSCut ) goto TheDipolesREDKJ
      call TopDecay(ATop(plus_,plus_,1),  DKJ_RE_TQ,MomExt(1:4,5:9),GluonHel=int((-1)**plus_),Gluon2Hel=int((-1)**plus_))
      call TopDecay(ATop(plus_,minus_,1), DKJ_RE_TQ,MomExt(1:4,5:9),GluonHel=int((-1)**plus_),Gluon2Hel=int((-1)**minus_))
      call TopDecay(ATop(minus_,plus_,1), DKJ_RE_TQ,MomExt(1:4,5:9),GluonHel=int((-1)**minus_),Gluon2Hel=int((-1)**plus_))
      call TopDecay(ATop(minus_,minus_,1),DKJ_RE_TQ,MomExt(1:4,5:9),GluonHel=int((-1)**minus_),Gluon2Hel=int((-1)**minus_))
      call TopDecay(Top(0,0,1),DK_LO,MomExt(1:4,10:12))
      ColCorrDK(1:2,1:2) = 1d0


ELSEIF( Contr.eq.4 ) THEN
      TopCorr = ATop_+Top_
      call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,3),yRnd(5:11),MomExt(1:4,5:8),PSWgt2)
      call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,4),yRnd(12:18),MomExt(1:4,9:12),PSWgt3)
      PSWgt = PSWgt1*PSWgt2*PSWgt3
      if( PSWgt.eq.0d0 ) cycle
      call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,12),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,9),MomExt(1:4,10),MomExt(1:4,11)/),applyPSCut,NBin)
!DEC$ IF(_FactCheck .EQ.1)
  if( (MomExt(1:4,5).dot.MomExt(1:4,8))/m_Top**2.lt.1d-3 .and. (MomExt(1:4,9).dot.MomExt(1:4,12))/m_Top**2.lt.1d-3 ) applyPSCut=.true.
!DEC$ ENDIF
      if( applyPSCut ) goto TheDipolesREDKJ
      call TopDecay(ATop(plus_,0,1), DKJ_LO_T,MomExt(1:4,5:8), GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,0,1),DKJ_LO_T,MomExt(1:4,5:8), GluonHel=int((-1)**minus_))
      call TopDecay(Top(plus_,0,1),  DKJ_LO_T,MomExt(1:4,9:12),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,0,1), DKJ_LO_T,MomExt(1:4,9:12),GluonHel=int((-1)**minus_))
      ColCorrDK(1:2,1:2) = 1d0



ELSEIF( Contr.eq.5 ) THEN
      if( TopDecays.ne.2 .and. TopDecays.ne.3 ) cycle
      TopCorr = ATop_+Top_
      call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,3),yRnd(5:11),MomExt(1:4,5:8),PSWgt2)
      call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,4),yRnd(12:18),MomExt(1:4,9:12),PSWgt3)
      PSWgt = PSWgt1*PSWgt2*PSWgt3
      if( PSWgt.eq.0d0 ) cycle
      call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,12),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,9),MomExt(1:4,10),MomExt(1:4,11)/),applyPSCut,NBin)
      if( applyPSCut ) goto TheDipolesREDKJ
      call TopDecay(ATop(plus_,0,1), DKJ_LO_T,MomExt(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,0,1),DKJ_LO_T,MomExt(1:4,5:8),GluonHel=int((-1)**minus_))
      call TopDecay(Top(plus_,0,1),  DKJ_LO_Q,MomExt(1:4,9:12),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,0,1), DKJ_LO_Q,MomExt(1:4,9:12),GluonHel=int((-1)**minus_))
      ColCorrDK(1:2,1:2) = 1d0



ELSEIF( Contr.eq.6 ) THEN
      if( TopDecays.ne.2 .and. TopDecays.ne.4 ) cycle
      TopCorr = ATop_+Top_
      call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,3),yRnd(5:11),MomExt(1:4,5:8),PSWgt2)
      call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,4),yRnd(12:18),MomExt(1:4,9:12),PSWgt3)
      PSWgt = PSWgt1*PSWgt2*PSWgt3
      if( PSWgt.eq.0d0 ) cycle
      call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,12),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,9),MomExt(1:4,10),MomExt(1:4,11)/),applyPSCut,NBin)
!DEC$ IF(_FactCheck .EQ.1)
  if( (MomExt(1:4,9).dot.MomExt(1:4,12))/m_Top**2.lt.1d-3 .and. ( (MomExt(1:4,6).dot.MomExt(1:4,8))/m_W**2.lt.1d-3 .or. (MomExt(1:4,7).dot.MomExt(1:4,8))/m_W**2.lt.1d-3) ) applyPSCut=.true.
!DEC$ ENDIF
      if( applyPSCut ) goto TheDipolesREDKJ
      call TopDecay(ATop(plus_,0,1), DKJ_LO_Q,MomExt(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,0,1),DKJ_LO_Q,MomExt(1:4,5:8),GluonHel=int((-1)**minus_))
      call TopDecay(Top(plus_,0,1),  DKJ_LO_T,MomExt(1:4,9:12),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,0,1), DKJ_LO_T,MomExt(1:4,9:12),GluonHel=int((-1)**minus_))
      ColCorrDK(1:2,1:2) = 1d0



ELSEIF( Contr.eq.7 ) THEN
      if( TopDecays.ne.2 ) cycle
      TopCorr = ATop_+Top_
      call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,3),yRnd(5:11),MomExt(1:4,5:8),PSWgt2)
      call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,4),yRnd(12:18),MomExt(1:4,9:12),PSWgt3)
      PSWgt = PSWgt1*PSWgt2*PSWgt3
      if( PSWgt.eq.0d0 ) cycle
      call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,12),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,9),MomExt(1:4,10),MomExt(1:4,11)/),applyPSCut,NBin)
      if( applyPSCut ) goto TheDipolesREDKJ
      call TopDecay(ATop(plus_,0,1),  DKJ_LO_Q,MomExt(1:4,5:8),GluonHel=int((-1)**plus_))
      call TopDecay(ATop(minus_,0,1), DKJ_LO_Q,MomExt(1:4,5:8),GluonHel=int((-1)**minus_))
      call TopDecay(Top(plus_,0,1),  DKJ_LO_Q,MomExt(1:4,9:12),GluonHel=int((-1)**plus_))
      call TopDecay(Top(minus_,0,1), DKJ_LO_Q,MomExt(1:4,9:12),GluonHel=int((-1)**minus_))
      ColCorrDK(1:2,1:2) = 1d0



ELSEIF( Contr.eq.8 ) THEN
      TopCorr = Top_
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,3),yRnd(5:8),MomExt(1:4,5:7),PSWgt2)
      call EvalPhasespace_TopDK(T_BGG_W, MomExt(1:4,4),yRnd(9:18),MomExt(1:4,8:12),PSWgt3); PSWgt3=PSWgt3*SymmFac
      PSWgt = PSWgt1*PSWgt2*PSWgt3
      if( PSWgt.eq.0d0 ) cycle
      call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,11),MomExt(1:4,12),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,10)/),applyPSCut,NBin)
      if( applyPSCut ) goto TheDipolesREDKJ

      call TopDecay(ATop(0,0,1),DK_LO,MomExt(1:4,5:7))
      do iPrimAmpDK=1,NumDKPrimAmp(Contr)
            call TopDecay(Top(plus_,plus_,iPrimAmpDK),  DKJ_RE_TT,MomExt(1:4,8:12),PartAmp=iPrimAmpDK,GluonHel=int((-1)**plus_),Gluon2Hel=int((-1)**plus_))
            call TopDecay(Top(plus_,minus_,iPrimAmpDK), DKJ_RE_TT,MomExt(1:4,8:12),PartAmp=iPrimAmpDK,GluonHel=int((-1)**plus_),Gluon2Hel=int((-1)**minus_))
            call TopDecay(Top(minus_,plus_,iPrimAmpDK), DKJ_RE_TT,MomExt(1:4,8:12),PartAmp=iPrimAmpDK,GluonHel=int((-1)**minus_),Gluon2Hel=int((-1)**plus_))
            call TopDecay(Top(minus_,minus_,iPrimAmpDK),DKJ_RE_TT,MomExt(1:4,8:12),PartAmp=iPrimAmpDK,GluonHel=int((-1)**minus_),Gluon2Hel=int((-1)**minus_))
      enddo
      ColCorrDK(1:2,1:2) = ColLO_ttbgg(1:2,1:2) * 1d0/3d0



ELSEIF( Contr.eq.9 ) THEN
      if( TopDecays.ne.2 .and. TopDecays.ne.3 ) cycle
      TopCorr = Top_
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,3),yRnd(5:8),MomExt(1:4,5:7),PSWgt2)
      call EvalPhasespace_TopDK(T_B_WGG, MomExt(1:4,4),yRnd(9:18),MomExt(1:4,8:12),PSWgt3); !PSWgt3=PSWgt3*SymmFac removed because Andreas has it inside
      PSWgt = PSWgt1*PSWgt2*PSWgt3
      if( PSWgt.eq.0d0 ) cycle
      call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,11),MomExt(1:4,12),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,10)/),applyPSCut,NBin)
      if( applyPSCut ) goto TheDipolesREDKJ

      call TopDecay(ATop(0,0,1),DK_LO,MomExt(1:4,5:7))
      do iPrimAmpDK=1,NumDKPrimAmp(Contr)
            call TopDecay(Top(plus_,plus_,iPrimAmpDK),  DKJ_RE_QQ,MomExt(1:4,8:12),PartAmp=iPrimAmpDK,GluonHel=int((-1)**plus_),Gluon2Hel=int((-1)**plus_))
            call TopDecay(Top(plus_,minus_,iPrimAmpDK), DKJ_RE_QQ,MomExt(1:4,8:12),PartAmp=iPrimAmpDK,GluonHel=int((-1)**plus_),Gluon2Hel=int((-1)**minus_))
            call TopDecay(Top(minus_,plus_,iPrimAmpDK), DKJ_RE_QQ,MomExt(1:4,8:12),PartAmp=iPrimAmpDK,GluonHel=int((-1)**minus_),Gluon2Hel=int((-1)**plus_))
            call TopDecay(Top(minus_,minus_,iPrimAmpDK),DKJ_RE_QQ,MomExt(1:4,8:12),PartAmp=iPrimAmpDK,GluonHel=int((-1)**minus_),Gluon2Hel=int((-1)**minus_))
      enddo
      ColCorrDK(1,1) = 16d0/3d0;    ColCorrDK(2,2) = 16d0/3d0;
      ColCorrDK(1,2) =- 2d0/3d0;    ColCorrDK(2,1) =- 2d0/3d0;

!        "check factor 2 with andreas!! symm. fact??


ELSEIF( Contr.eq.10 ) THEN
      if( TopDecays.ne.2 .and. TopDecays.ne.3 ) cycle
      TopCorr = Top_
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,3),yRnd(5:8),MomExt(1:4,5:7),PSWgt2)
      call EvalPhasespace_TopDK(T_BG_WG, MomExt(1:4,4),yRnd(9:18),MomExt(1:4,8:12),PSWgt3)
      PSWgt = PSWgt1*PSWgt2*PSWgt3
      if( PSWgt.eq.0d0 ) cycle
      call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,11),MomExt(1:4,12),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,10)/),applyPSCut,NBin)
      if( applyPSCut ) goto TheDipolesREDKJ

      call TopDecay(ATop(0,0,1),DK_LO,MomExt(1:4,5:7))
      call TopDecay(Top(plus_,plus_,1),  DKJ_RE_TQ,MomExt(1:4,8:12),GluonHel=int((-1)**plus_),Gluon2Hel=int((-1)**plus_))
      call TopDecay(Top(plus_,minus_,1), DKJ_RE_TQ,MomExt(1:4,8:12),GluonHel=int((-1)**plus_),Gluon2Hel=int((-1)**minus_))
      call TopDecay(Top(minus_,plus_,1), DKJ_RE_TQ,MomExt(1:4,8:12),GluonHel=int((-1)**minus_),Gluon2Hel=int((-1)**plus_))
      call TopDecay(Top(minus_,minus_,1),DKJ_RE_TQ,MomExt(1:4,8:12),GluonHel=int((-1)**minus_),Gluon2Hel=int((-1)**minus_))
      ColCorrDK(1:2,1:2) = 1d0



ELSEIF( Contr.eq.11 ) THEN!   qqb-splitting into u,d,c,s quarks
      TopCorr = ATop_
      call EvalPhasespace_TopDK(T_BGG_W,MomExt(1:4,3),yRnd(5:14),MomExt(1:4,5:9),PSWgt2);    PSWgt2 = PSWgt2 * 4d0! this accounts for 4 quark flavors
      call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,4),yRnd(15:18),MomExt(1:4,10:12),PSWgt3)
      PSWgt = PSWgt1*PSWgt2*PSWgt3
      if( PSWgt.eq.0d0 ) cycle
      call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
      if( applyPSCut ) goto TheDipolesREDKJ
      call TopDecay(ATop(plus_,0,1),  DKJ_RE2_TT,MomExt(1:4,5:9),GluonHel=int((-1)**plus_));  ATop(plus_, 1,1)%Pol(1:4) = (0d0,0d0)
      call TopDecay(ATop(minus_,0,1), DKJ_RE2_TT,MomExt(1:4,5:9),GluonHel=int((-1)**minus_)); ATop(minus_,1,1)%Pol(1:4) = (0d0,0d0)
      call TopDecay(Top(0,0,1),DK_LO,MomExt(1:4,10:12))
      ColCorrDK(1:2,1:2) = 0d0
      ColCorrDK(1,1) = ColLO_ttbqqb(1,1) * 1d0/3d0



ELSEIF( Contr.eq.12 ) THEN!   qqb-splitting into b quarks
      TopCorr = ATop_
      call EvalPhasespace_TopDK(T_BGG_W,MomExt(1:4,3),yRnd(5:14),MomExt(1:4,5:9),PSWgt2); PSWgt2=PSWgt2*SymmFac
      call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,4),yRnd(15:18),MomExt(1:4,10:12),PSWgt3)
      PSWgt = PSWgt1*PSWgt2*PSWgt3
      if( PSWgt.eq.0d0 ) cycle
      call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
      if( applyPSCut ) goto TheDipolesREDKJ

      iPrimAmpDK=1
      call TopDecay(ATop(plus_ ,0,iPrimAmpDK),DKJ_RE2_TT,MomExt(1:4,5:9),GluonHel=int((-1)**plus_));  ATop(plus_, 1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)
      call TopDecay(ATop(minus_,0,iPrimAmpDK),DKJ_RE2_TT,MomExt(1:4,5:9),GluonHel=int((-1)**minus_)); ATop(minus_,1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)

      iPrimAmpDK=2
      call swapMom(MomExt(1:4,5),MomExt(1:4,8))
      call TopDecay(ATop(plus_ ,0,iPrimAmpDK),DKJ_RE2_TT,MomExt(1:4,5:9),GluonHel=int((-1)**plus_));  ATop(plus_, 1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)
      call TopDecay(ATop(minus_,0,iPrimAmpDK),DKJ_RE2_TT,MomExt(1:4,5:9),GluonHel=int((-1)**minus_)); ATop(minus_,1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)
      call swapMom(MomExt(1:4,5),MomExt(1:4,8))
      ColCorrDK(1,1) = 8d0;       ColCorrDK(2,2) = 8d0
      ColCorrDK(1,2) = 8d0/3d0;   ColCorrDK(2,1) = 8d0/3d0
      ColCorrDK(1:2,1:2) = ColCorrDK(1:2,1:2)  *1d0/3d0

      call TopDecay(Top(0,0,1),DK_LO,MomExt(1:4,10:12))



ELSEIF( Contr.eq.13 ) THEN!   qqb-splitting into u,d,c,s quarks
      TopCorr = Top_
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,3),yRnd(5:8),MomExt(1:4,5:7),PSWgt2)
      call EvalPhasespace_TopDK(T_BGG_W, MomExt(1:4,4),yRnd(9:18),MomExt(1:4,8:12),PSWgt3);    PSWgt3 = PSWgt3 * 4d0! this accounts for 4 quark flavors
      PSWgt = PSWgt1*PSWgt2*PSWgt3
      if( PSWgt.eq.0d0 ) cycle
      call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,11),MomExt(1:4,12),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,10)/),applyPSCut,NBin)
      if( applyPSCut ) goto TheDipolesREDKJ

      call TopDecay(ATop(0,0,1),DK_LO,MomExt(1:4,5:7))
      call TopDecay(Top(plus_,0,1),  DKJ_RE2_TT,MomExt(1:4,8:12),GluonHel=int((-1)**plus_));  Top(plus_, 1,1)%Pol(1:4) = (0d0,0d0)
      call TopDecay(Top(minus_,0,1), DKJ_RE2_TT,MomExt(1:4,8:12),GluonHel=int((-1)**minus_)); Top(minus_,1,1)%Pol(1:4) = (0d0,0d0)
      ColCorrDK(1:2,1:2) = 0d0
      ColCorrDK(1,1) = ColLO_ttbqqb(1,1) * 1d0/3d0



ELSEIF( Contr.eq.14 ) THEN!   qqb-splitting into b quarks
      TopCorr = Top_
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,3),yRnd(5:8),MomExt(1:4,5:7),PSWgt2); PSWgt2=PSWgt2*SymmFac
      call EvalPhasespace_TopDK(T_BGG_W, MomExt(1:4,4),yRnd(9:18),MomExt(1:4,8:12),PSWgt3)
      PSWgt = PSWgt1*PSWgt2*PSWgt3
      if( PSWgt.eq.0d0 ) cycle
      call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,11),MomExt(1:4,12),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,10)/),applyPSCut,NBin)
      if( applyPSCut ) goto TheDipolesREDKJ

      call TopDecay(ATop(0,0,1),DK_LO,MomExt(1:4,5:7))

      iPrimAmpDK=1
      call TopDecay(Top(plus_,0,iPrimAmpDK), DKJ_RE2_TT,MomExt(1:4,8:12),GluonHel=int((-1)**plus_));  Top(plus_, 1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)
      call TopDecay(Top(minus_,0,iPrimAmpDK),DKJ_RE2_TT,MomExt(1:4,8:12),GluonHel=int((-1)**minus_)); Top(minus_,1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)

      iPrimAmpDK=2
      call swapMom(MomExt(1:4,8),MomExt(1:4,12))
      call TopDecay(Top(plus_,0,iPrimAmpDK), DKJ_RE2_TT,MomExt(1:4,8:12),GluonHel=int((-1)**plus_));  Top(plus_, 1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)
      call TopDecay(Top(minus_,0,iPrimAmpDK),DKJ_RE2_TT,MomExt(1:4,8:12),GluonHel=int((-1)**minus_)); Top(minus_,1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)
      call swapMom(MomExt(1:4,8),MomExt(1:4,12))
      ColCorrDK(1,1) = 8d0;       ColCorrDK(2,2) = 8d0
      ColCorrDK(1,2) = 8d0/3d0;   ColCorrDK(2,1) = 8d0/3d0
      ColCorrDK(1:2,1:2) = ColCorrDK(1:2,1:2)  *1d0/3d0



ELSEIF( Contr.eq.15 ) THEN!     W -> (ud+cc/ss/bb) and (cs+uu/dd/bb)
      if( TopDecays.ne.2 .and. TopDecays.ne.4 ) cycle
      TopCorr = ATop_
      call EvalPhasespace_TopDK(T_B_WGG,MomExt(1:4,3),yRnd(5:14),MomExt(1:4,5:9),PSWgt2);     PSWgt2 = PSWgt2 * 6d0! this accounts for 6 quark flavor combinations
      call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,4),yRnd(15:18),MomExt(1:4,10:12),PSWgt3)
      PSWgt = PSWgt1*PSWgt2*PSWgt3
      if( PSWgt.eq.0d0 ) cycle
      call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
      if( applyPSCut ) goto TheDipolesREDKJ
      call TopDecay(ATop(plus_, 0,1),DKJ_RE2_QQ,MomExt(1:4,5:9),GluonHel=int((-1)**plus_));  ATop(plus_, 1,1)%Pol(1:4) = (0d0,0d0)
      call TopDecay(ATop(minus_,0,1),DKJ_RE2_QQ,MomExt(1:4,5:9),GluonHel=int((-1)**minus_)); ATop(minus_,1,1)%Pol(1:4) = (0d0,0d0)
      ColCorrDK(1,1) = 2d0;
!       ColCorrDK(2,2) = 2d0! ??? why two?
!       ColCorrDK(1,2) = 2d0/3d0;  ColCorrDK(2,1) = 2d0/3d0????????/why is there a color matrix

      call TopDecay(Top(0,0,1),DK_LO,MomExt(1:4,10:12))


ELSEIF( Contr.eq.16 ) THEN!     W -> (ud+uu) and (cs+cc)
      if( TopDecays.ne.2 .and. TopDecays.ne.4 ) cycle
      TopCorr = ATop_
      call EvalPhasespace_TopDK(T_B_WGG,MomExt(1:4,3),yRnd(5:14),MomExt(1:4,5:9),PSWgt2);     PSWgt2 = PSWgt2 * 2d0 * SymmFac! this accounts for 2 quark flavor combinations and symm.fact.
      call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,4),yRnd(15:18),MomExt(1:4,10:12),PSWgt3)
      PSWgt = PSWgt1*PSWgt2*PSWgt3
      if( PSWgt.eq.0d0 ) cycle
      call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
      if( applyPSCut ) goto TheDipolesREDKJ

      iPrimAmpDK=1
      call TopDecay(ATop(plus_, 0,iPrimAmpDK),DKJ_RE2_QQ,MomExt(1:4,5:9),GluonHel=int((-1)**plus_));   ATop(plus_, 1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)
      call TopDecay(ATop(minus_,0,iPrimAmpDK),DKJ_RE2_QQ,MomExt(1:4,5:9),GluonHel=int((-1)**minus_)); ATop(minus_,1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)

      iPrimAmpDK=2
      call swapMom(MomExt(1:4,6),MomExt(1:4,9))
      call TopDecay(ATop(plus_, 0,iPrimAmpDK),DKJ_RE2_QQ,MomExt(1:4,5:9),GluonHel=int((-1)**plus_));   ATop(plus_, 1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)
      call TopDecay(ATop(minus_,0,iPrimAmpDK),DKJ_RE2_QQ,MomExt(1:4,5:9),GluonHel=int((-1)**minus_)); ATop(minus_,1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)
      call swapMom(MomExt(1:4,6),MomExt(1:4,9))
      ColCorrDK(1,1) = 2d0;      ColCorrDK(2,2) = 2d0
      ColCorrDK(1,2) = 2d0/3d0;  ColCorrDK(2,1) = 2d0/3d0

      call TopDecay(Top(0,0,1),DK_LO,MomExt(1:4,10:12))




ELSEIF( Contr.eq.17 ) THEN!     W -> (ud+dd) and (cs+ss)
      if( TopDecays.ne.2 .and. TopDecays.ne.4 ) cycle
      TopCorr = ATop_
      call EvalPhasespace_TopDK(T_B_WGG,MomExt(1:4,3),yRnd(5:14),MomExt(1:4,5:9),PSWgt2);     PSWgt2 = PSWgt2 * 2d0 * SymmFac! this accounts for 2 quark flavor combinations and symm.fact.
      call EvalPhasespace_TopDK(T_B_W, MomExt(1:4,4),yRnd(15:18),MomExt(1:4,10:12),PSWgt3)
      PSWgt = PSWgt1*PSWgt2*PSWgt3
      if( PSWgt.eq.0d0 ) cycle
      call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,10),MomExt(1:4,11),MomExt(1:4,12)/),applyPSCut,NBin)
      if( applyPSCut ) goto TheDipolesREDKJ

      iPrimAmpDK=1
      call TopDecay(ATop(plus_, 0,iPrimAmpDK),DKJ_RE2_QQ,MomExt(1:4,5:9),GluonHel=int((-1)**plus_));  ATop(plus_, 1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)
      call TopDecay(ATop(minus_,0,iPrimAmpDK),DKJ_RE2_QQ,MomExt(1:4,5:9),GluonHel=int((-1)**minus_)); ATop(minus_,1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)

      iPrimAmpDK=2
      call swapMom(MomExt(1:4,7),MomExt(1:4,8))
      call TopDecay(ATop(plus_, 0,iPrimAmpDK),DKJ_RE2_QQ,MomExt(1:4,5:9),GluonHel=int((-1)**plus_));  ATop(plus_, 1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)
      call TopDecay(ATop(minus_,0,iPrimAmpDK),DKJ_RE2_QQ,MomExt(1:4,5:9),GluonHel=int((-1)**minus_)); ATop(minus_,1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)
      call swapMom(MomExt(1:4,7),MomExt(1:4,8))
      ColCorrDK(1,1) = 2d0;      ColCorrDK(2,2) = 2d0
      ColCorrDK(1,2) = 2d0/3d0;  ColCorrDK(2,1) = 2d0/3d0

      call TopDecay(Top(0,0,1),DK_LO,MomExt(1:4,10:12))




ELSEIF( Contr.eq.18 ) THEN!     W -> (ud+cc/ss/bb) and (cs+uu/dd/bb)
      if( TopDecays.ne.2 .and. TopDecays.ne.3 ) cycle
      TopCorr = Top_
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,3),yRnd(5:8),MomExt(1:4,5:7),PSWgt2)
      call EvalPhasespace_TopDK(T_B_WGG, MomExt(1:4,4),yRnd(9:18),MomExt(1:4,8:12),PSWgt3);     PSWgt3 = PSWgt3 * 6d0! this accounts for 6 quark flavor combinations
      PSWgt = PSWgt1*PSWgt2*PSWgt3
      if( PSWgt.eq.0d0 ) cycle
      call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,11),MomExt(1:4,12),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,10)/),applyPSCut,NBin)
      if( applyPSCut ) goto TheDipolesREDKJ
      call TopDecay(ATop(0,0,1),DK_LO,MomExt(1:4,5:7))
      call TopDecay(Top(plus_, 0,1),DKJ_RE2_QQ,MomExt(1:4,8:12),GluonHel=int((-1)**plus_));  Top(plus_, 1,1)%Pol(1:4) = (0d0,0d0)
      call TopDecay(Top(minus_,0,1),DKJ_RE2_QQ,MomExt(1:4,8:12),GluonHel=int((-1)**minus_)); Top(minus_,1,1)%Pol(1:4) = (0d0,0d0)
      ColCorrDK(1,1) = 2d0;
!       ColCorrDK(2,2) = 2d0
!       ColCorrDK(1,2) = 2d0/3d0;  ColCorrDK(2,1) = 2d0/3d0




ELSEIF( Contr.eq.19 ) THEN!     W -> (ud+dd) and (cs+ss)
      if( TopDecays.ne.2 .and. TopDecays.ne.3 ) cycle
      TopCorr = Top_
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,3),yRnd(5:8),MomExt(1:4,5:7),PSWgt2)
      call EvalPhasespace_TopDK(T_B_WGG, MomExt(1:4,4),yRnd(9:18),MomExt(1:4,8:12),PSWgt3);     PSWgt3 = PSWgt3 * 2d0 * SymmFac! this accounts for 2 quark flavor combinations and symm.fact.
      PSWgt = PSWgt1*PSWgt2*PSWgt3
      if( PSWgt.eq.0d0 ) cycle
      call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,11),MomExt(1:4,12),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,10)/),applyPSCut,NBin)
      if( applyPSCut ) goto TheDipolesREDKJ
      call TopDecay(ATop(0,0,1),DK_LO,MomExt(1:4,5:7))


      iPrimAmpDK=1
      call TopDecay(Top(plus_, 0,iPrimAmpDK),DKJ_RE2_QQ,MomExt(1:4,8:12),GluonHel=int((-1)**plus_));  Top(plus_, 1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)
      call TopDecay(Top(minus_,0,iPrimAmpDK),DKJ_RE2_QQ,MomExt(1:4,8:12),GluonHel=int((-1)**minus_)); Top(minus_,1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)

      iPrimAmpDK=2
      call swapMom(MomExt(1:4,9),MomExt(1:4,12))
      call TopDecay(Top(plus_, 0,iPrimAmpDK),DKJ_RE2_QQ,MomExt(1:4,8:12),GluonHel=int((-1)**plus_));  Top(plus_, 1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)
      call TopDecay(Top(minus_,0,iPrimAmpDK),DKJ_RE2_QQ,MomExt(1:4,8:12),GluonHel=int((-1)**minus_)); Top(minus_,1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)
      call swapMom(MomExt(1:4,9),MomExt(1:4,12))
      ColCorrDK(1,1) = 2d0;      ColCorrDK(2,2) = 2d0
      ColCorrDK(1,2) = 2d0/3d0;  ColCorrDK(2,1) = 2d0/3d0



ELSEIF( Contr.eq.20 ) THEN!     W -> (ud+uu) and (cs+cc)
      if( TopDecays.ne.2 .and. TopDecays.ne.3 ) cycle
      TopCorr = Top_
      call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,3),yRnd(5:8),MomExt(1:4,5:7),PSWgt2)
      call EvalPhasespace_TopDK(T_B_WGG,MomExt(1:4,4),yRnd(9:18),MomExt(1:4,8:12),PSWgt3);     PSWgt3 = PSWgt3 * 2d0 * SymmFac! this accounts for 2 quark flavor combinations and symm.fact.
      PSWgt = PSWgt1*PSWgt2*PSWgt3
      if( PSWgt.eq.0d0 ) cycle
      call Kinematics_TTBARJET(1,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,11),MomExt(1:4,12),MomExt(1:4,3),MomExt(1:4,4)/),(/MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,10)/),applyPSCut,NBin)
      if( applyPSCut ) goto TheDipolesREDKJ
      call TopDecay(ATop(0,0,1),DK_LO,MomExt(1:4,5:7))

      iPrimAmpDK=1
      call TopDecay(Top(plus_, 0,iPrimAmpDK),DKJ_RE2_QQ,MomExt(1:4,8:12),GluonHel=int((-1)**plus_));  Top(plus_, 1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)
      call TopDecay(Top(minus_,0,iPrimAmpDK),DKJ_RE2_QQ,MomExt(1:4,8:12),GluonHel=int((-1)**minus_)); Top(minus_,1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)

      iPrimAmpDK=2
      call swapMom(MomExt(1:4,10),MomExt(1:4,11))
      call TopDecay(Top(plus_, 0,iPrimAmpDK),DKJ_RE2_QQ,MomExt(1:4,8:12),GluonHel=int((-1)**plus_));  Top(plus_, 1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)
      call TopDecay(Top(minus_,0,iPrimAmpDK),DKJ_RE2_QQ,MomExt(1:4,8:12),GluonHel=int((-1)**minus_)); Top(minus_,1,iPrimAmpDK)%Pol(1:4) = (0d0,0d0)
      call swapMom(MomExt(1:4,10),MomExt(1:4,11))
      ColCorrDK(1,1) = 2d0;      ColCorrDK(2,2) = 2d0
      ColCorrDK(1,2) = 2d0/3d0;  ColCorrDK(2,1) = 2d0/3d0



ENDIF! Contr



! this include files contains cross-check of the spinors against MadGraph
! #INCLUDE "./misc/CHECK_REDKJ_SPINORS.f90"



      PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
      NLO_Res_Unpol = (0d0,0d0)
      do iHel=1,NumHelicities ! loop over polarizations
      do GluHel1=plus_,minus_      ! loop over additional gluon polarization
      do GluHel2=plus_,minus_      ! loop over additional gluon polarization
            call HelCrossing(Helicities(iHel,1:NumExtParticles))
            call SetPolarizations()

            do iPrimAmpDK=1,NumDKPrimAmp(Contr)
                    if( TopCorr.eq.ATop_ ) then
                          ExtParticle(1)%Pol(1:4) = ATop(GluHel1,GluHel2,iPrimAmpDK)%Pol(1:4)
                          ExtParticle(2)%Pol(1:4) =  Top(0,0,1)%Pol(1:4)
                    elseif( TopCorr.eq.Top_ ) then
                          ExtParticle(1)%Pol(1:4) = ATop(0,0,1)%Pol(1:4)
                          ExtParticle(2)%Pol(1:4) =  Top(GluHel1,GluHel2,iPrimAmpDK)%Pol(1:4)
                    elseif( TopCorr.eq.ATop_+Top_ ) then
                          ExtParticle(1)%Pol(1:4) = ATop(GluHel1,0,1)%Pol(1:4)
                          ExtParticle(2)%Pol(1:4) =  Top(GluHel2,0,1)%Pol(1:4)
                    endif
                    do iPrimAmp=1,NumBornAmps
                      call EvalTree(BornAmps(iPrimAmp))
                      TreeResult(iPrimAmp,iPrimAmpDK-1) = BornAmps(iPrimAmp)%Result
                    enddo
            enddo

            HelCorr(1:2,1:2) = 1d0
            if( (Contr.eq.12 .or. Contr.eq.14).and.GluHel1.eq.minus_ ) then! reject 1-2 interference term for right-handed quarks
                  HelCorr(1,2) = 0d0
                  HelCorr(2,1) = 0d0
            endif

            NLO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
               do iPrimAmp=1,NumBornAmps
                  do iPrimAmpDK=1,NumDKPrimAmp(Contr)
                      do jPrimAmpDK=1,NumDKPrimAmp(Contr)
                          NLO_Res_Pol = NLO_Res_Pol + HelCorr(iPrimAmpDK,jPrimAmpDK) * ColCorrDK(iPrimAmpDK,jPrimAmpDK) * ColCorrLO(iPrimAmp,jPrimAmp)  &
                                                    * TreeResult(jPrimAmp,jPrimAmpDK-1) * dconjg(TreeResult(iPrimAmp,iPrimAmpDK-1))
                      enddo
                  enddo
               enddo
            enddo
            NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
      enddo!helicity loop
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac


      EvalCS_REDKJ_ttb = EvalCS_REDKJ_ttb + dble(NLO_Res_Unpol)
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),dreal(NLO_Res_UnPol))
      enddo

dummy(1) = dummy(1) + NLO_Res_Unpol
! print *, "Real",Contr,NLO_Res_Unpol
! pause

if( isnan(cdabs(NLO_Res_UnPol)) ) then
      print *, "Real:",NLO_Res_Unpol
      call printyrnd(yrnd(1:18))
      stop
endif




!------------------------------------------------------------
TheDipolesREDKJ continue!      dipole spinors               |
!------------------------------------------------------------
! dummy(2)=0d0
PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac


DO nDipole=1,NumDipoles(Contr);
          call DipolesDKJ(Contr,nDipole,ATop,Top,TopCorr,MomExt(1:4,1:12),Dipole,applyPSCut,NBin)
          if( applyPSCut ) cycle
          NLO_Res_Unpol = (0d0,0d0)
          do iHel=1,NumHelicities ! loop over polarizations
          do GluHel1=plus_,minus_      ! loop over additional gluon polarization
                call HelCrossing(Helicities(iHel,1:NumExtParticles))
                call SetPolarizations()
                if( TopCorr.eq.ATop_ ) then
                    ExtParticle(1)%Pol(1:4) = ATop(GluHel1,plus_,1)%Pol(1:4)!  hard jet
                    ExtParticle(2)%Pol(1:4) =  Top(1,1,1)%Pol(1:4)!            reduced process
                elseif( TopCorr.eq.Top_ ) then
                    ExtParticle(1)%Pol(1:4) = ATop(1,1,1)%Pol(1:4)!            reduced process
                    ExtParticle(2)%Pol(1:4) =  Top(GluHel1,plus_,1)%Pol(1:4)!  hard jet
                endif

                do iPrimAmp=1,NumBornAmps
                  call EvalTree(BornAmps(iPrimAmp))
                  TreeResult(iPrimAmp,GluHel1) = BornAmps(iPrimAmp)%Result
                enddo
          enddo!GluHel1 helicity loop
          NLO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,NumBornAmps
              do iPrimAmp=1,NumBornAmps
                      NLO_Res_Pol = NLO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * Dipole(plus_,plus_)   * TreeResult(jPrimAmp,plus_)*dconjg(TreeResult(iPrimAmp,plus_)) &
                                                + ColCorrLO(iPrimAmp,jPrimAmp) * Dipole(minus_,minus_) * TreeResult(jPrimAmp,minus_)*dconjg(TreeResult(iPrimAmp,minus_)) &
                                                + ColCorrLO(iPrimAmp,jPrimAmp) * Dipole(plus_,minus_)  * TreeResult(jPrimAmp,plus_)*dconjg(TreeResult(iPrimAmp,minus_)) &
                                                + ColCorrLO(iPrimAmp,jPrimAmp) * Dipole(minus_,plus_)  * TreeResult(jPrimAmp,minus_)*dconjg(TreeResult(iPrimAmp,plus_))
              enddo
          enddo
          NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
          enddo!helicity loop
          NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac

dummy(2) = dummy(2) + NLO_Res_Unpol
! print *, nDipole,dble(NLO_Res_Unpol)
! pause

          EvalCS_REDKJ_ttb = EvalCS_REDKJ_ttb + dble(NLO_Res_Unpol)
          do NHisto=1,NumHistograms
              call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_UnPol))
          enddo

if( isnan(cdabs(NLO_Res_UnPol)) ) then
      print *, "Dipole",NLO_Res_Unpol
      call printyrnd(yrnd(1:18))
      stop
endif
ENDDO! ndipole

ENDDO! Contr

! print *, "Real",dummy(1)
! print *, "Dipoles",dummy(2)
! print *, "ratio", dummy(1)/dummy(2), dummy(1)/dummy(2) +1d0
! print *, " Sing",(MomExt(1:4,5).dot.MomExt(1:4,8))/EHat**2
! print *, " Sing",(MomExt(1:4,5).dot.MomExt(1:4,9))/EHat**2
! print *, " Sing",(MomExt(1:4,8).dot.MomExt(1:4,9))/EHat**2
! pause



   EvalCS_REDKJ_ttb = EvalCS_REDKJ_ttb/VgsWgt
RETURN
END FUNCTION






END MODULE ModCrossSection_TTBJ
