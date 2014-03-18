MODULE ModCrossSection_TTBP
use ModTopDecay
implicit none

integer,private,parameter :: NumMaxHisto=45


integer,private,parameter :: nPhoRad1=1,nPhoRad2=2


contains




FUNCTION EvalCS_1L_ttbggp_MPI(yRnd,VgsWgt,res)
implicit none
integer :: EvalCS_1L_ttbggp_MPI
real(8) ::  yRnd(*),res(*),VgsWgt

res(1) = EvalCS_1L_ttbggp(yRnd,VgsWgt)
EvalCS_1L_ttbggp_MPI=0
RETURN
END FUNCTION






FUNCTION EvalCS_1L_ttbggp(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModIntDipoles_GGTTBGP
implicit none
real(8) ::  EvalCS_1L_ttbggp,yRnd(1:VegasMxDim),VgsWgt
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1),FermionLoopPartAmp(1:3,-2:1)
complex(8) :: BosonicPartAmp(1:3,-2:1),mydummy
integer :: iHel,jHel,kHel,iPrimAmp,jPrimAmp
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,ISFac
real(8) :: MomExt(1:4,1:12)
logical :: applyPSCut
real(8) :: MG_MOM(0:3,1:5)
real(8) :: MadGraph_tree
real(8),parameter :: Nc=3d0
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,AccPoles
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2),xE,HOp(1:3)
integer :: NBin(1:NumMaxHisto),NHisto,PhotonCouplCorr=2d0,nHel(1:2),NRndHel
! real(8) :: ThresholdCutOff = 1.0d0
include 'misc/global_import'
include 'vegas_common.f'


EvalCS_1L_ttbggp = 0d0
!print *, 'COMPARISON OF gg -> ttb+photon through (massive) fermions loops, in mod_CrosssSection_TTBP'

   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if( EHat.le.(2d0*m_Top+pT_pho_cut)  ) then
      EvalCS_1L_ttbggp = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to3(EHat,yRnd(3:7),MomExt(1:4,1:5),PSWgt)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))
   ISFac = MomCrossing(MomExt)

   NRndHel=8
IF( TOPDECAYS.NE.0 ) THEN
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomExt(1:4,6:8),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:15),.false.,MomExt(1:4,9:11),PSWgt3)
   PSWgt = PSWgt * PSWgt2*PSWgt3
   call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
   call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))
   NRndHel=16
ENDIF

   call Kinematics_TTBARPHOTON(0,MomExt(1:4,1:12),(/4,5,3,1,2,0,6,7,8,9,10,11/),applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_1L_ttbggp = 0d0
      return
   endif

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
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
   enddo!helicity loop


!------------ 1 LOOP --------------
ELSEIF( Correction.EQ.1 ) THEN

   do iHel=nHel(1),nHel(2)
!   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
!      print *, 'hel', iHel, Helicities(iHel,1:5)
      call SetPolarizations()
      do iPrimAmp=1,12
          call EvalTree(BornAmps(iPrimAmp))
      enddo

!       print *, "helicity ", ExtParticle(1)%Helicity,ExtParticle(2)%Helicity,ExtParticle(3)%Helicity,ExtParticle(4)%Helicity,ExtParticle(5)%Helicity
!       print *, "-----"
!       print *, "|Ta*Tb*M_0|^2:", PhotonCouplCorr* ColLO_ttbgg(1,1) * BornAmps(1)%Result*dconjg(BornAmps(1)%Result)
!       print *, "|Tb*Ta*M_0|^2:", PhotonCouplCorr* ColLO_ttbgg(2,2) * BornAmps(2)%Result*dconjg(BornAmps(2)%Result)


      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
      do iPrimAmp=1,NumBornAmps
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol

!------------ bosonic loops --------------
       do iPrimAmp=1,12
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
!           call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
!DEC$ IF (_QuadPrecImpr==1)
          AccPoles = CheckPoles(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv(1:2))


! print *, "This is for the check of ttb+Z(=>photon)"
! print *, "P LO",iHel,cdabs( BornAmps(1)%Result     )
! print *, "P VI",iHel,cdabs( PrimAmps(1)%Result(-2) )
! print *, "P VI",iHel,cdabs( PrimAmps(1)%Result(-1) )
! print *, "P VI",iHel,cdabs( PrimAmps(1)%Result(+0) )
! print *, "P VI",iHel,cdabs( PrimAmps(1)%Result(+1) )
! pause


          if( AccPoles.gt.1d-4 ) then
!               coeff4_128(:,:) = qcmplx( coeff4(:,:) )!   not re-calculating penta,quad contributions
!               coeff5_128(:,:) = qcmplx( coeff5(:,:) )
              call PentCut_128(PrimAmps(iPrimAmp))
              call QuadCut_128(PrimAmps(iPrimAmp))
              call TripCut_128(PrimAmps(iPrimAmp))
              call DoubCut_128(PrimAmps(iPrimAmp))
              call SingCut_128(PrimAmps(iPrimAmp))
              call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
              call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen*2)
              PrimAmps(iPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(iPrimAmp)%Result(-2:1)
              AccPoles = CheckPoles(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv(1:2))
              if( AccPoles.gt.1d-3 ) then
                  print *, "SKIP",AccPoles
!                   call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
                  EvalCS_1L_ttbggp = 0d0
                  SkipCounter = SkipCounter + 1
                  return
              endif
          endif
!DEC$ ENDIF
      enddo

      BosonicPartAmp(1,-2:1) =    Nc   *  PrimAmps(PrimAmp1_15234)%Result(-2:1) &
                             - 1d0/Nc *(  PrimAmps(PrimAmp1_15432)%Result(-2:1) &
                             +            PrimAmps(PrimAmp1_14532)%Result(-2:1) &
                             +            PrimAmps(PrimAmp1_14352)%Result(-2:1) )


      BosonicPartAmp(2,-2:1) =    Nc   *  PrimAmps(PrimAmp1_15243)%Result(-2:1) &
                             - 1d0/Nc *(  PrimAmps(PrimAmp1_15342)%Result(-2:1) &
                             +            PrimAmps(PrimAmp1_13542)%Result(-2:1) &
                             +            PrimAmps(PrimAmp1_13452)%Result(-2:1) )

      BosonicPartAmp(3,-2:1) = PrimAmps(PrimAmp1_15324)%Result(-2:1) &
                             + PrimAmps(PrimAmp1_13524)%Result(-2:1) &
                             + PrimAmps(PrimAmp1_15423)%Result(-2:1) &
                             + PrimAmps(PrimAmp1_14523)%Result(-2:1) &
                             + PrimAmps(PrimAmp1_15234)%Result(-2:1) &
                             + PrimAmps(PrimAmp1_15243)%Result(-2:1) &
                             + PrimAmps(PrimAmp1_15432)%Result(-2:1) &
                             + PrimAmps(PrimAmp1_14532)%Result(-2:1) &
                             + PrimAmps(PrimAmp1_14352)%Result(-2:1) &
                             + PrimAmps(PrimAmp1_15342)%Result(-2:1) &
                             + PrimAmps(PrimAmp1_13542)%Result(-2:1) &
                             + PrimAmps(PrimAmp1_13452)%Result(-2:1)

! print *, ""
! print *, "1-loop amplitudes:"
! print *, "bosonic loops"
! print *, "|Ta*Tb|^2*Re[M0*cong(M1)]:", PhotonCouplCorr*ColLO_ttbgg(1,1)*dreal( BornAmps(1)%Result * dconjg( BosonicPartAmp(1,-2:1) ))
! print *, "|Tb*Ta|^2*Re[M0*cong(M1)]:", PhotonCouplCorr*ColLO_ttbgg(2,2)*dreal( BornAmps(2)%Result * dconjg( BosonicPartAmp(2,-2:1) ))
! print *, "check",cdabs( BosonicPartAmp(3,-2:1) )
! stop

      NLO_Res_Pol(-2:1) = (0d0,0d0)
      do jPrimAmp=1,3
      do iPrimAmp=1,NumBornAmps
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + Col1L_ttbggp(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(BosonicPartAmp(jPrimAmp,-2:1)) )
      enddo
      enddo
      NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1) + NLO_Res_Pol(-2:1)
!      print *,'Bosonic loops'
!      print *, 'color-summed DP', NLO_Res_Pol(-2)/LO_Res_Pol
!      print *, 'color-summed SP', NLO_Res_Pol(-1)/LO_Res_Pol
!      print *, 'color-summed FIN',(NLO_Res_Pol(0)+NLO_Res_Pol(1))
!      pause



! ------------ fermionic loops --------------
      do iPrimAmp=13,28
!       do iPrimAmp=13,13
          call SetKirill(PrimAmps(iPrimAmp))
!           print *, ""
!           print *, 'iPrimAmp',iPrimAmp,PrimAmps(iPrimAmp)%ExtLine(:)
          call PentCut(PrimAmps(iPrimAmp))
          call QuadCut(PrimAmps(iPrimAmp))
          call TripCut(PrimAmps(iPrimAmp))
          call DoubCut(PrimAmps(iPrimAmp))
          call SingCut(PrimAmps(iPrimAmp))
          call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen**2)
          call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)

! print *, "This is for the check of ttb+Z(=>photon)"
 !print *, "P LO",iHel,cdabs( BornAmps(1)%Result     )
!print *, "P LO",iHel,BornAmps(1)%Result     
! print *, "P VI",iHel,cdabs( PrimAmps(iPrimAmp)%Result(-2) )
! print *, "P VI",iHel,cdabs( PrimAmps(iPrimAmp)%Result(-1) )
! print *, "P VI",iHel,cdabs( PrimAmps(iPrimAmp)%Result(+0) )
! print *, "P VI",iHel,cdabs( PrimAmps(iPrimAmp)%Result(+1) )
!
! print *, "P VI",iHel, PrimAmps(iPrimAmp)%Result(-2)/BornAmps(1)%Result 
! print *, "P VI",iHel, PrimAmps(iPrimAmp)%Result(-1) /BornAmps(1)%Result 
! print *, "P VI",iHel, PrimAmps(iPrimAmp)%Result(+0) /BornAmps(1)%Result 
! print *, "P VI",iHel, PrimAmps(iPrimAmp)%Result(+1) /BornAmps(1)%Result 
! pause

          PrimAmps(iPrimAmp)%Result(-2:1) = -(0d0,1d0)*PrimAmps(iPrimAmp)%Result(-2:1) !minus is from closed fermion loop
!           call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,rdiv(2),rdiv(1))
!           call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
      enddo
      FermionLoopPartAmp(1,-2:1) = nf_light * PrimAmps(PrimAmp2_15234)%Result(-2:1)  &
                                 +            PrimAmps(PrimAmp2m_15234)%Result(-2:1) &
                                 + sumQlight/Q_top *( PrimAmps(PrimAmp2_12534)%Result(-2:1) &
                                                   + PrimAmps(PrimAmp2_12354)%Result(-2:1) &
                                                   + PrimAmps(PrimAmp2_12345)%Result(-2:1) &
                                                   + PrimAmps(PrimAmp2_15234)%Result(-2:1) &
                                                   ) &
                                 + PrimAmps(PrimAmp2m_12534)%Result(-2:1) &
                                 + PrimAmps(PrimAmp2m_12354)%Result(-2:1) &
                                 + PrimAmps(PrimAmp2m_12345)%Result(-2:1) &
                                 + PrimAmps(PrimAmp2m_15234)%Result(-2:1)

      FermionLoopPartAmp(2,-2:1) = nf_light * PrimAmps(PrimAmp2_15243)%Result(-2:1)  &
                                 +            PrimAmps(PrimAmp2m_15243)%Result(-2:1) &

                                 + sumQlight/Q_top *( PrimAmps(PrimAmp2_12543)%Result(-2:1) &
                                                   + PrimAmps(PrimAmp2_12453)%Result(-2:1) &
                                                   + PrimAmps(PrimAmp2_12435)%Result(-2:1) &
                                                   + PrimAmps(PrimAmp2_15243)%Result(-2:1) &
                                                   ) &
                                 + PrimAmps(PrimAmp2m_12543)%Result(-2:1) &
                                 + PrimAmps(PrimAmp2m_12453)%Result(-2:1) &
                                 + PrimAmps(PrimAmp2m_12435)%Result(-2:1) &
                                 + PrimAmps(PrimAmp2m_15243)%Result(-2:1)


      FermionLoopPartAmp(3,-2:1) = -1d0/Nc * ( nf_light * PrimAmps(PrimAmp2_15234)%Result(-2:1) *0d0  &    ! cancels with  nf_light * PrimAmps(PrimAmp2_15243)  below
                                 +                        PrimAmps(PrimAmp2m_15234)%Result(-2:1)*0d0 &     ! cancels with  PrimAmps(PrimAmp2m_15243)  below
                                 + sumQlight/Q_top *( PrimAmps(PrimAmp2_12534)%Result(-2:1) &
                                                   + PrimAmps(PrimAmp2_12354)%Result(-2:1) &
                                                   + PrimAmps(PrimAmp2_12345)%Result(-2:1) &
                                                   + PrimAmps(PrimAmp2_15234)%Result(-2:1)*0d0 &       ! cancels with PrimAmps(PrimAmp2_15243)  below
                                                   ) &
                                 + PrimAmps(PrimAmp2m_12534)%Result(-2:1) &
                                 + PrimAmps(PrimAmp2m_12354)%Result(-2:1) &
                                 + PrimAmps(PrimAmp2m_12345)%Result(-2:1) &
                                 + PrimAmps(PrimAmp2m_15234)%Result(-2:1) &
                                 +             nf_light * PrimAmps(PrimAmp2_15243)%Result(-2:1)*0d0  &
                                 +                        PrimAmps(PrimAmp2m_15243)%Result(-2:1)*0d0 &
                                 + sumQlight/Q_top *( PrimAmps(PrimAmp2_12543)%Result(-2:1) &
                                                   + PrimAmps(PrimAmp2_12453)%Result(-2:1) &
                                                   + PrimAmps(PrimAmp2_12435)%Result(-2:1) &
                                                   + PrimAmps(PrimAmp2_15243)%Result(-2:1)*0d0 &
                                                   ) &
                                 + PrimAmps(PrimAmp2m_12543)%Result(-2:1) &
                                 + PrimAmps(PrimAmp2m_12453)%Result(-2:1) &
                                 + PrimAmps(PrimAmp2m_12435)%Result(-2:1) &
                                 + PrimAmps(PrimAmp2m_15243)%Result(-2:1) &
                                 )


! standarized output for cross check
! particle ordering: tbar,t,glu,glu,pho
! DUW kinpoint in all outgoing kinematics
! write (*,"(A,I5)") ExtParticle(1)%Helicity,ExtParticle(2)%Helicity,ExtParticle(3)%Helicity,ExtParticle(4)%Helicity,ExtParticle(5)%Helicity




! print *, ""
! print *, "fermion loops"
! print *, "Ta*Tb contribution:"
! print *, "photon couples only to external top quark line:"
! print *, "massive fermion loop ",PrimAmps(PrimAmp2m_15234)%Result(-2:1)
! print *, "massless fermion loop ",nf_light*PrimAmps(PrimAmp2_15234)%Result(-2:1)
! print *, "photon couples to fermion loop"
! print *, "massive fermion loop ",PrimAmps(PrimAmp2m_12534)%Result(-2:1) &
!                                  + PrimAmps(PrimAmp2m_12354)%Result(-2:1) &
!                                  + PrimAmps(PrimAmp2m_12345)%Result(-2:1) &
!                                  + PrimAmps(PrimAmp2m_15234)%Result(-2:1)
! print *, "massless fermion loop ",sumQlight/Q_up *( PrimAmps(PrimAmp2_12534)%Result(-2:1) &
!                                                    + PrimAmps(PrimAmp2_12354)%Result(-2:1) &
!                                                    + PrimAmps(PrimAmp2_12345)%Result(-2:1) &
!                                                    + PrimAmps(PrimAmp2_15234)%Result(-2:1) &
!                                                    )

! print *, ""
! print *, "sum of everything",FermionLoopPartAmp(7,-2:1)
!
! print *, "check with Andreas"
! print *, "** all massive fermion loops"
! FermionLoopPartAmp(7,-2:1) = 0D0 * nf_light * PrimAmps(PrimAmp2_15234)%Result(-2:1)  &
!                                  +            PrimAmps(PrimAmp2m_15234)%Result(-2:1) &
!
!                                  + 0D0 * sumQlight/Q_up *( PrimAmps(PrimAmp2_12534)%Result(-2:1) &
!                                                    + PrimAmps(PrimAmp2_12354)%Result(-2:1) &
!                                                    + PrimAmps(PrimAmp2_12345)%Result(-2:1) &
!                                                    + PrimAmps(PrimAmp2_15234)%Result(-2:1) &
!                                                    ) &
!                                  + PrimAmps(PrimAmp2m_12534)%Result(-2:1) &
!                                  + PrimAmps(PrimAmp2m_12354)%Result(-2:1) &
!                                  + PrimAmps(PrimAmp2m_12345)%Result(-2:1) &
!                                  + PrimAmps(PrimAmp2m_15234)%Result(-2:1)
! print *, "abs result: ", cdabs( FermionLoopPartAmp(7,0)+FermionLoopPartAmp(7,1) )
!
! print *, "** all massless fermion loops, no nf,Qf terms"
! FermionLoopPartAmp(7,-2:1) =                  PrimAmps(PrimAmp2_15234)%Result(-2:1)  &
!                                  +  0D0*      PrimAmps(PrimAmp2m_15234)%Result(-2:1) &
!
!                                  +                ( PrimAmps(PrimAmp2_12534)%Result(-2:1) &
!                                                    + PrimAmps(PrimAmp2_12354)%Result(-2:1) &
!                                                    + PrimAmps(PrimAmp2_12345)%Result(-2:1) &
!                                                    + PrimAmps(PrimAmp2_15234)%Result(-2:1) &
!                                                    ) &
!                                  + 0D0 * PrimAmps(PrimAmp2m_12534)%Result(-2:1) &
!                                  + 0D0 * PrimAmps(PrimAmp2m_12354)%Result(-2:1) &
!                                  + 0D0 * PrimAmps(PrimAmp2m_12345)%Result(-2:1) &
!                                  + 0D0 * PrimAmps(PrimAmp2m_15234)%Result(-2:1)
! print *, "result: ", cdabs( FermionLoopPartAmp(7,0)+FermionLoopPartAmp(7,1) )
!
! print *, "** massive fermion loops with photon coupling to external top line"
! FermionLoopPartAmp(7,-2:1) = 0D0 * nf_light * PrimAmps(PrimAmp2_15234)%Result(-2:1)  &
!                                  +            PrimAmps(PrimAmp2m_15234)%Result(-2:1) &
!
!                                  + 0D0 * sumQlight/Q_up *( PrimAmps(PrimAmp2_12534)%Result(-2:1) &
!                                                    + PrimAmps(PrimAmp2_12354)%Result(-2:1) &
!                                                    + PrimAmps(PrimAmp2_12345)%Result(-2:1) &
!                                                    + PrimAmps(PrimAmp2_15234)%Result(-2:1) &
!                                                    ) &
!                                  +0D0 *  PrimAmps(PrimAmp2m_12534)%Result(-2:1) &
!                                  +0D0 *  PrimAmps(PrimAmp2m_12354)%Result(-2:1) &
!                                  +0D0 *  PrimAmps(PrimAmp2m_12345)%Result(-2:1) &
!                                  +0D0 *  PrimAmps(PrimAmp2m_15234)%Result(-2:1)
! print *, "result: ", cdabs( FermionLoopPartAmp(7,0)+FermionLoopPartAmp(7,1) )
!
! print *, "** massless fermion loops, no nf,Qf terms with photon coupling to external top line"
! FermionLoopPartAmp(7,-2:1) =                  PrimAmps(PrimAmp2_15234)%Result(-2:1)  &
!                                  +  0D0*      PrimAmps(PrimAmp2m_15234)%Result(-2:1) &
!
!                                  +   0D0*        ( PrimAmps(PrimAmp2_12534)%Result(-2:1) &
!                                                    + PrimAmps(PrimAmp2_12354)%Result(-2:1) &
!                                                    + PrimAmps(PrimAmp2_12345)%Result(-2:1) &
!                                                    + PrimAmps(PrimAmp2_15234)%Result(-2:1) &
!                                                    ) &
!                                  + 0D0 * PrimAmps(PrimAmp2m_12534)%Result(-2:1) &
!                                  + 0D0 * PrimAmps(PrimAmp2m_12354)%Result(-2:1) &
!                                  + 0D0 * PrimAmps(PrimAmp2m_12345)%Result(-2:1) &
!                                  + 0D0 * PrimAmps(PrimAmp2m_15234)%Result(-2:1)
! print *, "result: ", cdabs( FermionLoopPartAmp(7,0)+FermionLoopPartAmp(7,1) )
!
! print *, ""
! print *, "** massless fermion loops, no nf,Qf terms with photon coupling to fermion loop"
! FermionLoopPartAmp(7,-2:1) =        0D0*          PrimAmps(PrimAmp2_15234)%Result(-2:1)  &
!                                  +  0D0*      PrimAmps(PrimAmp2m_15234)%Result(-2:1) &
!                                  +                 ( PrimAmps(PrimAmp2_12534)%Result(-2:1) &
!                                                    + PrimAmps(PrimAmp2_12354)%Result(-2:1) &
!                                                    + PrimAmps(PrimAmp2_12345)%Result(-2:1) &
!                                                    + PrimAmps(PrimAmp2_15234)%Result(-2:1) &
!                                                    ) &
!                                  + 0D0 * PrimAmps(PrimAmp2m_12534)%Result(-2:1) &
!                                  + 0D0 * PrimAmps(PrimAmp2m_12354)%Result(-2:1) &
!                                  + 0D0 * PrimAmps(PrimAmp2m_12345)%Result(-2:1) &
!                                  + 0D0 * PrimAmps(PrimAmp2m_15234)%Result(-2:1)
! print *, "result: ", cdabs( FermionLoopPartAmp(7,0)+FermionLoopPartAmp(7,1) )


      NLO_Res_Pol(-2:1) = (0d0,0d0)
      do jPrimAmp=1,3
      do iPrimAmp=1,NumBornAmps
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + Col1L_ttbggp(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result*dconjg(FermionLoopPartAmp(jPrimAmp,-2:1)) )
      enddo
      enddo
!      print *, 'fermionic loops:'
!      print *, 'color-summed DP',NLO_Res_Pol(-2)/LO_Res_Pol
!      print *, 'color-summed SP',NLO_Res_Pol(-1)/LO_Res_Pol
!      print *, 'color-summed FIN',(NLO_Res_Pol(0)+NLO_Res_Pol(1))
!      PAUSE
      NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) + NLO_Res_Pol(-2:1)
   enddo! helicity loop
ENDIF


IF( Correction.EQ.0 ) THEN
!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * Q_top**2*alpha4Pi*PhotonCouplCorr * WidthExpansion
   EvalCS_1L_ttbggp = LO_Res_Unpol * PreFac

ELSEIF( Correction.EQ.1 ) THEN
!  overall normalization: (4*Pi)^eps/Gamma(1-eps)
!  CT contributions
                         ! beta        !top WFRC
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
   LO_Res_Unpol = LO_Res_Unpol                         * ISFac * (alpha_s4Pi*RunFactor)**2                            * Q_top**2*alpha4Pi*PhotonCouplCorr
   NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1)           * ISFac * (alpha_s4Pi*RunFactor)**2 * alpha_sOver2Pi*RunFactor * Q_top**2*alpha4Pi*PhotonCouplCorr
   NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) * ISFac * (alpha_s4Pi*RunFactor)**2 * alpha_sOver2Pi*RunFactor * Q_top**2*alpha4Pi*PhotonCouplCorr
   EvalCS_1L_ttbggp = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac


ELSEIF( CORRECTION.EQ.3 ) THEN

! print *, "1-loop eps2:",(NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2) )* PreFac,  (NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2))/(alpha_sOver2Pi*RunFactor*LO_Res_Unpol)
! print *, "1-loop eps1:",(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1) )* PreFac,  (NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/(alpha_sOver2Pi*RunFactor*LO_Res_Unpol)
! print *, "tree virt",LO_Res_Unpol/RunFactor**2

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   IF( TOPDECAYS.GE.1 ) THEN
       xE = yRnd(16+HelSampling)
   ELSEIF( TOPDECAYS.EQ.0 ) THEN
       xE = yRnd(8+HelSampling)
   ENDIF
   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

   call EvalIntDipoles_GGTTBGP((/MomExt(1:4,4),MomExt(1:4,3),MomExt(1:4,5),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,6:11),xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
   EvalCS_1L_ttbggp = HOp(1)    * pdf(0,1)  * pdf(0,2)   &
                    + HOp(2)/xE * pdf_z(0,1)* pdf(0,2)   &
                    + HOp(3)/xE * pdf(0,1)  * pdf_z(0,2)

! print *, "real singularitites",EvalCS_1L_ttbggp
! pause
ENDIF



! !      MADGRAPH CHECK: gg->ttbp, mt=172, alpha_s=0.13
!        MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!        MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!        MG_MOM(0:3,3) = MomExt(1:4,5)*100d0
!        MG_MOM(0:3,4) = MomExt(1:4,4)*100d0
!        MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!        call coupsm(0)
!        call SGG_TTBA(MG_MOM,MadGraph_tree)
!        print *, ""
!        print *, "My tree:         ", LO_Res_Unpol/(100d0)**2
!        print *, "MadGraph hel.amp:", MadGraph_tree
!        print *, "MG/ME ratio: ", MadGraph_tree/(dble(LO_Res_Unpol)/(100d0)**2)
!        pause


   if( IsNan(EvalCS_1L_ttbggp) ) then
        print *, "NAN:",EvalCS_1L_ttbggp
        print *, yRnd(:)
        print *, NLO_Res_UnPol(0),NLO_Res_UnPol(1),NLO_Res_UnPol_Ferm(0),NLO_Res_UnPol_Ferm(1)
        print *, PSWgt , VgsWgt , PDFFac, sHatJacobi
        print *, eta1,eta2,MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,:)
        print *, "SKIP EVENT!!!!!"
        EvalCS_1L_ttbggp = 0d0
        return
   endif


   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_ttbggp)
   enddo


   EvalCS_1L_ttbggp = EvalCS_1L_ttbggp/VgsWgt

RETURN
END FUNCTION








FUNCTION EvalCS_1L_ttbqqbp(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModIntDipoles_QQBTTBGP
use ModIntDipoles_QGTTBQP
use ModIntDipoles_QBGTTBQBP
implicit none
real(8) ::  EvalCS_1L_ttbqqbp,yRnd(1:VegasMxDim),VgsWgt,xE
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1)
complex(8) :: BosonicPartAmp(1:2,-2:1),FermionPartAmp(1:2,-2:1),mydummy(1:2),LOPartAmp(1:2)
integer :: iHel,iPrimAmp,jPrimAmp
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,ISFac,AccPoles,HOp(1:2,1:3),pdf_z(-6:6,1:2)
real(8) :: MomExt(1:4,1:12)
logical :: applyPSCut
real(8) :: MG_MOM(0:3,1:NumExtParticles)
real(8) :: MadGraph_tree
real(8),parameter :: Nc=3d0
real(8) :: tau,eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac_a(1:2),PDFFac_b(1:2),PDFFac(1:2),pdf(-6:6,1:2)
integer :: NHisto,NBin(1:NumMaxHisto),npdf,ParityFlip=1,PhotonCouplCorr=2d0,nHel(1:2),NRndHel
integer,parameter :: up=1,dn=2
include 'misc/global_import'
include "vegas_common.f"



  EvalCS_1L_ttbqqbp = 0d0
  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top+pT_pho_cut ) then
      EvalCS_1L_ttbqqbp = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to3(EHat,yRnd(3:7),MomExt(1:4,1:5),PSWgt)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))

   NRndHel=8
IF( TOPDECAYS.NE.0 ) THEN
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomExt(1:4,6:8),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:15),.false.,MomExt(1:4,9:11),PSWgt3)
   PSWgt = PSWgt * PSWgt2*PSWgt3
   NRndHel=16
ENDIF

   call Kinematics_TTBARPHOTON(0,MomExt(1:4,1:12),(/4,5,3,1,2,0,6,7,8,9,10,11/),applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_1L_ttbqqbp = 0d0
      return
   endif

   call setPDFs(eta1,eta2,MuFac,pdf)
   IF( PROCESS.EQ.22 ) THEN
      PDFFac_a(up) = pdf(Up_,1)*pdf(AUp_,2) + pdf(Chm_,1)*pdf(AChm_,2)
      PDFFac_a(dn) = pdf(Dn_,1)*pdf(ADn_,2) + pdf(Str_,1)*pdf(AStr_,2) + pdf(Bot_,1)*pdf(ABot_,2)
      PDFFac_b(up) = pdf(Up_,2)*pdf(AUp_,1) + pdf(Chm_,2)*pdf(AChm_,1)
      PDFFac_b(dn) = pdf(Dn_,2)*pdf(ADn_,1) + pdf(Str_,2)*pdf(AStr_,1) + pdf(Bot_,2)*pdf(ABot_,1)
   ENDIF


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
        PDFFac(1:2) = PDFFac_a(1:2)
!         PDFFac(1:2) = PDFFac_a(1:2)+PDFFac_b(1:2)
    elseif(npdf.eq.2) then
        PDFFac(1:2) = PDFFac_b(1:2)
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
    endif
    ISFac = MomCrossing(MomExt)
    IF( TOPDECAYS.GE.1 ) THEN
          call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
          call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))
    ENDIF
!     call InitCurrCache()
    call SetPropagators()

    do iHel=nHel(1),nHel(2)
!     do iHel=27,27; print *, "helicity 27"
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
        call SetPolarizations()
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
! ./TOPAZ Collider=1 TopDK=0 Process=22 Correction=0 NLOParam=1 ObsSet=21  VegasNc0=1000 VegasNc1=1000
! if( ExtParticle(3)%Helicity*ExtParticle(4)%Helicity.eq.+1 ) cycle
! print *, ExtParticle(1)%Helicity,ExtParticle(2)%Helicity,ExtParticle(3)%Helicity,ExtParticle(4)%Helicity,ExtParticle(5)%Helicity
! print *, BornAmps(iPrimAmp)%ExtLine
! print *, "BornAmp",iPrimAmp,BornAmps(iPrimAmp)%Result * dsqrt(2d0)
! pause
        enddo

        LOPartAmp(up) = BornAmps(1)%Result + Q_up/Q_top * BornAmps(2)%Result
        LOPartAmp(dn) = BornAmps(1)%Result + Q_dn/Q_top * BornAmps(2)%Result
        LO_Res_Pol = ColLO_ttbqqb(1,1) * ( LOPartAmp(up)*dconjg(LOPartAmp(up))*PDFFac(up) + LOPartAmp(dn)*dconjg(LOPartAmp(dn))*PDFFac(dn))


!       print *, "helicity ", ExtParticle(1)%Helicity,ExtParticle(2)%Helicity,ExtParticle(3)%Helicity,ExtParticle(4)%Helicity,ExtParticle(5)%Helicity
!       print *, "-----"
!       mydummy(1) = cdabs(BornAmps(1)%Result)*dsqrt(2d0)
!       mydummy(2) = cdabs((-1.246106094083781d-002,2.384433430790109d-002) +(-0.236268898844133d0,-7.208990082860921d-002))
!       print *, "Markus:",mydummy(1)
!       print *, "Andreas:",mydummy(2)
!       print *, "ratio1",mydummy(1)/mydummy(2)
!       print *, ""
!       mydummy(1) = cdabs(BornAmps(2)%Result)*dsqrt(2d0)
!       mydummy(2) = cdabs((-3.449927694153093d-002,0.137660370291406d0) +(-6.553216668898655d-002,-8.502282966201087d-002) )
!       print *, "Markus:",mydummy(1)
!       print *, "Andreas:",mydummy(2)
!       print *, "ratio2",mydummy(1)/mydummy(2)
!
!       mydummy(1) = (-1.246106094083781d-002,2.384433430790109d-002) +(-0.236268898844133d0,-7.208990082860921d-002)
!       mydummy(2) = (-3.449927694153093d-002,0.137660370291406d0) +(-6.553216668898655d-002,-8.502282966201087d-002)
!       print *, "sum ratio",cdabs(BornAmps(1)%Result+BornAmps(2)%Result)*dsqrt(2d0)/cdabs(mydummy(1)+mydummy(2))
!       print *, "sum with colors",LO_Res_Pol*2d0/(0.243307611897825d0)
!       pause

        LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
    enddo!helicity loop
  enddo! npdf loop
  call swapMom(MomExt(1:4,1),MomExt(1:4,2))   ! swap back to original order, for ID below
! print *, "mom swap deactivated"

!------------ 1 LOOP --------------
ELSEIF( CORRECTION.EQ.1 ) THEN
  do npdf=1,2
    if(npdf.eq.1) then
        PDFFac(1:2) = PDFFac_a(1:2)
!         PDFFac(1:2) = PDFFac_a(1:2) + PDFFac_b(1:2)
    elseif(npdf.eq.2) then
        PDFFac(1:2) = PDFFac_b(1:2)
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
    endif
    ISFac = MomCrossing(MomExt)

IF( TOPDECAYS.GE.1 ) THEN
      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
      call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))
ENDIF
!     call InitCurrCache()
    call SetPropagators()

   do iHel=nHel(1),nHel(2)
!     do iHel=27,27; print *, "helicity 27"
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      do iPrimAmp=1,10
          call EvalTree(BornAmps(iPrimAmp))
! ./TOPAZ Collider=1 TopDK=0 Process=22 Correction=0 NLOParam=1 ObsSet=21  VegasNc0=1000 VegasNc1=1000
! if( ExtParticle(3)%Helicity*ExtParticle(4)%Helicity.eq.+1 ) cycle
! print *, ExtParticle(1)%Helicity,ExtParticle(2)%Helicity,ExtParticle(3)%Helicity,ExtParticle(4)%Helicity,ExtParticle(5)%Helicity
! print *, BornAmps(iPrimAmp)%ExtLine
! print *, "BornAmp",iPrimAmp,BornAmps(iPrimAmp)%Result * dsqrt(2d0)
! pause
      enddo


        LOPartAmp(up) = BornAmps(1)%Result + Q_up/Q_top * BornAmps(2)%Result
        LOPartAmp(dn) = BornAmps(1)%Result + Q_dn/Q_top * BornAmps(2)%Result
        LO_Res_Pol = ColLO_ttbqqb(1,1) * ( LOPartAmp(up)*dconjg(LOPartAmp(up))*PDFFac(up) + LOPartAmp(dn)*dconjg(LOPartAmp(dn))*PDFFac(dn))
        LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol

!  PDFFac(up)=0d0; print *, "pdffactor up deactivated"

! ------------ bosonic loops --------------
      do iPrimAmp=1,10
          call SetKirill(PrimAmps(iPrimAmp))
          call PentCut(PrimAmps(iPrimAmp))
          call QuadCut(PrimAmps(iPrimAmp))
          call TripCut(PrimAmps(iPrimAmp))
          call DoubCut(PrimAmps(iPrimAmp))
          call SingCut(PrimAmps(iPrimAmp))
          call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)

! if( ExtParticle(3)%Helicity*ExtParticle(4)%Helicity.eq.+1 ) cycle
! print *, "PrimAmpType",PrimAmps(iPrimAmp)%AmpType
          call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen**2)
          PrimAmps(iPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(iPrimAmp)%Result(-2:1)
          call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,3,rdiv(2),rdiv(1))
!           call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
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
              AccPoles = CheckPoles(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv(1:2))
              if( AccPoles.gt.1d-3 ) then
                  print *, "SKIP",AccPoles
!                 call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
                  EvalCS_1L_ttbqqbp = 0d0
                  SkipCounter = SkipCounter + 1
                  return
              endif
          endif
!DEC$ ENDIF
      enddo

      BosonicPartAmp(up,-2:1)=  &
                             +   Nc * ( PrimAmps(PrimAmp1_15234)%Result(-2:1) + Q_up/Q_top*PrimAmps(PrimAmp1_12354)%Result(-2:1) )&
                             - 2d0/Nc*( PrimAmps(PrimAmp1_15234)%Result(-2:1) + Q_up/Q_top*PrimAmps(PrimAmp1_12354)%Result(-2:1) &
                                      + PrimAmps(PrimAmp1_15243)%Result(-2:1) - Q_up/Q_top*PrimAmps(PrimAmp1_12453)%Result(-2:1) )  &
                             - 1d0/Nc*( PrimAmps(PrimAmp3_15432)%Result(-2:1) + PrimAmps(PrimAmp3_14532)%Result(-2:1) + PrimAmps(PrimAmp3_14352)%Result(-2:1) &
                                      - Q_up/Q_top*PrimAmps(PrimAmp3_14532)%Result(-2:1) + PrimAmps(PrimAmp4_15234)%Result(-2:1)  &
                                      - Q_up/Q_top*(PrimAmps(PrimAmp4_12534)%Result(-2:1)+PrimAmps(PrimAmp4_12345)%Result(-2:1)+PrimAmps(PrimAmp4_15234)%Result(-2:1)) )

      BosonicPartAmp(dn,-2:1)=  &
                             +   Nc * ( PrimAmps(PrimAmp1_15234)%Result(-2:1) + Q_dn/Q_top*PrimAmps(PrimAmp1_12354)%Result(-2:1) )&
                             - 2d0/Nc*( PrimAmps(PrimAmp1_15234)%Result(-2:1) + Q_dn/Q_top*PrimAmps(PrimAmp1_12354)%Result(-2:1) &
                                      + PrimAmps(PrimAmp1_15243)%Result(-2:1) - Q_dn/Q_top*PrimAmps(PrimAmp1_12453)%Result(-2:1) )  &
                             - 1d0/Nc*( PrimAmps(PrimAmp3_15432)%Result(-2:1) + PrimAmps(PrimAmp3_14532)%Result(-2:1) + PrimAmps(PrimAmp3_14352)%Result(-2:1) &
                                      - Q_dn/Q_top*PrimAmps(PrimAmp3_14532)%Result(-2:1) + PrimAmps(PrimAmp4_15234)%Result(-2:1) &
                                      - Q_dn/Q_top*(PrimAmps(PrimAmp4_12534)%Result(-2:1)+PrimAmps(PrimAmp4_12345)%Result(-2:1)+PrimAmps(PrimAmp4_15234)%Result(-2:1)) )


      NLO_Res_Pol(-2:1) = Col1L_ttbqqb(1,1) *( dreal(LOPartAmp(up)*dconjg(BosonicPartAmp(up,-2:1)))*PDFFac(up) &
                                             + dreal(LOPartAmp(dn)*dconjg(BosonicPartAmp(dn,-2:1)))*PDFFac(dn) )
      NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1) + NLO_Res_Pol(-2:1)



! ------------ fermionic loops --------------
      do iPrimAmp=11,14
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

!     massless closed fermion loops
      FermionPartAmp(up,-2:1) = Nf_light*( PrimAmps(PrimAmp2_15234)%Result(-2:1) + Q_up/Q_top*PrimAmps(PrimAmp2_12354)%Result(-2:1) ) &
                                         + PrimAmps(PrimAmp2m_15234)%Result(-2:1) + Q_up/Q_top*PrimAmps(PrimAmp2m_12354)%Result(-2:1)


      FermionPartAmp(dn,-2:1) = Nf_light*( PrimAmps(PrimAmp2_15234)%Result(-2:1) + Q_dn/Q_top*PrimAmps(PrimAmp2_12354)%Result(-2:1) ) &
                                         + PrimAmps(PrimAmp2m_15234)%Result(-2:1) + Q_dn/Q_top*PrimAmps(PrimAmp2m_12354)%Result(-2:1)


      NLO_Res_Pol(-2:1) = Col1L_ttbqqb(1,1) *( dreal(LOPartAmp(up)*dconjg(FermionPartAmp(up,-2:1)))*PDFFac(up) &
                                             + dreal(LOPartAmp(dn)*dconjg(FermionPartAmp(dn,-2:1)))*PDFFac(dn) )
      NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1) + NLO_Res_Pol(-2:1)

! print *, LOPartAmp(dn)*dsqrt(2d0)*Q_top
! print *, "",NLO_Res_UnPol(-2)*2d0*Q_top**2/PDFFac(dn)  *2d0
! print *, "",NLO_Res_UnPol(-1)*2d0*Q_top**2/PDFFac(dn)  *2d0
! print *, "",(NLO_Res_UnPol(0)+NLO_Res_UnPol(1))*2d0*Q_top**2/PDFFac(dn)  *2d0

! print *, "check",(BosonicPartAmp(dn,-2)+FermionPartAmp(dn,-2))/LOPartAmp(dn)! = -2*C_F
! pause

! !     check with Andreas
!       Q_in=Q_up
!       mydummy(1) = LOPartAmp(1)
!       print *, "1.) tree level:",dcmplx(mydummy(1)*dsqrt(2d0))
!
!       mydummy(1) =BosonicPartAmp(dn,-2)+FermionPartAmp(dn,-2)
!       print *, "2.) Summe aller 1-loop Beiträge:",dcmplx(mydummy(1)*dsqrt(2d0))*Q_top
!       mydummy(1) =BosonicPartAmp(dn,-1)+FermionPartAmp(dn,-1)
!       print *, "2.) Summe aller 1-loop Beiträge:",dcmplx(mydummy(1)*dsqrt(2d0))*Q_top
!       mydummy(1) =BosonicPartAmp(dn,0)+FermionPartAmp(dn,0)+BosonicPartAmp(dn,1)+FermionPartAmp(dn,1)
!       print *, "2.) Summe aller 1-loop Beiträge:",dcmplx(mydummy(1)*dsqrt(2d0))*Q_top
!
!       mydummy(1) =BosonicPartAmp(1,-2)
!       print *, "3.) 2.) ohne geschlossene fermion loops:",dcmplx(mydummy(1)*dsqrt(2d0))
!       mydummy(1) =BosonicPartAmp(1,-1)
!       print *, "3.) 2.) ohne geschlossene fermion loops:",dcmplx(mydummy(1)*dsqrt(2d0))
!       mydummy(1) =BosonicPartAmp(1,0)+BosonicPartAmp(1,1)
!       print *, "3.) 2.) ohne geschlossene fermion loops:",dcmplx(mydummy(1)*dsqrt(2d0))
!
!
!       mydummy(1) = FermionPartAmp(1,-2)
!       print *, "4.) alle fermion loops:",dcmplx(mydummy(1)*dsqrt(2d0))
!       mydummy(1) = FermionPartAmp(1,-1)
!       print *, "4.) alle fermion loops:",dcmplx(mydummy(1)*dsqrt(2d0))
!       mydummy(1) = FermionPartAmp(1,0)+FermionPartAmp(1,1)
!       print *, "4.) alle fermion loops:",dcmplx(mydummy(1)*dsqrt(2d0))
!
!
!       mydummy(1) = Nc * PrimAmps(PrimAmp1_15234)%Result(-2)
!       print *, "5.) leading color (LC), Q_top term:",dcmplx(mydummy(1)*dsqrt(2d0))
!       mydummy(1) = Nc * PrimAmps(PrimAmp1_15234)%Result(-1)
!       print *, "5.) leading color (LC), Q_top term:",dcmplx(mydummy(1)*dsqrt(2d0))
!       mydummy(1) = Nc * (PrimAmps(PrimAmp1_15234)%Result(0)+PrimAmps(PrimAmp1_15234)%Result(1))
!       print *, "5.) leading color (LC), Q_top term:",dcmplx(mydummy(1)*dsqrt(2d0))
!
!
!       mydummy(1) =Nc * ( Q_in/Q_top*PrimAmps(PrimAmp1_12354)%Result(-2) )
!       print *, "6.) LC, Q_q term:",dcmplx(mydummy(1)*dsqrt(2d0))
!       mydummy(1) =Nc * ( Q_in/Q_top*PrimAmps(PrimAmp1_12354)%Result(-1) )
!       print *, "6.) LC, Q_q term:",dcmplx(mydummy(1)*dsqrt(2d0))
!       mydummy(1) =Nc * Q_in/Q_top*(PrimAmps(PrimAmp1_12354)%Result(0)+PrimAmps(PrimAmp1_12354)%Result(1) )
!       print *, "6.) LC, Q_q term:",dcmplx(mydummy(1)*dsqrt(2d0))
!
!
!       mydummy(1) =  - 2d0/Nc*( PrimAmps(PrimAmp1_15234)%Result(-2)  &
!                    + PrimAmps(PrimAmp1_15243)%Result(-2)  )  &
!                    - 1d0/Nc * ( PrimAmps(PrimAmp3_15432)%Result(-2) + PrimAmps(PrimAmp3_14532)%Result(-2) + PrimAmps(PrimAmp3_14352)%Result(-2)   &
!                     + PrimAmps(PrimAmp4_15234)%Result(-2) )
!       print *, "7.) subLC, Q_top term:",dcmplx(mydummy(1)*dsqrt(2d0))! 2.111381630660386E-002
!
!       mydummy(1) =  - 2d0/Nc*( PrimAmps(PrimAmp1_15234)%Result(-1)  &
!                    + PrimAmps(PrimAmp1_15243)%Result(-1)  )  &
!                    - 1d0/Nc * ( PrimAmps(PrimAmp3_15432)%Result(-1) + PrimAmps(PrimAmp3_14532)%Result(-1) + PrimAmps(PrimAmp3_14352)%Result(-1)   &
!                     + PrimAmps(PrimAmp4_15234)%Result(-1) )
!       print *, "7.) subLC, Q_top term:",dcmplx(mydummy(1)*dsqrt(2d0)) ! 0.140413498945346
!
!       mydummy(1) =  - 2d0/Nc*( PrimAmps(PrimAmp1_15234)%Result(0)  &
!                    + PrimAmps(PrimAmp1_15243)%Result(0)  )  &
!                    - 1d0/Nc * ( PrimAmps(PrimAmp3_15432)%Result(0) + PrimAmps(PrimAmp3_14532)%Result(0) + PrimAmps(PrimAmp3_14352)%Result(0)   &
!                     + PrimAmps(PrimAmp4_15234)%Result(0) ) &
!                   - 2d0/Nc*( PrimAmps(PrimAmp1_15234)%Result(1)  &
!                    + PrimAmps(PrimAmp1_15243)%Result(1)  )  &
!                    - 1d0/Nc * ( PrimAmps(PrimAmp3_15432)%Result(1) + PrimAmps(PrimAmp3_14532)%Result(1) + PrimAmps(PrimAmp3_14352)%Result(1)   &
!                     + PrimAmps(PrimAmp4_15234)%Result(1) )
!       print *, "7.) subLC, Q_top term:",dcmplx(mydummy(1)*dsqrt(2d0))! 0.656123011732464
!
!
!       mydummy(1) =  - 2d0/Nc*(  + Q_in/Q_top*PrimAmps(PrimAmp1_12354)%Result(-2) &
!                     + (-1D0)*Q_in/Q_top*PrimAmps(PrimAmp1_12453)%Result(-2) )  &
!                    - 1d0/Nc * ( - Q_in/Q_top*PrimAmps(PrimAmp3_14532)%Result(-2)  &
!                   +(-1D0)* Q_in/Q_top*(PrimAmps(PrimAmp4_12534)%Result(-2)+PrimAmps(PrimAmp4_12345)%Result(-2)+PrimAmps(PrimAmp4_15234)%Result(-2)) )
!       print *, "8.) subLC, Q_q term:",dcmplx(mydummy(1)*dsqrt(2d0))! 9.419615601729865E-003! 8.335952209979481E-003 +I*( -4.386463135293729E-003 )
!
!       mydummy(1) =  - 2d0/Nc*(  + Q_in/Q_top*PrimAmps(PrimAmp1_12354)%Result(-1) &
!                     + (-1D0)*Q_in/Q_top*PrimAmps(PrimAmp1_12453)%Result(-1) )  &
!                    - 1d0/Nc * ( - Q_in/Q_top*PrimAmps(PrimAmp3_14532)%Result(-1)  &
!                   + (-1D0)*Q_in/Q_top*(PrimAmps(PrimAmp4_12534)%Result(-1)+PrimAmps(PrimAmp4_12345)%Result(-1)+PrimAmps(PrimAmp4_15234)%Result(-1)) )
!       print *, "8.) subLC, Q_q term:",dcmplx(mydummy(1)*dsqrt(2d0))!,mydummy(1)*dsqrt(2d0)! 6.264345095599894E-002! 4.397233889470472E-002 +I*(  4.461653684236241E-002 )
!
!
!       mydummy(1) =  - 2d0/Nc*(  + Q_in/Q_top*PrimAmps(PrimAmp1_12354)%Result(0) &
!                     + (-1D0)*Q_in/Q_top*PrimAmps(PrimAmp1_12453)%Result(0) )  &
!                    - 1d0/Nc * ( - Q_in/Q_top*PrimAmps(PrimAmp3_14532)%Result(0)  &
!                   + (-1D0)*Q_in/Q_top*(PrimAmps(PrimAmp4_12534)%Result(0)+PrimAmps(PrimAmp4_12345)%Result(0)+PrimAmps(PrimAmp4_15234)%Result(0)) ) &
!                 - 2d0/Nc*(  + Q_in/Q_top*PrimAmps(PrimAmp1_12354)%Result(1) &
!                     + (-1D0)*Q_in/Q_top*PrimAmps(PrimAmp1_12453)%Result(1) )  &
!                    - 1d0/Nc * ( - Q_in/Q_top*PrimAmps(PrimAmp3_14532)%Result(1)  &
!                   + (-1D0)*Q_in/Q_top*(PrimAmps(PrimAmp4_12534)%Result(1)+PrimAmps(PrimAmp4_12345)%Result(1)+PrimAmps(PrimAmp4_15234)%Result(1)) )
!       print *, "8.) subLC, Q_q term:",dcmplx(mydummy(1)*dsqrt(2d0))!,mydummy(1)*dsqrt(2d0)! 0.495491782923673! -9.830428212092179E-002 +I*(  0.485642229487480      )
!
!
!       mydummy(1) = PrimAmps(PrimAmp2m_15234)%Result(-2)
!       print *, "9a.) top loop, Q_top term:",dcmplx(mydummy(1)*dsqrt(2d0))
!       mydummy(1) = PrimAmps(PrimAmp2m_15234)%Result(-1)
!       print *, "9a.) top loop, Q_top term:",dcmplx(mydummy(1)*dsqrt(2d0))
!       mydummy(1) = PrimAmps(PrimAmp2m_15234)%Result(0)+PrimAmps(PrimAmp2m_15234)%Result(1)
!       print *, "9a.) top loop, Q_top term:",dcmplx(mydummy(1)*dsqrt(2d0))
!
!
!       mydummy(1) = PrimAmps(PrimAmp2m_12354)%Result(-2)
!       print *, "9b.) top loop, Q_q term:",dcmplx(mydummy(1)*dsqrt(2d0))
!       mydummy(1) = PrimAmps(PrimAmp2m_12354)%Result(-1)
!       print *, "9b.) top loop, Q_q term:",dcmplx(mydummy(1)*dsqrt(2d0))
!       mydummy(1) = PrimAmps(PrimAmp2m_12354)%Result(0)+PrimAmps(PrimAmp2m_12354)%Result(1)
!       print *, "9b.) top loop, Q_q term:",dcmplx(mydummy(1)*dsqrt(2d0))
!
!
! !       mydummy(1) = PrimAmps(PrimAmp2m_12534)%Result(-2)+PrimAmps(PrimAmp2m_12345)%Result(-2)+PrimAmps(PrimAmp2m_15234)%Result(-2)+PrimAmps(PrimAmp2m_12354)%Result(-2)
! !       print *, "9c.) top loop, Q_f term:",dcmplx(mydummy(1)*dsqrt(2d0))
! !       mydummy(1) = PrimAmps(PrimAmp2m_12534)%Result(-1)+PrimAmps(PrimAmp2m_12345)%Result(-1)+PrimAmps(PrimAmp2m_15234)%Result(-1)+PrimAmps(PrimAmp2m_12354)%Result(-1)
! !       print *, "9c.) top loop, Q_f term:",dcmplx(mydummy(1)*dsqrt(2d0))
! !       mydummy(1) = PrimAmps(PrimAmp2m_12534)%Result(0)+PrimAmps(PrimAmp2m_12345)%Result(0)+PrimAmps(PrimAmp2m_15234)%Result(0)+PrimAmps(PrimAmp2m_12354)%Result(0)+ PrimAmps(PrimAmp2m_12534)%Result(1)+PrimAmps(PrimAmp2m_12345)%Result(1)+PrimAmps(PrimAmp2m_15234)%Result(1)+PrimAmps(PrimAmp2m_12354)%Result(1)
! !       print *, "9c.) top loop, Q_f term:",dcmplx(mydummy(1)*dsqrt(2d0))
!
!
!
!       mydummy(1) =  PrimAmps(PrimAmp2_15234)%Result(-2)
!       print *, "10a.) light loops, Q_top term:",dcmplx(mydummy(1)*dsqrt(2d0)) * Nf_light
!       mydummy(1) =  PrimAmps(PrimAmp2_15234)%Result(-1)
!       print *, "10a.) light loops, Q_top term:",dcmplx(mydummy(1)*dsqrt(2d0)) * Nf_light
!       mydummy(1) =  PrimAmps(PrimAmp2_15234)%Result(0)+PrimAmps(PrimAmp2_15234)%Result(1)
!       print *, "10a.) light loops, Q_top term:",dcmplx(mydummy(1)*dsqrt(2d0)) * Nf_light
!
!
!       mydummy(1) =  PrimAmps(PrimAmp2_12354)%Result(-2)
!       print *, "10b.) light loops, Q_q term:",dcmplx(mydummy(1)*dsqrt(2d0)) * Nf_light
!       mydummy(1) =  PrimAmps(PrimAmp2_12354)%Result(-1)
!       print *, "10b.) light loops, Q_q term:",dcmplx(mydummy(1)*dsqrt(2d0)) * Nf_light
!       mydummy(1) =  PrimAmps(PrimAmp2_12354)%Result(0)+PrimAmps(PrimAmp2_12354)%Result(1)
!       print *, "10b.) light loops, Q_q term:",dcmplx(mydummy(1)*dsqrt(2d0)) * Nf_light
!
!
! !       mydummy(1) =  PrimAmps(PrimAmp2_12534)%Result(-2)+PrimAmps(PrimAmp2_12345)%Result(-2)+PrimAmps(PrimAmp2_15234)%Result(-2)+PrimAmps(PrimAmp2_12354)%Result(-2)
! !       print *, "10c.) light loops, summme über Q_f terms:",dcmplx(mydummy(1)*dsqrt(2d0))*sumQlight
! !       mydummy(1) =  PrimAmps(PrimAmp2_12534)%Result(-1)+PrimAmps(PrimAmp2_12345)%Result(-1)+PrimAmps(PrimAmp2_15234)%Result(-1)+PrimAmps(PrimAmp2_12354)%Result(-1)
! !       print *, "10c.) light loops, summme über Q_f terms:",dcmplx(mydummy(1)*dsqrt(2d0))*sumQlight
! !       mydummy(1) =  PrimAmps(PrimAmp2_12534)%Result(0)+PrimAmps(PrimAmp2_12345)%Result(0)+PrimAmps(PrimAmp2_15234)%Result(0)+PrimAmps(PrimAmp2_12354)%Result(0) &
! !                  +  PrimAmps(PrimAmp2_12534)%Result(1)+PrimAmps(PrimAmp2_12345)%Result(1)+PrimAmps(PrimAmp2_15234)%Result(1)+PrimAmps(PrimAmp2_12354)%Result(1)
! !       print *, "10c.) light loops, summme über Q_f terms:",dcmplx(mydummy(1)*dsqrt(2d0))*sumQlight
!       pause
   enddo!helicity loop
  enddo! npdf loop
  call swapMom(MomExt(1:4,1),MomExt(1:4,2))   ! swap back to original order, for ID below
! print *, "mom swap deactivated"
ENDIF




IF( CORRECTION.EQ.0 ) THEN
!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * Q_top**2*alpha4Pi*PhotonCouplCorr * WidthExpansion
   EvalCS_1L_ttbqqbp = LO_Res_Unpol * PreFac

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
   LO_Res_Unpol = LO_Res_Unpol                         * ISFac * (alpha_s4Pi*RunFactor)**2 * Q_top**2*alpha4Pi*PhotonCouplCorr
   NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1)           * ISFac * (alpha_s4Pi*RunFactor)**2 * Q_top**2*alpha4Pi*PhotonCouplCorr * alpha_sOver2Pi*RunFactor
   NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) * ISFac * (alpha_s4Pi*RunFactor)**2 * Q_top**2*alpha4Pi*PhotonCouplCorr * alpha_sOver2Pi*RunFactor


   EvalCS_1L_ttbqqbp = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac



ELSEIF( CORRECTION.EQ.3 ) THEN
! print *, "1-loop eps2:",(NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2) )* PreFac,  (NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2))/(alpha_sOver2Pi*RunFactor*LO_Res_Unpol)
! print *, "1-loop eps1:",(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1) )* PreFac,  (NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/(alpha_sOver2Pi*RunFactor*LO_Res_Unpol)
! print *, "tree virt",LO_Res_Unpol/RunFactor**2

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   IF( TOPDECAYS.GE.1 ) THEN
       xE = yRnd(16+HelSampling)
   ELSEIF( TOPDECAYS.EQ.0 ) THEN
       xE = yRnd(8+HelSampling)
   ENDIF

! xe=0.45d0;print *, "fixed xe!!"
   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

   IF( PROCESS.EQ.30 ) THEN
      call EvalIntDipoles_QQBTTBGP((/MomExt(1:4,4),MomExt(1:4,3),MomExt(1:4,5),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,6:11),xE,HOp(1:2,1:3))
      HOp(1:2,1:3) = HOp(1:2,1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ttbqqbp= HOp(1,1)    * (pdf(Up_,1)*pdf(AUp_,2)+pdf(Chm_,1)*pdf(AChm_,2) ) &
                        +HOp(1,2)/xE * (pdf_z(Up_,1)*pdf(AUp_,2)+pdf_z(Chm_,1)*pdf(AChm_,2) ) &
                        +HOp(1,3)/xE * (pdf(Up_,1)*pdf_z(AUp_,2)+pdf(Chm_,1)*pdf_z(AChm_,2) ) &
                        +HOp(2,1)    * (pdf(Dn_,1)*pdf(ADn_,2)+pdf(Str_,1)*pdf(AStr_,2)+pdf(Bot_,1)*pdf(ABot_,2) ) &
                        +HOp(2,2)/xE * (pdf_z(Dn_,1)*pdf(ADn_,2)+pdf_z(Str_,1)*pdf(AStr_,2)+pdf_z(Bot_,1)*pdf(ABot_,2) ) &
                        +HOp(2,3)/xE * (pdf(Dn_,1)*pdf_z(ADn_,2)+pdf(Str_,1)*pdf_z(AStr_,2)+pdf(Bot_,1)*pdf_z(ABot_,2) )

      call EvalIntDipoles_QQBTTBGP((/MomExt(1:4,4),MomExt(1:4,3),MomExt(1:4,5),-MomExt(1:4,2),-MomExt(1:4,1)/),MomExt(1:4,6:11),xE,HOp(1:2,1:3))
      HOp(1:2,1:3) = HOp(1:2,1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ttbqqbp= EvalCS_1L_ttbqqbp        &
                        +HOp(1,1)    * (pdf(Up_,2)*pdf(AUp_,1)+pdf(Chm_,2)*pdf(AChm_,1) ) &
                        +HOp(1,2)/xE * (pdf_z(Up_,2)*pdf(AUp_,1)+pdf_z(Chm_,2)*pdf(AChm_,1) ) &
                        +HOp(1,3)/xE * (pdf(Up_,2)*pdf_z(AUp_,1)+pdf(Chm_,2)*pdf_z(AChm_,1) ) &
                        +HOp(2,1)    * (pdf(Dn_,2)*pdf(ADn_,1)+pdf(Str_,2)*pdf(AStr_,1)+pdf(Bot_,2)*pdf(ABot_,1) ) &
                        +HOp(2,2)/xE * (pdf_z(Dn_,2)*pdf(ADn_,1)+pdf_z(Str_,2)*pdf(AStr_,1)+pdf_z(Bot_,2)*pdf(ABot_,1) ) &
                        +HOp(2,3)/xE * (pdf(Dn_,2)*pdf_z(ADn_,1)+pdf(Str_,2)*pdf_z(AStr_,1)+pdf(Bot_,2)*pdf_z(ABot_,1) )



   ELSEIF( PROCESS.EQ.24 ) THEN
      call EvalIntDipoles_QGTTBQP((/MomExt(1:4,4),MomExt(1:4,3),MomExt(1:4,5),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,6:11),xE,HOp(1:2,1:3))
      HOp(1:2,1:3) = HOp(1:2,1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ttbqqbp=HOp(1,1)    *  (pdf(Up_,1)*pdf(0,2)+pdf(Chm_,1)*pdf(0,2) ) &
                        +HOp(1,2)/xE * (pdf_z(Up_,1)*pdf(0,2)+pdf_z(Chm_,1)*pdf(0,2) ) &
                        +HOp(1,3)/xE * (pdf(Up_,1)*pdf_z(0,2)+pdf(Chm_,1)*pdf_z(0,2) ) &
                        +HOp(2,1)    * (pdf(Dn_,1)*pdf(0,2)+pdf(Str_,1)*pdf(0,2)+pdf(Bot_,1)*pdf(0,2) ) &
                        +HOp(2,2)/xE * (pdf_z(Dn_,1)*pdf(0,2)+pdf_z(Str_,1)*pdf(0,2)+pdf_z(Bot_,1)*pdf(0,2) ) &
                        +HOp(2,3)/xE * (pdf(Dn_,1)*pdf_z(0,2)+pdf(Str_,1)*pdf_z(0,2)+pdf(Bot_,1)*pdf_z(0,2) )

      call EvalIntDipoles_QGTTBQP((/MomExt(1:4,4),MomExt(1:4,3),MomExt(1:4,5),-MomExt(1:4,2),-MomExt(1:4,1)/),MomExt(1:4,6:11),xE,HOp(1:2,1:3))
      HOp(1:2,1:3) = HOp(1:2,1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ttbqqbp= EvalCS_1L_ttbqqbp  &
                        +HOp(1,1)    * (pdf(Up_,2)*pdf(0,1)+pdf(Chm_,2)*pdf(0,1) ) &
                        +HOp(1,2)/xE * (pdf_z(Up_,2)*pdf(0,1)+pdf_z(Chm_,2)*pdf(0,1) ) &
                        +HOp(1,3)/xE * (pdf(Up_,2)*pdf_z(0,1)+pdf(Chm_,2)*pdf_z(0,1) ) &
                        +HOp(2,1)    * (pdf(Dn_,2)*pdf(0,1)+pdf(Str_,2)*pdf(0,1)+pdf(Bot_,2)*pdf(0,1) ) &
                        +HOp(2,2)/xE * (pdf_z(Dn_,2)*pdf(0,1)+pdf_z(Str_,2)*pdf(0,1)+pdf_z(Bot_,2)*pdf(0,1) ) &
                        +HOp(2,3)/xE * (pdf(Dn_,2)*pdf_z(0,1)+pdf(Str_,2)*pdf_z(0,1)+pdf(Bot_,2)*pdf_z(0,1) )


   ELSEIF( PROCESS.EQ.26 ) THEN
      call EvalIntDipoles_QBGTTBQBP((/MomExt(1:4,4),MomExt(1:4,3),MomExt(1:4,5),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,6:11),xE,HOp(1:2,1:3))
      HOp(1:2,1:3) = HOp(1:2,1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ttbqqbp= HOp(1,1)    * (pdf(AUp_,1)*pdf(0,2)+pdf(AChm_,1)*pdf(0,2) ) &
                        +HOp(1,2)/xE * (pdf_z(AUp_,1)*pdf(0,2)+pdf_z(AChm_,1)*pdf(0,2) ) &
                        +HOp(1,3)/xE * (pdf(AUp_,1)*pdf_z(0,2)+pdf(AChm_,1)*pdf_z(0,2) ) &
                        +HOp(2,1)    * (pdf(ADn_,1)*pdf(0,2)+pdf(AStr_,1)*pdf(0,2)+pdf(ABot_,1)*pdf(0,2) ) &
                        +HOp(2,2)/xE * (pdf_z(ADn_,1)*pdf(0,2)+pdf_z(AStr_,1)*pdf(0,2)+pdf_z(ABot_,1)*pdf(0,2) ) &
                        +HOp(2,3)/xE * (pdf(ADn_,1)*pdf_z(0,2)+pdf(AStr_,1)*pdf_z(0,2)+pdf(ABot_,1)*pdf_z(0,2) )


      call EvalIntDipoles_QBGTTBQBP((/MomExt(1:4,4),MomExt(1:4,3),MomExt(1:4,5),-MomExt(1:4,2),-MomExt(1:4,1)/),MomExt(1:4,6:11),xE,HOp(1:2,1:3))
      HOp(1:2,1:3) = HOp(1:2,1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ttbqqbp= EvalCS_1L_ttbqqbp &
                        +HOp(1,1)    * (pdf(AUp_,2)*pdf(0,1)+pdf(AChm_,2)*pdf(0,1) ) &
                        +HOp(1,2)/xE * (pdf_z(AUp_,2)*pdf(0,1)+pdf_z(AChm_,2)*pdf(0,1) ) &
                        +HOp(1,3)/xE * (pdf(AUp_,2)*pdf_z(0,1)+pdf(AChm_,2)*pdf_z(0,1) ) &
                        +HOp(2,1)    * (pdf(ADn_,2)*pdf(0,1)+pdf(AStr_,2)*pdf(0,1)+pdf(ABot_,2)*pdf(0,1) ) &
                        +HOp(2,2)/xE * (pdf_z(ADn_,2)*pdf(0,1)+pdf_z(AStr_,2)*pdf(0,1)+pdf_z(ABot_,2)*pdf(0,1) ) &
                        +HOp(2,3)/xE * (pdf(ADn_,2)*pdf_z(0,1)+pdf(AStr_,2)*pdf_z(0,1)+pdf(ABot_,2)*pdf_z(0,1) )

   ENDIF

! print *, "real singularitites",EvalCS_1L_ttbqqbp
! print *, "real singularitites",EvalCS_1L_ttbqqbp/(alpha_sOver2Pi*RunFactor*LO_Res_Unpol* PreFac)
! pause

ENDIF



! !      MADGRAPH CHECK: uub->ttbp! mt=172, alpha_s=0.13
!        MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!        MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!        MG_MOM(0:3,3) = MomExt(1:4,5)*100d0
!        MG_MOM(0:3,4) = MomExt(1:4,4)*100d0
!        MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!        call coupsm(0)
!        call SUUB_TTBA(MG_MOM,MadGraph_tree)
!        print *, ""
!        print *, "My tree:         ", LO_Res_Unpol/PDFFac(up)/1d4
!        print *, "MadGraph hel.amp:", MadGraph_tree
!        print *, "MG/ME ratio: ", MadGraph_tree/dble(LO_Res_Unpol/PDFFac(up)/1d4)
!        pause


   if( IsNan(EvalCS_1L_ttbqqbp) ) then
        print *, "NAN:",EvalCS_1L_ttbqqbp
        print *, yRnd(:)
        print *, NLO_Res_UnPol(0),NLO_Res_UnPol(1),NLO_Res_UnPol_Ferm(0),NLO_Res_UnPol_Ferm(1)
        print *, PSWgt , VgsWgt , PDFFac_a,PDFFac_b, sHatJacobi
        print *, eta1,eta2,((1d0-eta1)*xE+eta1),((1d0-eta2)*xE+eta2),MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,:)
        print *, "SKIP EVENT!!!!!"
        EvalCS_1L_ttbqqbp = 0d0
        return
   endif


   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_ttbqqbp)
   enddo



   EvalCS_1L_ttbqqbp = EvalCS_1L_ttbqqbp/VgsWgt

RETURN
END FUNCTION








FUNCTION EvalCS_Real_ttbgggp(yRnd,VgsWgt)
use ModParameters
use ModKinematics
use ModAmplitudes
use ModMisc
use ModProcess
use ModDipoles_GGTTBGP
implicit none
real(8) ::  EvalCS_Real_ttbgggp,yRnd(1:VegasMxDim),VgsWgt,DipoleResult
complex(8) :: LO_Res_Pol,LO_Res_Unpol,PartAmp(1:4)
integer :: iHel,jPrimAmp,iPrimAmp,NHisto,NBin(1:NumMaxHisto)
real(8) :: EHat,PSWgt,PSWgt2,PSWgt3,ISFac,RunFactor,PreFac,sij
real(8) :: eta1,eta2,sHatJacobi,FluxFac,PDFFac
real(8) :: MomExt(1:4,1:12),pdf(-6:6,1:2)
real(8) :: MG_MOM(0:3,1:6)
real(8) :: MadGraph_tree
logical :: applyPSCut,applySingCut
real(8),parameter :: PhotonCouplCorr=2d0
include "vegas_common.f"


  EvalCS_Real_ttbgggp= 0d0
  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top+pT_pho_cut ) then
      EvalCS_Real_ttbgggp = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)


   call EvalPhaseSpace_2to4(EHat,yRnd(3:10),MomExt(1:4,1:6),PSWgt)
   call boost2Lab(eta1,eta2,6,MomExt(1:4,1:6))
   ISFac = MomCrossing(MomExt(1:4,1:6))

   PSWgt2 = 1d0
   PSWgt3 = 1d0
IF( TopDecays.GE.1 ) THEN
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(11:14),.false.,MomExt(1:4,7:9),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,6),yRnd(15:18),.false.,MomExt(1:4,10:12),PSWgt3)
   PSWgt = PSWgt * PSWgt2*PSWgt3
   call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,7:9))
   call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,10:12))
ENDIF

   call CheckSing(MomExt,applySingCut)
   if( applySingCut ) then
       EvalCS_Real_ttbgggp = 0d0
       return
   endif

   call Kinematics_TTBARPHOTON(1,MomExt(1:4,1:12),(/5,6,4,1,2,3,7,8,9,10,11,12/),applyPSCut,NBin)

   call setPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   RunFactor = RunAlphaS(NLOParam,MuRen)

   if( applyPSCut ) then
       EvalCS_Real_ttbgggp = 0d0
   else
        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities
          call HelCrossing(Helicities(iHel,1:NumExtParticles))
          call SetPolarizations()

          do iPrimAmp=1,NumBornAmps
              call EvalTree(BornAmps(iPrimAmp))
          enddo

          LO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,6
          do iPrimAmp=1,6
              LO_Res_Pol = LO_Res_Pol + ColLO_ttbggg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result * dconjg(BornAmps(jPrimAmp)%Result)
          enddo
          enddo
          LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop

        LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**3 * Q_top**2*alpha4Pi  *PhotonCouplCorr
        EvalCS_Real_ttbgggp = LO_Res_Unpol * PreFac

        do NHisto=1,NumHistograms
               call intoHisto(NHisto,NBin(NHisto),EvalCS_Real_ttbgggp)
        enddo
endif!applyPSCut


    PreFac = PreFac * ISFac * (alpha_s4Pi*RunFactor)**3 * Q_top**2*alpha4Pi*PhotonCouplCorr /PSWgt2/PSWgt3
    call EvalDipoles_GGTTBGP((/MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,6),-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,3)/),yRnd(11:18),PreFac,DipoleResult)

!      sij = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,3))
!      sij = MomExt(1,3)**2
!      print *,  sij/EHat**2,EvalCS_Real_ttbgggp,DipoleResult,(1d0+EvalCS_Real_ttbgggp/(DipoleResult))
!      pause

! !      MADGRAPH CHECK: gg->tbtgp
!        MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!        MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!        MG_MOM(0:3,3) = MomExt(1:4,5)*100d0
!        MG_MOM(0:3,4) = MomExt(1:4,6)*100d0
!        MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!        MG_MOM(0:3,6) = MomExt(1:4,4)*100d0
!        call coupsm(0)
!        call SGG_TBTGA(MG_MOM,MadGraph_tree)
!        LO_Res_Unpol = LO_Res_Unpol * 100d0**(-4)
!        print *, ""
!        print *, "My tree:         ", LO_Res_Unpol
!        print *, "MadGraph hel.amp:", MadGraph_tree
!        print *, "MG/ME ratios: ", MadGraph_tree/dble(LO_Res_Unpol)
!        pause

    EvalCS_Real_ttbgggp = (EvalCS_Real_ttbgggp + DipoleResult)/VgsWgt
RETURN
END FUNCTION







FUNCTION EvalCS_Real_ttbqqbgp(yRnd,VgsWgt)
use ModParameters
use ModKinematics
use ModAmplitudes
use ModProcess
use ModMisc
use ModDipoles_QQBTTBGP
use ModDipoles_QGTTBQP
use ModDipoles_QBGTTBQBP
implicit none
real(8) ::  EvalCS_Real_ttbqqbgp,EvalCS_Dips_ttbqqbgp,yRnd(1:VegasMxDim),VgsWgt,DipoleResult(1:2)
complex(8) :: LO_Res_Pol,LO_Res_Unpol,PartAmp(1:4)
integer :: iHel,jPrimAmp,iPrimAmp,NHisto,NBin(1:NumMaxHisto),NPDF
real(8) :: EHat,PSWgt,PSWgt2,PSWgt3,ISFac,RunFactor,PreFac,PreFacDip
real(8) :: eta1,eta2,sHatJacobi,FluxFac,PDFFac(1:2),PDFFac_a(1:2),PDFFac_b(1:2)
real(8) :: MomExt(1:4,1:12),sij,pdf(-6:6,1:2)
real(8) :: MG_MOM(0:3,1:6)
real(8) :: MadGraph_tree
logical :: applyPSCut,applySingCut
real(8),parameter :: PhotonCouplCorr=2d0
integer,parameter :: up=1, dn=2
include "vegas_common.f"


EvalCS_Real_ttbqqbgp= 0d0
EvalCS_Dips_ttbqqbgp= 0d0

   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if( EHat.le.2d0*m_Top+pT_pho_cut ) then
      EvalCS_Real_ttbqqbgp = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to4(EHat,yRnd(3:10),MomExt(1:4,1:6),PSWgt)
   call boost2Lab(eta1,eta2,6,MomExt(1:4,1:6))
   call CheckSing(MomExt,applySingCut)
   if( applySingCut ) then
       EvalCS_Real_ttbqqbgp = 0d0
       return
   endif

   PSWgt2 = 1d0
   PSWgt3 = 1d0
   IF( TopDecays.GE.1 ) THEN
       call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(11:14),.false.,MomExt(1:4,7:9),PSWgt2)
       call EvalPhasespace_TopDecay(MomExt(1:4,6),yRnd(15:18),.false.,MomExt(1:4,10:12),PSWgt3)
       PSWgt = PSWgt * PSWgt2*PSWgt3
   ENDIF
   call Kinematics_TTBARPHOTON(1,MomExt(1:4,1:12),(/5,6,4,1,2,3,7,8,9,10,11,12/),applyPSCut,NBin)

   call setPDFs(eta1,eta2,MuFac,pdf)
   IF( PROCESS.EQ.30 ) THEN
      PDFFac_a(up) = pdf(Up_,1)*pdf(AUp_,2) + pdf(Chm_,1)*pdf(AChm_,2)
      PDFFac_a(dn) = pdf(Dn_,1)*pdf(ADn_,2) + pdf(Str_,1)*pdf(AStr_,2) + pdf(Bot_,1)*pdf(ABot_,2)
      PDFFac_b(up) = pdf(Up_,2)*pdf(AUp_,1) + pdf(Chm_,2)*pdf(AChm_,1)
      PDFFac_b(dn) = pdf(Dn_,2)*pdf(ADn_,1) + pdf(Str_,2)*pdf(AStr_,1) + pdf(Bot_,2)*pdf(ABot_,1)
   ELSEIF( PROCESS.EQ.24 ) THEN
      PDFFac_a(up) = pdf(Up_,1)*pdf(0,2) + pdf(Chm_,1)*pdf(0,2)
      PDFFac_a(dn) = pdf(Dn_,1)*pdf(0,2) + pdf(Str_,1)*pdf(0,2) + pdf(Bot_,1)*pdf(0,2)
      PDFFac_b(up) = pdf(Up_,2)*pdf(0,1) + pdf(Chm_,2)*pdf(0,1)
      PDFFac_b(dn) = pdf(Dn_,2)*pdf(0,1) + pdf(Str_,2)*pdf(0,1) + pdf(Bot_,2)*pdf(0,1)
   ELSEIF( PROCESS.EQ.26 ) THEN
      PDFFac_a(up) = pdf(AUp_,1)*pdf(0,2) + pdf(AChm_,1)*pdf(0,2)
      PDFFac_a(dn) = pdf(ADn_,1)*pdf(0,2) + pdf(AStr_,1)*pdf(0,2) + pdf(ABot_,1)*pdf(0,2)
      PDFFac_b(up) = pdf(AUp_,2)*pdf(0,1) + pdf(AChm_,2)*pdf(0,1)
      PDFFac_b(dn) = pdf(ADn_,2)*pdf(0,1) + pdf(AStr_,2)*pdf(0,1) + pdf(ABot_,2)*pdf(0,1)
   ENDIF


! PDFFac_a(up) = 1d0; print *, "setting pdfs to one"
! PDFFac_a(dn) = 0d0
! PDFFac_b(up) = 0d0
! PDFFac_b(dn) = 0d0


   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   RunFactor = RunAlphaS(NLOParam,MuRen)
   DO NPDF=1,2
        if(npdf.eq.1) then
              PDFFac(1:2) = PDFFac_a(1:2)
!               PDFFac(1:2) = PDFFac_a(1:2) + PDFFac_b(1:2)
        elseif(npdf.eq.2) then
              PDFFac(1:2) = PDFFac_b(1:2)
              call swapMom(MomExt(1:4,1),MomExt(1:4,2))
        endif
        ISFac = MomCrossing(MomExt(1:4,1:6))
        IF( TopDecays.GE.1 ) THEN
          call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,7:9))
          call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,10:12))
        ENDIF


if( applyPSCut ) then
            EvalCS_Real_ttbqqbgp = 0d0
else
        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities
          call HelCrossing(Helicities(iHel,1:NumExtParticles))
          call SetPolarizations()

          do iPrimAmp=1,NumBornAmps
              call EvalTree(BornAmps(iPrimAmp))
          enddo

          Q_in = Q_up
          PartAmp(1) = BornAmps(PrimAmp1_162345)%Result + Q_in/Q_top*BornAmps(PrimAmp1_123645)%Result
          PartAmp(3) = BornAmps(PrimAmp1_162534)%Result + Q_in/Q_top*BornAmps(PrimAmp1_125364)%Result
          PartAmp(2) = BornAmps(PrimAmp1_156234)%Result + BornAmps(PrimAmp1_165234)%Result + Q_in/Q_top*BornAmps(PrimAmp1_152364)%Result
          PartAmp(4) = BornAmps(PrimAmp1_162354)%Result + Q_in/Q_top*BornAmps(PrimAmp1_123564)%Result + Q_in/Q_top*BornAmps(PrimAmp1_123654)%Result
          LO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,4
          do iPrimAmp=1,4
              LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqbg(iPrimAmp,jPrimAmp) * PartAmp(iPrimAmp) * dconjg(PartAmp(jPrimAmp))
          enddo
          enddo
          LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol * PDFFac(1)

          Q_in = Q_dn
          PartAmp(1) = BornAmps(PrimAmp1_162345)%Result + Q_in/Q_top*BornAmps(PrimAmp1_123645)%Result
          PartAmp(3) = BornAmps(PrimAmp1_162534)%Result + Q_in/Q_top*BornAmps(PrimAmp1_125364)%Result
          PartAmp(2) = BornAmps(PrimAmp1_156234)%Result + BornAmps(PrimAmp1_165234)%Result + Q_in/Q_top*BornAmps(PrimAmp1_152364)%Result
          PartAmp(4) = BornAmps(PrimAmp1_162354)%Result + Q_in/Q_top*BornAmps(PrimAmp1_123564)%Result + Q_in/Q_top*BornAmps(PrimAmp1_123654)%Result
          LO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,4
          do iPrimAmp=1,4
              LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqbg(iPrimAmp,jPrimAmp) * PartAmp(iPrimAmp) * dconjg(PartAmp(jPrimAmp))
          enddo
          enddo
          LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol * PDFFac(2)
        enddo!helicity loop


! print *, "real unpol p",LO_Res_UnPol* Q_top**2  *PhotonCouplCorr
! pause


        LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**3 * Q_top**2*alpha4Pi  *PhotonCouplCorr
        EvalCS_Real_ttbqqbgp = EvalCS_Real_ttbqqbgp + dble(LO_Res_Unpol*PreFac)

        do NHisto=1,NumHistograms
               call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol*PreFac))
        enddo
endif!applyPSCut

    PreFacDip = PreFac * ISFac * (alpha_s4Pi*RunFactor)**3 * Q_top**2*alpha4Pi*PhotonCouplCorr /PSWgt2/PSWgt3
    IF( PROCESS.EQ.30 ) THEN
        call EvalDipoles_QQBTTBGP((/MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,6),-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,3)/),yRnd(11:18),(/PreFacDip*PDFFac(1),PreFacDip*PDFFac(2)/),DipoleResult)
    ELSEIF( PROCESS.EQ.24 ) THEN
        call EvalDipoles_QGTTBQP((/MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,6),-MomExt(1:4,1),MomExt(1:4,3),-MomExt(1:4,2)/),yRnd(11:18),(/PreFacDip*PDFFac(1),PreFacDip*PDFFac(2)/),DipoleResult)
    ELSEIF( PROCESS.EQ.26 ) THEN
        call EvalDipoles_QBGTTBQBP((/MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,6),MomExt(1:4,3),-MomExt(1:4,1),-MomExt(1:4,2)/),yRnd(11:18),(/PreFacDip*PDFFac(1),PreFacDip*PDFFac(2)/),DipoleResult)
    ENDIF
    EvalCS_Dips_ttbqqbgp = EvalCS_Dips_ttbqqbgp + DipoleResult(1) + DipoleResult(2)

  ENDDO! loop over a<-->b pdfs


!      sij = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,3))
! !      sij = MomExt(1,3)**2
!      print *,  sij/EHat**2,EvalCS_Real_ttbqqbgp,EvalCS_Dips_ttbqqbgp,(1d0+EvalCS_Real_ttbqqbgp/EvalCS_Dips_ttbqqbgp)
!      pause


! !      MADGRAPH CHECK: ddb->ttbgp
!        MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!        MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!        MG_MOM(0:3,3) = MomExt(1:4,6)*100d0
!        MG_MOM(0:3,4) = MomExt(1:4,5)*100d0
!        MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!        MG_MOM(0:3,6) = MomExt(1:4,4)*100d0
!        call coupsm(0)
!        call SDDB_TTBGA(MG_MOM,MadGraph_tree)
!        LO_Res_Unpol = LO_Res_Unpol * 100d0**(-4)
!        print *, ""
!        print *, "My tree:         ", LO_Res_Unpol
!        print *, "MadGraph hel.amp:", MadGraph_tree
!        print *, "MG/ME ratio: ", MadGraph_tree/dble(LO_Res_Unpol)
!        pause


! !      MADGRAPH CHECK: ug->ttbup
!        MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!        MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!        MG_MOM(0:3,3) = MomExt(1:4,6)*100d0
!        MG_MOM(0:3,4) = MomExt(1:4,5)*100d0
!        MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!        MG_MOM(0:3,6) = MomExt(1:4,4)*100d0
!        call coupsm(0)
!        call SUG_TTBUA(MG_MOM,MadGraph_tree)
!        LO_Res_Unpol = LO_Res_Unpol * 100d0**(-4)
!        print *, ""
!        print *, "My tree:         ", LO_Res_Unpol
!        print *, "MadGraph hel.amp:", MadGraph_tree
!        print *, "MG/ME ratio: ", MadGraph_tree/dble(LO_Res_Unpol)
!        pause



! !      MADGRAPH CHECK: dbg->ttbdbp
!        MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!        MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!        MG_MOM(0:3,3) = MomExt(1:4,6)*100d0
!        MG_MOM(0:3,4) = MomExt(1:4,5)*100d0
!        MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!        MG_MOM(0:3,6) = MomExt(1:4,4)*100d0
!        call coupsm(0)
!        call SDBG_TTBDBA(MG_MOM,MadGraph_tree)
!        LO_Res_Unpol = LO_Res_Unpol * 100d0**(-4)
!        print *, ""
!        print *, "My tree:         ", LO_Res_Unpol
!        print *, "MadGraph hel.amp:", MadGraph_tree
!        print *, "MG/ME ratio: ", MadGraph_tree/dble(LO_Res_Unpol)
!        pause

    EvalCS_Real_ttbqqbgp = (EvalCS_Real_ttbqqbgp + EvalCS_Dips_ttbqqbgp) /VgsWgt

END FUNCTION








FUNCTION EvalCS_NLODK_ttbp(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModHadrWDecay
implicit none
real(8) ::  EvalCS_NLODK_ttbp,yRnd(1:VegasMxDim),VgsWgt
complex(8) :: LO_Res_Pol,LO_Res_Unpol,Dip_Res_Unpol,NLO_Res_Pol,NLO_Res_UnPol
complex(8) :: TreeResult(1:NumBornAmps),DKResult(1:NumBornAmps),LOPartAmp(1:2),NLOPartAmp(1:2)
integer :: iHel,jHel,kHel,GluHel,iPrimAmp,jPrimAmp,ndip
real(8) :: EHat,PSWgt1,PSWgt2,PSWgt3,ISFac,dip_res_w
real(8) :: MomExt(1:4,1:12),MomExtTd(1:4,1:12)
logical :: applyPSCut,applySingCut
real(8) :: tau,eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac_a(1:2),PDFFac_b(1:2),PDFFac(1:2),RunFactor
real(8) :: pdf(-6:6,1:2),ColCorrLO(1:NumBornAmps,1:NumBornAmps)
integer :: NBin(1:NumMaxHisto),NHisto,npdf
real(8) :: pbDpg,ptDpg,ptDpb,z,omz,Dipole,rsq,y
real(8), parameter :: CF=4d0/3d0,PhotonCouplCorr=2d0
real(8) :: MomBoost(1:4),MomLep1(1:4),MomLep2(1:4)
integer,parameter :: up=1,dn=2,glu=1
include "vegas_common.f"




  EvalCS_NLODK_ttbp = 0d0
  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top+pT_pho_cut ) then
      EvalCS_NLODK_ttbp = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

  call EvalPhaseSpace_2to3(EHat,yRnd(3:7),MomExt(1:4,1:5),PSWgt1)
  call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))

   call setPDFs(eta1,eta2,MuFac,pdf)
   IF( PROCESS.EQ.20 ) THEN
      PDFFac(glu)   = pdf(0,1) * pdf(0,2)
      PDFFac_a(glu) = PDFFac(glu)
      PDFFac_b(glu) = 0d0
      ColCorrLO(1:NumBornAmps,1:NumBornAmps) = ColLO_ttbgg(1:NumBornAmps,1:NumBornAmps)
   ELSEIF( PROCESS.EQ.22 ) THEN
      PDFFac_a(up) = pdf(Up_,1)*pdf(AUp_,2) + pdf(Chm_,1)*pdf(AChm_,2)
      PDFFac_a(dn) = pdf(Dn_,1)*pdf(ADn_,2) + pdf(Str_,1)*pdf(AStr_,2) + pdf(Bot_,1)*pdf(ABot_,2)
      PDFFac_b(up) = pdf(Up_,2)*pdf(AUp_,1) + pdf(Chm_,2)*pdf(AChm_,1)
      PDFFac_b(dn) = pdf(Dn_,2)*pdf(ADn_,1) + pdf(Str_,2)*pdf(AStr_,1) + pdf(Bot_,2)*pdf(ABot_,1)
      ColCorrLO(1:1,1:1) = ColLO_ttbqqb(1:1,1:1)
   ENDIF

IF( CORRECTION.EQ.4 ) THEN
!----------------------------------------
! one loop correction to Anti-top decay |
!----------------------------------------
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomExt(1:4,6:8),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:15),.false.,MomExt(1:4,9:11),PSWgt3)
   call Kinematics_TTBARPHOTON(0,MomExt(1:4,1:12),(/4,5,3,1,2,0,6,7,8,9,10,11/),applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_NLODK_ttbp = 0d0
      return
   endif

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt
   RunFactor = RunAlphaS(2,MuRen)

do npdf=1,2
    if(npdf.eq.1) then
        PDFFac(1:2) = PDFFac_a(1:2)
    elseif(npdf.eq.2) then
        if( Process.eq.20 ) cycle
        PDFFac(1:2) = PDFFac_b(1:2)
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
    endif
    ISFac = MomCrossing(MomExt)

   LO_Res_Unpol = (0d0,0d0)
   NLO_Res_UnPol= (0d0,0d0)
   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
      call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
          TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
      enddo

      call TopDecay(ExtParticle(1),DK_1L_T,MomExt(1:4,6:8))
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
          DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
      enddo

      NLO_Res_Pol = (0d0,0d0)
      if(PROCESS.EQ.20) then
          do jPrimAmp=1,NumBornAmps
          do iPrimAmp=1,NumBornAmps
              NLO_Res_Pol = NLO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(DKResult(jPrimAmp)) )*PDFFac(glu)
          enddo
          enddo
      elseif(PROCESS.EQ.22) then
         LOPartAmp(up) = TreeResult(1) + Q_up/Q_top*TreeResult(2)
         LOPartAmp(dn) = TreeResult(1) + Q_dn/Q_top*TreeResult(2)

         NLOPartAmp(up) = DKResult(1) + Q_up/Q_top*DKResult(2)
         NLOPartAmp(dn) = DKResult(1) + Q_dn/Q_top*DKResult(2)

         NLO_Res_Pol   = ColLO_ttbqqb(1,1) * ( dreal(LOPartAmp(up)*dconjg(NLOPartAmp(up)))*PDFFac(up)  +  dreal(LOPartAmp(dn)*dconjg(NLOPartAmp(dn)))*PDFFac(dn))
      endif
      NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
   enddo!helicity loop
!  normalization
   NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * Q_top**2*alpha4Pi*PhotonCouplCorr * PreFac
   EvalCS_NLODK_ttbp = EvalCS_NLODK_ttbp + dble(NLO_Res_Unpol)

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo

enddo! npdf loop
call swapMom(MomExt(1:4,1),MomExt(1:4,2))! swap back








!-------------------------------------
! one loop correction to top-decay   |
!-------------------------------------
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt
   RunFactor = RunAlphaS(2,MuRen)

do npdf=1,2
    if(npdf.eq.1) then
        PDFFac(1:2) = PDFFac_a(1:2)
    elseif(npdf.eq.2) then
        if( Process.eq.20 ) cycle
        PDFFac(1:2) = PDFFac_b(1:2)
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
    endif
    ISFac = MomCrossing(MomExt)

   LO_Res_Unpol = (0d0,0d0)
   NLO_Res_UnPol= (0d0,0d0)
   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
      call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
          TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
      enddo

      call TopDecay(ExtParticle(2),DK_1L_T,MomExt(1:4,9:11))
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
          DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
      enddo

      NLO_Res_Pol = (0d0,0d0)
      if(PROCESS.EQ.20) then
          do jPrimAmp=1,NumBornAmps
          do iPrimAmp=1,NumBornAmps
              NLO_Res_Pol = NLO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(DKResult(jPrimAmp)) )*PDFFac(glu)
          enddo
          enddo
      elseif(PROCESS.EQ.22) then
         LOPartAmp(up) = TreeResult(1) + Q_up/Q_top*TreeResult(2)
         LOPartAmp(dn) = TreeResult(1) + Q_dn/Q_top*TreeResult(2)

         NLOPartAmp(up) = DKResult(1) + Q_up/Q_top*DKResult(2)
         NLOPartAmp(dn) = DKResult(1) + Q_dn/Q_top*DKResult(2)

         NLO_Res_Pol   = ColLO_ttbqqb(1,1) * ( dreal(LOPartAmp(up)*dconjg(NLOPartAmp(up)))*PDFFac(up)  +  dreal(LOPartAmp(dn)*dconjg(NLOPartAmp(dn)))*PDFFac(dn))
      endif
      NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
   enddo!helicity loop
!  normalization
   NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * Q_top**2*alpha4Pi*PhotonCouplCorr * PreFac
   EvalCS_NLODK_ttbp = EvalCS_NLODK_ttbp + NLO_Res_Unpol


   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo

enddo! npdf loop
call swapMom(MomExt(1:4,1),MomExt(1:4,2))! swap back





ELSEIF( CORRECTION.EQ.5 ) THEN
if( DKRE_switch.eq.0 .or. DKRE_switch.eq.1 ) then
!----------------------------------------
! real gluon emission for Anti-top decay |
!----------------------------------------
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:14),.true.,MomExt(1:4,6:9),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(15:18),.false.,MomExt(1:4,10:12),PSWgt3)
   call CheckSing(MomExt(1:4,6:9),applySingCut)
   if( applySingCut) then
      EvalCS_NLODK_ttbp = 0d0
      goto 13
   endif
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt
   RunFactor = RunAlphaS(2,MuRen)

do npdf=1,2
    if(npdf.eq.1) then
        PDFFac(1:2) = PDFFac_a(1:2)
    elseif(npdf.eq.2) then
        if( Process.eq.20 ) cycle
        PDFFac(1:2) = PDFFac_b(1:2)
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
    endif
    ISFac = MomCrossing(MomExt)
    call Kinematics_TTBARPHOTON(1,MomExt(1:4,1:12),(/4,5,3,1,2,9,6,7,8,10,11,12/),applyPSCut,NBin)
    if( applyPSCut ) then
      goto 14
    endif


   LO_Res_Unpol = (0d0,0d0)
   do iHel=1,NumHelicities ! loop over initial state chiralities
   do GluHel=1,-1,-2 ! loop over additional gluon chiralities from decay
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_RE_T,MomExt(1:4,6:9),GluonHel=GluHel)
      call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,10:12))

      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo

      LO_Res_Pol = (0d0,0d0)
      if(PROCESS.EQ.20) then
          do jPrimAmp=1,NumBornAmps
          do iPrimAmp=1,NumBornAmps
              LO_Res_Pol = LO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)*PDFFac(glu)
          enddo
          enddo
      elseif(PROCESS.EQ.22) then
         LOPartAmp(up) = BornAmps(1)%Result + Q_up/Q_top * BornAmps(2)%Result
         LOPartAmp(dn) = BornAmps(1)%Result + Q_dn/Q_top * BornAmps(2)%Result
         LO_Res_Pol    = ColLO_ttbqqb(1,1) * ( LOPartAmp(up)*dconjg(LOPartAmp(up))*PDFFac(up) + LOPartAmp(dn)*dconjg(LOPartAmp(dn))*PDFFac(dn))
      endif
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
   enddo!helicity loop
   enddo!helicity loop

!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * Q_top**2*alpha4Pi*PhotonCouplCorr * PreFac
   EvalCS_NLODK_ttbp = EvalCS_NLODK_ttbp + dble(LO_Res_Unpol)

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol))
   enddo


14 continue
!-------------------------------------
! dipole subtraction for Atop-decay |
!-------------------------------------
  call WTransform(MomExt(1:4,6:9),MomExtTd(1:4,6:8),pbDpg,ptDpg,ptDpb)
  omz=ptDpg/(ptDpb+ptDpg-pbDpg)  !  for some reason this is not (1-z) as defined in the paper...
  rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
  z=1d0-omz
  y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2

  Dipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
  Dipole = Dipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )


   MomExtTd(1:4,1:5)   = MomExt(1:4,1:5)
   MomExtTd(1:4,10:12) = MomExt(1:4,10:12)
   call Kinematics_TTBARPHOTON(0,MomExtTd(1:4,1:12),(/4,5,3,1,2,0,6,7,8,10,11,12/),applyPSCut,NBin)
   if( applyPSCut ) then
      goto 13
   endif

   Dip_Res_Unpol= (0d0,0d0)
   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()

      call TopDecay(ExtParticle(1),DK_LO,MomExtTd(1:4,6:8))
      call TopDecay(ExtParticle(2),DK_LO,MomExtTd(1:4,10:12))

      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo

      LO_Res_Pol = (0d0,0d0)
      if(PROCESS.EQ.20) then
          do jPrimAmp=1,NumBornAmps
          do iPrimAmp=1,NumBornAmps
              LO_Res_Pol = LO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)*PDFFac(glu)
          enddo
          enddo
      elseif(PROCESS.EQ.22) then
         LOPartAmp(up) = BornAmps(1)%Result + Q_up/Q_top * BornAmps(2)%Result
         LOPartAmp(dn) = BornAmps(1)%Result + Q_dn/Q_top * BornAmps(2)%Result
         LO_Res_Pol    = ColLO_ttbqqb(1,1) * ( LOPartAmp(up)*dconjg(LOPartAmp(up))*PDFFac(up) + LOPartAmp(dn)*dconjg(LOPartAmp(dn))*PDFFac(dn))
      endif
      Dip_Res_UnPol = Dip_Res_UnPol + LO_Res_Pol
   enddo!helicity loop

!  normalization
   Dip_Res_Unpol = Dip_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * Dipole * Q_top**2*alpha4Pi*PhotonCouplCorr * PreFac
   EvalCS_NLODK_ttbp = EvalCS_NLODK_ttbp + Dip_Res_Unpol

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(Dip_Res_Unpol))
   enddo

enddo! npdf loop
call swapMom(MomExt(1:4,1),MomExt(1:4,2))! swap back

!             print *, MomDK(1,4)/EHat,EvalCS_NLODK_ttbp/dble(Dip_Res_Unpol)
!             print *, (MomDK(1:4,1).dot.MomDK(1:4,4))/EHat**2,EvalCS_NLODK_ttbp/dble(Dip_Res_Unpol)
!             EvalCS_NLODK_ttbp=EvalCS_NLODK_ttbp/Vgswgt
!             pause
!             return

endif! DKRE_switch

13 continue



if( DKRE_switch.eq.0 .or. DKRE_switch.eq.2 ) then
!-------------------------------------
! real gluon emission for top-decay  |
!-------------------------------------
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomExt(1:4,6:8),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:18),.true.,MomExt(1:4,9:12),PSWgt3)
   call CheckSing(MomExt(1:4,9:12),applySingCut)
   if( applySingCut ) then
      goto 17
   endif
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt
   RunFactor = RunAlphaS(2,MuRen)

do npdf=1,2
    if(npdf.eq.1) then
        PDFFac(1:2) = PDFFac_a(1:2)
    elseif(npdf.eq.2) then
        if( Process.eq.20 ) cycle
        PDFFac(1:2) = PDFFac_b(1:2)
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
    endif
    ISFac = MomCrossing(MomExt)
    call Kinematics_TTBARPHOTON(1,MomExt(1:4,1:12),(/4,5,3,1,2,12,6,7,8,9,10,11/),applyPSCut,NBin)
    if( applyPSCut ) then
      goto 15
    endif

   LO_Res_Unpol = (0d0,0d0)
   do iHel=1,NumHelicities ! loop over initial state chiralities
   do GluHel=1,-1,-2 ! loop over additional gluon chiralities from decay
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
      call TopDecay(ExtParticle(2),DK_RE_T,MomExt(1:4,9:12),GluonHel=GluHel)
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo

      LO_Res_Pol = (0d0,0d0)
      if(PROCESS.EQ.20) then
          do jPrimAmp=1,NumBornAmps
          do iPrimAmp=1,NumBornAmps
              LO_Res_Pol = LO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)*PDFFac(glu)
          enddo
          enddo
      elseif(PROCESS.EQ.22) then
         LOPartAmp(up) = BornAmps(1)%Result + Q_up/Q_top * BornAmps(2)%Result
         LOPartAmp(dn) = BornAmps(1)%Result + Q_dn/Q_top * BornAmps(2)%Result
         LO_Res_Pol    = ColLO_ttbqqb(1,1) * ( LOPartAmp(up)*dconjg(LOPartAmp(up))*PDFFac(up) + LOPartAmp(dn)*dconjg(LOPartAmp(dn))*PDFFac(dn))
      endif
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
   enddo!helicity loop
   enddo!helicity loop

!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * Q_top**2*alpha4Pi*PhotonCouplCorr * PreFac
   EvalCS_NLODK_ttbp = EvalCS_NLODK_ttbp + dble(LO_Res_Unpol)

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol))
   enddo


15 continue
!-------------------------------------
! dipole subtraction for top-decay   |
!-------------------------------------
  call WTransform(MomExt(1:4,9:12),MomExtTd(1:4,9:11),pbDpg,ptDpg,ptDpb)
  omz=ptDpg/(ptDpb+ptDpg-pbDpg)
  rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
  z=1d0-omz
  y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2
  Dipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
  Dipole = Dipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )

   MomExtTd(1:4,1:8) = MomExt(1:4,1:8)
   call Kinematics_TTBARPHOTON(0,MomExtTd(1:4,1:12),(/4,5,3,1,2,0,6,7,8,9,10,11/),applyPSCut,NBin)
   if( applyPSCut ) then
      goto 17
   endif

   Dip_Res_Unpol= (0d0,0d0)
   do iHel=1,NumHelicities ! loop over initial state chiralities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()

      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
      call TopDecay(ExtParticle(2),DK_LO,MomExtTd(1:4,9:11))

      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo

      LO_Res_Pol = (0d0,0d0)
      if(PROCESS.EQ.20) then
          do jPrimAmp=1,NumBornAmps
          do iPrimAmp=1,NumBornAmps
              LO_Res_Pol = LO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)*PDFFac(glu)
          enddo
          enddo
      elseif(PROCESS.EQ.22) then
         LOPartAmp(up) = BornAmps(1)%Result + Q_up/Q_top * BornAmps(2)%Result
         LOPartAmp(dn) = BornAmps(1)%Result + Q_dn/Q_top * BornAmps(2)%Result
         LO_Res_Pol    = ColLO_ttbqqb(1,1) * ( LOPartAmp(up)*dconjg(LOPartAmp(up))*PDFFac(up) + LOPartAmp(dn)*dconjg(LOPartAmp(dn))*PDFFac(dn))
      endif
      Dip_Res_UnPol = Dip_Res_UnPol + LO_Res_Pol
   enddo!helicity loop

!  normalization
   Dip_Res_Unpol = Dip_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * Dipole * Q_top**2*alpha4Pi*PhotonCouplCorr * PreFac
   EvalCS_NLODK_ttbp = EvalCS_NLODK_ttbp + dble(Dip_Res_Unpol)

!             print *, MomDK(1,7)/EHat,pbDpg/EHat**2,EvalCS_NLODK_ttbp/dble(Dip_Res_Unpol)
!             print *, dble(LO_Res_Unpol),dble(Dip_Res_Unpol)
!             pause

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(Dip_Res_Unpol))
   enddo


enddo! npdf loop
call swapMom(MomExt(1:4,1),MomExt(1:4,2))! swap back

endif! DKRE_switch
17 continue
ENDIF


   EvalCS_NLODK_ttbp = EvalCS_NLODK_ttbp/VgsWgt
return

END FUNCTION






FUNCTION EvalCS_NLODKP_ttb(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModTTBP_NLODK
implicit none
real(8) ::  EvalCS_NLODKP_ttb,yRnd(1:VegasMxDim),VgsWgt
complex(8) :: LO_Res_Pol,LO_Res_Unpol,Dip_Res_Unpol,NLO_Res_Pol,NLO_Res_UnPol
complex(8) :: TreeResult(1:NumBornAmps), DKResult(1:NumBornAmps)
integer :: iHel,PhoHel,iPrimAmp,jPrimAmp, GluHel
real(8) :: EHat,PSWgt1,PSWgt2,PSWgt3,ISFac
real(8) :: MomExt(1:4,1:12),MomExtTilde(1:4,1:12), MomExtBornTilde(1:4,1:4)
logical :: applyPSCut,applySingCut
integer, parameter :: ParityFlip=1
real(8) :: tau,eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,RunFactor
real(8) :: pdf(-6:6,1:2),ColCorrLO(1:NumBornAmps,1:NumBornAmps)
integer :: NBin(1:NumMaxHisto),NHisto,nPhoRad
real(8), parameter :: CF=4d0/3d0,PhotonCouplCorr=2d0
include "vegas_common.f"
! SCHARFS STUFF
integer :: Top_Atop,i, count_t1, count_t2
real(8) :: dipole, galle, soft, coll
complex(8) :: test(1:4)
logical :: nan_t
!Dipole Stuff for contribution II.b
real(8) :: pbDpg,ptDpg,ptDpb, omz, rsq,z,y


EvalCS_NLODKP_ttb = 0d0

  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top) then
      EvalCS_NLODKP_ttb = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to2(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt1)
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
   ISFac = MomCrossing(MomExt)
   RunFactor = RunAlphaS(NLOParam,MuRen)

   call setPDFs(eta1,eta2,MuFac,pdf)
IF( PROCESS.EQ.21 ) THEN
   ColCorrLO(1:NumBornAmps,1:NumBornAmps) = ColLO_ttbgg(1:NumBornAmps,1:NumBornAmps)
   PDFFac = pdf(0,1) * pdf(0,2)
ELSEIF( PROCESS.EQ.23 ) THEN
   ColCorrLO(1:NumBornAmps,1:NumBornAmps) = ColLO_ttbqqb(1:NumBornAmps,1:NumBornAmps)
   PDFFac =       ( pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
                  + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
                  + pdf(Bot_,1)*pdf(ABot_,2) ) +                         &
                  ( pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
                  + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
                  + pdf(Bot_,2)*pdf(ABot_,1) )
ENDIF



IF( CORRECTION.EQ.4 ) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!VIRTUAL CORRECTIONS  CORRECTION = 4
!----------------------------------
! photon emission off anti-top  with one loop correction   |
!----------------------------------
!goto 4444  !!!!!!!!!SWITCH IT OFF
!goto 6633
!goto 3366
      call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:11),.true.,MomExt(1:4,5:8),PSWgt2)
      call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(12:15),.false.,MomExt(1:4,9:11),PSWgt3)
      PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt * PDFFac

      NLO_Res_Unpol = (0d0,0d0)
      call Kinematics_TTBARPHOTON(0,MomExt(1:4,1:12),(/3,4,8,1,2,0,5,6,7,9,10,11/),applyPSCut,NBin)
      if( applyPSCut ) then
         goto 991
      endif


      do iHel=1,NumHelicities ! loop over initial state chiralities
         do PhoHel=1,-1,-2 ! loop over additional photon polarization
            call HelCrossing(Helicities(iHel,1:NumExtParticles))
            call SetPolarizations()
            !*******************
            ! 1st CALL FOR INTERFERENCE TERMS
            ! Photon from Anti-Top
            ! NO PHOTON FROM TOP
            call TopDecay(ExtParticle(1),DKP_LO_T,MomExt(1:4,5:8),PhotonHel=PhoHel)
            call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))
!            write(*,*) "Markus, ExtParticle(1)", ExtParticle(1)%Pol
!            call TopDecay(ExtParticle(1),DKP_1L_T,MomExt(1:4,5:8),PhotonHel=PhoHel)
!            write(*,*) "me, ExtParticle(1)", ExtParticle(1)%Pol
!            pause
            ! Production Process calculated with ExtParticle(1) and ExtParticle(2)
            do iPrimAmp=1,NumBornAmps
               call EvalTree(BornAmps(iPrimAmp))
               TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
            enddo
            !*******************

            !*******************
            ! 2nd CALL FOR INTERFERENCE TERMS
            ! NLO CALL FOR ANTI-TOP + Photon
            call TopDecay(ExtParticle(1),DKP_1L_T,MomExt(1:4,5:8),PhotonHel=PhoHel)
            ! LO CALL FOR TOP, NO PHOTON
!            write(*,*) "Spino-Atop", ExtParticle(1)%Pol(1:4)
!            write(*,*) "p1t_out", MomExt(1:4,3)
            call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))
            do iPrimAmp=1,NumBornAmps
               call EvalTree(BornAmps(iPrimAmp))
               DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
            enddo
            !*******************
            NLO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
               do iPrimAmp=1,NumBornAmps
                  NLO_Res_Pol = NLO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*dreal(TreeResult(jPrimAmp)* &
                  dconjg(DKResult(iPrimAmp)))
!                  NLO_Res_Pol = NLO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*TreeResult(iPrimAmp)* &
!                  dconjg(TreeResult(jPrimAmp))

!                  write(*,*) "NLO_RES_pol",NLO_Res_Pol
!                  write(*,*) "ColCorrLO(iPrimAmp,jPrimAmp)",  ColCorrLO(iPrimAmp,jPrimAmp), ParityFlip
!                  write(*,*) "i,TreeResult ",iPrimAmp, TreeResult(iPrimAmp)
!                  write(*,*) "j,DKResult", jPrimAmp, DKResult(jPrimAmp)
               enddo
            enddo
            NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
         enddo!helicity loop
      enddo!helicity loop
!      write(*,*) "WQ-ATop", NLO_RES_UNPOL
!      write(*,*) "*****************"
!      pause

!      write(*,*) "WQ-Anti-Top", NLO_RES_UNPOL
!      pause
!  normalization
   NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac


   EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb + dble(NLO_Res_Unpol)

!   galle =1.d0
!   NBIn(1:NumHistoGrams) = 1
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dreal(NLO_Res_UnPol))
!      call intoHisto(NHisto,NBin(NHisto),galle)
   enddo

!   nan_t = IsNan(EvalCS_NLODKP_ttb  )
!   write(*,*) "Hallo 1", EvalCS_NLODKP_ttb
!   if(nan_t) then
!      write(*,*) "1"
!   endif


991 continue
!goto 6789
!---------------------------------
! photon emission off top-decay    with one loop correction   |
!---------------------------------
   call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:8),.false.,MomExt(1:4,5:7),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(9:15),.true.,MomExt(1:4,8:11),PSWgt3)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt * PDFFac

   NLO_Res_Unpol = (0d0,0d0)
   call Kinematics_TTBARPHOTON(0,MomExt(1:4,1:12),(/3,4,11,1,2,0,5,6,7,8,9,10/),applyPSCut,NBin)
   if( applyPSCut ) then
!      EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb/VgsWgt
!      return
      goto 6633
   endif


   do iHel=1,NumHelicities ! loop over initial state chiralities
   do PhoHel=1,-1,-2 ! loop over additional photon chiralities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      !*******************
      ! 1st CALL FOR INTERFERENCE TERMS
      ! NO PHOTON FROM ANTI-TOP
      ! Photon from Top
      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,5:7))
      call TopDecay(ExtParticle(2),DKP_LO_T,MomExt(1:4,8:11),PhotonHel=PhoHel)
      ! Production Process calculated with ExtParticle(1) and ExtParticle(2)

!            TEST(1:4) = ExtParticle(2)%Pol(1:4)
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
          TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
      enddo
      !*******************
      !*******************
      ! 2nd CALL FOR INTERFERENCE TERMS
      ! LO CALL FOR ANTI-TOP, NO PHOTON
      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,5:7))
      ! NLO CALL FOR TOP + Photon
      call TopDecay(ExtParticle(2),DKP_1L_T,MomExt(1:4,8:11),PhotonHel=PhoHel)
      do iPrimAmp=1,NumBornAmps
         call EvalTree(BornAmps(iPrimAmp))
         DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
      enddo
      !*******************
!      write(*,*) "Markus-Spi", test
!      write(*,*) "Me-Spi/Markus-Spi", Extparticle(2)%Pol(1:4)/ test

      NLO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
         do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*dreal(TreeResult(iPrimAmp)* &
            dconjg(DKResult(jPrimAmp)))
         enddo
      enddo
      NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol

   enddo!helicity loop
   enddo!helicity loop

!  normalization
   NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac

   EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb + dble(NLO_Res_Unpol)
   galle =1.d0
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
!      call intoHisto(NHisto,NBin(NHisto),galle)
   enddo
!   nan_t = IsNan(EvalCS_NLODKP_ttb  )
!      write(*,*) "Halo 2"
!   if(nan_t) then
!      write(*,*) "2"
!   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CON IIa2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CON IIa2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CON IIa2

6633 continue
!----------------------------------
! photon emission off W-  with one loop correction  at anti-top |
!----------------------------------
  call EvalPhasespace_TopDecay2(MomExt(1:4,3),yRnd(5:11),.true.,MomExt(1:4,5:8),PSWgt2)
  call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(12:15),.false.,MomExt(1:4,9:11),PSWgt3)
  PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2 *PSWgt3* VgsWgt * PDFFac

   NLO_Res_Unpol = (0.d0,0.d0)
     call Kinematics_TTBARPHOTON(0,MomExt(1:4,1:12),(/3,4,8,1,2,0,5,6,7,9,10,11/),applyPSCut,NBin)
     if( applyPSCut ) then
        goto 3366
     endif

      do iHel=1,NumHelicities ! loop over initial state chiralities
         do PhoHel=1,-1,-2 ! loop over additional photon polarization
            call HelCrossing(Helicities(iHel,1:NumExtParticles))
            call SetPolarizations()
            !*******************
            ! 1st CALL FOR INTERFERENCE TERMS
            ! Photon from Anti-Top/W-
            call TopDecay(ExtParticle(1),DKP_LO_L,MomExt(1:4,5:8),PhotonHel=PhoHel)
            ! NO PHOTON FROM TOP
            call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))
            ! Production Process calculated with ExtParticle(1) and ExtParticle(2)
            do iPrimAmp=1,NumBornAmps
               call EvalTree(BornAmps(iPrimAmp))
               TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
            enddo
            !*******************

            !*******************
            ! 2nd CALL FOR INTERFERENCE TERMS
            ! Photon from W-
               call TopDecay(ExtParticle(1),DKP_1L_L,MomExt(1:4,5:8),PhotonHel=PhoHel)
            ! LO CALL FOR TOP, NO PHOTON
            call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))
            do iPrimAmp=1,NumBornAmps
               call EvalTree(BornAmps(iPrimAmp))
               DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
            enddo
            !*******************
            NLO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
               do iPrimAmp=1,NumBornAmps
                  NLO_Res_Pol = NLO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*dreal(TreeResult(jPrimAmp)* &
                  dconjg(DKResult(iPrimAmp)))
!                  NLO_Res_Pol = NLO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*TreeResult(iPrimAmp)* &
!                  dconjg(TreeResult(jPrimAmp))

!                  write(*,*) "NLO_RES_pol",NLO_Res_Pol
!                  write(*,*) "ColCorrLO(iPrimAmp,jPrimAmp)",  ColCorrLO(iPrimAmp,jPrimAmp), ParityFlip
!                  write(*,*) "i,TreeResult ",iPrimAmp, TreeResult(iPrimAmp)
!                  write(*,*) "j,DKResult", jPrimAmp, DKResult(jPrimAmp)
               enddo
            enddo
            NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
         enddo!helicity loop
      enddo!helicity loop
!      write(*,*) "WQ-ATop", NLO_RES_UNPOL
!      write(*,*) "*****************"
!      pause

!  normalization
   NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
   !   nan_t = IsNan(NLO_Res_Unpol )


   EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb + dble(NLO_Res_Unpol)

!   galle =1.d0
!   NBIn(1:NumHistoGrams) = 1
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dreal(NLO_Res_UnPol))
!      call intoHisto(NHisto,NBin(NHisto),galle)
   enddo
!----------------------------------
! photon emission off W+  with one loop correction  at top |
!----------------------------------
3366 continue
!goto 6789

    call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:8),.false.,MomExt(1:4,5:7),PSWgt2)
    call EvalPhasespace_TopDecay2(MomExt(1:4,4),yRnd(9:15),.true.,MomExt(1:4,8:11),PSWgt3)
    PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt * PDFFac

    NLO_Res_Unpol = (0.d0,0.d0)
    call Kinematics_TTBARPHOTON(0,MomExt(1:4,1:12),(/3,4,11,1,2,0,5,6,7,8,9,10/),applyPSCut,NBin)
    if( applyPSCut ) then
       goto 4444
    endif


   do iHel=1,NumHelicities ! loop over initial state chiralities
   do PhoHel=1,-1,-2 ! loop over additional photon chiralities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      !*******************
      ! 1st CALL FOR INTERFERENCE TERMS
      ! Photon from Anti-Top
      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,5:7))
      call TopDecay(ExtParticle(2),DKP_LO_L,MomExt(1:4,8:11),PhotonHel=PhoHel)
      ! Production Process calculated with ExtParticle(1) and ExtParticle(2)
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
          TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
      enddo
      !*******************
      !*******************
      ! 2nd CALL FOR INTERFERENCE TERMS
      ! NLO CALL FOR TOP NO Photon
      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,5:7))
      call TopDecay(ExtParticle(2),DKP_1L_L,MomExt(1:4,8:11),PhotonHel=PhoHel)
      do iPrimAmp=1,NumBornAmps
         call EvalTree(BornAmps(iPrimAmp))
         DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
      enddo
      !*******************
!      write(*,*) "Markus-Spi", test
!      write(*,*) "Me-Spi/Markus-Spi", Extparticle(2)%Pol(1:4)/ test

      NLO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
         do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*dreal(TreeResult(iPrimAmp)* &
            dconjg(DKResult(jPrimAmp)))
         enddo
      enddo
      NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol

   enddo!helicity loop
   enddo!helicity loop

!  normalization
   NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac

   EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb + dble(NLO_Res_Unpol)

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
!      call intoHisto(NHisto,NBin(NHisto),galle)
   enddo





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CON IIb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CON IIb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CON IIb

4444 continue
!goto 6789
!goto 4441
!----------------------------------
! photon emission off anti-top  with one loop correction  at top |
!----------------------------------
! tb -> bb W- gamma at LO including gamma-radiation from W; t -> W b at NLO

 do nPhoRad=nPhoRad1,nPhoRad2!   nPhoRad=1: photon radiation off top/bot/W,
 !                 nPhoRad=2: photon radiation off W/lep
    if( nPhoRad.eq.1 ) then
       call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:11),.true.,MomExt(1:4,5:8),PSWgt2)
    else
       call EvalPhasespace_TopDecay2(MomExt(1:4,3),yRnd(5:11),.true.,MomExt(1:4,5:8),PSWgt2)
    endif
    call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(12:15),.false.,MomExt(1:4,9:11),PSWgt3)
    PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt * PDFFac


    NLO_Res_Unpol = (0.d0,0.d0)
     call Kinematics_TTBARPHOTON(0,MomExt(1:4,1:12),(/3,4,8,1,2,0,5,6,7,9,10,11/),applyPSCut,NBin)
     if( applyPSCut ) then
        cycle
     endif


      do iHel=1,NumHelicities ! loop over initial state chiralities
         do PhoHel=1,-1,-2 ! loop over additional photon polarization
            call HelCrossing(Helicities(iHel,1:NumExtParticles))
            call SetPolarizations()
            !*******************
            ! 1st CALL FOR INTERFERENCE TERMS
            ! Photon from Anti-Top/W-
            if( nPhoRad.eq.1 ) then
               call TopDecay(ExtParticle(1),DKP_LO_T,MomExt(1:4,5:8),PhotonHel=PhoHel)
            else
               call TopDecay(ExtParticle(1),DKP_LO_L,MomExt(1:4,5:8),PhotonHel=PhoHel)
            endif
            ! NO PHOTON FROM TOP
            call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))
            ! Production Process calculated with ExtParticle(1) and ExtParticle(2)
            do iPrimAmp=1,NumBornAmps
               call EvalTree(BornAmps(iPrimAmp))
               TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
            enddo
            !*******************

            !*******************
            ! 2nd CALL FOR INTERFERENCE TERMS
            ! Photon from Anti-Top/W-
            if( nPhoRad.eq.1 ) then
               call TopDecay(ExtParticle(1),DKP_LO_T,MomExt(1:4,5:8),PhotonHel=PhoHel)
            else
               call TopDecay(ExtParticle(1),DKP_LO_L,MomExt(1:4,5:8),PhotonHel=PhoHel)
            endif
            ! NLO CALL FOR TOP, NO PHOTON
            call TopDecay(ExtParticle(2),DK_1L_T,MomExt(1:4,9:11))
            do iPrimAmp=1,NumBornAmps
               call EvalTree(BornAmps(iPrimAmp))
               DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
            enddo
            !*******************
            NLO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
               do iPrimAmp=1,NumBornAmps
                  NLO_Res_Pol = NLO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*dreal(TreeResult(jPrimAmp)* &
                  dconjg(DKResult(iPrimAmp)))
!                  NLO_Res_Pol = NLO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*TreeResult(iPrimAmp)* &
!                  dconjg(TreeResult(jPrimAmp))

!                  write(*,*) "NLO_RES_pol",NLO_Res_Pol
!                  write(*,*) "ColCorrLO(iPrimAmp,jPrimAmp)",  ColCorrLO(iPrimAmp,jPrimAmp), ParityFlip
!                  write(*,*) "i,TreeResult ",iPrimAmp, TreeResult(iPrimAmp)
!                  write(*,*) "j,DKResult", jPrimAmp, DKResult(jPrimAmp)
               enddo
            enddo
            NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
         enddo!helicity loop
      enddo!helicity loop
!      write(*,*) "WQ-ATop", NLO_RES_UNPOL
!      write(*,*) "*****************"
!      pause

!  normalization
   NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
   !   nan_t = IsNan(NLO_Res_Unpol )


   EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb + dble(NLO_Res_Unpol)

!   galle =1.d0
!   NBIn(1:NumHistoGrams) = 1
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dreal(NLO_Res_UnPol))
!      call intoHisto(NHisto,NBin(NHisto),galle)
   enddo
enddo ! nPhorad

!goto 6789
!4441 continue
!---------------------------------
! photon emission off top-decay  |
!---------------------------------
! t -> b W- gamma at LO including gamma-radiation from W; tb -> W bb at NLO
 do nPhoRad=nPhoRad1,nPhoRad2
    call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:8),.false.,MomExt(1:4,5:7),PSWgt2)
    if( nPhoRad.eq.1 ) then
       call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(9:15),.true.,MomExt(1:4,8:11),PSWgt3)
    else
       call EvalPhasespace_TopDecay2(MomExt(1:4,4),yRnd(9:15),.true.,MomExt(1:4,8:11),PSWgt3)
   endif
    PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt * PDFFac

    NLO_Res_Unpol = (0.d0,0.d0)
    call Kinematics_TTBARPHOTON(0,MomExt(1:4,1:12),(/3,4,11,1,2,0,5,6,7,8,9,10/),applyPSCut,NBin)
    if( applyPSCut ) then
       cycle
    endif


   do iHel=1,NumHelicities ! loop over initial state chiralities
   do PhoHel=1,-1,-2 ! loop over additional photon chiralities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      !*******************
      ! 1st CALL FOR INTERFERENCE TERMS
      ! Photon from Anti-Top
      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,5:7))
      if( nPhoRad .eq. 1) then
         call TopDecay(ExtParticle(2),DKP_LO_T,MomExt(1:4,8:11),PhotonHel=PhoHel)
      else
         call TopDecay(ExtParticle(2),DKP_LO_L,MomExt(1:4,8:11),PhotonHel=PhoHel)
      endif
      ! Production Process calculated with ExtParticle(1) and ExtParticle(2)
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
          TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
      enddo
      !*******************
      !*******************
      ! 2nd CALL FOR INTERFERENCE TERMS
      ! NLO CALL FOR TOP NO Photon
      call TopDecay(ExtParticle(1),DK_1L_T,MomExt(1:4,5:7))
      if( nPhoRad.eq.1 ) then
         call TopDecay(ExtParticle(2),DKP_LO_T,MomExt(1:4,8:11),PhotonHel=PhoHel)
      else
         call TopDecay(ExtParticle(2),DKP_LO_L,MomExt(1:4,8:11),PhotonHel=PhoHel)
      endif
      do iPrimAmp=1,NumBornAmps
         call EvalTree(BornAmps(iPrimAmp))
         DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
      enddo
      !*******************
!      write(*,*) "Markus-Spi", test
!      write(*,*) "Me-Spi/Markus-Spi", Extparticle(2)%Pol(1:4)/ test

      NLO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,NumBornAmps
         do iPrimAmp=1,NumBornAmps
            NLO_Res_Pol = NLO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*dreal(TreeResult(iPrimAmp)* &
            dconjg(DKResult(jPrimAmp)))
         enddo
      enddo
      NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol

   enddo!helicity loop
   enddo!helicity loop

!  normalization
   NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac

   EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb + dble(NLO_Res_Unpol)
!   galle =1.d0

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
!      call intoHisto(NHisto,NBin(NHisto),galle)
   enddo

enddo !nPhoRad





ELSEIF (CORRECTION .eq. 5) then
!----------------------------------
! photon emission off anti-top  with real correction  |
!----------------------------------
!goto 3333
!goto 1999
!goto 1881
!goto 1818

      call EvalPhasespace_TopDecay3(MomExt(1:4,3),yRnd(5:14),MomExt(1:4,5:9),PSWgt2)
      if(PSWgt2.eq.0.d0) then
         goto 1999
      endif
      call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(15:18),.false.,MomExt(1:4,10:12),PSWgt3)
      PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt * PDFFac
      if(PSWgt3.eq.0.d0) then
         goto 1999
      endif
!goto 191
!      count_t1=0
      call Kinematics_TTBARPHOTON(1,MomExt(1:4,1:12),(/3,4,9,1,2,8,5,6,7,10,11,12/),applyPSCut,NBin)
      NLO_Res_Unpol = (0.d0,0.d0)
      if( applyPSCut ) then
         goto 191
      endif

      count_t1 =1
      do iHel=1,NumHelicities ! loop over initial state chiralities
         do PhoHel=1,-1,-2 ! loop over additional photon chiralities
            do GluHel=1,-1,-2 ! loop over additional real gluon chiralities
               call HelCrossing(Helicities(iHel,1:NumExtParticles))
               call SetPolarizations()
               !*******************
               ! 1st CALL FOR INTERFERENCE TERMS
               ! Photon + Gluon  from Anti-Top
               ! NO PHOTON FROM TOP
               call TopDecay(ExtParticle(1),DKP_RE_T,MomExt(1:4,5:9),GluonHel=GluHel,PhotonHel=PhoHel)
               call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,10:12))
               ! Production Process calculated with ExtParticle(1) and ExtParticle(2)
               do iPrimAmp=1,NumBornAmps
                  call EvalTree(BornAmps(iPrimAmp))
                  DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
               enddo
               NLO_Res_Pol = (0.d0,0.d0)
               do jPrimAmp=1,NumBornAmps
                  do iPrimAmp=1,NumBornAmps
                     NLO_Res_Pol = NLO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*DKResult(iPrimAmp)* &
                          dconjg(DKResult(jPrimAmp))
                  enddo
               enddo
               NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
            enddo
         enddo
      enddo
      soft =  MomExt(1,8)/m_top
      coll = dsqrt(MomExt(1:4,8).dot.MomExt(1:4,5)/m_top**2)

      NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac

      EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb + dble(NLO_Res_Unpol)


      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
      enddo

! print *, (MomExt(1:4,8).dot.MomExt(1:4,5))/m_top**2
! print *, "real",dble(NLO_Res_Unpol)
! print *, "mom",MomExt(1:4,9).dot.(MomExt(1:4,5)+MomExt(1:4,8))

191 continue
!----------------------------------
! photon emission off anti-top dipole -contribution
!----------------------------------
      Dip_Res_Unpol = (0.d0,0.d0)
      Top_Atop =-1
      MomExtTilde = MomExt
      call EvalTTBP_DIPOLDK(Top_Atop,MomExt(1:4,5:9),MomExtTilde(1:4,5:9),dipole)
      call Kinematics_TTBARPHOTON(0,MomExtTilde(1:4,1:12),(/3,4,9,1,2,0,5,6,7,10,11,12/),applyPSCut,NBin)
      if (applyPScut) then
         goto 1999
      endif


!      write(*,*) "MomExtTilde-gluon",MomExtTilde(1:4,8)
      soft =  MomExt(1,8)/m_top
      coll = dsqrt(MomExt(1:4,8).dot.MomExt(1:4,5)/m_top**2)


      do iHel=1,NumHelicities ! loop over initial state chiralities
         do PhoHel=1,-1,-2 ! loop over additional photon chiralities
            call HelCrossing(Helicities(iHel,1:NumExtParticles))
            call SetPolarizations()
            !*******************
            ! HERE ONLY BORN-INTEFERENCE IS NEEDED
            ! SET MOMENTA FIRST SINCE MomExtTilde(1:4,8) is the gluon momentum and not-existing in reduced kinematics
            MomExtBornTilde(1:4,1) = MomExtTilde(1:4,5)
            MomExtBornTilde(1:4,2) = MomExtTilde(1:4,6)
            MomExtBornTilde(1:4,3) = MomExtTilde(1:4,7)
            MomExtBornTilde(1:4,4) = MomExtTilde(1:4,9)
            call TopDecay(ExtParticle(1),DKP_LO_T,MomExtBornTilde(1:4,1:4),PhotonHel=PhoHel)
            call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,10:12))
            do iPrimAmp=1,NumBornAmps
               call EvalTree(BornAmps(iPrimAmp))
               DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
            enddo
            LO_Res_Pol = (0.d0,0.d0)

            do jPrimAmp=1,NumBornAmps
               do iPrimAmp=1,NumBornAmps
                  LO_Res_Pol = LO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*DKResult(iPrimAmp)* &
                       dconjg(DKResult(jPrimAmp))
               enddo
            enddo
            Dip_Res_UnPol = Dip_Res_UnPol + LO_Res_Pol
         enddo
      enddo

      Dip_Res_UnPol = Dip_Res_UnPol* ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac * (-dipole)

      EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb + dble(Dip_Res_Unpol)
      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),dble(Dip_Res_Unpol))
      enddo

! if(dble(Dip_Res_Unpol).ne.0d0.and.abs(-dble(Dip_Res_Unpol)/dble(NLO_Res_Unpol)-1d0).gt.1d-2) then
!     print *, "dipol",dble(Dip_Res_Unpol)
!     print *, "real",dble(NLO_Res_Unpol)
!     print *, "ratio",-dble(Dip_Res_Unpol)/dble(NLO_Res_Unpol)-1d0
!     ! print *, "mom",MomExtBornTilde(1:4,1).dot.MomExtBornTilde(1:4,4)
! !     pause
! endif
!      soft =  MomExt(1,8)/m_top
!      coll = dsqrt(MomExt(1:4,8).dot.MomExt(1:4,5)/m_top**2)
!      if((soft .le. 1.d-3 .or. coll .lt. 1.d-3) .and. dble(Dip_Res_Unpol) .ne. 0.d0 .and.  dble(NLO_Res_UnPol) .eq. 0.d0) then
!         write(66,*) "Anti-Top soft/coll event dipole"
!         write(66,*) "softness", soft
!         write(66,*) "kollinear", coll
!         write(66,*) "Dipol", Dip_Res_UnPol
!         write(66,*) "Real", NLO_Res_UnPol
!         write(66,*) "Quotient", NLO_Res_Unpol/Dip_Res_Unpol
!         write(66,*) "yRnd(1) = ", yRnd(1)
!         write(66,*) "yRnd(2) = ", yRnd(2)
!         write(66,*) "yRnd(3) = ", yRnd(3)
!         write(66,*) "yRnd(4) = ", yRnd(4)
!         write(66,*) "yRnd(5) = ", yRnd(5)
!         write(66,*) "yRnd(6) = ", yRnd(6)
!         write(66,*) "yRnd(7) = ", yRnd(7)
!         write(66,*) "yRnd(8) = ", yRnd(8)
!         write(66,*) "yRnd(9) = ", yRnd(9)
!         write(66,*) "yRnd(10) = ", yRnd(10)
!         write(66,*) "yRnd(11) = ", yRnd(11)
!         write(66,*) "yRnd(12) = ", yRnd(12)
!         write(66,*) "yRnd(13) = ", yRnd(13)
!         write(66,*) "yRnd(14) = ", yRnd(14)
!         write(66,*) "yRnd(15) = ", yRnd(15)
!         write(66,*) "yRnd(16) = ", yRnd(16)
!         write(66,*) "yRnd(17) = ", yRnd(17)
!         write(66,*) "yRnd(18) = ", yRnd(18)
!         do i=1,12
!            write(66,*) "MomExt",i,MomExt(1:4,i)
!            write(66,*) "MomExtTilde",i,MomExtTilde(1:4,i)
!         enddo
!         stop
!      endif
!      if((soft .le. 1.d-3 .or. coll .lt. 1.d-3) .and. dble(Dip_Res_Unpol) .eq. 0.d0 .and.  dble(NLO_Res_UnPol) .ne. 0.d0) then
!         write(77,*) "Anti-Top soft/coll event dipole"
!         write(77,*) "softness", soft
!         write(77,*) "kollinear", coll
!         write(77,*) "Dipol", Dip_Res_UnPol
!         write(77,*) "Real", NLO_Res_UnPol
!         write(77,*) "Quotient", NLO_Res_Unpol/Dip_Res_Unpol
!      endif


1999 continue





!goto 6789
!----------------------------------
! photon emission off top  with real correction  |
!----------------------------------
      call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:8),.false.,MomExt(1:4,5:7),PSWgt2)
      call EvalPhasespace_TopDecay3(MomExt(1:4,4),yRnd(9:18),MomExt(1:4,8:12),PSWgt3)
      if(PSWgt3.eq.0.d0) then
         goto 1818
      endif

      PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt * PDFFac
      call Kinematics_TTBARPHOTON(1,MomExt(1:4,1:12),(/3,4,12,1,2,11,5,6,7,8,9,10/),applyPSCut,NBin)

      NLO_Res_Unpol = (0.d0,0.d0)
      if( applyPSCut ) then
         goto 1199
      endif


      do iHel=1,NumHelicities ! loop over initial state chiralities
         do PhoHel=1,-1,-2 ! loop over additional photon chiralities
            do GluHel=1,-1,-2 ! loop over additional real gluon chiralities
               call HelCrossing(Helicities(iHel,1:NumExtParticles))
               call SetPolarizations()
               !*******************
               ! 1st CALL FOR INTERFERENCE TERMS
               ! Photon + Gluon  from Anti-Top
               ! NO PHOTON FROM TOP
               call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,5:7))
               call TopDecay(ExtParticle(2),DKP_RE_T,MomExt(1:4,8:12),GluonHel=GluHel,PhotonHel=PhoHel)
               ! Production Process calculated with ExtParticle(1) and ExtParticle(2)
               do iPrimAmp=1,NumBornAmps
                  call EvalTree(BornAmps(iPrimAmp))
                  DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
               enddo
               NLO_Res_Pol = (0.d0,0.d0)
               do jPrimAmp=1,NumBornAmps
                  do iPrimAmp=1,NumBornAmps
                     NLO_Res_Pol = NLO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*DKResult(iPrimAmp)* &
                          dconjg(DKResult(jPrimAmp))
                  enddo
               enddo
               NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
            enddo
         enddo
      enddo
      NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac

      EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb + dble(NLO_Res_Unpol)

      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
      enddo
!      write(*,*) "RealBin", NBin




1199 continue
!----------------------------------
! photon emission off top dipole-contribution
!----------------------------------
      MomExtTilde = MomExt
      Dip_Res_Unpol = (0.d0,0.d0)
      Top_Atop =1
      call EvalTTBP_DIPOLDK(Top_Atop,MomExt(1:4,8:12),MomExtTilde(1:4,8:12),dipole)
      call Kinematics_TTBARPHOTON(0,MomExtTilde(1:4,1:12),(/3,4,12,1,2,0,5,6,7,8,9,10/),applyPSCut,NBin)
      if (applyPScut) then
!         EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb/VgsWgt
!         return
         goto 1818
      endif
      do iHel=1,NumHelicities ! loop over initial state chiralities
         do PhoHel=1,-1,-2 ! loop over additional photon chiralities
            call HelCrossing(Helicities(iHel,1:NumExtParticles))
            call SetPolarizations()
            !*******************
            ! HERE ONLY BORN-INTEFERENCE IS NEEDED
            ! SET MOMENTA FIRST SINCE MomExtTilde(1:4,11) is the gluon momentum and no-existing in reduced kinematics
            MomExtBornTilde(1:4,1) = MomExtTilde(1:4,8)
            MomExtBornTilde(1:4,2) = MomExtTilde(1:4,9)
            MomExtBornTilde(1:4,3) = MomExtTilde(1:4,10)
            MomExtBornTilde(1:4,4) = MomExtTilde(1:4,12)
            call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,5:7))
            call TopDecay(ExtParticle(2),DKP_LO_T,MomExtBornTilde(1:4,1:4),PhotonHel=PhoHel)
            do iPrimAmp=1,NumBornAmps
               call EvalTree(BornAmps(iPrimAmp))
               DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
            enddo
            LO_Res_Pol = (0.d0,0.d0)

            do jPrimAmp=1,NumBornAmps
               do iPrimAmp=1,NumBornAmps
                  LO_Res_Pol = LO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*DKResult(iPrimAmp)* &
                       dconjg(DKResult(jPrimAmp))
               enddo
            enddo
            Dip_Res_UnPol = Dip_Res_UnPol + LO_Res_Pol
         enddo
      enddo

      Dip_Res_UnPol = Dip_Res_UnPol* ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac * (-dipole)

      EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb + dble(Dip_Res_Unpol)


      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),dble(Dip_Res_Unpol))
      enddo

!      soft =  MomExt(1,11)/m_top
!      coll = dsqrt(MomExt(1:4,11).dot.MomExt(1:4,8)/m_top**2)
!      if((soft .le. 1.d-3 .or. coll .lt. 1.d-3) .and. dble(Dip_Res_Unpol) .ne. 0.d0 .and.  dble(NLO_Res_UnPol) .eq. 0.d0) then
!         write(66,*) "Top soft/coll event dipole"
!         write(66,*) "softness", soft
!         write(66,*) "kollinear", coll
!         write(66,*) "Dipol", Dip_Res_UnPol
!         write(66,*) "Real", NLO_Res_UnPol
!         write(66,*) "Quotient", NLO_Res_Unpol/Dip_Res_Unpol
!         write(66,*) "yRnd(1) = ", yRnd(1)
!         write(66,*) "yRnd(2) = ", yRnd(2)
!         write(66,*) "yRnd(3) = ", yRnd(3)
!         write(66,*) "yRnd(4) = ", yRnd(4)
!         write(66,*) "yRnd(5) = ", yRnd(5)
!         write(66,*) "yRnd(6) = ", yRnd(6)
!         write(66,*) "yRnd(7) = ", yRnd(7)
!         write(66,*) "yRnd(8) = ", yRnd(8)
!         write(66,*) "yRnd(9) = ", yRnd(9)
!         write(66,*) "yRnd(10) = ", yRnd(10)
!         write(66,*) "yRnd(11) = ", yRnd(11)
!         write(66,*) "yRnd(12) = ", yRnd(12)
!         write(66,*) "yRnd(13) = ", yRnd(13)
!         write(66,*) "yRnd(14) = ", yRnd(14)
!         write(66,*) "yRnd(15) = ", yRnd(15)
!         write(66,*) "yRnd(16) = ", yRnd(16)
!         write(66,*) "yRnd(17) = ", yRnd(17)
!         write(66,*) "yRnd(18) = ", yRnd(18)
!         do i=1,12
!            write(66,*) "MomExt",i,MomExt(1:4,i)
!            write(66,*) "MomExtTilde",i,MomExtTilde(1:4,i)
!         enddo
!         stop
!      endif
!      if((soft .le. 1.d-3 .or. coll .lt. 1.d-3) .and. dble(Dip_Res_Unpol) .eq. 0.d0 .and.  dble(NLO_Res_UnPol) .ne. 0.d0) then
!         write(77,*) "Top soft/coll event dipole"
!         write(77,*) "softness", soft
!         write(77,*) "kollinear", coll
!         write(77,*) "Dipol", Dip_Res_UnPol
!         write(77,*) "Real", NLO_Res_UnPol
!         write(77,*) "Quotient", NLO_Res_Unpol/Dip_Res_Unpol
!      endif



!      write(*,*) "softness", MomExt(1,11)/MomExt(1,4)
!      write(*,*) "kollinear", MomExt(1:4,11).dot.MomExt(1:4,8)/m_top**2
!      write(*,*) "Quotient", NLO_Res_Unpol/Dip_Res_Unpol+1.d0
!      write(*,*) "Dipol-Bin", NBin
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CON IIa2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CON IIa2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CON IIa2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CON IIa2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CON IIa2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CON IIa2
1818 continue
!----------------------------------
! photon emission off W-  with real correction  at anti-top |
!----------------------------------
  call EvalPhasespace_TopDecay4(MomExt(1:4,3),yRnd(5:14),MomExt(1:4,5:9),PSWgt2)
  call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(15:18),.false.,MomExt(1:4,10:12),PSWgt3)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt * PDFFac

!    call CheckSing(MomExt(1:4,9:12),applySingCut)
    if (PSwgt2 .eq. 0.d0 .or. PSwgt3 .eq. 0.d0) then
       goto 1881
    endif

    call Kinematics_TTBARPHOTON(1,MomExt(1:4,1:12),(/3,4,8,1,2,9,5,6,7,10,11,12/),applyPSCut,NBin)
    NLO_Res_Unpol = (0.d0,0.d0)

    if( applyPSCut ) then
       goto 811
    endif

      do iHel=1,NumHelicities ! loop over initial state chiralities
         do PhoHel=1,-1,-2 ! loop over additional photon chiralities
            do GluHel=1,-1,-2 ! loop over additional real gluon chiralities
               call HelCrossing(Helicities(iHel,1:NumExtParticles))
               call SetPolarizations()
               !*******************
               ! 1st CALL FOR INTERFERENCE TERMS
               ! Photon   from Anti-Top
               ! Gluon FROM TOP
               call TopDecay(ExtParticle(1),DKP_RE_L,MomExt(1:4,5:9),GluonHel=GluHel,PhotonHel=PhoHel)
               call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,10:12))
               ! Production Process calculated with ExtParticle(1) and ExtParticle(2)
               do iPrimAmp=1,NumBornAmps
                  call EvalTree(BornAmps(iPrimAmp))
                  DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
               enddo
               NLO_Res_Pol = (0.d0,0.d0)
               do jPrimAmp=1,NumBornAmps
                  do iPrimAmp=1,NumBornAmps
                     NLO_Res_Pol = NLO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*DKResult(iPrimAmp)* &
                          dconjg(DKResult(jPrimAmp))
                  enddo
               enddo
               NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
            enddo
         enddo
      enddo

      NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
      EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb + dble(NLO_Res_Unpol)

      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
      enddo

811  continue
!----------------------------------
! photon emission off W-, anti-top dipole-contribution for tb-> bb W- gluon
!----------------------------------
      MomExtTilde = MomExt
      MomExtTilde(1:4,9) = 0.d0
      call WTransform2(MomExt(1:4,5:9),MomExtTilde(1:4,5:8),pbDpg,ptDpg,ptDpb)
      omz=ptDpg/(ptDpb+ptDpg-pbDpg) !  for some reason this is not (1-z) as defined in the paper...
      rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
      z=1d0-omz
      y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2

      Dipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
      Dipole = Dipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )

      Dip_Res_UnPol =(0.d0,0.d0)
      call Kinematics_TTBARPHOTON(0,MomExtTilde(1:4,1:12),(/3,4,8,1,2,0,5,6,7,10,11,12/),applyPSCut,NBin)
      if( applyPSCut ) then
         goto 1881
      endif


      do iHel=1,NumHelicities ! loop over initial state chiralities
         do PhoHel=1,-1,-2 ! loop over additional photon chiralities
            call HelCrossing(Helicities(iHel,1:NumExtParticles))
            call SetPolarizations()
            ! Photon radiation from W-, Dipole from anti-top
            call TopDecay(ExtParticle(1),DKP_LO_L,MomExtTilde(1:4,5:8),PhotonHel=PhoHel)
            call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,10:12))
            do iPrimAmp=1,NumBornAmps
               call EvalTree(BornAmps(iPrimAmp))
               DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
            enddo
            LO_Res_Pol = (0.d0,0.d0)

            do jPrimAmp=1,NumBornAmps
               do iPrimAmp=1,NumBornAmps
                  LO_Res_Pol = LO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*DKResult(iPrimAmp)* &
                       dconjg(DKResult(jPrimAmp))
               enddo
            enddo
            Dip_Res_UnPol = Dip_Res_UnPol + LO_Res_Pol
         enddo
      enddo
      Dip_Res_UnPol = Dip_Res_UnPol* ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac * (dipole)

      EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb + dble(Dip_Res_Unpol)
!      soft =  MomExt(1,9)/m_top
!      coll = dsqrt(MomExt(1:4,9).dot.MomExt(1:4,5)/m_top**2)
!      if(soft .lt. 1.d-3 .or. coll .lt. 1.d-3) then
!         write(*,*) "soft", soft
!         write(*,*) "coll", coll
!         write(*,*) "real-atop",  NLO_Res_Unpol
!         write(*,*) "dip-atop",  dip_Res_Unpol
!         pause
!      endif
      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),dble(Dip_Res_Unpol))
      enddo

   if( IsNan(EvalCS_NLODKP_ttb )) then
        write(45,*) "ATOP-FacTest"
        write(45,*) "NAN:",EvalCS_NLODKP_ttb
         write(45,*) "yRnd(1) = ", yRnd(1)
         write(45,*) "yRnd(2) = ", yRnd(2)
         write(45,*) "yRnd(3) = ", yRnd(3)
         write(45,*) "yRnd(4) = ", yRnd(4)
         write(45,*) "yRnd(5) = ", yRnd(5)
         write(45,*) "yRnd(6) = ", yRnd(6)
         write(45,*) "yRnd(7) = ", yRnd(7)
         write(45,*) "yRnd(8) = ", yRnd(8)
         write(45,*) "yRnd(9) = ", yRnd(9)
         write(45,*) "yRnd(10) = ", yRnd(10)
         write(45,*) "yRnd(11) = ", yRnd(11)
         write(45,*) "yRnd(12) = ", yRnd(12)
         write(45,*) "yRnd(13) = ", yRnd(13)
         write(45,*) "yRnd(14) = ", yRnd(14)
         write(45,*) "yRnd(15) = ", yRnd(15)
         write(45,*) "yRnd(16) = ", yRnd(16)
         write(45,*) "yRnd(17) = ", yRnd(17)
         write(45,*) "yRnd(18) = ", yRnd(18)
        write(45,*) NLO_Res_UnPol,Dip_Res_Unpol
        write(45,*) PreFac , VgsWgt , PDFFac, sHatJacobi
        write(45,*) "Mom"
        write(45,*) MomExt(1:4,:)
        write(45,*) "SKIP EVENT!!!!!"
        EvalCS_NLODKP_ttb = 0d0
        return
   endif


1881 continue
!goto 6789
!----------------------------------
! photon emission off W+  with real correction  at the top quark|
!----------------------------------
  call EvalPhasespace_TopDecay4(MomExt(1:4,4),yRnd(5:14),MomExt(1:4,8:12),PSWgt2)
  call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(15:18),.false.,MomExt(1:4,5:7),PSWgt3)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt * PDFFac

!    call CheckSing(MomExt(1:4,5:8),applySingCut)
    if (PSwgt2 .eq. 0.d0 .or. PSwgt3 .eq. 0.d0) then
       goto 3333
    endif
!    write(*,*) "hall"
!   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt * PDFFac

    call Kinematics_TTBARPHOTON(1,MomExt(1:4,1:12),(/3,4,11,1,2,12,5,6,7,8,9,10/),applyPSCut,NBin)
    NLO_Res_Unpol = (0.d0,0.d0)
    if( applyPSCut ) then
      goto 188
!       write(*,*) "applypscut", applyPScut
!       stop
    endif


      do iHel=1,NumHelicities ! loop over initial state chiralities
         do PhoHel=1,-1,-2 ! loop over additional photon chiralities
            do GluHel=1,-1,-2 ! loop over additional real gluon chiralities
               call HelCrossing(Helicities(iHel,1:NumExtParticles))
               call SetPolarizations()
               !*******************
               ! 1st CALL FOR INTERFERENCE TERMS
               ! Gluon  from Top, photon from W+
               ! NO PHOTON FROM ANti-TOP
               call TopDecay(ExtParticle(2),DKP_RE_L,MomExt(1:4,8:12),GluonHel=GluHel,PhotonHel=PhoHel)
               call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,5:7))
               ! Production Process calculated with ExtParticle(1) and ExtParticle(2)
               do iPrimAmp=1,NumBornAmps
                  call EvalTree(BornAmps(iPrimAmp))
                  DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
               enddo
               NLO_Res_Pol = (0.d0,0.d0)
               do jPrimAmp=1,NumBornAmps
                  do iPrimAmp=1,NumBornAmps
                     NLO_Res_Pol = NLO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*DKResult(iPrimAmp)* &
                          dconjg(DKResult(jPrimAmp))
                  enddo
               enddo
               NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
            enddo
         enddo
      enddo
      NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
!      if(abs(dble(NLO_res_unpol)) .gt. 1.d-2) then
!         write(*,*) "real",NLO_Res_Unpol
!         write(*,*) "soft", MomExt(1,12)/m_top**2
!         write(*,*) "coll", (MomExt(1:4,12).dot.MomExt(1:4,8))/m_top**2
!         pause
!      endif
      EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb + dble(NLO_Res_Unpol)


      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
      enddo
!      write(*,*) "RealBin", NBin


188 continue
!----------------------------------
! photon emission off W+ with dipole-contribution  at the top quark
!----------------------------------
      MomExtTilde = MomExt
      MomExtTilde(1:4,12) = 0.d0
      call WTransform2(MomExt(1:4,8:12),MomExtTilde(1:4,8:11),pbDpg,ptDpg,ptDpb)
      omz=ptDpg/(ptDpb+ptDpg-pbDpg) !  for some reason this is not (1-z) as defined in the paper...
      rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
      z=1d0-omz
      y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2

      Dipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
      Dipole = Dipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )

      Dip_Res_UnPol =(0.d0,0.d0)
      call Kinematics_TTBARPHOTON(0,MomExtTilde(1:4,1:12),(/3,4,11,1,2,0,5,6,7,8,9,10/),applyPSCut,NBin)
      if (applyPScut) then
         goto 3333
      endif

      do iHel=1,NumHelicities ! loop over initial state chiralities
         do PhoHel=1,-1,-2 ! loop over additional photon chiralities
            call HelCrossing(Helicities(iHel,1:NumExtParticles))
            call SetPolarizations()
            !*******************
            call TopDecay(ExtParticle(2),DKP_LO_L,MomExtTilde(1:4,8:11),PhotonHel=PhoHel)
            call TopDecay(ExtParticle(1),DK_LO,MomExtTilde(1:4,5:7))

            do iPrimAmp=1,NumBornAmps
               call EvalTree(BornAmps(iPrimAmp))
               DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
            enddo
            LO_Res_Pol = (0.d0,0.d0)

            do jPrimAmp=1,NumBornAmps
               do iPrimAmp=1,NumBornAmps
                  LO_Res_Pol = LO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*DKResult(iPrimAmp)* &
                       dconjg(DKResult(jPrimAmp))
               enddo
            enddo
            Dip_Res_UnPol = Dip_Res_UnPol + LO_Res_Pol
         enddo
      enddo

      Dip_Res_UnPol = Dip_Res_UnPol* ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac * (dipole)
      EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb + dble(Dip_Res_Unpol)
!      soft =  MomExt(1,12)/m_top
!      coll = dsqrt(MomExt(1:4,12).dot.MomExt(1:4,8)/m_top**2)
!      if(soft .lt. 1.d-3 .or. coll .lt. 1.d-3) then
!         write(*,*) "soft", soft
!         write(*,*) "coll", coll
!         write(*,*) "real-top",  NLO_Res_Unpol
!         write(*,*) "dip-top",  dip_Res_Unpol
!         pause
!      endif

!      if(abs(dble(dip_res_unpol)) .gt. 1.d-2) then
!         write(*,*) "dip", Dip_res_unpol
!         pause
!      endif
      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),dble(Dip_Res_Unpol))
      enddo

      if( IsNan(EvalCS_NLODKP_ttb )) then
        write(55,*) "TOP-FacTest"
        write(55,*) "NAN:",EvalCS_NLODKP_ttb
         write(55,*) "yRnd(1) = ", yRnd(1)
         write(55,*) "yRnd(2) = ", yRnd(2)
         write(55,*) "yRnd(3) = ", yRnd(3)
         write(55,*) "yRnd(4) = ", yRnd(4)
         write(55,*) "yRnd(5) = ", yRnd(5)
         write(55,*) "yRnd(6) = ", yRnd(6)
         write(55,*) "yRnd(7) = ", yRnd(7)
         write(55,*) "yRnd(8) = ", yRnd(8)
         write(55,*) "yRnd(9) = ", yRnd(9)
         write(55,*) "yRnd(10) = ", yRnd(10)
         write(55,*) "yRnd(11) = ", yRnd(11)
         write(55,*) "yRnd(12) = ", yRnd(12)
         write(55,*) "yRnd(13) = ", yRnd(13)
         write(55,*) "yRnd(14) = ", yRnd(14)
         write(55,*) "yRnd(15) = ", yRnd(15)
         write(55,*) "yRnd(16) = ", yRnd(16)
         write(55,*) "yRnd(17) = ", yRnd(17)
         write(55,*) "yRnd(18) = ", yRnd(18)
        write(55,*) NLO_Res_UnPol,Dip_Res_Unpol
        write(55,*) PreFac, VgsWgt , PDFFac, sHatJacobi
        write(55,*) "Mom"
        write(55,*) MomExt(1:4,:)
        write(55,*) "SKIP EVENT!!!!!"
        EvalCS_NLODKP_ttb = 0d0
        return
   endif




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CON IIb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CON IIb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CON IIb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CON IIb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CON IIb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CON IIb

! !----------------------------------
! ! photon emission off anti-top    |
! !----------------------------------
! tb -> bb W- gamma at LO including gamma-radiation from W-, and t -> b W+  at NLO real
 3333 continue
!goto 4442
!goto 6789
 do nPhoRad=nPhoRad1,nPhoRad2!   nPhoRad=1: photon radiation off top/bot/W,
 !                 nPhoRad=2: photon radiation off W/lep
    if( nPhoRad.eq.1 ) then
       call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:11),.true.,MomExt(1:4,5:8),PSWgt2)
    else
       call EvalPhasespace_TopDecay2(MomExt(1:4,3),yRnd(5:11),.true.,MomExt(1:4,5:8),PSWgt2)
    endif
    call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(12:18),.true.,MomExt(1:4,9:12),PSWgt3)
    PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt * PDFFac

    call CheckSing(MomExt(1:4,9:12),applySingCut)
    if (applySingcut) then
       goto 7777
    endif
    if (PSWgt2 .eq. 0.d0 .or. PSWgt3 .eq. 0.d0) then
       goto 7777
    endif

    call Kinematics_TTBARPHOTON(1,MomExt(1:4,1:12),(/3,4,8,1,2,12,5,6,7,9,10,11/),applyPSCut,NBin)
      NLO_Res_Unpol = (0.d0,0.d0)

      if( applyPSCut ) then
         goto 177
      endif



      do iHel=1,NumHelicities ! loop over initial state chiralities
         do PhoHel=1,-1,-2 ! loop over additional photon chiralities
            do GluHel=1,-1,-2 ! loop over additional real gluon chiralities
               call HelCrossing(Helicities(iHel,1:NumExtParticles))
               call SetPolarizations()
               !*******************
               ! 1st CALL FOR INTERFERENCE TERMS
               ! Photon   from Anti-Top
               ! Gluon FROM TOP
               if(nPhoRad .eq. 1) THEN
                  call TopDecay(ExtParticle(1),DKP_LO_T,MomExt(1:4,5:8),PhotonHel=PhoHel)
               else
                  call TopDecay(ExtParticle(1),DKP_LO_L,MomExt(1:4,5:8),PhotonHel=PhoHel)
               endif
               call TopDecay(ExtParticle(2),DK_RE_T,MomExt(1:4,9:12),GluonHel=GluHel)
               ! Production Process calculated with ExtParticle(1) and ExtParticle(2)
               do iPrimAmp=1,NumBornAmps
                  call EvalTree(BornAmps(iPrimAmp))
                  DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
               enddo
               NLO_Res_Pol = (0.d0,0.d0)
               do jPrimAmp=1,NumBornAmps
                  do iPrimAmp=1,NumBornAmps
                     NLO_Res_Pol = NLO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*DKResult(iPrimAmp)* &
                          dconjg(DKResult(jPrimAmp))
                  enddo
               enddo
               NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
            enddo
         enddo
      enddo


      NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
      EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb + dble(NLO_Res_Unpol)
!      soft =  MomExt(1,8)/m_top
!      coll = dsqrt(MomExt(1:4,8).dot.MomExt(1:4,5)/m_top**2)

!      if(soft .le. 1.d-3 .or. coll .lt. 1.d-3) then
!      write(*,*) "Anti-Top soft/coll event dipole"
!         write(*,*) "softness", soft
!         write(*,*) "kollinear", coll
!         write(*,*) "Dipol", Dip_Res_UnPol
!         write(*,*) "Real", NLO_Res_UnPol
!         write(*,*) "Quotient", NLO_Res_Unpol/Dip_Res_Unpol
!         pause
!      endif


      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
      enddo




177  continue
!----------------------------------
! photon emission off anti-top, dipole-contribution for t-> b W+ gluon
!----------------------------------
      MomExtTilde = MomExt
      MomExtTilde(1:4,12) = 0.d0
      call WTransform(MomExt(1:4,9:12),MomExtTilde(1:4,9:11),pbDpg,ptDpg,ptDpb)
      omz=ptDpg/(ptDpb+ptDpg-pbDpg) !  for some reason this is not (1-z) as defined in the paper...
      rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
      z=1d0-omz
      y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2

      Dipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
      Dipole = Dipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )

      Dip_Res_UnPol =(0.d0,0.d0)
      call Kinematics_TTBARPHOTON(0,MomExtTilde(1:4,1:12),(/3,4,8,1,2,0,5,6,7,9,10,11/),applyPSCut,NBin)
      if( applyPSCut ) then
         cycle
      endif


      do iHel=1,NumHelicities ! loop over initial state chiralities
         do PhoHel=1,-1,-2 ! loop over additional photon chiralities
            call HelCrossing(Helicities(iHel,1:NumExtParticles))
            call SetPolarizations()
            !*******************
            ! HERE ONLY BORN-INTEFERENCE IS NEEDED
            ! SET MOMENTA FIRST SINCE MomExtTilde(1:4,12) is the gluon momentum and not-existing in reduced kinematics
            MomExtBornTilde(1:4,1) = MomExtTilde(1:4,9)
            MomExtBornTilde(1:4,2) = MomExtTilde(1:4,10)
            MomExtBornTilde(1:4,3) = MomExtTilde(1:4,11)
            ! Photon radiation from tbar/W-
            if (nPhoRad .eq. 1) then
               call TopDecay(ExtParticle(1),DKP_LO_T,MomExt(1:4,5:8),PhotonHel=PhoHel)
            else
               call TopDecay(ExtParticle(1),DKP_LO_L,MomExt(1:4,5:8),PhotonHel=PhoHel)
            endif
            ! Dipole from top
            call TopDecay(ExtParticle(2),DK_LO,MomExtBornTilde(1:4,1:3))
            do iPrimAmp=1,NumBornAmps
               call EvalTree(BornAmps(iPrimAmp))
               DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
            enddo
            LO_Res_Pol = (0.d0,0.d0)

            do jPrimAmp=1,NumBornAmps
               do iPrimAmp=1,NumBornAmps
                  LO_Res_Pol = LO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*DKResult(iPrimAmp)* &
                       dconjg(DKResult(jPrimAmp))
               enddo
            enddo
            Dip_Res_UnPol = Dip_Res_UnPol + LO_Res_Pol
         enddo
      enddo
      Dip_Res_UnPol = Dip_Res_UnPol* ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac * (dipole)

      EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb + dble(Dip_Res_Unpol)
      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),dble(Dip_Res_Unpol))
      enddo

!         soft =  MomExt(1,12)/m_top
!         coll = dsqrt(MomExt(1:4,12).dot.MomExt(1:4,9)/m_top**2)
!      if((soft .le. 1.d-3 .or. coll .lt. 1.d-3) .and. dble(Dip_Res_Unpol) .ne. 0.d0 .and.  dble(NLO_Res_UnPol) .eq. 0.d0) then
!         write(44,*) "IIb Top soft/coll event dipole"
!         write(44,*) "softness", soft
!         write(44,*) "kollinear", coll
!         write(44,*) "Dipol", Dip_Res_UnPol
!         write(44,*) "Real", NLO_Res_UnPol
!         write(44,*) "Quotient", NLO_Res_Unpol/Dip_Res_Unpol
!      endif
!      if((soft .le. 1.d-3 .or. coll .lt. 1.d-3) .and. dble(Dip_Res_Unpol) .eq. 0.d0 .and.  dble(NLO_Res_UnPol) .ne. 0.d0) then
!         write(55,*) "IIb Top soft/coll event dipole"
!         write(55,*) "softness", soft
!         write(55,*) "kollinear", coll
!         write(55,*) "Dipol", Dip_Res_UnPol
!         write(55,*) "Real", NLO_Res_UnPol
!         write(55,*) "Quotient", NLO_Res_Unpol/Dip_Res_Unpol
!      endif
!
!
   enddo!nPhoRad


7777 continue

!goto 6789
!   EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb/VgsWgt
!return
!4442 continue
!----------------------------------
! photon emission off top  with real correction  at the antitop quark|
!----------------------------------
 do nPhoRad=nPhoRad1,nPhoRad2!   nPhoRad=1: photon radiation off top/bot/W,
                             !   nPhoRad=2: photon radiation off W/lep
    if( nPhoRad.eq.1 ) then
       call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(12:18),.true.,MomExt(1:4,9:12),PSWgt2) ! Top-Momenta
    else
       call EvalPhasespace_TopDecay2(MomExt(1:4,4),yRnd(12:18),.true.,MomExt(1:4,9:12),PSWgt2)
    endif
    call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:11),.true.,MomExt(1:4,5:8),PSWgt3) ! A-Top-Momenta

    call CheckSing(MomExt(1:4,5:8),applySingCut)
    if (applySingcut) then
        EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb/VgsWgt
        return
    endif
    if (PSWgt2 .eq. 0.d0 .or. PSWgt3 .eq. 0.d0) then
        EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb/VgsWgt
        return
    endif

    PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt * PDFFac

    call Kinematics_TTBARPHOTON(1,MomExt(1:4,1:12),(/3,4,12,1,2,8,5,6,7,9,10,11/),applyPSCut,NBin)
    NLO_Res_Unpol = (0.d0,0.d0)
    if( applyPSCut ) then
      goto 7711
!       write(*,*) "applypscut", applyPScut
!       stop
    endif


      do iHel=1,NumHelicities ! loop over initial state chiralities
         do PhoHel=1,-1,-2 ! loop over additional photon chiralities
            do GluHel=1,-1,-2 ! loop over additional real gluon chiralities
               call HelCrossing(Helicities(iHel,1:NumExtParticles))
               call SetPolarizations()
               !*******************
               ! 1st CALL FOR INTERFERENCE TERMS
               ! Photon + Gluon  from Anti-Top
               ! NO PHOTON FROM TOP
               if(nPhoRad .eq. 1) then
                  call TopDecay(ExtParticle(2),DKP_LO_T,MomExt(1:4,9:12),PhotonHel=PhoHel)
               else
                  call TopDecay(ExtParticle(2),DKP_LO_L,MomExt(1:4,9:12),PhotonHel=PhoHel)
               endif
               call TopDecay(ExtParticle(1),DK_RE_T,MomExt(1:4,5:8),GluonHel=GluHel)
               ! Production Process calculated with ExtParticle(1) and ExtParticle(2)
               do iPrimAmp=1,NumBornAmps
                  call EvalTree(BornAmps(iPrimAmp))
                  DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
               enddo
               NLO_Res_Pol = (0.d0,0.d0)
               do jPrimAmp=1,NumBornAmps
                  do iPrimAmp=1,NumBornAmps
                     NLO_Res_Pol = NLO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*DKResult(iPrimAmp)* &
                          dconjg(DKResult(jPrimAmp))
                  enddo
               enddo
               NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
            enddo
         enddo
      enddo
      NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac
      EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb + dble(NLO_Res_Unpol)


      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
      enddo
!      write(*,*) "RealBin", NBin




7711 continue
!----------------------------------
! photon emission off top with dipole-contribution  at the antitop quark
!----------------------------------
      MomExtTilde = MomExt
      MomExtTilde(1:4,8) = 0.d0
      call WTransform(MomExt(1:4,5:8),MomExtTilde(1:4,5:7),pbDpg,ptDpg,ptDpb)
      omz=ptDpg/(ptDpb+ptDpg-pbDpg) !  for some reason this is not (1-z) as defined in the paper...
      rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
      z=1d0-omz
      y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2

      Dipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
      Dipole = Dipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )

      Dip_Res_UnPol =(0.d0,0.d0)

      call Kinematics_TTBARPHOTON(0,MomExtTilde(1:4,1:12),(/3,4,12,1,2,0,5,6,7,9,10,11/),applyPSCut,NBin)

      if (applyPScut) then
         cycle
      endif

      do iHel=1,NumHelicities ! loop over initial state chiralities
         do PhoHel=1,-1,-2 ! loop over additional photon chiralities
            call HelCrossing(Helicities(iHel,1:NumExtParticles))
            call SetPolarizations()
            !*******************
            ! HERE ONLY BORN-INTEFERENCE IS NEEDED
            ! SET MOMENTA FIRST SINCE MomExtTilde(1:4,11) is the gluon momentum and no-existing in reduced kinematics
            MomExtBornTilde(1:4,1) = MomExtTilde(1:4,5)
            MomExtBornTilde(1:4,2) = MomExtTilde(1:4,6)
            MomExtBornTilde(1:4,3) = MomExtTilde(1:4,7)
            if(nPhoRad .eq. 1) then
               call TopDecay(ExtParticle(2),DKP_LO_T,MomExt(1:4,9:12),PhotonHel=PhoHel)
            else
               call TopDecay(ExtParticle(2),DKP_LO_L,MomExt(1:4,9:12),PhotonHel=PhoHel)
            endif
            call TopDecay(ExtParticle(1),DK_LO,MomExtBornTilde(1:4,1:3))

            do iPrimAmp=1,NumBornAmps
               call EvalTree(BornAmps(iPrimAmp))
               DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
            enddo
            LO_Res_Pol = (0.d0,0.d0)

            do jPrimAmp=1,NumBornAmps
               do iPrimAmp=1,NumBornAmps
                  LO_Res_Pol = LO_Res_Pol + ColCorrLO(iPrimAmp,jPrimAmp) * ParityFlip*DKResult(iPrimAmp)* &
                       dconjg(DKResult(jPrimAmp))
               enddo
            enddo
            Dip_Res_UnPol = Dip_Res_UnPol + LO_Res_Pol
         enddo
      enddo

      Dip_Res_UnPol = Dip_Res_UnPol* ISFac * (alpha_s4Pi*RunFactor)**2 * PreFac * (dipole)
      EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb + dble(Dip_Res_Unpol)


      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),dble(Dip_Res_Unpol))
      enddo

!      soft =  MomExt(1,8)/m_top
!      coll = dsqrt(MomExt(1:4,5).dot.MomExt(1:4,8)/m_top**2)
!      if((soft .le. 1.d-3 .or. coll .lt. 1.d-3) .and. dble(Dip_Res_Unpol) .ne. 0.d0 .and.  dble(NLO_Res_UnPol) .eq. 0.d0) then
!         write(44,*) "IIb Anti-Top soft/coll event dipole"
!         write(44,*) "softness", soft
!         write(44,*) "kollinear", coll
!         write(44,*) "Dipol", Dip_Res_UnPol
!         write(44,*) "Real", NLO_Res_UnPol
!         write(44,*) "Quotient", NLO_Res_Unpol/Dip_Res_Unpol
!      endif
!      if((soft .le. 1.d-3 .or. coll .lt. 1.d-3) .and. dble(Dip_Res_Unpol) .eq. 0.d0 .and.  dble(NLO_Res_UnPol) .ne. 0.d0) then
!         write(44,*) "IIb Anti-Top soft/coll event dipole"
!         write(55,*) "softness", soft
!         write(55,*) "kollinear", coll
!         write(55,*) "Dipol", Dip_Res_UnPol
!         write(55,*) "Real", NLO_Res_UnPol
!         write(55,*) "Quotient", NLO_Res_Unpol/Dip_Res_Unpol
!      endif

!      if(soft .le. 1.d-3 .or. coll .lt. 1.d-3) then
!         write(*,*) "nPhoRad", nPhorad
!         write(*,*) "Top soft/coll event"
!         write(*,*) "softness", soft
!         write(*,*) "kollinear", coll
!         write(*,*) "Dipol", Dip_Res_UnPol
!         write(*,*) "Real", NLO_Res_UnPol
!         write(*,*) "Quotient", NLO_Res_Unpol/Dip_Res_Unpol
!
!         pause

!      endif



!      write(*,*) "softness", MomExt(1,11)/MomExt(1,4)
!      write(*,*) "kollinear", MomExt(1:4,11).dot.MomExt(1:4,8)/m_top**2
!      write(*,*) "Quotient", NLO_Res_Unpol/Dip_Res_Unpol+1.d0
!      write(*,*) "Dipol-Bin", NBin
!
   enddo!nPhoRad


ENDIF
6789 continue

   EvalCS_NLODKP_ttb = EvalCS_NLODKP_ttb/VgsWgt
RETURN
END FUNCTION








FUNCTION EvalCS_DKP_1L_ttbgg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModIntDipoles_DKP_GGTTBG
implicit none
real(8) ::  EvalCS_DKP_1L_ttbgg,EvalCS_DKP_1L_ttbgg_1,EvalCS_DKP_1L_ttbgg_2
real(8) ::  yRnd(1:VegasMxDim),VgsWgt,IOp(-2:0),HOp(1:3)
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1),FermionLoopPartAmp(7:8,-2:1)
integer :: iHel,jHel,kHel,iPrimAmp,jPrimAmp,nPhoRad,PhoHel
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,ISFac
real(8) :: MomExt(1:4,1:12),MomP(1:4,1:4)
logical :: applyPSCut
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,AccPoles
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2),xE
integer :: NBin(1:NumMaxHisto),NHisto,ParityFlip,nHel(1:2)
include 'misc/global_import'
include 'vegas_common.f'

  ParityFlip=1
  EvalCS_DKP_1L_ttbgg = 0d0
  EvalCS_DKP_1L_ttbgg_1 = 0d0
  EvalCS_DKP_1L_ttbgg_2 = 0d0

  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top ) then
      EvalCS_DKP_1L_ttbgg = 0d0
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

!----------------------------------
! photon emission off anti-top    |
!----------------------------------
do nPhoRad=nPhoRad1,nPhoRad2!   nPhoRad=1: photon radiation off top/bot/W, nPhoRad=2: photon radiation off W/lep
   LO_Res_Unpol             = (0d0,0d0)
   NLO_Res_Unpol(-2:1)      = (0d0,0d0)
   NLO_Res_Unpol_Ferm(-2:1) = (0d0,0d0)
   if( nPhoRad.eq.1 ) then
      call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:11),.true.,MomExt(1:4,5:8),PSWgt2)
   else
      call EvalPhasespace_TopDecay2(MomExt(1:4,3),yRnd(5:11),.true.,MomExt(1:4,5:8),PSWgt2)
   endif
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(12:15),.false.,MomExt(1:4,9:11),PSWgt3)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt*PSWgt2*PSWgt3 * VgsWgt * PDFFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))

   call Kinematics_TTBARPHOTON(0,MomExt(1:4,1:12),(/3,4,8,1,2,0,5,6,7,9,10,11/),applyPSCut,NBin)
   if( applyPSCut ) then
      cycle
   endif
!------------ LO --------------
IF( Correction.EQ.0 ) THEN
   do iHel=nHel(1),nHel(2)
   do PhoHel=1,-1,-2 ! loop over additional photon polarization
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      if( nPhoRad.eq.1 ) then
         call TopDecay(ExtParticle(1),DKP_LO_T,MomExt(1:4,5:8),PhotonHel=PhoHel)
      else
         call TopDecay(ExtParticle(1),DKP_LO_L,MomExt(1:4,5:8),PhotonHel=PhoHel)
      endif
      call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))

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

!------------ 1 LOOP --------------
ELSEIF( Correction.EQ.1 ) THEN
   do iHel=nHel(1),nHel(2)
   do PhoHel=1,-1,-2 ! loop over additional photon polarization
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      if( nPhoRad.eq.1 ) then
         call TopDecay(ExtParticle(1),DKP_LO_T,MomExt(1:4,5:8),PhotonHel=PhoHel)
      else
         call TopDecay(ExtParticle(1),DKP_LO_L,MomExt(1:4,5:8),PhotonHel=PhoHel)
      endif
      call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))

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
                  EvalCS_DKP_1L_ttbgg = 0d0
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
!           call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,rdiv(2),rdiv(1))
!           call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv)
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
   enddo! helicity loop
ENDIF


IF( Correction.EQ.0 ) THEN
!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * WidthExpansion
   EvalCS_DKP_1L_ttbgg = LO_Res_Unpol * PreFac

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

   EvalCS_DKP_1L_ttbgg = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac


ELSEIF( Correction.EQ.3 ) THEN
! print *, nPhoRad
! print *, "1-loop eps2:",(NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2) )* PreFac,  (NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2))/(alpha_sOver2Pi*RunFactor*LO_Res_Unpol)
! print *, "1-loop eps1:",(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1) )* PreFac,  (NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/(alpha_sOver2Pi*RunFactor*LO_Res_Unpol)
! print *, "tree virt",LO_Res_Unpol * alpha_sOver2Pi*RunFactor* PreFac


   MomP(1:4,1) = MomExt(1:4,3)
   MomP(1:4,2) = MomExt(1:4,4)
   MomP(1:4,3) =-MomExt(1:4,1)
   MomP(1:4,4) =-MomExt(1:4,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt*PSWgt2*PSWgt3 * VgsWgt
   xE = yRnd(16+HelSampling)
! xE=0.45d0

   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
   call EvalIntDipoles_DKP_GGTTBG(MomP(1:4,1:4),MomExt(1:4,5:11),ATop_,nPhoRad,xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac

   EvalCS_DKP_1L_ttbgg = HOp(1)    * pdf(0,1)  * pdf(0,2)   &
                      + HOp(2)/xE * pdf_z(0,1)* pdf(0,2)   &
                      + HOp(3)/xE * pdf(0,1)  * pdf_z(0,2)
ENDIF

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_DKP_1L_ttbgg)
   enddo

EvalCS_DKP_1L_ttbgg_1 = EvalCS_DKP_1L_ttbgg_1 + EvalCS_DKP_1L_ttbgg

! print *,EvalCS_DKP_1L_ttbgg
! pause
enddo! nPhoRad loop


! EvalCS_DKP_1L_ttbgg = (EvalCS_DKP_1L_ttbgg_1)/VgsWgt
! return
! -----------------------------------------------------------------------------


!----------------------------------
! photon emission off top         |
!----------------------------------
do nPhoRad=nPhoRad1,nPhoRad2!   nPhoRad=1: photon radiation off top/bot/W,nPhoRad=2: photon radiation off W/lep
   LO_Res_Unpol             = (0d0,0d0)
   NLO_Res_Unpol(-2:1)      = (0d0,0d0)
   NLO_Res_Unpol_Ferm(-2:1) = (0d0,0d0)
   call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:8),.false.,MomExt(1:4,5:7),PSWgt2)
   if( nPhoRad.eq.1 ) then
      call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(9:15),.true.,MomExt(1:4,8:11),PSWgt3)
   else
      call EvalPhasespace_TopDecay2(MomExt(1:4,4),yRnd(9:15),.true.,MomExt(1:4,8:11),PSWgt3)
   endif
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt*PSWgt2*PSWgt3 * VgsWgt * PDFFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))

   call Kinematics_TTBARPHOTON(0,MomExt(1:4,1:12),(/3,4,11,1,2,0,5,6,7,8,9,10/),applyPSCut,NBin)
   if( applyPSCut ) then
      cycle
   endif

!------------ LO --------------
IF( Correction.EQ.0 ) THEN
   do iHel=nHel(1),nHel(2)
   do PhoHel=1,-1,-2 ! loop over additional photon polarization
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,5:7))
      if( nPhoRad.eq.1 ) then
         call TopDecay(ExtParticle(2),DKP_LO_T,MomExt(1:4,8:11),PhotonHel=PhoHel)
      else
         call TopDecay(ExtParticle(2),DKP_LO_L,MomExt(1:4,8:11),PhotonHel=PhoHel)
      endif

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

!------------ 1 LOOP --------------
ELSEIF( Correction.EQ.1 ) THEN
   do iHel=nHel(1),nHel(2)
   do PhoHel=1,-1,-2 ! loop over additional photon polarization
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,5:7))
      if( nPhoRad.eq.1 ) then
         call TopDecay(ExtParticle(2),DKP_LO_T,MomExt(1:4,8:11),PhotonHel=PhoHel)
      else
         call TopDecay(ExtParticle(2),DKP_LO_L,MomExt(1:4,8:11),PhotonHel=PhoHel)
      endif

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
                  EvalCS_DKP_1L_ttbgg = 0d0
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
!           call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,rdiv(2),rdiv(1))
!           call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv)
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
   enddo! helicity loop
ENDIF


IF( Correction.EQ.0 ) THEN
!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * WidthExpansion
   EvalCS_DKP_1L_ttbgg = LO_Res_Unpol * PreFac

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

   EvalCS_DKP_1L_ttbgg = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac


ELSEIF( Correction.EQ.3 ) THEN
! print *, nPhoRad
! print *, "1-loop eps2:",(NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2) )* PreFac,  (NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2))/(alpha_sOver2Pi*RunFactor*LO_Res_Unpol)
! print *, "1-loop eps1:",(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1) )* PreFac,  (NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/(alpha_sOver2Pi*RunFactor*LO_Res_Unpol)
! print *, "tree virt",LO_Res_Unpol * alpha_sOver2Pi*RunFactor* PreFac


   MomP(1:4,1) = MomExt(1:4,3)
   MomP(1:4,2) = MomExt(1:4,4)
   MomP(1:4,3) =-MomExt(1:4,1)
   MomP(1:4,4) =-MomExt(1:4,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt*PSWgt2*PSWgt3 * VgsWgt

      xE = yRnd(16+HelSampling)
! xE=0.55d0
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

   call EvalIntDipoles_DKP_GGTTBG(MomP(1:4,1:4),MomExt(1:4,5:11),Top_,nPhoRad,xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac

   EvalCS_DKP_1L_ttbgg = HOp(1)    * pdf(0,1)  * pdf(0,2)   &
                      + HOp(2)/xE * pdf_z(0,1)* pdf(0,2)   &
                      + HOp(3)/xE * pdf(0,1)  * pdf_z(0,2)
ENDIF

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_DKP_1L_ttbgg)
   enddo

   EvalCS_DKP_1L_ttbgg_2 = EvalCS_DKP_1L_ttbgg_2 + EvalCS_DKP_1L_ttbgg

! print *,EvalCS_DKP_1L_ttbgg
! pause

enddo! nPhoRad loop



   EvalCS_DKP_1L_ttbgg = (EvalCS_DKP_1L_ttbgg_1+EvalCS_DKP_1L_ttbgg_2)/VgsWgt
return
END FUNCTION






FUNCTION EvalCS_DKP_1L_ttbqqb(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModIntDipoles_DKP_QQBTTBG
use ModIntDipoles_DKP_QGTTBQ
implicit none
real(8) :: EvalCS_DKP_1L_ttbqqb,EvalCS_DKP_1L_ttbqqb_1,EvalCS_DKP_1L_ttbqqb_2
real(8) :: yRnd(1:VegasMxDim),VgsWgt,xE,HOp(1:3)
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1)
complex(8) :: BosonicPartAmp(1:2,-2:1),FermionPartAmp(1:2,-2:1)
integer :: iHel,iPrimAmp,jPrimAmp,nPhoRad,PhoHel
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,ISFac
real(8) :: MomExt(1:4,1:12),MomP(1:4,1:4)
logical :: applyPSCut
real(8),parameter :: Nc=3d0
real(8) :: tau,eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac_a,PDFFac_b,PDFFac,AccPoles
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2)
integer :: NHisto,NBin(1:NumMaxHisto),ParityFlip,npdf,nHel(1:2),NRndHel
include 'misc/global_import'
include "vegas_common.f"


  ParityFlip=1
  EvalCS_DKP_1L_ttbqqb = 0d0
  EvalCS_DKP_1L_ttbqqb_1 = 0d0
  EvalCS_DKP_1L_ttbqqb_2 = 0d0

  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top ) then
      EvalCS_DKP_1L_ttbqqb = 0d0
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


!----------------------------------
! photon emission off anti-top    |
!----------------------------------
do nPhoRad=nPhoRad1,nPhoRad2!   nPhoRad=1: photon radiation off top/bot/W, nPhoRad=2: photon radiation off W/lep
   LO_Res_Unpol             = (0d0,0d0)
   NLO_Res_Unpol(-2:1)      = (0d0,0d0)
   NLO_Res_Unpol_Ferm(-2:1) = (0d0,0d0)
   if( nPhoRad.eq.1 ) then
      call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:11),.true.,MomExt(1:4,5:8),PSWgt2)
   else
      call EvalPhasespace_TopDecay2(MomExt(1:4,3),yRnd(5:11),.true.,MomExt(1:4,5:8),PSWgt2)
   endif
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(12:15),.false.,MomExt(1:4,9:11),PSWgt3)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt*PSWgt2*PSWgt3 * VgsWgt * dble(NumHelicities/(nHel(2)-nHel(1)+1))

   call Kinematics_TTBARPHOTON(0,MomExt(1:4,1:12),(/3,4,8,1,2,0,5,6,7,9,10,11/),applyPSCut,NBin)
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
    do PhoHel=1,-1,-2 ! loop over additional photon polarization
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
        call SetPolarizations()
        if( nPhoRad.eq.1 ) then
          call TopDecay(ExtParticle(1),DKP_LO_T,MomExt(1:4,5:8),PhotonHel=PhoHel)
        else
          call TopDecay(ExtParticle(1),DKP_LO_L,MomExt(1:4,5:8),PhotonHel=PhoHel)
        endif
        call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))

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
    do PhoHel=1,-1,-2 ! loop over additional photon polarization
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      if( nPhoRad.eq.1 ) then
         call TopDecay(ExtParticle(1),DKP_LO_T,MomExt(1:4,5:8),PhotonHel=PhoHel)
      else
         call TopDecay(ExtParticle(1),DKP_LO_L,MomExt(1:4,5:8),PhotonHel=PhoHel)
      endif
      call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))
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
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + ParityFlip*Col1L_ttbqqb(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(BosonicPartAmp(jPrimAmp,-2:1)) )
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
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + ParityFlip*Col1L_ttbqqb(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(FermionPartAmp(jPrimAmp,-2:1)) )
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
   EvalCS_DKP_1L_ttbqqb = LO_Res_Unpol * PreFac

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

   EvalCS_DKP_1L_ttbqqb = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac

ELSEIF( CORRECTION.EQ.3 ) THEN

! xE=0.23d0
! print *, "PhoRad",nPhoRad
! print *, "xE=",xE

   MomP(1:4,1) = MomExt(1:4,3)
   MomP(1:4,2) = MomExt(1:4,4)
   MomP(1:4,3) =-MomExt(1:4,1)
   MomP(1:4,4) =-MomExt(1:4,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt*PSWgt2*PSWgt3 * VgsWgt
   xE = yRnd(16+HelSampling)

! print *, "-2",(NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2))* PreFac
! print *, "-1",(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))* PreFac


IF( PROCESS.EQ.31 ) THEN
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

      call EvalIntDipoles_DKP_QQBTTBG(MomP(1:4,1:4),MomExt(1:4,5:11),ATop_,nPhoRad,xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_DKP_1L_ttbqqb = HOp(1) * (pdf(Up_,1)*pdf(AUp_,2)+pdf(Dn_,1)*pdf(ADn_,2)+pdf(Chm_,1)*pdf(AChm_,2)+pdf(Str_,1)*pdf(AStr_,2)+pdf(Bot_,1)*pdf(ABot_,2) ) &
                          + HOp(2)/xE * (pdf_z(Up_,1)*pdf(AUp_,2)+pdf_z(Dn_,1)*pdf(ADn_,2)+pdf_z(Chm_,1)*pdf(AChm_,2)+pdf_z(Str_,1)*pdf(AStr_,2)+pdf_z(Bot_,1)*pdf(ABot_,2) ) &
                          + HOp(3)/xE * (pdf(Up_,1)*pdf_z(AUp_,2)+pdf(Dn_,1)*pdf_z(ADn_,2)+pdf(Chm_,1)*pdf_z(AChm_,2)+pdf(Str_,1)*pdf_z(AStr_,2)+pdf(Bot_,1)*pdf_z(ABot_,2) )

      call swapMom(MomP(1:4,3),MomP(1:4,4))
      call EvalIntDipoles_DKP_QQBTTBG(MomP(1:4,1:4),MomExt(1:4,5:11),ATop_,nPhoRad,xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_DKP_1L_ttbqqb = EvalCS_DKP_1L_ttbqqb &
                           + HOp(1)    * (pdf(Up_,2)*pdf(AUp_,1)+pdf(Dn_,2)*pdf(ADn_,1)+pdf(Chm_,2)*pdf(AChm_,1)+pdf(Str_,2)*pdf(AStr_,1)+pdf(Bot_,2)*pdf(ABot_,1) ) &
                           + HOp(2)/xE * (pdf_z(Up_,2)*pdf(AUp_,1)+pdf_z(Dn_,2)*pdf(ADn_,1)+pdf_z(Chm_,2)*pdf(AChm_,1)+pdf_z(Str_,2)*pdf(AStr_,1)+pdf_z(Bot_,2)*pdf(ABot_,1) ) &
                           + HOp(3)/xE * (pdf(Up_,2)*pdf_z(AUp_,1)+pdf(Dn_,2)*pdf_z(ADn_,1)+pdf(Chm_,2)*pdf_z(AChm_,1)+pdf(Str_,2)*pdf_z(AStr_,1)+pdf(Bot_,2)*pdf_z(ABot_,1) )
      call swapMom(MomP(1:4,3),MomP(1:4,4))! swap back

ELSEIF( PROCESS.EQ.25 ) THEN
   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

   call EvalIntDipoles_DKP_QGTTBQ(MomP(1:4,1:4),MomExt(1:4,5:11),ATop_,nPhoRad,xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
   EvalCS_DKP_1L_ttbqqb = HOp(1) * (pdf(Up_,1)*pdf(0,2)+pdf(Dn_,1)*pdf(0,2)+pdf(Chm_,1)*pdf(0,2)+pdf(Str_,1)*pdf(0,2)+pdf(Bot_,1)*pdf(0,2) ) &
                       + HOp(2)/xE * (pdf_z(Up_,1)*pdf(0,2)+pdf_z(Dn_,1)*pdf(0,2)+pdf_z(Chm_,1)*pdf(0,2)+pdf_z(Str_,1)*pdf(0,2)+pdf_z(Bot_,1)*pdf(0,2) ) &
                       + HOp(3)/xE * (pdf(Up_,1)*pdf_z(0,2)+pdf(Dn_,1)*pdf_z(0,2)+pdf(Chm_,1)*pdf_z(0,2)+pdf(Str_,1)*pdf_z(0,2)+pdf(Bot_,1)*pdf_z(0,2) )

   call swapMom(MomP(1:4,3),MomP(1:4,4))
   call EvalIntDipoles_DKP_QGTTBQ(MomP(1:4,1:4),MomExt(1:4,5:11),ATop_,nPhoRad,xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
   EvalCS_DKP_1L_ttbqqb = EvalCS_DKP_1L_ttbqqb &
                       + HOp(1)    * (pdf(Up_,2)*pdf(0,1)+pdf(Dn_,2)*pdf(0,1)+pdf(Chm_,2)*pdf(0,1)+pdf(Str_,2)*pdf(0,1)+pdf(Bot_,2)*pdf(0,1) ) &
                       + HOp(2)/xE * (pdf_z(Up_,2)*pdf(0,1)+pdf_z(Dn_,2)*pdf(0,1)+pdf_z(Chm_,2)*pdf(0,1)+pdf_z(Str_,2)*pdf(0,1)+pdf_z(Bot_,2)*pdf(0,1) ) &
                       + HOp(3)/xE * (pdf(Up_,2)*pdf_z(0,1)+pdf(Dn_,2)*pdf_z(0,1)+pdf(Chm_,2)*pdf_z(0,1)+pdf(Str_,2)*pdf_z(0,1)+pdf(Bot_,2)*pdf_z(0,1) )
   call swapMom(MomP(1:4,3),MomP(1:4,4))! swap back

ELSEIF( PROCESS.EQ.27 ) THEN
   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

   call EvalIntDipoles_DKP_QGTTBQ(MomP(1:4,1:4),MomExt(1:4,5:11),ATop_,nPhoRad,xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
   EvalCS_DKP_1L_ttbqqb = HOp(1) * (pdf(AUp_,1)*pdf(0,2)+pdf(ADn_,1)*pdf(0,2)+pdf(AChm_,1)*pdf(0,2)+pdf(AStr_,1)*pdf(0,2)+pdf(ABot_,1)*pdf(0,2) ) &
                       + HOp(2)/xE * (pdf_z(AUp_,1)*pdf(0,2)+pdf_z(ADn_,1)*pdf(0,2)+pdf_z(AChm_,1)*pdf(0,2)+pdf_z(AStr_,1)*pdf(0,2)+pdf_z(ABot_,1)*pdf(0,2) ) &
                       + HOp(3)/xE * (pdf(AUp_,1)*pdf_z(0,2)+pdf(ADn_,1)*pdf_z(0,2)+pdf(AChm_,1)*pdf_z(0,2)+pdf(AStr_,1)*pdf_z(0,2)+pdf(ABot_,1)*pdf_z(0,2) )

   call swapMom(MomP(1:4,3),MomP(1:4,4))
   call EvalIntDipoles_DKP_QGTTBQ(MomP(1:4,1:4),MomExt(1:4,5:11),ATop_,nPhoRad,xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
   EvalCS_DKP_1L_ttbqqb = EvalCS_DKP_1L_ttbqqb &
                       + HOp(1)    * (pdf(AUp_,2)*pdf(0,1)+pdf(ADn_,2)*pdf(0,1)+pdf(AChm_,2)*pdf(0,1)+pdf(AStr_,2)*pdf(0,1)+pdf(ABot_,2)*pdf(0,1) ) &
                       + HOp(2)/xE * (pdf_z(AUp_,2)*pdf(0,1)+pdf_z(ADn_,2)*pdf(0,1)+pdf_z(AChm_,2)*pdf(0,1)+pdf_z(AStr_,2)*pdf(0,1)+pdf_z(ABot_,2)*pdf(0,1) ) &
                       + HOp(3)/xE * (pdf(AUp_,2)*pdf_z(0,1)+pdf(ADn_,2)*pdf_z(0,1)+pdf(AChm_,2)*pdf_z(0,1)+pdf(AStr_,2)*pdf_z(0,1)+pdf(ABot_,2)*pdf_z(0,1) )
   call swapMom(MomP(1:4,3),MomP(1:4,4))! swap back

ENDIF
ENDIF

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_DKP_1L_ttbqqb)
   enddo

EvalCS_DKP_1L_ttbqqb_1 = EvalCS_DKP_1L_ttbqqb_1 + EvalCS_DKP_1L_ttbqqb
enddo! nPhoRad loop


!    EvalCS_DKP_1L_ttbqqb = (EvalCS_DKP_1L_ttbqqb_1)/VgsWgt
!    return
!-----------------------------------------------------------------------------


!----------------------------------
! photon emission off top    |
!----------------------------------
do nPhoRad=nPhoRad1,nPhoRad2!   nPhoRad=1: photon radiation off top/bot/W, nPhoRad=2: photon radiation off W/lep
   LO_Res_Unpol             = (0d0,0d0)
   NLO_Res_Unpol(-2:1)      = (0d0,0d0)
   NLO_Res_Unpol_Ferm(-2:1) = (0d0,0d0)
   call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:8),.false.,MomExt(1:4,5:7),PSWgt2)
   if( nPhoRad.eq.1 ) then
      call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(9:15),.true.,MomExt(1:4,8:11),PSWgt3)
   else
      call EvalPhasespace_TopDecay2(MomExt(1:4,4),yRnd(9:15),.true.,MomExt(1:4,8:11),PSWgt3)
   endif
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt*PSWgt2*PSWgt3 * VgsWgt * dble(NumHelicities/(nHel(2)-nHel(1)+1))

   call Kinematics_TTBARPHOTON(0,MomExt(1:4,1:12),(/3,4,11,1,2,0,5,6,7,8,9,10/),applyPSCut,NBin)
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
    do PhoHel=1,-1,-2 ! loop over additional photon polarization
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
        call SetPolarizations()
        call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,5:7))
        if( nPhoRad.eq.1 ) then
          call TopDecay(ExtParticle(2),DKP_LO_T,MomExt(1:4,8:11),PhotonHel=PhoHel)
        else
          call TopDecay(ExtParticle(2),DKP_LO_L,MomExt(1:4,8:11),PhotonHel=PhoHel)
        endif

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
    do PhoHel=1,-1,-2 ! loop over additional photon polarization
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,5:7))
      if( nPhoRad.eq.1 ) then
          call TopDecay(ExtParticle(2),DKP_LO_T,MomExt(1:4,8:11),PhotonHel=PhoHel)
      else
          call TopDecay(ExtParticle(2),DKP_LO_L,MomExt(1:4,8:11),PhotonHel=PhoHel)
      endif

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
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + ParityFlip*Col1L_ttbqqb(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(BosonicPartAmp(jPrimAmp,-2:1)) )
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
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + ParityFlip*Col1L_ttbqqb(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(FermionPartAmp(jPrimAmp,-2:1)) )
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
   EvalCS_DKP_1L_ttbqqb = LO_Res_Unpol * PreFac

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

   EvalCS_DKP_1L_ttbqqb = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac


ELSEIF( CORRECTION.EQ.3 ) THEN
   MomP(1:4,1) = MomExt(1:4,3)
   MomP(1:4,2) = MomExt(1:4,4)
   MomP(1:4,3) =-MomExt(1:4,1)
   MomP(1:4,4) =-MomExt(1:4,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt*PSWgt2*PSWgt3 * VgsWgt
   xE = yRnd(16+HelSampling)

! xE=0.23d0
! print *, "PhoRad",nPhoRad
! print *, "xE=",xE
! print *, "-2",(NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2))* PreFac
! print *, "-1",(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))* PreFac

IF( PROCESS.EQ.31 ) THEN
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

      call EvalIntDipoles_DKP_QQBTTBG(MomP(1:4,1:4),MomExt(1:4,5:11),Top_,nPhoRad,xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_DKP_1L_ttbqqb = HOp(1) * (pdf(Up_,1)*pdf(AUp_,2)+pdf(Dn_,1)*pdf(ADn_,2)+pdf(Chm_,1)*pdf(AChm_,2)+pdf(Str_,1)*pdf(AStr_,2)+pdf(Bot_,1)*pdf(ABot_,2) ) &
                          + HOp(2)/xE * (pdf_z(Up_,1)*pdf(AUp_,2)+pdf_z(Dn_,1)*pdf(ADn_,2)+pdf_z(Chm_,1)*pdf(AChm_,2)+pdf_z(Str_,1)*pdf(AStr_,2)+pdf_z(Bot_,1)*pdf(ABot_,2) ) &
                          + HOp(3)/xE * (pdf(Up_,1)*pdf_z(AUp_,2)+pdf(Dn_,1)*pdf_z(ADn_,2)+pdf(Chm_,1)*pdf_z(AChm_,2)+pdf(Str_,1)*pdf_z(AStr_,2)+pdf(Bot_,1)*pdf_z(ABot_,2) )

      call swapMom(MomP(1:4,3),MomP(1:4,4))
      call EvalIntDipoles_DKP_QQBTTBG(MomP(1:4,1:4),MomExt(1:4,5:11),Top_,nPhoRad,xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_DKP_1L_ttbqqb = EvalCS_DKP_1L_ttbqqb &
                           + HOp(1)    * (pdf(Up_,2)*pdf(AUp_,1)+pdf(Dn_,2)*pdf(ADn_,1)+pdf(Chm_,2)*pdf(AChm_,1)+pdf(Str_,2)*pdf(AStr_,1)+pdf(Bot_,2)*pdf(ABot_,1) ) &
                           + HOp(2)/xE * (pdf_z(Up_,2)*pdf(AUp_,1)+pdf_z(Dn_,2)*pdf(ADn_,1)+pdf_z(Chm_,2)*pdf(AChm_,1)+pdf_z(Str_,2)*pdf(AStr_,1)+pdf_z(Bot_,2)*pdf(ABot_,1) ) &
                           + HOp(3)/xE * (pdf(Up_,2)*pdf_z(AUp_,1)+pdf(Dn_,2)*pdf_z(ADn_,1)+pdf(Chm_,2)*pdf_z(AChm_,1)+pdf(Str_,2)*pdf_z(AStr_,1)+pdf(Bot_,2)*pdf_z(ABot_,1) )
      call swapMom(MomP(1:4,3),MomP(1:4,4))! swap back

ELSEIF( PROCESS.EQ.25 ) THEN
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

      call EvalIntDipoles_DKP_QGTTBQ(MomP(1:4,1:4),MomExt(1:4,5:11),Top_,nPhoRad,xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_DKP_1L_ttbqqb = HOp(1) * (pdf(Up_,1)*pdf(0,2)+pdf(Dn_,1)*pdf(0,2)+pdf(Chm_,1)*pdf(0,2)+pdf(Str_,1)*pdf(0,2)+pdf(Bot_,1)*pdf(0,2) ) &
                          + HOp(2)/xE * (pdf_z(Up_,1)*pdf(0,2)+pdf_z(Dn_,1)*pdf(0,2)+pdf_z(Chm_,1)*pdf(0,2)+pdf_z(Str_,1)*pdf(0,2)+pdf_z(Bot_,1)*pdf(0,2) ) &
                          + HOp(3)/xE * (pdf(Up_,1)*pdf_z(0,2)+pdf(Dn_,1)*pdf_z(0,2)+pdf(Chm_,1)*pdf_z(0,2)+pdf(Str_,1)*pdf_z(0,2)+pdf(Bot_,1)*pdf_z(0,2) )

      call swapMom(MomP(1:4,3),MomP(1:4,4))
      call EvalIntDipoles_DKP_QGTTBQ(MomP(1:4,1:4),MomExt(1:4,5:11),Top_,nPhoRad,xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_DKP_1L_ttbqqb = EvalCS_DKP_1L_ttbqqb &
                          + HOp(1)    * (pdf(Up_,2)*pdf(0,1)+pdf(Dn_,2)*pdf(0,1)+pdf(Chm_,2)*pdf(0,1)+pdf(Str_,2)*pdf(0,1)+pdf(Bot_,2)*pdf(0,1) ) &
                          + HOp(2)/xE * (pdf_z(Up_,2)*pdf(0,1)+pdf_z(Dn_,2)*pdf(0,1)+pdf_z(Chm_,2)*pdf(0,1)+pdf_z(Str_,2)*pdf(0,1)+pdf_z(Bot_,2)*pdf(0,1) ) &
                          + HOp(3)/xE * (pdf(Up_,2)*pdf_z(0,1)+pdf(Dn_,2)*pdf_z(0,1)+pdf(Chm_,2)*pdf_z(0,1)+pdf(Str_,2)*pdf_z(0,1)+pdf(Bot_,2)*pdf_z(0,1) )
      call swapMom(MomP(1:4,3),MomP(1:4,4))! swap back

ELSEIF( PROCESS.EQ.27 ) THEN
      call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

      call EvalIntDipoles_DKP_QGTTBQ(MomP(1:4,1:4),MomExt(1:4,5:11),Top_,nPhoRad,xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_DKP_1L_ttbqqb = HOp(1) * (pdf(AUp_,1)*pdf(0,2)+pdf(ADn_,1)*pdf(0,2)+pdf(AChm_,1)*pdf(0,2)+pdf(AStr_,1)*pdf(0,2)+pdf(ABot_,1)*pdf(0,2) ) &
                          + HOp(2)/xE * (pdf_z(AUp_,1)*pdf(0,2)+pdf_z(ADn_,1)*pdf(0,2)+pdf_z(AChm_,1)*pdf(0,2)+pdf_z(AStr_,1)*pdf(0,2)+pdf_z(ABot_,1)*pdf(0,2) ) &
                          + HOp(3)/xE * (pdf(AUp_,1)*pdf_z(0,2)+pdf(ADn_,1)*pdf_z(0,2)+pdf(AChm_,1)*pdf_z(0,2)+pdf(AStr_,1)*pdf_z(0,2)+pdf(ABot_,1)*pdf_z(0,2) )

      call swapMom(MomP(1:4,3),MomP(1:4,4))
      call EvalIntDipoles_DKP_QGTTBQ(MomP(1:4,1:4),MomExt(1:4,5:11),Top_,nPhoRad,xE,HOp(1:3))
      HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
      EvalCS_DKP_1L_ttbqqb = EvalCS_DKP_1L_ttbqqb &
                          + HOp(1)    * (pdf(AUp_,2)*pdf(0,1)+pdf(ADn_,2)*pdf(0,1)+pdf(AChm_,2)*pdf(0,1)+pdf(AStr_,2)*pdf(0,1)+pdf(ABot_,2)*pdf(0,1) ) &
                          + HOp(2)/xE * (pdf_z(AUp_,2)*pdf(0,1)+pdf_z(ADn_,2)*pdf(0,1)+pdf_z(AChm_,2)*pdf(0,1)+pdf_z(AStr_,2)*pdf(0,1)+pdf_z(ABot_,2)*pdf(0,1) ) &
                          + HOp(3)/xE * (pdf(AUp_,2)*pdf_z(0,1)+pdf(ADn_,2)*pdf_z(0,1)+pdf(AChm_,2)*pdf_z(0,1)+pdf(AStr_,2)*pdf_z(0,1)+pdf(ABot_,2)*pdf_z(0,1) )
      call swapMom(MomP(1:4,3),MomP(1:4,4))! swap back

ENDIF
ENDIF

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_DKP_1L_ttbqqb)
   enddo

   EvalCS_DKP_1L_ttbqqb_2 = EvalCS_DKP_1L_ttbqqb_2 + EvalCS_DKP_1L_ttbqqb
enddo! nPhoRad loop

   EvalCS_DKP_1L_ttbqqb = (EvalCS_DKP_1L_ttbqqb_1+EvalCS_DKP_1L_ttbqqb_2)/VgsWgt

return
END FUNCTION








FUNCTION EvalCS_DKP_Real_ttbggg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModDipoles_DKP_GGTTBG
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModJPsiFrag
implicit none
real(8) :: EvalCS_DKP_Real_ttbggg,EvalCS_DKP_Real_ttbggg_1,EvalCS_DKP_Real_ttbggg_2
real(8) :: yRnd(1:VegasMxDim),VgsWgt
complex(8) :: LO_Res_Pol,LO_Res_Unpol
integer :: iHel,iPrimAmp,jPrimAmp,nPhoRad,PhoHel
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,DipoleResult,ISFac
real(8) :: MomExt(1:4,1:12)
logical :: applyPSCut,applySingCut
real(8) :: tau,eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac
real(8) :: pdf(-6:6,1:2)
integer :: NBin(1:NumMaxHisto),NHisto,ParityFlip
real(8) :: s13,s23,s12,s55
include "vegas_common.f"

! print *, "yrnd fixed"
! yrnd(1)=0.108672696449176d0
! yrnd(2)=0.319632736183532d0
! yrnd(3)=8.270438624671873d-002
! yrnd(4)=  0.122117752033341d0
! yrnd(5)=  0.370806539138224d0
! yrnd(6)=  0.112870671606097d0
! yrnd(7)=  0.222487953827944d0
! yrnd(8)=  0.329036836898437d0
! yrnd(9)=  0.124350515484973d0
! yrnd(10)=  0.157048880428564d0
! yrnd(11)=  0.433058424355955d0
! yrnd(12)=  0.123806856164619d0
! yrnd(13)=  0.291977992184450d0
! yrnd(14)=  0.454681325449926d0
! yrnd(15)=  6.672940988406983d-002
! yrnd(16)=  0.201397819771151d0
! yrnd(17)=  3.388324265083451d-002
!  yrnd(18)= 0.326361954829824d0



  ParityFlip=1
  EvalCS_DKP_Real_ttbggg = 0d0
  EvalCS_DKP_Real_ttbggg_1 = 0d0
  EvalCS_DKP_Real_ttbggg_2 = 0d0

  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top ) then
      EvalCS_DKP_Real_ttbggg = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to3(EHat,yRnd(3:7),MomExt(1:4,1:5),PSWgt)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))
   ISFac = MomCrossing(MomExt)


   call CheckSing(MomExt,applySingCut)
   if( applySingCut ) then
       EvalCS_DKP_Real_ttbggg = 0d0
       return
   endif

   call setPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
   RunFactor = RunAlphaS(2,MuRen)


!----------------------------------
! photon emission off anti-top    |
!----------------------------------
   EvalCS_DKP_Real_ttbggg = 0d0
do nPhoRad=nPhoRad1,nPhoRad2!   nPhoRad=1: photon radiation off top/bot/W,nPhoRad=2: photon radiation off W/lep
   if( nPhoRad.eq.1 ) then
      call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:14),.true.,MomExt(1:4,6:9),PSWgt2)
   else
      call EvalPhasespace_TopDecay2(MomExt(1:4,4),yRnd(8:14),.true.,MomExt(1:4,6:9),PSWgt2)
   endif
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(15:18),.false.,MomExt(1:4,10:12),PSWgt3)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt*PSWgt2*PSWgt3 * VgsWgt * PDFFac

   call Kinematics_TTBARPHOTON(1,MomExt(1:4,1:12),(/4,5,9,1,2,3,6,7,8,10,11,12/),applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_DKP_Real_ttbggg = 0d0
      goto 50
   endif

        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
        do PhoHel=1,-1,-2 ! loop over additional photon polarization
          call HelCrossing(Helicities(iHel,1:NumExtParticles))
          call SetPolarizations()
          if( nPhoRad.eq.1 ) then
            call TopDecay(ExtParticle(1),DKP_LO_T,MomExt(1:4,6:9),PhotonHel=PhoHel)
          else
            call TopDecay(ExtParticle(1),DKP_LO_L,MomExt(1:4,6:9),PhotonHel=PhoHel)
          endif
          call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,10:12))

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
      enddo!helicity loop
      EvalCS_DKP_Real_ttbggg = LO_Res_UnPol * PreFac * (alpha_s4Pi*RunFactor)**3 * ISFac

      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),EvalCS_DKP_Real_ttbggg)
      enddo


50 continue


      call EvalDipoles_DKP_GGTTBG((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:18),PreFac,ATop_,nPhoRad,DipoleResult)
      EvalCS_DKP_Real_ttbggg_1 = EvalCS_DKP_Real_ttbggg_1 + (EvalCS_DKP_Real_ttbggg+DipoleResult)


! ! print *, yrnd(1:18)
!   s12 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,2))
!   s13 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,3))
!   print *, s13/s12,EvalCS_DKP_Real_ttbggg,DipoleResult,EvalCS_DKP_Real_ttbggg/DipoleResult+1d0
!   pause
enddo! nPhoRad loop





!----------------------------------
! photon emission off top    |
!----------------------------------
   EvalCS_DKP_Real_ttbggg = 0d0
do nPhoRad=nPhoRad1,nPhoRad2!   nPhoRad=1: photon radiation off top/bot/W, nPhoRad=2: photon radiation off W/lep
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomExt(1:4,6:8),PSWgt2)
   if( nPhoRad.eq.1 ) then
      call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:18),.true.,MomExt(1:4,9:12),PSWgt3)
   else
      call EvalPhasespace_TopDecay2(MomExt(1:4,5),yRnd(12:18),.true.,MomExt(1:4,9:12),PSWgt3)
   endif
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt*PSWgt2*PSWgt3 * VgsWgt * PDFFac

   call Kinematics_TTBARPHOTON(1,MomExt(1:4,1:12),(/4,5,12,1,2,3,6,7,8,9,10,11/),applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_DKP_Real_ttbggg = 0d0
      goto 51
   endif

        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities! helicity summation
        do PhoHel=1,-1,-2 ! loop over additional photon polarization
          call HelCrossing(Helicities(iHel,1:NumExtParticles))
          call SetPolarizations()

          call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
          if( nPhoRad.eq.1 ) then
            call TopDecay(ExtParticle(2),DKP_LO_T,MomExt(1:4,9:12),PhotonHel=PhoHel)
          else
            call TopDecay(ExtParticle(2),DKP_LO_L,MomExt(1:4,9:12),PhotonHel=PhoHel)
          endif

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
      enddo!helicity loop
      EvalCS_DKP_Real_ttbggg = LO_Res_UnPol * PreFac * (alpha_s4Pi*RunFactor)**3 * ISFac

      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),EvalCS_DKP_Real_ttbggg)
      enddo


51 continue


    call EvalDipoles_DKP_GGTTBG((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:18),PreFac,Top_,nPhoRad,DipoleResult)
    EvalCS_DKP_Real_ttbggg_2 = EvalCS_DKP_Real_ttbggg_2 + (EvalCS_DKP_Real_ttbggg+DipoleResult)

! ! print *, yrnd(1:18)
!   s12 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,2))
!   s13 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,3))
!   print *, s13/s12,EvalCS_DKP_Real_ttbggg,DipoleResult,EvalCS_DKP_Real_ttbggg/DipoleResult+1d0
!   pause

enddo! nPhoRad loop




EvalCS_DKP_Real_ttbggg = (EvalCS_DKP_Real_ttbggg_1+EvalCS_DKP_Real_ttbggg_2)/VgsWgt
if( isNaN(EvalCS_DKP_Real_ttbggg) ) then
    print *, "NaN event"
    print *,  EvalCS_DKP_Real_ttbggg
    print *, EvalCS_DKP_Real_ttbggg_1,EvalCS_DKP_Real_ttbggg_2,VgsWgt
    print *, yRnd(:)
    print *, ""
    EvalCS_DKP_Real_ttbggg = 0d0
endif



RETURN
END FUNCTION








FUNCTION EvalCS_DKP_Real_ttbqqbg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModDipoles_DKP_QQBTTBG
use ModDipoles_DKP_QGTTBQ
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModJPsiFrag
implicit none
real(8) ::  EvalCS_DKP_Real_ttbqqbg,EvalCS_DKP_Real_ttbqqbg_1,EvalCS_DKP_Real_ttbqqbg_2
real(8) ::  EvalCS_DK_Dips_ttbqqbg,EvalCS_DK_Dips_ttbqqbg_1,EvalCS_DK_Dips_ttbqqbg_2
real(8) ::  yRnd(1:VegasMxDim),VgsWgt,DipoleResult
complex(8) :: LO_Res_Pol,LO_Res_Unpol
integer :: iHel,iPrimAmp,jPrimAmp,nPhoRad,PhoHel
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,s12,s13
real(8) :: PDFFac_a,PDFFac_b,PDFFac,pdf(-6:6,1:2),tau,eta1,eta2,FluxFac,sHatJacobi,ISFac,PreFac
real(8) :: MomExt(1:4,1:12)
integer :: NBin(1:NumMaxHisto),NHisto,ParityFlip,npdf
logical :: applyPSCut,applySingCut
include "vegas_common.f"

  ParityFlip=1
  EvalCS_DKP_Real_ttbqqbg = 0d0
  EvalCS_DKP_Real_ttbqqbg_1 = 0d0
  EvalCS_DKP_Real_ttbqqbg_2 = 0d0
  EvalCS_DK_Dips_ttbqqbg = 0d0
  EvalCS_DK_Dips_ttbqqbg_1 = 0d0
  EvalCS_DK_Dips_ttbqqbg_2 = 0d0

  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top ) then
      EvalCS_DKP_Real_ttbqqbg = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to3(EHat,yRnd(3:7),MomExt(1:4,1:5),PSWgt)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))

   call CheckSing(MomExt,applySingCut)
   if( applySingCut ) then
      EvalCS_DKP_Real_ttbqqbg = 0d0
      return
   endif


   call setPDFs(eta1,eta2,MuFac,pdf)
IF( PROCESS.EQ.31 ) THEN
   PDFFac_a = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
            + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
            + pdf(Bot_,1)*pdf(ABot_,2)
   PDFFac_b = pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
            + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
            + pdf(Bot_,2)*pdf(ABot_,1)
ELSEIF( PROCESS.EQ.25 ) THEN
   PDFFac_a = pdf(Up_,1) *pdf(0,2) + pdf(Dn_,1) *pdf(0,2)  &
            + pdf(Chm_,1)*pdf(0,2) + pdf(Str_,1)*pdf(0,2)  &
            + pdf(Bot_,1)*pdf(0,2)
   PDFFac_b = pdf(Up_,2) *pdf(0,1) + pdf(Dn_,2) *pdf(0,1)  &
            + pdf(Chm_,2)*pdf(0,1) + pdf(Str_,2)*pdf(0,1)  &
            + pdf(Bot_,2)*pdf(0,1)
ELSEIF( PROCESS.EQ.27 ) THEN
   PDFFac_a = pdf(AUp_,1) *pdf(0,2) + pdf(ADn_,1) *pdf(0,2)  &
            + pdf(AChm_,1)*pdf(0,2) + pdf(AStr_,1)*pdf(0,2)  &
            + pdf(ABot_,1)*pdf(0,2)
   PDFFac_b = pdf(AUp_,2) *pdf(0,1) + pdf(ADn_,2) *pdf(0,1)  &
            + pdf(AChm_,2)*pdf(0,1) + pdf(AStr_,2)*pdf(0,1)  &
            + pdf(ABot_,2)*pdf(0,1)
ENDIF
   RunFactor = RunAlphaS(2,MuRen)


!----------------------------------
! photon emission off anti-top    |
!----------------------------------
EvalCS_DKP_Real_ttbqqbg_1 = 0d0
do nPhoRad=nPhoRad1,nPhoRad2!   nPhoRad=1: photon radiation off top/bot/W, nPhoRad=2: photon radiation off W/lep
   if( nPhoRad.eq.1 ) then
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

        call Kinematics_TTBARPHOTON(1,MomExt(1:4,1:12),(/4,5,9,1,2,3,6,7,8,10,11,12/),applyPSCut,NBin)
        if( applyPSCut ) then
            EvalCS_DKP_Real_ttbqqbg = 0d0
            goto 60
        endif

        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities
        do PhoHel=1,-1,-2 ! loop over additional photon polarization
          call HelCrossing(Helicities(iHel,1:NumExtParticles))
          call SetPolarizations()

          if( nPhoRad.eq.1 ) then
            call TopDecay(ExtParticle(1),DKP_LO_T,MomExt(1:4,6:9),PhotonHel=PhoHel)
          else
            call TopDecay(ExtParticle(1),DKP_LO_L,MomExt(1:4,6:9),PhotonHel=PhoHel)
          endif
          call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,10:12))

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
      enddo!helicity loop
      EvalCS_DKP_Real_ttbqqbg = LO_Res_UnPol * PreFac * (alpha_s4Pi*RunFactor)**3 * ISFac

      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),EvalCS_DKP_Real_ttbqqbg)
      enddo


60    continue


      IF( PROCESS.EQ.31 ) THEN
          call EvalDipoles_DKP_QQBTTBG((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:18),PreFac,ATop_,nPhoRad,DipoleResult)
      ELSEIF( PROCESS.EQ.25 ) THEN
          call EvalDipoles_DKP_QGTTBQ((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:18),PreFac,ATop_,nPhoRad,DipoleResult)
      ELSEIF( PROCESS.EQ.27 ) THEN
          call EvalDipoles_DKP_QGTTBQ((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:18),PreFac,ATop_,nPhoRad,DipoleResult)
      ENDIF

      EvalCS_DKP_Real_ttbqqbg_1 = EvalCS_DKP_Real_ttbqqbg_1 + (EvalCS_DKP_Real_ttbqqbg + DipoleResult)

!   print *, nPdf,nPhoRad
!   s12 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,2))
!   s13 = 2d0*(MomExt(1:4,2).dot.MomExt(1:4,3))
!   print *, s13/s12,EvalCS_DKP_Real_ttbqqbg,DipoleResult,EvalCS_DKP_Real_ttbqqbg/DipoleResult+1d0
!   pause

enddo! nPdf loop
call swapMom(MomExt(1:4,1),MomExt(1:4,2))! swap back
enddo! nPhoRad loop




!----------------------------------
! photon emission off top         |
!----------------------------------
  EvalCS_DKP_Real_ttbqqbg_2 = 0d0
do nPhoRad=nPhoRad1,nPhoRad2!   nPhoRad=1: photon radiation off top/bot/W, nPhoRad=2: photon radiation off W/lep
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomExt(1:4,6:8),PSWgt2)
   if( nPhoRad.eq.1 ) then
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

        call Kinematics_TTBARPHOTON(1,MomExt(1:4,1:12),(/4,5,12,1,2,3,6,7,8,9,10,11/),applyPSCut,NBin)
        if( applyPSCut ) then
            EvalCS_DKP_Real_ttbqqbg = 0d0
            goto 61
        endif

        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities
        do PhoHel=1,-1,-2 ! loop over additional photon polarization
          call HelCrossing(Helicities(iHel,1:NumExtParticles))
          call SetPolarizations()
          call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
          if( nPhoRad.eq.1 ) then
            call TopDecay(ExtParticle(2),DKP_LO_T,MomExt(1:4,9:12),PhotonHel=PhoHel)
          else
            call TopDecay(ExtParticle(2),DKP_LO_L,MomExt(1:4,9:12),PhotonHel=PhoHel)
          endif

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
      enddo!helicity loop
      EvalCS_DKP_Real_ttbqqbg = LO_Res_UnPol * PreFac * (alpha_s4Pi*RunFactor)**3 * ISFac

      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),EvalCS_DKP_Real_ttbqqbg)
      enddo


61 continue


      IF( PROCESS.EQ.31 ) THEN
          call EvalDipoles_DKP_QQBTTBG((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:18),PreFac,Top_,nPhoRad,DipoleResult)
      ELSEIF( PROCESS.EQ.25 ) THEN
          call EvalDipoles_DKP_QGTTBQ((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:18),PreFac,Top_,nPhoRad,DipoleResult)
      ELSEIF( PROCESS.EQ.27 ) THEN
          call EvalDipoles_DKP_QGTTBQ((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/),(/0d0,0d0,m_Top**2,m_Top**2,0d0/),yRnd(8:18),PreFac,Top_,nPhoRad,DipoleResult)
      ENDIF

!       print *, nPdf,nPhoRad
!       s12 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,2))
!       s13 = 2d0*(MomExt(1:4,2).dot.MomExt(1:4,3))
!       print *, s13/s12,EvalCS_DKP_Real_ttbqqbg,DipoleResult,EvalCS_DKP_Real_ttbqqbg/DipoleResult+1d0
!       pause

      EvalCS_DKP_Real_ttbqqbg_2 = EvalCS_DKP_Real_ttbqqbg_2 + (EvalCS_DKP_Real_ttbqqbg + DipoleResult)

enddo! nPdf loop
call swapMom(MomExt(1:4,1),MomExt(1:4,2))! swap back
enddo! nPhoRad loop


EvalCS_DKP_Real_ttbqqbg = (EvalCS_DKP_Real_ttbqqbg_1+EvalCS_DKP_Real_ttbqqbg_2)/VgsWgt

if( isNaN(EvalCS_DKP_Real_ttbqqbg) ) then
    print *, "NaN event"
    print *, EvalCS_DKP_Real_ttbqqbg
    print *, EvalCS_DKP_Real_ttbqqbg_1,EvalCS_DKP_Real_ttbqqbg_2,VgsWgt
    print *, yRnd(:)
    print *, ""
    EvalCS_DKP_Real_ttbqqbg = 0d0
endif
RETURN
END FUNCTION





END MODULE ModCrossSection_TTBP


