  MODULE ModCrossSection_TTBZ
use ModTopDecay
implicit none

integer,private,parameter :: NumMaxHisto=45


contains



FUNCTION EvalCS_1L_ttbggZ(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_new
use ModUCuts_128
use ModUCuts_128_new
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModZDecay
use ModIntDipoles_GGTTBGZ
implicit none
real(8) ::  EvalCS_1L_ttbggZ,yRnd(1:VegasMxDim),VgsWgt
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1),FermionLoopPartAmp(1:3,-2:1)
complex(8) :: BosonicPartAmp(1:3,-2:1),mydummy,ZPolVec(1:4),BarSpi(1:4),Spi(1:4)
integer :: iHel,jHel,kHel,iPrimAmp,jPrimAmp,BPrimAmp,APrimAmp,ListPrimAmps(14)
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,PSWgt4,ISFac,ZDKMatel,Msq_T_BWENU
real(8) :: MomExt(1:4,1:14)
logical :: applyPSCut
real(8) :: Col1Lf_ttbggZ(2,2), Col1L_ttbggZ(2,3)
real(8) :: MG_MOM(0:3,1:5),tmpmom(1:4)
real(8) :: MadGraph_tree
real(8),parameter :: Nc=3d0, Cf=4d0/3d0
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,AccPoles
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2),xE,HOp(1:3),MZ_Inv
real(8) :: QPtol,DPtol
integer :: NBin(1:NumMaxHisto),NHisto,PhotonCouplCorr=2d0,nHel(1:2),NRndHel
integer :: ZQcoupl,jj,lastSister
integer :: QPredo(1:8,1:5)
include 'misc/global_import'
include 'vegas_common.f'



! RR : here you can set the coupling of the Z to the light quarks in the currents
! ZQcoupl=1           ! left- and right-handed (default)
! ZQcoupl=2           ! up and down
 ZQcoupl=3           !vector and axial-vector
DPtol=1d-4
QPtol=1d-4


EvalCS_1L_ttbggZ = 0d0

   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)


   if( EHat.le.2d0*m_Top+M_Z ) then
      EvalCS_1L_ttbggZ = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)

   if( ZDecays.le.10 ) then  ! decaying on-shell Z
      MZ_Inv = m_Z
   elseif( ZDecays.gt.10 ) then  ! decaying off-shell Z
      call Error("need to implement phase space for off-shell Z's")
      ! need to think about threshold cut: EHat.le.2d0*m_Top+M_Z   when z is off-shell! 
   endif

   call EvalPhaseSpace_2to3M(EHat,MZ_Inv,yRnd(3:7),MomExt(1:4,1:5),PSWgt)
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
IF( ZDECAYS.NE.0 .AND. ZDECAYS.NE.-2 ) THEN!  -2: Z=photon
   call EvalPhasespace_ZDecay(MZ_Inv,MomExt(1:4,3),yRnd(16:17),MomExt(1:4,12:13),PSWgt4)
   PSWgt = PSWgt * PSWgt4
ENDIF

IF( ZDECAYS.EQ.-1 ) THEN!   spin un-correlated decays
     ExtParticle(5)%Helicity=+1
     call ZDecay(ExtParticle(5),DK_LO,MomExt(1:4,12:13))
     call pol_massSR(ExtParticle(5)%Mom(1:4),M_Z,+1,ZPolVec(1:4))
     ZDKMatel = cdabs((ExtParticle(5)%Pol(1:4).dot.ZPolVec(1:4)))**2
     call pol_massSR(ExtParticle(5)%Mom(1:4),M_Z,-1,ZPolVec(1:4))
     ZDKMatel = ZDKMatel + cdabs((ExtParticle(5)%Pol(1:4).dot.ZPolVec(1:4)))**2
     call pol_massSR(ExtParticle(5)%Mom(1:4),M_Z,0,ZPolVec(1:4))
     ZDKMatel = ZDKMatel + cdabs((ExtParticle(5)%Pol(1:4).dot.ZPolVec(1:4)))**2

     ExtParticle(5)%Helicity=-1
     call ZDecay(ExtParticle(5),DK_LO,MomExt(1:4,12:13))
     call pol_massSR(ExtParticle(5)%Mom(1:4),M_Z,+1,ZPolVec(1:4))
     ZDKMatel = ZDKMatel + cdabs((ExtParticle(5)%Pol(1:4).dot.ZPolVec(1:4)))**2
     call pol_massSR(ExtParticle(5)%Mom(1:4),M_Z,-1,ZPolVec(1:4))
     ZDKMatel = ZDKMatel + cdabs((ExtParticle(5)%Pol(1:4).dot.ZPolVec(1:4)))**2
     call pol_massSR(ExtParticle(5)%Mom(1:4),M_Z,0,ZPolVec(1:4))
     ZDKMatel = ZDKMatel + cdabs((ExtParticle(5)%Pol(1:4).dot.ZPolVec(1:4)))**2

     PSWgt = PSWgt * ZDKMatel/3d0

      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
      Spi(1:4) = ExtParticle(1)%Pol(1:4)
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,-1,BarSpi(1:4))
      Msq_T_BWENU = cdabs(psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top))**2/2d0
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,+1,BarSpi(1:4))
      Msq_T_BWENU = Msq_T_BWENU + cdabs(psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top))**2/2d0

      PSWgt = PSWgt * Msq_T_BWENU


      call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))
      BarSpi(1:4) = ExtParticle(2)%Pol(1:4)
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,-1,Spi(1:4))
      Msq_T_BWENU = cdabs(psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top))**2/2d0
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,+1,Spi(1:4))
      Msq_T_BWENU = Msq_T_BWENU + cdabs(psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top))**2/2d0

      PSWgt = PSWgt * Msq_T_BWENU
ENDIF


   call Kinematics_TTBARZ(0,MomExt(1:4,1:14),(/4,5,3,1,2,0,6,7,8,9,10,11,12,13/),applyPSCut,NBin)
!   tmpmom(1:4)=ExtParticle(3)%Mom(1:4)
!   ExtParticle(3)%Mom(1:4)=ExtParticle(4)%Mom(1:4)
!   ExtParticle(4)%Mom(1:4)=tmpmom(1:4)
! this is a troublesome point
!ExtParticle(1)%Mom(1)  = (2.32342593857319d0,0d0)
!ExtParticle(1)%Mom(2) = (0.237295808622009d0,0.000000000000000d+000)
!ExtParticle(1)%Mom(3) = (5.101878197379041d-003,0.000000000000000d+000)
!ExtParticle(1)%Mom(4) = (-1.51786816360446d0,0.000000000000000d+000)
!ExtParticle(2)%Mom(1) = (2.02856223854126d0,0.000000000000000d+000)
!ExtParticle(2)%Mom(2) = (-0.262317688261449d0,0.000000000000000d+000)
!ExtParticle(2)%Mom(3) = (-6.934418703624884d-002,0.000000000000000d+000)
!ExtParticle(2)%Mom(4) = (-1.00169684525058d0,0.000000000000000d+000)
!ExtParticle(3)%Mom(1) = (-1.02489198372927d0,0.000000000000000d+000)
!ExtParticle(3)%Mom(2) = (0.000000000000000d+000,0.000000000000000d+000)
!ExtParticle(3)%Mom(3) = (0.000000000000000d+000,0.000000000000000d+000)
!ExtParticle(3)%Mom(4) = (-1.02489198372927d0,0.000000000000000d+000)
!ExtParticle(4)%Mom(1) = (-5.35957053786124d0,0.00000000000000d+000)
!ExtParticle(4)%Mom(2) = (0.000000000000000d+000,0.000000000000000d+000)
!ExtParticle(4)%Mom(3) = (0.00000000000000d+000,0.000000000000000d+000)
!ExtParticle(4)%Mom(4) = (5.35957053786124d0,0.000000000000000d+000)
!ExtParticle(5)%Mom(1) = (2.03247434447607d0,0.000000000000000d+000)
!ExtParticle(5)%Mom(2) = (2.502187963944023d-002,0.000000000000000d+000)
!ExtParticle(5)%Mom(3) = (6.424230883886980d-002,0.000000000000000d+000)
!ExtParticle(5)%Mom(4) = (-1.81511354527693d0,0.000000000000000d+000)

!   print *, 'p1', ExtParticle(1)%Mom(1:4)
!   print *, 'p2', ExtParticle(2)%Mom(1:4)
!   print *, 'p3', ExtParticle(3)%Mom(1:4)
!   print *, 'p4', ExtParticle(4)%Mom(1:4)
!   print *, 'p5', ExtParticle(5)%Mom(1:4)
!   print *, '========================================='
!   if (ZQcoupl .eq. 1) then
!      print * , 'using left and right handed couplings in the currents' 
!   elseif (ZQcoupl .eq. 2) then  
!      print * , 'using up or down  couplings in the currents' 
!   elseif (ZQcoupl .eq.3) then  
!      print *, 'using vector and then axial-vector couplings in the currents'
!   else
!      print *, 'Error: ZQcoupl not set' 
!      stop
!   endif
!   print *, '========================================='
!
   if( applyPSCut ) then
      EvalCS_1L_ttbggZ = 0d0
      return
   endif

   call SetPropagators()
   call SetPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
!   print *, yRnd(1:2)
!   pause
!   return
!   print *, 'ehat', Ehat
!   print *, 'fbGev2,FluxFac , sHatJacobi , PSWgt , VgsWgt , PDFFac',fbGev2,FluxFac , sHatJacobi , PSWgt , VgsWgt , PDFFac
!   pause
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   RunFactor = RunAlphaS(NLOParam,MuRen)
   nHel(1:2) = getHelicity(yrnd(NRndHel))
 !  PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))

   LO_Res_Unpol             = (0d0,0d0)
   NLO_Res_Unpol(-2:1)      = (0d0,0d0)
   NLO_Res_Unpol_Ferm(-2:1) = (0d0,0d0)
!------------ LO --------------
IF( Correction.EQ.0 ) THEN
!   do iHel=nHel(1),nHel(2)
    do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()

      if( ZDecays.gt.0 ) then
          if( ExtParticle(5)%Helicity.eq.0 ) cycle!   this can be done more elegantly in mod_process
          call ZDecay(ExtParticle(5),DK_LO,MomExt(1:4,12:13))

      endif
      do iPrimAmp=1,2!  set this explicitely to avoid a problem with warm-up run for Correction=1 
          call EvalTree(BornAmps(iPrimAmp))
      enddo
      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,2
      do iPrimAmp=1,2
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
   enddo!helicity loop


!------------ 1 LOOP --------------
ELSEIF( Correction.EQ.1 ) THEN

!   do iHel=nHel(1),nHel(2)
   do iHel=1,NumHelicities
!       useQP=0
!       pole_skipped=0
      QPredo=-1
!      print *, 'Helicity', iHel
!      print *, 'Born...'
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
!      print *, Helicities(iHel,1:5)
      call SetPolarizations()
      if( ZDecays.gt.0 ) then
          if( ExtParticle(5)%Helicity.eq.0 ) cycle!   this can be done more elegantly in mod_process
          call ZDecay(ExtParticle(5),DK_LO,MomExt(1:4,12:13))
!          print *, 'eZ', ExtParticle(5)%Pol(1:4)
!          pause
      endif

      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo



       LO_Res_Pol = (0d0,0d0)   ! MARKUS: has to be disabled for now because NumBornAmps counts more than just 1,2.
       do jPrimAmp=1,2
          do iPrimAmp=1,2
           LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
       enddo
       enddo
       LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol

! ----------------- bosonic loops --------------------

!      print *, 'virt...'
      ListPrimAmps=(/1,2,3,5,7,10,13,16,19,20,21,24,27,28/)
      do iPrimAmp=1,12
          call SetKirill(PrimAmps(iPrimAmp))
          call PentCut_new(PrimAmps(:),iPrimAmp)
       enddo

      do iPrimAmp=1,12
          call SetKirill(PrimAmps(iPrimAmp))
          call QuadCut_new(PrimAmps(:),iPrimAmp)
       enddo

      do iPrimAmp=1,12
          call SetKirill(PrimAmps(iPrimAmp))
          call TripCut_new(PrimAmps(:),iPrimAmp)
       enddo

      do iPrimAmp=1,12
          call SetKirill(PrimAmps(iPrimAmp))
          call DoubCut_new(PrimAmps(:),iPrimAmp)
       enddo

      do iPrimAmp=1,12
          call SetKirill(PrimAmps(iPrimAmp))
          call SingCut_new(PrimAmps(:),iPrimAmp)
       enddo

      do iPrimAmp=1,12
         call SetKirill(PrimAmps(iPrimAmp))
         call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
       enddo


! now combine into real prims:
       PrimAmps(3)%Result=PrimAmps(3)%Result+PrimAmps(4)%Result
       PrimAmps(5)%Result=PrimAmps(5)%Result+PrimAmps(6)%Result
       PrimAmps(7)%Result=PrimAmps(7)%Result+PrimAmps(8)%Result + PrimAmps(9)%Result
       PrimAmps(10)%Result=PrimAmps(10)%Result+PrimAmps(11)%Result + PrimAmps(12)%Result


! check on poles -- bosonic loops only
       do iPrimAmp=1,6
          APrimAmp=ListPrimAmps(iPrimAmp)
          call RenormalizeUV(PrimAmps(APrimAmp),BornAmps(APrimAmp),MuRen**2)
          PrimAmps(APrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(APrimAmp)%Result(-2:1)
          call OneLoopDiv(PrimAmps(APrimAmp),MuRen**2,3,rdiv(2),rdiv(1))

!          print *, 'CHANGING THE ANALYTIC POLES'
          rdiv(1)=rdiv(1)+1.5d0
          
!          print *, 'APrimAmp',APrimAmp,PrimAmps(APrimAmp)%Result(-2:1)
!          print *, 'BornAmp', APrimAmp,BornAmps(APrimAmp)%Result
!          print *, 'fin', PrimAmps(APrimAmp)%Result(0)+PrimAmps(APrimAmp)%Result(1)
!          print *, 'fin/Born', (PrimAmps(APrimAmp)%Result(0)+PrimAmps(APrimAmp)%Result(1))/BornAmps(APrimAmp)%Result
!          print *, 'DP',APrimAmp,PrimAmps(APrimAmp)%Result(-2)/BornAmps(APrimAmp)%Result
!          print *, 'SP',APrimAmp,PrimAmps(APrimAmp)%Result(-1)/BornAmps(APrimAmp)%Result
!          print *, 'CC',APrimAmp,PrimAmps(APrimAmp)%Result(0)/BornAmps(APrimAmp)%Result
!          print *, 'R',APrimAmp,PrimAmps(APrimAmp)%Result(1)/BornAmps(APrimAmp)%Result
!          print *, rdiv(1:2)
!          pause
!
          AccPoles = CheckPoles(PrimAmps(APrimAmp),BornAmps(APrimAmp),rdiv(1:2))
 !         print *, 'acc',AccPoles

! QP 
          if ( AccPoles .gt. DPtol ) then
             useQP=useQP+1
!              print *, 'QP'
!             QPredo(iPrimAmp,1)=APrimAmp
!             QPredo(iPrimAmp,2:PrimAmps(APrimAmp)%NumSisters+1)=PrimAmps(APrimAmp)%Sisters(1:PrimAmps(APrimAmp)%NumSisters)
             if (PrimAmps(APrimAmp)%NumSisters .gt. 0) then
                lastSister=PrimAmps(APrimAmp)%Sisters(PrimAmps(APrimAmp)%NumSisters)
             else
                lastSister=APrimAmp
             endif

             do jPrimAmp=APrimAmp,lastSister
                call SetKirill(PrimAmps(jPrimAmp))
                call PentCut_128_new(PrimAmps(:),jPrimAmp)
             enddo

             do jPrimAmp=APrimAmp,lastSister
                call SetKirill(PrimAmps(jPrimAmp))
                call QuadCut_128_new(PrimAmps(:),jPrimAmp)
             enddo

             do jPrimAmp=APrimAmp,lastSister
                call SetKirill(PrimAmps(jPrimAmp))
                call TripCut_128_new(PrimAmps(:),jPrimAmp)
             enddo

             do jPrimAmp=APrimAmp,lastSister
                call SetKirill(PrimAmps(jPrimAmp))
                call DoubCut_128_new(PrimAmps(:),jPrimAmp)
             enddo

             do jPrimAmp=APrimAmp,lastSister
                call SetKirill(PrimAmps(jPrimAmp))
                call SingCut_128_new(PrimAmps(:),jPrimAmp)
             enddo

             do jPrimAmp=APrimAmp,lastSister
                call SetKirill(PrimAmps(jPrimAmp))
                call EvalMasterIntegrals(PrimAmps(jPrimAmp),MuRen**2)
             enddo

!  sum up sisters to one prim
             do jPrimAmp=APrimAmp+1,lastSister
                PrimAmps(APrimAmp)%Result=PrimAmps(APrimAmp)%Result+PrimAmps(jPrimAmp)%Result
             enddo

             call RenormalizeUV(PrimAmps(APrimAmp),BornAmps(APrimAmp),MuRen**2)
             PrimAmps(APrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(APrimAmp)%Result(-2:1)

             AccPoles = CheckPoles(PrimAmps(APrimAmp),BornAmps(APrimAmp),rdiv(1:2))

!             print *, 'accpoles after QP', accpoles 
!          pause          
             if ( AccPoles .gt. QPtol) then
!                 print *, 'QP fails: ', AccPoles
                PrimAmps(APrimAmp)%Result=0d0
                pole_skipped=pole_skipped+1
                SkipCounter = SkipCounter + 1
                
                RETURN ! reject the whole event instead of just this primamp
                
             endif
          endif
!
       enddo



!------------ fermionic loops --------------




      if (ZQcoupl .eq. 1)  then
!         print *, 'left-handed loops'
         couplZQQ_left_dyn =one
         couplZQQ_right_dyn=zero
      elseif (ZQcoupl .eq. 2) then
!         print *, 'up in loop'
         couplZQQ_left_dyn =couplZUU_left
         couplZQQ_right_dyn=couplZUU_right
      elseif (ZQcoupl .eq. 3) then
!         print *, 'vector loops'
         couplZQQ_left_dyn=one/two
         couplZQQ_right_dyn=one/two
      endif


      do iPrimAmp=13,28
          call SetKirill(PrimAmps(iPrimAmp))
          call PentCut_new(PrimAmps(:),iPrimAmp)
       enddo

      do iPrimAmp=13,28
          call SetKirill(PrimAmps(iPrimAmp))
          call QuadCut_new(PrimAmps(:),iPrimAmp)
       enddo

      do iPrimAmp=13,28
          call SetKirill(PrimAmps(iPrimAmp))
          call TripCut_new(PrimAmps(:),iPrimAmp)
       enddo

      do iPrimAmp=13,28
          call SetKirill(PrimAmps(iPrimAmp))
          call DoubCut_new(PrimAmps(:),iPrimAmp)
       enddo

      do iPrimAmp=13,28
          call SetKirill(PrimAmps(iPrimAmp))
          call SingCut_new(PrimAmps(:),iPrimAmp)
       enddo

      do iPrimAmp=13,28
         call SetKirill(PrimAmps(iPrimAmp))
         call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
       enddo


! the fermion loops are combined into gauge invariant prims first
       PrimAmps(13)%Result = PrimAmps(13)%Result +PrimAmps(14)%Result +PrimAmps(15)%Result 
       PrimAmps(16)%Result = PrimAmps(16)%Result +PrimAmps(17)%Result +PrimAmps(18)%Result 
       PrimAmps(21)%Result = PrimAmps(21)%Result +PrimAmps(22)%Result +PrimAmps(23)%Result 
       PrimAmps(24)%Result = PrimAmps(24)%Result +PrimAmps(25)%Result +PrimAmps(26)%Result 
       
! check on poles
       do iPrimAmp=7,14
          APrimAmp=ListPrimAmps(iPrimAmp)
          call RenormalizeUV(PrimAmps(APrimAmp),BornAmps(APrimAmp),MuRen**2)
          PrimAmps(APrimAmp)%Result(-2:1) = -(0d0,1d0) * PrimAmps(APrimAmp)%Result(-2:1)
          call OneLoopDiv(PrimAmps(APrimAmp),MuRen**2,3,rdiv(2),rdiv(1))

!          print *, 'CHANGING THE ANALYTIC POLES'
          rdiv(1)=rdiv(1)+1.0d0/3.0d0
          
!          print *, 'APrimAmp',APrimAmp,PrimAmps(APrimAmp)%Result(-2:1)
!          print *, 'BornAmp', BornAmps(APrimAmp)%Result/Q_top/sqrt(2d0)
!          print *, 'fin', APrimAmp,abs( PrimAmps(APrimAmp)%Result(0)+PrimAmps(APrimAmp)%Result(1))
!          print *, 'fin/Born', (PrimAmps(APrimAmp)%Result(0)+PrimAmps(APrimAmp)%Result(1))/BornAmps(APrimAmp)%Result
!
!          print *, 'DP',APrimAmp,PrimAmps(APrimAmp)%Result(-2)/BornAmps(APrimAmp)%Result
!          print *, 'SP',APrimAmp,PrimAmps(APrimAmp)%Result(-1)/BornAmps(APrimAmp)%Result
!          print *, 'CC',APrimAmp,PrimAmps(APrimAmp)%Result(0)/BornAmps(APrimAmp)%Result
!          print *, 'R',APrimAmp,PrimAmps(APrimAmp)%Result(1)/BornAmps(APrimAmp)%Result
!          print *, 'DP',APrimAmp,PrimAmps(APrimAmp)%Result(-2)/sqrt(2d0)
!          print *, 'SP',APrimAmp,PrimAmps(APrimAmp)%Result(-1)/sqrt(2d0)
!          print *, 'CC',APrimAmp,PrimAmps(APrimAmp)%Result(0)/sqrt(2d0)
!          print *, 'R',APrimAmp,PrimAmps(APrimAmp)%Result(1)/sqrt(2d0)
!          print *, rdiv(1:2)
!          pause
!
          AccPoles = CheckPoles(PrimAmps(APrimAmp),BornAmps(APrimAmp),rdiv(1:2))
!          print *, 'acc',AccPoles

! QP 
          if ( AccPoles .gt. DPtol ) then

             useQP=useQP+1

             if (PrimAmps(APrimAmp)%NumSisters .gt. 0) then
                lastSister=PrimAmps(APrimAmp)%Sisters(PrimAmps(APrimAmp)%NumSisters)
             else
                lastSister=APrimAmp
             endif

             do jPrimAmp=APrimAmp,lastSister
                call SetKirill(PrimAmps(jPrimAmp))
                call PentCut_128_new(PrimAmps(:),jPrimAmp)
             enddo

             do jPrimAmp=APrimAmp,lastSister
                call SetKirill(PrimAmps(jPrimAmp))
                call QuadCut_128_new(PrimAmps(:),jPrimAmp)
             enddo

             do jPrimAmp=APrimAmp,lastSister
                call SetKirill(PrimAmps(jPrimAmp))
                call TripCut_128_new(PrimAmps(:),jPrimAmp)
             enddo

             do jPrimAmp=APrimAmp,lastSister
                call SetKirill(PrimAmps(jPrimAmp))
                call DoubCut_128_new(PrimAmps(:),jPrimAmp)
             enddo

             do jPrimAmp=APrimAmp,lastSister
                call SetKirill(PrimAmps(jPrimAmp))
                call SingCut_128_new(PrimAmps(:),jPrimAmp)
             enddo

             do jPrimAmp=APrimAmp,lastSister
                call SetKirill(PrimAmps(jPrimAmp))
                call EvalMasterIntegrals(PrimAmps(jPrimAmp),MuRen**2)
             enddo

!  sum up sisters to one prim
             do jPrimAmp=APrimAmp+1,lastSister
                PrimAmps(APrimAmp)%Result=PrimAmps(APrimAmp)%Result+PrimAmps(jPrimAmp)%Result
             enddo

             call RenormalizeUV(PrimAmps(APrimAmp),BornAmps(APrimAmp),MuRen**2)
             PrimAmps(APrimAmp)%Result(-2:1) = -(0d0,1d0) * PrimAmps(APrimAmp)%Result(-2:1)

             AccPoles = CheckPoles(PrimAmps(APrimAmp),BornAmps(APrimAmp),rdiv(1:2))

!             print *, 'accpoles after QP', accpoles 
!          pause
             if ( AccPoles .gt. QPtol) then
!                 print *, 'QP fails: ', AccPoles
                PrimAmps(APrimAmp)%Result=0d0
                pole_skipped=pole_skipped+1
                
                RETURN ! reject the whole event instead of just this primamp                
                
             endif
          endif
!
        enddo
       

! now combine into fermloops partials, with appropriate Z-coupls


       if (ZQcoupl .eq. 1) then 
          FermionLoopPartAmp(1,-2:1)=(2d0*couplZUU_left+3d0*couplZDD_left)*PrimAmps(PrimAmp2_12534)%Result + &
               & nf_light * PrimAmps(PrimAmp2_15234)%Result + &
               & PrimAmps(PrimAmp2m_12534)%Result + &
               & PrimAmps(PrimAmp2m_15234)%Result

          FermionLoopPartAmp(2,-2:1)=(2d0*couplZUU_left+3d0*couplZDD_left)*PrimAmps(PrimAmp2_12543)%Result + &
               & nf_light * PrimAmps(PrimAmp2_15243)%Result + &
               & PrimAmps(PrimAmp2m_12543)%Result + &
               & PrimAmps(PrimAmp2m_15243)%Result
          
       elseif (ZQcoupl .eq. 2) then

          FermionLoopPartAmp(1,-2:1)=2d0*PrimAmps(PrimAmp2_12534)%Result &
               & +nf_light * PrimAmps(PrimAmp2_15234)%Result + &
               & PrimAmps(PrimAmp2m_12534)%Result + &
               & PrimAmps(PrimAmp2m_15234)%Result
 
          FermionLoopPartAmp(2,-2:1)=2d0*PrimAmps(PrimAmp2_12543)%Result &
               & +nf_light * PrimAmps(PrimAmp2_15243)%Result + &
               & PrimAmps(PrimAmp2m_12543)%Result+ &
               & PrimAmps(PrimAmp2m_15243)%Result
          
       elseif (ZQcoupl .eq. 3) then
          FermionLoopPartAmp(1,-2:1)=&
               & (2d0*couplZUU_left+3d0*couplZDD_left+2d0*couplZUU_right+3d0*couplZDD_right)*PrimAmps(PrimAmp2_12534)%Result &
               & + nf_light * PrimAmps(PrimAmp2_15234)%Result  &
               & + PrimAmps(PrimAmp2m_12534)%Result &
               & + PrimAmps(PrimAmp2m_15234)%Result

          FermionLoopPartAmp(2,-2:1)=&
               & (2d0*couplZUU_left+3d0*couplZDD_left+2d0*couplZUU_right+3d0*couplZDD_right)*PrimAmps(PrimAmp2_12543)%Result  &
               & + nf_light * PrimAmps(PrimAmp2_15243)%Result &
               & + PrimAmps(PrimAmp2m_12543)%Result  &
               & + PrimAmps(PrimAmp2m_15243)%Result


       endif

!       print *, 'for comparison:'
!       print *, 'fin', (FermionLoopPartAmp(1,0)+FermionLoopPartAmp(1,1))/Q_top
!       print *, 'fin/Born', (FermionLoopPartAmp(1,0)+FermionLoopPartAmp(1,1))/BornAmps(1)%Result
!       pause
!       
!  -- other ferm loops

       if (ZQcoupl .eq. 1) then
!          print *, 'right-handed loops'
          couplZQQ_left_dyn =zero
          couplZQQ_right_dyn=one
       elseif (ZQcoupl .eq. 2) then
!          print *, 'down in loop'
          couplZQQ_left_dyn =couplZDD_left
          couplZQQ_right_dyn=couplZDD_right
       elseif (ZQcoupl .eq. 3) then
!          print *, 'axial-vector loops'
          couplZQQ_left_dyn=-one/two
          couplZQQ_right_dyn=one/two
       endif

      do iPrimAmp=13,18
          call SetKirill(PrimAmps(iPrimAmp))
          call PentCut_new(PrimAmps(:),iPrimAmp)
       enddo

      do iPrimAmp=13,18
          call SetKirill(PrimAmps(iPrimAmp))
          call QuadCut_new(PrimAmps(:),iPrimAmp)
       enddo

      do iPrimAmp=13,18
          call SetKirill(PrimAmps(iPrimAmp))
          call TripCut_new(PrimAmps(:),iPrimAmp)
       enddo

      do iPrimAmp=13,18
          call SetKirill(PrimAmps(iPrimAmp))
          call DoubCut_new(PrimAmps(:),iPrimAmp)
       enddo

      do iPrimAmp=13,18
          call SetKirill(PrimAmps(iPrimAmp))
          call SingCut_new(PrimAmps(:),iPrimAmp)
       enddo

      do iPrimAmp=13,18
         call SetKirill(PrimAmps(iPrimAmp))
         call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
       enddo


! the fermion loops are combined into gauge invariant prims first
       PrimAmps(13)%Result = PrimAmps(13)%Result +PrimAmps(14)%Result +PrimAmps(15)%Result 
       PrimAmps(16)%Result = PrimAmps(16)%Result +PrimAmps(17)%Result +PrimAmps(18)%Result 
       
! check on poles
       do iPrimAmp=7,8
          APrimAmp=ListPrimAmps(iPrimAmp)
          call RenormalizeUV(PrimAmps(APrimAmp),BornAmps(APrimAmp),MuRen**2)
          PrimAmps(APrimAmp)%Result(-2:1) = -(0d0,1d0) * PrimAmps(APrimAmp)%Result(-2:1)
          call OneLoopDiv(PrimAmps(APrimAmp),MuRen**2,3,rdiv(2),rdiv(1))

!          print *, 'CHANGING THE ANALYTIC POLES'
          rdiv(1)=rdiv(1)+1.0d0/3.0d0
          
!          print *, 'APrimAmp',APrimAmp,PrimAmps(APrimAmp)%Result(-2:1)
!          print *, 'BornAmp', APrimAmp,BornAmps(APrimAmp)%Result
!          print *, 'fin', PrimAmps(APrimAmp)%Result(0)+PrimAmps(APrimAmp)%Result(1)
!          print *, 'fin/Born', (PrimAmps(APrimAmp)%Result(0)+PrimAmps(APrimAmp)%Result(1))/BornAmps(APrimAmp)%Result
!          print *, 'DP',APrimAmp,PrimAmps(APrimAmp)%Result(-2)/BornAmps(APrimAmp)%Result
!          print *, 'SP',APrimAmp,PrimAmps(APrimAmp)%Result(-1)/BornAmps(APrimAmp)%Result
!          print *, 'CC',APrimAmp,PrimAmps(APrimAmp)%Result(0)/BornAmps(APrimAmp)%Result
!          print *, 'R',APrimAmp,PrimAmps(APrimAmp)%Result(1)/BornAmps(APrimAmp)%Result
!          print *, rdiv(1:2)
!          pause
!
          AccPoles = CheckPoles(PrimAmps(APrimAmp),BornAmps(APrimAmp),rdiv(1:2))
!          print *, 'acc',AccPoles

! QP 
          if ( AccPoles .gt. DPtol ) then
             useQP=useQP+1

             if (PrimAmps(APrimAmp)%NumSisters .gt. 0) then
                lastSister=PrimAmps(APrimAmp)%Sisters(PrimAmps(APrimAmp)%NumSisters)
             else
                lastSister=APrimAmp
             endif

             do jPrimAmp=APrimAmp,lastSister
                call SetKirill(PrimAmps(jPrimAmp))
                call PentCut_128_new(PrimAmps(:),jPrimAmp)
             enddo

             do jPrimAmp=APrimAmp,lastSister
                call SetKirill(PrimAmps(jPrimAmp))
                call QuadCut_128_new(PrimAmps(:),jPrimAmp)
             enddo

             do jPrimAmp=APrimAmp,lastSister
                call SetKirill(PrimAmps(jPrimAmp))
                call TripCut_128_new(PrimAmps(:),jPrimAmp)
             enddo

             do jPrimAmp=APrimAmp,lastSister
                call SetKirill(PrimAmps(jPrimAmp))
                call DoubCut_128_new(PrimAmps(:),jPrimAmp)
             enddo

             do jPrimAmp=APrimAmp,lastSister
                call SetKirill(PrimAmps(jPrimAmp))
                call SingCut_128_new(PrimAmps(:),jPrimAmp)
             enddo

             do jPrimAmp=APrimAmp,lastSister
                call SetKirill(PrimAmps(jPrimAmp))
                call EvalMasterIntegrals(PrimAmps(jPrimAmp),MuRen**2)
             enddo

!  sum up sisters to one prim
             do jPrimAmp=APrimAmp+1,lastSister
                PrimAmps(APrimAmp)%Result=PrimAmps(APrimAmp)%Result+PrimAmps(jPrimAmp)%Result
             enddo

             call RenormalizeUV(PrimAmps(APrimAmp),BornAmps(APrimAmp),MuRen**2)
             PrimAmps(APrimAmp)%Result(-2:1) = -(0d0,1d0) * PrimAmps(APrimAmp)%Result(-2:1)

             AccPoles = CheckPoles(PrimAmps(APrimAmp),BornAmps(APrimAmp),rdiv(1:2))

!             print *, 'accpoles after QP', accpoles 
             if ( AccPoles .gt. QPtol) then
!                 print *, 'QP fails: ', AccPoles
                PrimAmps(APrimAmp)%Result=0d0
                pole_skipped=pole_skipped+1
             endif
!             pause
          endif
!          pause
       enddo
       

! now combine into fermloops partials, with appropriate Z-coupls

       if (ZQcoupl .eq. 1) then
          FermionLoopPartAmp(1,-2:1)=FermionLoopPartAmp(1,-2:1) + (2d0*couplZUU_right+3d0*couplZDD_right)*(PrimAmps(PrimAmp2_12534)%Result)

          FermionLoopPartAmp(2,-2:1)=FermionLoopPartAmp(2,-2:1) + (2d0*couplZUU_right+3d0*couplZDD_right)*(PrimAmps(PrimAmp2_12543)%Result)

       elseif (ZQcoupl .eq. 2) then
          FermionLoopPartAmp(1,-2:1)=FermionLoopPartAmp(1,-2:1) + 3d0*PrimAmps(PrimAmp2_12534)%Result

          FermionLoopPartAmp(2,-2:1)=FermionLoopPartAmp(2,-2:1) + 3d0*PrimAmps(PrimAmp2_12543)%Result

       elseif (ZQcoupl .eq. 3) then
! axial-vector coupl
!          print *, ' axial-vector' ,-couplZUU_left-couplZDD_left+couplZUU_right+couplZDD_right
          FermionLoopPartAmp(1,-2:1)=FermionLoopPartAmp(1,-2:1)  + &
               & (-2d0*couplZUU_left-3d0*couplZDD_left+2d0*couplZUU_right+3d0*couplZDD_right)*(PrimAmps(PrimAmp2_12534)%Result)

          FermionLoopPartAmp(2,-2:1)=FermionLoopPartAmp(2,-2:1) + &
               & (-2d0*couplZUU_left-3d0*couplZDD_left+2d0*couplZUU_right+3d0*couplZDD_right)*(PrimAmps(PrimAmp2_12543)%Result)
       endif

!       print *, 'full ferm loops', FermionLoopPartAmp(1,-2:1)/Q_top
!       print *, 'full ferm loop/born', FermionLoopPartAmp(1,-2:1)/BornAmps(1)%Result
!       pause


!print *,"This is for the check against ttb+photon. Remember to change couplings to Q_up"
!!run ! ./TOPAZ Collider=1 TopDK=0 ZDK=-2 Process=81 Correction=1 NLOParam=0 ObsSet=51
!!run ! ./TOPAZ Collider=1 TopDK=0 Process=20 Correction=1 NLOParam=0 ObsSet=21
!print *, "Z LO",iHel,cdabs( BornAmps(2)%Result )/( Q_up *dsqrt(2d0) )
!print *, "Z VI",iHel,cdabs( PrimAmps(2)%Result(-2) )/( Q_up *dsqrt(2d0) )
!print *, "Z VI",iHel,cdabs( PrimAmps(2)%Result(-1) )/( Q_up *dsqrt(2d0) )
!print *, "Z VI",iHel,cdabs( PrimAmps(2)%Result(+0) )/( Q_up *dsqrt(2d0) )
!print *, "Z VI",iHel,cdabs( PrimAmps(2)%Result(+1) ) /( Q_up *dsqrt(2d0) )
!pause


       BosonicPartAmp(1,-2:1) =   PrimAmps(PrimAmp1_15234)%Result(-2:1) &
            & - 1d0/Nc**2 *(  PrimAmps(PrimAmp1_15432)%Result(-2:1) )
       BosonicPartAmp(2,-2:1) =   PrimAmps(PrimAmp1_15243)%Result(-2:1) &
            & - 1d0/Nc**2 *(  PrimAmps(PrimAmp1_15342)%Result(-2:1) )

       BosonicPartAmp(3,-2:1) =   PrimAmps(PrimAmp1_15234)%Result(-2:1) &
            + PrimAmps(PrimAmp1_15243)%Result(-2:1) &
            + PrimAmps(PrimAmp1_13524)%Result(-2:1) &
            + PrimAmps(PrimAmp1_14523)%Result(-2:1) &
            + PrimAmps(PrimAmp1_15342)%Result(-2:1) &
            + PrimAmps(PrimAmp1_15432)%Result(-2:1) 

!! print *, ""
!! print *, "1-loop amplitudes:"
!! print *, "bosonic loops"
!! print *, "|Ta*Tb|^2*Re[M0*cong(M1)]:", PhotonCouplCorr*ColLO_ttbgg(1,1)*dreal( BornAmps(1)%Result * dconjg( BosonicPartAmp(1,-2:1) ))
!! print *, "|Tb*Ta|^2*Re[M0*cong(M1)]:", PhotonCouplCorr*ColLO_ttbgg(2,2)*dreal( BornAmps(2)%Result * dconjg( BosonicPartAmp(2,-2:1) ))
!! print *, "check",cdabs( BosonicPartAmp(3,-2:1) )
!! stop
!
      NLO_Res_Pol(-2:1) = (0d0,0d0)
! this should really be set elsewhere...
      Col1L_ttbggZ = 0d0
      Col1L_ttbggZ(1,1)= 4d0 * Cf**2 * Nc**2
      Col1L_ttbggZ(2,2)= Col1L_ttbggZ(1,1)
      Col1L_ttbggZ(1,2)= - 2d0 * Cf * Nc
      Col1L_ttbggZ(2,1)= Col1L_ttbggZ(1,2)
      Col1L_ttbggZ(1,3)= 2d0 * Cf * Nc
      Col1L_ttbggZ(2,3)= Col1L_ttbggZ(1,3)
      
      do jPrimAmp=1,3
         do iPrimAmp=1,2
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + Col1L_ttbggZ(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(BosonicPartAmp(jPrimAmp,-2:1)) )
      enddo
   enddo
      NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1) + NLO_Res_Pol(-2:1)

!       print *, 'bosonic loops:'
!       print *, 'color-summed DP',NLO_Res_Pol(-2)/LO_Res_Pol
!       print *, 'color-summed SP',NLO_Res_Pol(-1)/LO_Res_Pol
!       print *, 'color-summed boson FIN',(NLO_Res_Pol(0)+NLO_Res_Pol(1))/2d0/Q_top**2
!      PAUSE
!      print *, 'color-summed SP', NLO_Res_Pol(-1)/LO_Res_Pol
!      pause

      Col1Lf_ttbggZ = 0d0
      Col1Lf_ttbggZ(1,1) = 4d0 * Cf**2 * Nc - 2d0*Cf
      Col1Lf_ttbggZ(1,2) = -4d0*Cf
      Col1Lf_ttbggZ(2,2) = Col1Lf_ttbggZ(1,1)
      Col1Lf_ttbggZ(2,1) = Col1Lf_ttbggZ(1,2)

      NLO_Res_Pol(-2:1) = (0d0,0d0)
      do jPrimAmp=1,2
         do iPrimAmp=1,2
          NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + Col1Lf_ttbggZ(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result*dconjg(FermionLoopPartAmp(jPrimAmp,-2:1)) )
! NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + ColLO_ttbgg(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result*dconjg(FermionLoopPartAmp(jPrimAmp,-2:1)) )
      enddo
      enddo
!       print *, 'fermionic loops:'
!       print *, 'color-summed DP',NLO_Res_Pol(-2)/LO_Res_Pol
!       print *, 'color-summed SP',NLO_Res_Pol(-1)/LO_Res_Pol
!       print *, 'color-summed ferm loop FIN',(NLO_Res_Pol(0)+NLO_Res_Pol(1))/2d0/Q_top**2
!       pause
      NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) + NLO_Res_Pol(-2:1)

!      print *, useQP, pole_skipped
!      pause
    enddo! helicity loop
!    print *, 'boson',NLO_Res_UnPol(0:1)
!    print *, 'ferm',NLO_Res_UnPol_Ferm(0:1)
!    print *, 'UV stuff:'
!    stop
    
ENDIF
!print *, 'bos DP', NLO_Res_UnPol(-2)
!print *, 'bos SP', NLO_Res_UnPol(-1)
!print *, 'bos CC', NLO_Res_UnPol(0)
!print *, 'bos rat', NLO_Res_UnPol(1)
!print *, 'ferm DP',  NLO_Res_UnPol_Ferm(-2)
!print *, 'ferm SP',  NLO_Res_UnPol_Ferm(-1)
!print *, 'ferm CC',  NLO_Res_UnPol_Ferm(0)
!print *, 'ferm rat', NLO_Res_UnPol_Ferm(1)
!print *, 'No. times QP used : ', useQP
!print *, 'No. times QP fails: ', pole_skipped
!stop

IF( Correction.EQ.0 ) THEN
!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * alpha4Pi * WidthExpansion
   EvalCS_1L_ttbggZ = LO_Res_Unpol * PreFac
!   print *, prefac
!   print *, 'ds', EvalCS_1L_ttbggZ

ELSEIF( Correction.EQ.1 ) THEN
!  overall normalization: (4*Pi)^eps/Gamma(1-eps)
!  CT contributions
                         ! beta        !top WFRC
   NLO_Res_UnPol(-1) = NLO_Res_UnPol(-1) + (-11d0/3d0*3d0 - 3d0*4d0/3d0 )*LO_Res_Unpol
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-3d0*4d0/3d0)*2d0*dlog(MuRen/m_top)*LO_Res_Unpol  ! finite log(mu) contrib. from  top WFRC
   NLO_Res_UnPol_Ferm(-1) = NLO_Res_UnPol_Ferm(-1) - (-2d0/3d0*Nf_light)*LO_Res_Unpol
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-5d0/2d0*8d0/3d0 )*LO_Res_Unpol   ! finite contribution from top WFRC's
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + LO_Res_Unpol  ! shift alpha_s^DR --> alpha_s^MSbar
!   print *, 'boson',NLO_Res_UnPol(0:1)
!   print *, 'ferm',NLO_Res_UnPol_Ferm(0:1)
!   print *,'coupl'
!  normalization
   LO_Res_Unpol = LO_Res_Unpol                         * ISFac * (alpha_s4Pi*RunFactor)**2                            * alpha4Pi
   NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1)           * ISFac * (alpha_s4Pi*RunFactor)**2 * alpha_sOver2Pi*RunFactor * alpha4Pi
   NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) * ISFac * (alpha_s4Pi*RunFactor)**2 * alpha_sOver2Pi*RunFactor * alpha4Pi
!   print *, 'boson',NLO_Res_UnPol(0:1)
!   print *, 'ferm',NLO_Res_UnPol_Ferm(0:1)
!    print *, 'PreFac', PreFac
   EvalCS_1L_ttbggZ = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac
!    print *, 'ds', EvalCS_1L_ttbggZ


ELSEIF( CORRECTION.EQ.3 ) THEN

! NLO_Res_UnPol      = NLO_Res_UnPol      * PreFac
! NLO_Res_UnPol_Ferm = NLO_Res_UnPol_Ferm * PreFac
! LO_Res_Unpol       = LO_Res_Unpol       * PreFac

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   IF( TOPDECAYS.GE.1 ) THEN
       xE = yRnd(18+HelSampling)
   ELSEIF( TOPDECAYS.EQ.0 ) THEN
       xE = yRnd(8+HelSampling)
   ENDIF

! xe=0.3d0

   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

  call EvalIntDipoles_GGTTBGZ((/MomExt(1:4,4),MomExt(1:4,3),MomExt(1:4,5),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,6:13),xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
   EvalCS_1L_ttbggZ = HOp(1)    * pdf(0,1)  * pdf(0,2)   &
                    + HOp(2)/xE * pdf_z(0,1)* pdf(0,2)   &
                    + HOp(3)/xE * pdf(0,1)  * pdf_z(0,2)

! print *, "1L check",NLO_Res_UnPol(-2)                           /(alpha_sOver2Pi*RunFactor)/(LO_Res_Unpol)
! print *, "1L check",( NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1) )/(alpha_sOver2Pi*RunFactor)/(LO_Res_Unpol)
! print *, "ID check", EvalCS_1L_ttbggZ/(alpha_sOver2Pi*RunFactor)/(LO_Res_Unpol)
! pause



ENDIF




!      careful, alphas=0.13 only for NLOParam=0 PDFSet=2
!      MADGRAPH CHECK: gg->ttbZ, mt=172, alpha_s=0.13 mZ=91.19
! if (ZDecays .eq. 0) then
!       MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!       MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!       MG_MOM(0:3,3) = MomExt(1:4,5)*100d0
!       MG_MOM(0:3,4) = MomExt(1:4,4)*100d0
!       MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!       call coupsm(0)
!       call SGG_TTBZ(MG_MOM,MadGraph_tree)
!       print *, ""
!       print *, alpha_s*RunFactor,m_top,m_z
!       print *, "My tree:         ", LO_Res_Unpol/(100d0)**2
!       print *, "MadGraph hel.amp:", MadGraph_tree
!       print *, "MG/ME ratio: ", MadGraph_tree/(dble(LO_Res_Unpol)/(100d0)**2)
!       pause
!    else
!       print *, "My tree:         ", LO_Res_Unpol/(100d0)**4
!       print *, "MG/ME ratio: ", 0.872680470745814d-13/(LO_Res_Unpol/(100d0)**4)
!    endif
!    stop
       



   if( IsNan(EvalCS_1L_ttbggZ) ) then
        print *, "NAN:",EvalCS_1L_ttbggZ
        print *, yRnd(:)
        print *, NLO_Res_UnPol(0),NLO_Res_UnPol(1),NLO_Res_UnPol_Ferm(0),NLO_Res_UnPol_Ferm(1)
        print *, PSWgt , VgsWgt , PDFFac, sHatJacobi
        print *, eta1,eta2,MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,:)
        print *, "SKIP EVENT!!!!!"
        EvalCS_1L_ttbggZ = 0d0
!         pause
        return
   endif


   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_ttbggZ)
   enddo
   EvalCounter = EvalCounter + 1 


   EvalCS_1L_ttbggZ = EvalCS_1L_ttbggZ/VgsWgt

RETURN
END FUNCTION







FUNCTION EvalCS_1L_ttbqqbZ(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModZDecay
!use ModIntDipoles_QQBTTBGP
!use ModIntDipoles_QGTTBQP
!use ModIntDipoles_QBGTTBQBP
implicit none
real(8) ::  EvalCS_1L_ttbqqbZ,yRnd(1:VegasMxDim),VgsWgt,xE
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1)
complex(8) :: BosonicPartAmp(1:2,-2:1),FermionPartAmp(1:2,-2:1),mydummy(1:2),LOPartAmp(1:2)
integer :: iHel,iPrimAmp,jPrimAmp,tmphel,APrimAmp,ListPrimAmps(4)
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,PSWgt4,ISFac,AccPoles,HOp(1:2,1:3),pdf_z(-6:6,1:2)
real(8) :: MomExt(1:4,1:14),MomZl(1:4), MomZa(1:4),pZsq,MZ_Inv,ZDKMatel,Msq_T_BWENU
complex(8) :: propZ,ZPolVec(1:4),BarSpi(1:4),Spi(1:4)
logical :: applyPSCut
real(8) :: couplZUU,couplZDD,couplZLL,couplGUU,couplGDD,couplGLL,couplZTT,couplGTT
real(8) :: MG_MOM(0:3,1:NumExtParticles)
real(8) :: MadGraph_tree
real(8),parameter :: Nc=3d0
real(8) :: tau,eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac_a(1:2),PDFFac_b(1:2),PDFFac(1:2),pdf(-6:6,1:2)
integer :: NHisto,NBin(1:NumMaxHisto),npdf,ParityFlip=1,PhotonCouplCorr=2d0,nHel(1:2),NRndHel
integer,parameter :: up=1,dn=2
include 'misc/global_import'
include "vegas_common.f"


  EvalCS_1L_ttbqqbZ = 0d0
  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top+M_Z ) then
      EvalCS_1L_ttbqqbZ = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

  if( ZDecays.le.10 ) then  ! decaying on-shell Z
      MZ_Inv = m_Z
   elseif( ZDecays.gt.10 ) then  ! decaying off-shell Z
      call Error("need to implement phase space for off-shell Z's")
   endif
  call EvalPhaseSpace_2to3M(EHat,MZ_Inv,yRnd(3:7),MomExt(1:4,1:5),PSWgt)
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

IF( ZDECAYS.NE.0 ) THEN
   call EvalPhasespace_ZDecay(MZ_Inv,MomExt(1:4,3),yRnd(16:17),MomExt(1:4,12:13),PSWgt4)
   PSWgt = PSWgt * PSWgt4
ENDIF

IF( ZDECAYS.EQ.-1 ) THEN
     ExtParticle(5)%Helicity=+1
     call ZDecay(ExtParticle(5),DK_LO,MomExt(1:4,12:13))
     call pol_massSR(ExtParticle(5)%Mom(1:4),M_Z,+1,ZPolVec(1:4))
     ZDKMatel = cdabs((ExtParticle(5)%Pol(1:4).dot.ZPolVec(1:4)))**2
     call pol_massSR(ExtParticle(5)%Mom(1:4),M_Z,-1,ZPolVec(1:4))
     ZDKMatel = ZDKMatel + cdabs((ExtParticle(5)%Pol(1:4).dot.ZPolVec(1:4)))**2
     call pol_massSR(ExtParticle(5)%Mom(1:4),M_Z,0,ZPolVec(1:4))
     ZDKMatel = ZDKMatel + cdabs((ExtParticle(5)%Pol(1:4).dot.ZPolVec(1:4)))**2

     ExtParticle(5)%Helicity=-1
     call ZDecay(ExtParticle(5),DK_LO,MomExt(1:4,12:13))
     call pol_massSR(ExtParticle(5)%Mom(1:4),M_Z,+1,ZPolVec(1:4))
     ZDKMatel = ZDKMatel + cdabs((ExtParticle(5)%Pol(1:4).dot.ZPolVec(1:4)))**2
     call pol_massSR(ExtParticle(5)%Mom(1:4),M_Z,-1,ZPolVec(1:4))
     ZDKMatel = ZDKMatel + cdabs((ExtParticle(5)%Pol(1:4).dot.ZPolVec(1:4)))**2
     call pol_massSR(ExtParticle(5)%Mom(1:4),M_Z,0,ZPolVec(1:4))
     ZDKMatel = ZDKMatel + cdabs((ExtParticle(5)%Pol(1:4).dot.ZPolVec(1:4)))**2

     PSWgt = PSWgt * ZDKMatel/3d0

      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
      Spi(1:4) = ExtParticle(1)%Pol(1:4)
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,-1,BarSpi(1:4))
      Msq_T_BWENU = cdabs(psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top))**2/2d0
      call vBarSpi(ExtParticle(1)%Mom(1:4),m_Top,+1,BarSpi(1:4))
      Msq_T_BWENU = Msq_T_BWENU + cdabs(psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top))**2/2d0

      PSWgt = PSWgt * Msq_T_BWENU


      call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))
      BarSpi(1:4) = ExtParticle(2)%Pol(1:4)
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,-1,Spi(1:4))
      Msq_T_BWENU = cdabs(psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top))**2/2d0
      call uSpi(ExtParticle(2)%Mom(1:4),m_Top,+1,Spi(1:4))
      Msq_T_BWENU = Msq_T_BWENU + cdabs(psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top))**2/2d0

      PSWgt = PSWgt * Msq_T_BWENU
ENDIF



   call Kinematics_TTBARZ(0,MomExt(1:4,1:14),(/4,5,3,1,2,0,6,7,8,9,10,11,12,13/),applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_1L_ttbqqbZ = 0d0
      return
   endif

! we still need this for the massless qqZ couplings
   pZsq=MomExt(1,3)*MomExt(1,3)-MomExt(2,3)*MomExt(2,3)-MomExt(3,3)*MomExt(3,3)-MomExt(4,3)*MomExt(4,3)

   if ( ZDecays .lt. 10) then
      propZ = (1d0,0d0)/dsqrt(2d0*Ga_Zexp*m_Z)
   elseif (ZDecays .gt. 10) then 
      propZ=cone/(pZsq-m_Z**2+ci*Ga_ZExp*m_Z)
   endif
   propZ=pZsq*propZ

   call setPDFs(eta1,eta2,MuFac,pdf)

IF( PROCESS.EQ.72 ) THEN
      PDFFac_a(up) = pdf(Up_,1)*pdf(AUp_,2) + pdf(Chm_,1)*pdf(AChm_,2)
      PDFFac_a(dn) = pdf(Dn_,1)*pdf(ADn_,2) + pdf(Str_,1)*pdf(AStr_,2) + pdf(Bot_,1)*pdf(ABot_,2)
      PDFFac_b(up) = pdf(Up_,2)*pdf(AUp_,1) + pdf(Chm_,2)*pdf(AChm_,1)
      PDFFac_b(dn) = pdf(Dn_,2)*pdf(ADn_,1) + pdf(Str_,2)*pdf(AStr_,1) + pdf(Bot_,2)*pdf(ABot_,1)
   ENDIF

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   RunFactor = RunAlphaS(NLOParam,MuRen)
   nHel(1:2) = getHelicity(yrnd(NRndHel))

!   PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))

   LO_Res_Unpol             = (0d0,0d0)
   NLO_Res_Unpol(-2:1)      = (0d0,0d0)
   NLO_Res_Unpol_Ferm(-2:1) = (0d0,0d0)
!------------ LO --------------
IF( CORRECTION.EQ.0 ) THEN
   do npdf=1,2
      if (npdf .eq. 1) then
         PDFFac(1:2)=PDFFac_a(1:2)
      elseif (npdf .eq. 2) then
         PDFFac(1:2)=PDFFac_b(1:2)
         call swapMom(MomExt(1:4,1),MomExt(1:4,2))
     endif
         
    ISFac = MomCrossing(MomExt)
    call SetPropagators()

!    do iHel=nHel(1),nHel(2)
    do iHel=1,NumHelicities
       call HelCrossing(Helicities(iHel,1:NumExtParticles))

        call SetPolarizations()
        
        if( ZDecays.ne.0 ) then
           if( ExtParticle(5)%Helicity.eq.0 ) cycle!   this can be more elegantly done in mod_process
           call ZDecay(ExtParticle(5),DK_LO,MomExt(1:4,12:13))
           call ZGamLCoupl(1,Helicities(iHel,5),couplZLL,couplGLL)  ! charged lept
        endif


       if (npdf .eq. 2) then
! change helicities of the massless quarks for the couplings to Z          
          Helicities(iHel,3)=-Helicities(iHel,3)
          Helicities(iHel,4)=-Helicities(iHel,4)
       endif

           call ZGamQcoupl(Up_,Helicities(iHel,3),couplZUU,couplGUU)
! one could here add a routine for the anomalous ttZ coupling, or modify this one. For now, use this for ttZ coupl
           call ZGamQcoupl(Dn_,Helicities(iHel,3),couplZDD,couplGDD)

      ! RR -- not sure about this -- recheck
      couplZQQ_left_dyn=one
      couplZQQ_right_dyn=one                  
      do iPrimAmp=1,NumBornAmps
         call EvalTree(BornAmps(iPrimAmp))
      enddo

! RR : this checks the f_4fV current
        print *,BornAmps(1)%Result
        print *, BornAmps(2)%Result
        print *, BornAmps(1)%Result+ BornAmps(2)%Result
        print *, BornAmps(3)%Result
        pause

        LO_Res_Pol = (0d0,0d0)

! Amp 1 has the Z coupl in the currents
! Amp 2 has the Z/gamma attached to the massless line, which could be up or down
!        BornAmps(2)%Result=BornAmps(2)%Result*( couplZUU*propZ*couplZLL + couplGUU*couplGLL)
!        BornAmps(2)%Result=BornAmps(2)%Result*( couplZDD*propZ*couplZLL + couplGDD*couplGLL)

! these are actually just different flavor structures, same color structure, so
! add them
!        BornAmps(1)%Result=BornAmps(1)%Result+BornAmps(2)%Result

        if (Zdecays .eq. 0) then
           LOPartAmp(up)=BornAmps(1)%Result+BornAmps(2)%Result*( couplZUU )
           LOPartAmp(dn)=BornAmps(1)%Result+BornAmps(2)%Result*( couplZDD )
        elseif( Zdecays.lt.10 ) then
           LOPartAmp(up)=BornAmps(1)%Result+BornAmps(2)%Result*( couplZUU*propZ*couplZLL )
           LOPartAmp(dn)=BornAmps(1)%Result+BornAmps(2)%Result*( couplZDD*propZ*couplZLL )
        elseif( Zdecays.gt.10 ) then
           LOPartAmp(up)=BornAmps(1)%Result+BornAmps(2)%Result*( couplZUU*propZ*couplZLL + couplGUU*couplGLL )
           LOPartAmp(dn)=BornAmps(1)%Result+BornAmps(2)%Result*( couplZDD*propZ*couplZLL + couplGDD*couplGLL )
        endif

        LO_Res_Pol =  ColLO_ttbqqb(1,1) * ( LOPartAmp(up)*dconjg(LOPartAmp(up))*PDFFac(up) + LOPartAmp(dn)*dconjg(LOPartAmp(dn))*PDFFac(dn))

      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
   enddo!helicity loop
 enddo! npdf loop
  call swapMom(MomExt(1:4,1),MomExt(1:4,2))   ! swap back to original order, for ID below

!------------ 1 LOOP --------------
ELSEIF( CORRECTION.EQ.1 ) THEN
!   stop " This correction not yet implemented for this process"
   do npdf=1,2
      if (npdf .eq. 1) then
         PDFFac(1:2)=PDFFac_a(1:2)
      elseif (npdf .eq. 2) then
         PDFFac(1:2)=PDFFac_b(1:2)
         call swapMom(MomExt(1:4,1),MomExt(1:4,2))
      endif
         
      !   do iHel=nHel(1),nHel(2)
      do iHel=1,NumHelicities
         call HelCrossing(Helicities(iHel,1:NumExtParticles))
         call SetPolarizations()
         if( ZDecays.gt.0 ) then
            if( ExtParticle(5)%Helicity.eq.0 ) cycle!   this can be done more elegantly in mod_process
            call ZDecay(ExtParticle(5),DK_LO,MomExt(1:4,12:13))
            call ZGamLCoupl(1,Helicities(iHel,5),couplZLL,couplGLL)  ! charged lept
      endif

      if (npdf .eq. 2) then
         ! change helicities of the massless quarks for the couplings to Z          
         Helicities(iHel,3)=-Helicities(iHel,3)
         Helicities(iHel,4)=-Helicities(iHel,4)
      endif

      call ZGamQcoupl(Up_,Helicities(iHel,3),couplZUU,couplGUU)
      ! one could here add a routine for the anomalous ttZ coupling, or modify this one. For now, use this for ttZ coupl
      call ZGamQcoupl(Dn_,Helicities(iHel,3),couplZDD,couplGDD)

      ! RR -- not sure about this -- see if works for LO
      couplZQQ_left_dyn=one
      couplZQQ_right_dyn=one

      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo

      LO_Res_Pol = (0d0,0d0)
               
      if (Zdecays .eq. 0) then
         LOPartAmp(up)=BornAmps(1)%Result+BornAmps(2)%Result*( couplZUU )
         LOPartAmp(dn)=BornAmps(1)%Result+BornAmps(2)%Result*( couplZDD )
      elseif( Zdecays.lt.10 ) then
         LOPartAmp(up)=BornAmps(1)%Result+BornAmps(2)%Result*( couplZUU*propZ*couplZLL )
         LOPartAmp(dn)=BornAmps(1)%Result+BornAmps(2)%Result*( couplZDD*propZ*couplZLL )
      elseif( Zdecays.gt.10 ) then
         LOPartAmp(up)=BornAmps(1)%Result+BornAmps(2)%Result*( couplZUU*propZ*couplZLL + couplGUU*couplGLL )
         LOPartAmp(dn)=BornAmps(1)%Result+BornAmps(2)%Result*( couplZDD*propZ*couplZLL + couplGDD*couplGLL )
      endif

      LO_Res_Pol =  ColLO_ttbqqb(1,1) * ( LOPartAmp(up)*dconjg(LOPartAmp(up))*PDFFac(up) + LOPartAmp(dn)*dconjg(LOPartAmp(dn))*PDFFac(dn))

      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol


!! ------------ bosonic loops --------------

! since there are no sisters, I think we can do this in the usual way...

      do iPrimAmp=1,NumPrimAmps
         call SetKirill(PrimAmps(iPrimAmp))
         call PentCut(PrimAmps(iPrimAmp))
         call QuadCut(PrimAmps(iPrimAmp))
         call TripCut(PrimAmps(iPrimAmp))
         call DoubCut(PrimAmps(iPrimAmp))
         call SingCut(PrimAmps(iPrimAmp))
         call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
      enddo
! the prims are nevertheless not true prims, since not EW gauge inv, so now combine them:
          PrimAmps(1)%Result=PrimAmps(1)%Result+PrimAmps(2)%Result
          PrimAmps(3)%Result=PrimAmps(3)%Result+PrimAmps(4)%Result
          PrimAmps(5)%Result=PrimAmps(5)%Result+PrimAmps(6)%Result
          PrimAmps(7)%Result=PrimAmps(7)%Result+PrimAmps(8)%Result

          ListPrimAmps=(/1,3,5,7/)
          do iPrimAmp=1,4
             APrimAmp=ListPrimAmps(iPrimAmp)
             call RenormalizeUV(PrimAmps(APrimAmp),BornAmps(APrimAmp),MuRen**2)
             PrimAmps(APrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(APrimAmp)%Result(-2:1)
             call OneLoopDiv(PrimAmps(APrimAmp),MuRen**2,3,rdiv(2),rdiv(1))


             print *, 'APrimAmp',APrimAmp,PrimAmps(APrimAmp)%Result(-2:1)
             print *, 'BornAmp', APrimAmp,BornAmps(APrimAmp)%Result
             print *, 'fin', PrimAmps(APrimAmp)%Result(0)+PrimAmps(APrimAmp)%Result(1)
             print *, 'fin/Born', (PrimAmps(APrimAmp)%Result(0)+PrimAmps(APrimAmp)%Result(1))/BornAmps(APrimAmp)%Result
             print *, 'DP',APrimAmp,PrimAmps(APrimAmp)%Result(-2)/BornAmps(APrimAmp)%Result
             print *, 'SP',APrimAmp,PrimAmps(APrimAmp)%Result(-1)/BornAmps(APrimAmp)%Result
             print *, 'CC',APrimAmp,PrimAmps(APrimAmp)%Result(0)/BornAmps(APrimAmp)%Result
             print *, 'R',APrimAmp,PrimAmps(APrimAmp)%Result(1)/BornAmps(APrimAmp)%Result
             print *, rdiv(1:2)
             pause

       enddo

!
!!           call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
!!DEC$ IF (_QuadPrecImpr==1)
!          AccPoles = CheckPoles(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv(1:2))
!          if( AccPoles.gt.1d-4 ) then
!!              print *, ""
!!               print *, "QP re-run",AccPoles
!!             coeff4_128(:,:) = qcmplx( coeff4(:,:) )
!!             coeff5_128(:,:) = qcmplx( coeff5(:,:) )
!              call PentCut_128(PrimAmps(iPrimAmp))
!              call QuadCut_128(PrimAmps(iPrimAmp))
!              call TripCut_128(PrimAmps(iPrimAmp))
!              call DoubCut_128(PrimAmps(iPrimAmp))
!              call SingCut_128(PrimAmps(iPrimAmp))
!              call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
!              call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen*2)
!              PrimAmps(iPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(iPrimAmp)%Result(-2:1)
!              AccPoles = CheckPoles(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv(1:2))
!              if( AccPoles.gt.1d-3 ) then
!                  print *, "SKIP",AccPoles
!!                 call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv,(/EHat/))
!                  EvalCS_1L_ttbqqbZ = 0d0
!                  SkipCounter = SkipCounter + 1
!                  return
!              endif
!          endif
!!DEC$ ENDIF
!      enddo
!
!      BosonicPartAmp(up,-2:1)=  &
!                             +   Nc * ( PrimAmps(PrimAmp1_15234)%Result(-2:1) + Q_up/Q_top*PrimAmps(PrimAmp1_12354)%Result(-2:1) )&
!                             - 2d0/Nc*( PrimAmps(PrimAmp1_15234)%Result(-2:1) + Q_up/Q_top*PrimAmps(PrimAmp1_12354)%Result(-2:1) &
!                                      + PrimAmps(PrimAmp1_15243)%Result(-2:1) - Q_up/Q_top*PrimAmps(PrimAmp1_12453)%Result(-2:1) )  &
!                             - 1d0/Nc*( PrimAmps(PrimAmp3_15432)%Result(-2:1) + PrimAmps(PrimAmp3_14532)%Result(-2:1) + PrimAmps(PrimAmp3_14352)%Result(-2:1) &
!                                      - Q_up/Q_top*PrimAmps(PrimAmp3_14532)%Result(-2:1) + PrimAmps(PrimAmp4_15234)%Result(-2:1)  &
!                                      - Q_up/Q_top*(PrimAmps(PrimAmp4_12534)%Result(-2:1)+PrimAmps(PrimAmp4_12345)%Result(-2:1)+PrimAmps(PrimAmp4_15234)%Result(-2:1)) )
!
!      BosonicPartAmp(dn,-2:1)=  &
!                             +   Nc * ( PrimAmps(PrimAmp1_15234)%Result(-2:1) + Q_dn/Q_top*PrimAmps(PrimAmp1_12354)%Result(-2:1) )&
!                             - 2d0/Nc*( PrimAmps(PrimAmp1_15234)%Result(-2:1) + Q_dn/Q_top*PrimAmps(PrimAmp1_12354)%Result(-2:1) &
!                                      + PrimAmps(PrimAmp1_15243)%Result(-2:1) - Q_dn/Q_top*PrimAmps(PrimAmp1_12453)%Result(-2:1) )  &
!                             - 1d0/Nc*( PrimAmps(PrimAmp3_15432)%Result(-2:1) + PrimAmps(PrimAmp3_14532)%Result(-2:1) + PrimAmps(PrimAmp3_14352)%Result(-2:1) &
!                                      - Q_dn/Q_top*PrimAmps(PrimAmp3_14532)%Result(-2:1) + PrimAmps(PrimAmp4_15234)%Result(-2:1) &
!                                      - Q_dn/Q_top*(PrimAmps(PrimAmp4_12534)%Result(-2:1)+PrimAmps(PrimAmp4_12345)%Result(-2:1)+PrimAmps(PrimAmp4_15234)%Result(-2:1)) )
!
!
!      NLO_Res_Pol(-2:1) = Col1L_ttbqqb(1,1) *( dreal(LOPartAmp(up)*dconjg(BosonicPartAmp(up,-2:1)))*PDFFac(up) &
!                                             + dreal(LOPartAmp(dn)*dconjg(BosonicPartAmp(dn,-2:1)))*PDFFac(dn) )
!      NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1) + NLO_Res_Pol(-2:1)
!
!
!
!! ------------ fermionic loops --------------
!      do iPrimAmp=11,14
!          call SetKirill(PrimAmps(iPrimAmp))
!          call PentCut(PrimAmps(iPrimAmp))
!          call QuadCut(PrimAmps(iPrimAmp))
!          call TripCut(PrimAmps(iPrimAmp))
!          call DoubCut(PrimAmps(iPrimAmp))
!          call SingCut(PrimAmps(iPrimAmp))
!          call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
!          PrimAmps(iPrimAmp)%Result(-2:1) = -(0d0,1d0)*PrimAmps(iPrimAmp)%Result(-2:1) !minus if from closed fermion loop
!!           call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,rdiv(2),rdiv(1))
!!           call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv)
!      enddo
!
!!     massless closed fermion loops
!      FermionPartAmp(up,-2:1) = Nf_light*( PrimAmps(PrimAmp2_15234)%Result(-2:1) + Q_up/Q_top*PrimAmps(PrimAmp2_12354)%Result(-2:1) ) &
!                                         + PrimAmps(PrimAmp2m_15234)%Result(-2:1) + Q_up/Q_top*PrimAmps(PrimAmp2m_12354)%Result(-2:1)
!
!
!      FermionPartAmp(dn,-2:1) = Nf_light*( PrimAmps(PrimAmp2_15234)%Result(-2:1) + Q_dn/Q_top*PrimAmps(PrimAmp2_12354)%Result(-2:1) ) &
!                                         + PrimAmps(PrimAmp2m_15234)%Result(-2:1) + Q_dn/Q_top*PrimAmps(PrimAmp2m_12354)%Result(-2:1)
!
!
!      NLO_Res_Pol(-2:1) = Col1L_ttbqqb(1,1) *( dreal(LOPartAmp(up)*dconjg(FermionPartAmp(up,-2:1)))*PDFFac(up) &
!                                             + dreal(LOPartAmp(dn)*dconjg(FermionPartAmp(dn,-2:1)))*PDFFac(dn) )
!      NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1) + NLO_Res_Pol(-2:1)
!
    enddo!helicity loop
 enddo ! npdf
  call swapMom(MomExt(1:4,1),MomExt(1:4,2))   ! swap back to original order, for ID below
!! print *, "mom swap deactivated"
ENDIF
stop



IF( CORRECTION.EQ.0 ) THEN
!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * alpha4Pi * WidthExpansion
   EvalCS_1L_ttbqqbZ = LO_Res_Unpol * PreFac


ELSEIF( CORRECTION.EQ.1 ) THEN
!  overall normalization: (4*Pi)^eps/Gamma(1-eps)
!  CT contributions                           ! beta           !top WFRC
   NLO_Res_UnPol(-1) = NLO_Res_UnPol(-1) + (-11d0/3d0*3d0 - 3d0*4d0/3d0)*LO_Res_Unpol
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-3d0*4d0/3d0)*2d0*dlog(MuRen/m_top)*LO_Res_Unpol  ! finite log(mu2) contrib. from  top WFRC

   NLO_Res_UnPol_Ferm(-1) = NLO_Res_UnPol_Ferm(-1) + (+2d0/3d0*Nf_light+2d0/3d0*Nf_heavy)*LO_Res_Unpol
   NLO_Res_UnPol_Ferm( 0) = NLO_Res_UnPol_Ferm( 0) + (2d0/3d0*Nf_heavy)*2d0*dlog(MuRen/m_top)*LO_Res_Unpol  ! finite log(mu2) contrib. from heavy flavor in alpha_s ren.

   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-5d0/2d0*8d0/3d0 )*LO_Res_Unpol   ! finite contribution from top WFRC's
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + LO_Res_Unpol ! shift alpha_s^DR --> alpha_s^MSbar

!  normalization
   LO_Res_Unpol = LO_Res_Unpol                         * ISFac * (alpha_s4Pi*RunFactor)**2 * Q_top**2*alpha4Pi*PhotonCouplCorr
   NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1)           * ISFac * (alpha_s4Pi*RunFactor)**2 * Q_top**2*alpha4Pi*PhotonCouplCorr * alpha_sOver2Pi*RunFactor
   NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) * ISFac * (alpha_s4Pi*RunFactor)**2 * Q_top**2*alpha4Pi*PhotonCouplCorr * alpha_sOver2Pi*RunFactor


   EvalCS_1L_ttbqqbZ = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac



ELSEIF( CORRECTION.EQ.3 ) THEN
! print *, "1-loop eps2:",(NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2) )* PreFac,  (NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2))/(alpha_sOver2Pi*RunFactor*LO_Res_Unpol)
! print *, "1-loop eps1:",(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1) )* PreFac,  (NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/(alpha_sOver2Pi*RunFactor*LO_Res_Unpol)
! print *, "tree virt",LO_Res_Unpol/RunFactor**2


   stop " This correction not yet implemented for this process"
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   IF( TOPDECAYS.GE.1 ) THEN
       xE = yRnd(16+HelSampling)
   ELSEIF( TOPDECAYS.EQ.0 ) THEN
       xE = yRnd(8+HelSampling)
   ENDIF

! xe=0.45d0;print *, "fixed xe!!"
   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

   IF( PROCESS.EQ.30 ) THEN
!      call EvalIntDipoles_QQBTTBGP((/MomExt(1:4,4),MomExt(1:4,3),MomExt(1:4,5),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,6:11),xE,HOp(1:2,1:3))
      HOp(1:2,1:3) = HOp(1:2,1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ttbqqbZ= HOp(1,1)    * (pdf(Up_,1)*pdf(AUp_,2)+pdf(Chm_,1)*pdf(AChm_,2) ) &
                        +HOp(1,2)/xE * (pdf_z(Up_,1)*pdf(AUp_,2)+pdf_z(Chm_,1)*pdf(AChm_,2) ) &
                        +HOp(1,3)/xE * (pdf(Up_,1)*pdf_z(AUp_,2)+pdf(Chm_,1)*pdf_z(AChm_,2) ) &
                        +HOp(2,1)    * (pdf(Dn_,1)*pdf(ADn_,2)+pdf(Str_,1)*pdf(AStr_,2)+pdf(Bot_,1)*pdf(ABot_,2) ) &
                        +HOp(2,2)/xE * (pdf_z(Dn_,1)*pdf(ADn_,2)+pdf_z(Str_,1)*pdf(AStr_,2)+pdf_z(Bot_,1)*pdf(ABot_,2) ) &
                        +HOp(2,3)/xE * (pdf(Dn_,1)*pdf_z(ADn_,2)+pdf(Str_,1)*pdf_z(AStr_,2)+pdf(Bot_,1)*pdf_z(ABot_,2) )

!      call EvalIntDipoles_QQBTTBGP((/MomExt(1:4,4),MomExt(1:4,3),MomExt(1:4,5),-MomExt(1:4,2),-MomExt(1:4,1)/),MomExt(1:4,6:11),xE,HOp(1:2,1:3))
      HOp(1:2,1:3) = HOp(1:2,1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ttbqqbZ= EvalCS_1L_ttbqqbZ        &
                        +HOp(1,1)    * (pdf(Up_,2)*pdf(AUp_,1)+pdf(Chm_,2)*pdf(AChm_,1) ) &
                        +HOp(1,2)/xE * (pdf_z(Up_,2)*pdf(AUp_,1)+pdf_z(Chm_,2)*pdf(AChm_,1) ) &
                        +HOp(1,3)/xE * (pdf(Up_,2)*pdf_z(AUp_,1)+pdf(Chm_,2)*pdf_z(AChm_,1) ) &
                        +HOp(2,1)    * (pdf(Dn_,2)*pdf(ADn_,1)+pdf(Str_,2)*pdf(AStr_,1)+pdf(Bot_,2)*pdf(ABot_,1) ) &
                        +HOp(2,2)/xE * (pdf_z(Dn_,2)*pdf(ADn_,1)+pdf_z(Str_,2)*pdf(AStr_,1)+pdf_z(Bot_,2)*pdf(ABot_,1) ) &
                        +HOp(2,3)/xE * (pdf(Dn_,2)*pdf_z(ADn_,1)+pdf(Str_,2)*pdf_z(AStr_,1)+pdf(Bot_,2)*pdf_z(ABot_,1) )



   ELSEIF( PROCESS.EQ.24 ) THEN
!      call EvalIntDipoles_QGTTBQP((/MomExt(1:4,4),MomExt(1:4,3),MomExt(1:4,5),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,6:11),xE,HOp(1:2,1:3))
      HOp(1:2,1:3) = HOp(1:2,1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ttbqqbZ=HOp(1,1)    *  (pdf(Up_,1)*pdf(0,2)+pdf(Chm_,1)*pdf(0,2) ) &
                        +HOp(1,2)/xE * (pdf_z(Up_,1)*pdf(0,2)+pdf_z(Chm_,1)*pdf(0,2) ) &
                        +HOp(1,3)/xE * (pdf(Up_,1)*pdf_z(0,2)+pdf(Chm_,1)*pdf_z(0,2) ) &
                        +HOp(2,1)    * (pdf(Dn_,1)*pdf(0,2)+pdf(Str_,1)*pdf(0,2)+pdf(Bot_,1)*pdf(0,2) ) &
                        +HOp(2,2)/xE * (pdf_z(Dn_,1)*pdf(0,2)+pdf_z(Str_,1)*pdf(0,2)+pdf_z(Bot_,1)*pdf(0,2) ) &
                        +HOp(2,3)/xE * (pdf(Dn_,1)*pdf_z(0,2)+pdf(Str_,1)*pdf_z(0,2)+pdf(Bot_,1)*pdf_z(0,2) )

!      call EvalIntDipoles_QGTTBQP((/MomExt(1:4,4),MomExt(1:4,3),MomExt(1:4,5),-MomExt(1:4,2),-MomExt(1:4,1)/),MomExt(1:4,6:11),xE,HOp(1:2,1:3))
      HOp(1:2,1:3) = HOp(1:2,1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ttbqqbZ= EvalCS_1L_ttbqqbZ  &
                        +HOp(1,1)    * (pdf(Up_,2)*pdf(0,1)+pdf(Chm_,2)*pdf(0,1) ) &
                        +HOp(1,2)/xE * (pdf_z(Up_,2)*pdf(0,1)+pdf_z(Chm_,2)*pdf(0,1) ) &
                        +HOp(1,3)/xE * (pdf(Up_,2)*pdf_z(0,1)+pdf(Chm_,2)*pdf_z(0,1) ) &
                        +HOp(2,1)    * (pdf(Dn_,2)*pdf(0,1)+pdf(Str_,2)*pdf(0,1)+pdf(Bot_,2)*pdf(0,1) ) &
                        +HOp(2,2)/xE * (pdf_z(Dn_,2)*pdf(0,1)+pdf_z(Str_,2)*pdf(0,1)+pdf_z(Bot_,2)*pdf(0,1) ) &
                        +HOp(2,3)/xE * (pdf(Dn_,2)*pdf_z(0,1)+pdf(Str_,2)*pdf_z(0,1)+pdf(Bot_,2)*pdf_z(0,1) )


   ELSEIF( PROCESS.EQ.26 ) THEN
!      call EvalIntDipoles_QBGTTBQBP((/MomExt(1:4,4),MomExt(1:4,3),MomExt(1:4,5),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,6:11),xE,HOp(1:2,1:3))
      HOp(1:2,1:3) = HOp(1:2,1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ttbqqbZ= HOp(1,1)    * (pdf(AUp_,1)*pdf(0,2)+pdf(AChm_,1)*pdf(0,2) ) &
                        +HOp(1,2)/xE * (pdf_z(AUp_,1)*pdf(0,2)+pdf_z(AChm_,1)*pdf(0,2) ) &
                        +HOp(1,3)/xE * (pdf(AUp_,1)*pdf_z(0,2)+pdf(AChm_,1)*pdf_z(0,2) ) &
                        +HOp(2,1)    * (pdf(ADn_,1)*pdf(0,2)+pdf(AStr_,1)*pdf(0,2)+pdf(ABot_,1)*pdf(0,2) ) &
                        +HOp(2,2)/xE * (pdf_z(ADn_,1)*pdf(0,2)+pdf_z(AStr_,1)*pdf(0,2)+pdf_z(ABot_,1)*pdf(0,2) ) &
                        +HOp(2,3)/xE * (pdf(ADn_,1)*pdf_z(0,2)+pdf(AStr_,1)*pdf_z(0,2)+pdf(ABot_,1)*pdf_z(0,2) )


!      call EvalIntDipoles_QBGTTBQBP((/MomExt(1:4,4),MomExt(1:4,3),MomExt(1:4,5),-MomExt(1:4,2),-MomExt(1:4,1)/),MomExt(1:4,6:11),xE,HOp(1:2,1:3))
      HOp(1:2,1:3) = HOp(1:2,1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ttbqqbZ= EvalCS_1L_ttbqqbZ &
                        +HOp(1,1)    * (pdf(AUp_,2)*pdf(0,1)+pdf(AChm_,2)*pdf(0,1) ) &
                        +HOp(1,2)/xE * (pdf_z(AUp_,2)*pdf(0,1)+pdf_z(AChm_,2)*pdf(0,1) ) &
                        +HOp(1,3)/xE * (pdf(AUp_,2)*pdf_z(0,1)+pdf(AChm_,2)*pdf_z(0,1) ) &
                        +HOp(2,1)    * (pdf(ADn_,2)*pdf(0,1)+pdf(AStr_,2)*pdf(0,1)+pdf(ABot_,2)*pdf(0,1) ) &
                        +HOp(2,2)/xE * (pdf_z(ADn_,2)*pdf(0,1)+pdf_z(AStr_,2)*pdf(0,1)+pdf_z(ABot_,2)*pdf(0,1) ) &
                        +HOp(2,3)/xE * (pdf(ADn_,2)*pdf_z(0,1)+pdf(AStr_,2)*pdf_z(0,1)+pdf(ABot_,2)*pdf_z(0,1) )

   ENDIF

! print *, "real singularitites",EvalCS_1L_ttbqqbZ
! print *, "real singularitites",EvalCS_1L_ttbqqbZ/(alpha_sOver2Pi*RunFactor*LO_Res_Unpol* PreFac)
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

!      careful, alphas=0.13 only for NLOParam=0 PDFSet=2
!      MADGRAPH CHECK: uub->ttbZ, mt=172, alpha_s=0.13 mZ=91.19
!if (ZDecays .eq. 0) then 
!   MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!   MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!   MG_MOM(0:3,3) = MomExt(1:4,5)*100d0
!   MG_MOM(0:3,4) = MomExt(1:4,4)*100d0
!   MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!   call coupsm(0)
!   print *, "My tree:         ", LO_Res_Unpol/(100d0)**2
!   call SUUB_TTBZ(MG_MOM,MadGraph_tree)
!   print *, "MadGraph hel.amp uub:", MadGraph_tree
!   print *, "MG/ME ratio: ", MadGraph_tree/(dble(LO_Res_Unpol)/(100d0)**2)
!   call SDDB_TTBZ(MG_MOM,MadGraph_tree)
!   print *, "MadGraph hel.amp ddb:", MadGraph_tree
!   print *, "MG/ME ratio: ", MadGraph_tree/(dble(LO_Res_Unpol)/(100d0)**2)
!else
!   print *, "My tree:         ", LO_Res_Unpol/(100d0)**4
!   print *, "MG/ME ratio uub: ", 0.287320393640445d-12/(LO_Res_Unpol/(100d0)**4)
!   print *, "MG/ME ratio ddb: ", 0.189323066770051d-12/(LO_Res_Unpol/(100d0)**4)
!   print *, "MG/ME ratio ubu: ", 0.184639487297707d-12/(LO_Res_Unpol/(100d0)**4)
!   print *, "MG/ME ratio dbd: ", 0.819198133023360d-13/(LO_Res_Unpol/(100d0)**4)
!   !        print *, "MG/ME ratio: ", MadGraph_tree/(dble(LO_Res_Unpol)/(100d0)**2)
!   
!
!   print *, ""
!endif
!   stop

   if( IsNan(EvalCS_1L_ttbqqbZ) ) then
        print *, "NAN:",EvalCS_1L_ttbqqbZ
        print *, yRnd(:)
        print *, NLO_Res_UnPol(0),NLO_Res_UnPol(1),NLO_Res_UnPol_Ferm(0),NLO_Res_UnPol_Ferm(1)
        print *, PSWgt , VgsWgt , PDFFac_a,PDFFac_b, sHatJacobi
        print *, eta1,eta2,((1d0-eta1)*xE+eta1),((1d0-eta2)*xE+eta2),MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,:)
        print *, "SKIP EVENT!!!!!"
        EvalCS_1L_ttbqqbZ = 0d0
        return
   endif


   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_ttbqqbZ)
   enddo



   EvalCS_1L_ttbqqbZ = EvalCS_1L_ttbqqbZ/VgsWgt

RETURN
END FUNCTION








FUNCTION EvalCS_Real_ttbgggZ(yRnd,VgsWgt)
use ModParameters
use ModKinematics
use ModAmplitudes
use ModMisc
use ModProcess
use ModDipoles_GGTTBGZ
use ModZDecay
implicit none
real(8) ::  EvalCS_Real_ttbgggZ,yRnd(1:VegasMxDim),VgsWgt,DipoleResult
complex(8) :: LO_Res_Pol,LO_Res_Unpol,PartAmp(1:4)
integer :: iHel,jPrimAmp,iPrimAmp,NHisto,NBin(1:NumMaxHisto)
real(8) :: EHat,PSWgt,PSWgt2,PSWgt3,PSWgt4,ISFac,RunFactor,PreFac,sij
real(8) :: eta1,eta2,sHatJacobi,FluxFac,PDFFac,MZ_Inv
real(8) :: MomExt(1:4,1:14),pdf(-6:6,1:2)
real(8) :: MG_MOM(0:3,1:6)
real(8) :: MadGraph_tree
logical :: applyPSCut,applySingCut
include "vegas_common.f"



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
! print *, "fixing yrnds"


  EvalCS_Real_ttbgggZ= 0d0
  DipoleResult = 0d0
  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top+M_Z) then
      EvalCS_Real_ttbgggZ = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

   if( ZDecays.le.10 ) then  ! decaying on-shell Z
      MZ_Inv = m_Z
   elseif( ZDecays.gt.10 ) then  ! decaying off-shell Z
      call Error("need to implement phase space for off-shell Z's")
      ! need to think about threshold cut: EHat.le.2d0*m_Top+M_Z   when z is off-shell! 
   endif

   call EvalPhaseSpace_2to4M(EHat,MZ_Inv,yRnd(3:7),MomExt(1:4,1:6),PSWgt)! gluon gluon gluon Z tb t
   call boost2Lab(eta1,eta2,6,MomExt(1:4,1:6))
   ISFac = MomCrossing(MomExt)


   PSWgt2 = 1d0
   PSWgt3 = 1d0
   PSWgt4 = 1d0
IF( TopDecays.GE.1 ) THEN
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(11:14),.false.,MomExt(1:4,7:9),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,6),yRnd(15:18),.false.,MomExt(1:4,10:12),PSWgt3)
   PSWgt = PSWgt * PSWgt2*PSWgt3
   call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,7:9))
   call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,10:12))
ENDIF
IF( ZDECAYS.NE.0 ) THEN
   call EvalPhasespace_ZDecay(MZ_Inv,MomExt(1:4,4),yRnd(19:20),MomExt(1:4,13:14),PSWgt4)
   PSWgt = PSWgt * PSWgt4
ENDIF

   call CheckSing(MomExt,applySingCut)
   if( applySingCut ) then
       EvalCS_Real_ttbgggZ = 0d0
       SkipCounter = SkipCounter + 1
       return
   endif

   call Kinematics_TTBARZ(1,MomExt(1:4,1:14),(/5,6,4,1,2,3,7,8,9,10,11,12,13,14/),applyPSCut,NBin)
   call setPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   RunFactor = RunAlphaS(NLOParam,MuRen)

   if( applyPSCut ) then
       EvalCS_Real_ttbgggZ = 0d0
   else
        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities
          call HelCrossing(Helicities(iHel,1:NumExtParticles))
          call SetPolarizations()
          if( ZDecays.ne.0 ) then
              if( ExtParticle(6)%Helicity.eq.0 ) cycle!   this can be more elegantly done in mod_process
              call ZDecay(ExtParticle(6),DK_LO,MomExt(1:4,13:14))
          endif
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

        LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**3 * alpha4Pi
        EvalCS_Real_ttbgggZ = LO_Res_Unpol * PreFac

        do NHisto=1,NumHistograms
               call intoHisto(NHisto,NBin(NHisto),EvalCS_Real_ttbgggZ)
        enddo
        EvalCounter = EvalCounter + 1
endif!applyPSCut


    PreFac = PreFac * ISFac * (alpha_s4Pi*RunFactor)**3 * alpha4Pi /PSWgt2/PSWgt3/PSWgt4
    call EvalDipoles_GGTTBGZ((/MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,6),-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,3)/),yRnd(11:20),PreFac,DipoleResult)
! DipoleResult=0d0
!      sij = 2d0*(MomExt(1:4,2).dot.MomExt(1:4,3))
! !      sij = MomExt(1,3)**2
!      print *,  sij/EHat**2,EvalCS_Real_ttbgggZ,DipoleResult,(1d0+EvalCS_Real_ttbgggZ/DipoleResult)
!      pause




!      careful, alphas=0.13 only for NLOParam=0 PDFSet=2
!      MADGRAPH CHECK: gg->ttbZ, mt=172, alpha_s=0.13 mZ=91.19
! if (ZDecays .eq. 0) then
!       MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!       MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!       MG_MOM(0:3,3) = MomExt(1:4,6)*100d0
!       MG_MOM(0:3,4) = MomExt(1:4,5)*100d0
!       MG_MOM(0:3,5) = MomExt(1:4,4)*100d0
!       MG_MOM(0:3,6) = MomExt(1:4,3)*100d0
!       call coupsm(0)
!       call SGG_TTBZG(MG_MOM,MadGraph_tree)
!       print *, ""
!       print *, alpha_s*RunFactor,m_top,m_z
!       print *, "My tree:         ", LO_Res_Unpol/(100d0)**2
!       print *, "MadGraph hel.amp:", MadGraph_tree
!       print *, "MG/ME ratio: ", MadGraph_tree/(dble(LO_Res_Unpol)/(100d0)**2)
!       pause
!    else
!       print *, "My tree:         ", LO_Res_Unpol/(100d0)**4
!       print *, "MG/ME ratio: ", 0.872680470745814d-13/(LO_Res_Unpol/(100d0)**4)
!    endif
!    stop
       

    EvalCS_Real_ttbgggZ = (EvalCS_Real_ttbgggZ + DipoleResult)/VgsWgt
RETURN
END FUNCTION







FUNCTION EvalCS_Real_ttbqqbgZ(yRnd,VgsWgt)
use ModParameters
use ModKinematics
use ModAmplitudes
use ModProcess
use ModMisc
use ModDipoles_QQBTTBGZ
use ModDipoles_QGTTBQZ
use ModDipoles_QBGTTBQBZ
implicit none
real(8) ::  EvalCS_Real_ttbqqbgZ,EvalCS_Dips_ttbqqbgZ,yRnd(1:VegasMxDim),VgsWgt,DipoleResult(1:2)
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



EvalCS_Real_ttbqqbgZ= 0d0
EvalCS_Dips_ttbqqbgZ= 0d0

   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if( EHat.le.2d0*m_Top+pT_pho_cut ) then
      EvalCS_Real_ttbqqbgZ = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to4(EHat,yRnd(3:10),MomExt(1:4,1:6),PSWgt)
   call boost2Lab(eta1,eta2,6,MomExt(1:4,1:6))
   call CheckSing(MomExt,applySingCut)
   if( applySingCut ) then
       EvalCS_Real_ttbqqbgZ = 0d0
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
            EvalCS_Real_ttbqqbgZ = 0d0
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

        LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**3 * Q_top**2*alpha4Pi  *PhotonCouplCorr
        EvalCS_Real_ttbqqbgZ = EvalCS_Real_ttbqqbgZ + dble(LO_Res_Unpol*PreFac)

        do NHisto=1,NumHistograms
               call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol*PreFac))
        enddo
endif!applyPSCut

    PreFacDip = PreFac * ISFac * (alpha_s4Pi*RunFactor)**3 * Q_top**2*alpha4Pi*PhotonCouplCorr /PSWgt2/PSWgt3
    IF( PROCESS.EQ.30 ) THEN
        call EvalDipoles_QQBTTBGZ((/MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,6),-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,3)/),yRnd(11:20),(/PreFacDip*PDFFac(1),PreFacDip*PDFFac(2)/),DipoleResult)
    ELSEIF( PROCESS.EQ.24 ) THEN
        call EvalDipoles_QGTTBQZ((/MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,6),-MomExt(1:4,1),MomExt(1:4,3),-MomExt(1:4,2)/),yRnd(11:20),(/PreFacDip*PDFFac(1),PreFacDip*PDFFac(2)/),DipoleResult)
    ELSEIF( PROCESS.EQ.26 ) THEN
        call EvalDipoles_QBGTTBQBZ((/MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,6),MomExt(1:4,3),-MomExt(1:4,1),-MomExt(1:4,2)/),yRnd(11:20),(/PreFacDip*PDFFac(1),PreFacDip*PDFFac(2)/),DipoleResult)
    ENDIF
    EvalCS_Dips_ttbqqbgZ = EvalCS_Dips_ttbqqbgZ + DipoleResult(1) + DipoleResult(2)

  ENDDO! loop over a<-->b pdfs


!      sij = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,3))
!      sij = MomExt(1,3)**2
!      print *,  sij/EHat**2,EvalCS_Real_ttbqqbgZ,EvalCS_Dips_ttbqqbgZ,(1d0+EvalCS_Real_ttbqqbgZ/EvalCS_Dips_ttbqqbgZ)
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

    EvalCS_Real_ttbqqbgZ = (EvalCS_Real_ttbqqbgZ + EvalCS_Dips_ttbqqbgZ) /VgsWgt

END FUNCTION








FUNCTION EvalCS_NLODK_ttbZ(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModHadrWDecay
implicit none
real(8) ::  EvalCS_NLODK_ttbZ,yRnd(1:VegasMxDim),VgsWgt
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




  EvalCS_NLODK_ttbZ = 0d0
  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top+pT_pho_cut ) then
      EvalCS_NLODK_ttbZ = 0d0
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
      EvalCS_NLODK_ttbZ = 0d0
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
   EvalCS_NLODK_ttbZ = EvalCS_NLODK_ttbZ + dble(NLO_Res_Unpol)

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
   EvalCS_NLODK_ttbZ = EvalCS_NLODK_ttbZ + NLO_Res_Unpol


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
      EvalCS_NLODK_ttbZ = 0d0
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
   EvalCS_NLODK_ttbZ = EvalCS_NLODK_ttbZ + dble(LO_Res_Unpol)

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
   EvalCS_NLODK_ttbZ = EvalCS_NLODK_ttbZ + Dip_Res_Unpol

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(Dip_Res_Unpol))
   enddo

enddo! npdf loop
call swapMom(MomExt(1:4,1),MomExt(1:4,2))! swap back

!             print *, MomDK(1,4)/EHat,EvalCS_NLODK_ttbZ/dble(Dip_Res_Unpol)
!             print *, (MomDK(1:4,1).dot.MomDK(1:4,4))/EHat**2,EvalCS_NLODK_ttbZ/dble(Dip_Res_Unpol)
!             EvalCS_NLODK_ttbZ=EvalCS_NLODK_ttbZ/Vgswgt
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
   EvalCS_NLODK_ttbZ = EvalCS_NLODK_ttbZ + dble(LO_Res_Unpol)

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
   EvalCS_NLODK_ttbZ = EvalCS_NLODK_ttbZ + dble(Dip_Res_Unpol)

!             print *, MomDK(1,7)/EHat,pbDpg/EHat**2,EvalCS_NLODK_ttbZ/dble(Dip_Res_Unpol)
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


   EvalCS_NLODK_ttbZ = EvalCS_NLODK_ttbZ/VgsWgt
return

END FUNCTION







END MODULE ModCrossSection_TTBZ


