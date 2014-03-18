MODULE ModCrossSection_ZprimeTTB
use ModTopDecay
implicit none

integer,private,parameter :: NumMaxHisto=45

contains






FUNCTION EvalCS_1L_Zprime_ttbqqb_MPI(yRnd,VgsWgt,res)
implicit none
integer :: EvalCS_1L_Zprime_ttbqqb_MPI
real(8) ::  yRnd(*),res(*),VgsWgt

res(1) = EvalCS_1L_Zprime_ttbqqb(yRnd,VgsWgt)
EvalCS_1L_Zprime_ttbqqb_MPI=0
RETURN
END FUNCTION





FUNCTION EvalCS_1L_Zprime_ttbqqb(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes_Zprime
use ModParameters
use ModIntDipoles_ZprimeTTB
implicit none
real(8) ::  EvalCS_1L_Zprime_ttbqqb,yRnd(1:VegasMxDim),VgsWgt
complex(8) :: LO_Res_Pol_Left, LO_Res_Pol_Right
real(8) :: LO_Res_Unpol_Left, LO_Res_Unpol_Right
complex(8) :: Virt_Res_Pol_Left, Virt_Res_Pol_Right
real(8) :: Virt_Res_Unpol_Left, Virt_Res_Unpol_Right
integer :: iHel
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,ISFac
real(8) :: MomExt(1:4,1:NumExtParticles),MomDK(1:4,1:6)
logical :: applyPSCut
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac_L, PDFFac_R, PDFFac_dip_L(3), PDFFac_dip_R(3)
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2), z
real(8) :: IDip(3), IDipAmp
integer :: NHisto,NBin(1:NumMaxHisto),npdf,nHel(1:2),NRndHel
real(8) :: MadGraphRes
real(8) :: xpdf
include 'misc/global_import'
include "vegas_common.f"

  EvalCS_1L_Zprime_ttbqqb = 0d0

!  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi) ! Flat mapping
  call PDFMapping(62,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi) ! BW mapping for Ehat

  if( EHat.le.2d0*m_Top) then
     EvalCS_1L_Zprime_ttbqqb = 0d0
     return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

  call EvalPhaseSpace_2to2(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt)
  call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))

  NRndHel=5
  IF( TOPDECAYS.NE.0 ) THEN
     call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:8),.false.,MomDK(1:4,1:3),PSWgt2)
     call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(9:12),.false.,MomDK(1:4,4:6),PSWgt3)
     PSWgt = PSWgt * PSWgt2*PSWgt3
     NRndHel=13
  ENDIF


  call Kinematics_TTBARZprime(.false.,MomExt,MomDK,applyPSCut,NBin)

  if( applyPSCut ) then
     EvalCS_1L_Zprime_ttbqqb = 0d0
     return
  endif


  call setPDFs(eta1,eta2,MuFac,pdf)

  call random_number(xpdf)

  if( xpdf.le.0.5d0 ) then
     PDFFac_L = gL_Zpr(up_)**2 * (pdf(up_,1)*pdf(aup_,2))
     PDFFac_L = PDFFac_L + gL_Zpr(dn_)**2 * (pdf(dn_,1)*pdf(adn_,2))
     PDFFac_L = PDFFac_L + gL_Zpr(chm_)**2 * (pdf(chm_,1)*pdf(achm_,2))
     PDFFac_L = PDFFac_L + gL_Zpr(str_)**2 * (pdf(str_,1)*pdf(astr_,2))
     PDFFac_L = PDFFac_L + gL_Zpr(bot_)**2 * (pdf(bot_,1)*pdf(abot_,2))
     
     PDFFac_R = gR_Zpr(up_)**2 * (pdf(up_,1)*pdf(aup_,2))
     PDFFac_R = PDFFac_R + gR_Zpr(dn_)**2 * (pdf(dn_,1)*pdf(adn_,2))
     PDFFac_R = PDFFac_R + gR_Zpr(chm_)**2 * (pdf(chm_,1)*pdf(achm_,2))
     PDFFac_R = PDFFac_R + gR_Zpr(str_)**2 * (pdf(str_,1)*pdf(astr_,2))
     PDFFac_R = PDFFac_R + gR_Zpr(bot_)**2 * (pdf(bot_,1)*pdf(abot_,2))
  else
     PDFFac_L = gL_Zpr(up_)**2 * (pdf(Up_,2)*pdf(aup_,1))
     PDFFac_L = PDFFac_L + gL_Zpr(dn_)**2 * (pdf(dn_,2)*pdf(adn_,1))
     PDFFac_L = PDFFac_L + gL_Zpr(chm_)**2 * (pdf(chm_,2)*pdf(achm_,1))
     PDFFac_L = PDFFac_L + gL_Zpr(str_)**2 * (pdf(str_,2)*pdf(astr_,1))
     PDFFac_L = PDFFac_L + gL_Zpr(bot_)**2 * (pdf(bot_,2)*pdf(abot_,1))
     
     PDFFac_R = gR_Zpr(up_)**2 * (pdf(Up_,2)*pdf(aup_,1))
     PDFFac_R = PDFFac_R + gR_Zpr(dn_)**2 * (pdf(dn_,2)*pdf(adn_,1))
     PDFFac_R = PDFFac_R + gR_Zpr(chm_)**2 * (pdf(chm_,2)*pdf(achm_,1))
     PDFFac_R = PDFFac_R + gR_Zpr(str_)**2 * (pdf(str_,2)*pdf(astr_,1))
     PDFFac_R = PDFFac_R + gR_Zpr(bot_)**2 * (pdf(bot_,2)*pdf(abot_,1))
     call swapMom(MomExt(1:4,1),MomExt(1:4,2))
  endif
  PDFFac_R = PDFFac_R * 2d0!  factor two for sampling over pdfs
  PDFFac_L = PDFFac_L * 2d0


  PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
  RunFactor = RunAlphaS(NLOParam,MuRen)
  nHel(1:2) = getHelicity(yrnd(NRndHel))
  PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))

  LO_Res_Unpol_Left = 0d0
  LO_Res_Unpol_Right = 0d0

  Virt_Res_Unpol_Left = 0d0
  Virt_Res_Unpol_Right = 0d0

  !------------ LO --------------
  IF( CORRECTION.EQ.0 ) THEN

     !Momext(1:4,4) = (/2.73325721108289d0, -0.86989643239260d0, -0.55008440614255d0, 1.85822020357279d0/)
     !Momext(1:4,3) = (/9.03205582647898d0, 0.86989643239260d0, 0.55008440614255d0, -8.80683369864915d0/)
     !Momext(1:4,1) = -(/-2.40834977124275d0, 0d0, 0d0, -2.40834977124275d0/)
     !Momext(1:4,2) = -(/-9.35696326631912d0, 0d0, 0d0, 9.35696326631912d0/)
     !EHat = dsqrt(2d0*(momext(1:4,1).dot.Momext(1:4,2)))

     ISFac = MomCrossing(MomExt)

     IF( TOPDECAYS.GE.1 ) THEN
        !     top decays with spin correlations
        call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
        call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
     ENDIF

     do iHel=nHel(1),nHel(2)
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
        call SetPolarizations()
        call Tree_Zprime_tbtqbq(LO_Res_Pol_Left, LO_Res_Pol_Right)
        LO_Res_UnPol_Left = LO_Res_UnPol_Left + dreal(LO_Res_Pol_Left*dconjg(LO_Res_Pol_Left))
        LO_Res_UnPol_Right = LO_Res_UnPol_Right + dreal(LO_Res_Pol_Right*dconjg(LO_Res_Pol_Right))
     enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! F- MadGraph check, Oct 13th
!
!     open(unit=66,file='/home/caola/tools/MG_ME_V4.5.2/zprime/SubProcesses/kin.dat', status='unknown')
!     do iHel = 1,4
!        write(66,*) MomExt(1,iHel) * 100d0
!        write(66,*) MomExt(2,iHel) * 100d0
!        write(66,*) MomExt(3,iHel) * 100d0
!        write(66,*) MomExt(4,iHel) * 100d0
!     enddo
!     close(66)
!     print *, 'MZpr', m_zpr * 100d0
!     print *, 'Ga_zpr', Ga_Zpr * 100d0
!     print *, 'gL_Zpr(up_), gL_Zpr(dn_)', gL_Zpr(up_), gL_Zpr(dn_)
!     print *, 'gR_Zpr(up_), gR_Zpr(dn_)', gR_Zpr(up_), gR_Zpr(dn_)
!     print *, 'uub->ttb ', ISFac * (gL_Zpr(up_)**2 * LO_Res_Unpol_Left + gR_Zpr(up_)**2 * LO_Res_Unpol_Right)*9d0
!     print *, 'ddb->ttb ', ISFac * (gL_Zpr(dn_)**2 * LO_Res_Unpol_Left + gR_Zpr(dn_)**2 * LO_Res_Unpol_Right)*9d0
!     pause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! M-F checked on May 26th

!     call coupsm(0)
!     call SUUB_ZPRIME_TTB((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,4),MomExt(1:4,3)/)*100d0,MadGraphRes)
!     print *, 'madgraph up', MadGraphRes
!     print *, 'us         ', ISFac * (gL_Zpr(up_)**2 * LO_Res_Unpol_Left + gR_Zpr(up_)**2 * LO_Res_Unpol_Right)*9d0
!     print *, ISFac*(gL_Zpr(up_)**2 * LO_Res_Unpol_Left + gR_Zpr(up_)**2 * LO_Res_Unpol_Right)*9d0/MadGraphRes
!     print *, '' 
!     call SDDB_ZPRIME_TTB((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,4),MomExt(1:4,3)/)*100d0,MadGraphRes)
!     print *, 'madgraph dn', MadGraphRes
!     print *, 'us         ', ISFac * (gL_Zpr(dn_)**2 * LO_Res_Unpol_Left + gR_Zpr(dn_)**2 * LO_Res_Unpol_Right)*9d0
!     print *, ISFac*(gL_Zpr(dn_)**2 * LO_Res_Unpol_Left + gR_Zpr(dn_)**2 * LO_Res_Unpol_Right)*9d0/MadGraphRes
!     pause

     EvalCS_1L_Zprime_ttbqqb = ISFac*PreFac* (PDFFac_L * LO_Res_Unpol_Left + PDFFac_R * LO_Res_Unpol_Right)
     EvalCS_1L_Zprime_ttbqqb = EvalCS_1L_Zprime_ttbqqb * 9d0 ! Sum over initial / final colors





! !     adding the photon ./TOPAZ_MPI Collider=2 TopDK=0 ObsSet=67 Process=62 NLOParam=1 Correction=0  MZpr=0.9119 GaZpr=0.02
!       if( xpdf.le.0.5d0 ) then
!         PDFFac_L = Q_up**2 * (pdf(up_,1)*pdf(aup_,2))
!         PDFFac_L = PDFFac_L + Q_dn**2 * (pdf(dn_,1)*pdf(adn_,2))
!         PDFFac_L = PDFFac_L + Q_up**2 * (pdf(chm_,1)*pdf(achm_,2))
!         PDFFac_L = PDFFac_L + Q_dn**2 * (pdf(str_,1)*pdf(astr_,2))
!         PDFFac_L = PDFFac_L + Q_dn**2 * (pdf(bot_,1)*pdf(abot_,2))
!         
!         PDFFac_R = Q_up**2 * (pdf(up_,1)*pdf(aup_,2))
!         PDFFac_R = PDFFac_R + Q_dn**2 * (pdf(dn_,1)*pdf(adn_,2))
!         PDFFac_R = PDFFac_R + Q_up**2 * (pdf(chm_,1)*pdf(achm_,2))
!         PDFFac_R = PDFFac_R + Q_dn**2 * (pdf(str_,1)*pdf(astr_,2))
!         PDFFac_R = PDFFac_R + Q_dn**2 * (pdf(bot_,1)*pdf(abot_,2))
!       else
!         PDFFac_L = Q_up**2 * (pdf(Up_,2)*pdf(aup_,1))
!         PDFFac_L = PDFFac_L + Q_dn**2 * (pdf(dn_,2)*pdf(adn_,1))
!         PDFFac_L = PDFFac_L + Q_up**2 * (pdf(chm_,2)*pdf(achm_,1))
!         PDFFac_L = PDFFac_L + Q_dn**2 * (pdf(str_,2)*pdf(astr_,1))
!         PDFFac_L = PDFFac_L + Q_dn**2 * (pdf(bot_,2)*pdf(abot_,1))
!         
!         PDFFac_R = Q_up**2 * (pdf(Up_,2)*pdf(aup_,1))
!         PDFFac_R = PDFFac_R + Q_dn**2 * (pdf(dn_,2)*pdf(adn_,1))
!         PDFFac_R = PDFFac_R + Q_up**2 * (pdf(chm_,2)*pdf(achm_,1))
!         PDFFac_R = PDFFac_R + Q_dn**2 * (pdf(str_,2)*pdf(astr_,1))
!         PDFFac_R = PDFFac_R + Q_dn**2 * (pdf(bot_,2)*pdf(abot_,1))
!         call swapMom(MomExt(1:4,1),MomExt(1:4,2))
!       endif
!       do iHel=nHel(1),nHel(2)
!         call HelCrossing(Helicities(iHel,1:NumExtParticles))
!         call SetPolarizations()
!         call Tree_photon_tbtqbq(LO_Res_Pol_Left, LO_Res_Pol_Right)
!         LO_Res_UnPol_Left = LO_Res_UnPol_Left + dreal(LO_Res_Pol_Left*dconjg(LO_Res_Pol_Left))
!         LO_Res_UnPol_Right = LO_Res_UnPol_Right + dreal(LO_Res_Pol_Right*dconjg(LO_Res_Pol_Right))
!       enddo
!       EvalCS_1L_Zprime_ttbqqb = EvalCS_1L_Zprime_ttbqqb  &
!                               + ISFac*PreFac*EL**4* (PDFFac_L * LO_Res_Unpol_Left + PDFFac_R * LO_Res_Unpol_Right)* 9d0 ! Sum over initial / final colors



  ELSEIF(CORRECTION.EQ.1) THEN


     !Momext(1:4,4) = (/2.73325721108289d0, -0.86989643239260d0, -0.55008440614255d0, 1.85822020357279d0/)
     !Momext(1:4,3) = (/9.03205582647898d0, 0.86989643239260d0, 0.55008440614255d0, -8.80683369864915d0/)
     !Momext(1:4,1) = -(/-2.40834977124275d0, 0d0, 0d0, -2.40834977124275d0/)
     !Momext(1:4,2) = -(/-9.35696326631912d0, 0d0, 0d0, 9.35696326631912d0/)
     !EHat = dsqrt(2d0*(momext(1:4,1).dot.Momext(1:4,2)))

     ISFac = MomCrossing(MomExt)

     IF( TOPDECAYS.GE.1 ) THEN
        !     top decays with spin correlations
        call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
        call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
     ENDIF


     do iHel=nHel(1),nHel(2)
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
        call SetPolarizations()
        call Tree_Zprime_tbtqbq(LO_Res_Pol_Left, LO_Res_Pol_Right)
        call Virt_Zprime_tbtqbq(Virt_Res_Pol_Left, Virt_Res_Pol_Right)
        LO_Res_UnPol_Left = LO_Res_UnPol_Left + dreal(LO_Res_Pol_Left*dconjg(LO_Res_Pol_Left))
        LO_Res_UnPol_Right = LO_Res_UnPol_Right + dreal(LO_Res_Pol_Right*dconjg(LO_Res_Pol_Right))
        Virt_Res_UnPol_Left = Virt_Res_UnPol_Left + 2d0 * dreal(LO_Res_Pol_Left * dconjg(Virt_Res_Pol_Left))
        Virt_Res_UnPol_Right = Virt_Res_UnPol_Right + 2d0 * dreal(LO_Res_Pol_Right * dconjg(Virt_Res_Pol_Right))
     enddo!helicity loop

     EvalCS_1L_Zprime_ttbqqb = ISFac*PreFac* (PDFFac_L * Virt_Res_Unpol_Left + PDFFac_R * Virt_Res_Unpol_Right)
     EvalCS_1L_Zprime_ttbqqb = EvalCS_1L_Zprime_ttbqqb * 9d0*(4d0*Pi*alpha_s*RunFactor) ! Color and alpha_s


  !!! Integrated dipoles !!!
  ELSEIF(CORRECTION.EQ.3) THEN

     IF ( TOPDECAYS.EQ.0 ) THEN 
        z = yRnd(5)
     ELSE
        z = yRnd(13)
     END IF

     call setPDFs(eta1/z,eta2/z,MuFac,pdf_z)

     ISFac = MomCrossing(MomExt)

     IF( TOPDECAYS.GE.1 ) THEN
        !     top decays with spin correlations
        call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
        call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
     ENDIF


     do iHel=nHel(1),nHel(2)
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
        call SetPolarizations()
        call Tree_Zprime_tbtqbq(LO_Res_Pol_Left, LO_Res_Pol_Right)
        LO_Res_UnPol_Left = LO_Res_UnPol_Left + dreal(LO_Res_Pol_Left*dconjg(LO_Res_Pol_Left))
        LO_Res_UnPol_Right = LO_Res_UnPol_Right + dreal(LO_Res_Pol_Right*dconjg(LO_Res_Pol_Right))
     enddo!helicity loop


     IF(PROCESS.EQ.63) THEN

        ! qg initial state !

        PDFFac_dip_L(1) = 0d0
        PDFFac_dip_L(2) = 0d0

        if( xpdf.le.0.5d0 ) then
            PDFFac_dip_L(3) =                   gL_Zpr(up_)**2 * (pdf_z(0,2)*pdf(up_,1))! here was a bug: 1<-->2; correct now
            PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(dn_)**2 * (pdf_z(0,2)*pdf(dn_,1))
            PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(chm_)**2 * (pdf_z(0,2)*pdf(chm_,1))
            PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(str_)**2 * (pdf_z(0,2)*pdf(str_,1))
            PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(bot_)**2 * (pdf_z(0,2)*pdf(bot_,1))
            PDFFac_dip_R(3) =                   gR_Zpr(up_)**2 * (pdf_z(0,2)*pdf(up_,1))
            PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(dn_)**2 * (pdf_z(0,2)*pdf(dn_,1))
            PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(chm_)**2 * (pdf_z(0,2)*pdf(chm_,1))
            PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(str_)**2 * (pdf_z(0,2)*pdf(str_,1))
            PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(bot_)**2 * (pdf_z(0,2)*pdf(bot_,1))
        else
            PDFFac_dip_L(3) =                   gL_Zpr(up_)**2 * (pdf_z(0,1)*pdf(up_,2))
            PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(dn_)**2 * (pdf_z(0,1)*pdf(dn_,2))
            PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(chm_)**2 * (pdf_z(0,1)*pdf(chm_,2))
            PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(str_)**2 * (pdf_z(0,1)*pdf(str_,2))
            PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(bot_)**2 * (pdf_z(0,1)*pdf(bot_,2))
            PDFFac_dip_R(3) =                   gR_Zpr(up_)**2 * (pdf_z(0,1)*pdf(up_,2))
            PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(dn_)**2 * (pdf_z(0,1)*pdf(dn_,2))
            PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(chm_)**2 * (pdf_z(0,1)*pdf(chm_,2))
            PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(str_)**2 * (pdf_z(0,1)*pdf(str_,2))
            PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(bot_)**2 * (pdf_z(0,1)*pdf(bot_,2))
        endif

        PDFFac_dip_R(3) = PDFFac_dip_R(3) * 2d0!  factor two for sampling over pdfs
        PDFFac_dip_L(3) = PDFFac_dip_L(3) * 2d0


!        PDFFac_dip_L(3) = gL_Zpr(up_)**2 * (pdf_z(0,1)*pdf(up_,2)+pdf_z(0,2)*pdf(up_,1))
!        PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(dn_)**2 * (pdf_z(0,1)*pdf(dn_,2)+pdf_z(0,2)*pdf(dn_,1))
!        PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(chm_)**2 * (pdf_z(0,1)*pdf(chm_,2)+pdf_z(0,2)*pdf(chm_,1))
!        PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(str_)**2 * (pdf_z(0,1)*pdf(str_,2)+pdf_z(0,2)*pdf(str_,1))
!        PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(bot_)**2 * (pdf_z(0,1)*pdf(bot_,2)+pdf_z(0,2)*pdf(bot_,1))

!        PDFFac_dip_R(3) = gR_Zpr(up_)**2 * (pdf_z(0,1)*pdf(up_,2)+pdf_z(0,2)*pdf(up_,1))
!        PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(dn_)**2 * (pdf_z(0,1)*pdf(dn_,2)+pdf_z(0,2)*pdf(dn_,1))
!        PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(chm_)**2 * (pdf_z(0,1)*pdf(chm_,2)+pdf_z(0,2)*pdf(chm_,1))
!        PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(str_)**2 * (pdf_z(0,1)*pdf(str_,2)+pdf_z(0,2)*pdf(str_,1))
!        PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(bot_)**2 * (pdf_z(0,1)*pdf(bot_,2)+pdf_z(0,2)*pdf(bot_,1))

        call IntDip_qg_Zprime_ttb(MomExt(1:4,1:4), z, IDip) ! Normalization: as/(2 Pi)

        IDip = IDip * (alpha_s*RunFactor)/(2d0*Pi) * PreFac

        IDipAmp = IDip(1) * (PDFFac_dip_L(1) * LO_Res_UnPol_Left + PDFFac_dip_R(1) * LO_Res_UnPol_Right)
        IDipAmp = IDipAmp + IDip(2)/z * (PDFFac_dip_L(2) * LO_Res_UnPol_Left + PDFFac_dip_R(2) * LO_Res_UnPol_Right)
        IDipAmp = IDipAmp + IDip(3)/z * (PDFFac_dip_L(3) * LO_Res_UnPol_Left + PDFFac_dip_R(3) * LO_Res_UnPol_Right)
        
        IDipAmp = IDipAmp * ISFac * 9d0 * (4d0/3d0) * 2d0!  9* CF/TR

        EvalCS_1L_Zprime_ttbqqb = EvalCS_1L_Zprime_ttbqqb + IDipAmp

     ELSEIF(PROCESS.EQ.64) THEN

        ! gqb initial state !

        PDFFac_dip_L(1) = 0d0
        PDFFac_dip_L(3) = 0d0

        if( xpdf.le.0.5d0 ) then
           PDFFac_dip_L(2) = gL_Zpr(up_)**2 * (pdf_z(0,1)*pdf(aup_,2))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(dn_)**2 * (pdf_z(0,1)*pdf(adn_,2))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(chm_)**2 * (pdf_z(0,1)*pdf(achm_,2))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(str_)**2 * (pdf_z(0,1)*pdf(astr_,2))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(bot_)**2 * (pdf_z(0,1)*pdf(abot_,2))
           PDFFac_dip_R(2) = gR_Zpr(up_)**2 * (pdf_z(0,1)*pdf(aup_,2))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(dn_)**2 * (pdf_z(0,1)*pdf(adn_,2))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(chm_)**2 * (pdf_z(0,1)*pdf(achm_,2))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(str_)**2 * (pdf_z(0,1)*pdf(astr_,2))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(bot_)**2 * (pdf_z(0,1)*pdf(abot_,2))
        else
           PDFFac_dip_L(2) = gL_Zpr(up_)**2 * (pdf_z(0,2)*pdf(aup_,1))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(dn_)**2 * (pdf_z(0,2)*pdf(adn_,1))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(chm_)**2 * (pdf_z(0,2)*pdf(achm_,1))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(str_)**2 * (pdf_z(0,2)*pdf(astr_,1))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(bot_)**2 * (pdf_z(0,2)*pdf(abot_,1))
           PDFFac_dip_R(2) = gR_Zpr(up_)**2 * (pdf_z(0,2)*pdf(aup_,1))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(dn_)**2 * (pdf_z(0,2)*pdf(adn_,1))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(chm_)**2 * (pdf_z(0,2)*pdf(achm_,1))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(str_)**2 * (pdf_z(0,2)*pdf(astr_,1))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(bot_)**2 * (pdf_z(0,2)*pdf(abot_,1))
        endif

        PDFFac_dip_R(2) = PDFFac_dip_R(2) * 2d0!  factor two for sampling over pdfs
        PDFFac_dip_L(2) = PDFFac_dip_L(2) * 2d0

!        PDFFac_dip_L(2) = gL_Zpr(up_)**2 * (pdf_z(0,1)*pdf(aup_,2)+pdf_z(0,2)*pdf(aup_,1))
!        PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(dn_)**2 * (pdf_z(0,1)*pdf(adn_,2)+pdf_z(0,2)*pdf(adn_,1))
!        PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(chm_)**2 * (pdf_z(0,1)*pdf(achm_,2)+pdf_z(0,2)*pdf(achm_,1))
!        PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(str_)**2 * (pdf_z(0,1)*pdf(astr_,2)+pdf_z(0,2)*pdf(astr_,1))
!        PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(bot_)**2 * (pdf_z(0,1)*pdf(abot_,2)+pdf_z(0,2)*pdf(abot_,1))
        
!        PDFFac_dip_R(2) = gR_Zpr(up_)**2 * (pdf_z(0,1)*pdf(aup_,2)+pdf_z(0,2)*pdf(aup_,1))
!        PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(dn_)**2 * (pdf_z(0,1)*pdf(adn_,2)+pdf_z(0,2)*pdf(adn_,1))
!        PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(chm_)**2 * (pdf_z(0,1)*pdf(achm_,2)+pdf_z(0,2)*pdf(achm_,1))
!        PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(str_)**2 * (pdf_z(0,1)*pdf(astr_,2)+pdf_z(0,2)*pdf(astr_,1))
!        PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(bot_)**2 * (pdf_z(0,1)*pdf(abot_,2)+pdf_z(0,2)*pdf(abot_,1))
        
        call IntDip_gqb_Zprime_ttb(MomExt(1:4,1:4), z, IDip) ! Normalization: as/(2 Pi)

        IDip = IDip * (alpha_s*RunFactor)/(2d0*Pi) * PreFac

        IDipAmp = IDip(1) * (PDFFac_dip_L(1) * LO_Res_UnPol_Left + PDFFac_dip_R(1) * LO_Res_UnPol_Right)
        IDipAmp = IDipAmp + IDip(2)/z * (PDFFac_dip_L(2) * LO_Res_UnPol_Left + PDFFac_dip_R(2) * LO_Res_UnPol_Right)
        IDipAmp = IDipAmp + IDip(3)/z * (PDFFac_dip_L(3) * LO_Res_UnPol_Left + PDFFac_dip_R(3) * LO_Res_UnPol_Right)
        
        IDipAmp = IDipAmp * ISFac * 9d0 * (4d0/3d0) * 2d0!  9* CF/TR
        
        EvalCS_1L_Zprime_ttbqqb = EvalCS_1L_Zprime_ttbqqb + IDipAmp


     ELSEIF(PROCESS.EQ.66) THEN
        
        ! qqb initial state !

        PDFFac_dip_L(1) = PDFFac_L
        PDFFac_dip_R(1) = PDFFac_R

        if( xpdf.le.0.5d0 ) then
           PDFFac_dip_L(2) = gL_Zpr(up_)**2 * (pdf_z(up_,1)*pdf(aup_,2))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(dn_)**2 * (pdf_z(dn_,1)*pdf(adn_,2))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(chm_)**2 * (pdf_z(chm_,1)*pdf(achm_,2))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(str_)**2 * (pdf_z(str_,1)*pdf(astr_,2))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(bot_)**2 * (pdf_z(bot_,1)*pdf(abot_,2))
           PDFFac_dip_L(3) = gL_Zpr(up_)**2 * (pdf(up_,1)*pdf_z(aup_,2))
           PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(dn_)**2 * (pdf(dn_,1)*pdf_z(adn_,2))
           PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(chm_)**2 * (pdf(chm_,1)*pdf_z(achm_,2))
           PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(str_)**2 * (pdf(str_,1)*pdf_z(astr_,2))
           PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(bot_)**2 * (pdf(bot_,1)*pdf_z(abot_,2))
           PDFFac_dip_R(2) = gR_Zpr(up_)**2 * (pdf_z(up_,1)*pdf(aup_,2))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(dn_)**2 * (pdf_z(dn_,1)*pdf(adn_,2))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(chm_)**2 * (pdf_z(chm_,1)*pdf(achm_,2))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(str_)**2 * (pdf_z(str_,1)*pdf(astr_,2))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(bot_)**2 * (pdf_z(bot_,1)*pdf(abot_,2))
           PDFFac_dip_R(3) = gR_Zpr(up_)**2 * (pdf(up_,1)*pdf_z(aup_,2))
           PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(dn_)**2 * (pdf(dn_,1)*pdf_z(adn_,2))
           PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(chm_)**2 * (pdf(chm_,1)*pdf_z(achm_,2))
           PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(str_)**2 * (pdf(str_,1)*pdf_z(astr_,2))
           PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(bot_)**2 * (pdf(bot_,1)*pdf_z(abot_,2))
        else
           PDFFac_dip_L(2) = gL_Zpr(up_)**2 * (pdf_z(up_,2)*pdf(aup_,1))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(dn_)**2 * (pdf_z(dn_,2)*pdf(adn_,1))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(chm_)**2 * (pdf_z(chm_,2)*pdf(achm_,1))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(str_)**2 * (pdf_z(str_,2)*pdf(astr_,1))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(bot_)**2 * (pdf_z(bot_,2)*pdf(abot_,1))
           PDFFac_dip_L(3) = gL_Zpr(up_)**2 * (pdf(up_,2)*pdf_z(aup_,1))
           PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(dn_)**2 * (pdf(dn_,2)*pdf_z(adn_,1))
           PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(chm_)**2 * (pdf(chm_,2)*pdf_z(achm_,1))
           PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(str_)**2 * (pdf(str_,2)*pdf_z(astr_,1))
           PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(bot_)**2 * (pdf(bot_,2)*pdf_z(abot_,1))
           PDFFac_dip_R(2) = gR_Zpr(up_)**2 * (pdf_z(up_,2)*pdf(aup_,1))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(dn_)**2 * (pdf_z(dn_,2)*pdf(adn_,1))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(chm_)**2 * (pdf_z(chm_,2)*pdf(achm_,1))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(str_)**2 * (pdf_z(str_,2)*pdf(astr_,1))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(bot_)**2 * (pdf_z(bot_,2)*pdf(abot_,1))
           PDFFac_dip_R(3) = gR_Zpr(up_)**2 * (pdf(up_,2)*pdf_z(aup_,1))
           PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(dn_)**2 * (pdf(dn_,2)*pdf_z(adn_,1))
           PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(chm_)**2 * (pdf(chm_,2)*pdf_z(achm_,1))
           PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(str_)**2 * (pdf(str_,2)*pdf_z(astr_,1))
           PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(bot_)**2 * (pdf(bot_,2)*pdf_z(abot_,1))
        endif

        PDFFac_dip_R(2) = PDFFac_dip_R(2) * 2d0!  factor two for sampling over pdfs
        PDFFac_dip_L(2) = PDFFac_dip_L(2) * 2d0

        PDFFac_dip_R(3) = PDFFac_dip_R(3) * 2d0
        PDFFac_dip_L(3) = PDFFac_dip_L(3) * 2d0

!        PDFFac_dip_L(2) = gL_Zpr(up_)**2 * (pdf_z(up_,1)*pdf(aup_,2)+pdf_z(up_,2)*pdf(aup_,1))
!        PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(dn_)**2 * (pdf_z(dn_,1)*pdf(adn_,2)+pdf_z(dn_,2)*pdf(adn_,1))
!        PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(chm_)**2 * (pdf_z(chm_,1)*pdf(achm_,2)+pdf_z(chm_,2)*pdf(achm_,1))
!        PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(str_)**2 * (pdf_z(str_,1)*pdf(astr_,2)+pdf_z(str_,2)*pdf(astr_,1))
!        PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(bot_)**2 * (pdf_z(bot_,1)*pdf(abot_,2)+pdf_z(bot_,2)*pdf(abot_,1))
        
!        PDFFac_dip_L(3) = gL_Zpr(up_)**2 * (pdf(up_,1)*pdf_z(aup_,2)+pdf(up_,2)*pdf_z(aup_,1))
!        PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(dn_)**2 * (pdf(dn_,1)*pdf_z(adn_,2)+pdf(dn_,2)*pdf_z(adn_,1))
!        PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(chm_)**2 * (pdf(chm_,1)*pdf_z(achm_,2)+pdf(chm_,2)*pdf_z(achm_,1))
!        PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(str_)**2 * (pdf(str_,1)*pdf_z(astr_,2)+pdf(str_,2)*pdf_z(astr_,1))
!        PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(bot_)**2 * (pdf(bot_,1)*pdf_z(abot_,2)+pdf(bot_,2)*pdf_z(abot_,1))
        
!        PDFFac_dip_R(1) = PDFFac_R
        
!        PDFFac_dip_R(2) = gR_Zpr(up_)**2 * (pdf_z(up_,1)*pdf(aup_,2)+pdf_z(up_,2)*pdf(aup_,1))
!        PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(dn_)**2 * (pdf_z(dn_,1)*pdf(adn_,2)+pdf_z(dn_,2)*pdf(adn_,1))
!        PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(chm_)**2 * (pdf_z(chm_,1)*pdf(achm_,2)+pdf_z(chm_,2)*pdf(achm_,1))
!        PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(str_)**2 * (pdf_z(str_,1)*pdf(astr_,2)+pdf_z(str_,2)*pdf(astr_,1))
!        PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(bot_)**2 * (pdf_z(bot_,1)*pdf(abot_,2)+pdf_z(bot_,2)*pdf(abot_,1))
        
!        PDFFac_dip_R(3) = gR_Zpr(up_)**2 * (pdf(up_,1)*pdf_z(aup_,2)+pdf(up_,2)*pdf_z(aup_,1))
!        PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(dn_)**2 * (pdf(dn_,1)*pdf_z(adn_,2)+pdf(dn_,2)*pdf_z(adn_,1))
!        PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(chm_)**2 * (pdf(chm_,1)*pdf_z(achm_,2)+pdf(chm_,2)*pdf_z(achm_,1))
!        PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(str_)**2 * (pdf(str_,1)*pdf_z(astr_,2)+pdf(str_,2)*pdf_z(astr_,1))
!        PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(bot_)**2 * (pdf(bot_,1)*pdf_z(abot_,2)+pdf(bot_,2)*pdf_z(abot_,1))
        
        call IntDip_qqb_Zprime_ttb(MomExt(1:4,1:4), z, IDip) ! Normalization: as/(2 Pi)

        IDip = IDip * (alpha_s*RunFactor)/(2d0*Pi) * PreFac

        IDipAmp = IDip(1) * (PDFFac_dip_L(1) * LO_Res_UnPol_Left + PDFFac_dip_R(1) * LO_Res_UnPol_Right)
        IDipAmp = IDipAmp + IDip(2)/z * (PDFFac_dip_L(2) * LO_Res_UnPol_Left + PDFFac_dip_R(2) * LO_Res_UnPol_Right)
        IDipAmp = IDipAmp + IDip(3)/z * (PDFFac_dip_L(3) * LO_Res_UnPol_Left + PDFFac_dip_R(3) * LO_Res_UnPol_Right)
        
        IDipAmp = IDipAmp * ISFac * 9d0 

        EvalCS_1L_Zprime_ttbqqb = EvalCS_1L_Zprime_ttbqqb + IDipAmp

     ENDIF
     

  ENDIF

!--F Oct 13: poles cancellation checks
!  print *, 'full result', EvalCS_1L_Zprime_ttbqqb
!  pause

  do NHisto=1,NumHistograms
     call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_Zprime_ttbqqb)
  enddo

  EvalCS_1L_Zprime_ttbqqb = EvalCS_1L_Zprime_ttbqqb/VgsWgt

  return

END FUNCTION EvalCS_1L_Zprime_ttbqqb



FUNCTION EvalCS_Real_Zprime_ttbqqb(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes_Zprime
use ModParameters
use ModDipoles_Zprime
implicit none
real(8) ::  EvalCS_Real_Zprime_ttbqqb,yRnd(1:VegasMxDim),VgsWgt,DipoleResult,EvalCS_Dips_Zprime_ttbqqb
complex(8) :: Res_Pol_L, Res_Pol_R, prim1_L, prim1_R, prim2_L, prim2_R
real(8) :: Res_UnPol_L, Res_UnPol_R
integer :: iHel
real(8) :: EHat,RunFactor,PSWgt,PS,PSWgt2,PSWgt3,xFrag, sHatJacobi
real(8) :: FluxFac, PreFac, ISFac
real(8) :: pdf(-6:6,1:2),eta1,eta2, PDFFac_L, PDFFac_R
real(8) :: MomExt(1:4,1:5),MomDK(1:4,1:6), MomExtTd(1:4,1:4), MomDKTd(1:4,1:6)
integer :: NBin(1:NumMaxHisto),NHisto
logical :: applyPSCut,applySingCut
real(8) :: colf
real(8),parameter :: ThresholdCutOff = 1.00d0
real(8),parameter :: Cf = 4d0/3d0
integer :: nDip
real(8) :: resdip, dipoles
real(8) :: resLO_L, resLO_R
complex(8) :: ampLO_L, ampLO_R
real(8) :: xpdf
include "vegas_common.f"
integer :: i1,i2,i3,i4,i5


  EvalCS_Real_Zprime_ttbqqb = 0d0
  EvalCS_Dips_Zprime_ttbqqb = 0d0


!  yrnd(1:15) = 0.4d0; yrnd(8)=0.6d0; yrnd(16)=0.6d0; print *, "fixed y,y16"
!  yRnd(1:2) = (/0.1d0,0.3d0/); print *, 'fixed yRnd(1:2)'

  if ((TOPDECAYS.EQ.0 .AND. yRnd(8).GT.0.5d0) .OR. (TOPDECAYS.NE.0 .AND. yRnd(16).GT.0.5d0)) then
     call PDFMapping(62,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi) ! BW mapping for Ehat
  else
     call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi) ! no mapping 
  endif

  if( EHat.le.2d0*m_Top * ThresholdCutOff ) then
      EvalCS_Real_Zprime_ttbqqb = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

  if ((TOPDECAYS.EQ.0 .AND. yRnd(8).LE.0.5d0) .OR. (TOPDECAYS.NE.0 .AND. yRnd(16).LE.0.5d0)) then
     call EvalPhasespaceBWMapp(EHat,(/m_Top,m_Top,0d0/),yRnd(3:7),MomExt(1:4,1:5),PSWgt)!   BW mapping for M_ttbar
     call swapmom(MomExt(1:4,3),MomExt(1:4,5))
  else
     call EvalPhaseSpace_2to3(EHat,yRnd(3:7),MomExt(1:4,1:5),PSWgt) ! no mapping
  endif

  call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))

  call CheckSing(MomExt,applySingCut)
  if( applySingCut ) then
     EvalCS_Real_Zprime_ttbqqb = 0d0
     SkipCounter = SkipCounter + 1
     return
  endif
  
  IF( TOPDECAYS.NE.0 ) THEN
     call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomDK(1:4,1:3),PSWgt2) ! anti-top
     call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:15),.false.,MomDK(1:4,4:6),PSWgt3)
     PSWgt = PSWgt * PSWgt2*PSWgt3
  ENDIF

  call Kinematics_TTBARZprime(.true.,MomExt,MomDK,applyPSCut,NBin)

! applyPSCut=.true.! this is for inverted alpha check

  call setPDFs(eta1,eta2,MuFac,pdf)


  IF (PROCESS.EQ.63 ) THEN
     ! gluon index is 0

     call random_number(xpdf)

     if( xpdf.le.0.5d0 ) then !   here was a bug: 1<-->2; corrected now
        PDFFac_L =            gL_Zpr(up_)**2 * (pdf(0,2)*pdf(up_,1))
        PDFFac_L = PDFFac_L + gL_Zpr(dn_)**2 * (pdf(0,2)*pdf(dn_,1))
        PDFFac_L = PDFFac_L + gL_Zpr(chm_)**2 * (pdf(0,2)*pdf(chm_,1))
        PDFFac_L = PDFFac_L + gL_Zpr(str_)**2 * (pdf(0,2)*pdf(str_,1))
        PDFFac_L = PDFFac_L + gL_Zpr(bot_)**2 * (pdf(0,2)*pdf(bot_,1))
        PDFFac_R =            gR_Zpr(up_)**2 * (pdf(0,2)*pdf(up_,1))
        PDFFac_R = PDFFac_R + gR_Zpr(dn_)**2 * (pdf(0,2)*pdf(dn_,1))
        PDFFac_R = PDFFac_R + gR_Zpr(chm_)**2 * (pdf(0,2)*pdf(chm_,1))
        PDFFac_R = PDFFac_R + gR_Zpr(str_)**2 * (pdf(0,2)*pdf(str_,1))
        PDFFac_R = PDFFac_R + gR_Zpr(bot_)**2 * (pdf(0,2)*pdf(bot_,1))
     else
        PDFFac_L =            gL_Zpr(up_)**2 * (pdf(0,1)*pdf(up_,2))
        PDFFac_L = PDFFac_L + gL_Zpr(dn_)**2 * (pdf(0,1)*pdf(dn_,2))
        PDFFac_L = PDFFac_L + gL_Zpr(chm_)**2 * (pdf(0,1)*pdf(chm_,2))
        PDFFac_L = PDFFac_L + gL_Zpr(str_)**2 * (pdf(0,1)*pdf(str_,2))
        PDFFac_L = PDFFac_L + gL_Zpr(bot_)**2 * (pdf(0,1)*pdf(bot_,2))
        PDFFac_R =            gR_Zpr(up_)**2 * (pdf(0,1)*pdf(up_,2))
        PDFFac_R = PDFFac_R + gR_Zpr(dn_)**2 * (pdf(0,1)*pdf(dn_,2))
        PDFFac_R = PDFFac_R + gR_Zpr(chm_)**2 * (pdf(0,1)*pdf(chm_,2))
        PDFFac_R = PDFFac_R + gR_Zpr(str_)**2 * (pdf(0,1)*pdf(str_,2))
        PDFFac_R = PDFFac_R + gR_Zpr(bot_)**2 * (pdf(0,1)*pdf(bot_,2))
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
     endif

     PDFFac_R = PDFFac_R * 2d0!  factor two for sampling over pdfs
     PDFFac_L = PDFFac_L * 2d0

!     PDFFac_L = gL_Zpr(up_)**2 * (pdf(0,1)*pdf(up_,2) + pdf(0,2)*pdf(up_,1))
!     PDFFac_L = PDFFac_L + gL_Zpr(dn_)**2 * (pdf(0,1)*pdf(dn_,2) + pdf(0,2)*pdf(dn_,1))
!     PDFFac_L = PDFFac_L + gL_Zpr(chm_)**2 * (pdf(0,1)*pdf(chm_,2) + pdf(0,2)*pdf(chm_,1))
!     PDFFac_L = PDFFac_L + gL_Zpr(str_)**2 * (pdf(0,1)*pdf(str_,2) + pdf(0,2)*pdf(str_,1))
!     PDFFac_L = PDFFac_L + gL_Zpr(bot_)**2 * (pdf(0,1)*pdf(bot_,2) + pdf(0,2)*pdf(bot_,1))

!     PDFFac_R = gR_Zpr(up_)**2 * (pdf(0,1)*pdf(up_,2) + pdf(0,2)*pdf(up_,1))
!     PDFFac_R = PDFFac_R + gR_Zpr(dn_)**2 * (pdf(0,1)*pdf(dn_,2) + pdf(0,2)*pdf(dn_,1))
!     PDFFac_R = PDFFac_R + gR_Zpr(chm_)**2 * (pdf(0,1)*pdf(chm_,2) + pdf(0,2)*pdf(chm_,1))
!     PDFFac_R = PDFFac_R + gR_Zpr(str_)**2 * (pdf(0,1)*pdf(str_,2) + pdf(0,2)*pdf(str_,1))
!     PDFFac_R = PDFFac_R + gR_Zpr(bot_)**2 * (pdf(0,1)*pdf(bot_,2) + pdf(0,2)*pdf(bot_,1))

  ELSEIF (PROCESS.EQ.64 ) THEN

     call random_number(xpdf)

    if( xpdf.le.0.5d0 ) then
        PDFFac_L = gL_Zpr(up_)**2 * (pdf(0,1)*pdf(aup_,2))
        PDFFac_L = PDFFac_L + gL_Zpr(dn_)**2 * (pdf(0,1)*pdf(adn_,2))
        PDFFac_L = PDFFac_L + gL_Zpr(chm_)**2 * (pdf(0,1)*pdf(achm_,2))
        PDFFac_L = PDFFac_L + gL_Zpr(str_)**2 * (pdf(0,1)*pdf(astr_,2))
        PDFFac_L = PDFFac_L + gL_Zpr(bot_)**2 * (pdf(0,1)*pdf(abot_,2))
        PDFFac_R = gR_Zpr(up_)**2 * (pdf(0,1)*pdf(aup_,2))
        PDFFac_R = PDFFac_R + gR_Zpr(dn_)**2 * (pdf(0,1)*pdf(adn_,2))
        PDFFac_R = PDFFac_R + gR_Zpr(chm_)**2 * (pdf(0,1)*pdf(achm_,2))
        PDFFac_R = PDFFac_R + gR_Zpr(str_)**2 * (pdf(0,1)*pdf(astr_,2))
        PDFFac_R = PDFFac_R + gR_Zpr(bot_)**2 * (pdf(0,1)*pdf(abot_,2))
     else
        PDFFac_L = gL_Zpr(up_)**2 * (pdf(0,2)*pdf(aup_,1))
        PDFFac_L = PDFFac_L + gL_Zpr(dn_)**2 * (pdf(0,2)*pdf(adn_,1))
        PDFFac_L = PDFFac_L + gL_Zpr(chm_)**2 * (pdf(0,2)*pdf(achm_,1))
        PDFFac_L = PDFFac_L + gL_Zpr(str_)**2 * (pdf(0,2)*pdf(astr_,1))
        PDFFac_L = PDFFac_L + gL_Zpr(bot_)**2 * (pdf(0,2)*pdf(abot_,1))
        PDFFac_R = gR_Zpr(up_)**2 * (pdf(0,2)*pdf(aup_,1))
        PDFFac_R = PDFFac_R + gR_Zpr(dn_)**2 * (pdf(0,2)*pdf(adn_,1))
        PDFFac_R = PDFFac_R + gR_Zpr(chm_)**2 * (pdf(0,2)*pdf(achm_,1))
        PDFFac_R = PDFFac_R + gR_Zpr(str_)**2 * (pdf(0,2)*pdf(astr_,1))
        PDFFac_R = PDFFac_R + gR_Zpr(bot_)**2 * (pdf(0,2)*pdf(abot_,1))
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
     endif
     
     PDFFac_R = PDFFac_R * 2d0!  factor two for sampling over pdfs
     PDFFac_L = PDFFac_L * 2d0

!     PDFFac_L = gL_Zpr(up_)**2 * (pdf(0,1)*pdf(aup_,2) + pdf(0,2)*pdf(aup_,1))
!     PDFFac_L = PDFFac_L + gL_Zpr(dn_)**2 * (pdf(0,1)*pdf(adn_,2) + pdf(0,2)*pdf(adn_,1))
!     PDFFac_L = PDFFac_L + gL_Zpr(chm_)**2 * (pdf(0,1)*pdf(achm_,2) + pdf(0,2)*pdf(achm_,1))
!     PDFFac_L = PDFFac_L + gL_Zpr(str_)**2 * (pdf(0,1)*pdf(astr_,2) + pdf(0,2)*pdf(astr_,1))
!     PDFFac_L = PDFFac_L + gL_Zpr(bot_)**2 * (pdf(0,1)*pdf(abot_,2) + pdf(0,2)*pdf(abot_,1))

!     PDFFac_R = gR_Zpr(up_)**2 * (pdf(0,1)*pdf(aup_,2) + pdf(0,2)*pdf(aup_,1))
!     PDFFac_R = PDFFac_R + gR_Zpr(dn_)**2 * (pdf(0,1)*pdf(adn_,2) + pdf(0,2)*pdf(adn_,1))
!     PDFFac_R = PDFFac_R + gR_Zpr(chm_)**2 * (pdf(0,1)*pdf(achm_,2) + pdf(0,2)*pdf(achm_,1))
!     PDFFac_R = PDFFac_R + gR_Zpr(str_)**2 * (pdf(0,1)*pdf(astr_,2) + pdf(0,2)*pdf(astr_,1))
!     PDFFac_R = PDFFac_R + gR_Zpr(bot_)**2 * (pdf(0,1)*pdf(abot_,2) + pdf(0,2)*pdf(abot_,1))

  ELSEIF( PROCESS.EQ.66 ) THEN

    call random_number(xpdf)

    if( xpdf.le.0.5d0 ) then
       PDFFac_L =            gL_Zpr(up_)**2 * (pdf(up_,1)*pdf(aup_,2))
       PDFFac_L = PDFFac_L + gL_Zpr(dn_)**2 * (pdf(dn_,1)*pdf(adn_,2))
       PDFFac_L = PDFFac_L + gL_Zpr(chm_)**2 * (pdf(chm_,1)*pdf(achm_,2))
       PDFFac_L = PDFFac_L + gL_Zpr(str_)**2 * (pdf(str_,1)*pdf(astr_,2))
       PDFFac_L = PDFFac_L + gL_Zpr(bot_)**2 * (pdf(bot_,1)*pdf(abot_,2))
       
       PDFFac_R =            gR_Zpr(up_)**2 * (pdf(up_,1)*pdf(aup_,2))
       PDFFac_R = PDFFac_R + gR_Zpr(dn_)**2 * (pdf(dn_,1)*pdf(adn_,2))
       PDFFac_R = PDFFac_R + gR_Zpr(chm_)**2 * (pdf(chm_,1)*pdf(achm_,2))
       PDFFac_R = PDFFac_R + gR_Zpr(str_)**2 * (pdf(str_,1)*pdf(astr_,2))
       PDFFac_R = PDFFac_R + gR_Zpr(bot_)**2 * (pdf(bot_,1)*pdf(abot_,2))
    else
       PDFFac_L =            gL_Zpr(up_)**2 * (pdf(Up_,2)*pdf(aup_,1))
       PDFFac_L = PDFFac_L + gL_Zpr(dn_)**2 * (pdf(dn_,2)*pdf(adn_,1))
       PDFFac_L = PDFFac_L + gL_Zpr(chm_)**2 * (pdf(chm_,2)*pdf(achm_,1))
       PDFFac_L = PDFFac_L + gL_Zpr(str_)**2 * (pdf(str_,2)*pdf(astr_,1))
       PDFFac_L = PDFFac_L + gL_Zpr(bot_)**2 * (pdf(bot_,2)*pdf(abot_,1))
       
       PDFFac_R =            gR_Zpr(up_)**2 * (pdf(Up_,2)*pdf(aup_,1))
       PDFFac_R = PDFFac_R + gR_Zpr(dn_)**2 * (pdf(dn_,2)*pdf(adn_,1))
       PDFFac_R = PDFFac_R + gR_Zpr(chm_)**2 * (pdf(chm_,2)*pdf(achm_,1))
       PDFFac_R = PDFFac_R + gR_Zpr(str_)**2 * (pdf(str_,2)*pdf(astr_,1))
       PDFFac_R = PDFFac_R + gR_Zpr(bot_)**2 * (pdf(bot_,2)*pdf(abot_,1))
       call swapMom(MomExt(1:4,1),MomExt(1:4,2))
    endif
  
  PDFFac_R = PDFFac_R * 2d0!  factor two for sampling over pdfs
  PDFFac_L = PDFFac_L * 2d0

!     PDFFac_L = gL_Zpr(up_)**2 * (pdf(up_,1)*pdf(aup_,2) + pdf(up_,2)*pdf(aup_,1))
!     PDFFac_L = PDFFac_L + gL_Zpr(dn_)**2 * (pdf(dn_,1)*pdf(adn_,2) + pdf(dn_,2)*pdf(adn_,1))
!     PDFFac_L = PDFFac_L + gL_Zpr(chm_)**2 * (pdf(chm_,1)*pdf(achm_,2) + pdf(chm_,2)*pdf(achm_,1))
!     PDFFac_L = PDFFac_L + gL_Zpr(str_)**2 * (pdf(str_,1)*pdf(astr_,2) + pdf(str_,2)*pdf(astr_,1))
!     PDFFac_L = PDFFac_L + gL_Zpr(bot_)**2 * (pdf(bot_,1)*pdf(abot_,2) + pdf(bot_,2)*pdf(abot_,1))

!     PDFFac_R = gR_Zpr(up_)**2 * (pdf(up_,1)*pdf(aup_,2) + pdf(up_,2)*pdf(aup_,1))
!     PDFFac_R = PDFFac_R + gR_Zpr(dn_)**2 * (pdf(dn_,1)*pdf(adn_,2) + pdf(dn_,2)*pdf(adn_,1))
!     PDFFac_R = PDFFac_R + gR_Zpr(chm_)**2 * (pdf(chm_,1)*pdf(achm_,2) + pdf(chm_,2)*pdf(achm_,1))
!     PDFFac_R = PDFFac_R + gR_Zpr(str_)**2 * (pdf(str_,1)*pdf(astr_,2) + pdf(str_,2)*pdf(astr_,1))
!     PDFFac_R = PDFFac_R + gR_Zpr(bot_)**2 * (pdf(bot_,1)*pdf(abot_,2) + pdf(bot_,2)*pdf(abot_,1))
  ENDIF

  PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
  colf = 2d0 * 9d0 * Cf

  RunFactor = RunAlphaS(NLOParam,MuRen)
  ISFac = MomCrossing(MomExt)

  IF( TOPDECAYS.GE.1 ) THEN
     call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
     call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
  ENDIF

  if( applyPSCut ) then
     EvalCS_Real_Zprime_ttbqqb = 0d0
  else

     ISFac = MomCrossing(MomExt)

     Res_UnPol_L = 0d0
     Res_UnPol_R = 0d0
     do iHel=1,NumHelicities
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
        call SetPolarizations()

        call Tree_Zprime_tbtqbqg_i(prim1_L, prim1_R)
        call Tree_Zprime_tbtqbqg_f(prim2_L, prim2_R)

        Res_UnPol_L = Res_UnPol_L +  dreal(prim1_L*dconjg(prim1_L) + prim2_L*dconjg(prim2_L))
        Res_UnPol_R = Res_UnPol_R +  dreal(prim1_R*dconjg(prim1_R) + prim2_R*dconjg(prim2_R))

     enddo!helicity loop

     EvalCS_Real_Zprime_ttbqqb = (PDFFac_R * Res_UnPol_R + PDFFac_L * Res_UnPol_L)
     EvalCS_Real_Zprime_ttbqqb = EvalCS_Real_Zprime_ttbqqb * (4d0*Pi*alpha_s*RunFactor) * PreFac * ISFac * colf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! F- MadGraph check, Oct 13th
!
!     open(unit=66,file='/home/caola/tools/MG_ME_V4.5.2/zprime/SubProcesses/kin.dat', status='unknown')
!     do iHel = 1,5
!        write(66,*) MomExt(1,iHel) * 100d0
!        write(66,*) MomExt(2,iHel) * 100d0
!        write(66,*) MomExt(3,iHel) * 100d0
!        write(66,*) MomExt(4,iHel) * 100d0
!     enddo
!     close(66)
!     print *, 'MZpr', m_zpr * 100d0
!     print *, 'Ga_zpr', Ga_Zpr * 100d0
!     print *, 'gL_Zpr(up_), gL_Zpr(dn_)', gL_Zpr(up_), gL_Zpr(dn_)
!     print *, 'gR_Zpr(up_), gR_Zpr(dn_)', gR_Zpr(up_), gR_Zpr(dn_)
!     print *, 'alpha_s', alpha_s*RunFactor
!     print *, 'gs', dsqrt(alpha_s*RunFactor*4d0*Pi)
!     print *, ''
!     print *, 'up',  ISFac * colf * (4d0*Pi*alpha_s*RunFactor) * (gL_Zpr(up_)**2 * Res_UnPol_L + gR_Zpr(up_)**2 * Res_UnPol_R)*GeV**2
!     print *, 'dn',  ISFac * colf * (4d0*Pi*alpha_s*RunFactor) * (gL_Zpr(dn_)**2 * Res_UnPol_L + gR_Zpr(dn_)**2 * Res_UnPol_R)*GeV**2
!     pause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     do NHisto=1,NumHistograms
        call intoHisto(NHisto,NBin(NHisto),EvalCS_Real_Zprime_ttbqqb)
     enddo
     EvalCounter = EvalCounter + 1

  endif!applyPSCut





  !!! Dipoles section
  IF ( PROCESS .EQ. 66 ) THEN ! qqb
     dipoles = 0d0
     do nDip = 1,4

        call Dipoles_qqb_Zprime_ttb(nDip, MomExt(1:4,1:5), MomExtTd(1:4,1:4), resdip)! here was a bug: MomExtTd was running from 1:5; correct now
 
        if (resdip .EQ. 0d0) cycle

        Crossing(:) = (/3,4,-1,-2,0/) ! For dipoles
        NumExtParticles = 4
        ISFac = MomCrossing(MomExtTd)
        Crossing(:) = (/4,5,-1,-2,3/) ! Original: this was wrong, corrected
        NumExtParticles = 5
        
        IF( TOPDECAYS.NE.0 ) THEN
           call EvalPhasespace_TopDecay(MomExtTd(1:4,3),yRnd(8:11),.false.,MomDKTd(1:4,1:3),PSWgt2)
           call EvalPhasespace_TopDecay(MomExtTd(1:4,4),yRnd(12:15),.false.,MomDKTd(1:4,4:6),PSWgt3)
           PSWgt = PSWgt * PSWgt2*PSWgt3! this is better never be used anywhere below
           
           call Kinematics_TTBARZprime(.false.,MomExtTd,MomDKTd,applyPSCut,NBin)

           if ( applyPSCut ) then
              resdip = 0d0
           else
              
              call TopDecay(ExtParticle(1),DK_LO,MomDKTd(1:4,1:3))
              call TopDecay(ExtParticle(2),DK_LO,MomDKTd(1:4,4:6))
              
              resLO_L = 0d0
              resLO_R = 0d0
              
              i1 = 0
              i2 = 0
              do i3 = -1,1,2
                 ExtParticle(1)%Helicity = i1
                 ExtParticle(2)%Helicity = i2
                 ExtParticle(3)%Helicity = i3
                 ExtParticle(4)%Helicity = -i3
                 call SetPolarizations()
                 call Tree_Zprime_tbtqbq(ampLO_L, ampLO_R)
                 resLO_L = resLO_L +  dreal(ampLO_L * dconjg(ampLO_L))
                 resLO_R = resLO_R +  dreal(ampLO_R * dconjg(ampLO_R))
              enddo
              
           endif
           
        ELSE
           
           call Kinematics_TTBARZprime(.false.,MomExtTd,MomDK,applyPSCut,NBin)
           
           if ( applyPSCut ) then
              resdip = 0d0
           else
              
              resLO_L = 0d0
              resLO_R = 0d0
              
              do i1 = -1,1,2
                 do i2 = -1,1,2
                    do i3 = -1,1,2
                       ExtParticle(1)%Helicity = i1
                       ExtParticle(2)%Helicity = i2
                       ExtParticle(3)%Helicity = i3
                       ExtParticle(4)%Helicity = -i3
                       call SetPolarizations()
                       call Tree_Zprime_tbtqbq(ampLO_L, ampLO_R)
                       resLO_L = resLO_L +  dreal(ampLO_L * dconjg(ampLO_L))
                       resLO_R = resLO_R +  dreal(ampLO_R * dconjg(ampLO_R))
                    enddo
                 enddo
              enddo
              
           endif
           
        ENDIF !TopDecays
        
        ! Check against Markus
        resdip = resdip * (PDFFac_R * resLO_R + PDFFac_L * resLO_L)
        resdip = resdip * (4d0*Pi*alpha_s*RunFactor) * PreFac * ISFac * colf
        
        do NHisto=1,NumHistograms
           call intoHisto(NHisto,NBin(NHisto),resdip)
        enddo
     
        dipoles = dipoles + resdip
     
     enddo
  

  ELSEIF ( PROCESS .EQ. 64) THEN
     dipoles = 0d0

     do nDip = 1, 1
   
        call Dipoles_gqb_Zprime_ttbqb(1, MomExt(1:4,1:5), MomExtTd(1:4,1:4), resdip)

        if (resdip .eq. 0d0) cycle
     
        Crossing(:) = (/3,4,-1,-2,0/) ! For dipoles
        NumExtParticles = 4
        ISFac = MomCrossing(MomExtTd)
        Crossing(:) = (/4,5,3,-2,-1/) ! Original: this was wrong, corrected
        NumExtParticles = 5
        
        IF( TOPDECAYS.NE.0 ) THEN
           call EvalPhasespace_TopDecay(MomExtTd(1:4,3),yRnd(8:11),.false.,MomDKTd(1:4,1:3),PSWgt2)
           call EvalPhasespace_TopDecay(MomExtTd(1:4,4),yRnd(12:15),.false.,MomDKTd(1:4,4:6),PSWgt3)
           PSWgt = PSWgt * PSWgt2*PSWgt3
           
           call Kinematics_TTBARZprime(.false.,MomExtTd,MomDKTd,applyPSCut,NBin)
           
           if ( applyPSCut) then
              resdip = 0d0
           else
              
              call TopDecay(ExtParticle(1),DK_LO,MomDKTd(1:4,1:3))
              call TopDecay(ExtParticle(2),DK_LO,MomDKTd(1:4,4:6))
              
              resLO_L = 0d0
              resLO_R = 0d0
              
              i1 = 0
              i2 = 0
              
              do i3 = -1,1,2
                 ExtParticle(1)%Helicity = i1
                 ExtParticle(2)%Helicity = i2
                 ExtParticle(3)%Helicity = i3
                 ExtParticle(4)%Helicity = -i3
                 call SetPolarizations()
                 call Tree_Zprime_tbtqbq(ampLO_L, ampLO_R)
                 resLO_L = resLO_L +  dreal(ampLO_L * dconjg(ampLO_L))
                 resLO_R = resLO_R +  dreal(ampLO_R * dconjg(ampLO_R))
              enddo
              
           endif
           
        ELSE

           call Kinematics_TTBARZprime(.false.,MomExtTd,MomDK,applyPSCut,NBin)
           
           if ( applyPSCut ) then
              resdip = 0d0
           else
              
              resLO_L = 0d0
              resLO_R = 0d0
              
              do i1 = -1,1,2
                 do i2 = -1,1,2
                    do i3 = -1,1,2
                       ExtParticle(1)%Helicity = i1
                       ExtParticle(2)%Helicity = i2
                       ExtParticle(3)%Helicity = i3
                       ExtParticle(4)%Helicity = -i3
                       call SetPolarizations()
                       call Tree_Zprime_tbtqbq(ampLO_L, ampLO_R)
                       resLO_L = resLO_L +  dreal(ampLO_L * dconjg(ampLO_L))
                       resLO_R = resLO_R +  dreal(ampLO_R * dconjg(ampLO_R))
                    enddo
                 enddo
              enddo
              
           endif
        
        ENDIF !TopDecays
        
        ! Check against Markus
        resdip = resdip * (PDFFac_R * resLO_R + PDFFac_L * resLO_L)
        resdip = resdip * (4d0*Pi*alpha_s*RunFactor) * PreFac * ISFac * colf
     
        do NHisto=1,NumHistograms
           call intoHisto(NHisto,NBin(NHisto),resdip)
        enddo
     
        dipoles = dipoles + resdip

     enddo
     
  ELSEIF ( PROCESS .EQ. 63) THEN
     dipoles = 0d0
   
     do nDip = 1,1

        call Dipoles_qg_Zprime_ttbq(1, MomExt(1:4,1:5), MomExtTd(1:4,1:4), resdip)
     
        if (resdip .EQ. 0d0 ) cycle

        Crossing(:) = (/3,4,-1,-2,0/) ! For dipoles
        NumExtParticles=4
        ISFac = MomCrossing(MomExtTd)
        Crossing(:) = (/4,5,-1,3,-2/) ! Original: this was wrong, corrected
        NumExtParticles=5
        
        IF( TOPDECAYS.NE.0 ) THEN
           call EvalPhasespace_TopDecay(MomExtTd(1:4,3),yRnd(8:11),.false.,MomDKTd(1:4,1:3),PSWgt2)
           call EvalPhasespace_TopDecay(MomExtTd(1:4,4),yRnd(12:15),.false.,MomDKTd(1:4,4:6),PSWgt3)
           PSWgt = PSWgt * PSWgt2*PSWgt3
           
           call Kinematics_TTBARZprime(.false.,MomExtTd,MomDKTd,applyPSCut,NBin)
           
           if ( applyPSCut) then
              resdip = 0d0
           else
              
              call TopDecay(ExtParticle(1),DK_LO,MomDKTd(1:4,1:3))
              call TopDecay(ExtParticle(2),DK_LO,MomDKTd(1:4,4:6))
              
              resLO_L = 0d0
              resLO_R = 0d0
              
              i1 = 0
              i2 = 0
              
              do i3 = -1,1,2
                 ExtParticle(1)%Helicity = i1
                 ExtParticle(2)%Helicity = i2
                 ExtParticle(3)%Helicity = i3
                 ExtParticle(4)%Helicity = -i3
                 call SetPolarizations()
                 call Tree_Zprime_tbtqbq(ampLO_L, ampLO_R)
                 resLO_L = resLO_L +  dreal(ampLO_L * dconjg(ampLO_L))
                 resLO_R = resLO_R +  dreal(ampLO_R * dconjg(ampLO_R))
              enddo
           
           endif
        
        ELSE

           call Kinematics_TTBARZprime(.false.,MomExtTd,MomDK,applyPSCut,NBin)
        
           if ( applyPSCut ) then
              resdip = 0d0
           else
              
              resLO_L = 0d0
              resLO_R = 0d0
              
              do i1 = -1,1,2
                 do i2 = -1,1,2
                    do i3 = -1,1,2
                       ExtParticle(1)%Helicity = i1
                       ExtParticle(2)%Helicity = i2
                       ExtParticle(3)%Helicity = i3
                       ExtParticle(4)%Helicity = -i3
                       call SetPolarizations()
                       call Tree_Zprime_tbtqbq(ampLO_L, ampLO_R)
                       resLO_L = resLO_L +  dreal(ampLO_L * dconjg(ampLO_L))
                       resLO_R = resLO_R +  dreal(ampLO_R * dconjg(ampLO_R))
                    enddo
                 enddo
              enddo
              
           endif
           
        ENDIF !TopDecays
        
        ! Check against Markus
        resdip = resdip * (PDFFac_R * resLO_R + PDFFac_L * resLO_L)
        resdip = resdip * (4d0*Pi*alpha_s*RunFactor) * PreFac * ISFac * colf
        
        do NHisto=1,NumHistograms
           call intoHisto(NHisto,NBin(NHisto),resdip)
        enddo
        
        dipoles = dipoles + resdip
        
     enddo
    
  ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--F Oct 13: dipoles check
!  print *, 'inv    ', (MomExt(1:4,1).dot.MomExt(1:4,3))/Ehat**2,(MomExt(1:4,2).dot.MomExt(1:4,3))/Ehat**2
!  print *, 'real   ', EvalCS_Real_Zprime_ttbqqb
!  print *, 'dipoles', dipoles
!  print *, 'dipoles/real+1', dipoles/EvalCS_Real_Zprime_ttbqqb + 1d0
!  pause 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  EvalCS_Real_Zprime_ttbqqb = (EvalCS_Real_Zprime_ttbqqb+dipoles)/VgsWgt

  RETURN

END FUNCTION EvalCS_Real_Zprime_ttbqqb



FUNCTION EvalCS_Virt_Zprime_Interf(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes_Zprime
use ModParameters
use ModJPsiFrag
use ModIntDipoles_ZprimeTTB
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
implicit none
real(8), parameter :: Nc = 3d0
real(8) ::  EvalCS_Virt_Zprime_Interf,yRnd(1:VegasMxDim),VgsWgt
integer :: iPrimAmp
complex(8) :: LO_Res_Pol_Left, LO_Res_Pol_Right, LO_Res_Pol_glue
complex(8) :: Virt_Res_Pol_Left, Virt_Res_Pol_Right
complex(8) :: BosonicPartAmp(4,-2:1), FermionPartAmp(2,-2:1), NLO_Res_Pol(2)
real(8) :: Virt_Res_Unpol_Left, Virt_Res_Unpol_Right, NLO_Res_Unpol_Left, NLO_Res_Unpol_Right
integer :: iHel, sig_t, sig_tb
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,ISFac,xFrag
real(8) :: MomExt(1:4,1:NumExtParticles),MomDK(1:4,1:6)
logical :: applyPSCut
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac_L, PDFFac_R, PDFFac_dip_L(3), PDFFac_dip_R(3)
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2), z,AccPoles
real(8) :: IDip(3), IDipAmp, LO_for_dip_Left, LO_for_dip_Right
integer :: NHisto,NBin(1:NumMaxHisto),npdf,nHel(1:2),NRndHel
real(8) :: xpdf,beta,AndreasResult
include 'misc/global_import'
include "vegas_common.f"

complex(8) :: rdiv(1:2)

  EvalCS_Virt_Zprime_Interf = 0d0

  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi) ! Flat mapping

  if( EHat.le.2d0*m_Top) then
     EvalCS_Virt_Zprime_Interf = 0d0
     return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

  call EvalPhaseSpace_2to2(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt)
  call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
!print *, "MaRKUS: boost off for checks"

  NRndHel=5
  IF( TOPDECAYS.NE.0 ) THEN
     call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:8),.false.,MomDK(1:4,1:3),PSWgt2)
     call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(9:12),.false.,MomDK(1:4,4:6),PSWgt3)
     PSWgt = PSWgt * PSWgt2*PSWgt3
     NRndHel=13
  ENDIF


  call Kinematics_TTBARZprime(.false.,MomExt,MomDK,applyPSCut,NBin)
  if( applyPSCut ) then
     EvalCS_Virt_Zprime_Interf = 0d0
     return
  endif


  ! Override masterprocess Helicities to avoid spin-flip symmetry

  IF( TOPDECAYS.GE.1 ) THEN
     NumHelicities = 2
     deallocate(Helicities)
     allocate(Helicities(1:NumHelicities,1:NumExtParticles))
     Helicities(1,1:4) = (/0,0,-1,+1/)
     Helicities(2,1:4) = (/0,0,+1,-1/)
  ELSE
     NumHelicities = 8
     deallocate(Helicities)
     allocate(Helicities(1:NumHelicities,1:NumExtParticles))
     sig_tb=+1; sig_t =+1;
     Helicities(1,1:NumExtParticles) = (/sig_tb,sig_t,+1,-1/)
     Helicities(2,1:NumExtParticles) = (/sig_tb,sig_t,-1,+1/)
     sig_tb=+1; sig_t =-1;
     Helicities(3,1:NumExtParticles) = (/sig_tb,sig_t,+1,-1/)
     Helicities(4,1:NumExtParticles) = (/sig_tb,sig_t,-1,+1/)
     sig_tb=-1; sig_t =+1;
     Helicities(5,1:NumExtParticles) = (/sig_tb,sig_t,+1,-1/)
     Helicities(6,1:NumExtParticles) = (/sig_tb,sig_t,-1,+1/)
     sig_tb=-1; sig_t =-1;
     Helicities(7,1:NumExtParticles) = (/sig_tb,sig_t,+1,-1/)
     Helicities(8,1:NumExtParticles) = (/sig_tb,sig_t,-1,+1/)
  ENDIF


  call setPDFs(eta1,eta2,MuFac,pdf)

  call random_number(xpdf)

!  xpdf = 0.4d0; print *, 'fixed xpdf'

  if( xpdf.le.0.5d0 ) then
     PDFFac_L = gL_Zpr(up_) * (pdf(up_,1)*pdf(aup_,2))
     PDFFac_L = PDFFac_L + gL_Zpr(dn_) * (pdf(dn_,1)*pdf(adn_,2))
     PDFFac_L = PDFFac_L + gL_Zpr(chm_) * (pdf(chm_,1)*pdf(achm_,2))
     PDFFac_L = PDFFac_L + gL_Zpr(str_) * (pdf(str_,1)*pdf(astr_,2))
     PDFFac_L = PDFFac_L + gL_Zpr(bot_) * (pdf(bot_,1)*pdf(abot_,2))
     
     PDFFac_R = gR_Zpr(up_) * (pdf(up_,1)*pdf(aup_,2))
     PDFFac_R = PDFFac_R + gR_Zpr(dn_) * (pdf(dn_,1)*pdf(adn_,2))
     PDFFac_R = PDFFac_R + gR_Zpr(chm_) * (pdf(chm_,1)*pdf(achm_,2))
     PDFFac_R = PDFFac_R + gR_Zpr(str_) * (pdf(str_,1)*pdf(astr_,2))
     PDFFac_R = PDFFac_R + gR_Zpr(bot_) * (pdf(bot_,1)*pdf(abot_,2))
  else
     PDFFac_L = gL_Zpr(up_) * (pdf(up_,2)*pdf(aup_,1))
     PDFFac_L = PDFFac_L + gL_Zpr(dn_) * (pdf(dn_,2)*pdf(adn_,1))
     PDFFac_L = PDFFac_L + gL_Zpr(chm_) * (pdf(chm_,2)*pdf(achm_,1))
     PDFFac_L = PDFFac_L + gL_Zpr(str_) * (pdf(str_,2)*pdf(astr_,1))
     PDFFac_L = PDFFac_L + gL_Zpr(bot_) * (pdf(bot_,2)*pdf(abot_,1))
     
     PDFFac_R = gR_Zpr(up_) * (pdf(up_,2)*pdf(aup_,1))
     PDFFac_R = PDFFac_R + gR_Zpr(dn_) * (pdf(dn_,2)*pdf(adn_,1))
     PDFFac_R = PDFFac_R + gR_Zpr(chm_) * (pdf(chm_,2)*pdf(achm_,1))
     PDFFac_R = PDFFac_R + gR_Zpr(str_) * (pdf(str_,2)*pdf(astr_,1))
     PDFFac_R = PDFFac_R + gR_Zpr(bot_) * (pdf(bot_,2)*pdf(abot_,1))
     call swapMom(MomExt(1:4,1),MomExt(1:4,2))
  endif
  PDFFac_R = PDFFac_R * 2d0!  factor two for sampling over pdfs
  PDFFac_L = PDFFac_L * 2d0

!  PDFFac_L = gL_Zpr(up_) * (pdf(up_,1)*pdf(aup_,2) + pdf(Up_,2)*pdf(aup_,1))
!  PDFFac_L = PDFFac_L + gL_Zpr(dn_) * (pdf(dn_,1)*pdf(adn_,2) + pdf(dn_,2)*pdf(adn_,1))
!  PDFFac_L = PDFFac_L + gL_Zpr(chm_) * (pdf(chm_,1)*pdf(achm_,2) + pdf(chm_,2)*pdf(achm_,1))
!  PDFFac_L = PDFFac_L + gL_Zpr(str_) * (pdf(str_,1)*pdf(astr_,2) + pdf(str_,2)*pdf(astr_,1))
!  PDFFac_L = PDFFac_L + gL_Zpr(bot_) * (pdf(bot_,1)*pdf(abot_,2) + pdf(bot_,2)*pdf(abot_,1))
  
!  PDFFac_R = gR_Zpr(up_) * (pdf(up_,1)*pdf(aup_,2) + pdf(Up_,2)*pdf(aup_,1))
!  PDFFac_R = PDFFac_R + gR_Zpr(dn_) * (pdf(dn_,1)*pdf(adn_,2) + pdf(dn_,2)*pdf(adn_,1))
!  PDFFac_R = PDFFac_R + gR_Zpr(chm_) * (pdf(chm_,1)*pdf(achm_,2) + pdf(chm_,2)*pdf(achm_,1))
!  PDFFac_R = PDFFac_R + gR_Zpr(str_) * (pdf(str_,1)*pdf(astr_,2) + pdf(str_,2)*pdf(astr_,1))
!  PDFFac_R = PDFFac_R + gR_Zpr(bot_) * (pdf(bot_,1)*pdf(abot_,2) + pdf(bot_,2)*pdf(abot_,1))

  
  PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
  RunFactor = RunAlphaS(NLOParam,MuRen)
  nHel(1:2) = getHelicity(yrnd(NRndHel))
  PreFac = PreFac * dble(NumHelicities/(nHel(2)-nHel(1)+1))
  

  Virt_Res_Unpol_Left = 0d0
  Virt_Res_Unpol_Right = 0d0
  NLO_Res_UnPol_Left = 0d0
  NLO_Res_UnPol_Right = 0d0

  LO_for_dip_left = 0d0
  LO_for_dip_right = 0d0
  
  ISFac = MomCrossing(MomExt)

  IF( TOPDECAYS.GE.1 ) THEN
     call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
     call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
  ENDIF
  
  IF(CORRECTION.EQ.1) THEN

     call InitCurrCache()
     call SetPropagators()


     do iHel=nHel(1),nHel(2)
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
        call SetPolarizations()
        
        ! Compute first the Zprime part
        call Tree_Zprime_tbtqbq(LO_Res_Pol_Left, LO_Res_Pol_Right)
        call Virt_Zprime_box(Virt_Res_Pol_Left, Virt_Res_Pol_Right)
        
        ! Now compute the gluon part
        do iPrimAmp=1,2
           call EvalTree(BornAmps(iPrimAmp)) ! Vertex: gs/sqrt(2)
        enddo
        
        ! ------------ bosonic loops --------------
        do iPrimAmp=1,2
           call SetKirill(PrimAmps(iPrimAmp))
           call PentCut(PrimAmps(iPrimAmp))
           call QuadCut(PrimAmps(iPrimAmp))
           call TripCut(PrimAmps(iPrimAmp))
           call DoubCut(PrimAmps(iPrimAmp))
           call SingCut(PrimAmps(iPrimAmp))
           call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
!           call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen**2)
           PrimAmps(iPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(iPrimAmp)%Result(-2:1)
           call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,2,rdiv(2),rdiv(1))
           AccPoles = CheckPoles(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv(1:2))
           !call WritePrimAmpResult(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv); print *, 'checking poles'
           if( AccPoles.gt.1d-4 ) then
!                 print *, "PrimAmp",iPrimAmp
!                 print *, "AccPoles",AccPoles
!                 call PrintYrnd(Yrnd(1:12))
!                 SkipCounter = SkipCounter + 1
!                 EvalCS_Virt_Zprime_Interf = 0d0
!                 return

!                 print *, "PrimAmp",iPrimAmp
!                 print *, "AccPoles",AccPoles
                 call PentCut_128(PrimAmps(iPrimAmp))
                 call QuadCut_128(PrimAmps(iPrimAmp))
                 call TripCut_128(PrimAmps(iPrimAmp))
                 call DoubCut_128(PrimAmps(iPrimAmp))
                 call SingCut_128(PrimAmps(iPrimAmp))
                 call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
                 PrimAmps(iPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(iPrimAmp)%Result(-2:1)
                 call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,2,rdiv(2),rdiv(1))
                 AccPoles = CheckPoles(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),rdiv(1:2))
                if( AccPoles.gt.1d-3 ) then
                   print *, "SKIP PrimAmp",iPrimAmp
                   print *, "AccPoles",AccPoles
                   call PrintYrnd(Yrnd(1:12))
                   SkipCounter = SkipCounter + 1
                   EvalCS_Virt_Zprime_Interf = 0d0
                   return
                endif 

           endif
        enddo


        BosonicPartAmp(1,-2:1) =  & ! Coefficient of d_tqbar d_qtbar
             + Nc * PrimAmps(PrimAmp1_1234)%Result(-2:1) &
             - 2d0/Nc*( PrimAmps(PrimAmp1_1234)%Result(-2:1) + PrimAmps(PrimAmp1_1243)%Result(-2:1) )
           
        BosonicPartAmp(2,-2:1) = PrimAmps(PrimAmp1_1243)%Result(-2:1)  & ! Coefficient of d_ttbar d_qqbar
             + 1d0/Nc**2 * ( PrimAmps(PrimAmp1_1234)%Result(-2:1) + PrimAmps(PrimAmp1_1243)%Result(-2:1) ) 
                   

        NLO_Res_Pol(1) = (BosonicPartAmp(1,0) + BosonicPartAmp(1,1)) ! cc + rational
        NLO_Res_Pol(2) = (BosonicPartAmp(2,0) + BosonicPartAmp(2,1)) ! cc + rational

        ! eps-cancellation check
!        NLO_Res_Pol(1) = (BosonicPartAmp(1,-1)); print *, 'eps cancellation check'
!        NLO_Res_Pol(2) = (BosonicPartAmp(2,-1)); print *, 'eps cancellation check'

        ! Z' virt, glue tree
        LO_res_pol_glue = BornAmps(1)%Result
        NLO_Res_Unpol_Left = NLO_Res_UnPol_Left + (nc**2-1d0) * dreal(LO_res_pol_glue*dconjg(Virt_Res_Pol_Left)) ! factor two taken care of by alpha_sOver2Pi
        NLO_Res_Unpol_Right = NLO_Res_UnPol_Right + (nc**2-1d0) * dreal(LO_res_pol_glue*dconjg(Virt_Res_Pol_Right))

        ! Z' tree, glue loop
        NLO_Res_UnPol_Left = NLO_Res_UnPol_Left + dreal( (nc * NLO_Res_Pol(1) + nc**2 * NLO_Res_Pol(2)) * dconjg(LO_Res_Pol_Left))
        NLO_Res_UnPol_Right = NLO_Res_UnPol_Right + dreal( (nc * NLO_Res_Pol(1) + nc**2 * NLO_Res_Pol(2)) * dconjg(LO_Res_Pol_Right))

     enddo!helicity loop

     EvalCS_Virt_Zprime_Interf = ISFac*PreFac* (PDFFac_L * NLO_Res_Unpol_Left + PDFFac_R * NLO_Res_Unpol_Right)
     EvalCS_Virt_Zprime_Interf = EvalCS_Virt_Zprime_Interf * (alpha_sOver2Pi*RunFactor) * (4d0*Pi*alpha_s*RunFactor)


  ELSEIF ( CORRECTION.EQ.3 ) THEN 

     !!! Integrated Dipoles !!!

!Ga_Zpr=0d0; print *, "MARKUS: ga_zpr to zero for checks"
!gR_Zpr(top_) = 1d0; gL_Zpr(top_) = 1d0! vector contribution
!gR_Zpr(top_) = -1d0; gL_Zpr(top_) = 1d0! axial contribution
!m_zpr=91d0*GeV

     do iHel=nHel(1),nHel(2)
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
        call SetPolarizations()
    
        ! Compute first the Zprime part
        call Tree_Zprime_tbtqbq(LO_Res_Pol_Left, LO_Res_Pol_Right)
        
        ! Now compute the gluon part
        call Tree_Zprime_tbtqbq_gluefordip(LO_Res_Pol_glue)

        LO_for_dip_Left  = LO_for_dip_Left  + dreal(LO_res_pol_Left*dconjg(LO_res_pol_glue))
        LO_for_dip_Right = LO_for_dip_Right + dreal(LO_res_pol_Right*dconjg(LO_res_pol_glue))


     enddo!helicity loop

!beta = dsqrt(1d0-4d0*m_top**2/Ehat**2)
!z = Get_CosTheta(Momext(1:4,3))
!print *, beta,z

! vector contribution
! AndreasResult = 2d0*ehat**2/(ehat**2-m_Zpr**2) *( 2d0-beta**2*(1-z**2) )! vector contribution
! print *, "TOPAZ: ", (LO_for_dip_Left+LO_for_dip_Right)
! print *, "Andreas: ",AndreasResult
! print *, "ratio vector", AndreasResult / (LO_for_dip_Left  + LO_for_dip_Right )! vector contribution

! axial contribution
!AndreasResult = 2d0*ehat**2/(ehat**2-m_Zpr**2) *( 2d0*beta*z )! axial  contribution
!print *, "TOPAZ: ", (-LO_for_dip_Left+LO_for_dip_Right)
!print *, "Andreas: ",AndreasResult
!print *, "ratio axial", AndreasResult / (-LO_for_dip_Left  + LO_for_dip_Right )! axial  contribution
!pause

     IF ( TOPDECAYS.EQ.0 ) THEN 
        z = yRnd(5)
     ELSE
        z = yRnd(13)
     END IF

     call setPDFs(eta1/z,eta2/z,MuFac,pdf_z)

     IF ( PROCESS.EQ.67 ) THEN

        ! qg initial state !
        
        PDFFac_dip_L(1) = 0d0
        PDFFac_dip_L(2) = 0d0

        if( xpdf.le.0.5d0 ) then
            PDFFac_dip_L(3) = gL_Zpr(up_) * (pdf_z(0,1)*pdf(up_,2))
            PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(dn_) * (pdf_z(0,1)*pdf(dn_,2))
            PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(chm_) * (pdf_z(0,1)*pdf(chm_,2))
            PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(str_) * (pdf_z(0,1)*pdf(str_,2))
            PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(bot_) * (pdf_z(0,1)*pdf(bot_,2))
            PDFFac_dip_R(3) = gR_Zpr(up_) * (pdf_z(0,1)*pdf(up_,2))
            PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(dn_) * (pdf_z(0,1)*pdf(dn_,2))
            PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(chm_) * (pdf_z(0,1)*pdf(chm_,2))
            PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(str_) * (pdf_z(0,1)*pdf(str_,2))
            PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(bot_) * (pdf_z(0,1)*pdf(bot_,2))
        else
            PDFFac_dip_L(3) = gL_Zpr(up_) * (pdf_z(0,2)*pdf(up_,1))
            PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(dn_) * (pdf_z(0,2)*pdf(dn_,1))
            PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(chm_) * (pdf_z(0,2)*pdf(chm_,1))
            PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(str_) * (pdf_z(0,2)*pdf(str_,1))
            PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(bot_) * (pdf_z(0,2)*pdf(bot_,1))
            PDFFac_dip_R(3) = gR_Zpr(up_) * (pdf_z(0,2)*pdf(up_,1))
            PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(dn_) * (pdf_z(0,2)*pdf(dn_,1))
            PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(chm_) * (pdf_z(0,2)*pdf(chm_,1))
            PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(str_) * (pdf_z(0,2)*pdf(str_,1))
            PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(bot_) * (pdf_z(0,2)*pdf(bot_,1))
        endif

        PDFFac_dip_R(3) = PDFFac_dip_R(3) * 2d0!  factor two for sampling over pdfs
        PDFFac_dip_L(3) = PDFFac_dip_L(3) * 2d0

        call IntDip_qg_ZprimeInt_ttb(MomExt(1:4,1:4), z, IDip) ! Normalization: as/(2 Pi)


     ELSEIF (PROCESS.EQ.68) THEN
        
        ! qbg initial state !

        PDFFac_dip_L(1) = 0d0
        PDFFac_dip_L(3) = 0d0

        if( xpdf.le.0.5d0 ) then
           PDFFac_dip_L(2) = gL_Zpr(up_) * (pdf_z(0,1)*pdf(aup_,2))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(dn_) * (pdf_z(0,1)*pdf(adn_,2))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(chm_) * (pdf_z(0,1)*pdf(achm_,2))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(str_) * (pdf_z(0,1)*pdf(astr_,2))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(bot_) * (pdf_z(0,1)*pdf(abot_,2))
           PDFFac_dip_R(2) = gR_Zpr(up_) * (pdf_z(0,1)*pdf(aup_,2))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(dn_) * (pdf_z(0,1)*pdf(adn_,2))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(chm_) * (pdf_z(0,1)*pdf(achm_,2))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(str_) * (pdf_z(0,1)*pdf(astr_,2))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(bot_) * (pdf_z(0,1)*pdf(abot_,2))
        else
           PDFFac_dip_L(2) = gL_Zpr(up_) * (pdf_z(0,2)*pdf(aup_,1))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(dn_) * (pdf_z(0,2)*pdf(adn_,1))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(chm_) * (pdf_z(0,2)*pdf(achm_,1))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(str_) * (pdf_z(0,2)*pdf(astr_,1))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(bot_) * (pdf_z(0,2)*pdf(abot_,1))
           PDFFac_dip_R(2) = gR_Zpr(up_) * (pdf_z(0,2)*pdf(aup_,1))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(dn_) * (pdf_z(0,2)*pdf(adn_,1))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(chm_) * (pdf_z(0,2)*pdf(achm_,1))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(str_) * (pdf_z(0,2)*pdf(astr_,1))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(bot_) * (pdf_z(0,2)*pdf(abot_,1))
        endif

        PDFFac_dip_R(2) = PDFFac_dip_R(2) * 2d0!  factor two for sampling over pdfs
        PDFFac_dip_L(2) = PDFFac_dip_L(2) * 2d0

        call IntDip_qbg_ZprimeInt_ttb(MomExt(1:4,1:4), z, IDip) ! Normalization: as/(2 Pi)

     ELSEIF ( PROCESS.EQ.69 ) THEN 
     
        ! qqb initial state !

        PDFFac_dip_L(1) = PDFFac_L
        PDFFac_dip_R(1) = PDFFac_R

        if( xpdf.le.0.5d0 ) then
           PDFFac_dip_L(2) = gL_Zpr(up_) * (pdf_z(up_,1)*pdf(aup_,2))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(dn_) * (pdf_z(dn_,1)*pdf(adn_,2))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(chm_) * (pdf_z(chm_,1)*pdf(achm_,2))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(str_) * (pdf_z(str_,1)*pdf(astr_,2))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(bot_) * (pdf_z(bot_,1)*pdf(abot_,2))
           PDFFac_dip_L(3) = gL_Zpr(up_) * (pdf(up_,1)*pdf_z(aup_,2))
           PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(dn_) * (pdf(dn_,1)*pdf_z(adn_,2))
           PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(chm_) * (pdf(chm_,1)*pdf_z(achm_,2))
           PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(str_) * (pdf(str_,1)*pdf_z(astr_,2))
           PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(bot_) * (pdf(bot_,1)*pdf_z(abot_,2))
           PDFFac_dip_R(2) = gR_Zpr(up_) * (pdf_z(up_,1)*pdf(aup_,2))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(dn_) * (pdf_z(dn_,1)*pdf(adn_,2))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(chm_) * (pdf_z(chm_,1)*pdf(achm_,2))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(str_) * (pdf_z(str_,1)*pdf(astr_,2))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(bot_) * (pdf_z(bot_,1)*pdf(abot_,2))
           PDFFac_dip_R(3) = gR_Zpr(up_) * (pdf(up_,1)*pdf_z(aup_,2))
           PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(dn_) * (pdf(dn_,1)*pdf_z(adn_,2))
           PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(chm_) * (pdf(chm_,1)*pdf_z(achm_,2))
           PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(str_) * (pdf(str_,1)*pdf_z(astr_,2))
           PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(bot_) * (pdf(bot_,1)*pdf_z(abot_,2))
        else
           PDFFac_dip_L(2) = gL_Zpr(up_) * (pdf_z(up_,2)*pdf(aup_,1))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(dn_) * (pdf_z(dn_,2)*pdf(adn_,1))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(chm_) * (pdf_z(chm_,2)*pdf(achm_,1))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(str_) * (pdf_z(str_,2)*pdf(astr_,1))
           PDFFac_dip_L(2) = PDFFac_dip_L(2) + gL_Zpr(bot_) * (pdf_z(bot_,2)*pdf(abot_,1))
           PDFFac_dip_L(3) = gL_Zpr(up_) * (pdf(up_,2)*pdf_z(aup_,1))
           PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(dn_) * (pdf(dn_,2)*pdf_z(adn_,1))
           PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(chm_) * (pdf(chm_,2)*pdf_z(achm_,1))
           PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(str_) * (pdf(str_,2)*pdf_z(astr_,1))
           PDFFac_dip_L(3) = PDFFac_dip_L(3) + gL_Zpr(bot_) * (pdf(bot_,2)*pdf_z(abot_,1))
           PDFFac_dip_R(2) = gR_Zpr(up_) * (pdf_z(up_,2)*pdf(aup_,1))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(dn_) * (pdf_z(dn_,2)*pdf(adn_,1))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(chm_) * (pdf_z(chm_,2)*pdf(achm_,1))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(str_) * (pdf_z(str_,2)*pdf(astr_,1))
           PDFFac_dip_R(2) = PDFFac_dip_R(2) + gR_Zpr(bot_) * (pdf_z(bot_,2)*pdf(abot_,1))
           PDFFac_dip_R(3) = gR_Zpr(up_) * (pdf(up_,2)*pdf_z(aup_,1))
           PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(dn_) * (pdf(dn_,2)*pdf_z(adn_,1))
           PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(chm_) * (pdf(chm_,2)*pdf_z(achm_,1))
           PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(str_) * (pdf(str_,2)*pdf_z(astr_,1))
           PDFFac_dip_R(3) = PDFFac_dip_R(3) + gR_Zpr(bot_) * (pdf(bot_,2)*pdf_z(abot_,1))
        endif

        PDFFac_dip_R(2) = PDFFac_dip_R(2) * 2d0!  factor two for sampling over pdfs
        PDFFac_dip_L(2) = PDFFac_dip_L(2) * 2d0

        PDFFac_dip_R(3) = PDFFac_dip_R(3) * 2d0
        PDFFac_dip_L(3) = PDFFac_dip_L(3) * 2d0
        
        call IntDip_qqb_ZprimeInt_ttb(MomExt(1:4,1:4), z, IDip) ! Normalization: as/(2 Pi)

        IDip = IDip * (nc**2-1d0)

     ENDIF ! Integrated Dipoles process
     
     IDip = IDip * (alpha_s*RunFactor)/(2d0*Pi) * PreFac * (4d0*Pi*alpha_s*RunFactor)
        
     IDipAmp = IDip(1) * (PDFFac_dip_L(1) * LO_for_dip_Left + PDFFac_dip_R(1) * LO_for_dip_Right)
     IDipAmp = IDipAmp + IDip(2)/z * (PDFFac_dip_L(2) * LO_for_dip_Left + PDFFac_dip_R(2) * LO_for_dip_Right)
     IDipAmp = IDipAmp + IDip(3)/z * (PDFFac_dip_L(3) * LO_for_dip_Left + PDFFac_dip_R(3) * LO_for_dip_Right)
     
     IDipAmp = IDipAmp * ISFac 

     EvalCS_Virt_Zprime_Interf = IDipAmp
     
  ENDIF ! Correction

!--F Oct 15: 1/ep poles cancellation checked
!     print *, 'virt', EvalCS_Virt_Zprime_Interf
!     print *, 'idip', IDipAmp
!     pause


  do NHisto=1,NumHistograms
     call intoHisto(NHisto,NBin(NHisto),EvalCS_Virt_Zprime_Interf)
  enddo


  EvalCS_Virt_Zprime_Interf = EvalCS_Virt_Zprime_Interf/VgsWgt

  EvalCounter = EvalCounter + 1


return

END FUNCTION EvalCS_Virt_Zprime_Interf



FUNCTION EvalCS_Real_Zprime_Interf(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes
use ModAmplitudes_Zprime
use ModParameters
use ModDipoles_Zprime
implicit none
real(8) ::  EvalCS_Real_Zprime_Interf,yRnd(1:VegasMxDim),VgsWgt,DipoleResult,EvalCS_Dips_Zprime_ttbqqb
complex(8) :: Res_Pol_L, Res_Pol_R, prim1_L, prim1_R, prim2_L, prim2_R
real(8) :: Res_UnPol_L, Res_UnPol_R
integer :: iHel
real(8) :: EHat,RunFactor,PSWgt,PS,PSWgt2,PSWgt3,xFrag, sHatJacobi
real(8) :: FluxFac, PreFac, ISFac
real(8) :: pdf(-6:6,1:2),eta1,eta2, PDFFac_L, PDFFac_R
real(8) :: MomExt(1:4,1:5),MomDK(1:4,1:6), MomExtTd(1:4,1:4), MomDKTd(1:4,1:6)
integer :: NBin(1:NumMaxHisto),NHisto
logical :: applyPSCut,applySingCut
real(8) :: colf
real(8),parameter :: ThresholdCutOff = 1.00d0
real(8),parameter :: Cf = 4d0/3d0
integer :: nDip
real(8) :: resdip, dipoles
real(8) :: resLO_L, resLO_R
complex(8) :: ampLO_L, ampLO_R, ampLO_QCD
real(8) :: xpdf
include "vegas_common.f"
integer :: i1,i2,i3,i4,i5, iPrimAmp

  EvalCS_Real_Zprime_interf = 0d0

!   yrnd(1:15) = 0.4d0; yrnd(8)=0.6d0; yrnd(16)=0.6d0; print *, "fixed y,y16"
  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi) ! no mapping 

  if( EHat.le.2d0*m_Top * ThresholdCutOff ) then
      EvalCS_Real_Zprime_interf = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

  call EvalPhaseSpace_2to3(EHat,yRnd(3:7),MomExt(1:4,1:5),PSWgt) ! no mapping
  
  call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))

  call CheckSing(MomExt,applySingCut)
  if( applySingCut ) then
     EvalCS_Real_Zprime_interf = 0d0
     SkipCounter = SkipCounter + 1
     return
  endif
  
  IF( TOPDECAYS.NE.0 ) THEN
     call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomDK(1:4,1:3),PSWgt2) ! anti-top
     call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:15),.false.,MomDK(1:4,4:6),PSWgt3)
     PSWgt = PSWgt * PSWgt2*PSWgt3
  ENDIF

  call Kinematics_TTBARZprime(.true.,MomExt,MomDK,applyPSCut,NBin)

  call setPDFs(eta1,eta2,MuFac,pdf)


  IF (PROCESS.EQ.67 ) THEN !-- qg channel

     call random_number(xpdf)

     if( xpdf.le.0.5d0 ) then 
        PDFFac_L =            gL_Zpr(up_) * (pdf(0,2)*pdf(up_,1))
        PDFFac_L = PDFFac_L + gL_Zpr(dn_) * (pdf(0,2)*pdf(dn_,1))
        PDFFac_L = PDFFac_L + gL_Zpr(chm_) * (pdf(0,2)*pdf(chm_,1))
        PDFFac_L = PDFFac_L + gL_Zpr(str_) * (pdf(0,2)*pdf(str_,1))
        PDFFac_L = PDFFac_L + gL_Zpr(bot_) * (pdf(0,2)*pdf(bot_,1))
        PDFFac_R =            gR_Zpr(up_) * (pdf(0,2)*pdf(up_,1))
        PDFFac_R = PDFFac_R + gR_Zpr(dn_) * (pdf(0,2)*pdf(dn_,1))
        PDFFac_R = PDFFac_R + gR_Zpr(chm_) * (pdf(0,2)*pdf(chm_,1))
        PDFFac_R = PDFFac_R + gR_Zpr(str_) * (pdf(0,2)*pdf(str_,1))
        PDFFac_R = PDFFac_R + gR_Zpr(bot_) * (pdf(0,2)*pdf(bot_,1))
     else
        PDFFac_L =            gL_Zpr(up_) * (pdf(0,1)*pdf(up_,2))
        PDFFac_L = PDFFac_L + gL_Zpr(dn_) * (pdf(0,1)*pdf(dn_,2))
        PDFFac_L = PDFFac_L + gL_Zpr(chm_) * (pdf(0,1)*pdf(chm_,2))
        PDFFac_L = PDFFac_L + gL_Zpr(str_) * (pdf(0,1)*pdf(str_,2))
        PDFFac_L = PDFFac_L + gL_Zpr(bot_) * (pdf(0,1)*pdf(bot_,2))
        PDFFac_R =            gR_Zpr(up_) * (pdf(0,1)*pdf(up_,2))
        PDFFac_R = PDFFac_R + gR_Zpr(dn_) * (pdf(0,1)*pdf(dn_,2))
        PDFFac_R = PDFFac_R + gR_Zpr(chm_) * (pdf(0,1)*pdf(chm_,2))
        PDFFac_R = PDFFac_R + gR_Zpr(str_) * (pdf(0,1)*pdf(str_,2))
        PDFFac_R = PDFFac_R + gR_Zpr(bot_) * (pdf(0,1)*pdf(bot_,2))
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
     endif

     PDFFac_R = PDFFac_R * 2d0!  factor two for sampling over pdfs
     PDFFac_L = PDFFac_L * 2d0

  ELSEIF (PROCESS.EQ.68 ) THEN !-- gqb channel

     call random_number(xpdf)

    if( xpdf.le.0.5d0 ) then
        PDFFac_L = gL_Zpr(up_) * (pdf(0,1)*pdf(aup_,2))
        PDFFac_L = PDFFac_L + gL_Zpr(dn_) * (pdf(0,1)*pdf(adn_,2))
        PDFFac_L = PDFFac_L + gL_Zpr(chm_) * (pdf(0,1)*pdf(achm_,2))
        PDFFac_L = PDFFac_L + gL_Zpr(str_) * (pdf(0,1)*pdf(astr_,2))
        PDFFac_L = PDFFac_L + gL_Zpr(bot_) * (pdf(0,1)*pdf(abot_,2))
        PDFFac_R = gR_Zpr(up_) * (pdf(0,1)*pdf(aup_,2))
        PDFFac_R = PDFFac_R + gR_Zpr(dn_) * (pdf(0,1)*pdf(adn_,2))
        PDFFac_R = PDFFac_R + gR_Zpr(chm_) * (pdf(0,1)*pdf(achm_,2))
        PDFFac_R = PDFFac_R + gR_Zpr(str_) * (pdf(0,1)*pdf(astr_,2))
        PDFFac_R = PDFFac_R + gR_Zpr(bot_) * (pdf(0,1)*pdf(abot_,2))
     else
        PDFFac_L = gL_Zpr(up_) * (pdf(0,2)*pdf(aup_,1))
        PDFFac_L = PDFFac_L + gL_Zpr(dn_) * (pdf(0,2)*pdf(adn_,1))
        PDFFac_L = PDFFac_L + gL_Zpr(chm_) * (pdf(0,2)*pdf(achm_,1))
        PDFFac_L = PDFFac_L + gL_Zpr(str_) * (pdf(0,2)*pdf(astr_,1))
        PDFFac_L = PDFFac_L + gL_Zpr(bot_) * (pdf(0,2)*pdf(abot_,1))
        PDFFac_R = gR_Zpr(up_) * (pdf(0,2)*pdf(aup_,1))
        PDFFac_R = PDFFac_R + gR_Zpr(dn_) * (pdf(0,2)*pdf(adn_,1))
        PDFFac_R = PDFFac_R + gR_Zpr(chm_) * (pdf(0,2)*pdf(achm_,1))
        PDFFac_R = PDFFac_R + gR_Zpr(str_) * (pdf(0,2)*pdf(astr_,1))
        PDFFac_R = PDFFac_R + gR_Zpr(bot_) * (pdf(0,2)*pdf(abot_,1))
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
     endif
     
     PDFFac_R = PDFFac_R * 2d0!  factor two for sampling over pdfs
     PDFFac_L = PDFFac_L * 2d0

  ELSEIF( PROCESS.EQ.69 ) THEN !-- qqb channel

    call random_number(xpdf)

!    xpdf = 0.4d0; print *, 'fixed xpdf'

    if( xpdf.le.0.5d0 ) then
       PDFFac_L =            gL_Zpr(up_) * (pdf(up_,1)*pdf(aup_,2))
       PDFFac_L = PDFFac_L + gL_Zpr(dn_) * (pdf(dn_,1)*pdf(adn_,2))
       PDFFac_L = PDFFac_L + gL_Zpr(chm_) * (pdf(chm_,1)*pdf(achm_,2))
       PDFFac_L = PDFFac_L + gL_Zpr(str_) * (pdf(str_,1)*pdf(astr_,2))
       PDFFac_L = PDFFac_L + gL_Zpr(bot_) * (pdf(bot_,1)*pdf(abot_,2))
       
       PDFFac_R =            gR_Zpr(up_) * (pdf(up_,1)*pdf(aup_,2))
       PDFFac_R = PDFFac_R + gR_Zpr(dn_) * (pdf(dn_,1)*pdf(adn_,2))
       PDFFac_R = PDFFac_R + gR_Zpr(chm_) * (pdf(chm_,1)*pdf(achm_,2))
       PDFFac_R = PDFFac_R + gR_Zpr(str_) * (pdf(str_,1)*pdf(astr_,2))
       PDFFac_R = PDFFac_R + gR_Zpr(bot_) * (pdf(bot_,1)*pdf(abot_,2))
    else
       PDFFac_L =            gL_Zpr(up_) * (pdf(up_,2)*pdf(aup_,1))
       PDFFac_L = PDFFac_L + gL_Zpr(dn_) * (pdf(dn_,2)*pdf(adn_,1))
       PDFFac_L = PDFFac_L + gL_Zpr(chm_) * (pdf(chm_,2)*pdf(achm_,1))
       PDFFac_L = PDFFac_L + gL_Zpr(str_) * (pdf(str_,2)*pdf(astr_,1))
       PDFFac_L = PDFFac_L + gL_Zpr(bot_) * (pdf(bot_,2)*pdf(abot_,1))
       
       PDFFac_R =            gR_Zpr(up_) * (pdf(up_,2)*pdf(aup_,1))
       PDFFac_R = PDFFac_R + gR_Zpr(dn_) * (pdf(dn_,2)*pdf(adn_,1))
       PDFFac_R = PDFFac_R + gR_Zpr(chm_) * (pdf(chm_,2)*pdf(achm_,1))
       PDFFac_R = PDFFac_R + gR_Zpr(str_) * (pdf(str_,2)*pdf(astr_,1))
       PDFFac_R = PDFFac_R + gR_Zpr(bot_) * (pdf(bot_,2)*pdf(abot_,1))
       call swapMom(MomExt(1:4,1),MomExt(1:4,2))
    endif
  
  PDFFac_R = PDFFac_R * 2d0!  factor two for sampling over pdfs
  PDFFac_L = PDFFac_L * 2d0

  ENDIF

  PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt

  RunFactor = RunAlphaS(NLOParam,MuRen)
  ISFac = MomCrossing(MomExt)

  IF( TOPDECAYS.GE.1 ) THEN
     call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
     call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
  ENDIF

!if (1.eq.1) then !!! INVERTED DIPOLES CHECK 
  if( applyPSCut ) then
     EvalCS_Real_Zprime_interf = 0d0
  else

     PrimAmp1_12345 = 1
     PrimAmp1_15234 = 2
     PrimAmp1_12534 = 3
     PrimAmp1_12354 = 4

     colf = 8d0 ! NC^2 - 1
     Res_UnPol_L = 0d0
     Res_UnPol_R = 0d0
     do iHel=1,NumHelicities
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
        call SetPolarizations()
        
        ! Z' amplitudes
        call Tree_Zprime_tbtqbqg_i(prim1_L, prim1_R)
        call Tree_Zprime_tbtqbqg_f(prim2_L, prim2_R)

        ! QCD amplitudes
        do iPrimAmp=1,NumPrimAmps
           call EvalTree(BornAmps(iPrimAmp))
        enddo

        !--F Oct 14: sign fixed by comparing against MadGraph
        Res_UnPol_L = Res_UnPol_L - ( + 2d0 * dreal(BornAmps(PrimAmp1_12345)%Result * dconjg(prim1_L)) + &
             2d0 * dreal(BornAmps(PrimAmp1_12345)%Result * dconjg(prim2_L)) + &
             2d0 * dreal(BornAmps(PrimAmp1_15234)%Result * dconjg(prim2_L)) + &
             2d0 * dreal(BornAmps(PrimAmp1_12534)%Result * dconjg(prim1_L)) + &
             2d0 * dreal(BornAmps(PrimAmp1_12534)%Result * dconjg(prim2_L)) + &
             2d0 * dreal(BornAmps(PrimAmp1_12354)%Result * dconjg(prim1_L)))

        !--F Oct 14: sign fixed by comparing against MadGraph
        Res_UnPol_R = Res_UnPol_R - ( 2d0 * dreal(BornAmps(PrimAmp1_12345)%Result * dconjg(prim1_R)) + &
             2d0 * dreal(BornAmps(PrimAmp1_12345)%Result * dconjg(prim2_R)) + &
             2d0 * dreal(BornAmps(PrimAmp1_15234)%Result * dconjg(prim2_R)) + &
             2d0 * dreal(BornAmps(PrimAmp1_12534)%Result * dconjg(prim1_R)) + &
             2d0 * dreal(BornAmps(PrimAmp1_12534)%Result * dconjg(prim2_R)) + &
             2d0 * dreal(BornAmps(PrimAmp1_12354)%Result * dconjg(prim1_R)) )

     enddo!helicity loop

     EvalCS_Real_Zprime_interf = (PDFFac_R * Res_UnPol_R + PDFFac_L * Res_UnPol_L)
     EvalCS_Real_Zprime_interf = EvalCS_Real_Zprime_interf * (4d0*Pi*alpha_s*RunFactor)**2 * PreFac * ISFac * colf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--F Oct 14 MadGraph Check 
!     open(unit=66,file='/home/caola/tools/MG_ME_V4.5.2/zprime/SubProcesses/kin.dat', status='unknown')
!     do iHel = 1,5
!        write(66,*) MomExt(1,iHel) * 100d0
!        write(66,*) MomExt(2,iHel) * 100d0
!        write(66,*) MomExt(3,iHel) * 100d0
!        write(66,*) MomExt(4,iHel) * 100d0
!     enddo
!     close(66)
!     print *, 'MZpr', m_zpr * 100d0
!     print *, 'Ga_zpr', Ga_Zpr * 100d0
!     print *, 'gL_Zpr(up_), gL_Zpr(dn_)', gL_Zpr(up_), gL_Zpr(dn_)
!     print *, 'gR_Zpr(up_), gR_Zpr(dn_)', gR_Zpr(up_), gR_Zpr(dn_)
!     print *, 'alpha_s', alpha_s*RunFactor
!     print *, 'gs', dsqrt(alpha_s*RunFactor*4d0*Pi)
!     print *, ''
!     print *, 'up',  ISFac * colf * (4d0*Pi*alpha_s*RunFactor)**2 * (gL_Zpr(up_) * Res_UnPol_L + gR_Zpr(up_) * Res_UnPol_R)*GeV**2
!     print *, 'dn',  ISFac * colf * (4d0*Pi*alpha_s*RunFactor)**2 * (gL_Zpr(dn_) * Res_UnPol_L + gR_Zpr(dn_) * Res_UnPol_R)*GeV**2
!     pause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     do NHisto=1,NumHistograms
        call intoHisto(NHisto,NBin(NHisto),EvalCS_Real_Zprime_interf)
     enddo

  endif!applyPSCut


  !!! Dipoles section
  colf = 8d0 ! NC^2 - 1
  
  dipoles = 0d0

  IF ( PROCESS .EQ. 69 ) THEN ! qqb

     do nDip = 1,8

        call Dipoles_qqb_Zprime_interf(nDip, MomExt(1:4,1:5), MomExtTd(1:4,1:5), resdip)

        if (resdip .EQ. 0d0) cycle

        Crossing(:) = (/3,4,-1,-2,0/) ! For dipoles
        NumExtParticles = 4
        ISFac = MomCrossing(MomExtTd)
        Crossing(:) = (/4,5,-1,-2,3/) ! Original: this was wrong, corrected
        NumExtParticles = 5
        
        IF( TOPDECAYS.NE.0 ) THEN
           call EvalPhasespace_TopDecay(MomExtTd(1:4,3),yRnd(8:11),.false.,MomDKTd(1:4,1:3),PSWgt2)
           call EvalPhasespace_TopDecay(MomExtTd(1:4,4),yRnd(12:15),.false.,MomDKTd(1:4,4:6),PSWgt3)
!           PSWgt = PSWgt * PSWgt2*PSWgt3! this is better never be used anywhere below
           
           call Kinematics_TTBARZprime(.false.,MomExtTd,MomDKTd,applyPSCut,NBin)

           if ( applyPSCut ) then
              resdip = 0d0
           else
              
              call TopDecay(ExtParticle(1),DK_LO,MomDKTd(1:4,1:3))
              call TopDecay(ExtParticle(2),DK_LO,MomDKTd(1:4,4:6))
              
              resLO_L = 0d0
              resLO_R = 0d0
              
              i1 = 0
              i2 = 0
              do i3 = -1,1,2
                 ExtParticle(1)%Helicity = i1
                 ExtParticle(2)%Helicity = i2
                 ExtParticle(3)%Helicity = i3
                 ExtParticle(4)%Helicity = -i3
                 call SetPolarizations()
                 call Tree_Zprime_tbtqbq_gluefordip(ampLO_QCD)
                 call Tree_Zprime_tbtqbq(ampLO_L, ampLO_R)
                 resLO_L = resLO_L +  2d0 * dreal(ampLO_L * dconjg(ampLO_QCD))
                 resLO_R = resLO_R +  2d0 * dreal(ampLO_R * dconjg(ampLO_QCD))
              enddo
              
           endif
           
        ELSE
           
           call Kinematics_TTBARZprime(.false.,MomExtTd,MomDK,applyPSCut,NBin)
           
           if ( applyPSCut ) then
              resdip = 0d0
           else
              
              resLO_L = 0d0
              resLO_R = 0d0
              
              do i1 = -1,1,2
                 do i2 = -1,1,2
                    do i3 = -1,1,2
                       ExtParticle(1)%Helicity = i1
                       ExtParticle(2)%Helicity = i2
                       ExtParticle(3)%Helicity = i3
                       ExtParticle(4)%Helicity = -i3
                       call SetPolarizations()
                       call Tree_Zprime_tbtqbq_gluefordip(ampLO_QCD)
                       call Tree_Zprime_tbtqbq(ampLO_L, ampLO_R)
                       resLO_L = resLO_L +  2d0 * dreal(ampLO_L * dconjg(ampLO_QCD))
                       resLO_R = resLO_R +  2d0 * dreal(ampLO_R * dconjg(ampLO_QCD))
                    enddo
                 enddo
              enddo
              
           endif
           
        ENDIF !TopDecays
        
        ! Check against Markus

        resdip = resdip * (PDFFac_R * resLO_R + PDFFac_L * resLO_L)
        resdip = resdip * (4d0*Pi*alpha_s*RunFactor)**2 * PreFac * ISFac * colf
        !-- FM Nov 7, sign changed to agree with MadGraph BEFORE filling histogram!
        resdip = -resdip

        do NHisto=1,NumHistograms
           call intoHisto(NHisto,NBin(NHisto),resdip)
        enddo

        dipoles = dipoles + resdip

     enddo

  ENDIF

!  print *, 'inv    ', MomExt(1,3)/Ehat !(MomExt(1:4,1).dot.MomExt(1:4,3))/Ehat**2,(MomExt(1:4,2).dot.MomExt(1:4,3))/Ehat**2
!  print *, 'real   ', EvalCS_Real_Zprime_interf
!  print *, 'dipoles', dipoles
!  print *, 'dipoles/real+1', dipoles/EvalCS_Real_Zprime_interf + 1d0
!  pause 

  EvalCS_Real_Zprime_interf = (EvalCS_Real_Zprime_interf+dipoles)/VgsWgt

  RETURN

END FUNCTION EvalCS_Real_Zprime_Interf



FUNCTION EvalCS_NLODK_Zprime_ttb(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes_Zprime
use ModParameters
use ModHadrWDecay
implicit none
#DEFINE TheReal_Top 2233
#DEFINE TheDipole_Top 2244
#DEFINE TheDipole_ATop 2255
#DEFINE TheWm_decay 2260
#DEFINE TheReal_Top_Wdecay 2266
#DEFINE TheDipole_Top_Wdecay 2277
#DEFINE TheWp_decay 2280
#DEFINE TheDipole_ATop_Wdecay 2288
#DEFINE TheEnd 2299
real(8) ::  EvalCS_NLODK_Zprime_ttb, yRnd(1:VegasMxDim),VgsWgt
complex(8) :: LO_Res_Pol_Left, LO_Res_Pol_Right
complex(8) :: NLO_Res_Pol_Left, NLO_Res_Pol_Right, NLO_Res_Unpol_Left, NLO_Res_Unpol_Right, NLO_Res_Unpol
complex(8) :: Dip_Res_Unpol
integer :: iHel,jHel,kHel,GluHel,iPrimAmp,jPrimAmp,ndip
real(8) :: EHat,PSWgt1,PSWgt2,PSWgt3,ISFac,dip_res_w
real(8) :: MomExt(1:4,1:NumExtParticles),MomDK(1:4,1:7),MomDKTd(1:4,1:6),MomDKx(1:4,1:7)
logical :: applyPSCut,applySingCut
real(8) :: tau,eta1,eta2,sHatJacobi,PreFac,FluxFac,RunFactor
real(8) :: pdf(-6:6,1:2), PDFFac_L, PDFFac_R,xpdf
integer :: NBin(1:NumMaxHisto),NHisto,npdf
real(8) :: pbDpg,ptDpg,ptDpb,z,omz,Dipole,rsq,y
real(8), parameter :: CF=4d0/3d0,PhotonCouplCorr=2d0
real(8) :: MomBoost(1:4),MomLep1(1:4),MomLep2(1:4)
integer,parameter :: up=1,dn=2,glu=1
include "vegas_common.f"




  EvalCS_NLODK_Zprime_ttb = 0d0

  call PDFMapping(62,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top) then
      EvalCS_NLODK_Zprime_ttb = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

  call EvalPhaseSpace_2to2(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt1)
  call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))

  call setPDFs(eta1,eta2,MuFac,pdf)
  call random_number(xpdf)
  if( xpdf.le.0.5d0 ) then
     PDFFac_L = gL_Zpr(up_)**2 * (pdf(up_,1)*pdf(aup_,2))
     PDFFac_L = PDFFac_L + gL_Zpr(dn_)**2 * (pdf(dn_,1)*pdf(adn_,2))
     PDFFac_L = PDFFac_L + gL_Zpr(chm_)**2 * (pdf(chm_,1)*pdf(achm_,2))
     PDFFac_L = PDFFac_L + gL_Zpr(str_)**2 * (pdf(str_,1)*pdf(astr_,2))
     PDFFac_L = PDFFac_L + gL_Zpr(bot_)**2 * (pdf(bot_,1)*pdf(abot_,2))
     PDFFac_R = gR_Zpr(up_)**2 * (pdf(up_,1)*pdf(aup_,2))
     PDFFac_R = PDFFac_R + gR_Zpr(dn_)**2 * (pdf(dn_,1)*pdf(adn_,2))
     PDFFac_R = PDFFac_R + gR_Zpr(chm_)**2 * (pdf(chm_,1)*pdf(achm_,2))
     PDFFac_R = PDFFac_R + gR_Zpr(str_)**2 * (pdf(str_,1)*pdf(astr_,2))
     PDFFac_R = PDFFac_R + gR_Zpr(bot_)**2 * (pdf(bot_,1)*pdf(abot_,2))
  else
     PDFFac_L = gL_Zpr(up_)**2 * (pdf(Up_,2)*pdf(aup_,1))
     PDFFac_L = PDFFac_L + gL_Zpr(dn_)**2 * (pdf(dn_,2)*pdf(adn_,1))
     PDFFac_L = PDFFac_L + gL_Zpr(chm_)**2 * (pdf(chm_,2)*pdf(achm_,1))
     PDFFac_L = PDFFac_L + gL_Zpr(str_)**2 * (pdf(str_,2)*pdf(astr_,1))
     PDFFac_L = PDFFac_L + gL_Zpr(bot_)**2 * (pdf(bot_,2)*pdf(abot_,1))
     PDFFac_R = gR_Zpr(up_)**2 * (pdf(Up_,2)*pdf(aup_,1))
     PDFFac_R = PDFFac_R + gR_Zpr(dn_)**2 * (pdf(dn_,2)*pdf(adn_,1))
     PDFFac_R = PDFFac_R + gR_Zpr(chm_)**2 * (pdf(chm_,2)*pdf(achm_,1))
     PDFFac_R = PDFFac_R + gR_Zpr(str_)**2 * (pdf(str_,2)*pdf(astr_,1))
     PDFFac_R = PDFFac_R + gR_Zpr(bot_)**2 * (pdf(bot_,2)*pdf(abot_,1))
     call swapMom(MomExt(1:4,1),MomExt(1:4,2))
  endif
  PDFFac_R = PDFFac_R * 2d0!  factor two for sampling over pdfs
  PDFFac_L = PDFFac_L * 2d0

!   PDFFac_L = gL_Zpr(up_)**2 * (pdf(up_,1)*pdf(aup_,2) + pdf(Up_,2)*pdf(aup_,1))
!   PDFFac_L = PDFFac_L + gL_Zpr(dn_)**2 * (pdf(dn_,1)*pdf(adn_,2) + pdf(dn_,2)*pdf(adn_,1))
!   PDFFac_L = PDFFac_L + gL_Zpr(chm_)**2 * (pdf(chm_,1)*pdf(achm_,2) + pdf(chm_,2)*pdf(achm_,1))
!   PDFFac_L = PDFFac_L + gL_Zpr(str_)**2 * (pdf(str_,1)*pdf(astr_,2) + pdf(str_,2)*pdf(astr_,1))
!   PDFFac_L = PDFFac_L + gL_Zpr(bot_)**2 * (pdf(bot_,1)*pdf(abot_,2) + pdf(bot_,2)*pdf(abot_,1))
!   PDFFac_R = gR_Zpr(up_)**2 * (pdf(up_,1)*pdf(aup_,2) + pdf(Up_,2)*pdf(aup_,1))
!   PDFFac_R = PDFFac_R + gR_Zpr(dn_)**2 * (pdf(dn_,1)*pdf(adn_,2) + pdf(dn_,2)*pdf(adn_,1))
!   PDFFac_R = PDFFac_R + gR_Zpr(chm_)**2 * (pdf(chm_,1)*pdf(achm_,2) + pdf(chm_,2)*pdf(achm_,1))
!   PDFFac_R = PDFFac_R + gR_Zpr(str_)**2 * (pdf(str_,1)*pdf(astr_,2) + pdf(str_,2)*pdf(astr_,1))
!   PDFFac_R = PDFFac_R + gR_Zpr(bot_)**2 * (pdf(bot_,1)*pdf(abot_,2) + pdf(bot_,2)*pdf(abot_,1))


IF( CORRECTION.EQ.4 ) THEN
   call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,3),yRnd(5:8) ,MomDK(1:4,1:3),PSWgt2)
   call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,4),yRnd(9:12),MomDK(1:4,4:6),PSWgt3)

   call Kinematics_TTBARZprime(.false.,MomExt,MomDK,applyPSCut,NBin)!,xJPsiFrag=xFrag)
   if( applyPSCut ) then
      EvalCS_NLODK_Zprime_ttb = 0d0
      return
   endif

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt
   RunFactor = RunAlphaS(2,MuRen)
   ISFac = MomCrossing(MomExt)

!----------------------------------------
! one loop correction to Anti-top decay |
!----------------------------------------
   NLO_Res_UnPol_Left = (0d0,0d0)
   NLO_Res_UnPol_Right = (0d0,0d0)
   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
      call Tree_Zprime_tbtqbq(LO_Res_Pol_Left, LO_Res_Pol_Right)

      call TopDecay(ExtParticle(1),DK_1L_T,MomDK(1:4,1:3))
      call Tree_Zprime_tbtqbq(NLO_Res_Pol_Left, NLO_Res_Pol_Right)

      !--F By looking at Markus' implementation I assume that the factor 2 in 2*dreal(Alo*conjg(Avirt)) is already taken care of
      NLO_Res_UnPol_Left = NLO_Res_UnPol_Left + dreal(LO_Res_Pol_Left*dconjg(NLO_Res_Pol_Left))
      NLO_Res_UnPol_Right = NLO_Res_UnPol_Right + dreal(LO_Res_Pol_Right*dconjg(NLO_Res_Pol_Right))
   enddo!helicity loop

   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * (PDFFac_L * NLO_Res_Unpol_Left + PDFFac_R * NLO_Res_Unpol_Right)

   NLO_Res_Unpol = NLO_Res_Unpol * 9d0 ! color sum
   EvalCS_NLODK_Zprime_ttb = EvalCS_NLODK_Zprime_ttb + dble(NLO_Res_Unpol)
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo



!-------------------------------------
! one loop correction to top-decay   |
!-------------------------------------
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt
   RunFactor = RunAlphaS(2,MuRen)
   ISFac = MomCrossing(MomExt)
   NLO_Res_UnPol_Left = (0d0,0d0)
   NLO_Res_UnPol_Right = (0d0,0d0)
   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
      call Tree_Zprime_tbtqbq(LO_Res_Pol_Left,LO_Res_Pol_Right)

      call TopDecay(ExtParticle(2),DK_1L_T,MomDK(1:4,4:6))
      call Tree_Zprime_tbtqbq(NLO_Res_Pol_Left,NLO_Res_Pol_Right)

      !--F By looking at Markus' implementation I assume that the factor 2 in 2*dreal(Alo*conjg(Avirt)) is already taken care of
      NLO_Res_UnPol_Left = NLO_Res_UnPol_Left + dreal(LO_Res_Pol_Left*dconjg(NLO_Res_Pol_Left))
      NLO_Res_UnPol_Right = NLO_Res_UnPol_Right + dreal(LO_Res_Pol_Right*dconjg(NLO_Res_Pol_Right))
   enddo!helicity loop

   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * (PDFFac_L * NLO_Res_Unpol_Left + PDFFac_R * NLO_Res_Unpol_Right)
   NLO_Res_Unpol = NLO_Res_Unpol * 9d0 ! color sum
   EvalCS_NLODK_Zprime_ttb = EvalCS_NLODK_Zprime_ttb + dble(NLO_Res_Unpol)
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo




IF( TOPDECAYS.EQ.2 .OR. TOPDECAYS.EQ.4 ) THEN
!----------------------------------------
! one loop correction to W- decay |
!----------------------------------------
   NLO_Res_UnPol_Left = (0d0,0d0)
   NLO_Res_UnPol_Right = (0d0,0d0)
   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
      call Tree_Zprime_tbtqbq(LO_Res_Pol_Left, LO_Res_Pol_Right)

      call TopDecay(ExtParticle(1),DK_1L_Q,MomDK(1:4,1:3))
      call Tree_Zprime_tbtqbq(NLO_Res_Pol_Left, NLO_Res_Pol_Right)

      NLO_Res_UnPol_Left = NLO_Res_UnPol_Left + dreal(LO_Res_Pol_Left*dconjg(NLO_Res_Pol_Left))
      NLO_Res_UnPol_Right = NLO_Res_UnPol_Right + dreal(LO_Res_Pol_Right*dconjg(NLO_Res_Pol_Right))
   enddo!helicity loop

   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * (PDFFac_L * NLO_Res_Unpol_Left + PDFFac_R * NLO_Res_Unpol_Right)
   NLO_Res_Unpol = NLO_Res_Unpol * 9d0 ! color sum
   EvalCS_NLODK_Zprime_ttb = EvalCS_NLODK_Zprime_ttb + dble(NLO_Res_Unpol)
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt
   RunFactor = RunAlphaS(2,MuRen)
   ISFac = MomCrossing(MomExt)
ENDIF


IF( TOPDECAYS.EQ.2 .OR. TOPDECAYS.EQ.3 ) THEN
!-------------------------------------
! one loop correction to W+ decay   |
!-------------------------------------
   NLO_Res_UnPol_Left = (0d0,0d0)
   NLO_Res_UnPol_Right = (0d0,0d0)
   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
      call Tree_Zprime_tbtqbq(LO_Res_Pol_Left,LO_Res_Pol_Right)

      call TopDecay(ExtParticle(2),DK_1L_Q,MomDK(1:4,4:6))
      call Tree_Zprime_tbtqbq(NLO_Res_Pol_Left,NLO_Res_Pol_Right)

      NLO_Res_UnPol_Left = NLO_Res_UnPol_Left + dreal(LO_Res_Pol_Left*dconjg(NLO_Res_Pol_Left))
      NLO_Res_UnPol_Right = NLO_Res_UnPol_Right + dreal(LO_Res_Pol_Right*dconjg(NLO_Res_Pol_Right))
   enddo!helicity loop

   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * (PDFFac_L * NLO_Res_Unpol_Left + PDFFac_R * NLO_Res_Unpol_Right)
   NLO_Res_Unpol = NLO_Res_Unpol * 9d0 ! color sum
   EvalCS_NLODK_Zprime_ttb = EvalCS_NLODK_Zprime_ttb + dble(NLO_Res_Unpol)
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo
ENDIF



ELSEIF( CORRECTION.EQ.5 ) THEN
!------------------------------------
! real correction to Anti-top decay |
!------------------------------------
   call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,3),yRnd(5:11),MomDK(1:4,1:4),PSWgt2) 
   call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,4),yRnd(12:15),MomDK(1:4,5:7),PSWgt3) 
!   write(6,*) dsqrt(dabs(MomExt(1,4)**2-MomExt(2,4)**2-MomExt(3,4)**2-MomExt(4,4)**2))
!   pause

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt
   RunFactor = RunAlphaS(2,MuRen)
   ISFac = MomCrossing(MomExt)


   MomDKx(1:4,1:3) = MomDK(1:4,1:3) 
   MomDKx(1:4,4:6) = MomDK(1:4,5:7) 
   MomDKx(1:4,7) = MomDK(1:4,4) 
   call Kinematics_TTBARZprime(.true.,MomExt,MomDKx,applyPSCut,NBin)
   if( applyPSCut ) goto TheDipole_ATop

   NLO_Res_UnPol_Left = (0d0,0d0)
   NLO_Res_UnPol_Right = (0d0,0d0)
   do iHel=1,NumHelicities
   do GluHel=-1,1,2 
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_RE_T,MomDK(1:4,1:4),GluonHel=GluHel) 
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,5:7)) 
      call Tree_Zprime_tbtqbq(NLO_Res_Pol_Left,NLO_Res_Pol_Right)

      NLO_Res_UnPol_Left = NLO_Res_UnPol_Left + NLO_Res_Pol_Left*dconjg(NLO_Res_Pol_Left)
      NLO_Res_UnPol_Right = NLO_Res_UnPol_Right + NLO_Res_Pol_Right*dconjg(NLO_Res_Pol_Right)
   enddo!helicity loop 
   enddo!helicity loop

   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * (PDFFac_L * NLO_Res_Unpol_Left + PDFFac_R * NLO_Res_Unpol_Right) 
   NLO_Res_Unpol = NLO_Res_Unpol * 9d0 ! color sum
   EvalCS_NLODK_Zprime_ttb = EvalCS_NLODK_Zprime_ttb + dble(NLO_Res_Unpol)
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo
!print *, "real  ", dble(NLO_Res_Unpol)


TheDipole_ATop continue

   call WTransform(MomDK(1:4,1:4),MomDKTd(1:4,1:3),pbDpg,ptDpg,ptDpb)
   MomDKTd(1:4,4:6) = MomDK(1:4,5:7)
   omz=ptDpg/(ptDpb+ptDpg-pbDpg)
   rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
   z=1d0-omz
   y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2
   Dipole = - alpha_s*4d0*Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
   Dipole = Dipole * (1d0 - StepFunc(1d0-alpha_DK-z) * StepFunc(y-alpha_DK*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )

   call Kinematics_TTBARZprime(.false.,MomExt,MomDKTd,applyPSCut,NBin)
   if( applyPSCut ) goto TheReal_Top
   NLO_Res_UnPol_Left = (0d0,0d0)
   NLO_Res_UnPol_Right = (0d0,0d0)
   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomDKTd(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_LO,MomDKTd(1:4,4:6))
      call Tree_Zprime_tbtqbq(NLO_Res_Pol_Left,NLO_Res_Pol_Right)

      NLO_Res_UnPol_Left = NLO_Res_UnPol_Left + NLO_Res_Pol_Left*dconjg(NLO_Res_Pol_Left)
      NLO_Res_UnPol_Right = NLO_Res_UnPol_Right + NLO_Res_Pol_Right*dconjg(NLO_Res_Pol_Right)
   enddo!helicity loop

   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * (PDFFac_L * NLO_Res_Unpol_Left + PDFFac_R * NLO_Res_Unpol_Right) * Dipole
   NLO_Res_Unpol = NLO_Res_Unpol * 9d0 ! color sum
   EvalCS_NLODK_Zprime_ttb = EvalCS_NLODK_Zprime_ttb + dble(NLO_Res_Unpol) 
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo
! print *, "dipole", dble(NLO_Res_Unpol)
! print *, "sing", (MomDK(1:4,1).dot.MomDK(1:4,4))/m_Top**2
! pause

TheReal_Top continue

!------------------------------------
! real correction to Top decay |
!------------------------------------
   call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,3),yRnd(5:8),MomDK(1:4,1:3),PSWgt2)
   call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,4),yRnd(9:15),MomDK(1:4,4:7),PSWgt3)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt
   RunFactor = RunAlphaS(2,MuRen)
   ISFac = MomCrossing(MomExt)
   call Kinematics_TTBARZprime(.true.,MomExt,MomDK(1:4,1:7),applyPSCut,NBin)
   if( applyPSCut ) goto TheDipole_Top

   NLO_Res_UnPol_Left = (0d0,0d0)
   NLO_Res_UnPol_Right = (0d0,0d0)
   do iHel=1,NumHelicities
   do GluHel=-1,1,2
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_RE_T,MomDK(1:4,4:7),GluonHel=GluHel)
      call Tree_Zprime_tbtqbq(NLO_Res_Pol_Left,NLO_Res_Pol_Right)

      NLO_Res_UnPol_Left = NLO_Res_UnPol_Left + NLO_Res_Pol_Left*dconjg(NLO_Res_Pol_Left)
      NLO_Res_UnPol_Right = NLO_Res_UnPol_Right + NLO_Res_Pol_Right*dconjg(NLO_Res_Pol_Right)
   enddo!helicity loop
   enddo!helicity loop


   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * (PDFFac_L * NLO_Res_Unpol_Left + PDFFac_R * NLO_Res_Unpol_Right)
   NLO_Res_Unpol = NLO_Res_Unpol * 9d0 ! color sum
   EvalCS_NLODK_Zprime_ttb = EvalCS_NLODK_Zprime_ttb + dble(NLO_Res_Unpol) 
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo
!   print *, 'real  ', dble(NLO_Res_Unpol)

TheDipole_Top continue

   call WTransform(MomDK(1:4,4:7),MomDKTd(1:4,4:6),pbDpg,ptDpg,ptDpb)
   MomDKTd(1:4,1:3) = MomDK(1:4,1:3)
   omz=ptDpg/(ptDpb+ptDpg-pbDpg)
   rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
   z=1d0-omz
   y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2
   Dipole = - alpha_s*4d0*Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
   Dipole = Dipole * (1d0 - StepFunc(1d0-alpha_DK-z) * StepFunc(y-alpha_DK*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )

   call Kinematics_TTBARZprime(.false.,MomExt,MomDKTd,applyPSCut,NBin)
   if( applyPSCut ) goto TheWm_decay
   NLO_Res_UnPol_Left = (0d0,0d0)
   NLO_Res_UnPol_Right = (0d0,0d0)
   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomDKTd(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_LO,MomDKTd(1:4,4:6))
      call Tree_Zprime_tbtqbq(NLO_Res_Pol_Left,NLO_Res_Pol_Right)

      NLO_Res_UnPol_Left = NLO_Res_UnPol_Left + NLO_Res_Pol_Left*dconjg(NLO_Res_Pol_Left)
      NLO_Res_UnPol_Right = NLO_Res_UnPol_Right + NLO_Res_Pol_Right*dconjg(NLO_Res_Pol_Right)
   enddo!helicity loop

   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * (PDFFac_L * NLO_Res_Unpol_Left + PDFFac_R * NLO_Res_Unpol_Right) * Dipole
   NLO_Res_Unpol = NLO_Res_Unpol * 9d0 ! color sum
   EvalCS_NLODK_Zprime_ttb = EvalCS_NLODK_Zprime_ttb + dble(NLO_Res_Unpol) 
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo
! print *, "dipole", dble(NLO_Res_Unpol)
! print *, "sing", (MomDK(1:4,4).dot.MomDK(1:4,7))/m_Top**2
! pause

TheWm_decay continue

IF (TOPDECAYS.EQ.2 .OR. TOPDECAYS.EQ.4) THEN
!------------------------------
! real correction to W- decay |
!------------------------------
   call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,3),yRnd(5:11),MomDK(1:4,1:4),PSWgt2) 
   call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,4),yRnd(12:15),MomDK(1:4,5:7),PSWgt3) 
!   write(6,*) dsqrt(dabs(MomExt(1,4)**2-MomExt(2,4)**2-MomExt(3,4)**2-MomExt(4,4)**2))
!   pause

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt
   RunFactor = RunAlphaS(2,MuRen)
   ISFac = MomCrossing(MomExt)

   MomDKx(1:4,1:3) = MomDK(1:4,1:3) 
   MomDKx(1:4,4:6) = MomDK(1:4,5:7) 
   MomDKx(1:4,7) = MomDK(1:4,4) 
   call Kinematics_TTBARZprime(.true.,MomExt,MomDKx,applyPSCut,NBin)
   if( applyPSCut ) goto TheDipole_ATop_Wdecay

   NLO_Res_UnPol_Left = (0d0,0d0)
   NLO_Res_UnPol_Right = (0d0,0d0)
   do iHel=1,NumHelicities
   do GluHel=-1,1,2 
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_RE_Q,MomDK(1:4,1:4),GluonHel=GluHel) 
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,5:7)) 
      call Tree_Zprime_tbtqbq(NLO_Res_Pol_Left,NLO_Res_Pol_Right)

      NLO_Res_UnPol_Left = NLO_Res_UnPol_Left + NLO_Res_Pol_Left*dconjg(NLO_Res_Pol_Left)
      NLO_Res_UnPol_Right = NLO_Res_UnPol_Right + NLO_Res_Pol_Right*dconjg(NLO_Res_Pol_Right)
   enddo!helicity loop 
   enddo!helicity loop

   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * (PDFFac_L * NLO_Res_Unpol_Left + PDFFac_R * NLO_Res_Unpol_Right) 
   NLO_Res_Unpol = NLO_Res_Unpol * 9d0 ! color sum
   EvalCS_NLODK_Zprime_ttb = EvalCS_NLODK_Zprime_ttb + dble(NLO_Res_Unpol)

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo

!print *, "real              ", dble(NLO_Res_Unpol)

TheDipole_ATop_Wdecay continue

do ndip=1,2
   call wdec_trans(ndip,MomDK(1:4,1:4),MomDKTd(1:4,1:3),alpha_DKWff,dip_res_w)
   if( dip_res_w.eq.0d0 ) goto 20
   Dipole = - alpha_s4Pi*RunFactor * CF * dip_res_w

   MomDKTd(1:4,4:6) = MomDK(1:4,5:7)

   call Kinematics_TTBARZprime(.false.,MomExt,MomDKTd,applyPSCut,NBin)
   if( applyPSCut ) then
      goto 20
      return
   endif

   NLO_Res_UnPol_Left = (0d0,0d0)
   NLO_Res_UnPol_Right = (0d0,0d0)
   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomDKTd(1:4,1:3)) 
      call TopDecay(ExtParticle(2),DK_LO,MomDKTd(1:4,4:6)) 
      call Tree_Zprime_tbtqbq(NLO_Res_Pol_Left,NLO_Res_Pol_Right)

      NLO_Res_UnPol_Left = NLO_Res_UnPol_Left + NLO_Res_Pol_Left*dconjg(NLO_Res_Pol_Left)
      NLO_Res_UnPol_Right = NLO_Res_UnPol_Right + NLO_Res_Pol_Right*dconjg(NLO_Res_Pol_Right)
   enddo!helicity loop

   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * (PDFFac_L * NLO_Res_Unpol_Left + PDFFac_R * NLO_Res_Unpol_Right) 
   NLO_Res_Unpol = NLO_Res_Unpol * 9d0 ! color sum
   NLO_Res_Unpol = NLO_Res_Unpol * Dipole
   EvalCS_NLODK_Zprime_ttb = EvalCS_NLODK_Zprime_ttb + dble(NLO_Res_Unpol)

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo

!   print *, "dipole", ndip, dble(NLO_Res_Unpol)
!   print *, "sing", (MomDK(1:4,2).dot.MomDK(1:4,4))/m_Top**2

20 continue
enddo ! dipole loop

ENDIF ! Wm decay real correction


TheWp_decay continue

IF (TOPDECAYS.EQ.2 .OR. TOPDECAYS.EQ.3) THEN
!------------------------------
! real correction to W+ decay |
!------------------------------
   call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,3),yRnd(5:8),MomDK(1:4,1:3),PSWgt2)
   call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,4),yRnd(9:15),MomDK(1:4,4:7),PSWgt3)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt
   RunFactor = RunAlphaS(2,MuRen)
   ISFac = MomCrossing(MomExt)
   call Kinematics_TTBARZprime(.true.,MomExt,MomDK(1:4,1:7),applyPSCut,NBin)
   if( applyPSCut ) goto TheDipole_Top_Wdecay

   NLO_Res_UnPol_Left = (0d0,0d0)
   NLO_Res_UnPol_Right = (0d0,0d0)
   do iHel=1,NumHelicities
   do GluHel=-1,1,2
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_RE_Q,MomDK(1:4,4:7),GluonHel=GluHel)
      call Tree_Zprime_tbtqbq(NLO_Res_Pol_Left,NLO_Res_Pol_Right)

      NLO_Res_UnPol_Left = NLO_Res_UnPol_Left + NLO_Res_Pol_Left*dconjg(NLO_Res_Pol_Left)
      NLO_Res_UnPol_Right = NLO_Res_UnPol_Right + NLO_Res_Pol_Right*dconjg(NLO_Res_Pol_Right)
   enddo!helicity loop
   enddo!helicity loop


   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * (PDFFac_L * NLO_Res_Unpol_Left + PDFFac_R * NLO_Res_Unpol_Right)
   NLO_Res_Unpol = NLO_Res_Unpol * 9d0 ! color sum
   EvalCS_NLODK_Zprime_ttb = EvalCS_NLODK_Zprime_ttb + dble(NLO_Res_Unpol) 
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo

!print *, "real              ", dble(NLO_Res_Unpol)

TheDipole_Top_Wdecay continue

do ndip=1,2
   call wdec_trans(ndip,MomDK(1:4,4:7),MomDKTd(1:4,4:6),alpha_DKWff,dip_res_w)
   if( dip_res_w.eq.0d0 ) goto 22
   Dipole = - alpha_s4Pi*RunFactor * CF * dip_res_w

   MomDKTd(1:4,1:3) = MomDK(1:4,1:3)

   call Kinematics_TTBARZprime(.false.,MomExt,MomDKTd,applyPSCut,NBin)
   if( applyPSCut ) then
      goto 22
      return
   endif

   NLO_Res_UnPol_Left = (0d0,0d0)
   NLO_Res_UnPol_Right = (0d0,0d0)
   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomDKTd(1:4,1:3)) 
      call TopDecay(ExtParticle(2),DK_LO,MomDKTd(1:4,4:6)) 
      call Tree_Zprime_tbtqbq(NLO_Res_Pol_Left,NLO_Res_Pol_Right)

      NLO_Res_UnPol_Left = NLO_Res_UnPol_Left + NLO_Res_Pol_Left*dconjg(NLO_Res_Pol_Left)
      NLO_Res_UnPol_Right = NLO_Res_UnPol_Right + NLO_Res_Pol_Right*dconjg(NLO_Res_Pol_Right)
   enddo!helicity loop

   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * (PDFFac_L * NLO_Res_Unpol_Left + PDFFac_R * NLO_Res_Unpol_Right) 
   NLO_Res_Unpol = NLO_Res_Unpol * 9d0 ! color sum
   NLO_Res_Unpol = NLO_Res_Unpol * Dipole
   EvalCS_NLODK_Zprime_ttb = EvalCS_NLODK_Zprime_ttb + dble(NLO_Res_Unpol)

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo

! print *, "dipole", ndip, dble(NLO_Res_Unpol)
! print *, "sing", (MomDK(1:4,4).dot.MomDK(1:4,7))/m_Top**2

22 continue
enddo ! dipole loop

! pause

ENDIF ! Wp decay real correction

TheEnd continue


ENDIF

   EvalCS_NLODK_Zprime_ttb = EvalCS_NLODK_Zprime_ttb/VgsWgt
RETURN
END FUNCTION


END MODULE ModCrossSection_ZprimeTTB


