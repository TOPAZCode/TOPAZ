module modIntDipoles_QGSTSTBQ
use ModTopdecay
implicit none

public:: EvalIntDipoles_QGSTSTBQ

integer, parameter,private  :: dp = selected_real_kind(15)
real(dp),parameter,private :: zero = 0.0_dp

double precision, private :: MomExt(1:4,1:14)


contains






      SUBROUTINE EvalIntDipoles_QGSTSTBQ(p,pDK,z,res)
      use modParameters
      use modMisc
      use ModIntDipoles
      implicit none
      real(dp), intent(in) :: p(4,4),pDK(1:4,1:10)
      real(dp), intent(out) :: res(1:3)
      real(dp) :: dipsoft,dipfini,dipplus,mtrsq,AP(1:3),epcorr
      real(dp) :: Tree_21,Tree_13,Tree_14,Tree_23,Tree_24,Tree_12,z
      real(dp) :: CA,CF,TR,L,Q2,tmp
      integer  :: n,emi,in1,in2
      complex(dp) :: TreeMom(1:4,1:4)

       MomExt(1:4,5:14) = pDK(1:4,1:10)

       res(1:3) = zero
       CF=4d0/3d0
       CA=3d0
       TR=0.5d0


!  tree momenta for g g -> t tb  &  q qb -> t tb
   TreeMom(1:4,1) =-dcmplx( p(1:4,3) )
   TreeMom(1:4,2) =-dcmplx( p(1:4,4) )
   TreeMom(1:4,3) = dcmplx( p(1:4,2) )
   TreeMom(1:4,4) = dcmplx( p(1:4,1) )

   Tree_12 = Tree_GG_TTb_00(0,TreeMom,(/0d0,0d0,m_stop**2,m_stop**2/))
   Tree_21 = Tree_UUb_TTb_00(0,TreeMom,(/0d0,0d0,m_stop**2,m_stop**2/))

   do n=5,6   ! implementation for ii-contributions only
      if(n.eq.5) then
        dipsoft =ii_qg(zero,zero,p,3,4,z,1)
        dipfini =ii_qg(zero,zero,p,3,4,z,2)
        dipplus =ii_qg(zero,zero,p,3,4,z,3)
        emi = 1 ! mom #3 is emitting
        mtrsq = Tree_12*(-CF)
      endif
      if(n.eq.6) then
        dipsoft =ii_gq(zero,zero,p,4,3,z,1)
        dipfini =ii_gq(zero,zero,p,4,3,z,2)
        dipplus =ii_gq(zero,zero,p,4,3,z,3)
        emi = 2 ! mom #4 is emitting
        mtrsq = Tree_21   !*(-1d0)  ! removed this minus sign!  ! TR is included in splitting functions
      endif

      if(emi.eq.1) then
        res(1) = res(1) + (dipsoft-dipplus)*mtrsq
        res(2) = res(2) + (dipfini+dipplus)*mtrsq
      endif
      if(emi.eq.2) then
        res(1) = res(1) + (dipsoft-dipplus)*mtrsq
        res(3) = res(3) + (dipfini+dipplus)*mtrsq
      endif

   enddo
   res(1:3) = alpha_sOver2Pi * res(1:3)!    removed a minus sign here


! !        epcorr=epinv+2d0*dlog(renscale/facscale)
       mtrsq = Tree_21
       epcorr=epinv
       AP(1)= 0d0
       AP(2)= TR*(z**2+(1d0-z)**2) * mtrsq * alpha_sOver2Pi *epcorr
       AP(3)= 0d0
       res(1) = res(1) + (AP(1)-AP(3))
       res(3) = res(3) + (AP(2)+AP(3))

       mtrsq = Tree_12
       AP(1)= 0d0
       AP(2)= CF * (1d0+(1d0-z)**2)/z * mtrsq * alpha_sOver2Pi *epcorr
       AP(3)= 0d0
       res(1) = res(1) + (AP(1)-AP(3))
       res(2) = res(2) + (AP(2)+AP(3))


  RETURN
  END SUBROUTINE












      FUNCTION Tree_UUb_TTb_00(icorr,MomTd,Mass2Td) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      use modKinematics
      use ModExoticDecay
      implicit none
      double complex SqAmp,Amp(1,0:1),Res(0:3)
      integer col,colP,iHel,icorr,ngl,A0barHel,A0Hel,SecHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      double precision :: PSWgt2,PSWgt3,PSWgt4,PSWgt5
      type(Particle),target :: TopQuark(1:2)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:1),Quarks(1:4)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(9)
      real(8)::ColCorr(1:1,1:1)
!       real(8),parameter::ColCorr1(1:1,1:1)=(/8.D0/3.D0/)
!       real(8),parameter::ColCorr2(1:1,1:1)=(/-56.D0/3.D0/)
!       real(8),parameter::ColCorr3(1:1,1:1)=(/-16.D0/3.D0/)
!       real(8),parameter::ColCorr4(1:1,1:1)=(/8.D0/3.D0/)
!       real(8),parameter::ColCorr5(1:1,1:1)=(/-16.D0/3.D0/)
!       real(8),parameter::ColCorr6(1:1,1:1)=(/-56.D0/3.D0/)
!       real(8),parameter::ColCorr7(1:1,1:1)=(/-56.D0/3.D0/)
!       real(8),parameter::ColCorr8(1:1,1:1)=(/-16.D0/3.D0/)
!       real(8),parameter::ColCorr9(1:1,1:1)=(/8.D0/3.D0/)
!       real(8),parameter::ColCorr10(1:1,1:1)=(/-16.D0/3.D0/)
!       real(8),parameter::ColCorr11(1:1,1:1)=(/-56.D0/3.D0/)
!       real(8),parameter::ColCorr12(1:1,1:1)=(/8.D0/3.D0/)
      integer,target :: HelList(1:16,1:4)
      double complex :: MomTmp(0:3)
      logical :: applyPSCut

      SqAmp = 0d0
!       if(icorr.eq.1) ColCorr(1:1,1:1)=ColCorr1(1:1,1:1)
!       if(icorr.eq.2) ColCorr(1:1,1:1)=ColCorr2(1:1,1:1)
!       if(icorr.eq.3) ColCorr(1:1,1:1)=ColCorr3(1:1,1:1)
!       if(icorr.eq.4) ColCorr(1:1,1:1)=ColCorr4(1:1,1:1)
!       if(icorr.eq.5) ColCorr(1:1,1:1)=ColCorr5(1:1,1:1)
!       if(icorr.eq.6) ColCorr(1:1,1:1)=ColCorr6(1:1,1:1)
!       if(icorr.eq.7) ColCorr(1:1,1:1)=ColCorr7(1:1,1:1)
!       if(icorr.eq.8) ColCorr(1:1,1:1)=ColCorr8(1:1,1:1)
!       if(icorr.eq.9) ColCorr(1:1,1:1)=ColCorr9(1:1,1:1)
!       if(icorr.eq.10) ColCorr(1:1,1:1)=ColCorr10(1:1,1:1)
!       if(icorr.eq.11) ColCorr(1:1,1:1)=ColCorr11(1:1,1:1)
!       if(icorr.eq.12) ColCorr(1:1,1:1)=ColCorr12(1:1,1:1)

      ColCorr(1:1,1:1)=8d0



      HelList(1,1:4)=(/0,0,0,0/)
      HelList(2,1:4)=(/0,1,0,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

      MomExt(1:4,1)=MomTd(0:3,1)
      MomExt(1:4,2)=MomTd(0:3,2)
      MomExt(1:4,3)=MomTd(0:3,4)
      MomExt(1:4,4)=MomTd(0:3,3)

      TopQuark(1)%PartType = STop_
      TopQuark(1)%Mom(1:4) = MomTd(0:3,3)
      TopQuark(2)%PartType = ASTop_
      TopQuark(2)%Mom(1:4) = MomTd(0:3,4)


!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)



!      set particles
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),+1,PolV(0:3,plus,1))
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),-1,PolV(0:3,minus,1))
      Quarks(1)%Mom => MomTd(:,1)
      Quarks(1)%PartType => AUp_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(1)
      Quarks(1)%Mass  => MassTd(1)

      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),+1,PolV(0:3,plus,2))
      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),-1,PolV(0:3,minus,2))
      Quarks(2)%Mom => MomTd(:,2)
      Quarks(2)%PartType => Up_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(2)
      Quarks(2)%Mass  => MassTd(2)

      Quarks(3)%Mom => MomTd(:,3)
      Quarks(3)%PartType => STop_
      Quarks(3)%ExtRef => ExtRef
      Quarks(3)%Mass2 => Mass2Td(3)
      Quarks(3)%Mass  => MassTd(3)

      Quarks(4)%Mom => MomTd(:,4)
      Quarks(4)%PartType => ASTop_
      Quarks(4)%ExtRef => ExtRef
      Quarks(4)%Mass2 => Mass2Td(4)
      Quarks(4)%Mass  => MassTd(4)


!      sum over helicities
      if( XTopDecays.eq.0 ) then
            SecHel=1
      else
            SecHel=-1
      endif


      SqAmp = (0d0,0d0)
      do A0barHel=1,SecHel,-2
      do A0Hel=1,SecHel,-2


      call STopDecay(TopQuark(2),DKX_STChi0_LO,A0barHel,MomExt(1:4,5:9))
      call STopDecay(TopQuark(1),DKX_STChi0_LO,A0Hel,MomExt(1:4,10:14))


      Quarks(3)%Pol => TopQuark(1)%Pol
      Quarks(4)%Pol => TopQuark(2)%Pol


      do iHel=1,2
      Quarks(2)%Pol => PolV(:,HelList(iHel,2),2)
      Quarks(2)%Helicity => HelList(iHel,2)


!      calc currents
      Res(0:3) = cur_f_ffss((/Gluons(1)/),Quarks(3:4),Quarks(2:2),(/0,0,0,0,0/))
      Amp(1,plus)  = psp1_(Res(0:3),PolV(0:3,plus,1))
      Amp(1,minus) = psp1_(Res(0:3),PolV(0:3,minus,1))



      do col =1,1
      do colP=1,1
        SqAmp = SqAmp + ColCorr(colP,col) * ( dconjg(Amp(colP,plus))*Amp(col,plus) + dconjg(Amp(colP,minus))*Amp(col,minus) )
      enddo
      enddo

      enddo
      enddo
      enddo

      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

      return
      END FUNCTION




      FUNCTION Tree_GG_TTb_00(icorr,MomTd,Mass2Td) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      use modKinematics
      use ModExoticDecay
      implicit none
      double complex SqAmp,Amp(2,0:1),Res(0:3)
      integer col,colP,iHel,icorr,ngl,A0barHel,A0Hel,SecHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      double precision :: PSWgt2,PSWgt3,PSWgt4,PSWgt5
      type(Particle),target :: TopQuark(1:2)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      real(8)::ColCorr(1:2,1:2)
      real(8),parameter::ColCorr0(1:2,1:2)=(/ 64.D0/3.D0, -8.D0/3.D0, -8.D0/3.D0, 64.D0/3.D0 /)
!       real(8),parameter::ColCorr1(1:2,1:2)=(/-72.D0,0.D0,0.D0,-72.D0/)
!       real(8),parameter::ColCorr2(1:2,1:2)=(/-64.D0,8.D0,8.D0,8.D0/)
!       real(8),parameter::ColCorr3(1:2,1:2)=(/8.D0,8.D0,8.D0,-64.D0/)
!       real(8),parameter::ColCorr4(1:2,1:2)=(/-72.D0,0.D0,0.D0,-72.D0/)
!       real(8),parameter::ColCorr5(1:2,1:2)=(/8.D0,8.D0,8.D0,-64.D0/)
!       real(8),parameter::ColCorr6(1:2,1:2)=(/-64.D0,8.D0,8.D0,8.D0/)
!       real(8),parameter::ColCorr7(1:2,1:2)=(/-64.D0,8.D0,8.D0,8.D0/)
!       real(8),parameter::ColCorr8(1:2,1:2)=(/8.D0,8.D0,8.D0,-64.D0/)
!       real(8),parameter::ColCorr9(1:2,1:2)=(/-8.D0/9.D0,-80.D0/9.D0,-80.D0/9.D0,-8.D0/9.D0/)
!       real(8),parameter::ColCorr10(1:2,1:2)=(/8.D0,8.D0,8.D0,-64.D0/)
!       real(8),parameter::ColCorr11(1:2,1:2)=(/-64.D0,8.D0,8.D0,8.D0/)
!       real(8),parameter::ColCorr12(1:2,1:2)=(/-8.D0/9.D0,-80.D0/9.D0,-80.D0/9.D0,-8.D0/9.D0/)
      integer,target :: HelList(1:4,1:3)
      double complex :: MomTmp(0:3)
      logical :: applyPSCut

      SqAmp = 0d0
!       if(icorr.eq.1) ColCorr(1:2,1:2)=ColCorr1(1:2,1:2)
!       if(icorr.eq.2) ColCorr(1:2,1:2)=ColCorr2(1:2,1:2)
!       if(icorr.eq.3) ColCorr(1:2,1:2)=ColCorr3(1:2,1:2)
!       if(icorr.eq.4) ColCorr(1:2,1:2)=ColCorr4(1:2,1:2)
!       if(icorr.eq.5) ColCorr(1:2,1:2)=ColCorr5(1:2,1:2)
!       if(icorr.eq.6) ColCorr(1:2,1:2)=ColCorr6(1:2,1:2)
!       if(icorr.eq.7) ColCorr(1:2,1:2)=ColCorr7(1:2,1:2)
!       if(icorr.eq.8) ColCorr(1:2,1:2)=ColCorr8(1:2,1:2)
!       if(icorr.eq.9) ColCorr(1:2,1:2)=ColCorr9(1:2,1:2)
!       if(icorr.eq.10) ColCorr(1:2,1:2)=ColCorr10(1:2,1:2)
!       if(icorr.eq.11) ColCorr(1:2,1:2)=ColCorr11(1:2,1:2)
!       if(icorr.eq.12) ColCorr(1:2,1:2)=ColCorr12(1:2,1:2)

      ColCorr(1:2,1:2)=ColCorr0(1:2,1:2)


      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))


      MomExt(1:4,1)=MomTd(0:3,1)
      MomExt(1:4,2)=MomTd(0:3,2)
      MomExt(1:4,3)=MomTd(0:3,4)
      MomExt(1:4,4)=MomTd(0:3,3)

      TopQuark(1)%PartType = STop_
      TopQuark(1)%Mom(1:4) = MomTd(0:3,3)
      TopQuark(2)%PartType = ASTop_
      TopQuark(2)%Mom(1:4) = MomTd(0:3,4)


!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

      if( icorr.eq.4 .or. icorr.eq.5 .or. icorr.eq.6 .or. icorr.eq.10 .or. icorr.eq.11.or. icorr.eq.12 ) then
        ngl=1
      else
        ngl=2
      endif


!      set particles
      call pol_mless(MomTd(0:3,1),+1,PolV(0:3,plus,1))
      call pol_mless(MomTd(0:3,1),-1,PolV(0:3,minus,1))

      call pol_mless(MomTd(0:3,2),+1,PolV(0:3,plus,2))
      call pol_mless(MomTd(0:3,2),-1,PolV(0:3,minus,2))
      Gluons(ngl)%Mom => MomTd(:,ngl)
      Gluons(ngl)%PartType => Glu_
      Gluons(ngl)%ExtRef => ExtRef
      Gluons(ngl)%Mass2 => Mass2Td(ngl)
      Gluons(ngl)%Mass  => MassTd(ngl)

      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => STop_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ASTop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)


!      sum over helicities
      if( XTopDecays.eq.0 ) then
            SecHel=1
      else
            SecHel=-1
      endif

      SqAmp = (0d0,0d0)
      do A0barHel=1,SecHel,-2
      do A0Hel=1,SecHel,-2

            call STopDecay(TopQuark(2),DKX_STChi0_LO,A0barHel,MomExt(1:4,5:9))
            call STopDecay(TopQuark(1),DKX_STChi0_LO,A0Hel,MomExt(1:4,10:14))


      Quarks(1)%Pol => TopQuark(1)%Pol
      Quarks(2)%Pol => TopQuark(2)%Pol

      do iHel=1,2
      Gluons(ngl)%Pol => PolV(:,HelList(iHel,1),ngl)
      Gluons(ngl)%Helicity => HelList(iHel,1)

!      calc currents
      Res(0:3) = cur_g_2s( (/Gluons(ngl)/),(/Quarks(2),Quarks(1)/),(/2,1,0,0/) )
      Amp(3-ngl,plus)  = Res(0:3).dot.PolV(0:3,plus,3-ngl)
      Amp(3-ngl,minus) = Res(0:3).dot.PolV(0:3,minus,3-ngl)

      Res(0:3) = cur_g_2s( (/Gluons(ngl)/),(/Quarks(2),Quarks(1)/),(/2,0,0,1/) )
      Amp(ngl,plus)  = Res(0:3).dot.PolV(0:3,plus,3-ngl)
      Amp(ngl,minus) = Res(0:3).dot.PolV(0:3,minus,3-ngl)

!       EpsDotP = PolV(0:3,plus,3-ngl).dot.Split_V(0:3)
      SpinCorr_plus_plus  =-1d0
      SpinCorr_plus_minus = 0d0
      do col =1,2
      do colP=1,2
        SqAmp = SqAmp + ColCorr(colP,col) *(       SpinCorr_plus_plus * ( dconjg(Amp(colP,plus)) *Amp(col,plus)     &
                                                                        + dconjg(Amp(colP,minus))*Amp(col,minus) )  &
                                           +        SpinCorr_plus_minus * dconjg(Amp(colP,plus)) *Amp(col,minus)    &
                                           + dconjg(SpinCorr_plus_minus)* dconjg(Amp(colP,minus))*Amp(col,plus)  )
      enddo
      enddo

      enddo
      enddo
      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg *SqAmp



      return
      END FUNCTION











END MODULE
