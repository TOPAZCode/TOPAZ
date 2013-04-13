module modIntDipoles_GGHTHTBG
use ModTopDecay
implicit none
private

public:: EvalIntDipoles_GGHTHTBG,Tree_GG_TTb_ij

integer, parameter,private  :: dp = selected_real_kind(15)
real(dp),parameter,private :: zero = 0.0_dp

double precision, private :: MomDK(1:4,5:14)

contains






      SUBROUTINE EvalIntDipoles_GGHTHTBG(p,pDK,z,res)
      use modParameters
      use modMisc
      use ModIntDipoles
      implicit none
      real(dp), intent(in) :: p(4,4),pDK(1:4,1:10)
      real(dp), intent(out) :: res(1:3)
      real(dp) :: dipsoft,dipfini,dipplus,mtrsq,AP(1:3),epcorr
      real(dp) :: Tree_ij(0:6),Tree_12,Tree_13,Tree_14,Tree_23,Tree_24,Tree_34,z
      real(dp) :: CA,L,Q2
      integer  :: n,emi,in1,in2
      complex(dp) :: TreeMom(1:4,1:4)

       res(1:3) = zero
       CA=3d0

       MomDK(1:4,5:14) = pDK(1:4,1:10)

!  tree momenta for g g -> t tb
    TreeMom(1:4,1) =-dcmplx( p(1:4,3) )
    TreeMom(1:4,2) =-dcmplx( p(1:4,4) )
    TreeMom(1:4,3) = dcmplx( p(1:4,2) )
    TreeMom(1:4,4) = dcmplx( p(1:4,1) )


    Tree_ij= Tree_GG_TTb_ij(TreeMom,(/0d0,0d0,m_HTop**2,m_HTop**2/))
    Tree_12  =  Tree_ij(1) * 0.5d0
    Tree_13  =  Tree_ij(2) * 0.5d0
    Tree_14  =  Tree_ij(3) * 0.5d0
    Tree_23  =  Tree_ij(4) * 0.5d0
    Tree_24  =  Tree_ij(5) * 0.5d0
    Tree_34  =  Tree_ij(6) * 0.5d0

   do n=1,12
      if(n.eq.1) then
!       {if}=initial emitter, final spectator, {gg}=gluon parent, gluon not emitted, {zero,m_HTop}=parent,spectator mass, {3,1}=parent, spectator in p
        dipsoft =if_gg(zero,m_HTop,p,3,1,z,1)
        dipfini =if_gg(zero,m_HTop,p,3,1,z,2)
        dipplus =if_gg(zero,m_HTop,p,3,1,z,3)
        emi = 1 ! mom #3 is emitting
        mtrsq = Tree_14   ! (1-4)-color corr. for (1-3)-dipole because t and tbar are exchanged
      endif
      if(n.eq.2) then
        dipsoft =if_gg(zero,m_HTop,p,3,2,z,1)
        dipfini =if_gg(zero,m_HTop,p,3,2,z,2)
        dipplus =if_gg(zero,m_HTop,p,3,2,z,3)
        emi = 1 ! mom #3 is emitting
        mtrsq = Tree_13
      endif
      if(n.eq.3) then
        dipsoft =if_gg(zero,m_HTop,p,4,1,z,1)
        dipfini =if_gg(zero,m_HTop,p,4,1,z,2)
        dipplus =if_gg(zero,m_HTop,p,4,1,z,3)
        emi = 2 ! mom #4 is emitting
        mtrsq = Tree_24
      endif
      if(n.eq.4) then
        dipsoft =if_gg(zero,m_HTop,p,4,2,z,1)
        dipfini =if_gg(zero,m_HTop,p,4,2,z,2)
        dipplus =if_gg(zero,m_HTop,p,4,2,z,3)
        emi = 2 ! mom #4 is emitting
        mtrsq = Tree_23
      endif
      if(n.eq.5) then
        dipsoft =ii_gg(zero,zero,p,3,4,z,1)
        dipfini =ii_gg(zero,zero,p,3,4,z,2)
        dipplus =ii_gg(zero,zero,p,3,4,z,3)
        emi = 1 ! mom #3 is emitting
        mtrsq = Tree_12
      endif
      if(n.eq.6) then
        dipsoft =ii_gg(zero,zero,p,4,3,z,1)
        dipfini =ii_gg(zero,zero,p,4,3,z,2)
        dipplus =ii_gg(zero,zero,p,4,3,z,3)
        emi = 2 ! mom #4 is emitting
        mtrsq = Tree_12
      endif
      if(n.eq.7) then
        dipsoft =fi_qq(m_HTop,zero,p,1,3,z,1)
        dipfini =fi_qq(m_HTop,zero,p,1,3,z,2)
        dipplus =fi_qq(m_HTop,zero,p,1,3,z,3)
        emi = 1
        mtrsq = Tree_14
      endif
      if(n.eq.8) then
        dipsoft =fi_qq(m_HTop,zero,p,1,4,z,1)
        dipfini =fi_qq(m_HTop,zero,p,1,4,z,2)
        dipplus =fi_qq(m_HTop,zero,p,1,4,z,3)
        emi = 2
        mtrsq = Tree_24
      endif
      if(n.eq.9) then
        dipsoft =fi_qq(m_HTop,zero,p,2,3,z,1)
        dipfini =fi_qq(m_HTop,zero,p,2,3,z,2)
        dipplus =fi_qq(m_HTop,zero,p,2,3,z,3)
        emi = 1
        mtrsq = Tree_13
      endif
      if(n.eq.10) then
        dipsoft =fi_qq(m_HTop,zero,p,2,4,z,1)
        dipfini =fi_qq(m_HTop,zero,p,2,4,z,2)
        dipplus =fi_qq(m_HTop,zero,p,2,4,z,3)
        emi = 2
        mtrsq = Tree_23
      endif
      if(n.eq.11) then
        dipsoft =ff_qq(m_HTop,m_HTop,p,1,2,z,1)
        dipfini =ff_qq(m_HTop,m_HTop,p,1,2,z,2)
        dipplus =ff_qq(m_HTop,m_HTop,p,1,2,z,3)
        emi = 3
        mtrsq = Tree_34
      endif
      if(n.eq.12) then
        dipsoft =ff_qq(m_HTop,m_HTop,p,2,1,z,1)
        dipfini =ff_qq(m_HTop,m_HTop,p,2,1,z,2)
        dipplus =ff_qq(m_HTop,m_HTop,p,2,1,z,3)
        emi = 3
        mtrsq = Tree_34
      endif


      if(emi.eq.1) then
        res(1) = res(1) + (dipsoft-dipplus)*mtrsq
        res(2) = res(2) + (dipfini+dipplus)*mtrsq
      endif
      if(emi.eq.2) then
        res(1) = res(1) + (dipsoft-dipplus)*mtrsq
        res(3) = res(3) + (dipfini+dipplus)*mtrsq
      endif
      if(emi.eq.3) then
        res(1) = res(1) + (dipsoft-dipplus)*mtrsq
        res(2) = res(2) + (dipfini+dipplus)*0.5_dp*mtrsq
        res(3) = res(3) + (dipfini+dipplus)*0.5_dp*mtrsq
      endif
! res(1) = res(1) + dipsoft*mtrsq; print *, "eps check"


   enddo
   res(1:3) = -alpha_sOver2Pi * res(1:3)

! print *,"1", res(2:3)


       mtrsq = Tree_ij(0)
! !        epcorr=epinv+2d0*dlog(renscale/facscale)
       epcorr=epinv
       AP(1)= (11.0_dp/6.0_dp*CA -   5d0/3.0_dp)  ! 5d0=nf
       AP(2)= CA * two*(one/z+z*(one-z)-two)
       AP(3)= CA * two/(one-z)
       AP(1:3) = AP(1:3) * alpha_sOver2Pi * epcorr * mtrsq

       res(1) = res(1) + 2d0*(AP(1)-AP(3))
       res(2) = res(2) + (AP(2) + AP(3))
       res(3) = res(3) + (AP(2) + AP(3))
! res(1) = res(1)+2*ap(1); print *, "eps check"

! print *,"2", res(2:3)
! pause

! 
! print *, "LO ID", Tree_ij(0)
! print *, "ID check",res(1)/alpha_sOver2Pi/Tree_ij(0)
! print *, "ID check",res(2)/alpha_sOver2Pi/Tree_ij(0)
! print *, "ID check",res(3)/alpha_sOver2Pi/Tree_ij(0)
! pause



  RETURN 
  END SUBROUTINE



      FUNCTION Tree_GG_TTb_ij(MomTd,Mass2Td) Result(SqAmp)! this needs to be modified
      use modProcess
      use modParameters
      use modMyRecurrence
      use modKinematics
      use ModExoticDecay
      implicit none
      double precision,parameter :: Split_A=1d0
      double complex SqAmp(0:6),Amp(2,0:1),Res(0:3)
      integer col,colP,iHel,icorr,SecHel,A0barHel,A0Hel,HTHel,HTbarHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
!       double precision :: PSWgt2,PSWgt3,MomDK_ATop(0:3,1:3),MomDK_Top(0:3,1:3)
      type(Particle),target :: TopQuark(1:2)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      double precision :: ColCorr(1:2,1:2)
      real(8),parameter::ColCorr0(1:2,1:2)=(/64.D0/3.D0,-8.D0/3.D0,-8.D0/3.D0,64.D0/3.D0/)
      real(8),parameter::ColCorr1(1:2,1:2)=(/-72.D0,0.D0,0.D0,-72.D0/)
      real(8),parameter::ColCorr2(1:2,1:2)=(/-64.D0,8.D0,8.D0,8.D0/)
      real(8),parameter::ColCorr3(1:2,1:2)=(/8.D0,8.D0,8.D0,-64.D0/)
      real(8),parameter::ColCorr4(1:2,1:2)=(/8.D0,8.D0,8.D0,-64.D0/)
      real(8),parameter::ColCorr5(1:2,1:2)=(/-64.D0,8.D0,8.D0,8.D0/)
      real(8),parameter::ColCorr6(1:2,1:2)=(/-8.D0/9.D0,-80.D0/9.D0,-80.D0/9.D0,-8.D0/9.D0/)
      integer,target :: HelList(1:4,1:3)

      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
!       HelList(3,1:3)=(/0,1,0/)
!       HelList(4,1:3)=(/1,1,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)


              TopQuark(1)%PartType = Top_
              TopQuark(1)%Mom(1:4) = MomTd(0:3,3)
              TopQuark(2)%PartType = ATop_
              TopQuark(2)%Mom(1:4) = MomTd(0:3,4)


!      set particles
      call pol_mless(MomTd(0:3,1),+1,PolV(0:3,plus,1))
      call pol_mless(MomTd(0:3,1),-1,PolV(0:3,minus,1))

      call pol_mless(MomTd(0:3,2),+1,PolV(0:3,plus,2))
      call pol_mless(MomTd(0:3,2),-1,PolV(0:3,minus,2))
      Gluons(2)%Mom => MomTd(:,2)
      Gluons(2)%PartType => Glu_
      Gluons(2)%ExtRef => ExtRef
      Gluons(2)%Mass2 => Mass2Td(2)
      Gluons(2)%Mass  => MassTd(2)

      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => Top_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ATop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)

!      sum over helicities
      if( XTopDecays.eq.0 ) then
            SecHel=2
      elseif( XTopDecays.eq.1 ) then
            SecHel=1
      else
          call Error("Error in dipole Tree_GG_TTb")
      endif


      SqAmp(0:6) = (0d0,0d0)
do icorr=0,6
      if(icorr.eq.0) ColCorr(1:2,1:2)=ColCorr0(1:2,1:2)
      if(icorr.eq.1) ColCorr(1:2,1:2)=ColCorr1(1:2,1:2)
      if(icorr.eq.2) ColCorr(1:2,1:2)=ColCorr2(1:2,1:2)
      if(icorr.eq.3) ColCorr(1:2,1:2)=ColCorr3(1:2,1:2)
      if(icorr.eq.4) ColCorr(1:2,1:2)=ColCorr4(1:2,1:2)
      if(icorr.eq.5) ColCorr(1:2,1:2)=ColCorr5(1:2,1:2)
      if(icorr.eq.6) ColCorr(1:2,1:2)=ColCorr6(1:2,1:2)


      do A0barHel=-1,+1,SecHel
      do A0Hel=-1,+1,SecHel
        IF( XTOPDECAYS.EQ.1 ) THEN
         call HTopBHDecay(TopQuark(2),DKX_HTBH_LO,A0barHel,MomDK(1:4,5:9))
         call HTopBHDecay(TopQuark(1),DKX_HTBH_LO,A0Hel,MomDK(1:4,10:14))
        ELSEIF(XTOPDECAYS.EQ.0 ) THEN
         TopQuark(1)%Helicity = A0Hel
         TopQuark(2)%Helicity = A0barHel
         TopQuark(1)%Mass = MassTd(3)
         TopQuark(2)%Mass = MassTd(4)
         TopQuark(1)%Mass2= Mass2Td(3)
         TopQuark(2)%Mass2= Mass2Td(4)
         call HTopBHDecay(TopQuark(2),DKX_HTBH_LO,A0barHel,MomDK(1:4,5:9))
         call HTopBHDecay(TopQuark(1),DKX_HTBH_LO,A0Hel,MomDK(1:4,10:14))
        ENDIF

        Quarks(1)%Pol => TopQuark(1)%Pol
        Quarks(2)%Pol => TopQuark(2)%Pol
        Quarks(1)%Helicity => TopQuark(1)%Helicity
        Quarks(2)%Helicity => TopQuark(2)%Helicity

      do iHel=1,2

      Gluons(2)%Pol => PolV(:,HelList(iHel,1),2)
      Gluons(2)%Helicity => HelList(iHel,1)
!       Quarks(1)%Pol => PolV(:,HelList(iHel,2),3)
!       Quarks(1)%Helicity => HelList(iHel,2)
!       Quarks(2)%Pol => PolV(:,HelList(iHel,3),4)
!       Quarks(2)%Helicity => HelList(iHel,3)

!      calc currents
      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,1,0,0/) )
      Amp(1,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(1,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,0,0,1/) )
      Amp(2,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(2,minus) = Res(0:3).dot.PolV(0:3,minus,1)


      do col =1,2
      do colP=1,2

        SqAmp(icorr) = SqAmp(icorr) + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus)) *Amp(col,plus)     &
                                                                   + dconjg(Amp(colP,minus))*Amp(col,minus) )
      enddo
      enddo

      enddo

      enddo
      enddo
      SqAmp(icorr) = alpha_s4Pi**2 / SpinAvg / ColAvg *SqAmp(icorr)
enddo

!      momentum crossing backwards
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

      return
      END FUNCTION





END MODULE
