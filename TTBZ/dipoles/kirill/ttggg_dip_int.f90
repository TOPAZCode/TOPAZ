module ttgggg_dipint
implicit none
private

public:: ttgggg_dip_int

integer, parameter,private  :: dp = selected_real_kind(15)
real(dp),parameter,private :: zero = 0.0_dp


contains






      SUBROUTINE ttgggg_dip_int(p,z,fx1,fx2,fx1z,fx2z,res)
      use modParameters
      use modMisc
      use mod_dipoles2
      implicit none
      real(dp), intent(in) :: p(4,4),fx1,fx2,fx1z,fx2z
      real(dp), intent(out) :: res
      real(dp) :: dipsoft,dipfini,dipplus,mtrsq,AP(1:3),epcorr
      real(dp) :: Tree_12,Tree_13,Tree_14,Tree_23,Tree_24,Tree_34,z
      real(dp) :: CA,L,Q2
      integer  :: n,emi,in1,in2
      complex(dp) :: TreeMom(1:4,1:4)

       res = zero
       CA=3d0


!  tree momenta for g g -> t tb
   TreeMom(1:4,1) =-dcmplx( p(1:4,3) )
   TreeMom(1:4,2) =-dcmplx( p(1:4,4) )
   TreeMom(1:4,3) = dcmplx( p(1:4,2) )
   TreeMom(1:4,4) = dcmplx( p(1:4,1) )

   Tree_12  = Tree_GG_TTb_12(TreeMom,(/0d0,0d0,m_Top**2,m_Top**2/)) * 0.5d0
   Tree_13  = Tree_GG_TTb_13(TreeMom,(/0d0,0d0,m_Top**2,m_Top**2/)) * 0.5d0
   Tree_14  = Tree_GG_TTb_14(TreeMom,(/0d0,0d0,m_Top**2,m_Top**2/)) * 0.5d0
   Tree_23  = Tree_GG_TTb_23(TreeMom,(/0d0,0d0,m_Top**2,m_Top**2/)) * 0.5d0
   Tree_24  = Tree_GG_TTb_24(TreeMom,(/0d0,0d0,m_Top**2,m_Top**2/)) * 0.5d0
   Tree_34  = Tree_GG_TTb_34(TreeMom,(/0d0,0d0,m_Top**2,m_Top**2/)) * 0.5d0

   do n=1,12
      if(n.eq.1) then
!       {if}=initial emitter, final spectator, {gg}=gluon parent, gluon emitted, {zero,m_Top}=parent,spectator mass, {3,1}=parent, spectator in p
        dipsoft =if_gg(zero,m_Top,p,3,1,z,1)
        dipfini =if_gg(zero,m_Top,p,3,1,z,2)
        dipplus =if_gg(zero,m_Top,p,3,1,z,3)
        emi = 1 ! mom #3 is emitting
        mtrsq = Tree_14   ! (1-4)-color corr. for (1-3)-dipole because t and tbar are exchanged
      endif
      if(n.eq.2) then
        dipsoft =if_gg(zero,m_Top,p,3,2,z,1)
        dipfini =if_gg(zero,m_Top,p,3,2,z,2)
        dipplus =if_gg(zero,m_Top,p,3,2,z,3)
        emi = 1 ! mom #3 is emitting
        mtrsq = Tree_13
      endif
      if(n.eq.3) then
        dipsoft =if_gg(zero,m_Top,p,4,1,z,1)
        dipfini =if_gg(zero,m_Top,p,4,1,z,2)
        dipplus =if_gg(zero,m_Top,p,4,1,z,3)
        emi = 2 ! mom #4 is emitting
        mtrsq = Tree_24
      endif
      if(n.eq.4) then
        dipsoft =if_gg(zero,m_Top,p,4,2,z,1)
        dipfini =if_gg(zero,m_Top,p,4,2,z,2)
        dipplus =if_gg(zero,m_Top,p,4,2,z,3)
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
        dipsoft =fi_qq(m_Top,zero,p,1,3,z,1)
        dipfini =fi_qq(m_Top,zero,p,1,3,z,2)
        dipplus =fi_qq(m_Top,zero,p,1,3,z,3)
        emi = 3
        mtrsq = Tree_14
      endif
      if(n.eq.8) then
        dipsoft =fi_qq(m_Top,zero,p,1,4,z,1)
        dipfini =fi_qq(m_Top,zero,p,1,4,z,2)
        dipplus =fi_qq(m_Top,zero,p,1,4,z,3)
        emi = 3
        mtrsq = Tree_24
      endif
      if(n.eq.9) then
        dipsoft =fi_qq(m_Top,zero,p,2,3,z,1)
        dipfini =fi_qq(m_Top,zero,p,2,3,z,2)
        dipplus =fi_qq(m_Top,zero,p,2,3,z,3)
        emi = 3
        mtrsq = Tree_13
      endif
      if(n.eq.10) then
        dipsoft =fi_qq(m_Top,zero,p,2,4,z,1)
        dipfini =fi_qq(m_Top,zero,p,2,4,z,2)
        dipplus =fi_qq(m_Top,zero,p,2,4,z,3)
        emi = 3
        mtrsq = Tree_23
      endif
      if(n.eq.11) then
        dipsoft =ff_qq(m_Top,m_Top,p,1,2,z,1)
        dipfini =ff_qq(m_Top,m_Top,p,1,2,z,2)
        dipplus =ff_qq(m_Top,m_Top,p,1,2,z,3)
        emi = 3
        mtrsq = Tree_34
      endif
      if(n.eq.12) then
        dipsoft =ff_qq(m_Top,m_Top,p,2,1,z,1)
        dipfini =ff_qq(m_Top,m_Top,p,2,1,z,2)
        dipplus =ff_qq(m_Top,m_Top,p,2,1,z,3)
        emi = 3
        mtrsq = Tree_34
      endif


      if(emi.eq.1) then
      res = res + (dipsoft-dipplus)*mtrsq*fx1 *fx2    &
                + (dipfini+dipplus)*mtrsq*fx1z*fx2/z
      endif
      if(emi.eq.2) then
      res = res + (dipsoft-dipplus)*mtrsq*fx1*fx2     &
                + (dipfini+dipplus)*mtrsq*fx1*fx2z/z
      endif
      if(emi.eq.3) then
      res = res + (dipsoft-dipplus)*mtrsq*fx1*fx2             &
                + (dipfini+dipplus)*0.5_dp*mtrsq*fx1 *fx2z/z  &
                + (dipfini+dipplus)*0.5_dp*mtrsq*fx1z*fx2/z
      endif

!       res= res + (dipfini+dipplus)*mtrsq
!       res= res + dipsoft*mtrsq
   enddo
   res = -alpha_sOver2Pi * res

       mtrsq = Tree_GG_TTb_00(TreeMom,(/0d0,0d0,m_Top**2,m_Top**2/))
! !        epcorr=epinv+2d0*dlog(renscale/facscale)
       epcorr=epinv
       AP(1)= two*(11.0_dp/6.0_dp*CA - 5d0/3.0_dp) *epcorr  ! 5d0=nf
       AP(2)= CA * four*(one/z+z*(one-z)-two)*epcorr
       AP(3)= CA * four/(one-z)*epcorr
       AP(1:3) = AP(1:3) * alpha_sOver2Pi
       res = res   + mtrsq*(AP(1) - AP(3))*fx1*fx2    &
                   + mtrsq*(AP(2) + AP(3))*fx1z*fx2/z &
                   + mtrsq*(AP(1) - AP(3))*fx1*fx2    &
                   + mtrsq*(AP(2) + AP(3))*fx1*fx2z/z


!   checks
!        epcorr=epinv  ! this assumes mu_f=mu_r
!        AP(1)= two*(11.0_dp/6.0_dp*CA - 5d0/3.0_dp) *epcorr  ! 5d0=nf    ! this contribution is checked vs. I-operator
!        AP(2)= CA * four*(one/z+z*(one-z)-two)*epcorr
!        AP(3)= CA * four/(one-z)*epcorr
!        AP(1:3) = AP(1:3) * alpha_sOver2Pi
!        res = res + mtrsq*(AP(2)+AP(3))
!        res = res + mtrsq*AP(1)

   z=0.323d0
   do n=1,3
        print *, n
!         Q2 = (2d0*(p(1:4,3).dot.p(1:4,4))) ! agreement with MCFM
!         L  = dlog(Q2/(2d0*m_Top)**2)
!         print *, ii_gg(zero,zero,p,3,4,z,n)
!         print *, ii_mgg(z,L,0d0,n)

!         Q2 = (2d0*(-p(1:4,3).dot.p(1:4,2))) ! agreement with MCFM
!         L  = dlog(Q2/(2d0*m_Top)**2)
!         print *, if_gg(zero,m_Top,p,3,2,z,n)
!         print *, if_mgg(z,L,m_Top/dsqrt(Q2),n)

        Q2 = (2d0*(-p(1:4,3).dot.p(1:4,2)))  ! no direct agreement with MCFM because of different formulae
        L  = dlog(Q2/(2d0*m_Top)**2)
        print *, fi_qq(m_Top,zero,p,2,3,z,n)
        print *, fi_mqq(z,L,m_Top/dsqrt(Q2),n)

!         Q2 = (2d0*(p(1:4,1).dot.p(1:4,2)))+2d0*m_Top**2 ! agreement with MCFM
!         L  = dlog(Q2/(2d0*m_Top)**2)
!         print *, ff_qq(m_Top,m_Top,p,1,2,z,n)
!         print *, ff_mqq(z,L,m_Top/dsqrt(Q2),n)
   enddo
   print *, fi_qq(m_Top,zero,p,2,3,z,2)+fi_qq(m_Top,zero,p,2,3,z,3), fi_mqq(z,L,m_Top/dsqrt(Q2),2)+fi_mqq(z,L,m_Top/dsqrt(Q2),3)
   stop


  RETURN
  END SUBROUTINE ttgggg_dip_int












      FUNCTION Tree_GG_TTb_12(MomTd,Mass2Td) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      implicit none
      double precision,parameter :: Split_A=1d0
      double complex SqAmp,Amp(2,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      real(8),parameter::ColCorr(1:2,1:2)=(/-72.D0,0.D0,0.D0,-72.D0/)
      integer,target :: HelList(1:4,1:3)

      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
      HelList(3,1:3)=(/0,1,0/)
      HelList(4,1:3)=(/1,1,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

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

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => Top_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ATop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,4
      Gluons(2)%Pol => PolV(:,HelList(iHel,1),2)
      Gluons(2)%Helicity => HelList(iHel,1)
      Quarks(1)%Pol => PolV(:,HelList(iHel,2),3)
      Quarks(1)%Helicity => HelList(iHel,2)
      Quarks(2)%Pol => PolV(:,HelList(iHel,3),4)
      Quarks(2)%Helicity => HelList(iHel,3)

!      calc currents
      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,1,0,0/) )
      Amp(1,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(1,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,0,0,1/) )
      Amp(2,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(2,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      do col =1,2
      do colP=1,2
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus)) *Amp(col,plus)     &
                                                    + dconjg(Amp(colP,minus))*Amp(col,minus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * 2d0*SqAmp

!      momentum crossing backwards
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

      return
      END FUNCTION






      FUNCTION Tree_GG_TTb_13(MomTd,Mass2Td) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      implicit none
      double precision,parameter :: Split_A=1d0
      double complex SqAmp,Amp(2,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      real(8),parameter::ColCorr(1:2,1:2)=(/-64.D0,8.D0,8.D0,8.D0/)
      integer,target :: HelList(1:4,1:3)

      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
      HelList(3,1:3)=(/0,1,0/)
      HelList(4,1:3)=(/1,1,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

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

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => Top_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ATop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,4
      Gluons(2)%Pol => PolV(:,HelList(iHel,1),2)
      Gluons(2)%Helicity => HelList(iHel,1)
      Quarks(1)%Pol => PolV(:,HelList(iHel,2),3)
      Quarks(1)%Helicity => HelList(iHel,2)
      Quarks(2)%Pol => PolV(:,HelList(iHel,3),4)
      Quarks(2)%Helicity => HelList(iHel,3)

!      calc currents
      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,1,0,0/) )
      Amp(1,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(1,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,0,0,1/) )
      Amp(2,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(2,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      do col =1,2
      do colP=1,2
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus)) *Amp(col,plus)     &
                                                    + dconjg(Amp(colP,minus))*Amp(col,minus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * 2d0*SqAmp

!      momentum crossing backwards
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

      return
      END FUNCTION






      FUNCTION Tree_GG_TTb_14(MomTd,Mass2Td) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      implicit none
      double precision,parameter :: Split_A=1d0
      double complex SqAmp,Amp(2,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      real(8),parameter::ColCorr(1:2,1:2)=(/8.D0,8.D0,8.D0,-64.D0/)
      integer,target :: HelList(1:4,1:3)

      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
      HelList(3,1:3)=(/0,1,0/)
      HelList(4,1:3)=(/1,1,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

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

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => Top_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ATop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,4
      Gluons(2)%Pol => PolV(:,HelList(iHel,1),2)
      Gluons(2)%Helicity => HelList(iHel,1)
      Quarks(1)%Pol => PolV(:,HelList(iHel,2),3)
      Quarks(1)%Helicity => HelList(iHel,2)
      Quarks(2)%Pol => PolV(:,HelList(iHel,3),4)
      Quarks(2)%Helicity => HelList(iHel,3)

!      calc currents
      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,1,0,0/) )
      Amp(1,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(1,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,0,0,1/) )
      Amp(2,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(2,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      do col =1,2
      do colP=1,2
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus)) *Amp(col,plus)     &
                                                    + dconjg(Amp(colP,minus))*Amp(col,minus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * 2d0*SqAmp

!      momentum crossing backwards
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

      return
      END FUNCTION






      FUNCTION Tree_GG_TTb_23(MomTd,Mass2Td) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      implicit none
      double precision,parameter :: Split_A=1d0
      double complex SqAmp,Amp(2,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      real(8),parameter::ColCorr(1:2,1:2)=(/8.D0,8.D0,8.D0,-64.D0/)
      integer,target :: HelList(1:4,1:3)

      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
      HelList(3,1:3)=(/0,1,0/)
      HelList(4,1:3)=(/1,1,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

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

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => Top_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ATop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,4
      Gluons(2)%Pol => PolV(:,HelList(iHel,1),2)
      Gluons(2)%Helicity => HelList(iHel,1)
      Quarks(1)%Pol => PolV(:,HelList(iHel,2),3)
      Quarks(1)%Helicity => HelList(iHel,2)
      Quarks(2)%Pol => PolV(:,HelList(iHel,3),4)
      Quarks(2)%Helicity => HelList(iHel,3)

!      calc currents
      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,1,0,0/) )
      Amp(1,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(1,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,0,0,1/) )
      Amp(2,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(2,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      do col =1,2
      do colP=1,2
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus)) *Amp(col,plus)     &
                                                    + dconjg(Amp(colP,minus))*Amp(col,minus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * 2d0*SqAmp

!      momentum crossing backwards
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

      return
      END FUNCTION






      FUNCTION Tree_GG_TTb_24(MomTd,Mass2Td) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      implicit none
      double precision,parameter :: Split_A=1d0
      double complex SqAmp,Amp(2,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      real(8),parameter::ColCorr(1:2,1:2)=(/-64.D0,8.D0,8.D0,8.D0/)
      integer,target :: HelList(1:4,1:3)

      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
      HelList(3,1:3)=(/0,1,0/)
      HelList(4,1:3)=(/1,1,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

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

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => Top_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ATop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,4
      Gluons(2)%Pol => PolV(:,HelList(iHel,1),2)
      Gluons(2)%Helicity => HelList(iHel,1)
      Quarks(1)%Pol => PolV(:,HelList(iHel,2),3)
      Quarks(1)%Helicity => HelList(iHel,2)
      Quarks(2)%Pol => PolV(:,HelList(iHel,3),4)
      Quarks(2)%Helicity => HelList(iHel,3)

!      calc currents
      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,1,0,0/) )
      Amp(1,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(1,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,0,0,1/) )
      Amp(2,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(2,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      do col =1,2
      do colP=1,2
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus)) *Amp(col,plus)     &
                                                    + dconjg(Amp(colP,minus))*Amp(col,minus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * 2d0*SqAmp

!      momentum crossing backwards
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

      return
      END FUNCTION






      FUNCTION Tree_GG_TTb_34(MomTd,Mass2Td) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      implicit none
      double precision,parameter :: Split_A=1d0
      double complex SqAmp,Amp(2,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      real(8),parameter::ColCorr(1:2,1:2)=(/-8.D0/9.D0,-80.D0/9.D0,-80.D0/9.D0,-8.D0/9.D0/)
      integer,target :: HelList(1:4,1:3)

      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
      HelList(3,1:3)=(/0,1,0/)
      HelList(4,1:3)=(/1,1,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

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

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => Top_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ATop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,4
      Gluons(2)%Pol => PolV(:,HelList(iHel,1),2)
      Gluons(2)%Helicity => HelList(iHel,1)
      Quarks(1)%Pol => PolV(:,HelList(iHel,2),3)
      Quarks(1)%Helicity => HelList(iHel,2)
      Quarks(2)%Pol => PolV(:,HelList(iHel,3),4)
      Quarks(2)%Helicity => HelList(iHel,3)

!      calc currents
      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,1,0,0/) )
      Amp(1,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(1,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,0,0,1/) )
      Amp(2,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(2,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      do col =1,2
      do colP=1,2
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus)) *Amp(col,plus)     &
                                                    + dconjg(Amp(colP,minus))*Amp(col,minus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * 2d0*SqAmp

!      momentum crossing backwards
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

      return
      END FUNCTION






      FUNCTION Tree_GG_TTb_00(MomTd,Mass2Td) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      implicit none
      double precision,parameter :: Split_A=1d0
      double complex SqAmp,Amp(2,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      real(8),parameter::ColCorr(1:2,1:2)=(/64.D0/3.D0,-8.D0/3.D0,-8.D0/3.D0,64.D0/3.D0/)
      integer,target :: HelList(1:4,1:3)

      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
      HelList(3,1:3)=(/0,1,0/)
      HelList(4,1:3)=(/1,1,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

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

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => Top_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ATop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,4
      Gluons(2)%Pol => PolV(:,HelList(iHel,1),2)
      Gluons(2)%Helicity => HelList(iHel,1)
      Quarks(1)%Pol => PolV(:,HelList(iHel,2),3)
      Quarks(1)%Helicity => HelList(iHel,2)
      Quarks(2)%Pol => PolV(:,HelList(iHel,3),4)
      Quarks(2)%Helicity => HelList(iHel,3)

!      calc currents
      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,1,0,0/) )
      Amp(1,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(1,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,0,0,1/) )
      Amp(2,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(2,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      do col =1,2
      do colP=1,2
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus)) *Amp(col,plus)     &
                                                    + dconjg(Amp(colP,minus))*Amp(col,minus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * 2d0*SqAmp

!      momentum crossing backwards
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

      return
      END FUNCTION




END MODULE ttgggg_dipint
