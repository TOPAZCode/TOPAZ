module modIntDipoles_QQBTTBG_noDK
implicit none
private

public:: EvalIntDipoles_QQBTTBG_noDK

integer, parameter,private  :: dp = selected_real_kind(15)
real(dp),parameter,private :: zero = 0.0_dp


contains






      SUBROUTINE EvalIntDipoles_QQBTTBG_noDK(p,z,res)
      use modParameters
      use modMisc
      use ModIntDipoles
      implicit none
      real(dp), intent(in) :: p(4,4)
      real(dp), intent(out) :: res(1:3)
      real(dp) :: dipsoft,dipfini,dipplus,mtrsq,AP(1:3),epcorr
      real(dp) :: Tree_12,Tree_13,Tree_14,Tree_23,Tree_24,Tree_34,z
      real(dp) :: CF,L,Q2
      integer  :: n,emi,in1,in2
      complex(dp) :: TreeMom(1:4,1:4)

       res(1:3) = zero
       CF=4d0/3d0


!  tree momenta for g g -> t tb
   TreeMom(1:4,1) =-dcmplx( p(1:4,3) )
   TreeMom(1:4,2) =-dcmplx( p(1:4,4) )
   TreeMom(1:4,3) = dcmplx( p(1:4,2) )
   TreeMom(1:4,4) = dcmplx( p(1:4,1) )

   Tree_12 = Tree_UUb_TTb_12(TreeMom,(/0d0,0d0,m_Top**2,m_Top**2/)) * 0.5d0
   Tree_13 = Tree_UUb_TTb_13(TreeMom,(/0d0,0d0,m_Top**2,m_Top**2/)) * 0.5d0
   Tree_14 = Tree_UUb_TTb_14(TreeMom,(/0d0,0d0,m_Top**2,m_Top**2/)) * 0.5d0
   Tree_23 = Tree_UUb_TTb_23(TreeMom,(/0d0,0d0,m_Top**2,m_Top**2/)) * 0.5d0
   Tree_24 = Tree_UUb_TTb_24(TreeMom,(/0d0,0d0,m_Top**2,m_Top**2/)) * 0.5d0
   Tree_34 = Tree_UUb_TTb_34(TreeMom,(/0d0,0d0,m_Top**2,m_Top**2/)) * 0.5d0


   do n=1,12
      if(n.eq.1) then
        dipsoft =if_qq(zero,m_Top,p,3,1,z,1)
        dipfini =if_qq(zero,m_Top,p,3,1,z,2)
        dipplus =if_qq(zero,m_Top,p,3,1,z,3)
        emi = 1 ! mom #3 is emitting
        mtrsq = Tree_14   ! (1-4)-color corr. for (1-3)-dipole because t and tbar are exchanged
      endif
      if(n.eq.2) then
        dipsoft =if_qq(zero,m_Top,p,3,2,z,1)
        dipfini =if_qq(zero,m_Top,p,3,2,z,2)
        dipplus =if_qq(zero,m_Top,p,3,2,z,3)
        emi = 1 ! mom #3 is emitting
        mtrsq = Tree_13
      endif
      if(n.eq.3) then
        dipsoft =if_qq(zero,m_Top,p,4,1,z,1)
        dipfini =if_qq(zero,m_Top,p,4,1,z,2)
        dipplus =if_qq(zero,m_Top,p,4,1,z,3)
        emi = 2 ! mom #4 is emitting
        mtrsq = Tree_24
      endif
      if(n.eq.4) then
        dipsoft =if_qq(zero,m_Top,p,4,2,z,1)
        dipfini =if_qq(zero,m_Top,p,4,2,z,2)
        dipplus =if_qq(zero,m_Top,p,4,2,z,3)
        emi = 2 ! mom #4 is emitting
        mtrsq = Tree_23
      endif
      if(n.eq.5) then
        dipsoft =ii_qq(zero,zero,p,3,4,z,1)
        dipfini =ii_qq(zero,zero,p,3,4,z,2)
        dipplus =ii_qq(zero,zero,p,3,4,z,3)
        emi = 1 ! mom #3 is emitting
        mtrsq = Tree_34
      endif
      if(n.eq.6) then
        dipsoft =ii_qq(zero,zero,p,4,3,z,1)
        dipfini =ii_qq(zero,zero,p,4,3,z,2)
        dipplus =ii_qq(zero,zero,p,4,3,z,3)
        emi = 2 ! mom #4 is emitting
        mtrsq = Tree_34
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
        mtrsq = Tree_12
      endif
      if(n.eq.12) then
        dipsoft =ff_qq(m_Top,m_Top,p,2,1,z,1)
        dipfini =ff_qq(m_Top,m_Top,p,2,1,z,2)
        dipplus =ff_qq(m_Top,m_Top,p,2,1,z,3)
        emi = 3
        mtrsq = Tree_12
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

! res(1) = res(1) +  dipsoft*mtrsq  ; print *, "for delta-fct. check
   enddo
   res(1:3) = -alpha_sOver2Pi * res(1:3)

       mtrsq = Tree_UUb_TTb_00(TreeMom,(/0d0,0d0,m_Top**2,m_Top**2/))
! !        epcorr=epinv+2d0*dlog(renscale/facscale)
       epcorr=epinv
       AP(1)= 3d0/2d0*CF
       AP(2)= (-1d0-z)*CF
       AP(3)= 2d0*CF/(1d0-z)
       AP(1:3) = AP(1:3) * alpha_sOver2Pi *epcorr * mtrsq

       res(1) = res(1) + 2d0*(AP(1)-AP(3))
       res(2) = res(2) + (AP(2) + AP(3))
       res(3) = res(3) + (AP(2) + AP(3))
! res(1) = res(1) + 2d0*AP(1)  ! for delta-fct. check


!    z=0.323d0
!    do n=1,3
!         print *, n
!         Q2 = (2d0*(p(1:4,3).dot.p(1:4,4))) ! agreement with MCFM
!         L  = dlog(Q2/(2d0*m_Top)**2)
!         print *, ii_qq(zero,zero,p,3,4,z,n)
!         print *, ii_mqq(z,L,0d0,n)

!         Q2 = (2d0*(-p(1:4,3).dot.p(1:4,2))) ! NO agreement with MCFM
!         L  = dlog(Q2/(2d0*m_Top)**2)
!         print *, if_gg(zero,m_Top,p,3,2,z,n)
!         print *, if_mgg(z,L,m_Top/dsqrt(Q2),n)
!
!         Q2 = (2d0*(-p(1:4,3).dot.p(1:4,2)))  ! agreement with MCFM
!         L  = dlog(Q2/(2d0*m_Top)**2)
!         print *, fi_qq(m_Top,zero,p,2,3,z,n)
!         print *, fi_mqq(z,L,m_Top/dsqrt(Q2),n)
!
!         Q2 = (2d0*(p(1:4,1).dot.p(1:4,2)))+2d0*m_Top**2 ! agreement with MCFM
!         L  = dlog(Q2/(2d0*m_Top)**2)
!         print *, ff_qq(m_Top,m_Top,p,1,2,z,n)
!         print *, ff_mqq(z,L,m_Top/dsqrt(Q2),n)
!    enddo
!    stop


  RETURN
  END SUBROUTINE




      FUNCTION Tree_UUb_TTb_12(MomTd,Mass2Td) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      implicit none
      double precision,parameter :: Split_A=1d0
      double complex SqAmp,Amp(1,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:1),Quarks(1:4)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(9)
      real(8),parameter::ColCorr(1:1,1:1)=(/8.D0/3.D0/)
      integer,target :: HelList(1:16,1:4)

      HelList(1,1:4)=(/0,0,0,0/)
      HelList(2,1:4)=(/1,0,0,0/)
      HelList(3,1:4)=(/0,1,0,0/)
      HelList(4,1:4)=(/1,1,0,0/)
      HelList(5,1:4)=(/0,0,1,0/)
      HelList(6,1:4)=(/1,0,1,0/)
      HelList(7,1:4)=(/0,1,1,0/)
      HelList(8,1:4)=(/1,1,1,0/)
      HelList(9,1:4)=(/0,0,0,1/)
      HelList(10,1:4)=(/1,0,0,1/)
      HelList(11,1:4)=(/0,1,0,1/)
      HelList(12,1:4)=(/1,1,0,1/)
      HelList(13,1:4)=(/0,0,1,1/)
      HelList(14,1:4)=(/1,0,1,1/)
      HelList(15,1:4)=(/0,1,1,1/)
      HelList(16,1:4)=(/1,1,1,1/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

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

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(3)%Mom => MomTd(:,3)
      Quarks(3)%PartType => Top_
      Quarks(3)%ExtRef => ExtRef
      Quarks(3)%Mass2 => Mass2Td(3)
      Quarks(3)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(4)%Mom => MomTd(:,4)
      Quarks(4)%PartType => ATop_
      Quarks(4)%ExtRef => ExtRef
      Quarks(4)%Mass2 => Mass2Td(4)
      Quarks(4)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,16
      Quarks(1)%Pol => PolV(:,HelList(iHel,1),1)
      Quarks(1)%Helicity => HelList(iHel,1)
      Quarks(2)%Pol => PolV(:,HelList(iHel,2),2)
      Quarks(2)%Helicity => HelList(iHel,2)
      Quarks(3)%Pol => PolV(:,HelList(iHel,3),3)
      Quarks(3)%Helicity => HelList(iHel,3)
      Quarks(4)%Pol => PolV(:,HelList(iHel,4),4)
      Quarks(4)%Helicity => HelList(iHel,4)

!      calc currents
      Res(0:3) = cur_f_4f((/Gluons(1)/),(/Quarks(2),Quarks(3),Quarks(4)/),Quarks(1)%PartType,(/0,0,0,0,0/),0)
      Amp(1,plus)  = psp1_(Res(0:3),Quarks(1)%Pol(1:4))

      do col =1,1
      do colP=1,1
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus)) *Amp(col,plus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

!      momentum crossing backwards
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

      return
      END FUNCTION






      FUNCTION Tree_UUb_TTb_13(MomTd,Mass2Td) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      implicit none
      double precision,parameter :: Split_A=1d0
      double complex SqAmp,Amp(1,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:1),Quarks(1:4)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(9)
      real(8),parameter::ColCorr(1:1,1:1)=(/-56.D0/3.D0/)
      integer,target :: HelList(1:16,1:4)

      HelList(1,1:4)=(/0,0,0,0/)
      HelList(2,1:4)=(/1,0,0,0/)
      HelList(3,1:4)=(/0,1,0,0/)
      HelList(4,1:4)=(/1,1,0,0/)
      HelList(5,1:4)=(/0,0,1,0/)
      HelList(6,1:4)=(/1,0,1,0/)
      HelList(7,1:4)=(/0,1,1,0/)
      HelList(8,1:4)=(/1,1,1,0/)
      HelList(9,1:4)=(/0,0,0,1/)
      HelList(10,1:4)=(/1,0,0,1/)
      HelList(11,1:4)=(/0,1,0,1/)
      HelList(12,1:4)=(/1,1,0,1/)
      HelList(13,1:4)=(/0,0,1,1/)
      HelList(14,1:4)=(/1,0,1,1/)
      HelList(15,1:4)=(/0,1,1,1/)
      HelList(16,1:4)=(/1,1,1,1/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

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

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(3)%Mom => MomTd(:,3)
      Quarks(3)%PartType => Top_
      Quarks(3)%ExtRef => ExtRef
      Quarks(3)%Mass2 => Mass2Td(3)
      Quarks(3)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(4)%Mom => MomTd(:,4)
      Quarks(4)%PartType => ATop_
      Quarks(4)%ExtRef => ExtRef
      Quarks(4)%Mass2 => Mass2Td(4)
      Quarks(4)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,16
      Quarks(1)%Pol => PolV(:,HelList(iHel,1),1)
      Quarks(1)%Helicity => HelList(iHel,1)
      Quarks(2)%Pol => PolV(:,HelList(iHel,2),2)
      Quarks(2)%Helicity => HelList(iHel,2)
      Quarks(3)%Pol => PolV(:,HelList(iHel,3),3)
      Quarks(3)%Helicity => HelList(iHel,3)
      Quarks(4)%Pol => PolV(:,HelList(iHel,4),4)
      Quarks(4)%Helicity => HelList(iHel,4)

!      calc currents
      Res(0:3) = cur_f_4f((/Gluons(1)/),(/Quarks(2),Quarks(3),Quarks(4)/),Quarks(1)%PartType,(/0,0,0,0,0/),0)
      Amp(1,plus)  = psp1_(Res(0:3),Quarks(1)%Pol(1:4))

      do col =1,1
      do colP=1,1
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus)) *Amp(col,plus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

!      momentum crossing backwards
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

      return
      END FUNCTION






      FUNCTION Tree_UUb_TTb_14(MomTd,Mass2Td) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      implicit none
      double precision,parameter :: Split_A=1d0
      double complex SqAmp,Amp(1,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:1),Quarks(1:4)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(9)
      real(8),parameter::ColCorr(1:1,1:1)=(/-16.D0/3.D0/)
      integer,target :: HelList(1:16,1:4)

      HelList(1,1:4)=(/0,0,0,0/)
      HelList(2,1:4)=(/1,0,0,0/)
      HelList(3,1:4)=(/0,1,0,0/)
      HelList(4,1:4)=(/1,1,0,0/)
      HelList(5,1:4)=(/0,0,1,0/)
      HelList(6,1:4)=(/1,0,1,0/)
      HelList(7,1:4)=(/0,1,1,0/)
      HelList(8,1:4)=(/1,1,1,0/)
      HelList(9,1:4)=(/0,0,0,1/)
      HelList(10,1:4)=(/1,0,0,1/)
      HelList(11,1:4)=(/0,1,0,1/)
      HelList(12,1:4)=(/1,1,0,1/)
      HelList(13,1:4)=(/0,0,1,1/)
      HelList(14,1:4)=(/1,0,1,1/)
      HelList(15,1:4)=(/0,1,1,1/)
      HelList(16,1:4)=(/1,1,1,1/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

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

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(3)%Mom => MomTd(:,3)
      Quarks(3)%PartType => Top_
      Quarks(3)%ExtRef => ExtRef
      Quarks(3)%Mass2 => Mass2Td(3)
      Quarks(3)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(4)%Mom => MomTd(:,4)
      Quarks(4)%PartType => ATop_
      Quarks(4)%ExtRef => ExtRef
      Quarks(4)%Mass2 => Mass2Td(4)
      Quarks(4)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,16
      Quarks(1)%Pol => PolV(:,HelList(iHel,1),1)
      Quarks(1)%Helicity => HelList(iHel,1)
      Quarks(2)%Pol => PolV(:,HelList(iHel,2),2)
      Quarks(2)%Helicity => HelList(iHel,2)
      Quarks(3)%Pol => PolV(:,HelList(iHel,3),3)
      Quarks(3)%Helicity => HelList(iHel,3)
      Quarks(4)%Pol => PolV(:,HelList(iHel,4),4)
      Quarks(4)%Helicity => HelList(iHel,4)

!      calc currents
      Res(0:3) = cur_f_4f((/Gluons(1)/),(/Quarks(2),Quarks(3),Quarks(4)/),Quarks(1)%PartType,(/0,0,0,0,0/),0)
      Amp(1,plus)  = psp1_(Res(0:3),Quarks(1)%Pol(1:4))

      do col =1,1
      do colP=1,1
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus)) *Amp(col,plus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

!      momentum crossing backwards
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

      return
      END FUNCTION






      FUNCTION Tree_UUb_TTb_23(MomTd,Mass2Td) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      implicit none
      double precision,parameter :: Split_A=1d0
      double complex SqAmp,Amp(1,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:1),Quarks(1:4)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(9)
      real(8),parameter::ColCorr(1:1,1:1)=(/-16.D0/3.D0/)
      integer,target :: HelList(1:16,1:4)

      HelList(1,1:4)=(/0,0,0,0/)
      HelList(2,1:4)=(/1,0,0,0/)
      HelList(3,1:4)=(/0,1,0,0/)
      HelList(4,1:4)=(/1,1,0,0/)
      HelList(5,1:4)=(/0,0,1,0/)
      HelList(6,1:4)=(/1,0,1,0/)
      HelList(7,1:4)=(/0,1,1,0/)
      HelList(8,1:4)=(/1,1,1,0/)
      HelList(9,1:4)=(/0,0,0,1/)
      HelList(10,1:4)=(/1,0,0,1/)
      HelList(11,1:4)=(/0,1,0,1/)
      HelList(12,1:4)=(/1,1,0,1/)
      HelList(13,1:4)=(/0,0,1,1/)
      HelList(14,1:4)=(/1,0,1,1/)
      HelList(15,1:4)=(/0,1,1,1/)
      HelList(16,1:4)=(/1,1,1,1/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

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

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(3)%Mom => MomTd(:,3)
      Quarks(3)%PartType => Top_
      Quarks(3)%ExtRef => ExtRef
      Quarks(3)%Mass2 => Mass2Td(3)
      Quarks(3)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(4)%Mom => MomTd(:,4)
      Quarks(4)%PartType => ATop_
      Quarks(4)%ExtRef => ExtRef
      Quarks(4)%Mass2 => Mass2Td(4)
      Quarks(4)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,16
      Quarks(1)%Pol => PolV(:,HelList(iHel,1),1)
      Quarks(1)%Helicity => HelList(iHel,1)
      Quarks(2)%Pol => PolV(:,HelList(iHel,2),2)
      Quarks(2)%Helicity => HelList(iHel,2)
      Quarks(3)%Pol => PolV(:,HelList(iHel,3),3)
      Quarks(3)%Helicity => HelList(iHel,3)
      Quarks(4)%Pol => PolV(:,HelList(iHel,4),4)
      Quarks(4)%Helicity => HelList(iHel,4)

!      calc currents
      Res(0:3) = cur_f_4f((/Gluons(1)/),(/Quarks(2),Quarks(3),Quarks(4)/),Quarks(1)%PartType,(/0,0,0,0,0/),0)
      Amp(1,plus)  = psp1_(Res(0:3),Quarks(1)%Pol(1:4))

      do col =1,1
      do colP=1,1
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus)) *Amp(col,plus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

!      momentum crossing backwards
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

      return
      END FUNCTION






      FUNCTION Tree_UUb_TTb_24(MomTd,Mass2Td) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      implicit none
      double precision,parameter :: Split_A=1d0
      double complex SqAmp,Amp(1,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:1),Quarks(1:4)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(9)
      real(8),parameter::ColCorr(1:1,1:1)=(/-56.D0/3.D0/)
      integer,target :: HelList(1:16,1:4)

      HelList(1,1:4)=(/0,0,0,0/)
      HelList(2,1:4)=(/1,0,0,0/)
      HelList(3,1:4)=(/0,1,0,0/)
      HelList(4,1:4)=(/1,1,0,0/)
      HelList(5,1:4)=(/0,0,1,0/)
      HelList(6,1:4)=(/1,0,1,0/)
      HelList(7,1:4)=(/0,1,1,0/)
      HelList(8,1:4)=(/1,1,1,0/)
      HelList(9,1:4)=(/0,0,0,1/)
      HelList(10,1:4)=(/1,0,0,1/)
      HelList(11,1:4)=(/0,1,0,1/)
      HelList(12,1:4)=(/1,1,0,1/)
      HelList(13,1:4)=(/0,0,1,1/)
      HelList(14,1:4)=(/1,0,1,1/)
      HelList(15,1:4)=(/0,1,1,1/)
      HelList(16,1:4)=(/1,1,1,1/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

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

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(3)%Mom => MomTd(:,3)
      Quarks(3)%PartType => Top_
      Quarks(3)%ExtRef => ExtRef
      Quarks(3)%Mass2 => Mass2Td(3)
      Quarks(3)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(4)%Mom => MomTd(:,4)
      Quarks(4)%PartType => ATop_
      Quarks(4)%ExtRef => ExtRef
      Quarks(4)%Mass2 => Mass2Td(4)
      Quarks(4)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,16
      Quarks(1)%Pol => PolV(:,HelList(iHel,1),1)
      Quarks(1)%Helicity => HelList(iHel,1)
      Quarks(2)%Pol => PolV(:,HelList(iHel,2),2)
      Quarks(2)%Helicity => HelList(iHel,2)
      Quarks(3)%Pol => PolV(:,HelList(iHel,3),3)
      Quarks(3)%Helicity => HelList(iHel,3)
      Quarks(4)%Pol => PolV(:,HelList(iHel,4),4)
      Quarks(4)%Helicity => HelList(iHel,4)

!      calc currents
      Res(0:3) = cur_f_4f((/Gluons(1)/),(/Quarks(2),Quarks(3),Quarks(4)/),Quarks(1)%PartType,(/0,0,0,0,0/),0)
      Amp(1,plus)  = psp1_(Res(0:3),Quarks(1)%Pol(1:4))

      do col =1,1
      do colP=1,1
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus)) *Amp(col,plus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

!      momentum crossing backwards
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

      return
      END FUNCTION






      FUNCTION Tree_UUb_TTb_34(MomTd,Mass2Td) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      implicit none
      double precision,parameter :: Split_A=1d0
      double complex SqAmp,Amp(1,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:1),Quarks(1:4)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(9)
      real(8),parameter::ColCorr(1:1,1:1)=(/8.D0/3.D0/)
      integer,target :: HelList(1:16,1:4)

      HelList(1,1:4)=(/0,0,0,0/)
      HelList(2,1:4)=(/1,0,0,0/)
      HelList(3,1:4)=(/0,1,0,0/)
      HelList(4,1:4)=(/1,1,0,0/)
      HelList(5,1:4)=(/0,0,1,0/)
      HelList(6,1:4)=(/1,0,1,0/)
      HelList(7,1:4)=(/0,1,1,0/)
      HelList(8,1:4)=(/1,1,1,0/)
      HelList(9,1:4)=(/0,0,0,1/)
      HelList(10,1:4)=(/1,0,0,1/)
      HelList(11,1:4)=(/0,1,0,1/)
      HelList(12,1:4)=(/1,1,0,1/)
      HelList(13,1:4)=(/0,0,1,1/)
      HelList(14,1:4)=(/1,0,1,1/)
      HelList(15,1:4)=(/0,1,1,1/)
      HelList(16,1:4)=(/1,1,1,1/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

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

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(3)%Mom => MomTd(:,3)
      Quarks(3)%PartType => Top_
      Quarks(3)%ExtRef => ExtRef
      Quarks(3)%Mass2 => Mass2Td(3)
      Quarks(3)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(4)%Mom => MomTd(:,4)
      Quarks(4)%PartType => ATop_
      Quarks(4)%ExtRef => ExtRef
      Quarks(4)%Mass2 => Mass2Td(4)
      Quarks(4)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,16
      Quarks(1)%Pol => PolV(:,HelList(iHel,1),1)
      Quarks(1)%Helicity => HelList(iHel,1)
      Quarks(2)%Pol => PolV(:,HelList(iHel,2),2)
      Quarks(2)%Helicity => HelList(iHel,2)
      Quarks(3)%Pol => PolV(:,HelList(iHel,3),3)
      Quarks(3)%Helicity => HelList(iHel,3)
      Quarks(4)%Pol => PolV(:,HelList(iHel,4),4)
      Quarks(4)%Helicity => HelList(iHel,4)

!      calc currents
      Res(0:3) = cur_f_4f((/Gluons(1)/),(/Quarks(2),Quarks(3),Quarks(4)/),Quarks(1)%PartType,(/0,0,0,0,0/),0)
      Amp(1,plus)  = psp1_(Res(0:3),Quarks(1)%Pol(1:4))

      do col =1,1
      do colP=1,1
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus)) *Amp(col,plus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

!      momentum crossing backwards
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

      return
      END FUNCTION






      FUNCTION Tree_UUb_TTb_00(MomTd,Mass2Td) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      implicit none
      double precision,parameter :: Split_A=1d0
      double complex SqAmp,Amp(1,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:1),Quarks(1:4)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(9)
      real(8),parameter::ColCorr(1:1,1:1)=(/8.D0/)
      integer,target :: HelList(1:16,1:4)

      HelList(1,1:4)=(/0,0,0,0/)
      HelList(2,1:4)=(/1,0,0,0/)
      HelList(3,1:4)=(/0,1,0,0/)
      HelList(4,1:4)=(/1,1,0,0/)
      HelList(5,1:4)=(/0,0,1,0/)
      HelList(6,1:4)=(/1,0,1,0/)
      HelList(7,1:4)=(/0,1,1,0/)
      HelList(8,1:4)=(/1,1,1,0/)
      HelList(9,1:4)=(/0,0,0,1/)
      HelList(10,1:4)=(/1,0,0,1/)
      HelList(11,1:4)=(/0,1,0,1/)
      HelList(12,1:4)=(/1,1,0,1/)
      HelList(13,1:4)=(/0,0,1,1/)
      HelList(14,1:4)=(/1,0,1,1/)
      HelList(15,1:4)=(/0,1,1,1/)
      HelList(16,1:4)=(/1,1,1,1/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

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

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(3)%Mom => MomTd(:,3)
      Quarks(3)%PartType => Top_
      Quarks(3)%ExtRef => ExtRef
      Quarks(3)%Mass2 => Mass2Td(3)
      Quarks(3)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(4)%Mom => MomTd(:,4)
      Quarks(4)%PartType => ATop_
      Quarks(4)%ExtRef => ExtRef
      Quarks(4)%Mass2 => Mass2Td(4)
      Quarks(4)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,16
      Quarks(1)%Pol => PolV(:,HelList(iHel,1),1)
      Quarks(1)%Helicity => HelList(iHel,1)
      Quarks(2)%Pol => PolV(:,HelList(iHel,2),2)
      Quarks(2)%Helicity => HelList(iHel,2)
      Quarks(3)%Pol => PolV(:,HelList(iHel,3),3)
      Quarks(3)%Helicity => HelList(iHel,3)
      Quarks(4)%Pol => PolV(:,HelList(iHel,4),4)
      Quarks(4)%Helicity => HelList(iHel,4)

!      calc currents
      Res(0:3) = cur_f_4f((/Gluons(1)/),(/Quarks(2),Quarks(3),Quarks(4)/),Quarks(1)%PartType,(/0,0,0,0,0/),0)
      Amp(1,plus)  = psp1_(Res(0:3),Quarks(1)%Pol(1:4))

      do col =1,1
      do colP=1,1
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus)) *Amp(col,plus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

!      momentum crossing backwards
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

      return
      END FUNCTION







END MODULE
