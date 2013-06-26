module ModIntDipoles_GGTTBGZ
use ModParameters
use ModProcess
use ModMisc
implicit none
private

public:: EvalIntDipoles_GGTTBGZ
integer, parameter  :: dp = selected_real_kind(15)


!       double precision, private, parameter :: NCol=3d0
!       double precision, private, parameter :: TR=0.5d0
!       double precision, private, parameter :: CA=2d0*TR*NCol
!       double precision, private, parameter :: CF=TR*(NCol**2-1d0)/NCol

contains


!-------------------------------------------------------------
!    ordering of the momenta in the p-array:
!     outgoing convention:
!     bar t = 1, Z =2, t=3 ,  incoming glue=4,
!     5=incoming glue
!-----------------------------------------------------------

      subroutine EvalIntDipoles_GGTTBGZ(p,MomDK,z,res)
      use ModIntDipoles
      use ModAmplitudes
      implicit none
      real(dp), intent(in) :: p(4,5),MomDK(1:4,1:8),z
      real(dp), intent(out) :: res(1:3)
      real(dp) :: dipsoft,dipfini,dipplus,AP(1:3),CA
      type(Particle),save :: ExtParticles(1:5)
      type(TreeProcess),save :: TreeAmpsDip(1:2)
      complex(dp) :: Bm(-1:1,1:2,-1:1,-1:1,-1:1,-1:1),mtrsq
      integer :: iTree,i1,i2,i3,i4,i5,i,j,i6,n
      integer :: hel(1:5),emi , Nmax(5),N2Jump
      character(2) :: diptype
      real(dp) :: C(1:2,1:2),epcorr
      logical, save :: first_time = .true.

      res(1:3) = zero
      CA=3d0

       if( first_time ) then
             call InitTrees(2,2,2,TreeAmpsDip,NumBoson=1)
             call InitProcess_TbTGGZ(ExtParticles(1:5))
             TreeAmpsDip(1)%PartRef(1:5) = (/1,5,2,3,4/)
             TreeAmpsDip(2)%PartRef(1:5) = (/1,5,2,4,3/)
             do iTree=1,2
                call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles(1:5))
             enddo
             first_time=.false.
      endif

      BM = (zero,zero) ! first, the ``amplitude matrix'' is set to zero

!--- after momentum mapping -- sum over colors and polarizations

        Nmax = 1
        if (TopDecays.ge.1) then ! Top decays
              Nmax(3) = -1
              Nmax(1) = -1
        endif
        N2Jump = 1
        if (ZDecays.ge.1) then ! Z decays
            N2jump=2
        endif

       do i1=-1,Nmax(1),2
         do i2 = -1,Nmax(2),N2Jump! Z boson
           do i3 = -1,Nmax(3),2
             do i4 = -1,Nmax(4),2
               do i5 = -1,Nmax(5),2


           hel(1) = i1
           hel(2) = i2
           hel(3) = i3
           hel(4) = i4
           hel(5) = i5


    call SetPolarization((/p(1:4,1),p(1:4,3),p(1:4,4),p(1:4,5),p(1:4,2)/),momDK(1:4,1:8),(/hel(1),hel(3),hel(4),hel(5),hel(2)/),ExtParticles(1:5))


      do i6 = 1,2
         call EvalTree2(TreeAmpsDip(i6),Bm(i1,i6,i2,i3,i4,i5))
      enddo


      enddo
      enddo
      enddo
      enddo
      enddo




   do n=1,12

     mtrsq = (zero,zero)
     call cc_gg_ttgamg_soft(n,C)

       do i1=-1,Nmax(1),2
         do i2 = -1,Nmax(2),N2Jump! Z boson
           do i3 = -1,Nmax(3),2
             do i4 = -1,Nmax(4),2
               do i5 = -1,Nmax(5),2

     do i=1,2
     do j=1,2
          mtrsq = mtrsq + C(i,j)*Bm(i1,i,i2,i3,i4,i5)*dconjg(Bm(i1,j,i2,i3,i4,i5))
     enddo
     enddo

      enddo
      enddo
      enddo
      enddo
      enddo


! spin & color avg., color convention adjustment
      mtrsq = mtrsq * alpha_s4Pi**2 * alpha4Pi /4d0/64d0


      if(n.eq.1) then
      dipsoft =ff_qq(m_TOP,m_TOP,p,1,3,z,1)
      dipfini =ff_qq(m_TOP,m_TOP,p,1,3,z,2)
      dipplus =ff_qq(m_TOP,m_TOP,p,1,3,z,3)
      emi = 3 ! final state dipole
      endif
      if(n.eq.2) then
      dipsoft =fi_qq(m_TOP,zero,p,1,4,z,1)
      dipfini =fi_qq(m_TOP,zero,p,1,4,z,2)
      dipplus =fi_qq(m_TOP,zero,p,1,4,z,3)
!      emi = 3 ! final state dipole
      emi = 1
      endif
      if(n.eq.3) then
      dipsoft =fi_qq(m_TOP,zero,p,1,5,z,1)
      dipfini =fi_qq(m_TOP,zero,p,1,5,z,2)
      dipplus =fi_qq(m_TOP,zero,p,1,5,z,3)
!      emi = 3 ! final state dipole
      emi = 2
      endif
      if(n.eq.4) then
      dipsoft =ff_qq(m_TOP,m_TOP,p,3,1,z,1)
      dipfini =ff_qq(m_TOP,m_TOP,p,3,1,z,2)
      dipplus =ff_qq(m_TOP,m_TOP,p,3,1,z,3)
      emi = 3 ! final state dipole
      endif
      if(n.eq.5) then
      dipsoft =fi_qq(m_TOP,zero,p,3,4,z,1)
      dipfini =fi_qq(m_TOP,zero,p,3,4,z,2)
      dipplus =fi_qq(m_TOP,zero,p,3,4,z,3)
!      emi = 3 ! final state dipole
      emi = 1
      endif
      if(n.eq.6) then
      dipsoft =fi_qq(m_TOP,zero,p,3,5,z,1)
      dipfini =fi_qq(m_TOP,zero,p,3,5,z,2)
      dipplus =fi_qq(m_TOP,zero,p,3,5,z,3)
!      emi = 3 ! final state dipole
      emi = 2
      endif
      if(n.eq.7) then
      dipsoft =if_gg(zero,m_TOP,p,4,1,z,1)
      dipfini =if_gg(zero,m_TOP,p,4,1,z,2)
      dipplus =if_gg(zero,m_TOP,p,4,1,z,3)
     emi = 1 ! pdf1 is emitting
      endif
      if(n.eq.8) then
      dipsoft =if_gg(zero,m_TOP,p,4,3,z,1)
      dipfini =if_gg(zero,m_TOP,p,4,3,z,2)
      dipplus =if_gg(zero,m_TOP,p,4,3,z,3)
     emi = 1 ! pdf1 is emitting
      endif
      if(n.eq.9) then
      dipsoft =ii_gg(zero,zero,p,4,5,z,1)
      dipfini =ii_gg(zero,zero,p,4,5,z,2)
      dipplus =ii_gg(zero,zero,p,4,5,z,3)
     emi = 1 ! pdf1 is emitting
      endif
      if(n.eq.10) then
      dipsoft =if_gg(zero,m_TOP,p,5,1,z,1)
      dipfini =if_gg(zero,m_TOP,p,5,1,z,2)
      dipplus =if_gg(zero,m_TOP,p,5,1,z,3)
     emi = 2 ! pdf2 is emitting
      endif
      if(n.eq.11) then
      dipsoft =if_gg(zero,m_TOP,p,5,3,z,1)
      dipfini =if_gg(zero,m_TOP,p,5,3,z,2)
      dipplus =if_gg(zero,m_TOP,p,5,3,z,3)
     emi = 2 ! pdf2 is emitting
      endif
      if(n.eq.12) then
      dipsoft =ii_gg(zero,zero,p,5,4,z,1)
      dipfini =ii_gg(zero,zero,p,5,4,z,2)
      dipplus =ii_gg(zero,zero,p,5,4,z,3)
     emi = 2 ! pdf2 is emitting
      endif


! print *, "comment here"
!       RES(1) = RES(1) +  DIPSOFT*MTRSQ  ! FOR DELTA-FCT. CHECK

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


   enddo
   res(1:3)=alpha_sOver2Pi*0.5d0*res(1:3)  ! overall normalization: alpha_s/4pi


!    eval tree
     mtrsq = (zero,zero)
     call  cc_gg_ttgam(C)
       do i1=-1,Nmax(1),2
         do i2 = -1,Nmax(2),N2Jump! Z boson
           do i3 = -1,Nmax(3),2
             do i4 = -1,Nmax(4),2
               do i5 = -1,Nmax(5),2

     do i=1,2
     do j=1,2
          mtrsq = mtrsq + C(i,j)*Bm(i1,i,i2,i3,i4,i5)*dconjg(Bm(i1,j,i2,i3,i4,i5))
     enddo
     enddo

      enddo
      enddo
      enddo
      enddo
      enddo
     mtrsq = mtrsq * alpha_s4Pi**2 * alpha4Pi /4d0/64d0
 

! !        epcorr=epinv+2d0*dlog(renscale/facscale)
     epcorr=epinv
     AP(1)= (11.0_dp/6.0_dp*CA - nf_light/3.0_dp )
     AP(2)= CA * two*(one/z+z*(one-z)-two)
     AP(3)= CA * two/(one-z)
     AP(1:3) = AP(1:3) * alpha_sOver2Pi * epcorr * mtrsq
! print *, "comment here"
!        RES(1) = RES(1) + 2D0*AP(1)  ! FOR DELTA-FCT. CHECK
     res(1) = res(1) + 2d0*(AP(1)-AP(3))
     res(2) = res(2) + (AP(2) + AP(3))
     res(3) = res(3) + (AP(2) + AP(3))


! print *, "tree idip",mtrsq
! print *, "RES(1)/mtrsq",RES(1)/(mtrsq*alpha_sOver2Pi)
! pause

  return
  end subroutine





      subroutine cc_gg_ttgamg_soft(n,C) ! soft cc matrix, gg-> tt+gamma+g
      integer, intent(in) :: n
      real(dp), intent(out)  :: C(2,2)

      if(n.eq.1) then
      C(1,1)=8.0_dp/9.0_dp
      C(1,2)=80.0_dp/9.0_dp
      C(2,1)=80.0_dp/9.0_dp
      C(2,2)=8.0_dp/9.0_dp
      endif
      if(n.eq.2) then
      C(1,1)=-8.0_dp
      C(1,2)=-8.0_dp
      C(2,1)=-8.0_dp
      C(2,2)=64.0_dp
      endif
      if(n.eq.3) then
      C(1,1)=64.0_dp
      C(1,2)=-8.0_dp
      C(2,1)=-8.0_dp
      C(2,2)=-8.0_dp
      endif
      if(n.eq.4) then
      C(1,1)=8.0_dp/9.0_dp
      C(1,2)=80.0_dp/9.0_dp
      C(2,1)=80.0_dp/9.0_dp
      C(2,2)=8.0_dp/9.0_dp
      endif
      if(n.eq.5) then
      C(1,1)=64.0_dp
      C(1,2)=-8.0_dp
      C(2,1)=-8.0_dp
      C(2,2)=-8.0_dp
      endif
      if(n.eq.6) then
      C(1,1)=-8.0_dp
      C(1,2)=-8.0_dp
      C(2,1)=-8.0_dp
      C(2,2)=64.0_dp
      endif
      if(n.eq.7) then
      C(1,1)=-8.0_dp
      C(1,2)=-8.0_dp
      C(2,1)=-8.0_dp
      C(2,2)=64.0_dp
      endif
      if(n.eq.8) then
      C(1,1)=64.0_dp
      C(1,2)=-8.0_dp
      C(2,1)=-8.0_dp
      C(2,2)=-8.0_dp
      endif
      if(n.eq.9) then
      C(1,1)=72.0_dp
      C(1,2)=0.0_dp
      C(2,1)=0.0_dp
      C(2,2)=72.0_dp
      endif
      if(n.eq.10) then
      C(1,1)=64.0_dp
      C(1,2)=-8.0_dp
      C(2,1)=-8.0_dp
      C(2,2)=-8.0_dp
      endif
      if(n.eq.11) then
      C(1,1)=-8.0_dp
      C(1,2)=-8.0_dp
      C(2,1)=-8.0_dp
      C(2,2)=64.0_dp
      endif
      if(n.eq.12) then
      C(1,1)=72.0_dp
      C(1,2)=0.0_dp
      C(2,1)=0.0_dp
      C(2,2)=72.0_dp
      endif


      end subroutine cc_gg_ttgamg_soft


      subroutine cc_gg_ttgam(C)
      real(dp), intent(out) :: C(2,2)


       C(1,1) = 64.0_dp/3.0_dp
       C(1,2) =  - 8.0_dp/3
       C(2,1) =  - 8.0_dp/3.0_dp
       C(2,2) = 64.0_dp/3.0_dp

      end subroutine cc_gg_ttgam



SUBROUTINE InitProcess_TbTGGZ(ExtParticles)
use ModProcess
implicit none
type(Particle) :: ExtParticles(:)

  ExtParticles(1)%PartType = ATop_
  ExtParticles(1)%ExtRef   = 1
  ExtParticles(1)%Mass = m_Top
  ExtParticles(1)%Mass2= ExtParticles(1)%Mass**2

  ExtParticles(2)%PartType = Top_
  ExtParticles(2)%ExtRef   = 2
  ExtParticles(2)%Mass = m_Top
  ExtParticles(2)%Mass2= ExtParticles(2)%Mass**2

  ExtParticles(3)%PartType = Glu_
  ExtParticles(3)%ExtRef   = 3
  ExtParticles(3)%Mass = 0d0
  ExtParticles(3)%Mass2= 0d0

  ExtParticles(4)%PartType = Glu_
  ExtParticles(4)%ExtRef   = 4
  ExtParticles(4)%Mass = 0d0
  ExtParticles(4)%Mass2= 0d0

  ExtParticles(5)%PartType = Z0_
  ExtParticles(5)%ExtRef   = 5
  ExtParticles(5)%Mass = m_Z
  ExtParticles(5)%Mass2= ExtParticles(5)%Mass**2


RETURN
END SUBROUTINE InitProcess_TbTGGZ



SUBROUTINE SetPolarization(Mom,MomDK,Hel,ExtParticles)
use ModMisc
use ModProcess
use ModTopDecay
use ModZDecay
implicit none
type(Particle) :: ExtParticles(1:5)
real(8) :: Mom(1:4,1:5),MomDK(1:4,1:8)
integer :: Hel(1:5)

     ExtParticles(1)%Mom(1:4) = dcmplx(Mom(1:4,1))   ! HERE WAS A BUG: this was inside the (TopDecays.ge.1) condition
     ExtParticles(2)%Mom(1:4) = dcmplx(Mom(1:4,2))
     if (TopDecays.ge.1) then
        call TopDecay(ExtParticles(1),DK_LO,MomDK(1:4,1:3))
        call TopDecay(ExtParticles(2),DK_LO,MomDK(1:4,4:6))
     else
        call vSpi(ExtParticles(1)%Mom(1:4),ExtParticles(1)%Mass,Hel(1),ExtParticles(1)%Pol(1:4))
        call ubarSpi(ExtParticles(2)%Mom(1:4),ExtParticles(2)%Mass,Hel(2),ExtParticles(2)%Pol(1:4))
    endif


     ExtParticles(5)%Mom(1:4) = dcmplx(Mom(1:4,5))
     ExtParticles(5)%Helicity = Hel(5)
     if (ZDecays.ge.1) then
          call ZDecay(ExtParticles(5),DK_LO,MomDK(1:4,7:8))
     else
          call pol_massSR(ExtParticles(5)%Mom(1:4),ExtParticles(5)%Mass,ExtParticles(5)%Helicity,ExtParticles(5)%Pol(1:4))
    endif

    ExtParticles(3)%Mom(1:4) = dcmplx(Mom(1:4,3))
    call pol_mless(ExtParticles(3)%Mom(1:4),Hel(3),ExtParticles(3)%Pol(1:4))

    ExtParticles(4)%Mom(1:4) = dcmplx(Mom(1:4,4))
    call pol_mless(ExtParticles(4)%Mom(1:4),Hel(4),ExtParticles(4)%Pol(1:4))


! ExtParticles(3)%Pol(1:4)=ExtParticles(3)%Mom(1:4); print *, "check gauge inv. in dipoles"

RETURN
END SUBROUTINE




end  module
