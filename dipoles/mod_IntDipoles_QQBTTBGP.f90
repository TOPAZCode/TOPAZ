module ModIntDipoles_QQBTTBGP
use ModParameters
use ModProcess
use ModMisc
implicit none
private

public:: EvalIntDipoles_QQBTTBGP
integer, parameter  :: dp = selected_real_kind(15)


!       double precision, private, parameter :: NCol=3d0
!       double precision, private, parameter :: TR=0.5d0
!       double precision, private, parameter :: CA=2d0*TR*NCol
!       double precision, private, parameter :: CF=TR*(NCol**2-1d0)/NCol

contains


!-------------------------------------------------------------
!    ordering of the momenta in the p-array:
!     outgoing convention:
!     bar t = 1, photon =2, t=3 ,  incoming quark=4,
!     5=incoming antiquark
!-----------------------------------------------------------

      subroutine EvalIntDipoles_QQBTTBGP(p,MomDK,z,res)
      use ModIntDipoles
      use ModAmplitudes
      implicit none
      real(dp), intent(in) :: p(4,5),MomDK(1:4,1:6),z
      real(dp), intent(out) :: res(1:2,1:3)
      real(dp) :: dipsoft,dipfini,dipplus,AP(1:2,1:3),CF
      type(Particle),save :: ExtParticles(1:5)
      type(TreeProcess),save :: TreeAmpsDip(1:2)
      complex(dp) :: Bm(1:2),mtrsq(2)
      complex(dp) :: AM(1:2,1:2)
      integer :: iTree,i1,i2,i3,i4,i5,i,j,i6,n
      integer :: hel(1:5),emi , Nmax(5)
      character(2) :: diptype
      real(dp) :: Cup(1:2,1:2),Cdn(1:2,1:2),epcorr
      logical, save :: first_time = .true.
      real,parameter :: PhotonCouplCorr=2d0

      res(1:2,1:3) = zero
      CF=4d0/3d0

       if( first_time ) then
              call InitTrees(4,1,2,TreeAmpsDip)
              call InitProcess_TbTQbQG(ExtParticles(1:5))
             TreeAmpsDip(1)%PartRef(1:5) = (/1,5,2,3,4/)
             TreeAmpsDip(2)%PartRef(1:5) = (/1,2,3,5,4/)
             do iTree=1,2
             call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles(1:5))
             enddo
             first_time=.false.
      endif

      AM = (zero,zero) ! first, the ``amplitude matrix'' is set to zero

!--- after momentum mapping -- sum over colors and polarizations

         Nmax = 1

        if (TopDecays.ge.1) then
            Nmax(1) = -1
            Nmax(3) = -1
        endif

       do i1=-1,Nmax(1),2
         do i2 = -1,Nmax(2),2
           do i3 = -1,Nmax(3),2
             do i4 = -1,Nmax(4),2
               do i5 = -1,Nmax(5),2

          hel(1) = i1
          hel(2) = i2
          hel(3) = i3
          hel(4) = i4
          hel(5) = i5

  call GenerateEventttqqg2((/p(1:4,1),p(1:4,3),p(1:4,4),p(1:4,5),p(1:4,2)/), &
  momDK(1:4,1:6),(/hel(1),hel(3),hel(4),hel(5),hel(2)/),ExtParticles(1:5))

          do i6 = 1,2
              call EvalTree2(TreeAmpsDip(i6),Bm(i6))
          enddo

          do i=1,2
          do j=1,2
              AM(i,j) = AM(i,j) + Bm(i)*dconjg(Bm(j))
          enddo
          enddo

      enddo
      enddo
      enddo
      enddo
      enddo


   do n=1,12

     mtrsq = (zero,zero)

        call cc_qq_ttgamg_soft(n,'up',Cup)
        call cc_qq_ttgamg_soft(n,'dn',Cdn)



     do i=1,2
     do j=1,2
          mtrsq(1) = mtrsq(1) + Cup(i,j)*real(AM(i,j),dp)
          mtrsq(2) = mtrsq(2) + Cdn(i,j)*real(AM(i,j),dp)
     enddo
     enddo

! spin & color avg., color convention adjustment

     mtrsq = mtrsq * alpha_s4Pi**2 * alpha4Pi * Q_top**2 /4d0/9d0 * PhotonCouplCorr


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
      dipsoft =if_qq(zero,m_TOP,p,4,1,z,1)
      dipfini =if_qq(zero,m_TOP,p,4,1,z,2)
      dipplus =if_qq(zero,m_TOP,p,4,1,z,3)
      emi = 1 ! pdf1 is emitting
      endif
      if(n.eq.8) then
      dipsoft =if_qq(zero,m_TOP,p,4,3,z,1)
      dipfini =if_qq(zero,m_TOP,p,4,3,z,2)
      dipplus =if_qq(zero,m_TOP,p,4,3,z,3)
     emi = 1 ! pdf1 is emitting
      endif
      if(n.eq.9) then
      dipsoft =ii_qq(zero,zero,p,4,5,z,1)
      dipfini =ii_qq(zero,zero,p,4,5,z,2)
      dipplus =ii_qq(zero,zero,p,4,5,z,3)
     emi = 1 ! pdf1 is emitting
      endif
      if(n.eq.10) then
      dipsoft =if_qq(zero,m_TOP,p,5,1,z,1)
      dipfini =if_qq(zero,m_TOP,p,5,1,z,2)
      dipplus =if_qq(zero,m_TOP,p,5,1,z,3)
     emi = 2 ! pdf2 is emitting
      endif
      if(n.eq.11) then
      dipsoft =if_qq(zero,m_TOP,p,5,3,z,1)
      dipfini =if_qq(zero,m_TOP,p,5,3,z,2)
      dipplus =if_qq(zero,m_TOP,p,5,3,z,3)
     emi = 2 ! pdf2 is emitting
      endif
      if(n.eq.12) then
      dipsoft =ii_qq(zero,zero,p,5,4,z,1)
      dipfini =ii_qq(zero,zero,p,5,4,z,2)
      dipplus =ii_qq(zero,zero,p,5,4,z,3)
     emi = 2 ! pdf2 is emitting
      endif


! print *, "comment here";  RES(:,1) = RES(:,1) +  DIPSOFT*MTRSQ(:)  ! FOR DELTA-FCT. CHECK
! print *, "dipsoft=0";   dipsoft=0d0

      if(emi.eq.1) then
        res(:,1) = res(:,1) + (dipsoft-dipplus)*mtrsq(:)
        res(:,2) = res(:,2) + (dipfini+dipplus)*mtrsq(:)
      endif
      if(emi.eq.2) then
        res(:,1) = res(:,1) + (dipsoft-dipplus)*mtrsq(:)
        res(:,3) = res(:,3) + (dipfini+dipplus)*mtrsq(:)
      endif
      if(emi.eq.3) then
        res(:,1) = res(:,1) + (dipsoft-dipplus)*mtrsq(:)
        res(:,2) = res(:,2) + (dipfini+dipplus)*0.5_dp*mtrsq(:)
        res(:,3) = res(:,3) + (dipfini+dipplus)*0.5_dp*mtrsq(:)
      endif

   enddo

! overall normalization: alpha_s/4pi
   res(1:2,1:3)=alpha_sOver2Pi*0.5d0*res(1:2,1:3)

! print *, "1",res(1:2,1:3)

!    eval tree
     mtrsq = (zero,zero)
     call  cc_qq_ttgam('up',Cup)
     call  cc_qq_ttgam('dn',Cdn)
     do i=1,2
     do j=1,2
          mtrsq(1) = mtrsq(1) + Cup(i,j)*real(AM(i,j),dp)
          mtrsq(2) = mtrsq(2) + Cdn(i,j)*real(AM(i,j),dp)
     enddo
     enddo
     mtrsq = mtrsq * alpha_s4Pi**2 * alpha4Pi * Q_top**2 /4d0/9d0 * PhotonCouplCorr


! !        epcorr=epinv+2d0*dlog(renscale/facscale)
       epcorr=epinv
       AP(:,1)= 3d0/2d0*CF
       AP(:,2)= (-1d0-z)*CF
       AP(:,3)= 2d0*CF/(1d0-z)
       AP(1,1:3) = AP(1,1:3) * alpha_sOver2Pi *epcorr * mtrsq(1)
       AP(2,1:3) = AP(2,1:3) * alpha_sOver2Pi *epcorr * mtrsq(2)


! print *, "comment here";  RES(:,1) = RES(:,1) + 2D0*AP(:,1)  ! FOR DELTA-FCT. CHECK
! print *, "AP(1)=0";   AP(:,1)=0d0
       res(1:2,1) = res(1:2,1) + 2d0*(AP(1:2,1)-AP(1:2,3)) ! here was a bug 2d0*AP(1:2,1)-AP(1:2,3)
       res(1:2,2) = res(1:2,2) + AP(1:2,2) + AP(1:2,3)
       res(1:2,3) = res(1:2,3) + AP(1:2,2) + AP(1:2,3)



! print *, "2",res(1:2,1:3)

! print *, "RES(1,1)/mtrsq",RES(1,1)/(mtrsq(1)*alpha_sOver2Pi)
! print *, "RES(2,1)/mtrsq",RES(2,1)/(mtrsq(2)*alpha_sOver2Pi)
! pause
! RES(1:2,1)=mtrsq(1:2)

  return
  end subroutine




      subroutine cc_qq_ttgamg_soft(n,ixq,C)  ! ixq labels the power of the
      integer, intent(in) :: n               ! electric charge of the
      real(dp), intent(out) :: C(2,2)        ! emitting quark
      character, intent(in) :: ixq*2

      C = 0.0_dp

      if (ixq.eq.'up') then



      if(n.eq.1) then
      C(1,1)=-8.0_dp/3.0_dp
      C(1,2)=-8.0_dp/3.0_dp*Q_up/Q_top
      C(2,1)=-8.0_dp/3.0_dp*Q_up/Q_top
      C(2,2)=-8.0_dp/3.0_dp*Q_up**2/Q_top**2
      endif
      if(n.eq.2) then
      C(1,1)=16.0_dp/3.0_dp
      C(1,2)=16.0_dp/3.0_dp*Q_up/Q_top
      C(2,1)=16.0_dp/3.0_dp*Q_up/Q_top
      C(2,2)=16.0_dp/3.0_dp*Q_up**2/Q_top**2
      endif
      if(n.eq.3) then
      C(1,1)=56.0_dp/3.0_dp
      C(1,2)=56.0_dp/3.0_dp*Q_up/Q_top
      C(2,1)=56.0_dp/3.0_dp*Q_up/Q_top
      C(2,2)=56.0_dp/3.0_dp*Q_up**2/Q_top**2
      endif
      if(n.eq.4) then
      C(1,1)=-8.0_dp/3.0_dp
      C(1,2)=-8.0_dp/3.0_dp*Q_up/Q_top
      C(2,1)=-8.0_dp/3.0_dp*Q_up/Q_top
      C(2,2)=-8.0_dp/3.0_dp*Q_up**2/Q_top**2
      endif
      if(n.eq.5) then
      C(1,1)=56.0_dp/3.0_dp
      C(1,2)=56.0_dp/3.0_dp*Q_up/Q_top
      C(2,1)=56.0_dp/3.0_dp*Q_up/Q_top
      C(2,2)=56.0_dp/3.0_dp*Q_up**2/Q_top**2
      endif
      if(n.eq.6) then
      C(1,1)=16.0_dp/3.0_dp
      C(1,2)=16.0_dp/3.0_dp*Q_up/Q_top
      C(2,1)=16.0_dp/3.0_dp*Q_up/Q_top
      C(2,2)=16.0_dp/3.0_dp*Q_up**2/Q_top**2
      endif
      if(n.eq.7) then
      C(1,1)=16.0_dp/3.0_dp
      C(1,2)=16.0_dp/3.0_dp*Q_up/Q_top
      C(2,1)=16.0_dp/3.0_dp*Q_up/Q_top
      C(2,2)=16.0_dp/3.0_dp*Q_up**2/Q_top**2
      endif
      if(n.eq.8) then
      C(1,1)=56.0_dp/3.0_dp
      C(1,2)=56.0_dp/3.0_dp*Q_up/Q_top
      C(2,1)=56.0_dp/3.0_dp*Q_up/Q_top
      C(2,2)=56.0_dp/3.0_dp*Q_up**2/Q_top**2
      endif
      if(n.eq.9) then
      C(1,1)=-8.0_dp/3.0_dp
      C(1,2)=-8.0_dp/3.0_dp*Q_up/Q_top
      C(2,1)=-8.0_dp/3.0_dp*Q_up/Q_top
      C(2,2)=-8.0_dp/3.0_dp*Q_up**2/Q_top**2
      endif
      if(n.eq.10) then
      C(1,1)=56.0_dp/3.0_dp
      C(1,2)=56.0_dp/3.0_dp*Q_up/Q_top
      C(2,1)=56.0_dp/3.0_dp*Q_up/Q_top
      C(2,2)=56.0_dp/3.0_dp*Q_up**2/Q_top**2
      endif
      if(n.eq.11) then
      C(1,1)=16.0_dp/3.0_dp
      C(1,2)=16.0_dp/3.0_dp*Q_up/Q_top
      C(2,1)=16.0_dp/3.0_dp*Q_up/Q_top
      C(2,2)=16.0_dp/3.0_dp*Q_up**2/Q_top**2
      endif
      if(n.eq.12) then
      C(1,1)=-8.0_dp/3.0_dp
      C(1,2)=-8.0_dp/3.0_dp*Q_up/Q_top
      C(2,1)=-8.0_dp/3.0_dp*Q_up/Q_top
      C(2,2)=-8.0_dp/3.0_dp*Q_up**2/Q_top**2
      endif



      elseif(ixq.eq.'dn') then

      if(n.eq.1) then
      C(1,1)=-8.0_dp/3.0_dp
      C(1,2)=-8.0_dp/3.0_dp*Q_dn/Q_top
      C(2,1)=-8.0_dp/3.0_dp*Q_dn/Q_top
      C(2,2)=-8.0_dp/3.0_dp*Q_dn**2/Q_top**2
      endif
      if(n.eq.2) then
      C(1,1)=16.0_dp/3.0_dp
      C(1,2)=16.0_dp/3.0_dp*Q_dn/Q_top
      C(2,1)=16.0_dp/3.0_dp*Q_dn/Q_top
      C(2,2)=16.0_dp/3.0_dp*Q_dn**2/Q_top**2
      endif
      if(n.eq.3) then
      C(1,1)=56.0_dp/3.0_dp
      C(1,2)=56.0_dp/3.0_dp*Q_dn/Q_top
      C(2,1)=56.0_dp/3.0_dp*Q_dn/Q_top
      C(2,2)=56.0_dp/3.0_dp*Q_dn**2/Q_top**2
      endif
      if(n.eq.4) then
      C(1,1)=-8.0_dp/3.0_dp
      C(1,2)=-8.0_dp/3.0_dp*Q_dn/Q_top
      C(2,1)=-8.0_dp/3.0_dp*Q_dn/Q_top
      C(2,2)=-8.0_dp/3.0_dp*Q_dn**2/Q_top**2
      endif
      if(n.eq.5) then
      C(1,1)=56.0_dp/3.0_dp
      C(1,2)=56.0_dp/3.0_dp*Q_dn/Q_top
      C(2,1)=56.0_dp/3.0_dp*Q_dn/Q_top
      C(2,2)=56.0_dp/3.0_dp*Q_dn**2/Q_top**2
      endif
      if(n.eq.6) then
      C(1,1)=16.0_dp/3.0_dp
      C(1,2)=16.0_dp/3.0_dp*Q_dn/Q_top
      C(2,1)=16.0_dp/3.0_dp*Q_dn/Q_top
      C(2,2)=16.0_dp/3.0_dp*Q_dn**2/Q_top**2
      endif
      if(n.eq.7) then
      C(1,1)=16.0_dp/3.0_dp
      C(1,2)=16.0_dp/3.0_dp*Q_dn/Q_top
      C(2,1)=16.0_dp/3.0_dp*Q_dn/Q_top
      C(2,2)=16.0_dp/3.0_dp*Q_dn**2/Q_top**2
      endif
      if(n.eq.8) then
      C(1,1)=56.0_dp/3.0_dp
      C(1,2)=56.0_dp/3.0_dp*Q_dn/Q_top
      C(2,1)=56.0_dp/3.0_dp*Q_dn/Q_top
      C(2,2)=56.0_dp/3.0_dp*Q_dn**2/Q_top**2
      endif
      if(n.eq.9) then
      C(1,1)=-8.0_dp/3.0_dp
      C(1,2)=-8.0_dp/3.0_dp*Q_dn/Q_top
      C(2,1)=-8.0_dp/3.0_dp*Q_dn/Q_top
      C(2,2)=-8.0_dp/3.0_dp*Q_dn**2/Q_top**2
      endif
      if(n.eq.10) then
      C(1,1)=56.0_dp/3.0_dp
      C(1,2)=56.0_dp/3.0_dp*Q_dn/Q_top
      C(2,1)=56.0_dp/3.0_dp*Q_dn/Q_top
      C(2,2)=56.0_dp/3.0_dp*Q_dn**2/Q_top**2
      endif
      if(n.eq.11) then
      C(1,1)=16.0_dp/3.0_dp
      C(1,2)=16.0_dp/3.0_dp*Q_dn/Q_top
      C(2,1)=16.0_dp/3.0_dp*Q_dn/Q_top
      C(2,2)=16.0_dp/3.0_dp*Q_dn**2/Q_top**2
      endif
      if(n.eq.12) then
      C(1,1)=-8.0_dp/3.0_dp
      C(1,2)=-8.0_dp/3.0_dp*Q_dn/Q_top
      C(2,1)=-8.0_dp/3.0_dp*Q_dn/Q_top
      C(2,2)=-8.0_dp/3.0_dp*Q_dn**2/Q_top**2
      endif


!       if(n.eq.1) then
!       C(1,1)=-8.0_dp/3.0_dp
!       C(1,2)=4.0_dp/3.0_dp
!       C(2,1)=4.0_dp/3.0_dp
!       C(2,2)=-2.0_dp/3.0_dp
!       endif
!       if(n.eq.2) then
!       C(1,1)=16.0_dp/3.0_dp
!       C(1,2)=-8.0_dp/3.0_dp
!       C(2,1)=-8.0_dp/3.0_dp
!       C(2,2)=4.0_dp/3.0_dp
!       endif
!       if(n.eq.3) then
!       C(1,1)=56.0_dp/3.0_dp
!       C(1,2)=-28.0_dp/3.0_dp
!       C(2,1)=-28.0_dp/3.0_dp
!       C(2,2)=14.0_dp/3.0_dp
!       endif
!       if(n.eq.4) then
!       C(1,1)=-8.0_dp/3.0_dp
!       C(1,2)=4.0_dp/3.0_dp
!       C(2,1)=4.0_dp/3.0_dp
!       C(2,2)=-2.0_dp/3.0_dp
!       endif
!       if(n.eq.5) then
!       C(1,1)=56.0_dp/3.0_dp
!       C(1,2)=-28.0_dp/3.0_dp
!       C(2,1)=-28.0_dp/3.0_dp
!       C(2,2)=14.0_dp/3.0_dp
!       endif
!       if(n.eq.6) then
!       C(1,1)=16.0_dp/3.0_dp
!       C(1,2)=-8.0_dp/3.0_dp
!       C(2,1)=-8.0_dp/3.0_dp
!       C(2,2)=4.0_dp/3.0_dp
!       endif
!       if(n.eq.7) then
!       C(1,1)=16.0_dp/3.0_dp
!       C(1,2)=-8.0_dp/3.0_dp
!       C(2,1)=-8.0_dp/3.0_dp
!       C(2,2)=4.0_dp/3.0_dp
!       endif
!       if(n.eq.8) then
!       C(1,1)=56.0_dp/3.0_dp
!       C(1,2)=-28.0_dp/3.0_dp
!       C(2,1)=-28.0_dp/3.0_dp
!       C(2,2)=14.0_dp/3.0_dp
!       endif
!       if(n.eq.9) then
!       C(1,1)=-8.0_dp/3.0_dp
!       C(1,2)=4.0_dp/3.0_dp
!       C(2,1)=4.0_dp/3.0_dp
!       C(2,2)=-2.0_dp/3.0_dp
!       endif
!       if(n.eq.10) then
!       C(1,1)=56.0_dp/3.0_dp
!       C(1,2)=-28.0_dp/3.0_dp
!       C(2,1)=-28.0_dp/3.0_dp
!       C(2,2)=14.0_dp/3.0_dp
!       endif
!       if(n.eq.11) then
!       C(1,1)=16.0_dp/3.0_dp
!       C(1,2)=-8.0_dp/3.0_dp
!       C(2,1)=-8.0_dp/3.0_dp
!       C(2,2)=4.0_dp/3.0_dp
!       endif
!       if(n.eq.12) then
!       C(1,1)=-8.0_dp/3.0_dp
!       C(1,2)=4.0_dp/3.0_dp
!       C(2,1)=4.0_dp/3.0_dp
!       C(2,2)=-2.0_dp/3.0_dp
!       endif

      endif


      end subroutine


      subroutine cc_qq_ttgam(ixq,C)
      character, intent(in) :: ixq*2
      real(dp), intent(out) :: C(2,2)

      C = 0.0_dp

      if (ixq.eq.'up') then

       C(1,1)  =  8.0_dp
       C(1,2)  =  8.0_dp*Q_up/Q_top
       C(2,1)  =  8.0_dp*Q_up/Q_top
       C(2,2)  =  8.0_dp*Q_up**2/Q_top**2

       elseif(ixq.eq.'dn') then

       C(1,1)  =  8.0_dp
       C(1,2)  =  8.0_dp*Q_dn/Q_top
       C(2,1)  =  8.0_dp*Q_dn/Q_top
       C(2,2)  =  8.0_dp*Q_dn**2/Q_top**2

       endif

      end subroutine




end  module
