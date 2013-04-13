module gg_ttqqgg_dipint
use types
use consts_dp
use mod_dipoles
use colorcorr
use ModAmplitudes
use ModMisc
implicit none
private

public:: gg_ttqqgg_dip_int


contains

      subroutine MomTrans(n,p,q)   ! subroutine to transform momenta for required integrated dipoles
        real(dp), intent(in) :: p(4,5)
        integer, intent(in) :: n
        real(dp), intent(out) :: q(4,5)

          if (n.eq.2) then
             q = p
             q(:,4) = p(:,5)
             q(:,5) = p(:,4)
          endif

          if (n.eq.3) then
             q = p
             q(:,3) = p(:,4)
             q(:,4) = p(:,5)
             q(:,5) = p(:,3)
          endif

          if (n.eq.4) then
             q = p
             q(:,3) = p(:,5)
             q(:,4) = p(:,3)
             q(:,5) = p(:,4)
          endif

          if (n.eq.5) then
             q = p
             q(:,3) = p(:,5)
             q(:,5) = p(:,3)

          endif

      end subroutine MomTrans


      subroutine EvalIntDipoles_GGTTBQQB(p,res)
      real(dp), intent(in) :: p(4,5)
      real(dp), intent(out) :: res
      complex(dp) :: cres
      type(Particle) :: ExtParticles(1:5)
      type(TreeProcess) :: TreeAmpsDip(1:4),TreeAmpsDip1(1:6)
      complex(dp) :: Bm(1:4),Bm1(1:6)
      complex(dp) :: AM(1:4,1:4),AM1(1:6,1:6)
      integer :: i1,i2,i3,i4,hel(5),i,j,i5,i6,emi,n, in1, in2,iTree
      real(dp) :: dipsoft, dipfini, dipplus
      real(dp) :: z, C(4,4), C1(6,6),fx1(-5:5), fx2(-5:5),fx1z(-5:5),fx2z(-5:5)
      real(dp) :: mtrs, q(4,5)
      character :: diptype*2

!---- only gluon PDF contribues here



      fx1  = 1.0_dp
      fx2  = 1.0_dp
      fx1z = 1.0_dp
      fx2z = 1.0_dp


      fx1(0) = 0.0_dp
      fx2(0) = 0.0_dp
      fx1z(0) = 0.0_dp
      fx2z(0) = 0.0_dp


      res = zero
      cres = (zero,zero)
      z = 0.5_dp

!    in this case different underyling born structures for different dipoles


   do n=1,5

      if(n.eq.1) then


        call cc_ttggg(C1)  ! this is the color correlation matrix for ttggg


      call InitTrees(2,3,6,TreeAmpsDip1)

      call InitProcess_TbTGGG(ExtParticles(1:5))

      TreeAmpsDip1(1)%PartRef(1:5) = (/1,2,3,4,5/)
      TreeAmpsDip1(2)%PartRef(1:5) = (/1,2,3,5,4/)
      TreeAmpsDip1(3)%PartRef(1:5) = (/1,2,4,3,5/)
      TreeAmpsDip1(4)%PartRef(1:5) = (/1,2,4,5,3/)
      TreeAmpsDip1(5)%PartRef(1:5) = (/1,2,5,3,4/)
      TreeAmpsDip1(6)%PartRef(1:5) = (/1,2,5,4,3/)


      do iTree=1,6
      call LinkTreeParticles(TreeAmpsDip1(iTree),ExtParticles(1:5))
      enddo


      AM1 = (zero,zero) ! first, the ``amplitude matrix'' is set to zero

      do i1=-1,1,2
         do i2=-1,1,2
            do i3 = -1,1,2
               do i4 = -1,1,2
                  do i5 = -1,1,2

                 hel(1) = i1
                 hel(2) = i2
                 hel(3) = i3
                 hel(4) = i4
                 hel(5) = i5

      call GenerateEvent5(p,hel,ExtParticles(1:5))

      do i6 = 1,6
      call EvalTree(TreeAmpsDip1(i6),Bm1(i6))
      enddo

       do i=1,6
          do j=1,6
       AM1(i,j) = AM1(i,j) + Bm1(i)*conjg(Bm1(j)) !color correlation matrix
         enddo
       enddo


                  enddo
                enddo
              enddo
            enddo
          enddo


        mtrs = zero

      do i=1,6
      do j=1,6
    mtrs = mtrs + C1(i,j)*real(AM1(i,j),dp)
      enddo
      enddo


      dipsoft =fi_qg(zero,zero,p,5,3,z,1)
      dipfini =fi_qg(zero,zero,p,5,3,z,2)
      dipplus =fi_qg(zero,zero,p,5,3,z,3)
      diptype ='fi'
      emi = 3 ! final state dipole


      else   ! if n <> 1, then the born structure is ttqqg


        call cc_ttqqg(C)

! init external particles
      call InitProcess_TbTQbQG(ExtParticles(1:5))
! init tree processes for 0-> tb t qq g
       call InitTrees(4,1,4,TreeAmpsDip)

      TreeAmpsDip(1)%PartRef(1:5) = (/1,2,3,4,5/)
      TreeAmpsDip(2)%PartRef(1:5) = (/1,2,3,5,4/)
      TreeAmpsDip(3)%PartRef(1:5) = (/1,2,4,3,5/)
      TreeAmpsDip(4)%PartRef(1:5) = (/1,2,4,5,3/)


      do iTree=1,4
      call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles(1:5))
      enddo



      AM = (zero,zero) ! first, the ``amplitude matrix'' is set to zero

      do i1=-1,1,2
         do i2=-1,1,2
            do i3 = -1,1,2
               do i4 = -1,1,2
                  do i5 = -1,1,2

                 hel(1) = i1
                 hel(2) = i2
                 hel(3) = i3
                 hel(4) = i4
                 hel(5) = i5

      call Momtrans(n,p,q)  ! this is the momentum transform required
                            !to order momenta in a way that is
                            ! consistent with amplitude calculus

      call GenerateEventttqqg(q,hel,ExtParticles(1:5))

      do i6 = 1,4
      call EvalTree(TreeAmpsDip(i6),Bm(i6))
      enddo

       do i=1,4
          do j=1,4
       AM(i,j) = AM(i,j) + Bm(i)*conjg(Bm(j)) !color correlation matrix
         enddo
       enddo


                  enddo
                enddo
              enddo
            enddo
          enddo


        mtrs = zero

      do i=1,4
      do j=1,4
    mtrs = mtrs + C(i,j)*real(AM(i,j),dp)
      enddo
      enddo

      if(n.eq.2) then
      dipsoft =ii_gq(zero,zero,p,3,4,z,1)
      dipfini =ii_gq(zero,zero,p,3,4,z,2)
      dipplus =ii_gq(zero,zero,p,3,4,z,3)
      diptype ='ii'
     emi = 1 ! mom #3 is emitting
      endif
      if(n.eq.3) then
      dipsoft =ii_gq(zero,zero,p,4,3,z,1)
      dipfini =ii_gq(zero,zero,p,4,3,z,2)
      dipplus =ii_gq(zero,zero,p,4,3,z,3)
      diptype ='ii'
     emi = 2 ! mom #4 is emitting
      endif
      if(n.eq.4) then
      dipsoft =ii_gq(zero,zero,p,3,4,z,1)
      dipfini =ii_gq(zero,zero,p,3,4,z,2)
      dipplus =ii_gq(zero,zero,p,3,4,z,3)
      diptype ='ii'
     emi = 1 ! mom #3 is emitting
      endif
      if(n.eq.5) then
      dipsoft =ii_gq(zero,zero,p,4,3,z,1)
      dipfini =ii_gq(zero,zero,p,4,3,z,2)
      dipplus =ii_gq(zero,zero,p,4,3,z,3)
      diptype ='ii'
     emi = 2 ! mom #4 is emitting
      endif


      endif ! continuation of n = 1 condition


      if(emi.eq.1) then
      res = res + (dipsoft-dipplus)*mtrs*fx1(in1)*fx2(in2) &
      +mtrs*(dipfini+dipplus)*fx1z(in1)*fx2(in2)/z
      endif
      if(emi.eq.2) then
      res = res + (dipsoft-dipplus)*mtrs*fx1(in1)*fx2(in2) &
      +mtrs*(dipfini+dipplus)*fx1(in1)*fx2z(in2)/z
      endif
      if(emi.eq.3) then
      res = res + (dipsoft-dipplus)*mtrs*fx1(in1)*fx2(in2) &
      + 0.5_dp*mtrs*(dipfini+dipplus)*fx1(in1)*fx2z(in2)/z &
      + 0.5_dp*mtrs*(dipfini+dipplus)*fx1z(in1)*fx2(in2)/z
      endif



       in1 = 0
       in2 = 0


      if(emi.eq.1) then
      res = res + (dipsoft-dipplus)*mtrs*fx1(in1)*fx2(in2) &
      +mtrs*(dipfini+dipplus)*fx1z(in1)*fx2(in2)/z
      endif
      if(emi.eq.2) then
      res = res + (dipsoft-dipplus)*mtrs*fx1(in1)*fx2(in2) &
      +mtrs*(dipfini+dipplus)*fx1(in1)*fx2z(in2)/z
      endif
      if(emi.eq.3) then
      res = res + (dipsoft-dipplus)*mtrs*fx1(in1)*fx2(in2) &
      + 0.5_dp*mtrs*(dipfini+dipplus)*fx1(in1)*fx2z(in2)/z &
      + 0.5_dp*mtrs*(dipfini+dipplus)*fx1z(in1)*fx2(in2)/z
      endif


   enddo


  end subroutine gg_ttqqgg_dip_int

end  module gg_ttqqgg_dipint


