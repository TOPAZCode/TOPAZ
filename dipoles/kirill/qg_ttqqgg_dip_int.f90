module qg_ttqqgg_dipint
use types
use consts_dp
use mod_dipoles
use colorcorr
use ModAmplitudes
use ModMisc
implicit none
private

public:: qg_ttqqgg_dip_int


contains


!---  the momentum p that goes into that subroutine has to be order as
!     p(1,:) -- tbar
!     p(2,:) -- t
!     p(3,:) -- initial
!     p(4,:) -- final
!     p(5,:) -- initial
!----------------------------

      subroutine qg_ttqqgg_dip_int(p,res)
      real(dp), intent(in) :: p(4,5)
      real(dp), intent(out) :: res
      complex(dp) :: cres
      type(Particle) :: ExtParticles1(1:5),ExtParticles2(1:5),ExtParticles3(1:5)
      type(TreeProcess) :: TreeAmpsDip(1:4),TreeAmpsDipg(1:4),TreeAmpsDipq(1:6)
      complex(dp) :: Bm(1:4),Bq(1:6)
      complex(dp) :: AM(1:4,1:4)
      integer :: i1,i2,i3,i4,hel(5),i,j,i5,i6,emi,n, in1, in2,iTree
      real(dp) :: dipsoft, dipfini, dipplus
      real(dp) :: z, C(4,4), fx1(-5:5), fx2(-5:5),fx1z(-5:5),fx2z(-5:5)
      real(dp) :: Cq(6,6)
      real(dp) :: mtrs
      real(dp) :: q(4,5)
      character :: diptype*2


      res = zero
      cres = (zero,zero)

! init external particles
      call InitProcess_TbTQbQG(ExtParticles1(1:5))
! init tree processes for 0-> tb t bq q g g
       call InitTrees(4,1,4,TreeAmpsDip)

      TreeAmpsDip(1)%PartRef(1:5) = (/1,2,3,4,5/)
      TreeAmpsDip(2)%PartRef(1:5) = (/1,5,2,3,4/)
      TreeAmpsDip(3)%PartRef(1:5) = (/1,2,5,3,4/)
      TreeAmpsDip(4)%PartRef(1:5) = (/1,2,3,5,4/)


      do iTree=1,4
      call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles1(1:5))
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

      call GenerateEventttqqg(p,hel,ExtParticles1(1:5))

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


   do n=1,20



        call cc_qg_ttqqgg(n,C)

        mtrs = zero

      do i=1,4
      do j=1,4
    mtrs = mtrs + C(i,j)*real(AM(i,j),dp)
      enddo
      enddo


      if(n.eq.1) then
      dipsoft =if_qq(zero,mt,p,3,1,z,1)
      dipfini =if_qq(zero,mt,p,3,1,z,2)
      dipplus =if_qq(zero,mt,p,3,1,z,3)
      diptype ='if'
     emi = 1 ! mom #3 is emitting
      endif
      if(n.eq.2) then
      dipsoft =if_qq(zero,mt,p,3,2,z,1)
      dipfini =if_qq(zero,mt,p,3,2,z,2)
      dipplus =if_qq(zero,mt,p,3,2,z,3)
      diptype ='if'
     emi = 1 ! mom #3 is emitting
      endif
      if(n.eq.3) then
      dipsoft =if_qq(zero,zero,p,3,4,z,1)
      dipfini =if_qq(zero,zero,p,3,4,z,2)
      dipplus =if_qq(zero,zero,p,3,4,z,3)
      diptype ='if'
     emi = 1 ! mom #3 is emitting
      endif
      if(n.eq.4) then
      dipsoft =ii_qq(zero,zero,p,3,5,z,1)
      dipfini =ii_qq(zero,zero,p,3,5,z,2)
      dipplus =ii_qq(zero,zero,p,3,5,z,3)
      diptype ='ii'
     emi = 1 ! mom #3 is emitting
      endif
      if(n.eq.5) then
      dipsoft =ff_qq(zero,mt,p,4,1,z,1)
      dipfini =ff_qq(zero,mt,p,4,1,z,2)
      dipplus =ff_qq(zero,mt,p,4,1,z,3)
      diptype ='ff'
      emi = 3 ! final state dipole
      endif
      if(n.eq.6) then
      dipsoft =ff_qq(zero,mt,p,4,2,z,1)
      dipfini =ff_qq(zero,mt,p,4,2,z,2)
      dipplus =ff_qq(zero,mt,p,4,2,z,3)
      diptype ='ff'
      emi = 3 ! final state dipole
      endif
      if(n.eq.7) then
      dipsoft =fi_qq(zero,zero,p,4,3,z,1)
      dipfini =fi_qq(zero,zero,p,4,3,z,2)
      dipplus =fi_qq(zero,zero,p,4,3,z,3)
      diptype ='fi'
      emi = 3 ! final state dipole
      endif
      if(n.eq.8) then
      dipsoft =fi_qq(zero,zero,p,4,5,z,1)
      dipfini =fi_qq(zero,zero,p,4,5,z,2)
      dipplus =fi_qq(zero,zero,p,4,5,z,3)
      diptype ='fi'
      emi = 3 ! final state dipole
      endif
      if(n.eq.9) then
      dipsoft =ff_qq(mt,mt,p,1,2,z,1)
      dipfini =ff_qq(mt,mt,p,1,2,z,2)
      dipplus =ff_qq(mt,mt,p,1,2,z,3)
      diptype ='ff'
      emi = 3 ! final state dipole
      endif
      if(n.eq.10) then
      dipsoft =fi_qq(mt,zero,p,1,3,z,1)
      dipfini =fi_qq(mt,zero,p,1,3,z,2)
      dipplus =fi_qq(mt,zero,p,1,3,z,3)
      diptype ='fi'
      emi = 3 ! final state dipole
      endif
      if(n.eq.11) then
      dipsoft =ff_qq(mt,zero,p,1,4,z,1)
      dipfini =ff_qq(mt,zero,p,1,4,z,2)
      dipplus =ff_qq(mt,zero,p,1,4,z,3)
      diptype ='ff'
      emi = 3 ! final state dipole
      endif
      if(n.eq.12) then
      dipsoft =fi_qq(mt,zero,p,1,5,z,1)
      dipfini =fi_qq(mt,zero,p,1,5,z,2)
      dipplus =fi_qq(mt,zero,p,1,5,z,3)
      diptype ='fi'
      emi = 3 ! final state dipole
      endif
      if(n.eq.13) then
      dipsoft =ff_qq(mt,mt,p,2,1,z,1)
      dipfini =ff_qq(mt,mt,p,2,1,z,2)
      dipplus =ff_qq(mt,mt,p,2,1,z,3)
      diptype ='ff'
      emi = 3 ! final state dipole
      endif
      if(n.eq.14) then
      dipsoft =fi_qq(mt,zero,p,2,3,z,1)
      dipfini =fi_qq(mt,zero,p,2,3,z,2)
      dipplus =fi_qq(mt,zero,p,2,3,z,3)
      diptype ='fi'
      emi = 3 ! final state dipole
      endif
      if(n.eq.15) then
      dipsoft =ff_qq(mt,zero,p,2,4,z,1)
      dipfini =ff_qq(mt,zero,p,2,4,z,2)
      dipplus =ff_qq(mt,zero,p,2,4,z,3)
      diptype ='ff'
      emi = 3 ! final state dipole
      endif
      if(n.eq.16) then
      dipsoft =fi_qq(mt,zero,p,2,5,z,1)
      dipfini =fi_qq(mt,zero,p,2,5,z,2)
      dipplus =fi_qq(mt,zero,p,2,5,z,3)
      diptype ='fi'
      emi = 3 ! final state dipole
      endif
      if(n.eq.17) then
      dipsoft =if_gg(zero,mt,p,5,1,z,1)
      dipfini =if_gg(zero,mt,p,5,1,z,2)
      dipplus =if_gg(zero,mt,p,5,1,z,3)
      diptype ='if'
     emi = 2 ! mom #4 is emitting
      endif
      if(n.eq.18) then
      dipsoft =if_gg(zero,mt,p,5,2,z,1)
      dipfini =if_gg(zero,mt,p,5,2,z,2)
      dipplus =if_gg(zero,mt,p,5,2,z,3)
      diptype ='if'
     emi = 2 ! mom #4 is emitting
      endif
      if(n.eq.19) then
      dipsoft =ii_gg(zero,zero,p,5,3,z,1)
      dipfini =ii_gg(zero,zero,p,5,3,z,2)
      dipplus =ii_gg(zero,zero,p,5,3,z,3)
      diptype ='ii'
     emi = 2 ! mom #4 is emitting
      endif
      if(n.eq.20) then
      dipsoft =if_gg(zero,zero,p,5,4,z,1)
      dipfini =if_gg(zero,zero,p,5,4,z,2)
      dipplus =if_gg(zero,zero,p,5,4,z,3)
      diptype ='if'
     emi = 2 ! mom #4 is emitting
      endif


      do in1 = 1,5
         do in2 = -5,-1

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
       enddo

   enddo  ! end of loop over 20(!) dipoles




      ! now add additional two dipoles for collinear singularities associated with  collinear quark

       ! dipole # 21

      call InitTrees(2,3,6,TreeAmpsDipq)

      call InitProcess_TbTGGG(ExtParticles2(1:5))
      TreeAmpsDipq(1)%PartRef(1:5) = (/1,2,3,4,5/)
      TreeAmpsDipq(2)%PartRef(1:5) = (/1,2,3,5,4/)
      TreeAmpsDipq(3)%PartRef(1:5) = (/1,2,4,3,5/)
      TreeAmpsDipq(4)%PartRef(1:5) = (/1,2,4,5,3/)
      TreeAmpsDipq(5)%PartRef(1:5) = (/1,2,5,3,4/)
      TreeAmpsDipq(6)%PartRef(1:5) = (/1,2,5,4,3/)



      do iTree=1,6
      call LinkTreeParticles(TreeAmpsDipq(iTree),ExtParticles2(1:5))
      enddo

!----- color factor

        call cc_ttggg(Cq)

!---- we need to rearrange momenta here

        q = p
        q(:,4) = p(:,5)
        q(:,5) = p(:,4)

      call GenerateEvent5(q,hel,ExtParticles2(1:5))


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

      do i6 = 1,6
      call EvalTree(TreeAmpsDipq(i6),Bq(i6))
      enddo


        mtrs = zero

       do i=1,6
       do j=1,6
         mtrs = mtrs + Cq(i,j)*real(Bq(i)*conjg(Bq(j)),dp)
       enddo
       enddo

                 enddo
               enddo
             enddo
           enddo
         enddo


      dipsoft =ii_qg(zero,zero,p,3,5,z,1)
      dipfini =ii_qg(zero,zero,p,3,5,z,2)
      dipplus =ii_qg(zero,zero,p,3,5,z,3)
      diptype ='ii'
      emi = 1 ! mom #3 is emitting



      do in1 = 1,5
         do in2 = -5,-1

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
       enddo

!------ the dipole number 22

        call cc_ttqqg(C)

! init external particles
      call InitProcess_TbTQbQG(ExtParticles3(1:5))
! init tree processes for 0-> tb t qq g

       call InitTrees(4,1,4,TreeAmpsDipg)

      TreeAmpsDipg(1)%PartRef(1:5) = (/1,2,3,4,5/)
      TreeAmpsDipg(2)%PartRef(1:5) = (/1,5,2,3,4/)
      TreeAmpsDipg(3)%PartRef(1:5) = (/1,2,5,3,4/)
      TreeAmpsDipg(4)%PartRef(1:5) = (/1,2,3,5,4/)

      do iTree=1,4
      call LinkTreeParticles(TreeAmpsDipg(iTree),ExtParticles3(1:5))
      enddo



      Bm = (zero,zero) ! first, the ``amplitude matrix'' is set to zero

      q = p
      q(:,4) = p(:,5)
      q(:,5) = p(:,4)

        mtrs = zero

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


      call GenerateEventttqqg(q,hel,ExtParticles3(1:5))

      do i6 = 1,4
      call EvalTree(TreeAmpsDipg(i6),Bm(i6))
      enddo

       do i=1,4
          do j=1,4
       mtrs  = mtrs + C(i,j)*real(Bm(i)*conjg(Bm(j)),dp) !color correlation matrix
         enddo
       enddo


                  enddo
                enddo
              enddo
            enddo
          enddo

      dipsoft =ii_gq(zero,zero,p,5,3,z,1)
      dipfini =ii_gq(zero,zero,p,5,3,z,2)
      dipplus =ii_gq(zero,zero,p,5,3,z,3)
      diptype ='ii'
      emi = 2 ! mom #4 is emitting


      do in1 = 1,5
         do in2 = -5,-1

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
       enddo

  end subroutine

end  module


