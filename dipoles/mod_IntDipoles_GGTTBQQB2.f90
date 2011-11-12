module ModIntDipoles_GGTTBQQB
use ModParameters
use ModProcess
use ModMisc
implicit none


public:: EvalIntDipoles_GGTTBQQB
integer, parameter  :: dp = selected_real_kind(15)


contains


      subroutine EvalIntDipoles_GGTTBQQB(p,MomDK,z,res)
      use ModIntDipoles
      use ModAmplitudes
      implicit none
      real(dp), intent(in) :: p(4,5),MomDK(1:4,1:6),z
      real(dp), intent(out) :: res(1:3)
      type(Particle),save :: ExtParticles(1:5),ExtParticles2(1:5)
      type(TreeProcess),save :: TreeAmpsDip(1:4),TreeAmpsDip1(1:6)
      complex(dp) :: Bm(1:4),Bm1(1:6)
      complex(dp) :: AM(1:4,1:4),AM1(1:6,1:6)
      integer :: i1,i2,i3,i4,hel(5),i,j,i5,i6,emi,n, in1, in2,iTree
      real(dp) :: dipsoft, dipfini, dipplus,AP(1:3),TR
      real(dp) :: C(4,4), C1(6,6)
      real(dp) :: mtrs,mtrs0,mtrs1,mtrs2,epcorr, q(4,5)
      character :: diptype*2
      integer :: Nmax(2)
      logical, save :: first_time = .true.
      logical, save :: first_time2= .true.
      real(8) :: MadGraph_tree,MG_MOM(0:3,1:5)

      TR=0.5d0
      res(1:3) = zero



   do n=1,5

      if(n.eq.1) then

      if( first_time ) then
          call InitProcess_TbTGGG(ExtParticles(1:5))
          call InitTrees(2,3,6,TreeAmpsDip1)

          TreeAmpsDip1(1)%PartRef(1:5) = (/1,2,3,4,5/)
          TreeAmpsDip1(2)%PartRef(1:5) = (/1,2,3,5,4/)
          TreeAmpsDip1(3)%PartRef(1:5) = (/1,2,4,3,5/)
          TreeAmpsDip1(4)%PartRef(1:5) = (/1,2,4,5,3/)
          TreeAmpsDip1(5)%PartRef(1:5) = (/1,2,5,3,4/)
          TreeAmpsDip1(6)%PartRef(1:5) = (/1,2,5,4,3/)

          do iTree=1,6
            call LinkTreeParticles(TreeAmpsDip1(iTree),ExtParticles(1:5))
          enddo

          first_time=.false.
      endif


      AM1 = (zero,zero) ! first, the ``amplitude matrix'' is set to zero

        Nmax = 1
        if (TopDecays.ge.1) then
            Nmax(1) = -1
            Nmax(2) = -1
        endif


      do i1=-1,Nmax(1),2
         do i2=-1,Nmax(2),2
            do i3 = -1,1,2
               do i4 = -1,1,2
                 do i5 = -1,1,2
                 hel(1) = i1
                 hel(2) = i2
                 hel(3) = i3
                 hel(4) = i4
                 hel(5) = i5

                  call GenerateEvent52(p,momDK(1:4,1:6),hel,ExtParticles(1:5))

                  do i6 = 1,6
                  call EvalTree2(TreeAmpsDip1(i6),Bm1(i6))
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
      call cc_ttggg(C1)  ! this is the color correlation matrix for ttggg
      do i=1,6
      do j=1,6
          mtrs = mtrs + C1(i,j)*real(AM1(i,j),dp)
      enddo
      enddo
     mtrs = mtrs * alpha_s4Pi**3 /4d0/64d0   ! spin & color avg.
     mtrs0 = mtrs

      dipsoft =fi_qg(zero,zero,p,5,3,z,1)    *nf_light
      dipfini =fi_qg(zero,zero,p,5,3,z,2)    *nf_light
      dipplus =fi_qg(zero,zero,p,5,3,z,3)    *nf_light
      diptype ='fi'
      emi = 1 ! final state dipole


      else   ! if n <> 1, then the born structure is ttqqg

      if( first_time2 ) then
! init external particles
      call InitProcess_TbTQbQG(ExtParticles2(1:5))
! init tree processes for 0-> tb t qq g
      call InitTrees(4,1,4,TreeAmpsDip)

      TreeAmpsDip(1)%PartRef(1:5) = (/1,2,3,4,5/)
      TreeAmpsDip(2)%PartRef(1:5) = (/1,5,2,3,4/)
      TreeAmpsDip(3)%PartRef(1:5) = (/1,2,5,3,4/)
      TreeAmpsDip(4)%PartRef(1:5) = (/1,2,3,5,4/)


      do iTree=1,4
      call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles2(1:5))
      enddo

          first_time2=.false.
      endif

         Nmax = 1

        if (TopDecays.ge.1) then
            Nmax(1) = -1
            Nmax(2) = -1
        endif


      AM = (zero,zero) ! first, the ``amplitude matrix'' is set to zero

      do i1=-1,Nmax(1),2
         do i2=-1,Nmax(2),2
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

      call GenerateEventttqqg2(q,momDK(1:4,1:6),hel,ExtParticles2(1:5))

      do i6 = 1,4
      call EvalTree2(TreeAmpsDip(i6),Bm(i6))
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

        call cc_ttqqg(C)
        mtrs = zero
      do i=1,4
      do j=1,4
    mtrs = mtrs + C(i,j)*real(AM(i,j),dp)
      enddo
      enddo
     mtrs = mtrs * alpha_s4Pi**3 /4d0/24d0  ! spin & color avg.


!       if(n.eq.2) then
!           MG_MOM(0:3,1) =-p(1:4,3)*100d0
!           MG_MOM(0:3,2) =-p(1:4,4)*100d0
!           MG_MOM(0:3,3) = p(1:4,2)*100d0
!           MG_MOM(0:3,4) = p(1:4,1)*100d0
!           MG_MOM(0:3,5) = p(1:4,5)*100d0
!           call coupsm(0)
!           call SUG_TTBU(MG_MOM,MadGraph_tree)
!           MadGraph_tree=MadGraph_tree*100d2
!           mtrs=mtrs*RunAlphaS(NLOParam,MuRen)**3
!           print *, n,mtrs/MadGraph_tree
!           pause
!       elseif(n.eq.3) then
!           MG_MOM(0:3,1) =-p(1:4,3)*100d0
!           MG_MOM(0:3,2) =-p(1:4,4)*100d0
!           MG_MOM(0:3,3) = p(1:4,2)*100d0
!           MG_MOM(0:3,4) = p(1:4,1)*100d0
!           MG_MOM(0:3,5) = p(1:4,5)*100d0
!           call coupsm(0)
!           call SGU_TTBU(MG_MOM,MadGraph_tree)
!           MadGraph_tree=MadGraph_tree*100d2
!           mtrs=mtrs*RunAlphaS(NLOParam,MuRen)**3
!           print *, n,mtrs/MadGraph_tree
!           pause
!       elseif(n.eq.4) then
!           MG_MOM(0:3,1) =-p(1:4,3)*100d0
!           MG_MOM(0:3,2) =-p(1:4,4)*100d0
!           MG_MOM(0:3,3) = p(1:4,2)*100d0
!           MG_MOM(0:3,4) = p(1:4,1)*100d0
!           MG_MOM(0:3,5) = p(1:4,5)*100d0
!           call coupsm(0)
!           call SUBG_TTBUB(MG_MOM,MadGraph_tree)
!           MadGraph_tree=MadGraph_tree*100d2
!           mtrs=mtrs*RunAlphaS(NLOParam,MuRen)**3
!           print *, n,mtrs/MadGraph_tree
!           pause
!       elseif(n.eq.5) then
!           MG_MOM(0:3,1) =-p(1:4,3)*100d0
!           MG_MOM(0:3,2) =-p(1:4,4)*100d0
!           MG_MOM(0:3,3) = p(1:4,2)*100d0
!           MG_MOM(0:3,4) = p(1:4,1)*100d0
!           MG_MOM(0:3,5) = p(1:4,5)*100d0
!           call coupsm(0)
!           call SGUB_TTBUB(MG_MOM,MadGraph_tree)
!           MadGraph_tree=MadGraph_tree*100d2
!           mtrs=mtrs*RunAlphaS(NLOParam,MuRen)**3
!           print *, n,mtrs/MadGraph_tree
!           pause
!       endif


      if(n.eq.2) then
      dipsoft =ii_gq(zero,zero,p,3,4,z,1)    *nf_light *2d0  ! factor 2d0 fixes normalization for ii_gq dipole
      dipfini =ii_gq(zero,zero,p,3,4,z,2)    *nf_light *2d0
      dipplus =ii_gq(zero,zero,p,3,4,z,3)    *nf_light *2d0
      diptype ='ii'
      emi = 1 ! mom #3 is emitting
      endif
      if(n.eq.3) then
      dipsoft =ii_gq(zero,zero,p,4,3,z,1)    *nf_light *2d0
      dipfini =ii_gq(zero,zero,p,4,3,z,2)    *nf_light *2d0
      dipplus =ii_gq(zero,zero,p,4,3,z,3)    *nf_light *2d0
      diptype ='ii'
      emi = 2 ! mom #4 is emitting
      endif
      if(n.eq.4) then
      dipsoft =ii_gq(zero,zero,p,3,4,z,1)    *nf_light *2d0
      dipfini =ii_gq(zero,zero,p,3,4,z,2)    *nf_light *2d0
      dipplus =ii_gq(zero,zero,p,3,4,z,3)    *nf_light *2d0
      diptype ='ii'
      emi = 1 ! mom #3 is emitting
      endif
      if(n.eq.5) then
      dipsoft =ii_gq(zero,zero,p,4,3,z,1)    *nf_light *2d0
      dipfini =ii_gq(zero,zero,p,4,3,z,2)    *nf_light *2d0
      dipplus =ii_gq(zero,zero,p,4,3,z,3)    *nf_light *2d0
      diptype ='ii'
      emi = 2 ! mom #4 is emitting
      endif

      endif ! continuation of n = 1 condition



      if(emi.eq.1) then
        res(1) = res(1) + (dipsoft-dipplus)*mtrs
        res(2) = res(2) + (dipfini+dipplus)*mtrs
        mtrs1=mtrs
      endif
      if(emi.eq.2) then
        res(1) = res(1) + (dipsoft-dipplus)*mtrs
        res(3) = res(3) + (dipfini+dipplus)*mtrs
        mtrs2=mtrs
      endif
      if(emi.eq.3) then
        res(1) = res(1) + (dipsoft-dipplus)*mtrs
        res(2) = res(2) + (dipfini+dipplus)*0.5_dp*mtrs
        res(3) = res(3) + (dipfini+dipplus)*0.5_dp*mtrs
      endif
!       RES(1) = RES(1) +  DIPSOFT*MTRS  ! FOR DELTA-FCT. CHECK

   enddo
   res(1:3) = alpha_sOver2Pi*0.5d0 * res(1:3)  ! overall normalization: alpha_s/4pi
!  print * , "eval only 1 dipole + eval only delta part"

! !        epcorr=epinv+2d0*dlog(renscale/facscale)
     epcorr=epinv
     AP(1)= 0d0
     AP(2)= TR*(z**2+(1d0-z)**2) *2d0 ! factor two is for q and qb contribution
     AP(3)= 0d0
     AP(1:3) = AP(1:3) * alpha_sOver2Pi * epcorr  * nf_light

        res(1) = res(1) + 2d0*(AP(1)-AP(3))
        res(2) = res(2) + (AP(2) + AP(3))* mtrs1
        res(3) = res(3) + (AP(2) + AP(3))* mtrs2
!         RES(1) = RES(1) + 2D0*AP(1)  ! FOR DELTA-FCT. CHECK


  end subroutine







      subroutine MomTrans(n,p,q)   ! subroutine to transform momenta for required integrated dipoles
        real(dp), intent(in) :: p(4,5)
        integer, intent(in) :: n
        real(dp), intent(out) :: q(4,5)

          if (n.eq.2) then   ! qbar g q
             q = p
             q(:,4) = p(:,5)
             q(:,5) = p(:,4)
          endif

          if (n.eq.3) then   ! q g qbar
             q = p
             q(:,3) = p(:,4)
             q(:,4) = p(:,5)
             q(:,5) = p(:,3)
          endif

          if (n.eq.4) then   ! g qbar q
             q = p
             q(:,3) = p(:,5)
             q(:,4) = p(:,3)
             q(:,5) = p(:,4)
          endif

          if (n.eq.5) then   ! g q qbar
             q = p
             q(:,3) = p(:,5)
             q(:,5) = p(:,3)
          endif

      end subroutine MomTrans







end  module


