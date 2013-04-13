module ModIntDipoles_QBGTTBQBG
use ModParameters
use ModProcess
use ModMisc
use ModIntDipoles
implicit none

public:: EvalIntDipoles_QBGTTBQBG
integer, parameter  :: dp = selected_real_kind(15)

contains


!----------------------------
!---  the momentum p that goes into that subroutine has to be order as
!     p(1,:) -- tbar
!     p(2,:) -- t
!     p(3,:) -- final
!     p(4,:) -- initial
!     p(5,:) -- initial  gluon
!----------------------------


      subroutine EvalIntDipoles_QBGTTBQBG(p,MomDK,z,res)
      use ModAmplitudes
      implicit none
      real(dp), intent(in) :: p(4,5),MomDK(1:4,1:6)
      real(dp), intent(out) :: res(1:3)
      type(Particle),save :: ExtParticles1(1:5),ExtParticles2(1:5),ExtParticles3(1:5)
      type(TreeProcess),save :: TreeAmpsDip(1:4),TreeAmpsDipg(1:4),TreeAmpsDipq(1:6)
      complex(dp) :: Bm(1:4),Bq(1:6)
      complex(dp) :: AM(1:4,1:4)
      integer :: i1,i2,i3,i4,hel(5),i,j,i5,i6,emi,n, in1, in2,iTree
      real(dp) :: dipsoft, dipfini, dipplus
      real(dp) :: z, C(4,4)
      real(dp) :: Cq(6,6)
      real(dp) :: mtrs, Nmax(2)
      real(dp) :: q(4,5),epcorr,AP(1:3),CF,CA,TR
      character :: diptype*2
      logical, save :: first_time1= .true.
      logical, save :: first_time2= .true.
      logical, save :: first_time3= .true.
      real(8) :: MadGraph_tree,MG_MOM(0:3,1:5)


      CF=4d0/3d0
      CA=3d0
      TR=0.5d0

         Nmax = 1
        if (TopDecays.ge.1) then
            Nmax(1) = -1
            Nmax(2) = -1
	 else
            Nmax(1) = +1
            Nmax(2) = +1
        endif
      res(1:3) = zero


      if( first_time1 ) then
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

      first_time1=.false.
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

      call GenerateEventttqqg2(p,momDK(1:4,1:6),hel,ExtParticles1(1:5))

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


   do n=1,20

      call cc_qg_ttqqgg(n,C)
      mtrs = zero
      do i=1,4
      do j=1,4
          mtrs = mtrs + C(i,j)*real(AM(i,j),dp)
      enddo
      enddo
      mtrs = mtrs * alpha_s4Pi**3 /4d0/3d0/8d0   ! spin & color avg.      qg-->t tb q

      if(n.eq.1) then
      dipsoft =ff_qq(zero,m_top,p,3,1,z,1)
      dipfini =ff_qq(zero,m_top,p,3,1,z,2)
      dipplus =ff_qq(zero,m_top,p,3,1,z,3)
      diptype ='ff'
      emi = 3 ! final state dipole
      endif
      if(n.eq.2) then
      dipsoft =ff_qq(zero,m_top,p,3,2,z,1)
      dipfini =ff_qq(zero,m_top,p,3,2,z,2)
      dipplus =ff_qq(zero,m_top,p,3,2,z,3)
      diptype ='ff'
      emi = 3 ! final state dipole
      endif
      if(n.eq.3) then
      dipsoft =fi_qq(zero,zero,p,3,4,z,1)
      dipfini =fi_qq(zero,zero,p,3,4,z,2)
      dipplus =fi_qq(zero,zero,p,3,4,z,3)
      diptype ='fi'
      emi = 1 ! final state dipole
      endif
      if(n.eq.4) then
      dipsoft =fi_qq(zero,zero,p,3,5,z,1)
      dipfini =fi_qq(zero,zero,p,3,5,z,2)
      dipplus =fi_qq(zero,zero,p,3,5,z,3)
      diptype ='fi'
      emi = 2 ! final state dipole
      endif
      if(n.eq.5) then
      dipsoft =if_qq(zero,m_top,p,4,1,z,1)
      dipfini =if_qq(zero,m_top,p,4,1,z,2)
      dipplus =if_qq(zero,m_top,p,4,1,z,3)
      diptype ='if'
     emi = 1 ! mom #3 is emitting
      endif
      if(n.eq.6) then
      dipsoft =if_qq(zero,m_top,p,4,2,z,1)
      dipfini =if_qq(zero,m_top,p,4,2,z,2)
      dipplus =if_qq(zero,m_top,p,4,2,z,3)
      diptype ='if'
     emi = 1 ! mom #3 is emitting
      endif
      if(n.eq.7) then
      dipsoft =if_qq(zero,zero,p,4,3,z,1)
      dipfini =if_qq(zero,zero,p,4,3,z,2)
      dipplus =if_qq(zero,zero,p,4,3,z,3)
      diptype ='if'
     emi = 1 ! mom #3 is emitting
      endif
      if(n.eq.8) then
      dipsoft =ii_qq(zero,zero,p,4,5,z,1)
      dipfini =ii_qq(zero,zero,p,4,5,z,2)
      dipplus =ii_qq(zero,zero,p,4,5,z,3)
      diptype ='ii'
     emi = 1 ! mom #3 is emitting
      endif
      if(n.eq.9) then
      dipsoft =ff_qq(m_top,m_top,p,1,2,z,1)
      dipfini =ff_qq(m_top,m_top,p,1,2,z,2)
      dipplus =ff_qq(m_top,m_top,p,1,2,z,3)
      diptype ='ff'
      emi = 3 ! final state dipole
      endif
      if(n.eq.10) then
      dipsoft =ff_qq(m_top,zero,p,1,3,z,1)
      dipfini =ff_qq(m_top,zero,p,1,3,z,2)
      dipplus =ff_qq(m_top,zero,p,1,3,z,3)
      diptype ='ff'
      emi = 3 ! final state dipole
      endif
      if(n.eq.11) then
      dipsoft =fi_qq(m_top,zero,p,1,4,z,1)
      dipfini =fi_qq(m_top,zero,p,1,4,z,2)
      dipplus =fi_qq(m_top,zero,p,1,4,z,3)
      diptype ='fi'
      emi = 1 ! final state dipole  ! changed 2 to 1 here, KM.
      endif
      if(n.eq.12) then
      dipsoft =fi_qq(m_top,zero,p,1,5,z,1)
      dipfini =fi_qq(m_top,zero,p,1,5,z,2)
      dipplus =fi_qq(m_top,zero,p,1,5,z,3)
      diptype ='fi'
      emi = 2 ! final state dipole
      endif
      if(n.eq.13) then
      dipsoft =ff_qq(m_top,m_top,p,2,1,z,1)
      dipfini =ff_qq(m_top,m_top,p,2,1,z,2)
      dipplus =ff_qq(m_top,m_top,p,2,1,z,3)
      diptype ='ff'
      emi = 3 ! final state dipole
      endif
      if(n.eq.14) then
      dipsoft =ff_qq(m_top,zero,p,2,3,z,1)
      dipfini =ff_qq(m_top,zero,p,2,3,z,2)
      dipplus =ff_qq(m_top,zero,p,2,3,z,3)
      diptype ='ff'
      emi = 3 ! final state dipole
      endif
      if(n.eq.15) then
      dipsoft =fi_qq(m_top,zero,p,2,4,z,1)
      dipfini =fi_qq(m_top,zero,p,2,4,z,2)
      dipplus =fi_qq(m_top,zero,p,2,4,z,3)
      diptype ='fi'
      emi = 1 ! final state dipole
      endif
      if(n.eq.16) then
      dipsoft =fi_qq(m_top,zero,p,2,5,z,1)
      dipfini =fi_qq(m_top,zero,p,2,5,z,2)
      dipplus =fi_qq(m_top,zero,p,2,5,z,3)
      diptype ='fi'
      emi = 2 ! final state dipole
      endif
      if(n.eq.17) then
      dipsoft =if_gg(zero,m_top,p,5,1,z,1)
      dipfini =if_gg(zero,m_top,p,5,1,z,2)
      dipplus =if_gg(zero,m_top,p,5,1,z,3)
      diptype ='if'
     emi = 2 ! mom #4 is emitting
      endif
      if(n.eq.18) then
      dipsoft =if_gg(zero,m_top,p,5,2,z,1)
      dipfini =if_gg(zero,m_top,p,5,2,z,2)
      dipplus =if_gg(zero,m_top,p,5,2,z,3)
      diptype ='if'
     emi = 2 ! mom #4 is emitting
      endif
      if(n.eq.19) then
      dipsoft =if_gg(zero,zero,p,5,3,z,1)
      dipfini =if_gg(zero,zero,p,5,3,z,2)
      dipplus =if_gg(zero,zero,p,5,3,z,3)
      diptype ='if'
     emi = 2 ! mom #4 is emitting
      endif
      if(n.eq.20) then
      dipsoft =ii_gg(zero,zero,p,5,4,z,1)
      dipfini =ii_gg(zero,zero,p,5,4,z,2)
      dipplus =ii_gg(zero,zero,p,5,4,z,3)
      diptype ='ii'
     emi = 2 ! mom #4 is emitting
      endif



      if(emi.eq.1) then
        res(1) = res(1) + (dipsoft-dipplus)*mtrs
        res(2) = res(2) + (dipfini+dipplus)*mtrs
      endif
      if(emi.eq.2) then
        res(1) = res(1) + (dipsoft-dipplus)*mtrs
        res(3) = res(3) + (dipfini+dipplus)*mtrs
      endif
      if(emi.eq.3) then
        res(1) = res(1) + (dipsoft-dipplus)*mtrs
        res(2) = res(2) + (dipfini+dipplus)*0.5_dp*mtrs
        res(3) = res(3) + (dipfini+dipplus)*0.5_dp*mtrs
      endif
!       RES(1) = RES(1) +  DIPSOFT*MTRS  ! FOR DELTA-FCT. CHECK
!       print *, n,DIPSOFT,MTRS

   enddo
   res(1:3) = res(1:3) * alpha_sOver2Pi*0.5d0   ! overall normalization: alpha_s/4pi

        if( epinv.ne.0d0) then
      ! !        epcorr=epinv+2d0*dlog(renscale/facscale)
              epcorr=epinv
              call cc_ttqqg(C)
              mtrs = zero
              do i=1,4
              do j=1,4
                  mtrs = mtrs + C(i,j)*real(AM(i,j),dp)
              enddo
              enddo
              mtrs = mtrs * alpha_s4Pi**3 /4d0/3d0/8d0   ! spin & color avg.  ! check this against MG!!!

              AP(1)= 3d0/2d0*CF
              AP(2)= (-1d0-z)*CF
              AP(3)= 2d0*CF/(1d0-z)
              AP(1:3) = AP(1:3) * mtrs * alpha_sOver2Pi *epcorr
              res(1) = res(1) + (AP(1)-AP(3))
              res(2) = res(2) + (AP(2)+AP(3))
        !        RES(1) = RES(1) +  AP(1)

              AP(1)= (11d0/6d0*CA - nf_light/3.0)
              AP(2)= CA * two*(one/z+z*(one-z)-two)
              AP(3)= CA * two/(one-z)
              AP(1:3) = AP(1:3) * alpha_sOver2Pi * epcorr * mtrs
              res(1) = res(1) + (AP(1)-AP(3))
              res(3) = res(3) + (AP(2)+AP(3))
        !        RES(1) = RES(1) +  AP(1)
        endif


      ! now add additional two dipoles for collinear singularities associated with  collinear quark

       ! dipole # 21
      if( first_time2 ) then
      call InitTrees(2,3,6,TreeAmpsDipq)  ! gg -> ttb g

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

      first_time2=.false.
      endif

        call cc_ttggg(Cq)

!---- we need to rearrange momenta here
        q = p
        q(:,3) = p(:,4)
        q(:,4) = p(:,5)
        q(:,5) = p(:,3)


      mtrs = zero
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
                 call GenerateEvent52(q,momDK(1:4,1:6),hel,ExtParticles2(1:5))

      do i6 = 1,6
      call EvalTree2(TreeAmpsDipq(i6),Bq(i6))
      enddo

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
       mtrs = mtrs * alpha_s4Pi**3 /4d0/8d0/8d0   ! spin & color avg.

!           MG_MOM(0:3,1) =-p(1:4,4)*100d0
!           MG_MOM(0:3,2) =-p(1:4,5)*100d0
!           MG_MOM(0:3,3) = p(1:4,2)*100d0
!           MG_MOM(0:3,4) = p(1:4,1)*100d0
!           MG_MOM(0:3,5) = p(1:4,3)*100d0
!           call coupsm(0)
!           call SGG_TTBG(MG_MOM,MadGraph_tree)
!           MadGraph_tree=MadGraph_tree*100d2
!           mtrs=mtrs*RunAlphaS(NLOParam,MuRen)**3
!           print *, "MG ratio1",mtrs/MadGraph_tree

      dipsoft =ii_qg(zero,zero,p,4,5,z,1)  *2d0*CF
      dipfini =ii_qg(zero,zero,p,4,5,z,2)  *2d0*CF
      dipplus =ii_qg(zero,zero,p,4,5,z,3)  *2d0*CF
      diptype ='ii'
      emi = 1 ! mom #3 is emitting
      if(emi.eq.1) then
        res(1) = res(1) + (dipsoft-dipplus)*mtrs * alpha_sOver2Pi*0.5d0   ! overall normalization: alpha_s/4pi
        res(2) = res(2) + (dipfini+dipplus)*mtrs * alpha_sOver2Pi*0.5d0   ! overall normalization: alpha_s/4pi
      endif

      if( epinv.ne.0d0) then
          AP(1)= 0d0
          AP(2)= CF * (1d0+(1d0-z)**2)/z * mtrs * alpha_sOver2Pi *epcorr
          AP(3)= 0d0
          res(1) = res(1) + (AP(1)-AP(3))
          res(2) = res(2) + (AP(2)+AP(3))
      endif


!------ the dipole number 22

        call cc_ttqqg(C)

      if( first_time3 ) then
! init external particles
      call InitProcess_TbTQbQG(ExtParticles3(1:5))   ! q qb -> t tb g
! init tree processes for 0-> tb t qq g

       call InitTrees(4,1,4,TreeAmpsDipg)
      TreeAmpsDipg(1)%PartRef(1:5) = (/1,2,3,4,5/)
      TreeAmpsDipg(2)%PartRef(1:5) = (/1,5,2,3,4/)
      TreeAmpsDipg(3)%PartRef(1:5) = (/1,2,5,3,4/)
      TreeAmpsDipg(4)%PartRef(1:5) = (/1,2,3,5,4/)

      do iTree=1,4
      call LinkTreeParticles(TreeAmpsDipg(iTree),ExtParticles3(1:5))
      enddo

      first_time3=.false.
      endif


      Bm = (zero,zero) ! first, the ``amplitude matrix'' is set to zero

      q = p
      q(:,3) = p(:,5)
      q(:,5) = p(:,3)


      mtrs = zero
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


      call GenerateEventttqqg2(q,momDK(1:4,1:6),hel,ExtParticles3(1:5))

      do i6 = 1,4
      call EvalTree2(TreeAmpsDipg(i6),Bm(i6))
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
       mtrs = mtrs * alpha_s4Pi**3 /4d0/3d0/3d0   ! spin & color avg.

!           MG_MOM(0:3,1) =-p(1:4,5)*100d0
!           MG_MOM(0:3,2) =-p(1:4,4)*100d0
!           MG_MOM(0:3,3) = p(1:4,2)*100d0
!           MG_MOM(0:3,4) = p(1:4,1)*100d0
!           MG_MOM(0:3,5) = p(1:4,3)*100d0
!           call coupsm(0)
!           call SUUB_TTBG(MG_MOM,MadGraph_tree)
!           MadGraph_tree=MadGraph_tree*100d2
!           mtrs=mtrs*RunAlphaS(NLOParam,MuRen)**3
!           print *, "MG ratio2",mtrs/MadGraph_tree
!           pause


      dipsoft =ii_gq(zero,zero,p,5,4,z,1) *2d0  ! factor 2d0 fixes normalization for ii_gq dipole
      dipfini =ii_gq(zero,zero,p,5,4,z,2) *2d0
      dipplus =ii_gq(zero,zero,p,5,4,z,3) *2d0
      diptype ='ii'
      emi = 2 ! mom #4 is emitting
      if(emi.eq.2) then
        res(1) = res(1) + (dipsoft-dipplus)*mtrs * alpha_sOver2Pi*0.5d0   ! overall normalization: alpha_s/4pi
        res(3) = res(3) + (dipfini+dipplus)*mtrs * alpha_sOver2Pi*0.5d0   ! overall normalization: alpha_s/4pi
      endif

       if( epinv.ne.0d0) then
          AP(1)= 0d0
          AP(2)= TR*(z**2+(1d0-z)**2) * mtrs * alpha_sOver2Pi *epcorr
          AP(3)= 0d0
          res(1) = res(1) + (AP(1)-AP(3))
          res(3) = res(3) + (AP(2)+AP(3))
       endif


!       print *, "res(2)",res(2)
!       print *, "res(3)",res(3)
!       stop


  end subroutine





end  module


