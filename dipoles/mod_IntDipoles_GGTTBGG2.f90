module ModIntDipoles_GGTTBGG
use ModParameters
use ModProcess
use ModMisc
implicit none


public:: EvalIntDipoles_GGTTBGG
integer, parameter  :: dp = selected_real_kind(15)


!       double precision, private, parameter :: NCol=3d0
!       double precision, private, parameter :: TR=0.5d0
!       double precision, private, parameter :: CA=2d0*TR*NCol
!       double precision, private, parameter :: CF=TR*(NCol**2-1d0)/NCol

contains


!-------------------------------------------------------------
!    ordering of the momenta in the p-array:
!     outgoing convention: bar t, t 3= incoming, 4=incoming, 5=final state

      subroutine EvalIntDipoles_GGTTBGG(p,MomDK,z,res)
      use ModIntDipoles
      use ModAmplitudes
      implicit none
      real(dp), intent(in) :: p(4,5),MomDK(1:4,1:6),z
      real(dp), intent(out) :: res(1:3)
      real(dp) :: dipsoft,dipfini,dipplus,AP(1:3),CA
      type(Particle),save :: ExtParticles(1:5)
      type(TreeProcess),save :: TreeAmpsDip(1:6)
      complex(dp) :: Bm(1:6),mtrsq
      complex(dp) :: AM(1:6,1:6)
      integer :: iTree,i1,i2,i3,i4,i5,i,j,i6,n
      integer :: hel(1:5),emi , Nmax(5)
      character(2) :: diptype
      real(dp) :: C(1:6,1:6),epcorr
      logical, save :: first_time = .true.

      res(1:3) = zero
      CA=3d0


      if (first_time) then
! init tree processes for 0-> tb t g g g
      call InitTrees(2,3,6,TreeAmpsDip)
! init external particles
      call InitProcess_TbTGGG(ExtParticles(1:5))

      TreeAmpsDip(1)%PartRef(1:5) = (/1,2,3,4,5/)
      TreeAmpsDip(2)%PartRef(1:5) = (/1,2,3,5,4/)
      TreeAmpsDip(3)%PartRef(1:5) = (/1,2,4,3,5/)
      TreeAmpsDip(4)%PartRef(1:5) = (/1,2,4,5,3/)
      TreeAmpsDip(5)%PartRef(1:5) = (/1,2,5,3,4/)
      TreeAmpsDip(6)%PartRef(1:5) = (/1,2,5,4,3/)

      do iTree=1,6
          call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles(1:5))
      enddo

              first_time=.false.
      endif

      AM = (zero,zero) ! first, the ``amplitude matrix'' is set to zero

!--- after momentum mapping -- sum over colors and polarizations

         Nmax = 1

        if (TopDecays.ge.1) then
            Nmax(1) = -1
            Nmax(2) = -1
        endif

       do i1=-1,Nmax(1),2
         do i2 = -1,Nmax(2),2
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
              call EvalTree2(TreeAmpsDip(i6),Bm(i6))
          enddo

          do i=1,6
          do j=1,6
              AM(i,j) = AM(i,j) + Bm(i)*dconjg(Bm(j))
          enddo
          enddo

      enddo
      enddo
      enddo
      enddo
      enddo


   do n=1,36

     mtrsq = (zero,zero)
     call cc_gg_ttgggg(n,C)
     do i=1,6
     do j=1,6
          mtrsq = mtrsq + C(i,j)*real(AM(i,j),dp)
     enddo
     enddo
     mtrsq = mtrsq * alpha_s4Pi**3 /4d0/64d0   *0.5d0  ! spin & color avg., color convention adjustment


     if(n.eq.1) then
     dipsoft =if_gg(zero,m_Top,p,3,1,z,1)
     dipfini =if_gg(zero,m_Top,p,3,1,z,2)
     dipplus =if_gg(zero,m_Top,p,3,1,z,3)
     diptype ='if'
     emi = 1 ! mom #3 is emitting
     endif
     if(n.eq.2) then
     dipsoft =if_gg(zero,m_Top,p,3,2,z,1)
     dipfini =if_gg(zero,m_Top,p,3,2,z,2)
     dipplus =if_gg(zero,m_Top,p,3,2,z,3)
     diptype ='if'
     emi = 1 ! mom #3 is emitting
     endif
     if(n.eq.3) then
     dipsoft =ii_gg(zero,zero,p,3,4,z,1)
     dipfini =ii_gg(zero,zero,p,3,4,z,2)
     dipplus =ii_gg(zero,zero,p,3,4,z,3)
     diptype ='ii'
     emi = 1 ! mom #3 is emitting
     endif
     if(n.eq.4) then
     dipsoft =if_gg(zero,zero,p,3,5,z,1)
     dipfini =if_gg(zero,zero,p,3,5,z,2)
     dipplus =if_gg(zero,zero,p,3,5,z,3)
     diptype ='if'
     emi = 1 ! mom #3 is emitting
     endif
     if(n.eq.5) then
     dipsoft =if_gg(zero,m_Top,p,4,1,z,1)
     dipfini =if_gg(zero,m_Top,p,4,1,z,2)
     dipplus =if_gg(zero,m_Top,p,4,1,z,3)
     diptype ='if'
     emi = 2 ! mom #4 is emitting
     endif
     if(n.eq.6) then
     dipsoft =if_gg(zero,m_Top,p,4,2,z,1)
     dipfini =if_gg(zero,m_Top,p,4,2,z,2)
     dipplus =if_gg(zero,m_Top,p,4,2,z,3)
     diptype ='if'
     emi = 2 ! mom #4 is emitting
     endif
     if(n.eq.7) then
     dipsoft =ii_gg(zero,zero,p,4,3,z,1)
     dipfini =ii_gg(zero,zero,p,4,3,z,2)
     dipplus =ii_gg(zero,zero,p,4,3,z,3)
     diptype ='ii'
     emi = 2 ! mom #4 is emitting
     endif
     if(n.eq.8) then
     dipsoft =if_gg(zero,zero,p,4,5,z,1)
     dipfini =if_gg(zero,zero,p,4,5,z,2)
     dipplus =if_gg(zero,zero,p,4,5,z,3)
     diptype ='if'
     emi = 2 ! mom #4 is emitting
     endif
     if(n.eq.9) then
     dipsoft =ff_qq(m_Top,m_Top,p,1,2,z,1)
     dipfini =ff_qq(m_Top,m_Top,p,1,2,z,2)
     dipplus =ff_qq(m_Top,m_Top,p,1,2,z,3)
     diptype ='ff'
     emi = 3 ! final state is emitting
     endif
     if(n.eq.10) then
     dipsoft =fi_qq(m_Top,zero,p,1,3,z,1)
     dipfini =fi_qq(m_Top,zero,p,1,3,z,2)
     dipplus =fi_qq(m_Top,zero,p,1,3,z,3)
     diptype ='fi'
     emi = 1 ! final state is emitting
     endif
     if(n.eq.11) then
     dipsoft =fi_qq(m_Top,zero,p,1,4,z,1)
     dipfini =fi_qq(m_Top,zero,p,1,4,z,2)
     dipplus =fi_qq(m_Top,zero,p,1,4,z,3)
     diptype ='fi'
     emi = 2 ! final state is emitting
     endif
     if(n.eq.12) then
     dipsoft =ff_qq(m_Top,zero,p,1,5,z,1)
     dipfini =ff_qq(m_Top,zero,p,1,5,z,2)
     dipplus =ff_qq(m_Top,zero,p,1,5,z,3)
     diptype ='ff'
     emi = 3 ! final state is emitting
     endif
     if(n.eq.13) then
     dipsoft =ff_qq(m_Top,m_Top,p,2,1,z,1)
     dipfini =ff_qq(m_Top,m_Top,p,2,1,z,2)
     dipplus =ff_qq(m_Top,m_Top,p,2,1,z,3)
     diptype ='ff'
     emi = 3 ! final state is emitting
     endif
     if(n.eq.14) then
     dipsoft =fi_qq(m_Top,zero,p,2,3,z,1)
     dipfini =fi_qq(m_Top,zero,p,2,3,z,2)
     dipplus =fi_qq(m_Top,zero,p,2,3,z,3)
     diptype ='fi'
     emi = 1 ! final state is emitting
     endif
     if(n.eq.15) then
     dipsoft =fi_qq(m_Top,zero,p,2,4,z,1)
     dipfini =fi_qq(m_Top,zero,p,2,4,z,2)
     dipplus =fi_qq(m_Top,zero,p,2,4,z,3)
     diptype ='fi'
     emi = 2 ! final state is emitting
     endif
     if(n.eq.16) then
     dipsoft =ff_qq(m_Top,zero,p,2,5,z,1)
     dipfini =ff_qq(m_Top,zero,p,2,5,z,2)
     dipplus =ff_qq(m_Top,zero,p,2,5,z,3)
     diptype ='ff'
     emi = 3 ! final state is emitting
     endif
     if(n.eq.17) then
     dipsoft =ff_gg(zero,m_Top,p,5,1,z,1)
     dipfini =ff_gg(zero,m_Top,p,5,1,z,2)
     dipplus =ff_gg(zero,m_Top,p,5,1,z,3)
     diptype ='ff'
     emi = 3 ! final state is emitting
     endif
     if(n.eq.18) then
     dipsoft =ff_gg(zero,m_Top,p,5,2,z,1)
     dipfini =ff_gg(zero,m_Top,p,5,2,z,2)
     dipplus =ff_gg(zero,m_Top,p,5,2,z,3)
     diptype ='ff'
     emi = 3 ! final state is emitting
     endif
     if(n.eq.19) then
     dipsoft =fi_gg(zero,zero,p,5,3,z,1)
     dipfini =fi_gg(zero,zero,p,5,3,z,2)
     dipplus =fi_gg(zero,zero,p,5,3,z,3)
     diptype ='fi'
     emi = 1 ! final state is emitting
     endif
     if(n.eq.20) then
     dipsoft =fi_gg(zero,zero,p,5,4,z,1)
     dipfini =fi_gg(zero,zero,p,5,4,z,2)
     dipplus =fi_gg(zero,zero,p,5,4,z,3)
     diptype ='fi'
     emi = 2 ! final state is emitting
     endif
     if(n.eq.21) then
     dipsoft =if_gg(zero,m_Top,p,3,1,z,1)
     dipfini =if_gg(zero,m_Top,p,3,1,z,2)
     dipplus =if_gg(zero,m_Top,p,3,1,z,3)
     diptype ='if'
     emi = 1 ! mom #3 is emitting
     endif
     if(n.eq.22) then
     dipsoft =if_gg(zero,m_Top,p,3,2,z,1)
     dipfini =if_gg(zero,m_Top,p,3,2,z,2)
     dipplus =if_gg(zero,m_Top,p,3,2,z,3)
     diptype ='if'
     emi = 1 ! mom #3 is emitting
     endif
     if(n.eq.23) then
     dipsoft =ii_gg(zero,zero,p,3,4,z,1)
     dipfini =ii_gg(zero,zero,p,3,4,z,2)
     dipplus =ii_gg(zero,zero,p,3,4,z,3)
     diptype ='ii'
     emi = 1 ! mom #3 is emitting
     endif
     if(n.eq.24) then
     dipsoft =if_gg(zero,zero,p,3,5,z,1)
     dipfini =if_gg(zero,zero,p,3,5,z,2)
     dipplus =if_gg(zero,zero,p,3,5,z,3)
     diptype ='if'
    emi = 1 ! mom #3 is emitting
     endif
     if(n.eq.25) then
     dipsoft =if_gg(zero,m_Top,p,4,1,z,1)
     dipfini =if_gg(zero,m_Top,p,4,1,z,2)
     dipplus =if_gg(zero,m_Top,p,4,1,z,3)
     diptype ='if'
     emi = 2 ! mom #4 is emitting
     endif
     if(n.eq.26) then
     dipsoft =if_gg(zero,m_Top,p,4,2,z,1)
     dipfini =if_gg(zero,m_Top,p,4,2,z,2)
     dipplus =if_gg(zero,m_Top,p,4,2,z,3)
     diptype ='if'
     emi = 2 ! mom #4 is emitting
     endif
     if(n.eq.27) then
     dipsoft =ii_gg(zero,zero,p,4,3,z,1)
     dipfini =ii_gg(zero,zero,p,4,3,z,2)
     dipplus =ii_gg(zero,zero,p,4,3,z,3)
     diptype ='ii'
     emi = 2 ! mom #4 is emitting
     endif
     if(n.eq.28) then
     dipsoft =if_gg(zero,zero,p,4,5,z,1)
     dipfini =if_gg(zero,zero,p,4,5,z,2)
     dipplus =if_gg(zero,zero,p,4,5,z,3)
     diptype ='if'
     emi = 2 ! mom #4 is emitting
     endif
     if(n.eq.29) then
     dipsoft =ff_qq(m_Top,m_Top,p,1,2,z,1)
     dipfini =ff_qq(m_Top,m_Top,p,1,2,z,2)
     dipplus =ff_qq(m_Top,m_Top,p,1,2,z,3)
     diptype ='ff'
     emi = 3 ! final state is emitting
     endif
     if(n.eq.30) then
     dipsoft =fi_qq(m_Top,zero,p,1,3,z,1)
     dipfini =fi_qq(m_Top,zero,p,1,3,z,2)
     dipplus =fi_qq(m_Top,zero,p,1,3,z,3)
     diptype ='fi'
     emi = 1 ! final state is emitting
     endif
     if(n.eq.31) then
     dipsoft =fi_qq(m_Top,zero,p,1,4,z,1)
     dipfini =fi_qq(m_Top,zero,p,1,4,z,2)
     dipplus =fi_qq(m_Top,zero,p,1,4,z,3)
     diptype ='fi'
     emi = 2 ! final state is emitting
     endif
     if(n.eq.32) then
     dipsoft =ff_qq(m_Top,zero,p,1,5,z,1)
     dipfini =ff_qq(m_Top,zero,p,1,5,z,2)
     dipplus =ff_qq(m_Top,zero,p,1,5,z,3)
     diptype ='ff'
     emi = 3 ! final state is emitting
     endif
     if(n.eq.33) then
     dipsoft =ff_qq(m_Top,m_Top,p,2,1,z,1)
     dipfini =ff_qq(m_Top,m_Top,p,2,1,z,2)
     dipplus =ff_qq(m_Top,m_Top,p,2,1,z,3)
     diptype ='ff'
     emi = 3 ! final state is emitting
     endif
     if(n.eq.34) then
     dipsoft =fi_qq(m_Top,zero,p,2,3,z,1)
     dipfini =fi_qq(m_Top,zero,p,2,3,z,2)
     dipplus =fi_qq(m_Top,zero,p,2,3,z,3)
     diptype ='fi'
     emi = 1 ! final state is emitting
     endif
     if(n.eq.35) then
     dipsoft =fi_qq(m_Top,zero,p,2,4,z,1)
     dipfini =fi_qq(m_Top,zero,p,2,4,z,2)
     dipplus =fi_qq(m_Top,zero,p,2,4,z,3)
     diptype ='fi'
     emi = 2 ! final state is emitting
     endif
     if(n.eq.36) then
     dipsoft =ff_qq(m_Top,zero,p,2,5,z,1)
     dipfini =ff_qq(m_Top,zero,p,2,5,z,2)
     dipplus =ff_qq(m_Top,zero,p,2,5,z,3)
     diptype ='ff'
     emi = 3 ! final state is emitting
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

!       RES(1) = RES(1) +  DIPSOFT*MTRSQ  ! FOR DELTA-FCT. CHECK
   enddo
   res(1:3) = alpha_sOver2Pi*0.5d0 * res(1:3)  ! overall normalization: alpha_s/4pi

!    eval tree
     mtrsq = (zero,zero)
     call cc_gg_ttgggg(0,C)
     do i=1,6
     do j=1,6
          mtrsq = mtrsq + C(i,j)*AM(i,j)
     enddo
     enddo
     mtrsq = mtrsq * alpha_s4Pi**3 /4d0/64d0


! !        epcorr=epinv+2d0*dlog(renscale/facscale)
     epcorr=epinv
     AP(1)= (11.0_dp/6.0_dp*CA - nf_light/3.0_dp )
     AP(2)= CA * two*(one/z+z*(one-z)-two)
     AP(3)= CA * two/(one-z)
     AP(1:3) = AP(1:3) * alpha_sOver2Pi * epcorr * mtrsq

        res(1) = res(1) + 2d0*(AP(1)-AP(3))
        res(2) = res(2) + (AP(2) + AP(3))
        res(3) = res(3) + (AP(2) + AP(3))
!          RES(1) = RES(1) + 2D0*AP(1)  ! FOR DELTA-FCT. CHECK

!     print *, "eval only delta part"

!       res(1) = res(1)   - nf_light/3.0_dp* alpha_sOver2Pi * epcorr * mtrsq
!       print *, "xx",  - nf_light/3.0_dp* alpha_sOver2Pi * epcorr * mtrsq

! !    CHECK: color correlation
!      mtrsq = (zero,zero)
!      call cc_gg_ttgggg(0,C)
!      do i=1,6
!      do j=1,6
!           mtrsq = mtrsq + C(i,j)*real(AM(i,j),dp)
!      enddo
!      enddo
!      mtrsq = mtrsq * alpha_s4Pi**3 /4d0/64d0
!      print *, "uncorr. tree",mtrsq
! ! ---------
!
!      mtrsq = (zero,zero)
!      call cc_gg_ttgggg(12,C)
!      do i=1,6
!      do j=1,6
!           mtrsq =  + C(i,j)*real(AM(i,j),dp)
!      enddo
!      enddo
!      print *, mtrsq* alpha_s4Pi**3 /4d0/64d0
!      mtrsq = (zero,zero)
!      call cc_gg_ttgggg(16,C)
!      do i=1,6
!      do j=1,6
!           mtrsq =  + C(i,j)*real(AM(i,j),dp)
!      enddo
!      enddo
!      print *, mtrsq* alpha_s4Pi**3 /4d0/64d0
!      mtrsq = (zero,zero)
!      call cc_gg_ttgggg(19,C)
!      do i=1,6
!      do j=1,6
!           mtrsq =  + C(i,j)*real(AM(i,j),dp)
!      enddo
!      enddo
!      print *, mtrsq* alpha_s4Pi**3 /4d0/64d0
!      mtrsq = (zero,zero)
!      call cc_gg_ttgggg(20,C)
!      do i=1,6
!      do j=1,6
!           mtrsq =  + C(i,j)*real(AM(i,j),dp)
!      enddo
!      enddo
!      print *, mtrsq* alpha_s4Pi**3 /4d0/64d0
!      mtrsq = mtrsq * alpha_s4Pi**3 /4d0/64d0
!      print *, "color corr. tree",mtrsq
!      STOP
! !---------

!      mtrsq = (zero,zero)
!      call cc_gg_ttgggg(9,C)
!      do i=1,6
!      do j=1,6
!           mtrsq = mtrsq + C(i,j)*real(AM(i,j),dp)
!      enddo
!      enddo
!      call cc_gg_ttgggg(1,C)
!      do i=1,6
!      do j=1,6
!           mtrsq = mtrsq + C(i,j)*real(AM(i,j),dp)
!      enddo
!      enddo
!      call cc_gg_ttgggg(11,C)
!      do i=1,6
!      do j=1,6
!           mtrsq = mtrsq + C(i,j)*real(AM(i,j),dp)
!      enddo
!      enddo
!      call cc_gg_ttgggg(12,C)
!      do i=1,6
!      do j=1,6
!           mtrsq = mtrsq + C(i,j)*real(AM(i,j),dp)
!      enddo
!      enddo
!      mtrsq = mtrsq * alpha_s4Pi**3 /4d0/64d0
!      print *, "color corr. tree",mtrsq

!      STOP


  return
  end subroutine






end  module
