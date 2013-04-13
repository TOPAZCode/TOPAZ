module modIntDipoles_QQBTTBGG
use ModParameters
use ModProcess
use ModMisc
use ModIntDipoles
implicit none

public:: EvalIntDipoles_QQBTTBGG
integer, parameter  :: dp = selected_real_kind(15)

contains


      subroutine EvalIntDipoles_QQBTTBGG(p,MomDK,z,res)
      use ModAmplitudes
      implicit none
      real(dp), intent(in) :: p(4,5),MomDK(1:4,1:6),z
      real(dp), intent(out) :: res(1:3)
      type(Particle),save :: ExtParticles(1:5)
      type(TreeProcess),save :: TreeAmpsDip(1:4)
      complex(dp) :: Bm(1:4)
      complex(dp) :: AM(1:4,1:4)
      integer :: i1,i2,i3,i4,hel(5),i,j,i5,i6,emi,n, in1, in2,iTree
      integer ::  Nmax(5)
      real(dp) :: dipsoft, dipfini, dipplus
      real(dp) :: C(4,4)
      real(dp) :: mtrs,epcorr,AP(1:3),CF
      character :: diptype*2
      logical, save :: first_time = .true.

      CF=4d0/3d0
      res(1:3) = zero



        if( first_time ) then
! init external particles
              call InitTrees(4,1,4,TreeAmpsDip)
              call InitProcess_TbTQbQG(ExtParticles(1:5))
              TreeAmpsDip(1)%PartRef(1:5) = (/1,2,3,4,5/)
              TreeAmpsDip(2)%PartRef(1:5) = (/1,5,2,3,4/)
              TreeAmpsDip(3)%PartRef(1:5) = (/1,2,5,3,4/)
              TreeAmpsDip(4)%PartRef(1:5) = (/1,2,3,5,4/)
              do iTree=1,4
                  call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles(1:5))
              enddo
              first_time=.false.
      endif


      AM = (zero,zero) ! first, the ``amplitude matrix'' is set to zero


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

      call GenerateEventttqqg2(p,momDK(1:4,1:6),hel,ExtParticles(1:5))

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


   do n=1,36

        call cc_qq_ttqqgg(n,C)

        mtrs = zero
      do i=1,4
      do j=1,4
          mtrs = mtrs + C(i,j)*real(AM(i,j),dp)
      enddo
      enddo
      mtrs = mtrs * alpha_s4Pi**3 /4d0/9d0   ! spin & color avg.


!       if( n.eq.9 .or. n.eq.13 .or. n.eq.29 .or. n.eq.33 ) then
!       if( n.eq.12 .or. n.eq.16 .or. n.eq.32 .or. n.eq.36 ) then
!       if( n.eq.17 .or. n.eq.18 ) then
!           alpha_ff=alpha_ff
!       else
!           alpha_tmp=alpha_ff
!           alpha_ff= 1d0
!       endif


      if(n.eq.1) then
      dipsoft =if_qq(zero,m_Top,p,3,1,z,1)
      dipfini =if_qq(zero,m_Top,p,3,1,z,2)
      dipplus =if_qq(zero,m_Top,p,3,1,z,3)
      diptype ='if'
     emi = 1 ! mom #3 is emitting
      endif
      if(n.eq.2) then
      dipsoft =if_qq(zero,m_Top,p,3,2,z,1)
      dipfini =if_qq(zero,m_Top,p,3,2,z,2)
      dipplus =if_qq(zero,m_Top,p,3,2,z,3)
      diptype ='if'
     emi = 1 ! mom #3 is emitting
      endif
      if(n.eq.3) then
      dipsoft =ii_qq(zero,zero,p,3,4,z,1)
      dipfini =ii_qq(zero,zero,p,3,4,z,2)
      dipplus =ii_qq(zero,zero,p,3,4,z,3)
      diptype ='ii'
     emi = 1 ! mom #3 is emitting
      endif
      if(n.eq.4) then
      dipsoft =if_qq(zero,zero,p,3,5,z,1)
      dipfini =if_qq(zero,zero,p,3,5,z,2)
      dipplus =if_qq(zero,zero,p,3,5,z,3)
      diptype ='if'
     emi = 1 ! mom #3 is emitting
      endif
      if(n.eq.5) then
      dipsoft =if_qq(zero,m_Top,p,4,1,z,1)
      dipfini =if_qq(zero,m_Top,p,4,1,z,2)
      dipplus =if_qq(zero,m_Top,p,4,1,z,3)
      diptype ='if'
     emi = 2 ! mom #4 is emitting
      endif
      if(n.eq.6) then
      dipsoft =if_qq(zero,m_Top,p,4,2,z,1)
      dipfini =if_qq(zero,m_Top,p,4,2,z,2)
      dipplus =if_qq(zero,m_Top,p,4,2,z,3)
      diptype ='if'
     emi = 2 ! mom #4 is emitting
      endif
      if(n.eq.7) then
      dipsoft =ii_qq(zero,zero,p,4,3,z,1)
      dipfini =ii_qq(zero,zero,p,4,3,z,2)
      dipplus =ii_qq(zero,zero,p,4,3,z,3)
      diptype ='ii'
     emi = 2 ! mom #4 is emitting
      endif
      if(n.eq.8) then
      dipsoft =if_qq(zero,zero,p,4,5,z,1)
      dipfini =if_qq(zero,zero,p,4,5,z,2)
      dipplus =if_qq(zero,zero,p,4,5,z,3)
      diptype ='if'
     emi = 2 ! mom #4 is emitting
      endif
      if(n.eq.9) then
      dipsoft =ff_qq(m_Top,m_Top,p,1,2,z,1)
      dipfini =ff_qq(m_Top,m_Top,p,1,2,z,2)
      dipplus =ff_qq(m_Top,m_Top,p,1,2,z,3)
      diptype ='ff'
      emi = 3 ! final state dipole
      endif
      if(n.eq.10) then
      dipsoft =fi_qq(m_Top,zero,p,1,3,z,1)
      dipfini =fi_qq(m_Top,zero,p,1,3,z,2)
      dipplus =fi_qq(m_Top,zero,p,1,3,z,3)
      diptype ='fi'
      emi = 1 ! final state dipole
      endif
      if(n.eq.11) then
      dipsoft =fi_qq(m_Top,zero,p,1,4,z,1)
      dipfini =fi_qq(m_Top,zero,p,1,4,z,2)
      dipplus =fi_qq(m_Top,zero,p,1,4,z,3)
      diptype ='fi'
      emi = 2 ! final state dipole
      endif
      if(n.eq.12) then
      dipsoft =ff_qq(m_Top,zero,p,1,5,z,1)
      dipfini =ff_qq(m_Top,zero,p,1,5,z,2)
      dipplus =ff_qq(m_Top,zero,p,1,5,z,3)
      diptype ='ff'
      emi = 3 ! final state dipole
      endif
      if(n.eq.13) then
      dipsoft =ff_qq(m_Top,m_Top,p,2,1,z,1)
      dipfini =ff_qq(m_Top,m_Top,p,2,1,z,2)
      dipplus =ff_qq(m_Top,m_Top,p,2,1,z,3)
      diptype ='ff'
      emi = 3 ! final state dipole
      endif
      if(n.eq.14) then
      dipsoft =fi_qq(m_Top,zero,p,2,3,z,1)
      dipfini =fi_qq(m_Top,zero,p,2,3,z,2)
      dipplus =fi_qq(m_Top,zero,p,2,3,z,3)
      diptype ='fi'
      emi = 1 ! final state dipole
      endif
      if(n.eq.15) then
      dipsoft =fi_qq(m_Top,zero,p,2,4,z,1)
      dipfini =fi_qq(m_Top,zero,p,2,4,z,2)
      dipplus =fi_qq(m_Top,zero,p,2,4,z,3)
      diptype ='fi'
      emi = 2 ! final state dipole
      endif
      if(n.eq.16) then
      dipsoft =ff_qq(m_Top,zero,p,2,5,z,1)
      dipfini =ff_qq(m_Top,zero,p,2,5,z,2)
      dipplus =ff_qq(m_Top,zero,p,2,5,z,3)
      diptype ='ff'
      emi = 3 ! final state dipole
      endif
      if(n.eq.17) then
      dipsoft =ff_gg(zero,m_Top,p,5,1,z,1)
      dipfini =ff_gg(zero,m_Top,p,5,1,z,2)
      dipplus =ff_gg(zero,m_Top,p,5,1,z,3)
      diptype ='ff'
      emi = 3 ! final state dipole
      endif
      if(n.eq.18) then
      dipsoft =ff_gg(zero,m_Top,p,5,2,z,1)
      dipfini =ff_gg(zero,m_Top,p,5,2,z,2)
      dipplus =ff_gg(zero,m_Top,p,5,2,z,3)
      diptype ='ff'
      emi = 3 ! final state dipole
      endif
      if(n.eq.19) then
      dipsoft =fi_gg(zero,zero,p,5,3,z,1)
      dipfini =fi_gg(zero,zero,p,5,3,z,2)
      dipplus =fi_gg(zero,zero,p,5,3,z,3)
      diptype ='fi'
      emi = 1 ! final state dipole
      endif
      if(n.eq.20) then
      dipsoft =fi_gg(zero,zero,p,5,4,z,1)
      dipfini =fi_gg(zero,zero,p,5,4,z,2)
      dipplus =fi_gg(zero,zero,p,5,4,z,3)
      diptype ='fi'
      emi = 2 ! final state dipole
      endif
      if(n.eq.21) then
      dipsoft =if_qq(zero,m_Top,p,3,1,z,1)
      dipfini =if_qq(zero,m_Top,p,3,1,z,2)
      dipplus =if_qq(zero,m_Top,p,3,1,z,3)
      diptype ='if'
     emi = 1 ! mom #3 is emitting
      endif
      if(n.eq.22) then
      dipsoft =if_qq(zero,m_Top,p,3,2,z,1)
      dipfini =if_qq(zero,m_Top,p,3,2,z,2)
      dipplus =if_qq(zero,m_Top,p,3,2,z,3)
      diptype ='if'
     emi = 1 ! mom #3 is emitting
      endif
      if(n.eq.23) then
      dipsoft =ii_qq(zero,zero,p,3,4,z,1)
      dipfini =ii_qq(zero,zero,p,3,4,z,2)
      dipplus =ii_qq(zero,zero,p,3,4,z,3)
      diptype ='ii'
     emi = 1 ! mom #3 is emitting
      endif
      if(n.eq.24) then
      dipsoft =if_qq(zero,zero,p,3,5,z,1)
      dipfini =if_qq(zero,zero,p,3,5,z,2)
      dipplus =if_qq(zero,zero,p,3,5,z,3)
      diptype ='if'
     emi = 1 ! mom #3 is emitting
      endif
      if(n.eq.25) then
      dipsoft =if_qq(zero,m_Top,p,4,1,z,1)
      dipfini =if_qq(zero,m_Top,p,4,1,z,2)
      dipplus =if_qq(zero,m_Top,p,4,1,z,3)
      diptype ='if'
     emi = 2 ! mom #4 is emitting
      endif
      if(n.eq.26) then
      dipsoft =if_qq(zero,m_Top,p,4,2,z,1)
      dipfini =if_qq(zero,m_Top,p,4,2,z,2)
      dipplus =if_qq(zero,m_Top,p,4,2,z,3)
      diptype ='if'
     emi = 2 ! mom #4 is emitting
      endif
      if(n.eq.27) then
      dipsoft =ii_qq(zero,zero,p,4,3,z,1)
      dipfini =ii_qq(zero,zero,p,4,3,z,2)
      dipplus =ii_qq(zero,zero,p,4,3,z,3)
      diptype ='ii'
     emi = 2 ! mom #4 is emitting
      endif
      if(n.eq.28) then
      dipsoft =if_qq(zero,zero,p,4,5,z,1)
      dipfini =if_qq(zero,zero,p,4,5,z,2)
      dipplus =if_qq(zero,zero,p,4,5,z,3)
      diptype ='if'
     emi = 2 ! mom #4 is emitting
      endif
      if(n.eq.29) then
      dipsoft =ff_qq(m_Top,m_Top,p,1,2,z,1)
      dipfini =ff_qq(m_Top,m_Top,p,1,2,z,2)
      dipplus =ff_qq(m_Top,m_Top,p,1,2,z,3)
      diptype ='ff'
      emi = 3 ! final state dipole
      endif
      if(n.eq.30) then
      dipsoft =fi_qq(m_Top,zero,p,1,3,z,1)
      dipfini =fi_qq(m_Top,zero,p,1,3,z,2)
      dipplus =fi_qq(m_Top,zero,p,1,3,z,3)
      diptype ='fi'
      emi = 1 ! final state dipole
      endif
      if(n.eq.31) then
      dipsoft =fi_qq(m_Top,zero,p,1,4,z,1)
      dipfini =fi_qq(m_Top,zero,p,1,4,z,2)
      dipplus =fi_qq(m_Top,zero,p,1,4,z,3)
      diptype ='fi'
      emi = 2 ! final state dipole
      endif
      if(n.eq.32) then
      dipsoft =ff_qq(m_Top,zero,p,1,5,z,1)
      dipfini =ff_qq(m_Top,zero,p,1,5,z,2)
      dipplus =ff_qq(m_Top,zero,p,1,5,z,3)
      diptype ='ff'
      emi = 3 ! final state dipole
      endif
      if(n.eq.33) then
      dipsoft =ff_qq(m_Top,m_Top,p,2,1,z,1)
      dipfini =ff_qq(m_Top,m_Top,p,2,1,z,2)
      dipplus =ff_qq(m_Top,m_Top,p,2,1,z,3)
      diptype ='ff'
      emi = 3 ! final state dipole
      endif
      if(n.eq.34) then
      dipsoft =fi_qq(m_Top,zero,p,2,3,z,1)
      dipfini =fi_qq(m_Top,zero,p,2,3,z,2)
      dipplus =fi_qq(m_Top,zero,p,2,3,z,3)
      diptype ='fi'
      emi = 1 ! final state dipole
      endif
      if(n.eq.35) then
      dipsoft =fi_qq(m_Top,zero,p,2,4,z,1)
      dipfini =fi_qq(m_Top,zero,p,2,4,z,2)
      dipplus =fi_qq(m_Top,zero,p,2,4,z,3)
      diptype ='fi'
      emi = 2 ! final state dipole
      endif
      if(n.eq.36) then
      dipsoft =ff_qq(m_Top,zero,p,2,5,z,1)
      dipfini =ff_qq(m_Top,zero,p,2,5,z,2)
      dipplus =ff_qq(m_Top,zero,p,2,5,z,3)
      diptype ='ff'
      emi = 3 ! final state dipole
      endif

!       if( n.eq.9 .or. n.eq.13 .or. n.eq.29 .or. n.eq.33 ) then
!       if( n.eq.12 .or. n.eq.16 .or. n.eq.32 .or. n.eq.36 ) then
!       if( n.eq.17 .or. n.eq.18 ) then
!           alpha_ff=alpha_ff
!       else
!           alpha_ff= alpha_tmp
!       endif

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
!      RES(1) = RES(1) +  DIPSOFT*MTRS  ! FOR DELTA-FCT. CHECK
!         print *, n,z,dipfini,dipplus,dipsoft
!         pause
   enddo
   res(1:3) = alpha_sOver2Pi*0.5d0 *0.5d0 * res(1:3)  ! overall normalization: alpha_s/4pi


      call cc_ttqqg(C)
      mtrs = zero
      do i=1,4
      do j=1,4
          mtrs = mtrs + C(i,j)*real(AM(i,j),dp)
      enddo
      enddo
      mtrs = mtrs * alpha_s4Pi**3 /4d0/9d0   ! spin & color avg.  ! check this against MG!!!


! !        epcorr=epinv+2d0*dlog(renscale/facscale)
       epcorr=epinv
       AP(1)= 3d0/2d0*CF
       AP(2)= (-1d0-z)*CF
       AP(3)= 2d0*CF/(1d0-z)
       AP(1:3) = AP(1:3) * alpha_sOver2Pi *epcorr * mtrs

        res(1) = res(1) + 2d0*(AP(1)-AP(3))
        res(2) = res(2) + (AP(2) + AP(3))
        res(3) = res(3) + (AP(2) + AP(3))
!      res(1) = res(1) + 2d0*AP(1)  ! for delta-fct. check

!    print *, "eval only delta part"

  end subroutine

end  module


