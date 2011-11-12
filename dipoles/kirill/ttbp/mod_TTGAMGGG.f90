! module with bar t gamma t  ggg amplitude squared
      module Mod_TGAMTBGGG
      use ModAmplitudes
      use ModProcess
      use ModParameters
      use ModMisc
      use ModKinematics
      use ModIntDipoles2
      use colorcorr_ttgam
      implicit none
      private

      integer, parameter  :: dp = selected_real_kind(15)
      logical, parameter :: invert_alphaCut = .false.

      public :: TGAMTBGGG, TGAMTBGG
      contains


!---------------------------------------------------------------------------
!     The matrix element squared for amplitude
!      0 -> bar t(p1)  + t (p2) + g(3) + g(4) + g(5) + gam(6)
!--------------------------------------------------------------------------

      subroutine TGAMTBGGG(p,yRnDk,Wgt,res)
      real(dp), intent(out) ::  res
      real(dp), intent(in) :: p(4,6)
      real(dp), intent(in) :: yRnDK(1:8),Wgt
      real(dp), save  :: C(6,6)
      type(Particle),save :: ExtParticles(1:6)
      type(TreeProcess),save :: TreeAmpsReal(1:6)
      integer :: Njet,  Nhisto,iTree,i,j
      integer :: Nmax(6),i1,i2,i3,i4,i5,i6,hel(6)
      logical, save :: first_time = .true.
      real(dp) :: MomDK(1:4,1:6),PSWgt1,PSWgt2
      complex(dp) :: A(6),Ac(6)

      res = 0d0


       if( first_time ) then

      call cc_gg_ttgamg(C)  ! the color matrix

             call InitTrees(2,4,6,TreeAmpsReal)
             call InitProcess_TbTGGGG(ExtParticles(1:6))  ! 6 is a photon


             TreeAmpsReal(1)%PartRef(1:6) = (/1,6,2,3,4,5/)
             TreeAmpsReal(2)%PartRef(1:6) = (/1,6,2,3,5,4/)
             TreeAmpsReal(3)%PartRef(1:6) = (/1,6,2,4,3,5/)
             TreeAmpsReal(4)%PartRef(1:6) = (/1,6,2,4,5,3/)
             TreeAmpsReal(5)%PartRef(1:6) = (/1,6,2,5,3,4/)
             TreeAmpsReal(6)%PartRef(1:6) = (/1,6,2,5,4,3/)
             do iTree=1,6
             call LinkTreeParticles(TreeAmpsReal(iTree),ExtParticles(1:6))
             enddo
             first_time=.false.
      endif


       if (TopDecays.ge.1) then
          call EvalPhasespace_TopDecay(p(1:4,1),yRnDk(1:4),.false., &
          MomDK(1:4,1:3),PSWgt1)
          call EvalPhasespace_TopDecay(p(1:4,2),yRnDk(5:8),.false., &
          MomDK(1:4,4:6),PSWgt2)
       elseif(TopDecays.eq.0) then
         MomDK = 0d0
         PSWgt1 = one
         PSWgt2 = one
       endif

      Nmax = 1

        if (TopDecays.ge.1) then ! Top decays
            Nmax(1) =  - 1
            Nmax(2) =  - 1
       endif


       do i1 = -1,Nmax(1),2
          hel(1) = i1
       do i2 = -1,Nmax(2),2
          hel(2) = i2
       do i3 = -1,Nmax(3),2
          hel(3) = i3
       do i4 = -1,Nmax(4),2
          hel(4) = i4
       do i5 = -1,Nmax(5),2
          hel(5) = i5
       do i6 = -1,Nmax(6),2
          hel(6) = i6

   call GenerateEventttgggg((/p(1:4,1),p(1:4,2),p(1:4,3),p(1:4,4), &
     p(1:4,5),p(1:4,6)/), MomDK, hel, ExtParticles(1:6))





       do i=1,6
          call EvalTree2(TreeAmpsReal(i),A(i))
          Ac(i) = conjg(A(i))
        enddo


          do i=1,6
             do j=1,6
               res = res + real(C(i,j)*A(i)*Ac(j),dp)*PSWgt1*PSWgt2
             enddo
           enddo

       enddo
       enddo
       enddo
       enddo
       enddo
       enddo

       end subroutine

!---------------------------------------------------------------------------
!     The matrix element squared for amplitude
!      0 -> bar t(p1) + t (p2) + g(3) + g(4)+gam(5)
!--------------------------------------------------------------------------

      subroutine TGAMTBGG(p,yRnDk,Wgt,res)
      real(dp), intent(out) ::  res
      real(dp), intent(in) :: p(4,5)
      real(dp), intent(in) :: yRnDK(1:8),Wgt
      real(dp), save  :: C(2,2)
      type(Particle),save :: ExtParticles(1:5)
      type(TreeProcess),save :: TreeAmpsReal(1:2)
      integer :: Njet,  Nhisto,iTree,i,j
      integer :: Nmax(5),i1,i2,i3,i4,i5,i6,hel(5)
      logical, save :: first_time = .true.
      real(dp) :: MomDK(1:4,1:6),PSWgt1,PSWgt2
      complex(dp) :: A(2),Ac(2)

      res = 0d0


       if( first_time ) then

      call cc_gg_ttgam(C)  ! the color matrix


             call InitTrees(2,3,2,TreeAmpsReal)
             call InitProcess_TbTGGG(ExtParticles(1:5)) ! 5 is a photon
             TreeAmpsReal(1)%PartRef(1:5) = (/1,5,2,3,4/)
             TreeAmpsReal(2)%PartRef(1:5) = (/1,5,2,4,3/)
             do iTree=1,2
             call LinkTreeParticles(TreeAmpsReal(iTree),ExtParticles(1:5))
             enddo
             first_time=.false.
      endif




       if (TopDecays.ge.1) then
      call EvalPhasespace_TopDecay(p(1:4,1),yRnDk(1:4),.false., &
       MomDK(1:4,1:3),PSWgt1)
     call EvalPhasespace_TopDecay(p(1:4,2),yRnDk(5:8),.false., &
       MomDK(1:4,4:6),PSWgt2)
       elseif(TopDecays.eq.0) then
         MomDK = 0d0
         PSWgt1 = one
         PSWgt2 = one
       endif

      Nmax = 1

        if (TopDecays.ge.1) then ! Top decays
            Nmax(1) =  - 1
            Nmax(2) =  - 1
       endif


       do i1 = -1,Nmax(1),2
          hel(1) = i1
       do i2 = -1,Nmax(2),2
          hel(2) = i2
       do i3 = -1,Nmax(3),2
          hel(3) = i3
       do i4 = -1,Nmax(4),2
          hel(4) = i4
       do i5 = -1,Nmax(5),2
          hel(5) = i5

  call GenerateEvent52((/p(1:4,1),p(1:4,2),p(1:4,3),p(1:4,4),p(1:4,5)/), &
  momDK,hel,ExtParticles(1:5))

       do i=1,2
          call EvalTree2(TreeAmpsReal(i),A(i))
          Ac(i) = conjg(A(i))
!          print *, i, ';', i1, i2, i3, i4, i5, A(i)
       enddo


          do i=1,2
             do j=1,2
               res = res   &
       + real(C(i,j)*A(i)*Ac(j),dp)
             enddo
           enddo

       enddo
       enddo
       enddo
       enddo
       enddo


      end subroutine TGAMTBGG


      function delta1(i1,i2,i3,i4,i5,j1,j2,j3,j4,j5)
        implicit none
        real(dp) :: delta1
        integer, intent(in) :: i1,i2,i3,i4,i5,j1,j2,j3,j4,j5

        delta1 = 1d0


         if (i3.ne.j3.or.i4.ne.j4) then
!         if (abs(i3).ne.1.or.i4.ne.j4) then
        delta1 = 0d0
         endif

      end function delta1


      end module






