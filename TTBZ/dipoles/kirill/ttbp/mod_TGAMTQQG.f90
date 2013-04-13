! this is the file to compute the matrix element squared
! for qq -> tt + gamma + g
      module Mod_TGAMTQQG
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

      public :: TGAMTQQG


      contains

!---------------------------------------------------------------------------
!     The matrix element squared for amplitude
!     0 -> bar t(p1) + t (p2) + bar q(3) + q(4) + g(5) + gamma(6)
!     We return res(1) for up, res(2) for dn, because of the different
!     electric charges ; the square of the top quark charge is pulled
!     out from the results
!--------------------------------------------------------------------------

      subroutine TGAMTQQG(p,yRnDk,Wgt,res)
      real(dp), intent(out) ::  res(2)
      real(dp), intent(in) :: p(4,6)
      real(dp), intent(in) :: yRnDK(1:8), Wgt
      real(dp), save  :: Cup(10,10),Cdn(10,10)
      type(Particle),save :: ExtParticles(1:6)
      type(TreeProcess),save :: TreeAmpsReal(1:10)
      integer :: Njet, Nhisto,iTree,i,j
      integer :: Nmax(6),i1,i2,i3,i4,i5,i6,hel(6)
      logical, save :: first_time = .true.
      real(dp) :: MomDK(1:4,1:6),PSWgt1,PSWgt2
      complex(dp) :: A(10),Ac(10)
      character :: ixq*2

      res = 0d0



       if( first_time ) then

      call cc_qq_ttgamg('up',Cup)
      call cc_qq_ttgamg('dn',Cdn)


        call InitTrees(4,2,10,TreeAmpsReal)
        call InitProcess_TbTQbQGG(ExtParticles(1:6))  ! 5gluon,6 is a photon
              TreeAmpsReal(1)%PartRef(1:6)= (/1,6,2,3,4,5/)
              TreeAmpsReal(2)%PartRef(1:6)= (/1,6,2,3,5,4/)
              TreeAmpsReal(3)%PartRef(1:6)= (/1,6,2,5,3,4/)
              TreeAmpsReal(4)%PartRef(1:6)= (/1,6,5,2,3,4/)
              TreeAmpsReal(5)%PartRef(1:6)= (/1,2,3,6,4,5/)
              TreeAmpsReal(6)%PartRef(1:6)= (/1,2,3,5,6,4/)
              TreeAmpsReal(7)%PartRef(1:6)= (/1,2,5,3,6,4/)
              TreeAmpsReal(8)%PartRef(1:6)= (/1,2,3,6,5,4/)
              TreeAmpsReal(9)%PartRef(1:6)= (/1,5,6,2,3,4/)
             TreeAmpsReal(10)%PartRef(1:6)= (/1,5,2,3,6,4/)

             do iTree=1,10
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

   call GenerateEventttqqgg((/p(1:4,1),p(1:4,2),p(1:4,3),p(1:4,4), &
     p(1:4,5),p(1:4,6)/), MomDK,hel, ExtParticles(1:6))

       do i=1,10
          call EvalTree2(TreeAmpsReal(i),A(i))
          Ac(i) = dconjg(A(i))
       enddo


          do i=1,10
             do j=1,10
               res(1) = res(1)  + dreal(Cup(i,j)*A(i)*Ac(j),dp)*PSWgt1*PSWgt2
               res(2) = res(2)  + dreal(Cdn(i,j)*A(i)*Ac(j),dp)*PSWgt2*PSWgt2
             enddo
           enddo

       enddo
       enddo
       enddo
       enddo
       enddo
       enddo


      end subroutine

      end module
