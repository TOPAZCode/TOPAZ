module ttqqgg_mod
use types 
use consts_dp
use ModAmplitudes
use ModMisc
implicit none
private 


public :: ttqqgg, momtrans

contains

  subroutine momtrans(pin,i1,i2,j1,j2,pout)
  integer, intent(in) :: i1,i2,j1,j2
  real(dp), intent(in) :: pin(4,6) 
  real(dp), intent(out) :: pout(4,6) 

     pout = pin

     pout(:,i2) = pin(:,i1) 
     pout(:,i1) = pin(:,i2)


     pout(:,j2) = pin(:,j1) 
     pout(:,j1) = pin(:,j2)


  end subroutine momtrans



!------------------------------------------------------------
! the matrix element squared for 0 -> bar t_1 + t_2 + bar q_3 + q_4 
! + g_5 + g_6
!  in the all-outgoing convention 
! The momentum assignement is such that bar t and t are 1 & 2
! initial state particles are 3 and 4 and final state particles are 5 and 6
!------------------------------------------------------------

     subroutine ttqqgg(p,res)
     real(dp),  intent(in) :: p(4,6)
     real(dp), intent(out) :: res 
     real(dp), save  :: C(12,12) ! color matrix 
     logical, save :: first_time = .true.
     integer :: hel(6),i1,i2,i3,i4,i5,i6,i,j, iTree
     type(Particle) :: ExtParticles(1:6)
     type(TreeProcess) :: TreeAmpsReal(1:12)
     complex(dp) :: A(12), Ac(12), cres
!------------------------------------------------------------------

      res = zero
      cres = (zero,zero)


      if (first_time) then

      C(1,1)=64.0_dp 
      C(1,2)=-8.0_dp 
      C(1,3)=8.0_dp 
      C(1,4)=8.0_dp 
      C(1,5)=0.0_dp 
      C(1,6)=0.0_dp 
      C(1,7)=8.0_dp/9.0_dp 
      C(1,8)=-64.0_dp/9.0_dp 
      C(1,9)=8.0_dp/9.0_dp 
      C(1,10)=-64.0_dp/9.0_dp 
      C(1,11)=8.0_dp/9.0_dp 
      C(1,12)=-64.0_dp/9.0_dp 
      C(2,1)=-8.0_dp 
      C(2,2)=64.0_dp 
      C(2,3)=8.0_dp 
      C(2,4)=8.0_dp 
      C(2,5)=0.0_dp 
      C(2,6)=0.0_dp 
      C(2,7)=-64.0_dp/9.0_dp 
      C(2,8)=8.0_dp/9.0_dp 
      C(2,9)=-64.0_dp/9.0_dp 
      C(2,10)=8.0_dp/9.0_dp 
      C(2,11)=-64.0_dp/9.0_dp 
      C(2,12)=8.0_dp/9.0_dp 
      C(3,1)=8.0_dp 
      C(3,2)=8.0_dp 
      C(3,3)=64.0_dp 
      C(3,4)=-8.0_dp 
      C(3,5)=0.0_dp 
      C(3,6)=0.0_dp 
      C(3,7)=8.0_dp/9.0_dp 
      C(3,8)=-64.0_dp/9.0_dp 
      C(3,9)=-64.0_dp/9.0_dp 
      C(3,10)=8.0_dp/9.0_dp 
      C(3,11)=8.0_dp/9.0_dp 
      C(3,12)=-64.0_dp/9.0_dp 
      C(4,1)=8.0_dp 
      C(4,2)=8.0_dp 
      C(4,3)=-8.0_dp 
      C(4,4)=64.0_dp 
      C(4,5)=0.0_dp 
      C(4,6)=0.0_dp 
      C(4,7)=-64.0_dp/9.0_dp 
      C(4,8)=8.0_dp/9.0_dp 
      C(4,9)=8.0_dp/9.0_dp 
      C(4,10)=-64.0_dp/9.0_dp 
      C(4,11)=-64.0_dp/9.0_dp 
      C(4,12)=8.0_dp/9.0_dp 
      C(5,1)=0.0_dp 
      C(5,2)=0.0_dp 
      C(5,3)=0.0_dp 
      C(5,4)=0.0_dp 
      C(5,5)=64.0_dp 
      C(5,6)=8.0_dp 
      C(5,7)=8.0_dp/9.0_dp 
      C(5,8)=-64.0_dp/9.0_dp 
      C(5,9)=-64.0_dp/9.0_dp 
      C(5,10)=-64.0_dp/9.0_dp 
      C(5,11)=-64.0_dp/9.0_dp 
      C(5,12)=8.0_dp/9.0_dp 
      C(6,1)=0.0_dp 
      C(6,2)=0.0_dp 
      C(6,3)=0.0_dp 
      C(6,4)=0.0_dp 
      C(6,5)=8.0_dp 
      C(6,6)=64.0_dp 
      C(6,7)=-64.0_dp/9.0_dp 
      C(6,8)=8.0_dp/9.0_dp 
      C(6,9)=-64.0_dp/9.0_dp 
      C(6,10)=-64.0_dp/9.0_dp 
      C(6,11)=8.0_dp/9.0_dp 
      C(6,12)=-64.0_dp/9.0_dp 
      C(7,1)=8.0_dp/9.0_dp 
      C(7,2)=-64.0_dp/9.0_dp 
      C(7,3)=8.0_dp/9.0_dp 
      C(7,4)=-64.0_dp/9.0_dp 
      C(7,5)=8.0_dp/9.0_dp 
      C(7,6)=-64.0_dp/9.0_dp 
      C(7,7)=64.0_dp/9.0_dp 
      C(7,8)=-8.0_dp/9.0_dp 
      C(7,9)=0.0_dp 
      C(7,10)=0.0_dp 
      C(7,11)=8.0_dp/9.0_dp 
      C(7,12)=8.0_dp/9.0_dp 
      C(8,1)=-64.0_dp/9.0_dp 
      C(8,2)=8.0_dp/9.0_dp 
      C(8,3)=-64.0_dp/9.0_dp 
      C(8,4)=8.0_dp/9.0_dp 
      C(8,5)=-64.0_dp/9.0_dp 
      C(8,6)=8.0_dp/9.0_dp 
      C(8,7)=-8.0_dp/9.0_dp 
      C(8,8)=64.0_dp/9.0_dp 
      C(8,9)=0.0_dp 
      C(8,10)=0.0_dp 
      C(8,11)=8.0_dp/9.0_dp 
      C(8,12)=8.0_dp/9.0_dp 
      C(9,1)=8.0_dp/9.0_dp 
      C(9,2)=-64.0_dp/9.0_dp 
      C(9,3)=-64.0_dp/9.0_dp 
      C(9,4)=8.0_dp/9.0_dp 
      C(9,5)=-64.0_dp/9.0_dp 
      C(9,6)=-64.0_dp/9.0_dp 
      C(9,7)=0.0_dp 
      C(9,8)=0.0_dp 
      C(9,9)=64.0_dp/9.0_dp 
      C(9,10)=8.0_dp/9.0_dp 
      C(9,11)=0.0_dp 
      C(9,12)=0.0_dp 
      C(10,1)=-64.0_dp/9.0_dp 
      C(10,2)=8.0_dp/9.0_dp 
      C(10,3)=8.0_dp/9.0_dp 
      C(10,4)=-64.0_dp/9.0_dp 
      C(10,5)=-64.0_dp/9.0_dp 
      C(10,6)=-64.0_dp/9.0_dp 
      C(10,7)=0.0_dp 
      C(10,8)=0.0_dp 
      C(10,9)=8.0_dp/9.0_dp 
      C(10,10)=64.0_dp/9.0_dp 
      C(10,11)=0.0_dp 
      C(10,12)=0.0_dp 
      C(11,1)=8.0_dp/9.0_dp 
      C(11,2)=-64.0_dp/9.0_dp 
      C(11,3)=8.0_dp/9.0_dp 
      C(11,4)=-64.0_dp/9.0_dp 
      C(11,5)=-64.0_dp/9.0_dp 
      C(11,6)=8.0_dp/9.0_dp 
      C(11,7)=8.0_dp/9.0_dp 
      C(11,8)=8.0_dp/9.0_dp 
      C(11,9)=0.0_dp 
      C(11,10)=0.0_dp 
      C(11,11)=64.0_dp/9.0_dp 
      C(11,12)=-8.0_dp/9.0_dp 
      C(12,1)=-64.0_dp/9.0_dp 
      C(12,2)=8.0_dp/9.0_dp 
      C(12,3)=-64.0_dp/9.0_dp 
      C(12,4)=8.0_dp/9.0_dp 
      C(12,5)=8.0_dp/9.0_dp 
      C(12,6)=-64.0_dp/9.0_dp 
      C(12,7)=8.0_dp/9.0_dp 
      C(12,8)=8.0_dp/9.0_dp 
      C(12,9)=0.0_dp 
      C(12,10)=0.0_dp 
      C(12,11)=-8.0_dp/9.0_dp 
      C(12,12)=64.0_dp/9.0_dp 

       endif 

       if (first_time) first_time = .false.

! init external particles
call InitProcess_TbTQbQGG(ExtParticles(1:6))
! init tree processes for 0-> tb t qb q g g
call InitTrees(4,2,12,TreeAmpsReal)

!--- ordering of particles in a particular primitive amplitude; 
!--- has to be correlated with the form program that calculates 
!--- color correlation matrix


  TreeAmpsReal( 1)%PartRef(1:6) =  (/1, 2, 3, 4, 5, 6/)
  TreeAmpsReal( 2)%PartRef(1:6) =  (/1, 2, 3, 4, 6, 5/)
  TreeAmpsReal( 3)%PartRef(1:6) =  (/1, 2, 5, 6, 3, 4/)
  TreeAmpsReal( 4)%PartRef(1:6) =  (/1, 2, 6, 5, 3, 4/)
  TreeAmpsReal( 5)%PartRef(1:6) =  (/1, 2, 5, 3, 4, 6/)
  TreeAmpsReal( 6)%PartRef(1:6) =  (/1, 2, 6, 3, 4, 5/)
  TreeAmpsReal( 7)%PartRef(1:6) =  (/1, 5, 6, 2, 3, 4/)
  TreeAmpsReal( 8)%PartRef(1:6) =  (/1, 6, 5, 2, 3, 4/)
  TreeAmpsReal( 9)%PartRef(1:6) =  (/1, 5, 2, 3, 6, 4/)
  TreeAmpsReal( 10)%PartRef(1:6) = (/1, 6, 2, 3, 5, 4/)
  TreeAmpsReal( 11)%PartRef(1:6) = (/1, 2, 3, 5, 6, 4/)
  TreeAmpsReal( 12)%PartRef(1:6) = (/1, 2, 3, 6, 5, 4/)

  do iTree=1,12
      call LinkTreeParticles(TreeAmpsReal(iTree),ExtParticles(1:6))
  enddo

         do i1 =-1,1,2
             hel(1) = i1
        do i2 = -1,1,2
            hel(2) = i2
        do i3 = -1,1,2
            hel(3) = i3
        do i4 = -1,1,2
            hel(4) = i4
        do i5 = -1,1,2
            hel(5) = i5
        do i6 = -1,1,2
            hel(6) = i6

             call GenerateEventttqqgg(p,hel,ExtParticles(1:6))

             do i = 1,12
      call EvalTree(TreeAmpsReal(i),A(i))
        Ac(i) = conjg(A(i))
             enddo 


       do i=1,12
          do j=1,12
           cres = cres + A(i)*Ac(j)*C(i,j)
       enddo 
       enddo



!       do i=4,4
!           cres = cres + A(i)*Ac(i)*C(i,i)
!           cres = A(i)*Ac(j)*C(i,j)
!           if (abs(cres).gt.20.0_dp) then 
!           print *, i, j, cres
!           endif
!       enddo


       

      enddo 
      enddo 
      enddo 
      enddo 
      enddo 
      enddo 
      
      res = real(cres,dp)

      end subroutine ttqqgg




      end module ttqqgg_mod


