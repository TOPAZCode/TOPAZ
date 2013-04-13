! color correlation file 
      module colorcorr_ttgam
      implicit none
      private 

      public :: cc_gg_ttgamg, cc_gg_ttgamg_soft, &
       cc_qq_ttgamg_soft, cc_qq_ttgamg, cc_gg_ttgam, & 
       cc_qq_ttgam


      integer, parameter  :: dp = selected_real_kind(15)

      contains


      subroutine cc_qq_ttgam(ixq,C) 
      character, intent(in) :: ixq*2 
      real(dp), intent(out) :: C(2,2) 

      C = 0.0_dp

      if (ixq.eq.'up') then 

       C(1,1)  =  8.0_dp
       C(1,2)  =  8.0_dp
       C(2,1)  =  8.0_dp
       C(2,2)  =  8.0_dp

       elseif(ixq.eq.'dn') then 

       C(1,1)  =  8.0_dp
       C(1,2)  =  -4.0_dp
       C(2,1)  =  -4.0_dp
       C(2,2)  =   2.0_dp

       endif 

      end subroutine 




      subroutine cc_qq_ttgamg(ixq,C) 
      character, intent(in) :: ixq*2 
      real(dp), intent(out) :: C(10,10) 

      C = 0.0_dp

      if (ixq.eq.'up') then 

       C(1,1) = 24.0_dp
       C(1,2) = 8.0_dp/3.0_dp 
       C(1,4) = 8.0_dp/3.0_dp 
       C(1,5) = 24.0_dp 
       C(1,6) = 8.0_dp/3.0_dp 
       C(1,8) = 8.0_dp/3.0_dp 
       C(1,9) = 8.0_dp/3.0_dp 
       C(1,10) = 8.0_dp/3.0_dp 
       C(2,1) = 8.0_dp/3.0_dp 
       C(2,2) = 8.0_dp/3.0_dp 
       C(2,3) = 8.0_dp/3.0_dp 
       C(2,5) = 8.0_dp/3.0_dp 
       C(2,6) = 8.0_dp/3.0_dp 
       C(2,7) = 8.0_dp/3.0_dp 
       C(2,8) = 8.0_dp/3.0_dp 
       C(3,2) = 8.0_dp/3.0_dp 
       C(3,3) = 24.0_dp 
       C(3,4) = 8.0_dp/3.0_dp 
       C(3,6) = 8.0_dp/3.0_dp 
       C(3,7) = 24.0_dp 
       C(3,8) = 8.0_dp/3.0_dp 
       C(3,9) = 8.0_dp/3.0_dp 
       C(3,10) = 8.0_dp/3.0_dp 
       C(4,1) = 8.0_dp/3.0_dp 
       C(4,3) = 8.0_dp/3.0_dp 
       C(4,4) = 8.0_dp/3.0_dp 
       C(4,5) = 8.0_dp/3.0_dp 
       C(4,7) = 8.0_dp/3.0_dp 
       C(4,9) = 8.0_dp/3.0_dp 
       C(4,10) = 8.0_dp/3.0_dp 
       C(5,1) = 24.0_dp 
       C(5,2) = 8.0_dp/3.0_dp 
       C(5,4) = 8.0_dp/3.0_dp 
       C(5,5) = 24.0_dp 
       C(5,6) = 8.0_dp/3.0_dp 
       C(5,8) = 8.0_dp/3.0_dp 
       C(5,9) = 8.0_dp/3.0_dp 
       C(5,10) = 8.0_dp/3.0_dp 
       C(6,1) = 8.0_dp/3.0_dp 
       C(6,2) = 8.0_dp/3.0_dp 
       C(6,3) = 8.0_dp/3.0_dp 
       C(6,5) = 8.0_dp/3.0_dp 
       C(6,6) = 8.0_dp/3.0_dp 
       C(6,7) = 8.0_dp/3.0_dp 
       C(6,8) = 8.0_dp/3.0_dp 
       C(7,2) = 8.0_dp/3.0_dp 
       C(7,3) = 24.0_dp 
       C(7,4) = 8.0_dp/3.0_dp 
       C(7,6) = 8.0_dp/3.0_dp 
       C(7,7) = 24.0_dp 
       C(7,8) = 8.0_dp/3.0_dp 
       C(7,9) = 8.0_dp/3.0_dp 
       C(7,10) = 8.0_dp/3.0_dp 
       C(8,1) = 8.0_dp/3.0_dp 
       C(8,2) = 8.0_dp/3.0_dp 
       C(8,3) = 8.0_dp/3.0_dp 
       C(8,5) = 8.0_dp/3.0_dp 
       C(8,6) = 8.0_dp/3.0_dp 
       C(8,7) = 8.0_dp/3.0_dp 
       C(8,8) = 8.0_dp/3.0_dp 
       C(9,1) = 8.0_dp/3.0_dp 
       C(9,3) = 8.0_dp/3.0_dp 
       C(9,4) = 8.0_dp/3.0_dp 
       C(9,5) = 8.0_dp/3.0_dp 
       C(9,7) = 8.0_dp/3.0_dp 
       C(9,9) = 8.0_dp/3.0_dp 
       C(9,10) = 8.0_dp/3.0_dp 
       C(10,1) = 8.0_dp/3.0_dp 
       C(10,3) = 8.0_dp/3.0_dp 
       C(10,4) = 8.0_dp/3.0_dp 
       C(10,5) = 8.0_dp/3.0_dp 
       C(10,7) = 8.0_dp/3.0_dp 
       C(10,9) = 8.0_dp/3.0_dp 
       C(10,10) = 8.0_dp/3.0_dp 

       elseif(ixq.eq.'dn') then 


        C(1,1) = ( 24.0_dp )
        C(1,2) = ( 8.0_dp/3.0_dp )
        C(1,4) = ( 8.0_dp/3.0_dp )
        C(1,5) = (  - 12.0_dp )
        C(1,6) = (  - 4.0_dp/3.0_dp )
        C(1,8) = (  - 4.0_dp/3.0_dp )
        C(1,9) = ( 8.0_dp/3.0_dp )
        C(1,10) = (  - 4.0_dp/3.0_dp )
        C(2,1) = ( 8.0_dp/3.0_dp )
        C(2,2) = ( 8.0_dp/3.0_dp )
        C(2,3) = ( 8.0_dp/3.0_dp )
        C(2,5) = (  - 4.0_dp/3.0_dp )
        C(2,6) = (  - 4.0_dp/3.0_dp )
        C(2,7) = (  - 4.0_dp/3.0_dp )
        C(2,8) = (  - 4.0_dp/3.0_dp )
        C(3,2) = ( 8.0_dp/3.0_dp )
        C(3,3) = ( 24.0_dp )
        C(3,4) = ( 8.0_dp/3.0_dp )
        C(3,6) = (  - 4.0_dp/3.0_dp )
        C(3,7) = (  - 12.0_dp )
        C(3,8) = (  - 4.0_dp/3.0_dp )
        C(3,9) = ( 8.0_dp/3.0_dp )
        C(3,10) = (  - 4.0_dp/3.0_dp )
        C(4,1) = ( 8.0_dp/3.0_dp )
        C(4,3) = ( 8.0_dp/3.0_dp )
        C(4,4) = ( 8.0_dp/3.0_dp )
        C(4,5) = (  - 4.0_dp/3.0_dp )
        C(4,7) = (  - 4.0_dp/3.0_dp )
        C(4,9) = ( 8.0_dp/3.0_dp )
        C(4,10) = (  - 4.0_dp/3.0_dp )
        C(5,1) = (  - 12.0_dp )
        C(5,2) = (  - 4.0_dp/3.0_dp )
        C(5,4) = (  - 4.0_dp/3.0_dp )
        C(5,5) = ( 6 )
        C(5,6) = ( 2.0_dp/3.0_dp )
        C(5,8) = ( 2.0_dp/3.0_dp )
        C(5,9) = (  - 4.0_dp/3.0_dp )
        C(5,10) = ( 2.0_dp/3.0_dp )
        C(6,1) = (  - 4.0_dp/3.0_dp )
        C(6,2) = (  - 4.0_dp/3.0_dp )
        C(6,3) = (  - 4.0_dp/3.0_dp )
        C(6,5) = ( 2.0_dp/3.0_dp )
        C(6,6) = ( 2.0_dp/3.0_dp )
        C(6,7) = ( 2.0_dp/3.0_dp )
        C(6,8) = ( 2.0_dp/3.0_dp )
        C(7,2) = (  - 4.0_dp/3.0_dp )
        C(7,3) = (  - 12.0_dp )
        C(7,4) = (  - 4.0_dp/3.0_dp )
        C(7,6) = ( 2.0_dp/3.0_dp )
        C(7,7) = ( 6 )
        C(7,8) = ( 2.0_dp/3.0_dp )
        C(7,9) = (  - 4.0_dp/3.0_dp )
        C(7,10) = ( 2.0_dp/3.0_dp )
        C(8,1) = (  - 4.0_dp/3.0_dp )
        C(8,2) = (  - 4.0_dp/3.0_dp )
        C(8,3) = (  - 4.0_dp/3.0_dp )
        C(8,5) = ( 2.0_dp/3.0_dp )
        C(8,6) = ( 2.0_dp/3.0_dp )
        C(8,7) = ( 2.0_dp/3.0_dp )
        C(8,8) = ( 2.0_dp/3.0_dp )
        C(9,1) = ( 8.0_dp/3.0_dp )
        C(9,3) = ( 8.0_dp/3.0_dp )
        C(9,4) = ( 8.0_dp/3.0_dp )
        C(9,5) = (  - 4.0_dp/3.0_dp )
        C(9,7) = (  - 4.0_dp/3.0_dp )
        C(9,9) = ( 8.0_dp/3.0_dp )
        C(9,10) = (  - 4.0_dp/3.0_dp )
        C(10,1) = (  - 4.0_dp/3.0_dp )
        C(10,3) = (  - 4.0_dp/3.0_dp )
        C(10,4) = (  - 4.0_dp/3.0_dp )
        C(10,5) = ( 2.0_dp/3.0_dp )
        C(10,7) = ( 2.0_dp/3.0_dp )
        C(10,9) = (  - 4.0_dp/3.0_dp )
        C(10,10) = ( 2.0_dp/3.0_dp )

       endif 

      end subroutine 


      subroutine cc_qq_ttgamg_soft(n,ixq,C)  ! ixq labels the power of the 
      integer, intent(in) :: n           ! electric charge of the 
      real(dp), intent(out) :: C(2,2)        ! emitting quark
      character, intent(in) :: ixq*2 

      C = 0.0_dp

      if (ixq.eq.'up') then 



      if(n.eq.1) then  
      C(1,1)=-8.0_dp/3.0_dp 
      C(1,2)=-8.0_dp/3.0_dp 
      C(2,1)=-8.0_dp/3.0_dp 
      C(2,2)=-8.0_dp/3.0_dp 
      endif 
      if(n.eq.2) then  
      C(1,1)=16.0_dp/3.0_dp 
      C(1,2)=16.0_dp/3.0_dp 
      C(2,1)=16.0_dp/3.0_dp 
      C(2,2)=16.0_dp/3.0_dp 
      endif 
      if(n.eq.3) then  
      C(1,1)=56.0_dp/3.0_dp 
      C(1,2)=56.0_dp/3.0_dp 
      C(2,1)=56.0_dp/3.0_dp 
      C(2,2)=56.0_dp/3.0_dp 
      endif 
      if(n.eq.4) then  
      C(1,1)=-8.0_dp/3.0_dp 
      C(1,2)=-8.0_dp/3.0_dp 
      C(2,1)=-8.0_dp/3.0_dp 
      C(2,2)=-8.0_dp/3.0_dp 
      endif 
      if(n.eq.5) then  
      C(1,1)=56.0_dp/3.0_dp 
      C(1,2)=56.0_dp/3.0_dp 
      C(2,1)=56.0_dp/3.0_dp 
      C(2,2)=56.0_dp/3.0_dp 
      endif 
      if(n.eq.6) then  
      C(1,1)=16.0_dp/3.0_dp 
      C(1,2)=16.0_dp/3.0_dp 
      C(2,1)=16.0_dp/3.0_dp 
      C(2,2)=16.0_dp/3.0_dp 
      endif 
      if(n.eq.7) then  
      C(1,1)=16.0_dp/3.0_dp 
      C(1,2)=16.0_dp/3.0_dp 
      C(2,1)=16.0_dp/3.0_dp 
      C(2,2)=16.0_dp/3.0_dp 
      endif 
      if(n.eq.8) then  
      C(1,1)=56.0_dp/3.0_dp 
      C(1,2)=56.0_dp/3.0_dp 
      C(2,1)=56.0_dp/3.0_dp 
      C(2,2)=56.0_dp/3.0_dp 
      endif 
      if(n.eq.9) then  
      C(1,1)=-8.0_dp/3.0_dp 
      C(1,2)=-8.0_dp/3.0_dp 
      C(2,1)=-8.0_dp/3.0_dp 
      C(2,2)=-8.0_dp/3.0_dp 
      endif 
      if(n.eq.10) then  
      C(1,1)=56.0_dp/3.0_dp 
      C(1,2)=56.0_dp/3.0_dp 
      C(2,1)=56.0_dp/3.0_dp 
      C(2,2)=56.0_dp/3.0_dp 
      endif 
      if(n.eq.11) then  
      C(1,1)=16.0_dp/3.0_dp 
      C(1,2)=16.0_dp/3.0_dp 
      C(2,1)=16.0_dp/3.0_dp 
      C(2,2)=16.0_dp/3.0_dp 
      endif 
      if(n.eq.12) then  
      C(1,1)=-8.0_dp/3.0_dp 
      C(1,2)=-8.0_dp/3.0_dp 
      C(2,1)=-8.0_dp/3.0_dp 
      C(2,2)=-8.0_dp/3.0_dp 
      endif 



      elseif(ixq.eq.'dn') then 

      if(n.eq.1) then  
      C(1,1)=-8.0_dp/3.0_dp 
      C(1,2)=4.0_dp/3.0_dp 
      C(2,1)=4.0_dp/3.0_dp 
      C(2,2)=-2.0_dp/3.0_dp 
      endif 
      if(n.eq.2) then  
      C(1,1)=16.0_dp/3.0_dp 
      C(1,2)=-8.0_dp/3.0_dp 
      C(2,1)=-8.0_dp/3.0_dp 
      C(2,2)=4.0_dp/3.0_dp 
      endif 
      if(n.eq.3) then  
      C(1,1)=56.0_dp/3.0_dp 
      C(1,2)=-28.0_dp/3.0_dp 
      C(2,1)=-28.0_dp/3.0_dp 
      C(2,2)=14.0_dp/3.0_dp 
      endif 
      if(n.eq.4) then  
      C(1,1)=-8.0_dp/3.0_dp 
      C(1,2)=4.0_dp/3.0_dp 
      C(2,1)=4.0_dp/3.0_dp 
      C(2,2)=-2.0_dp/3.0_dp 
      endif 
      if(n.eq.5) then  
      C(1,1)=56.0_dp/3.0_dp 
      C(1,2)=-28.0_dp/3.0_dp 
      C(2,1)=-28.0_dp/3.0_dp 
      C(2,2)=14.0_dp/3.0_dp 
      endif 
      if(n.eq.6) then  
      C(1,1)=16.0_dp/3.0_dp 
      C(1,2)=-8.0_dp/3.0_dp 
      C(2,1)=-8.0_dp/3.0_dp 
      C(2,2)=4.0_dp/3.0_dp 
      endif 
      if(n.eq.7) then  
      C(1,1)=16.0_dp/3.0_dp 
      C(1,2)=-8.0_dp/3.0_dp 
      C(2,1)=-8.0_dp/3.0_dp 
      C(2,2)=4.0_dp/3.0_dp 
      endif 
      if(n.eq.8) then  
      C(1,1)=56.0_dp/3.0_dp 
      C(1,2)=-28.0_dp/3.0_dp 
      C(2,1)=-28.0_dp/3.0_dp 
      C(2,2)=14.0_dp/3.0_dp 
      endif 
      if(n.eq.9) then  
      C(1,1)=-8.0_dp/3.0_dp 
      C(1,2)=4.0_dp/3.0_dp 
      C(2,1)=4.0_dp/3.0_dp 
      C(2,2)=-2.0_dp/3.0_dp 
      endif 
      if(n.eq.10) then  
      C(1,1)=56.0_dp/3.0_dp 
      C(1,2)=-28.0_dp/3.0_dp 
      C(2,1)=-28.0_dp/3.0_dp 
      C(2,2)=14.0_dp/3.0_dp 
      endif 
      if(n.eq.11) then  
      C(1,1)=16.0_dp/3.0_dp 
      C(1,2)=-8.0_dp/3.0_dp 
      C(2,1)=-8.0_dp/3.0_dp 
      C(2,2)=4.0_dp/3.0_dp 
      endif 
      if(n.eq.12) then  
      C(1,1)=-8.0_dp/3.0_dp 
      C(1,2)=4.0_dp/3.0_dp 
      C(2,1)=4.0_dp/3.0_dp 
      C(2,2)=-2.0_dp/3.0_dp 
      endif 

      endif   
      

      end subroutine 


      subroutine cc_gg_ttgamg(C)   
      real(dp), intent(out) :: C(6,6) 


       C(1,1) = 512.0_dp/9.0_dp 
       C(1,2) =  - 64.0_dp/9 
       C(1,3) =  - 64.0_dp/9.0_dp 
       C(1,4) = 8.0_dp/9.0_dp 
       C(1,5) = 8.0_dp/9.0_dp 
       C(1,6) = 80.0_dp/9.0_dp
       C(2,1) =  - 64.0_dp/9.0_dp 
       C(2,2) = 512.0_dp/9.0_dp 
       C(2,3) = 8.0_dp/9.0_dp 
       C(2,4) = 80.0_dp/9.0_dp 
       C(2,5) =  - 64.0_dp/9.0_dp 
       C(2,6) = 8.0_dp/9.0_dp 
       C(3,1) =  - 64.0_dp/9.0_dp 
       C(3,2) = 8.0_dp/9.0_dp 
       C(3,3) = 512.0_dp/9.0_dp 
       C(3,4) =  - 64.0_dp/9.0_dp 
       C(3,5) = 80.0_dp/9.0_dp 
       C(3,6) = 8.0_dp/9.0_dp 
       C(4,1) = 8.0_dp/9.0_dp 
       C(4,2) = 80.0_dp/9.0_dp 
       C(4,3) =  - 64.0_dp/9.0_dp 
       C(4,4) = 512.0_dp/9.0_dp 
       C(4,5) = 8.0_dp/9.0_dp 
       C(4,6) =  - 64.0_dp/9.0_dp 
       C(5,1) = 8.0_dp/9.0_dp 
       C(5,2) =  - 64.0_dp/9.0_dp 
       C(5,3) = 80.0_dp/9.0_dp 
       C(5,4) = 8.0_dp/9.0_dp 
       C(5,5) = 512.0_dp/9.0_dp 
       C(5,6) =  - 64.0_dp/9.0_dp 
       C(6,1) = 80.0_dp/9.0_dp 
       C(6,2) = 8.0_dp/9.0_dp 
       C(6,3) = 8.0_dp/9.0_dp 
       C(6,4) =  - 64.0_dp/9.0_dp 
       C(6,5) =  - 64.0_dp/9.0_dp 
       C(6,6) = 512.0_dp/9.0_dp 


      end subroutine cc_gg_ttgamg


      subroutine cc_gg_ttgamg_soft(n,C) ! soft cc matrix, gg-> tt+gamma+g
      integer, intent(in) :: n
      real(dp), intent(out)  :: C(2,2)

      if(n.eq.1) then  
      C(1,1)=8.0_dp/9.0_dp 
      C(1,2)=80.0_dp/9.0_dp 
      C(2,1)=80.0_dp/9.0_dp 
      C(2,2)=8.0_dp/9.0_dp 
      endif 
      if(n.eq.2) then  
      C(1,1)=-8.0_dp 
      C(1,2)=-8.0_dp 
      C(2,1)=-8.0_dp 
      C(2,2)=64.0_dp 
      endif 
      if(n.eq.3) then  
      C(1,1)=64.0_dp 
      C(1,2)=-8.0_dp 
      C(2,1)=-8.0_dp 
      C(2,2)=-8.0_dp 
      endif 
      if(n.eq.4) then  
      C(1,1)=8.0_dp/9.0_dp 
      C(1,2)=80.0_dp/9.0_dp 
      C(2,1)=80.0_dp/9.0_dp 
      C(2,2)=8.0_dp/9.0_dp 
      endif 
      if(n.eq.5) then  
      C(1,1)=64.0_dp 
      C(1,2)=-8.0_dp 
      C(2,1)=-8.0_dp 
      C(2,2)=-8.0_dp 
      endif 
      if(n.eq.6) then  
      C(1,1)=-8.0_dp 
      C(1,2)=-8.0_dp 
      C(2,1)=-8.0_dp 
      C(2,2)=64.0_dp 
      endif 
      if(n.eq.7) then  
      C(1,1)=-8.0_dp 
      C(1,2)=-8.0_dp 
      C(2,1)=-8.0_dp 
      C(2,2)=64.0_dp
      endif 
      if(n.eq.8) then 
      C(1,1)=64.0_dp 
      C(1,2)=-8.0_dp 
      C(2,1)=-8.0_dp 
      C(2,2)=-8.0_dp
      endif 
      if(n.eq.9) then 
      C(1,1)=72.0_dp 
      C(1,2)=0.0_dp 
      C(2,1)=0.0_dp 
      C(2,2)=72.0_dp
      endif 
      if(n.eq.10) then  
      C(1,1)=64.0_dp 
      C(1,2)=-8.0_dp 
      C(2,1)=-8.0_dp 
      C(2,2)=-8.0_dp 
      endif 
      if(n.eq.11) then  
      C(1,1)=-8.0_dp 
      C(1,2)=-8.0_dp 
      C(2,1)=-8.0_dp 
      C(2,2)=64.0_dp 
      endif 
      if(n.eq.12) then  
      C(1,1)=72.0_dp 
      C(1,2)=0.0_dp 
      C(2,1)=0.0_dp 
      C(2,2)=72.0_dp 
      endif 


      end subroutine cc_gg_ttgamg_soft


      subroutine cc_gg_ttgam(C)   
      real(dp), intent(out) :: C(2,2) 


       C(1,1) = 64.0_dp/3.0_dp 
       C(1,2) =  - 8.0_dp/3 
       C(2,1) =  - 8.0_dp/3.0_dp 
       C(2,2) = 64.0_dp/3.0_dp 

      end subroutine cc_gg_ttgam



      end module colorcorr_ttgam

