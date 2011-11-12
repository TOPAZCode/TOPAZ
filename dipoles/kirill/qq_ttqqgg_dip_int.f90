module qq_ttqqgg_dipint
use types
use consts_dp
use mod_dipoles
use colorcorr
use ModAmplitudes
use ModMisc
implicit none
private

public:: qq_ttqqgg_dip_int


contains

      subroutine qq_ttqqgg_dip_int(p,res) 
      real(dp), intent(in) :: p(4,5) 
      real(dp), intent(out) :: res
      complex(dp) :: cres 
      type(Particle) :: ExtParticles(1:5)
      type(TreeProcess) :: TreeAmpsDip(1:4)
      complex(dp) :: Bm(1:4)
      complex(dp) :: AM(1:4,1:4)
      integer :: i1,i2,i3,i4,hel(5),i,j,i5,i6,emi,n, in1, in2,iTree
      real(dp) :: dipsoft, dipfini, dipplus
      real(dp) :: z, C(4,4), fx1(-5:5), fx2(-5:5),fx1z(-5:5),fx2z(-5:5)
      real(dp) :: mtrs
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

! init external particles 
      call InitProcess_TbTQbQG(ExtParticles(1:5))           
! init tree processes for 0-> tb t g g g g
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

      call GenerateEventttqqg(p,hel,ExtParticles(1:5))

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


   do n=1,37 

        call cc_qq_ttqqgg(n,C)

        mtrs = zero

      do i=1,4
      do j=1,4
    mtrs = mtrs + C(i,j)*real(AM(i,j),dp)
      enddo 
      enddo 


        z = 0.5_dp

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
      dipsoft =if_qq(zero,mt,p,4,1,z,1) 
      dipfini =if_qq(zero,mt,p,4,1,z,2) 
      dipplus =if_qq(zero,mt,p,4,1,z,3) 
      diptype ='if' 
     emi = 2 ! mom #4 is emitting 
      endif 
      if(n.eq.6) then  
      dipsoft =if_qq(zero,mt,p,4,2,z,1) 
      dipfini =if_qq(zero,mt,p,4,2,z,2) 
      dipplus =if_qq(zero,mt,p,4,2,z,3) 
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
      dipsoft =fi_qq(mt,zero,p,1,4,z,1) 
      dipfini =fi_qq(mt,zero,p,1,4,z,2) 
      dipplus =fi_qq(mt,zero,p,1,4,z,3) 
      diptype ='fi' 
      emi = 3 ! final state dipole 
      endif 
      if(n.eq.12) then  
      dipsoft =ff_qq(mt,zero,p,1,5,z,1) 
      dipfini =ff_qq(mt,zero,p,1,5,z,2) 
      dipplus =ff_qq(mt,zero,p,1,5,z,3) 
      diptype ='ff' 
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
      dipsoft =fi_qq(mt,zero,p,2,4,z,1) 
      dipfini =fi_qq(mt,zero,p,2,4,z,2) 
      dipplus =fi_qq(mt,zero,p,2,4,z,3) 
      diptype ='fi' 
      emi = 3 ! final state dipole 
      endif 
      if(n.eq.16) then  
      dipsoft =ff_qq(mt,zero,p,2,5,z,1) 
      dipfini =ff_qq(mt,zero,p,2,5,z,2) 
      dipplus =ff_qq(mt,zero,p,2,5,z,3) 
      diptype ='ff' 
      emi = 3 ! final state dipole 
      endif 
      if(n.eq.17) then  
      dipsoft =ff_gg(zero,mt,p,5,1,z,1) 
      dipfini =ff_gg(zero,mt,p,5,1,z,2) 
      dipplus =ff_gg(zero,mt,p,5,1,z,3) 
      diptype ='ff' 
      emi = 3 ! final state dipole 
      endif 
      if(n.eq.18) then  
      dipsoft =ff_gg(zero,mt,p,5,2,z,1) 
      dipfini =ff_gg(zero,mt,p,5,2,z,2) 
      dipplus =ff_gg(zero,mt,p,5,2,z,3) 
      diptype ='ff' 
      emi = 3 ! final state dipole 
      endif 
      if(n.eq.19) then  
      dipsoft =fi_gg(zero,zero,p,5,3,z,1) 
      dipfini =fi_gg(zero,zero,p,5,3,z,2) 
      dipplus =fi_gg(zero,zero,p,5,3,z,3) 
      diptype ='fi' 
      emi = 3 ! final state dipole 
      endif 
      if(n.eq.20) then  
      dipsoft =fi_gg(zero,zero,p,5,4,z,1) 
      dipfini =fi_gg(zero,zero,p,5,4,z,2) 
      dipplus =fi_gg(zero,zero,p,5,4,z,3) 
      diptype ='fi' 
      emi = 3 ! final state dipole 
      endif 
      if(n.eq.21) then  
      dipsoft =if_qq(zero,mt,p,3,1,z,1) 
      dipfini =if_qq(zero,mt,p,3,1,z,2) 
      dipplus =if_qq(zero,mt,p,3,1,z,3) 
      diptype ='if' 
     emi = 1 ! mom #3 is emitting 
      endif 
      if(n.eq.22) then  
      dipsoft =if_qq(zero,mt,p,3,2,z,1) 
      dipfini =if_qq(zero,mt,p,3,2,z,2) 
      dipplus =if_qq(zero,mt,p,3,2,z,3) 
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
      dipsoft =if_qq(zero,mt,p,4,1,z,1) 
      dipfini =if_qq(zero,mt,p,4,1,z,2) 
      dipplus =if_qq(zero,mt,p,4,1,z,3) 
      diptype ='if' 
     emi = 2 ! mom #4 is emitting 
      endif 
      if(n.eq.26) then  
      dipsoft =if_qq(zero,mt,p,4,2,z,1) 
      dipfini =if_qq(zero,mt,p,4,2,z,2) 
      dipplus =if_qq(zero,mt,p,4,2,z,3) 
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
      dipsoft =ff_qq(mt,mt,p,1,2,z,1) 
      dipfini =ff_qq(mt,mt,p,1,2,z,2) 
      dipplus =ff_qq(mt,mt,p,1,2,z,3) 
      diptype ='ff' 
      emi = 3 ! final state dipole 
      endif 
      if(n.eq.30) then  
      dipsoft =fi_qq(mt,zero,p,1,3,z,1) 
      dipfini =fi_qq(mt,zero,p,1,3,z,2) 
      dipplus =fi_qq(mt,zero,p,1,3,z,3) 
      diptype ='fi' 
      emi = 3 ! final state dipole 
      endif 
      if(n.eq.31) then  
      dipsoft =fi_qq(mt,zero,p,1,4,z,1) 
      dipfini =fi_qq(mt,zero,p,1,4,z,2) 
      dipplus =fi_qq(mt,zero,p,1,4,z,3) 
      diptype ='fi' 
      emi = 3 ! final state dipole 
      endif 
      if(n.eq.32) then  
      dipsoft =ff_qq(mt,zero,p,1,5,z,1) 
      dipfini =ff_qq(mt,zero,p,1,5,z,2) 
      dipplus =ff_qq(mt,zero,p,1,5,z,3) 
      diptype ='ff' 
      emi = 3 ! final state dipole 
      endif 
      if(n.eq.33) then  
      dipsoft =ff_qq(mt,mt,p,2,1,z,1) 
      dipfini =ff_qq(mt,mt,p,2,1,z,2) 
      dipplus =ff_qq(mt,mt,p,2,1,z,3) 
      diptype ='ff' 
      emi = 3 ! final state dipole 
      endif 
      if(n.eq.34) then  
      dipsoft =fi_qq(mt,zero,p,2,3,z,1) 
      dipfini =fi_qq(mt,zero,p,2,3,z,2) 
      dipplus =fi_qq(mt,zero,p,2,3,z,3) 
      diptype ='fi' 
      emi = 3 ! final state dipole 
      endif 
      if(n.eq.35) then  
      dipsoft =fi_qq(mt,zero,p,2,4,z,1) 
      dipfini =fi_qq(mt,zero,p,2,4,z,2) 
      dipplus =fi_qq(mt,zero,p,2,4,z,3) 
      diptype ='fi' 
      emi = 3 ! final state dipole 
      endif 
      if(n.eq.36) then  
      dipsoft =ff_qq(mt,zero,p,2,5,z,1) 
      dipfini =ff_qq(mt,zero,p,2,5,z,2) 
      dipplus =ff_qq(mt,zero,p,2,5,z,3) 
      diptype ='ff' 
      emi = 3 ! final state dipole 
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

   enddo


  end subroutine qq_ttqqgg_dip_int

end  module qq_ttqqgg_dipint


