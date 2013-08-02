module ModIntDipoles_QQBTTBGZ
use ModParameters
use ModProcess
use ModMisc
implicit none
private

public:: EvalIntDipoles_QQBTTBGZ
integer, parameter  :: dp = selected_real_kind(15)


!       double precision, private, parameter :: NCol=3d0
!       double precision, private, parameter :: TR=0.5d0
!       double precision, private, parameter :: CA=2d0*TR*NCol
!       double precision, private, parameter :: CF=TR*(NCol**2-1d0)/NCol

contains


!-------------------------------------------------------------
!    ordering of the momenta in the p-array:
!     outgoing convention:
!     bar t = 1, photon =2, t=3 ,  incoming quark=4,
!     5=incoming antiquark
!-----------------------------------------------------------

      subroutine EvalIntDipoles_QQBTTBGZ(p,MomDK,z,res)
      use ModIntDipoles
      use ModAmplitudes
      use ModZDecay
      implicit none
      real(dp), intent(in) :: p(4,5),MomDK(1:4,1:8),z
      real(dp), intent(out) :: res(1:2,1:3)
      real(dp) :: dipsoft,dipfini,dipplus,AP(1:2,1:3),CF
      type(Particle),save :: ExtParticles(1:5)
      type(TreeProcess),save :: TreeAmpsDip(1:2)
      complex(dp) :: B(1:2),Bm(1:2,1:2),mtrsq(2)
      complex(dp) :: AM(1:2,1:2,1:2),propZ,MZ_Inv
      integer :: iTree,i1,i2,i3,i4,i5,i,j,i6,n
      integer :: hel(1:5),emi , Nmax(5),N2Jump
      character(2) :: diptype
      real(dp) :: colcorr(1:2,1:2),epcorr,couplZLL,couplGLL,couplZUU,couplGUU,couplZDD,couplGDD
      logical, save :: first_time = .true.

      res(1:2,1:3) = zero
      CF=4d0/3d0

       if( first_time ) then
             call InitTrees(4,0,2,TreeAmpsDip,NumBoson=1)
             call InitProcess_TbTQbQZ(ExtParticles(1:5))
             TreeAmpsDip(1)%PartRef(1:5) = (/1,5,2,3,4/)
             TreeAmpsDip(2)%PartRef(1:5) = (/1,2,3,5,4/)
             do iTree=1,2
              call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles(1:5))
             enddo
             first_time=.false.
      endif


      AM = (zero,zero) ! first, the ``amplitude matrix'' is set to zero

!--- after momentum mapping -- sum over colors and polarizations

         Nmax = 1
        if (TopDecays.ge.1) then
            Nmax(1) = -1
            Nmax(3) = -1
        endif


       N2Jump = 1
       if (ZDecays.ge.1 .or. ZDecays.eq.-2) then ! Z decays
            N2jump=2
       endif
       do i1=-1,Nmax(1),2!   top
       do i2 = -1,Nmax(2),N2Jump! Z boson
       do i3 = -1,Nmax(3),2 !  top
       do i4 = -1,Nmax(4),2
       do i5 = -1,Nmax(5),2


          hel(1) = i1
          hel(2) = i2
          hel(3) = i3
          hel(4) = i4
          hel(5) = i5

          call SetPolarization((/p(1:4,1),p(1:4,3),p(1:4,4),p(1:4,5),p(1:4,2)/),momDK(1:4,1:8),(/hel(1),hel(3),hel(4),hel(5),hel(2)/),ExtParticles(1:5))
          if( ZDecays.gt.0 ) call ZGamLCoupl(1,hel(2),couplZLL,couplGLL)  ! charged lept
          call ZGamQcoupl(Up_,ExtParticles(3)%Helicity,couplZUU,couplGUU)
          call ZGamQcoupl(Dn_,ExtParticles(3)%Helicity,couplZDD,couplGDD)

          couplZQQ_left_dyn=one
          couplZQQ_right_dyn=one
          if( ZDecays.le.10 ) then  ! decaying on-shell Z
                MZ_Inv = m_Z
          elseif( ZDecays.gt.10 ) then  ! decaying off-shell Z
                call Error("need to implement phase space for off-shell Z's")
                ! need to think about threshold cut: EHat.le.2d0*m_Top+M_Z   when z is off-shell!
          endif
          if ( ZDecays.lt.10) then
              propZ = (1d0,0d0)/dsqrt(2d0*Ga_Zexp*m_Z)  * MZ_Inv**2
          elseif (ZDecays .gt. 10) then 
              propZ=cone/(MZ_Inv**2-m_Z**2+ci*Ga_ZExp*m_Z)   * MZ_Inv**2
          endif

          do i6 = 1,2
              call EvalTree2(TreeAmpsDip(i6),B(i6))
              if( i6.eq.1 ) then
                  Bm(up_,i6)=B(i6)
                  Bm(dn_,i6)=B(i6)
                  if( Zdecays.gt.0 .and. Zdecays.lt.10 ) then
                      Bm(up_,i6)=Bm(up_,i6)*propZ*couplZLL
                      Bm(dn_,i6)=Bm(dn_,i6)*propZ*couplZLL
                  endif
              else
                  Bm(up_,i6)=B(i6)*couplZUU
                  Bm(dn_,i6)=B(i6)*couplZDD
                  if( Zdecays.gt.0 .and. Zdecays.lt.10 ) then
                      Bm(up_,i6)=Bm(up_,i6)*propZ*couplZLL
                      Bm(dn_,i6)=Bm(dn_,i6)*propZ*couplZLL
                  elseif( Zdecays.gt.10 ) then
!                       Bm(up_,i6)=Bm(up_,i6)*propZ*( couplZUU*propZ*couplZLL + couplGUU*couplGLL )  /couplZUU
!                       Bm(dn_,i6)=Bm(dn_,i6)*propZ*( couplZUU*propZ*couplZLL + couplGUU*couplGLL )  /couplZDD
                      call Error("Zdecays.gt.10 not yet implemented")
                  endif
              endif
          enddo



          do i=1,2
          do j=1,2
              AM(up_,i,j) = AM(up_,i,j) + Bm(up_,i)*dconjg(Bm(up_,j))
              AM(dn_,i,j) = AM(dn_,i,j) + Bm(dn_,i)*dconjg(Bm(dn_,j))
          enddo
          enddo

      enddo
      enddo
      enddo
      enddo
      enddo


   do n=1,12

     mtrsq = (zero,zero)
     call cc_qq_ttgamg_soft(n,colcorr)
     do i=1,2
     do j=1,2
          mtrsq(up_) = mtrsq(up_) + colcorr(i,j)*real(AM(up_,i,j),dp)
          mtrsq(dn_) = mtrsq(dn_) + colcorr(i,j)*real(AM(dn_,i,j),dp)
     enddo
     enddo

! spin & color avg., color convention adjustment

     mtrsq = mtrsq * alpha_s4Pi**2 * alpha4Pi  /4d0/9d0


      if(n.eq.1) then
      dipsoft =ff_qq(m_TOP,m_TOP,p,1,3,z,1)
      dipfini =ff_qq(m_TOP,m_TOP,p,1,3,z,2)
      dipplus =ff_qq(m_TOP,m_TOP,p,1,3,z,3)
      emi = 3 ! final state dipole
      endif
      if(n.eq.2) then
      dipsoft =fi_qq(m_TOP,zero,p,1,4,z,1)
      dipfini =fi_qq(m_TOP,zero,p,1,4,z,2)
      dipplus =fi_qq(m_TOP,zero,p,1,4,z,3)
!      emi = 3 ! final state dipole
      emi = 1
      endif
      if(n.eq.3) then
      dipsoft =fi_qq(m_TOP,zero,p,1,5,z,1)
      dipfini =fi_qq(m_TOP,zero,p,1,5,z,2)
      dipplus =fi_qq(m_TOP,zero,p,1,5,z,3)
!      emi = 3 ! final state dipole
      emi = 2
      endif
      if(n.eq.4) then
      dipsoft =ff_qq(m_TOP,m_TOP,p,3,1,z,1)
      dipfini =ff_qq(m_TOP,m_TOP,p,3,1,z,2)
      dipplus =ff_qq(m_TOP,m_TOP,p,3,1,z,3)
      emi = 3 ! final state dipole
      endif
      if(n.eq.5) then
      dipsoft =fi_qq(m_TOP,zero,p,3,4,z,1)
      dipfini =fi_qq(m_TOP,zero,p,3,4,z,2)
      dipplus =fi_qq(m_TOP,zero,p,3,4,z,3)
!      emi = 3 ! final state dipole
      emi = 1
      endif
      if(n.eq.6) then
      dipsoft =fi_qq(m_TOP,zero,p,3,5,z,1)
      dipfini =fi_qq(m_TOP,zero,p,3,5,z,2)
      dipplus =fi_qq(m_TOP,zero,p,3,5,z,3)
!      emi = 3 ! final state dipole
      emi = 2
      endif
      if(n.eq.7) then
      dipsoft =if_qq(zero,m_TOP,p,4,1,z,1)
      dipfini =if_qq(zero,m_TOP,p,4,1,z,2)
      dipplus =if_qq(zero,m_TOP,p,4,1,z,3)
      emi = 1 ! pdf1 is emitting
      endif
      if(n.eq.8) then
      dipsoft =if_qq(zero,m_TOP,p,4,3,z,1)
      dipfini =if_qq(zero,m_TOP,p,4,3,z,2)
      dipplus =if_qq(zero,m_TOP,p,4,3,z,3)
     emi = 1 ! pdf1 is emitting
      endif
      if(n.eq.9) then
      dipsoft =ii_qq(zero,zero,p,4,5,z,1)
      dipfini =ii_qq(zero,zero,p,4,5,z,2)
      dipplus =ii_qq(zero,zero,p,4,5,z,3)
     emi = 1 ! pdf1 is emitting
      endif
      if(n.eq.10) then
      dipsoft =if_qq(zero,m_TOP,p,5,1,z,1)
      dipfini =if_qq(zero,m_TOP,p,5,1,z,2)
      dipplus =if_qq(zero,m_TOP,p,5,1,z,3)
     emi = 2 ! pdf2 is emitting
      endif
      if(n.eq.11) then
      dipsoft =if_qq(zero,m_TOP,p,5,3,z,1)
      dipfini =if_qq(zero,m_TOP,p,5,3,z,2)
      dipplus =if_qq(zero,m_TOP,p,5,3,z,3)
     emi = 2 ! pdf2 is emitting
      endif
      if(n.eq.12) then
      dipsoft =ii_qq(zero,zero,p,5,4,z,1)
      dipfini =ii_qq(zero,zero,p,5,4,z,2)
      dipplus =ii_qq(zero,zero,p,5,4,z,3)
     emi = 2 ! pdf2 is emitting
      endif


! print *, "comment here";  RES(:,1) = RES(:,1) +  DIPSOFT*MTRSQ(:)  ! FOR DELTA-FCT. CHECK
! print *, "dipsoft=0";   dipsoft=0d0

      if(emi.eq.1) then
        res(:,1) = res(:,1) + (dipsoft-dipplus)*mtrsq(:)
        res(:,2) = res(:,2) + (dipfini+dipplus)*mtrsq(:)
      endif
      if(emi.eq.2) then
        res(:,1) = res(:,1) + (dipsoft-dipplus)*mtrsq(:)
        res(:,3) = res(:,3) + (dipfini+dipplus)*mtrsq(:)
      endif
      if(emi.eq.3) then
        res(:,1) = res(:,1) + (dipsoft-dipplus)*mtrsq(:)
        res(:,2) = res(:,2) + (dipfini+dipplus)*0.5_dp*mtrsq(:)
        res(:,3) = res(:,3) + (dipfini+dipplus)*0.5_dp*mtrsq(:)
      endif

   enddo

! overall normalization: alpha_s/4pi
   res(1:2,1:3)=alpha_sOver2Pi*0.5d0*res(1:2,1:3)

! print *, "1",res(1:2,1)/alpha_sOver2Pi

!    eval tree
     mtrsq = (zero,zero)
     call cc_qq_ttgamg_soft(0,colcorr)
     do i=1,2
     do j=1,2
          mtrsq(1) = mtrsq(1) + colcorr(i,j)*real(AM(up_,i,j),dp)
          mtrsq(2) = mtrsq(2) + colcorr(i,j)*real(AM(dn_,i,j),dp)
     enddo
     enddo
     mtrsq = mtrsq * alpha_s4Pi**2 * alpha4Pi  /4d0/9d0
! print *, "ID tree",mtrsq(1),mtrsq(2)



! !        epcorr=epinv+2d0*dlog(renscale/facscale)
       epcorr=epinv
       AP(:,1)= 3d0/2d0*CF
       AP(:,2)= (-1d0-z)*CF
       AP(:,3)= 2d0*CF/(1d0-z)
       AP(1,1:3) = AP(1,1:3) * alpha_sOver2Pi *epcorr * mtrsq(1)
       AP(2,1:3) = AP(2,1:3) * alpha_sOver2Pi *epcorr * mtrsq(2)

! print *, "comment here";  RES(:,1) = RES(:,1) + 2D0*AP(:,1)  ! FOR DELTA-FCT. CHECK
! print *, "AP(1)=0";   AP(:,1)=0d0
       res(1:2,1) = res(1:2,1) + 2d0*(AP(1:2,1)-AP(1:2,3)) ! here was a bug 2d0*AP(1:2,1)-AP(1:2,3)
       res(1:2,2) = res(1:2,2) + AP(1:2,2) + AP(1:2,3)
       res(1:2,3) = res(1:2,3) + AP(1:2,2) + AP(1:2,3)



! print *, "2",res(1:2,1:3)

! print *, "RES(1,1)/mtrsq",RES(1,1)/(mtrsq(1)*alpha_sOver2Pi)
! print *, "RES(2,1)/mtrsq",RES(2,1)/(mtrsq(2)*alpha_sOver2Pi)
! pause
! RES(1:2,1)=mtrsq(1:2)

  return
  end subroutine




      subroutine cc_qq_ttgamg_soft(n,C)  ! ixq labels the power of the
      integer, intent(in) :: n               ! electric charge of the
      real(dp), intent(out) :: C(2,2)        ! emitting quark


      C = 0.0_dp


      if(n.eq.0) then! 0= tree level correlation
      C(1,1)=8.0_dp
      C(1,2)=8.0_dp
      C(2,1)=8.0_dp
      C(2,2)=8.0_dp
      endif
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

      end subroutine





SUBROUTINE SetPolarization(Mom,MomDK,Hel,ExtParticles)
use ModMisc
use ModProcess
use ModTopDecay
use ModZDecay
implicit none
type(Particle) :: ExtParticles(1:5)
real(8) :: Mom(1:4,1:5),MomDK(1:4,1:8)
integer :: Hel(1:5)

     ExtParticles(1)%Mom(1:4) = dcmplx(Mom(1:4,1))   ! HERE WAS A BUG: this was inside the (TopDecays.ge.1) condition
     ExtParticles(2)%Mom(1:4) = dcmplx(Mom(1:4,2))
     if (TopDecays.ge.1) then
        call TopDecay(ExtParticles(1),DK_LO,MomDK(1:4,1:3))
        call TopDecay(ExtParticles(2),DK_LO,MomDK(1:4,4:6))
     else
        call vSpi(ExtParticles(1)%Mom(1:4),ExtParticles(1)%Mass,Hel(1),ExtParticles(1)%Pol(1:4))
        call ubarSpi(ExtParticles(2)%Mom(1:4),ExtParticles(2)%Mass,Hel(2),ExtParticles(2)%Pol(1:4))
    endif


     ExtParticles(5)%Mom(1:4) = dcmplx(Mom(1:4,5))
     ExtParticles(5)%Helicity = Hel(5)
     if (ZDecays.ge.1) then
          call ZDecay(ExtParticles(5),DK_LO,MomDK(1:4,7:8))
     elseif (ZDecays.eq.0) then
          call pol_massSR(ExtParticles(5)%Mom(1:4),ExtParticles(5)%Mass,ExtParticles(5)%Helicity,ExtParticles(5)%Pol(1:4))
     elseif (ZDecays.eq.-2) then
         call pol_mless(ExtParticles(5)%Mom(1:4),ExtParticles(5)%Helicity,ExtParticles(5)%Pol(1:4))
    endif


    ExtParticles(3)%Mom(1:4) = dcmplx(Mom(1:4,3))
    ExtParticles(3)%Helicity = Hel(3)
    call vSpi(ExtParticles(3)%Mom(1:4),ExtParticles(3)%Mass,Hel(3),ExtParticles(3)%Pol(1:4))

    ExtParticles(4)%Mom(1:4) = dcmplx(Mom(1:4,4))
    ExtParticles(4)%Helicity = Hel(4)
    call ubarSpi(ExtParticles(4)%Mom(1:4),ExtParticles(4)%Mass,Hel(4),ExtParticles(4)%Pol(1:4))


RETURN
END SUBROUTINE





SUBROUTINE InitProcess_TbTQbQZ(ExtParticles)
use ModProcess
implicit none
type(Particle) :: ExtParticles(:)

  ExtParticles(1)%PartType = ATop_
  ExtParticles(1)%ExtRef   = 1
  ExtParticles(1)%Mass = m_Top
  ExtParticles(1)%Mass2= ExtParticles(1)%Mass**2

  ExtParticles(2)%PartType = Top_
  ExtParticles(2)%ExtRef   = 2
  ExtParticles(2)%Mass = m_Top
  ExtParticles(2)%Mass2= ExtParticles(2)%Mass**2

  ExtParticles(3)%PartType = AStr_
  ExtParticles(3)%ExtRef   = 3
  ExtParticles(3)%Mass = 0d0
  ExtParticles(3)%Mass2= 0d0

  ExtParticles(4)%PartType = Str_
  ExtParticles(4)%ExtRef   = 4
  ExtParticles(4)%Mass = 0d0
  ExtParticles(4)%Mass2= 0d0

  ExtParticles(5)%PartType = Z0_
  if( Process.eq.86 ) ExtParticles(5)%PartType = Pho_
  ExtParticles(5)%ExtRef   = 5
  ExtParticles(5)%Mass = m_Z
  ExtParticles(5)%Mass2= ExtParticles(5)%Mass**2

RETURN
END SUBROUTINE




end  module
