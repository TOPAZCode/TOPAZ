module ModIntDipoles_QGTTBQZ
use ModParameters
use ModProcess
use ModMisc
implicit none
private

public:: EvalIntDipoles_QGTTBQZ
integer, parameter  :: dp = selected_real_kind(15)



 contains


!-------------------------------------------------------------
!    ordering of the momenta in the p-array:
!     outgoing convention:
!     bar t = 1, Z =2, t=3 ,  incoming quark=4,
!     5=incoming gluon
!-----------------------------------------------------------

      subroutine EvalIntDipoles_QGTTBQZ(p,MomDK,z,npdf,res)
      use ModIntDipoles
      use ModAmplitudes
      use ModZDecay
      implicit none
      real(dp), intent(in) :: p(4,5),MomDK(1:4,1:8),z
      real(dp), intent(out) :: res(1:2,1:3)
      real(dp) :: dipsoft,dipfini,dipplus,AP(1:3),CF,TR
      type(Particle),save :: ExtParticles1(1:5)
      type(TreeProcess),save :: TreeAmpsDip1(1:2)
      type(Particle),save :: ExtParticles2(1:5)
      type(TreeProcess),save :: TreeAmpsDip2(1:2)
      complex(dp) :: B(1:2),Bm(1:2,1:2),mtrsq(2),MZ_Inv,propZ
      complex(dp) :: AM(1:2,1:2,1:2)
      integer :: iTree,i1,i2,i3,i4,i5,i,j,i6,n,npdf
      integer :: hel(1:5),emi , Nmax(5),N2jump
      character(2) :: diptype
      real(dp) :: C(1:2,1:2),epcorr,couplZLL,couplGLL,couplZUU,couplGUU,couplZDD,couplGDD
      logical, save :: first_time1 = .true.
      logical, save :: first_time2 = .true.
      real,parameter :: PhotonCouplCorr=2d0


      res(1:2,1:3) = zero
      CF=4d0/3d0
      TR = 0.5d0

      if( npdf.eq.2 ) then
            call swapMom(p(1:4,4),p(1:4,5))
      endif


         Nmax = 1

        if (TopDecays.ge.1) then
            Nmax(1) = -1
            Nmax(3) = -1
        endif

       if( first_time1 ) then
             call InitTrees(2,2,2,TreeAmpsDip1,NumBoson=1)
             call InitProcess_TbTGGZ(ExtParticles1(1:5))
             TreeAmpsDip1(1)%PartRef(1:5) = (/1,5,2,3,4/)
             TreeAmpsDip1(2)%PartRef(1:5) = (/1,5,2,4,3/)
             do iTree=1,2
                call LinkTreeParticles(TreeAmpsDip1(iTree),ExtParticles1(1:5))
             enddo
             first_time1=.false.
      endif
      call cc_gg_ttgam(C)


!--- after momentum mapping -- sum over colors and polarizations


        mtrsq = 0d0

       N2Jump = 1
       if (ZDecays.ge.1 .or. ZDecays.eq.-2) then ! Z decays
            N2jump=2
       endif
       do i1=-1,Nmax(1),2
          do i2 = -1,Nmax(2),N2jump!  Z boson
           do i3 = -1,Nmax(3),2
             do i4 = -1,Nmax(4),2
               do i5 = -1,Nmax(5),2

          hel(1) = i1
          hel(2) = i2
          hel(3) = i3
          hel(4) = i4
          hel(5) = i5

    call SetPolarization_GG((/p(1:4,1),p(1:4,3),p(1:4,4),p(1:4,5),p(1:4,2)/),momDK(1:4,1:8),(/hel(1),hel(3),hel(4),hel(5),hel(2)/),ExtParticles1(1:5))



          do i6 = 1,2
              call EvalTree2(TreeAmpsDip1(i6),Bm(1,i6))
          enddo

          do i=1,2
          do j=1,2
              mtrsq(1) =  mtrsq(1) + C(i,j)*Bm(1,i)*dconjg(Bm(1,j))
          enddo
          enddo

      enddo
      enddo
      enddo
      enddo
      enddo

! spin & color avg., color convention adjustment
      mtrsq(1) = mtrsq(1) * alpha_s4Pi**2 * alpha4Pi /4d0/64d0
      mtrsq(2) = mtrsq(1) ! here was a bug: this line was missing, q-> g splitting for down quarks



      res = 0d0

      dipsoft =ii_qg(zero,zero,p,4,5,z,1)*2d0*CF
      dipfini =ii_qg(zero,zero,p,4,5,z,2)*2d0*CF
      dipplus =ii_qg(zero,zero,p,4,5,z,3)*2d0*CF
      emi = 1 ! mom #3 is emitting

!       if(emi.eq.1) then
        res(:,1) = res(:,1) + (dipsoft-dipplus)*mtrsq(:)*alpha_sOver2Pi*0.5d0! overall normalization: alpha_s/4pi
        res(:,2) = res(:,2) + (dipfini+dipplus)*mtrsq(:)*alpha_sOver2Pi*0.5d0
!       endif
!       if(emi.eq.2) then
!         res(:,1) = res(:,1) + (dipsoft-dipplus)*mtrsq(:)
!         res(:,3) = res(:,3) + (dipfini+dipplus)*mtrsq(:)
!       endif
!       if(emi.eq.3) then
!         res(:,1) = res(:,1) + (dipsoft-dipplus)*mtrsq(:)
!         res(:,2) = res(:,2) + (dipfini+dipplus)*0.5_dp*mtrsq(:)
!         res(:,3) = res(:,3) + (dipfini+dipplus)*0.5_dp*mtrsq(:)
!       endif


! print *, "1u",res(1,1:3)
! print *, "1d",res(2,1:3)


     epcorr = epinv
     AP(1)= 0d0
     AP(2)= CF*(1d0+(1d0-z)**2)/z*mtrsq(1)*alpha_sOver2Pi*epinv
     AP(3)= 0d0

     res(1:2,1) = res(1:2,1) + (AP(1)-AP(3)) ! here was a bug:  q-> g splitting for down quarks was missing
     res(1:2,2) = res(1:2,2) + (AP(2)+AP(3))
     res(1:2,3) = 0d0

! print *, "2u",res(1,1:3)
! print *, "2d",res(2,1:3)






!-- ---- NEXT DIPOLE

      AM = (zero,zero) ! first, the ``amplitude matrix'' is set to zero

       if( first_time2 ) then
             call InitTrees(4,0,2,TreeAmpsDip2,NumBoson=1)
             call InitProcess_TbTQbQZ(ExtParticles2(1:5))
             TreeAmpsDip2(1)%PartRef(1:5) = (/1,5,2,3,4/)
             TreeAmpsDip2(2)%PartRef(1:5) = (/1,2,3,5,4/)
             do iTree=1,2
              call LinkTreeParticles(TreeAmpsDip2(iTree),ExtParticles2(1:5))
             enddo
             first_time2=.false.
      endif

      call cc_qq_tt(C)


        mtrsq = 0d0


       N2Jump = 1
       if (ZDecays.ge.1 .or. ZDecays.eq.-2) then ! Z decays
            N2jump=2
       endif
       do i1=-1,Nmax(1),2
       do i2 = -1,Nmax(2),N2Jump! Z boson
           do i3 = -1,Nmax(3),2
             do i4 = -1,Nmax(4),2
               do i5 = -1,Nmax(5),2

          hel(1) = i1
          hel(2) = i2
          hel(3) = i3
          hel(4) = i4
          hel(5) = i5

          call SetPolarization_QQ((/p(1:4,1),p(1:4,3),p(1:4,4),p(1:4,5),p(1:4,2)/),momDK(1:4,1:8),(/hel(1),hel(3),hel(4),hel(5),hel(2)/),ExtParticles2(1:5))
          if( ZDecays.gt.0 ) call ZGamLCoupl(1,hel(2),couplZLL,couplGLL)  ! charged lept
          if( npdf.eq.1 ) then
              call ZGamQcoupl(Up_,hel(4),couplZUU,couplGUU)
              call ZGamQcoupl(Dn_,hel(4),couplZDD,couplGDD)
          else
              call ZGamQcoupl(Up_,-hel(4),couplZUU,couplGUU)
              call ZGamQcoupl(Dn_,-hel(4),couplZDD,couplGDD)
          endif

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
              call EvalTree2(TreeAmpsDip2(i6),B(i6))
              if( i6.eq.1 ) then
                  Bm(up_,i6)=B(i6)
                  Bm(dn_,i6)=B(i6)
                  if( Zdecays.gt.0 .and. Zdecays.lt.10 ) then
!                       Bm(up_,i6)=Bm(up_,i6)*propZ*couplZLL
!                       Bm(dn_,i6)=Bm(dn_,i6)*propZ*couplZLL
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

      mtrsq = (zero,zero)
      do i=1,2
      do j=1,2
            mtrsq(up_) = mtrsq(up_) + C(i,j)*real(AM(up_,i,j),dp)
            mtrsq(dn_) = mtrsq(dn_) + C(i,j)*real(AM(dn_,i,j),dp)
      enddo
      enddo


! spin & color avg., color convention adjustment
      mtrsq = mtrsq * alpha_s4Pi**2 * alpha4Pi /4d0/9d0



      dipsoft =ii_gq(zero,zero,p,5,4,z,1) *2d0  ! factor 2d0 fixes normalization for ii_gq dipole
      dipfini =ii_gq(zero,zero,p,5,4,z,2) *2d0
      dipplus =ii_gq(zero,zero,p,5,4,z,3) *2d0
      diptype ='ii'
      emi = 2 !

!       if(emi.eq.2) then
         res(:,1) = res(:,1) + (dipsoft-dipplus)*mtrsq(:) * alpha_sOver2Pi*0.5d0   ! overall normalization: alpha_s/4pi
         res(:,3) = res(:,3) + (dipfini+dipplus)*mtrsq(:) * alpha_sOver2Pi*0.5d0   ! overall normalization: alpha_s/4pi
!       endif


! print *, "3u",res(1,1:3)
! print *, "3d",res(2,1:3)

! !        epcorr=epinv+2d0*dlog(renscale/facscale)
          epcorr=epinv
          AP(1)= 0d0
          AP(2)= TR*(z**2+(1d0-z)**2)* alpha_sOver2Pi *epcorr
          AP(3)= 0d0

          res(1:2,1) = res(1:2,1) + (AP(1)-AP(3))*mtrsq(1:2)
          res(1:2,2) = res(1:2,2)
          res(1:2,3) = res(1:2,3) + (AP(2)+AP(3))*mtrsq(1:2)


! print *, "4u",res(1,1:3)
! print *, "4d",res(2,1:3)
! pause

  return
  end subroutine




SUBROUTINE SetPolarization_QQ(Mom,MomDK,Hel,ExtParticles)
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





SUBROUTINE InitProcess_TbTGGZ(ExtParticles)
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

  ExtParticles(3)%PartType = Glu_
  ExtParticles(3)%ExtRef   = 3
  ExtParticles(3)%Mass = 0d0
  ExtParticles(3)%Mass2= 0d0

  ExtParticles(4)%PartType = Glu_
  ExtParticles(4)%ExtRef   = 4
  ExtParticles(4)%Mass = 0d0
  ExtParticles(4)%Mass2= 0d0

  ExtParticles(5)%PartType = Z0_
  ExtParticles(5)%ExtRef   = 5
  ExtParticles(5)%Mass = m_Z
  ExtParticles(5)%Mass2= ExtParticles(5)%Mass**2


RETURN
END SUBROUTINE InitProcess_TbTGGZ



SUBROUTINE SetPolarization_GG(Mom,MomDK,Hel,ExtParticles)
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
     else
          call pol_massSR(ExtParticles(5)%Mom(1:4),ExtParticles(5)%Mass,ExtParticles(5)%Helicity,ExtParticles(5)%Pol(1:4))
    endif

    ExtParticles(3)%Mom(1:4) = dcmplx(Mom(1:4,3))
    call pol_mless(ExtParticles(3)%Mom(1:4),Hel(3),ExtParticles(3)%Pol(1:4))

    ExtParticles(4)%Mom(1:4) = dcmplx(Mom(1:4,4))
    call pol_mless(ExtParticles(4)%Mom(1:4),Hel(4),ExtParticles(4)%Pol(1:4))


! ExtParticles(3)%Pol(1:4)=ExtParticles(3)%Mom(1:4); print *, "check gauge inv. in dipoles"

RETURN
END SUBROUTINE




      subroutine cc_gg_ttgam(C)
      real(dp), intent(out) :: C(2,2)


       C(1,1) = 64.0_dp/3.0_dp
       C(1,2) =  - 8.0_dp/3
       C(2,1) =  - 8.0_dp/3.0_dp
       C(2,2) = 64.0_dp/3.0_dp

      end subroutine cc_gg_ttgam



      subroutine cc_qq_tt(C)
      real(dp), intent(out) :: C(2,2)

      C = 0.0_dp


       C(1,1)  =  8.0_dp
       C(1,2)  =  8.0_dp
       C(2,1)  =  8.0_dp
       C(2,2)  =  8.0_dp

      end subroutine



end  module
