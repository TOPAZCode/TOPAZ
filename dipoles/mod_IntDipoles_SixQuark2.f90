module ModSixFermionIntDip
use ModAmplitudes
use ModParameters
use ModProcess
use ModMisc
use ModIntDipoles

implicit none

integer, parameter  :: dp = selected_real_kind(15)

private

public:: sixquark_intdip

contains

!  momentum input for all these  subroutines in all-outgoing convention
!  should be :
!  tbar     1
!  t        2
!  initial  3
!  initial  4
!  final    5


  subroutine sixquark_intdip(p,MomDK,z,in1,in2,res)
  implicit none
  real(dp), intent(in) :: p(4,5),MomDK(1:4,1:6),z
  integer, intent(in) :: in1, in2
  real(dp), intent(out) :: res(1:3)
  complex(dp) :: cres
  type(Particle) :: ExtParticles1(1:5),ExtParticles2(1:5),ExtParticles3(1:5)
  type(TreeProcess) :: TreeAmpsDip(1:4),TreeAmpsDipg(1:4),TreeAmpsDipq(1:6)
  complex(dp) :: Bm(1:4),Bq(1:6)
  complex(dp) :: AM(1:4,1:4)
  integer :: i1,i2,i3,i4,hel(5),i,j,i5,i6,emi,n,iTree
  real(dp) :: dipsoft_ii, dipfini_ii, dipplus_ii
  real(dp) :: dipsoft_fi, dipfini_fi, dipplus_fi
  real(dp) :: C(4,4),AP(1:3)  !, fx1(-5:5), fx2(-5:5),fx1z(-5:5),fx2z(-5:5)
  real(dp) :: Cq(6,6)
  real(dp) :: mtrs(2)
  real(dp) :: q(4,5)
  real(dp) :: p1(4,5),CF,epcorr
  character :: diptype*2


       res(1:3) = zero
       CF=4d0/3d0

      dipsoft_ii =ii_qg(zero,zero,p,3,4,z,1)
      dipfini_ii =ii_qg(zero,zero,p,3,4,z,2)
      dipplus_ii =ii_qg(zero,zero,p,3,4,z,3)



  if (abs(in1).ne.abs(in2)) then
!------------------------------------------------------------------
       if ((in1.gt.0).and.(in2.gt.0)) then    !  Q QP

          call qqq_ttqqqq_dip_int(p,MomDK,mtrs)
          res(1) = res(1) + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2)) * alpha_sOver2Pi*CF
          res(2) = res(2) + (dipfini_ii+dipplus_ii)*mtrs(1)           * alpha_sOver2Pi*CF
          res(3) = res(3) + (dipfini_ii+dipplus_ii)*mtrs(2)           * alpha_sOver2Pi*CF


!         epcorr=epinv+2d0*dlog(renscale/facscale)
          epcorr=epinv
          AP(1)= 0d0
          AP(2)= CF * (1d0+(1d0-z)**2)/z
          AP(3)= 0d0
          AP(1:3) = AP(1:3) * alpha_sOver2Pi * epcorr * mtrs(1)
          res(1) = res(1) + (AP(1)-AP(3))
          res(2) = res(2) + (AP(2)+AP(3))

          AP(1)= 0d0
          AP(2)= CF * (1d0+(1d0-z)**2)/z
          AP(3)= 0d0
          AP(1:3) = AP(1:3) * alpha_sOver2Pi * epcorr * mtrs(2)
          res(1) = res(1) + (AP(1)-AP(3))
          res(3) = res(3) + (AP(2)+AP(3))

!           print *, res(2)
!           print *, res(3)
!           pause


!------------------------------------------------------------------
       elseif ((in1.gt.0).and.(in2.lt.0)) then    !  Q QbarP

          call qqb_ttqqqq_dip_12_int(p,MomDK,mtrs)
          res(1) = res(1) + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2)) * alpha_sOver2Pi*CF
          res(2) = res(2) + (dipfini_ii+dipplus_ii)*mtrs(1)           * alpha_sOver2Pi*CF
          res(3) = res(3) + (dipfini_ii+dipplus_ii)*mtrs(2)           * alpha_sOver2Pi*CF

!           print *, res(2)
!           print *, res(3)

!         epcorr=epinv+2d0*dlog(renscale/facscale)
          epcorr=epinv
          AP(1)= 0d0
          AP(2)= CF * (1d0+(1d0-z)**2)/z
          AP(3)= 0d0
          AP(1:3) = AP(1:3) * alpha_sOver2Pi * epcorr * mtrs(1)
          res(1) = res(1) + (AP(1)-AP(3))
          res(2) = res(2) + (AP(2)+AP(3))

          AP(1)= 0d0
          AP(2)= CF * (1d0+(1d0-z)**2)/z
          AP(3)= 0d0
          AP(1:3) = AP(1:3) * alpha_sOver2Pi * epcorr * mtrs(2)
          res(1) = res(1) + (AP(1)-AP(3))
          res(3) = res(3) + (AP(2)+AP(3))

!           print *, res(2)
!           print *, res(3)
!           pause


!------------------------------------------------------------------
       elseif ((in1.lt.0).and.(in2.gt.0)) then    !  Qbar QP

          p1 = p
          p1(:,3) = p(:,4)
          p1(:,4) = p(:,3)
          call qqb_ttqqqq_dip_12_int(p1,MomDK,mtrs)
          res(1) = res(1) + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2)) * alpha_sOver2Pi*CF
          res(2) = res(2) + (dipfini_ii+dipplus_ii)*mtrs(1)           * alpha_sOver2Pi*CF
          res(3) = res(3) + (dipfini_ii+dipplus_ii)*mtrs(2)           * alpha_sOver2Pi*CF

!           print *, res(2)
!           print *, res(3)

!         epcorr=epinv+2d0*dlog(renscale/facscale)
          epcorr=epinv
          AP(1)= 0d0
          AP(2)= CF * (1d0+(1d0-z)**2)/z
          AP(3)= 0d0
          AP(1:3) = AP(1:3) * alpha_sOver2Pi * epcorr * mtrs(1)
          res(1) = res(1) + (AP(1)-AP(3))
          res(2) = res(2) + (AP(2)+AP(3))

          AP(1)= 0d0
          AP(2)= CF * (1d0+(1d0-z)**2)/z
          AP(3)= 0d0
          AP(1:3) = AP(1:3) * alpha_sOver2Pi * epcorr * mtrs(2)
          res(1) = res(1) + (AP(1)-AP(3))
          res(3) = res(3) + (AP(2)+AP(3))

!           print *, res(2)
!           print *, res(3)
!           pause



!------------------------------------------------------------------
       elseif ((in1.lt.0).and.(in2.lt.0)) then    !  Qbar QbarP

          call qbb_ttqqqq_dip_int(p,MomDK,mtrs)
          res(1) = res(1) + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2)) * alpha_sOver2Pi*CF
          res(2) = res(2) + (dipfini_ii+dipplus_ii)*mtrs(1)           * alpha_sOver2Pi*CF
          res(3) = res(3) + (dipfini_ii+dipplus_ii)*mtrs(2)           * alpha_sOver2Pi*CF

!           print *, res(2)
!           print *, res(3)

!         epcorr=epinv+2d0*dlog(renscale/facscale)
          epcorr=epinv
          AP(1)= 0d0
          AP(2)= CF * (1d0+(1d0-z)**2)/z
          AP(3)= 0d0
          AP(1:3) = AP(1:3) * alpha_sOver2Pi * epcorr * mtrs(1)
          res(1) = res(1) + (AP(1)-AP(3))
          res(2) = res(2) + (AP(2)+AP(3))

          AP(1)= 0d0
          AP(2)= CF * (1d0+(1d0-z)**2)/z
          AP(3)= 0d0
          AP(1:3) = AP(1:3) * alpha_sOver2Pi * epcorr * mtrs(2)
          res(1) = res(1) + (AP(1)-AP(3))
          res(3) = res(3) + (AP(2)+AP(3))

!           print *, res(2)
!           print *, res(3)
!           pause


       endif
  endif  ! for not same/ anti-same flavors



  if (abs(in1).eq.abs(in2)) then


!------------------------------------------------------------------
       if ((in1.gt.0).and.(in2.gt.0)) then    !  Q Q

          call qqq_ttqqqq_dip_int(p,MomDK,mtrs)
          mtrs(1:2) = mtrs(1:2) * 0.5d0 ! symmetry factor
          res(1) = res(1) + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2)) * alpha_sOver2Pi*CF
          res(2) = res(2) + (dipfini_ii+dipplus_ii)*mtrs(1)           * alpha_sOver2Pi*CF
          res(3) = res(3) + (dipfini_ii+dipplus_ii)*mtrs(2)           * alpha_sOver2Pi*CF

!           print *, res(2)
!           print *, res(3)

!         epcorr=epinv+2d0*dlog(renscale/facscale)
          epcorr=epinv
          AP(1)= 0d0
          AP(2)= CF * (1d0+(1d0-z)**2)/z
          AP(3)= 0d0
          AP(1:3) = AP(1:3) * alpha_sOver2Pi * epcorr * mtrs(1)
          res(1) = res(1) + (AP(1)-AP(3))
          res(2) = res(2) + (AP(2)+AP(3))

          AP(1)= 0d0
          AP(2)= CF * (1d0+(1d0-z)**2)/z
          AP(3)= 0d0
          AP(1:3) = AP(1:3) * alpha_sOver2Pi * epcorr * mtrs(2)
          res(1) = res(1) + (AP(1)-AP(3))
          res(3) = res(3) + (AP(2)+AP(3))

!           print *, res(2)
!           print *, res(3)
!           pause

          res(1:3) = two*res(1:3)   ! this factor of two is because there
                                    ! are twice as many dipoles
                                    ! but we also need the symmetry factor

!------------------------------------------------------------------
       elseif ((in1.gt.0).and.(in2.lt.0)) then    !  Q Qbar

          call qqb_ttqqqq_dip_12_int(p,MomDK,mtrs)
          res(1) = res(1) + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2)) * alpha_sOver2Pi*CF
          res(2) = res(2) + (dipfini_ii+dipplus_ii)*mtrs(1)           * alpha_sOver2Pi*CF
          res(3) = res(3) + (dipfini_ii+dipplus_ii)*mtrs(2)           * alpha_sOver2Pi*CF
!       RES(1) = RES(1) +  DIPSOFT_ii*(mtrs(1)+mtrs(2))  ! FOR DELTA-FCT. CHECK


!         epcorr=epinv+2d0*dlog(renscale/facscale)
          epcorr=epinv
          AP(1)= 0d0
          AP(2)= CF * (1d0+(1d0-z)**2)/z
          AP(3)= 0d0
          AP(1:3) = AP(1:3) * alpha_sOver2Pi * epcorr * mtrs(1)
          res(1) = res(1) + (AP(1)-AP(3))
          res(2) = res(2) + (AP(2)+AP(3))

          AP(1)= 0d0
          AP(2)= CF * (1d0+(1d0-z)**2)/z
          AP(3)= 0d0
          AP(1:3) = AP(1:3) * alpha_sOver2Pi * epcorr * mtrs(2)
          res(1) = res(1) + (AP(1)-AP(3))
          res(3) = res(3) + (AP(2)+AP(3))

          call qqb_ttqqqq_dip_11_int(p,MomDK,mtrs)   ! mtrs(2)=0

          dipsoft_fi = fi_qg(zero,zero,p,3,5,z,1)    *nf_light
          dipfini_fi = fi_qg(zero,zero,p,3,5,z,2)    *nf_light
          dipplus_fi = fi_qg(zero,zero,p,3,5,z,3)    *nf_light
          res(1) = res(1) + (dipsoft_fi-dipplus_fi)*(mtrs(1)+mtrs(2)) * alpha_sOver2Pi*0.5d0
          res(2) = res(2) + (dipfini_fi+dipplus_fi)*mtrs(1)           * alpha_sOver2Pi*0.5d0
          res(3) = res(3) + (dipfini_fi+dipplus_fi)*mtrs(2)           * alpha_sOver2Pi*0.5d0
!       RES(1) = RES(1) +  DIPSOFT_fi*(mtrs(1)+mtrs(2))  ! FOR DELTA-FCT. CHECK
!       print *, "eval only delta part for nf only"

!------------------------------------------------------------------
       elseif ((in1.lt.0).and.(in2.gt.0)) then    !  Qbar Q

          p1 = p
          p1(:,3) = p(:,4)
          p1(:,4) = p(:,3)
          call qqb_ttqqqq_dip_12_int(p1,MomDK,mtrs)
          res(1) = res(1) + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2)) * alpha_sOver2Pi*CF
          res(2) = res(2) + (dipfini_ii+dipplus_ii)*mtrs(1)           * alpha_sOver2Pi*CF
          res(3) = res(3) + (dipfini_ii+dipplus_ii)*mtrs(2)           * alpha_sOver2Pi*CF


!         epcorr=epinv+2d0*dlog(renscale/facscale)
          epcorr=epinv
          AP(1)= 0d0
          AP(2)= CF * (1d0+(1d0-z)**2)/z
          AP(3)= 0d0
          AP(1:3) = AP(1:3) * alpha_sOver2Pi * epcorr * mtrs(1)
          res(1) = res(1) + (AP(1)-AP(3))
          res(2) = res(2) + (AP(2)+AP(3))

          AP(1)= 0d0
          AP(2)= CF * (1d0+(1d0-z)**2)/z
          AP(3)= 0d0
          AP(1:3) = AP(1:3) * alpha_sOver2Pi * epcorr * mtrs(2)
          res(1) = res(1) + (AP(1)-AP(3))
          res(3) = res(3) + (AP(2)+AP(3))

          call qqb_ttqqqq_dip_11_int(p1,MomDK,mtrs)
          dipsoft_fi = fi_qg(zero,zero,p,4,5,z,1)    *nf_light
          dipfini_fi = fi_qg(zero,zero,p,4,5,z,2)    *nf_light
          dipplus_fi = fi_qg(zero,zero,p,4,5,z,3)    *nf_light
          res(1) = res(1) + (dipsoft_fi-dipplus_fi)*(mtrs(1)+mtrs(2)) * alpha_sOver2Pi*0.5d0
          res(2) = res(2) + (dipfini_fi+dipplus_fi)*mtrs(1)           * alpha_sOver2Pi*0.5d0
          res(3) = res(3) + (dipfini_fi+dipplus_fi)*mtrs(2)           * alpha_sOver2Pi*0.5d0

!------------------------------------------------------------------
       elseif ((in1.lt.0).and.(in2.lt.0)) then    !  Qbar Qbar

          call qbb_ttqqqq_dip_int(p,MomDK,mtrs)
          mtrs(1:2) = mtrs(1:2) * 0.5d0 ! symmetry factor
          res(1) = res(1) + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2)) * alpha_sOver2Pi*CF
          res(2) = res(2) + (dipfini_ii+dipplus_ii)*mtrs(1)           * alpha_sOver2Pi*CF
          res(3) = res(3) + (dipfini_ii+dipplus_ii)*mtrs(2)           * alpha_sOver2Pi*CF

!           print *, res(2)
!           print *, res(3)

!         epcorr=epinv+2d0*dlog(renscale/facscale)
          epcorr=epinv
          AP(1)= 0d0
          AP(2)= CF * (1d0+(1d0-z)**2)/z
          AP(3)= 0d0
          AP(1:3) = AP(1:3) * alpha_sOver2Pi * epcorr * mtrs(1)
          res(1) = res(1) + (AP(1)-AP(3))
          res(2) = res(2) + (AP(2)+AP(3))

          AP(1)= 0d0
          AP(2)= CF * (1d0+(1d0-z)**2)/z
          AP(3)= 0d0
          AP(1:3) = AP(1:3) * alpha_sOver2Pi * epcorr * mtrs(2)
          res(1) = res(1) + (AP(1)-AP(3))
          res(3) = res(3) + (AP(2)+AP(3))

!           print *, res(2)
!           print *, res(3)
!           pause

          res(1:3) = two*res(1:3)   ! this factor of two is because there
                                    ! are twice as many dipoles
                                    ! but we also need the symmetry factor
       endif
  endif  ! for not same/ anti-same flavors


  end subroutine





  subroutine qqb_ttqqqq_dip_12_int(p,MomDK,res)
  real(dp), intent(in) :: p(4,5), MomDK(1:4,1:6)
  real(dp), intent(out) :: res(2)
  complex(dp) :: cres
  type(Particle), save :: ExtParticles1(1:5),ExtParticles2(1:5),ExtParticles3(1:5)
  type(TreeProcess), save :: TreeAmpsDip(1:4),TreeAmpsDipg(1:4),TreeAmpsDipq(1:6)
  complex(dp) :: Bm(1:4),Bq(1:6)
  complex(dp) :: AM(1:4,1:4)
  integer :: i1,i2,i3,i4,hel(5),i,j,i5,i6,emi,n, in1, in2,iTree
  real(dp) :: dipsoft, dipfini, dipplus
  real(dp) :: z, C(4,4), fx1(-5:5), fx2(-5:5),fx1z(-5:5),fx2z(-5:5)
  real(dp) :: Cq(6,6)
  real(dp) :: mtrs
  real(dp) :: q(4,5)
  real(dp) :: pd1(4,5), pd2(4,5)
  integer :: Nmax(5)
  character :: diptype*2
  logical, save :: first_time = .true.


      res = zero
      cres = (zero,zero)

!---- there are two dipoles here
       pd1 = p
       pd1(:,3) = p(:,5)
       pd1(:,5) = p(:,3)

       pd2 = p
       pd2(:,4) = p(:,5)
       pd2(:,5) = p(:,4)

      if (first_time) then
  !      init external particles
        call InitProcess_TbTQbQG(ExtParticles1(1:5))
  !      init tree processes for 0-> tb t bq q g g
        call InitTrees(4,1,4,TreeAmpsDip)

        TreeAmpsDip(1)%PartRef(1:5) = (/1,2,3,4,5/)
        TreeAmpsDip(2)%PartRef(1:5) = (/1,5,2,3,4/)
        TreeAmpsDip(3)%PartRef(1:5) = (/1,2,5,3,4/)
        TreeAmpsDip(4)%PartRef(1:5) = (/1,2,3,5,4/)

        do iTree=1,4
        call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles1(1:5))
        enddo
        first_time = .false.
      endif


      do n = 1,2

         if (n.eq.1) then
            q = pd1
         elseif(n.eq.2) then
            q = pd2
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

      call GenerateEventttqqg2(q,MomDK,hel,ExtParticles1(1:5))

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


        call cc_ttqqg(C)

        mtrs = zero

      do i=1,4
      do j=1,4
    mtrs = mtrs + C(i,j)*real(AM(i,j),dp)
      enddo
      enddo


      res(n) = mtrs * alpha_s4Pi**3 /4d0/3d0/8d0  ! spin & color avg., color convention adjustment

       enddo


end subroutine qqb_ttqqqq_dip_12_int



  subroutine qqb_ttqqqq_dip_11_int(p,MomDK,res)
  real(dp), intent(in) :: p(4,5),  MomDK(1:4,1:6)
  real(dp), intent(out) :: res(2)  !but in reality only one dip here
  complex(dp) :: cres
  type(Particle), save :: ExtParticles1(1:5),ExtParticles2(1:5),ExtParticles3(1:5)
  type(TreeProcess), save :: TreeAmpsDip(1:4),TreeAmpsDipg(1:4),TreeAmpsDipq(1:6)
  complex(dp) :: Bm(1:4),Bq(1:6)
  complex(dp) :: AM(1:4,1:4)
  integer :: i1,i2,i3,i4,hel(5),i,j,i5,i6,emi,n, in1, in2,iTree
  real(dp) :: dipsoft, dipfini, dipplus
  real(dp) :: z, C(4,4), fx1(-5:5), fx2(-5:5),fx1z(-5:5),fx2z(-5:5)
  real(dp) :: Cq(6,6)
  real(dp) :: mtrs
  integer :: Nmax(5)
  real(dp) :: q(4,5)
  real(dp) :: pd1(4,5), pd2(4,5)
  character :: diptype*2
  logical, save :: first_time = .true.

      res = zero
      cres = (zero,zero)

!---- there is only one dipole here

       pd1 = p

      if (first_time) then
        !      init external particles
              call InitProcess_TbTQbQG(ExtParticles1(1:5))
        !      init tree processes for 0-> tb t bq q g g
              call InitTrees(4,1,4,TreeAmpsDip)

              TreeAmpsDip(1)%PartRef(1:5) = (/1,2,3,4,5/)
              TreeAmpsDip(2)%PartRef(1:5) = (/1,5,2,3,4/)
              TreeAmpsDip(3)%PartRef(1:5) = (/1,2,5,3,4/)
              TreeAmpsDip(4)%PartRef(1:5) = (/1,2,3,5,4/)

              do iTree=1,4
              call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles1(1:5))
              enddo
              first_time = .false.
      endif

      do n = 1,1

         q = pd1

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

      call GenerateEventttqqg2(q,MomDK(1:4,1:6),hel,ExtParticles1(1:5))

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


        call cc_ttqqg(C)

        mtrs = zero

      do i=1,4
      do j=1,4
          mtrs = mtrs + C(i,j)*real(AM(i,j),dp)
      enddo
      enddo


      res(n) = mtrs * alpha_s4Pi**3 /4d0/9d0  ! spin & color avg., color convention adjustment

       enddo


end subroutine qqb_ttqqqq_dip_11_int




  subroutine qqq_ttqqqq_dip_int(p,MomDK,res)
  real(dp), intent(in) :: p(4,5), MomDK(1:4,1:6)
  real(dp), intent(out) :: res(2)
  complex(dp) :: cres
  type(Particle), save :: ExtParticles1(1:5),ExtParticles2(1:5),ExtParticles3(1:5)
  type(TreeProcess), save :: TreeAmpsDip(1:4),TreeAmpsDipg(1:4),TreeAmpsDipq(1:6)
  complex(dp) :: Bm(1:4),Bq(1:6)
  complex(dp) :: AM(1:4,1:4)
  integer :: i1,i2,i3,i4,hel(5),i,j,i5,i6,emi,n, in1, in2,iTree
  real(dp) :: dipsoft, dipfini, dipplus
  real(dp) :: z, C(4,4), fx1(-5:5), fx2(-5:5),fx1z(-5:5),fx2z(-5:5)
  real(dp) :: Cq(6,6)
  real(dp) :: mtrs
  real(dp) :: q(4,5)
  integer :: Nmax(5)
  real(dp) :: pd1(4,5), pd2(4,5)
  character :: diptype*2
  logical, save :: first_time = .true.


      res = zero
      cres = (zero,zero)


!---- there are two dipoles here

      pd1 = p
      pd1(:,3) = p(:,4)
      pd1(:,4) = p(:,5)
      pd1(:,5) = p(:,3)


      pd2 = p
      pd2(:,4) = p(:,5)
      pd2(:,5) = p(:,4)



        if (first_time) then
      !      init external particles
            call InitProcess_TbTQbQG(ExtParticles1(1:5))
      !      init tree processes for 0-> tb t bq q g g
            call InitTrees(4,1,4,TreeAmpsDip)

            TreeAmpsDip(1)%PartRef(1:5) = (/1,2,3,4,5/)
            TreeAmpsDip(2)%PartRef(1:5) = (/1,5,2,3,4/)
            TreeAmpsDip(3)%PartRef(1:5) = (/1,2,5,3,4/)
            TreeAmpsDip(4)%PartRef(1:5) = (/1,2,3,5,4/)

            do iTree=1,4
            call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles1(1:5))
            enddo
            first_time = .false.
      endif


      do n = 1,2

         if (n.eq.1) then
            q = pd1
         elseif(n.eq.2) then
            q = pd2
         endif


      AM = (zero,zero) ! first, the ``amplitude matrix'' is set to zero


         Nmax = 1

        if (TopDecays.ge.1) then
            Nmax(1) = -1
            Nmax(2) = -1
        endif



      do i1= -1,Nmax(1),2
         do i2 = -1,Nmax(2),2
            do i3 = -1,1,2
               do i4 = -1,1,2
                  do i5 = -1,1,2

                 hel(1) = i1
                 hel(2) = i2
                 hel(3) = i3
                 hel(4) = i4
                 hel(5) = i5

      call GenerateEventttqqg2(q,MomDK,hel,ExtParticles1(1:5))

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


        call cc_ttqqg(C)

        mtrs = zero

      do i=1,4
      do j=1,4
    mtrs = mtrs + C(i,j)*real(AM(i,j),dp)
      enddo
      enddo


      res(n) = mtrs* alpha_s4Pi**3 /4d0/3d0/8d0


       enddo


end subroutine qqq_ttqqqq_dip_int



subroutine qbb_ttqqqq_dip_int(p,MomDK,res)
  real(dp), intent(in) :: p(4,5), MomDK(1:4,1:6)
  real(dp), intent(out) :: res(2)
  complex(dp) :: cres
  type(Particle) , save :: ExtParticles1(1:5),ExtParticles2(1:5),ExtParticles3(1:5)
  type(TreeProcess), save :: TreeAmpsDip(1:4),TreeAmpsDipg(1:4),TreeAmpsDipq(1:6)
  complex(dp) :: Bm(1:4),Bq(1:6)
  complex(dp) :: AM(1:4,1:4)
  integer :: i1,i2,i3,i4,hel(5),i,j,i5,i6,emi,n, in1, in2,iTree
  real(dp) :: dipsoft, dipfini, dipplus
  real(dp) :: z, C(4,4), fx1(-5:5), fx2(-5:5),fx1z(-5:5),fx2z(-5:5)
  real(dp) :: Cq(6,6)
  real(dp) :: mtrs
  real(dp) :: q(4,5)
  integer :: Nmax(5)
  real(dp) :: pd1(4,5), pd2(4,5)
  character :: diptype*2
  logical, save :: first_time = .true.

      res = zero
      cres = (zero,zero)
!---- there are two dipoles here

      pd1 = p
      pd1(:,3) = p(:,5)
      pd1(:,5) = p(:,3)


      pd2 = p
      pd2(:,3) = p(:,5)
      pd2(:,4) = p(:,3)
      pd2(:,5) = p(:,4)


      if (first_time) then
      !      init external particles
            call InitProcess_TbTQbQG(ExtParticles1(1:5))
      !      init tree processes for 0-> tb t bq q g g
            call InitTrees(4,1,4,TreeAmpsDip)

            TreeAmpsDip(1)%PartRef(1:5) = (/1,2,3,4,5/)
            TreeAmpsDip(2)%PartRef(1:5) = (/1,5,2,3,4/)
            TreeAmpsDip(3)%PartRef(1:5) = (/1,2,5,3,4/)
            TreeAmpsDip(4)%PartRef(1:5) = (/1,2,3,5,4/)


            do iTree=1,4
            call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles1(1:5))
            enddo
            first_time = .false.
      endif


      do n = 1,2

         if (n.eq.1) then
            q = pd1
         elseif(n.eq.2) then
            q = pd2
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

      call GenerateEventttqqg2(q,MomDK,hel,ExtParticles1(1:5))

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


        call cc_ttqqg(C)

        mtrs = zero

      do i=1,4
      do j=1,4
    mtrs = mtrs + C(i,j)*real(AM(i,j),dp)
      enddo
      enddo


      res(n) = mtrs* alpha_s4Pi**3 /4d0/3d0/8d0


       enddo


end subroutine qbb_ttqqqq_dip_int






end module ModSixFermionIntDip
