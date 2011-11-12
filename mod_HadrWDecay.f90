! this is the file with subroutine for wdecay amplitudes
      module ModHadrWDecay
      implicit none
      private

!----- notation for subroutines
      public ::  ampl_w_qbq,ampl_w_qbqg,wdec_trans,ampl_w_qbqffb


integer, parameter,private  :: dp = selected_real_kind(15)
  real(dp), private, parameter :: pi =&
       & 3.141592653589793238462643383279502884197_dp
  real(dp), private, parameter :: twopi =&
       & 6.283185307179586476925286766559005768394_dp
  real(dp), parameter, private :: pisq =&
       & 9.869604401089358618834490999876151135314_dp
  real(dp), parameter, private :: eulergamma =&
       & 0.577215664901532860606512090082402431042_dp
  real(dp), private, parameter :: half  = 0.5_dp, two = 2.0_dp
  real(dp), private, parameter :: zero  = 0.0_dp, one = 1.0_dp
  real(dp), private, parameter :: mone = -1._dp
  real(dp), private, parameter :: three = 3._dp
  real(dp), private, parameter :: four  = 4._dp
  real(dp), private, parameter :: sqrt2 = &
       &1.4142135623730950488016887242096980785696718753769_dp
  real(dp), private, parameter :: msqrt2 = &
       &-1.4142135623730950488016887242096980785696718753769_dp
  real(dp), private, parameter :: sqrt3 = &
       &1.7320508075688772935274463415058723669428_dp
  real(dp), private, parameter :: CF = four/three
  complex(dp), private, parameter :: ci = dcmplx(zero,1.0_dp)
  complex(dp), private, parameter :: cone = ci
  complex(dp), private, parameter :: czero = dcmplx(0.0_dp,0.0_dp)
!------------------------------------- alpha-parameter
! !   real(dp), private, parameter :: al_ff = 1.0_dp
!   real(dp), private, parameter :: epinv2 = 1.0_dp
!   real(dp), private, parameter :: epinv  = 1.0_dp
  real(dp), private, parameter :: tol = 0.00000010_dp
  complex(dp), private, parameter :: wpol(4)=(/czero,cone,czero,cone/)

      contains

! !----- a subroutine for w->qb + q ; all outgoing convention
! !       0-> w(p1) + qbar(p2) + q(p3) ; p1 = -p2 - p3
!       subroutine w_qbq(p,sum)
!       implicit none
!       real(dp), intent(out) ::  sum
!       real(dp), intent(in) :: p(4,2)
!       integer ::  iw
!       complex(dp) :: A(4), caux
!
!       sum = 0.0_dp
!
!
!           call ampl_w_qbq(p,A)              ! tree-level current
!           caux = sc(wpol,A)
!           sum = sum + caux*conjg(caux)      ! iw is the helicity of the W-boson
!
! !   with the normalization that we have, we get 2*MW^2 for sumed, squared,
! !   no color
!       end subroutine w_qbq


!
! !----- a subroutine for w->qb + q ; all outgoing convention
! !       0-> qbar(p2) + q(p3) ; momentum of the w is computed
! !     virtual corrections + integrated dipoles
!
!       subroutine w_qbq_v_intdip(p,sum)
!       implicit none
!       real(dp), intent(out) ::  sum
!       real(dp), intent(in) :: p(4,2)
!       real(dp), parameter :: ep2 = one/1d-2
!       real(dp), parameter :: ep1 = one/1d-1
!       real(dp) :: v_factor,  int_dip_factor
!       integer ::  iw
!       complex(dp) :: A(4),caux
!
!       sum = 0.0_dp
!
!
!           call ampl_w_qbq(p,A)            ! current
!
!           caux = sc(wpol,A)
!           sum = sum + caux*conjg(caux)
!
!
!       v_factor = as/two/pi*CF*(-2.0_dp*ep2-3.0_dp*ep1-8.0_dp+pisq)
!
! !---- first factor ``2'' is the number of dipoles
!
!       int_dip_factor = two*as/two/pi*CF*(one*ep2          &
!      +three/two*ep1+5.0_dp-pisq/two                       &
!      +three/two*(al_ff - one -log(al_ff))-log(al_ff)**2)
!
!       sum = sum*(v_factor+int_dip_factor)
!
!       end subroutine w_qbq_v_intdip



      subroutine ampl_w_qbq(p,A)
      implicit none
      complex(dp), intent(out) ::  A(4)
      real(dp), intent(in) :: p(4,2)
      integer ::  i1, i2, i3
      complex(dp) :: sp2(4),sp3(4),sp3aux(4)
      complex(dp) :: v(4,4)

      v = czero
      v(1,1) = cone
      v(2,2) = cone
      v(3,3) = cone
      v(4,4) = cone

      sp3 = ubar0(dcmplx(p(:,2)),-1)       ! q polarization/ left handed
      sp2 = v0(dcmplx(p(:,1)),1)           ! qbar polarization/left handed

      do i1=1,4
      sp3aux = spb2(sp3,v(i1,:))
      if (i1.eq.1) then
      A(i1) = -psp1(sp3aux,sp2)
      else
      A(i1) =  psp1(sp3aux,sp2)
      endif
      enddo

      end subroutine ampl_w_qbq


! !----- a subroutine for w->qb + q + g ; all outgoing convention
! !       0-> w(p1) + qbar(p2) + q(p3)+g(p4)
!       subroutine w_qbqg(p,sum)
!       implicit none
!       real(dp), intent(out) ::  sum
!       real(dp), intent(in) :: p(4,4)
!       integer ::  iw,ig
!       complex(dp) :: ew(4),sp2(4),sp3(4),eg(4)
!       real(dp) :: p13(4), p34(4)
!       complex(dp) :: A(4),caux
!
!
!       sum = 0.0_dp
!
!           do ig = 1,2
!         call ampl_w_qbqg(p,ig,A)  !this subroutine returns helicity amplitude
!         caux=sc(wpol,A)
!         sum = sum + caux*conjg(caux)
!          enddo
!
!        sum = gs**2*CF*sum
!
!       end subroutine w_qbqg



      subroutine ampl_w_qbqg(p,ig,A)
      implicit none
      complex(dp), intent(out) ::  A(4)
      real(dp), intent(in) :: p(4,3)
      integer, intent(in) :: ig
      complex(dp) ::  ew(4)
      integer ::  i1, i2, i3
      complex(dp) :: sp2(4),sp3(4),eg(4), res1, res2
      complex(dp) :: sp30(4)
      real(dp) :: p13(4), p34(4)
      complex(dp) :: v(4,4)

      v = czero
      v(1,1) = cone
      v(2,2) = cone
      v(3,3) = cone
      v(4,4) = cone

        sp30 = ubar0(dcmplx(p(:,2)),-1)  ! q polarization/ left handed
        sp2 = v0(dcmplx(p(:,1)),1)     ! qbar polarization/left handed
        eg = pol_mless2(dcmplx(p(:,3)),ig,'out')  ! gluon

        do i1=1,4
           if (i1.eq.1) then
           ew = -1d0*v(i1,:)
           else
           ew = v(i1,:)
           endif

            sp3 = spb2(sp30,eg)
            p34 = p(:,2)+p(:,3)
            sp3 = spb2(sp3,dcmplx(p34))/scr(p34,p34)
            sp3 = spb2(sp3,ew)
            res1  = psp1(sp3,sp2)

            sp3 = spb2(sp30,ew)
            p13 = -p(:,3)-p(:,1)
            sp3 = spb2(sp3,dcmplx(p13))/scr(p13,p13)
            sp3 = spb2(sp3,eg)
            res2 =   psp1(sp3,sp2)

            A(i1) = res1 + res2

          enddo

      end subroutine ampl_w_qbqg




      subroutine ampl_w_qbqffb(p,ig,A)
      implicit none
      complex(dp), intent(out) ::  A(4)
      real(dp), intent(in) :: p(4,4)
      integer, intent(in) :: ig
      complex(dp) ::  ew(4)
      integer ::  i1, i2, i3
      complex(dp) :: sp2(4),sp3(4),eg(4), res1, res2,sp4(4),sp5(4)
      complex(dp) :: sp30(4)
      real(dp) :: p13(4), p34(4)
      complex(dp) :: v(4,4)

      v = czero
      v(1,1) = cone
      v(2,2) = cone
      v(3,3) = cone
      v(4,4) = cone

        sp30 = ubar0(dcmplx(p(:,2)),-1)  ! q polarization/ left handed
        sp2 = v0(dcmplx(p(:,1)),1)     ! qbar polarization/left handed

        sp4 = ubar0(dcmplx(p(:,3)),-ig)
        sp5 = v0(dcmplx(p(:,4)),ig)
        eg = vbqq_weyl(4,sp4,sp5)/(2d0*scr(p(:,3),p(:,4)))


        do i1=1,4

           if (i1.eq.1) then
           ew = -1d0*v(i1,:)
           else
           ew = v(i1,:)
           endif

            sp3 = spb2(sp30,eg)
            p34 = p(:,2)+p(:,3)+p(:,4)
            sp3 = spb2(sp3,dcmplx(p34))/scr(p34,p34)
            sp3 = spb2(sp3,ew)
            res1  = psp1(sp3,sp2)


            sp3 = spb2(sp30,ew)
            p13 = -p(:,3)-p(:,4)-p(:,1)
            sp3 = spb2(sp3,dcmplx(p13))/scr(p13,p13)
            sp3 = spb2(sp3,eg)
            res2 =   psp1(sp3,sp2)

            A(i1) = res1 + res2

            enddo

      end subroutine ampl_w_qbqffb







       subroutine wdec_trans(n,p,ptilde,alpha_DKWff,res)! momentum order: p(1:4,i) i=1 bot, i=2,3 qua, i=4 glu
       implicit none
       integer,  intent(in) :: n
       real(dp), intent(in) :: p(4,4)
       real(dp), intent(out) :: ptilde(4,3),res
       real(dp) ::  yijk,zi,zj
       real(dp) :: pi(4),pj(4),pk(4),pij(4),pij2
       real(dp) :: pijt(4), pkt(4),q(1:4,3), paux(4)
       real(dp) :: weight, diag,alpha_DKWff
       integer :: i1,i,j,k
       complex(dp) :: A(4), caux
       integer ::  dip(2,3)
       data dip(1,1)/4/ , dip(1,2)/3/, dip(1,3)/2/
       data dip(2,1)/4/ , dip(2,2)/2/, dip(2,3)/3/

        res = 0.0_dp

!       momentum mapping

         i=dip(n,1)  ! emitted
         j=dip(n,2)  ! emittor
         k=dip(n,3)  ! spectator

        pi=p(:,i)
        pj=p(:,j)
        pk=p(:,k)

        pij = pi+pj
        pij2  = scr(pij,pij)

        yijk = scr(pi,pj)/(scr(pi,pj)+scr(pi,pk)+scr(pj,pk))

        if (yijk.gt.alpha_DKWff) return

        zi = scr(pi,pk)/(scr(pi,pk)+scr(pj,pk))
        zj = 1d0-zi

        pkt = pk/(one-yijk)
        pijt = pi+pj-yijk*pkt

        paux = zi*pi - zj*pj


!---------------------------------------momenta
        q(:,1) = p(:,1)
       if (n.eq.1) then
          q(:,3) = pijt
          q(:,2) = pkt
       elseif(n.eq.2) then
          q(:,3) = pkt
          q(:,2) = pijt
       endif

       ptilde = q   ! tilde momentum


!------------------------------------------------------ weight

       res = (two/(one-zj*(one-yijk))-(one + zj) ) * two*one  /(pij2)

!       res = gs**2*CF**res ---> normalization

      end subroutine wdec_trans





!---- THESE ARE POLARIZATION/OTHER ROUTINES

  ! -- massless vector polarization subroutine
  function pol_mless(p,i,outgoing)
    complex(dp), intent(in)    :: p(4)
    integer, intent(in)          :: i
    logical, intent(in),optional :: outgoing
    ! -------------------------------
    integer :: pol
    real(dp) :: p0,px,py,pz
    real(dp) :: pv,ct,st,cphi,sphi
    complex(dp) :: pol_mless(4)

!^^^IFmp
!    p0=(p(1)+conjg(p(1)))/two
!    px=(p(2)+conjg(p(2)))/two
!    py=(p(3)+conjg(p(3)))/two
!    pz=(p(4)+conjg(p(4)))/two
!^^^ELSE
    p0=real(p(1),dp)
    px=real(p(2),dp)
    py=real(p(3),dp)
    pz=real(p(4),dp)
!^^^END


    pv=sqrt(abs(p0**2))
    ct=pz/pv
    st=sqrt(abs(1.0_dp-ct**2))

    if (st < tol) then
       cphi=1.0_dp
       sphi=0.0_dp
    else
       cphi= px/pv/st
       sphi= py/pv/st
    endif


    ! -- distinguish between positive and negative energies
    if ( p0 > 0.0_dp) then
       pol=i
    else
       pol=-i
    endif

    ! -- take complex conjugate for outgoing;   removed to be in agreement with same routine in mod_Kinematics
!     if (present(outgoing)) then
!        if (outgoing) pol = -pol
!     endif

    pol_mless(1)=czero
    pol_mless(2)=ct*cphi/sqrt2 - ci*pol*sphi/sqrt2
    pol_mless(3)=ct*sphi/sqrt2 + ci*pol*cphi/sqrt2
    pol_mless(4)=-st/sqrt2

  end function pol_mless


  function pol_mless2(p,i,out)
    integer, intent(in)       :: i
    complex(dp), intent(in) :: p(4)
    character(len=*), intent(in):: out
    complex(dp)             :: pol_mless2(4)
    ! -------------------------------------

    if (out == 'out') then
       pol_mless2 = pol_mless(p,i,outgoing=.true.)
    else
       pol_mless2 = pol_mless(p,i,outgoing=.false.)
    endif
  end function pol_mless2


  function pol_dk2mom(plepton,antilepton,i,outgoing)
    integer, intent(in) :: i
    integer :: j
    complex(dp), intent(in) :: plepton(:),antilepton(:)
    logical, intent(in),optional :: outgoing
    complex(dp) :: pol_dk2mom(4),Ub(4),V(4),q(4),qsq


    q=plepton+antilepton
    qsq=q(1)**2-q(2)**2-q(3)**2-q(4)**2

    Ub(:)=ubar0(plepton,i)
    V(:)=v0(antilepton,-i)

    !---Now return in Kirill's notation  1=E,2=px,3=py,4=pz
    !   This is an expression for (-i)/qsq* (-i) Ub(+/-)) Gamma^\mu V(-/+)
    pol_dk2mom(1)=-(Ub(2)*V(4)+V(2)*Ub(4)+Ub(1)*V(3)+V(1)*Ub(3))
    pol_dk2mom(2)=-(-Ub(1)*V(4)+V(1)*Ub(4)-Ub(2)*V(3)+V(2)*Ub(3))
    pol_dk2mom(3)=-ci*(Ub(1)*V(4)+V(1)*Ub(4)-Ub(2)*V(3)-V(2)*Ub(3))
    pol_dk2mom(4)=-(Ub(2)*V(4)-V(2)*Ub(4)-Ub(1)*V(3)+V(1)*Ub(3))

    do j=1,4
       pol_dk2mom(j)=pol_dk2mom(j)/qsq
    enddo

    ! -- do nothing in this case
    if (present(outgoing)) then
       !if (outgoing) pol_dk2mom = conjg(pol_dk2mom)
    endif
  end function pol_dk2mom


   !     ubar spinor, massless

  function ubar0(p,i)
    complex(dp), intent(in) :: p(4)
    integer, intent(in)       :: i
    ! -------------------------------
    complex(dp) :: ubar0(4)
    complex(dp) :: fc, fc2
    real(dp)    :: p0,px,py,pz

!^^^IFmp
!    p0=(p(1)+conjg(p(1)))/two
!    px=(p(2)+conjg(p(2)))/two
!    py=(p(3)+conjg(p(3)))/two
!    pz=(p(4)+conjg(p(4)))/two
!^^^ELSE
    p0=real(p(1),dp)
    px=real(p(2),dp)
    py=real(p(3),dp)
    pz=real(p(4),dp)
!^^^END

    fc2 = p0 + pz
    fc=sqrt(fc2)

    if (abs(fc2).gt. tol) then

       if (i.eq.1) then
          ubar0(1)=czero
          ubar0(2)=czero
          ubar0(3)=fc
          ubar0(4)=(px-ci*py)/fc
       elseif (i.eq.-1) then
          ubar0(1)=(px+ci*py)/fc
          ubar0(2)=-fc
          ubar0(3)=czero
          ubar0(4)=czero
       else
          stop 'ubar0: i out of range'
       endif

    else
       if (i.eq.1) then
          ubar0(1) = czero
          ubar0(2) = czero
          ubar0(3) = czero
          ubar0(4) = sqrt(cone*two*p0)
       elseif (i.eq.-1) then
          ubar0(1) = sqrt(cone*(two*p0))
          ubar0(2) = czero
          ubar0(3) = czero
          ubar0(4) = czero
       else
          stop 'ubar0: i out of range'
       endif
    endif


  end function ubar0



  ! -- v0  spinor, massless
  function v0(p,i)
    complex(dp), intent(in) :: p(4)
    integer, intent(in)       :: i
    ! -------------------------------
    complex(dp) :: v0(4)
    complex(dp) :: fc2, fc
    real(dp)    :: p0,px,py,pz

!^^^IFmp
!    p0=(p(1)+conjg(p(1)))/two
!    px=(p(2)+conjg(p(2)))/two
!    py=(p(3)+conjg(p(3)))/two
!    pz=(p(4)+conjg(p(4)))/two
!^^^ELSE
    p0=real(p(1),dp)
    px=real(p(2),dp)
    py=real(p(3),dp)
    pz=real(p(4),dp)
!^^^END

    fc2 = p0 + pz
    fc=sqrt(fc2)



    if (abs(fc2).gt. tol) then

       if (i.eq.1) then
          v0(1)=czero
          v0(2)=czero
          v0(3)=(px-ci*py)/fc
          v0(4)=-fc
       elseif (i.eq.-1) then
          v0(1)=fc
          v0(2)=(px+ci*py)/fc
          v0(3)=czero
          v0(4)=czero
       else
          stop 'v0: i out of range'
       endif

    else

       if (i.eq.1) then
          v0(1)=czero
          v0(2)=czero
          v0(3)=sqrt(cone*two*p0)
          v0(4)=czero
       elseif (i.eq.-1) then
          v0(1)=czero
          v0(2)=sqrt(cone*two*p0)
          v0(3)=czero
          v0(4)=czero
       else
          stop 'v0: i out of range'
       endif

    endif

  end function v0



        FUNCTION vbqq_Weyl(Dv,sp1,sp2)
        implicit none
        complex(8), intent(in) :: sp1(:), sp2(:)
        integer, intent(in) ::  Dv
        integer :: i
        complex(8) :: vbqq_Weyl(Dv)
        complex(8) :: rr, va(Dv),sp1a(size(sp1))
        real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

            va=(0d0,0d0)
            vbqq_Weyl=(0d0,0d0)

            do i=1,Dv
              if (i.eq.1) then
                va(1)=(1d0,0d0)
              else
                va(i)=(-1d0,0d0)
              endif
              sp1a=spb2_Weyl(sp1,va)

              rr=(0d0,-1d0)/sqrt2*psp1(sp1a,sp2)
              if (i.eq.1) then
                    vbqq_Weyl = vbqq_Weyl + rr*va
                else
                    vbqq_Weyl = vbqq_Weyl - rr*va
              endif
              va(i)=(0d0,0d0)
            enddo

        END FUNCTION




      function spb2_Weyl(sp,v)
      implicit none
      complex(8), intent(in) :: sp(:),v(:)
      complex (8) :: spb2_Weyl(size(sp))
      complex(8) :: x0(4,4),xx(4,4),xy(4,4)
      complex(8) :: xz(4,4),x5(4,4)
      complex(8) :: y1,y2,y3,y4,bp,bm,cp,cm
      integer :: i,i1,i2,i3,Dv,Ds,imax

      Ds = size(sp)

      if (Ds == 4) Dv = 4
      if (Ds == 8) Dv = 6
      if (Ds == 16) Dv = 8

      imax = Ds/4

           do i=1,imax
           i1= 1+4*(i-1)
           i2=i1+3

           y1=sp(i1)
           y2=sp(i1+1)
           y3=sp(i1+2)
           y4=sp(i1+3)

           x0(1,i)=y3
           x0(2,i)=y4
           x0(3,i)=y1
           x0(4,i)=y2

           xx(1,i) = y4
           xx(2,i) = y3
           xx(3,i) = -y2
           xx(4,i) = -y1

           xy(1,i)=(0d0,1d0)*y4
           xy(2,i)=-(0d0,1d0)*y3
           xy(3,i)=-(0d0,1d0)*y2
           xy(4,i)=(0d0,1d0)*y1

           xz(1,i)=y3
           xz(2,i)=-y4
           xz(3,i)=-y1
           xz(4,i)=y2

           x5(1,i)=y1
           x5(2,i)=y2
           x5(3,i)=-y3
           x5(4,i)=-y4

           enddo

           if (Dv.eq.4) then

           do i=1,4

           spb2_Weyl(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)

           enddo

           endif

           if (Dv.eq.6) then
           bp = (v(5)+(0d0,1d0)*v(6))
           bm=(v(5)-(0d0,1d0)*v(6))

           do i=1,4

           spb2_Weyl(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)+bm*x5(i,2)

            i1 = i+4

            spb2_Weyl(i1)= v(1)*x0(i,2)-v(2)*xx(i,2)-v(3)*xy(i,2)-v(4)*xz(i,2) -bp*x5(i,1)


            enddo

           endif

           if (Dv.eq.8) then
           bp=(v(5)+(0d0,1d0)*v(6))
           bm=(v(5)-(0d0,1d0)*v(6))
           cp=(v(7)+(0d0,1d0)*v(8))
           cm=(v(7)-(0d0,1d0)*v(8))

           do i=1,4

       spb2_Weyl(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)+bm*x5(i,2)-cm*x5(i,3)

            i1 = i+4

       spb2_Weyl(i1) = v(1)*x0(i,2)-v(2)*xx(i,2) -v(3)*xy(i,2)-v(4)*xz(i,2) -bp*x5(i,1)+cm*x5(i,4)

             i2 = i1+4

       spb2_Weyl(i2)=v(1)*x0(i,3)-v(2)*xx(i,3)  -v(3)*xy(i,3)-v(4)*xz(i,3)+bm*x5(i,4)+cp*x5(i,1)

              i3=i2+4

       spb2_Weyl(i3)=v(1)*x0(i,4)-v(2)*xx(i,4)  -v(3)*xy(i,4)-v(4)*xz(i,4) -bp*x5(i,3)-cp*x5(i,2)

              enddo

              endif

               end function




!--------massive vector boson polarization routine

  function pol_mass(p,m,i,outgoing)
    integer, intent(in)       :: i
    complex(dp), intent(in) :: p(4)
    real(dp),  intent(in)   :: m
    logical, intent(in),optional :: outgoing
    complex(dp)             :: pol_mass(4)
    ! -------------------------------------
    real(dp) :: p0,px,py,pz, pv
    real(dp) :: ct,st,cphi,sphi
    integer :: pol

!^^^IFmp
!    p0=(p(1)+conjg(p(1)))/two
!    px=(p(2)+conjg(p(2)))/two
!    py=(p(3)+conjg(p(3)))/two
!    pz=(p(4)+conjg(p(4)))/two
!^^^ELSE
    p0=real(p(1),dp)
    px=real(p(2),dp)
    py=real(p(3),dp)
    pz=real(p(4),dp)
!^^^END


    ! i=0 is longitudinal polarization

    ! -- distinguish between positive and negative energies
    if ( p0 > zero) then
       pol=i
    else
       pol=-i
    endif

    ! -- take complex conjugate for outgoing
    if (present(outgoing)) then
       if (outgoing) pol = -pol
    endif



    pv= sqrt(abs(p0**2 - m**2))

    if (pv.gt.1d-8) then  ! for ``moving W'

    ct= pz/pv
    st= sqrt(abs(one-ct**2))

    if (st < tol) then
       cphi=one
       sphi=zero
    else
       cphi= px/pv/st
       sphi= py/pv/st
    endif


    if(pol == -1.or.pol == 1) then
       pol_mass(1)=czero
       pol_mass(2)=ct*cphi/sqrt2 - ci*pol*sphi/sqrt2
       pol_mass(3)=ct*sphi/sqrt2 + ci*pol*cphi/sqrt2
       pol_mass(4)=-st/sqrt2
    elseif (pol == 0) then
       pol_mass(1)= pv/m
       pol_mass(2)= p0/m/pv*px
       pol_mass(3)= p0/m/pv*py
       pol_mass(4)= p0/m/pv*pz
       else
          stop 'pol_mass: pol out of range'
       endif

       else

          pol_mass = czero

          if (pol == -1.or.pol == 1) then

        pol_mass(2)=one/sqrt2
        pol_mass(3)=+ ci*pol/sqrt2

          elseif(pol == 0) then

             pol_mass(4) = (1.0_dp,zero)

          endif

       endif

  end function pol_mass

  function pol_mass2(p,m,i,out)
    integer, intent(in)       :: i
    complex(dp), intent(in) :: p(4)
    real(dp),  intent(in)   :: m
    character(len=*), intent(in):: out
    complex(dp)             :: pol_mass2(4)
    ! -------------------------------------

    if (out == 'out') then
       pol_mass2 = pol_mass(p,m,i,outgoing=.true.)
    else
       pol_mass2 = pol_mass(p,m,i,outgoing=.false.)
    endif
  end function pol_mass2



    function sc(p1,p2)
      complex(dp), intent(in) :: p1(:)
      complex(dp), intent(in) :: p2(:)
      complex(dp)             :: sc
      integer :: sizemin

      sizemin=min(size(p1),size(p2))

      sc = p1(1)*p2(1)
      sc = sc - sum(p1(2:sizemin)*p2(2:sizemin))

    end function sc

    double precision function scr(p1,p2)   !scalar product of real vectors
    real(dp), intent(in) :: p1(4), p2(4)
    scr = p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)-p1(4)*p2(4)
    end function scr



    function spb2(sp,v)
      complex(dp), intent(in) :: sp(:),v(:)
      complex(dp) :: spb2(size(sp))
      complex(dp) :: x0(4,4),xx(4,4),xy(4,4)
      complex(dp) :: xz(4,4),x5(4,4)
      complex(dp) :: y1,y2,y3,y4,bp,bm,cp,cm
      integer :: i,i1,i2,i3,Dv,Ds,imax

      Ds = size(sp)

      if (Ds == 4) then
         Dv = 4
      elseif (Ds == 8) then
         Dv = 6
      elseif (Ds == 16) then
         Dv = 8
      else
         stop 'spb2:Dv not allowed'
      endif

      imax = Ds/4

      do i=1,imax
         i1= 1+4*(i-1)
         i2=i1+3

         y1=sp(i1)
         y2=sp(i1+1)
         y3=sp(i1+2)
         y4=sp(i1+3)

         x0(1,i)=y3
         x0(2,i)=y4
         x0(3,i)=y1
         x0(4,i)=y2

         xx(1,i) = y4
         xx(2,i) = y3
         xx(3,i) = -y2
         xx(4,i) = -y1

         xy(1,i)=ci*y4
         xy(2,i)=-ci*y3
         xy(3,i)=-ci*y2
         xy(4,i)=ci*y1

         xz(1,i)=y3
         xz(2,i)=-y4
         xz(3,i)=-y1
         xz(4,i)=y2

         x5(1,i)=y1
         x5(2,i)=y2
         x5(3,i)=-y3
         x5(4,i)=-y4
      enddo

      if (Dv.eq.4) then

         do i=1,4
            spb2(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)
         enddo

      elseif (Dv.eq.6) then
         bp = (v(5)+ci*v(6))
         bm=(v(5)-ci*v(6))

         do i=1,4

            spb2(i)=v(1)*x0(i,1)-v(2)*xx(i,1) &
                 &          -v(3)*xy(i,1)-v(4)*xz(i,1)+bm*x5(i,2)

            i1 = i+4
            spb2(i1)= v(1)*x0(i,2)-v(2)*xx(i,2) &
                 &            -v(3)*xy(i,2)-v(4)*xz(i,2)-bp*x5(i,1)
         enddo
      elseif (Dv.eq.8) then
         bp=(v(5)+ci*v(6))
         bm=(v(5)-ci*v(6))
         cp=(v(7)+ci*v(8))
         cm=(v(7)-ci*v(8))

         do i=1,4

            spb2(i)=v(1)*x0(i,1)-v(2)*xx(i,1) &
            &          -v(3)*xy(i,1)-v(4)*xz(i,1) &
            &         +bm*x5(i,2)-cm*x5(i,3)

            i1 = i+4

            spb2(i1) = v(1)*x0(i,2)-v(2)*xx(i,2) &
            &            -v(3)*xy(i,2)-v(4)*xz(i,2) &
            &            -bp*x5(i,1)+cm*x5(i,4)

            i2 = i1+4

            spb2(i2)=v(1)*x0(i,3)-v(2)*xx(i,3) &
            &             -v(3)*xy(i,3)-v(4)*xz(i,3) &
            &             +bm*x5(i,4)+cp*x5(i,1)

            i3=i2+4

            spb2(i3)=v(1)*x0(i,4)-v(2)*xx(i,4) &
            &             -v(3)*xy(i,4)-v(4)*xz(i,4) &
            &             -bp*x5(i,3)-cp*x5(i,2)

         enddo
      else
         stop 'spb2: Dv out of bound'
      endif

    end function spb2



    function spi2(v,sp)
      complex(dp), intent(in) :: sp(:),v(:)
      complex(dp) :: spi2(size(sp))
      complex(dp) :: x0(4,4),xx(4,4),xy(4,4)
      complex(dp) :: xz(4,4),x5(4,4)
      complex(dp) ::  y1,y2,y3,y4,bp,bm,cp,cm
      integer :: i,i1,i2,i3,imax,Dv,Ds

      Ds = size(sp)

      if (Ds == 4) then
         Dv = 4
      elseif (Ds == 8) then
         Dv = 6
      elseif (Ds == 16) then
         Dv = 8
      else
         stop 'spi2:Dv not allowed'
      endif



      imax = Ds/4

      do i=1,imax
         i1= 1+4*(i-1)
         i2=i1+3

         y1=sp(i1)
         y2=sp(i1+1)
         y3=sp(i1+2)
         y4=sp(i1+3)

         x0(1,i)=y3
         x0(2,i)=y4
         x0(3,i)=y1
         x0(4,i)=y2


         xx(1,i) = -y4
         xx(2,i) = -y3
         xx(3,i) = y2
         xx(4,i) = y1


         xy(1,i)=ci*y4
         xy(2,i)=-ci*y3
         xy(3,i)=-ci*y2
         xy(4,i)=ci*y1

         xz(1,i)=-y3
         xz(2,i)=y4
         xz(3,i)=y1
         xz(4,i)=-y2

         x5(1,i)=y1
         x5(2,i)=y2
         x5(3,i)=-y3
         x5(4,i)=-y4

      enddo

      if(Dv.eq.4) then

         do i=1,4

            spi2(i)=v(1)*x0(i,1)-v(2)*xx(i,1) &
            &           -v(3)*xy(i,1)-v(4)*xz(i,1)
         enddo

      elseif (Dv.eq.6) then
         bp = (v(5)+ci*v(6))
         bm=(v(5)-ci*v(6))


         do i=1,4

            spi2(i)=v(1)*x0(i,1)-v(2)*xx(i,1) &
            &           -v(3)*xy(i,1)-v(4)*xz(i,1) &
            &           -bp*x5(i,2)

            i1=i+4

            spi2(i1)=v(1)*x0(i,2)-v(2)*xx(i,2) &
            &             -v(3)*xy(i,2)-v(4)*xz(i,2) &
            &             +bm*x5(i,1)

         enddo

      elseif (Dv.eq.8) then

         bp = (v(5)+ci*v(6))
         bm=(v(5)-ci*v(6))
         cp=(v(7)+ci*v(8))
         cm=(v(7)-ci*v(8))

         do i=1,4

            spi2(i)=v(1)*x0(i,1)-v(2)*xx(i,1)&
            &           -v(3)*xy(i,1)-v(4)*xz(i,1)&
            &           -bp*x5(i,2)+ cp*x5(i,3)

            i1=i+4

            spi2(i1)=v(1)*x0(i,2)-v(2)*xx(i,2)&
            &             -v(3)*xy(i,2)-v(4)*xz(i,2)&
            &             +bm*x5(i,1)-cp*x5(i,4)

            i2=i1+4

            spi2(i2)=v(1)*x0(i,3)-v(2)*xx(i,3)&
            &           -v(3)*xy(i,3)-v(4)*xz(i,3)&
            &          -bp*x5(i,4)-cm*x5(i,1)

            i3=i2+4

            spi2(i3)=v(1)*x0(i,4)-v(2)*xx(i,4)&
            &           -v(3)*xy(i,4)-v(4)*xz(i,4)&
            &           +bm*x5(i,3)+cm*x5(i,2)


         enddo

      else
         stop 'spi2: Dv out of bounds'
      end if

    end function spi2


    function  psp1(sp1,sp2)
      complex(dp), intent(in) :: sp1(:)
      complex(dp), intent(in) :: sp2(:)
      complex(dp) :: psp1

      psp1 = sum(sp1(1:)*sp2(1:))

    end function psp1





       end module


