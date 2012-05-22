! this module is the same as mod_IntDipoles (in mod_IntDipoles.f90) apart from the size of the momentum arrays p(:,:) !!
! CHANGES:
! log(x) --> dlog(x)
! line 367: 1/2 --> 1d0/2d0
! ff_gg, ff_qq replaced
! closed bracket in front of 11d0/6d0 in line 359
! replaced kappa by global parameter kappa_ff
! sqrt(x) --> dsqrt(x)


module ModIntDipoles
use ModParameters
use ModMisc
implicit none

public :: ff_qq, ff_gg, fi_qq, fi_gg, fi_qg,if_gg, ii_gg
public :: if_qq, if_qg, if_gq, ii_qq, ii_qg, ii_gq
public :: ii_mgg, if_mgg, fi_mqq, ff_mqq, ii_mgq, ii_mqg, ii_mqq

integer, parameter,private  :: dp = selected_real_kind(15)
real(dp),parameter,private :: pisqo6=pi**2/6d0
!real(dp),private :: MuRen**2   = MuRen**2  ! assumes that MuRes==MuFac!!
real(dp),public,parameter :: epinv2 = 0d0    ! R(eps2,eps)=A/eps2 + B/eps+C  --> A=R(2,1)-R(1,1),  B=R(0,2)-R(0,1)
real(dp),public,parameter :: epinv  = 0d0
character,private,parameter :: scheme*(4)='fdh'  ! this parameter should match the scheme of the virtual corrections
!  character,parameter :: scheme*(4)='fdh', 'tH-V'


contains




! REFERENCES:
! CDST:  hep-ph/0201036
! CS:    hep-ph/9605323
! BCPPW: 0907.4723[hep-ph]
! CET:   hep-ph/0408158
! CT:    hep-ph/0506289
! MCFM:  see ../Need/dipoles_mass.f




!---------------------------------------------------------------------------
!    final-final dipoles
!---------------------------------------------------------------------------
  double precision function  ff_qq(mq1,mq2,p,i,j,z,s)
   real(dp), intent(in)  :: mq1,mq2,z
   real(dp), intent(in)  :: p(:,:)
   integer, intent(in) :: i,j,s
   real(dp) :: mu1, mu2, mu1s, mu2s, q2, Ieik, Icoll,L
   real(dp) ::  r, vijk, r1, r2
   real(dp) :: x,yp,xm,xp,a,b,c,d, Ieik_add, Icoll_add
   integer :: check

   Ieik = zero
   Icoll = zero
   Ieik_add = zero
   Icoll_add = zero

   if (s.eq.1) then

   check = 0

   q2 = two*scr(p(:,i),p(:,j))+mq1**2 + mq2**2
   L = dlog(q2/MuRen**2)


   mu1 = mq1/dsqrt(q2)
   mu2 =  mq2/dsqrt(q2)
   mu1s = mu1**2
   mu2s = mu2**2

   vijk = dsqrt(vl(1.0_dp,mu1s,mu2s))/(1.0_dp - mu1s - mu2s)
   r = dsqrt((one - vijk)/(one + vijk))
   r1 = dsqrt((one - vijk + &
   two*mu1s/(one - mu1s - mu2s))/(one + vijk + two*mu1s/(one-mu1s-mu2s)))
   r2 = dsqrt((one - vijk + &
   two*mu2s/(one - mu1s - mu2s))/(one + vijk + two*mu2s/(one-mu1s-mu2s)))


   if (mq1.eq.mq2.and.mq1.gt.zero) then ! massive case, no scheme dependence

      check = 1

   Ieik = 0.5_dp/vijk*dlog(r)*(epinv-L)      ! CDST eq.(5.34)

   Icoll = epinv -L      ! CDST eq.(5.35)

   Ieik = Ieik + one/vijk*(-dlog(r)*dlog(one - (mu1 + mu2)**2) &      ! CDST eq.(5.34)
   - 0.5_dp*dlog(r1)**2   - 0.5_dp*dlog(r2)**2 &
   + pisq/6.0_dp+2.0_dp*DLi2(-r) -2.0_dp*DLi2(one - r) &
   -0.5_dp*DLi2(one-r1**2) -0.5_dp*DLi2(one-r2**2) )

   Icoll = Icoll + dlog(mu1)-two-two*dlog((one-mu2)**2-mu1s)+dlog(one-mu2) &      ! CDST eq.(5.35)
   - two*mu1s/(one-mu1s-mu2s)*dlog(mu1/(one-mu2)) &
   +5.0_dp - mu2/(one-mu2) - two*mu2*(one-two*mu2)/(one-mu1s-mu2s)


  if (alpha_ff.ne.1d0) then   ! relative to Czakon, we redefine alpha_ff(CZ) = alpha_ff*yp

   yp= one - two*mu2*(one - mu2)/(one - mu1s-mu2s)

  x = yp*(one - alpha_ff) + dsqrt(yp*(one - alpha_ff)*(one/yp - yp*alpha_ff &
 + 4.0_dp*mu1s*mu2s/(mu1s - (one - mu2)**2)/(one - mu1s-mu2s)))
!   xp = ((one - mu2)**2 - mu1s)/(one - mu1s - mu2s) +vijk       ! HERE WAS A BUG: for some reasons values with xp < x appeared
!   xm = ((one - mu2)**2 - mu1s)/(one - mu1s - mu2s) -vijk
  xp = 1d0+(2d0*mu2s-2d0*mu2+dsqrt(vl(1.0_dp,mu1s,mu2s)))/(one - mu1s - mu2s)
  xm = 1d0+(2d0*mu2s-2d0*mu2-dsqrt(vl(1.0_dp,mu1s,mu2s)))/(one - mu1s - mu2s)

  a = two*mu2/(one - mu1s-mu2s)
  b = two*(one - mu2)/(one - mu1s - mu2s)
  c = two*mu2*(one - mu2)/(one - mu1s - mu2s)
  d = one/two*(one - mu1s - mu2s)

  if( xp.lt.x .or. xm.gt.x ) then
    print *, "error xp<x or xm>x",x,xm,xp
    ff_qq = 0d0
    return
  endif

 Ieik_add = 1d0/vijk*(-DLi2((a+x)/(a+xp)) + DLi2(a/(a+xp)) &      ! BCPPW eq.(A.9)
  +DLi2((xp-x)/(xp-b)) - DLi2(xp/(xp-b)) + DLi2((c+x)/(c+xp)) &
 - DLi2(c/(c+xp)) + DLi2((xm - x)/(xm +a)) - DLi2(xm/(xm+a)) &
 -DLi2((b-x)/(b-xm)) + DLi2(b/(b-xm))-DLi2((xm-x)/(xm+c)) &
 + DLi2(xm/(xm+c)) + DLi2((b-x)/(b+a)) - DLi2(b/(b+a)) &
 -DLi2((c+x)/(c-a)) + DLi2(c/(c-a))  &
 + dlog(c+x)*dlog((a-c)*(xp-x)/(a+x)/(c+xp)) - dlog(c)*dlog((a-c)*xp/a/(c+xp)) &
 +dlog(b-x)*dlog((a+x)*(xm-b)/(a+b)/(xm-x)) - dlog(b)*dlog((a*(xm-b)/(a+b)/xm)) &
 - dlog((a+x)*(b-xp))*dlog(xp-x)+dlog(a*(b-xp))*dlog(xp) &
 +dlog(d)*dlog((a+x)*xp*xm/a/(xp-x)/(xm-x))  &
 + dlog((xm-x)/xm)*dlog((c+xm)/(a+xm))  &
 + 1d0/2d0*dlog((a+x)/a)*dlog(a*(a+x)*(a+xp)**2))

 Icoll_add = three/two*(one + alpha_ff*yp) + one/(one - mu2) &      ! BCPPW eq.(A.20)
   -two*(two-two*mu1s-mu2)/(one - mu1s-mu2s)  &
   + (one-alpha_ff*yp)*mu1s/two/(mu1s+alpha_ff*yp*(one - mu1s - mu2s)) &
   -two*dlog(alpha_ff*yp*(one-mu1s-mu2s)/((one-mu2)**2-mu1s)) &
   +(one + mu1s-mu2s)/two/(one - mu1s - mu2s)*dlog((mu1s  &
   +alpha_ff*yp*(one - mu1s-mu2s))/(one-mu2)**2)


 Ieik = Ieik+Ieik_add
 Icoll = Icoll+Icoll_add

  endif

   elseif ((mq1.ne.zero).and.(mq2.eq.zero)) then ! massive case,no scheme dependence

      check = 1

     Ieik = dlog(mu1)*(epinv - L)     ! CET eq.(A.22)

     Icoll = epinv - L     ! CET eq.(A.23)

    Ieik = Ieik -dlog(mu1)**2-two*dlog(one-mu1s)*dlog(mu1) &     ! CET eq.(A.22) 1st line
   - two*DLi2(one - mu1s)                                       ! this is simplified version of the function f(j,k)

   Icoll = Icoll + 3d0*dlog(mu1)-two*dlog(one -mu1s) &     ! CET eq.(A.23) 1st line
   - two/(one-mu1s)*dlog(mu1)+3.0_dp                       ! HERE WERE SEVERAL BUGS

     if (alpha_ff.lt.1d0) then
   Ieik = Ieik -dlog(alpha_ff)*two*dlog(mu1) &     ! CET eq.(A.22) 2nd line
   + DLi2(alpha_ff*(mu1s-one)/mu1s) - DLi2((mu1s-one)/mu1s)

   Icoll = Icoll +one/two*(3.0_dp*alpha_ff -2.0_dp &     ! CET eq.(A.23) 2nd+3rd line
          - (3.0_dp-mu1s)/(one-mu1s)*dlog(alpha_ff+(one-alpha_ff)*mu1s) &
          -alpha_ff/(alpha_ff + (one-alpha_ff)*mu1s) )  &
          - two*dlog(alpha_ff) + two*dlog(alpha_ff+(one-alpha_ff)*mu1s)/(one-mu1s)
        endif
  endif


  if (mq1.eq.zero.and.mq2.ne.zero) then   ! scheme dependence

     check = 1

       yp = (one - mu2)/(one+mu2)
       xp = yp*(one-alpha_ff)+dsqrt((one-alpha_ff)*(one-alpha_ff*yp**2))


     Ieik = epinv*(epinv2/two - dlog(one-mu2s)) &     ! CET eq.(A.29)
       -L*(epinv/two-dlog(one-mu2s)) + L**2/4.0_dp &
      +DLi2(one-mu2s)-5.0_dp/12.0_dp*pisq + dlog(one-mu2s)**2  &
      +one/two*dlog((one-yp**2+two*xp*yp)/(one+yp-xp)/(one-yp+xp))**2  &      ! HERE WAS A BUG: dlog(...)**2
      -dlog((one+yp-xp)/(one+yp))**2  &
      +two*(dlog((one+yp)/two)*dlog((one-yp+xp)/(one-yp))  &
       +dlog((one+yp)/two/yp)*dlog((one-yp**2+two*xp*yp)/(one-yp**2)) &
       +DLi2((one-yp)/(one+yp))-DLi2((one-yp**2+two*xp*yp)/(one+yp)**2) &
       +DLi2((one-yp+xp)/two)-DLi2((one-yp)/two) )

      Icoll = three/two*(epinv-L) -three*dlog(one-mu2) + 5.0_dp - mu2/(one-mu2) &     ! CET eq.(A.31)
        -two*mu2*(one-two*mu2)/(one-mu2**2)   &
       - three/two*(dlog(alpha_ff)+yp*(one-alpha_ff))


      if (scheme.eq.'fdh') Icoll = Icoll - one/two


       endif


   if (check.eq.0) then
      print *, 'required ff_qq dipole is not implemented'
!        pause
     ff_qq = 0d0
     return
   endif

   endif  ! for s=1


       ff_qq = two*Ieik  + Icoll



   end function ff_qq



   double precision  function ff_gg(mq1,mq2,p,i,j,z,s)
   real(dp), intent(in)  :: mq1,mq2,z
   real(dp), intent(in)  :: p(:,:)
   integer, intent(in) :: i,j,s
   real(dp) :: mu1, mu2, mu1s, mu2s, q2, Ieik, Icoll, L
   real(dp) :: yp,x,xp,xm,a,b,c,d,Icoll_add,Ieik_add
   integer :: check


   Ieik = zero
   Icoll = zero


   if (s.eq.1) then

   check = 0


   q2 = two*scr(p(:,i),p(:,j))+mq1**2+mq2**2
   mu1 = mq1/dsqrt(q2)
   mu2 = mq2/dsqrt(q2)
   mu1s = mu1**2
   mu2s = mu2**2

   L = dlog(q2/MuRen**2)

   if ((mq1.eq.zero).and.(mq2.ne.zero)) then


     check = 1

       yp = (one - mu2)/(one+mu2)
       xp = yp*(one-alpha_ff)+dsqrt((one-alpha_ff)*(one-alpha_ff*yp**2))

     Ieik = epinv*(epinv2/two - dlog(one-mu2s)) &      ! same as above, BCPPW eq.(A.9)
       -L*(epinv/two-dlog(one-mu2s)) + L**2/4.0_dp &
      +DLi2(one-mu2s)-5.0_dp/12.0_dp*pisq + dlog(one-mu2s)**2  &
      +one/two*dlog((one-yp**2+two*xp*yp)/(one+yp-xp)/(one-yp+xp))**2  &     ! HERE WAS A BUG: dlog(...)**2
      -dlog((one+yp-xp)/(one+yp))**2  &
      +two*(dlog((one+yp)/two)*dlog((one-yp+xp)/(one-yp))  &
       +dlog((one+yp)/two/yp)*dlog((one-yp**2+two*xp*yp)/(one-yp**2)) &
       +DLi2((one-yp)/(one+yp))-DLi2((one-yp**2+two*xp*yp)/(one+yp)**2) &
       +DLi2((one-yp+xp)/two)-DLi2((one-yp)/two) )


   Icoll = 11.0_dp/6.0_dp*(epinv-L)      ! CDST eq.(5.37)
   Icoll = Icoll + 50.0_dp/9.0_dp - 11.0_dp/3.0_dp*(mu2/(one+mu2) &
   + dlog(one - mu2) ) + (two-three*kappa_ff)*mu2s/three/(one-mu2s)*dlog(two*mu2/(one+mu2))


    if (scheme.eq.'fdh') Icoll = Icoll - 1.0_dp/6.0_dp

     Icoll_add = -11.0_dp/6.0_dp*(   &      ! BCPPW eq.(A.22)
     (one - mu2-alpha_ff*yp*(one+mu2))/(one+mu2) &
   + dlog(alpha_ff*yp*(one+mu2)/(one - mu2))  ) &
   -(kappa_ff-two/three)*mu2s/(one-mu2s)*dlog((one-alpha_ff*yp)*(one+mu2)/two/mu2)

     Icoll = Icoll + Icoll_add

   endif

   if (check.eq.0) then
      print *, 'required ff_gg dipole is not implemented'
!        pause
   endif


   endif

       ff_gg = two*(two*Ieik + Icoll)

   end function  ff_gg




!-------------------------------------------------------------------------
!    final-initial dipoles
!------------------------------------------------------------------------

    double precision function  fi_qq(mq1,mq2,p,i,j,x,s)
    real(dp), intent(in)  :: mq1,mq2,x
    real(dp), intent(in)  :: p(:,:)
    integer, intent(in) :: i,j,s
    real(dp) :: mu1, mu2, q2, Ieik, Icoll,lx,L
    real(dp) :: mu1s,mu2s
    integer :: check

    check = 0

    q2 = -two*scr(p(:,i),p(:,j))  ! the minus sign here because of all-outgoing
                                  ! convention
    mu1 = mq1/dsqrt(q2)
    mu2 = mq2/dsqrt(q2)
    mu1s = mu1**2
    mu2s = mu2**2

    L  = dlog(q2/MuRen**2)  !  HERE WAS A BUG: /MuRen**2 missing
    lx = dlog(x)

    if (mq2.ne.zero) then
       print *, 'fi_qq dipole is not implemented properly :mq2'
       stop
    endif

 if ((mq1.ne.zero).and.(mq2.eq.zero)) then !massive quark -- no scheme dependence

       check = 1

       if (s.eq.1) then
      fi_qq = (epinv-L)*(one - dlog((one+mu1s)/mu1s)) &      ! CET eq.(A.13)
        +two-two*DLi2(-mu1s)+one/two*dlog(mu1s)**2 &
        +one/two*dlog(one+mu1s)**2 - pisq/three &
        -two*dlog(mu1s)*dlog(1d0+mu1s)+dlog(mu1s) &
        + two*dlog(alpha_fi)*(dlog((1.0_dp+mu1s)/mu1s) -1.0_dp) ! alpha_fi was missed here
       endif

       if (s.eq.2) then
        if (x .gt. 1d0-alpha_fi) then
       fi_qq  =(one - x)/two/(one - x + mu1s*x)**2 &      ! CET eq.(A.14)
   +two/(one-x)*dlog(mu1s/(one+mu1s)*(two-x+mu1s*x)/(one-x+mu1s*x))
       else
         fi_qq=0d0
       endif
       endif

       if (s.eq.3) then
        if (x .gt. 1d0-alpha_fi) then
          fi_qq = two/(one-x)*(dlog((one+mu1s)/mu1s) - one)      ! CET eq.(A.14)
        else
          fi_qq = 0d0
        endif
       endif

   elseif((mq1.eq.zero).and.(mq2.eq.zero)) then


      if (s.eq.1) then

         fi_qq = epinv*(epinv2-L+lx) + L**2/two+three/two*(epinv-L) &      ! CET eq.(A.17)
              +7.0_dp/2.0_dp-pisq/two-dlog(alpha_fi)*(three/two+dlog(alpha_fi))

       if (scheme.eq.'fdh') fi_qq = fi_qq - one/two

         elseif(s.eq.2) then

            if(x.gt.1d0-alpha_fi) then
               fi_qq = two*dlog(two-x)/(one-x)      ! CET eq.(A.17)
            else
               fi_qq = zero
            endif

         elseif(s.eq.3) then

              if (x.gt.1d0-alpha_fi) then!   here was a bug: alpha --> alpha_fi
               fi_qq = -one/(one-x)*(two*dlog(one-x)+three/two)      ! CET eq.(A.17)
              else
               fi_qq = zero
              endif

      endif

      else
       print *, 'fi_qq dipole is not implemented properly'
       fi_qq = 0d0
       return
    endif



    end function  fi_qq





    double precision function  fi_gg(mq1,mq2,p,i,j,x,s)
    real(dp), intent(in)  :: mq1,mq2,x
    real(dp), intent(in)  :: p(:,:)
    integer, intent(in) :: i,j,s
    real(dp) :: mu1, mu2, q2, Ieik, Icoll, L, lx
    real(dp) :: mu1s,mu2s
    real(dp) :: theta
    integer :: check


    check = 0

    q2 = -two*scr(p(:,i),p(:,j))
    mu1 = mq1/dsqrt(q2)
    mu2 = mq2/dsqrt(q2)
    mu1s = mu1**2
    mu2s = mu2**2

    L  = dlog(q2/MuRen**2)


    if (mq1.ne.zero.or.mq2.ne.zero) then
       print *, 'fi_gg dipole is not implemented properly'
       stop
    endif

! so we will assume that all the masses are zero

    if (mq1.eq.zero.and.mq2.eq.zero) then

       check = 1

    theta = 0d0
    if (x.gt.1.0_dp-alpha_fi) theta = 1.0_dp

    if (s.eq.1) then

       fi_gg = two*(epinv*(epinv2 + 11.0_dp/6.0_dp-L) &      ! CDST eq.(5.68)
   +L**2/two - 11.0_dp/6.0_dp*L &                            ! HERE WAS A BUG: -L**2/two
   + 67.0_dp/18.0_dp - pisq/two  &
    -dlog(alpha_fi)*(11.0_dp/6.0_dp+dlog(alpha_fi)))         ! MCFM

       if (scheme.eq.'fdh') fi_gg = fi_gg - two*one/6.0_dp

    endif

    if (s.eq.2) then

       fi_gg = two*( two/(one-x)*dlog(two-x) )*theta      ! CDST eq.(5.67)

    endif


    if (s.eq.3) then

       fi_gg = two*(-two/(one-x)*dlog(one-x) - 11.0_dp/6.0_dp/(one-x))*theta      ! CDST eq.(5.67)


    endif



    endif


    if (check.eq.0) then
       print *, 'required fi_gg dipole is not implemented'
!        pause
      fi_gg = 0d0
      return
    endif



    end function  fi_gg



!------------------------------------------------------------

    double precision function  fi_qg(mq1,mq2,p,i,j,x,s)
    real(dp), intent(in)  :: mq1,mq2,x
    real(dp), intent(in)  :: p(:,:)
    integer, intent(in) :: i,j,s
    real(dp) :: mu1, mu2, q2, Ieik, Icoll, L, lx
    real(dp) :: mu1s,mu2s
    real(dp) :: theta
    integer :: check


    check = 0

    q2 = -two*scr(p(:,i),p(:,j))
    mu1 = mq1/dsqrt(q2)
    mu2 = mq2/dsqrt(q2)
    mu1s = mu1**2
    mu2s = mu2**2

    L  = dlog(q2/MuRen**2)



    if (mq1.ne.zero.or.mq2.ne.zero) then
       print *, 'fi_qg dipole is not implemented properly'
       stop
    endif

! so we will assume that all the masses are zero

    if (mq1.eq.zero.and.mq2.eq.zero) then

       check = 1

    theta = 0d0
    if (x.gt.1.0_dp-alpha_fi) theta = 1.0_dp

    if (s.eq.1) then

       fi_qg = -two/three*(epinv-L-dlog(alpha_fi)) - 10.0_dp/9.0_dp      ! CS eq.(5.58)+MCFM

    endif

    if (s.eq.2) then

       fi_qg = zero

    endif


    if (s.eq.3) then

       fi_qg = two/three/(one-x)*theta      ! CS eq.(5.58)

    endif

    endif

    if (check.eq.0) then
       print *, 'required fi_qg dipole is not implemented'
!        pause
      fi_qg = 0d0
      return
    endif


    end function  fi_qg



!------------------------------------------------------------
    double precision function  if_gg(mq1,mq2,p,i,j,x,s)
    real(dp), intent(in)  :: mq1,mq2,x
    real(dp), intent(in)  :: p(:,:)
    integer, intent(in) :: i,j,s
    real(dp) :: mu1, mu2, q2, Ieik, Icoll,L,lx,zp
    real(dp) :: mu1s,mu2s, Pggreg
    integer :: check

    if (mq1.ne.zero) then
       print *, 'if_gg dipole is not implemented correctly'
    endif

    q2 = -two*scr(p(:,i),p(:,j)) !the minus sing due to all outgoing convention
    mu1 = mq1/dsqrt(q2)
    mu2 = mq2/dsqrt(q2)
    mu1s = mu1**2
    mu2s = mu2**2

    L = dlog(q2/MuRen**2)
    lx = dlog(x)

    check = 0
    Pggreg = two*((one-x)/x-one+x*(one-x))

    if (mq2.ne.zero) then
       check = 1
       if (s.eq.1) then
         if_gg = epinv*(epinv2-L) +L**2/two + pisq/6.0_dp &      ! CT eq.(A.11) with DLi2 simplifications
         + (epinv-L)*dlog(one+mu2s)+one/two*dlog(one+mu2s)**2  &
         + two*DLi2(one/(one+mu2s))-pisq/three

         if (scheme.eq.'fdh') if_gg = if_gg - one/6.0_dp
       endif

       if (s.eq.2) then
        if_gg = -(epinv-L+lx)*Pggreg+Pggreg*dlog(one-x) &      ! CT eq.(A.11)
    - two*dlog(two-x)/(one-x) - two*dlog(x)/(one-x) &
    + Pggreg*dlog((one-x)/(one-x+mu2s*x)) &
    + two*mu2s*dlog(mu2s*x/(one-x + mu2s*x)) &
    + two/(one-x)*dlog((two-x)*(one+mu2s)/(two - x + mu2s*x))

       zp=(1d0-x)/(1d0-x+x*mu2s)
       if (alpha_if .lt. zp) then
      if_gg = if_gg - (2d0/(1d0-x)*(dlog(zp*(1d0-x+alpha_if)/(alpha_if*(1d0-x+zp)))) &      ! CT eq.(A.11)
      + Pggreg*dlog(zp/alpha_if)+2d0*mu2s*dlog((1d0-zp)/(1d0-alpha_if)))
       endif

       endif

       if (s.eq.3) then
     if_gg = four*dlog(one-x)/(one-x) &      ! CT eq.(A.11)
              -(epinv-L)*two/(one-x) &
              -two/(one-x)*dlog(one+mu2s)
       endif


    endif

    if (mq2.eq.zero) then

       check = 1

       if (s.eq.1) then
         if_gg = epinv*(epinv2-L) +L**2/two + pisq/6.0_dp      ! CT eq.(A.11) with mu=0
         if (scheme.eq.'fdh') if_gg = if_gg - one/6.0_dp
       endif

       if (s.eq.2) if_gg = -(epinv-L+lx)*Pggreg+Pggreg*(dlog(one-x)+dlog(alpha_if) ) &      ! CT eq.(A.11) with mu=0
                           - two*dlog((one-x+alpha_if)/alpha_if)/(one-x) &
                           - two*dlog(x)/(one-x)

       if (s.eq.3) if_gg = four*dlog(one-x)/(one-x) &      ! CT eq.(A.11) with mu=0
                            -(epinv-L)*two/(one-x)

    endif

    end function  if_gg



!------------------------------------------------------------
    double precision function  if_qq(mq1,mq2,p,i,j,x,s)
    real(dp), intent(in)  :: mq1,mq2,x
    real(dp), intent(in)  :: p(:,:)
    integer, intent(in) :: i,j,s
    real(dp) :: mu1, mu2, q2, Ieik, Icoll,L,lx,zp
    real(dp) :: mu1s,mu2s, Pqqreg, Pqq_pr
    integer :: check

    if (mq1.ne.zero) then
       print *, 'if_qq dipole is not implemented correctly'
    endif

    q2 = -two*scr(p(:,i),p(:,j)) !the minus sing due to all outgoing convention
    mu1 = mq1/dsqrt(q2)
    mu2 = mq2/dsqrt(q2)
    mu1s = mu1**2
    mu2s = mu2**2

    L = dlog(q2/MuRen**2)
    lx = dlog(x)


    check = 0

    Pqqreg = -one-x
    Pqq_pr = one-x

    if (mq2.ne.zero) then

       check = 1

       if (s.eq.1) then      ! CET eq.(A.9)   with DLi2 simplifications
         if_qq = epinv*(epinv2-L) +L**2/two + pisq/6.0_dp &
         + (epinv-L)*dlog(one+mu2s)+one/two*dlog(one+mu2s)**2  &
         + two*DLi2(one/(one+mu2s))-pisq/three
         if (scheme.eq.'fdh') if_qq = if_qq - one/two
       endif


       if (s.eq.2) then
        if_qq = -(epinv-L+lx)*Pqqreg+Pqqreg*dlog(one-x) &
    - two*dlog(two-x)/(one-x) - two*dlog(x)/(one-x) + Pqq_pr &
    + Pqqreg*dlog((one-x)/(one-x+mu2s*x)) &
    + two/(one-x)*dlog((two-x)*(one+mu2s)/(two - x + mu2s*x))

         zp=(1d0-x)/(1d0-x+x*mu2s)
       if (alpha_if .lt. zp) then
        if_qq=if_qq-(2d0/(1d0-x)*(dlog(zp*(1d0-x+alpha_if)/(alpha_if*(1d0-x+zp)))) &
         +Pqqreg*dlog(zp/alpha_if))
         endif

       endif

       if (s.eq.3) then
     if_qq = four*dlog(one-x)/(one-x) &
              -(epinv-L)*two/(one-x) &
              -two/(one-x)*dlog(one+mu2s)
       endif

    endif



    if (mq2.eq.zero) then

       check = 1

       if (s.eq.1) then      ! CET eq.(A.11)
         if_qq = epinv*(epinv2-L) +L**2/two + pisq/6.0_dp
         if (scheme.eq.'fdh') if_qq = if_qq - one/two
       endif

   if (s.eq.2)  if_qq = -(epinv-L+lx-dlog(alpha_if))*Pqqreg+Pqqreg*dlog(one-x) &      ! CET eq.(A.11)
    - two*dlog((one-x+alpha_if)/alpha_if)/(one-x) - two*dlog(x)/(one-x) &
     + Pqq_pr

       if (s.eq.3) if_qq = four*dlog(one-x)/(one-x) &      ! CET eq.(A.11)
                            -(epinv-L)*two/(one-x)

    endif

    end function  if_qq


!------------------------------------------------------------
!         ----q--/---g-->hard scattering

    double precision function  if_qg(mq1,mq2,p,i,j,x,s)
    real(dp), intent(in)  :: mq1,mq2,x
    real(dp), intent(in)  :: p(:,:)
    integer, intent(in) :: i,j,s
    real(dp) :: mu1, mu2, q2, Ieik, Icoll,L,lx,zp
    real(dp) :: mu1s,mu2s, Pqgreg, Pqg_prim
    integer :: check

    if (mq1.ne.zero) then
       print *, 'if_qg dipole is not implemented correctly'
    endif

    q2 = -two*scr(p(:,i),p(:,j)) !the minus sing due to all outgoing convention
    mu1 = mq1/dsqrt(q2)
    mu2 = mq2/dsqrt(q2)
    mu1s = mu1**2
    mu2s = mu2**2

    L = dlog(q2/MuRen**2)
    lx = dlog(x)


    check = 0

    Pqgreg = (one + (one-x)**2)/x
    Pqg_prim = x

       check = 1

       if (s.eq.1) if_qg = zero

       if (s.eq.2) then
          if_qg = -(epinv-L+lx)*Pqgreg+Pqgreg*dlog(one-x) &      ! CT eq.(A.10)
             +Pqg_prim

         if (mq2.ne.zero) if_qg = if_qg + Pqgreg*dlog((one-x)/(one-x+mu2s*x)) &
             + two*mu2s*dlog(mu2s*x/(one-x + mu2s*x))

         zp = (one - x)/(one - x + mu2s*x)

         if (zp.gt.alpha_if) if_qg = if_qg - Pqgreg*dlog(zp/alpha_if)      ! CT eq.(A.10)
         if (zp.gt.alpha_if.and.mq2.ne.zero)  if_qg = if_qg &
                              -two*mu2s*dlog((one-zp)/(one-alpha_if))

        endif

       if (s.eq.3) if_qg = zero


     end function  if_qg




!         ----g--/---q-->hard scattering

    double precision function  if_gq(mq1,mq2,p,i,j,x,s)
    real(dp), intent(in)  :: mq1,mq2,x
    real(dp), intent(in)  :: p(:,:)
    integer, intent(in) :: i,j,s
    real(dp) :: mu1, mu2, q2, Ieik, Icoll,L,lx
    real(dp) :: mu1s,mu2s, Pgqreg, Pgq_prim
    integer :: check
    real(dp) :: Tr,zp

    if (mq1.ne.zero) then
       print *, 'if_gq dipole is not implemented correctly'
    endif

    q2 = -two*scr(p(:,i),p(:,j)) !the minus sing due to all outgoing convention
    mu1 = mq1/dsqrt(q2)
    mu2 = mq2/dsqrt(q2)
    mu1s = mu1**2
    mu2s = mu2**2

    L = dlog(q2/MuRen**2)
    lx = dlog(x)

     Tr = one/two    ! casimir

    check = 0

    Pgqreg = Tr*(x**2 + (one-x)**2)
    Pgq_prim = two*Tr*x*(one-x)

       check = 1

       if (s.eq.1) if_gq = zero

       if (s.eq.2) then
          if_gq = -(epinv-L+lx)*Pgqreg+Pgqreg*dlog((one-x)**2/(one-x+mu2s*x)) &      ! CT eq.(A.8)
           +Pgq_prim

         zp = (one - x)/(one - x + mu2s*x)

         if (zp.gt.alpha_if) if_gq = if_gq - Pgqreg*dlog(zp/alpha_if)      ! CT eq.(A.8)

        endif

       if (s.eq.3) if_gq = zero

    end function  if_gq


!--------------------------------------------------------------------------
!    initial initial dipoles
!--------------------------------------------------------------------------

    double precision function  ii_gg(mq1,mq2,p,i,j,x,s)      ! includes alpha paramter
    real(dp), intent(in)  :: mq1,mq2,x
    real(dp), intent(in)  :: p(:,:)
    integer, intent(in) :: i,j,s
    real(dp) :: mu1, mu2, q2, Ieik, Icoll,L,lx
    real(dp) :: mu1s,mu2s, Pggreg
    integer :: check

    if (mq1.ne.zero.or.mq2.ne.zero) then
       print *, 'ii_gg dipole is not implemented correctly'
    endif

    q2 = two*scr(p(:,i),p(:,j))

    L = dlog(q2/MuRen**2)
    lx = dlog(x)


    check = 0

    Pggreg = two*((one-x)/x-one+x*(one-x))


       if (s.eq.1) then
       ii_gg = epinv*(epinv2-L) +L**2/two - pisq/6.0_dp      ! CT eq.(A.5)
         if (scheme.eq.'fdh') ii_gg = ii_gg - one/6.0_dp
       endif


       if (s.eq.2) then
            ii_gg = -(epinv-L+lx)*Pggreg + two*Pggreg*dlog(one-x) &      ! CT eq.(A.5)
                            - two*dlog(x)/(one-x)
            if (alpha_ii/(1d0-x) .lt. 1d0) ii_gg=ii_gg +(2d0/(1d0-x)+Pggreg)*dlog(alpha_ii/(1d0-x))   ! NEW
       endif

       if (s.eq.3) ii_gg = four*dlog(one-x)/(one-x) &      ! CT eq.(A.5)
                            -(epinv-L)*two/(one-x)

    end function  ii_gg




!         ----q--/---q-->hard scattering
    double precision function  ii_qq(mq1,mq2,p,i,j,x,s)
    real(dp), intent(in)  :: mq1,mq2,x
    real(dp), intent(in)  :: p(:,:)
    integer, intent(in) :: i,j,s
    real(dp) :: mu1, mu2, q2, Ieik, Icoll,L,lx
    real(dp) :: mu1s,mu2s, Pqqreg
    integer :: check

    if (mq1.ne.zero.or.mq2.ne.zero) then
       print *, 'ii_qq dipole is not implemented correctly'
    endif

    q2 = two*scr(p(:,i),p(:,j))
    L = dlog(q2/MuRen**2)
    lx = dlog(x)


    check = 0

    Pqqreg = -one-x


       if (s.eq.1) then
       ii_qq = epinv*(epinv2-L) +L**2/two - pisq/6.0_dp      ! CET eq.(A.4)
       if (scheme.eq.'fdh') ii_qq = ii_qq - one/two
       endif

       if (s.eq.2) then
          ii_qq = -(epinv-L+lx)*Pqqreg+two*Pqqreg*dlog(one-x) - two*dlog(x)/(one-x)+one-x      ! CET eq.(A.4)
          if (alpha_ii/(1d0-x) .lt. 1d0) ii_qq=ii_qq + (2d0/(1d0-x)+Pqqreg)*dlog(alpha_ii/(1d0-x))   ! NEW
       endif

       if (s.eq.3) ii_qq = four*dlog(one-x)/(one-x) &      ! CET eq.(A.4)
                            -(epinv-L)*two/(one-x)


    end function  ii_qq




!         ----q--/---g-->hard scattering

    double precision function  ii_qg(mq1,mq2,p,i,j,x,s)
    real(dp), intent(in)  :: mq1,mq2,x
    real(dp), intent(in)  :: p(:,:)
    integer, intent(in) :: i,j,s
    real(dp) :: mu1, mu2, q2, Ieik, Icoll,L,lx
    real(dp) :: mu1s,mu2s, Pqgreg, Pqg_prim
    integer :: check

    if (mq1.ne.zero.or.mq2.ne.zero) then
       print *, 'ii_qg dipole is not implemented correctly'
    endif

    q2 = two*scr(p(:,i),p(:,j))

    L = dlog(q2/MuRen**2)
    lx = dlog(x)


    check = 0

    Pqgreg = (one + (one-x)**2)/x
    Pqg_prim = x


       if (s.eq.1) ii_qg = zero

       if (s.eq.2) then
          ii_qg = -(epinv-L+lx)*Pqgreg+two*Pqgreg*dlog(one-x)+Pqg_prim ! HERE WAS A BUG: Pgq_prim missing
          if (alpha_ii/(1d0-x) .lt. 1d0) ii_qg = ii_qg + Pqgreg*dlog(alpha_ii/(1d0-x))      ! MCFM
       endif

       if (s.eq.3) ii_qg = zero


    end function  ii_qg




!         ----g--/---q-->hard scattering
    double precision function  ii_gq(mq1,mq2,p,i,j,x,s)   ! needs a correction factor of 2 to fix normalization
    real(dp), intent(in)  :: mq1,mq2,x
    real(dp), intent(in)  :: p(:,:)
    integer, intent(in) :: i,j,s
    real(dp) :: mu1, mu2, q2, Ieik, Icoll,L,lx
    real(dp) :: mu1s,mu2s, Pgqreg, Pgq_prim
    integer :: check
    real(dp) :: Tr

    if (mq1.ne.zero.or.mq2.ne.zero) then
       print *, 'ii_gq dipole is not implemented correctly'
    endif

    q2 = two*scr(p(:,i),p(:,j))
    L = dlog(q2/MuRen**2)
    lx = dlog(x)

!     print *, q2,l,lx
!     pause

    Tr = one/two  !casimir

    check = 0

    Pgqreg = Tr*(x**2 + (one-x)**2)
    Pgq_prim = two*Tr*x*(one-x)


       if (s.eq.1) ii_gq = zero

       if (s.eq.2) then
          ii_gq = -(epinv-L+lx)*Pgqreg+two*Pgqreg*dlog(one-x)+ Pgq_prim ! HERE WAS A BUG: Pgq_prim missing      ! CET eq.(A.5)
          if (alpha_ii/(1d0-x) .lt. 1d0) ii_gq=ii_gq + Pgqreg*dlog(alpha_ii/(1d0-x))
       endif

       if (s.eq.3) ii_gq = zero


    end function  ii_gq














! ************************************************************************
! *     Author: J. M. Campbell                                           *
! *     April 15th, 2004                                                 *
! *                                                                      *
! *     Routines which return various pieces of the integrated           *
! *     MASSIVE subtraction terms, used in both _v and _z routines       *
! *                                                                      *
! *     The formulae implemented here are derived in the FORM programs   *
! *     testif.frm, testif_gg.frm and testfi.frm                         *
! *                                                                      *
! *     Other final-initial and all final-final dipoles are untested     *
! ************************************************************************
!
! ************************************************************************
! *                                                                      *
! *     The labelling of the routines is as follows:                     *
! *     The collinear pair is assumed to be incoming,                    *
! *     so a reversal has to be made for the final state cases           *
! *                                                                      *
! *              -------->------------>--------                          *
! *                  j        /         i                                *
! *                          /                                           *
! *                         /                                            *
! *                                                                      *
! *                represented by {ii/if}_ij                             *
! *                                                                      *
! ************************************************************************


! ***********************************************************************
! *************************** Initial-INITIAL ***************************
! ***********************************************************************
! ***************************** Quark-Quark *****************************
      double precision function ii_mqq(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,mbar,Pqqreg,alfax
!--- returns the integral of the subtraction term for an
!--- initial-initial quark-gluon antenna, either
!--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)

      if (vorz .eq. 1) then
        ii_mqq=epinv*(epinv2-L)+0.5d0*L**2-pisqo6
        if (scheme .eq. 'tH-V') then
           return
        elseif (scheme .eq. 'dred') then
           ii_mqq=ii_mqq-half
           return
        endif
      endif

      omx=one-x
      lomx=dlog(omx)
      lx=dlog(x)

      if (vorz .eq. 2) then
        Pqqreg=-one-x
        ii_mqq=omx+Pqqreg*(two*lomx+L-epinv)-(one+x**2)/omx*lx
        alfax=alpha_ii/omx
        if (alfax .lt. 1d0) ii_mqq=ii_mqq+(two/omx+Pqqreg)*dlog(alfax)
        return
      endif

      ii_mqq=two/omx*(two*lomx+L-epinv)

      return
      end function ii_mqq



!***************************** Quark-Gluon *****************************
      double precision function ii_mqg(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,Pqgreg,alfax,mbar
!--- returns the integral of the subtraction term for an
!--- initial-initial gluon-quark antenna, either
      ii_mqg=0d0
      if ((vorz .eq. 1) .or. (vorz .eq. 3)) return

      omx=one-x
      lomx=dlog(omx)
      lx=dlog(x)

      if (vorz .eq. 2) then
        Pqgreg=one-two*x*omx
        ii_mqg=Pqgreg*(two*lomx-lx+L-epinv)+two*x*omx
        alfax=alpha_ii/omx
        if (alfax .lt. 1d0) ii_mqg=ii_mqg+Pqgreg*dlog(alfax)
      endif
      return
      end function ii_mqg



!***************************** Gluon-Quark *****************************
      double precision function ii_mgq(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,Pgqreg,alfax,mbar
!--- returns the integral of the subtraction term for an
!--- initial-initial quark-quark (--> gluon) antenna, either
!--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)

      ii_mgq=0d0
      if ((vorz .eq. 1) .or. (vorz .eq. 3)) return

      omx=one-x
      lomx=dlog(omx)
      lx=dlog(x)

      if (vorz .eq. 2) then
        Pgqreg=(one+omx**2)/x
        ii_mgq=Pgqreg*(two*lomx-lx+L-epinv)+x
        alfax=alpha_ii/omx
        if (alfax .lt. 1d0) ii_mgq=ii_mgq+Pgqreg*dlog(alfax)
        return
      endif

      return
      end function ii_mgq



!***************************** Gluon-Gluon *****************************
      double precision function ii_mgg(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,lomx,Pggreg,alfax,mbar
!--- returns the integral of the subtraction term for an
!--- initial-initial gluon-gluon antenna, either
!--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)

      if (vorz .eq. 1) then
        ii_mgg=epinv*(epinv2-L)+half*L**2-pisqo6
        if (scheme .eq. 'tH-V') then
          return
        elseif (scheme .eq. 'dred') then
          ii_mgg=ii_mgg-1d0/6d0
          return
        endif
      endif

      omx=one-x
      lomx=dlog(omx)

      if (vorz .eq. 2) then
        Pggreg=omx/x+x*omx-one
        lx=dlog(x)
        ii_mgg=two*Pggreg*(two*lomx-lx+L-epinv)-two*lx/omx
        alfax=alpha_ii/omx
        if (alfax .lt. 1d0) ii_mgg=ii_mgg+two*(one/omx+Pggreg)*dlog(alfax)
        return
      endif

      ii_mgg=two*(two*lomx+L-epinv)/omx

      return
      end function ii_mgg



!***********************************************************************
!**************************** INITIAL-FINAL ****************************
!***********************************************************************
!***************************** Quark-Quark *****************************
      double precision function if_mqq(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,mbarsq,omx,lx,lomx,Pqqreg,zp
!--- returns the integral of the subtraction term for an
!--- initial-final quark-gluon antenna, either
!--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)

      mbarsq=mbar**2
      if_mqq=0d0
      if (vorz .eq. 1) then
          if_mqq=(epinv+dlog(1d0+mbarsq))*(epinv-L)+half*L**2-half*dlog(1d0+mbarsq)**2+2d0*dlog(mbarsq)*dlog(1d0+mbarsq)+2d0*DLi2(-mbarsq)+pisqo6
          if (scheme .eq. 'tH-V') then
            return
          elseif (scheme .eq. 'dred') then
            if_mqq=if_mqq-half
            return
          endif
      endif
      omx=one-x
      lomx=dlog(omx)
      zp=omx/(omx+x*mbarsq)
      if (vorz .eq. 2) then
         Pqqreg=-(1d0+x)
         lx=dlog(x)
         if_mqq=Pqqreg*(-(epinv-L)+2d0*lomx-lx-dlog(x*mbarsq+omx)) +omx-2d0/omx*(lx+dlog((1d0+x*mbarsq+omx)/(1d0+mbarsq)))
         if (alpha_if .lt. zp) then
            if_mqq=if_mqq-(two/omx*(dlog(zp*(omx+alpha_if)/(alpha_if*(omx+zp)))) +Pqqreg*dlog(zp/alpha_if))
         endif
      elseif (vorz .eq. 3) then
         if_mqq=2d0/omx*(-(epinv-L)+2d0*lomx-dlog(1d0+mbarsq))
      endif
      return
      end function



!***************************** Gluon-Gluon *****************************
      double precision function if_mgg(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,omx,lx,mbar,mbarsq,Pggreg,zp

      mbarsq=mbar**2
      if_mgg=0d0

      if (vorz .eq. 1) then  ! this smells like a BUG: singular part should be:  epinv*( (epinv2-L) + dlog(1d0+mbarsq) )
         if_mgg=(epinv+dlog(1d0+mbarsq))*(epinv-L)+half*L**2 -half*dlog(1d0+mbarsq)**2+2d0*dlog(mbarsq)*dlog(1d0+mbarsq)+2d0*DLi2(-mbarsq)+pisqo6

         if (scheme .eq. 'tH-V') then
           return
         elseif (scheme .eq. 'dred') then
           if_mgg=if_mgg-1d0/6d0
           return
         endif
         return
      endif
      omx=one-x
      zp=omx/(omx+x*mbarsq)
      if (vorz .eq. 2) then
         Pggreg=2d0*(omx/x-1d0+x*omx)
         lx=dlog(x)
         if_mgg=Pggreg*(-(epinv-L)+2d0*dlog(omx)-lx-dlog(x*mbarsq+omx)) +2d0*mbarsq*dlog(x*mbarsq/(x*mbarsq+omx)) -2d0/omx*(lx+dlog((1d0+x*mbarsq+omx)/(1d0+mbarsq)))
         if (alpha_if .lt. zp) then
           if (alpha_if .eq. 1d0) then
             write(6,*) 'zp > 1 in dipoles_mass.f - this is forbidden'
             stop
           endif
           if_mgg=if_mgg-(two/omx*(dlog(zp*(omx+alpha_if)/(alpha_if*(omx+zp)))) +Pggreg*dlog(zp/alpha_if)+2d0*mbarsq*dlog((1d0-zp)/(1d0-alpha_if)))
         endif
         return
      elseif (vorz .eq. 3) then
         if_mgg=2d0/omx*(-(epinv-L)+2d0*dlog(omx)-dlog(1d0+mbarsq))
         return
      endif
      return
      end function



!***************************** Quark-Gluon *****************************
!--- Not necessary because for off-diagonal (no soft singularity)
!--- we always choose to use the initial spectator

!***************************** Gluon-Quark *****************************
!--- Not necessary because for off-diagonal (no soft singularity)
!--- we always choose to use the initial spectator


!***********************************************************************
!**************************** FINAL-INITIAL ****************************
!***********************************************************************

!***************************** Quark-Quark *****************************
      double precision function fi_mqq(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,mbarsq,omx
!--- returns the integral of the subtraction term for an
!--- final-initial quark-gluon antenna, either
!--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)

      mbarsq=mbar**2
      if (vorz .eq. 1) then
        fi_mqq=(1d0+dlog(mbarsq/(1d0+mbarsq)))*(epinv-L) +dlog(mbarsq)+half*dlog(mbarsq)**2 +half*dlog(1d0+mbarsq)**2 -2d0*dlog(mbarsq)*dlog(1d0+mbarsq) -2d0*DLi2(-mbarsq)+2d0-pisq/3d0 +2d0*dlog(alpha_fi)*(dlog((1d0+mbarsq)/mbarsq)-1d0)


        return
      endif

      omx=one-x

      if (vorz .eq. 2) then
        if (x .gt. 1d0-alpha_fi) then
          fi_mqq=+omx/2d0/(x*mbarsq+omx)**2 +2d0/omx*(dlog((1d0+x*mbarsq+omx)*mbarsq/((1d0+mbarsq)*(omx+x*mbarsq))))
        else
          fi_mqq=0d0
        endif
        return
      endif

      if (vorz .eq. 3) then
        if (x .gt. 1d0-alpha_fi) then
          fi_mqq=2d0/omx*(dlog((1d0+mbarsq)/(mbarsq))-1d0)
        else
          fi_mqq=0d0
        endif
      endif

      return
      end function




!************************************************************************
!************************************************************************
!*              BELOW HERE FUNCTIONS HAVE NOT BEEN CHECKED              *
!************************************************************************
!************************************************************************


!***************************** Gluon-Gluon *****************************
      double precision function fi_mgg(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,omx
!--- returns the integral of the subtraction term for an
!--- final-initial gluon-gluon antenna, either
!--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)
!--MSbar
!-- Id,aqg=(-2/3*(epinv-L)-10/9)*[delta(1-x)]
!--  +0
!--  +2/3/[1-xp]
!-- Id,agg=
!--  (2*epinv*(epinv-L)+L^2+(epinv-L)*11/3+67/9-[pi]^2)
!--  *[delta(1-x)]
!--  +4*[ln(2-x)]/[1-x]
!--  +2*(-2*[ln(1-x)/(1-xp)]-11/6/[1-xp])


      if (vorz .eq. 1) then
        fi_mgg=two*epinv*(epinv2-L)+L**2+67d0/9d0-pisq+11d0*(epinv-L)/3d0
      endif

      omx=one-x

      if (vorz .eq. 2) then
        fi_mgg=four*dlog(two-x)/omx
        return
      endif

      fi_mgg=-(four*dlog(omx)+11d0/3d0)/omx
      return
      end function





!***************************** Quark-Gluon *****************************
      double precision function fi_mqg(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,mbarsq,omx,rt,JaS,JaNS
!--- returns the integral of the subtraction term for an
!--- final-initial gluon-quark antenna, either
!--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)
!--MSbar

      mbarsq=mbar**2
      if (1d0-4d0*mbarsq .lt. 0d0) then
      write(6,*) 'error in fi_mqg,(1d0-4d0*mbarsq .lt. 0d0)'
      stop
      else
      rt=dsqrt(1d0-4d0*mbarsq)
      endif
      if (vorz .eq. 1) then
!CDTS 5.63
      JaS=-2d0/3d0*dlog(mbarsq)-10d0/9d0
!CDTS 5.64
      JaNS=10d0/9d0*(1d0-rt)-8d0/9d0*mbarsq*rt +4d0/3d0*dlog(half*(1d0+rt))
      fi_mqg=JaS+JaNS

      elseif (vorz .eq. 2) then
!C regular
        fi_mqg=0d0
      elseif (vorz .eq. 3) then
!C plus at x=x+
      omx=1d0-x
!CDTS 5.62
      fi_mqg=2d0/3d0*(omx+2d0*mbarsq)/omx**2*dsqrt(1d0-4d0*mbarsq/omx)
      endif
      return
      end function




!***********************************************************************
!***************************** FINAL-FINAL *****************************
!***********************************************************************
!***************************** Quark-Quark *****************************
      double precision function ff_1mqq(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,mbarsq,afftmp,Icolla,Ieika,Icollb,Ieikb,arg,ommsq,logm,logomm,xp,yp,arg1,arg2,arg3,ypp,ypm
!       include 'phi.f'
      double precision phi

      ff_1mqq=0d0
      phi=1d0


      if (vorz .eq. 1) then
        if (mbar .gt. 1d0) then
            write(6,*) 'Problem with mbar in ff_1mqq, mbar=',mbar
            stop
        endif


      mbarsq=mbar**2
      ommsq=1d0-mbarsq
      logm=dlog(mbarsq)
      logomm=dlog(ommsq)

!----radiation from massive line with massless spectator
      afftmp=alpha_ff
      arg=afftmp+(1d0-afftmp)*mbarsq
      Ieika=half*logm*(epinv-L)-2d0*DLi2(ommsq)-logm*logomm-0.25d0*logm**2-dlog(afftmp)*logm-DLi2(-ommsq/mbarsq)+DLi2(-afftmp*ommsq/mbarsq)
      Icolla=epinv-L+phi+2d0+(1d0+half*phi)*logm+(phi-2d0)*logm/(ommsq)-2d0*logomm+half*phi*(3d0*afftmp-2d0-(3d0-mbarsq)/ommsq*dlog(arg)-afftmp/arg)-2d0*dlog(afftmp)+2d0*dlog(arg)/ommsq


!----radiation from massless line with massive spectator
      yp=(1d0-mbar)/(1d0+mbar)
      afftmp=alpha_ff*yp
      xp=(yp-afftmp)+dsqrt((yp-afftmp)*(1d0/yp-afftmp))
      arg1=0.25d0*(1d0-yp**2+2*xp*yp)
      arg2=half*(1d0+YP-xp)
      arg3=half*(1d0-YP+xp)
      ypp=half*(1d0+yp)
      ypm=half*(1d0-yp)

      Ieikb=half*epinv*epinv2-half*epinv*L+0.25d0*L**2-logomm*(epinv-L)+DLi2(ommsq)-2.5d0*pisqo6+logomm**2+half*dlog(arg1/(arg2*arg3))**2-dLOG(arg2/ypp)**2+2d0*(dLOG(ypp)*dLOG(arg3/ypm)+dLOG(ypp/yp)*dLOG(arg1/(ypp*ypm))+DLi2(ypm/ypp)-DLi2(arg1/ypp**2)+DLi2(arg3)-DLi2(ypm))

      Icollb=1.5d0*(epinv-L)-3d0*dlog(1d0-mbar)+5d0-mbar/(1d0-mbar)-2d0*mbar*(1d0-2d0*mbar)/ommsq+1.5d0*(dLOG(YP/afftmp)-YP+afftmp)

!--- Note: extra factor of half because we include this term once for each
!---  leg, but this is the sum of both legs
      ff_1mqq=half*(2d0*(Ieika+Ieikb)+Icolla+Icollb)


        if (scheme .eq. 'tH-V') then
           return
        elseif (scheme .eq. 'dred') then
!--- Note: extra factor of half because we include this term once for each
!---  leg, but this is the sum of both legs (as above)
           ff_1mqq=ff_1mqq-half*half
           return
        endif
      endif

      return
      end function





!***************************** Quark-Quark *****************************
      double precision function ff_mqq(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,mbarsq,Lro,ro,vtijk,arg
!C     mbarsq=mass**2/Qsq
!C     L=dlog(Qsq/MuRen**2)
!--- returns the integral of the subtraction term for an
!--- final-initial quark-gluon antenna, either
!--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)
!C --MSbar
!c Id,aqq=epinv*(epinv-L)+1/2*L^2+3/2*(epinv-L)+5-[pi]^2/2

!c g zipffqq=colfac*del(omx)*(
!c  1/vtijk*((epinv-L)*Lro-1/2*Lro**2-2*Lro*dlog(1-4*mbarsq)
!c           +4*DLi2(-ro)-6*DLi2(1-ro)+pisq/3)
!c +epinv-L+3-2*dlog(1-2*mbar)+dlog(1-mbar)+1/2*dlog(mbarsq)
!c  -2*mbar/(1-2*mbarsq)*(1-2*mbar)-mbar/(1-mbar)
!c  -2*mbarsq/(1-2*mbarsq)*dlog(mbar/(1-mu))
!c  )-ffqq;

      ff_mqq=0d0
      if (vorz .eq. 1) then
        mbarsq=mbar**2
        arg=1d0-4d0*mbarsq
        if (arg .lt. 0d0) then
            write(6,*) 'Threshold problem in ff_mqq'
            stop
        endif
        vtijk=dsqrt(arg)/(1d0-2d0*mbarsq)
        ro=dsqrt((1d0-vtijk)/(1d0+vtijk))
        Lro=dlog(ro)
        ff_mqq=1d0/vtijk*((epinv-L)*Lro-half*Lro**2-2d0*Lro*dlog(arg)+4d0*DLi2(-ro)-6d0*DLi2(1d0-ro)+pisq/3d0) +epinv-L+3d0 -2d0*dlog(1d0-2d0*mbar)+dlog(1d0-mbar)+half*dlog(mbarsq)-2d0*mbar/(1d0-2d0*mbarsq)*(1d0-2d0*mbar)-mbar/(1d0-mbar)-2d0*mbarsq/(1d0-2d0*mbarsq)*dlog(mbar/(1d0-mbar))
!c        Ieik=1d0/vtijk*(half*(epinv-L)*Lro-Lro*dlog(arg)
!c     .  -dlog(roj)**2+pisq/6d0
!c     .  +2d0*DLi2(-ro)-2d0*DLi2(1d0-ro)-DLi2(1d0-roj**2))
!c        Icoll=(epinv-L)+half*Lmbarsq-2d0-2d0*dlog((1d0-mbar)**2-mbarsq)
!c     . +dlog(1d0-mbar)-2d0*mbarsq/(1d0-2d0*mbarsq)*dlog(mbar/(1d0-mbar))
!c     . +5d0-mbar/(1d0-mbar)-2d0*mbar*(1d0-2d0*mbar)/(1d0-2d0*mbarsq)
!c        ff_mqq=2d0*Ieik+Icoll
        return
      endif
      return
      end function





!***************************** Quark-Gluon *****************************
      double precision function ff_mqg(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,mbarsq,ro,arg
!-----L=dlog(Qsq/MuRen**2)

!--- returns the integral of the subtraction term for an
!--- final-initial gluon-quark antenna, either
!--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)
!C --MSbar
!c Id,aqg=-2/3*(epinv-L)-16/9
!c
      ff_mqg=0d0
      if (vorz .eq. 1) then
        mbarsq=mbar**2
        arg=1d0-4d0*mbarsq
        if (arg .lt. 0d0) then
        write(6,*) 'Threshold problem in ff_mqg'
        stop
        else
        ro=dsqrt(arg)
        ff_mqg=-2d0/3d0*(2d0*dlog(mbarsq)-2d0*dlog(half*(1d0+ro))+2d0/3d0*ro*(3d0+ro**2))
        return
        endif
      endif
      return
      end function




!***************************** Gluon-Gluon *****************************
      double precision function ff_mgg(x,L,mbar,vorz)
      implicit none
      integer vorz
      double precision x,L,mbar,Ieik,Icoll,xn
!--- returns the integral of the subtraction term for an
!--- final-initial gluon-gluon antenna, either
!--- divergent for _v (vorz=1) or finite for _z (vorz=2,3 for reg,plus)
!C --MSbar
!c Id,aqg=-2/3*(epinv-L)-16/9
!c Id,agg=2*epinv*(epinv-L)+L^2+11/3*(epinv-L)+100/9-[pi]^2

      xn=3d0

      ff_mgg=0d0
!CDST Eq.5.32
      if (vorz .eq. 1) then
      Ieik=half*epinv*(epinv2-L)+half*L**2-pisq/4d0
      Icoll=11d0/6d0*(epinv-L)+50d0/9d0
!CDST Eq.5.36 needs to be added too - need to fix this
      Icoll=Icoll+half*(5d0)/xn*(-2d0/3d0*(epinv-L)-16d0/9d0)  ! 5d0=nf
!       Icoll=Icoll+half*(Nf_light)/xn*(-2d0/3d0*(epinv-L)-16d0/9d0)
      ff_mgg=2d0*(2d0*Ieik+Icoll)
      return
      endif
      return
      end function






! * Now write these expressions in a neater form
! g zipffqq=colfac*del(omx)*(
!  1/vtijk*((epinv-L)*Lro-1/2*Lro**2-2*Lro*dlog(1-4*mbarsq)
!           +4*DLi2(-ro)-6*DLi2(1-ro)+pisq/3)
! +epinv-L+3-2*dlog(1-2*mbar)+dlog(1-mbar)+1/2*dlog(mbarsq)
!  -2*mbar/(1-2*mbarsq)*(1-2*mbar)-mbar/(1-mbar)
!  -2*mbarsq/(1-2*mbarsq)*dlog(mbar/(1-mu))
!  )-ffqq;

! g zipifqq=colfac*(
!  del(omx)*((epinv+dlog(1+mbarsq))*(epinv-L)+1/2*L**2
!            +1/2*dlog(1+mbarsq)**2+2*DLi2(1/(1+MuRen**2))-pisq/6)
! +Preg(q,q)/colfac*(-(epinv-L)+2*dlog(omx)-dlog(x)-dlog(x*mbarsq+omx))
!  +omx-2/omx*(dlog(x)+dlog(1+x*mbarsq+omx,1+mbarsq))
! +2/omxp*(-[epinv-L]+2*dlog(omx)-dlog(1+mbarsq))
!  )-ifqq;

! g zipifgg=colfac*(
!  del(omx)*((epinv+dlog(1+mbarsq))*(epinv-L)+1/2*L**2
!            +1/2*dlog(1+mbarsq)**2+2*DLi2(1/(1+MuRen**2))-pisq/6)
! +Preg(g,g)/colfac*(-(epinv-L)+2*dlog(omx)-dlog(x)-dlog(x*mbarsq+omx))
!  +2*mbarsq*dlog(x*mbarsq,x*mbarsq+omx)
!  -2/omx*(dlog(x)+dlog(1+x*mbarsq+omx,1+mbarsq))
! +2/omxp*(-[epinv-L]+2*dlog(omx)-dlog(1+mbarsq))
!  )-ifgg;

! g zipfiqq=colfac*(
!  del(omx)*((1+dlog(mbarsq,1+mbarsq))*(epinv-L)
!            +1/2*dlog(mbarsq)-1/2*dlog(mbarsq)**2
!            +1/2*dlog(1+mbarsq)+1/2*dlog(1+mbarsq)**2
!            -2*dlog(mbarsq)*dlog(1+mbarsq)
!            -4*DLi2(-mbarsq)+mbarsq/2/(1+mbarsq)
!            +3/2-2/3*pisq)
! +2/omx*(dlog(1+x*mbarsq+omx,1+mbarsq))
! +omx/2/(x*mbarsq+omx)**2
! +2/omxp*(dlog(1+mbarsq,omx+x*mbarsq)-1)
!  )-fiqq;












SUBROUTINE InitProcess_TbTGGG(ExtParticles)
use ModMisc
use ModProcess
implicit none
type(Particle) :: ExtParticles(:)

  ExtParticles(1)%PartType = -5
  ExtParticles(1)%ExtRef   = 1
  ExtParticles(1)%Mass = m_Top
  ExtParticles(1)%Mass2= ExtParticles(1)%Mass**2

  ExtParticles(2)%PartType = 5
  ExtParticles(2)%ExtRef   = 2
  ExtParticles(2)%Mass = m_Top
  ExtParticles(2)%Mass2= ExtParticles(2)%Mass**2

  ExtParticles(3)%PartType = 10
  ExtParticles(3)%ExtRef   = 3
  ExtParticles(3)%Mass = 0d0
  ExtParticles(3)%Mass2= 0d0

  ExtParticles(4)%PartType = 10
  ExtParticles(4)%ExtRef   = 4
  ExtParticles(4)%Mass = 0d0
  ExtParticles(4)%Mass2= 0d0

  ExtParticles(5)%PartType = 10
  ExtParticles(5)%ExtRef   = 5
  ExtParticles(5)%Mass = 0d0
  ExtParticles(5)%Mass2= 0d0


RETURN
END SUBROUTINE InitProcess_TbTGGG



SUBROUTINE InitProcess_TbTGG(ExtParticles)
use ModMisc
use ModProcess
implicit none
type(Particle) :: ExtParticles(:)

  ExtParticles(1)%PartType = -5
  ExtParticles(1)%ExtRef   = 1
  ExtParticles(1)%Mass = m_Top
  ExtParticles(1)%Mass2= ExtParticles(1)%Mass**2

  ExtParticles(2)%PartType = 5
  ExtParticles(2)%ExtRef   = 2
  ExtParticles(2)%Mass = m_Top
  ExtParticles(2)%Mass2= ExtParticles(2)%Mass**2

  ExtParticles(3)%PartType = 10
  ExtParticles(3)%ExtRef   = 3
  ExtParticles(3)%Mass = 0d0
  ExtParticles(3)%Mass2= 0d0

  ExtParticles(4)%PartType = 10
  ExtParticles(4)%ExtRef   = 4
  ExtParticles(4)%Mass = 0d0
  ExtParticles(4)%Mass2= 0d0

RETURN
END SUBROUTINE InitProcess_TbTGG


    SUBROUTINE InitProcess_TbTQbQG(ExtParticles) ! NEW
    use ModMisc
    use ModProcess
    implicit none
    type(Particle) :: ExtParticles(:)

      ExtParticles(1)%PartType = -5
      ExtParticles(1)%ExtRef   = 1
      ExtParticles(1)%Mass = m_Top
      ExtParticles(1)%Mass2= ExtParticles(1)%Mass**2

      ExtParticles(2)%PartType = 5
      ExtParticles(2)%ExtRef   = 2
      ExtParticles(2)%Mass = m_Top
      ExtParticles(2)%Mass2= ExtParticles(2)%Mass**2

      ExtParticles(3)%PartType = -4
      ExtParticles(3)%ExtRef   = 3
      ExtParticles(3)%Mass = 0d0
      ExtParticles(3)%Mass2= 0d0

      ExtParticles(4)%PartType = 4
      ExtParticles(4)%ExtRef   = 4
      ExtParticles(4)%Mass = 0d0
      ExtParticles(4)%Mass2= 0d0

      ExtParticles(5)%PartType = 10
      ExtParticles(5)%ExtRef   = 5
      ExtParticles(5)%Mass = 0d0
      ExtParticles(5)%Mass2= 0d0

    RETURN
    END SUBROUTINE InitProcess_TbTQbQG




    SUBROUTINE InitProcess_TbTQbQ(ExtParticles)
    use ModMisc
    use ModProcess
    implicit none
    type(Particle) :: ExtParticles(:)

      ExtParticles(1)%PartType = -5
      ExtParticles(1)%ExtRef   = 1
      ExtParticles(1)%Mass = m_Top
      ExtParticles(1)%Mass2= ExtParticles(1)%Mass**2

      ExtParticles(2)%PartType = 5
      ExtParticles(2)%ExtRef   = 2
      ExtParticles(2)%Mass = m_Top
      ExtParticles(2)%Mass2= ExtParticles(2)%Mass**2

      ExtParticles(3)%PartType = -4
      ExtParticles(3)%ExtRef   = 3
      ExtParticles(3)%Mass = 0d0
      ExtParticles(3)%Mass2= 0d0

      ExtParticles(4)%PartType = 4
      ExtParticles(4)%ExtRef   = 4
      ExtParticles(4)%Mass = 0d0
      ExtParticles(4)%Mass2= 0d0

    RETURN
    END SUBROUTINE InitProcess_TbTQbQ




SUBROUTINE GenerateEventTTQQG2(Mom,MomDK,Hel,ExtParticles)
use ModMisc
use ModProcess
use ModTopDecay
implicit none
type(Particle) :: ExtParticles(:)
real(8) :: Mom(1:4,1:5),MomDK(1:4,1:6)
integer :: Hel(1:5)

     ExtParticles(1)%Mom(1:4) = dcmplx(Mom(1:4,1))
     ExtParticles(2)%Mom(1:4) = dcmplx(Mom(1:4,2))
     if (TopDecays.ge.1) then
        call TopDecay(ExtParticles(1),DK_LO,MomDK(1:4,1:3))
        call TopDecay(ExtParticles(2),DK_LO,MomDK(1:4,4:6))
     else
        call vSpi(ExtParticles(1)%Mom(1:4),ExtParticles(1)%Mass,Hel(1),ExtParticles(1)%Pol(1:4))
        call ubarSpi(ExtParticles(2)%Mom(1:4),ExtParticles(2)%Mass,Hel(2),ExtParticles(2)%Pol(1:4))
    endif

    ExtParticles(3)%Mom(1:4) = dcmplx(Mom(1:4,3))
    call vSpi(ExtParticles(3)%Mom(1:4),ExtParticles(3)%Mass,Hel(3),ExtParticles(3)%Pol(1:4))      ! THERE WAS A BUG HERE: POL_MLESS

    ExtParticles(4)%Mom(1:4) = dcmplx(Mom(1:4,4))
    call ubarSpi(ExtParticles(4)%Mom(1:4),ExtParticles(4)%Mass,Hel(4),ExtParticles(4)%Pol(1:4))   ! THERE WAS A BUG HERE: POL_MLESS

    ExtParticles(5)%Mom(1:4) = dcmplx(Mom(1:4,5))
    call pol_mless(ExtParticles(5)%Mom(1:4),Hel(5),ExtParticles(5)%Pol(1:4))

RETURN
END SUBROUTINE


SUBROUTINE GenerateEvent52(Mom,MomDK,Hel,ExtParticles)
use ModMisc
use ModProcess
use ModTopDecay
implicit none
type(Particle) :: ExtParticles(:)
real(8) :: Mom(1:4,1:5),MomDK(1:4,1:6)
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

    ExtParticles(3)%Mom(1:4) = dcmplx(Mom(1:4,3))
    call pol_mless(ExtParticles(3)%Mom(1:4),Hel(3),ExtParticles(3)%Pol(1:4))

    ExtParticles(4)%Mom(1:4) = dcmplx(Mom(1:4,4))
    call pol_mless(ExtParticles(4)%Mom(1:4),Hel(4),ExtParticles(4)%Pol(1:4))

    ExtParticles(5)%Mom(1:4) = dcmplx(Mom(1:4,5))
    call pol_mless(ExtParticles(5)%Mom(1:4),Hel(5),ExtParticles(5)%Pol(1:4))

RETURN
END SUBROUTINE



SUBROUTINE GenerateEventTTQQG(Mom,Hel,ExtParticles)
use ModMisc
use ModProcess
implicit none
type(Particle) :: ExtParticles(:)
real(8) :: Mom(1:4,1:5)
integer :: Hel(1:5)


    ExtParticles(1)%Mom(1:4) = dcmplx(Mom(1:4,1))
    call vSpi(ExtParticles(1)%Mom(1:4),ExtParticles(1)%Mass,Hel(1),ExtParticles(1)%Pol(1:4))

    ExtParticles(2)%Mom(1:4) = dcmplx(Mom(1:4,2))
    call ubarSpi(ExtParticles(2)%Mom(1:4),ExtParticles(2)%Mass,Hel(2),ExtParticles(2)%Pol(1:4))

    ExtParticles(3)%Mom(1:4) = dcmplx(Mom(1:4,3))
    call vSpi(ExtParticles(3)%Mom(1:4),ExtParticles(3)%Mass,Hel(3),ExtParticles(3)%Pol(1:4))      ! THERE WAS A BUG HERE: POL_MLESS

    ExtParticles(4)%Mom(1:4) = dcmplx(Mom(1:4,4))
    call ubarSpi(ExtParticles(4)%Mom(1:4),ExtParticles(4)%Mass,Hel(4),ExtParticles(4)%Pol(1:4))   ! THERE WAS A BUG HERE: POL_MLESS

    ExtParticles(5)%Mom(1:4) = dcmplx(Mom(1:4,5))
    call pol_mless(ExtParticles(5)%Mom(1:4),Hel(5),ExtParticles(5)%Pol(1:4))

RETURN
END SUBROUTINE


SUBROUTINE GenerateEvent5(Mom,Hel,ExtParticles)
use ModMisc
use ModProcess
implicit none
type(Particle) :: ExtParticles(:)
real(8) :: Mom(1:4,1:5),MomDK(1:4,1:6)
integer :: Hel(1:5)



    ExtParticles(1)%Mom(1:4) = dcmplx(Mom(1:4,1))
    call vSpi(ExtParticles(1)%Mom(1:4),ExtParticles(1)%Mass,Hel(1),ExtParticles(1)%Pol(1:4))

    ExtParticles(2)%Mom(1:4) = dcmplx(Mom(1:4,2))
    call ubarSpi(ExtParticles(2)%Mom(1:4),ExtParticles(2)%Mass,Hel(2),ExtParticles(2)%Pol(1:4))

    ExtParticles(3)%Mom(1:4) = dcmplx(Mom(1:4,3))
    call pol_mless(ExtParticles(3)%Mom(1:4),Hel(3),ExtParticles(3)%Pol(1:4))

    ExtParticles(4)%Mom(1:4) = dcmplx(Mom(1:4,4))
    call pol_mless(ExtParticles(4)%Mom(1:4),Hel(4),ExtParticles(4)%Pol(1:4))

    ExtParticles(5)%Mom(1:4) = dcmplx(Mom(1:4,5))
    call pol_mless(ExtParticles(5)%Mom(1:4),Hel(5),ExtParticles(5)%Pol(1:4))

RETURN
END SUBROUTINE GenerateEvent5


      subroutine cc_ttqqg(C) ! color correlation matrix for ttqqg
      real(dp), intent(out)  :: C(4,4)

        C(1,1)=24.0_dp
        C(1,2)=8.0_dp/3.0_dp
        C(1,3)=0.0_dp
        C(1,4)=8.0_dp/3.0_dp
        C(2,1)=8.0_dp/3.0_dp
        C(2,2)=8.0_dp/3.0_dp
        C(2,3)=8.0_dp/3.0_dp
        C(2,4)=0.0_dp
        C(3,1)=0.0_dp
        C(3,2)=8.0_dp/3.0_dp
        C(3,3)=24.0_dp
        C(3,4)=8.0_dp/3.0_dp
        C(4,1)=8.0_dp/3.0_dp
        C(4,2)=0.0_dp
        C(4,3)=8.0_dp/3.0_dp
        C(4,4)=8.0_dp/3.0_dp
      end subroutine cc_ttqqg


      subroutine cc_ttggg(C) ! color correlation matrix for ttggg
      real(dp), intent(out)  :: C(6,6)

      C(1,1)=512.0_dp/9.0_dp
      C(1,2)=-64.0_dp/9.0_dp
      C(1,3)=-64.0_dp/9.0_dp
      C(1,4)=8.0_dp/9.0_dp
      C(1,5)=8.0_dp/9.0_dp
      C(1,6)=80.0_dp/9.0_dp
      C(2,1)=-64.0_dp/9.0_dp
      C(2,2)=512.0_dp/9.0_dp
      C(2,3)=8.0_dp/9.0_dp
      C(2,4)=80.0_dp/9.0_dp
      C(2,5)=-64.0_dp/9.0_dp
      C(2,6)=8.0_dp/9.0_dp
      C(3,1)=-64.0_dp/9.0_dp
      C(3,2)=8.0_dp/9.0_dp
      C(3,3)=512.0_dp/9.0_dp
      C(3,4)=-64.0_dp/9.0_dp
      C(3,5)=80.0_dp/9.0_dp
      C(3,6)=8.0_dp/9.0_dp
      C(4,1)=8.0_dp/9.0_dp
      C(4,2)=80.0_dp/9.0_dp
      C(4,3)=-64.0_dp/9.0_dp
      C(4,4)=512.0_dp/9.0_dp
      C(4,5)=8.0_dp/9.0_dp
      C(4,6)=-64.0_dp/9.0_dp
      C(5,1)=8.0_dp/9.0_dp
      C(5,2)=-64.0_dp/9.0_dp
      C(5,3)=80.0_dp/9.0_dp
      C(5,4)=8.0_dp/9.0_dp
      C(5,5)=512.0_dp/9.0_dp
      C(5,6)=-64.0_dp/9.0_dp
      C(6,1)=80.0_dp/9.0_dp
      C(6,2)=8.0_dp/9.0_dp
      C(6,3)=8.0_dp/9.0_dp
      C(6,4)=-64.0_dp/9.0_dp
      C(6,5)=-64.0_dp/9.0_dp
      C(6,6)=512.0_dp/9.0_dp


      end subroutine cc_ttggg


      subroutine cc_gg_ttgggg(n,C)
      integer, intent(in)  :: n
      real(dp), intent(out)  :: C(6,6)

      if(n.eq.0) then
      C(1,1)=512.0_dp/9.0_dp
      C(1,2)=-64.0_dp/9.0_dp
      C(1,3)=-64.0_dp/9.0_dp
      C(1,4)=8.0_dp/9.0_dp
      C(1,5)=8.0_dp/9.0_dp
      C(1,6)=80.0_dp/9.0_dp
      C(2,1)=-64.0_dp/9.0_dp
      C(2,2)=512.0_dp/9.0_dp
      C(2,3)=8.0_dp/9.0_dp
      C(2,4)=80.0_dp/9.0_dp
      C(2,5)=-64.0_dp/9.0_dp
      C(2,6)=8.0_dp/9.0_dp
      C(3,1)=-64.0_dp/9.0_dp
      C(3,2)=8.0_dp/9.0_dp
      C(3,3)=512.0_dp/9.0_dp
      C(3,4)=-64.0_dp/9.0_dp
      C(3,5)=80.0_dp/9.0_dp
      C(3,6)=8.0_dp/9.0_dp
      C(4,1)=8.0_dp/9.0_dp
      C(4,2)=80.0_dp/9.0_dp
      C(4,3)=-64.0_dp/9.0_dp
      C(4,4)=512.0_dp/9.0_dp
      C(4,5)=8.0_dp/9.0_dp
      C(4,6)=-64.0_dp/9.0_dp
      C(5,1)=8.0_dp/9.0_dp
      C(5,2)=-64.0_dp/9.0_dp
      C(5,3)=80.0_dp/9.0_dp
      C(5,4)=8.0_dp/9.0_dp
      C(5,5)=512.0_dp/9.0_dp
      C(5,6)=-64.0_dp/9.0_dp
      C(6,1)=80.0_dp/9.0_dp
      C(6,2)=8.0_dp/9.0_dp
      C(6,3)=8.0_dp/9.0_dp
      C(6,4)=-64.0_dp/9.0_dp
      C(6,5)=-64.0_dp/9.0_dp
      C(6,6)=512.0_dp/9.0_dp
      endif

      if(n.eq.1) then
      C(1,1)=8.0_dp/3.0_dp;
      C(1,2)=80.0_dp/3.0_dp;
      C(1,3)=8.0_dp/3.0_dp;
      C(1,4)=8.0_dp/3.0_dp;
      C(1,5)=-64.0_dp/3.0_dp;
      C(1,6)=80.0_dp/3.0_dp;
      C(2,1)=80.0_dp/3.0_dp;
      C(2,2)=8.0_dp/3.0_dp;
      C(2,3)=-64.0_dp/3.0_dp;
      C(2,4)=80.0_dp/3.0_dp;
      C(2,5)=8.0_dp/3.0_dp;
      C(2,6)=8.0_dp/3.0_dp;
      C(3,1)=8.0_dp/3.0_dp;
      C(3,2)=-64.0_dp/3.0_dp;
      C(3,3)=-64.0_dp/3.0_dp;
      C(3,4)=-64.0_dp/3.0_dp;
      C(3,5)=-64.0_dp/3.0_dp;
      C(3,6)=8.0_dp/3.0_dp;
      C(4,1)=8.0_dp/3.0_dp;
      C(4,2)=80.0_dp/3.0_dp;
      C(4,3)=-64.0_dp/3.0_dp;
      C(4,4)=512.0_dp/3.0_dp;
      C(4,5)=8.0_dp/3.0_dp;
      C(4,6)=-64.0_dp/3.0_dp;
      C(5,1)=-64.0_dp/3.0_dp;
      C(5,2)=8.0_dp/3.0_dp;
      C(5,3)=-64.0_dp/3.0_dp;
      C(5,4)=8.0_dp/3.0_dp;
      C(5,5)=-64.0_dp/3.0_dp;
      C(5,6)=-64.0_dp/3.0_dp;
      C(6,1)=80.0_dp/3.0_dp;
      C(6,2)=8.0_dp/3.0_dp;
      C(6,3)=8.0_dp/3.0_dp;
      C(6,4)=-64.0_dp/3.0_dp;
      C(6,5)=-64.0_dp/3.0_dp;
      C(6,6)=512.0_dp/3.0_dp;
      endif
      if(n.eq.2) then
      C(1,1)=512.0_dp/3.0_dp;
      C(1,2)=-64.0_dp/3.0_dp;
      C(1,3)=-64.0_dp/3.0_dp;
      C(1,4)=8.0_dp/3.0_dp;
      C(1,5)=8.0_dp/3.0_dp;
      C(1,6)=80.0_dp/3.0_dp;
      C(2,1)=-64.0_dp/3.0_dp;
      C(2,2)=512.0_dp/3.0_dp;
      C(2,3)=8.0_dp/3.0_dp;
      C(2,4)=80.0_dp/3.0_dp;
      C(2,5)=-64.0_dp/3.0_dp;
      C(2,6)=8.0_dp/3.0_dp;
      C(3,1)=-64.0_dp/3.0_dp;
      C(3,2)=8.0_dp/3.0_dp;
      C(3,3)=-64.0_dp/3.0_dp;
      C(3,4)=8.0_dp/3.0_dp;
      C(3,5)=-64.0_dp/3.0_dp;
      C(3,6)=-64.0_dp/3.0_dp;
      C(4,1)=8.0_dp/3.0_dp;
      C(4,2)=80.0_dp/3.0_dp;
      C(4,3)=8.0_dp/3.0_dp;
      C(4,4)=8.0_dp/3.0_dp;
      C(4,5)=-64.0_dp/3.0_dp;
      C(4,6)=80.0_dp/3.0_dp;
      C(5,1)=8.0_dp/3.0_dp;
      C(5,2)=-64.0_dp/3.0_dp;
      C(5,3)=-64.0_dp/3.0_dp;
      C(5,4)=-64.0_dp/3.0_dp;
      C(5,5)=-64.0_dp/3.0_dp;
      C(5,6)=8.0_dp/3.0_dp;
      C(6,1)=80.0_dp/3.0_dp;
      C(6,2)=8.0_dp/3.0_dp;
      C(6,3)=-64.0_dp/3.0_dp;
      C(6,4)=80.0_dp/3.0_dp;
      C(6,5)=8.0_dp/3.0_dp;
      C(6,6)=8.0_dp/3.0_dp;
      endif
      if(n.eq.3) then
      C(1,1)=192.0_dp;
      C(1,2)=-24.0_dp;
      C(1,3)=0.0_dp;
      C(1,4)=0.0_dp;
      C(1,5)=24.0_dp;
      C(1,6)=48.0_dp;
      C(2,1)=-24.0_dp;
      C(2,2)=-24.0_dp;
      C(2,3)=0.0_dp;
      C(2,4)=-48.0_dp;
      C(2,5)=-24.0_dp;
      C(2,6)=0.0_dp;
      C(3,1)=0.0_dp;
      C(3,2)=0.0_dp;
      C(3,3)=192.0_dp;
      C(3,4)=-24.0_dp;
      C(3,5)=48.0_dp;
      C(3,6)=24.0_dp;
      C(4,1)=0.0_dp;
      C(4,2)=-48.0_dp;
      C(4,3)=-24.0_dp;
      C(4,4)=-24.0_dp;
      C(4,5)=0.0_dp;
      C(4,6)=-24.0_dp;
      C(5,1)=24.0_dp;
      C(5,2)=-24.0_dp;
      C(5,3)=48.0_dp;
      C(5,4)=0.0_dp;
      C(5,5)=192.0_dp;
      C(5,6)=0.0_dp;
      C(6,1)=48.0_dp;
      C(6,2)=0.0_dp;
      C(6,3)=24.0_dp;
      C(6,4)=-24.0_dp;
      C(6,5)=0.0_dp;
      C(6,6)=192.0_dp;
      endif
      if(n.eq.4) then
      C(1,1)=-24.0_dp;
      C(1,2)=-24.0_dp;
      C(1,3)=-24.0_dp;
      C(1,4)=0.0_dp;
      C(1,5)=0.0_dp;
      C(1,6)=-48.0_dp;
      C(2,1)=-24.0_dp;
      C(2,2)=192.0_dp;
      C(2,3)=24.0_dp;
      C(2,4)=48.0_dp;
      C(2,5)=0.0_dp;
      C(2,6)=0.0_dp;
      C(3,1)=-24.0_dp;
      C(3,2)=24.0_dp;
      C(3,3)=192.0_dp;
      C(3,4)=0.0_dp;
      C(3,5)=48.0_dp;
      C(3,6)=0.0_dp;
      C(4,1)=0.0_dp;
      C(4,2)=48.0_dp;
      C(4,3)=0.0_dp;
      C(4,4)=192.0_dp;
      C(4,5)=24.0_dp;
      C(4,6)=-24.0_dp;
      C(5,1)=0.0_dp;
      C(5,2)=0.0_dp;
      C(5,3)=48.0_dp;
      C(5,4)=24.0_dp;
      C(5,5)=192.0_dp;
      C(5,6)=-24.0_dp;
      C(6,1)=-48.0_dp;
      C(6,2)=0.0_dp;
      C(6,3)=0.0_dp;
      C(6,4)=-24.0_dp;
      C(6,5)=-24.0_dp;
      C(6,6)=-24.0_dp;
      endif
      if(n.eq.5) then
      C(1,1)=-64.0_dp/3.0_dp;
      C(1,2)=-64.0_dp/3.0_dp;
      C(1,3)=8.0_dp/3.0_dp;
      C(1,4)=-64.0_dp/3.0_dp;
      C(1,5)=8.0_dp/3.0_dp;
      C(1,6)=-64.0_dp/3.0_dp;
      C(2,1)=-64.0_dp/3.0_dp;
      C(2,2)=512.0_dp/3.0_dp;
      C(2,3)=8.0_dp/3.0_dp;
      C(2,4)=80.0_dp/3.0_dp;
      C(2,5)=-64.0_dp/3.0_dp;
      C(2,6)=8.0_dp/3.0_dp;
      C(3,1)=8.0_dp/3.0_dp;
      C(3,2)=8.0_dp/3.0_dp;
      C(3,3)=8.0_dp/3.0_dp;
      C(3,4)=80.0_dp/3.0_dp;
      C(3,5)=80.0_dp/3.0_dp;
      C(3,6)=-64.0_dp/3.0_dp;
      C(4,1)=-64.0_dp/3.0_dp;
      C(4,2)=80.0_dp/3.0_dp;
      C(4,3)=80.0_dp/3.0_dp;
      C(4,4)=8.0_dp/3.0_dp;
      C(4,5)=8.0_dp/3.0_dp;
      C(4,6)=8.0_dp/3.0_dp;
      C(5,1)=8.0_dp/3.0_dp;
      C(5,2)=-64.0_dp/3.0_dp;
      C(5,3)=80.0_dp/3.0_dp;
      C(5,4)=8.0_dp/3.0_dp;
      C(5,5)=512.0_dp/3.0_dp;
      C(5,6)=-64.0_dp/3.0_dp;
      C(6,1)=-64.0_dp/3.0_dp;
      C(6,2)=8.0_dp/3.0_dp;
      C(6,3)=-64.0_dp/3.0_dp;
      C(6,4)=8.0_dp/3.0_dp;
      C(6,5)=-64.0_dp/3.0_dp;
      C(6,6)=-64.0_dp/3.0_dp;
      endif
      if(n.eq.6) then
      C(1,1)=-64.0_dp/3.0_dp;
      C(1,2)=8.0_dp/3.0_dp;
      C(1,3)=-64.0_dp/3.0_dp;
      C(1,4)=8.0_dp/3.0_dp;
      C(1,5)=-64.0_dp/3.0_dp;
      C(1,6)=-64.0_dp/3.0_dp;
      C(2,1)=8.0_dp/3.0_dp;
      C(2,2)=8.0_dp/3.0_dp;
      C(2,3)=8.0_dp/3.0_dp;
      C(2,4)=80.0_dp/3.0_dp;
      C(2,5)=80.0_dp/3.0_dp;
      C(2,6)=-64.0_dp/3.0_dp;
      C(3,1)=-64.0_dp/3.0_dp;
      C(3,2)=8.0_dp/3.0_dp;
      C(3,3)=512.0_dp/3.0_dp;
      C(3,4)=-64.0_dp/3.0_dp;
      C(3,5)=80.0_dp/3.0_dp;
      C(3,6)=8.0_dp/3.0_dp;
      C(4,1)=8.0_dp/3.0_dp;
      C(4,2)=80.0_dp/3.0_dp;
      C(4,3)=-64.0_dp/3.0_dp;
      C(4,4)=512.0_dp/3.0_dp;
      C(4,5)=8.0_dp/3.0_dp;
      C(4,6)=-64.0_dp/3.0_dp;
      C(5,1)=-64.0_dp/3.0_dp;
      C(5,2)=80.0_dp/3.0_dp;
      C(5,3)=80.0_dp/3.0_dp;
      C(5,4)=8.0_dp/3.0_dp;
      C(5,5)=8.0_dp/3.0_dp;
      C(5,6)=8.0_dp/3.0_dp;
      C(6,1)=-64.0_dp/3.0_dp;
      C(6,2)=-64.0_dp/3.0_dp;
      C(6,3)=8.0_dp/3.0_dp;
      C(6,4)=-64.0_dp/3.0_dp;
      C(6,5)=8.0_dp/3.0_dp;
      C(6,6)=-64.0_dp/3.0_dp;
      endif
      if(n.eq.7) then
      C(1,1)=192.0_dp;
      C(1,2)=-24.0_dp;
      C(1,3)=0.0_dp;
      C(1,4)=0.0_dp;
      C(1,5)=24.0_dp;
      C(1,6)=48.0_dp;
      C(2,1)=-24.0_dp;
      C(2,2)=-24.0_dp;
      C(2,3)=0.0_dp;
      C(2,4)=-48.0_dp;
      C(2,5)=-24.0_dp;
      C(2,6)=0.0_dp;
      C(3,1)=0.0_dp;
      C(3,2)=0.0_dp;
      C(3,3)=192.0_dp;
      C(3,4)=-24.0_dp;
      C(3,5)=48.0_dp;
      C(3,6)=24.0_dp;
      C(4,1)=0.0_dp;
      C(4,2)=-48.0_dp;
      C(4,3)=-24.0_dp;
      C(4,4)=-24.0_dp;
      C(4,5)=0.0_dp;
      C(4,6)=-24.0_dp;
      C(5,1)=24.0_dp;
      C(5,2)=-24.0_dp;
      C(5,3)=48.0_dp;
      C(5,4)=0.0_dp;
      C(5,5)=192.0_dp;
      C(5,6)=0.0_dp;
      C(6,1)=48.0_dp;
      C(6,2)=0.0_dp;
      C(6,3)=24.0_dp;
      C(6,4)=-24.0_dp;
      C(6,5)=0.0_dp;
      C(6,6)=192.0_dp;
      endif
      if(n.eq.8) then
      C(1,1)=192.0_dp;
      C(1,2)=0.0_dp;
      C(1,3)=-24.0_dp;
      C(1,4)=24.0_dp;
      C(1,5)=0.0_dp;
      C(1,6)=48.0_dp;
      C(2,1)=0.0_dp;
      C(2,2)=192.0_dp;
      C(2,3)=0.0_dp;
      C(2,4)=48.0_dp;
      C(2,5)=-24.0_dp;
      C(2,6)=24.0_dp;
      C(3,1)=-24.0_dp;
      C(3,2)=0.0_dp;
      C(3,3)=-24.0_dp;
      C(3,4)=-24.0_dp;
      C(3,5)=-48.0_dp;
      C(3,6)=0.0_dp;
      C(4,1)=24.0_dp;
      C(4,2)=48.0_dp;
      C(4,3)=-24.0_dp;
      C(4,4)=192.0_dp;
      C(4,5)=0.0_dp;
      C(4,6)=0.0_dp;
      C(5,1)=0.0_dp;
      C(5,2)=-24.0_dp;
      C(5,3)=-48.0_dp;
      C(5,4)=0.0_dp;
      C(5,5)=-24.0_dp;
      C(5,6)=-24.0_dp;
      C(6,1)=48.0_dp;
      C(6,2)=24.0_dp;
      C(6,3)=0.0_dp;
      C(6,4)=0.0_dp;
      C(6,5)=-24.0_dp;
      C(6,6)=192.0_dp;
      endif
      if(n.eq.9) then
      C(1,1)=-8.0_dp/27.0_dp;
      C(1,2)=-80.0_dp/27.0_dp;
      C(1,3)=-80.0_dp/27.0_dp;
      C(1,4)=496.0_dp/27.0_dp;
      C(1,5)=496.0_dp/27.0_dp;
      C(1,6)=-224.0_dp/27.0_dp;
      C(2,1)=-80.0_dp/27.0_dp;
      C(2,2)=-8.0_dp/27.0_dp;
      C(2,3)=496.0_dp/27.0_dp;
      C(2,4)=-224.0_dp/27.0_dp;
      C(2,5)=-80.0_dp/27.0_dp;
      C(2,6)=496.0_dp/27.0_dp;
      C(3,1)=-80.0_dp/27.0_dp;
      C(3,2)=496.0_dp/27.0_dp;
      C(3,3)=-8.0_dp/27.0_dp;
      C(3,4)=-80.0_dp/27.0_dp;
      C(3,5)=-224.0_dp/27.0_dp;
      C(3,6)=496.0_dp/27.0_dp;
      C(4,1)=496.0_dp/27.0_dp;
      C(4,2)=-224.0_dp/27.0_dp;
      C(4,3)=-80.0_dp/27.0_dp;
      C(4,4)=-8.0_dp/27.0_dp;
      C(4,5)=496.0_dp/27.0_dp;
      C(4,6)=-80.0_dp/27.0_dp;
      C(5,1)=496.0_dp/27.0_dp;
      C(5,2)=-80.0_dp/27.0_dp;
      C(5,3)=-224.0_dp/27.0_dp;
      C(5,4)=496.0_dp/27.0_dp;
      C(5,5)=-8.0_dp/27.0_dp;
      C(5,6)=-80.0_dp/27.0_dp;
      C(6,1)=-224.0_dp/27.0_dp;
      C(6,2)=496.0_dp/27.0_dp;
      C(6,3)=496.0_dp/27.0_dp;
      C(6,4)=-80.0_dp/27.0_dp;
      C(6,5)=-80.0_dp/27.0_dp;
      C(6,6)=-8.0_dp/27.0_dp;
      endif
      if(n.eq.10) then
      C(1,1)=8.0_dp/3.0_dp;
      C(1,2)=80.0_dp/3.0_dp;
      C(1,3)=8.0_dp/3.0_dp;
      C(1,4)=8.0_dp/3.0_dp;
      C(1,5)=-64.0_dp/3.0_dp;
      C(1,6)=80.0_dp/3.0_dp;
      C(2,1)=80.0_dp/3.0_dp;
      C(2,2)=8.0_dp/3.0_dp;
      C(2,3)=-64.0_dp/3.0_dp;
      C(2,4)=80.0_dp/3.0_dp;
      C(2,5)=8.0_dp/3.0_dp;
      C(2,6)=8.0_dp/3.0_dp;
      C(3,1)=8.0_dp/3.0_dp;
      C(3,2)=-64.0_dp/3.0_dp;
      C(3,3)=-64.0_dp/3.0_dp;
      C(3,4)=-64.0_dp/3.0_dp;
      C(3,5)=-64.0_dp/3.0_dp;
      C(3,6)=8.0_dp/3.0_dp;
      C(4,1)=8.0_dp/3.0_dp;
      C(4,2)=80.0_dp/3.0_dp;
      C(4,3)=-64.0_dp/3.0_dp;
      C(4,4)=512.0_dp/3.0_dp;
      C(4,5)=8.0_dp/3.0_dp;
      C(4,6)=-64.0_dp/3.0_dp;
      C(5,1)=-64.0_dp/3.0_dp;
      C(5,2)=8.0_dp/3.0_dp;
      C(5,3)=-64.0_dp/3.0_dp;
      C(5,4)=8.0_dp/3.0_dp;
      C(5,5)=-64.0_dp/3.0_dp;
      C(5,6)=-64.0_dp/3.0_dp;
      C(6,1)=80.0_dp/3.0_dp;
      C(6,2)=8.0_dp/3.0_dp;
      C(6,3)=8.0_dp/3.0_dp;
      C(6,4)=-64.0_dp/3.0_dp;
      C(6,5)=-64.0_dp/3.0_dp;
      C(6,6)=512.0_dp/3.0_dp;
      endif
      if(n.eq.11) then
      C(1,1)=-64.0_dp/3.0_dp;
      C(1,2)=-64.0_dp/3.0_dp;
      C(1,3)=8.0_dp/3.0_dp;
      C(1,4)=-64.0_dp/3.0_dp;
      C(1,5)=8.0_dp/3.0_dp;
      C(1,6)=-64.0_dp/3.0_dp;
      C(2,1)=-64.0_dp/3.0_dp;
      C(2,2)=512.0_dp/3.0_dp;
      C(2,3)=8.0_dp/3.0_dp;
      C(2,4)=80.0_dp/3.0_dp;
      C(2,5)=-64.0_dp/3.0_dp;
      C(2,6)=8.0_dp/3.0_dp;
      C(3,1)=8.0_dp/3.0_dp;
      C(3,2)=8.0_dp/3.0_dp;
      C(3,3)=8.0_dp/3.0_dp;
      C(3,4)=80.0_dp/3.0_dp;
      C(3,5)=80.0_dp/3.0_dp;
      C(3,6)=-64.0_dp/3.0_dp;
      C(4,1)=-64.0_dp/3.0_dp;
      C(4,2)=80.0_dp/3.0_dp;
      C(4,3)=80.0_dp/3.0_dp;
      C(4,4)=8.0_dp/3.0_dp;
      C(4,5)=8.0_dp/3.0_dp;
      C(4,6)=8.0_dp/3.0_dp;
      C(5,1)=8.0_dp/3.0_dp;
      C(5,2)=-64.0_dp/3.0_dp;
      C(5,3)=80.0_dp/3.0_dp;
      C(5,4)=8.0_dp/3.0_dp;
      C(5,5)=512.0_dp/3.0_dp;
      C(5,6)=-64.0_dp/3.0_dp;
      C(6,1)=-64.0_dp/3.0_dp;
      C(6,2)=8.0_dp/3.0_dp;
      C(6,3)=-64.0_dp/3.0_dp;
      C(6,4)=8.0_dp/3.0_dp;
      C(6,5)=-64.0_dp/3.0_dp;
      C(6,6)=-64.0_dp/3.0_dp;
      endif
      if(n.eq.12) then
      C(1,1)=512.0_dp/3.0_dp;
      C(1,2)=-64.0_dp/3.0_dp;
      C(1,3)=-64.0_dp/3.0_dp;
      C(1,4)=8.0_dp/3.0_dp;
      C(1,5)=8.0_dp/3.0_dp;
      C(1,6)=80.0_dp/3.0_dp;
      C(2,1)=-64.0_dp/3.0_dp;
      C(2,2)=-64.0_dp/3.0_dp;
      C(2,3)=8.0_dp/3.0_dp;
      C(2,4)=-64.0_dp/3.0_dp;
      C(2,5)=8.0_dp/3.0_dp;
      C(2,6)=-64.0_dp/3.0_dp;
      C(3,1)=-64.0_dp/3.0_dp;
      C(3,2)=8.0_dp/3.0_dp;
      C(3,3)=512.0_dp/3.0_dp;
      C(3,4)=-64.0_dp/3.0_dp;
      C(3,5)=80.0_dp/3.0_dp;
      C(3,6)=8.0_dp/3.0_dp;
      C(4,1)=8.0_dp/3.0_dp;
      C(4,2)=-64.0_dp/3.0_dp;
      C(4,3)=-64.0_dp/3.0_dp;
      C(4,4)=-64.0_dp/3.0_dp;
      C(4,5)=-64.0_dp/3.0_dp;
      C(4,6)=8.0_dp/3.0_dp;
      C(5,1)=8.0_dp/3.0_dp;
      C(5,2)=8.0_dp/3.0_dp;
      C(5,3)=80.0_dp/3.0_dp;
      C(5,4)=-64.0_dp/3.0_dp;
      C(5,5)=8.0_dp/3.0_dp;
      C(5,6)=80.0_dp/3.0_dp;
      C(6,1)=80.0_dp/3.0_dp;
      C(6,2)=-64.0_dp/3.0_dp;
      C(6,3)=8.0_dp/3.0_dp;
      C(6,4)=8.0_dp/3.0_dp;
      C(6,5)=80.0_dp/3.0_dp;
      C(6,6)=8.0_dp/3.0_dp;
      endif
      if(n.eq.13) then
      C(1,1)=-8.0_dp/27.0_dp;
      C(1,2)=-80.0_dp/27.0_dp;
      C(1,3)=-80.0_dp/27.0_dp;
      C(1,4)=496.0_dp/27.0_dp;
      C(1,5)=496.0_dp/27.0_dp;
      C(1,6)=-224.0_dp/27.0_dp;
      C(2,1)=-80.0_dp/27.0_dp;
      C(2,2)=-8.0_dp/27.0_dp;
      C(2,3)=496.0_dp/27.0_dp;
      C(2,4)=-224.0_dp/27.0_dp;
      C(2,5)=-80.0_dp/27.0_dp;
      C(2,6)=496.0_dp/27.0_dp;
      C(3,1)=-80.0_dp/27.0_dp;
      C(3,2)=496.0_dp/27.0_dp;
      C(3,3)=-8.0_dp/27.0_dp;
      C(3,4)=-80.0_dp/27.0_dp;
      C(3,5)=-224.0_dp/27.0_dp;
      C(3,6)=496.0_dp/27.0_dp;
      C(4,1)=496.0_dp/27.0_dp;
      C(4,2)=-224.0_dp/27.0_dp;
      C(4,3)=-80.0_dp/27.0_dp;
      C(4,4)=-8.0_dp/27.0_dp;
      C(4,5)=496.0_dp/27.0_dp;
      C(4,6)=-80.0_dp/27.0_dp;
      C(5,1)=496.0_dp/27.0_dp;
      C(5,2)=-80.0_dp/27.0_dp;
      C(5,3)=-224.0_dp/27.0_dp;
      C(5,4)=496.0_dp/27.0_dp;
      C(5,5)=-8.0_dp/27.0_dp;
      C(5,6)=-80.0_dp/27.0_dp;
      C(6,1)=-224.0_dp/27.0_dp;
      C(6,2)=496.0_dp/27.0_dp;
      C(6,3)=496.0_dp/27.0_dp;
      C(6,4)=-80.0_dp/27.0_dp;
      C(6,5)=-80.0_dp/27.0_dp;
      C(6,6)=-8.0_dp/27.0_dp;
      endif
      if(n.eq.14) then
      C(1,1)=512.0_dp/3.0_dp;
      C(1,2)=-64.0_dp/3.0_dp;
      C(1,3)=-64.0_dp/3.0_dp;
      C(1,4)=8.0_dp/3.0_dp;
      C(1,5)=8.0_dp/3.0_dp;
      C(1,6)=80.0_dp/3.0_dp;
      C(2,1)=-64.0_dp/3.0_dp;
      C(2,2)=512.0_dp/3.0_dp;
      C(2,3)=8.0_dp/3.0_dp;
      C(2,4)=80.0_dp/3.0_dp;
      C(2,5)=-64.0_dp/3.0_dp;
      C(2,6)=8.0_dp/3.0_dp;
      C(3,1)=-64.0_dp/3.0_dp;
      C(3,2)=8.0_dp/3.0_dp;
      C(3,3)=-64.0_dp/3.0_dp;
      C(3,4)=8.0_dp/3.0_dp;
      C(3,5)=-64.0_dp/3.0_dp;
      C(3,6)=-64.0_dp/3.0_dp;
      C(4,1)=8.0_dp/3.0_dp;
      C(4,2)=80.0_dp/3.0_dp;
      C(4,3)=8.0_dp/3.0_dp;
      C(4,4)=8.0_dp/3.0_dp;
      C(4,5)=-64.0_dp/3.0_dp;
      C(4,6)=80.0_dp/3.0_dp;
      C(5,1)=8.0_dp/3.0_dp;
      C(5,2)=-64.0_dp/3.0_dp;
      C(5,3)=-64.0_dp/3.0_dp;
      C(5,4)=-64.0_dp/3.0_dp;
      C(5,5)=-64.0_dp/3.0_dp;
      C(5,6)=8.0_dp/3.0_dp;
      C(6,1)=80.0_dp/3.0_dp;
      C(6,2)=8.0_dp/3.0_dp;
      C(6,3)=-64.0_dp/3.0_dp;
      C(6,4)=80.0_dp/3.0_dp;
      C(6,5)=8.0_dp/3.0_dp;
      C(6,6)=8.0_dp/3.0_dp;
      endif
      if(n.eq.15) then
      C(1,1)=-64.0_dp/3.0_dp;
      C(1,2)=8.0_dp/3.0_dp;
      C(1,3)=-64.0_dp/3.0_dp;
      C(1,4)=8.0_dp/3.0_dp;
      C(1,5)=-64.0_dp/3.0_dp;
      C(1,6)=-64.0_dp/3.0_dp;
      C(2,1)=8.0_dp/3.0_dp;
      C(2,2)=8.0_dp/3.0_dp;
      C(2,3)=8.0_dp/3.0_dp;
      C(2,4)=80.0_dp/3.0_dp;
      C(2,5)=80.0_dp/3.0_dp;
      C(2,6)=-64.0_dp/3.0_dp;
      C(3,1)=-64.0_dp/3.0_dp;
      C(3,2)=8.0_dp/3.0_dp;
      C(3,3)=512.0_dp/3.0_dp;
      C(3,4)=-64.0_dp/3.0_dp;
      C(3,5)=80.0_dp/3.0_dp;
      C(3,6)=8.0_dp/3.0_dp;
      C(4,1)=8.0_dp/3.0_dp;
      C(4,2)=80.0_dp/3.0_dp;
      C(4,3)=-64.0_dp/3.0_dp;
      C(4,4)=512.0_dp/3.0_dp;
      C(4,5)=8.0_dp/3.0_dp;
      C(4,6)=-64.0_dp/3.0_dp;
      C(5,1)=-64.0_dp/3.0_dp;
      C(5,2)=80.0_dp/3.0_dp;
      C(5,3)=80.0_dp/3.0_dp;
      C(5,4)=8.0_dp/3.0_dp;
      C(5,5)=8.0_dp/3.0_dp;
      C(5,6)=8.0_dp/3.0_dp;
      C(6,1)=-64.0_dp/3.0_dp;
      C(6,2)=-64.0_dp/3.0_dp;
      C(6,3)=8.0_dp/3.0_dp;
      C(6,4)=-64.0_dp/3.0_dp;
      C(6,5)=8.0_dp/3.0_dp;
      C(6,6)=-64.0_dp/3.0_dp;
      endif
      if(n.eq.16) then
      C(1,1)=8.0_dp/3.0_dp;
      C(1,2)=8.0_dp/3.0_dp;
      C(1,3)=80.0_dp/3.0_dp;
      C(1,4)=-64.0_dp/3.0_dp;
      C(1,5)=8.0_dp/3.0_dp;
      C(1,6)=80.0_dp/3.0_dp;
      C(2,1)=8.0_dp/3.0_dp;
      C(2,2)=-64.0_dp/3.0_dp;
      C(2,3)=-64.0_dp/3.0_dp;
      C(2,4)=-64.0_dp/3.0_dp;
      C(2,5)=-64.0_dp/3.0_dp;
      C(2,6)=8.0_dp/3.0_dp;
      C(3,1)=80.0_dp/3.0_dp;
      C(3,2)=-64.0_dp/3.0_dp;
      C(3,3)=8.0_dp/3.0_dp;
      C(3,4)=8.0_dp/3.0_dp;
      C(3,5)=80.0_dp/3.0_dp;
      C(3,6)=8.0_dp/3.0_dp;
      C(4,1)=-64.0_dp/3.0_dp;
      C(4,2)=-64.0_dp/3.0_dp;
      C(4,3)=8.0_dp/3.0_dp;
      C(4,4)=-64.0_dp/3.0_dp;
      C(4,5)=8.0_dp/3.0_dp;
      C(4,6)=-64.0_dp/3.0_dp;
      C(5,1)=8.0_dp/3.0_dp;
      C(5,2)=-64.0_dp/3.0_dp;
      C(5,3)=80.0_dp/3.0_dp;
      C(5,4)=8.0_dp/3.0_dp;
      C(5,5)=512.0_dp/3.0_dp;
      C(5,6)=-64.0_dp/3.0_dp;
      C(6,1)=80.0_dp/3.0_dp;
      C(6,2)=8.0_dp/3.0_dp;
      C(6,3)=8.0_dp/3.0_dp;
      C(6,4)=-64.0_dp/3.0_dp;
      C(6,5)=-64.0_dp/3.0_dp;
      C(6,6)=512.0_dp/3.0_dp;
      endif
      if(n.eq.17) then
      C(1,1)=512.0_dp/3.0_dp;
      C(1,2)=-64.0_dp/3.0_dp;
      C(1,3)=-64.0_dp/3.0_dp;
      C(1,4)=8.0_dp/3.0_dp;
      C(1,5)=8.0_dp/3.0_dp;
      C(1,6)=80.0_dp/3.0_dp;
      C(2,1)=-64.0_dp/3.0_dp;
      C(2,2)=-64.0_dp/3.0_dp;
      C(2,3)=8.0_dp/3.0_dp;
      C(2,4)=-64.0_dp/3.0_dp;
      C(2,5)=8.0_dp/3.0_dp;
      C(2,6)=-64.0_dp/3.0_dp;
      C(3,1)=-64.0_dp/3.0_dp;
      C(3,2)=8.0_dp/3.0_dp;
      C(3,3)=512.0_dp/3.0_dp;
      C(3,4)=-64.0_dp/3.0_dp;
      C(3,5)=80.0_dp/3.0_dp;
      C(3,6)=8.0_dp/3.0_dp;
      C(4,1)=8.0_dp/3.0_dp;
      C(4,2)=-64.0_dp/3.0_dp;
      C(4,3)=-64.0_dp/3.0_dp;
      C(4,4)=-64.0_dp/3.0_dp;
      C(4,5)=-64.0_dp/3.0_dp;
      C(4,6)=8.0_dp/3.0_dp;
      C(5,1)=8.0_dp/3.0_dp;
      C(5,2)=8.0_dp/3.0_dp;
      C(5,3)=80.0_dp/3.0_dp;
      C(5,4)=-64.0_dp/3.0_dp;
      C(5,5)=8.0_dp/3.0_dp;
      C(5,6)=80.0_dp/3.0_dp;
      C(6,1)=80.0_dp/3.0_dp;
      C(6,2)=-64.0_dp/3.0_dp;
      C(6,3)=8.0_dp/3.0_dp;
      C(6,4)=8.0_dp/3.0_dp;
      C(6,5)=80.0_dp/3.0_dp;
      C(6,6)=8.0_dp/3.0_dp;
      endif
      if(n.eq.18) then
      C(1,1)=8.0_dp/3.0_dp;
      C(1,2)=8.0_dp/3.0_dp;
      C(1,3)=80.0_dp/3.0_dp;
      C(1,4)=-64.0_dp/3.0_dp;
      C(1,5)=8.0_dp/3.0_dp;
      C(1,6)=80.0_dp/3.0_dp;
      C(2,1)=8.0_dp/3.0_dp;
      C(2,2)=-64.0_dp/3.0_dp;
      C(2,3)=-64.0_dp/3.0_dp;
      C(2,4)=-64.0_dp/3.0_dp;
      C(2,5)=-64.0_dp/3.0_dp;
      C(2,6)=8.0_dp/3.0_dp;
      C(3,1)=80.0_dp/3.0_dp;
      C(3,2)=-64.0_dp/3.0_dp;
      C(3,3)=8.0_dp/3.0_dp;
      C(3,4)=8.0_dp/3.0_dp;
      C(3,5)=80.0_dp/3.0_dp;
      C(3,6)=8.0_dp/3.0_dp;
      C(4,1)=-64.0_dp/3.0_dp;
      C(4,2)=-64.0_dp/3.0_dp;
      C(4,3)=8.0_dp/3.0_dp;
      C(4,4)=-64.0_dp/3.0_dp;
      C(4,5)=8.0_dp/3.0_dp;
      C(4,6)=-64.0_dp/3.0_dp;
      C(5,1)=8.0_dp/3.0_dp;
      C(5,2)=-64.0_dp/3.0_dp;
      C(5,3)=80.0_dp/3.0_dp;
      C(5,4)=8.0_dp/3.0_dp;
      C(5,5)=512.0_dp/3.0_dp;
      C(5,6)=-64.0_dp/3.0_dp;
      C(6,1)=80.0_dp/3.0_dp;
      C(6,2)=8.0_dp/3.0_dp;
      C(6,3)=8.0_dp/3.0_dp;
      C(6,4)=-64.0_dp/3.0_dp;
      C(6,5)=-64.0_dp/3.0_dp;
      C(6,6)=512.0_dp/3.0_dp;
      endif
      if(n.eq.19) then
      C(1,1)=-24.0_dp;
      C(1,2)=-24.0_dp;
      C(1,3)=-24.0_dp;
      C(1,4)=0.0_dp;
      C(1,5)=0.0_dp;
      C(1,6)=-48.0_dp;
      C(2,1)=-24.0_dp;
      C(2,2)=192.0_dp;
      C(2,3)=24.0_dp;
      C(2,4)=48.0_dp;
      C(2,5)=0.0_dp;
      C(2,6)=0.0_dp;
      C(3,1)=-24.0_dp;
      C(3,2)=24.0_dp;
      C(3,3)=192.0_dp;
      C(3,4)=0.0_dp;
      C(3,5)=48.0_dp;
      C(3,6)=0.0_dp;
      C(4,1)=0.0_dp;
      C(4,2)=48.0_dp;
      C(4,3)=0.0_dp;
      C(4,4)=192.0_dp;
      C(4,5)=24.0_dp;
      C(4,6)=-24.0_dp;
      C(5,1)=0.0_dp;
      C(5,2)=0.0_dp;
      C(5,3)=48.0_dp;
      C(5,4)=24.0_dp;
      C(5,5)=192.0_dp;
      C(5,6)=-24.0_dp;
      C(6,1)=-48.0_dp;
      C(6,2)=0.0_dp;
      C(6,3)=0.0_dp;
      C(6,4)=-24.0_dp;
      C(6,5)=-24.0_dp;
      C(6,6)=-24.0_dp;
      endif
      if(n.eq.20) then
      C(1,1)=192.0_dp;
      C(1,2)=0.0_dp;
      C(1,3)=-24.0_dp;
      C(1,4)=24.0_dp;
      C(1,5)=0.0_dp;
      C(1,6)=48.0_dp;
      C(2,1)=0.0_dp;
      C(2,2)=192.0_dp;
      C(2,3)=0.0_dp;
      C(2,4)=48.0_dp;
      C(2,5)=-24.0_dp;
      C(2,6)=24.0_dp;
      C(3,1)=-24.0_dp;
      C(3,2)=0.0_dp;
      C(3,3)=-24.0_dp;
      C(3,4)=-24.0_dp;
      C(3,5)=-48.0_dp;
      C(3,6)=0.0_dp;
      C(4,1)=24.0_dp;
      C(4,2)=48.0_dp;
      C(4,3)=-24.0_dp;
      C(4,4)=192.0_dp;
      C(4,5)=0.0_dp;
      C(4,6)=0.0_dp;
      C(5,1)=0.0_dp;
      C(5,2)=-24.0_dp;
      C(5,3)=-48.0_dp;
      C(5,4)=0.0_dp;
      C(5,5)=-24.0_dp;
      C(5,6)=-24.0_dp;
      C(6,1)=48.0_dp;
      C(6,2)=24.0_dp;
      C(6,3)=0.0_dp;
      C(6,4)=0.0_dp;
      C(6,5)=-24.0_dp;
      C(6,6)=192.0_dp;
      endif
      if(n.eq.21) then
      C(1,1)=8.0_dp/3.0_dp;
      C(1,2)=80.0_dp/3.0_dp;
      C(1,3)=8.0_dp/3.0_dp;
      C(1,4)=8.0_dp/3.0_dp;
      C(1,5)=-64.0_dp/3.0_dp;
      C(1,6)=80.0_dp/3.0_dp;
      C(2,1)=80.0_dp/3.0_dp;
      C(2,2)=8.0_dp/3.0_dp;
      C(2,3)=-64.0_dp/3.0_dp;
      C(2,4)=80.0_dp/3.0_dp;
      C(2,5)=8.0_dp/3.0_dp;
      C(2,6)=8.0_dp/3.0_dp;
      C(3,1)=8.0_dp/3.0_dp;
      C(3,2)=-64.0_dp/3.0_dp;
      C(3,3)=-64.0_dp/3.0_dp;
      C(3,4)=-64.0_dp/3.0_dp;
      C(3,5)=-64.0_dp/3.0_dp;
      C(3,6)=8.0_dp/3.0_dp;
      C(4,1)=8.0_dp/3.0_dp;
      C(4,2)=80.0_dp/3.0_dp;
      C(4,3)=-64.0_dp/3.0_dp;
      C(4,4)=512.0_dp/3.0_dp;
      C(4,5)=8.0_dp/3.0_dp;
      C(4,6)=-64.0_dp/3.0_dp;
      C(5,1)=-64.0_dp/3.0_dp;
      C(5,2)=8.0_dp/3.0_dp;
      C(5,3)=-64.0_dp/3.0_dp;
      C(5,4)=8.0_dp/3.0_dp;
      C(5,5)=-64.0_dp/3.0_dp;
      C(5,6)=-64.0_dp/3.0_dp;
      C(6,1)=80.0_dp/3.0_dp;
      C(6,2)=8.0_dp/3.0_dp;
      C(6,3)=8.0_dp/3.0_dp;
      C(6,4)=-64.0_dp/3.0_dp;
      C(6,5)=-64.0_dp/3.0_dp;
      C(6,6)=512.0_dp/3.0_dp;
      endif
      if(n.eq.22) then
      C(1,1)=512.0_dp/3.0_dp;
      C(1,2)=-64.0_dp/3.0_dp;
      C(1,3)=-64.0_dp/3.0_dp;
      C(1,4)=8.0_dp/3.0_dp;
      C(1,5)=8.0_dp/3.0_dp;
      C(1,6)=80.0_dp/3.0_dp;
      C(2,1)=-64.0_dp/3.0_dp;
      C(2,2)=512.0_dp/3.0_dp;
      C(2,3)=8.0_dp/3.0_dp;
      C(2,4)=80.0_dp/3.0_dp;
      C(2,5)=-64.0_dp/3.0_dp;
      C(2,6)=8.0_dp/3.0_dp;
      C(3,1)=-64.0_dp/3.0_dp;
      C(3,2)=8.0_dp/3.0_dp;
      C(3,3)=-64.0_dp/3.0_dp;
      C(3,4)=8.0_dp/3.0_dp;
      C(3,5)=-64.0_dp/3.0_dp;
      C(3,6)=-64.0_dp/3.0_dp;
      C(4,1)=8.0_dp/3.0_dp;
      C(4,2)=80.0_dp/3.0_dp;
      C(4,3)=8.0_dp/3.0_dp;
      C(4,4)=8.0_dp/3.0_dp;
      C(4,5)=-64.0_dp/3.0_dp;
      C(4,6)=80.0_dp/3.0_dp;
      C(5,1)=8.0_dp/3.0_dp;
      C(5,2)=-64.0_dp/3.0_dp;
      C(5,3)=-64.0_dp/3.0_dp;
      C(5,4)=-64.0_dp/3.0_dp;
      C(5,5)=-64.0_dp/3.0_dp;
      C(5,6)=8.0_dp/3.0_dp;
      C(6,1)=80.0_dp/3.0_dp;
      C(6,2)=8.0_dp/3.0_dp;
      C(6,3)=-64.0_dp/3.0_dp;
      C(6,4)=80.0_dp/3.0_dp;
      C(6,5)=8.0_dp/3.0_dp;
      C(6,6)=8.0_dp/3.0_dp;
      endif
      if(n.eq.23) then
      C(1,1)=192.0_dp;
      C(1,2)=-24.0_dp;
      C(1,3)=0.0_dp;
      C(1,4)=0.0_dp;
      C(1,5)=24.0_dp;
      C(1,6)=48.0_dp;
      C(2,1)=-24.0_dp;
      C(2,2)=-24.0_dp;
      C(2,3)=0.0_dp;
      C(2,4)=-48.0_dp;
      C(2,5)=-24.0_dp;
      C(2,6)=0.0_dp;
      C(3,1)=0.0_dp;
      C(3,2)=0.0_dp;
      C(3,3)=192.0_dp;
      C(3,4)=-24.0_dp;
      C(3,5)=48.0_dp;
      C(3,6)=24.0_dp;
      C(4,1)=0.0_dp;
      C(4,2)=-48.0_dp;
      C(4,3)=-24.0_dp;
      C(4,4)=-24.0_dp;
      C(4,5)=0.0_dp;
      C(4,6)=-24.0_dp;
      C(5,1)=24.0_dp;
      C(5,2)=-24.0_dp;
      C(5,3)=48.0_dp;
      C(5,4)=0.0_dp;
      C(5,5)=192.0_dp;
      C(5,6)=0.0_dp;
      C(6,1)=48.0_dp;
      C(6,2)=0.0_dp;
      C(6,3)=24.0_dp;
      C(6,4)=-24.0_dp;
      C(6,5)=0.0_dp;
      C(6,6)=192.0_dp;
      endif
      if(n.eq.24) then
      C(1,1)=-24.0_dp;
      C(1,2)=-24.0_dp;
      C(1,3)=-24.0_dp;
      C(1,4)=0.0_dp;
      C(1,5)=0.0_dp;
      C(1,6)=-48.0_dp;
      C(2,1)=-24.0_dp;
      C(2,2)=192.0_dp;
      C(2,3)=24.0_dp;
      C(2,4)=48.0_dp;
      C(2,5)=0.0_dp;
      C(2,6)=0.0_dp;
      C(3,1)=-24.0_dp;
      C(3,2)=24.0_dp;
      C(3,3)=192.0_dp;
      C(3,4)=0.0_dp;
      C(3,5)=48.0_dp;
      C(3,6)=0.0_dp;
      C(4,1)=0.0_dp;
      C(4,2)=48.0_dp;
      C(4,3)=0.0_dp;
      C(4,4)=192.0_dp;
      C(4,5)=24.0_dp;
      C(4,6)=-24.0_dp;
      C(5,1)=0.0_dp;
      C(5,2)=0.0_dp;
      C(5,3)=48.0_dp;
      C(5,4)=24.0_dp;
      C(5,5)=192.0_dp;
      C(5,6)=-24.0_dp;
      C(6,1)=-48.0_dp;
      C(6,2)=0.0_dp;
      C(6,3)=0.0_dp;
      C(6,4)=-24.0_dp;
      C(6,5)=-24.0_dp;
      C(6,6)=-24.0_dp;
      endif
      if(n.eq.25) then
      C(1,1)=-64.0_dp/3.0_dp;
      C(1,2)=-64.0_dp/3.0_dp;
      C(1,3)=8.0_dp/3.0_dp;
      C(1,4)=-64.0_dp/3.0_dp;
      C(1,5)=8.0_dp/3.0_dp;
      C(1,6)=-64.0_dp/3.0_dp;
      C(2,1)=-64.0_dp/3.0_dp;
      C(2,2)=512.0_dp/3.0_dp;
      C(2,3)=8.0_dp/3.0_dp;
      C(2,4)=80.0_dp/3.0_dp;
      C(2,5)=-64.0_dp/3.0_dp;
      C(2,6)=8.0_dp/3.0_dp;
      C(3,1)=8.0_dp/3.0_dp;
      C(3,2)=8.0_dp/3.0_dp;
      C(3,3)=8.0_dp/3.0_dp;
      C(3,4)=80.0_dp/3.0_dp;
      C(3,5)=80.0_dp/3.0_dp;
      C(3,6)=-64.0_dp/3.0_dp;
      C(4,1)=-64.0_dp/3.0_dp;
      C(4,2)=80.0_dp/3.0_dp;
      C(4,3)=80.0_dp/3.0_dp;
      C(4,4)=8.0_dp/3.0_dp;
      C(4,5)=8.0_dp/3.0_dp;
      C(4,6)=8.0_dp/3.0_dp;
      C(5,1)=8.0_dp/3.0_dp;
      C(5,2)=-64.0_dp/3.0_dp;
      C(5,3)=80.0_dp/3.0_dp;
      C(5,4)=8.0_dp/3.0_dp;
      C(5,5)=512.0_dp/3.0_dp;
      C(5,6)=-64.0_dp/3.0_dp;
      C(6,1)=-64.0_dp/3.0_dp;
      C(6,2)=8.0_dp/3.0_dp;
      C(6,3)=-64.0_dp/3.0_dp;
      C(6,4)=8.0_dp/3.0_dp;
      C(6,5)=-64.0_dp/3.0_dp;
      C(6,6)=-64.0_dp/3.0_dp;
      endif
      if(n.eq.26) then
      C(1,1)=-64.0_dp/3.0_dp;
      C(1,2)=8.0_dp/3.0_dp;
      C(1,3)=-64.0_dp/3.0_dp;
      C(1,4)=8.0_dp/3.0_dp;
      C(1,5)=-64.0_dp/3.0_dp;
      C(1,6)=-64.0_dp/3.0_dp;
      C(2,1)=8.0_dp/3.0_dp;
      C(2,2)=8.0_dp/3.0_dp;
      C(2,3)=8.0_dp/3.0_dp;
      C(2,4)=80.0_dp/3.0_dp;
      C(2,5)=80.0_dp/3.0_dp;
      C(2,6)=-64.0_dp/3.0_dp;
      C(3,1)=-64.0_dp/3.0_dp;
      C(3,2)=8.0_dp/3.0_dp;
      C(3,3)=512.0_dp/3.0_dp;
      C(3,4)=-64.0_dp/3.0_dp;
      C(3,5)=80.0_dp/3.0_dp;
      C(3,6)=8.0_dp/3.0_dp;
      C(4,1)=8.0_dp/3.0_dp;
      C(4,2)=80.0_dp/3.0_dp;
      C(4,3)=-64.0_dp/3.0_dp;
      C(4,4)=512.0_dp/3.0_dp;
      C(4,5)=8.0_dp/3.0_dp;
      C(4,6)=-64.0_dp/3.0_dp;
      C(5,1)=-64.0_dp/3.0_dp;
      C(5,2)=80.0_dp/3.0_dp;
      C(5,3)=80.0_dp/3.0_dp;
      C(5,4)=8.0_dp/3.0_dp;
      C(5,5)=8.0_dp/3.0_dp;
      C(5,6)=8.0_dp/3.0_dp;
      C(6,1)=-64.0_dp/3.0_dp;
      C(6,2)=-64.0_dp/3.0_dp;
      C(6,3)=8.0_dp/3.0_dp;
      C(6,4)=-64.0_dp/3.0_dp;
      C(6,5)=8.0_dp/3.0_dp;
      C(6,6)=-64.0_dp/3.0_dp;
      endif
      if(n.eq.27) then
      C(1,1)=192.0_dp;
      C(1,2)=-24.0_dp;
      C(1,3)=0.0_dp;
      C(1,4)=0.0_dp;
      C(1,5)=24.0_dp;
      C(1,6)=48.0_dp;
      C(2,1)=-24.0_dp;
      C(2,2)=-24.0_dp;
      C(2,3)=0.0_dp;
      C(2,4)=-48.0_dp;
      C(2,5)=-24.0_dp;
      C(2,6)=0.0_dp;
      C(3,1)=0.0_dp;
      C(3,2)=0.0_dp;
      C(3,3)=192.0_dp;
      C(3,4)=-24.0_dp;
      C(3,5)=48.0_dp;
      C(3,6)=24.0_dp;
      C(4,1)=0.0_dp;
      C(4,2)=-48.0_dp;
      C(4,3)=-24.0_dp;
      C(4,4)=-24.0_dp;
      C(4,5)=0.0_dp;
      C(4,6)=-24.0_dp;
      C(5,1)=24.0_dp;
      C(5,2)=-24.0_dp;
      C(5,3)=48.0_dp;
      C(5,4)=0.0_dp;
      C(5,5)=192.0_dp;
      C(5,6)=0.0_dp;
      C(6,1)=48.0_dp;
      C(6,2)=0.0_dp;
      C(6,3)=24.0_dp;
      C(6,4)=-24.0_dp;
      C(6,5)=0.0_dp;
      C(6,6)=192.0_dp;
      endif
      if(n.eq.28) then
      C(1,1)=192.0_dp;
      C(1,2)=0.0_dp;
      C(1,3)=-24.0_dp;
      C(1,4)=24.0_dp;
      C(1,5)=0.0_dp;
      C(1,6)=48.0_dp;
      C(2,1)=0.0_dp;
      C(2,2)=192.0_dp;
      C(2,3)=0.0_dp;
      C(2,4)=48.0_dp;
      C(2,5)=-24.0_dp;
      C(2,6)=24.0_dp;
      C(3,1)=-24.0_dp;
      C(3,2)=0.0_dp;
      C(3,3)=-24.0_dp;
      C(3,4)=-24.0_dp;
      C(3,5)=-48.0_dp;
      C(3,6)=0.0_dp;
      C(4,1)=24.0_dp;
      C(4,2)=48.0_dp;
      C(4,3)=-24.0_dp;
      C(4,4)=192.0_dp;
      C(4,5)=0.0_dp;
      C(4,6)=0.0_dp;
      C(5,1)=0.0_dp;
      C(5,2)=-24.0_dp;
      C(5,3)=-48.0_dp;
      C(5,4)=0.0_dp;
      C(5,5)=-24.0_dp;
      C(5,6)=-24.0_dp;
      C(6,1)=48.0_dp;
      C(6,2)=24.0_dp;
      C(6,3)=0.0_dp;
      C(6,4)=0.0_dp;
      C(6,5)=-24.0_dp;
      C(6,6)=192.0_dp;
      endif
      if(n.eq.29) then
      C(1,1)=-8.0_dp/27.0_dp;
      C(1,2)=-80.0_dp/27.0_dp;
      C(1,3)=-80.0_dp/27.0_dp;
      C(1,4)=496.0_dp/27.0_dp;
      C(1,5)=496.0_dp/27.0_dp;
      C(1,6)=-224.0_dp/27.0_dp;
      C(2,1)=-80.0_dp/27.0_dp;
      C(2,2)=-8.0_dp/27.0_dp;
      C(2,3)=496.0_dp/27.0_dp;
      C(2,4)=-224.0_dp/27.0_dp;
      C(2,5)=-80.0_dp/27.0_dp;
      C(2,6)=496.0_dp/27.0_dp;
      C(3,1)=-80.0_dp/27.0_dp;
      C(3,2)=496.0_dp/27.0_dp;
      C(3,3)=-8.0_dp/27.0_dp;
      C(3,4)=-80.0_dp/27.0_dp;
      C(3,5)=-224.0_dp/27.0_dp;
      C(3,6)=496.0_dp/27.0_dp;
      C(4,1)=496.0_dp/27.0_dp;
      C(4,2)=-224.0_dp/27.0_dp;
      C(4,3)=-80.0_dp/27.0_dp;
      C(4,4)=-8.0_dp/27.0_dp;
      C(4,5)=496.0_dp/27.0_dp;
      C(4,6)=-80.0_dp/27.0_dp;
      C(5,1)=496.0_dp/27.0_dp;
      C(5,2)=-80.0_dp/27.0_dp;
      C(5,3)=-224.0_dp/27.0_dp;
      C(5,4)=496.0_dp/27.0_dp;
      C(5,5)=-8.0_dp/27.0_dp;
      C(5,6)=-80.0_dp/27.0_dp;
      C(6,1)=-224.0_dp/27.0_dp;
      C(6,2)=496.0_dp/27.0_dp;
      C(6,3)=496.0_dp/27.0_dp;
      C(6,4)=-80.0_dp/27.0_dp;
      C(6,5)=-80.0_dp/27.0_dp;
      C(6,6)=-8.0_dp/27.0_dp;
      endif
      if(n.eq.30) then
      C(1,1)=8.0_dp/3.0_dp;
      C(1,2)=80.0_dp/3.0_dp;
      C(1,3)=8.0_dp/3.0_dp;
      C(1,4)=8.0_dp/3.0_dp;
      C(1,5)=-64.0_dp/3.0_dp;
      C(1,6)=80.0_dp/3.0_dp;
      C(2,1)=80.0_dp/3.0_dp;
      C(2,2)=8.0_dp/3.0_dp;
      C(2,3)=-64.0_dp/3.0_dp;
      C(2,4)=80.0_dp/3.0_dp;
      C(2,5)=8.0_dp/3.0_dp;
      C(2,6)=8.0_dp/3.0_dp;
      C(3,1)=8.0_dp/3.0_dp;
      C(3,2)=-64.0_dp/3.0_dp;
      C(3,3)=-64.0_dp/3.0_dp;
      C(3,4)=-64.0_dp/3.0_dp;
      C(3,5)=-64.0_dp/3.0_dp;
      C(3,6)=8.0_dp/3.0_dp;
      C(4,1)=8.0_dp/3.0_dp;
      C(4,2)=80.0_dp/3.0_dp;
      C(4,3)=-64.0_dp/3.0_dp;
      C(4,4)=512.0_dp/3.0_dp;
      C(4,5)=8.0_dp/3.0_dp;
      C(4,6)=-64.0_dp/3.0_dp;
      C(5,1)=-64.0_dp/3.0_dp;
      C(5,2)=8.0_dp/3.0_dp;
      C(5,3)=-64.0_dp/3.0_dp;
      C(5,4)=8.0_dp/3.0_dp;
      C(5,5)=-64.0_dp/3.0_dp;
      C(5,6)=-64.0_dp/3.0_dp;
      C(6,1)=80.0_dp/3.0_dp;
      C(6,2)=8.0_dp/3.0_dp;
      C(6,3)=8.0_dp/3.0_dp;
      C(6,4)=-64.0_dp/3.0_dp;
      C(6,5)=-64.0_dp/3.0_dp;
      C(6,6)=512.0_dp/3.0_dp;
      endif
      if(n.eq.31) then
      C(1,1)=-64.0_dp/3.0_dp;
      C(1,2)=-64.0_dp/3.0_dp;
      C(1,3)=8.0_dp/3.0_dp;
      C(1,4)=-64.0_dp/3.0_dp;
      C(1,5)=8.0_dp/3.0_dp;
      C(1,6)=-64.0_dp/3.0_dp;
      C(2,1)=-64.0_dp/3.0_dp;
      C(2,2)=512.0_dp/3.0_dp;
      C(2,3)=8.0_dp/3.0_dp;
      C(2,4)=80.0_dp/3.0_dp;
      C(2,5)=-64.0_dp/3.0_dp;
      C(2,6)=8.0_dp/3.0_dp;
      C(3,1)=8.0_dp/3.0_dp;
      C(3,2)=8.0_dp/3.0_dp;
      C(3,3)=8.0_dp/3.0_dp;
      C(3,4)=80.0_dp/3.0_dp;
      C(3,5)=80.0_dp/3.0_dp;
      C(3,6)=-64.0_dp/3.0_dp;
      C(4,1)=-64.0_dp/3.0_dp;
      C(4,2)=80.0_dp/3.0_dp;
      C(4,3)=80.0_dp/3.0_dp;
      C(4,4)=8.0_dp/3.0_dp;
      C(4,5)=8.0_dp/3.0_dp;
      C(4,6)=8.0_dp/3.0_dp;
      C(5,1)=8.0_dp/3.0_dp;
      C(5,2)=-64.0_dp/3.0_dp;
      C(5,3)=80.0_dp/3.0_dp;
      C(5,4)=8.0_dp/3.0_dp;
      C(5,5)=512.0_dp/3.0_dp;
      C(5,6)=-64.0_dp/3.0_dp;
      C(6,1)=-64.0_dp/3.0_dp;
      C(6,2)=8.0_dp/3.0_dp;
      C(6,3)=-64.0_dp/3.0_dp;
      C(6,4)=8.0_dp/3.0_dp;
      C(6,5)=-64.0_dp/3.0_dp;
      C(6,6)=-64.0_dp/3.0_dp;
      endif
      if(n.eq.32) then
      C(1,1)=512.0_dp/3.0_dp;
      C(1,2)=-64.0_dp/3.0_dp;
      C(1,3)=-64.0_dp/3.0_dp;
      C(1,4)=8.0_dp/3.0_dp;
      C(1,5)=8.0_dp/3.0_dp;
      C(1,6)=80.0_dp/3.0_dp;
      C(2,1)=-64.0_dp/3.0_dp;
      C(2,2)=-64.0_dp/3.0_dp;
      C(2,3)=8.0_dp/3.0_dp;
      C(2,4)=-64.0_dp/3.0_dp;
      C(2,5)=8.0_dp/3.0_dp;
      C(2,6)=-64.0_dp/3.0_dp;
      C(3,1)=-64.0_dp/3.0_dp;
      C(3,2)=8.0_dp/3.0_dp;
      C(3,3)=512.0_dp/3.0_dp;
      C(3,4)=-64.0_dp/3.0_dp;
      C(3,5)=80.0_dp/3.0_dp;
      C(3,6)=8.0_dp/3.0_dp;
      C(4,1)=8.0_dp/3.0_dp;
      C(4,2)=-64.0_dp/3.0_dp;
      C(4,3)=-64.0_dp/3.0_dp;
      C(4,4)=-64.0_dp/3.0_dp;
      C(4,5)=-64.0_dp/3.0_dp;
      C(4,6)=8.0_dp/3.0_dp;
      C(5,1)=8.0_dp/3.0_dp;
      C(5,2)=8.0_dp/3.0_dp;
      C(5,3)=80.0_dp/3.0_dp;
      C(5,4)=-64.0_dp/3.0_dp;
      C(5,5)=8.0_dp/3.0_dp;
      C(5,6)=80.0_dp/3.0_dp;
      C(6,1)=80.0_dp/3.0_dp;
      C(6,2)=-64.0_dp/3.0_dp;
      C(6,3)=8.0_dp/3.0_dp;
      C(6,4)=8.0_dp/3.0_dp;
      C(6,5)=80.0_dp/3.0_dp;
      C(6,6)=8.0_dp/3.0_dp;
      endif
      if(n.eq.33) then
      C(1,1)=-8.0_dp/27.0_dp;
      C(1,2)=-80.0_dp/27.0_dp;
      C(1,3)=-80.0_dp/27.0_dp;
      C(1,4)=496.0_dp/27.0_dp;
      C(1,5)=496.0_dp/27.0_dp;
      C(1,6)=-224.0_dp/27.0_dp;
      C(2,1)=-80.0_dp/27.0_dp;
      C(2,2)=-8.0_dp/27.0_dp;
      C(2,3)=496.0_dp/27.0_dp;
      C(2,4)=-224.0_dp/27.0_dp;
      C(2,5)=-80.0_dp/27.0_dp;
      C(2,6)=496.0_dp/27.0_dp;
      C(3,1)=-80.0_dp/27.0_dp;
      C(3,2)=496.0_dp/27.0_dp;
      C(3,3)=-8.0_dp/27.0_dp;
      C(3,4)=-80.0_dp/27.0_dp;
      C(3,5)=-224.0_dp/27.0_dp;
      C(3,6)=496.0_dp/27.0_dp;
      C(4,1)=496.0_dp/27.0_dp;
      C(4,2)=-224.0_dp/27.0_dp;
      C(4,3)=-80.0_dp/27.0_dp;
      C(4,4)=-8.0_dp/27.0_dp;
      C(4,5)=496.0_dp/27.0_dp;
      C(4,6)=-80.0_dp/27.0_dp;
      C(5,1)=496.0_dp/27.0_dp;
      C(5,2)=-80.0_dp/27.0_dp;
      C(5,3)=-224.0_dp/27.0_dp;
      C(5,4)=496.0_dp/27.0_dp;
      C(5,5)=-8.0_dp/27.0_dp;
      C(5,6)=-80.0_dp/27.0_dp;
      C(6,1)=-224.0_dp/27.0_dp;
      C(6,2)=496.0_dp/27.0_dp;
      C(6,3)=496.0_dp/27.0_dp;
      C(6,4)=-80.0_dp/27.0_dp;
      C(6,5)=-80.0_dp/27.0_dp;
      C(6,6)=-8.0_dp/27.0_dp;
      endif
      if(n.eq.34) then
      C(1,1)=512.0_dp/3.0_dp;
      C(1,2)=-64.0_dp/3.0_dp;
      C(1,3)=-64.0_dp/3.0_dp;
      C(1,4)=8.0_dp/3.0_dp;
      C(1,5)=8.0_dp/3.0_dp;
      C(1,6)=80.0_dp/3.0_dp;
      C(2,1)=-64.0_dp/3.0_dp;
      C(2,2)=512.0_dp/3.0_dp;
      C(2,3)=8.0_dp/3.0_dp;
      C(2,4)=80.0_dp/3.0_dp;
      C(2,5)=-64.0_dp/3.0_dp;
      C(2,6)=8.0_dp/3.0_dp;
      C(3,1)=-64.0_dp/3.0_dp;
      C(3,2)=8.0_dp/3.0_dp;
      C(3,3)=-64.0_dp/3.0_dp;
      C(3,4)=8.0_dp/3.0_dp;
      C(3,5)=-64.0_dp/3.0_dp;
      C(3,6)=-64.0_dp/3.0_dp;
      C(4,1)=8.0_dp/3.0_dp;
      C(4,2)=80.0_dp/3.0_dp;
      C(4,3)=8.0_dp/3.0_dp;
      C(4,4)=8.0_dp/3.0_dp;
      C(4,5)=-64.0_dp/3.0_dp;
      C(4,6)=80.0_dp/3.0_dp;
      C(5,1)=8.0_dp/3.0_dp;
      C(5,2)=-64.0_dp/3.0_dp;
      C(5,3)=-64.0_dp/3.0_dp;
      C(5,4)=-64.0_dp/3.0_dp;
      C(5,5)=-64.0_dp/3.0_dp;
      C(5,6)=8.0_dp/3.0_dp;
      C(6,1)=80.0_dp/3.0_dp;
      C(6,2)=8.0_dp/3.0_dp;
      C(6,3)=-64.0_dp/3.0_dp;
      C(6,4)=80.0_dp/3.0_dp;
      C(6,5)=8.0_dp/3.0_dp;
      C(6,6)=8.0_dp/3.0_dp;
      endif
      if(n.eq.35) then
      C(1,1)=-64.0_dp/3.0_dp;
      C(1,2)=8.0_dp/3.0_dp;
      C(1,3)=-64.0_dp/3.0_dp;
      C(1,4)=8.0_dp/3.0_dp;
      C(1,5)=-64.0_dp/3.0_dp;
      C(1,6)=-64.0_dp/3.0_dp;
      C(2,1)=8.0_dp/3.0_dp;
      C(2,2)=8.0_dp/3.0_dp;
      C(2,3)=8.0_dp/3.0_dp;
      C(2,4)=80.0_dp/3.0_dp;
      C(2,5)=80.0_dp/3.0_dp;
      C(2,6)=-64.0_dp/3.0_dp;
      C(3,1)=-64.0_dp/3.0_dp;
      C(3,2)=8.0_dp/3.0_dp;
      C(3,3)=512.0_dp/3.0_dp;
      C(3,4)=-64.0_dp/3.0_dp;
      C(3,5)=80.0_dp/3.0_dp;
      C(3,6)=8.0_dp/3.0_dp;
      C(4,1)=8.0_dp/3.0_dp;
      C(4,2)=80.0_dp/3.0_dp;
      C(4,3)=-64.0_dp/3.0_dp;
      C(4,4)=512.0_dp/3.0_dp;
      C(4,5)=8.0_dp/3.0_dp;
      C(4,6)=-64.0_dp/3.0_dp;
      C(5,1)=-64.0_dp/3.0_dp;
      C(5,2)=80.0_dp/3.0_dp;
      C(5,3)=80.0_dp/3.0_dp;
      C(5,4)=8.0_dp/3.0_dp;
      C(5,5)=8.0_dp/3.0_dp;
      C(5,6)=8.0_dp/3.0_dp;
      C(6,1)=-64.0_dp/3.0_dp;
      C(6,2)=-64.0_dp/3.0_dp;
      C(6,3)=8.0_dp/3.0_dp;
      C(6,4)=-64.0_dp/3.0_dp;
      C(6,5)=8.0_dp/3.0_dp;
      C(6,6)=-64.0_dp/3.0_dp;
      endif
      if(n.eq.36) then
      C(1,1)=8.0_dp/3.0_dp;
      C(1,2)=8.0_dp/3.0_dp;
      C(1,3)=80.0_dp/3.0_dp;
      C(1,4)=-64.0_dp/3.0_dp;
      C(1,5)=8.0_dp/3.0_dp;
      C(1,6)=80.0_dp/3.0_dp;
      C(2,1)=8.0_dp/3.0_dp;
      C(2,2)=-64.0_dp/3.0_dp;
      C(2,3)=-64.0_dp/3.0_dp;
      C(2,4)=-64.0_dp/3.0_dp;
      C(2,5)=-64.0_dp/3.0_dp;
      C(2,6)=8.0_dp/3.0_dp;
      C(3,1)=80.0_dp/3.0_dp;
      C(3,2)=-64.0_dp/3.0_dp;
      C(3,3)=8.0_dp/3.0_dp;
      C(3,4)=8.0_dp/3.0_dp;
      C(3,5)=80.0_dp/3.0_dp;
      C(3,6)=8.0_dp/3.0_dp;
      C(4,1)=-64.0_dp/3.0_dp;
      C(4,2)=-64.0_dp/3.0_dp;
      C(4,3)=8.0_dp/3.0_dp;
      C(4,4)=-64.0_dp/3.0_dp;
      C(4,5)=8.0_dp/3.0_dp;
      C(4,6)=-64.0_dp/3.0_dp;
      C(5,1)=8.0_dp/3.0_dp;
      C(5,2)=-64.0_dp/3.0_dp;
      C(5,3)=80.0_dp/3.0_dp;
      C(5,4)=8.0_dp/3.0_dp;
      C(5,5)=512.0_dp/3.0_dp;
      C(5,6)=-64.0_dp/3.0_dp;
      C(6,1)=80.0_dp/3.0_dp;
      C(6,2)=8.0_dp/3.0_dp;
      C(6,3)=8.0_dp/3.0_dp;
      C(6,4)=-64.0_dp/3.0_dp;
      C(6,5)=-64.0_dp/3.0_dp;
      C(6,6)=512.0_dp/3.0_dp;
      endif
      if(n.eq.37) then
      C(1,1)=512.0_dp/3.0_dp;
      C(1,2)=-64.0_dp/3.0_dp;
      C(1,3)=-64.0_dp/3.0_dp;
      C(1,4)=8.0_dp/3.0_dp;
      C(1,5)=8.0_dp/3.0_dp;
      C(1,6)=80.0_dp/3.0_dp;
      C(2,1)=-64.0_dp/3.0_dp;
      C(2,2)=-64.0_dp/3.0_dp;
      C(2,3)=8.0_dp/3.0_dp;
      C(2,4)=-64.0_dp/3.0_dp;
      C(2,5)=8.0_dp/3.0_dp;
      C(2,6)=-64.0_dp/3.0_dp;
      C(3,1)=-64.0_dp/3.0_dp;
      C(3,2)=8.0_dp/3.0_dp;
      C(3,3)=512.0_dp/3.0_dp;
      C(3,4)=-64.0_dp/3.0_dp;
      C(3,5)=80.0_dp/3.0_dp;
      C(3,6)=8.0_dp/3.0_dp;
      C(4,1)=8.0_dp/3.0_dp;
      C(4,2)=-64.0_dp/3.0_dp;
      C(4,3)=-64.0_dp/3.0_dp;
      C(4,4)=-64.0_dp/3.0_dp;
      C(4,5)=-64.0_dp/3.0_dp;
      C(4,6)=8.0_dp/3.0_dp;
      C(5,1)=8.0_dp/3.0_dp;
      C(5,2)=8.0_dp/3.0_dp;
      C(5,3)=80.0_dp/3.0_dp;
      C(5,4)=-64.0_dp/3.0_dp;
      C(5,5)=8.0_dp/3.0_dp;
      C(5,6)=80.0_dp/3.0_dp;
      C(6,1)=80.0_dp/3.0_dp;
      C(6,2)=-64.0_dp/3.0_dp;
      C(6,3)=8.0_dp/3.0_dp;
      C(6,4)=8.0_dp/3.0_dp;
      C(6,5)=80.0_dp/3.0_dp;
      C(6,6)=8.0_dp/3.0_dp;
      endif
      if(n.eq.38) then
      C(1,1)=8.0_dp/3.0_dp;
      C(1,2)=8.0_dp/3.0_dp;
      C(1,3)=80.0_dp/3.0_dp;
      C(1,4)=-64.0_dp/3.0_dp;
      C(1,5)=8.0_dp/3.0_dp;
      C(1,6)=80.0_dp/3.0_dp;
      C(2,1)=8.0_dp/3.0_dp;
      C(2,2)=-64.0_dp/3.0_dp;
      C(2,3)=-64.0_dp/3.0_dp;
      C(2,4)=-64.0_dp/3.0_dp;
      C(2,5)=-64.0_dp/3.0_dp;
      C(2,6)=8.0_dp/3.0_dp;
      C(3,1)=80.0_dp/3.0_dp;
      C(3,2)=-64.0_dp/3.0_dp;
      C(3,3)=8.0_dp/3.0_dp;
      C(3,4)=8.0_dp/3.0_dp;
      C(3,5)=80.0_dp/3.0_dp;
      C(3,6)=8.0_dp/3.0_dp;
      C(4,1)=-64.0_dp/3.0_dp;
      C(4,2)=-64.0_dp/3.0_dp;
      C(4,3)=8.0_dp/3.0_dp;
      C(4,4)=-64.0_dp/3.0_dp;
      C(4,5)=8.0_dp/3.0_dp;
      C(4,6)=-64.0_dp/3.0_dp;
      C(5,1)=8.0_dp/3.0_dp;
      C(5,2)=-64.0_dp/3.0_dp;
      C(5,3)=80.0_dp/3.0_dp;
      C(5,4)=8.0_dp/3.0_dp;
      C(5,5)=512.0_dp/3.0_dp;
      C(5,6)=-64.0_dp/3.0_dp;
      C(6,1)=80.0_dp/3.0_dp;
      C(6,2)=8.0_dp/3.0_dp;
      C(6,3)=8.0_dp/3.0_dp;
      C(6,4)=-64.0_dp/3.0_dp;
      C(6,5)=-64.0_dp/3.0_dp;
      C(6,6)=512.0_dp/3.0_dp;
      endif
      if(n.eq.39) then
      C(1,1)=-24.0_dp;
      C(1,2)=-24.0_dp;
      C(1,3)=-24.0_dp;
      C(1,4)=0.0_dp;
      C(1,5)=0.0_dp;
      C(1,6)=-48.0_dp;
      C(2,1)=-24.0_dp;
      C(2,2)=192.0_dp;
      C(2,3)=24.0_dp;
      C(2,4)=48.0_dp;
      C(2,5)=0.0_dp;
      C(2,6)=0.0_dp;
      C(3,1)=-24.0_dp;
      C(3,2)=24.0_dp;
      C(3,3)=192.0_dp;
      C(3,4)=0.0_dp;
      C(3,5)=48.0_dp;
      C(3,6)=0.0_dp;
      C(4,1)=0.0_dp;
      C(4,2)=48.0_dp;
      C(4,3)=0.0_dp;
      C(4,4)=192.0_dp;
      C(4,5)=24.0_dp;
      C(4,6)=-24.0_dp;
      C(5,1)=0.0_dp;
      C(5,2)=0.0_dp;
      C(5,3)=48.0_dp;
      C(5,4)=24.0_dp;
      C(5,5)=192.0_dp;
      C(5,6)=-24.0_dp;
      C(6,1)=-48.0_dp;
      C(6,2)=0.0_dp;
      C(6,3)=0.0_dp;
      C(6,4)=-24.0_dp;
      C(6,5)=-24.0_dp;
      C(6,6)=-24.0_dp;
      endif
      if(n.eq.40) then
      C(1,1)=192.0_dp;
      C(1,2)=0.0_dp;
      C(1,3)=-24.0_dp;
      C(1,4)=24.0_dp;
      C(1,5)=0.0_dp;
      C(1,6)=48.0_dp;
      C(2,1)=0.0_dp;
      C(2,2)=192.0_dp;
      C(2,3)=0.0_dp;
      C(2,4)=48.0_dp;
      C(2,5)=-24.0_dp;
      C(2,6)=24.0_dp;
      C(3,1)=-24.0_dp;
      C(3,2)=0.0_dp;
      C(3,3)=-24.0_dp;
      C(3,4)=-24.0_dp;
      C(3,5)=-48.0_dp;
      C(3,6)=0.0_dp;
      C(4,1)=24.0_dp;
      C(4,2)=48.0_dp;
      C(4,3)=-24.0_dp;
      C(4,4)=192.0_dp;
      C(4,5)=0.0_dp;
      C(4,6)=0.0_dp;
      C(5,1)=0.0_dp;
      C(5,2)=-24.0_dp;
      C(5,3)=-48.0_dp;
      C(5,4)=0.0_dp;
      C(5,5)=-24.0_dp;
      C(5,6)=-24.0_dp;
      C(6,1)=48.0_dp;
      C(6,2)=24.0_dp;
      C(6,3)=0.0_dp;
      C(6,4)=0.0_dp;
      C(6,5)=-24.0_dp;
      C(6,6)=192.0_dp;
      endif


      end subroutine cc_gg_ttgggg



      subroutine cc_qq_ttqqgg(n,C)
      integer, intent(in)  :: n
      real(dp), intent(out)  :: C(4,4)


      if(n.eq.1) then
      C(1,1)=0.0_dp
      C(1,2)=8.0_dp/9.0_dp
      C(1,3)=-8.0_dp
      C(1,4)=-64.0_dp/9.0_dp
      C(2,1)=8.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=-64.0_dp/9.0_dp
      C(2,4)=-8.0_dp/9.0_dp
      C(3,1)=-8.0_dp
      C(3,2)=-64.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=8.0_dp/9.0_dp
      C(4,1)=-64.0_dp/9.0_dp
      C(4,2)=-8.0_dp/9.0_dp
      C(4,3)=8.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.2) then
      C(1,1)=64.0_dp
      C(1,2)=64.0_dp/9.0_dp
      C(1,3)=0.0_dp
      C(1,4)=64.0_dp/9.0_dp
      C(2,1)=64.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=-8.0_dp/9.0_dp
      C(2,4)=8.0_dp/9.0_dp
      C(3,1)=0.0_dp
      C(3,2)=-8.0_dp/9.0_dp
      C(3,3)=-8.0_dp
      C(3,4)=-8.0_dp/9.0_dp
      C(4,1)=64.0_dp/9.0_dp
      C(4,2)=8.0_dp/9.0_dp
      C(4,3)=-8.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.3) then
      C(1,1)=0.0_dp
      C(1,2)=64.0_dp/9.0_dp
      C(1,3)=8.0_dp
      C(1,4)=-8.0_dp/9.0_dp
      C(2,1)=64.0_dp/9.0_dp
      C(2,2)=64.0_dp/9.0_dp
      C(2,3)=64.0_dp/9.0_dp
      C(2,4)=0.0_dp
      C(3,1)=8.0_dp
      C(3,2)=64.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=-8.0_dp/9.0_dp
      C(4,1)=-8.0_dp/9.0_dp
      C(4,2)=0.0_dp
      C(4,3)=-8.0_dp/9.0_dp
      C(4,4)=-8.0_dp/9.0_dp
      endif
      if(n.eq.4) then
      C(1,1)=0.0_dp
      C(1,2)=-8.0_dp
      C(1,3)=0.0_dp
      C(1,4)=8.0_dp
      C(2,1)=-8.0_dp
      C(2,2)=0.0_dp
      C(2,3)=8.0_dp
      C(2,4)=0.0_dp
      C(3,1)=0.0_dp
      C(3,2)=8.0_dp
      C(3,3)=72.0_dp
      C(3,4)=8.0_dp
      C(4,1)=8.0_dp
      C(4,2)=0.0_dp
      C(4,3)=8.0_dp
      C(4,4)=8.0_dp
      endif
      if(n.eq.5) then
      C(1,1)=-8.0_dp
      C(1,2)=-8.0_dp/9.0_dp
      C(1,3)=0.0_dp
      C(1,4)=-8.0_dp/9.0_dp
      C(2,1)=-8.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=64.0_dp/9.0_dp
      C(2,4)=8.0_dp/9.0_dp
      C(3,1)=0.0_dp
      C(3,2)=64.0_dp/9.0_dp
      C(3,3)=64.0_dp
      C(3,4)=64.0_dp/9.0_dp
      C(4,1)=-8.0_dp/9.0_dp
      C(4,2)=8.0_dp/9.0_dp
      C(4,3)=64.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.6) then
      C(1,1)=0.0_dp
      C(1,2)=-64.0_dp/9.0_dp
      C(1,3)=-8.0_dp
      C(1,4)=8.0_dp/9.0_dp
      C(2,1)=-64.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=8.0_dp/9.0_dp
      C(2,4)=-8.0_dp/9.0_dp
      C(3,1)=-8.0_dp
      C(3,2)=8.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=-64.0_dp/9.0_dp
      C(4,1)=8.0_dp/9.0_dp
      C(4,2)=-8.0_dp/9.0_dp
      C(4,3)=-64.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.7) then
      C(1,1)=0.0_dp
      C(1,2)=64.0_dp/9.0_dp
      C(1,3)=8.0_dp
      C(1,4)=-8.0_dp/9.0_dp
      C(2,1)=64.0_dp/9.0_dp
      C(2,2)=64.0_dp/9.0_dp
      C(2,3)=64.0_dp/9.0_dp
      C(2,4)=0.0_dp
      C(3,1)=8.0_dp
      C(3,2)=64.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=-8.0_dp/9.0_dp
      C(4,1)=-8.0_dp/9.0_dp
      C(4,2)=0.0_dp
      C(4,3)=-8.0_dp/9.0_dp
      C(4,4)=-8.0_dp/9.0_dp
      endif
      if(n.eq.8) then
      C(1,1)=72.0_dp
      C(1,2)=8.0_dp
      C(1,3)=0.0_dp
      C(1,4)=8.0_dp
      C(2,1)=8.0_dp
      C(2,2)=0.0_dp
      C(2,3)=-8.0_dp
      C(2,4)=0.0_dp
      C(3,1)=0.0_dp
      C(3,2)=-8.0_dp
      C(3,3)=0.0_dp
      C(3,4)=8.0_dp
      C(4,1)=8.0_dp
      C(4,2)=0.0_dp
      C(4,3)=8.0_dp
      C(4,4)=8.0_dp
      endif
      if(n.eq.9) then
      C(1,1)=0.0_dp
      C(1,2)=-8.0_dp/9.0_dp
      C(1,3)=8.0_dp
      C(1,4)=64.0_dp/9.0_dp
      C(2,1)=-8.0_dp/9.0_dp
      C(2,2)=-8.0_dp/9.0_dp
      C(2,3)=-8.0_dp/9.0_dp
      C(2,4)=0.0_dp
      C(3,1)=8.0_dp
      C(3,2)=-8.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=64.0_dp/9.0_dp
      C(4,1)=64.0_dp/9.0_dp
      C(4,2)=0.0_dp
      C(4,3)=64.0_dp/9.0_dp
      C(4,4)=64.0_dp/9.0_dp
      endif
      if(n.eq.10) then
      C(1,1)=0.0_dp
      C(1,2)=8.0_dp/9.0_dp
      C(1,3)=-8.0_dp
      C(1,4)=-64.0_dp/9.0_dp
      C(2,1)=8.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=-64.0_dp/9.0_dp
      C(2,4)=-8.0_dp/9.0_dp
      C(3,1)=-8.0_dp
      C(3,2)=-64.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=8.0_dp/9.0_dp
      C(4,1)=-64.0_dp/9.0_dp
      C(4,2)=-8.0_dp/9.0_dp
      C(4,3)=8.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.11) then
      C(1,1)=-8.0_dp
      C(1,2)=-8.0_dp/9.0_dp
      C(1,3)=0.0_dp
      C(1,4)=-8.0_dp/9.0_dp
      C(2,1)=-8.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=64.0_dp/9.0_dp
      C(2,4)=8.0_dp/9.0_dp
      C(3,1)=0.0_dp
      C(3,2)=64.0_dp/9.0_dp
      C(3,3)=64.0_dp
      C(3,4)=64.0_dp/9.0_dp
      C(4,1)=-8.0_dp/9.0_dp
      C(4,2)=8.0_dp/9.0_dp
      C(4,3)=64.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.12) then
      C(1,1)=72.0_dp
      C(1,2)=8.0_dp
      C(1,3)=0.0_dp
      C(1,4)=8.0_dp
      C(2,1)=8.0_dp
      C(2,2)=8.0_dp
      C(2,3)=8.0_dp
      C(2,4)=0.0_dp
      C(3,1)=0.0_dp
      C(3,2)=8.0_dp
      C(3,3)=0.0_dp
      C(3,4)=-8.0_dp
      C(4,1)=8.0_dp
      C(4,2)=0.0_dp
      C(4,3)=-8.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.13) then
      C(1,1)=0.0_dp
      C(1,2)=-8.0_dp/9.0_dp
      C(1,3)=8.0_dp
      C(1,4)=64.0_dp/9.0_dp
      C(2,1)=-8.0_dp/9.0_dp
      C(2,2)=-8.0_dp/9.0_dp
      C(2,3)=-8.0_dp/9.0_dp
      C(2,4)=0.0_dp
      C(3,1)=8.0_dp
      C(3,2)=-8.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=64.0_dp/9.0_dp
      C(4,1)=64.0_dp/9.0_dp
      C(4,2)=0.0_dp
      C(4,3)=64.0_dp/9.0_dp
      C(4,4)=64.0_dp/9.0_dp
      endif
      if(n.eq.14) then
      C(1,1)=64.0_dp
      C(1,2)=64.0_dp/9.0_dp
      C(1,3)=0.0_dp
      C(1,4)=64.0_dp/9.0_dp
      C(2,1)=64.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=-8.0_dp/9.0_dp
      C(2,4)=8.0_dp/9.0_dp
      C(3,1)=0.0_dp
      C(3,2)=-8.0_dp/9.0_dp
      C(3,3)=-8.0_dp
      C(3,4)=-8.0_dp/9.0_dp
      C(4,1)=64.0_dp/9.0_dp
      C(4,2)=8.0_dp/9.0_dp
      C(4,3)=-8.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.15) then
      C(1,1)=0.0_dp
      C(1,2)=-64.0_dp/9.0_dp
      C(1,3)=-8.0_dp
      C(1,4)=8.0_dp/9.0_dp
      C(2,1)=-64.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=8.0_dp/9.0_dp
      C(2,4)=-8.0_dp/9.0_dp
      C(3,1)=-8.0_dp
      C(3,2)=8.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=-64.0_dp/9.0_dp
      C(4,1)=8.0_dp/9.0_dp
      C(4,2)=-8.0_dp/9.0_dp
      C(4,3)=-64.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.16) then
      C(1,1)=0.0_dp
      C(1,2)=8.0_dp
      C(1,3)=0.0_dp
      C(1,4)=-8.0_dp
      C(2,1)=8.0_dp
      C(2,2)=8.0_dp
      C(2,3)=8.0_dp
      C(2,4)=0.0_dp
      C(3,1)=0.0_dp
      C(3,2)=8.0_dp
      C(3,3)=72.0_dp
      C(3,4)=8.0_dp
      C(4,1)=-8.0_dp
      C(4,2)=0.0_dp
      C(4,3)=8.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.17) then
      C(1,1)=72.0_dp
      C(1,2)=8.0_dp
      C(1,3)=0.0_dp
      C(1,4)=8.0_dp
      C(2,1)=8.0_dp
      C(2,2)=8.0_dp
      C(2,3)=8.0_dp
      C(2,4)=0.0_dp
      C(3,1)=0.0_dp
      C(3,2)=8.0_dp
      C(3,3)=0.0_dp
      C(3,4)=-8.0_dp
      C(4,1)=8.0_dp
      C(4,2)=0.0_dp
      C(4,3)=-8.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.18) then
      C(1,1)=0.0_dp
      C(1,2)=8.0_dp
      C(1,3)=0.0_dp
      C(1,4)=-8.0_dp
      C(2,1)=8.0_dp
      C(2,2)=8.0_dp
      C(2,3)=8.0_dp
      C(2,4)=0.0_dp
      C(3,1)=0.0_dp
      C(3,2)=8.0_dp
      C(3,3)=72.0_dp
      C(3,4)=8.0_dp
      C(4,1)=-8.0_dp
      C(4,2)=0.0_dp
      C(4,3)=8.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.19) then
      C(1,1)=0.0_dp
      C(1,2)=-8.0_dp
      C(1,3)=0.0_dp
      C(1,4)=8.0_dp
      C(2,1)=-8.0_dp
      C(2,2)=0.0_dp
      C(2,3)=8.0_dp
      C(2,4)=0.0_dp
      C(3,1)=0.0_dp
      C(3,2)=8.0_dp
      C(3,3)=72.0_dp
      C(3,4)=8.0_dp
      C(4,1)=8.0_dp
      C(4,2)=0.0_dp
      C(4,3)=8.0_dp
      C(4,4)=8.0_dp
      endif
      if(n.eq.20) then
      C(1,1)=72.0_dp
      C(1,2)=8.0_dp
      C(1,3)=0.0_dp
      C(1,4)=8.0_dp
      C(2,1)=8.0_dp
      C(2,2)=0.0_dp
      C(2,3)=-8.0_dp
      C(2,4)=0.0_dp
      C(3,1)=0.0_dp
      C(3,2)=-8.0_dp
      C(3,3)=0.0_dp
      C(3,4)=8.0_dp
      C(4,1)=8.0_dp
      C(4,2)=0.0_dp
      C(4,3)=8.0_dp
      C(4,4)=8.0_dp
      endif
      if(n.eq.21) then
      C(1,1)=0.0_dp
      C(1,2)=8.0_dp/9.0_dp
      C(1,3)=-8.0_dp
      C(1,4)=-64.0_dp/9.0_dp
      C(2,1)=8.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=-64.0_dp/9.0_dp
      C(2,4)=-8.0_dp/9.0_dp
      C(3,1)=-8.0_dp
      C(3,2)=-64.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=8.0_dp/9.0_dp
      C(4,1)=-64.0_dp/9.0_dp
      C(4,2)=-8.0_dp/9.0_dp
      C(4,3)=8.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.22) then
      C(1,1)=64.0_dp
      C(1,2)=64.0_dp/9.0_dp
      C(1,3)=0.0_dp
      C(1,4)=64.0_dp/9.0_dp
      C(2,1)=64.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=-8.0_dp/9.0_dp
      C(2,4)=8.0_dp/9.0_dp
      C(3,1)=0.0_dp
      C(3,2)=-8.0_dp/9.0_dp
      C(3,3)=-8.0_dp
      C(3,4)=-8.0_dp/9.0_dp
      C(4,1)=64.0_dp/9.0_dp
      C(4,2)=8.0_dp/9.0_dp
      C(4,3)=-8.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.23) then
      C(1,1)=0.0_dp
      C(1,2)=64.0_dp/9.0_dp
      C(1,3)=8.0_dp
      C(1,4)=-8.0_dp/9.0_dp
      C(2,1)=64.0_dp/9.0_dp
      C(2,2)=64.0_dp/9.0_dp
      C(2,3)=64.0_dp/9.0_dp
      C(2,4)=0.0_dp
      C(3,1)=8.0_dp
      C(3,2)=64.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=-8.0_dp/9.0_dp
      C(4,1)=-8.0_dp/9.0_dp
      C(4,2)=0.0_dp
      C(4,3)=-8.0_dp/9.0_dp
      C(4,4)=-8.0_dp/9.0_dp
      endif
      if(n.eq.24) then
      C(1,1)=0.0_dp
      C(1,2)=-8.0_dp
      C(1,3)=0.0_dp
      C(1,4)=8.0_dp
      C(2,1)=-8.0_dp
      C(2,2)=0.0_dp
      C(2,3)=8.0_dp
      C(2,4)=0.0_dp
      C(3,1)=0.0_dp
      C(3,2)=8.0_dp
      C(3,3)=72.0_dp
      C(3,4)=8.0_dp
      C(4,1)=8.0_dp
      C(4,2)=0.0_dp
      C(4,3)=8.0_dp
      C(4,4)=8.0_dp
      endif
      if(n.eq.25) then
      C(1,1)=-8.0_dp
      C(1,2)=-8.0_dp/9.0_dp
      C(1,3)=0.0_dp
      C(1,4)=-8.0_dp/9.0_dp
      C(2,1)=-8.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=64.0_dp/9.0_dp
      C(2,4)=8.0_dp/9.0_dp
      C(3,1)=0.0_dp
      C(3,2)=64.0_dp/9.0_dp
      C(3,3)=64.0_dp
      C(3,4)=64.0_dp/9.0_dp
      C(4,1)=-8.0_dp/9.0_dp
      C(4,2)=8.0_dp/9.0_dp
      C(4,3)=64.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.26) then
      C(1,1)=0.0_dp
      C(1,2)=-64.0_dp/9.0_dp
      C(1,3)=-8.0_dp
      C(1,4)=8.0_dp/9.0_dp
      C(2,1)=-64.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=8.0_dp/9.0_dp
      C(2,4)=-8.0_dp/9.0_dp
      C(3,1)=-8.0_dp
      C(3,2)=8.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=-64.0_dp/9.0_dp
      C(4,1)=8.0_dp/9.0_dp
      C(4,2)=-8.0_dp/9.0_dp
      C(4,3)=-64.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.27) then
      C(1,1)=0.0_dp
      C(1,2)=64.0_dp/9.0_dp
      C(1,3)=8.0_dp
      C(1,4)=-8.0_dp/9.0_dp
      C(2,1)=64.0_dp/9.0_dp
      C(2,2)=64.0_dp/9.0_dp
      C(2,3)=64.0_dp/9.0_dp
      C(2,4)=0.0_dp
      C(3,1)=8.0_dp
      C(3,2)=64.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=-8.0_dp/9.0_dp
      C(4,1)=-8.0_dp/9.0_dp
      C(4,2)=0.0_dp
      C(4,3)=-8.0_dp/9.0_dp
      C(4,4)=-8.0_dp/9.0_dp
      endif
      if(n.eq.28) then
      C(1,1)=72.0_dp
      C(1,2)=8.0_dp
      C(1,3)=0.0_dp
      C(1,4)=8.0_dp
      C(2,1)=8.0_dp
      C(2,2)=0.0_dp
      C(2,3)=-8.0_dp
      C(2,4)=0.0_dp
      C(3,1)=0.0_dp
      C(3,2)=-8.0_dp
      C(3,3)=0.0_dp
      C(3,4)=8.0_dp
      C(4,1)=8.0_dp
      C(4,2)=0.0_dp
      C(4,3)=8.0_dp
      C(4,4)=8.0_dp
      endif
      if(n.eq.29) then
      C(1,1)=0.0_dp
      C(1,2)=-8.0_dp/9.0_dp
      C(1,3)=8.0_dp
      C(1,4)=64.0_dp/9.0_dp
      C(2,1)=-8.0_dp/9.0_dp
      C(2,2)=-8.0_dp/9.0_dp
      C(2,3)=-8.0_dp/9.0_dp
      C(2,4)=0.0_dp
      C(3,1)=8.0_dp
      C(3,2)=-8.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=64.0_dp/9.0_dp
      C(4,1)=64.0_dp/9.0_dp
      C(4,2)=0.0_dp
      C(4,3)=64.0_dp/9.0_dp
      C(4,4)=64.0_dp/9.0_dp
      endif
      if(n.eq.30) then
      C(1,1)=0.0_dp
      C(1,2)=8.0_dp/9.0_dp
      C(1,3)=-8.0_dp
      C(1,4)=-64.0_dp/9.0_dp
      C(2,1)=8.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=-64.0_dp/9.0_dp
      C(2,4)=-8.0_dp/9.0_dp
      C(3,1)=-8.0_dp
      C(3,2)=-64.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=8.0_dp/9.0_dp
      C(4,1)=-64.0_dp/9.0_dp
      C(4,2)=-8.0_dp/9.0_dp
      C(4,3)=8.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.31) then
      C(1,1)=-8.0_dp
      C(1,2)=-8.0_dp/9.0_dp
      C(1,3)=0.0_dp
      C(1,4)=-8.0_dp/9.0_dp
      C(2,1)=-8.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=64.0_dp/9.0_dp
      C(2,4)=8.0_dp/9.0_dp
      C(3,1)=0.0_dp
      C(3,2)=64.0_dp/9.0_dp
      C(3,3)=64.0_dp
      C(3,4)=64.0_dp/9.0_dp
      C(4,1)=-8.0_dp/9.0_dp
      C(4,2)=8.0_dp/9.0_dp
      C(4,3)=64.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.32) then
      C(1,1)=72.0_dp
      C(1,2)=8.0_dp
      C(1,3)=0.0_dp
      C(1,4)=8.0_dp
      C(2,1)=8.0_dp
      C(2,2)=8.0_dp
      C(2,3)=8.0_dp
      C(2,4)=0.0_dp
      C(3,1)=0.0_dp
      C(3,2)=8.0_dp
      C(3,3)=0.0_dp
      C(3,4)=-8.0_dp
      C(4,1)=8.0_dp
      C(4,2)=0.0_dp
      C(4,3)=-8.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.33) then
      C(1,1)=0.0_dp
      C(1,2)=-8.0_dp/9.0_dp
      C(1,3)=8.0_dp
      C(1,4)=64.0_dp/9.0_dp
      C(2,1)=-8.0_dp/9.0_dp
      C(2,2)=-8.0_dp/9.0_dp
      C(2,3)=-8.0_dp/9.0_dp
      C(2,4)=0.0_dp
      C(3,1)=8.0_dp
      C(3,2)=-8.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=64.0_dp/9.0_dp
      C(4,1)=64.0_dp/9.0_dp
      C(4,2)=0.0_dp
      C(4,3)=64.0_dp/9.0_dp
      C(4,4)=64.0_dp/9.0_dp
      endif
      if(n.eq.34) then
      C(1,1)=64.0_dp
      C(1,2)=64.0_dp/9.0_dp
      C(1,3)=0.0_dp
      C(1,4)=64.0_dp/9.0_dp
      C(2,1)=64.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=-8.0_dp/9.0_dp
      C(2,4)=8.0_dp/9.0_dp
      C(3,1)=0.0_dp
      C(3,2)=-8.0_dp/9.0_dp
      C(3,3)=-8.0_dp
      C(3,4)=-8.0_dp/9.0_dp
      C(4,1)=64.0_dp/9.0_dp
      C(4,2)=8.0_dp/9.0_dp
      C(4,3)=-8.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.35) then
      C(1,1)=0.0_dp
      C(1,2)=-64.0_dp/9.0_dp
      C(1,3)=-8.0_dp
      C(1,4)=8.0_dp/9.0_dp
      C(2,1)=-64.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=8.0_dp/9.0_dp
      C(2,4)=-8.0_dp/9.0_dp
      C(3,1)=-8.0_dp
      C(3,2)=8.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=-64.0_dp/9.0_dp
      C(4,1)=8.0_dp/9.0_dp
      C(4,2)=-8.0_dp/9.0_dp
      C(4,3)=-64.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.36) then
      C(1,1)=0.0_dp
      C(1,2)=8.0_dp
      C(1,3)=0.0_dp
      C(1,4)=-8.0_dp
      C(2,1)=8.0_dp
      C(2,2)=8.0_dp
      C(2,3)=8.0_dp
      C(2,4)=0.0_dp
      C(3,1)=0.0_dp
      C(3,2)=8.0_dp
      C(3,3)=72.0_dp
      C(3,4)=8.0_dp
      C(4,1)=-8.0_dp
      C(4,2)=0.0_dp
      C(4,3)=8.0_dp
      C(4,4)=0.0_dp
      endif


      end subroutine cc_qq_ttqqgg

      subroutine cc_qg_ttqqgg(n,C)
      integer, intent(in)  :: n
      real(dp), intent(out)  :: C(4,4)


      if(n.eq.1) then
      C(1,1)=0.0_dp
      C(1,2)=8.0_dp/9.0_dp
      C(1,3)=-8.0_dp
      C(1,4)=-64.0_dp/9.0_dp
      C(2,1)=8.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=-64.0_dp/9.0_dp
      C(2,4)=-8.0_dp/9.0_dp
      C(3,1)=-8.0_dp
      C(3,2)=-64.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=8.0_dp/9.0_dp
      C(4,1)=-64.0_dp/9.0_dp
      C(4,2)=-8.0_dp/9.0_dp
      C(4,3)=8.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.2) then
      C(1,1)=64.0_dp
      C(1,2)=64.0_dp/9.0_dp
      C(1,3)=0.0_dp
      C(1,4)=64.0_dp/9.0_dp
      C(2,1)=64.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=-8.0_dp/9.0_dp
      C(2,4)=8.0_dp/9.0_dp
      C(3,1)=0.0_dp
      C(3,2)=-8.0_dp/9.0_dp
      C(3,3)=-8.0_dp
      C(3,4)=-8.0_dp/9.0_dp
      C(4,1)=64.0_dp/9.0_dp
      C(4,2)=8.0_dp/9.0_dp
      C(4,3)=-8.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.3) then
      C(1,1)=0.0_dp
      C(1,2)=64.0_dp/9.0_dp
      C(1,3)=8.0_dp
      C(1,4)=-8.0_dp/9.0_dp
      C(2,1)=64.0_dp/9.0_dp
      C(2,2)=64.0_dp/9.0_dp
      C(2,3)=64.0_dp/9.0_dp
      C(2,4)=0.0_dp
      C(3,1)=8.0_dp
      C(3,2)=64.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=-8.0_dp/9.0_dp
      C(4,1)=-8.0_dp/9.0_dp
      C(4,2)=0.0_dp
      C(4,3)=-8.0_dp/9.0_dp
      C(4,4)=-8.0_dp/9.0_dp
      endif
      if(n.eq.4) then
      C(1,1)=0.0_dp
      C(1,2)=-8.0_dp
      C(1,3)=0.0_dp
      C(1,4)=8.0_dp
      C(2,1)=-8.0_dp
      C(2,2)=0.0_dp
      C(2,3)=8.0_dp
      C(2,4)=0.0_dp
      C(3,1)=0.0_dp
      C(3,2)=8.0_dp
      C(3,3)=72.0_dp
      C(3,4)=8.0_dp
      C(4,1)=8.0_dp
      C(4,2)=0.0_dp
      C(4,3)=8.0_dp
      C(4,4)=8.0_dp
      endif
      if(n.eq.5) then
      C(1,1)=-8.0_dp
      C(1,2)=-8.0_dp/9.0_dp
      C(1,3)=0.0_dp
      C(1,4)=-8.0_dp/9.0_dp
      C(2,1)=-8.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=64.0_dp/9.0_dp
      C(2,4)=8.0_dp/9.0_dp
      C(3,1)=0.0_dp
      C(3,2)=64.0_dp/9.0_dp
      C(3,3)=64.0_dp
      C(3,4)=64.0_dp/9.0_dp
      C(4,1)=-8.0_dp/9.0_dp
      C(4,2)=8.0_dp/9.0_dp
      C(4,3)=64.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.6) then
      C(1,1)=0.0_dp
      C(1,2)=-64.0_dp/9.0_dp
      C(1,3)=-8.0_dp
      C(1,4)=8.0_dp/9.0_dp
      C(2,1)=-64.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=8.0_dp/9.0_dp
      C(2,4)=-8.0_dp/9.0_dp
      C(3,1)=-8.0_dp
      C(3,2)=8.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=-64.0_dp/9.0_dp
      C(4,1)=8.0_dp/9.0_dp
      C(4,2)=-8.0_dp/9.0_dp
      C(4,3)=-64.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.7) then
      C(1,1)=0.0_dp
      C(1,2)=64.0_dp/9.0_dp
      C(1,3)=8.0_dp
      C(1,4)=-8.0_dp/9.0_dp
      C(2,1)=64.0_dp/9.0_dp
      C(2,2)=64.0_dp/9.0_dp
      C(2,3)=64.0_dp/9.0_dp
      C(2,4)=0.0_dp
      C(3,1)=8.0_dp
      C(3,2)=64.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=-8.0_dp/9.0_dp
      C(4,1)=-8.0_dp/9.0_dp
      C(4,2)=0.0_dp
      C(4,3)=-8.0_dp/9.0_dp
      C(4,4)=-8.0_dp/9.0_dp
      endif
      if(n.eq.8) then
      C(1,1)=72.0_dp
      C(1,2)=8.0_dp
      C(1,3)=0.0_dp
      C(1,4)=8.0_dp
      C(2,1)=8.0_dp
      C(2,2)=0.0_dp
      C(2,3)=-8.0_dp
      C(2,4)=0.0_dp
      C(3,1)=0.0_dp
      C(3,2)=-8.0_dp
      C(3,3)=0.0_dp
      C(3,4)=8.0_dp
      C(4,1)=8.0_dp
      C(4,2)=0.0_dp
      C(4,3)=8.0_dp
      C(4,4)=8.0_dp
      endif
      if(n.eq.9) then
      C(1,1)=0.0_dp
      C(1,2)=-8.0_dp/9.0_dp
      C(1,3)=8.0_dp
      C(1,4)=64.0_dp/9.0_dp
      C(2,1)=-8.0_dp/9.0_dp
      C(2,2)=-8.0_dp/9.0_dp
      C(2,3)=-8.0_dp/9.0_dp
      C(2,4)=0.0_dp
      C(3,1)=8.0_dp
      C(3,2)=-8.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=64.0_dp/9.0_dp
      C(4,1)=64.0_dp/9.0_dp
      C(4,2)=0.0_dp
      C(4,3)=64.0_dp/9.0_dp
      C(4,4)=64.0_dp/9.0_dp
      endif
      if(n.eq.10) then
      C(1,1)=0.0_dp
      C(1,2)=8.0_dp/9.0_dp
      C(1,3)=-8.0_dp
      C(1,4)=-64.0_dp/9.0_dp
      C(2,1)=8.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=-64.0_dp/9.0_dp
      C(2,4)=-8.0_dp/9.0_dp
      C(3,1)=-8.0_dp
      C(3,2)=-64.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=8.0_dp/9.0_dp
      C(4,1)=-64.0_dp/9.0_dp
      C(4,2)=-8.0_dp/9.0_dp
      C(4,3)=8.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.11) then
      C(1,1)=-8.0_dp
      C(1,2)=-8.0_dp/9.0_dp
      C(1,3)=0.0_dp
      C(1,4)=-8.0_dp/9.0_dp
      C(2,1)=-8.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=64.0_dp/9.0_dp
      C(2,4)=8.0_dp/9.0_dp
      C(3,1)=0.0_dp
      C(3,2)=64.0_dp/9.0_dp
      C(3,3)=64.0_dp
      C(3,4)=64.0_dp/9.0_dp
      C(4,1)=-8.0_dp/9.0_dp
      C(4,2)=8.0_dp/9.0_dp
      C(4,3)=64.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.12) then
      C(1,1)=72.0_dp
      C(1,2)=8.0_dp
      C(1,3)=0.0_dp
      C(1,4)=8.0_dp
      C(2,1)=8.0_dp
      C(2,2)=8.0_dp
      C(2,3)=8.0_dp
      C(2,4)=0.0_dp
      C(3,1)=0.0_dp
      C(3,2)=8.0_dp
      C(3,3)=0.0_dp
      C(3,4)=-8.0_dp
      C(4,1)=8.0_dp
      C(4,2)=0.0_dp
      C(4,3)=-8.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.13) then
      C(1,1)=0.0_dp
      C(1,2)=-8.0_dp/9.0_dp
      C(1,3)=8.0_dp
      C(1,4)=64.0_dp/9.0_dp
      C(2,1)=-8.0_dp/9.0_dp
      C(2,2)=-8.0_dp/9.0_dp
      C(2,3)=-8.0_dp/9.0_dp
      C(2,4)=0.0_dp
      C(3,1)=8.0_dp
      C(3,2)=-8.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=64.0_dp/9.0_dp
      C(4,1)=64.0_dp/9.0_dp
      C(4,2)=0.0_dp
      C(4,3)=64.0_dp/9.0_dp
      C(4,4)=64.0_dp/9.0_dp
      endif
      if(n.eq.14) then
      C(1,1)=64.0_dp
      C(1,2)=64.0_dp/9.0_dp
      C(1,3)=0.0_dp
      C(1,4)=64.0_dp/9.0_dp
      C(2,1)=64.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=-8.0_dp/9.0_dp
      C(2,4)=8.0_dp/9.0_dp
      C(3,1)=0.0_dp
      C(3,2)=-8.0_dp/9.0_dp
      C(3,3)=-8.0_dp
      C(3,4)=-8.0_dp/9.0_dp
      C(4,1)=64.0_dp/9.0_dp
      C(4,2)=8.0_dp/9.0_dp
      C(4,3)=-8.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.15) then
      C(1,1)=0.0_dp
      C(1,2)=-64.0_dp/9.0_dp
      C(1,3)=-8.0_dp
      C(1,4)=8.0_dp/9.0_dp
      C(2,1)=-64.0_dp/9.0_dp
      C(2,2)=0.0_dp
      C(2,3)=8.0_dp/9.0_dp
      C(2,4)=-8.0_dp/9.0_dp
      C(3,1)=-8.0_dp
      C(3,2)=8.0_dp/9.0_dp
      C(3,3)=0.0_dp
      C(3,4)=-64.0_dp/9.0_dp
      C(4,1)=8.0_dp/9.0_dp
      C(4,2)=-8.0_dp/9.0_dp
      C(4,3)=-64.0_dp/9.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.16) then
      C(1,1)=0.0_dp
      C(1,2)=8.0_dp
      C(1,3)=0.0_dp
      C(1,4)=-8.0_dp
      C(2,1)=8.0_dp
      C(2,2)=8.0_dp
      C(2,3)=8.0_dp
      C(2,4)=0.0_dp
      C(3,1)=0.0_dp
      C(3,2)=8.0_dp
      C(3,3)=72.0_dp
      C(3,4)=8.0_dp
      C(4,1)=-8.0_dp
      C(4,2)=0.0_dp
      C(4,3)=8.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.17) then
      C(1,1)=72.0_dp
      C(1,2)=8.0_dp
      C(1,3)=0.0_dp
      C(1,4)=8.0_dp
      C(2,1)=8.0_dp
      C(2,2)=8.0_dp
      C(2,3)=8.0_dp
      C(2,4)=0.0_dp
      C(3,1)=0.0_dp
      C(3,2)=8.0_dp
      C(3,3)=0.0_dp
      C(3,4)=-8.0_dp
      C(4,1)=8.0_dp
      C(4,2)=0.0_dp
      C(4,3)=-8.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.18) then
      C(1,1)=0.0_dp
      C(1,2)=8.0_dp
      C(1,3)=0.0_dp
      C(1,4)=-8.0_dp
      C(2,1)=8.0_dp
      C(2,2)=8.0_dp
      C(2,3)=8.0_dp
      C(2,4)=0.0_dp
      C(3,1)=0.0_dp
      C(3,2)=8.0_dp
      C(3,3)=72.0_dp
      C(3,4)=8.0_dp
      C(4,1)=-8.0_dp
      C(4,2)=0.0_dp
      C(4,3)=8.0_dp
      C(4,4)=0.0_dp
      endif
      if(n.eq.19) then
      C(1,1)=0.0_dp
      C(1,2)=-8.0_dp
      C(1,3)=0.0_dp
      C(1,4)=8.0_dp
      C(2,1)=-8.0_dp
      C(2,2)=0.0_dp
      C(2,3)=8.0_dp
      C(2,4)=0.0_dp
      C(3,1)=0.0_dp
      C(3,2)=8.0_dp
      C(3,3)=72.0_dp
      C(3,4)=8.0_dp
      C(4,1)=8.0_dp
      C(4,2)=0.0_dp
      C(4,3)=8.0_dp
      C(4,4)=8.0_dp
      endif
      if(n.eq.20) then
      C(1,1)=72.0_dp
      C(1,2)=8.0_dp
      C(1,3)=0.0_dp
      C(1,4)=8.0_dp
      C(2,1)=8.0_dp
      C(2,2)=0.0_dp
      C(2,3)=-8.0_dp
      C(2,4)=0.0_dp
      C(3,1)=0.0_dp
      C(3,2)=-8.0_dp
      C(3,3)=0.0_dp
      C(3,4)=8.0_dp
      C(4,1)=8.0_dp
      C(4,2)=0.0_dp
      C(4,3)=8.0_dp
      C(4,4)=8.0_dp
      endif

      end subroutine cc_qg_ttqqgg










end module
