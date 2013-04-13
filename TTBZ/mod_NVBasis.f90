MODULE ModNVBasis
use ModMisc
implicit none

public :: GramOrthog,NVBasis4,NVBasis3,NVBasis2,NVBasis1
public :: NVBasis1_Light,give1to4vect,give1to4vect_light

contains




SUBROUTINE GramOrthog(N,VIn,VOut)  ! N is the number of vectors that are returned
implicit none
integer :: N
double precision :: VIn(4-N,1:4),V1(1:4),V2(1:4),V3(1:4),V4(1:4)
double complex :: VOut(1:N,1:4),VAux(1:4,1:4),Norm
integer :: j


   if(N.eq.1) then
      V1(1:4) = VIn(1,1:4)
      V2(1:4) = VIn(2,1:4)
      V3(1:4) = VIn(3,1:4)

      V4(1)  = VIn(1,1) * (1.8d0)
      V4(2)  = VIn(1,2) *(-2.3d0)
      V4(3)  = VIn(1,3) * (1.3d0)
      V4(4)  = VIn(1,4) *(-0.9d0)

      V1(1)=1d0
      V1(2)=2d0
      V1(3)=0d0
      V1(4)=0d0

      V2(1)=1d0
      V2(2)=-1d0
      V2(3)=0d0
      V2(4)=0d0

      V3(1)=0d0
      V3(2)=1d0
      V3(3)=0d0
      V3(4)=1d0

      V4(1)=0d0
      V4(2)=1d0
      V4(3)=1d0
      V4(4)=0d0


   elseif(N.eq.2) then
      V1(1:4) = VIn(1,1:4)
      V2(1:4) = VIn(2,1:4)

      V3(1)  = (VIn(1,1)+VIn(2,1)) * 1.2d0
      V3(2)  = (VIn(1,2)+VIn(2,2)) * 1.2d0
      V3(3)  = (VIn(1,3)+VIn(2,3)) * 1.2d0
      V3(4)  = (VIn(1,4)+VIn(2,4)) * 1.2d0

      V4(1)  = (VIn(1,1)-VIn(2,1)) * 1d0
      V4(2)  = (VIn(1,2)-VIn(2,2)) * 1d0
      V4(3)  = (VIn(1,3)-VIn(2,3)) * 1d0
      V4(4)  = (VIn(1,4)-VIn(2,4)) * 1d0

   elseif(N.eq.3) then
      V1(1:4) = VIn(1,1:4)

      V2(1)  = VIn(1,1) * 1.8d0
      V2(2)  = VIn(1,2) *(-2.3d0)
      V2(3)  = VIn(1,3) * 1.3d0
      V2(4)  = VIn(1,4) *(-0.9d0)

      V3(1)  = VIn(1,1) *(-1.3d0)
      V3(2)  = VIn(1,2) *(-0.4d0)
      V3(3)  = VIn(1,3) * 1.0d0
      V3(4)  = VIn(1,4) * 1.2d0

      V4(1)  = VIn(1,1) * 1.1d0
      V4(2)  = VIn(1,2) * 1.4d0
      V4(3)  = VIn(1,3) *(-0.4d0)
      V4(4)  = VIn(1,4) * 1.2d0

   else
      call Error("SUBROUTINE GramOrthog")
   endif

   VAux(1,1:4) = V1(1:4)

   VAux(2,1:4) = V2(1:4)
   Norm = VAux(1,1:4).dot.VAux(1,1:4)
   if( cdabs(Norm) .lt. 1d-6 ) call Error("Too small Norm in GramOrthog, 1")
   VAux(2,1:4) = VAux(2,1:4) - (VAux(2,1:4).dot.VAux(1,1:4))*VAux(1,1:4)/Norm

   VAux(3,1:4) = V3(1:4)
   do j=1,2
        Norm = VAux(j,1:4).dot.VAux(j,1:4)
        if( cdabs(Norm) .lt. 1d-6 ) call Error("Too small Norm in GramOrthog, 2")
        VAux(3,1:4) =  VAux(3,1:4) - (VAux(3,1:4).dot.VAux(j,1:4)) * VAux(j,1:4)/Norm
   enddo

   VAux(4,1:4) = V4(1:4)
   do j=1,3
        Norm = VAux(j,1:4).dot.VAux(j,1:4)
        if( cdabs(Norm) .lt. 1d-6 ) call Error("Too small Norm in GramOrthog, 3")
        VAux(4,1:4) =  VAux(4,1:4) - (VAux(4,1:4).dot.VAux(j,1:4)) * VAux(j,1:4)/Norm
   enddo

   do j=1,N
      Norm = cdsqrt((VAux(4-j+1,1:4).dot.VAux(4-j+1,1:4)))
      if( cdabs(Norm) .lt. 1d-6 ) call Error("Too small Norm in GramOrthog, 4")
      VOut(j,1:4) = VAux(4-j+1,1:4)/Norm
   enddo

!    print *, VOut(1,1:4).dot.VOut(1,1:4)
!    print *, V1(1:4).dot.VOut(1,1:4)
!    print *, V2(1:4).dot.VOut(1,1:4)
!    print *, V3(1:4).dot.VOut(1,1:4)
!    print *, V4(1:4).dot.VOut(1,1:4)
!    stop

END SUBROUTINE





! SUBROUTINE GaussOrthog(N,VIn,VOut)
! use ModMisc
! implicit none
! integer :: N
! double precision :: VIn(1:4-N,1:4),VWork(1:4-N,1:4)
! double complex :: VOut(1:N,1:4),SolVec(1:4)
! double precision :: abslargest,tPr
! integer :: i,j,iLargest,jLargest,iPivo,jPivo
!
!    print *, "0"
!    write (*,"(F8.2,F8.2,F8.2,F8.2)") VIn(1,1:4)
!    write (*,"(F8.2,F8.2,F8.2,F8.2)") VIn(2,1:4)
!    print *, "---"
!
! ! search for largest element
!    abslargest = 0d0
!    do i=1,4-N
!       do j=1,4
!          if ( dabs(VIn(i,j)) .gt. abslargest ) then
!             abslargest = dabs(VIn(i,j))
!             iLargest = i
!             jLargest = j
!          endif
!       enddo
!    enddo
!
!
! !  full pivoting
!    do i=1,4-N
!       if( i.eq.1 ) then
!          iPivo=iLargest
!       elseif( i.eq.iLargest ) then
!          iPivo=1
!       else
!          iPivo=i
!       endif
!       do j=1,4
!          if( j.eq.1 ) then
!             jPivo=jLargest
!          elseif( j.eq.jLargest ) then
!             jPivo=1
!          else
!             jPivo=j
!          endif
!          VWork(i,j) = VIn(iPivo,jPivo)
!          if( jPivo.ne.1 ) VWork(i,j) = -VWork(i,j)   ! minus sign from metric
!       enddo
!    enddo
!
!    print *, "1"
!    write (*,"(F8.2,F8.2,F8.2,F8.2)") VWork(1,1:4)
!    write (*,"(F8.2,F8.2,F8.2,F8.2)") VWork(2,1:4)
!    print *, "---"
!
! !  Gauss elimination:
! !  1. generate leading 1
!    do j=2,4
!       VWork(1,j) = VWork(1,j)/VWork(1,1)
!    enddo
!    VWork(1,1) = 1d0
!
! !  2. generate zeros in 2nd column
!    do i=2,4-N
!       do j=2,4
!          VWork(i,j) = VWork(i,j) - VWork(1,j)*VWork(i,1)
!       enddo
!       VWork(i,1) = 0d0
!    enddo
!
! !  3. generate leading 1
!    do j=3,4
!       VWork(2,j) = VWork(2,j)/VWork(2,2)
!    enddo
!    VWork(2,2) = 1d0
!
!    print *, "2"
!    write (*,"(F8.2,F8.2,F8.2,F8.2)") VWork(1,1:4)
!    write (*,"(F8.2,F8.2,F8.2,F8.2)") VWork(2,1:4)
!    print *, "---"
!
! ! solution for 1st vector
!    SolVec(1) = VWork(1,2)*( VWork(2,3)+VWork(2,4) ) - VWork(1,3) - VWork(1,4)
!    SolVec(2) = -VWork(2,3) - VWork(2,4)
!    SolVec(3) = 1d0
!    SolVec(4) = 1d0
!
! !  re-ordering wrt. to Pivot element
! !    VOut(1,1) = SolVec(jLargest)
! !    do j=2,4
! !       if( j.eq.jLargest ) then
! !          jPivo = 1
! !       else
! !          jPivo = j
! !       endif
! !       VOut(1,j) = SolVec(jPivo)
! !    enddo
!    VOut(1,1:4)=SolVec(1:4)
!
! ! solution for 2nd vector
!    tPr =(-1d0 + VWork(1,3)**2 + VWork(1,3)*(VWork(1,4) - VWork(1,2)*(2d0*VWork(2,3) + VWork(2,4))) -  &
!     VWork(2,3)*(VWork(2,3) + VWork(2,4) + VWork(1,2)*(VWork(1,4) - VWork(1,2)*(VWork(2,3) + VWork(2,4)))))/   &
!   (-1d0 - VWork(2,4)*(VWork(2,3) + VWork(2,4)) + (VWork(1,4) - VWork(1,2)*VWork(2,4))*(VWork(1,3) + VWork(1,4) -   & VWork(1,2)*(VWork(2,3) + VWork(2,4))))
!
!
!    SolVec(1) = VWork(1,2)*(-VWork(2,3)+ tPr*VWork(2,4) ) + VWork(1,3) - tPr*VWork(1,4)
!    SolVec(2) = VWork(2,3) - tPr*VWork(2,4)
!    SolVec(3) = -1d0
!    SolVec(4) = tPr
!
! !  re-ordering wrt. to Pivot element
! !    VOut(2,1) = SolVec(jLargest)
! !    do j=2,4
! !       if( j.eq.jLargest ) then
! !          jPivo = 1
! !       else
! !          jPivo = j
! !       endif
! !       VOut(2,j) = SolVec(jPivo)
! !    enddo
!    VOut(2,1:4)=SolVec(1:4)
!
! ! re-pivoting zerstört othogon. für vout 1 und 2!
! ! lösung: 4 versch. lösungen für tpr implementieren!
!
!    print *, VOut(1,1:4).dot.VOut(2,1:4)
!    print *, VIn(1,1:4).dot.VOut(1,1:4)
!    print *, VIn(1,1:4).dot.VOut(2,1:4)
!    print *, VIn(2,1:4).dot.VOut(1,1:4)
!    print *, VIn(2,1:4).dot.VOut(2,1:4)
! stop
! END SUBROUTINE






SUBROUTINE NVBasis4(KMom,VMom)
implicit none
double precision :: KMom(1:4,1:4),VMom(1:4,1:4)
double precision :: GramDet4
double precision :: s11,s22,s33,s44,s12,s13,s14,s23,s24,s34

   s11 = KMom(1,1:4) .dot. KMom(1,1:4)
   s22 = KMom(2,1:4) .dot. KMom(2,1:4)
   s33 = KMom(3,1:4) .dot. KMom(3,1:4)
   s44 = KMom(4,1:4) .dot. KMom(4,1:4)
   s12 = KMom(1,1:4) .dot. KMom(2,1:4)
   s13 = KMom(1,1:4) .dot. KMom(3,1:4)
   s14 = KMom(1,1:4) .dot. KMom(4,1:4)
   s23 = KMom(2,1:4) .dot. KMom(3,1:4)
   s24 = KMom(2,1:4) .dot. KMom(4,1:4)
   s34 = KMom(3,1:4) .dot. KMom(4,1:4)

   GramDet4 = s14**2*(s23**2 - s22*s33) + 2d0*s14*(-(s13*s23*s24) + s12*s24*s33 + s13*s22*s34 - s12*s23*s34)     &
            + s13**2*(s24**2 - s22*s44) + 2d0*s12*s13*(-(s24*s34) + s23*s44) + s12**2*(s34**2 - s33*s44)         &
            - s11*(s24**2*s33 - 2d0*s23*s24*s34 + s23**2*s44 + s22*(s34**2 - s33*s44))


   VMom(1,1:4) = KMom(4,1:4)*(-(s13*s23*s24) + s12*s24*s33 + s14*(s23**2 - s22*s33) + s13*s22*s34 - s12*s23*s34)/GramDet4     &
               + KMom(3,1:4)*(-(s14*s23*s24) + s13*s24**2 + s14*s22*s34 -s12*s24*s34 - s13*s22*s44 + s12*s23*s44)/GramDet4    &
               + KMom(2,1:4)*(s14*s24*s33 - s14*s23*s34 - s13*s24*s34 + s12*s34**2 + s13*s23*s44 - s12*s33*s44)/GramDet4      &
               + KMom(1,1:4)*(-(s24**2*s33) + 2d0*s23*s24*s34 - s22*s34**2 - s23**2*s44 + s22*s33*s44)/GramDet4

   VMom(2,1:4) = KMom(4,1:4)*(s13**2*s24 + s12*s14*s33 - s11*s24*s33 + s11*s23*s34 - s13*(s14*s23 + s12*s34))/GramDet4        &
               + KMom(3,1:4)*(s14**2*s23 + s11*s24*s34 - s14*(s13*s24 + s12*s34) + s12*s13*s44 - s11*s23*s44)/GramDet4        &
               + KMom(2,1:4)*(-(s14**2*s33) + 2d0*s13*s14*s34 - s11*s34**2 - s13**2*s44 + s11*s33*s44)/GramDet4               &
               + KMom(1,1:4)*(s14*s24*s33 - s14*s23*s34 - s13*s24*s34 + s12*s34**2 + s13*s23*s44 - s12*s33*s44)/GramDet4

   VMom(3,1:4) = KMom(4,1:4)*(s13*s14*s22 - s12*s14*s23 - s12*s13*s24 + s11*s23*s24 + s12**2*s34 - s11*s22*s34)/GramDet4      &
               + KMom(3,1:4)*(-(s14**2*s22) + 2d0*s12*s14*s24 - s11*s24**2 - s12**2*s44 + s11*s22*s44)/GramDet4               &
               + KMom(2,1:4)*(s14**2*s23 + s11*s24*s34 - s14*(s13*s24 + s12*s34) + s12*s13*s44 - s11*s23*s44)/GramDet4        &
               + KMom(1,1:4)*(-(s14*s23*s24) + s13*s24**2 + s14*s22*s34 - s12*s24*s34 - s13*s22*s44 + s12*s23*s44)/GramDet4

   VMom(4,1:4) = KMom(4,1:4)*(-(s13**2*s22) + 2d0*s12*s13*s23 - s11*s23**2 - s12**2*s33 + s11*s22*s33)/GramDet4               &
               + KMom(3,1:4)*(s13*s14*s22 - s12*s14*s23 - s12*s13*s24 + s11*s23*s24 + s12**2*s34 - s11*s22*s34)/GramDet4      &
               + KMom(2,1:4)*(s13**2*s24 + s12*s14*s33 - s11*s24*s33 + s11*s23*s34 - s13*(s14*s23 + s12*s34))/GramDet4        &
               + KMom(1,1:4)*(-(s13*s23*s24) + s12*s24*s33 + s14*(s23**2 - s22*s33) + s13*s22*s34 - s12*s23*s34)/GramDet4
END SUBROUTINE





SUBROUTINE NVBasis3(KMom,VMom,NMom)
implicit none
double precision :: KMom(1:3,1:4),VMom(1:3,1:4)
double complex :: NMom(1,1:4)
double precision :: GramDet3
double precision :: s11,s22,s33,s12,s13,s23

double precision :: det,ko1,ko2,ko3
double complex :: VAux(1:4),VAux2(1:4),VAux3(1:4),Norm,a1v,a2v,a3v


   s11 = KMom(1,1:4) .dot. KMom(1,1:4)
   s22 = KMom(2,1:4) .dot. KMom(2,1:4)
   s33 = KMom(3,1:4) .dot. KMom(3,1:4)
   s12 = KMom(1,1:4) .dot. KMom(2,1:4)
   s13 = KMom(1,1:4) .dot. KMom(3,1:4)
   s23 = KMom(2,1:4) .dot. KMom(3,1:4)

   GramDet3 = s23*(2d0*s12*s13-s11*s23)-s22*s13**2+s33*(s11*s22-s12**2)

   VMom(1,1:4) = KMom(1,1:4)*(-s23**2+s22*s33)/GramDet3      + KMom(2,1:4)*(s23*s13-s12*s33)/GramDet3   + KMom(3,1:4)*(-s22*s13+s12*s23)/GramDet3
   VMom(2,1:4) = KMom(1,1:4)*(s13*s23 - s12*s33)/GramDet3    + KMom(2,1:4)*(-s13**2 + s11*s33)/GramDet3 + KMom(3,1:4)*(s12*s13 - s11*s23)/GramDet3
   VMom(3,1:4) = KMom(1,1:4)*(-(s13*s22) + s12*s23)/GramDet3 + KMom(2,1:4)*(s12*s13 - s11*s23)/GramDet3 + KMom(3,1:4)*(-s12**2 + s11*s22)/GramDet3

!    my algorithm
!    call GramOrthog2(1,VMom(1:3,1:4),NMom(1:1,1:4))

!       Kirill's algorithm
        det = s12**2 -s11*s22
        ko1 = s22/det
        ko2 = s11/det
        ko3 =-s12/det

        VAux(1:4) = KMom(1,1:4)*( ko1*s13+ko3*s23 ) + KMom(2,1:4)*( ko2*s23+ko3*s13 ) + KMom(3,1:4)
        Norm = cdsqrt( VAux.dot.VAux )
        VAux(1:4) = Vaux(1:4)/Norm   ! orthog. to k1 and k2


        VAux2(1)=dcmplx(0.5d0)
        VAux2(2)=dcmplx(2.7d0)
        VAux2(3)=dcmplx(3.2d0)
        VAux2(4)=dcmplx(4.3d0)
!         VAux2(1)=dcmplx(KMom(1,1)*(+1.5d0))
!         VAux2(2)=dcmplx(KMom(1,1)*(-0.8d0))
!         VAux2(3)=dcmplx(KMom(1,1)*(-1.1d0))
!         VAux2(4)=dcmplx(KMom(1,1)*(+0.5d0))

        a1v = KMom(1,1:4).dot.VAux2
        a2v = KMom(2,1:4).dot.VAux2
        VAux3(1:4) = VAux2(1:4) + KMom(1,1:4)*( ko1*a1v+ ko3*a2v ) + KMom(2,1:4)*( ko2*a2v+ko3*a1v )  ! orthog. to k1 and k2

        a3v = VAux3.dot.VAux
        VAux2(1:4) = VAux3(1:4) - a3v*VAux(1:4)
        Norm = cdsqrt( VAux2.dot.VAux2 )

!!!!
        if( dreal(Norm).eq.0d0 .and. dimag(Norm).lt.0d0 ) Norm=Norm*(-1d0,0d0)        ! CONVERT TO KIRILL's F77-CONVENTION
!!!!
        NMom(1,1:4) = VAux2(1:4)/Norm

END SUBROUTINE




! SUBROUTINE NVBasis3_Q(KMom,NMom)
! use ModMisc
! implicit none
! double precision :: KMom(1:3,1:4)
! complex(16) :: NMom(1,1:4),a1v,a2v,a3v,VAux2(1:4),VAux3(1:4), VAux(1:4),Norm
! real(16) :: GramDet3
! real(16) :: s11,s22,s33,s12,s13,s23
!
! real(16) :: det,ko1,ko2,ko3
!
!
!    s11 = qext(KMom(1,1:4) .dot. KMom(1,1:4))
!    s22 = qext(KMom(2,1:4) .dot. KMom(2,1:4))
!    s33 = qext(KMom(3,1:4) .dot. KMom(3,1:4))
!    s12 = qext(KMom(1,1:4) .dot. KMom(2,1:4))
!    s13 = qext(KMom(1,1:4) .dot. KMom(3,1:4))
!    s23 = qext(KMom(2,1:4) .dot. KMom(3,1:4))
!
!    GramDet3 = s23*(2q0*s12*s13-s11*s23)-s22*s13**2+s33*(s11*s22-s12**2)
!
! !    my algorithm
! !    call GramOrthog2(1,VMom(1:3,1:4),NMom(1:1,1:4))
!
! !       Kirill's algorithm
!         det = s12**2 -s11*s22
!         ko1 = s22/det
!         ko2 = s11/det
!         ko3 =-s12/det
!
!         VAux(1:4) = KMom(1,1:4)*( ko1*s13+ko3*s23 ) + KMom(2,1:4)*( ko2*s23+ko3*s13 ) + KMom(3,1:4)
!         Norm = cqsqrt( VAux.dot.VAux )
!         VAux(1:4) = Vaux(1:4)/Norm   ! orthog. to k1 and k2
!
!
!         VAux2(1)=(0.5q0)
!         VAux2(2)=(2.7q0)
!         VAux2(3)=(3.2q0)
!         VAux2(4)=(4.3q0)
! !         VAux2(1)=dcmplx(KMom(1,1)*(+1.5d0))
! !         VAux2(2)=dcmplx(KMom(1,1)*(-0.8d0))
! !         VAux2(3)=dcmplx(KMom(1,1)*(-1.1d0))
! !         VAux2(4)=dcmplx(KMom(1,1)*(+0.5d0))
!
!         a1v = qcmplx(KMom(1,1:4)).dot.VAux2
!         a2v = qcmplx(KMom(2,1:4)).dot.VAux2
!         VAux3(1:4) = VAux2(1:4) + qcmplx(KMom(1,1:4),0q0)*( ko1*a1v+ ko3*a2v ) + qcmplx(KMom(2,1:4),0q0)*( ko2*a2v+ko3*a1v )  ! orthog. to k1 and k2
!
!         a3v = VAux3.dot.VAux
!         VAux2(1:4) = VAux3(1:4) - a3v*VAux(1:4)
!         Norm = cqsqrt( VAux2.dot.VAux2 )
!
! !!!!
!         if( qreal(Norm).eq.0d0 .and. qimag(Norm).lt.0d0 ) Norm=Norm*(-1q0,0q0)        ! CONVERT TO KIRILL's F77-CONVENTION
! !!!!
!         NMom(1,1:4) = VAux2(1:4)/Norm
!
! END SUBROUTINE







SUBROUTINE NVBasis2(KMom,VMom,NMom)
implicit none
double precision :: KMom(1:2,1:4),VMom(1:2,1:4)
double complex :: NMom(1:2,1:4)
double precision :: GramDet2
double precision :: s11,s22,s12

double precision :: det,ko1,ko2,ko3
double complex :: VAux(1:4),VAux2(1:4),VAux3(1:4),Norm,a1v,a2v,a3v

   s11 = KMom(1,1:4) .dot. KMom(1,1:4)
   s22 = KMom(2,1:4) .dot. KMom(2,1:4)
   s12 = KMom(1,1:4) .dot. KMom(2,1:4)

   GramDet2 = -s12**2 + s11*s22

   VMom(1,1:4) =  KMom(1,1:4)*s22/GramDet2 - KMom(2,1:4)*s12/GramDet2
   VMom(2,1:4) = -KMom(1,1:4)*s12/GramDet2 + KMom(2,1:4)*s11/GramDet2

!     call GaussOrthog(2,VMom(1:2,1:4),NMom(1:2,1:4))



!       use Kirill's algorithm
        det = -GramDet2
        ko1 = s22/det
        ko2 = s11/det
        ko3 =-s12/det

        VAux(1)=dcmplx(1.3d0)
        VAux(2)=dcmplx(1.7d0)
        VAux(3)=dcmplx(2.4d0)
        VAux(4)=dcmplx(3.5d0)
!         VAux(1)=dcmplx(KMom(1,1)*(+1.5d0))
!         VAux(2)=dcmplx(KMom(1,1)*(-0.8d0))
!         VAux(3)=dcmplx(KMom(1,1)*(-1.1d0))
!         VAux(4)=dcmplx(KMom(1,1)*(+0.5d0))

        a1v = KMom(1,1:4).dot.VAux(1:4)
        a2v = KMom(2,1:4).dot.VAux(1:4)
        VAux(1:4) = VAux(1:4) + KMom(1,1:4)*( ko1*a1v+ ko3*a2v ) + KMom(2,1:4)*( ko2*a2v+ko3*a1v )
        Norm = cdsqrt( VAux(1:4).dot.VAux(1:4) )

!!!!
        if( dreal(Norm).eq.0d0 .and. dimag(Norm).lt.0d0 ) Norm=Norm*(-1d0,0d0)        ! CONVERT TO KIRILL's F77-CONVENTION
!!!!
        NMom(1,1:4) = VAux(1:4)/Norm

        Vaux(1)=dcmplx(2.1d0)
        Vaux(2)=dcmplx(1.2d0)
        Vaux(3)=dcmplx(3.4d0)
        Vaux(4)=dcmplx(0.5d0)
!         VAux(1)=dcmplx(KMom(2,1)*(-1.2d0))
!         VAux(2)=dcmplx(KMom(2,1)*(+0.9d0))
!         VAux(3)=dcmplx(KMom(2,1)*(-1.7d0))
!         VAux(4)=dcmplx(KMom(2,1)*(-0.3d0))

        a1v = KMom(1,1:4).dot.VAux(1:4)
        a2v = KMom(2,1:4).dot.VAux(1:4)
        VAux(1:4) = VAux(1:4) + KMom(1,1:4)*( ko1*a1v+ ko3*a2v ) + KMom(2,1:4)*( ko2*a2v+ko3*a1v )

        a3v = VAux(1:4).dot.NMom(1,1:4)
        Norm = NMom(1,1:4).dot.NMom(1,1:4)

        VAux(1:4) = VAux(1:4) - a3v*NMom(1,1:4)/Norm
        Norm = cdsqrt( VAux.dot.VAux )


!!!!
        if( dreal(Norm).eq.0d0 .and. dimag(Norm).lt.0d0 ) Norm=Norm*(-1d0,0d0)        ! CONVERT TO KIRILL's F77-CONVENTION
!!!!
        NMom(2,1:4) = VAux(1:4)/Norm
END SUBROUTINE




SUBROUTINE NVBasis1(KMom,VMom,NMom)
implicit none
double precision :: KMom(1,1:4),VMom(1,1:4)
double complex :: NMom(1:3,1:4)
double precision :: GramDet1

double precision :: det,ko1,ko2,ko3
double complex :: VAux(1:4),VAux2(1:4),VAux3(1:4),Norm,a1v,a2v,a3v,s11,s12,s22


   GramDet1 = KMom(1,1:4) .dot. KMom(1,1:4)

   VMom(1,1:4) =  KMom(1,1:4)/GramDet1

!    call GramOrthog(3,VMom(1,1:4),NMom(1:3,1:4))



        VAux(1)=dcmplx(1.3d0)
        VAux(2)=dcmplx(1.7d0)
        VAux(3)=dcmplx(2.3d0)
        VAux(4)=dcmplx(3.4d0)
!         VAux(1)=dcmplx(KMom(1,1)*(+1.5d0))
!         VAux(2)=dcmplx(KMom(1,1)*(-0.8d0))
!         VAux(3)=dcmplx(KMom(1,1)*(-1.1d0))
!         VAux(4)=dcmplx(KMom(1,1)*(+0.5d0))

        a1v = VMom(1,1:4).dot.VAux(1:4)
        VAux(1:4) = VAux(1:4) - KMom(1,1:4)*a1v
        Norm = cdsqrt( VAux(1:4).dot.VAux(1:4) )

!!!!
        if( dreal(Norm).eq.0d0 .and. dimag(Norm).lt.0d0 ) Norm=Norm*(-1d0,0d0)        ! CONVERT TO KIRILL's F77-CONVENTION
!!!!
        NMom(1,1:4) = VAux(1:4)/Norm

        s11 = VMom(1,1:4).dot.VMom(1,1:4)
        s12 = NMom(1,1:4).dot.VMom(1,1:4)
        s22 = NMom(1,1:4).dot.NMom(1,1:4)

        det = s12**2 - s11*s22
        ko1 = s22/det
        ko2 = s11/det
        ko3 =-s12/det

        Vaux(1)=dcmplx(1.5d0)
        Vaux(2)=dcmplx(13.2d0)
        Vaux(3)=dcmplx(14.1d0)
        Vaux(4)=dcmplx(9.1d0)
!         VAux(1)=dcmplx(KMom(1,1)*(-1.2d0))
!         VAux(2)=dcmplx(KMom(1,1)*(+0.9d0))
!         VAux(3)=dcmplx(KMom(1,1)*(-1.7d0))
!         VAux(4)=dcmplx(KMom(1,1)*(-0.3d0))

        a1v = VMom(1,1:4).dot.VAux
        a2v = NMom(1,1:4).dot.VAux
        VAux(1:4) = VAux(1:4) + VMom(1,1:4)*( ko1*a1v+ ko3*a2v ) + NMom(1,1:4)*( ko2*a2v+ko3*a1v )
        Norm = cdsqrt( VAux.dot.VAux )

!!!!
        if( dreal(Norm).eq.0d0 .and. dimag(Norm).lt.0d0 ) Norm=Norm*(-1d0,0d0)        ! CONVERT TO KIRILL's F77-CONVENTION
!!!!
        NMom(2,1:4) = VAux(1:4)/Norm

        Vaux2(1)=dcmplx(2.4d0)
        Vaux2(2)=dcmplx(1.3d0)
        Vaux2(3)=dcmplx(3.45d0)
        Vaux2(4)=dcmplx(0.56d0)
!         VAux2(1)=dcmplx(KMom(1,1)*(+1.2d0))
!         VAux2(2)=dcmplx(KMom(1,1)*(-0.9d0))
!         VAux2(3)=dcmplx(KMom(1,1)*(+0.6d0))
!         VAux2(4)=dcmplx(KMom(1,1)*(-0.3d0))

        a1v = VMom(1,1:4).dot.VAux2
        a2v = NMom(1,1:4).dot.VAux2
        VAux(1:4) = VAux2(1:4) + VMom(1,1:4)*( ko1*a1v+ ko3*a2v ) + NMom(1,1:4)*( ko2*a2v+ko3*a1v )

        a3v = VAux(1:4).dot.NMom(2,1:4)
        Norm = NMom(2,1:4).dot.NMom(2,1:4)

        VAux(1:4) = VAux(1:4) - a3v*NMom(2,1:4)/Norm
        Norm = cdsqrt( VAux.dot.VAux )

!!!!
        if( dreal(Norm).eq.0d0 .and. dimag(Norm).lt.0d0 ) Norm=Norm*(-1d0,0d0)        ! CONVERT TO KIRILL's F77-CONVENTION
!!!!
        NMom(3,1:4) = VAux(1:4)/Norm

END SUBROUTINE




SUBROUTINE NVBasis1_Light(K1,VMom,NMom)
implicit none
integer i
double precision K1(1:4),VMom(1,1:4)
double complex NMom(1:3,1:4)
double complex v3(4),v4(4)
double complex k1d(4)
double complex vaux(4),vaux1(4)
double complex d2,ko3
double complex dnorm
double complex a1v,a2v,a33,av3

      k1d(1)   = dcmplx( k1(1) )
      k1d(2:4) = dcmplx(-k1(2:4) )

      call sc(4,k1d,dcmplx(k1),d2)

      if (zabs(d2).lt.1D-8) then
      k1d(1)=dcmplx(1d0,0d0)
      k1d(2)=dcmplx(1d0/dsqrt(3d0),0d0)
      k1d(3)=dcmplx(1d0/dsqrt(3d0),0d0)
      k1d(4)=dcmplx(1d0/dsqrt(3d0),0d0)
      endif

      call sc(4,k1d,dcmplx(k1),d2)

      if (zabs(d2).lt.1D-8) then
      k1d(1)=dcmplx(1d0,0d0)
      k1d(2)=dcmplx(1d0/dsqrt(3d0),0d0)
      k1d(3)=dcmplx(dsqrt(2d0)/dsqrt(3d0),0d0)
      k1d(4)=dcmplx(0d0,0d0)
      endif

      call sc(4,k1d,dcmplx(k1),d2)

      if (zabs(d2).lt.1D-8) then
         call Error("can not choose a dual vector")
      endif


      ko3=-1D0/d2

      vaux(1)=dcmplx(1.5d0,0d0)
      vaux(2)=dcmplx(13.2d0,0d0)
      vaux(3)=dcmplx(14.1d0,0d0)
      vaux(4)=dcmplx(9.1d0,0d0)

      call sc(4,dcmplx(k1),vaux,a1v)
      call sc(4,k1d,vaux,a2v)

      vaux(1:4)=vaux(1:4)+ko3*(a1v*k1d(1:4)+a2v*k1(1:4))

      call sc(4,vaux,vaux,dnorm)
      dnorm=zsqrt(dnorm)
      v3(1:4)=vaux(1:4)/dnorm


      !obtaining the forth vector
      vaux1(1)=dcmplx(2.4d0,0d0)
      vaux1(2)=dcmplx(1.35d0,0d0)
      vaux1(3)=dcmplx(3.61d0,0d0)
      vaux1(4)=dcmplx(0.54d0,0d0)

      call sc(4,dcmplx(k1),vaux1,a1v)
      call sc(4,k1d,vaux1,a2v)

      vaux(1:4)=vaux1(1:4)+ko3*(a1v*k1d(1:4)+a2v*k1(1:4))
      call sc(4,vaux,v3,av3)
      call sc(4,v3,v3,a33)

      vaux(1:4)=vaux(1:4)-av3/a33*v3(1:4)
      call sc(4,vaux,vaux,dnorm)
      dnorm=zsqrt(dnorm)

      v4(1:4)=vaux(1:4)/dnorm

      VMom(1,1:4)=k1(1:4)
      NMom(1,1:4)=k1d(1:4)
      NMom(2,1:4)=v3(1:4)
      NMom(3,1:4)=v4(1:4)

return
END SUBROUTINE






!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! the procedure returns 4 vectors from a not-light like vector
! the first vector, v1, is parallel to k1 and it is not normalized
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        SUBROUTINE give1to4vect(k1,v)
        implicit none
        integer i
        complex(8) k1(*),v(4,4),d1
        double complex v1(4),v2(4),v3(4),v4(4)
        double complex vaux(4),vaux1(4)
        double complex v2aux(4)
        double complex d2,ko1,ko2,ko3
        double complex a11,a12,a22,dnorm
        double complex a1v,a2v,a3v,a33

        call sc(4,k1,k1,d1)

        do i=1,4
        v1(i)=k1(i)/d1
        enddo

        v2aux(1)=dcmplx(1.3d0,0d0)
        v2aux(2)=dcmplx(1.7d0,0d0)
        v2aux(3)=dcmplx(2.3d0,0d0)
        v2aux(4)=dcmplx(3.4d0,0d0)

        call sc(4,v2aux,v1,a1v)

         do i=1,4
         v2aux(i)=v2aux(i)-a1v*k1(i)
         enddo

         call sc(4,v2aux,v2aux,dnorm)
         dnorm=zsqrt(dnorm)

         do i=1,4
         v2(i) = v2aux(i)/dnorm
         enddo

         call sc(4,v1,v2,a12)
         call sc(4,v1,v1,a11)
         call sc(4,v2,v2,a22)

         d2=a12**2 - a11*a22
         ko1= a22/d2
         ko2= a11/d2
         ko3=-a12/d2

         vaux(1)=dcmplx(1.5d0,0d0)
         vaux(2)=dcmplx(13.2d0,0d0)
         vaux(3)=dcmplx(14.1d0,0d0)
         vaux(4)=dcmplx(9.1d0,0d0)

         call sc(4,v1,vaux,a1v)
         call sc(4,v2,vaux,a2v)

         do i=1,4
         vaux(i)=vaux(i)+ko1*a1v*v1(i)+ko2*a2v*v2(i)+ ko3*(a1v*v2(i)+a2v*v1(i))
         enddo

          call sc(4,vaux,vaux,dnorm)
          dnorm=zsqrt(dnorm)



          do i=1,4
          v3(i)=vaux(i)/dnorm
          enddo

          vaux1(1)=dcmplx(2.4d0,0d0)
          vaux1(2)=dcmplx(1.3d0,0d0)
          vaux1(3)=dcmplx(3.45d0,0d0)
          vaux1(4)=dcmplx(0.56d0,0d0)

          call sc(4,v1,vaux1,a1v)
          call sc(4,v2,vaux1,a2v)

          do i=1,4
          vaux(i)=vaux1(i)+ko1*a1v*v1(i)+ko2*a2v*v2(i) + ko3*(a1v*v2(i)+a2v*v1(i))
          enddo

          call sc(4,vaux,v3,a3v)
          call sc(4,v3,v3,a33)

          do i=1,4
          vaux(i)=vaux(i) - a3v/a33*v3(i)
          enddo

          call sc(4,vaux,vaux,dnorm)
          dnorm=zsqrt(dnorm)

          do i=1,4
          v4(i)=vaux(i)/dnorm
          enddo

          do i=1,4
          v(1,i) = v1(i)
          v(2,i) = v2(i)
          v(3,i) = v3(i)
          v(4,i) = v4(i)
          enddo

          return
          END SUBROUTINE


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! for a light-like vector, this procedure returns the ``dual'' light-like
! vector -- k1d and two unit vectors (v3,v4) in the plane transverse to k1 & k1d
! the first vector that is returned is k1
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        SUBROUTINE give1to4vect_light(k1,v)
        implicit none
        integer i
        complex(8) k1(1:4),k1d(4),v(4,4),d2
        complex(8) v3(4),v4(4)
        complex(8) vaux(4),vaux1(4)
        complex(8) ko3
        complex(8) dnorm
        complex(8) a1v,a2v,a33,av3

        k1d(1)=k1(1)
        k1d(2:4) = -k1(2:4)

        call sc(4,k1d,k1,d2)

        if (zabs(d2).lt.1D-8) then
          k1d(1)=dcmplx(1d0,0d0)
          k1d(2)=dcmplx(1d0/dsqrt(3d0),0d0)
          k1d(3)=dcmplx(1d0/dsqrt(3d0),0d0)
          k1d(4)=dcmplx(1d0/dsqrt(3d0),0d0)
        endif

        call sc(4,k1d,k1,d2)

        if (zabs(d2).lt.1D-8) then
          k1d(1)=dcmplx(1d0,0d0)
          k1d(2)=dcmplx(1d0/dsqrt(3d0),0d0)
          k1d(3)=dcmplx(dsqrt(2d0)/dsqrt(3d0),0d0)
          k1d(4)=dcmplx(0d0,0d0)
        endif

        call sc(4,k1d,k1,d2)

        if (zabs(d2).lt.1D-8) then
            call Error("can not choose a dual vector")
        endif


        ko3=-1D0/d2

        vaux(1)=dcmplx(1.5d0,0d0)
        vaux(2)=dcmplx(13.2d0,0d0)
        vaux(3)=dcmplx(14.1d0,0d0)
        vaux(4)=dcmplx(9.1d0,0d0)

        call sc(4,k1,vaux,a1v)
        call sc(4,k1d,vaux,a2v)

        vaux(1:4)=vaux(1:4)+ko3*(a1v*k1d(1:4)+a2v*k1(1:4))

        call sc(4,vaux,vaux,dnorm)
        dnorm=zsqrt(1d0*dnorm)
        v3(1:4)=vaux(1:4)/dnorm


!        obtaining the forth vector
         vaux1(1)=dcmplx(2.4d0,0d0)
         vaux1(2)=dcmplx(1.35d0,0d0)
         vaux1(3)=dcmplx(3.61d0,0d0)
         vaux1(4)=dcmplx(0.54d0,0d0)

         call sc(4,k1,vaux1,a1v)
         call sc(4,k1d,vaux1,a2v)

         vaux(1:4)=vaux1(1:4)+ko3*(a1v*k1d(1:4)+a2v*k1(1:4))
         call sc(4,vaux,v3,av3)
         call sc(4,v3,v3,a33)

         vaux(1:4)=vaux(1:4)-av3/a33*v3(1:4)
         call sc(4,vaux,vaux,dnorm)
         dnorm=zsqrt(1d0*dnorm)

         v4(1:4)=vaux(1:4)/dnorm

         v(1,1:4)=k1(1:4)
         v(2,1:4)=k1d(1:4)
         v(3,1:4)=v3(1:4)
         v(4,1:4)=v4(1:4)

        return
        END SUBROUTINE













END MODULE
