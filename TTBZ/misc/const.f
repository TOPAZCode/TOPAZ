      module consts_dp
      use types
      implicit none
      private


c-----algebraic constants

      double precision, public, parameter :: pi =
     ,3.141592653589793238462643383279502884197d0

       double precision, public, parameter :: one = 1.0d0,
     , mone = -1.0d0
       double precision, public, parameter :: half  = 0.5d0,
     , two = 2.0d0
       double precision, public, parameter :: zero  = 0.0d0
       double precision, public, parameter :: three = 3.d0
       double precision, public, parameter :: sqrt2 =
     ,1.41421356237309504880168872420969807856967d0
       double precision, public, parameter :: msqrt2=
     ,-1.41421356237309504880168872420969807856967d0

       double precision, public, parameter :: sqrt3 =
     ,1.7320508075688772935274463415058723669428d0

       double precision, public, parameter :: propcut = 1D-10

       double complex, parameter, public :: czero = 0.0d0
       double complex, parameter, public :: cone = 1.0d0
       double complex, parameter, public :: ci=(0.0d0,1.0d0)
       double complex, parameter, public :: ne=(0.0d0,1.0d0)

c------masses of various particles
       double precision, parameter, public :: mt=1.75d0
       double precision, parameter, public :: mb=0d0
!        double precision, parameter, public :: mb=0.0d0
!        double precision, parameter, public :: mw=0.8d0


c------common data to be used to construct spinors

!-----dirac spinors

!       double complex, public :: u(8,16),uz(8,16),ux(8,16)
!
!         data u(1,1)/ one/, u(2,2)/ one/, u(3,5)/ one/,
!      ,u(4,6)/ one/, u(5,9)/ one/, u(6,10)/ one/,
!      ,u(7,13)/ one/, u(8,14)/ one/
!
!       data uz(1,1)/ one/, uz(1,3)/ one/,
!      ,     uz(2,2)/ one/, uz(2,4)/ mone/,
!      ,     uz(3,5)/ one/, uz(3,7)/ one/,
!      ,     uz(4,6)/ one/, uz(4,8)/ mone/,
!      ,     uz(5,9)/ one/, uz(5,11)/ one/,
!      ,     uz(6,10)/ one/, uz(6,12)/ mone/,
!      ,     uz(7,13)/ one/, uz(7,15)/ one/,
!      ,     uz(8,14)/ one/, uz(8,16)/ mone/
!
!       data ux(1,1)/ sqrt2/, ux(1,2)/ sqrt2/,
!      ,     ux(1,3)/ sqrt2/, ux(1,4)/ sqrt2/,
!      ,     ux(2,1)/ msqrt2/, ux(2,2)/ sqrt2/,
!      ,     ux(2,3)/ sqrt2/,  ux(2,4)/ msqrt2/,
!      ,     ux(3,5)/ sqrt2/,  ux(3,6)/ sqrt2/,
!      ,     ux(3,7)/ sqrt2/,  ux(3,8)/ sqrt2/,
!      ,     ux(4,5)/ msqrt2/, ux(4,6)/ sqrt2/,
!      ,     ux(4,7)/ sqrt2/,  ux(4,8)/ msqrt2/,
!      ,     ux(5,9)/ sqrt2/, ux(5,10)/ sqrt2/,
!      ,     ux(5,11)/ sqrt2/, ux(5,12)/ sqrt2/,
!      ,     ux(6,9)/ msqrt2/, ux(6,10)/ sqrt2/,
!      ,     ux(6,11)/ sqrt2/, ux(6,12)/ msqrt2/,
!      ,     ux(7,13)/ sqrt2/, ux(7,14)/ sqrt2/,
!      ,     ux(7,15)/ sqrt2/, ux(7,16)/ sqrt2/,
!      ,     ux(8,13)/ msqrt2/, ux(8,14)/ sqrt2/,
!      ,     ux(8,15)/ sqrt2/, ux(8,16)/ msqrt2/
!
!       double precision, public :: hc(16)
!
!              data hc(1)/ one/, hc(2)/ one/,
!      ,hc(3)/ mone/, hc(4)/ mone/, hc(5)/ one/, hc(6)/ one/,
!      ,hc(7)/ mone/, hc(8)/ mone/, hc(9)/ one/, hc(10)/ one/,
!      ,hc(11)/ mone/, hc(12)/ mone/, hc(13)/ one/, hc(14)/ one/,
!      ,hc(15)/ mone/, hc(16)/ mone/


!-----Weyl spinors

!       double complex, public :: wz(8,16),wx(8,16), bwz(8,16), bwx(8,16)
!
!
!       data wz(1,1)/ one/,
!      ,     wz(2,4)/ mone/,
!      ,     wz(3,5)/ one/,
!      ,     wz(4,8)/ mone/,
!      ,     wz(5,9)/ one/,
!      ,     wz(6,12)/ mone/,
!      ,     wz(7,13)/ one/,
!      ,     wz(8,16)/ mone/
!
!       data wx(1,1)/ one/, wx(1,2)/ one/,
!      ,     wx(2,3)/ one/,  wx(2,4)/ mone/,
!      ,     wx(3,5)/ one/,  wx(3,6)/ one/,
!      ,     wx(4,7)/ one/,  wx(4,8)/ mone/,
!      ,     wx(5,9)/ one/, wx(5,10)/ one/,
!      ,     wx(6,11)/ one/, wx(6,12)/ mone/,
!      ,     wx(7,13)/ one/, wx(7,14)/ one/,
!      ,     wx(8,15)/ one/, wx(8,16)/ mone/


!       data bwz(1,3)/ one/,
!      ,     bwz(2,2)/ mone/,
!      ,     bwz(3,7)/ one/,
!      ,     bwz(4,6)/ mone/,
!      ,     bwz(5,11)/ one/,
!      ,     bwz(6,10)/ mone/,
!      ,     bwz(7,15)/ one/,
!      ,     bwz(8,14)/ mone/
!
!       data bwx(1,3)/ one/, bwx(1,4)/ one/,
!      ,     bwx(2,1)/ one/,  bwx(2,2)/ mone/,
!      ,     bwx(3,7)/ one/,  bwx(3,8)/ one/,
!      ,     bwx(4,5)/ one/,  bwx(4,6)/ mone/,
!      ,     bwx(5,11)/ one/, bwx(5,12)/ one/,
!      ,     bwx(6,9)/ one/, bwx(6,10)/ mone/,
!      ,     bwx(7,15)/ one/, bwx(7,16)/ one/,
!      ,     bwx(8,13)/ one/, bwx(8,14)/ mone/




       end module constsd0
