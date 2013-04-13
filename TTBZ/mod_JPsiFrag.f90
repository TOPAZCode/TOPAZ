! this is the file with subroutine for topdecay amplitudes
      module ModJPsiFrag
      use ModParameters
      integer, parameter,private  :: dp = selected_real_kind(15)

      private

!----- notation for subroutines
      public ::  FF, psi, psi1, Dnp, fitFF !, D_np


!------------------------------------- for Mellin transform
real(dp), private, parameter :: theta = 0.8d0*pi/2.d0
real(dp), private, parameter :: mco = cos(theta)
real(dp), private, parameter :: msi = sin(theta)
real(dp), private, parameter :: r0 = 1.0_dp

!-------------------------------------- to fit fragmentation
integer, public, parameter :: Nfitpoint = 50  ! maximum is 2000!
real(dp), public :: A(1:1000), B(1:1000), C(1:1000)

real(dp), private, parameter :: tol = 0.00000010_dp
real(dp), private, parameter :: mb = 4.5_dp /100d0 !bottom quark mass

!------------------------------------- factorization scale for top decay
real(dp), private, parameter :: b0 = (33d0-2d0*5d0)/12d0/pi
real(dp), private, parameter :: mu0 = mb
real(dp), private :: muf_t,b1
real(dp), private  :: gs,as0 ! strong coupling at MZ
real(dp), private  :: as ! as at mt
real(dp), private  :: as_f
real(dp), private :: asmu0
real(dp), private, parameter :: CF = 4d0/3d0
!------------------------------------- alpha-parameter
! real(dp), private, parameter :: al_ff = 1_dp

real(dp), private, parameter :: epinv2 = 0.0_dp
real(dp), private, parameter :: epinv  = 0.0_dp

      contains

!------ various auxiliary functions


!--- this is the staff for the integration / to obtain fragment function



    function IntF(x,mu)
     implicit none
     real(8) :: IntF
     real(8):: y,x,mu, Lx,t, h
     integer :: i
     integer(8), parameter  :: Ni=300000

     IntF = 0d0

     if (x.eq.1d0) return
     if (x.eq.0d0) return

     h = one/real(Ni,dp)

       do i=1,Ni
          y = h*i
          IntF = IntF + h*F(y,x,mu)
       enddo

    end function IntF


     function F(y,x,mu)  ! y is the integration variable   ! changed
     implicit none
     complex(8) :: z, z1, z2
     real(8):: y,x,mu, Lx,t
     real(dp) :: F

       Lx = log(one/x)

       t =one/Lx/mco*log(one/y)

       z = dcmplx(r0-mco*t,t*msi)
       z1=dcmplx(r0,t*msi)
       z2 =dcmplx(-mco,msi)

       if (abs(z).lt.200d0) then
        F =two*real(z2/two/pi/(0d0,1d0)*exp(z1*Lx)*Dnp(z,mu)/Lx/mco,dp)
       else
        F = 0d0
       endif

     end function F



     subroutine fitFF(mu)
     implicit none
     real(8):: y,x,mu, Lx,t, h, f_3, f_2, f_1,porder
     integer :: i, i1

        muf_t = mu
        if(NLOParam.le.1) then
          porder=0d0
        else
          porder=1d0
        endif
        b1= porder*(153d0-19d0*5d0)/24d0/pi**2
        as0 = 0.1300_dp - porder*(0.130_dp - 0.118_dp)
        Lambda_QCD = 0.167d0/100d0  + porder*(0.226235d0 - 0.167d0 )/100d0
        as   = one/b0/log(muren**2/Lambda_QCD**2)*( one - b1/b0**2*log(log(muren**2/Lambda_QCD**2))/log(muren**2/Lambda_QCD**2))
        as_f = one/b0/log(muf_t**2/Lambda_QCD**2)*( one - b1/b0**2*log(log(muf_t**2/Lambda_QCD**2))/log(muf_t**2/Lambda_QCD**2))
        asmu0= one/b0/log(mu0**2/Lambda_QCD**2)  *( one - b1/b0**2*log(log(mu0**2/Lambda_QCD**2))/log(mu0**2/Lambda_QCD**2) )
        gs = sqrt(4.0_dp*pi*as)

!         print *, "# porder=",porder

       A=0d0
       B=0d0
       C=0d0
!        call GetFitValues(mu)


       if( all(A.eq.0d0) .and. all(B.eq.0d0) .and. all(C.eq.0d0) ) then
              continue
              print *, "# fitting fragmentation function at mu_f=",mu*100d0," GeV"
       else
              return
       endif



     h = one/real(Nfitpoint,8)

     do i1=1,Nfitpoint/2
          i = 2*i1-1
          x = h*i

          f_1 = IntF(x-h,mu)
          f_2 = IntF(x,mu)
          f_3 = IntF(x+h,mu)

          C(i1) = one/two/h**2*(f_3+f_1-two*f_2)
          B(i1) = one/h*(f_3-f_2)-two*x*C(i1)-C(i1)*h
          A(i1) = f_2 - B(i1)*x-C(i1)*x**2
!           print *, "# fit point ",i1,"/",Nfitpoint/2
     enddo

!      do i1=1,Nfitpoint/2
!         write(*,"(A2,I3,A2,1F24.16,A)") "A(",i1,")=",A(i1),"d0"
!      enddo
!      do i1=1,Nfitpoint/2
!         write(*,"(A2,I3,A2,1F24.16,A)") "B(",i1,")=",B(i1),"d0"
!      enddo
!      do i1=1,Nfitpoint/2
!         write(*,"(A2,I3,A2,1F24.16,A)") "C(",i1,")=",C(i1),"d0"
!      enddo


     end subroutine fitFF


!--- The fitted function

     function FFfit(x)
     implicit none
     real(8) :: FFfit
     real(8):: y,x, x1,x2,x3,mu, Lx,t, h
     integer :: i,i1

     FFfit = 0d0

     h = one/real(Nfitpoint,dp)

       do i1 = 1,Nfitpoint/2
       i = 2*i1-1
       x2 = h*i-h
       x3 = h*i+h
       if (x.le.x3.and.x.ge.x2) then
          FFfit = A(i1)+B(i1)*x + C(i1)*x**2
          return
       endif
       enddo
     end function FFfit



        function FF(x)
     implicit none
        real(dp) :: FF
        real(dp) :: x
        real(dp) :: Beta, a ,b

        FF = zero

        if (x.gt.1d0) return


        FF = FFfit(x)

!        FF = D_np(x) !+ as_f/2.0_dp/pi*4.0_dp/3.0_dp*( &
!        !(one+y**2)/(one-y)*log(muf_t**2/mb**2)-two*log(one-y) - one) &
!        !*(D_np(x/y)/y - D_np(x))

        end function FF




!------ lnGamma(z)

        function lnGamma(z) ! checked against Mathematica
     implicit none
        complex(8) :: z
        complex(dp) :: lnGamma
        real(dp) :: x
        real(dp), parameter :: gg = 7.0_dp
        real(dp) ::  C(1:9)
        complex(dp) ::  zk, z2, z1, zsq, rem, t
        integer(dp) :: k,i1

    C(1)=0.99999999999980993d0
    C(2)=676.5203681218851d0
    C(3)=-1259.1392167224028d0
    C(4)=771.32342877765313d0
    C(5)=-176.61502916214059d0
    C(6)=12.507343278686905d0
    C(7)=-0.13857109526572012d0
    C(8)=9.9843695780195716*1d-6
    C(9)=1.5056327351493116*1d-7

        x = real(z,dp)
        rem = (0.0_dp,0.0_dp)

         if (x.lt.0.5_dp) then
            z1 = one -z
            z2 = -z

           rem = C(1)

           do i1=2,9
              rem = rem + C(i1)/(z2+real(i1,dp)-one)
           enddo

           t = z2 + gg + 0.5_dp

   lnGamma = log(pi/sin(pi*z)/(sqrt(two*pi)*rem*t**(z2+0.5_dp)*exp(-t)))

           else
           z2 = z-one
           rem = C(1)

           do i1=2,9
              rem = rem + C(i1)/(z2+real(i1,dp)-one)
           enddo

           t = z2 + gg + 0.5_dp


           lnGamma = log(sqrt(two*pi)*rem*t**(z2+0.5_dp)*exp(-t))


         endif


        end function lnGamma




!----- digamma function Psi(z) ! checked against Mathematica

        function psi(z)
     implicit none
        complex(8) :: z
        complex(dp) :: psi
        real(dp) :: x
        real(dp), parameter :: x0 = 7.0_dp
        real(dp) ::  B(1:7)
        complex(dp) ::  zk, z2, z1, zsq, rem
        integer(dp) :: k,i1



        B(1)=1d0/6d0
        B(2)=-1d0/30d0
        B(3)=1d0/42d0
        B(4)=-1d0/30d0
        B(5)=5d0/66d0
        B(6)=-691d0/2730d0
        B(7)=7d0/6d0


        x = real(z,dp)
        rem = (0.0_dp,0.0_dp)

         if (x.lt.0.0_dp) then
            z1 = -z
            rem = -one/z - pi/tan(pi*z)
         else
            z1 = z
         endif


         x = real(z1,dp)

         if (x.lt.0.0_dp) then
            print *, 'something is wrong'
            pause
         elseif(x.gt.0.0_dp.and.x.lt.x0) then
            z2 = z1+8.0_dp
            do i1 = 0,7
            rem = rem - one/(z1+i1)
            enddo
         else
            z2 = z1
         endif

           psi = log(z2) - one/two/z2  + rem
           zk = (one,0.0_dp)
           zsq = z2**2

           do k=1,7
              zk = zk*zsq
              psi = psi-B(k)/two/real(k,dp)/zk
           enddo



        end function psi



!-------- trigamma function

        function psi1(z) ! checked against Mathematica
     implicit none
        complex(8) :: z
        complex(dp) :: psi1
        real(dp) :: x
        real(dp), parameter :: x0 = 7.0_dp
        real(dp) ::  B(1:7)
        complex(dp) ::  zk, z2, z1, zsq, rem
        integer(dp) :: k,i1
        real(dp) :: puh


        puh = 1.0_dp


        B(1)=1d0/6d0
        B(2)=-1d0/30d0
        B(3)=1d0/42d0
        B(4)=-1d0/30d0
        B(5)=5d0/66d0
        B(6)=-691d0/2730d0
        B(7)=7d0/6d0

        x = real(z,dp)
        rem = (0.0_dp,0.0_dp)

         if (x.lt.0.0_dp) then
            z1 = -z
            if (abs(sin(pi*z)).gt.1d16) then
            rem = one/z**2
            else
            rem = one/z**2  +  pi**2/sin(pi*z)**2
            endif
            puh = -1.0_dp
         else
            z1 = z
         endif


         x = real(z1,dp)

         if (x.lt.0.0_dp) then
            print *, 'something is wrong'
            pause
         elseif(x.gt.0.0_dp.and.x.lt.x0) then
            z2 = z1+8.0_dp
            do i1 = 0,7
            rem = rem + puh*(one/(z1+i1)**2)
            enddo
         else
            z2 = z1
         endif

           psi1 = puh*(one/z2 + one/two/z2**2)  + rem
           zk = (one,0.0_dp)
           zsq = z2**2


           do k=1,7
              zk = zk*zsq
              psi1 = psi1+puh*B(k)/zk/z2
           enddo


        end function psi1





!-------- four-gamma function

        function psi2(z) ! checked against Mathematica
     implicit none
        complex(8) :: z
        complex(dp) :: psi2
        real(dp) :: x
        real(dp), parameter :: x0 = 7.0_dp
        real(dp) ::  B(1:7)
        complex(dp) ::  zk, z2, z1, zsq, rem
        integer(dp) :: k,i1
        real(dp) :: puh


        puh = 1.0_dp


        B(1)=1d0/6d0
        B(2)=-1d0/30d0
        B(3)=1d0/42d0
        B(4)=-1d0/30d0
        B(5)=5d0/66d0
        B(6)=-691d0/2730d0
        B(7)=7d0/6d0

        x = real(z,dp)
        rem = (0.0_dp,0.0_dp)

         if (x.lt.0.0_dp) then
            z1 = -z
            if (abs(sin(pi*z)).gt.1d16) then
            rem = -two/z**3
            else
            rem = -two/z**3 - two*pi**3*cos(pi*z)/sin(pi*z)**3
            endif
         else
            z1 = z
         endif


         x = real(z1,dp)

         if (x.lt.0.0_dp) then
            print *, 'something is wrong'
            pause
         elseif(x.gt.0.0_dp.and.x.lt.x0) then
            z2 = z1+8.0_dp
            do i1 = 0,7
            rem = rem - puh*two/(z1+i1)**3
            enddo
         else
            z2 = z1
         endif


           psi2 = (-1d0)*(one/z2**2 + one/z2**3)  + rem
           zk = (one,0.0_dp)
           zsq = z2**2


           do k=1,7
              zk = zk*zsq
              psi2 = psi2+(-1d0)*puh*B(k)*(two*k+1)/zk/z2**2
           enddo



        end function psi2





!       non-perturbative fragmentation function
!       / power like and  Kartelishvili
!       all Mellin space
        function Dnp_i(z)      !  / initial np function
        implicit none
        complex(8) :: Dnp_i, z, tmp

        if (fragm_func_type.eq.1) then ! double power
        tmp =   lnGamma(beta_frag+z)    &
              - lnGamma(dcmplx(beta_frag+one,0.0_dp))  &
              + lnGamma(dcmplx(beta_frag+alpha_frag+two,0_dp)) &
              - lnGamma(beta_frag+alpha_frag+z+one)

        Dnp_i = exp(tmp)
        elseif(fragm_func_type.eq.2) then ! karteshvili
        Dnp_i = (one+delta_frag)*(two+delta_frag)/ &
              (z+delta_frag)/(z+delta_frag+one)
        endif

        end function Dnp_i






!------------------------------------------------------------------
        function Dbi(z)  ! initial condition in the mellin space
        implicit none
        complex(8) :: Dbi,z,z1,S1N,S2N, z3
        real(8), parameter :: gammaE = -0.5772156649_dp

        Dbi = one
        if (NLOParam.LE.1) return

        z1 = z+one
        S1N = psi(z1)-gammaE

!------------ Mellin of (1+x^2)/(1-x)_+, essentially
        Dbi = one + asmu0/two/pi*CF*(                      &
        (log(mu0**2/mb**2)-one)*(three/two + one/z/z1 - two*S1N))


!------------subtract  Mellin of 4*log(-x)/(1-x)_+
        Dbi = Dbi - asmu0/two/pi*CF* &
         two*((psi(z)-gammaE)**2 -psi1(z)+pisq/6d0)


!----------- add the Mellin of two*(1+x)*log(1-x)

        Dbi = Dbi + asmu0/two/pi*CF*two*( &
        -(two*z+one)/(z+one)/z*psi(z)  &
        +(two*z+one)/(z+one)/z*gammaE        &
        -(three*z+one+three*z**2)/z**2/(z+one)**2 + 7d0/4d0)

        end function Dbi






        function Dnp(z,mu)   ! solution of the RG equation
        implicit none
        complex(8) :: Dnp,z
        real(8) :: mu,b1
        real(8) :: asmu, t
        complex(8) :: P0N, PFN, DelN, PGN, PNFN, P1N
        complex(8) :: z1, S1N,S2N,S3N
        complex(8) :: phase
        real(8), parameter :: gammaE = -0.5772156649_dp
        real(8), parameter :: zeta3 = 1.202056903_dp
        real(8), parameter :: CA=3d0,CF=4d0/3d0

        if( NLOParam.LE.1 ) then
            b1=0d0
        else
            b1= (153d0-19d0*5d0)/24d0/pi**2
        endif


        asmu = one/b0/log(mu**2/Lambda_QCD**2)*( one &
  - b1/b0**2*log(log(mu**2/Lambda_QCD**2))/log(mu**2/Lambda_QCD**2)  )

        t = one/two/pi/b0*log(asmu0/asmu)

        z1 = z + one
        S1N = psi(z1)-gammaE
        S2N = -psi1(z1)+pi**2/6.0_dp
        S3N = one/two*(psi2(z1)+two*zeta3)

!--------- terms for the LL evolution
        P0N = CF*(three/two + one/z/z1 - two*S1N)

!--------- terms for the NLL evolution

        PFN = (two*S1N-one/z/z1)*(two*S2N-two*pi**2/6.0_dp) &
           -two*(two*z+one)/z**2/z1**2*S1N &
           +4.0_dp*S3N - 3.0_dp*S2N + 3.0_dp*pi**2/6.0_dp &
           +(3.0_dp*z**3+z**2-one)/z**3/z1**3 - 23.0_dp/8.0_dp


        PNFN = 20.0_dp/9.0_dp*S1N - 4.0_dp/3.0_dp*S2N - one/6.0_dp  &
         - two*(11.0_dp*z**2 + 5.0_dp*z - three)/(9.0_dp*z**2*z1**2)

        PGN = -PFN + S1N*(-134.0_dp/9.0_dp-two*(two*z+one)/z**2/z1**2) &
         + 4.0_dp*S1N*S2N+S2N*(13.0_dp/3.0_dp-two/z/z1)  &
         +43.0_dp/24.0_dp + (151.0_dp*z**4   &
         + 263.0_dp*z**3 + 97.0_dp*z**2  &
         + 3.0_dp*z+9.0_dp)/9.0_dp/z**3/z1**3

        DelN = two*(-two*S1N+three/two+one/z/z1)*(two*S2N - pi**2/3.0_dp &
        - (two*z+one)/z**2/z1**2)

        P1N = CF**2*(PFN+DelN)+one/two*CF*CA*PGN + one/two*CF*5d0*PNFN


         if ( NLOParam.LE.1 ) then
            phase = P0N*t
         else
            phase = P0N*t +one/four/pi**2/b0*(asmu0-asmu)*(P1N - two*pi*b1/b0*P0N)
         endif

      Dnp = Dnp_i(z)*Dbi(z)*exp(phase)

      end function Dnp


!------ non-perturbative fragmentation function in x-space; input

!         function D_np(x)
!         implicit none
!         real(dp) :: D_np
!         real(dp) :: x
!         real(dp) :: Beta, a ,b, tmp
!
!         D_np = zero
!
!         if (x.gt.1d0) return
!
!
!         if (fragm_func_type.eq.1) then ! double power
!
! !        Beta = 0.01167667555_dp
!         tmp = real(  &
!           lnGamma(dcmplx(alpha_frag+one)) &
!         + lnGamma(dcmplx(beta_frag+one))  &
!         - lnGamma(dcmplx(alpha_frag++beta_frag + two)),dp)
!         Beta = exp(tmp)
!     D_np = one/Beta*dexp(alpha_frag*dlog(dabs(one-x)) &
!         + beta_frag*dlog(dabs(x)))
!
!         elseif(fragm_func_type.eq.2) then
!
!     D_np = one*(one+delta_frag)*(two+delta_frag)* &
!        (one-x)*exp(delta_frag*log(x))
!
!         endif
!
!
!         end function D_np







SUBROUTINE GetFitValues(mu)
implicit none
real(8) :: mu


! if( Fragm_Func_Type.ne.2 .or. mb.ne.4.5d0/100d0 .or. NFitPoint.ne.50 .or. &
!     alpha_frag.ne.0.66d0 .or. beta_frag.ne.12.39d0 .or. delta_frag.ne.14.97d0 &
!   ) then
! !       print *, Fragm_Func_Type,mb,NFitPoint,alpha_frag,beta_frag,delta_frag
!       return
! endif

    print *, "fit values for mt=173, porder=nlo"
    A(  1)=      0.0000000000000001d0
    A(  2)=      0.3406426445705983d0
    A(  3)=      0.3666618865143286d0
    A(  4)=      0.3801760607699525d0
    A(  5)=      0.3906925748794750d0
    A(  6)=      0.4022738417364065d0
    A(  7)=      0.4184002853419221d0
    A(  8)=      0.4434441690721482d0
    A(  9)=      0.4838141849980891d0
    A( 10)=      0.5491636154050319d0
    A( 11)=      0.6531129356673160d0
    A( 12)=      0.8113120035223162d0
    A( 13)=      1.0314690864532368d0
    A( 14)=      1.2846654501193930d0
    A( 15)=      1.4404691793271858d0
    A( 16)=      1.1438652997260550d0
    A( 17)=     -0.3791450457285965d0
    A( 18)=     -4.5525166412548570d0
    A( 19)=    -13.5416917331326978d0
    A( 20)=    -29.6858501050806893d0
    A( 21)=    -52.9849303830990479d0
    A( 22)=    -73.7022301051200941d0
    A( 23)=    -53.4544326304085970d0
    A( 24)=    121.5827655144488091d0
    A( 25)=    676.6910858484186519d0
    B(  1)=     26.0176453138009194d0
    B(  2)=      1.4963820257374110d0
    B(  3)=      0.8056260953254859d0
    B(  4)=      0.5751067710509010d0
    B(  5)=      0.4424760004283435d0
    B(  6)=      0.3264249799199037d0
    B(  7)=      0.1921326790785383d0
    B(  8)=      0.0135344678546285d0
    B(  9)=     -0.2383201989906134d0
    B( 10)=     -0.6007327931097486d0
    B( 11)=     -1.1196735803784326d0
    B( 12)=     -1.8379189906479527d0
    B( 13)=     -2.7547635283422016d0
    B( 14)=     -3.7293968351073072d0
    B( 15)=     -4.2896637153437993d0
    B( 16)=     -3.3107091980300218d0
    B( 17)=      1.4290728819281564d0
    B( 18)=     13.6699882674705524d0
    B( 19)=     38.5909251598450282d0
    B( 20)=     81.0208740706118249d0
    B( 21)=    139.2458082828557622d0
    B( 22)=    188.6832033907831772d0
    B( 23)=    143.1490118767898707d0
    B( 24)=   -235.9379589547010596d0
    B( 25)=  -1390.3053723630932836d0
    C(  1)=   -404.6182646823959317d0
    C(  2)=     -4.4883353374314261d0
    C(  3)=      0.0806072390097806d0
    C(  4)=      1.0631172846574388d0
    C(  5)=      1.4812582686452025d0
    C(  6)=      1.7719816997641136d0
    C(  7)=      2.0515599740081325d0
    C(  8)=      2.3699734359018674d0
    C(  9)=      2.7627808330164818d0
    C( 10)=      3.2652430511950770d0
    C( 11)=      3.9129117677275116d0
    C( 12)=      4.7281437546262133d0
    C( 13)=      5.6826936468237772d0
    C( 14)=      6.6206113829617159d0
    C( 15)=      7.1242648387162539d0
    C( 16)=      6.3165736420853218d0
    C( 17)=      2.6289511183583136d0
    C( 18)=     -6.3469374340693889d0
    C( 19)=    -23.6191200548868672d0
    C( 20)=    -51.4975875068168776d0
    C( 21)=    -87.8739423377181197d0
    C( 22)=   -117.3667316923637856d0
    C( 23)=    -91.7697649340959458d0
    C( 24)=    113.4793076322015963d0
    C( 25)=    713.6142865146745180d0




if(     mu.eq.172d0/100d0 .and. NLOParam.le.1 ) then
!     A(  1)=      0.0000000000000000d0
!     A(  2)=      0.5203747085569792d0
!     A(  3)=      0.4894932856381513d0
!     A(  4)=      0.4718179688628537d0
!     A(  5)=      0.4607547319810033d0
!     A(  6)=      0.4548936672124516d0
!     A(  7)=      0.4546807422132016d0
!     A(  8)=      0.4620561646426216d0
!     A(  9)=      0.4809828757127768d0
!     A( 10)=      0.5187850294018997d0
!     A( 11)=      0.5885470455475232d0
!     A( 12)=      0.7126382952418419d0
!     A( 13)=      0.9263369929119925d0
!     A( 14)=      1.2775792801420467d0
!     A( 15)=      1.8126454668558134d0
!     A( 16)=      2.5266273893807099d0
!     A( 17)=      3.2415221361869913d0
!     A( 18)=      3.3581030722789276d0
!     A( 19)=      1.4276285539057920d0
!     A( 20)=     -5.4428464528270180d0
!     A( 21)=    -22.0594324632942893d0
!     A( 22)=    -54.2772788827742616d0
!     A( 23)=   -102.5162781808016064d0
!     A( 24)=   -131.5626114161593421d0
!     A( 25)=    178.2500728961877599d0
!     B(  1)=     38.6188923143889440d0
!     B(  2)=     -1.4119482341718959d0
!     B(  3)=     -0.5970212246312889d0
!     B(  4)=     -0.2952807765261217d0
!     B(  5)=     -0.1548106250053375d0
!     B(  6)=     -0.0952530937166352d0
!     B(  7)=     -0.0929203544066087d0
!     B(  8)=     -0.1451536373113682d0
!     B(  9)=     -0.2629713963767017d0
!     B( 10)=     -0.4723778838015700d0
!     B( 11)=     -0.8203397877864338d0
!     B( 12)=     -1.3831815563675689d0
!     B( 13)=     -2.2719388394360305d0
!     B( 14)=     -3.6208400763372817d0
!     B( 15)=     -5.5298398551630168d0
!     B( 16)=     -7.9091635071765580d0
!     B( 17)=    -10.1465913808837129d0
!     B( 18)=    -10.5016625265856085d0
!     B( 19)=     -5.1676893842294351d0
!     B( 20)=     12.8583799003890054d0
!     B( 21)=     54.3134046145481335d0
!     B( 22)=    130.9162816566980609d0
!     B( 23)=    240.5091468939564265d0
!     B( 24)=    304.0545459766992735d0
!     B( 25)=   -336.7728458935459344d0
!     C(  1)=   -666.7896232646955923d0
!     C(  2)=      8.7471976012132995d0
!     C(  3)=      3.3858323130225854d0
!     C(  4)=      2.0987811326529657d0
!     C(  5)=      1.6530003763453456d0
!     C(  6)=      1.5017393391156264d0
!     C(  7)=      1.4957162065608320d0
!     C(  8)=      1.5881893798250224d0
!     C(  9)=      1.7715387141097050d0
!     C( 10)=      2.0615401167627114d0
!     C( 11)=      2.4954322758147240d0
!     C( 12)=      3.1336509146648317d0
!     C( 13)=      4.0577168790863203d0
!     C( 14)=      5.3527847635485850d0
!     C( 15)=      7.0555069874919667d0
!     C( 16)=      9.0377632893898223d0
!     C( 17)=     10.7883958391122299d0
!     C( 18)=     11.0584378870909319d0
!     C( 19)=      7.3740510350012167d0
!     C( 20)=     -4.4495944801947740d0
!     C( 21)=    -30.3049597315385739d0
!     C( 22)=    -75.8386477926581364d0
!     C( 23)=   -138.0838984522139583d0
!     C( 24)=   -172.8374829521733886d0
!     C( 25)=    158.5227729973581745d0


elseif(     mu.eq.172d0/100d0 .and. NLOParam.gt.1 ) then

!     A(  1)=      0.0000000000000000d0
!     A(  2)=      0.3655316605239591d0
!     A(  3)=      0.3927112310645687d0
!     A(  4)=      0.4064646912969148d0
!     A(  5)=      0.4169081418041090d0
!     A(  6)=      0.4282713858470729d0
!     A(  7)=      0.4440924977290778d0
!     A(  8)=      0.4686421827707030d0
!     A(  9)=      0.5076502147522954d0
!     A( 10)=      0.5683441005043444d0
!     A( 11)=      0.6577436857972475d0
!     A( 12)=      0.7769125091827811d0
!     A( 13)=      0.9072024445543958d0
!     A( 14)=      0.9826561412215356d0
!     A( 15)=      0.8414835127472120d0
!     A( 16)=      0.1509229986727800d0
!     A( 17)=     -1.6917726019354578d0
!     A( 18)=     -5.6546613913644936d0
!     A( 19)=    -13.0733116932768958d0
!     A( 20)=    -25.3219628402052486d0
!     A( 21)=    -42.5443167024170563d0
!     A( 22)=    -59.7854472791486202d0
!     A( 23)=    -55.1138336279069136d0
!     A( 24)=     58.4278835036341988d0
!     A( 25)=    976.1274783971858824d0
!     B(  1)=     27.8788041083089198d0
!     B(  2)=      1.5554597267661783d0
!     B(  3)=      0.8333750732569714d0
!     B(  4)=      0.5986457695596192d0
!     B(  5)=      0.4668959903750797d0
!     B(  6)=      0.3530234911967134d0
!     B(  7)=      0.2212767575374974d0
!     B(  8)=      0.0461951634954327d0
!     B(  9)=     -0.1972080495716209d0
!     B( 10)=     -0.5339293254640829d0
!     B( 11)=     -0.9805089274985606d0
!     B( 12)=     -1.5220718600239547d0
!     B( 13)=     -2.0656229361051714d0
!     B( 14)=     -2.3581152038068658d0
!     B( 15)=     -1.8590710638494317d0
!     B( 16)=      0.4330921259051599d0
!     B( 17)=      6.1752565461297539d0
!     B( 18)=     17.8063149652527173d0
!     B( 19)=     38.3810997767582904d0
!     B( 20)=     70.5792936502810022d0
!     B( 21)=    113.6168501093207084d0
!     B( 22)=    154.7199993708096031d0
!     B( 23)=    144.3832234063433759d0
!     B( 24)=   -101.2867914004572754d0
!     B( 25)=  -2001.9566567761560236d0
!     C(  1)=   -434.2704177612354783d0
!     C(  2)=     -4.6440960501416271d0
!     C(  3)=      0.1351542217532115d0
!     C(  4)=      1.1361303475404538d0
!     C(  5)=      1.5516191820065528d0
!     C(  6)=      1.8369005768242852d0
!     C(  7)=      2.1111732191203547d0
!     C(  8)=      2.4233308886376093d0
!     C(  9)=      2.8030281171519125d0
!     C( 10)=      3.2700479749392386d0
!     C( 11)=      3.8277495719447883d0
!     C( 12)=      4.4430329754863163d0
!     C( 13)=      5.0099365394939976d0
!     C( 14)=      5.2933769334945202d0
!     C( 15)=      4.8523945039605216d0
!     C( 16)=      2.9503461712429568d0
!     C( 17)=     -1.5230171791855152d0
!     C( 18)=    -10.0572536315526868d0
!     C( 19)=    -24.3226756268682962d0
!     C( 20)=    -45.4826067157687604d0
!     C( 21)=    -72.3696243798624437d0
!     C( 22)=    -96.8672360620040394d0
!     C( 23)=    -91.1534586892140339d0
!     C( 24)=     41.7325248773212749d0
!     C( 25)=   1025.8291783789702549d0



elseif( mu.eq.4.5d0/100d0 .and. NLOParam.le.1 ) then
!     print *, "reading fragmentation function fit at mu_f=",4.5d0," GeV"
!     A(  1)=      0.0000000000000000d0
!     A(  2)=      0.0000446494856518d0
!     A(  3)=      0.0000286472658737d0
!     A(  4)=      0.0000214481082062d0
!     A(  5)=      0.0000172700081663d0
!     A(  6)=      0.0000169939490540d0
!     A(  7)=      0.0000436670241398d0
!     A(  8)=      0.0002587071830038d0
!     A(  9)=      0.0014979373150616d0
!     A( 10)=      0.0072316795234153d0
!     A( 11)=      0.0295221076514480d0
!     A( 12)=      0.1045005598652120d0
!     A( 13)=      0.3272759208362784d0
!     A( 14)=      0.9204549312962717d0
!     A( 15)=      2.3490975197504045d0
!     A( 16)=      5.4743629689000590d0
!     A( 17)=     11.6699448195029429d0
!     A( 18)=     22.6561044554412732d0
!     A( 19)=     39.4417406971653577d0
!     A( 20)=     59.0494280214439442d0
!     A( 21)=     66.3833183185712130d0
!     A( 22)=     16.3146871930257724d0
!     A( 23)=   -202.8515795733746074d0
!     A( 24)=   -844.3618566081320296d0
!     A( 25)=  -2403.0284183672729341d0
!     B(  1)=      0.0033969902306001d0
!     B(  2)=     -0.0006961611582260d0
!     B(  3)=     -0.0002695229798316d0
!     B(  4)=     -0.0001462963215707d0
!     B(  5)=     -0.0000932086526715d0
!     B(  6)=     -0.0000895189126224d0
!     B(  7)=     -0.0003074266384521d0
!     B(  8)=     -0.0018227823468191d0
!     B(  9)=     -0.0094881173514899d0
!     B( 10)=     -0.0410824699119433d0
!     B( 11)=     -0.1518019620976921d0
!     B( 12)=     -0.4907821470748135d0
!     B( 13)=     -1.4148985442654887d0
!     B( 14)=     -3.6879725627255149d0
!     B( 15)=     -8.7747410932256233d0
!     B( 16)=    -19.1663258158283902d0
!     B( 17)=    -38.4891812879164519d0
!     B( 18)=    -70.7546696085575917d0
!     B( 19)=   -117.3457776452267325d0
!     B( 20)=   -168.9739746344259856d0
!     B( 21)=   -187.5234891846965581d0
!     B( 22)=    -68.9584162530029516d0
!     B( 23)=    427.6074381119585723d0
!     B( 24)=   1818.7686105335919819d0
!     B( 25)=   5059.3248081988913327d0
!     C(  1)=     -0.0707122554556846d0
!     C(  2)=      0.0037106007324919d0
!     C(  3)=      0.0008779703428966d0
!     C(  4)=      0.0003510230287455d0
!     C(  5)=      0.0001824321309336d0
!     C(  6)=      0.0001708849084955d0
!     C(  7)=      0.0006157595458826d0
!     C(  8)=      0.0032848850494441d0
!     C(  9)=      0.0151372001806631d0
!     C( 10)=      0.0586574525989469d0
!     C( 11)=      0.1961410072631143d0
!     C( 12)=      0.5792651249086174d0
!     C( 13)=      1.5376006703967142d0
!     C( 14)=      3.7151875015328839d0
!     C( 15)=      8.2430822356716842d0
!     C( 16)=     16.8810971923717013d0
!     C( 17)=     31.9471266150608493d0
!     C( 18)=     55.6373783929536856d0
!     C( 19)=     87.9675511257424461d0
!     C( 20)=    121.9525621844558287d0
!     C( 21)=    133.6802517830326735d0
!     C( 22)=     63.4902997747033311d0
!     C( 23)=   -217.7747894230830354d0
!     C( 24)=   -971.9784775057216848d0
!     C( 25)=  -2656.2963898316183986d0

endif


return
END SUBROUTINE



       end module ModJPsiFrag


