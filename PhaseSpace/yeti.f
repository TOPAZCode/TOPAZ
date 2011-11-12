
       subroutine Yeti3(Ecm,xRnd,Mass,Mom,Wgt)
c---------------------------------------------------------------------------------------------
c                           Yeti  updated version
c              3-particle phase space generator
c
c   integration  parametrization:  E_1     --> x2            (energy of particle 1)
c                                  E_3     --> x1            (energy of particle 3)
c                                  theta_3 --> x3            (polar angle of particle 3)
c                                  phi_3   --> x5            (azimuth angle of particle 3)
c                                  eta     --> x4            (azimuth angle around particle 3)
c
c           INPUT:
c            Ecm           partonic center of mass energy
c            xRnd          integration variables, range[0..1]
c            Mass(i)       mass of particle i, i=1,2,3
c
c           OUTPUT:
c            Wgt           integration weight
c            Mom(0:3,i)    4-momentum vector (P_mu)=(E,px,py,pz) for particle i=1,2,3
c
c
c     parameter 'check' allows for checking on-shell relations and momentum conservation
c---------------------------------------------------------------------------------------------
       implicit none
       double precision Ecm,Wgt,Norm
       double precision xRnd(1:5),Mass(1:3),Mom(0:3,1:3)
       double precision E3,E1
       double precision phi3,sinP3,cosP3,cosT3,Eta
       double precision sinT3,cosEta,sinEta,cosXi,sinXi
       double precision Jac1,Jac2
       double precision E1max,E1min,E3max,E3min
       double precision k1abs,k3abs
       double precision sig,tau,mp,mm
       double precision dummy1,dummy2
       double precision OS1,OS2,OS3,EC,MC1,MC2,MC3
       double precision twoPi
       parameter ( twoPi=6.2831853071795864d0 )
       logical check
       parameter (check=.false.)

c      1st integration variable: E_3 with E3min<=E3<=E3max
       E3min = Mass(3)
       E3max = (Ecm -((Mass(1)+Mass(2))**2-Mass(3)**2)/Ecm)/2d0
       E3    = (E3max-E3min)*xRnd(1) + E3min
       Jac1  = (E3max-E3min)



c      2nd integration variable: E_1 with E1min<=E1<=E1max
       k3abs= dsqrt(E3**2 - Mass(3)**2)
       sig  = Ecm - E3
       tau  = sig**2 - k3abs**2
       mp   = Mass(1) + Mass(2)
       mm   = Mass(1) - Mass(2)
       dummy1 = 0.5d0/tau * sig*(tau+mp*mm)
       dummy2 = 0.5d0/tau * k3abs*dsqrt((tau-mp**2)*(tau-mm**2))
       E1min  = dummy1 - dummy2
       E1max  = dummy1 + dummy2
       Jac2   = (E1max-E1min)
       E1     = Jac2*xRnd(2) + E1min


c      3rd integration variable: cos(theta_3) with -1<=cos(theta_3)<=+1
       cosT3 = 2d0*xRnd(3) - 1d0
       sinT3 = dsqrt(1d0-cosT3**2)

c      4th integration variable: eta with 0<=eta<=2*Pi
       eta  = twoPi * xRnd(4)

c      5th integration variable: eta with 0<=phi_3<=2*Pi
       Phi3  = twoPi * xRnd(5)
       sinP3 = dsin(Phi3)
       cosP3 = dcos(Phi3)

c      set p3
       Mom(0,3) = E3
       Mom(1,3) = k3abs * sinT3 * cosP3
       Mom(2,3) = k3abs * sinT3 * sinP3
       Mom(3,3) = k3abs * cosT3

c      set p1
       k1abs  = dsqrt(E1**2 - Mass(1)**2)
       cosEta = dcos(eta)
!        sinEta = dsin(eta)
       sinEta = dsqrt(1d0-cosEta**2)

c      usual formula for cosXi:
!        cosXi  = ( Mass(1)**2-Mass(2)**2+Mass(3)**2
!      #             +2d0*E1*E3+Ecm**2-2d0*Ecm*(E1+E3)
!      #            ) / 2d0 / k1abs / k3abs
c      more accurate evaluation for cosXi:
       cosXi  = ( Mass(1)**2-Mass(2)**2+Mass(3)**2
     #           +2d0*Ecm*(Ecm/2d0-E1-E3) )/(2d0*k1abs*k3abs)
     # + dsqrt((1d0+(Mass(1)/k1abs)**2)*(1d0+(Mass(3)/k3abs)**2))

c      safety exit for |cosXi|=1d0 + 1d-16
       if( cosXi .gt.  1d0 ) cosXi= 1d0
       if( cosXi .lt. -1d0 ) cosXi=-1d0
       sinXi  = dsqrt(1d0-cosXi**2)

       Mom(0,1) = E1
       Mom(1,1) = k1abs * (cosP3*(cosT3*cosEta*sinXi+sinT3*cosXi)
     #                             -sinP3*sinXi*sinEta)
       Mom(2,1) = k1abs * (cosP3*(sinEta*sinXi)+sinP3*(cosXi*sinT3
     #                             +cosEta*cosT3*sinXi))
       Mom(3,1) = k1abs * (-sinT3*cosEta*sinXi+cosT3*cosXi )

c      set p2
       Mom(0,2) = Ecm - E1 - E3
       Mom(1,2) = - Mom(1,1) - Mom(1,3)
       Mom(2,2) = - Mom(2,1) - Mom(2,3)
       Mom(3,2) = - Mom(3,1) - Mom(3,3)

c      set integration weight
       Norm = 0.25d0 * (twoPi)**(-3)
       Wgt  = Norm * Jac1*Jac2


c      check on-shell relations and energy-momentum conservation
       if( check ) then
        OS1=(Mass(1)**2
     #     -(Mom(0,1)**2-Mom(1,1)**2-Mom(2,1)**2-Mom(3,1)**2 )
     #      )/Ecm**2
        OS2=(Mass(2)**2
     #     -(Mom(0,2)**2-Mom(1,2)**2-Mom(2,2)**2-Mom(3,2)**2 )
     #      )/Ecm**2
        OS3=(Mass(3)**2
     #     -(Mom(0,3)**2-Mom(1,3)**2-Mom(2,3)**2-Mom(3,3)**2 )
     #      )/Ecm**2
        EC =(Ecm-Mom(0,1)-Mom(0,2)-Mom(0,3))/Ecm
        MC1=(Mom(1,1)+Mom(1,2)+Mom(1,3))/Ecm
        MC2=(Mom(2,1)+Mom(2,2)+Mom(2,3))/Ecm
        MC3=(Mom(3,1)+Mom(3,2)+Mom(3,3))/Ecm
        if( dabs(OS1).gt.1d-15 ) print *,'OS1=',OS1
        if( dabs(OS2).gt.1d-15 ) print *,'OS2=',OS2
        if( dabs(OS3).gt.1d-15 ) print *,'OS3=',OS3
        if( dabs( EC).gt.1d-15 ) print *,' EC=', EC
        if( dabs(MC1).gt.1d-15 ) print *,'MC1=',MC1
        if( dabs(MC2).gt.1d-15 ) print *,'MC3=',MC2
        if( dabs(MC3).gt.1d-15 ) print *,'MC3=',MC3
       endif

       return
       end







!        subroutine Yeti3PPS(Ecm,x1,x2,x3,x4,Wgt)
! c---------------------------------------------------------------------------------------------
! c                           Yeti   older version
! c              3-particle phase space generator
! c
! c   integration  parametrization:  E_1     --> x2            (energy of particle 1)
! c                                  E_3     --> x1            (energy of particle 3)
! c                                  theta_3 --> x3            (polar angle of particle 3)
! c                                  eta     --> x4            (azimuth angle around particle 3)
! c
! c           INPUT:
! c            Ecm           partonic center of mass energy
! c          yetiM(i)        mass for particle i, i=1,2,3
! c         x1,x2,x3,x4      integration variables, range[0..1]
! c
! c           OUTPUT:
! c            Wgt           integration weight
! c        yetiP(i,mu)       4-momentum (P_mu)=(P_0,P_1,P_2,P_3)=(E,px,py,pz) for particle i
! c
! c     'commonB.yeti.f' includes:
! c     double precision yetiM(1:3),yetiP(1:3,0:3)
! c     common/yeti/ yetiM,yetiP
! c
! c     parameter 'check' allows for checking on-shell relations and momentum conservation
! c---------------------------------------------------------------------------------------------
!        implicit none
!        include 'commonB.yeti.f'
!        double precision Ecm,Wgt,Norm
!        double precision x1,x2,x3,x4
!        double precision E3,E1
!        double precision cosT3,Eta,sinT3,cosEta,sinEta,cosXi,sinXi
!        double precision Jac1,Jac2,Jac3,Jac4
!        double precision E1max,E1min,E3max,E3min
!        double precision k1abs,k3abs
!        double precision sig,tau,mp,mm
!        double precision dummy1,dummy2
!        double precision OS1,OS2,OS3,EC,MC1,MC2,MC3
!        double precision dblPi
!        parameter ( dblPi=3.1415926535897932d0 )
!        logical check
!        parameter (check=.false.)
!
! c      1st integration variable: E_3 with E3min<=E3<=E3max
!        E3min = yetiM(3)
!        E3max = (Ecm -((yetiM(1)+yetiM(2))**2-yetiM(3)**2)/Ecm)/2d0
!        E3    = (E3max-E3min)*x1 + E3min
!        Jac1  = (E3max-E3min)
!
!
!
! c      2nd integration variable: E_1 with E1min<=E1<=E1max
!        k3abs= dsqrt(E3**2 - yetiM(3)**2)
!        sig  = Ecm - E3
!        tau  = sig**2 - k3abs**2
!        mp   = yetiM(1) + yetiM(2)
!        mm   = yetiM(1) - yetiM(2)
!        dummy1 = 0.5d0/tau * sig*(tau+mp*mm)
!        dummy2 = 0.5d0/tau * k3abs*dsqrt((tau-mp**2)*(tau-mm**2))
!        E1min  = dummy1 - dummy2
!        E1max  = dummy1 + dummy2
!        E1     = (E1max-E1min)*x2 + E1min
!        Jac2   = (E1max-E1min)
!
!
! c      3rd integration variable: cos(theta_3) with -1<=cos(theta_3)<=+1
!        cosT3 = 2d0*x3 - 1d0
!        sinT3 = dsqrt(1d0-cosT3**2)
!        Jac3  = 2d0
!
! c      4th integration variable: eta with 0<=eta<=2*Pi
!        eta  = 2d0*dblPi * x4
!        Jac4 = 2d0*dblPi
!
! c      set p3
!        yetiP(3,0) = E3
!        yetiP(3,1) = k3abs * sinT3
!        yetiP(3,2) = 0d0
!        yetiP(3,3) = k3abs * cosT3
!
! c      set p1
!        k1abs  = dsqrt(E1**2 - yetiM(1)**2)
!        cosEta = dcos(eta)
! !        sinEta = dsin(eta)
!        sinEta = dsqrt(1d0-cosEta**2)
!
! c      usual formula for cosXi:
! !        cosXi  = ( yetiM(1)**2-yetiM(2)**2+yetiM(3)**2
! !      #             +2d0*E1*E3+Ecm**2-2d0*Ecm*(E1+E3)
! !      #            ) / 2d0 / k1abs / k3abs
! c      more accurate evaluation for cosXi:
!        cosXi  = ( yetiM(1)**2-yetiM(2)**2+yetiM(3)**2
!      #           +2d0*Ecm*(Ecm/2d0-E1-E3) )/(2d0*k1abs*k3abs)
!      # + dsqrt((1d0+(yetiM(1)/k1abs)**2)*(1d0+(yetiM(3)/k3abs)**2))
!
! c      safety exit for |cosXi|=1d0 + 1d-16
!        if( cosXi .gt.  1d0 ) cosXi= 1d0
!        if( cosXi .lt. -1d0 ) cosXi=-1d0
!        sinXi  = dsqrt(1d0-cosXi**2)
!
!        yetiP(1,0) = E1
!        yetiP(1,1) = k1abs * ( cosT3*cosEta*sinXi + sinT3*cosXi )
!        yetiP(1,2) = k1abs * ( sinEta*sinXi )
!        yetiP(1,3) = k1abs * (-sinT3*cosEta*sinXi + cosT3*cosXi )
!
! c      set p2
!        yetiP(2,0) = Ecm - E1 - E3
!        yetiP(2,1) = - yetiP(1,1) - yetiP(3,1)
!        yetiP(2,2) = - yetiP(1,2) - yetiP(3,2)
!        yetiP(2,3) = - yetiP(1,3) - yetiP(3,3)
!
! c      set integration weight
!        Norm = 1d0/8d0 * (2d0*dblPi)**(-4)
!        Wgt  = Norm * Jac1*Jac2*Jac3*Jac4
!
!
! c      check on-shell relations and energy-momentum conservation
!        if( check ) then
!         OS1=(yetiM(1)**2
!      #     -(yetiP(1,0)**2-yetiP(1,1)**2-yetiP(1,2)**2-yetiP(1,3)**2 )
!      #      )/Ecm**2
!         OS2=(yetiM(2)**2
!      #     -(yetiP(2,0)**2-yetiP(2,1)**2-yetiP(2,2)**2-yetiP(2,3)**2 )
!      #      )/Ecm**2
!         OS3=(yetiM(3)**2
!      #     -(yetiP(3,0)**2-yetiP(3,1)**2-yetiP(3,2)**2-yetiP(3,3)**2 )
!      #      )/Ecm**2
!         EC =(Ecm-yetiP(1,0)-yetiP(2,0)-yetiP(3,0))/Ecm
!         MC1=(yetiP(1,1)+yetiP(2,1)+yetiP(3,1))/Ecm
!         MC2=(yetiP(1,2)+yetiP(2,2)+yetiP(3,2))/Ecm
!         MC3=(yetiP(1,3)+yetiP(2,3)+yetiP(3,3))/Ecm
!         if( dabs(OS1).gt.5d-16 ) print *,'OS1=',OS1
!         if( dabs(OS2).gt.5d-16 ) print *,'OS2=',OS2
!         if( dabs(OS3).gt.5d-16 ) print *,'OS3=',OS3
!         if( dabs( EC).gt.1d-16 ) print *,' EC=', EC
!         if( dabs(MC1).gt.1d-16 ) print *,'MC1=',MC1
!         if( dabs(MC2).gt.1d-16 ) print *,'MC3=',MC2
!         if( dabs(MC3).gt.1d-16 ) print *,'MC3=',MC3
!        endif
!
!        return
!        end
