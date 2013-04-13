MODULE ModSpinCorrel
use ModTopDecay
implicit none

save

contains




SUBROUTINE setLOPDFs(x1,x2,MuPDF,pdf)
use ModParameters
implicit none
real(8) :: x1,x2,PDFScale,MuPDF
real(8) :: upv(1:2),dnv(1:2),usea(1:2),dsea(1:2),str(1:2),chm(1:2),bot(1:2),glu(1:2)
integer,parameter :: swPDF_u=1, swPDF_d=1, swPDF_c=1, swPDF_s=1, swPDF_b=1, swPDF_g=1
real(8) :: pdf(-6:6,1:2)

        PDFScale=MuPDF*100d0
        call mrstlo(x1,PDFScale,1,upv(1),dnv(1),usea(1),dsea(1),str(1),chm(1),bot(1),glu(1))
        call mrstlo(x2,PDFScale,1,upv(2),dnv(2),usea(2),dsea(2),str(2),chm(2),bot(2),glu(2))

        print *, "these PDFs should be updated"; stop;


IF( COLLIDER.EQ.1 ) THEN
!       PROTON CONTENT
        pdf(Up_,1)   = (upv(1) + usea(1))  * swPDF_u / x1
        pdf(AUp_,1)  = usea(1)             * swPDF_u / x1
        pdf(Dn_,1)   = (dnv(1) + dsea(1))  * swPDF_d / x1
        pdf(ADn_,1)  = dsea(1)             * swPDF_d / x1
        pdf(Chm_,1)  = chm(1)              * swPDF_c / x1
        pdf(AChm_,1) = chm(1)              * swPDF_c / x1
        pdf(Str_,1)  = str(1)              * swPDF_s / x1
        pdf(AStr_,1) = str(1)              * swPDF_s / x1
        pdf(Bot_,1)  = bot(1)              * swPDF_b / x1
        pdf(ABot_,1) = bot(1)              * swPDF_b / x1
        pdf(0,1)     = glu(1)              * swPDF_g / x1

!       PROTON CONTENT
        pdf(Up_,2)   = (upv(2) + usea(2))  * swPDF_u / x2
        pdf(AUp_,2)  = usea(2)             * swPDF_u / x2
        pdf(Dn_,2)   = (dnv(2) + dsea(2))  * swPDF_d / x2
        pdf(ADn_,2)  = dsea(2)             * swPDF_d / x2
        pdf(Chm_,2)  = chm(2)              * swPDF_c / x2
        pdf(AChm_,2) = chm(2)              * swPDF_c / x2
        pdf(Str_,2)  = str(2)              * swPDF_s / x2
        pdf(AStr_,2) = str(2)              * swPDF_s / x2
        pdf(Bot_,2)  = bot(2)              * swPDF_b / x2
        pdf(ABot_,2) = bot(2)              * swPDF_b / x2
        pdf(0,2)     = glu(2)              * swPDF_g / x2

ELSEIF( COLLIDER.EQ.2 ) THEN
!       PROTON CONTENT
        pdf(Up_,1)   = (upv(1) + usea(1))  * swPDF_u / x1
        pdf(AUp_,1)  = usea(1)             * swPDF_u / x1
        pdf(Dn_,1)   = (dnv(1) + dsea(1))  * swPDF_d / x1
        pdf(ADn_,1)  = dsea(1)             * swPDF_d / x1
        pdf(Chm_,1)  = chm(1)              * swPDF_c / x1
        pdf(AChm_,1) = chm(1)              * swPDF_c / x1
        pdf(Str_,1)  = str(1)              * swPDF_s / x1
        pdf(AStr_,1) = str(1)              * swPDF_s / x1
        pdf(Bot_,1)  = bot(1)              * swPDF_b / x1
        pdf(ABot_,1) = bot(1)              * swPDF_b / x1
        pdf(0,1)     = glu(1)              * swPDF_g / x1

!       ANTI-PROTON CONTENT
        pdf(Up_,2)   = usea(2)             * swPDF_u / x2
        pdf(AUp_,2)  = (upv(2)+usea(2))    * swPDF_u / x2
        pdf(Dn_,2)   = dsea(2)             * swPDF_d / x2
        pdf(ADn_,2)  = (dnv(2) + dsea(2))  * swPDF_d / x2
        pdf(Chm_,2)  = chm(2)              * swPDF_c / x2
        pdf(AChm_,2) = chm(2)              * swPDF_c / x2
        pdf(Str_,2)  = str(2)              * swPDF_s / x2
        pdf(AStr_,2) = str(2)              * swPDF_s / x2
        pdf(Bot_,2)  = bot(2)              * swPDF_b / x2
        pdf(ABot_,2) = bot(2)              * swPDF_b / x2
        pdf(0,2)     = glu(2)              * swPDF_g / x2
ENDIF

RETURN
END SUBROUTINE







    SUBROUTINE GetAllSolutions(MomDK,kount,Jacobian,MomREC)
    use ModMisc
    use ModParameters
    implicit none
    real(8), intent(in) :: MomDK(1:4,6)
    real(8), intent(out) :: MomREC(4,1:4,6)
    integer, intent(out) :: kount
    integer :: i,j
    real(8) :: ptmiss(4), mom_e1(4), mom_e2(4), mom_b1(4), mom_b2(4)
    real(8) :: ga_1,ga_2,al_1,al_2, seb1, seb2, f1perp, f2perp
    real(8) :: v1(4,4),v2(4,4),C(4)
    real(8) :: a1(4), b1(4), a2(4), b2(4)
    real(8) :: xx(6), deta1, kap1,kap2,kap3
    integer, parameter :: Nscan = 100000
    real(8) :: scan(1:4), h , phi2, phi2_check, phi1_check
    real(8) :: ptest2(4),ptest1(4)
    real(8) :: sol_sin(4),sol_cos(4),J_ne(1:2),J_nebar(1:2),Jacobian(1:4)
    real(8) :: cosphi1,sinphi1,cosphi2,sinphi2,rx,ry,sol1,sol2,sol3,sol4
    logical,parameter :: checkSol=.false.

!------------- missing transverse momentum

    ptmiss(1) = 0d0
    ptmiss(4) = 0d0
    ptmiss(2) = MomDK(2,3)+MomDK(2,6)
    ptmiss(3) = MomDK(3,3)+MomDK(3,6)



!------------ now reconstructing solutions


    ! pick neg.lep. and choose smalles inv.mass with b-jet
    mom_e1(:) = MomDK(:,5)
    mom_e2(:) = MomDK(:,2)
    if( get_MInv2(MomDK(1:4,2)+MomDK(1:4,1)) .lt. get_MInv2(MomDK(1:4,2)+MomDK(1:4,4)) ) then
        mom_b1(:) = MomDK(:,4)
        mom_b2(:) = MomDK(:,1)
    else
        mom_b1(:) = MomDK(:,1)
        mom_b2(:) = MomDK(:,4)
    endif

    seb1 = two*dot(mom_b1,mom_e1)
    seb2 = two*dot(mom_b2,mom_e2)
    al_1 = (m_top**2 - m_W**2)/seb1 - 1d0
    al_2 = (m_top**2 - m_W**2)/seb2 - 1d0

    if(al_1.lt.0d0 .or. al_2.lt.0d0) then! this rejects unphysical solutions from "incorrect" pairings
        kount=0
        return
    endif

    ga_1 = m_W**2/seb1
    ga_2 = m_W**2/seb2

    f1perp = dsqrt(dabs(al_1*ga_1*seb1))
    f2perp = dsqrt(dabs(al_2*ga_2*seb2))


    call give2to4vect(mom_e1,mom_b1,v1)
    call give2to4vect(mom_e2,mom_b2,v2)


!-- after that, the neutrino momentum becomes
!--- nu1 = al_1*mom_e1 + ga_1*mom_b1 + ftr1*(cos(phi)*v1(3,:) +sin(phi)*v1(4,:))
!--- nu2 = al_2*mom_e2 + ga_2*mom_b2 + ftr2*(cos(phi)*v2(3,:) +sin(phi)*v2(4,:))



!--- get at least one value for phi2

!     phi1_check = atan(scr(momDK(:,6),v1(4,:))/scr(momDK(:,6),v1(3,:)))
!    print *, phi1_check
!    pause
!    phi2_check = pi+atan(scr(momDK(:,3),v2(4,:))/scr(momDK(:,3),v2(3,:)))
!    print *, phi2_check
!    pause


!    ptest2(:) =  al_2*mom_e2(:)+ ga_2*mom_b2(:) &
!     + f2perp*v2(3,:)*cos(phi2_check) + f2perp*v2(4,:)*sin(phi2_check)

!    ptest1(:) =  al_1*mom_e1(:)+ ga_1*mom_b1(:) &
!     + f1perp*v1(3,:)*cos(phi1_check) + f1perp*v1(4,:)*sin(phi1_check)


!    print *, scr(ptest1,v1(3,:)), scr(momDK(:,6),v1(3,:))
!    print *, ptest1 - momDK(:,6)
!    pause


!--- calculating the set of equations from constraints of the missing
!--- transverse momentum

    C= 0d0
    do i=2,3
        C(i) = 1d0/m_top*( ptmiss(i) - (al_1*mom_e1(i)+ga_1*mom_b1(i) + al_2*mom_e2(i)+ga_2*mom_b2(i)) )
        a1(i) = f1perp*v1(3,i)/m_top
        b1(i) = f1perp*v1(4,i)/m_top
        a2(i) = f2perp*v2(3,i)/m_top
        b2(i) = f2perp*v2(4,i)/m_top
    enddo

    deta1 = a1(2)*b1(3)-a1(3)*b1(2)
    kap1 = b1(3)**2 +a1(3)**2
    kap2 = b1(2)**2 +a1(2)**2
    kap3 = - 2d0*(b1(3)*b1(2) + a1(2)*a1(3))

    xx(1)= kap1*a2(2)**2 + kap2*a2(3)**2+kap3*a2(2)*a2(3)
    xx(2)= kap1*b2(2)**2 + kap2*b2(3)**2+kap3*b2(2)*b2(3)
    xx(3)= kap1*2d0*a2(2)*b2(2)+kap2*2d0*a2(3)*b2(3)  + kap3*(a2(2)*b2(3)+b2(2)*a2(3))
    xx(4)= kap1*2d0*(-C(2)*a2(2))+kap2*2d0*(-C(3)*a2(3))  +kap3*(-C(2)*a2(3)-C(3)*a2(2))
    xx(5)= kap1*2d0*(-C(2)*b2(2))+kap2*2d0*(-C(3)*b2(3))  +kap3*(-C(2)*b2(3)-C(3)*b2(2))
    xx(6)= kap1*C(2)**2+kap2*C(3)**2+kap3*C(2)*C(3) -deta1**2


!-- the equation is---------------------------------------------------
!   0 = xx1*cos(phi2)**2 + xx2*sin(phi2)**2 + xx3*cos(phi2)*sin(phi2)+ xx4*cos(phi2)+xx5*sin(phi2)+x6
!-------------------------------------------------------------------
!     h = 2d0/dble(Nscan)
!     kount = 0
!
!     phi2 = -pi/2d0
!     scan(2) = xx(1)*dcos(phi2)**2 + xx(2)*dsin(phi2)**2 + xx(3)*dcos(phi2)*dsin(phi2) + xx(4)*dcos(phi2)+xx(5)*dsin(phi2) +xx(6)
!     scan(4) = xx(1)*dcos(phi2)**2 + xx(2)*dsin(phi2)**2 - xx(3)*dcos(phi2)*dsin(phi2) - xx(4)*dcos(phi2)+xx(5)*dsin(phi2) +xx(6)
!     do i=2,Nscan
!        sinphi2 = -1d0+(i-1)*h
!        cosphi2 = dsqrt(1d0-sinphi2**2)
!        scan(1) = (xx(2)-xx(1))*sinphi2**2 + xx(3)*cosphi2*sinphi2 + xx(4)*cosphi2+xx(5)*sinphi2 +(xx(6)+xx(1))
!        if (scan(1)*scan(2).lt.0d0) then
!           kount = kount + 1
!           sol_sin(kount) = sinphi2
!           sol_cos(kount) = cosphi2
!        endif
!        scan(2) = scan(1)
! !------
!        cosphi2 = -cosphi2
!        scan(3) = (xx(2)-xx(1))*sinphi2**2 + xx(3)*cosphi2*sinphi2 + xx(4)*cosphi2+xx(5)*sinphi2 +(xx(6)+xx(1))
!        if (scan(3)*scan(4).lt.0d0) then
!           kount = kount + 1
!           sol_sin(kount) = sinphi2
!           sol_cos(kount) = cosphi2
!        endif
!        scan(4) = scan(3)
!     enddo


    call SolveForPhi(xx(1:6),kount,sol_sin(1:4))
    sol_cos(1:4) = (-xx(6)-xx(5)*sol_sin(1:4)-xx(1)*(1d0-sol_sin(1:4)**2)-xx(2)*sol_sin(1:4)**2)/(xx(3)*sol_sin(1:4)+xx(4))

!--- now, getting all the solutions,storing momenta
    Jacobian(1:4) = 0d0
    do j=1,kount
        cosphi2 = sol_cos(j)
        sinphi2 = sol_sin(j)
        rx=  C(2)-a2(2)*cosphi2-b2(2)*sinphi2
        ry = C(3)-a2(3)*cosphi2-b2(3)*sinphi2
        cosphi1 = 1d0/deta1*(+b1(3)*rx-b1(2)*ry )
        sinphi1 = 1d0/deta1*(-a1(3)*rx+a1(2)*ry)

        momREC(j,:,1) = mom_b2(:)
        momREC(j,:,2) = mom_e2(:)
        momREC(j,:,3) = al_2*mom_e2(:)+ ga_2*mom_b2(:) + f2perp*v2(3,:)*cosphi2 + f2perp*v2(4,:)*sinphi2
        momREC(j,:,4) = mom_b1(:)
        momREC(j,:,5) = mom_e1(:)
        momREC(j,:,6) = al_1*mom_e1(:)+ ga_1*mom_b1(:)  + f1perp*v1(3,:)*cosphi1 + f1perp*v1(4,:)*sinphi1

        J_ne(1) = (momREC(j,1:4,3).dot.v1(4,1:4))*v1(3,2) - (momREC(j,1:4,3).dot.v1(3,1:4))*v1(4,2) ! = Jx_nu
        J_ne(2) = (momREC(j,1:4,3).dot.v1(4,1:4))*v1(3,3) - (momREC(j,1:4,3).dot.v1(3,1:4))*v1(4,3) ! = Jy_nu
        J_nebar(1) = (momREC(j,1:4,6).dot.v2(4,1:4))*v2(3,2) - (momREC(j,1:4,6).dot.v2(3,1:4))*v2(4,2) ! = Jx_nubar
        J_nebar(2) = (momREC(j,1:4,6).dot.v2(4,1:4))*v2(3,3) - (momREC(j,1:4,6).dot.v2(3,1:4))*v2(4,3) ! = Jy_nubar
        Jacobian(j) = 1d0/dabs( J_ne(1)*J_nebar(2) - J_ne(2)*J_nebar(1) )
!         if(Jacobian(j).gt.1d6) print *, "Jacobian:",Jacobian(j)
!         if(Jacobian(j).gt.MaxValue) MaxValue=Jacobian(j)
!         if(Jacobian(j).lt.MinValue) MinValue=Jacobian(j)

        if( checkSol ) then
           if( dabs((xx(1)*cosphi2**2 + xx(2)*sinphi2**2 + xx(3)*cosphi2*sinphi2 + xx(4)*cosphi2 + xx(5)*sinphi2 +xx(6))/(xx(1)+xx(2)+xx(3)+xx(4)+xx(5)+xx(6))).gt.1d-4 ) then
                print *, "solution ",j,(xx(1)*cosphi2**2 + xx(2)*sinphi2**2 + xx(3)*cosphi2*sinphi2 + xx(4)*cosphi2 + xx(5)*sinphi2 +xx(6))/(xx(1)+xx(2)+xx(3)+xx(4)+xx(5)+xx(6))
!                 print *, xx(1)+xx(2)+xx(3)+xx(4)+xx(5)+xx(6)
!             print *, "check1",sol_sin(j)**2+sol_cos(j)**2-1d0
           endif
           if( dabs(sol_sin(j)**2+sol_cos(j)**2-1d0).gt.1d-4 ) print *, "sin^2 relation ",sol_sin(j)**2+sol_cos(j)**2-1d0
           if( kount.ne.2 .and. kount.ne.4 ) print *, "num. sol ",kount
        endif

        if(dabs((xx(1)*cosphi2**2 + xx(2)*sinphi2**2 + xx(3)*cosphi2*sinphi2 + xx(4)*cosphi2 + xx(5)*sinphi2 +xx(6))/(xx(1)+xx(2)+xx(3)+xx(4)+xx(5)+xx(6))).gt.1d-4) then
!              print *, "skip solution"
             Jacobian(j) = 0d0
        endif
        if( abs(momREC(j,1:4,6).dot.momREC(j,1:4,6)).gt.1d-4 ) then
!              print *, "skip solution 2"
             Jacobian(j) = 0d0
             kount=0
             return
        endif
        if( abs(momREC(j,1:4,3).dot.momREC(j,1:4,3)).gt.1d-4 ) then
!              print *, "skip solution 1"
             Jacobian(j) = 0d0
             kount=0
             return
        endif
    enddo

    if(all(Jacobian(:).eq.0d0).and.kount.ne.0) then
          kount=0
          return
    endif

    END SUBROUTINE




        SUBROUTINE SolveForPhi(x,NumSol,sinP)
        implicit none
        real(8) :: x(1:6)
        real(8) :: sinP(1:4)
        real(8),parameter :: third=1d0/3d0
        complex(8) :: s1,s2,s3
        complex(8) :: solS,solC,t1,t2,t3,t4
        integer :: NumSol


        NumSol=0
        t1=(x(3)*x(4)+(-x(1)+x(2))*x(5))**2/((x(1)-x(2))**2+x(3)**2)**2+(2d0*(x(3)**2-x(4)**2-x(5)**2+2d0*(x(1)-x(2))*(x(1)+x(6))))/(3d0*((x(1)-x(2))**2+x(3)**2))
        t2=((108d0,0d0)*(x(3)*x(4)+(-x(1)+x(2))*x(5))**2*(x(1)-x(4)+x(6))*(x(1)+x(4)+x(6))-72d0*((x(1)-x(2))**2+x(3)**2)*(x(1)-x(4)+x(6))*(x(1)+x(4)+x(6))*(-x(3)**2+x(4)**2+x(5)**2-2d0*(x(1)-x(2))*(x(1)+x(6)))+2d0*(-x(3)**2+x(4)**2+x(5)**2-2d0*(x(1)-x(2))*(x(1)+x(6)))**3+36d0*(x(3)*x(4)+(-x(1)+x(2))*x(5))*(-x(3)**2+x(4)**2+x(5)**2-2d0*(x(1)-x(2))*(x(1)+x(6)))*(x(3)*x(4)-x(5)*(x(1)+x(6)))+108d0*((x(1)-x(2))**2+x(3)**2)*(x(3)*x(4)-x(5)*(x(1)+x(6)))**2+cdsqrt(-(4d0,0d0)*(12d0*((x(1)-x(2))**2+x(3)**2)*(x(1)-x(4)+x(6))*(x(1)+x(4)+x(6))+(-x(3)**2+x(4)**2+x(5)**2-2d0*(x(1)-x(2))*(x(1)+x(6)))**2+12d0*(x(3)*x(4)+(-x(1)+x(2))*x(5))*(x(3)*x(4)-x(5)*(x(1)+x(6))))**3+4d0*(54d0*(x(3)*x(4)+(-x(1)+x(2))*x(5))**2*(x(1)-x(4)+x(6))*(x(1)+x(4)+x(6))-36d0*((x(1)-x(2))**2+x(3)**2)*(x(1)-x(4)+x(6))*(x(1)+x(4)+x(6))*(-x(3)**2+x(4)**2+x(5)**2-2d0*(x(1)-x(2))*(x(1)+x(6)))+(-x(3)**2+x(4)**2+x(5)**2-2d0*(x(1)-x(2))*(x(1)+x(6)))**3+18d0*(x(3)*x(4)+(-x(1)+x(2))*x(5))*(-x(3)**2+x(4)**2+x(5)**2-2d0*(x(1)-x(2))*(x(1)+x(6)))*(x(3)*x(4)-x(5)*(x(1)+x(6)))+54d0*((x(1)-x(2))**2+x(3)**2)*(x(3)*x(4)-x(5)*(x(1)+x(6)))**2)**2))**third
        t3=t1+t2/(3d0*2**third*((x(1)-x(2))**2+x(3)**2))+(2**third*(12d0*((x(1)-x(2))**2+x(3)**2)*(x(1)-x(4)+x(6))*(x(1)+x(4)+x(6))+(-x(3)**2+x(4)**2+x(5)**2-2d0*(x(1)-x(2))*(x(1)+x(6)))**2+12d0*(x(3)*x(4)+(-x(1)+x(2))*x(5))*(x(3)*x(4)-x(5)*(x(1)+x(6)))))/(3d0*t2*((x(1)-x(2))**2+x(3)**2))
        t4=(8d0*((-(x(3)*x(4))+(x(1)-x(2))*x(5))**3-((x(1)-x(2))**2+x(3)**2)*(-(x(3)*x(4))+(x(1)-x(2))*x(5))*(-x(3)**2+x(4)**2+x(5)**2-2d0*(x(1)-x(2))*(x(1)+x(6)))-2d0*((x(1)-x(2))**2+x(3)**2)**2*(-(x(3)*x(4))+x(5)*(x(1)+x(6)))))/((x(1)-x(2))**2+x(3)**2)**3

        s1 = 0.5d0*cdSqrt(t3)
        s2 = 0.5d0*cdSqrt(3d0*t1-t3-t4/(4d0*cdSqrt(t3)))
        s3 = (-(x(3)*x(4))+x(1)*x(5)-x(2)*x(5))/(2d0*(x(1)**2-2d0*x(1)*x(2)+x(2)**2+x(3)**2))

        solS=-s1 - s2 + s3
        if( abs(dimag(solS)).lt.1d-6 ) then
            NumSol = NumSol+1
            sinP(NumSol) = dble(solS)
        endif

        solS=-s1 + s2 + s3
        if( abs(dimag(solS)).lt.1d-6 ) then
            NumSol = NumSol+1
            sinP(NumSol) = dble(solS)
        endif

        s2 = 0.5d0*cdSqrt(3d0*t1-t3+t4/(4d0*cdSqrt(t3)))
        solS=+s1 - s2 + s3
        if( abs(dimag(solS)).lt.1d-6 ) then
            NumSol = NumSol+1
            sinP(NumSol) = dble(solS)
        endif

        solS=+s1 + s2 + s3
        if( abs(dimag(solS)).lt.1d-6 ) then
            NumSol = NumSol+1
            sinP(NumSol) = dble(solS)
        endif

! if(NumSol.eq.0) then
!     print *, "no real solutions"
!     print *, -s1 - s2 + s3
!     print *, -s1 + s2 + s3
!     print *, +s1 - s2 + s3
!     print *, +s1 + s2 + s3
! pause
! endif

        END SUBROUTINE



FUNCTION calc_rgg(MomProd,MomDK)
use ModProcess
use ModParameters
use ModAmplitudes
use ModMyRecurrence
use ModIntDipoles
implicit none
real(8) :: calc_rgg
logical,save :: first_time=.true.
type(Particle),save :: ExtParticles_tbtgg(1:4)
type(TreeProcess),save :: TreeAmps_tbtgg(1:2)
integer :: iTree,hel1,hel2,hel3,hel4,kcount,NumSol,iPrimAmp,jPrimAmp
real(8) :: MomProd(1:4,1:4),MomExt(1:4,1:4),MomDK(1:4,1:6),MomDKInt(1:4,1:4,1:6),NuJac(1:4),X_c,X_u
real(8) :: RunFactor,eta1,eta2,zBoostMom(1:4),pdf(-6:6,1:2),PDFFac
complex(8) :: Spi(1:4),BarSpi(1:4),Msq_T_BW,Msq_W_ENU,Msq_DK,M_T_BW
complex(8) :: Amp_tbtgg(1:2),LO_Res_Unpol,LO_Res_Pol


     if( first_time ) then
          call InitTrees(2,2,2,TreeAmps_tbtgg)
          call InitProcess_TbTGG(ExtParticles_tbtgg(1:4))
          TreeAmps_tbtgg(1)%PartRef(1:4) = (/1,2,3,4/)
          TreeAmps_tbtgg(2)%PartRef(1:4) = (/1,2,4,3/)
          TreeAmps_tbtgg(1)%NumQua = 2
          TreeAmps_tbtgg(2)%NumQua = 2

          do iTree=1,2
               call LinkTreeParticles(TreeAmps_tbtgg(iTree),ExtParticles_tbtgg(1:4))
          enddo
          ExtParticles_tbtgg(1)%Helicity=0
          ExtParticles_tbtgg(2)%Helicity=0
          first_time=.false.
     endif

     call GetAllSolutions(MomDK,NumSol,NuJac,MomDKInt)
     if( NumSol.eq.0) then
          calc_rgg = -1d0
          return
     endif


! calculation of spin uncorrelated decay
! spherical decay of ATop
!   Msq_T_BW = dsqrt(2d0)*GF*m_top**4*(1d0-m_w**2/m_top**2)*(1d0+2d0*m_w**2/m_top**2)  ! sq.mat.el. for T -> b W
!   Msq_W_ENU = dsqrt(2d0)*GF*m_W**4/3d0*4d0  ! sq.mat.el. for W -> l nu
!   Msq_DK = Msq_T_BW/(2d0*Ga_Top(0)*m_Top)
!   Msq_DK = Msq_DK * Msq_W_ENU/(2d0*Ga_W*m_W)

! spherical decay of Top
!   Msq_DK = Msq_DK * Msq_T_BW/(2d0*Ga_Top(0)*m_Top)
!   Msq_DK = Msq_DK * Msq_W_ENU/(2d0*Ga_W*m_W)

   X_c=0d0; X_u=0d0
   RunFactor = RunAlphaS(NLOParam,MuRen)
   do kcount=1,NumSol
      MomExt(1:4,3) = MomDKInt(kcount,1:4,1)+MomDKInt(kcount,1:4,2)+MomDKInt(kcount,1:4,3)
      MomExt(1:4,4) = MomDKInt(kcount,1:4,4)+MomDKInt(kcount,1:4,5)+MomDKInt(kcount,1:4,6)
      zBoostMom(1:4) = MomExt(1:4,3)+MomExt(1:4,4)
      eta1 = (zBoostMom(1)+zBoostMom(4))/Collider_Energy
      eta2 = (zBoostMom(1)-zBoostMom(4))/Collider_Energy
      if(eta1.ge.1d0 .or. eta2.ge.1d0) cycle
      if(eta1.lt.0d0 .or. eta2.lt.0d0) cycle

      MomExt(1,1) = Collider_Energy/2d0*eta1
      MomExt(2,1) = 0d0
      MomExt(3,1) = 0d0
      MomExt(4,1) = Collider_Energy/2d0*eta1
      MomExt(1,2) = Collider_Energy/2d0*eta2
      MomExt(2,2) = 0d0
      MomExt(3,2) = 0d0
      MomExt(4,2) =-Collider_Energy/2d0*eta2

      call SetLOPDFs(eta1,eta2,m_top,pdf)
      PDFFac = pdf(0,1) * pdf(0,2)
      ExtParticles_tbtgg(1)%Mom(1:4) = MomExt(1:4,3)
      ExtParticles_tbtgg(2)%Mom(1:4) = MomExt(1:4,4)
      ExtParticles_tbtgg(3)%Mom(1:4) =-MomExt(1:4,1)
      ExtParticles_tbtgg(4)%Mom(1:4) =-MomExt(1:4,2)

!     spin correlated sq. matrix element
      call TopDecay(ExtParticles_tbtgg(1),DK_LO,MomDKInt(kcount,1:4,1:3))
      call TopDecay(ExtParticles_tbtgg(2),DK_LO,MomDKInt(kcount,1:4,4:6))

!     ATop
      Spi(1:4) = ExtParticles_tbtgg(1)%Pol(1:4)
      call vBarSpi(ExtParticles_tbtgg(1)%Mom(1:4),m_Top,-1,BarSpi(1:4))
      M_T_BW = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BW = M_T_BW*dconjg(M_T_BW)/2d0

      call vBarSpi(ExtParticles_tbtgg(1)%Mom(1:4),m_Top,+1,BarSpi(1:4))
      M_T_BW = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BW = Msq_T_BW + M_T_BW*dconjg(M_T_BW)/2d0
      Msq_DK = Msq_T_BW

!     Top
      BarSpi(1:4) = ExtParticles_tbtgg(2)%Pol(1:4)
      call uSpi(ExtParticles_tbtgg(2)%Mom(1:4),m_Top,-1,Spi(1:4))
      M_T_BW = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BW = M_T_BW*dconjg(M_T_BW)/2d0

      call uSpi(ExtParticles_tbtgg(2)%Mom(1:4),m_Top,+1,Spi(1:4))
      M_T_BW = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BW = Msq_T_BW + M_T_BW*dconjg(M_T_BW)/2d0
      Msq_DK = Msq_DK * Msq_T_BW


      LO_Res_Unpol = (0d0,0d0)
      do hel3=-1,1,2
      do hel4=-1,1,2
          call pol_mless(ExtParticles_tbtgg(3)%Mom(1:4),hel3,ExtParticles_tbtgg(3)%Pol(1:4))
          call pol_mless(ExtParticles_tbtgg(4)%Mom(1:4),hel4,ExtParticles_tbtgg(4)%Pol(1:4))
          ExtParticles_tbtgg(3)%Helicity=hel3
          ExtParticles_tbtgg(4)%Helicity=hel4
          do iPrimAmp=1,2
              call EvalTree2(TreeAmps_tbtgg(iPrimAmp),Amp_tbtgg(iPrimAmp))
          enddo
          LO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,2
          do iPrimAmp=1,2
              LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * Amp_tbtgg(iPrimAmp)*dconjg(Amp_tbtgg(jPrimAmp))
          enddo
          enddo
          LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol/8d0/8d0/4d0 * (alpha_s4Pi*RunFactor)**2
      enddo!helicity loop
      enddo!helicity loop
      X_c = X_c + LO_Res_Unpol * NuJac(kcount) * PDFFac


!     spin uncorrelated sq. matrix element
      LO_Res_Unpol = (0d0,0d0)
      do hel1=-1,1,2
      do hel2=-1,1,2
      do hel3=-1,1,2
      do hel4=-1,1,2
          call vSpi(ExtParticles_tbtgg(1)%Mom(1:4),m_top,hel1,ExtParticles_tbtgg(1)%Pol(1:4))
          call ubarSpi(ExtParticles_tbtgg(2)%Mom(1:4),m_top,hel2,ExtParticles_tbtgg(2)%Pol(1:4))
          call pol_mless(ExtParticles_tbtgg(3)%Mom(1:4),hel3,ExtParticles_tbtgg(3)%Pol(1:4))
          call pol_mless(ExtParticles_tbtgg(4)%Mom(1:4),hel4,ExtParticles_tbtgg(4)%Pol(1:4))
          ExtParticles_tbtgg(1)%Helicity=hel1
          ExtParticles_tbtgg(2)%Helicity=hel2
          ExtParticles_tbtgg(3)%Helicity=hel3
          ExtParticles_tbtgg(4)%Helicity=hel4

          do iPrimAmp=1,2
              call EvalTree2(TreeAmps_tbtgg(iPrimAmp),Amp_tbtgg(iPrimAmp))
          enddo
          LO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,2
          do iPrimAmp=1,2
              LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * Amp_tbtgg(iPrimAmp)*dconjg(Amp_tbtgg(jPrimAmp))
          enddo
          enddo
          LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol/8d0/8d0/4d0 * (alpha_s4Pi*RunFactor)**2 * Msq_DK
      enddo!helicity loop
      enddo!helicity loop
      enddo!helicity loop
      enddo!helicity loop
      X_u = X_u + LO_Res_Unpol * NuJac(kcount) * PDFFac

! if(kcount.eq.1) then
! write (*,"(A,F20.10,A)") "R=",X_c/(X_c+X_u)
! write (*,"(A,F20.10,A)") "ratio=",X_c/X_u
! write (*,"(A,4F20.10,A)") "Mom(1:4,glu)=(/",MomExt(1:4,1)*100d0,"/)"
! write (*,"(A,4F20.10,A)") "Mom(1:4,glu)=(/",MomExt(1:4,2)*100d0,"/)"
! write (*,"(A,4F20.10,A)") "Mom(1:4,tbar)=(/",MomExt(1:4,3)*100d0,"/)"
! write (*,"(A,4F20.10,A)") "Mom(1:4,t)=(/",MomExt(1:4,4)*100d0,"/)"
! write (*,"(A,4F20.10,A)") "Mom(1:4,bbar)=(/",MomDKInt(kcount,1:4,1)*100d0,"/)"
! write (*,"(A,4F20.10,A)") "Mom(1:4,lep-)=(/",MomDKInt(kcount,1:4,2)*100d0,"/)"
! write (*,"(A,4F20.10,A)") "Mom(1:4,nubar)=(/",MomDKInt(kcount,1:4,3)*100d0,"/)"
! write (*,"(A,4F20.10,A)") "Mom(1:4,b)=(/",MomDKInt(kcount,1:4,4)*100d0,"/)"
! write (*,"(A,4F20.10,A)") "Mom(1:4,lep+)=(/",MomDKInt(kcount,1:4,5)*100d0,"/)"
! write (*,"(A,4F20.10,A)") "Mom(1:4,nu)=(/",MomDKInt(kcount,1:4,6)*100d0,"/)"
! pause
! endif

   enddo! solutions loop


   calc_rgg = X_c/(X_c+X_u)


return
END FUNCTION






FUNCTION calc_rqq(MomProd,MomDK)
use ModProcess
use ModParameters
use ModAmplitudes
use ModMyRecurrence
use ModIntDipoles
implicit none
real(8) :: calc_rqq
logical,save :: first_time=.true.
type(Particle),save :: ExtParticles_tbtqbq(1:4)
type(TreeProcess),save :: TreeAmps_tbtqbq(1:1)
integer :: iTree,hel1,hel2,hel3,hel4,kcount,NumSol,iPrimAmp,jPrimAmp
real(8) :: MomProd(1:4,1:4),MomExt(1:4,1:4),MomDK(1:4,1:6),MomDKInt(1:4,1:4,1:6),NuJac(1:4),X_c,X_u
real(8) :: RunFactor,zBoostMom(1:4),eta1,eta2,pdf(-6:6,1:2),PDFFac
complex(8) :: Amp_tbtqbq(1:2),LO_Res_Unpol,LO_Res_Pol
complex(8) :: Spi(1:4),BarSpi(1:4),Msq_T_BW,Msq_W_ENU,Msq_DK,M_T_BW


     if( first_time ) then
          call InitTrees(4,0,1,TreeAmps_tbtqbq)
          call InitProcess_TbTQbQ(ExtParticles_tbtqbq(1:4))
          TreeAmps_tbtqbq(1)%PartRef(1:4) = (/1,2,3,4/)
          TreeAmps_tbtqbq(1)%NumQua = 4

          do iTree=1,1
               call LinkTreeParticles(TreeAmps_tbtqbq(iTree),ExtParticles_tbtqbq(1:4))
          enddo
          ExtParticles_tbtqbq(1)%Helicity=0
          ExtParticles_tbtqbq(2)%Helicity=0

          first_time=.false.
     endif

     call GetAllSolutions(MomDK,NumSol,NuJac,MomDKInt)
     if( NumSol.eq.0) then
          calc_rqq = -1d0
          return
     endif

   X_c=0d0; X_u=0d0
   RunFactor = RunAlphaS(NLOParam,MuRen)
! calculation of spin uncorrelated process
! spherical decay of ATop
!   Msq_T_BW = dsqrt(2d0)*GF*m_top**4*(1d0-m_w**2/m_top**2)*(1d0+2d0*m_w**2/m_top**2)  ! sq.mat.el. for T -> b W
!   Msq_W_ENU = dsqrt(2d0)*GF*m_W**4/3d0*4d0  ! sq.mat.el. for W -> l nu

!   Msq_DK = Msq_T_BW/(2d0*Ga_Top(0)*m_Top)
!   Msq_DK = Msq_DK * Msq_W_ENU/(2d0*Ga_W*m_W)
! spherical decay of Top
!   Msq_DK = Msq_DK * Msq_T_BW/(2d0*Ga_Top(0)*m_Top)
!   Msq_DK = Msq_DK * Msq_W_ENU/(2d0*Ga_W*m_W)

   do kcount=1,NumSol
      MomExt(1:4,3) = MomDKInt(kcount,1:4,1)+MomDKInt(kcount,1:4,2)+MomDKInt(kcount,1:4,3)
      MomExt(1:4,4) = MomDKInt(kcount,1:4,4)+MomDKInt(kcount,1:4,5)+MomDKInt(kcount,1:4,6)
      zBoostMom(1:4) = MomExt(1:4,3)+MomExt(1:4,4)
      eta1 = (zBoostMom(1)+zBoostMom(4))/Collider_Energy
      eta2 = (zBoostMom(1)-zBoostMom(4))/Collider_Energy
      if(eta1.ge.1d0 .or. eta2.ge.1d0) cycle

      MomExt(1,1) = Collider_Energy/2d0*eta1
      MomExt(2,1) = 0d0
      MomExt(3,1) = 0d0
      MomExt(4,1) = Collider_Energy/2d0*eta1
      MomExt(1,2) = Collider_Energy/2d0*eta2
      MomExt(2,2) = 0d0
      MomExt(3,2) = 0d0
      MomExt(4,2) =-Collider_Energy/2d0*eta2
      call SetLOPDFs(eta1,eta2,m_top,pdf)
      PDFFac = pdf(Up_,1) *pdf(AUp_,2) + pdf(Dn_,1) *pdf(ADn_,2)
      ExtParticles_tbtqbq(1)%Mom(1:4) = MomExt(1:4,3)
      ExtParticles_tbtqbq(2)%Mom(1:4) = MomExt(1:4,4)
      ExtParticles_tbtqbq(3)%Mom(1:4) =-MomExt(1:4,1)
      ExtParticles_tbtqbq(4)%Mom(1:4) =-MomExt(1:4,2)

!     spin correlated sq. matrix element
      call TopDecay(ExtParticles_tbtqbq(1),DK_LO,MomDKInt(kcount,1:4,1:3))
      call TopDecay(ExtParticles_tbtqbq(2),DK_LO,MomDKInt(kcount,1:4,4:6))


!     ATop
      Spi(1:4) = ExtParticles_tbtqbq(1)%Pol(1:4)
      call vBarSpi(ExtParticles_tbtqbq(1)%Mom(1:4),m_Top,-1,BarSpi(1:4))
      M_T_BW = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BW = M_T_BW*dconjg(M_T_BW)/2d0

      call vBarSpi(ExtParticles_tbtqbq(1)%Mom(1:4),m_Top,+1,BarSpi(1:4))
      M_T_BW = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BW = Msq_T_BW + M_T_BW*dconjg(M_T_BW)/2d0
      Msq_DK = Msq_T_BW

!     Top
      BarSpi(1:4) = ExtParticles_tbtqbq(2)%Pol(1:4)
      call uSpi(ExtParticles_tbtqbq(2)%Mom(1:4),m_Top,-1,Spi(1:4))
      M_T_BW = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BW = M_T_BW*dconjg(M_T_BW)/2d0

      call uSpi(ExtParticles_tbtqbq(2)%Mom(1:4),m_Top,+1,Spi(1:4))
      M_T_BW = psp1_(Spi(1:4),BarSpi(1:4))/(2d0*m_Top)
      Msq_T_BW = Msq_T_BW + M_T_BW*dconjg(M_T_BW)/2d0
      Msq_DK = Msq_DK * Msq_T_BW

      LO_Res_Unpol = (0d0,0d0)
      do hel3=-1,1,2
      do hel4=-1,1,2
          call vSpi(ExtParticles_tbtqbq(3)%Mom(1:4),0d0,hel3,ExtParticles_tbtqbq(3)%Pol(1:4))
          call ubarSpi(ExtParticles_tbtqbq(4)%Mom(1:4),0d0,hel4,ExtParticles_tbtqbq(4)%Pol(1:4))
          ExtParticles_tbtqbq(3)%Helicity=hel3
          ExtParticles_tbtqbq(4)%Helicity=hel4

          do iPrimAmp=1,1
              call EvalTree2(TreeAmps_tbtqbq(iPrimAmp),Amp_tbtqbq(iPrimAmp))
          enddo
          LO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,1
          do iPrimAmp=1,1
              LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * Amp_tbtqbq(iPrimAmp)*dconjg(Amp_tbtqbq(jPrimAmp))
          enddo
          enddo
          LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol/3d0/3d0/4d0 * (alpha_s4Pi*RunFactor)**2
      enddo!helicity loop
      enddo!helicity loop
      X_c = X_c + LO_Res_Unpol * NuJac(kcount) * PDFFac


!     spin uncorrelated sq. matrix element
      LO_Res_Unpol = (0d0,0d0)
      do hel1=-1,1,2
      do hel2=-1,1,2
      do hel3=-1,1,2
      do hel4=-1,1,2
          call vSpi(ExtParticles_tbtqbq(1)%Mom(1:4),m_top,hel1,ExtParticles_tbtqbq(1)%Pol(1:4))
          call ubarSpi(ExtParticles_tbtqbq(2)%Mom(1:4),m_top,hel2,ExtParticles_tbtqbq(2)%Pol(1:4))
          call vSpi(ExtParticles_tbtqbq(3)%Mom(1:4),0d0,hel3,ExtParticles_tbtqbq(3)%Pol(1:4))
          call ubarSpi(ExtParticles_tbtqbq(4)%Mom(1:4),0d0,hel4,ExtParticles_tbtqbq(4)%Pol(1:4))
          ExtParticles_tbtqbq(1)%Helicity=hel1
          ExtParticles_tbtqbq(2)%Helicity=hel2
          ExtParticles_tbtqbq(3)%Helicity=hel3
          ExtParticles_tbtqbq(4)%Helicity=hel4

          do iPrimAmp=1,1
              call EvalTree2(TreeAmps_tbtqbq(iPrimAmp),Amp_tbtqbq(iPrimAmp))
          enddo
          LO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,1
          do iPrimAmp=1,1
              LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqb(iPrimAmp,jPrimAmp) * Amp_tbtqbq(iPrimAmp)*dconjg(Amp_tbtqbq(jPrimAmp))
          enddo
          enddo
          LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol/3d0/3d0/4d0 * (alpha_s4Pi*RunFactor)**2 * Msq_DK
      enddo!helicity loop
      enddo!helicity loop
      enddo!helicity loop
      enddo!helicity loop
      X_u = X_u + LO_Res_Unpol * NuJac(kcount) * PDFFac


! if(kcount.eq.1) then
! write (*,"(A,F20.10,A)") "R=",X_c/(X_c+X_u)
! write (*,"(A,F20.10,A)") "ratio=",X_c/X_u
! write (*,"(A,4F20.10,A)") "Mom(1:4,glu)=(/",MomExt(1:4,1)*100d0,"/)"
! write (*,"(A,4F20.10,A)") "Mom(1:4,glu)=(/",MomExt(1:4,2)*100d0,"/)"
! write (*,"(A,4F20.10,A)") "Mom(1:4,tbar)=(/",MomExt(1:4,3)*100d0,"/)"
! write (*,"(A,4F20.10,A)") "Mom(1:4,t)=(/",MomExt(1:4,4)*100d0,"/)"
! write (*,"(A,4F20.10,A)") "Mom(1:4,bbar)=(/",MomDKInt(kcount,1:4,1)*100d0,"/)"
! write (*,"(A,4F20.10,A)") "Mom(1:4,lep-)=(/",MomDKInt(kcount,1:4,2)*100d0,"/)"
! write (*,"(A,4F20.10,A)") "Mom(1:4,nubar)=(/",MomDKInt(kcount,1:4,3)*100d0,"/)"
! write (*,"(A,4F20.10,A)") "Mom(1:4,b)=(/",MomDKInt(kcount,1:4,4)*100d0,"/)"
! write (*,"(A,4F20.10,A)") "Mom(1:4,lep+)=(/",MomDKInt(kcount,1:4,5)*100d0,"/)"
! write (*,"(A,4F20.10,A)") "Mom(1:4,nu)=(/",MomDKInt(kcount,1:4,6)*100d0,"/)"
! pause
! endif


   enddo! solutions loop

   calc_rqq = X_c/(X_c+X_u)
return
END FUNCTION








  subroutine give2vect(k2,k1,v)
    use ModMisc
    implicit none
    real(8), intent(in)  :: k1(4), k2(4)
    real(8), intent(out) :: v(4)

    v = dot(k2,k2)*k1-dot(k1,k2)*k2

  end subroutine give2vect



  subroutine give2to4vect(k1,k2,v)
    use ModMisc
    implicit none
    real(8),  intent(in) :: k1(4),k2(4)
    real(8),  intent(out) ::  v(4,4)
    ! ---------------------------------------------
    real(8) :: v1(4),v2(4),v3(4),v4(4)
    real(8) :: vaux(4)
    real(8) :: d2,ko1,ko2,ko3
    real(8) :: a11,a12,a22,dnorm
    real(8) :: a1v,a2v,a3v,a33

    call give2vect(k2,k1,v1)
    call give2vect(k1,k2,v2)

    d2 = dot(k1,v1)

    v1=v1/d2
    v2=v2/d2


    a11 = dot(k1,k1)
    a12 = dot(k1,k2)
    a22 = dot(k2,k2)

    d2=a12**2-a11*a22
    ko1=a22/d2
    ko2=a11/d2
    ko3=-a12/d2

    vaux(1)=1.3d0
    vaux(2)=1.7d0
    vaux(3)=2.4d0
    vaux(4)=3.5d0

    a1v = dot(k1,vaux)
    a2v = dot(k2,vaux)

    vaux=vaux+ko1*a1v*k1+ko2*a2v*k2+ko3*(a1v*k2 +a2v*k1)


    dnorm = scr(vaux,vaux)
    dnorm = dsqrt(dabs(dnorm))

    v3=vaux/dnorm



    vaux(1)=2.1d0
    vaux(2)=1.2d0
    vaux(3)=3.4d0
    vaux(4)=0.5d0

    a1v = dot(k1,vaux)
    a2v = dot(k2,vaux)

    vaux=vaux+ko1*a1v*k1+ko2*a2v*k2+ko3*(a1v*k2+a2v*k1)

    a3v =  dot(vaux,v3)
    a33 = dot(v3,v3)
    vaux = vaux - a3v/a33*v3

    dnorm = dot(vaux,vaux)
    dnorm=dsqrt(dabs(dnorm))

    v4 = vaux/dnorm

    v(1,:) = v1
    v(2,:) = v2
    v(3,:) = v3
    v(4,:) = v4

  end subroutine give2to4vect








END MODULE