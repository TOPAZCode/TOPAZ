MODULE ModIntegrals
implicit none

private :: tBub,vBub

contains






!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     The result in terms of master integrals
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE EvalMasterIntegrals(ThePrimAmp,Mu2)
      use ModAmplitudes
      use ModParameters
      use ModMisc
      use ModProcess
      implicit none
      type(PrimitiveAmplitude),target :: ThePrimAmp
      include 'misc/global_import'
      integer i,j,j1,xe,ia,ib
      complex(8),pointer :: res(:)
      complex(8) d1,d5,c1,c8
      complex(8) b1,b10,b2,b5,pr1,pr2
      complex(8) a1,r1
      complex(8) k1d(16)
      complex(8) k1(16),k2(16)
      complex(8) p1(4),p2(4),p3(4),p4(4)
      complex(8) p12(4),p23(4)
      complex(8) cs11,cs22,cs33,cs44,cs12,cs23
      real(8) s11,s22,s33,s44,s12,s23
      real(8)  m12,m22,m32,m42
      real(8) argk(5,6),argm(5,4)
      real(8) mu2,MassSq
      complex(8)  qlI4,qlI3,qlI2,pr1dc,pr2dc
      complex(8) d(5,1)
!       complex(8) tBub,vBub

      call qlinit

      res => ThePrimAmp%Result

      res(-2)=zero
      res(-1)=zero
      res(0)=zero
      res(1)=zero


! !     5 cut
!       do i=1,N5
!          call penttobox(i,d,argk,argm)
!          do j=1,5
!             s11=argk(j,1)
!             s22=argk(j,2)
!             s33=argk(j,3)
!             s44=argk(j,4)
!             s12=argk(j,5)
!             s23=argk(j,6)
!             m12=argm(j,1)
!             m22=argm(j,2)
!             m32=argm(j,3)
!             m42=argm(j,4)
!
!          do xe=-2,0
!             res(xe) = res(xe)+ d(j,1)*qlI4(s11,s22,s33,s44,s12,s23,m12,m22,m32,m42,mu2,xe)
!          enddo
!
!          enddo
!        enddo



!     4 cut
      do i=1,N4
         if (ThePrimAmp%UCuts(4)%skip(i)) cycle
      d1=ThePrimAmp%UCuts(4)%Coeff(i,0)
      d5=ThePrimAmp%UCuts(4)%Coeff(i,4)

!------ finding the momenta

        do j=1,4
          p1(j)=zero
          p2(j)=zero
          p3(j)=zero
          p4(j)=zero
         enddo

         ia = ThePrimAmp%UCuts(4)%CutProp(i,1)
         ib = ThePrimAmp%UCuts(4)%CutProp(i,2)-1
         do j=ia,ib
         do j1=1,4
           p1(j1) = p1(j1) + mom(j,j1)
         enddo
         enddo


         ia = ThePrimAmp%UCuts(4)%CutProp(i,2)
         ib = ThePrimAmp%UCuts(4)%CutProp(i,3)-1
         do j=ia,ib
         do j1=1,4
           p2(j1) = p2(j1) + mom(j,j1)
         enddo
         enddo


         ia = ThePrimAmp%UCuts(4)%CutProp(i,3)
         ib = ThePrimAmp%UCuts(4)%CutProp(i,4)-1
         do j=ia,ib
         do j1=1,4
           p3(j1) = p3(j1) + mom(j,j1)
         enddo
         enddo


         do j1=1,4
            p4(j1) = -p1(j1)-p2(j1)-p3(j1)
         enddo

         do j=1,4
          p12(j)=p1(j)+p2(j)
          p23(j)=p2(j)+p3(j)
         enddo

         call sc(4,p1,p1,cs11)
         call sc(4,p2,p2,cs22)
         call sc(4,p3,p3,cs33)
         call sc(4,p4,p4,cs44)

         call sc(4,p12,p12,cs12)
         call sc(4,p23,p23,cs23)

         s11=dreal(cs11)
         s22=dreal(cs22)
         s33=dreal(cs33)
         s44=dreal(cs44)

         s12=dreal(cs12)
         s23=dreal(cs23)

         m12=ThePrimAmp%IntPart(ThePrimAmp%UCuts(4)%CutProp(i,1))%Mass2
         m22=ThePrimAmp%IntPart(ThePrimAmp%UCuts(4)%CutProp(i,2))%Mass2
         m32=ThePrimAmp%IntPart(ThePrimAmp%UCuts(4)%CutProp(i,3))%Mass2
         m42=ThePrimAmp%IntPart(ThePrimAmp%UCuts(4)%CutProp(i,4))%Mass2


         MassSq= ExtParticle(1)%Mass2! this assumes that the first particle is always a top or stop


         if (dabs(s11-MassSq).lt.1D-6) then
         s11=MassSq
         endif
         if (dabs(s22-MassSq).lt.1D-6) then
         s22=MassSq
         endif
         if (dabs(s33-MassSq).lt.1D-6) then
         s33=MassSq
         endif
         if (dabs(s44-MassSq).lt.1D-6) then
         s44=MassSq
         endif
         if (dabs(s11).lt.1D-6) then
         s11=0d0
         endif
         if (dabs(s22).lt.1D-6) then
         s22=0d0
         endif
         if (dabs(s33).lt.1D-6) then
         s33=0d0
         endif
         if (dabs(s44).lt.1D-6) then
         s44=0d0
         endif


         do xe=-2,0
            res(xe) = res(xe)+ d1*qlI4(s11,s22,s33,s44,s12,s23,m12,m22,m32,m42,mu2,xe)
!             if(xe.eq.-2) print *, "D0: ",i,d1,qlI4(s11,s22,s33,s44,s12,s23,m12,m22,m32,m42,mu2,xe)
         enddo
         res(1) = res(1)+d5*dcmplx(-1d0/6d0,0d0)
!  	  print *, qlI4(s11,s22,s33,s44,s12,s23,m12,m22,m32,m42,mu2,0)
         enddo
!         print *, 'res4', res(-1)

!     3 cut
      do i=1,N3
         if (ThePrimAmp%UCuts(3)%skip(i)) cycle
      c1=ThePrimAmp%UCuts(3)%Coeff(i,0)
      c8=ThePrimAmp%UCuts(3)%Coeff(i,7)

!------finding the momenta

        do j=1,4
          p1(j)=zero
          p2(j)=zero
          p3(j)=zero
         enddo

         ia = ThePrimAmp%UCuts(3)%CutProp(i,1)
         ib = ThePrimAmp%UCuts(3)%CutProp(i,2)-1

         do j=ia,ib
         do j1=1,4
           p1(j1) = p1(j1) + mom(j,j1)
         enddo
         enddo


         ia = ThePrimAmp%UCuts(3)%CutProp(i,2)
         ib = ThePrimAmp%UCuts(3)%CutProp(i,3)-1

         do j=ia,ib
         do j1=1,4
           p2(j1) = p2(j1) + mom(j,j1)
         enddo
         enddo


         do j1=1,4
            p3(j1) = -p1(j1)-p2(j1)
         enddo

         call sc(4,p1,p1,cs11)
         call sc(4,p2,p2,cs22)
         call sc(4,p3,p3,cs33)



         s11=dreal(cs11)
         s22=dreal(cs22)
         s33=dreal(cs33)

         m12=ThePrimAmp%IntPart(ThePrimAmp%UCuts(3)%CutProp(i,1))%Mass2
         m22=ThePrimAmp%IntPart(ThePrimAmp%UCuts(3)%CutProp(i,2))%Mass2
         m32=ThePrimAmp%IntPart(ThePrimAmp%UCuts(3)%CutProp(i,3))%Mass2


         MassSq= ExtParticle(1)%Mass2! this assumes that the first particle is always a top or stop



         if (dabs(s11-MassSq).lt.1D-6) then
         s11=MassSq
         endif

         if (dabs(s22-MassSq).lt.1D-6) then
         s22=MassSq
         endif

         if (dabs(s33-MassSq).lt.1D-6) then
         s33=MassSq
         endif

         if (dabs(s11).lt.1D-6) then
         s11=0d0
         endif

         if (dabs(s22).lt.1D-6) then
         s22=0d0
         endif

         if (dabs(s33).lt.1D-6) then
         s33=0d0
         endif


         do xe=-2,0
          res(xe) = res(xe) + c1*qlI3(s11,s22,s33,m12,m22,m32,mu2,xe)
!            if(xe.eq.-2) print *, "C0: ",i,c1,qlI3(s11,s22,s33,m12,m22,m32,mu2,xe)
         enddo
         res(1) = res(1) + c8*dcmplx(-0.5d0,0d0)

         enddo

!      print *, 'res3', res(-1)
!     2 cut
      do i=1,N2
         if (ThePrimAmp%UCuts(2)%skip(i)) cycle

!      if (tagdcut(i,1).eq.666) then
         if (ThePrimAmp%UCuts(2)%tagcuts(i) .eq. 666) then
      !-----buble off the light cone

      b1 =ThePrimAmp%UCuts(2)%Coeff(i,0)
      b10=ThePrimAmp%UCuts(2)%Coeff(i,9)

!      finding the momenta
        do j=1,4
          p1(j)=zero
          p2(j)=zero
        enddo

         ia = ThePrimAmp%UCuts(2)%CutProp(i,1)
         ib = ThePrimAmp%UCuts(2)%CutProp(i,2)-1
         do j=ia,ib
         do j1=1,4
           p1(j1) = p1(j1) + mom(j,j1)
         enddo
         enddo

         do j1=1,4
            p2(j1) = -p1(j1)
         enddo

         call sc(4,p1,p1,cs11)

         s11=dreal(cs11)

         m12=ThePrimAmp%IntPart(ThePrimAmp%UCuts(2)%CutProp(i,1))%Mass2
         m22=ThePrimAmp%IntPart(ThePrimAmp%UCuts(2)%CutProp(i,2))%Mass2


         do xe=-2,0
            res(xe) = res(xe) + b1*qlI2(s11,m12,m22,mu2,xe)
!             if(xe.eq.-1) print *, "B0: ",i,b1,qlI2(s11,m12,m22,mu2,xe)
         enddo
         endif


!      if (tagdcut(i,1).eq.999) then
         if (ThePrimAmp%UCuts(2)%tagcuts(i) .eq. 999) then
      !-----bubble on a light-cone

      b1=ThePrimAmp%UCuts(2)%Coeff(i,0)
      b2=ThePrimAmp%UCuts(2)%Coeff(i,1)
      b5=ThePrimAmp%UCuts(2)%Coeff(i,4)
      b10=ThePrimAmp%UCuts(2)%Coeff(i,9)

       do j=1,4
       k1(j)=propv2(i,j)
       k2(j)=propv2(i,4+j)
       k1d(j)=refvect2(i,j)
       enddo


!      finding the momenta

        do j=1,4
          p1(j)=zero
          p2(j)=zero
         enddo

         ia = ThePrimAmp%UCuts(2)%CutProp(i,1)
         ib = ThePrimAmp%UCuts(2)%CutProp(i,2)-1
         do j=ia,ib
         do j1=1,4
           p1(j1) = p1(j1) + mom(j,j1)
         enddo
       enddo

         do j1=1,4
            p2(j1) = -p1(j1)
         enddo

         call sc(4,p1,p1,cs11)
         s11=dreal(cs11)

         m12=ThePrimAmp%IntPart(ThePrimAmp%UCuts(2)%CutProp(i,1))%Mass2
         m22=ThePrimAmp%IntPart(ThePrimAmp%UCuts(2)%CutProp(i,2))%Mass2


         call sc(4,k1,k1d,pr1)
         call sc(4,k2,k1d,pr2)
         pr1dc=-pr1
         pr2dc=-pr2


         do xe=-2,0
            res(xe) = res(xe)+ b1*qlI2(s11,m12,m22,mu2,xe)+ b2*vBub(pr1dc,pr2dc,m12,m22,mu2,xe)+ b5*tBub(pr1dc,pr2dc,m12,m22,mu2,xe)
!             if(xe.eq.-1) print *, "B0(lightcone): ",i,b1,qlI2(s11,m12,m22,mu2,xe)
!             if(xe.eq.-1) print *, "B0(lightcone): ",i,b2,vBub(pr1dc,pr2dc,m12,m22,mu2,xe)
!             if(xe.eq.-1) print *, "B0(lightcone): ",i,b5,tBub(pr1dc,pr2dc,m12,m22,mu2,xe)
         enddo

         endif

!--------contribution to the rational part

         res(1) = res(1) +b10*dcmplx(-0.5d0*(m12+m22-s11/3d0),0d0)
      enddo
!      print *, 'res2',res(-1)

! print *, "checker",qlI2(m_stop**2,0d0,m_stop**2,m_stop**2,-1)
! print *, "checker",qlI2(m_stop**2,0d0,m_stop**2,m_stop**2,0)

! print *, "checker",qlI2(m_stop**2,1d8*m_stop**2,m_top**2,m_stop**2,-1)
! print *, "checker",qlI2(m_stop**2,(1d4*m_stop)**2,m_top**2,m_stop**2,0)
! print *,  1d0 - dlog((1d4*m_stop)**2/m_stop**2) 

! call ffxdb0(res(0),res(1),m_stop**2,1d8*m_stop**2,m_top**2,xe)! m*dB
! print *, "checker",2d0*res(0)*( m_stop**2-1d8*m_stop**2-m_top**2 )

! call ffxdb0(res(0),res(1),m_stop**2,1d-10,m_stop**2,xe)! m*dB
! print *, "checker",res(1)


! print *, "checker",qlI3(m_stop**2,m_stop**2,(1000d0*GeV)**2,m_top**2,1d12*m_stop**2,m_top**2,MuRen**2,-2)
! print *, "checker",qlI3(m_stop**2,m_stop**2,(1000d0*GeV)**2,m_top**2,1d12*m_stop**2,m_top**2,MuRen**2,-1)
! print *, "checker",qlI3(m_stop**2,m_stop**2,(1000d0*GeV)**2,m_top**2,1d12*m_stop**2,m_top**2,MuRen**2,0)
! pause


!     1 cut
      do i=1,N1
         if (ThePrimAmp%UCuts(1)%skip(i)) cycle
         a1=ThePrimAmp%UCuts(1)%Coeff(i,0)
         m12=ThePrimAmp%IntPart(ThePrimAmp%UCuts(1)%CutProp(i,1))%Mass2
! print *, "A0: ",i,a1,dcmplx(m12,0d0)

         res(-1)=res(-1) + a1*dcmplx(m12,0d0)
         res(0)=res(0) + a1*dcmplx(m12*(1d0-dlog(m12/mu2)),0d0)
      enddo

       return
       END SUBROUTINE EvalMasterIntegrals





!-----auxiliary master integrals (vector & tensor bubbles on a light-cone)
!---- currently, for equal masses only
      FUNCTION vBub(pr1,pr2,m12,m22,mu2,xe)
      implicit none
      integer xe
      complex(8)  pr1,pr2,vBub
      real(8) m12,m22,mu2

      if (xe.eq.-2) then
         vBub=dcmplx(0d0,0d0)
      endif

      if (xe.eq.-1) then
      vBub =-0.5d0*pr2-0.5d0*pr1
      endif

      if (xe.eq.0) then
      vBub =0.5d0*(pr2+pr1)*dlog(m12/mu2)
      endif

      if(m12.ne.m22) print *,"There may be an error in vBub",m12,m22

      return
      END FUNCTION vBub





      FUNCTION tBub(pr1,pr2,m12,m22,mu2,xe)
      implicit none
      integer xe
      complex(8)  pr1,pr2,tBub
      real(8) m12,m22,mu2


      if (xe.eq.-2) then
         tBub=dcmplx(0d0,0d0)
      endif

      if (xe.eq.-1) then
      tBub =-1d0/3d0*(pr1-pr2)**2-pr2*(pr1-pr2)-pr2**2
      endif

      if (xe.eq.0) then
      tBub = (pr2**2+ pr2*(pr1-pr2)+1d0/3d0*(pr1-pr2)**2)*dlog(m12/mu2)
      endif

      if(m12.ne.m22) print *,"There may be an error in tBub"

      return
      END FUNCTION tBub








!         SUBROUTINE  mymatch45a(l4cut,l5cut,ns)
!         implicit none
!         include 'global_import'
!         integer i,j1,j2,ns,xtot
!         integer l5cut(5),l4cut(4)
!
!          ns=0
!
!          xtot = 0
!
!             do j1=1,4
!             do j2=1,5
!                    if (l4cut(j1).eq.l5cut(j2)) then
!                     xtot= xtot + 1
!                     endif
!               enddo
!             enddo
!
!          if (xtot.eq.4) then
!            ns = ns+1
!          endif
!         return
!         END SUBROUTINE


!        SUBROUTINE penttobox(i,d,argk,argm)
!        use ModNVBasis
!        use ModMisc
!        use ModProcess
!        use ModParameters
!        implicit none
!        include 'global_import'
!        integer i,j,pos1,pos2,pos3,pos4,j1,ns
!        integer nops,ia,ib,j2
!        integer l5cut(5),l4cut(4)
!        complex(8) res(-2:1)
!        complex(8) summ,sump
!        complex(8) vprop(4),propa,ltr,lv4(4)
!        complex(8) v1(4),v2(4),v3(4),v4(4),v(1:4,1:4)
!        real(8) KMom(1:3,1:4),VMom(1:3,1:4)
!        complex(8) NMom(1,1:4)
!        complex(8) V123(4)
!        real(8) m(4),argm(5,4),argk(5,6)
!        complex(8) p1(4),p2(4),p3(4),p4(4)
!        complex(8) p12(4),p23(4)
!        complex(8) cs11,cs22,cs33,cs44
!        complex(8) cs12,cs23
!        real(8) m12,m22,m32,m42
!        real(8) s11,s22,s33,s44,s12,s23
!        complex(8) e(1),d(5,1),r1,r2,r3,x1,x2,x3
!
!
!           d(1:5,1) = (0d0,0d0)
!           do j1=1,5
!             l5cut(j1)=Lc5(i,j1)
!           enddo
!
!           e(1)=coeff5(i,1)
!           nops=0
!
!
!
!           do j=1,N4
!
!           do j1=1,4
!             l4cut(j1)=Lc4(j,j1)
!           enddo
!
!
!           call mymatch45a(l4cut,l5cut,ns)
!           if (ns.eq.0) then
!               goto 111
!           else
!
!
!          nops = nops+1
!          call mismatch4_5(l4cut,i,pos1)
!
!          KMom(1,1:4)=momline(l4cut(2),1:4)-momline(l4cut(1),1:4)
!          KMom(2,1:4)=momline(l4cut(3),1:4)-momline(l4cut(1),1:4)
!          KMom(3,1:4)=momline(l4cut(4),1:4)-momline(l4cut(1),1:4)
!
!
!
!       do j1=1,4
!             if (Lab_in(l4cut(j1)).eq.'glu') then
!                m(j1)=0d0
!             else
!                m(j1)=m_Top
!             endif
!       enddo
!
!       call NVBasis3(KMom(1:3,1:4),VMom(1:3,1:4),NMom(1,1:4))
!
!
!       v1(1:4)=VMom(1,1:4)
!       v2(1:4)=VMom(2,1:4)
!       v3(1:4)=VMom(3,1:4)
!       v4(1:4)=NMom(1,1:4)
!
!       r1=KMom(1,1:4).dot.KMom(1,1:4)
!       r2=KMom(2,1:4).dot.KMom(2,1:4)
!       r3=KMom(3,1:4).dot.KMom(3,1:4)
!
!       x1=-0.5d0*(r1+dcmplx(m(1)**2-m(2)**2,0d0))
!       x2=-0.5d0*(r2+dcmplx(m(1)**2-m(3)**2,0d0))
!       x3=-0.5d0*(r3+dcmplx(m(1)**2-m(4)**2,0d0))
!
!       do j1=1,4
!          V123(j1)=x1*v1(j1)+x2*v2(j1)+x3*v3(j1)
!       enddo
!
!       call sc(4,V123,V123,r1)
!       ltr=zsqrt(dcmplx(m(1)**2)-r1)
!
!
! !     4-dim coefficients
!        do j1=1,4
!           lv4(j1)=ltr*v4(j1) + V123(j1)
!        enddo
!        do j1=1,4
!           vprop(j1)=momline(pos1,j1)-momline(l4cut(1),j1)+lv4(j1)
!        enddo
!        call sc(4,vprop,vprop,r1)
!        propa = r1 - dcmplx(mass5(i,pos1)**2,0d0)
!        sump = e(1)/propa
!
!
! !         print *, "Prop",propa
! !         if( cdabs(propa).lt.1d-10 ) then
! !             print *, "Prop",cdabs(propa)
! !             stop
! !         endif
!
!
!        do j1=1,4
!           lv4(j1)= -ltr*v4(j1) + V123(j1)
!        enddo
!        do j1=1,4
!           vprop(j1)=momline(pos1,j1)-momline(l4cut(1),j1)+lv4(j1)
!        enddo
!        call sc(4,vprop,vprop,r1)
!        propa = r1 - dcmplx(mass5(i,pos1)**2,0d0)
!        summ = e(1)/propa
!
!        d(nops,1)=0.5d0*(sump +summ)
!
!
!
! !-----finding the momenta
!
!         do j1=1,4
!           p1(j1)=zero
!           p2(j1)=zero
!           p3(j1)=zero
!           p4(j1)=zero
!          enddo
!
!          ia = Lc4(j,1)
!          ib = Lc4(j,2)-1
!
!          do j1=ia,ib
!          do j2=1,4
!            p1(j2) = p1(j2) + mom(j1,j2)
!          enddo
!        enddo
!
!
!          ia = Lc4(j,2)
!          ib = Lc4(j,3)-1
!
!          do j1=ia,ib
!          do j2=1,4
!            p2(j2) = p2(j2) + mom(j1,j2)
!          enddo
!        enddo
!
!
!          ia = Lc4(j,3)
!          ib = Lc4(j,4)-1
!
!
!
!          do j1=ia,ib
!          do j2=1,4
!            p3(j2) = p3(j2) + mom(j1,j2)
!          enddo
!        enddo
!
!
!          do j1=1,4
!             p4(j1) = -p1(j1)-p2(j1)-p3(j1)
!          enddo
!
!
!          do j1=1,4
!           p12(j1)=p1(j1)+p2(j1)
!           p23(j1)=p2(j1)+p3(j1)
!          enddo
!
!          call sc(4,p1,p1,cs11)
!          call sc(4,p2,p2,cs22)
!          call sc(4,p3,p3,cs33)
!          call sc(4,p4,p4,cs44)
!
!          call sc(4,p12,p12,cs12)
!          call sc(4,p23,p23,cs23)
!
!          s11=dreal(cs11)
!          s22=dreal(cs22)
!          s33=dreal(cs33)
!          s44=dreal(cs44)
!
!          s12=dreal(cs12)
!          s23=dreal(cs23)
!
!          m12=mass4(j,1)**2
!          m22=mass4(j,2)**2
!          m32=mass4(j,3)**2
!          m42=mass4(j,4)**2
!
!
!          if (dabs(s11-m_Top**2).lt.1D-6) then
!          s11=m_Top**2
!          endif
!
!          if (dabs(s22-m_Top**2).lt.1D-6) then
!          s22=m_Top**2
!          endif
!
!          if (dabs(s33-m_Top**2).lt.1D-6) then
!          s33=m_Top**2
!          endif
!
!          if (dabs(s44-m_Top**2).lt.1D-6) then
!          s44=m_Top**2
!          endif
!
!
!          if (dabs(s11).lt.1D-6) then
!          s11=0d0
!          endif
!
!          if (dabs(s22).lt.1D-6) then
!          s22=0d0
!          endif
!
!          if (dabs(s33).lt.1D-6) then
!          s33=0d0
!          endif
!
!          if (dabs(s44).lt.1D-6) then
!          s44=0d0
!          endif
!
!          argk(nops,1)=s11
!          argk(nops,2)=s22
!          argk(nops,3)=s33
!          argk(nops,4)=s44
!          argk(nops,5)=s12
!          argk(nops,6)=s23
!
!          argm(nops,1)=m12
!          argm(nops,2)=m22
!          argm(nops,3)=m32
!          argm(nops,4)=m42
!
!             endif
!  111        continue
!
!          enddo
!
!        return
!        END SUBROUTINE penttobox






END MODULE ModIntegrals

