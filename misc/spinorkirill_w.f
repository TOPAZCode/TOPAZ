c-----subroutine for u spinor
      subroutine uspinor(p,m,i,f)
C     Weyl representation
C   Gamma0= 
c            [ 0  0   1    0 ]
c            [               ]
c            [ 0  0   0    1 ]
c            [               ]
c            [ 1  0   0   0  ]
c            [               ]
c            [ 0  1   0   0  ]

c  Gamma1=
c            [  0    0   0  -1 ]
c            [                 ]
c            [  0    0   -1  0 ]
c            [                 ]
c            [  0   +1  0  0   ]
c            [                 ]
c            [ +1   0   0  0   ]
C  Gamma2=
c            [  0    0   0   +%i]
c            [                  ]
c            [  0    0   -%i   0]
c            [                  ]
c            [  0    -%i  0    0]
c            [                  ]
c            [ + %i  0   0    0 ]
c Gamma3= 
c            [  0   0  -1  0 ]
c            [               ]
c            [  0   0  0   1 ]
c            [               ]
c            [  1  0  0   0  ]
c            [               ]
c            [  0  -1  0   0 ]
c 
      implicit none
      integer i
      double precision m
      double complex p(4),f(4),fc,czip
      double precision  px,py,pz
      parameter(czip=(0d0,0d0))

      px=+dreal(p(2))
      py=+dreal(p(3))
      pz=+dreal(p(4))

      fc=sqrt(p(1)+p(4))

      if (i.eq.1) then 
      f(1)=fc
      f(2)=dcmplx(px,py)/fc
      f(3)=dcmplx(m,0d0)/fc
      f(4)=czip
      endif 

      if (i.eq.-1) then 
      f(1)=czip
      f(2)=dcmplx(-m,0d0)/fc
      f(3)=dcmplx(px,-py)/fc
      f(4)=-fc
      endif 
         
      return
      end

c-----subroutine for ubar spinor

      subroutine ubarspinor(p,m,i,f)
C  As usual we modifiy this solution by exchanging (px <--> pz) 
C  and changing the sign of py, (to remain in right handed system)
      implicit none
      integer i
      double precision m
      double complex p(4),f(4),fc,czip
      double precision  p0,px,py,pz
      parameter(czip=(0d0,0d0))

      px=+dreal(p(2))
      py=+dreal(p(3))
      pz=+dreal(p(4))

      fc=sqrt(p(1)+p(4))

      if (i.eq.1) then 
      f(1)=dcmplx(m,0d0)/fc
      f(2)=czip
      f(3)=fc
      f(4)=dcmplx(px,-py)/fc
      endif 

      if (i.eq.-1) then 
      f(1)=dcmplx(px,py)/fc
      f(2)=-fc
      f(3)=czip
      f(4)=-dcmplx(m,0d0)/fc
      endif 

      return
      end


c-----subroutine for anti-quark spinor
C-----these formula are directly copied from BjD formula 3.7
C-----(after modification of overall normalization)
C-----However in the light of BjD (3.16), RKE would prefer to 
C-----exchange (i.eq.1) and (i.eq. -1).
C-----If we sum over i it doesn't matter, but could be of importance to 
C-----identify particular values in table.

      subroutine vspinor(p,m,i,f)
      implicit none
      integer i
      double precision m
      double complex p(4),f(4),fc,czip
      double precision p0,px,py,pz
      parameter(czip=(0d0,0d0))

      px=+dreal(p(2))
      py=+dreal(p(3))
      pz=+dreal(p(4))

      fc=sqrt(p(1)+p(4))

      if (i.eq.-1) then 
      f(1)=czip
      f(2)=dcmplx(m,0d0)/fc
      f(3)=dcmplx(px,-py)/fc
      f(4)=-fc
      endif

      if (i.eq.1) then
      f(1)=fc
      f(2)=dcmplx(px,py)/fc
      f(3)=dcmplx(-m,0d0)/fc
      f(4)=czip
      endif
      return
      end

c-----subroutine for anti-quark spinor
      subroutine Vbspinor(p,m,i,f)
      implicit none
      integer i
      double precision m
      double complex p(4),f(4),fc,czip
      double precision px,py,pz
      parameter(czip=(0d0,0d0))

      px=+dreal(p(2))
      py=+dreal(p(3))
      pz=+dreal(p(4))

      fc=sqrt(p(1)+p(4))

      if (i.eq.-1) then 
      f(1)=dcmplx(px,py)/fc
      f(2)=-fc
      f(3)=czip
      f(4)=dcmplx(m,0d0)/fc
      endif

      if (i.eq.1) then
      f(1)=dcmplx(-m,0d0)/fc
      f(2)=czip
      f(3)=fc
      f(4)=dcmplx(px,-py)/fc
      endif
      return
      end
