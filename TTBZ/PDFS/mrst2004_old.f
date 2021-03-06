      subroutine mrst2004(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
C***************************************************************C
C                                                               C
C  This is a package for the new MRST 2004 'physical gluon' NLO C
C  and NNLO parton distributions.                               C
C  Reference: A.D. Martin, R.G. Roberts, W.J. Stirling and      C
C  R.S. Thorne, hep-ph/0410230                                  C
C                                                               C
C  There are 2 pdf sets corresponding to mode = 1, 2            C
C                                                               C
C  Mode=1 gives the NLO set with Lambda(4) = 347 MeV            C
C  This set reads a grid called mrst2004nlo.dat                 C
C  whose first number is 0.00910                                C
C                                                               C
C  Mode=2 gives the NNLO set with Lambda(4) = 251 MeV           C
C  This set reads a grid called mrst2004nnlo.dat                C
C  whose first number is 0.00673                                C
C                                                               C
C  These fits use a new, physically motivated parametrisation   C
C  for the gluon at the starting scale, Q_0^2 = 1 GeV^2         C
C                                                               C
C         Comments to : W.J.Stirling@durham.ac.uk               C
C                                                               C
C***************************************************************C
      implicit real*8(a-h,o-z)
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      q2=q*q
      if(q2.lt.qsqmin.or.q2.gt.qsqmax) print 99,q2
      if(x.lt.xmin.or.x.gt.xmax)       print 98,x
          if(mode.eq.1) then
        call mrst1(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
      elseif(mode.eq.2) then
        call mrst2(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
      endif
  99  format('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ','q2= ',f10.5)
  98  format('  WARNING:   X  VALUE IS OUT OF RANGE   ','x= ',f10.5)
      return
      end

      subroutine mrst1(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
       open(unit=33,file='/home/schulze/lib/PDFS/mrst2004nlo.dat'
     .       ,status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe1(nx,nq,xxl,qql,f1,cc1)
      call jeppe1(nx,nq,xxl,qql,f2,cc2)
      call jeppe1(nx,nq,xxl,qql,f3,cc3)
      call jeppe1(nx,nq,xxl,qql,f4,cc4)
      call jeppe1(nx,nq,xxl,qql,f6,cc6)
      call jeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue

      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then
      call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then
      call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end

      subroutine mrst2(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='mrst2004nnlo.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe1(nx,nq,xxl,qql,f1,cc1)
      call jeppe1(nx,nq,xxl,qql,f2,cc2)
      call jeppe1(nx,nq,xxl,qql,f3,cc3)
      call jeppe1(nx,nq,xxl,qql,f4,cc4)
      call jeppe1(nx,nq,xxl,qql,f6,cc6)
      call jeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue

      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then
      call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then
      call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end


!       subroutine jeppe2(x,y,nx,my,xx,yy,cc,z)
! C--   G.W. 02/07/2007 Allow extrapolation to small x and large q.
!       implicit real*8(a-h,o-z)
!       dimension xx(nx),yy(my),cc(nx,my,4,4)
!
!       n=locx(xx,nx,x)
!       m=locx(yy,my,y)
!
!       if (n.gt.0.and.n.lt.nx.and.m.gt.0.and.m.lt.my) then
! C--   Do usual interpolation.
!          t=(x-xx(n))/(xx(n+1)-xx(n))
!          u=(y-yy(m))/(yy(m+1)-yy(m))
!          z=0.d0
!          do l=4,1,-1
!             z=t*z+((cc(n,m,l,4)*u+cc(n,m,l,3))*u
!      &           +cc(n,m,l,2))*u+cc(n,m,l,1)
!          enddo
!
!       else if (n.eq.0.and.m.gt.0.and.m.lt.my) then
! C--   Extrapolate to small x.
!          call jeppe3(xx(1),y,nx,my,xx,yy,cc,f0)
!          call jeppe3(xx(2),y,nx,my,xx,yy,cc,f1)
!          if (f0.gt.0.d0.and.f1.gt.0.d0) then
!             z = exp(log(f0)+(log(f1)-log(f0))/(xx(2)-xx(1))*(x-xx(1)))
!          else
!             z = f0+(f1-f0)/(xx(2)-xx(1))*(x-xx(1))
!          end if
!
!       else if (n.gt.0.and.m.eq.my) then
! C--   Extrapolate to large q.
!          call jeppe3(x,yy(my),nx,my,xx,yy,cc,f0)
!          call jeppe3(x,yy(my-1),nx,my,xx,yy,cc,f1)
!          if (f0.gt.0.d0.and.f1.gt.0.d0) then
!             z = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))*
!      &           (y-yy(my)))
!          else
!             z = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
!          end if
!
!       else if (n.eq.0.and.m.eq.my) then
! C--   Extrapolate to small x AND large q.
!          call jeppe3(xx(1),yy(my),nx,my,xx,yy,cc,f0)
!          call jeppe3(xx(1),yy(my-1),nx,my,xx,yy,cc,f1)
!          if (f0.gt.0.d0.and.f1.gt.0.d0) then
!             z0 = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))*
!      &           (y-yy(my)))
!          else
!             z0 = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
!          end if
!          call jeppe3(xx(2),yy(my),nx,my,xx,yy,cc,f0)
!          call jeppe3(xx(2),yy(my-1),nx,my,xx,yy,cc,f1)
!          if (f0.gt.0.d0.and.f1.gt.0.d0) then
!             z1 = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))*
!      &           (y-yy(my)))
!          else
!             z1 = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
!          end if
!          if (z0.gt.0.d0.and.z1.gt.0.d0) then
!             z = exp(log(z0)+(log(z1)-log(z0))/(xx(2)-xx(1))*(x-xx(1)))
!          else
!             z = z0+(z1-z0)/(xx(2)-xx(1))*(x-xx(1))
!          end if
!
!       else
! C--   Set parton distribution to zero otherwise.
!          z = 0.d0
!
!       end if
!
!       return
!       end

C--   G.W. 02/07/2007 Copy of the original jeppe2,
C--   only used for extrapolation.
!       subroutine jeppe3(x,y,nx,my,xx,yy,cc,z)
!       implicit real*8(a-h,o-z)
!       dimension xx(nx),yy(my),cc(nx,my,4,4)
!       n=locx(xx,nx,x)
!       m=locx(yy,my,y)
!       t=(x-xx(n))/(xx(n+1)-xx(n))
!       u=(y-yy(m))/(yy(m+1)-yy(m))
!       z=0.d0
!       do l=4,1,-1
!          z=t*z+((cc(n,m,l,4)*u+cc(n,m,l,3))*u
!      &        +cc(n,m,l,2))*u+cc(n,m,l,1)
!       enddo
!       return
!       end
!
!       integer function locx(xx,nx,x)
!       implicit real*8(a-h,o-z)
!       dimension xx(nx)
! c$$$      if(x.le.xx(1)) then
!       if(x.eq.xx(1)) then ! G.W. 02/07/2007
!       locx=1
!       return
!       endif
! c$$$      if(x.ge.xx(nx)) then
!       if(x.eq.xx(nx)) then ! G.W. 02/07/2007
!       locx=nx-1
!       return
!       endif
!       ju=nx+1
!       jl=0
!     1 if((ju-jl).le.1) go to 2
!       jm=(ju+jl)/2
!       if(x.ge.xx(jm)) then
!       jl=jm
!       else
!       ju=jm
!       endif
!       go to 1
!     2 locx=jl
!       return
!       end

!       real*8 function  polderiv(x1,x2,x3,y1,y2,y3)
!       implicit real*8(a-h,o-z)
!       polderiv=(x3*x3*(y1-y2)-2.0*x2*(x3*(y1-y2)+x1*
!      .(y2-y3))+x2*x2*(y1-y3)+x1*x1*(y2-y3))/((x1-x2)*(x1-x3)*(x2-x3))
!       return
!       end
