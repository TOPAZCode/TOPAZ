! this is the file to subtract dipoles for ttgggg amplitudes
      module ModDipoles_GGTTBGG
      use ModAmplitudes
      use ModProcess
      use ModParameters
      use ModMisc
      use ModKinematics
      use ModIntDipoles
      implicit none

      public :: EvalDipoles_GGTTBGG
      integer, parameter  :: dp = selected_real_kind(15)
      real(dp), private :: yRnDK(1:8), Wgt_ext


      logical, parameter :: invert_alphaCut = .false.

      contains

!  we have a list of dipoles that we need to go through
!  We label things as 0-> bar t(p1)+t(p2) + g(p3)+g(p4)+g(p5) + g(p6)
!  and we  assume that gluons in the initial state have momenta p3 and p4,
!  so only p5 or p6 are in the final state

      subroutine EvalDipoles_GGTTBGG(p,yRnDk1,Wgt,sum_dip)
      real(dp), intent(out) ::  sum_dip
      real(dp), intent(in) :: p(4,6)
      real(dp), intent(in) :: yRnDK1(1:8), Wgt
      integer, parameter :: ndip = 36
      integer, parameter :: in1 = 3
      integer, parameter :: in2 = 4
      real(dp) :: res
      integer ::  dip(ndip,3)
      data dip(1,1)/5/ , dip(1,2)/3/, dip(1,3)/1/
      data dip(2,1)/5/ , dip(2,2)/3/, dip(2,3)/2/
      data dip(3,1)/5/ , dip(3,2)/3/, dip(3,3)/4/
      data dip(4,1)/5/ , dip(4,2)/3/, dip(4,3)/6/
      data dip(5,1)/5/ , dip(5,2)/4/, dip(5,3)/1/
      data dip(6,1)/5/ , dip(6,2)/4/, dip(6,3)/2/
      data dip(7,1)/5/ , dip(7,2)/4/, dip(7,3)/3/
      data dip(8,1)/5/ , dip(8,2)/4/, dip(8,3)/6/
      data dip(9,1)/5/ ,  dip(9,2)/1/, dip(9,3)/2/
      data dip(10,1)/5/, dip(10,2)/1/, dip(10,3)/3/
      data dip(11,1)/5/, dip(11,2)/1/, dip(11,3)/4/
      data dip(12,1)/5/, dip(12,2)/1/, dip(12,3)/6/
      data dip(13,1)/5/ , dip(13,2)/2/, dip(13,3)/1/
      data dip(14,1)/5/ , dip(14,2)/2/, dip(14,3)/3/
      data dip(15,1)/5/ , dip(15,2)/2/, dip(15,3)/4/
      data dip(16,1)/5/ , dip(16,2)/2/, dip(16,3)/6/
      data dip(17,1)/5/ , dip(17,2)/6/, dip(17,3)/1/
      data dip(18,1)/5/ , dip(18,2)/6/, dip(18,3)/2/
      data dip(19,1)/5/ , dip(19,2)/6/, dip(19,3)/3/
      data dip(20,1)/5/ , dip(20,2)/6/, dip(20,3)/4/
      data dip(21,1)/6/ , dip(21,2)/3/, dip(21,3)/1/
      data dip(22,1)/6/ , dip(22,2)/3/, dip(22,3)/2/
      data dip(23,1)/6/ , dip(23,2)/3/, dip(23,3)/4/
      data dip(24,1)/6/ , dip(24,2)/3/, dip(24,3)/5/
      data dip(25,1)/6/ , dip(25,2)/4/, dip(25,3)/1/
      data dip(26,1)/6/ , dip(26,2)/4/, dip(26,3)/2/
      data dip(27,1)/6/ , dip(27,2)/4/, dip(27,3)/3/
      data dip(28,1)/6/ , dip(28,2)/4/, dip(28,3)/5/
      data dip(29,1)/6/, dip(29,2)/1/, dip(29,3)/2/
      data dip(30,1)/6/, dip(30,2)/1/, dip(30,3)/3/
      data dip(31,1)/6/, dip(31,2)/1/, dip(31,3)/4/
      data dip(32,1)/6/, dip(32,2)/1/, dip(32,3)/5/
      data dip(33,1)/6/ , dip(33,2)/2/, dip(33,3)/1/
      data dip(34,1)/6/ , dip(34,2)/2/, dip(34,3)/3/
      data dip(35,1)/6/ , dip(35,2)/2/, dip(35,3)/4/
      data dip(36,1)/6/ , dip(36,2)/2/, dip(36,3)/5/
!      data dip(37,1)/6/ , dip(37,2)/5/, dip(37,3)/1/
!      data dip(38,1)/6/ , dip(38,2)/5/, dip(38,3)/2/
!      data dip(39,1)/6/ , dip(39,2)/5/, dip(39,3)/3/
!      data dip(40,1)/6/ , dip(40,2)/5/, dip(40,3)/4/
      real(dp) :: mass(6)
!-----flavor list qm -- massive quark, qu -- massless quark
!----- gl -- gluon
      character :: fl(6)*2
      integer :: n,i,j,k

       yRnDK = yRnDK1
       Wgt_ext = Wgt


        fl =    'gl'
        fl(1) = 'qm'
        fl(2) = 'qm'

!       print *, "SP"
!       print *, p(1:4,1).dot.p(1:4,1)
!       print *, p(1:4,2).dot.p(1:4,2)
!       print *, p(1:4,3).dot.p(1:4,3)
!       print *, p(1:4,4).dot.p(1:4,4)
!       print *, p(1:4,5).dot.p(1:4,5)
!       print *, p(1:4,6).dot.p(1:4,6)
!       print *, "EMC"
!       print *, p(1,1)+p(1,2)+p(1,3)+p(1,4)+p(1,5)+p(1,6)
!       print *, p(2,1)+p(2,2)+p(2,3)+p(2,4)+p(2,5)+p(2,6)
!       print *, p(3,1)+p(3,2)+p(3,3)+p(3,4)+p(3,5)+p(3,6)
!       print *, p(4,1)+p(4,2)+p(4,3)+p(4,4)+p(4,5)+p(4,6)


!     mass assignment
      mass = zero
      mass(1) = m_top
      mass(2) = m_top

      sum_dip = zero

            do n=1,ndip
!             do n=1,4

         i=dip(n,1)  ! emitted
         j=dip(n,2)  ! emittor
         k=dip(n,3)  ! spectator

         res = zero

!      if( i.eq.6 .and. j.eq.5 .and. k.eq.1 ) then
!         alpha_ff=1d-2
!      else
!         alpha_ff=1d0
!      endif


      if (j.ne.in1.and.j.ne.in2.and.k.ne.in1.and.k.ne.in2) then
        call dipff(n,i,j,k,mass(i),mass(j),mass(k),fl(i),fl(j),p,res)
!         print *, "ff"
      endif

      if ((j.ne.in1.and.j.ne.in2).and.(k.eq.in1.or.k.eq.in2)) then
        call dipfi(n,i,j,k,mass(i),mass(j),mass(k),fl(i),fl(j),p,res)
!         print *, "fi"
      endif


      if ( (k.ne.in1.and.k.ne.in2).and.(j.eq.in1.or.j.eq.in2)) then
        call dipif(n,i,j,k,mass(i),mass(j),mass(k),fl(i),fl(j),p,res)
!         print *, "if"
      endif

      if ( (j.eq.in1.and.k.eq.in2).or.(j.eq.in2.and.k.eq.in1)) then
        call dipii(n,i,j,k,mass(i),mass(j),mass(k),fl(i),fl(j),p,res)
!         print *, "ii"
      endif

        sum_dip = sum_dip + res
!         print *, n,res
        enddo

      end subroutine






!      dipole subroutines

       subroutine dipff(n,i,j,k,mi,mj,mk,fl1,fl2,p,res)
       implicit none
       integer, intent(in) :: n,i,j,k
       real(dp), intent(in) :: mi,mj,mk
       character, intent(in) :: fl1*2,fl2*2
       real(dp), intent(in) :: p(4,6)
       real(dp), intent(out) :: res
       complex(dp) :: cres
       integer :: oh(6,5)
       real(dp) :: mij
       real(dp) :: C(6,6)
       real(dp) :: mui,muj,muk,q2sqrt,yp
       real(dp) ::  pi(4), pj(4), pk(4), pkt(4), pij(4), qq(4), qqij(4)
       real(dp) :: q(4,5), paux(4), pijk(4)
       real(dp) :: xx1, xx2, zi, zj, yijk, vijk, zim
       real(dp) :: pij2,pijk2,tvijk,q2,tpij2, qij2, zjm
       real(dp) ::  pijt(4)
       real(dp) :: zp, zm, xxx, viji
       integer :: i1, i2, i3, i4, i5, i6, hel(5), i8
       character :: fl(5)*2
       integer :: pos, iTree
       real(dp) :: diag, offdiag
       complex(dp) :: xp, xm
       complex(dp) :: HH(-1:1,-1:1)
       type(Particle),save :: ExtParticles(1:5)
       type(TreeProcess),save :: TreeAmpsDip(1:6)
       complex(8)        :: ResAmpsDip_mt(-1:1,1:6)
       integer :: c1, c2, j1
       complex(8)        :: Bm(-1:1,1:6,-1:1,-1:1,-1:1,-1:1)
       complex(8)        :: Am(-1:1,-1:1,1:6,1:6)
       complex(dp) :: POL1(-1:1,4)
       real(dp) ::  pjetout(4,5), weight
       logical :: Not_Passed_Cuts
       integer :: NBin(1:NumHistograms)
       real(dp) :: MomDK(1:4,1:6),PSWgt1,PSWgt2
       integer :: Njet, Nmax(5), Nhisto
       logical, save :: first_time = .true.


       res = zero
       cres = (0d0,0d0)
       Bm =   (0d0,0d0)
       weight = zero    ! this weight is the 0 if counterevent fails to pass
                        ! jet cuts


        if( first_time ) then
              call InitTrees(2,3,6,TreeAmpsDip)
              call InitProcess_TbTGGG(ExtParticles(1:5))
              TreeAmpsDip(1)%PartRef(1:5) = (/1,2,3,4,5/)
              TreeAmpsDip(2)%PartRef(1:5) = (/1,2,3,5,4/)
              TreeAmpsDip(3)%PartRef(1:5) = (/1,2,4,3,5/)
              TreeAmpsDip(4)%PartRef(1:5) = (/1,2,4,5,3/)
              TreeAmpsDip(5)%PartRef(1:5) = (/1,2,5,3,4/)
              TreeAmpsDip(6)%PartRef(1:5) = (/1,2,5,4,3/)

              do iTree=1,6
              call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles(1:5))
              enddo
              first_time=.false.
      endif

!       momentum mapping

        pi=p(:,i)
        pj=p(:,j)
        pk=p(:,k)

        pijk = pi+pj+pk
        pij = pi+pj

        pijk2 = scr(pijk,pijk)
        pij2  = scr(pij,pij)

        mij = mj

      xx1=sqrt(vl(pijk2,mij**2,mk**2))
      xx2=sqrt(vl(pijk2,pij2,mk**2))

      pkt = xx1/xx2*(pk - scr(pijk,pk)/pijk2*pijk) &
      +(pijk2+mk**2-mij**2)/two/pijk2*pijk
      pijt = pijk - pkt


        zi = scr(pi,pk)/(scr(pi,pk)+scr(pj,pk))
        zj = 1d0-zi
        yijk = scr(pi,pj)/(scr(pi,pj)+scr(pi,pk)+scr(pj,pk))

        q2sqrt=dsqrt(pijk2)
        muk=mk/q2sqrt
        mui=mi/q2sqrt
        muj=mj/q2sqrt
        yp = 1d0-2d0*muk*(1d0-muk)/(1d0-mui**2-muj**2-muk**2)
        if( .not. invert_alphacut ) then
            if (yijk .gt. alpha_ff*yp) then
                res = 0d0
                return
            endif
        else
            if (yijk .lt. alpha_ff*yp) then
                res = 0d0
                return
            endif
        endif

        vijk  = vel(pij,pk)
        tvijk = vel(pijt,pkt)
        viji  = vel(pij,pi)

       zim = zi - one/two*(one-vijk)
       zjm = zj - one/two*(one-vijk)

       xxx = (pijk2-mi**2-mj**2-mk**2)*yijk
       xxx = (two*mi**2 + xxx )/two*(mi**2+mj**2 + xxx )
       zp = xxx*(one + viji*vijk)
       zm = xxx*(one - viji*vijk)

       paux = zim*pi - zjm*pj

          q(:,1) = p(:,1)
          fl(1) = 'qm'
          q(:,2) = p(:,2)
          fl(2) = 'qm'
          q(:,3) = p(:,3)
          fl(3) = 'gl'
          q(:,4) = p(:,4)
          fl(4) = 'gl'
          fl(5) = 'gl'

       if (i.eq.6) then
          q(:,5) = p(:,5)
          q(:,j) = pijt(:)
          q(:,k) = pkt(:)
          pos = j
       endif

       if (i.eq.5) then

          q(:,5) = p(:,6)

            if (j.eq.6) then
            q(:,5) = pijt(:)
            pos = 5
            else
            q(:,j) = pijt(:)
            pos = j
            endif

            if(k.eq.6) then
            q(:,5) = pkt(:)
             else
            q(:,k) = pkt(:)
            endif
       endif

!---- now we can check that the subtraction kinematics passes jet cuts
!    call jetktalg(q,5,pjetout,Njet,weight)
!       weight=one
!     if(weight.eq.one) then
!----- now generate the top decay products

       if (TopDecays.ge.1) then
         call EvalPhasespace_TopDecay(q(1:4,1),yRnDk(1:4),.false., &
  MomDK(1:4,1:3),PSWgt1)
         call EvalPhasespace_TopDecay(q(1:4,2),yRnDk(5:8),.false., &
  MomDK(1:4,4:6),PSWgt2)
       elseif(TopDecays.eq.0) then
         PSWgt1 = one
         PSWgt2 = one
       else
          print *, 'topdecays are not defined'
          pause
       endif


!-----------------     initial   initial     final   top      top
   call Kinematics_TTBARJET(0,(/-q(1:4,3),-q(1:4,4),q(1:4,5),q(1:4,1),q(1:4,2)/) &
   ,MomDK(1:4,1:6),Not_Passed_Cuts,NBin(1:NumHistograms))


     if(Not_Passed_Cuts.eq..false.) then

        call cc_gg_ttgggg(n,C)

!--- after momentum mapping -- sum over colors and polarizations

         Nmax = 1

        if (TopDecays.ge.1) then
            Nmax(2) = -1
           if (pos.eq.1.or.pos.eq.2) then
            Nmax(1) = -1
           elseif(pos.ge.3) then
            Nmax(pos) = -1
            endif
        endif


       do i1=-1,Nmax(1),2
         do i2 = -1,Nmax(2),2
           do i3 = -1,Nmax(3),2
             do i4 = -1,Nmax(4),2
               do i5 = -1,Nmax(5),2


           hel(1) = i1
           hel(2) = i2
           hel(3) = i3
           hel(4) = i4
           hel(5) = i5

           if (pos.eq.1) then
             hel(1) = i1
          endif

          if (pos.eq.2) then
             hel(1) = i2
             hel(2) = i1
          endif


           if (pos.eq.3) then
              hel(3) = i1
              hel(1) = i3
           endif

           if (pos.eq.4) then
              hel(4) = i1
              hel(1) = i4
           endif

           if (pos.eq.5) then
              hel(5) = i1
              hel(1) =i5
           endif

      call GenerateEvent52(q,momDK(1:4,1:6),hel,ExtParticles(1:5))


   if (pos.eq.1) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%quarks(1)%pol(1:4)
   endif


   if (pos.eq.2) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%quarks(2)%pol(1:4)
   endif


   if (pos.eq.3) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%gluons(1)%pol(1:4)
   endif

   if (pos.eq.4) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%gluons(2)%pol(1:4)
   endif

   if (pos.eq.5) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%gluons(3)%pol(1:4)
   endif


      do i6 = 1,6
      call EvalTree2(TreeAmpsDip(i6),Bm(i1,i6,i2,i3,i4,i5))
      enddo


      enddo
      enddo
      enddo
      enddo
      enddo


      do i1=-1,1,2
         do j1 = -1,1,2
            do c1 = 1,6
               do c2 = 1,6

       Am(i1,j1,c1,c2) = (0.0_dp,0.0_dp)

      do i2=-1,1,2
      do i3=-1,1,2
      do i4=-1,1,2
      do i5=-1,1,2

         Am(i1,j1,c1,c2) = Am(i1,j1,c1,c2) + &
  Bm(i1,c1,i2,i3,i4,i5)*dconjg(Bm(j1,c2,i2,i3,i4,i5))

      enddo
      enddo
      enddo
      enddo


      enddo
      enddo
      enddo
      enddo


!     now build the helicity matrix

       HH = zero

      if ( (fl1.eq.'gl').and.(fl2.eq.'gl') ) then

      diag= two*(one/(one-zi*(one-yijk)) &
      +one/(one-zj*(one-yijk)) &
       -(two - kappa_ff*zp*zm)/vijk )

       offdiag = two/vijk/scr(pi,pj)


       xm = scc(POL1(-1,:),paux)
       xp = scc(POL1(1,:),paux)

      HH(-1,1)=offdiag*dconjg(xm)*xp
      HH(1,-1)=offdiag*dconjg(xp)*xm
      HH(-1,-1)=diag + offdiag*dconjg(xm)*xm
      HH(1,1) = diag + offdiag*dconjg(xp)*xp



      endif


      if ( (fl1.eq.'gl').and.(fl2.eq.'qm') ) then

       diag = two/(one-zj*(one-yijk)) &
     -tvijk/vijk*(one + zj+mj**2/scr(pi,pj))


       HH(-1,1) = cmplx(0.0_dp,0.0_dp)
       HH(1,-1) = cmplx(0.0_dp,0.0_dp)
       HH(1,1) =  cmplx(diag,0.0_dp)
       HH(-1,-1) = cmplx(diag,0.0_dp)

      endif


        do i3=1,6
           do i4=1,6
      do i1 = -1,1,2
         do i2 = -1,1,2

      cres = cres + C(i3,i4)*Am(i1,i2,i3,i4)*HH(i1,i2)

            enddo
          enddo

           enddo
         enddo

          res = one/(pij2 -mij**2)*real(cres,dp)

!-------- account for all the weights, change the sing
          res = (-1d0)*res*PSWgt1*PSwgt2*Wgt_ext

!------- fill in histograms

         do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),res)
         enddo

         endif   ! endif for passed cuts

       end subroutine dipff


       subroutine dipfi(n,i,j,a,mi,mj,ma,fl1,fl2,p,res)
       implicit none
       integer, intent(in) :: n,i,j,a
       real(dp), intent(in) :: mi,mj,ma
       character, intent(in) :: fl1*2,fl2*2
       real(dp), intent(in) :: p(4,6)
       real(dp), intent(out) :: res
       complex(dp) :: cres
       integer :: oh(6,5)
       real(dp) :: mij
       real(dp) :: C(6,6)
       real(dp) ::  pi(4), pj(4), pk(4), pkt(4), pij(4), qq(4), qqij(4)
       real(dp) :: q(4,5), paux(4), pijk(4), pa(4),  pat(4)
       real(dp) :: xx1, xx2, zi, zj, yijk, vijk, zim
       real(dp) :: pij2,pijk2,tvijk,q2,tpij2, qij2, zjm
       real(dp) ::  pijt(4), pija(4), pija2, xija
       real(dp) :: zp, zm, xxx, viji
       integer :: i1, i2, i3, i4, i5, i6, hel(5)
       character :: fl(5)*2
       integer :: pos, iTree, c1, c2, j1
       real(dp) :: diag, offdiag
       complex(dp) :: xp, xm
       complex(dp) :: HH(-1:1,-1:1)
       type(Particle),save :: ExtParticles(1:5)
       type(TreeProcess),save :: TreeAmpsDip(1:6)
       complex(8)        :: ResAmpsDip_mt(-1:1,1:6)
       complex(8)        :: Bm(-1:1,1:6,-1:1,-1:1,-1:1,-1:1)
       complex(8)        :: Am(-1:1,-1:1,1:6,1:6)
       complex(dp) :: POL1(-1:1,4)
       real(dp) ::  pjetout(4,5), weight
       integer :: Njet
       logical, save :: first_time = .true.
       logical :: Not_Passed_Cuts
       integer :: NBin(1:NumHistograms), Nmax(5), Nhisto
       real(dp) :: MomDK(1:4,1:6),PSWgt1,PSWgt2

       res = zero
       cres = (0d0,0d0)
       Bm =   (0d0,0d0)
       weight = zero

        if( first_time ) then
                call InitTrees(2,3,6,TreeAmpsDip)
                call InitProcess_TbTGGG(ExtParticles(1:5))
                TreeAmpsDip(1)%PartRef(1:5) = (/1,2,3,4,5/)
                TreeAmpsDip(2)%PartRef(1:5) = (/1,2,3,5,4/)
                TreeAmpsDip(3)%PartRef(1:5) = (/1,2,4,3,5/)
                TreeAmpsDip(4)%PartRef(1:5) = (/1,2,4,5,3/)
                TreeAmpsDip(5)%PartRef(1:5) = (/1,2,5,3,4/)
                TreeAmpsDip(6)%PartRef(1:5) = (/1,2,5,4,3/)
                do iTree=1,6
                    call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles(1:5))
                enddo
                first_time = .false.
        endif

!       momentum mapping
        pi=p(:,i)
        pj=p(:,j)
        pa= -p(:,a)  ! the `-' sign accounts for the all outgoing convention

        pija = pi+pj-pa
        pij = pi+pj

        pija2 = scr(pija,pija)
        pij2  = scr(pij,pij)

        mij = mj

        xxx = scr(pa,pi)+scr(pa,pj)
        xija = one -(scr(pi,pj)+half*(mij**2-mi**2-mj**2))/xxx
        zi = scr(pa,pi)/xxx
        zj = scr(pa,pj)/xxx

        pat = xija*pa
        pijt = pi + pj - (one - xija)*pa
        if( .not. invert_alphacut ) then
            if( alpha_fi.lt.one-xija ) then
                res = 0d0
                return
            endif
        else
            if( alpha_fi.gt.one-xija ) then
                res = 0d0
                return
            endif
        endif

          fl(1) = 'qm'
          fl(2) = 'qm'
          fl(3) = 'gl'
          fl(4) = 'gl'
          fl(5) = 'gl'

       if (i.eq.6) then
          q(:,1) = p(:,1)
          q(:,2) = p(:,2)
          q(:,3) = p(:,3)
          q(:,4) = p(:,4)
          q(:,5) = p(:,5)

          q(:,j) = pijt(:)
          q(:,a) = -pat(:)
          pos = j
       endif

       if (i.eq.5) then

          q(:,1) = p(:,1)
          q(:,2) = p(:,2)
          q(:,3) = p(:,3)
          q(:,4) = p(:,4)
          q(:,5) = p(:,6)

          q(:,a) = -pat(:)

            if (j.eq.6) then
            q(:,5) = pijt(:)
            pos = 5
            else
            q(:,j) = pijt(:)
            pos = j
            endif
       endif


!    call jetktalg(q,5,pjetout,Njet,weight)
!       weight=one
!     if (weight.eq.one) then


!----- now generate the top decay products

       if (TopDecays.ge.1) then
         call EvalPhasespace_TopDecay(q(1:4,1),yRnDk(1:4),.false., &
  MomDK(1:4,1:3),PSWgt1)
         call EvalPhasespace_TopDecay(q(1:4,2),yRnDk(5:8),.false., &
  MomDK(1:4,4:6),PSWgt2)
       elseif(TopDecays.eq.0) then
         PSWgt1 = one
         PSWgt2 = one
       else
          print *, 'topdecays are not defined'
          pause
       endif


!-----------------     initial   initial       final  top      top
   call Kinematics_TTBARJET(0,(/-q(1:4,3),-q(1:4,4),q(1:4,5),q(1:4,1),q(1:4,2)/) &
   ,MomDK(1:4,1:6),Not_Passed_Cuts,NBin(1:NumHistograms))

     if(Not_Passed_Cuts.eq..false.) then

        call cc_gg_ttgggg(n,C)

         Nmax = 1

        if (TopDecays.ge.1) then
            Nmax(2) = -1
           if (pos.eq.1.or.pos.eq.2) then
            Nmax(1) = -1
           elseif(pos.ge.3) then
            Nmax(pos) = -1
            endif
        endif


!--- after momentum mapping -- sum over colors and polarizations

       do i1=-1,Nmax(1),2
          do i2 = -1,Nmax(2),2
             do i3 = -1,Nmax(3),2
               do i4 = -1,Nmax(4),2
                 do i5=-1,Nmax(5),2

           hel(1) = i1
           hel(2) = i2
           hel(3) = i3
           hel(4) = i4
           hel(5) = i5

           if (pos.eq.1) then
             hel(1) = i1
          endif

          if (pos.eq.2) then
             hel(1) = i2
             hel(2) = i1
          endif

           if (pos.eq.3) then
              hel(3) = i1
              hel(1) = i3
           endif

           if (pos.eq.4) then
              hel(4) = i1
              hel(1) = i4
           endif

           if (pos.eq.5) then
              hel(5) = i1
              hel(1) =i5
           endif

      call GenerateEvent52(q,momDK(1:4,1:6),hel,ExtParticles(1:5))

      if (pos.eq.1) then
       if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(1)%pOL(1:4)
       if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(1)%pOL(1:4)
      endif


      if (pos.eq.2) then
       if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(2)%pOL(1:4)
       if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(2)%pOL(1:4)
      endif


      if (pos.eq.3) then
       if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%gLUONS(1)%pOL(1:4)
       if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%gLUONS(1)%pOL(1:4)
      endif

      if (pos.eq.4) then
       if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%gLUONS(2)%pOL(1:4)
       if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%gLUONS(2)%pOL(1:4)
      endif

      if (pos.eq.5) then
       if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%gLUONS(3)%pOL(1:4)
       if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%gLUONS(3)%pOL(1:4)
      endif

      do i6 = 1,6
      call EvalTree2(TreeAmpsDip(i6),Bm(i1,i6,i2,i3,i4,i5))
      enddo

      enddo
      enddo
      enddo
      enddo
      enddo

      do i1=-1,1,2
         do j1 = -1,1,2
            do c1 = 1,6
               do c2 = 1,6

       Am(i1,j1,c1,c2) = (0.0_dp,0.0_dp)

      do i2=-1,1,2
      do i3=-1,1,2
      do i4=-1,1,2
      do i5=-1,1,2

         Am(i1,j1,c1,c2) = Am(i1,j1,c1,c2) + &
  Bm(i1,c1,i2,i3,i4,i5)*dconjg(Bm(j1,c2,i2,i3,i4,i5))

      enddo
      enddo
      enddo
      enddo


      enddo
      enddo
      enddo
      enddo

!     now build the helicity matrix

       HH = zero

      if ( (fl1.eq.'gl').and.(fl2.eq.'gl') ) then

       paux = zi*pi - zj*pj

      diag= two*(one/(one-zi + one - xija)+one/(one-zj+ one-xija) - two)
      offdiag = two/scr(pi,pj)

       xm = scc(POL1(-1,:),paux)
       xp = scc(POL1(1,:),paux)

      HH(-1,1)=offdiag*dconjg(xm)*xp
      HH(1,-1)=offdiag*dconjg(xp)*xm
      HH(-1,-1)=diag + offdiag*dconjg(xm)*xm
      HH(1,1) = diag + offdiag*dconjg(xp)*xp


      endif


      if ( (fl1.eq.'gl').and.(fl2.eq.'qm') ) then

    diag = two/(two - xija - zj) - one - zj -mj**2/scr(pi,pj)

       HH(-1,1) = cmplx(0.0_dp,0.0_dp)
       HH(1,-1) = cmplx(0.0_dp,0.0_dp)
       HH(1,1) =  cmplx(diag,0.0_dp)
       HH(-1,-1) = cmplx(diag,0.0_dp)
      endif


      do i1 = -1,1,2
         do i2 = -1,1,2

        do i3=1,6
           do i4=1,6


!      cres = cres + C(i3,i4)*Am(i1,i3)*dconjg(Am(i2,i4))*HH(i1,i2)
      cres = cres + C(i3,i4)*Am(i1,i2,i3,i4)*HH(i1,i2)


            enddo
          enddo
         enddo
         enddo

         res = one/xija/(pij2-mij**2)*real(cres,dp)


!-------- account for all the weights, change the sing
         res = (-1d0)*res*PSWgt1*PSwgt2*Wgt_ext


!------- fill in histograms

         do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),res)
         enddo


         endif ! for passed cuts

        end subroutine dipfi


       subroutine dipif(n,i,a,j,mi,ma,mj,fl1,fl2,p,res)
       implicit none
       integer, intent(in) :: n,i,j,a
       real(dp), intent(in) :: mi,mj,ma
       character, intent(in) :: fl1*2,fl2*2
       real(dp), intent(in) :: p(4,6)
       real(dp), intent(out) :: res
       complex(dp) :: cres
       real(dp) :: mij
       real(dp) :: C(6,6)
       real(dp) ::  pi(4), pj(4), pk(4), pkt(4), pij(4), qq(4), qqij(4)
       real(dp) :: q(4,5), paux(4), pijk(4), pa(4),  pat(4)
       real(dp) :: xx1, xx2, zi, zj, yijk, vijk, zim, u
       real(dp) :: pij2,pijk2,tvijk,q2,tpij2, qij2, zjm
       real(dp) ::  pijt(4), pija(4), pija2, xija
       real(dp) :: zp, zm, xxx, viji, pait(4), pjt(4)
       integer :: i1, i2, i3, i4, i5, i6, hel(5)
       character :: fl(5)*2
       integer :: pos, iTree
       real(dp) :: diag, offdiag
       complex(dp) :: xp, xm
       complex(dp) :: HH(-1:1,-1:1)
       type(Particle), save :: ExtParticles(1:5)
       type(TreeProcess), save :: TreeAmpsDip(1:6)
       complex(8)        :: ResAmpsDip_mt(-1:1,1:6)
       integer :: c1, c2, j1
       complex(8)        :: Bm(-1:1,1:6,-1:1,-1:1,-1:1,-1:1)
       complex(8)        :: Am(-1:1,-1:1,1:6,1:6)
       complex(dp) :: POL1(-1:1,4)
       real(dp) ::  pjetout(4,5), weight
       integer :: Njet
       logical, save :: first_time = .true.
       logical :: Not_Passed_Cuts
       integer :: NBin(1:NumHistograms), Nmax(5), Nhisto
       real(dp) :: MomDK(1:4,1:6),PSWgt1,PSWgt2


       res = zero
       cres = (0d0,0d0)
       Bm = (0d0,0d0)
       weight = zero


        if( first_time ) then
              call InitTrees(2,3,6,TreeAmpsDip)
              call InitProcess_TbTGGG(ExtParticles(1:5))
              TreeAmpsDip(1)%PartRef(1:5) = (/1,2,3,4,5/)
              TreeAmpsDip(2)%PartRef(1:5) = (/1,2,3,5,4/)
              TreeAmpsDip(3)%PartRef(1:5) = (/1,2,4,3,5/)
              TreeAmpsDip(4)%PartRef(1:5) = (/1,2,4,5,3/)
              TreeAmpsDip(5)%PartRef(1:5) = (/1,2,5,3,4/)
              TreeAmpsDip(6)%PartRef(1:5) = (/1,2,5,4,3/)
              do iTree=1,6
                  call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles(1:5))
              enddo
              first_time = .false.
      endif


!       momentum mapping
        pi =  p(:,i)
        pa = -p(:,a) ! the `-' sign accounts for the all outgoing convention
        pj =  p(:,j)

        pija = pi+pj-pa
        pij  = pi+pj

        pija2 = scr(pija,pija)
        pij2  = scr(pij,pij)

        mij = mj

        xxx = scr(pa,pi)+scr(pa,pj)
        xija = one -scr(pi,pj)/xxx
        zi = scr(pa,pi)/xxx
        zj = scr(pa,pj)/xxx

        u = scr(pa,pi)/xxx
        if( .not. invert_alphacut ) then
            if( alpha_if.lt.u ) then
                res = 0d0
                return
            endif
        else
            if( alpha_if.gt.u ) then
                res = 0d0
                return
            endif
        endif

        pait = xija*pa
        pjt =  pi + pj - (one - xija)*pa


          fl(1) = 'qm'
          fl(2) = 'qm'
          fl(3) = 'gl'
          fl(4) = 'gl'
          fl(5) = 'gl'


       if (i.eq.6) then
          q(:,1) = p(:,1)
          q(:,2) = p(:,2)
          q(:,3) = p(:,3)
          q(:,4) = p(:,4)
          q(:,5) = p(:,5)

          q(:,a) = -pait(:)

          q(:,j) = pjt(:)
          pos = a
       endif


       if (i.eq.5) then

          q(:,1) = p(:,1)
          q(:,2) = p(:,2)
          q(:,3) = p(:,3)
          q(:,4) = p(:,4)
          q(:,5) = p(:,6)

          if (j.eq.6) then
          q(:,5) = pjt(:)
          else
          q(:,j) = pjt(:)
          endif


          q(:,a) = -pait(:)

          pos = a

        endif

!    call jetktalg(q,5,pjetout,Njet,weight)
!     if(weight.eq.one) then


       if (TopDecays.ge.1) then
         call EvalPhasespace_TopDecay(q(1:4,1),yRnDk(1:4),.false., &
  MomDK(1:4,1:3),PSWgt1)
         call EvalPhasespace_TopDecay(q(1:4,2),yRnDk(5:8),.false., &
  MomDK(1:4,4:6),PSWgt2)
       elseif(TopDecays.eq.0) then
         PSWgt1 = one
         PSWgt2 = one
       else
          print *, 'topdecays are not defined'
          pause
       endif

!----------------- initial   initial        final top      top
   call Kinematics_TTBARJET(0,(/-q(1:4,3),-q(1:4,4),q(1:4,5),q(1:4,1),q(1:4,2)/) &
   ,MomDK(1:4,1:6),Not_Passed_Cuts,NBin(1:NumHistograms))


     if(Not_Passed_Cuts.eq..false.) then

        call cc_gg_ttgggg(n,C)


!--- after momentum mapping -- sum over colors and polarizations

         Nmax = 1

        if (TopDecays.ge.1) then
            Nmax(2) = -1
           if (pos.eq.1.or.pos.eq.2) then
            Nmax(1) = -1
           elseif(pos.ge.3) then
            Nmax(pos) = -1
            endif
        endif

       do i1 = -1,Nmax(1),2
          do i2 = -1,Nmax(2),2
             do i3 = -1,Nmax(3),2
                do i4 = -1,Nmax(4),2
                   do i5 = -1,Nmax(5),2


           hel(1) = i1
           hel(2) = i2
           hel(3) = i3
           hel(4) = i4
           hel(5) = i5

           if (pos.eq.1) then
             hel(1) = i1
          endif

          if (pos.eq.2) then
             hel(1) = i2
             hel(2) = i1
          endif


           if (pos.eq.3) then
              hel(3) = i1
              hel(1) = i3
           endif

           if (pos.eq.4) then
              hel(4) = i1
              hel(1) = i4
           endif

           if (pos.eq.5) then
              hel(5) = i1
              hel(1) =i5
           endif

      call GenerateEvent52(q,momDK(1:4,1:6),hel,ExtParticles(1:5))

!       print *, "if dipole"
!       print *, ExtParticles(1)%Mom(1:4).dot.ExtParticles(1)%Mom(1:4)
!       print *, ExtParticles(2)%Mom(1:4).dot.ExtParticles(2)%Mom(1:4)
!       print *, ExtParticles(3)%Mom(1:4).dot.ExtParticles(3)%Mom(1:4)
!       print *, ExtParticles(4)%Mom(1:4).dot.ExtParticles(4)%Mom(1:4)
!       print *, ExtParticles(5)%Mom(1:4).dot.ExtParticles(5)%Mom(1:4)



       if (pos.eq.1) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(1)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(1)%pOL(1:4)
      endif


       if (pos.eq.2) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(2)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(2)%pOL(1:4)
      endif


       if (pos.eq.3) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%gLUONS(1)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%gLUONS(1)%pOL(1:4)
      endif

      if (pos.eq.4) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%gLUONS(2)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%gLUONS(2)%pOL(1:4)
      endif

      if (pos.eq.5) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%gLUONS(3)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%gLUONS(3)%pOL(1:4)
      endif

      do i6 = 1,6
      call EvalTree2(TreeAmpsDip(i6),Bm(i1,i6,i2,i3,i4,i5))
      enddo

      enddo
      enddo
      enddo
      enddo
      enddo


      do i1=-1,1,2
         do j1 = -1,1,2
            do c1 = 1,6
               do c2 = 1,6

       Am(i1,j1,c1,c2) = (0.0_dp,0.0_dp)

      do i2=-1,1,2
      do i3=-1,1,2
      do i4=-1,1,2
      do i5=-1,1,2

         Am(i1,j1,c1,c2) = Am(i1,j1,c1,c2) + &
  Bm(i1,c1,i2,i3,i4,i5)*dconjg(Bm(j1,c2,i2,i3,i4,i5))

      enddo
      enddo
      enddo
      enddo


      enddo
      enddo
      enddo
      enddo



!     now build the helicity matrix

       HH = zero

      if ( (fl1.eq.'gl').and.(fl2.eq.'gl') ) then

       paux = pi/zi - pj/zj

      diag= two*(one/(two-xija - zj)-one + xija*(one-xija))
      offdiag = two*(one-xija)/xija*zi*zj/scr(pi,pj)

       xm = scc(POL1(-1,:),paux)
       xp = scc(POL1(1,:),paux)

      HH(-1,1)=offdiag*dconjg(xm)*xp
      HH(1,-1)=offdiag*dconjg(xp)*xm
      HH(-1,-1)=diag + offdiag*dconjg(xm)*xm
      HH(1,1) = diag + offdiag*dconjg(xp)*xp


      endif




      do i1 = -1,1,2
         do i2 = -1,1,2

        do i3=1,6
           do i4=1,6

      cres = cres + C(i3,i4)*Am(i1,i2,i3,i4)*HH(i1,i2)
!      cres = cres + C(i3,i4)*Am(i1,i3)*dconjg(Am(i2,i4))*HH(i1,i2)

            enddo
           enddo
          enddo
         enddo

         res = one/xija/two/scr(pa,pi)*real(cres,dp)

!-------- account for all the weights, change the sing
         res = (-1d0)*res*PSWgt1*PSwgt2*Wgt_ext

!------- fill in histograms

         do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),res)
         enddo

         endif ! -- ! passed cuts

        end subroutine dipif



       subroutine dipii(n,i,a,b,mi,ma,mb,fl1,fl2,p,res)
       implicit none
       integer, intent(in) :: n,i,a,b
       real(dp), intent(in) :: mi,ma,mb
       character, intent(in) :: fl1*2,fl2*2
       real(dp), intent(in) :: p(4,6)
       real(dp), intent(out) :: res
       complex(dp) :: cres
       real(dp) :: mij
       real(dp) :: C(6,6)
       real(dp) ::  pi(4), pj(4), pk(4), pkt(4), pij(4), qq(4), qqij(4)
       real(dp) :: q(4,5), paux(4), pijk(4), pa(4),  pat(4)
       real(dp) :: xx1, xx2, zi, zj, yijk, vijk, zim
       real(dp) :: pij2,pijk2,tvijk,q2,tpij2, qij2, zjm
       real(dp) ::  pijt(4), pija(4), pb(4),pija2, xija
       real(dp) :: pait(4), xiab
       real(dp) ::  K(4), tK(4),K2,tK2, K2tK2, k1(4)
       real(dp) :: zp, zm, xxx, viji,v
       integer :: i1, i2, i3, i4, i5, i6, hel(5), j
       character :: fl(5)*2
       integer :: pos, iTree
       real(dp) :: diag, offdiag
       complex(dp) :: xp, xm
       complex(dp) :: HH(-1:1,-1:1)
       type(Particle),save :: ExtParticles(1:5)
       type(TreeProcess),save :: TreeAmpsDip(1:6)
       integer :: c1, c2, j1
       complex(8)        :: Bm(-1:1,1:6,-1:1,-1:1,-1:1,-1:1)
       complex(8)        :: Am(-1:1,-1:1,1:6,1:6)
       complex(8)        :: ResAmpsDip_mt(-1:1,1:6)
       complex(dp) :: POL1(-1:1,4)
       real(dp) ::  pjetout(4,5), weight
       integer :: Njet
       logical, save :: first_time = .true.
       logical :: Not_Passed_Cuts
       integer :: NBin(1:NumHistograms), Nmax(5), Nhisto
       real(dp) :: MomDK(1:4,1:6),PSWgt1,PSWgt2


       res = zero
       cres = (0d0,0d0)
       Bm = (0d0,0d0)

        if( first_time ) then
              call InitTrees(2,3,6,TreeAmpsDip)
              call InitProcess_TbTGGG(ExtParticles(1:5))
              TreeAmpsDip(1)%PartRef(1:5) = (/1,2,3,4,5/)
              TreeAmpsDip(2)%PartRef(1:5) = (/1,2,3,5,4/)
              TreeAmpsDip(3)%PartRef(1:5) = (/1,2,4,3,5/)
              TreeAmpsDip(4)%PartRef(1:5) = (/1,2,4,5,3/)
              TreeAmpsDip(5)%PartRef(1:5) = (/1,2,5,3,4/)
              TreeAmpsDip(6)%PartRef(1:5) = (/1,2,5,4,3/)
              do iTree=1,6
                  call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles(1:5))
              enddo
              first_time=.false.
      endif

!       momentum mapping
        pi =  p(:,i)
        pa = -p(:,a) ! the `-' sign accounts for the all outgoing convention
        pb = -p(:,b) ! the `-' sing accounts for the all outgoing convention

        xiab = one - (scr(pi,pa)+scr(pi,pb))/scr(pa,pb)
        pait = xiab*pa

        v = scr(pi,pa)/scr(pa,pb)
        if( .not. invert_alphacut ) then
            if( alpha_ii.lt.v ) then
                res = 0d0
                return
            endif
        else
            if( alpha_ii.gt.v ) then
                res = 0d0
                return
            endif
        endif

          fl(1) = 'qm'
          fl(2) = 'qm'
          fl(3) = 'gl'
          fl(4) = 'gl'
          fl(5) = 'gl'

       if (i.eq.6) then

          q(:,1) = p(:,1)
          q(:,2) = p(:,2)
          q(:,3) = p(:,3)
          q(:,4) = p(:,4)
          q(:,5) = p(:,5)

          q(:,a) = -pait(:)
          q(:,b) = -pb(:)
          pos = a
       endif

       if (i.eq.5) then

          q(:,1) = p(:,1)
          q(:,2) = p(:,2)
          q(:,3) = p(:,3)
          q(:,4) = p(:,4)
          q(:,5) = p(:,6)

          q(:,a) = -pait(:)
          q(:,b) = -pb(:)
          pos = a

        endif


!-- now Lorentz transform

        K = pa + pb - pi
        tK = pait + pb
        K2 = scr(K,K)
        tK2 = scr(tK,tK)
        K2tK2 = scr(tK,K)
        xxx = K2+two*K2tK2 + tK2


        do j=1,5
     if (j.ne.3.and.j.ne.4) then
     k1 = q(:,j)
     k1 = k1 - two*(scr(k1,K)+scr(k1,tK))/xxx*(K+tK)+two*scr(k1,K)/K2*tK
     q(:,j) = k1
     endif
        enddo

!    call jetktalg(q,5,pjetout,Njet,weight)
!       weight=one
!     if(weight.eq.one) then


       if (TopDecays.ge.1) then
         call EvalPhasespace_TopDecay(q(1:4,1),yRnDk(1:4),.false., &
  MomDK(1:4,1:3),PSWgt1)
         call EvalPhasespace_TopDecay(q(1:4,2),yRnDk(5:8),.false., &
  MomDK(1:4,4:6),PSWgt2)
       elseif(TopDecays.eq.0) then
         PSWgt1 = one
         PSWgt2 = one
       else
          print *, 'topdecays are not defined'
          pause
       endif

!-----------------     initial   initial       final  top      top
   call Kinematics_TTBARJET(0,(/-q(1:4,3),-q(1:4,4),q(1:4,5),q(1:4,1),q(1:4,2)/) &
   ,MomDK(1:4,1:6),Not_Passed_Cuts,NBin(1:NumHistograms))


     if(Not_Passed_Cuts.eq..false.) then

        call cc_gg_ttgggg(n,C)


!--- after momentum mapping -- sum over colors and polarizations

         Nmax = 1

        if (TopDecays.ge.1) then
            Nmax(2) = -1
           if (pos.eq.1.or.pos.eq.2) then
            Nmax(1) = -1
           elseif(pos.ge.3) then
            Nmax(pos) = -1
            endif
        endif



       do i1=-1,Nmax(1),2
          do i2 = -1,Nmax(2),2
             do i3 = -1,Nmax(3),2
                do i4 = -1,Nmax(4),2
                   do i5=-1,Nmax(5),2

           hel(1) = i1
           hel(2) = i2
           hel(3) = i3
           hel(4) = i4
           hel(5) = i5

           if (pos.eq.1) then
             hel(1) = i1
          endif

          if (pos.eq.2) then
             hel(1) = i2
             hel(2) = i1
          endif


           if (pos.eq.3) then
              hel(3) = i1
              hel(1) = i3
           endif

           if (pos.eq.4) then
              hel(4) = i1
              hel(1) = i4
           endif

           if (pos.eq.5) then
              hel(5) = i1
              hel(1) =i5
           endif

      call GenerateEvent52(q,momDK(1:4,1:6),hel,ExtParticles(1:5))

       if (pos.eq.1) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(1)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(1)%pOL(1:4)
      endif


       if (pos.eq.2) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(2)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(2)%pOL(1:4)
      endif


       if (pos.eq.3) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%gLUONS(1)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%gLUONS(1)%pOL(1:4)
      endif

      if (pos.eq.4) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%gLUONS(2)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%gLUONS(2)%pOL(1:4)
      endif

      if (pos.eq.5) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%gLUONS(3)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%gLUONS(3)%pOL(1:4)
      endif


      do i6 = 1,6
      call EvalTree2(TreeAmpsDip(i6),Bm(i1,i6,i2,i3,i4,i5))
      enddo

      enddo
      enddo
      enddo
      enddo
      enddo


      do i1=-1,1,2
         do j1 = -1,1,2
            do c1 = 1,6
               do c2 = 1,6

       Am(i1,j1,c1,c2) = (0.0_dp,0.0_dp)

      do i2=-1,1,2
      do i3=-1,1,2
      do i4=-1,1,2
      do i5=-1,1,2

         Am(i1,j1,c1,c2) = Am(i1,j1,c1,c2) + &
  Bm(i1,c1,i2,i3,i4,i5)*dconjg(Bm(j1,c2,i2,i3,i4,i5))

      enddo
      enddo
      enddo
      enddo


      enddo
      enddo
      enddo
      enddo


!     now build the helicity matrix

       HH = zero

      if ( (fl1.eq.'gl').and.(fl2.eq.'gl') ) then

       paux = pi - scr(pi,pa)/scr(pi,pb)*pb

!       diag= two*xiab/(one-xiab)+xiab*(one-xiab)   ! HERE WAS A BUG, compare to line below
      diag= two*(xiab/(one-xiab)+xiab*(one-xiab))
      offdiag = two*(one-xiab)/xiab*scr(pa,pb)/scr(pi,pa)/scr(pi,pb)

       xm = scc(POL1(-1,:),paux)
       xp = scc(POL1(1,:),paux)

      HH(-1,1)=offdiag*dconjg(xm)*xp
      HH(1,-1)=offdiag*dconjg(xp)*xm
      HH(-1,-1)=diag + offdiag*dconjg(xm)*xm
      HH(1,1) = diag + offdiag*dconjg(xp)*xp

      endif




      do i1 = -1,1,2
         do i2 = -1,1,2

        do i3=1,6
           do i4=1,6


!      cres = cres + C(i3,i4)*Am(i1,i3)*dconjg(Am(i2,i4))*HH(i1,i2)
      cres = cres + C(i3,i4)*Am(i1,i2,i3,i4)*HH(i1,i2)

            enddo
          enddo
         enddo
         enddo


         res = one/xiab/two/scr(pa,pi)*real(cres,dp)


!-------- account for all the weights, change the sing
         res = (-1d0)*res*PSWgt1*PSwgt2*Wgt_ext


!------- fill in histograms

         do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),res)
         enddo

         endif !  for passed cuts

        end subroutine dipii










!       subroutine cc_gg_ttgggg(n,C)
!       integer, intent(in)  :: n
!       real(dp), intent(out)  :: C(6,6)
!
!       if(n.eq.1) then
!       C(1,1)=8.0_dp/3.0_dp;
!       C(1,2)=80.0_dp/3.0_dp;
!       C(1,3)=8.0_dp/3.0_dp;
!       C(1,4)=8.0_dp/3.0_dp;
!       C(1,5)=-64.0_dp/3.0_dp;
!       C(1,6)=80.0_dp/3.0_dp;
!       C(2,1)=80.0_dp/3.0_dp;
!       C(2,2)=8.0_dp/3.0_dp;
!       C(2,3)=-64.0_dp/3.0_dp;
!       C(2,4)=80.0_dp/3.0_dp;
!       C(2,5)=8.0_dp/3.0_dp;
!       C(2,6)=8.0_dp/3.0_dp;
!       C(3,1)=8.0_dp/3.0_dp;
!       C(3,2)=-64.0_dp/3.0_dp;
!       C(3,3)=-64.0_dp/3.0_dp;
!       C(3,4)=-64.0_dp/3.0_dp;
!       C(3,5)=-64.0_dp/3.0_dp;
!       C(3,6)=8.0_dp/3.0_dp;
!       C(4,1)=8.0_dp/3.0_dp;
!       C(4,2)=80.0_dp/3.0_dp;
!       C(4,3)=-64.0_dp/3.0_dp;
!       C(4,4)=512.0_dp/3.0_dp;
!       C(4,5)=8.0_dp/3.0_dp;
!       C(4,6)=-64.0_dp/3.0_dp;
!       C(5,1)=-64.0_dp/3.0_dp;
!       C(5,2)=8.0_dp/3.0_dp;
!       C(5,3)=-64.0_dp/3.0_dp;
!       C(5,4)=8.0_dp/3.0_dp;
!       C(5,5)=-64.0_dp/3.0_dp;
!       C(5,6)=-64.0_dp/3.0_dp;
!       C(6,1)=80.0_dp/3.0_dp;
!       C(6,2)=8.0_dp/3.0_dp;
!       C(6,3)=8.0_dp/3.0_dp;
!       C(6,4)=-64.0_dp/3.0_dp;
!       C(6,5)=-64.0_dp/3.0_dp;
!       C(6,6)=512.0_dp/3.0_dp;
!       endif
!       if(n.eq.2) then
!       C(1,1)=512.0_dp/3.0_dp;
!       C(1,2)=-64.0_dp/3.0_dp;
!       C(1,3)=-64.0_dp/3.0_dp;
!       C(1,4)=8.0_dp/3.0_dp;
!       C(1,5)=8.0_dp/3.0_dp;
!       C(1,6)=80.0_dp/3.0_dp;
!       C(2,1)=-64.0_dp/3.0_dp;
!       C(2,2)=512.0_dp/3.0_dp;
!       C(2,3)=8.0_dp/3.0_dp;
!       C(2,4)=80.0_dp/3.0_dp;
!       C(2,5)=-64.0_dp/3.0_dp;
!       C(2,6)=8.0_dp/3.0_dp;
!       C(3,1)=-64.0_dp/3.0_dp;
!       C(3,2)=8.0_dp/3.0_dp;
!       C(3,3)=-64.0_dp/3.0_dp;
!       C(3,4)=8.0_dp/3.0_dp;
!       C(3,5)=-64.0_dp/3.0_dp;
!       C(3,6)=-64.0_dp/3.0_dp;
!       C(4,1)=8.0_dp/3.0_dp;
!       C(4,2)=80.0_dp/3.0_dp;
!       C(4,3)=8.0_dp/3.0_dp;
!       C(4,4)=8.0_dp/3.0_dp;
!       C(4,5)=-64.0_dp/3.0_dp;
!       C(4,6)=80.0_dp/3.0_dp;
!       C(5,1)=8.0_dp/3.0_dp;
!       C(5,2)=-64.0_dp/3.0_dp;
!       C(5,3)=-64.0_dp/3.0_dp;
!       C(5,4)=-64.0_dp/3.0_dp;
!       C(5,5)=-64.0_dp/3.0_dp;
!       C(5,6)=8.0_dp/3.0_dp;
!       C(6,1)=80.0_dp/3.0_dp;
!       C(6,2)=8.0_dp/3.0_dp;
!       C(6,3)=-64.0_dp/3.0_dp;
!       C(6,4)=80.0_dp/3.0_dp;
!       C(6,5)=8.0_dp/3.0_dp;
!       C(6,6)=8.0_dp/3.0_dp;
!       endif
!       if(n.eq.3) then
!       C(1,1)=192.0_dp;
!       C(1,2)=-24.0_dp;
!       C(1,3)=0.0_dp;
!       C(1,4)=0.0_dp;
!       C(1,5)=24.0_dp;
!       C(1,6)=48.0_dp;
!       C(2,1)=-24.0_dp;
!       C(2,2)=-24.0_dp;
!       C(2,3)=0.0_dp;
!       C(2,4)=-48.0_dp;
!       C(2,5)=-24.0_dp;
!       C(2,6)=0.0_dp;
!       C(3,1)=0.0_dp;
!       C(3,2)=0.0_dp;
!       C(3,3)=192.0_dp;
!       C(3,4)=-24.0_dp;
!       C(3,5)=48.0_dp;
!       C(3,6)=24.0_dp;
!       C(4,1)=0.0_dp;
!       C(4,2)=-48.0_dp;
!       C(4,3)=-24.0_dp;
!       C(4,4)=-24.0_dp;
!       C(4,5)=0.0_dp;
!       C(4,6)=-24.0_dp;
!       C(5,1)=24.0_dp;
!       C(5,2)=-24.0_dp;
!       C(5,3)=48.0_dp;
!       C(5,4)=0.0_dp;
!       C(5,5)=192.0_dp;
!       C(5,6)=0.0_dp;
!       C(6,1)=48.0_dp;
!       C(6,2)=0.0_dp;
!       C(6,3)=24.0_dp;
!       C(6,4)=-24.0_dp;
!       C(6,5)=0.0_dp;
!       C(6,6)=192.0_dp;
!       endif
!       if(n.eq.4) then
!       C(1,1)=-24.0_dp;
!       C(1,2)=-24.0_dp;
!       C(1,3)=-24.0_dp;
!       C(1,4)=0.0_dp;
!       C(1,5)=0.0_dp;
!       C(1,6)=-48.0_dp;
!       C(2,1)=-24.0_dp;
!       C(2,2)=192.0_dp;
!       C(2,3)=24.0_dp;
!       C(2,4)=48.0_dp;
!       C(2,5)=0.0_dp;
!       C(2,6)=0.0_dp;
!       C(3,1)=-24.0_dp;
!       C(3,2)=24.0_dp;
!       C(3,3)=192.0_dp;
!       C(3,4)=0.0_dp;
!       C(3,5)=48.0_dp;
!       C(3,6)=0.0_dp;
!       C(4,1)=0.0_dp;
!       C(4,2)=48.0_dp;
!       C(4,3)=0.0_dp;
!       C(4,4)=192.0_dp;
!       C(4,5)=24.0_dp;
!       C(4,6)=-24.0_dp;
!       C(5,1)=0.0_dp;
!       C(5,2)=0.0_dp;
!       C(5,3)=48.0_dp;
!       C(5,4)=24.0_dp;
!       C(5,5)=192.0_dp;
!       C(5,6)=-24.0_dp;
!       C(6,1)=-48.0_dp;
!       C(6,2)=0.0_dp;
!       C(6,3)=0.0_dp;
!       C(6,4)=-24.0_dp;
!       C(6,5)=-24.0_dp;
!       C(6,6)=-24.0_dp;
!       endif
!       if(n.eq.5) then
!       C(1,1)=-64.0_dp/3.0_dp;
!       C(1,2)=-64.0_dp/3.0_dp;
!       C(1,3)=8.0_dp/3.0_dp;
!       C(1,4)=-64.0_dp/3.0_dp;
!       C(1,5)=8.0_dp/3.0_dp;
!       C(1,6)=-64.0_dp/3.0_dp;
!       C(2,1)=-64.0_dp/3.0_dp;
!       C(2,2)=512.0_dp/3.0_dp;
!       C(2,3)=8.0_dp/3.0_dp;
!       C(2,4)=80.0_dp/3.0_dp;
!       C(2,5)=-64.0_dp/3.0_dp;
!       C(2,6)=8.0_dp/3.0_dp;
!       C(3,1)=8.0_dp/3.0_dp;
!       C(3,2)=8.0_dp/3.0_dp;
!       C(3,3)=8.0_dp/3.0_dp;
!       C(3,4)=80.0_dp/3.0_dp;
!       C(3,5)=80.0_dp/3.0_dp;
!       C(3,6)=-64.0_dp/3.0_dp;
!       C(4,1)=-64.0_dp/3.0_dp;
!       C(4,2)=80.0_dp/3.0_dp;
!       C(4,3)=80.0_dp/3.0_dp;
!       C(4,4)=8.0_dp/3.0_dp;
!       C(4,5)=8.0_dp/3.0_dp;
!       C(4,6)=8.0_dp/3.0_dp;
!       C(5,1)=8.0_dp/3.0_dp;
!       C(5,2)=-64.0_dp/3.0_dp;
!       C(5,3)=80.0_dp/3.0_dp;
!       C(5,4)=8.0_dp/3.0_dp;
!       C(5,5)=512.0_dp/3.0_dp;
!       C(5,6)=-64.0_dp/3.0_dp;
!       C(6,1)=-64.0_dp/3.0_dp;
!       C(6,2)=8.0_dp/3.0_dp;
!       C(6,3)=-64.0_dp/3.0_dp;
!       C(6,4)=8.0_dp/3.0_dp;
!       C(6,5)=-64.0_dp/3.0_dp;
!       C(6,6)=-64.0_dp/3.0_dp;
!       endif
!       if(n.eq.6) then
!       C(1,1)=-64.0_dp/3.0_dp;
!       C(1,2)=8.0_dp/3.0_dp;
!       C(1,3)=-64.0_dp/3.0_dp;
!       C(1,4)=8.0_dp/3.0_dp;
!       C(1,5)=-64.0_dp/3.0_dp;
!       C(1,6)=-64.0_dp/3.0_dp;
!       C(2,1)=8.0_dp/3.0_dp;
!       C(2,2)=8.0_dp/3.0_dp;
!       C(2,3)=8.0_dp/3.0_dp;
!       C(2,4)=80.0_dp/3.0_dp;
!       C(2,5)=80.0_dp/3.0_dp;
!       C(2,6)=-64.0_dp/3.0_dp;
!       C(3,1)=-64.0_dp/3.0_dp;
!       C(3,2)=8.0_dp/3.0_dp;
!       C(3,3)=512.0_dp/3.0_dp;
!       C(3,4)=-64.0_dp/3.0_dp;
!       C(3,5)=80.0_dp/3.0_dp;
!       C(3,6)=8.0_dp/3.0_dp;
!       C(4,1)=8.0_dp/3.0_dp;
!       C(4,2)=80.0_dp/3.0_dp;
!       C(4,3)=-64.0_dp/3.0_dp;
!       C(4,4)=512.0_dp/3.0_dp;
!       C(4,5)=8.0_dp/3.0_dp;
!       C(4,6)=-64.0_dp/3.0_dp;
!       C(5,1)=-64.0_dp/3.0_dp;
!       C(5,2)=80.0_dp/3.0_dp;
!       C(5,3)=80.0_dp/3.0_dp;
!       C(5,4)=8.0_dp/3.0_dp;
!       C(5,5)=8.0_dp/3.0_dp;
!       C(5,6)=8.0_dp/3.0_dp;
!       C(6,1)=-64.0_dp/3.0_dp;
!       C(6,2)=-64.0_dp/3.0_dp;
!       C(6,3)=8.0_dp/3.0_dp;
!       C(6,4)=-64.0_dp/3.0_dp;
!       C(6,5)=8.0_dp/3.0_dp;
!       C(6,6)=-64.0_dp/3.0_dp;
!       endif
!       if(n.eq.7) then
!       C(1,1)=192.0_dp;
!       C(1,2)=-24.0_dp;
!       C(1,3)=0.0_dp;
!       C(1,4)=0.0_dp;
!       C(1,5)=24.0_dp;
!       C(1,6)=48.0_dp;
!       C(2,1)=-24.0_dp;
!       C(2,2)=-24.0_dp;
!       C(2,3)=0.0_dp;
!       C(2,4)=-48.0_dp;
!       C(2,5)=-24.0_dp;
!       C(2,6)=0.0_dp;
!       C(3,1)=0.0_dp;
!       C(3,2)=0.0_dp;
!       C(3,3)=192.0_dp;
!       C(3,4)=-24.0_dp;
!       C(3,5)=48.0_dp;
!       C(3,6)=24.0_dp;
!       C(4,1)=0.0_dp;
!       C(4,2)=-48.0_dp;
!       C(4,3)=-24.0_dp;
!       C(4,4)=-24.0_dp;
!       C(4,5)=0.0_dp;
!       C(4,6)=-24.0_dp;
!       C(5,1)=24.0_dp;
!       C(5,2)=-24.0_dp;
!       C(5,3)=48.0_dp;
!       C(5,4)=0.0_dp;
!       C(5,5)=192.0_dp;
!       C(5,6)=0.0_dp;
!       C(6,1)=48.0_dp;
!       C(6,2)=0.0_dp;
!       C(6,3)=24.0_dp;
!       C(6,4)=-24.0_dp;
!       C(6,5)=0.0_dp;
!       C(6,6)=192.0_dp;
!       endif
!       if(n.eq.8) then
!       C(1,1)=192.0_dp;
!       C(1,2)=0.0_dp;
!       C(1,3)=-24.0_dp;
!       C(1,4)=24.0_dp;
!       C(1,5)=0.0_dp;
!       C(1,6)=48.0_dp;
!       C(2,1)=0.0_dp;
!       C(2,2)=192.0_dp;
!       C(2,3)=0.0_dp;
!       C(2,4)=48.0_dp;
!       C(2,5)=-24.0_dp;
!       C(2,6)=24.0_dp;
!       C(3,1)=-24.0_dp;
!       C(3,2)=0.0_dp;
!       C(3,3)=-24.0_dp;
!       C(3,4)=-24.0_dp;
!       C(3,5)=-48.0_dp;
!       C(3,6)=0.0_dp;
!       C(4,1)=24.0_dp;
!       C(4,2)=48.0_dp;
!       C(4,3)=-24.0_dp;
!       C(4,4)=192.0_dp;
!       C(4,5)=0.0_dp;
!       C(4,6)=0.0_dp;
!       C(5,1)=0.0_dp;
!       C(5,2)=-24.0_dp;
!       C(5,3)=-48.0_dp;
!       C(5,4)=0.0_dp;
!       C(5,5)=-24.0_dp;
!       C(5,6)=-24.0_dp;
!       C(6,1)=48.0_dp;
!       C(6,2)=24.0_dp;
!       C(6,3)=0.0_dp;
!       C(6,4)=0.0_dp;
!       C(6,5)=-24.0_dp;
!       C(6,6)=192.0_dp;
!       endif
!       if(n.eq.9) then
!       C(1,1)=-8.0_dp/27.0_dp;
!       C(1,2)=-80.0_dp/27.0_dp;
!       C(1,3)=-80.0_dp/27.0_dp;
!       C(1,4)=496.0_dp/27.0_dp;
!       C(1,5)=496.0_dp/27.0_dp;
!       C(1,6)=-224.0_dp/27.0_dp;
!       C(2,1)=-80.0_dp/27.0_dp;
!       C(2,2)=-8.0_dp/27.0_dp;
!       C(2,3)=496.0_dp/27.0_dp;
!       C(2,4)=-224.0_dp/27.0_dp;
!       C(2,5)=-80.0_dp/27.0_dp;
!       C(2,6)=496.0_dp/27.0_dp;
!       C(3,1)=-80.0_dp/27.0_dp;
!       C(3,2)=496.0_dp/27.0_dp;
!       C(3,3)=-8.0_dp/27.0_dp;
!       C(3,4)=-80.0_dp/27.0_dp;
!       C(3,5)=-224.0_dp/27.0_dp;
!       C(3,6)=496.0_dp/27.0_dp;
!       C(4,1)=496.0_dp/27.0_dp;
!       C(4,2)=-224.0_dp/27.0_dp;
!       C(4,3)=-80.0_dp/27.0_dp;
!       C(4,4)=-8.0_dp/27.0_dp;
!       C(4,5)=496.0_dp/27.0_dp;
!       C(4,6)=-80.0_dp/27.0_dp;
!       C(5,1)=496.0_dp/27.0_dp;
!       C(5,2)=-80.0_dp/27.0_dp;
!       C(5,3)=-224.0_dp/27.0_dp;
!       C(5,4)=496.0_dp/27.0_dp;
!       C(5,5)=-8.0_dp/27.0_dp;
!       C(5,6)=-80.0_dp/27.0_dp;
!       C(6,1)=-224.0_dp/27.0_dp;
!       C(6,2)=496.0_dp/27.0_dp;
!       C(6,3)=496.0_dp/27.0_dp;
!       C(6,4)=-80.0_dp/27.0_dp;
!       C(6,5)=-80.0_dp/27.0_dp;
!       C(6,6)=-8.0_dp/27.0_dp;
!       endif
!       if(n.eq.10) then
!       C(1,1)=8.0_dp/3.0_dp;
!       C(1,2)=80.0_dp/3.0_dp;
!       C(1,3)=8.0_dp/3.0_dp;
!       C(1,4)=8.0_dp/3.0_dp;
!       C(1,5)=-64.0_dp/3.0_dp;
!       C(1,6)=80.0_dp/3.0_dp;
!       C(2,1)=80.0_dp/3.0_dp;
!       C(2,2)=8.0_dp/3.0_dp;
!       C(2,3)=-64.0_dp/3.0_dp;
!       C(2,4)=80.0_dp/3.0_dp;
!       C(2,5)=8.0_dp/3.0_dp;
!       C(2,6)=8.0_dp/3.0_dp;
!       C(3,1)=8.0_dp/3.0_dp;
!       C(3,2)=-64.0_dp/3.0_dp;
!       C(3,3)=-64.0_dp/3.0_dp;
!       C(3,4)=-64.0_dp/3.0_dp;
!       C(3,5)=-64.0_dp/3.0_dp;
!       C(3,6)=8.0_dp/3.0_dp;
!       C(4,1)=8.0_dp/3.0_dp;
!       C(4,2)=80.0_dp/3.0_dp;
!       C(4,3)=-64.0_dp/3.0_dp;
!       C(4,4)=512.0_dp/3.0_dp;
!       C(4,5)=8.0_dp/3.0_dp;
!       C(4,6)=-64.0_dp/3.0_dp;
!       C(5,1)=-64.0_dp/3.0_dp;
!       C(5,2)=8.0_dp/3.0_dp;
!       C(5,3)=-64.0_dp/3.0_dp;
!       C(5,4)=8.0_dp/3.0_dp;
!       C(5,5)=-64.0_dp/3.0_dp;
!       C(5,6)=-64.0_dp/3.0_dp;
!       C(6,1)=80.0_dp/3.0_dp;
!       C(6,2)=8.0_dp/3.0_dp;
!       C(6,3)=8.0_dp/3.0_dp;
!       C(6,4)=-64.0_dp/3.0_dp;
!       C(6,5)=-64.0_dp/3.0_dp;
!       C(6,6)=512.0_dp/3.0_dp;
!       endif
!       if(n.eq.11) then
!       C(1,1)=-64.0_dp/3.0_dp;
!       C(1,2)=-64.0_dp/3.0_dp;
!       C(1,3)=8.0_dp/3.0_dp;
!       C(1,4)=-64.0_dp/3.0_dp;
!       C(1,5)=8.0_dp/3.0_dp;
!       C(1,6)=-64.0_dp/3.0_dp;
!       C(2,1)=-64.0_dp/3.0_dp;
!       C(2,2)=512.0_dp/3.0_dp;
!       C(2,3)=8.0_dp/3.0_dp;
!       C(2,4)=80.0_dp/3.0_dp;
!       C(2,5)=-64.0_dp/3.0_dp;
!       C(2,6)=8.0_dp/3.0_dp;
!       C(3,1)=8.0_dp/3.0_dp;
!       C(3,2)=8.0_dp/3.0_dp;
!       C(3,3)=8.0_dp/3.0_dp;
!       C(3,4)=80.0_dp/3.0_dp;
!       C(3,5)=80.0_dp/3.0_dp;
!       C(3,6)=-64.0_dp/3.0_dp;
!       C(4,1)=-64.0_dp/3.0_dp;
!       C(4,2)=80.0_dp/3.0_dp;
!       C(4,3)=80.0_dp/3.0_dp;
!       C(4,4)=8.0_dp/3.0_dp;
!       C(4,5)=8.0_dp/3.0_dp;
!       C(4,6)=8.0_dp/3.0_dp;
!       C(5,1)=8.0_dp/3.0_dp;
!       C(5,2)=-64.0_dp/3.0_dp;
!       C(5,3)=80.0_dp/3.0_dp;
!       C(5,4)=8.0_dp/3.0_dp;
!       C(5,5)=512.0_dp/3.0_dp;
!       C(5,6)=-64.0_dp/3.0_dp;
!       C(6,1)=-64.0_dp/3.0_dp;
!       C(6,2)=8.0_dp/3.0_dp;
!       C(6,3)=-64.0_dp/3.0_dp;
!       C(6,4)=8.0_dp/3.0_dp;
!       C(6,5)=-64.0_dp/3.0_dp;
!       C(6,6)=-64.0_dp/3.0_dp;
!       endif
!       if(n.eq.12) then
!       C(1,1)=512.0_dp/3.0_dp;
!       C(1,2)=-64.0_dp/3.0_dp;
!       C(1,3)=-64.0_dp/3.0_dp;
!       C(1,4)=8.0_dp/3.0_dp;
!       C(1,5)=8.0_dp/3.0_dp;
!       C(1,6)=80.0_dp/3.0_dp;
!       C(2,1)=-64.0_dp/3.0_dp;
!       C(2,2)=-64.0_dp/3.0_dp;
!       C(2,3)=8.0_dp/3.0_dp;
!       C(2,4)=-64.0_dp/3.0_dp;
!       C(2,5)=8.0_dp/3.0_dp;
!       C(2,6)=-64.0_dp/3.0_dp;
!       C(3,1)=-64.0_dp/3.0_dp;
!       C(3,2)=8.0_dp/3.0_dp;
!       C(3,3)=512.0_dp/3.0_dp;
!       C(3,4)=-64.0_dp/3.0_dp;
!       C(3,5)=80.0_dp/3.0_dp;
!       C(3,6)=8.0_dp/3.0_dp;
!       C(4,1)=8.0_dp/3.0_dp;
!       C(4,2)=-64.0_dp/3.0_dp;
!       C(4,3)=-64.0_dp/3.0_dp;
!       C(4,4)=-64.0_dp/3.0_dp;
!       C(4,5)=-64.0_dp/3.0_dp;
!       C(4,6)=8.0_dp/3.0_dp;
!       C(5,1)=8.0_dp/3.0_dp;
!       C(5,2)=8.0_dp/3.0_dp;
!       C(5,3)=80.0_dp/3.0_dp;
!       C(5,4)=-64.0_dp/3.0_dp;
!       C(5,5)=8.0_dp/3.0_dp;
!       C(5,6)=80.0_dp/3.0_dp;
!       C(6,1)=80.0_dp/3.0_dp;
!       C(6,2)=-64.0_dp/3.0_dp;
!       C(6,3)=8.0_dp/3.0_dp;
!       C(6,4)=8.0_dp/3.0_dp;
!       C(6,5)=80.0_dp/3.0_dp;
!       C(6,6)=8.0_dp/3.0_dp;
!       endif
!       if(n.eq.13) then
!       C(1,1)=-8.0_dp/27.0_dp;
!       C(1,2)=-80.0_dp/27.0_dp;
!       C(1,3)=-80.0_dp/27.0_dp;
!       C(1,4)=496.0_dp/27.0_dp;
!       C(1,5)=496.0_dp/27.0_dp;
!       C(1,6)=-224.0_dp/27.0_dp;
!       C(2,1)=-80.0_dp/27.0_dp;
!       C(2,2)=-8.0_dp/27.0_dp;
!       C(2,3)=496.0_dp/27.0_dp;
!       C(2,4)=-224.0_dp/27.0_dp;
!       C(2,5)=-80.0_dp/27.0_dp;
!       C(2,6)=496.0_dp/27.0_dp;
!       C(3,1)=-80.0_dp/27.0_dp;
!       C(3,2)=496.0_dp/27.0_dp;
!       C(3,3)=-8.0_dp/27.0_dp;
!       C(3,4)=-80.0_dp/27.0_dp;
!       C(3,5)=-224.0_dp/27.0_dp;
!       C(3,6)=496.0_dp/27.0_dp;
!       C(4,1)=496.0_dp/27.0_dp;
!       C(4,2)=-224.0_dp/27.0_dp;
!       C(4,3)=-80.0_dp/27.0_dp;
!       C(4,4)=-8.0_dp/27.0_dp;
!       C(4,5)=496.0_dp/27.0_dp;
!       C(4,6)=-80.0_dp/27.0_dp;
!       C(5,1)=496.0_dp/27.0_dp;
!       C(5,2)=-80.0_dp/27.0_dp;
!       C(5,3)=-224.0_dp/27.0_dp;
!       C(5,4)=496.0_dp/27.0_dp;
!       C(5,5)=-8.0_dp/27.0_dp;
!       C(5,6)=-80.0_dp/27.0_dp;
!       C(6,1)=-224.0_dp/27.0_dp;
!       C(6,2)=496.0_dp/27.0_dp;
!       C(6,3)=496.0_dp/27.0_dp;
!       C(6,4)=-80.0_dp/27.0_dp;
!       C(6,5)=-80.0_dp/27.0_dp;
!       C(6,6)=-8.0_dp/27.0_dp;
!       endif
!       if(n.eq.14) then
!       C(1,1)=512.0_dp/3.0_dp;
!       C(1,2)=-64.0_dp/3.0_dp;
!       C(1,3)=-64.0_dp/3.0_dp;
!       C(1,4)=8.0_dp/3.0_dp;
!       C(1,5)=8.0_dp/3.0_dp;
!       C(1,6)=80.0_dp/3.0_dp;
!       C(2,1)=-64.0_dp/3.0_dp;
!       C(2,2)=512.0_dp/3.0_dp;
!       C(2,3)=8.0_dp/3.0_dp;
!       C(2,4)=80.0_dp/3.0_dp;
!       C(2,5)=-64.0_dp/3.0_dp;
!       C(2,6)=8.0_dp/3.0_dp;
!       C(3,1)=-64.0_dp/3.0_dp;
!       C(3,2)=8.0_dp/3.0_dp;
!       C(3,3)=-64.0_dp/3.0_dp;
!       C(3,4)=8.0_dp/3.0_dp;
!       C(3,5)=-64.0_dp/3.0_dp;
!       C(3,6)=-64.0_dp/3.0_dp;
!       C(4,1)=8.0_dp/3.0_dp;
!       C(4,2)=80.0_dp/3.0_dp;
!       C(4,3)=8.0_dp/3.0_dp;
!       C(4,4)=8.0_dp/3.0_dp;
!       C(4,5)=-64.0_dp/3.0_dp;
!       C(4,6)=80.0_dp/3.0_dp;
!       C(5,1)=8.0_dp/3.0_dp;
!       C(5,2)=-64.0_dp/3.0_dp;
!       C(5,3)=-64.0_dp/3.0_dp;
!       C(5,4)=-64.0_dp/3.0_dp;
!       C(5,5)=-64.0_dp/3.0_dp;
!       C(5,6)=8.0_dp/3.0_dp;
!       C(6,1)=80.0_dp/3.0_dp;
!       C(6,2)=8.0_dp/3.0_dp;
!       C(6,3)=-64.0_dp/3.0_dp;
!       C(6,4)=80.0_dp/3.0_dp;
!       C(6,5)=8.0_dp/3.0_dp;
!       C(6,6)=8.0_dp/3.0_dp;
!       endif
!       if(n.eq.15) then
!       C(1,1)=-64.0_dp/3.0_dp;
!       C(1,2)=8.0_dp/3.0_dp;
!       C(1,3)=-64.0_dp/3.0_dp;
!       C(1,4)=8.0_dp/3.0_dp;
!       C(1,5)=-64.0_dp/3.0_dp;
!       C(1,6)=-64.0_dp/3.0_dp;
!       C(2,1)=8.0_dp/3.0_dp;
!       C(2,2)=8.0_dp/3.0_dp;
!       C(2,3)=8.0_dp/3.0_dp;
!       C(2,4)=80.0_dp/3.0_dp;
!       C(2,5)=80.0_dp/3.0_dp;
!       C(2,6)=-64.0_dp/3.0_dp;
!       C(3,1)=-64.0_dp/3.0_dp;
!       C(3,2)=8.0_dp/3.0_dp;
!       C(3,3)=512.0_dp/3.0_dp;
!       C(3,4)=-64.0_dp/3.0_dp;
!       C(3,5)=80.0_dp/3.0_dp;
!       C(3,6)=8.0_dp/3.0_dp;
!       C(4,1)=8.0_dp/3.0_dp;
!       C(4,2)=80.0_dp/3.0_dp;
!       C(4,3)=-64.0_dp/3.0_dp;
!       C(4,4)=512.0_dp/3.0_dp;
!       C(4,5)=8.0_dp/3.0_dp;
!       C(4,6)=-64.0_dp/3.0_dp;
!       C(5,1)=-64.0_dp/3.0_dp;
!       C(5,2)=80.0_dp/3.0_dp;
!       C(5,3)=80.0_dp/3.0_dp;
!       C(5,4)=8.0_dp/3.0_dp;
!       C(5,5)=8.0_dp/3.0_dp;
!       C(5,6)=8.0_dp/3.0_dp;
!       C(6,1)=-64.0_dp/3.0_dp;
!       C(6,2)=-64.0_dp/3.0_dp;
!       C(6,3)=8.0_dp/3.0_dp;
!       C(6,4)=-64.0_dp/3.0_dp;
!       C(6,5)=8.0_dp/3.0_dp;
!       C(6,6)=-64.0_dp/3.0_dp;
!       endif
!       if(n.eq.16) then
!       C(1,1)=8.0_dp/3.0_dp;
!       C(1,2)=8.0_dp/3.0_dp;
!       C(1,3)=80.0_dp/3.0_dp;
!       C(1,4)=-64.0_dp/3.0_dp;
!       C(1,5)=8.0_dp/3.0_dp;
!       C(1,6)=80.0_dp/3.0_dp;
!       C(2,1)=8.0_dp/3.0_dp;
!       C(2,2)=-64.0_dp/3.0_dp;
!       C(2,3)=-64.0_dp/3.0_dp;
!       C(2,4)=-64.0_dp/3.0_dp;
!       C(2,5)=-64.0_dp/3.0_dp;
!       C(2,6)=8.0_dp/3.0_dp;
!       C(3,1)=80.0_dp/3.0_dp;
!       C(3,2)=-64.0_dp/3.0_dp;
!       C(3,3)=8.0_dp/3.0_dp;
!       C(3,4)=8.0_dp/3.0_dp;
!       C(3,5)=80.0_dp/3.0_dp;
!       C(3,6)=8.0_dp/3.0_dp;
!       C(4,1)=-64.0_dp/3.0_dp;
!       C(4,2)=-64.0_dp/3.0_dp;
!       C(4,3)=8.0_dp/3.0_dp;
!       C(4,4)=-64.0_dp/3.0_dp;
!       C(4,5)=8.0_dp/3.0_dp;
!       C(4,6)=-64.0_dp/3.0_dp;
!       C(5,1)=8.0_dp/3.0_dp;
!       C(5,2)=-64.0_dp/3.0_dp;
!       C(5,3)=80.0_dp/3.0_dp;
!       C(5,4)=8.0_dp/3.0_dp;
!       C(5,5)=512.0_dp/3.0_dp;
!       C(5,6)=-64.0_dp/3.0_dp;
!       C(6,1)=80.0_dp/3.0_dp;
!       C(6,2)=8.0_dp/3.0_dp;
!       C(6,3)=8.0_dp/3.0_dp;
!       C(6,4)=-64.0_dp/3.0_dp;
!       C(6,5)=-64.0_dp/3.0_dp;
!       C(6,6)=512.0_dp/3.0_dp;
!       endif
!       if(n.eq.17) then
!       C(1,1)=512.0_dp/3.0_dp;
!       C(1,2)=-64.0_dp/3.0_dp;
!       C(1,3)=-64.0_dp/3.0_dp;
!       C(1,4)=8.0_dp/3.0_dp;
!       C(1,5)=8.0_dp/3.0_dp;
!       C(1,6)=80.0_dp/3.0_dp;
!       C(2,1)=-64.0_dp/3.0_dp;
!       C(2,2)=-64.0_dp/3.0_dp;
!       C(2,3)=8.0_dp/3.0_dp;
!       C(2,4)=-64.0_dp/3.0_dp;
!       C(2,5)=8.0_dp/3.0_dp;
!       C(2,6)=-64.0_dp/3.0_dp;
!       C(3,1)=-64.0_dp/3.0_dp;
!       C(3,2)=8.0_dp/3.0_dp;
!       C(3,3)=512.0_dp/3.0_dp;
!       C(3,4)=-64.0_dp/3.0_dp;
!       C(3,5)=80.0_dp/3.0_dp;
!       C(3,6)=8.0_dp/3.0_dp;
!       C(4,1)=8.0_dp/3.0_dp;
!       C(4,2)=-64.0_dp/3.0_dp;
!       C(4,3)=-64.0_dp/3.0_dp;
!       C(4,4)=-64.0_dp/3.0_dp;
!       C(4,5)=-64.0_dp/3.0_dp;
!       C(4,6)=8.0_dp/3.0_dp;
!       C(5,1)=8.0_dp/3.0_dp;
!       C(5,2)=8.0_dp/3.0_dp;
!       C(5,3)=80.0_dp/3.0_dp;
!       C(5,4)=-64.0_dp/3.0_dp;
!       C(5,5)=8.0_dp/3.0_dp;
!       C(5,6)=80.0_dp/3.0_dp;
!       C(6,1)=80.0_dp/3.0_dp;
!       C(6,2)=-64.0_dp/3.0_dp;
!       C(6,3)=8.0_dp/3.0_dp;
!       C(6,4)=8.0_dp/3.0_dp;
!       C(6,5)=80.0_dp/3.0_dp;
!       C(6,6)=8.0_dp/3.0_dp;
!       endif
!       if(n.eq.18) then
!       C(1,1)=8.0_dp/3.0_dp;
!       C(1,2)=8.0_dp/3.0_dp;
!       C(1,3)=80.0_dp/3.0_dp;
!       C(1,4)=-64.0_dp/3.0_dp;
!       C(1,5)=8.0_dp/3.0_dp;
!       C(1,6)=80.0_dp/3.0_dp;
!       C(2,1)=8.0_dp/3.0_dp;
!       C(2,2)=-64.0_dp/3.0_dp;
!       C(2,3)=-64.0_dp/3.0_dp;
!       C(2,4)=-64.0_dp/3.0_dp;
!       C(2,5)=-64.0_dp/3.0_dp;
!       C(2,6)=8.0_dp/3.0_dp;
!       C(3,1)=80.0_dp/3.0_dp;
!       C(3,2)=-64.0_dp/3.0_dp;
!       C(3,3)=8.0_dp/3.0_dp;
!       C(3,4)=8.0_dp/3.0_dp;
!       C(3,5)=80.0_dp/3.0_dp;
!       C(3,6)=8.0_dp/3.0_dp;
!       C(4,1)=-64.0_dp/3.0_dp;
!       C(4,2)=-64.0_dp/3.0_dp;
!       C(4,3)=8.0_dp/3.0_dp;
!       C(4,4)=-64.0_dp/3.0_dp;
!       C(4,5)=8.0_dp/3.0_dp;
!       C(4,6)=-64.0_dp/3.0_dp;
!       C(5,1)=8.0_dp/3.0_dp;
!       C(5,2)=-64.0_dp/3.0_dp;
!       C(5,3)=80.0_dp/3.0_dp;
!       C(5,4)=8.0_dp/3.0_dp;
!       C(5,5)=512.0_dp/3.0_dp;
!       C(5,6)=-64.0_dp/3.0_dp;
!       C(6,1)=80.0_dp/3.0_dp;
!       C(6,2)=8.0_dp/3.0_dp;
!       C(6,3)=8.0_dp/3.0_dp;
!       C(6,4)=-64.0_dp/3.0_dp;
!       C(6,5)=-64.0_dp/3.0_dp;
!       C(6,6)=512.0_dp/3.0_dp;
!       endif
!       if(n.eq.19) then
!       C(1,1)=-24.0_dp;
!       C(1,2)=-24.0_dp;
!       C(1,3)=-24.0_dp;
!       C(1,4)=0.0_dp;
!       C(1,5)=0.0_dp;
!       C(1,6)=-48.0_dp;
!       C(2,1)=-24.0_dp;
!       C(2,2)=192.0_dp;
!       C(2,3)=24.0_dp;
!       C(2,4)=48.0_dp;
!       C(2,5)=0.0_dp;
!       C(2,6)=0.0_dp;
!       C(3,1)=-24.0_dp;
!       C(3,2)=24.0_dp;
!       C(3,3)=192.0_dp;
!       C(3,4)=0.0_dp;
!       C(3,5)=48.0_dp;
!       C(3,6)=0.0_dp;
!       C(4,1)=0.0_dp;
!       C(4,2)=48.0_dp;
!       C(4,3)=0.0_dp;
!       C(4,4)=192.0_dp;
!       C(4,5)=24.0_dp;
!       C(4,6)=-24.0_dp;
!       C(5,1)=0.0_dp;
!       C(5,2)=0.0_dp;
!       C(5,3)=48.0_dp;
!       C(5,4)=24.0_dp;
!       C(5,5)=192.0_dp;
!       C(5,6)=-24.0_dp;
!       C(6,1)=-48.0_dp;
!       C(6,2)=0.0_dp;
!       C(6,3)=0.0_dp;
!       C(6,4)=-24.0_dp;
!       C(6,5)=-24.0_dp;
!       C(6,6)=-24.0_dp;
!       endif
!       if(n.eq.20) then
!       C(1,1)=192.0_dp;
!       C(1,2)=0.0_dp;
!       C(1,3)=-24.0_dp;
!       C(1,4)=24.0_dp;
!       C(1,5)=0.0_dp;
!       C(1,6)=48.0_dp;
!       C(2,1)=0.0_dp;
!       C(2,2)=192.0_dp;
!       C(2,3)=0.0_dp;
!       C(2,4)=48.0_dp;
!       C(2,5)=-24.0_dp;
!       C(2,6)=24.0_dp;
!       C(3,1)=-24.0_dp;
!       C(3,2)=0.0_dp;
!       C(3,3)=-24.0_dp;
!       C(3,4)=-24.0_dp;
!       C(3,5)=-48.0_dp;
!       C(3,6)=0.0_dp;
!       C(4,1)=24.0_dp;
!       C(4,2)=48.0_dp;
!       C(4,3)=-24.0_dp;
!       C(4,4)=192.0_dp;
!       C(4,5)=0.0_dp;
!       C(4,6)=0.0_dp;
!       C(5,1)=0.0_dp;
!       C(5,2)=-24.0_dp;
!       C(5,3)=-48.0_dp;
!       C(5,4)=0.0_dp;
!       C(5,5)=-24.0_dp;
!       C(5,6)=-24.0_dp;
!       C(6,1)=48.0_dp;
!       C(6,2)=24.0_dp;
!       C(6,3)=0.0_dp;
!       C(6,4)=0.0_dp;
!       C(6,5)=-24.0_dp;
!       C(6,6)=192.0_dp;
!       endif
!       if(n.eq.21) then
!       C(1,1)=8.0_dp/3.0_dp;
!       C(1,2)=80.0_dp/3.0_dp;
!       C(1,3)=8.0_dp/3.0_dp;
!       C(1,4)=8.0_dp/3.0_dp;
!       C(1,5)=-64.0_dp/3.0_dp;
!       C(1,6)=80.0_dp/3.0_dp;
!       C(2,1)=80.0_dp/3.0_dp;
!       C(2,2)=8.0_dp/3.0_dp;
!       C(2,3)=-64.0_dp/3.0_dp;
!       C(2,4)=80.0_dp/3.0_dp;
!       C(2,5)=8.0_dp/3.0_dp;
!       C(2,6)=8.0_dp/3.0_dp;
!       C(3,1)=8.0_dp/3.0_dp;
!       C(3,2)=-64.0_dp/3.0_dp;
!       C(3,3)=-64.0_dp/3.0_dp;
!       C(3,4)=-64.0_dp/3.0_dp;
!       C(3,5)=-64.0_dp/3.0_dp;
!       C(3,6)=8.0_dp/3.0_dp;
!       C(4,1)=8.0_dp/3.0_dp;
!       C(4,2)=80.0_dp/3.0_dp;
!       C(4,3)=-64.0_dp/3.0_dp;
!       C(4,4)=512.0_dp/3.0_dp;
!       C(4,5)=8.0_dp/3.0_dp;
!       C(4,6)=-64.0_dp/3.0_dp;
!       C(5,1)=-64.0_dp/3.0_dp;
!       C(5,2)=8.0_dp/3.0_dp;
!       C(5,3)=-64.0_dp/3.0_dp;
!       C(5,4)=8.0_dp/3.0_dp;
!       C(5,5)=-64.0_dp/3.0_dp;
!       C(5,6)=-64.0_dp/3.0_dp;
!       C(6,1)=80.0_dp/3.0_dp;
!       C(6,2)=8.0_dp/3.0_dp;
!       C(6,3)=8.0_dp/3.0_dp;
!       C(6,4)=-64.0_dp/3.0_dp;
!       C(6,5)=-64.0_dp/3.0_dp;
!       C(6,6)=512.0_dp/3.0_dp;
!       endif
!       if(n.eq.22) then
!       C(1,1)=512.0_dp/3.0_dp;
!       C(1,2)=-64.0_dp/3.0_dp;
!       C(1,3)=-64.0_dp/3.0_dp;
!       C(1,4)=8.0_dp/3.0_dp;
!       C(1,5)=8.0_dp/3.0_dp;
!       C(1,6)=80.0_dp/3.0_dp;
!       C(2,1)=-64.0_dp/3.0_dp;
!       C(2,2)=512.0_dp/3.0_dp;
!       C(2,3)=8.0_dp/3.0_dp;
!       C(2,4)=80.0_dp/3.0_dp;
!       C(2,5)=-64.0_dp/3.0_dp;
!       C(2,6)=8.0_dp/3.0_dp;
!       C(3,1)=-64.0_dp/3.0_dp;
!       C(3,2)=8.0_dp/3.0_dp;
!       C(3,3)=-64.0_dp/3.0_dp;
!       C(3,4)=8.0_dp/3.0_dp;
!       C(3,5)=-64.0_dp/3.0_dp;
!       C(3,6)=-64.0_dp/3.0_dp;
!       C(4,1)=8.0_dp/3.0_dp;
!       C(4,2)=80.0_dp/3.0_dp;
!       C(4,3)=8.0_dp/3.0_dp;
!       C(4,4)=8.0_dp/3.0_dp;
!       C(4,5)=-64.0_dp/3.0_dp;
!       C(4,6)=80.0_dp/3.0_dp;
!       C(5,1)=8.0_dp/3.0_dp;
!       C(5,2)=-64.0_dp/3.0_dp;
!       C(5,3)=-64.0_dp/3.0_dp;
!       C(5,4)=-64.0_dp/3.0_dp;
!       C(5,5)=-64.0_dp/3.0_dp;
!       C(5,6)=8.0_dp/3.0_dp;
!       C(6,1)=80.0_dp/3.0_dp;
!       C(6,2)=8.0_dp/3.0_dp;
!       C(6,3)=-64.0_dp/3.0_dp;
!       C(6,4)=80.0_dp/3.0_dp;
!       C(6,5)=8.0_dp/3.0_dp;
!       C(6,6)=8.0_dp/3.0_dp;
!       endif
!       if(n.eq.23) then
!       C(1,1)=192.0_dp;
!       C(1,2)=-24.0_dp;
!       C(1,3)=0.0_dp;
!       C(1,4)=0.0_dp;
!       C(1,5)=24.0_dp;
!       C(1,6)=48.0_dp;
!       C(2,1)=-24.0_dp;
!       C(2,2)=-24.0_dp;
!       C(2,3)=0.0_dp;
!       C(2,4)=-48.0_dp;
!       C(2,5)=-24.0_dp;
!       C(2,6)=0.0_dp;
!       C(3,1)=0.0_dp;
!       C(3,2)=0.0_dp;
!       C(3,3)=192.0_dp;
!       C(3,4)=-24.0_dp;
!       C(3,5)=48.0_dp;
!       C(3,6)=24.0_dp;
!       C(4,1)=0.0_dp;
!       C(4,2)=-48.0_dp;
!       C(4,3)=-24.0_dp;
!       C(4,4)=-24.0_dp;
!       C(4,5)=0.0_dp;
!       C(4,6)=-24.0_dp;
!       C(5,1)=24.0_dp;
!       C(5,2)=-24.0_dp;
!       C(5,3)=48.0_dp;
!       C(5,4)=0.0_dp;
!       C(5,5)=192.0_dp;
!       C(5,6)=0.0_dp;
!       C(6,1)=48.0_dp;
!       C(6,2)=0.0_dp;
!       C(6,3)=24.0_dp;
!       C(6,4)=-24.0_dp;
!       C(6,5)=0.0_dp;
!       C(6,6)=192.0_dp;
!       endif
!       if(n.eq.24) then
!       C(1,1)=-24.0_dp;
!       C(1,2)=-24.0_dp;
!       C(1,3)=-24.0_dp;
!       C(1,4)=0.0_dp;
!       C(1,5)=0.0_dp;
!       C(1,6)=-48.0_dp;
!       C(2,1)=-24.0_dp;
!       C(2,2)=192.0_dp;
!       C(2,3)=24.0_dp;
!       C(2,4)=48.0_dp;
!       C(2,5)=0.0_dp;
!       C(2,6)=0.0_dp;
!       C(3,1)=-24.0_dp;
!       C(3,2)=24.0_dp;
!       C(3,3)=192.0_dp;
!       C(3,4)=0.0_dp;
!       C(3,5)=48.0_dp;
!       C(3,6)=0.0_dp;
!       C(4,1)=0.0_dp;
!       C(4,2)=48.0_dp;
!       C(4,3)=0.0_dp;
!       C(4,4)=192.0_dp;
!       C(4,5)=24.0_dp;
!       C(4,6)=-24.0_dp;
!       C(5,1)=0.0_dp;
!       C(5,2)=0.0_dp;
!       C(5,3)=48.0_dp;
!       C(5,4)=24.0_dp;
!       C(5,5)=192.0_dp;
!       C(5,6)=-24.0_dp;
!       C(6,1)=-48.0_dp;
!       C(6,2)=0.0_dp;
!       C(6,3)=0.0_dp;
!       C(6,4)=-24.0_dp;
!       C(6,5)=-24.0_dp;
!       C(6,6)=-24.0_dp;
!       endif
!       if(n.eq.25) then
!       C(1,1)=-64.0_dp/3.0_dp;
!       C(1,2)=-64.0_dp/3.0_dp;
!       C(1,3)=8.0_dp/3.0_dp;
!       C(1,4)=-64.0_dp/3.0_dp;
!       C(1,5)=8.0_dp/3.0_dp;
!       C(1,6)=-64.0_dp/3.0_dp;
!       C(2,1)=-64.0_dp/3.0_dp;
!       C(2,2)=512.0_dp/3.0_dp;
!       C(2,3)=8.0_dp/3.0_dp;
!       C(2,4)=80.0_dp/3.0_dp;
!       C(2,5)=-64.0_dp/3.0_dp;
!       C(2,6)=8.0_dp/3.0_dp;
!       C(3,1)=8.0_dp/3.0_dp;
!       C(3,2)=8.0_dp/3.0_dp;
!       C(3,3)=8.0_dp/3.0_dp;
!       C(3,4)=80.0_dp/3.0_dp;
!       C(3,5)=80.0_dp/3.0_dp;
!       C(3,6)=-64.0_dp/3.0_dp;
!       C(4,1)=-64.0_dp/3.0_dp;
!       C(4,2)=80.0_dp/3.0_dp;
!       C(4,3)=80.0_dp/3.0_dp;
!       C(4,4)=8.0_dp/3.0_dp;
!       C(4,5)=8.0_dp/3.0_dp;
!       C(4,6)=8.0_dp/3.0_dp;
!       C(5,1)=8.0_dp/3.0_dp;
!       C(5,2)=-64.0_dp/3.0_dp;
!       C(5,3)=80.0_dp/3.0_dp;
!       C(5,4)=8.0_dp/3.0_dp;
!       C(5,5)=512.0_dp/3.0_dp;
!       C(5,6)=-64.0_dp/3.0_dp;
!       C(6,1)=-64.0_dp/3.0_dp;
!       C(6,2)=8.0_dp/3.0_dp;
!       C(6,3)=-64.0_dp/3.0_dp;
!       C(6,4)=8.0_dp/3.0_dp;
!       C(6,5)=-64.0_dp/3.0_dp;
!       C(6,6)=-64.0_dp/3.0_dp;
!       endif
!       if(n.eq.26) then
!       C(1,1)=-64.0_dp/3.0_dp;
!       C(1,2)=8.0_dp/3.0_dp;
!       C(1,3)=-64.0_dp/3.0_dp;
!       C(1,4)=8.0_dp/3.0_dp;
!       C(1,5)=-64.0_dp/3.0_dp;
!       C(1,6)=-64.0_dp/3.0_dp;
!       C(2,1)=8.0_dp/3.0_dp;
!       C(2,2)=8.0_dp/3.0_dp;
!       C(2,3)=8.0_dp/3.0_dp;
!       C(2,4)=80.0_dp/3.0_dp;
!       C(2,5)=80.0_dp/3.0_dp;
!       C(2,6)=-64.0_dp/3.0_dp;
!       C(3,1)=-64.0_dp/3.0_dp;
!       C(3,2)=8.0_dp/3.0_dp;
!       C(3,3)=512.0_dp/3.0_dp;
!       C(3,4)=-64.0_dp/3.0_dp;
!       C(3,5)=80.0_dp/3.0_dp;
!       C(3,6)=8.0_dp/3.0_dp;
!       C(4,1)=8.0_dp/3.0_dp;
!       C(4,2)=80.0_dp/3.0_dp;
!       C(4,3)=-64.0_dp/3.0_dp;
!       C(4,4)=512.0_dp/3.0_dp;
!       C(4,5)=8.0_dp/3.0_dp;
!       C(4,6)=-64.0_dp/3.0_dp;
!       C(5,1)=-64.0_dp/3.0_dp;
!       C(5,2)=80.0_dp/3.0_dp;
!       C(5,3)=80.0_dp/3.0_dp;
!       C(5,4)=8.0_dp/3.0_dp;
!       C(5,5)=8.0_dp/3.0_dp;
!       C(5,6)=8.0_dp/3.0_dp;
!       C(6,1)=-64.0_dp/3.0_dp;
!       C(6,2)=-64.0_dp/3.0_dp;
!       C(6,3)=8.0_dp/3.0_dp;
!       C(6,4)=-64.0_dp/3.0_dp;
!       C(6,5)=8.0_dp/3.0_dp;
!       C(6,6)=-64.0_dp/3.0_dp;
!       endif
!       if(n.eq.27) then
!       C(1,1)=192.0_dp;
!       C(1,2)=-24.0_dp;
!       C(1,3)=0.0_dp;
!       C(1,4)=0.0_dp;
!       C(1,5)=24.0_dp;
!       C(1,6)=48.0_dp;
!       C(2,1)=-24.0_dp;
!       C(2,2)=-24.0_dp;
!       C(2,3)=0.0_dp;
!       C(2,4)=-48.0_dp;
!       C(2,5)=-24.0_dp;
!       C(2,6)=0.0_dp;
!       C(3,1)=0.0_dp;
!       C(3,2)=0.0_dp;
!       C(3,3)=192.0_dp;
!       C(3,4)=-24.0_dp;
!       C(3,5)=48.0_dp;
!       C(3,6)=24.0_dp;
!       C(4,1)=0.0_dp;
!       C(4,2)=-48.0_dp;
!       C(4,3)=-24.0_dp;
!       C(4,4)=-24.0_dp;
!       C(4,5)=0.0_dp;
!       C(4,6)=-24.0_dp;
!       C(5,1)=24.0_dp;
!       C(5,2)=-24.0_dp;
!       C(5,3)=48.0_dp;
!       C(5,4)=0.0_dp;
!       C(5,5)=192.0_dp;
!       C(5,6)=0.0_dp;
!       C(6,1)=48.0_dp;
!       C(6,2)=0.0_dp;
!       C(6,3)=24.0_dp;
!       C(6,4)=-24.0_dp;
!       C(6,5)=0.0_dp;
!       C(6,6)=192.0_dp;
!       endif
!       if(n.eq.28) then
!       C(1,1)=192.0_dp;
!       C(1,2)=0.0_dp;
!       C(1,3)=-24.0_dp;
!       C(1,4)=24.0_dp;
!       C(1,5)=0.0_dp;
!       C(1,6)=48.0_dp;
!       C(2,1)=0.0_dp;
!       C(2,2)=192.0_dp;
!       C(2,3)=0.0_dp;
!       C(2,4)=48.0_dp;
!       C(2,5)=-24.0_dp;
!       C(2,6)=24.0_dp;
!       C(3,1)=-24.0_dp;
!       C(3,2)=0.0_dp;
!       C(3,3)=-24.0_dp;
!       C(3,4)=-24.0_dp;
!       C(3,5)=-48.0_dp;
!       C(3,6)=0.0_dp;
!       C(4,1)=24.0_dp;
!       C(4,2)=48.0_dp;
!       C(4,3)=-24.0_dp;
!       C(4,4)=192.0_dp;
!       C(4,5)=0.0_dp;
!       C(4,6)=0.0_dp;
!       C(5,1)=0.0_dp;
!       C(5,2)=-24.0_dp;
!       C(5,3)=-48.0_dp;
!       C(5,4)=0.0_dp;
!       C(5,5)=-24.0_dp;
!       C(5,6)=-24.0_dp;
!       C(6,1)=48.0_dp;
!       C(6,2)=24.0_dp;
!       C(6,3)=0.0_dp;
!       C(6,4)=0.0_dp;
!       C(6,5)=-24.0_dp;
!       C(6,6)=192.0_dp;
!       endif
!       if(n.eq.29) then
!       C(1,1)=-8.0_dp/27.0_dp;
!       C(1,2)=-80.0_dp/27.0_dp;
!       C(1,3)=-80.0_dp/27.0_dp;
!       C(1,4)=496.0_dp/27.0_dp;
!       C(1,5)=496.0_dp/27.0_dp;
!       C(1,6)=-224.0_dp/27.0_dp;
!       C(2,1)=-80.0_dp/27.0_dp;
!       C(2,2)=-8.0_dp/27.0_dp;
!       C(2,3)=496.0_dp/27.0_dp;
!       C(2,4)=-224.0_dp/27.0_dp;
!       C(2,5)=-80.0_dp/27.0_dp;
!       C(2,6)=496.0_dp/27.0_dp;
!       C(3,1)=-80.0_dp/27.0_dp;
!       C(3,2)=496.0_dp/27.0_dp;
!       C(3,3)=-8.0_dp/27.0_dp;
!       C(3,4)=-80.0_dp/27.0_dp;
!       C(3,5)=-224.0_dp/27.0_dp;
!       C(3,6)=496.0_dp/27.0_dp;
!       C(4,1)=496.0_dp/27.0_dp;
!       C(4,2)=-224.0_dp/27.0_dp;
!       C(4,3)=-80.0_dp/27.0_dp;
!       C(4,4)=-8.0_dp/27.0_dp;
!       C(4,5)=496.0_dp/27.0_dp;
!       C(4,6)=-80.0_dp/27.0_dp;
!       C(5,1)=496.0_dp/27.0_dp;
!       C(5,2)=-80.0_dp/27.0_dp;
!       C(5,3)=-224.0_dp/27.0_dp;
!       C(5,4)=496.0_dp/27.0_dp;
!       C(5,5)=-8.0_dp/27.0_dp;
!       C(5,6)=-80.0_dp/27.0_dp;
!       C(6,1)=-224.0_dp/27.0_dp;
!       C(6,2)=496.0_dp/27.0_dp;
!       C(6,3)=496.0_dp/27.0_dp;
!       C(6,4)=-80.0_dp/27.0_dp;
!       C(6,5)=-80.0_dp/27.0_dp;
!       C(6,6)=-8.0_dp/27.0_dp;
!       endif
!       if(n.eq.30) then
!       C(1,1)=8.0_dp/3.0_dp;
!       C(1,2)=80.0_dp/3.0_dp;
!       C(1,3)=8.0_dp/3.0_dp;
!       C(1,4)=8.0_dp/3.0_dp;
!       C(1,5)=-64.0_dp/3.0_dp;
!       C(1,6)=80.0_dp/3.0_dp;
!       C(2,1)=80.0_dp/3.0_dp;
!       C(2,2)=8.0_dp/3.0_dp;
!       C(2,3)=-64.0_dp/3.0_dp;
!       C(2,4)=80.0_dp/3.0_dp;
!       C(2,5)=8.0_dp/3.0_dp;
!       C(2,6)=8.0_dp/3.0_dp;
!       C(3,1)=8.0_dp/3.0_dp;
!       C(3,2)=-64.0_dp/3.0_dp;
!       C(3,3)=-64.0_dp/3.0_dp;
!       C(3,4)=-64.0_dp/3.0_dp;
!       C(3,5)=-64.0_dp/3.0_dp;
!       C(3,6)=8.0_dp/3.0_dp;
!       C(4,1)=8.0_dp/3.0_dp;
!       C(4,2)=80.0_dp/3.0_dp;
!       C(4,3)=-64.0_dp/3.0_dp;
!       C(4,4)=512.0_dp/3.0_dp;
!       C(4,5)=8.0_dp/3.0_dp;
!       C(4,6)=-64.0_dp/3.0_dp;
!       C(5,1)=-64.0_dp/3.0_dp;
!       C(5,2)=8.0_dp/3.0_dp;
!       C(5,3)=-64.0_dp/3.0_dp;
!       C(5,4)=8.0_dp/3.0_dp;
!       C(5,5)=-64.0_dp/3.0_dp;
!       C(5,6)=-64.0_dp/3.0_dp;
!       C(6,1)=80.0_dp/3.0_dp;
!       C(6,2)=8.0_dp/3.0_dp;
!       C(6,3)=8.0_dp/3.0_dp;
!       C(6,4)=-64.0_dp/3.0_dp;
!       C(6,5)=-64.0_dp/3.0_dp;
!       C(6,6)=512.0_dp/3.0_dp;
!       endif
!       if(n.eq.31) then
!       C(1,1)=-64.0_dp/3.0_dp;
!       C(1,2)=-64.0_dp/3.0_dp;
!       C(1,3)=8.0_dp/3.0_dp;
!       C(1,4)=-64.0_dp/3.0_dp;
!       C(1,5)=8.0_dp/3.0_dp;
!       C(1,6)=-64.0_dp/3.0_dp;
!       C(2,1)=-64.0_dp/3.0_dp;
!       C(2,2)=512.0_dp/3.0_dp;
!       C(2,3)=8.0_dp/3.0_dp;
!       C(2,4)=80.0_dp/3.0_dp;
!       C(2,5)=-64.0_dp/3.0_dp;
!       C(2,6)=8.0_dp/3.0_dp;
!       C(3,1)=8.0_dp/3.0_dp;
!       C(3,2)=8.0_dp/3.0_dp;
!       C(3,3)=8.0_dp/3.0_dp;
!       C(3,4)=80.0_dp/3.0_dp;
!       C(3,5)=80.0_dp/3.0_dp;
!       C(3,6)=-64.0_dp/3.0_dp;
!       C(4,1)=-64.0_dp/3.0_dp;
!       C(4,2)=80.0_dp/3.0_dp;
!       C(4,3)=80.0_dp/3.0_dp;
!       C(4,4)=8.0_dp/3.0_dp;
!       C(4,5)=8.0_dp/3.0_dp;
!       C(4,6)=8.0_dp/3.0_dp;
!       C(5,1)=8.0_dp/3.0_dp;
!       C(5,2)=-64.0_dp/3.0_dp;
!       C(5,3)=80.0_dp/3.0_dp;
!       C(5,4)=8.0_dp/3.0_dp;
!       C(5,5)=512.0_dp/3.0_dp;
!       C(5,6)=-64.0_dp/3.0_dp;
!       C(6,1)=-64.0_dp/3.0_dp;
!       C(6,2)=8.0_dp/3.0_dp;
!       C(6,3)=-64.0_dp/3.0_dp;
!       C(6,4)=8.0_dp/3.0_dp;
!       C(6,5)=-64.0_dp/3.0_dp;
!       C(6,6)=-64.0_dp/3.0_dp;
!       endif
!       if(n.eq.32) then
!       C(1,1)=512.0_dp/3.0_dp;
!       C(1,2)=-64.0_dp/3.0_dp;
!       C(1,3)=-64.0_dp/3.0_dp;
!       C(1,4)=8.0_dp/3.0_dp;
!       C(1,5)=8.0_dp/3.0_dp;
!       C(1,6)=80.0_dp/3.0_dp;
!       C(2,1)=-64.0_dp/3.0_dp;
!       C(2,2)=-64.0_dp/3.0_dp;
!       C(2,3)=8.0_dp/3.0_dp;
!       C(2,4)=-64.0_dp/3.0_dp;
!       C(2,5)=8.0_dp/3.0_dp;
!       C(2,6)=-64.0_dp/3.0_dp;
!       C(3,1)=-64.0_dp/3.0_dp;
!       C(3,2)=8.0_dp/3.0_dp;
!       C(3,3)=512.0_dp/3.0_dp;
!       C(3,4)=-64.0_dp/3.0_dp;
!       C(3,5)=80.0_dp/3.0_dp;
!       C(3,6)=8.0_dp/3.0_dp;
!       C(4,1)=8.0_dp/3.0_dp;
!       C(4,2)=-64.0_dp/3.0_dp;
!       C(4,3)=-64.0_dp/3.0_dp;
!       C(4,4)=-64.0_dp/3.0_dp;
!       C(4,5)=-64.0_dp/3.0_dp;
!       C(4,6)=8.0_dp/3.0_dp;
!       C(5,1)=8.0_dp/3.0_dp;
!       C(5,2)=8.0_dp/3.0_dp;
!       C(5,3)=80.0_dp/3.0_dp;
!       C(5,4)=-64.0_dp/3.0_dp;
!       C(5,5)=8.0_dp/3.0_dp;
!       C(5,6)=80.0_dp/3.0_dp;
!       C(6,1)=80.0_dp/3.0_dp;
!       C(6,2)=-64.0_dp/3.0_dp;
!       C(6,3)=8.0_dp/3.0_dp;
!       C(6,4)=8.0_dp/3.0_dp;
!       C(6,5)=80.0_dp/3.0_dp;
!       C(6,6)=8.0_dp/3.0_dp;
!       endif
!       if(n.eq.33) then
!       C(1,1)=-8.0_dp/27.0_dp;
!       C(1,2)=-80.0_dp/27.0_dp;
!       C(1,3)=-80.0_dp/27.0_dp;
!       C(1,4)=496.0_dp/27.0_dp;
!       C(1,5)=496.0_dp/27.0_dp;
!       C(1,6)=-224.0_dp/27.0_dp;
!       C(2,1)=-80.0_dp/27.0_dp;
!       C(2,2)=-8.0_dp/27.0_dp;
!       C(2,3)=496.0_dp/27.0_dp;
!       C(2,4)=-224.0_dp/27.0_dp;
!       C(2,5)=-80.0_dp/27.0_dp;
!       C(2,6)=496.0_dp/27.0_dp;
!       C(3,1)=-80.0_dp/27.0_dp;
!       C(3,2)=496.0_dp/27.0_dp;
!       C(3,3)=-8.0_dp/27.0_dp;
!       C(3,4)=-80.0_dp/27.0_dp;
!       C(3,5)=-224.0_dp/27.0_dp;
!       C(3,6)=496.0_dp/27.0_dp;
!       C(4,1)=496.0_dp/27.0_dp;
!       C(4,2)=-224.0_dp/27.0_dp;
!       C(4,3)=-80.0_dp/27.0_dp;
!       C(4,4)=-8.0_dp/27.0_dp;
!       C(4,5)=496.0_dp/27.0_dp;
!       C(4,6)=-80.0_dp/27.0_dp;
!       C(5,1)=496.0_dp/27.0_dp;
!       C(5,2)=-80.0_dp/27.0_dp;
!       C(5,3)=-224.0_dp/27.0_dp;
!       C(5,4)=496.0_dp/27.0_dp;
!       C(5,5)=-8.0_dp/27.0_dp;
!       C(5,6)=-80.0_dp/27.0_dp;
!       C(6,1)=-224.0_dp/27.0_dp;
!       C(6,2)=496.0_dp/27.0_dp;
!       C(6,3)=496.0_dp/27.0_dp;
!       C(6,4)=-80.0_dp/27.0_dp;
!       C(6,5)=-80.0_dp/27.0_dp;
!       C(6,6)=-8.0_dp/27.0_dp;
!       endif
!       if(n.eq.34) then
!       C(1,1)=512.0_dp/3.0_dp;
!       C(1,2)=-64.0_dp/3.0_dp;
!       C(1,3)=-64.0_dp/3.0_dp;
!       C(1,4)=8.0_dp/3.0_dp;
!       C(1,5)=8.0_dp/3.0_dp;
!       C(1,6)=80.0_dp/3.0_dp;
!       C(2,1)=-64.0_dp/3.0_dp;
!       C(2,2)=512.0_dp/3.0_dp;
!       C(2,3)=8.0_dp/3.0_dp;
!       C(2,4)=80.0_dp/3.0_dp;
!       C(2,5)=-64.0_dp/3.0_dp;
!       C(2,6)=8.0_dp/3.0_dp;
!       C(3,1)=-64.0_dp/3.0_dp;
!       C(3,2)=8.0_dp/3.0_dp;
!       C(3,3)=-64.0_dp/3.0_dp;
!       C(3,4)=8.0_dp/3.0_dp;
!       C(3,5)=-64.0_dp/3.0_dp;
!       C(3,6)=-64.0_dp/3.0_dp;
!       C(4,1)=8.0_dp/3.0_dp;
!       C(4,2)=80.0_dp/3.0_dp;
!       C(4,3)=8.0_dp/3.0_dp;
!       C(4,4)=8.0_dp/3.0_dp;
!       C(4,5)=-64.0_dp/3.0_dp;
!       C(4,6)=80.0_dp/3.0_dp;
!       C(5,1)=8.0_dp/3.0_dp;
!       C(5,2)=-64.0_dp/3.0_dp;
!       C(5,3)=-64.0_dp/3.0_dp;
!       C(5,4)=-64.0_dp/3.0_dp;
!       C(5,5)=-64.0_dp/3.0_dp;
!       C(5,6)=8.0_dp/3.0_dp;
!       C(6,1)=80.0_dp/3.0_dp;
!       C(6,2)=8.0_dp/3.0_dp;
!       C(6,3)=-64.0_dp/3.0_dp;
!       C(6,4)=80.0_dp/3.0_dp;
!       C(6,5)=8.0_dp/3.0_dp;
!       C(6,6)=8.0_dp/3.0_dp;
!       endif
!       if(n.eq.35) then
!       C(1,1)=-64.0_dp/3.0_dp;
!       C(1,2)=8.0_dp/3.0_dp;
!       C(1,3)=-64.0_dp/3.0_dp;
!       C(1,4)=8.0_dp/3.0_dp;
!       C(1,5)=-64.0_dp/3.0_dp;
!       C(1,6)=-64.0_dp/3.0_dp;
!       C(2,1)=8.0_dp/3.0_dp;
!       C(2,2)=8.0_dp/3.0_dp;
!       C(2,3)=8.0_dp/3.0_dp;
!       C(2,4)=80.0_dp/3.0_dp;
!       C(2,5)=80.0_dp/3.0_dp;
!       C(2,6)=-64.0_dp/3.0_dp;
!       C(3,1)=-64.0_dp/3.0_dp;
!       C(3,2)=8.0_dp/3.0_dp;
!       C(3,3)=512.0_dp/3.0_dp;
!       C(3,4)=-64.0_dp/3.0_dp;
!       C(3,5)=80.0_dp/3.0_dp;
!       C(3,6)=8.0_dp/3.0_dp;
!       C(4,1)=8.0_dp/3.0_dp;
!       C(4,2)=80.0_dp/3.0_dp;
!       C(4,3)=-64.0_dp/3.0_dp;
!       C(4,4)=512.0_dp/3.0_dp;
!       C(4,5)=8.0_dp/3.0_dp;
!       C(4,6)=-64.0_dp/3.0_dp;
!       C(5,1)=-64.0_dp/3.0_dp;
!       C(5,2)=80.0_dp/3.0_dp;
!       C(5,3)=80.0_dp/3.0_dp;
!       C(5,4)=8.0_dp/3.0_dp;
!       C(5,5)=8.0_dp/3.0_dp;
!       C(5,6)=8.0_dp/3.0_dp;
!       C(6,1)=-64.0_dp/3.0_dp;
!       C(6,2)=-64.0_dp/3.0_dp;
!       C(6,3)=8.0_dp/3.0_dp;
!       C(6,4)=-64.0_dp/3.0_dp;
!       C(6,5)=8.0_dp/3.0_dp;
!       C(6,6)=-64.0_dp/3.0_dp;
!       endif
!       if(n.eq.36) then
!       C(1,1)=8.0_dp/3.0_dp;
!       C(1,2)=8.0_dp/3.0_dp;
!       C(1,3)=80.0_dp/3.0_dp;
!       C(1,4)=-64.0_dp/3.0_dp;
!       C(1,5)=8.0_dp/3.0_dp;
!       C(1,6)=80.0_dp/3.0_dp;
!       C(2,1)=8.0_dp/3.0_dp;
!       C(2,2)=-64.0_dp/3.0_dp;
!       C(2,3)=-64.0_dp/3.0_dp;
!       C(2,4)=-64.0_dp/3.0_dp;
!       C(2,5)=-64.0_dp/3.0_dp;
!       C(2,6)=8.0_dp/3.0_dp;
!       C(3,1)=80.0_dp/3.0_dp;
!       C(3,2)=-64.0_dp/3.0_dp;
!       C(3,3)=8.0_dp/3.0_dp;
!       C(3,4)=8.0_dp/3.0_dp;
!       C(3,5)=80.0_dp/3.0_dp;
!       C(3,6)=8.0_dp/3.0_dp;
!       C(4,1)=-64.0_dp/3.0_dp;
!       C(4,2)=-64.0_dp/3.0_dp;
!       C(4,3)=8.0_dp/3.0_dp;
!       C(4,4)=-64.0_dp/3.0_dp;
!       C(4,5)=8.0_dp/3.0_dp;
!       C(4,6)=-64.0_dp/3.0_dp;
!       C(5,1)=8.0_dp/3.0_dp;
!       C(5,2)=-64.0_dp/3.0_dp;
!       C(5,3)=80.0_dp/3.0_dp;
!       C(5,4)=8.0_dp/3.0_dp;
!       C(5,5)=512.0_dp/3.0_dp;
!       C(5,6)=-64.0_dp/3.0_dp;
!       C(6,1)=80.0_dp/3.0_dp;
!       C(6,2)=8.0_dp/3.0_dp;
!       C(6,3)=8.0_dp/3.0_dp;
!       C(6,4)=-64.0_dp/3.0_dp;
!       C(6,5)=-64.0_dp/3.0_dp;
!       C(6,6)=512.0_dp/3.0_dp;
!       endif
!       if(n.eq.37) then
!       C(1,1)=512.0_dp/3.0_dp;
!       C(1,2)=-64.0_dp/3.0_dp;
!       C(1,3)=-64.0_dp/3.0_dp;
!       C(1,4)=8.0_dp/3.0_dp;
!       C(1,5)=8.0_dp/3.0_dp;
!       C(1,6)=80.0_dp/3.0_dp;
!       C(2,1)=-64.0_dp/3.0_dp;
!       C(2,2)=-64.0_dp/3.0_dp;
!       C(2,3)=8.0_dp/3.0_dp;
!       C(2,4)=-64.0_dp/3.0_dp;
!       C(2,5)=8.0_dp/3.0_dp;
!       C(2,6)=-64.0_dp/3.0_dp;
!       C(3,1)=-64.0_dp/3.0_dp;
!       C(3,2)=8.0_dp/3.0_dp;
!       C(3,3)=512.0_dp/3.0_dp;
!       C(3,4)=-64.0_dp/3.0_dp;
!       C(3,5)=80.0_dp/3.0_dp;
!       C(3,6)=8.0_dp/3.0_dp;
!       C(4,1)=8.0_dp/3.0_dp;
!       C(4,2)=-64.0_dp/3.0_dp;
!       C(4,3)=-64.0_dp/3.0_dp;
!       C(4,4)=-64.0_dp/3.0_dp;
!       C(4,5)=-64.0_dp/3.0_dp;
!       C(4,6)=8.0_dp/3.0_dp;
!       C(5,1)=8.0_dp/3.0_dp;
!       C(5,2)=8.0_dp/3.0_dp;
!       C(5,3)=80.0_dp/3.0_dp;
!       C(5,4)=-64.0_dp/3.0_dp;
!       C(5,5)=8.0_dp/3.0_dp;
!       C(5,6)=80.0_dp/3.0_dp;
!       C(6,1)=80.0_dp/3.0_dp;
!       C(6,2)=-64.0_dp/3.0_dp;
!       C(6,3)=8.0_dp/3.0_dp;
!       C(6,4)=8.0_dp/3.0_dp;
!       C(6,5)=80.0_dp/3.0_dp;
!       C(6,6)=8.0_dp/3.0_dp;
!       endif
!       if(n.eq.38) then
!       C(1,1)=8.0_dp/3.0_dp;
!       C(1,2)=8.0_dp/3.0_dp;
!       C(1,3)=80.0_dp/3.0_dp;
!       C(1,4)=-64.0_dp/3.0_dp;
!       C(1,5)=8.0_dp/3.0_dp;
!       C(1,6)=80.0_dp/3.0_dp;
!       C(2,1)=8.0_dp/3.0_dp;
!       C(2,2)=-64.0_dp/3.0_dp;
!       C(2,3)=-64.0_dp/3.0_dp;
!       C(2,4)=-64.0_dp/3.0_dp;
!       C(2,5)=-64.0_dp/3.0_dp;
!       C(2,6)=8.0_dp/3.0_dp;
!       C(3,1)=80.0_dp/3.0_dp;
!       C(3,2)=-64.0_dp/3.0_dp;
!       C(3,3)=8.0_dp/3.0_dp;
!       C(3,4)=8.0_dp/3.0_dp;
!       C(3,5)=80.0_dp/3.0_dp;
!       C(3,6)=8.0_dp/3.0_dp;
!       C(4,1)=-64.0_dp/3.0_dp;
!       C(4,2)=-64.0_dp/3.0_dp;
!       C(4,3)=8.0_dp/3.0_dp;
!       C(4,4)=-64.0_dp/3.0_dp;
!       C(4,5)=8.0_dp/3.0_dp;
!       C(4,6)=-64.0_dp/3.0_dp;
!       C(5,1)=8.0_dp/3.0_dp;
!       C(5,2)=-64.0_dp/3.0_dp;
!       C(5,3)=80.0_dp/3.0_dp;
!       C(5,4)=8.0_dp/3.0_dp;
!       C(5,5)=512.0_dp/3.0_dp;
!       C(5,6)=-64.0_dp/3.0_dp;
!       C(6,1)=80.0_dp/3.0_dp;
!       C(6,2)=8.0_dp/3.0_dp;
!       C(6,3)=8.0_dp/3.0_dp;
!       C(6,4)=-64.0_dp/3.0_dp;
!       C(6,5)=-64.0_dp/3.0_dp;
!       C(6,6)=512.0_dp/3.0_dp;
!       endif
!       if(n.eq.39) then
!       C(1,1)=-24.0_dp;
!       C(1,2)=-24.0_dp;
!       C(1,3)=-24.0_dp;
!       C(1,4)=0.0_dp;
!       C(1,5)=0.0_dp;
!       C(1,6)=-48.0_dp;
!       C(2,1)=-24.0_dp;
!       C(2,2)=192.0_dp;
!       C(2,3)=24.0_dp;
!       C(2,4)=48.0_dp;
!       C(2,5)=0.0_dp;
!       C(2,6)=0.0_dp;
!       C(3,1)=-24.0_dp;
!       C(3,2)=24.0_dp;
!       C(3,3)=192.0_dp;
!       C(3,4)=0.0_dp;
!       C(3,5)=48.0_dp;
!       C(3,6)=0.0_dp;
!       C(4,1)=0.0_dp;
!       C(4,2)=48.0_dp;
!       C(4,3)=0.0_dp;
!       C(4,4)=192.0_dp;
!       C(4,5)=24.0_dp;
!       C(4,6)=-24.0_dp;
!       C(5,1)=0.0_dp;
!       C(5,2)=0.0_dp;
!       C(5,3)=48.0_dp;
!       C(5,4)=24.0_dp;
!       C(5,5)=192.0_dp;
!       C(5,6)=-24.0_dp;
!       C(6,1)=-48.0_dp;
!       C(6,2)=0.0_dp;
!       C(6,3)=0.0_dp;
!       C(6,4)=-24.0_dp;
!       C(6,5)=-24.0_dp;
!       C(6,6)=-24.0_dp;
!       endif
!       if(n.eq.40) then
!       C(1,1)=192.0_dp;
!       C(1,2)=0.0_dp;
!       C(1,3)=-24.0_dp;
!       C(1,4)=24.0_dp;
!       C(1,5)=0.0_dp;
!       C(1,6)=48.0_dp;
!       C(2,1)=0.0_dp;
!       C(2,2)=192.0_dp;
!       C(2,3)=0.0_dp;
!       C(2,4)=48.0_dp;
!       C(2,5)=-24.0_dp;
!       C(2,6)=24.0_dp;
!       C(3,1)=-24.0_dp;
!       C(3,2)=0.0_dp;
!       C(3,3)=-24.0_dp;
!       C(3,4)=-24.0_dp;
!       C(3,5)=-48.0_dp;
!       C(3,6)=0.0_dp;
!       C(4,1)=24.0_dp;
!       C(4,2)=48.0_dp;
!       C(4,3)=-24.0_dp;
!       C(4,4)=192.0_dp;
!       C(4,5)=0.0_dp;
!       C(4,6)=0.0_dp;
!       C(5,1)=0.0_dp;
!       C(5,2)=-24.0_dp;
!       C(5,3)=-48.0_dp;
!       C(5,4)=0.0_dp;
!       C(5,5)=-24.0_dp;
!       C(5,6)=-24.0_dp;
!       C(6,1)=48.0_dp;
!       C(6,2)=24.0_dp;
!       C(6,3)=0.0_dp;
!       C(6,4)=0.0_dp;
!       C(6,5)=-24.0_dp;
!       C(6,6)=192.0_dp;
!       endif
!
!       end subroutine cc_gg_ttgggg
!
!
!
!
!
!
! SUBROUTINE LinkTreeParticles(TheTreeAmp,TheParticles) ! NEW
! implicit none
! type(TreeProcess) :: TheTreeAmp
! type(Particle),target :: TheParticles(:)
! integer :: iPart,PartRef,PartType,ig,iq,NPart,counter,QuarkPos(1:4)
!
!
!             counter = 0
!             do NPart=1,TheTreeAmp%NumPart
!                   TheTreeAmp%PartType(NPart) = TheParticles( TheTreeAmp%PartRef(NPart) )%PartType
!                   if( TheTreeAmp%PartType(NPart).ne.10 ) then  ! not a gluon
! !                      TheTreeAmp%NumQua = TheTreeAmp%NumQua + 1
!                      counter = counter + 1
!                      QuarkPos(counter) = NPart
!                   endif
!             enddo
!
! !           set number of gluons between quark lines
!             if( TheTreeAmp%PartType(1).ne.10 ) then ! not a gluon
!               if( TheTreeAmp%NumQua .eq. 2 ) then
!                     TheTreeAmp%NumGlu(1) = QuarkPos(2) - QuarkPos(1) - 1
!                     TheTreeAmp%NumGlu(2) = TheTreeAmp%NumPart - QuarkPos(2)
!               endif
!               if( TheTreeAmp%NumQua .eq. 4 ) then
!                     TheTreeAmp%NumGlu(1) = QuarkPos(2) - QuarkPos(1) - 1
!                     TheTreeAmp%NumGlu(2) = QuarkPos(3) - QuarkPos(2) - 1
!                     TheTreeAmp%NumGlu(3) = QuarkPos(4) - QuarkPos(3) - 1
!                     TheTreeAmp%NumGlu(4) = TheTreeAmp%NumPart - QuarkPos(4)
!               endif
!             elseif( TheTreeAmp%PartType(1).eq.0 ) then
!               if( TheTreeAmp%NumQua .eq. 2 ) then
!                     TheTreeAmp%NumGlu(1) = QuarkPos(1) - 2
!                     TheTreeAmp%NumGlu(2) = QuarkPos(2) - QuarkPos(1) - 1
!                     TheTreeAmp%NumGlu(3) = TheTreeAmp%NumPart - QuarkPos(2)
!               endif
!               if( TheTreeAmp%NumQua .eq. 4 ) then
!                     TheTreeAmp%NumGlu(1) = QuarkPos(1) - 2
!                     TheTreeAmp%NumGlu(2) = QuarkPos(2) - QuarkPos(1) - 1
!                     TheTreeAmp%NumGlu(3) = QuarkPos(3) - QuarkPos(2) - 1
!                     TheTreeAmp%NumGlu(4) = QuarkPos(4) - QuarkPos(3) - 1
!                     TheTreeAmp%NumGlu(5) = TheTreeAmp%NumPart - QuarkPos(4)
!               endif
!             endif
!
!
!    ig=0; iq=0;
!    do iPart=1,TheTreeAmp%NumPart
!       PartRef = TheTreeAmp%PartRef(iPart)
!       PartType= TheParticles(PartRef)%PartType
!       if( PartType.eq.10 ) then  ! PartType==Gluon
!             ig=ig+1
!             TheTreeAmp%Gluons(ig)%PartType => TheParticles(PartRef)%PartType
!             TheTreeAmp%Gluons(ig)%ExtRef   => TheParticles(PartRef)%ExtRef
!             TheTreeAmp%Gluons(ig)%Mass     => TheParticles(PartRef)%Mass
!             TheTreeAmp%Gluons(ig)%Mass2    => TheParticles(PartRef)%Mass2
!             TheTreeAmp%Gluons(ig)%Helicity => TheParticles(PartRef)%Helicity
!             TheTreeAmp%Gluons(ig)%Mom      => TheParticles(PartRef)%Mom
!             TheTreeAmp%Gluons(ig)%Pol      => TheParticles(PartRef)%Pol
!             if( PartType.ne.TheParticles(PartRef)%PartType ) print *,"Error1 in LinkTreeParticles"
!       else ! PartType==Quark
!             iq=iq+1
!             TheTreeAmp%Quarks(iq)%PartType => TheParticles(PartRef)%PartType
!             TheTreeAmp%Quarks(iq)%ExtRef   => TheParticles(PartRef)%ExtRef
!             TheTreeAmp%Quarks(iq)%Mass     => TheParticles(PartRef)%Mass
!             TheTreeAmp%Quarks(iq)%Mass2    => TheParticles(PartRef)%Mass2
!             TheTreeAmp%Quarks(iq)%Helicity => TheParticles(PartRef)%Helicity
!             TheTreeAmp%Quarks(iq)%Mom      => TheParticles(PartRef)%Mom
!             TheTreeAmp%Quarks(iq)%Pol      => TheParticles(PartRef)%Pol
!             if( PartType.ne.TheParticles(PartRef)%PartType ) print *,"Error2 in LinkTreeParticles"
!       endif
!    enddo
!    if( ig.ne.TheTreeAmp%NumGlu(0) .OR. iq.ne.TheTreeAmp%NumQua ) print *,"Error3 in LinkTreeParticles"
!
! return
! END SUBROUTINE
!
!
!
!
!
! SUBROUTINE InitProcess_TbTGGG(ExtParticles)
! use ModMisc
! implicit none
! type(Particle) :: ExtParticles(:)
!
!   ExtParticles(1)%PartType = -5
!   ExtParticles(1)%ExtRef   = 1
!   ExtParticles(1)%Mass = m_Top
!   ExtParticles(1)%Mass2= ExtParticles(1)%Mass**2
!
!   ExtParticles(2)%PartType = 5
!   ExtParticles(2)%ExtRef   = 2
!   ExtParticles(2)%Mass = m_Top
!   ExtParticles(2)%Mass2= ExtParticles(2)%Mass**2
!
!   ExtParticles(3)%PartType = 10
!   ExtParticles(3)%ExtRef   = 3
!   ExtParticles(3)%Mass = 0d0
!   ExtParticles(3)%Mass2= 0d0
!
!   ExtParticles(4)%PartType = 10
!   ExtParticles(4)%ExtRef   = 4
!   ExtParticles(4)%Mass = 0d0
!   ExtParticles(4)%Mass2= 0d0
!
!   ExtParticles(5)%PartType = 10
!   ExtParticles(5)%ExtRef   = 5
!   ExtParticles(5)%Mass = 0d0
!   ExtParticles(5)%Mass2= 0d0
!
!
! RETURN
! END SUBROUTINE InitProcess_TbTGGG
!
!
!
!
! SUBROUTINE InitTrees(NumQuarks,NumGluons,NumTrees,TreeAmpsDip) ! NEW
! use ModMisc
! implicit none
! integer :: NumQuarks,NumGluons,NumTrees
! type(TreeProcess) :: TreeAmpsDip(:)
! integer :: NumParticles,iTree
!
!   NumParticles=NumQuarks+NumGluons
!   do iTree=1,NumTrees
!       TreeAmpsDip(iTree)%NumPart=NumParticles
!       TreeAmpsDip(iTree)%NumQua=NumQuarks
!       allocate( TreeAmpsDip(iTree)%NumGlu(0:NumQuarks)     )
!       allocate( TreeAmpsDip(iTree)%PartRef(1:NumParticles) )
!       allocate( TreeAmpsDip(iTree)%PartType(1:NumParticles))
!       allocate( TreeAmpsDip(iTree)%Quarks(1:NumQuarks)     )
!       allocate( TreeAmpsDip(iTree)%Gluons(1:NumGluons)     )
!       TreeAmpsDip(iTree)%NumGlu(0) = NumGluons
!   enddo
!
! RETURN
! END SUBROUTINE InitTrees
!
!
!
!
!
! SUBROUTINE GenerateEvent5(Mom,Hel,ExtParticles)
! use ModMisc
! implicit none
! type(Particle) :: ExtParticles(:)
! real(8) :: Mom(1:4,1:5)
! integer :: Hel(1:5)
!
!     ExtParticles(1)%Mom(1:4) = dcmplx(Mom(1:4,1))
!     call vSpi(ExtParticles(1)%Mom(1:4),ExtParticles(1)%Mass,Hel(1),ExtParticles(1)%Pol(1:4))
!
!     ExtParticles(2)%Mom(1:4) = dcmplx(Mom(1:4,2))
!     call ubarSpi(ExtParticles(2)%Mom(1:4),ExtParticles(2)%Mass,Hel(2),ExtParticles(2)%Pol(1:4))
!
!     ExtParticles(3)%Mom(1:4) = dcmplx(Mom(1:4,3))
!     call pol_mless(ExtParticles(3)%Mom(1:4),Hel(3),ExtParticles(3)%Pol(1:4))
!
!     ExtParticles(4)%Mom(1:4) = dcmplx(Mom(1:4,4))
!     call pol_mless(ExtParticles(4)%Mom(1:4),Hel(4),ExtParticles(4)%Pol(1:4))
!
!     ExtParticles(5)%Mom(1:4) = dcmplx(Mom(1:4,5))
!     call pol_mless(ExtParticles(5)%Mom(1:4),Hel(5),ExtParticles(5)%Pol(1:4))
!
! RETURN
! END SUBROUTINE GenerateEvent5



! SUBROUTINE jetktalg(q,NPart,pjetout,Njet,weight)
! use ModMisc
! implicit none
! real(8) :: q(1:4,1:5),pjetout(1:4,1:5),weight
! integer :: NPart,Njet
!
! !     print *, "E4", q(1,5)/q(1,3)
!     weight=1d0
!     if( dabs(q(1,5)/q(1,3)).lt.1d-5 ) weight=0d0
!     if( dabs((q(1:4,5).dot.q(1:4,1))/q(1,3)**2).lt.1d-2 ) weight=0d0
!     if( dabs((q(1:4,5).dot.q(1:4,1))/q(1,4)**2).lt.1d-2 ) weight=0d0
!
! END SUBROUTINE


!
!    subroutine jetktalg(p,Np,pjetout,Njet,weight)
!    implicit none
!    real(dp), intent(in) :: p(1:4,Np)
!    integer, intent(in) :: Np
!    real(dp), intent(out) :: pjetout(4,Np),weight
!    integer, intent(out) :: Njet
!    real(dp) :: kt5, kt6, kt56
!    real(dp) :: p56(4), d55, d66, d56
!    real(dp) :: eta5,eta6,phi5,phi6, R56,ptjet,deltar,deltaphi
!
!
!     ptjet=20d0/100d0
!     deltar=0.4d0
!
!    if (Np.eq.5) then
!       kt5=sqrt(p(2,5)**2+p(3,5)**2)
!       if (kt5.lt.ptjet) then
!         weight = zero
!         Njet = 0
!         pjetout = zero
!       else
!         weight = one
!         Njet = 5   !total # of reconstructed particles
!         pjetout = p
!       endif
!     endif
!
!    if (Np.eq.6) then
!
! !  kinematics
!
!       kt5 = sqrt(p(2,5)**2+p(3,5)**2)
!       kt6 = sqrt(p(2,6)**2+p(3,6)**2)
!
!
!       eta5 = half*log((p(1,5)+p(4,5))/(p(1,5)-p(4,5)))
!       eta6 = half*log((p(1,6)+p(4,6))/(p(1,6)-p(4,6)))
!
!       phi5 = acos(p(2,5)/kt5)
!       phi6 = acos(p(2,6)/kt6)
!
!       deltaphi = dacos( (p(2,5)*p(2,6)+p(3,5)*p(3,6))/dsqrt((p(2,5)**2+p(3,5)**2)*(p(2,6)**2+p(3,6)**2))  )
! !       R56 = (eta5-eta6)**2 + (phi5-phi6)**2
!       R56 = (eta5-eta6)**2 + deltaphi**2
!
!        d55 = kt5**2
!        d66 = kt6**2
!        d56 = min(kt5,kt6)**2*R56/deltaR**2
!
!
!
!        if (d56.lt.min(d55,d66)) then
!             p56 = p(:,5) + p(:,6)                 ! combine
!             kt5=sqrt(p56(2)**2+p56(3)**2)
!             if (kt5.lt.ptjet) then
!             weight = zero
!             Njet = 0
!             pjetout = zero
!             else
!             weight = one
!             Njet = 5   !total # of reconstructed particles
!             pjetout(:,1) = p(:,1)
!             pjetout(:,2) = p(:,2)
!             pjetout(:,3) = p(:,3)
!             pjetout(:,4) = p(:,4)
!             pjetout(:,5) = p56(:)
!             endif
!         endif
!
!        if (d56.ge.min(d55,d66)) then
!
!           if (kt5.gt.ptjet.and.kt6.gt.ptjet) then
!              weight = one
!              Njet = 6
!              pjetout = p
!           endif
!
!           if (kt5.gt.ptjet.and.kt6.lt.ptjet) then
!              weight = one
!              Njet = 5
!             pjetout(:,1) = p(:,1)
!             pjetout(:,2) = p(:,2)
!             pjetout(:,3) = p(:,3)
!             pjetout(:,4) = p(:,4)
!             pjetout(:,5) = p(:,5)
!           endif
!
!           if (kt5.lt.ptjet.and.kt6.gt.ptjet) then
!              weight = one
!              Njet = 5
!             pjetout(:,1) = p(:,1)
!             pjetout(:,2) = p(:,2)
!             pjetout(:,3) = p(:,3)
!             pjetout(:,4) = p(:,4)
!             pjetout(:,5) = p(:,6)
!           endif
!
!           if (kt5.lt.ptjet.and.kt6.lt.ptjet) then
!              weight = zero
!              Njet = 0
!              pjetout = zero
!           endif
!
!           endif
!      endif
!
!
!    end subroutine jetktalg



!
!        double precision function vl(x1,x2,x3)
!        real(dp), intent(in) ::  x1,x2,x3
!        vl = x1**2+x2**2+x3**2-2.0_dp*x1*x2-2.0_dp*x1*x3-2.0_dp*x2*x3
!        end function vl
!
!
!        double precision function vel(p1,p2)
!        real(dp), intent(in) ::  p1(4),p2(4)
!        real(dp) :: p12(4),x,y,z
!        p12 = p1+p2
!        x = scr(p12,p12)
!        y = scr(p1,p1)
!        z = scr(p2,p2)
!        vel = sqrt(vl(x,y,z))/(x-y-z)
!        end function vel
!
!
!        double precision function scr(p1,p2)
!        real(dp), intent(in) :: p1(4), p2(4)
!        scr = p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)-p1(4)*p2(4)
!        end function scr
!
!
!        double complex function scc(cp1,p2)
!        complex(dp), intent(in) :: cp1(4)
!        real(dp), intent(in) :: p2(4)
!        scc = cp1(1)*p2(1)-cp1(2)*p2(2)-cp1(3)*p2(3)-cp1(4)*p2(4)
!        end function scc



       end module
