! this is the file to subtract dipoles for ttqqqq amplitudes
! for qbqb initial state
      module ModDipoles_QBQBTTBQpBQpB
      use ModAmplitudes
      use ModProcess
      use ModParameters
      use ModMisc
      use ModIntDipoles
      implicit none
      private

!----- notation for subroutines
      public :: EvalDipoles_QBQBTTBQpBQpB
      private:: jetktalg
      integer, parameter  :: dp = selected_real_kind(15)

      contains

!  so, we assume that there is a list of dipoles that we need to go through
!  I will label things as 0-> bar t(p1)+t(p2) + bar q(p3)+q (p4)
!   + bar q1(p5) + q1 (p6)
!  and here, therefore, 4 and 6 will refer to the inital state
!  and 3 and 5 are in the final state


      subroutine EvalDipoles_QBQBTTBQpBQpB(p,sum_dip)
      real(dp), intent(out) ::  sum_dip
      real(dp), intent(in) :: p(4,6)
      integer, parameter :: ndip = 2
      integer, parameter :: in1 = 4
      integer, parameter :: in2 = 6
      real(dp) :: res
      integer ::  dip(ndip,3)
      data dip(1,1)/3/ , dip(1,2)/4/, dip(1,3)/6/
      data dip(2,1)/5/ , dip(2,2)/6/, dip(2,3)/4/
      real(dp) :: mass(6)
!-----flavor list qm -- massive quark, qu -- massless quark
!----- gl -- gluon
      character :: fl(6)*2
      integer :: n,i,j,k

      fl =    'qu'

      fl(1) = 'qm'
      fl(2) = 'qm'

!     mass assignment

      mass = zero

      mass(1) = m_Top
      mass(2) = m_Top


      sum_dip = zero

            do n=1,ndip

         i=dip(n,1)  ! emitted
         j=dip(n,2)  ! emittor
         k=dip(n,3)  ! spectator

         res = zero

      if (j.ne.in1.and.j.ne.in2.and.k.ne.in1.and.k.ne.in2) &
    call dipff(n,i,j,k,mass(i),mass(j),mass(k),fl(i),fl(j),p,res)

      if ((j.ne.in1.and.j.ne.in2).and.(k.eq.in1.or.k.eq.in2)) &
    call dipfi(n,i,j,k,mass(i),mass(j),mass(k),fl(i),fl(j),p,res)

      if ( (k.ne.in1.and.k.ne.in2).and.(j.eq.in1.or.j.eq.in2)) &
     call dipif(n,i,j,k,mass(i),mass(j),mass(k),fl(i),fl(j),p,res)

      if ( (j.eq.in1.and.k.eq.in2).or.(j.eq.in2.and.k.eq.in1)) &
    call dipii(n,i,j,k,mass(i),mass(j),mass(k),fl(i),fl(j),p,res)


        sum_dip = sum_dip + res


        enddo

      end subroutine


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


!      dipole subroutines

       subroutine dipff(n,i,j,k,mi,mj,mk,fl1,fl2,p,res)
       implicit none
       integer, intent(in) :: n,i,j,k
       real(dp), intent(in) :: mi,mj,mk
       character, intent(in) :: fl1*2,fl2*2
       real(dp), intent(in) :: p(4,6)
       real(dp), intent(out) :: res
       complex(dp) :: cres
       real(dp) :: mij
       real(dp) :: C(4,4)
       real(dp) ::  pi(4), pj(4), pk(4), pkt(4), pij(4), qq(4), qqij(4)
       real(dp) :: q(4,5), paux(4), pijk(4)
       real(dp) :: xx1, xx2, zi, zj, yijk, vijk, zim, kap
       real(dp) :: pij2,pijk2,tvijk,q2,tpij2, qij2, zjm
       real(dp) ::  pijt(4)
       real(dp) :: zp, zm, xxx, viji
       integer :: i1, i2, i3, i4, i5, i6, hel(5), i8
       character :: fl(5)*2
       integer :: pos, iTree
       real(dp) :: diag, offdiag
       complex(dp) :: xp, xm
       complex(dp) :: HH(-1:1,-1:1)
       type(Particle) :: ExtParticles(1:5)
       type(TreeProcess) :: TreeAmpsDip(1:4)
       complex(8)        :: ResAmpsDip_mt(-1:1,1:4)
       integer :: c1, c2, j1
       complex(8)        :: Bm(-1:1,1:4,-1:1,-1:1,-1:1,-1:1)
       complex(8)        :: Am(-1:1,-1:1,1:4,1:4)
       complex(dp) :: POL1(-1:1,4)
       real(dp) ::  pjetout(4,5), weight
       integer :: Njet

       res = zero
       cres = (0d0,0d0)
       weight = zero    ! this weight is the 0 if counterevent fails to pass
                        ! jet cuts



      call InitTrees(4,1,4,TreeAmpsDip)
      call InitProcess_TbTQbQG(ExtParticles(1:5))
      TreeAmpsDip(1)%PartRef(1:5) = (/1,2,3,4,5/)
      TreeAmpsDip(2)%PartRef(1:5) = (/1,5,2,3,4/)
      TreeAmpsDip(3)%PartRef(1:5) = (/1,2,5,3,4/)
      TreeAmpsDip(4)%PartRef(1:5) = (/1,2,3,5,4/)

      do iTree=1,4
      call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles(1:5))
      enddo
!-----

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
        zj = 1-zi
        yijk = scr(pi,pj)/(scr(pi,pj)+scr(pi,pk)+scr(pj,pk))

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
          fl(3) = 'qu'
          q(:,4) = p(:,4)
          fl(4) = 'qu'
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
     call jetktalg(q,5,pjetout,Njet,weight)



     if(weight.eq.one) then

        call cc_qq_ttqqgg(n,C)


!--- after momentum mapping -- sum over colors and polarizations


       do i1=-1,1,2
          do i2 = -1,1,2
             do i3 = -1,1,2
                do i4 = -1,1,2
                   do i5=-1,1,2

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

      call GenerateEventttqqg(q,hel,ExtParticles(1:5))


   if (pos.eq.1) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(1)%pOL(1:4)
   endif


   if (pos.eq.2) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(2)%pOL(1:4)
   endif


   if (pos.eq.3) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(3)%pOL(1:4)
   endif

   if (pos.eq.4) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(4)%pOL(1:4)
   endif

   if (pos.eq.5) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%gLUONS(1)%pOL(1:4)
   endif


      do i6 = 1,4
      call EvalTree2(TreeAmpsDip(i6),Bm(i1,i6,i2,i3,i4,i5))
      enddo


      enddo
      enddo
      enddo
      enddo
      enddo


      do i1=-1,1,2
         do j1 = -1,1,2
            do c1 = 1,4
               do c2 = 1,4

       Am(i1,j1,c1,c2) = (0.0_dp,0.0_dp)

      do i2=-1,1,2
      do i3=-1,1,2
      do i4=-1,1,2
      do i5=-1,1,2

         Am(i1,j1,c1,c2) = Am(i1,j1,c1,c2) + &
  Bm(i1,c1,i2,i3,i4,i5)*conjg(Bm(j1,c2,i2,i3,i4,i5))

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

         kap = zero

      diag= two*(one/(one-zi*(one-yijk)) &
      +one/(one-zj*(one-yijk)) &
       -(two - kap*zp*zm)/vijk  )

       offdiag = two/vijk/scr(pi,pj)


       xm = scc(POL1(-1,:),paux)
       xp = scc(POL1(1,:),paux)

      HH(-1,1)=offdiag*conjg(xm)*xp
      HH(1,-1)=offdiag*conjg(xp)*xm
      HH(-1,-1)=diag + offdiag*conjg(xm)*xm
      HH(1,1) = diag + offdiag*conjg(xp)*xp


      endif


      if ( (fl1.eq.'gl').and.(fl2.eq.'qm') ) then

       diag = two/(one-zj*(one-yijk)) &
     -tvijk/vijk*(one + zj+mj**2/scr(pi,pj))


       HH(-1,1) = cmplx(0.0_dp,0.0_dp)
       HH(1,-1) = cmplx(0.0_dp,0.0_dp)
       HH(1,1) =  cmplx(diag,0.0_dp)
       HH(-1,-1) = cmplx(diag,0.0_dp)

      endif


        do i3=1,4
           do i4=1,4
      do i1 = -1,1,2
         do i2 = -1,1,2

      cres = cres + C(i3,i4)*Am(i1,i2,i3,i4)*HH(i1,i2)

            enddo
          enddo

           enddo
         enddo

         res = one/(pij2 -mij**2)*real(cres,dp)


         endif   ! endif for weight = 1 condition

       end subroutine dipff




       subroutine dipfi(n,i,j,a,mi,mj,ma,fl1,fl2,p,res)
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
       real(dp) :: xx1, xx2, zi, zj, yijk, vijk, zim, kap
       real(dp) :: pij2,pijk2,tvijk,q2,tpij2, qij2, zjm
       real(dp) ::  pijt(4), pija(4), pija2, xija
       real(dp) :: zp, zm, xxx, viji
       integer :: i1, i2, i3, i4, i5, i6, hel(5)
       character :: fl(5)*2
       integer :: pos, iTree, c1, c2, j1
       real(dp) :: diag, offdiag
       complex(dp) :: xp, xm
       complex(dp) :: HH(-1:1,-1:1)
       type(Particle) :: ExtParticles(1:5)
       type(TreeProcess) :: TreeAmpsDip(1:6)
       complex(8)        :: ResAmpsDip_mt(-1:1,1:6)
       complex(8)        :: Bm(-1:1,1:6,-1:1,-1:1,-1:1,-1:1)
       complex(8)        :: Am(-1:1,-1:1,1:6,1:6)
       complex(dp) :: POL1(-1:1,4)
       real(dp) ::  pjetout(4,5), weight
       integer :: Njet
       real(dp), parameter :: TR = half

       res = zero
       cres = (0d0,0d0)
       weight = zero


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


          fl(1) = 'qm'
          fl(2) = 'qm'
          fl(3) = 'gl'
          fl(4) = 'gl'
          fl(5) = 'gl'

          q(:,1) = p(:,1)
          q(:,2) = p(:,2)
          q(:,3) = -pat(:)
          q(:,4) = p(:,6)
          q(:,5) = pijt(:)

           pos = 5


     call jetktalg(q,5,pjetout,Njet,weight)


     if (weight.eq.one) then

        call cc_ttggg(C)  ! this is the color correlation matrix for ttggg

!--- after momentum mapping -- sum over colors and polarizations


       do i1=-1,1,2
          do i2 = -1,1,2
             do i3 = -1,1,2
                do i4 = -1,1,2
                   do i5=-1,1,2

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
              hel(1) = i5
           endif



!-------------- generates ttggg

      call GenerateEvent5(q,hel,ExtParticles(1:5))


       if (pos.eq.1) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(1)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(1)%pOL(1:4)
      endif


       if (pos.eq.2) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(2)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(2)%pOL(1:4)
      endif


       if (pos.eq.3) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(3)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(3)%pOL(1:4)
      endif

      if (pos.eq.4) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(4)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(4)%pOL(1:4)
      endif

      if (pos.eq.5) then ! third gluon in the list of gluons is what we need
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
  Bm(i1,c1,i2,i3,i4,i5)*conjg(Bm(j1,c2,i2,i3,i4,i5))

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


      if ( (fl1.eq.'qu').and.(fl2.eq.'qu') ) then

       paux = zi*pi - zj*pj

       diag= two*TR
       offdiag = two*(-TR*two/scr(pi,pj))


       xm = scc(POL1(-1,:),paux)
       xp = scc(POL1(1,:),paux)



      HH(-1,1)=offdiag*conjg(xm)*xp
      HH(1,-1)=offdiag*conjg(xp)*xm
      HH(-1,-1)=diag + offdiag*conjg(xm)*xm
      HH(1,1) = diag + offdiag*conjg(xp)*xp


      endif



      do i1 = -1,1,2
         do i2 = -1,1,2

        do i3=1,6
           do i4=1,6

      cres = cres + C(i3,i4)*Am(i1,i2,i3,i4)*HH(i1,i2)


            enddo
          enddo
         enddo
         enddo

         res = one/xija/(pij2-mij**2)*real(cres,dp)

         endif ! for weight =1

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
       real(dp) :: C(4,4)
       real(dp) ::  pi(4), pj(4), pk(4), pkt(4), pij(4), qq(4), qqij(4)
       real(dp) :: q(4,5), paux(4), pijk(4), pa(4),  pat(4)
       real(dp) :: xx1, xx2, zi, zj, yijk, vijk, zim, kap
       real(dp) :: pij2,pijk2,tvijk,q2,tpij2, qij2, zjm
       real(dp) ::  pijt(4), pija(4), pija2, xija
       real(dp) :: zp, zm, xxx, viji, pait(4), pjt(4)
       integer :: i1, i2, i3, i4, i5, i6, hel(5)
       character :: fl(5)*2
       integer :: pos, iTree
       real(dp) :: diag, offdiag
       complex(dp) :: xp, xm
       complex(dp) :: HH(-1:1,-1:1)
       type(Particle) :: ExtParticles(1:5)
       type(TreeProcess) :: TreeAmpsDip(1:4)
       complex(8)        :: ResAmpsDip_mt(-1:1,1:4)
       integer :: c1, c2, j1
       complex(8)        :: Bm(-1:1,1:4,-1:1,-1:1,-1:1,-1:1)
       complex(8)        :: Am(-1:1,-1:1,1:4,1:4)
       complex(dp) :: POL1(-1:1,4)
       real(dp) ::  pjetout(4,5), weight
       integer :: Njet


       res = zero
       cres = (0d0,0d0)
       weight = zero

      call InitTrees(4,1,4,TreeAmpsDip)

      call InitProcess_TbTQbQG(ExtParticles(1:5))
      TreeAmpsDip(1)%PartRef(1:5) = (/1,2,3,4,5/)
      TreeAmpsDip(2)%PartRef(1:5) = (/1,5,2,3,4/)
      TreeAmpsDip(3)%PartRef(1:5) = (/1,2,5,3,4/)
      TreeAmpsDip(4)%PartRef(1:5) = (/1,2,3,5,4/)




      do iTree=1,4
      call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles(1:5))
      enddo




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

        pait = xija*pa
        pjt =  pi + pj - (one - xija)*pa


          fl(1) = 'qm'
          fl(2) = 'qm'
          fl(3) = 'qu'
          fl(4) = 'qu'
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


     call jetktalg(q,5,pjetout,Njet,weight)

     if(weight.eq.one) then


        call cc_qq_ttqqgg(n,C)


!--- after momentum mapping -- sum over colors and polarizations


       do i1=-1,1,2

          do i2 = -1,1,2
             do i3 = -1,1,2
                do i4 = -1,1,2
                   do i5=-1,1,2

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


      call GenerateEventttqqg(q,hel,ExtParticles(1:5))

       if (pos.eq.1) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(1)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(1)%pOL(1:4)
      endif


       if (pos.eq.2) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(2)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(2)%pOL(1:4)
      endif


       if (pos.eq.3) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(3)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(3)%pOL(1:4)
      endif

      if (pos.eq.4) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(4)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(4)%pOL(1:4)
      endif

      if (pos.eq.5) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%gLUONS(1)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%gLUONS(1)%pOL(1:4)
      endif

      do i6 = 1,4
      call EvalTree2(TreeAmpsDip(i6),Bm(i1,i6,i2,i3,i4,i5))
      enddo

      enddo
      enddo
      enddo
      enddo
      enddo


      do i1=-1,1,2
         do j1 = -1,1,2
            do c1 = 1,4
               do c2 = 1,4

       Am(i1,j1,c1,c2) = (0.0_dp,0.0_dp)

      do i2=-1,1,2
      do i3=-1,1,2
      do i4=-1,1,2
      do i5=-1,1,2

         Am(i1,j1,c1,c2) = Am(i1,j1,c1,c2) + &
  Bm(i1,c1,i2,i3,i4,i5)*conjg(Bm(j1,c2,i2,i3,i4,i5))

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

      HH(-1,1)=offdiag*conjg(xm)*xp
      HH(1,-1)=offdiag*conjg(xp)*xm
      HH(-1,-1)=diag + offdiag*conjg(xm)*xm
      HH(1,1) = diag + offdiag*conjg(xp)*xp


      endif



      if ( (fl1.eq.'gl').and.(fl2.eq.'qu') ) then

      diag= two/(two-xija - zj)-one - xija
      offdiag = zero

       xm = scc(POL1(-1,:),paux)
       xp = scc(POL1(1,:),paux)

      HH(-1,1)=offdiag*conjg(xm)*xp
      HH(1,-1)=offdiag*conjg(xp)*xm
      HH(-1,-1)=diag + offdiag*conjg(xm)*xm
      HH(1,1) = diag + offdiag*conjg(xp)*xp


      endif





      do i1 = -1,1,2
         do i2 = -1,1,2

        do i3=1,4
           do i4=1,4

      cres = cres + C(i3,i4)*Am(i1,i2,i3,i4)*HH(i1,i2)
!      cres = cres + C(i3,i4)*Am(i1,i3)*conjg(Am(i2,i4))*HH(i1,i2)

            enddo
          enddo
         enddo
         enddo

         res = one/xija/two/scr(pa,pi)*real(cres,dp)

         endif ! -- ! for weight = 1

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
       real(dp) :: C(4,4)
       real(dp) ::  pi(4), pj(4), pk(4), pkt(4), pij(4), qq(4), qqij(4)
       real(dp) :: q(1:4,5), paux(4), pijk(4), pa(4),  pat(4)
       real(dp) :: xx1, xx2, zi, zj, yijk, vijk, zim, kap
       real(dp) :: pij2,pijk2,tvijk,q2,tpij2, qij2, zjm
       real(dp) ::  pijt(4), pija(4), pb(4),pija2, xija
       real(dp) :: pait(4), xiab
       real(dp) ::  K(4), tK(4),K2,tK2, K2tK2, k1(4)
       real(dp) :: zp, zm, xxx, viji
       integer :: i1, i2, i3, i4, i5, i6, hel(5), j
       character :: fl(5)*2
       integer :: pos, iTree
       real(dp) :: diag, offdiag
       complex(dp) :: xp, xm
       complex(dp) :: HH(-1:1,-1:1)
       type(Particle) :: ExtParticles(1:5)
       type(TreeProcess) :: TreeAmpsDip(1:4)
       integer :: c1, c2, j1
       complex(8)        :: Bm(-1:1,1:4,-1:1,-1:1,-1:1,-1:1)
       complex(8)        :: Am(-1:1,-1:1,1:4,1:4)
       complex(8)        :: ResAmpsDip_mt(-1:1,1:4)
       complex(dp) :: POL1(-1:1,4)
       real(dp) ::  pjetout(4,5), weight
       integer :: Njet, in1,in2
       real(dp), parameter :: TR=half
       real(dp), parameter :: CF=four/three
       real(dp) :: q1(1:4,5)

       res = zero
       cres = (0d0,0d0)

      call InitTrees(4,1,4,TreeAmpsDip)

      call InitProcess_TbTQbQG(ExtParticles(1:5))
      TreeAmpsDip(1)%PartRef(1:5) = (/1,2,3,4,5/)
      TreeAmpsDip(2)%PartRef(1:5) = (/1,5,2,3,4/)
      TreeAmpsDip(3)%PartRef(1:5) = (/1,2,5,3,4/)
      TreeAmpsDip(4)%PartRef(1:5) = (/1,2,3,5,4/)


      do iTree=1,4
      call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles(1:5))
      enddo


!       momentum mapping
        pi =  p(:,i)
        pa = -p(:,a) ! the `-' sign accounts for the all outgoing convention
        pb = -p(:,b) ! the `-' sing accounts for the all outgoing convention

        xiab = one - (scr(pi,pa)+scr(pi,pb))/scr(pa,pb)
        pait = xiab*pa


          fl(1) = 'qm'
          fl(2) = 'qm'
          fl(3) = 'qu'
          fl(4) = 'qu'
          fl(5) = 'gl'

       if (i.eq.3) then

          q(:,1) =  p(:,1)
          q(:,2) =  p(:,2)
          q(:,3) =  p(:,5)
          q(:,4) =  -pb(:)
          q(:,5) = -pait(:)

          in1 = 4   !  those in1 and in2 refer to the position of initial
          in2 = 5   !  state vectors in the q-array

          pos = 5

       endif

       if (i.eq.5) then

          q(:,1) = p(:,1)
          q(:,2) = p(:,2)
          q(:,3) =  p(:,3)
          q(:,4) = -pb(:)
          q(:,5) = -pait(:)

          in1 = 4   ! those in1 and in2 refer to the position of initial
          in2 = 5   ! state vectors in the q-array

          pos = 5

        endif


!-- now Lorentz transform

        K = pa + pb - pi
        tK = pait + pb
        K2 = scr(K,K)
        tK2 = scr(tK,tK)
        K2tK2 = scr(tK,K)
        xxx = K2+two*K2tK2 + tK2


        do j=1,5
     if (j.ne.in1.and.j.ne.in2) then
     k1 = q(:,j)
     k1 = k1 - two*(scr(k1,K)+scr(k1,tK))/xxx*(K+tK)+two*scr(k1,K)/K2*tK
     q(:,j) = k1
     endif
        enddo


        if (i.eq.3) then     ! this is needed for jet staff, to put
          q1=q               ! pre-jet momenta in the right order
          q1(:,3) = q(:,4)
          q1(:,4) = q(:,5)
          q1(:,5) = q(:,3)
        endif

        if (i.eq.5) then
          q1=q               ! pre-jet momenta in the right order
          q1(:,3) = q(:,4)
          q1(:,4) = q(:,5)
          q1(:,5) = q(:,3)
        endif

!        print *, i
!        print *, q1(:,1)
!        print *, q1(:,2)
!        print *, q1(:,3)
!        print *, q1(:,4)
!        print *, q1(:,5)
!        pause


     call jetktalg(q1,5,pjetout,Njet,weight)




     if(weight.eq.one) then



!------ color matrix for ttqqg
        call cc_ttqqg(C)


!--- after momentum mapping -- sum over colors and polarizations


       do i1=-1,1,2
          do i2 = -1,1,2
             do i3 = -1,1,2
                do i4 = -1,1,2
                   do i5=-1,1,2

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

           if (pos.eq.5) then   ! this basically means that the gluon
              hel(5) = i1       ! whose polarization correlation matrix
              hel(1) = i5       ! we need is in position 5
           endif

      call GenerateEventttqqg(q,hel,ExtParticles(1:5))

       if (pos.eq.1) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(1)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(1)%pOL(1:4)
      endif


       if (pos.eq.2) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(2)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(2)%pOL(1:4)
      endif

       if (pos.eq.3) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(3)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(3)%pOL(1:4)
      endif

      if (pos.eq.4) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(4)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(4)%pOL(1:4)
      endif

      if (pos.eq.5) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%gLUONS(1)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%gLUONS(1)%pOL(1:4)
      endif


      do i6 = 1,4
      call EvalTree2(TreeAmpsDip(i6),Bm(i1,i6,i2,i3,i4,i5))
      enddo

      enddo
      enddo
      enddo
      enddo
      enddo


      do i1=-1,1,2
         do j1 = -1,1,2
            do c1 = 1,4
               do c2 = 1,4

       Am(i1,j1,c1,c2) = (0.0_dp,0.0_dp)

      do i2=-1,1,2
      do i3=-1,1,2
      do i4=-1,1,2
      do i5=-1,1,2

         Am(i1,j1,c1,c2) = Am(i1,j1,c1,c2) + &
  Bm(i1,c1,i2,i3,i4,i5)*conjg(Bm(j1,c2,i2,i3,i4,i5))


      enddo
      enddo
      enddo
      enddo


      enddo
      enddo
      enddo
      enddo

      pause

!     now build the helicity matrix

       HH = zero

      if ( (fl1.eq.'gl').and.(fl2.eq.'gl') ) then

       paux = pi - scr(pi,pa)/scr(pi,pb)*pb

      diag= two*xiab/(one-xiab)+xiab*(one-xiab)
      offdiag = two*(one-xiab)/xiab*scr(pa,pb)/scr(pi,pa)/scr(pi,pb)

       xm = scc(POL1(-1,:),paux)
       xp = scc(POL1(1,:),paux)

      HH(-1,1)=offdiag*conjg(xm)*xp
      HH(1,-1)=offdiag*conjg(xp)*xm
      HH(-1,-1)=diag + offdiag*conjg(xm)*xm
      HH(1,1) = diag + offdiag*conjg(xp)*xp

      endif


      if ( (fl1.eq.'qu').and.(fl2.eq.'gl') ) then

      diag= two*CF*(one - two*xiab*(one-xiab))   ! here we change color
                                                 ! factor because we do not
                                                 ! do color averaging
      HH(-1,1)=dcmplx(0.0_dp,0.0_dp)
      HH(1,-1)=dcmplx(0.0_dp,0.0_dp)
      HH(-1,-1)=dcmplx(diag,0.0_dp)
      HH(1,1) = dcmplx(diag,0.0_dp)

      endif


      if ( (fl1.eq.'qu').and.(fl2.eq.'qu') ) then !this is the relevant piece

       paux = pi - scr(pi,pa)/scr(pi,pb)*pb

      diag= xiab
      offdiag = two*(one-xiab)/xiab*scr(pa,pb)/scr(pi,pa)/scr(pi,pb)

       xm = scc(POL1(-1,:),paux)
       xp = scc(POL1(1,:),paux)

      HH(-1,1)=offdiag*conjg(xm)*xp
      HH(1,-1)=offdiag*conjg(xp)*xm
      HH(-1,-1)=diag + offdiag*conjg(xm)*xm
      HH(1,1) = diag + offdiag*conjg(xp)*xp

      endif



      do i1 = -1,1,2
         do i2 = -1,1,2
        do i3=1,4
            do i4=1,4

      cres = cres + C(i3,i4)*Am(i1,i2,i3,i4)*HH(i1,i2)

            enddo
          enddo
         enddo
         enddo


         res = one/xiab/two/scr(pa,pi)*real(cres,dp)

         endif !  for weight = 1

        end subroutine dipii




   subroutine jetktalg(p,Np,pjetout,Njet,weight)
   implicit none
   real(dp), intent(in) :: p(1:4,Np)
   integer, intent(in) :: Np
   real(dp), intent(out) :: pjetout(4,Np),weight
   integer, intent(out) :: Njet
   real(dp) :: kt5, kt6, kt56
   real(dp) :: p56(4), d55, d66, d56
   real(dp) :: eta5,eta6,phi5,phi6, R56,ptjet,deltar,deltaphi


    ptjet=20d0/100d0
    deltar=0.4d0

   if (Np.eq.5) then
      kt5=sqrt(p(2,5)**2+p(3,5)**2)
      if (kt5.lt.ptjet) then
        weight = zero
        Njet = 0
        pjetout = zero
      else
        weight = one
        Njet = 5   !total # of reconstructed particles
        pjetout = p
      endif
    endif

   if (Np.eq.6) then

!  kinematics

      kt5 = sqrt(p(2,5)**2+p(3,5)**2)
      kt6 = sqrt(p(2,6)**2+p(3,6)**2)


      eta5 = half*log((p(1,5)+p(4,5))/(p(1,5)-p(4,5)))
      eta6 = half*log((p(1,6)+p(4,6))/(p(1,6)-p(4,6)))

      phi5 = acos(p(2,5)/kt5)
      phi6 = acos(p(2,6)/kt6)

      deltaphi = dacos( (p(2,5)*p(2,6)+p(3,5)*p(3,6))/dsqrt((p(2,5)**2+p(3,5)**2)*(p(2,6)**2+p(3,6)**2))  )
!       R56 = (eta5-eta6)**2 + (phi5-phi6)**2
      R56 = (eta5-eta6)**2 + deltaphi**2

       d55 = kt5**2
       d66 = kt6**2
       d56 = min(kt5,kt6)**2*R56/deltaR**2



       if (d56.lt.min(d55,d66)) then
            p56 = p(:,5) + p(:,6)                 ! combine
            kt5=sqrt(p56(2)**2+p56(3)**2)
            if (kt5.lt.ptjet) then
            weight = zero
            Njet = 0
            pjetout = zero
            else
            weight = one
            Njet = 5   !total # of reconstructed particles
            pjetout(:,1) = p(:,1)
            pjetout(:,2) = p(:,2)
            pjetout(:,3) = p(:,3)
            pjetout(:,4) = p(:,4)
            pjetout(:,5) = p56(:)
            endif
        endif

       if (d56.ge.min(d55,d66)) then

          if (kt5.gt.ptjet.and.kt6.gt.ptjet) then
             weight = one
             Njet = 6
             pjetout = p
          endif

          if (kt5.gt.ptjet.and.kt6.lt.ptjet) then
             weight = one
             Njet = 5
            pjetout(:,1) = p(:,1)
            pjetout(:,2) = p(:,2)
            pjetout(:,3) = p(:,3)
            pjetout(:,4) = p(:,4)
            pjetout(:,5) = p(:,5)
          endif

          if (kt5.lt.ptjet.and.kt6.gt.ptjet) then
             weight = one
             Njet = 5
            pjetout(:,1) = p(:,1)
            pjetout(:,2) = p(:,2)
            pjetout(:,3) = p(:,3)
            pjetout(:,4) = p(:,4)
            pjetout(:,5) = p(:,6)
          endif

          if (kt5.lt.ptjet.and.kt6.lt.ptjet) then
             weight = zero
             Njet = 0
             pjetout = zero
          endif

          endif
     endif


   end subroutine jetktalg





       end module
