! this is the file to subtract dipoles for ttggg+gamma amplitudes
      module ModDipoles_GGTTBGP
      use ModAmplitudes
      use ModProcess
      use ModParameters
      use ModMisc
      use ModIntDipoles
      use ModKinematics
      implicit none
      private

      public :: EvalDipoles_GGTTBGP
      integer, parameter, private   :: dp = selected_real_kind(15)
      real(dp), private :: yRnDK(1:8), Wgt_ext

      logical, parameter, private  :: invert_alphaCut = .false.
      real(dp), parameter, private  :: MomZero(1:4)=0d0

      contains

!  we have a list of dipoles that we need to go through
!  We label things as 0-> bar t(p1)+ gamma(2) + t(3) + g(p4)+g(p5)+g(p6)
!  and we  assume that gluons in the initial state have momenta p4 and p5,


      subroutine EvalDipoles_GGTTBGP(p,yRnDk1,Wgt,sum_dip)
      real(dp), intent(out) ::  sum_dip
      real(dp), intent(in) :: p(4,6)
      real(dp), intent(in) :: yRnDK1(1:8), Wgt
      integer, parameter :: ndip = 12
      integer, parameter :: in1 = 4
      integer, parameter :: in2 = 5
      real(dp) :: res
      integer ::  dip(ndip,3)
      data dip(1,1)/6/, dip(1,2)/1/, dip(1,3)/3/
      data dip(2,1)/6/, dip(2,2)/1/, dip(2,3)/4/
      data dip(3,1)/6/, dip(3,2)/1/, dip(3,3)/5/
      data dip(4,1)/6/, dip(4,2)/3/, dip(4,3)/1/
      data dip(5,1)/6/, dip(5,2)/3/, dip(5,3)/4/
      data dip(6,1)/6/, dip(6,2)/3/, dip(6,3)/5/
      data dip(7,1)/6/, dip(7,2)/4/, dip(7,3)/1/
      data dip(8,1)/6/, dip(8,2)/4/, dip(8,3)/3/
      data dip(9,1)/6/, dip(9,2)/4/, dip(9,3)/5/
      data dip(10,1)/6/, dip(10,2)/5/, dip(10,3)/1/
      data dip(11,1)/6/, dip(11,2)/5/, dip(11,3)/3/
      data dip(12,1)/6/, dip(12,2)/5/, dip(12,3)/4/
      real(dp) :: mass(6)
!-----flavor list qm -- massive quark, qu -- massless quark
!----- gl -- gluon
      character :: fl(6)*2
      integer :: n,i,j,k

       yRnDK = yRnDK1
       Wgt_ext = Wgt


        fl =    'gl'
        fl(1) = 'qm'
        fl(2) = 'gm'  ! gm is the photon
        fl(3) = 'qm'


!     mass assignment
      mass = zero
      mass(1) = m_top
      mass(3) = m_top

      sum_dip = zero

            do n=1,ndip

!               print *, '-----------------------------------------------'
!               print *, 'begin dipole', n

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
      endif

      if ((j.ne.in1.and.j.ne.in2).and.(k.eq.in1.or.k.eq.in2)) then
        call dipfi(n,i,j,k,mass(i),mass(j),mass(k),fl(i),fl(j),p,res)
      endif


      if ( (k.ne.in1.and.k.ne.in2).and.(j.eq.in1.or.j.eq.in2)) then
        call dipif(n,i,j,k,mass(i),mass(j),mass(k),fl(i),fl(j),p,res)
      endif

      if ( (j.eq.in1.and.k.eq.in2).or.(j.eq.in2.and.k.eq.in1)) then
        call dipii(n,i,j,k,mass(i),mass(j),mass(k),fl(i),fl(j),p,res)
      endif

        sum_dip = sum_dip + res
!         print *, 'res_dipole', n,";", i,j,k, res
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
       real(dp) :: C(2,2)
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
       type(TreeProcess),save :: TreeAmpsDip(1:2)
       complex(8)        :: ResAmpsDip_mt(-1:1,1:2)
       integer :: c1, c2, j1
       complex(8)        :: Bm(-1:1,1:2,-1:1,-1:1,-1:1,-1:1)
       complex(8)        :: Am(-1:1,-1:1,1:2,1:2)
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
             call InitTrees(2,3,2,TreeAmpsDip)
             call InitProcess_TbTGGG(ExtParticles(1:5))
             TreeAmpsDip(1)%PartRef(1:5) = (/1,5,2,3,4/)
             TreeAmpsDip(2)%PartRef(1:5) = (/1,5,2,4,3/)
             do iTree=1,2
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
          fl(2) = 'gm'
          q(:,3) = p(:,3)
          fl(3) = 'qm'
          q(:,4) = p(:,4)
          fl(4) = 'gl'
          q(:,5) = p(:,5)
          fl(5) = 'gl'

       if (i.eq.6) then
          q(:,j) = pijt(:)
          q(:,k) = pkt(:)
          pos = j
       endif


!---- now we can check that the subtraction kinematics passes jet cuts
!    call jetktalg(q,5,pjetout,Njet,weight)
!       weight=one
!     if(weight.eq.one) then
!----- now generate the top decay products

       if (TopDecays.ge.1) then
         call EvalPhasespace_TopDecay(q(1:4,1),yRnDk(1:4),.false., &
  MomDK(1:4,1:3),PSWgt1)
         call EvalPhasespace_TopDecay(q(1:4,3),yRnDk(5:8),.false., &
  MomDK(1:4,4:6),PSWgt2)
       elseif(TopDecays.eq.0) then
         MomDK = 0d0
         PSWgt1 = one
         PSWgt2 = one
       else
          print *, 'topdecays are not defined'
          pause
       endif


!----------------- ini   ini     final   top   top, real


  call Kinematics_TTBARPHOTON(0, &
 (/-q(1:4,4),-q(1:4,5),q(1:4,2),q(1:4,1),q(1:4,3), Momzero, &
    MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3), &
    MomDK(1:4,4),MomDK(1:4,5),MomDK(1:4,6)/), &
    (/4,5,3,1,2,6,7,8,9,10,11,12/), &
    Not_Passed_Cuts,NBin(1:NumHistograms) &
    )



     if(Not_Passed_Cuts.eq..false.) then

        call cc_gg_ttgamg_soft(n,C)
!--- after momentum mapping -- sum over colors and polarizations

        Nmax = 1

        if (TopDecays.ge.1) then ! Top decays
               Nmax(3) = -1
               Nmax(1) = -1
       if (pos.eq.4.or.pos.eq.5) then  ! this is needed because of swappping
               Nmax(pos) = -1               ! helicity index below
               Nmax(1) = 1
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
              hel(1) = i5
           endif

  call GenerateEvent52((/q(1:4,1),q(1:4,3),q(1:4,4),q(1:4,5),q(1:4,2)/), &
  momDK(1:4,1:6),(/hel(1),hel(3),hel(4),hel(5),hel(2)/),ExtParticles(1:5))




   if (pos.eq.1) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%quarks(1)%pol(1:4)
   endif

   if (pos.eq.2) then    ! photons
      print *, 'error, pos = 2, photon'
      stop
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%gluons(1)%pol(1:4)
   endif

   if (pos.eq.3) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%quarks(2)%pol(1:4)
   endif


   if (pos.eq.4) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%gluons(2)%pol(1:4)
   endif

   if (pos.eq.5) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%gluons(3)%pol(1:4)
   endif


      do i6 = 1,2
      call EvalTree2(TreeAmpsDip(i6),Bm(i1,i6,i2,i3,i4,i5))
      enddo


      enddo
      enddo
      enddo
      enddo
      enddo


      do i1=-1,1,2
         do j1 = -1,1,2
            do c1 = 1,2
               do c2 = 1,2

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

      diag= two*(one/(one-zi*(one-yijk)) &
      +one/(one-zj*(one-yijk)) &
       -(two - kappa_ff*zp*zm)/vijk )

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


           do i3=1,2
           do i4=1,2
           do i1 = -1,1,2
           do i2 = -1,1,2

      cres = cres + C(i3,i4)*Am(i1,i2,i3,i4)*HH(i1,i2)

          enddo
          enddo
          enddo
          enddo

          res = one/(pij2 -mij**2)*real(cres,dp)

!-------- account for all the weights, change the sign
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
       real(dp) :: C(2,2)
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
       type(TreeProcess),save :: TreeAmpsDip(1:2)
       complex(8)        :: ResAmpsDip_mt(-1:1,1:2)
       complex(8)        :: Bm(-1:1,1:2,-1:1,-1:1,-1:1,-1:1)
       complex(8)        :: Am(-1:1,-1:1,1:2,1:2)
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
             call InitTrees(2,3,2,TreeAmpsDip)
             call InitProcess_TbTGGG(ExtParticles(1:5))
             TreeAmpsDip(1)%PartRef(1:5) = (/1,5,2,3,4/)
             TreeAmpsDip(2)%PartRef(1:5) = (/1,5,2,4,3/)
             do iTree=1,2
             call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles(1:5))
             enddo
             first_time=.false.
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
          fl(2) = 'gm'   ! photon
          fl(3) = 'qm'
          fl(4) = 'gl'
          fl(5) = 'gl'


          q(:,1) = p(:,1)
          q(:,2) = p(:,2)
          q(:,3) = p(:,3)
          q(:,4) = p(:,4)
          q(:,5) = p(:,5)

          q(:,j) = pijt(:)
          q(:,a) = -pat(:)
          pos = j



!    call jetktalg(q,5,pjetout,Njet,weight)
!       weight=one
!     if (weight.eq.one) then


!----- now generate the top decay products

       if (TopDecays.ge.1) then
         call EvalPhasespace_TopDecay(q(1:4,1),yRnDk(1:4),.false., &
  MomDK(1:4,1:3),PSWgt1)
         call EvalPhasespace_TopDecay(q(1:4,3),yRnDk(5:8),.false., &
  MomDK(1:4,4:6),PSWgt2)
       elseif(TopDecays.eq.0) then
         MomDK = 0d0
         PSWgt1 = one
         PSWgt2 = one
       else
          print *, 'topdecays are not defined'
          pause
       endif


!-----------------     initial   initial       final  top      top

  call Kinematics_TTBARPHOTON(0, &
 (/-q(1:4,4),-q(1:4,5),q(1:4,2),q(1:4,1),q(1:4,3), Momzero, &
    MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3), &
    MomDK(1:4,4),MomDK(1:4,5),MomDK(1:4,6)/), &
    (/4,5,3,1,2,6,7,8,9,10,11,12/), &
    Not_Passed_Cuts,NBin(1:NumHistograms) &
    )

     if(Not_Passed_Cuts.eq..false.) then

        call cc_gg_ttgamg_soft(n,C)

         Nmax = 1

        if (TopDecays.ge.1) then ! Top decays
               Nmax(3) = -1
               Nmax(1) = -1
            if (pos.eq.4.or.pos.eq.5) then  ! this is needed because of swappping
               Nmax(pos) = -1               ! helicity index below
               Nmax(1) = 1
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

  call GenerateEvent52((/q(1:4,1),q(1:4,3),q(1:4,4),q(1:4,5),q(1:4,2)/), &
  momDK(1:4,1:6),(/hel(1),hel(3),hel(4),hel(5),hel(2)/),ExtParticles(1:5))

   if (pos.eq.1) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%quarks(1)%pol(1:4)
   endif

   if (pos.eq.2) then    ! photons
      print *, 'error, pos = 2, photon'
      stop
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%gluons(1)%pol(1:4)
   endif

   if (pos.eq.3) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%quarks(2)%pol(1:4)
   endif


   if (pos.eq.4) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%gluons(2)%pol(1:4)
   endif

   if (pos.eq.5) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%gluons(3)%pol(1:4)
   endif


      do i6 = 1,2
      call EvalTree2(TreeAmpsDip(i6),Bm(i1,i6,i2,i3,i4,i5))
      enddo

      enddo
      enddo
      enddo
      enddo
      enddo

      do i1=-1,1,2
         do j1 = -1,1,2
            do c1 = 1,2
               do c2 = 1,2

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

       paux = zi*pi - zj*pj

      diag= two*(one/(one-zi + one - xija)+one/(one-zj+ one-xija) - two)
      offdiag = two/scr(pi,pj)

       xm = scc(POL1(-1,:),paux)
       xp = scc(POL1(1,:),paux)

      HH(-1,1)=offdiag*conjg(xm)*xp
      HH(1,-1)=offdiag*conjg(xp)*xm
      HH(-1,-1)=diag + offdiag*conjg(xm)*xm
      HH(1,1) = diag + offdiag*conjg(xp)*xp


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

        do i3=1,2
           do i4=1,2


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
       real(dp) :: C(2,2)
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
       type(TreeProcess), save :: TreeAmpsDip(1:2)
       complex(8)        :: ResAmpsDip_mt(-1:1,1:2)
       integer :: c1, c2, j1
       complex(8)        :: Bm(-1:1,1:2,-1:1,-1:1,-1:1,-1:1)
       complex(8)        :: Am(-1:1,-1:1,1:2,1:2)
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
             call InitTrees(2,3,2,TreeAmpsDip)
             call InitProcess_TbTGGG(ExtParticles(1:5))
             TreeAmpsDip(1)%PartRef(1:5) = (/1,5,2,3,4/)
             TreeAmpsDip(2)%PartRef(1:5) = (/1,5,2,4,3/)
             do iTree=1,2
             call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles(1:5))
             enddo
             first_time=.false.
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
          fl(2) = 'gm'
          fl(3) = 'qm'
          fl(4) = 'gl'
          fl(5) = 'gl'





          q(:,1) = p(:,1)
          q(:,2) = p(:,2)
          q(:,3) = p(:,3)
          q(:,4) = p(:,4)
          q(:,5) = p(:,5)

          q(:,a) = -pait(:)

          q(:,j) = pjt(:)
          pos = a



!    call jetktalg(q,5,pjetout,Njet,weight)
!     if(weight.eq.one) then

       if (TopDecays.ge.1) then
         call EvalPhasespace_TopDecay(q(1:4,1),yRnDk(1:4),.false., &
  MomDK(1:4,1:3),PSWgt1)
         call EvalPhasespace_TopDecay(q(1:4,3),yRnDk(5:8),.false., &
  MomDK(1:4,4:6),PSWgt2)
       elseif(TopDecays.eq.0) then
         PSWgt1 = one
         PSWgt2 = one
       else
          print *, 'topdecays are not defined'
          pause
       endif

!----------------- initial   initial        final top      top

  call Kinematics_TTBARPHOTON(0, &
 (/-q(1:4,4),-q(1:4,5),q(1:4,2),q(1:4,1),q(1:4,3), Momzero, &
    MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3), &
    MomDK(1:4,4),MomDK(1:4,5),MomDK(1:4,6)/), &
    (/4,5,3,1,2,6,7,8,9,10,11,12/), &
    Not_Passed_Cuts,NBin(1:NumHistograms) &
    )


     if(Not_Passed_Cuts.eq..false.) then

        call cc_gg_ttgamg_soft(n,C)

!--- after momentum mapping -- sum over colors and polarizations

         Nmax = 1

        if (TopDecays.ge.1) then ! Top decays
               Nmax(3) = -1
               Nmax(1) = -1
        if (pos.eq.4.or.pos.eq.5) then  ! this is needed because of swappping
               Nmax(pos) = -1               ! helicity index below
               Nmax(1) = 1
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


  call GenerateEvent52((/q(1:4,1),q(1:4,3),q(1:4,4),q(1:4,5),q(1:4,2)/), &
  momDK(1:4,1:6),(/hel(1),hel(3),hel(4),hel(5),hel(2)/),ExtParticles(1:5))


!       print *, "if dipole"
!       print *, ExtParticles(1)%Mom(1:4).dot.ExtParticles(1)%Mom(1:4)
!       print *, ExtParticles(2)%Mom(1:4).dot.ExtParticles(2)%Mom(1:4)
!       print *, ExtParticles(3)%Mom(1:4).dot.ExtParticles(3)%Mom(1:4)
!       print *, ExtParticles(4)%Mom(1:4).dot.ExtParticles(4)%Mom(1:4)
!       print *, ExtParticles(5)%Mom(1:4).dot.ExtParticles(5)%Mom(1:4)


   if (pos.eq.1) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%quarks(1)%pol(1:4)
   endif

   if (pos.eq.2) then    ! photons
      print *, 'error, pos = 2, photon'
      stop
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%gluons(1)%pol(1:4)
   endif

   if (pos.eq.3) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%quarks(2)%pol(1:4)
   endif


   if (pos.eq.4) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%gluons(2)%pol(1:4)
   endif

   if (pos.eq.5) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%gluons(3)%pol(1:4)
   endif



      do i6 = 1,2
      call EvalTree2(TreeAmpsDip(i6),Bm(i1,i6,i2,i3,i4,i5))
      enddo

      enddo
      enddo
      enddo
      enddo
      enddo



      do i1=-1,1,2
         do j1 = -1,1,2
            do c1 = 1,2
               do c2 = 1,2

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


      do i1 = -1,1,2
         do i2 = -1,1,2
           do i3=1,2
             do i4=1,2

      cres = cres + C(i3,i4)*Am(i1,i2,i3,i4)*HH(i1,i2)

            enddo
           enddo
          enddo
         enddo



         res = one/xija/two/scr(pa,pi)*real(cres,dp)


!-------- account for all the weights, change the sign
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
       real(dp) :: C(2,2)
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
       type(TreeProcess),save :: TreeAmpsDip(1:2)
       integer :: c1, c2, j1
       complex(8)        :: Bm(-1:1,1:2,-1:1,-1:1,-1:1,-1:1)
       complex(8)        :: Am(-1:1,-1:1,1:2,1:2)
       complex(8)        :: ResAmpsDip_mt(-1:1,1:2)
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
             call InitTrees(2,3,2,TreeAmpsDip)
             call InitProcess_TbTGGG(ExtParticles(1:5))
             TreeAmpsDip(1)%PartRef(1:5) = (/1,5,2,3,4/)
             TreeAmpsDip(2)%PartRef(1:5) = (/1,5,2,4,3/)
             do iTree=1,2
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
          fl(2) = 'gm'
          fl(3) = 'qm'
          fl(4) = 'gl'
          fl(5) = 'gl'


          q(:,1) = p(:,1)
          q(:,2) = p(:,2)
          q(:,3) = p(:,3)
          q(:,4) = p(:,4)
          q(:,5) = p(:,5)

          q(:,a) = -pait(:)
          q(:,b) = -pb(:)
          pos = a


!-- now Lorentz transform

        K = pa + pb - pi
        tK = pait + pb
        K2 = scr(K,K)
        tK2 = scr(tK,tK)
        K2tK2 = scr(tK,K)
        xxx = K2+two*K2tK2 + tK2


        do j=1,5
     if (j.ne.a.and.j.ne.b) then
     k1 = q(:,j)
     k1 = k1 - two*(scr(k1,K)+scr(k1,tK))/xxx*(K+tK)+two*scr(k1,K)/K2*tK
     q(:,j) = k1
     endif
        enddo

       if (TopDecays.ge.1) then
         call EvalPhasespace_TopDecay(q(1:4,1),yRnDk(1:4),.false., &
  MomDK(1:4,1:3),PSWgt1)
         call EvalPhasespace_TopDecay(q(1:4,3),yRnDk(5:8),.false., &
  MomDK(1:4,4:6),PSWgt2)
       elseif(TopDecays.eq.0) then
         PSWgt1 = one
         PSWgt2 = one
       else
          print *, 'topdecays are not defined'
          pause
       endif

!-----------------     initial   initial       final  top      top

  call Kinematics_TTBARPHOTON(0, &
 (/-q(1:4,4),-q(1:4,5),q(1:4,2),q(1:4,1),q(1:4,3), Momzero, &
    MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3), &
    MomDK(1:4,4),MomDK(1:4,5),MomDK(1:4,6)/), &
    (/4,5,3,1,2,6,7,8,9,10,11,12/), &
    Not_Passed_Cuts,NBin(1:NumHistograms) &
    )


     if(Not_Passed_Cuts.eq..false.) then

        call cc_gg_ttgamg_soft(n,C)

!--- after momentum mapping -- sum over colors and polarizations

         Nmax = 1

        if (TopDecays.ge.1) then ! Top decays
               Nmax(3) = -1
               Nmax(1) = -1
        if (pos.eq.4.or.pos.eq.5) then  ! this is needed because of swappping
               Nmax(pos) = -1           ! helicity index below
               Nmax(1) = 1
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

  call GenerateEvent52((/q(1:4,1),q(1:4,3),q(1:4,4),q(1:4,5),q(1:4,2)/), &
  momDK(1:4,1:6),(/hel(1),hel(3),hel(4),hel(5),hel(2)/),ExtParticles(1:5))


   if (pos.eq.1) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%quarks(1)%pol(1:4)
   endif

   if (pos.eq.2) then    ! photons
      print *, 'error, pos = 2, photon'
      stop
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%gluons(1)%pol(1:4)
   endif

   if (pos.eq.3) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%quarks(2)%pol(1:4)
   endif


   if (pos.eq.4) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%gluons(2)%pol(1:4)
   endif

   if (pos.eq.5) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%gluons(3)%pol(1:4)
   endif


      do i6 = 1,2
      call EvalTree2(TreeAmpsDip(i6),Bm(i1,i6,i2,i3,i4,i5))
      enddo

      enddo
      enddo
      enddo
      enddo
      enddo


      do i1=-1,1,2
         do j1 = -1,1,2
             do c1 = 1,2
               do c2 = 1,2

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

      paux = pi - scr(pi,pa)/scr(pi,pb)*pb
      diag= two*(xiab/(one-xiab)+xiab*(one-xiab))
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
           do i3=1,2
              do i4=1,2

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



      subroutine cc_gg_ttgamg_soft(n,C) ! soft cc matrix, gg-> tt+gamma+g
      integer, intent(in) :: n
      real(dp), intent(out)  :: C(2,2)

      if(n.eq.1) then
      C(1,1)=8.0_dp/9.0_dp
      C(1,2)=80.0_dp/9.0_dp
      C(2,1)=80.0_dp/9.0_dp
      C(2,2)=8.0_dp/9.0_dp
      endif
      if(n.eq.2) then
      C(1,1)=-8.0_dp
      C(1,2)=-8.0_dp
      C(2,1)=-8.0_dp
      C(2,2)=64.0_dp
      endif
      if(n.eq.3) then
      C(1,1)=64.0_dp
      C(1,2)=-8.0_dp
      C(2,1)=-8.0_dp
      C(2,2)=-8.0_dp
      endif
      if(n.eq.4) then
      C(1,1)=8.0_dp/9.0_dp
      C(1,2)=80.0_dp/9.0_dp
      C(2,1)=80.0_dp/9.0_dp
      C(2,2)=8.0_dp/9.0_dp
      endif
      if(n.eq.5) then
      C(1,1)=64.0_dp
      C(1,2)=-8.0_dp
      C(2,1)=-8.0_dp
      C(2,2)=-8.0_dp
      endif
      if(n.eq.6) then
      C(1,1)=-8.0_dp
      C(1,2)=-8.0_dp
      C(2,1)=-8.0_dp
      C(2,2)=64.0_dp
      endif
      if(n.eq.7) then
      C(1,1)=-8.0_dp
      C(1,2)=-8.0_dp
      C(2,1)=-8.0_dp
      C(2,2)=64.0_dp
      endif
      if(n.eq.8) then
      C(1,1)=64.0_dp
      C(1,2)=-8.0_dp
      C(2,1)=-8.0_dp
      C(2,2)=-8.0_dp
      endif
      if(n.eq.9) then
      C(1,1)=72.0_dp
      C(1,2)=0.0_dp
      C(2,1)=0.0_dp
      C(2,2)=72.0_dp
      endif
      if(n.eq.10) then
      C(1,1)=64.0_dp
      C(1,2)=-8.0_dp
      C(2,1)=-8.0_dp
      C(2,2)=-8.0_dp
      endif
      if(n.eq.11) then
      C(1,1)=-8.0_dp
      C(1,2)=-8.0_dp
      C(2,1)=-8.0_dp
      C(2,2)=64.0_dp
      endif
      if(n.eq.12) then
      C(1,1)=72.0_dp
      C(1,2)=0.0_dp
      C(2,1)=0.0_dp
      C(2,2)=72.0_dp
      endif

      end subroutine cc_gg_ttgamg_soft



       end module
