! this is the file to subtract dipoles for qq->tt+g+Z amplitudes
      module ModDipoles_QQBTTBGZ
      use ModAmplitudes
      use ModProcess
      use ModParameters
      use ModMisc
      use ModKinematics
      use ModIntDipoles
      use ModZDecay
      implicit none
      private

      public :: EvalDipoles_QQBTTBGZ
      integer, parameter  :: dp = selected_real_kind(15)
      real(dp), private :: yRnDK(1:10), Wgt_ext(1:2)

      real(dp), parameter, private  :: Momzero(1:4)=0d0
      logical, parameter :: invert_alphaCut = .false.

      contains

!  we have a list of dipoles that we need to go through
!  We label things as 0-> bar t(p1)+ Z(p2) + t(p3) + qb(p4)+q(p5)+g(p6)
!  and we  assume that quarks in the initial state have momenta p4 and p5,


!--- I pass TWO external weights to the subroutine, one for `up' quarks
!-------------------------------------------------, the other for 'dn' quarks

      subroutine EvalDipoles_QQBTTBGZ(p,yRnDk1,Wgt,sum_dip)
      real(dp), intent(out) ::  sum_dip(2)
      real(dp), intent(in) :: p(4,6)
      real(dp), intent(in) :: yRnDK1(1:10), Wgt(1:2)
      integer, parameter :: ndip = 12
      integer, parameter :: in1 = 4
      integer, parameter :: in2 = 5
      real(dp) :: res(2)
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
        fl(4) = 'qu'
        fl(5) = 'qu'

!     mass assignment
      mass = zero
      mass(1) = m_top
      mass(3) = m_top

      sum_dip = zero

            do n=1,ndip

         i=dip(n,1)  ! emitted
         j=dip(n,2)  ! emittor
         k=dip(n,3)  ! spectator

         res = zero

      if (j.ne.in1.and.j.ne.in2.and.k.ne.in1.and.k.ne.in2) then
        call dipff(n,i,j,k,mass(i),mass(j),mass(k),fl(i),fl(j),p,res)
! print *, "ff",n,res
      endif

      if ((j.ne.in1.and.j.ne.in2).and.(k.eq.in1.or.k.eq.in2)) then
        call dipfi(n,i,j,k,mass(i),mass(j),mass(k),fl(i),fl(j),p,res)
! print *, "fi",n,res
      endif


      if ( (k.ne.in1.and.k.ne.in2).and.(j.eq.in1.or.j.eq.in2)) then
        call dipif(n,i,j,k,mass(i),mass(j),mass(k),fl(i),fl(j),p,res)
! print *, "if",n,res
      endif

      if ( (j.eq.in1.and.k.eq.in2).or.(j.eq.in2.and.k.eq.in1)) then
        call dipii(n,i,j,k,mass(i),mass(j),mass(k),fl(i),fl(j),p,res)
! print *, "ii",n,res
      endif

      sum_dip = sum_dip + res
      enddo
! pause

      end subroutine



!      dipole subroutines
       subroutine dipff(n,i,j,k,mi,mj,mk,fl1,fl2,p,res)
       implicit none
       integer, intent(in) :: n,i,j,k
       real(dp), intent(in) :: mi,mj,mk
       character, intent(in) :: fl1*2,fl2*2
       real(dp), intent(in) :: p(4,6)
       real(dp), intent(out) :: res(2)
       complex(dp) :: cres(2)
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
       complex(8)        :: Bm(-1:1,1:2,-1:1,-1:1,-1:1,-1:1), Bm_up(-1:1,1:2,-1:1,-1:1,-1:1,-1:1), Bm_dn(-1:1,1:2,-1:1,-1:1,-1:1,-1:1)
       complex(8)        :: Am_up(-1:1,-1:1,1:2,1:2),Am_dn(-1:1,-1:1,1:2,1:2)
       complex(dp) :: POL1(-1:1,4),propZ
       real(dp) ::  pjetout(4,5), weight,couplZLL,couplGLL,couplZUU,couplGUU,couplZDD,couplGDD
       logical :: Not_Passed_Cuts
       integer :: NBin(1:NumHistograms)
       real(dp) :: MomDK(1:4,1:8),PSWgt1,PSWgt2,PSWgt3,MZ_Inv
       integer :: Njet, Nmax(5), Nhisto,N2jump
       logical, save :: first_time = .true.


       res = zero
       cres = (0d0,0d0)
       Bm =   (0d0,0d0); Bm_up =   (0d0,0d0); Bm_dn =   (0d0,0d0)
       weight = zero    ! this weight is the 0 if counterevent fails to pass
                        ! jet cuts

       if( first_time ) then
             call InitTrees(4,0,2,TreeAmpsDip,NumBoson=1)
             call InitProcess_TbTQbQZ(ExtParticles(1:5))
             TreeAmpsDip(1)%PartRef(1:5) = (/1,5,2,3,4/)
             TreeAmpsDip(2)%PartRef(1:5) = (/1,2,3,5,4/)
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
          fl(2) = 'gm'   ! z boson
          q(:,3) = p(:,3)
          fl(3) = 'qm'
          q(:,4) = p(:,4)
          fl(4) = 'qu'
          q(:,5) = p(:,5)
          fl(5) = 'qu'

       if (i.eq.6) then
          q(:,j) = pijt(:)
          q(:,k) = pkt(:)
          pos = j
       endif


       if (TopDecays.ge.1) then
         call EvalPhasespace_TopDecay(q(1:4,1),yRnDk(1:4),.false.,MomDK(1:4,1:3),PSWgt1)
         call EvalPhasespace_TopDecay(q(1:4,3),yRnDk(5:8),.false.,MomDK(1:4,4:6),PSWgt2)
       elseif(TopDecays.eq.0) then
         MomDK = 0d0
         PSWgt1 = one
         PSWgt2 = one
       else
          print *, 'topdecays are not defined'
          pause
       endif

       if( ZDecays.le.10 ) then  ! decaying on-shell Z
            MZ_Inv = m_Z
       elseif( ZDecays.gt.10 ) then  ! decaying off-shell Z
            call Error("need to implement phase space for off-shell Z's")
            ! need to think about threshold cut: EHat.le.2d0*m_Top+M_Z   when z is off-shell!
       endif
       IF( ZDECAYS.GT.0 ) THEN
          call EvalPhasespace_ZDecay(MZ_Inv,q(1:4,2),yRnDk(9:10),MomDK(1:4,7:8),PSWgt3)
          PSWgt1 = PSWgt1 * PSWgt3
       ENDIF




!----------------- ini   ini     final   top   top, real


  call Kinematics_TTBARZ(0,(/-q(1:4,4),-q(1:4,5),q(1:4,2),q(1:4,1),q(1:4,3), Momzero,MomDK(1:4,1:8)/), (/4,5,3,1,2,0,7,8,9,10,11,12,13,14/), Not_Passed_Cuts,NBin(1:NumHistograms)   )



     if(Not_Passed_Cuts.eq..false.) then

        call cc_qq_ttgamg_soft(n,C)
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
        N2Jump = 1
        if (ZDecays.ge.1 .or. ZDecays.eq.-2) then ! Z decays
            N2jump=2
        endif

       do i1=-1,Nmax(1),2
         do i2 = -1,Nmax(2),N2Jump! Z boson
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

    call SetPolarization((/q(1:4,1),q(1:4,3),q(1:4,4),q(1:4,5),q(1:4,2)/),momDK(1:4,1:8),(/hel(1),hel(3),hel(4),hel(5),hel(2)/),ExtParticles(1:5))
    if( ZDecays.gt.0 ) call ZGamLCoupl(1,hel(2),couplZLL,couplGLL)  ! charged lept


   if (pos.eq.1) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%quarks(1)%pol(1:4)
   endif

   if (pos.eq.2) then    ! z
      print *, 'error, pos = 2, z boson'
      stop
   endif

   if (pos.eq.3) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%quarks(2)%pol(1:4)
   endif


   if (pos.eq.4) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%quarks(3)%pol(1:4)
   endif

   if (pos.eq.5) then
   if (i1.eq.-1.or.i1.eq.1) POL1(i1,:)= TreeAmpsDip(1)%quarks(4)%pol(1:4)
   endif


      call ZGamQcoupl(Up_,ExtParticles(3)%Helicity,couplZUU,couplGUU)
      call ZGamQcoupl(Dn_,ExtParticles(3)%Helicity,couplZDD,couplGDD)
      couplZQQ_left_dyn=one
      couplZQQ_right_dyn=one

      if ( ZDecays .lt. 10) then
          propZ = (1d0,0d0)/dsqrt(2d0*Ga_Zexp*m_Z)  * MZ_Inv**2
      elseif (ZDecays .gt. 10) then 
          propZ=cone/(MZ_Inv**2-m_Z**2+ci*Ga_ZExp*m_Z)   * MZ_Inv**2
      endif

      do i6 = 1,2
          call EvalTree2(TreeAmpsDip(i6),Bm(i1,i6,i2,i3,i4,i5))
          if( i6.eq.1 ) then 
              Bm_up(i1,i6,i2,i3,i4,i5)=Bm(i1,i6,i2,i3,i4,i5)
              Bm_dn(i1,i6,i2,i3,i4,i5)=Bm(i1,i6,i2,i3,i4,i5)
          else
              Bm_up(i1,i6,i2,i3,i4,i5)=Bm(i1,i6,i2,i3,i4,i5)*couplZUU
              Bm_dn(i1,i6,i2,i3,i4,i5)=Bm(i1,i6,i2,i3,i4,i5)*couplZDD

              if( Zdecays.gt.0 .and. Zdecays.lt.10 ) then
                  Bm_up(i1,i6,i2,i3,i4,i5)=Bm_up(i1,i6,i2,i3,i4,i5)*propZ*couplZLL
                  Bm_dn(i1,i6,i2,i3,i4,i5)=Bm_dn(i1,i6,i2,i3,i4,i5)*propZ*couplZLL
              elseif( Zdecays.gt.10 ) then
    !             LOPartAmp(up)=BornAmps(1)%Result+BornAmps(2)%Result*( couplZUU*propZ*couplZLL + couplGUU*couplGLL )
    !             LOPartAmp(dn)=BornAmps(1)%Result+BornAmps(2)%Result*( couplZDD*propZ*couplZLL + couplGDD*couplGLL )
                  call Error("Zdecays.gt.10 not yet implemented")
              endif
          endif

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

       Am_up(i1,j1,c1,c2) = (0.0_dp,0.0_dp); Am_dn(i1,j1,c1,c2) = (0.0_dp,0.0_dp)

      do i2=-1,1,N2jump!  Z boson
      do i3=-1,1,2
      do i4=-1,1,2
      do i5=-1,1,2
         Am_up(i1,j1,c1,c2) = Am_up(i1,j1,c1,c2) + Bm_up(i1,c1,i2,i3,i4,i5)*conjg(Bm_up(j1,c2,i2,i3,i4,i5))
         Am_dn(i1,j1,c1,c2) = Am_dn(i1,j1,c1,c2) + Bm_dn(i1,c1,i2,i3,i4,i5)*conjg(Bm_dn(j1,c2,i2,i3,i4,i5))
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
              cres(1) = cres(1) + C(i3,i4)*Am_up(i1,i2,i3,i4)*HH(i1,i2)
              cres(2) = cres(2) + C(i3,i4)*Am_dn(i1,i2,i3,i4)*HH(i1,i2)
          enddo
          enddo
          enddo
          enddo

          res = one/(pij2 -mij**2)*real(cres,dp)


!---------- this part requires getting the proper combination of the
!---------- phase-space weights

!-------- account for all the weights, change the sign
          res = (-1d0)*res*PSWgt1*PSwgt2*Wgt_ext

!------- fill in histograms --- this  is not going to

         do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),res(1)+res(2))
         enddo


         endif   ! endif for passed cuts

       end subroutine dipff


       subroutine dipfi(n,i,j,a,mi,mj,ma,fl1,fl2,p,res)
       implicit none
       integer, intent(in) :: n,i,j,a
       real(dp), intent(in) :: mi,mj,ma
       character, intent(in) :: fl1*2,fl2*2
       real(dp), intent(in) :: p(4,6)
       real(dp), intent(out) :: res(2)
       complex(dp) :: cres(2)
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
       complex(8)        :: Bm(-1:1,1:2,-1:1,-1:1,-1:1,-1:1), Bm_up(-1:1,1:2,-1:1,-1:1,-1:1,-1:1), Bm_dn(-1:1,1:2,-1:1,-1:1,-1:1,-1:1)
       complex(8)        :: Am_up(-1:1,-1:1,1:2,1:2),Am_dn(-1:1,-1:1,1:2,1:2)
       complex(dp) :: POL1(-1:1,4),propZ
       real(dp) ::  pjetout(4,5), weight,couplZLL,couplGLL,couplZUU,couplGUU,couplZDD,couplGDD
       integer :: Njet
       logical, save :: first_time = .true.
       logical :: Not_Passed_Cuts
       integer :: NBin(1:NumHistograms), Nmax(5), Nhisto,N2jump
       real(dp) :: MomDK(1:4,1:8),PSWgt1,PSWgt2,PSWgt3,MZ_Inv

       res = zero
       cres = (0d0,0d0)
       Bm =   (0d0,0d0); Bm_up =   (0d0,0d0); Bm_dn =   (0d0,0d0)
       weight = zero


       if( first_time ) then
             call InitTrees(4,0,2,TreeAmpsDip,NumBoson=1)
             call InitProcess_TbTQbQZ(ExtParticles(1:5))
             TreeAmpsDip(1)%PartRef(1:5) = (/1,5,2,3,4/)
             TreeAmpsDip(2)%PartRef(1:5) = (/1,2,3,5,4/)
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
          fl(2) = 'gm'   ! z boson
          fl(3) = 'qm'
          fl(4) = 'qu'
          fl(5) = 'qu'


          q(:,1) = p(:,1)
          q(:,2) = p(:,2)
          q(:,3) = p(:,3)
          q(:,4) = p(:,4)
          q(:,5) = p(:,5)

          q(:,j) = pijt(:)
          q(:,a) = -pat(:)
          pos = j


!----- now generate the top decay products

       if (TopDecays.ge.1) then
         call EvalPhasespace_TopDecay(q(1:4,1),yRnDk(1:4),.false.,MomDK(1:4,1:3),PSWgt1)
         call EvalPhasespace_TopDecay(q(1:4,3),yRnDk(5:8),.false.,MomDK(1:4,4:6),PSWgt2)
       elseif(TopDecays.eq.0) then
         MomDK = 0d0
         PSWgt1 = one
         PSWgt2 = one
       else
          print *, 'topdecays are not defined'
          pause
       endif

       if( ZDecays.le.10 ) then  ! decaying on-shell Z
            MZ_Inv = m_Z
       elseif( ZDecays.gt.10 ) then  ! decaying off-shell Z
            call Error("need to implement phase space for off-shell Z's")
            ! need to think about threshold cut: EHat.le.2d0*m_Top+M_Z   when z is off-shell!
       endif
       IF( ZDECAYS.GT.0 ) THEN
          call EvalPhasespace_ZDecay(MZ_Inv,q(1:4,2),yRnDk(9:10),MomDK(1:4,7:8),PSWgt3)
          PSWgt1 = PSWgt1 * PSWgt3
       ENDIF





!----------------- ini   ini     final   top   top, real


  call Kinematics_TTBARZ(0,(/-q(1:4,4),-q(1:4,5),q(1:4,2),q(1:4,1),q(1:4,3), Momzero,MomDK(1:4,1:8)/), (/4,5,3,1,2,0,7,8,9,10,11,12,13,14/), Not_Passed_Cuts,NBin(1:NumHistograms)   )



     if(Not_Passed_Cuts.eq..false.) then


        call cc_qq_ttgamg_soft(n,C)

        Nmax = 1
        if (TopDecays.ge.1) then ! Top decays
              Nmax(3) = -1
              Nmax(1) = -1
              if (pos.eq.4.or.pos.eq.5) then  ! this is needed because of swappping
                  Nmax(pos) = -1               ! helicity index below
                  Nmax(1) = 1
              endif
        endif
        N2Jump = 1
        if (ZDecays.ge.1 .or. ZDecays.eq.-2) then ! Z decays
            N2jump=2
        endif

       do i1=-1,Nmax(1),2
         do i2 = -1,Nmax(2),N2Jump! Z boson
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

    call SetPolarization((/q(1:4,1),q(1:4,3),q(1:4,4),q(1:4,5),q(1:4,2)/),momDK(1:4,1:8),(/hel(1),hel(3),hel(4),hel(5),hel(2)/),ExtParticles(1:5))
    if( ZDecays.gt.0 ) call ZGamLCoupl(1,hel(2),couplZLL,couplGLL)  ! charged lept

      if (pos.eq.1) then
       if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(1)%pOL(1:4)
       if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(1)%pOL(1:4)
      endif

      if (pos.eq.2) then
         print *, 'pos eq 2, something is wrong'
         stop
      endif


      if (pos.eq.3) then
       if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(2)%pOL(1:4)
       if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(2)%pOL(1:4)
      endif

      if (pos.eq.4) then
       if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(3)%pOL(1:4)
       if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(3)%pOL(1:4)
      endif

      if (pos.eq.5) then
       if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(4)%pOL(1:4)
       if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(4)%pOL(1:4)
      endif


      call ZGamQcoupl(Up_,ExtParticles(3)%Helicity,couplZUU,couplGUU)
      call ZGamQcoupl(Dn_,ExtParticles(3)%Helicity,couplZDD,couplGDD)
      couplZQQ_left_dyn=one
      couplZQQ_right_dyn=one

      if ( ZDecays .lt. 10) then
          propZ = (1d0,0d0)/dsqrt(2d0*Ga_Zexp*m_Z)  * MZ_Inv**2
      elseif (ZDecays .gt. 10) then 
          propZ=cone/(MZ_Inv**2-m_Z**2+ci*Ga_ZExp*m_Z)   * MZ_Inv**2
      endif

      do i6 = 1,2
          call EvalTree2(TreeAmpsDip(i6),Bm(i1,i6,i2,i3,i4,i5))
          if( i6.eq.1 ) then 
              Bm_up(i1,i6,i2,i3,i4,i5)=Bm(i1,i6,i2,i3,i4,i5)
              Bm_dn(i1,i6,i2,i3,i4,i5)=Bm(i1,i6,i2,i3,i4,i5)
          else
              Bm_up(i1,i6,i2,i3,i4,i5)=Bm(i1,i6,i2,i3,i4,i5)*couplZUU
              Bm_dn(i1,i6,i2,i3,i4,i5)=Bm(i1,i6,i2,i3,i4,i5)*couplZDD

              if( Zdecays.gt.0 .and. Zdecays.lt.10 ) then
                  Bm_up(i1,i6,i2,i3,i4,i5)=Bm_up(i1,i6,i2,i3,i4,i5)*propZ*couplZLL
                  Bm_dn(i1,i6,i2,i3,i4,i5)=Bm_dn(i1,i6,i2,i3,i4,i5)*propZ*couplZLL
              elseif( Zdecays.gt.10 ) then
    !             LOPartAmp(up)=BornAmps(1)%Result+BornAmps(2)%Result*( couplZUU*propZ*couplZLL + couplGUU*couplGLL )
    !             LOPartAmp(dn)=BornAmps(1)%Result+BornAmps(2)%Result*( couplZDD*propZ*couplZLL + couplGDD*couplGLL )
                  call Error("Zdecays.gt.10 not yet implemented")
              endif
          endif

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

       Am_up(i1,j1,c1,c2) = (0.0_dp,0.0_dp); Am_dn(i1,j1,c1,c2) = (0.0_dp,0.0_dp)

      do i2=-1,1,N2jump!  Z boson
      do i3=-1,1,2
      do i4=-1,1,2
      do i5=-1,1,2
         Am_up(i1,j1,c1,c2) = Am_up(i1,j1,c1,c2) + Bm_up(i1,c1,i2,i3,i4,i5)*conjg(Bm_up(j1,c2,i2,i3,i4,i5))
         Am_dn(i1,j1,c1,c2) = Am_dn(i1,j1,c1,c2) + Bm_dn(i1,c1,i2,i3,i4,i5)*conjg(Bm_dn(j1,c2,i2,i3,i4,i5))
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
              cres(1) = cres(1) + C(i3,i4)*Am_up(i1,i2,i3,i4)*HH(i1,i2)
              cres(2) = cres(2) + C(i3,i4)*Am_dn(i1,i2,i3,i4)*HH(i1,i2)
             enddo
           enddo
           enddo
          enddo

         res = one/xija/(pij2-mij**2)*real(cres,dp)


!-------- account for all the weights, change the sing
         res = (-1d0)*res*PSWgt1*PSwgt2*Wgt_ext


!------- fill in histograms

         do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),res(1)+res(2))
         enddo


         endif ! for passed cuts

        end subroutine dipfi


       subroutine dipif(n,i,a,j,mi,ma,mj,fl1,fl2,p,res)
       implicit none
       integer, intent(in) :: n,i,j,a
       real(dp), intent(in) :: mi,mj,ma
       character, intent(in) :: fl1*2,fl2*2
       real(dp), intent(in) :: p(4,6)
       real(dp), intent(out) :: res(2)
       complex(dp) :: cres(2)
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
       complex(8)        :: Bm(-1:1,1:2,-1:1,-1:1,-1:1,-1:1), Bm_up(-1:1,1:2,-1:1,-1:1,-1:1,-1:1), Bm_dn(-1:1,1:2,-1:1,-1:1,-1:1,-1:1)
       complex(8)        :: Am_up(-1:1,-1:1,1:2,1:2),Am_dn(-1:1,-1:1,1:2,1:2)
       complex(dp) :: POL1(-1:1,4),propZ
       real(dp) ::  pjetout(4,5), weight,couplZLL,couplGLL,couplZUU,couplGUU,couplZDD,couplGDD
       integer :: Njet
       logical, save :: first_time = .true.
       logical :: Not_Passed_Cuts
       integer :: NBin(1:NumHistograms), Nmax(5), Nhisto,N2jump
       real(dp) :: MomDK(1:4,1:8),PSWgt1,PSWgt2,PSWgt3,MZ_Inv

       res = zero
       cres = (0d0,0d0)
       Bm =   (0d0,0d0); Bm_up =   (0d0,0d0); Bm_dn =   (0d0,0d0)
       weight = zero


       if( first_time ) then
             call InitTrees(4,0,2,TreeAmpsDip,NumBoson=1)
             call InitProcess_TbTQbQZ(ExtParticles(1:5))
             TreeAmpsDip(1)%PartRef(1:5) = (/1,5,2,3,4/)
             TreeAmpsDip(2)%PartRef(1:5) = (/1,2,3,5,4/)
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
          fl(4) = 'qu'
          fl(5) = 'qu'


          q(:,1) = p(:,1)
          q(:,2) = p(:,2)
          q(:,3) = p(:,3)
          q(:,4) = p(:,4)
          q(:,5) = p(:,5)

          q(:,a) = -pait(:)

          q(:,j) = pjt(:)
          pos = a



       if (TopDecays.ge.1) then
         call EvalPhasespace_TopDecay(q(1:4,1),yRnDk(1:4),.false.,MomDK(1:4,1:3),PSWgt1)
         call EvalPhasespace_TopDecay(q(1:4,3),yRnDk(5:8),.false.,MomDK(1:4,4:6),PSWgt2)
       elseif(TopDecays.eq.0) then
         MomDK = 0d0
         PSWgt1 = one
         PSWgt2 = one
       else
          print *, 'topdecays are not defined'
          pause
       endif

       if( ZDecays.le.10 ) then  ! decaying on-shell Z
            MZ_Inv = m_Z
       elseif( ZDecays.gt.10 ) then  ! decaying off-shell Z
            call Error("need to implement phase space for off-shell Z's")
            ! need to think about threshold cut: EHat.le.2d0*m_Top+M_Z   when z is off-shell!
       endif
       IF( ZDECAYS.GT.0 ) THEN
          call EvalPhasespace_ZDecay(MZ_Inv,q(1:4,2),yRnDk(9:10),MomDK(1:4,7:8),PSWgt3)
          PSWgt1 = PSWgt1 * PSWgt3
       ENDIF





!----------------- ini   ini     final   top   top, real


  call Kinematics_TTBARZ(0,(/-q(1:4,4),-q(1:4,5),q(1:4,2),q(1:4,1),q(1:4,3), Momzero,MomDK(1:4,1:8)/), (/4,5,3,1,2,0,7,8,9,10,11,12,13,14/), Not_Passed_Cuts,NBin(1:NumHistograms)   )



     if(Not_Passed_Cuts.eq..false.) then

        call cc_qq_ttgamg_soft(n,C)

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
        N2Jump = 1
        if (ZDecays.ge.1 .or. ZDecays.eq.-2) then ! Z decays
            N2jump=2
        endif


! print *, "Dip",n
! print *, "TB",TreeAmpsDip(1)%QUARKS(1)%Mom(1:4)
! print *, "T ",TreeAmpsDip(1)%QUARKS(2)%Mom(1:4)
! print *, "QB",TreeAmpsDip(1)%QUARKS(3)%Mom(1:4)
! print *, "Q ",TreeAmpsDip(1)%QUARKS(4)%Mom(1:4)
! print *, "Z ",TreeAmpsDip(1)%Boson%Mom(1:4)
! pause

       do i1=-1,Nmax(1),2
         do i2 = -1,Nmax(2),N2Jump! Z boson
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



       call SetPolarization((/q(1:4,1),q(1:4,3),q(1:4,4),q(1:4,5),q(1:4,2)/),momDK(1:4,1:8),(/hel(1),hel(3),hel(4),hel(5),hel(2)/),ExtParticles(1:5))
       if( ZDecays.gt.0 ) call ZGamLCoupl(1,hel(2),couplZLL,couplGLL)  ! charged lept



         if (pos.eq.1) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(1)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(1)%pOL(1:4)
         endif


       if (pos.eq.2) then
          print *, 'priehali, pos = 2'
          stop
      endif



       if (pos.eq.3) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(2)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(2)%pOL(1:4)
      endif


      if (pos.eq.4) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(3)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(3)%pOL(1:4)
      endif

      if (pos.eq.5) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(4)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(4)%pOL(1:4)
      endif

      call ZGamQcoupl(Up_,ExtParticles(3)%Helicity,couplZUU,couplGUU)
      call ZGamQcoupl(Dn_,ExtParticles(3)%Helicity,couplZDD,couplGDD)
      couplZQQ_left_dyn=one
      couplZQQ_right_dyn=one

      if ( ZDecays .lt. 10) then
          propZ = (1d0,0d0)/dsqrt(2d0*Ga_Zexp*m_Z)    * MZ_Inv**2
      elseif (ZDecays .gt. 10) then 
          propZ=cone/(MZ_Inv**2-m_Z**2+ci*Ga_ZExp*m_Z)   * MZ_Inv**2
      endif

      do i6 = 1,2
          call EvalTree2(TreeAmpsDip(i6),Bm(i1,i6,i2,i3,i4,i5))
          if( i6.eq.1 ) then 
              Bm_up(i1,i6,i2,i3,i4,i5)=Bm(i1,i6,i2,i3,i4,i5)
              Bm_dn(i1,i6,i2,i3,i4,i5)=Bm(i1,i6,i2,i3,i4,i5)
          else
              Bm_up(i1,i6,i2,i3,i4,i5)=Bm(i1,i6,i2,i3,i4,i5)*couplZUU
              Bm_dn(i1,i6,i2,i3,i4,i5)=Bm(i1,i6,i2,i3,i4,i5)*couplZDD

              if( Zdecays.gt.0 .and. Zdecays.lt.10 ) then
                  Bm_up(i1,i6,i2,i3,i4,i5)=Bm_up(i1,i6,i2,i3,i4,i5)*propZ*couplZLL
                  Bm_dn(i1,i6,i2,i3,i4,i5)=Bm_dn(i1,i6,i2,i3,i4,i5)*propZ*couplZLL
              elseif( Zdecays.gt.10 ) then
    !             LOPartAmp(up)=BornAmps(1)%Result+BornAmps(2)%Result*( couplZUU*propZ*couplZLL + couplGUU*couplGLL )
    !             LOPartAmp(dn)=BornAmps(1)%Result+BornAmps(2)%Result*( couplZDD*propZ*couplZLL + couplGDD*couplGLL )
                  call Error("Zdecays.gt.10 not yet implemented")
              endif
          endif
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

       Am_up(i1,j1,c1,c2) = (0.0_dp,0.0_dp)
       Am_dn(i1,j1,c1,c2) = (0.0_dp,0.0_dp)

      do i2=-1,1,N2jump!  Z boson
      do i3=-1,1,2
      do i4=-1,1,2
      do i5=-1,1,2
         Am_up(i1,j1,c1,c2) = Am_up(i1,j1,c1,c2) +  Bm_up(i1,c1,i2,i3,i4,i5)*conjg(Bm_up(j1,c2,i2,i3,i4,i5))
         Am_dn(i1,j1,c1,c2) = Am_dn(i1,j1,c1,c2) +  Bm_dn(i1,c1,i2,i3,i4,i5)*conjg(Bm_dn(j1,c2,i2,i3,i4,i5))
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
      diag=  two/(two-xija - zj)-one - xija
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
           do i3=1,2
             do i4=1,2
                cres(1) = cres(1) + C(i3,i4)*Am_up(i1,i2,i3,i4)*HH(i1,i2)
                cres(2) = cres(2) + C(i3,i4)*Am_dn(i1,i2,i3,i4)*HH(i1,i2)
            enddo
           enddo
          enddo
         enddo

         res = one/xija/two/scr(pa,pi)*real(cres,dp)

!-------- account for all the weights, change the sing
         res = (-1d0)*res*PSWgt1*PSwgt2*Wgt_ext

!------- fill in histograms

         do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),res(1)+res(2))
         enddo

         endif ! -- ! passed cuts

        end subroutine dipif



       subroutine dipii(n,i,a,b,mi,ma,mb,fl1,fl2,p,res)
       implicit none
       integer, intent(in) :: n,i,a,b
       real(dp), intent(in) :: mi,ma,mb
       character, intent(in) :: fl1*2,fl2*2
       real(dp), intent(in) :: p(4,6)
       real(dp), intent(out) :: res(2)
       complex(dp) :: cres(2)
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
       complex(8)        :: Bm(-1:1,1:2,-1:1,-1:1,-1:1,-1:1), Bm_up(-1:1,1:2,-1:1,-1:1,-1:1,-1:1), Bm_dn(-1:1,1:2,-1:1,-1:1,-1:1,-1:1)
       complex(8)        :: Am_up(-1:1,-1:1,1:2,1:2),Am_dn(-1:1,-1:1,1:2,1:2)
       complex(8)        :: ResAmpsDip_mt(-1:1,1:2)
       complex(dp) :: POL1(-1:1,4),propZ
       real(dp) ::  pjetout(4,5), weight,couplZLL,couplGLL,couplZUU,couplGUU,couplZDD,couplGDD
       integer :: Njet
       logical, save :: first_time = .true.
       logical :: Not_Passed_Cuts
       integer :: NBin(1:NumHistograms), Nmax(5), Nhisto,N2jump
       real(dp) :: MomDK(1:4,1:8),PSWgt1,PSWgt2,PSWgt3,MZ_Inv


       res = zero
       cres = (0d0,0d0)
       Bm =   (0d0,0d0); Bm_up =   (0d0,0d0); Bm_dn =   (0d0,0d0)


       if( first_time ) then
             call InitTrees(4,0,2,TreeAmpsDip,NumBoson=1)
             call InitProcess_TbTQbQZ(ExtParticles(1:5))
             TreeAmpsDip(1)%PartRef(1:5) = (/1,5,2,3,4/)
             TreeAmpsDip(2)%PartRef(1:5) = (/1,2,3,5,4/)
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
          fl(4) = 'qu'
          fl(5) = 'qu'


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
         call EvalPhasespace_TopDecay(q(1:4,1),yRnDk(1:4),.false.,MomDK(1:4,1:3),PSWgt1)
         call EvalPhasespace_TopDecay(q(1:4,3),yRnDk(5:8),.false.,MomDK(1:4,4:6),PSWgt2)
       elseif(TopDecays.eq.0) then
         MomDK = 0d0
         PSWgt1 = one
         PSWgt2 = one
       else
          print *, 'topdecays are not defined'
          pause
       endif

       if( ZDecays.le.10 ) then  ! decaying on-shell Z
            MZ_Inv = m_Z
       elseif( ZDecays.gt.10 ) then  ! decaying off-shell Z
            call Error("need to implement phase space for off-shell Z's")
            ! need to think about threshold cut: EHat.le.2d0*m_Top+M_Z   when z is off-shell!
       endif
       IF( ZDECAYS.GT.0 ) THEN
          call EvalPhasespace_ZDecay(MZ_Inv,q(1:4,2),yRnDk(9:10),MomDK(1:4,7:8),PSWgt3)
          PSWgt1 = PSWgt1 * PSWgt3
       ENDIF





!----------------- ini   ini     final   top   top, real


  call Kinematics_TTBARZ(0,(/-q(1:4,4),-q(1:4,5),q(1:4,2),q(1:4,1),q(1:4,3), Momzero,MomDK(1:4,1:8)/), (/4,5,3,1,2,0,7,8,9,10,11,12,13,14/), Not_Passed_Cuts,NBin(1:NumHistograms)   )


     if(Not_Passed_Cuts.eq..false.) then

        call cc_qq_ttgamg_soft(n,C)

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
        N2Jump = 1
        if (ZDecays.ge.1 .or. ZDecays.eq.-2) then ! Z decays
            N2jump=2
        endif

! print *, "TB",dreal( TreeAmpsDip(1)%QUARKS(1)%Mom(1:4) )
! print *, "T ",dreal( TreeAmpsDip(1)%QUARKS(2)%Mom(1:4) )
! print *, "QB",dreal( TreeAmpsDip(1)%QUARKS(3)%Mom(1:4) )
! print *, "Q ",dreal( TreeAmpsDip(1)%QUARKS(4)%Mom(1:4) )
! print *, "Z ",dreal( TreeAmpsDip(1)%Boson%Mom(1:4) )
! pause

       do i1=-1,Nmax(1),2
         do i2 = -1,Nmax(2),N2Jump! Z boson
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

    call SetPolarization((/q(1:4,1),q(1:4,3),q(1:4,4),q(1:4,5),q(1:4,2)/),momDK(1:4,1:8),(/hel(1),hel(3),hel(4),hel(5),hel(2)/),ExtParticles(1:5))
    if( ZDecays.gt.0 ) call ZGamLCoupl(1,hel(2),couplZLL,couplGLL)  ! charged lept

! print *, hel(4),hel(5)
! print *, TreeAmpsDip(1)%QUARKS(3)%pol(1:4)
! print *, TreeAmpsDip(1)%QUARKS(4)%pol(1:4)
! ! print *, TreeAmpsDip(1)%boson%pol(1:4)
! pause


     if (pos.eq.1) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(1)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(1)%pOL(1:4)
     endif


       if (pos.eq.3) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(2)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(2)%pOL(1:4)
      endif


       if (pos.eq.2) then
          print *, 'pos = 2, something is wrong'
          stop
      endif

      if (pos.eq.4) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(3)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(3)%pOL(1:4)
      endif

      if (pos.eq.5) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%QUARKS(4)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%QUARKS(4)%pOL(1:4)
      endif


      call ZGamQcoupl(Up_,ExtParticles(3)%Helicity,couplZUU,couplGUU)
      call ZGamQcoupl(Dn_,ExtParticles(3)%Helicity,couplZDD,couplGDD)
      couplZQQ_left_dyn=one
      couplZQQ_right_dyn=one

      if ( ZDecays .lt. 10) then
          propZ = (1d0,0d0)/dsqrt(2d0*Ga_Zexp*m_Z)    * MZ_Inv**2
      elseif (ZDecays .gt. 10) then 
          propZ=cone/(MZ_Inv**2-m_Z**2+ci*Ga_ZExp*m_Z)   * MZ_Inv**2
      endif

      do i6 = 1,2
          call EvalTree2(TreeAmpsDip(i6),Bm(i1,i6,i2,i3,i4,i5))
          if( i6.eq.1 ) then 
              Bm_up(i1,i6,i2,i3,i4,i5)=Bm(i1,i6,i2,i3,i4,i5)
              Bm_dn(i1,i6,i2,i3,i4,i5)=Bm(i1,i6,i2,i3,i4,i5)
          else
              Bm_up(i1,i6,i2,i3,i4,i5)=Bm(i1,i6,i2,i3,i4,i5)*couplZUU
              Bm_dn(i1,i6,i2,i3,i4,i5)=Bm(i1,i6,i2,i3,i4,i5)*couplZDD

              if( Zdecays.gt.0 .and. Zdecays.lt.10 ) then
                  Bm_up(i1,i6,i2,i3,i4,i5)=Bm_up(i1,i6,i2,i3,i4,i5)*propZ*couplZLL
                  Bm_dn(i1,i6,i2,i3,i4,i5)=Bm_dn(i1,i6,i2,i3,i4,i5)*propZ*couplZLL
              elseif( Zdecays.gt.10 ) then
    !             LOPartAmp(up)=BornAmps(1)%Result+BornAmps(2)%Result*( couplZUU*propZ*couplZLL + couplGUU*couplGLL )
    !             LOPartAmp(dn)=BornAmps(1)%Result+BornAmps(2)%Result*( couplZDD*propZ*couplZLL + couplGDD*couplGLL )
                  call Error("Zdecays.gt.10 not yet implemented")
              endif
          endif

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

       Am_up(i1,j1,c1,c2) = (0.0_dp,0.0_dp); Am_dn(i1,j1,c1,c2) = (0.0_dp,0.0_dp)

      do i2=-1,1,N2jump!  Z boson
      do i3=-1,1,2
      do i4=-1,1,2
      do i5=-1,1,2
         Am_up(i1,j1,c1,c2) = Am_up(i1,j1,c1,c2) + Bm_up(i1,c1,i2,i3,i4,i5)*conjg(Bm_up(j1,c2,i2,i3,i4,i5))
         Am_dn(i1,j1,c1,c2) = Am_dn(i1,j1,c1,c2) + Bm_dn(i1,c1,i2,i3,i4,i5)*conjg(Bm_dn(j1,c2,i2,i3,i4,i5))
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


      if ( (fl1.eq.'gl').and.(fl2.eq.'qu') ) then

      diag= two/(one-xiab) - (one+xiab)
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
           do i3=1,2
             do i4=1,2
              cres(1) = cres(1) + C(i3,i4)*Am_up(i1,i2,i3,i4)*HH(i1,i2)
              cres(2) = cres(2) + C(i3,i4)*Am_dn(i1,i2,i3,i4)*HH(i1,i2)
            enddo
          enddo
         enddo
        enddo


         res = one/xiab/two/scr(pa,pi)*real(cres,dp)

!-------- account for all the weights, change the sing
         res = (-1d0)*res*PSWgt1*PSwgt2*Wgt_ext

!------- fill in histograms

         do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),res(1)+res(2))
         enddo

         endif !  for passed cuts

        end subroutine dipii







      subroutine cc_qq_ttgamg_soft(n,C)  
      integer, intent(in) :: n           
      real(dp), intent(out) :: C(2,2)    

      C = 0.0_dp


      if(n.eq.0) then! 0= tree level correlation
      C(1,1)=8.0_dp
      C(1,2)=8.0_dp
      C(2,1)=8.0_dp
      C(2,2)=8.0_dp
      endif
      if(n.eq.1) then
      C(1,1)=-8.0_dp/3.0_dp
      C(1,2)=-8.0_dp/3.0_dp
      C(2,1)=-8.0_dp/3.0_dp
      C(2,2)=-8.0_dp/3.0_dp
      endif
      if(n.eq.2) then
      C(1,1)=16.0_dp/3.0_dp
      C(1,2)=16.0_dp/3.0_dp
      C(2,1)=16.0_dp/3.0_dp
      C(2,2)=16.0_dp/3.0_dp
      endif
      if(n.eq.3) then
      C(1,1)=56.0_dp/3.0_dp
      C(1,2)=56.0_dp/3.0_dp
      C(2,1)=56.0_dp/3.0_dp
      C(2,2)=56.0_dp/3.0_dp
      endif
      if(n.eq.4) then
      C(1,1)=-8.0_dp/3.0_dp
      C(1,2)=-8.0_dp/3.0_dp
      C(2,1)=-8.0_dp/3.0_dp
      C(2,2)=-8.0_dp/3.0_dp
      endif
      if(n.eq.5) then
      C(1,1)=56.0_dp/3.0_dp
      C(1,2)=56.0_dp/3.0_dp
      C(2,1)=56.0_dp/3.0_dp
      C(2,2)=56.0_dp/3.0_dp
      endif
      if(n.eq.6) then
      C(1,1)=16.0_dp/3.0_dp
      C(1,2)=16.0_dp/3.0_dp
      C(2,1)=16.0_dp/3.0_dp
      C(2,2)=16.0_dp/3.0_dp
      endif
      if(n.eq.7) then
      C(1,1)=16.0_dp/3.0_dp
      C(1,2)=16.0_dp/3.0_dp
      C(2,1)=16.0_dp/3.0_dp
      C(2,2)=16.0_dp/3.0_dp
      endif
      if(n.eq.8) then
      C(1,1)=56.0_dp/3.0_dp
      C(1,2)=56.0_dp/3.0_dp
      C(2,1)=56.0_dp/3.0_dp
      C(2,2)=56.0_dp/3.0_dp
      endif
      if(n.eq.9) then
      C(1,1)=-8.0_dp/3.0_dp
      C(1,2)=-8.0_dp/3.0_dp
      C(2,1)=-8.0_dp/3.0_dp
      C(2,2)=-8.0_dp/3.0_dp
      endif
      if(n.eq.10) then
      C(1,1)=56.0_dp/3.0_dp
      C(1,2)=56.0_dp/3.0_dp
      C(2,1)=56.0_dp/3.0_dp
      C(2,2)=56.0_dp/3.0_dp
      endif
      if(n.eq.11) then
      C(1,1)=16.0_dp/3.0_dp
      C(1,2)=16.0_dp/3.0_dp
      C(2,1)=16.0_dp/3.0_dp
      C(2,2)=16.0_dp/3.0_dp
      endif
      if(n.eq.12) then
      C(1,1)=-8.0_dp/3.0_dp
      C(1,2)=-8.0_dp/3.0_dp
      C(2,1)=-8.0_dp/3.0_dp
      C(2,2)=-8.0_dp/3.0_dp
      endif




      end subroutine





SUBROUTINE SetPolarization(Mom,MomDK,Hel,ExtParticles)
use ModMisc
use ModProcess
use ModTopDecay
use ModZDecay
implicit none
type(Particle) :: ExtParticles(1:5)
real(8) :: Mom(1:4,1:5),MomDK(1:4,1:8)
integer :: Hel(1:5)

     ExtParticles(1)%Mom(1:4) = dcmplx(Mom(1:4,1))   ! HERE WAS A BUG: this was inside the (TopDecays.ge.1) condition
     ExtParticles(2)%Mom(1:4) = dcmplx(Mom(1:4,2))
     if (TopDecays.ge.1) then
        call TopDecay(ExtParticles(1),DK_LO,MomDK(1:4,1:3))
        call TopDecay(ExtParticles(2),DK_LO,MomDK(1:4,4:6))
     else
        call vSpi(ExtParticles(1)%Mom(1:4),ExtParticles(1)%Mass,Hel(1),ExtParticles(1)%Pol(1:4))
        call ubarSpi(ExtParticles(2)%Mom(1:4),ExtParticles(2)%Mass,Hel(2),ExtParticles(2)%Pol(1:4))
    endif


     ExtParticles(5)%Mom(1:4) = dcmplx(Mom(1:4,5))
     ExtParticles(5)%Helicity = Hel(5)
     if (ZDecays.ge.1) then
          call ZDecay(ExtParticles(5),DK_LO,MomDK(1:4,7:8))
     elseif (ZDecays.eq.0) then
          call pol_massSR(ExtParticles(5)%Mom(1:4),ExtParticles(5)%Mass,ExtParticles(5)%Helicity,ExtParticles(5)%Pol(1:4))
     elseif (ZDecays.eq.-2) then
         call pol_mless(ExtParticles(5)%Mom(1:4),ExtParticles(5)%Helicity,ExtParticles(5)%Pol(1:4))
    endif




    ExtParticles(3)%Mom(1:4) = dcmplx(Mom(1:4,3))
    ExtParticles(3)%Helicity = Hel(3)
    call vSpi(ExtParticles(3)%Mom(1:4),ExtParticles(3)%Mass,Hel(3),ExtParticles(3)%Pol(1:4))

    ExtParticles(4)%Mom(1:4) = dcmplx(Mom(1:4,4))
    ExtParticles(4)%Helicity = Hel(4)
    call ubarSpi(ExtParticles(4)%Mom(1:4),ExtParticles(4)%Mass,Hel(4),ExtParticles(4)%Pol(1:4))


RETURN
END SUBROUTINE






SUBROUTINE InitProcess_TbTQbQZ(ExtParticles)
use ModProcess
implicit none
type(Particle) :: ExtParticles(:)

  ExtParticles(1)%PartType = ATop_
  ExtParticles(1)%ExtRef   = 1
  ExtParticles(1)%Mass = m_Top
  ExtParticles(1)%Mass2= ExtParticles(1)%Mass**2

  ExtParticles(2)%PartType = Top_
  ExtParticles(2)%ExtRef   = 2
  ExtParticles(2)%Mass = m_Top
  ExtParticles(2)%Mass2= ExtParticles(2)%Mass**2

  ExtParticles(3)%PartType = AStr_
  ExtParticles(3)%ExtRef   = 3
  ExtParticles(3)%Mass = 0d0
  ExtParticles(3)%Mass2= 0d0

  ExtParticles(4)%PartType = Str_
  ExtParticles(4)%ExtRef   = 4
  ExtParticles(4)%Mass = 0d0
  ExtParticles(4)%Mass2= 0d0

  ExtParticles(5)%PartType = Z0_
  if( Process.eq.86 ) ExtParticles(5)%PartType = Pho_
  ExtParticles(5)%ExtRef   = 5
  ExtParticles(5)%Mass = m_Z
  ExtParticles(5)%Mass2= ExtParticles(5)%Mass**2

RETURN
END SUBROUTINE



       end module
