! this is the file to subtract dipoles for qq->tt+g+gamma amplitudes
      module ModDipoles_QGTTBQZ
      use ModAmplitudes
      use ModProcess
      use ModParameters
      use ModMisc
      use ModKinematics
      use ModIntDipoles
      use ModZDecay
      implicit none
      private

      integer, parameter  :: dp = selected_real_kind(15)
      real(dp), private :: yRnDK(1:10), Wgt_ext(2)
      real(dp), parameter, private  :: Momzero(1:4)=0d0

      public :: EvalDipoles_QGTTBQZ

      logical, parameter :: invert_alphaCut = .false.

      contains

!  we have a list of dipoles that we need to go through
!  We label things as 0-> bar t(p1)+ gamma(2) + t(3) + q(p4)+q(p5)+g(p6)
!  and we  assume that quark in the initial state has momentum p4 and
! the  gluon in the initial state has momentum  p6;
! final state (massless) quark has momentum p5


!--- I pass TWO external weights to the subroutine, one for `up' quarks
!-------------------------------------------------, the other for 'dn' quarks

      subroutine EvalDipoles_QGTTBQZ(p,yRnDk1,Wgt,sum_dip)
      real(dp), intent(out) ::  sum_dip(2)
      real(dp), intent(in) :: p(4,6)
      real(dp), intent(in) :: yRnDK1(1:10), Wgt(2)
      integer, parameter :: ndip = 2
      integer, parameter :: in1 = 4
      integer, parameter :: in2 = 6
      real(dp) :: res(2)
      integer ::  dip(ndip,3)
      data dip(1,1)/5/, dip(1,2)/4/, dip(1,3)/6/
      data dip(2,1)/5/, dip(2,2)/6/, dip(2,3)/4/
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

      if ( (j.eq.in1.and.k.eq.in2).or.(j.eq.in2.and.k.eq.in1)) then
         if (n.eq.1) then
!------------------------------- here quark emits
        call dipii_q(n,i,j,k,mass(i),mass(j),mass(k),fl(i),fl(j),p,res)
         elseif(n.eq.2) then
!------------------------------- here gluon emits
        call dipii_g(n,i,j,k,mass(i),mass(j),mass(k),fl(i),fl(j),p,res)
         endif
      endif

        sum_dip = sum_dip + res

print *, "dip",n,res

        enddo

      end subroutine






!      dipole subroutines

       subroutine dipii_q(n,i,a,b,mi,ma,mb,fl1,fl2,p,res)
       implicit none
       integer, intent(in) :: n,i,a,b
       real(dp), intent(in) :: mi,ma,mb
       character, intent(in) :: fl1*2,fl2*2
       real(dp), intent(in) :: p(4,6)
       real(dp), intent(out) :: res(2)
       complex(dp) :: cres(2)
       real(dp) :: mij
       real(dp) ::   C(2,2)
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
       real(dp), parameter :: CF=4.0_dp/3.0_dp
       complex(dp) :: xp, xm
       complex(dp) :: HH(-1:1,-1:1)
       type(Particle),save :: ExtParticles(1:5)
       type(TreeProcess),save :: TreeAmpsDip(1:2)
       integer :: c1, c2, j1, in1, in2
       complex(8)        :: Bm(-1:1,1:2,-1:1,-1:1,-1:1,-1:1)
       complex(8)        :: Am(-1:1,-1:1,1:2,1:2)
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
       Bm = (0d0,0d0)


       if( first_time ) then
             call InitTrees(2,2,2,TreeAmpsDip,NumBoson=1)
             call InitProcess_TbTGGZ(ExtParticles(1:5))
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
          q(:,5) = -pb(:)         !-- this is specific for this dipole
          pos = a

          in1 = a
          in2 = 5


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

       if (TopDecays.ge.1) then
         call EvalPhasespace_TopDecay(q(1:4,1),yRnDk(1:4),.false., MomDK(1:4,1:3),PSWgt1)
         call EvalPhasespace_TopDecay(q(1:4,3),yRnDk(5:8),.false., MomDK(1:4,4:6),PSWgt2)
       elseif(TopDecays.eq.0) then
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
       IF( ZDECAYS.NE.0 ) THEN
          call EvalPhasespace_ZDecay(MZ_Inv,q(1:4,2),yRnDk(9:10),MomDK(1:4,7:8),PSWgt3)
          PSWgt1 = PSWgt1 * PSWgt3
       ENDIF

!-----------------     initial   initial       final  top      top


  call Kinematics_TTBARZ(0,(/-q(1:4,4),-q(1:4,5),q(1:4,2),q(1:4,1),q(1:4,3), Momzero,MomDK(1:4,1:8)/), (/4,5,3,1,2,0,7,8,9,10,11,12,13,14/), Not_Passed_Cuts,NBin(1:NumHistograms)   )


     if(Not_Passed_Cuts.eq..false.) then

        call cc_gg_ttgam(C)


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

        N2Jump = 1
        if (ZDecays.ge.1) then ! Z decays
            N2jump=2
        endif

       do i1=-1,Nmax(1),2
          do i2 = -1,Nmax(2),N2jump!  Z boson
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
  call SetPolarization_GG((/q(1:4,1),q(1:4,3),q(1:4,4),q(1:4,5),q(1:4,2)/),momDK(1:4,1:8),(/hel(1),hel(3),hel(4),hel(5),hel(2)/),ExtParticles(1:5))

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
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%GLUONS(2)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%GLUONS(2)%pOL(1:4)
      endif

      if (pos.eq.5) then
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%GLUONS(3)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%GLUONS(3)%pOL(1:4)
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


      if ( (fl1.eq.'qu').and.(fl2.eq.'qu') ) then

       paux = pi - scr(pi,pa)/scr(pi,pb)*pb

       diag= xiab
       offdiag = two*(one-xiab)/xiab*scr(pa,pb)/scr(pi,pa)/scr(pi,pb)

       xm = scc(POL1(-1,:),paux)
       xp = scc(POL1(+1,:),paux)

      HH(-1,1)=offdiag*conjg(xm)*xp
      HH(1,-1)=offdiag*conjg(xp)*xm
      HH(-1,-1)=diag + offdiag*conjg(xm)*xm
      HH(1,1) = diag + offdiag*conjg(xp)*xp

      HH = HH

      endif



      do i1 = -1,1,2
         do i2 = -1,1,2
           do i3=1,2
             do i4=1,2

       cres(1) = cres(1) + C(i3,i4)*Am(i1,i2,i3,i4)*HH(i1,i2)
       cres(2) = cres(2) + C(i3,i4)*Am(i1,i2,i3,i4)*HH(i1,i2)

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

        end subroutine dipii_q



       subroutine dipii_g(n,i,a,b,mi,ma,mb,fl1,fl2,p,res)
       implicit none
       integer, intent(in) :: n,i,a,b
       real(dp), intent(in) :: mi,ma,mb
       character, intent(in) :: fl1*2,fl2*2
       real(dp), intent(in) :: p(4,6)
       real(dp), intent(out) :: res(2)
       real(dp) ::   C(2,2)
       complex(dp) :: cres(2)
       real(dp) :: mij
       real(dp) :: Cup(2,2),Cdn(2,2)
       real(dp) ::  pi(4), pj(4), pk(4), pkt(4), pij(4), qq(4), qqij(4)
       real(dp) :: q(4,6), paux(4), pijk(4), pa(4),  pat(4)
       real(dp) :: xx1, xx2, zi, zj, yijk, vijk, zim
       real(dp) :: pij2,pijk2,tvijk,q2,tpij2, qij2, zjm
       real(dp) ::  pijt(4), pija(4), pb(4),pija2, xija
       real(dp) :: pait(4), xiab
       real(dp) ::  K(4), tK(4),K2,tK2, K2tK2, k1(4)
       real(dp) :: zp, zm, xxx, viji,v
       integer :: i1, i2, i3, i4, i5, i6, hel(5), j
       character :: fl(5)*2
       integer :: pos, iTree, in1, in2
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
       real(dp), parameter :: CF=4.0_dp/3.0_dp
       integer :: Njet
       logical, save :: first_time = .true.
       logical :: Not_Passed_Cuts
       integer :: NBin(1:NumHistograms), Nmax(5), Nhisto,N2Jump
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
          fl(4) = 'gl'
          fl(5) = 'gl'


          q(:,1) = p(:,1)
          q(:,2) = p(:,2)
          q(:,3) = p(:,3)
          q(:,4) = p(:,4)
          q(:,5) = p(:,5)


          q(:,a) = -pait(:)
          q(:,5) = -pb(:)         !-- this is specific for this dipole
          pos = a
          in1 = a
          in2 = 5

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


print *, "ou1",-q(1:4,4)
print *, "ou1",-q(1:4,5)
!-----------------     initial   initial       final  top      top
  call Kinematics_TTBARZ(0,(/-q(1:4,4),-q(1:4,5),q(1:4,2),q(1:4,1),q(1:4,3), Momzero,MomDK(1:4,1:8)/), (/4,5,3,1,2,0,7,8,9,10,11,12,13,14/), Not_Passed_Cuts,NBin(1:NumHistograms)   )



     if(Not_Passed_Cuts.eq..false.) then

        call cc_qq_tt(C)

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


    call SetPolarization_QQB((/q(1:4,1),q(1:4,3),q(1:4,4),q(1:4,5),q(1:4,2)/),momDK(1:4,1:8),(/hel(1),hel(3),hel(4),hel(5),hel(2)/),ExtParticles(1:5))
    if( ZDecays.gt.0 ) call ZGamLCoupl(1,hel(2),couplZLL,couplGLL)  ! charged lept


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
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsDip(1)%gLUONS(1)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsDip(1)%gLUONS(1)%pOL(1:4)
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



      if ( (fl1.eq.'qu').and.(fl2.eq.'gl') ) then

      diag= two*CF*(one - two*xiab*(one-xiab))   ! here we change color
                                                 ! factor because we do not
                                                 ! do color averaging
      HH(-1,1)=dcmplx(0.0_dp,0.0_dp)
      HH(1,-1)=dcmplx(0.0_dp,0.0_dp)
      HH(-1,-1)=dcmplx(diag,0.0_dp)
      HH(1,1) = dcmplx(diag,0.0_dp)

      HH = HH


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

        end subroutine dipii_g


      subroutine cc_gg_ttgam(C)
      real(dp), intent(out) :: C(2,2)


       C(1,1) = 64.0_dp/3.0_dp
       C(1,2) =  - 8.0_dp/3
       C(2,1) =  - 8.0_dp/3.0_dp
       C(2,2) = 64.0_dp/3.0_dp

      end subroutine cc_gg_ttgam



      subroutine cc_qq_tt(C)
      real(dp), intent(out) :: C(2,2)

      C = 0.0_dp


       C(1,1)  =  8.0_dp
       C(1,2)  =  8.0_dp
       C(2,1)  =  8.0_dp
       C(2,2)  =  8.0_dp

      end subroutine




SUBROUTINE InitProcess_TbTGGZ(ExtParticles)
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

  ExtParticles(3)%PartType = Glu_
  ExtParticles(3)%ExtRef   = 3
  ExtParticles(3)%Mass = 0d0
  ExtParticles(3)%Mass2= 0d0

  ExtParticles(4)%PartType = Glu_
  ExtParticles(4)%ExtRef   = 4
  ExtParticles(4)%Mass = 0d0
  ExtParticles(4)%Mass2= 0d0

  ExtParticles(5)%PartType = Z0_
  ExtParticles(5)%ExtRef   = 5
  ExtParticles(5)%Mass = m_Z
  ExtParticles(5)%Mass2= ExtParticles(5)%Mass**2


RETURN
END SUBROUTINE InitProcess_TbTGGZ



SUBROUTINE SetPolarization_GG(Mom,MomDK,Hel,ExtParticles)
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
     else
          call pol_massSR(ExtParticles(5)%Mom(1:4),ExtParticles(5)%Mass,ExtParticles(5)%Helicity,ExtParticles(5)%Pol(1:4))
    endif

    ExtParticles(3)%Mom(1:4) = dcmplx(Mom(1:4,3))
    call pol_mless(ExtParticles(3)%Mom(1:4),Hel(3),ExtParticles(3)%Pol(1:4))

    ExtParticles(4)%Mom(1:4) = dcmplx(Mom(1:4,4))
    call pol_mless(ExtParticles(4)%Mom(1:4),Hel(4),ExtParticles(4)%Pol(1:4))


! ExtParticles(3)%Pol(1:4)=ExtParticles(3)%Mom(1:4); print *, "check gauge inv. in dipoles"

RETURN
END SUBROUTINE





SUBROUTINE SetPolarization_QQB(Mom,MomDK,Hel,ExtParticles)
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
