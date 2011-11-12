module ModSixFermionProcs
use ModAmplitudes
use ModProcess
use ModParameters
use ModIntDipoles
use ModKinematics
use mod_qqb_ttqqqq_dip
use mod_qqq_ttqqqq_dip
use mod_qbb_ttqqqq_dip
use ModTopdecay
implicit none


      integer, parameter  :: dp = selected_real_kind(15)
      real(dp), private :: yRnDK(1:8), Wgt_ext

contains



!------------------------------------------------------------
! the six-quark matrix element squared
! momenta: outgoing convention
! 0 -> bar t_1 + t_2 + 3_in(fl1) + 4_in(fl2) + 5_out + 6_out
!  in the all-outgoing convention
! fl1 and fl2 are flavors of colliding partons with
! mcmf convention d - 1 u -2  s -3  c -4  b-5
! and antiquarks with the minus sign
!------------------------------------------------------------
      subroutine SixQuark(pin,yRnDK1,Wgt,fl1,fl2,subtr,res,dip)
      real(dp),  intent(in) :: pin(4,6)
      integer,  intent(in) :: fl1, fl2
      logical, intent(in) :: subtr
      real(dp), intent(in) :: yRnDK1(1:8), Wgt
      real(dp), intent(out) :: res,dip
      real(dp) :: p(4,6)
      real(dp), parameter :: nf = Nf_light
      real(dp) :: res1, res2,res_subtr
      real(dp) ::  pjetout(4,6), weight
      logical :: Not_Passed_Cuts
      integer :: NBin(1:NumHistograms)
      real(dp) :: MomDK(1:4,1:6),PSWgt1,PSWgt2
      integer :: Njet, Nmax(5), Nhisto
      logical, save :: first_time = .true.


        res_subtr = zero
        yRnDK = yRnDK1               ! HERE WAS A BUG: yRnDK = yRnDK
        Wgt_ext = Wgt

       if (TopDecays.ge.1) then
         call EvalPhasespace_TopDecay(pin(1:4,1),yRnDk(1:4),.false., MomDK(1:4,1:3),PSWgt1)
         call EvalPhasespace_TopDecay(pin(1:4,2),yRnDk(5:8),.false., MomDK(1:4,4:6),PSWgt2)
       elseif(TopDecays.eq.0) then
         PSWgt1 = one
         PSWgt2 = one
       else
          print *, 'topdecays are not defined'
          pause
       endif

        dip = zero





      if (abs(fl1).ne.abs(fl2)) then   ! different quarks/antiquark/etc

!-----------------------------------------------------------------
         if ((fl1.gt.0).and.(fl2.gt.0)) then   !checked
         p = pin
         p(:,5) = pin(:,4)
         p(:,4) = pin(:,5)

!         call jetktalg(pin,6,pjetout,Njet,weight)

!-----------------     initial   initial   top      top     final      !  HERE WAS A BUG: NPart
     call Kinematics_TTBARJET(1,(/-pin(1:4,3),-pin(1:4,4),pin(1:4,5),pin(1:4,6),pin(1:4,1),pin(1:4,2)/),MomDK(1:4,1:6),Not_Passed_Cuts,NBin(1:NumHistograms))

     if(Not_Passed_Cuts.eq..false.) then
         call ttqqqq(p,res)
         res = Wgt*PSWgt1*PSWgt2*res    ! re-weight
         do NHisto=1,NumHistograms
              call intoHisto(NHisto,NBin(NHisto),res)
         enddo
     else
         res = zero
     endif

          if (subtr) then
           call qqq_ttqqqq_dip(p,yRnDK1,Wgt,res_subtr)

           res = res + res_subtr   ! for subtr, reweightning done inside
           dip = dip + res_subtr
          endif


!-----------------------------------------------------------------
         elseif((fl1.gt.0).and.(fl2.lt.0)) then  ! checked

         p = pin
         p(:,6) = pin(:,4)
         p(:,4) = pin(:,6)

!         call jetktalg(pin,6,pjetout,Njet,weight)
!         if( weight.eq.1d0) then

      !  HERE WAS A BUG: NPart
       call Kinematics_TTBARJET(1,(/-pin(1:4,3),-pin(1:4,4),pin(1:4,5),pin(1:4,6),pin(1:4,1),pin(1:4,2)/),MomDK(1:4,1:6),Not_Passed_Cuts,NBin(1:NumHistograms))
       if(Not_Passed_Cuts.eq..false.) then
            call ttqqqq(p,res)
            res = Wgt*PSWgt1*PSWgt2*res    ! re-weight
            do NHisto=1,NumHistograms
                  call intoHisto(NHisto,NBin(NHisto),res)
            enddo
       else
            res = zero
       endif

         if (subtr) then
         call qqb_ttqqqq_dip_12(p,yRnDK1,Wgt,res_subtr)
!           print *, (p(1:4,3).dot.p(1:4,4))/(p(1:4,3).dot.p(1:4,6)),(p(1:4,3).dot.p(1:4,4))/(p(1:4,3).dot.p(1:4,6))*res,res_subtr*(p(1:4,3).dot.p(1:4,4))/(p(1:4,3).dot.p(1:4,6))
!           pause
          res = res + res_subtr
          dip = dip + res_subtr
         endif


!-----------------------------------------------------------------
         elseif((fl1.lt.0).and.(fl2.gt.0)) then  ! checked

         p = pin
         p(:,6) = pin(:,3)
         p(:,3) = pin(:,4)
         p(:,4) = pin(:,6)

!         call jetktalg(pin,6,pjetout,Njet,weight)
      !  HERE WAS A BUG: NPart
         call Kinematics_TTBARJET(1,(/-pin(1:4,3),-pin(1:4,4),pin(1:4,5),pin(1:4,6),pin(1:4,1),pin(1:4,2)/),MomDK(1:4,1:6),Not_Passed_Cuts,NBin(1:NumHistograms))


         if(Not_Passed_Cuts.eq..false.) then
            call ttqqqq(p,res)
            res = Wgt*PSWgt1*PSWgt2*res    ! re-weight
            do NHisto=1,NumHistograms
                  call intoHisto(NHisto,NBin(NHisto),res)
            enddo
         else
            res = zero
         endif

         if (subtr) then
         call qqb_ttqqqq_dip_12(p,yRnDK1,Wgt,res_subtr)
          res = res + res_subtr
          dip = dip + res_subtr
         endif


!-----------------------------------------------------------------
         elseif((fl1.lt.0).and.(fl2.lt.0)) then !checked

          p = pin
          p(:,4) = pin(:,3)
          p(:,6) = pin(:,4)
          p(:,3) = pin(:,6)

!         call jetktalg(pin,6,pjetout,Njet,weight)
!         if( weight.eq.1d0) then

      !  HERE WAS A BUG: NPart
   call Kinematics_TTBARJET(1,(/-pin(1:4,3),-pin(1:4,4),pin(1:4,5),pin(1:4,6),pin(1:4,1),pin(1:4,2)/),MomDK(1:4,1:6),Not_Passed_Cuts,NBin(1:NumHistograms))


       if(Not_Passed_Cuts.eq..false.) then
           call ttqqqq(p,res)
            res = Wgt*PSWgt1*PSWgt2*res    ! re-weight
            do NHisto=1,NumHistograms
                  call intoHisto(NHisto,NBin(NHisto),res)
            enddo
         else
            res = zero
         endif

         if (subtr) then
         call qbb_ttqqqq_dip(p,yRnDK1,Wgt,res_subtr)
          res = res + res_subtr
          dip = dip + res_subtr
         endif

         endif

      endif   ! for non-id quarks/antiquarks/etc.





     if (abs(fl1).eq.abs(fl2)) then ! same or anti-same flavor initial state

!-----------------------------------------------------------------
         if ((fl1.gt.0).and.(fl2.gt.0)) then !checked

         p = pin
         p(:,5) = pin(:,4)
         p(:,4) = pin(:,5)

!         call jetktalg(pin,6,pjetout,Njet,weight)
!         if( weight.eq.1d0) then
      !  HERE WAS A BUG: NPart
   call Kinematics_TTBARJET(1,(/-pin(1:4,3),-pin(1:4,4),pin(1:4,5),pin(1:4,6),pin(1:4,1),pin(1:4,2)/),MomDK(1:4,1:6),Not_Passed_Cuts,NBin(1:NumHistograms))

       if(Not_Passed_Cuts.eq..false.) then
            call ttqqqq_id(p,res)
            res = Wgt*PSWgt1*PSWgt2*res    ! re-weight
            do NHisto=1,NumHistograms
                  call intoHisto(NHisto,NBin(NHisto),res)
            enddo
         else
            res = zero
         endif


         if (subtr) then

         call qqq_ttqqqq_dip(p,yRnDK1,Wgt,res_subtr)
         res = res + res_subtr
         dip = dip + res_subtr


         p(:,4) = pin(:,6)
         p(:,6) = pin(:,5)


         call qqq_ttqqqq_dip(p,yRnDK1,Wgt,res_subtr)
         res = res + res_subtr
         dip = dip + res_subtr

         endif !

!-----------------------------------------------------------------
         elseif((fl1.gt.0).and.(fl2.lt.0)) then  ! checked

         p = pin

!         call jetktalg(pin,6,pjetout,Njet,weight)
!         if( weight.eq.1d0) then
      !  HERE WAS A BUG: NPart
       call Kinematics_TTBARJET(1,(/-pin(1:4,3),-pin(1:4,4),pin(1:4,5),pin(1:4,6),pin(1:4,1),pin(1:4,2)/),MomDK(1:4,1:6),Not_Passed_Cuts,NBin(1:NumHistograms))

       if(Not_Passed_Cuts.eq..false.) then
            call ttqqqq_id(p,res1)
            call ttqqqq(p,res2)
            res = res1 + dble(nf-1)*res2
            res = Wgt*PSWgt1*PSWgt2*res    ! re-weight
            do NHisto=1,NumHistograms
                  call intoHisto(NHisto,NBin(NHisto),res)
            enddo
       else
            res = zero
       endif

         if (subtr) then
                call qqb_ttqqqq_dip_11(p,yRnDK1,dble(nf)*Wgt,res_subtr)

                  res = res + res_subtr
                  dip = dip + res_subtr

                  p(:,6) = pin(:,4)
                  p(:,4) = pin(:,6)

                call qqb_ttqqqq_dip_12(p,yRnDK1,Wgt,res_subtr)

                res = res + res_subtr
                dip = dip + res_subtr
         endif

!-----------------------------------------------------------------
         elseif((fl1.lt.0).and.(fl2.gt.0)) then ! checked

         p = pin
         p(:,4) = pin(:,3)
         p(:,3) = pin(:,4)


!         call jetktalg(pin,6,pjetout,Njet,weight)
!         if( weight.eq.1d0) then
      !  HERE WAS A BUG: NPart
       call Kinematics_TTBARJET(1,(/-pin(1:4,3),-pin(1:4,4),pin(1:4,5),pin(1:4,6),pin(1:4,1),pin(1:4,2)/),MomDK(1:4,1:6),Not_Passed_Cuts,NBin(1:NumHistograms))


       if(Not_Passed_Cuts.eq..false.) then
            call ttqqqq_id(p,res1)
            call ttqqqq(p,res2)
            res = res1 + dble(nf-1)*res2
            res = Wgt*PSWgt1*PSWgt2*res    ! re-weight
            do NHisto=1,NumHistograms
                  call intoHisto(NHisto,NBin(NHisto),res)
            enddo
         else
            res = zero
         endif
!	print *, "mat.el.",res


         if (subtr) then

         call qqb_ttqqqq_dip_11(p,yRnDK1,dble(nf)*Wgt,res_subtr)

         res = res + res_subtr
         dip = dip + res_subtr

          p(:,6) = pin(:,3)
          p(:,3) = pin(:,4)
          p(:,4) = pin(:,6)

         call qqb_ttqqqq_dip_12(p,yRnDK1,Wgt,res_subtr)

         res = res + res_subtr
         dip = dip + res_subtr

         endif

!-----------------------------------------------------------------
         elseif((fl1.lt.0).and.(fl2.lt.0)) then  ! checked

          p = pin
          p(:,4) = pin(:,3)
          p(:,6) = pin(:,4)
          p(:,3) = pin(:,6)

!         call jetktalg(pin,6,pjetout,Njet,weight)
!         if( weight.eq.1d0) then
      !  HERE WAS A BUG: NPart
   call Kinematics_TTBARJET(1,(/-pin(1:4,3),-pin(1:4,4),pin(1:4,5),pin(1:4,6),pin(1:4,1),pin(1:4,2)/),MomDK(1:4,1:6),Not_Passed_Cuts,NBin(1:NumHistograms))

       if(Not_Passed_Cuts.eq..false.) then
            call ttqqqq_id(p,res)
            res = Wgt*PSWgt1*PSWgt2*res    ! re-weight
            do NHisto=1,NumHistograms
                  call intoHisto(NHisto,NBin(NHisto),res)
            enddo
         else
            res = zero
         endif

         if (subtr) then
         call qbb_ttqqqq_dip(p,yRnDK1,Wgt,res_subtr)
         res = res + res_subtr
         dip = dip + res_subtr

         p(:,3) = pin(:,5)
         p(:,5) = pin(:,6)

         call qbb_ttqqqq_dip(p,yRnDK1,Wgt,res_subtr)
         res = res + res_subtr
         dip = dip + res_subtr

         endif

         endif

     endif
      end subroutine sixquark



     subroutine ttqqqq(p,res)
     use ModProcess
     use ModMisc
     implicit none
     real(8),  intent(in) :: p(4,6)
     real(8), intent(out) :: res
     real(8), save  :: C(5,5) ! color matrix
     logical, save :: first_time = .true.
     integer :: hel(6),i1,i2,i3,i4,i5,i6,i,j, iTree
     type(Particle),save :: ExtParticles(1:6)
     type(TreeProcess),save :: TreeAmpsReal(1:5)  ! there are 4 independent ampls
     complex(8) :: A(5), Ac(5), cres
     integer(8) :: Nmax
     real(dp) :: MomDK(1:4,1:6),PSWgt1,PSWgt2

      res = 0d0
      cres = (0d0,0d0)

      if (first_time) then

       C = 0.0d0
       C(1,1) = 27.0d0
       C(1,2) = 3.0d0
       C(1,3) =  - 17.0d0/3.0d0
       C(1,4) =  - 17.0d0/3.0d0
       C(1,5) =  - 17.0d0/3.0d0
       C(2,1) = 3.0d0
       C(2,2) = 27.0d0
       C(2,3) =  - 17.0d0/3.0d0
       C(2,4) =  - 17.0d0/3.0d0
       C(2,5) =  - 17.0d0/3.0d0
       C(3,1) =  - 17.0d0/3.0d0
       C(3,2) =  - 17.0d0/3.0d0
       C(3,3) = 17.0d0/3.0d0
       C(3,4) = 3.0d0
       C(3,5) = 3.0d0
       C(4,1) =  - 17.0d0/3.0d0
       C(4,2) =  - 17.0d0/3.0d0
       C(4,3) = 3.0d0
       C(4,4) = 17.0d0/3.0d0
       C(4,5) = 3.0d0
       C(5,1) =  - 17.0d0/3.0d0
       C(5,2) =  - 17.0d0/3.0d0
       C(5,3) = 3.0d0
       C(5,4) = 3.0d0
       C(5,5) = 17.0d0/3.0d0


      ! init external particles
      call InitProcess_TbTqbqqbq(ExtParticles(1:6))
      ! init tree processes for 0-> tb t qb q Qb Q
      call InitTrees2(6,0,5,TreeAmpsReal)

      !--- ordering of particles in a particular primitive amplitude;
      !--- has to be correlated with the form program that calculates
      !--- color correlation matrix
        TreeAmpsReal( 1)%PartRef(1:6) =  (/1, 2, 3, 4, 5, 6/)
        TreeAmpsReal( 2)%PartRef(1:6) =  (/1, 2, 5, 6, 3, 4/)
        TreeAmpsReal( 3)%PartRef(1:6) =  (/1, 6, 3, 4, 5, 2/)
        TreeAmpsReal( 4)%PartRef(1:6) =  (/1, 6, 5, 2, 3, 4/)
        TreeAmpsReal( 5)%PartRef(1:6) =  (/1, 4, 5, 6, 3, 2/)

        do iTree=1,5
            call LinkTreeParticles2(TreeAmpsReal(iTree),ExtParticles(1:6))
        enddo

        first_time = .false.
      endif


       if (TopDecays.ge.1) then
         call EvalPhasespace_TopDecay(p(1:4,1),yRnDk(1:4),.false.,MomDK(1:4,1:3),PSWgt1)
         call EvalPhasespace_TopDecay(p(1:4,2),yRnDk(5:8),.false.,MomDK(1:4,4:6),PSWgt2)
       elseif(TopDecays.eq.0) then
         PSWgt1 = one
         PSWgt2 = one
       else
          print *, 'topdecays are not defined'
          pause
       endif

       if (TopDecays.ge.1) then
          Nmax = -1
        elseif(TopDecays.eq.0) then
           Nmax = 1
       endif


         do i1 =-1,Nmax,2
             hel(1) = i1
        do i2 = -1,Nmax,2
            hel(2) = i2
        do i3 = -1,1,2
            hel(3) = i3
        do i4 = -1,1,2
            hel(4) = i4
        do i5 = -1,1,2
            hel(5) = i5
        do i6 = -1,1,2
            hel(6) = i6


      call GenerateEventttqqqq(p,MomDK(1:4,1:6),hel,ExtParticles(1:6))

             do i = 1,5
                call EvalTree2(TreeAmpsReal(i),A(i))
                Ac(i) = conjg(A(i))
             enddo

       do i=1,5
          do j=1,5
           cres = cres + A(i)*Ac(j)*C(i,j)
       enddo
       enddo


      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      res = real(cres)

      end subroutine ttqqqq




     subroutine ttqqqq_id(p,res)
     use ModProcess
     implicit none
     real(8),  intent(in) :: p(4,6)
     real(8), intent(out) :: res
     real(8), save  :: C(10,10) ! color matrix
     logical, save :: first_time = .true.
     integer :: hel(6),i1,i2,i3,i4,i5,i6,i,j, iTree
     type(Particle), save :: ExtParticles(1:6), ExtParticles1(1:6)
     type(TreeProcess), save  :: TreeAmpsReal(1:5)  ! there are 5 ind ampls
     type(TreeProcess), save :: TreeAmpsReal1(1:5) ! there are 5 other ind ampls
     complex(8) :: A(10), Ac(10), cres
     integer(8) :: Nmax
     real(dp) :: MomDK(1:4,1:6),PSWgt1,PSWgt2
!------------------------------------------------------------------

      res = 0d0
      cres = (0d0,0d0)

      if (first_time) then

       C = 0.0d0

      C(1,1)= 27.0d0
      C(1,2)= 3.0d0
      C(1,3)=  - 17.0d0/3.0d0
      C(1,4)=  - 17.0d0/3.0d0
      C(1,5)=  - 17.0d0/3.0d0
      C(1,6)=  - 9.0d0
      C(1,7)=  - 9.0d0
      C(1,8)= 1.0d0
      C(1,9)= 9.0d0
      C(1,10)= 9.0d0
      C(2,1)= 3.0d0
      C(2,2)= 27.0d0
      C(2,3)=  - 17.0d0/3.0d0
      C(2,4)=  - 17.0d0/3.0d0
      C(2,5)=  - 17.0d0/3.0d0
      C(2,6)=  - 9.0d0
      C(2,7)=  - 9.0d0
      C(2,8)= 9.0d0
      C(2,9)= 9.0d0
      C(2,10)= 1.0d0
      C(3,1)=  - 17.0d0/3.0d0
      C(3,2)=  - 17.0d0/3.0d0
      C(3,3)= 17.0d0/3.0d0
      C(3,4)= 3.0d0
      C(3,5)= 3.0d0
      C(3,6)= 1.0d0
      C(3,7)= 9.0d0
      C(3,8)=  - 17.0d0/9.0d0
      C(3,9)=  - 25.0d0/9.0d0
      C(3,10)=  - 17.0d0/9.0d0
      C(4,1)=  - 17.0d0/3.0d0
      C(4,2)=  - 17.0d0/3.0d0
      C(4,3)= 3.0d0
      C(4,4)= 17.0d0/3.0d0
      C(4,5)= 3.0d0
      C(4,6)= 9.0d0
      C(4,7)= 9.0d0
      C(4,8)=  - 25.0d0/9.0d0
      C(4,9)=  - 11.0d0/3.0d0
      C(4,10)=  - 25.0d0/9.0d0
      C(5,1)=  - 17.0d0/3.0d0
      C(5,2)=  - 17.0d0/3.0d0
      C(5,3)= 3.0d0
      C(5,4)= 3.0d0
      C(5,5)= 17.0d0/3.0d0
      C(5,6)= 9.0d0
      C(5,7)= 1.0d0
      C(5,8)=  - 17.0d0/9.0d0
      C(5,9)=  - 25.0d0/9.0d0
      C(5,10)=  - 17.0d0/9.0d0
      C(6,1)=  - 9.0d0
      C(6,2)=  - 9.0d0
      C(6,3)= 1.0d0
      C(6,4)= 9.0d0
      C(6,5)= 9.0d0
      C(6,6)= 27.0d0
      C(6,7)= 3.0d0
      C(6,8)=  - 17.0d0/3.0d0
      C(6,9)=  - 17.0d0/3.0d0
      C(6,10)=  - 17.0d0/3.0d0
      C(7,1)=  - 9.0d0
      C(7,2)=  - 9.0d0
      C(7,3)= 9.0d0
      C(7,4)= 9.0d0
      C(7,5)= 1.0d0
      C(7,6)= 3.0d0
      C(7,7)= 27.0d0
      C(7,8)=  - 17.0d0/3.0d0
      C(7,9)=  - 17.0d0/3.0d0
      C(7,10)=  - 17.0d0/3.0d0
      C(8,1)= 1.0d0
      C(8,2)= 9.0d0
      C(8,3)=  - 17.0d0/9.0d0
      C(8,4)=  - 25.0d0/9.0d0
      C(8,5)=  - 17.0d0/9.0d0
      C(8,6)=  - 17.0d0/3.0d0
      C(8,7)=  - 17.0d0/3.0d0
      C(8,8)= 17.0d0/3.0d0
      C(8,9)= 3.0d0
      C(8,10)= 3.0d0
      C(9,1)= 9.0d0
      C(9,2)= 9.0d0
      C(9,3)=  - 25.0d0/9.0d0
      C(9,4)=  - 11.0d0/3.0d0
      C(9,5)=  - 25.0d0/9.0d0
      C(9,6)=  - 17.0d0/3.0d0
      C(9,7)=  - 17.0d0/3.0d0
      C(9,8)= 3.0d0
      C(9,9)= 17.0d0/3.0d0
      C(9,10)= 3.0d0
      C(10,1)= 9.0d0
      C(10,2)= 1.0d0
      C(10,3)=  - 17.0d0/9.0d0
      C(10,4)=  - 25.0d0/9.0d0
      C(10,5)=  - 17.0d0/9.0d0
      C(10,6)=  - 17.0d0/3.0d0
      C(10,7)=  - 17.0d0/3.0d0
      C(10,8)= 3.0d0
      C(10,9)= 3.0d0
      C(10,10)= 17.0d0/3.0d0



      ! init external particles
      call InitProcess_TbTqbqqbq(ExtParticles(1:6))
      ! init tree processes for 0-> tb t qb q Qb Q
      call InitTrees2(6,0,5,TreeAmpsReal)


      call InitProcess_TbTqbqqbq1(ExtParticles1(1:6))
      ! init tree processes for 0-> tb t qb q Qb Q
      call InitTrees2(6,0,5,TreeAmpsReal1)


      !--- ordering of particles in a particular primitive amplitude;
      !--- has to be correlated with the form program that calculates
      !--- color correlation matrix


        TreeAmpsReal(1)%PartRef(1:6) =  (/1, 2, 3, 4, 5, 6/)
        TreeAmpsReal(2)%PartRef(1:6) =  (/1, 2, 5, 6, 3, 4/)
        TreeAmpsReal(3)%PartRef(1:6) =  (/1, 6, 3, 4, 5, 2/)
        TreeAmpsReal(4)%PartRef(1:6) =  (/1, 6, 5, 2, 3, 4/)
        TreeAmpsReal(5)%PartRef(1:6) =  (/1, 4, 5, 6, 3, 2/)

        TreeAmpsReal1(1)%PartRef(1:6) =  (/1, 2, 3, 6, 5, 4/)
        TreeAmpsReal1(2)%PartRef(1:6) =  (/1, 2, 5, 4, 3, 6/)
        TreeAmpsReal1(3)%PartRef(1:6) =  (/1, 4, 3, 6, 5, 2/)
        TreeAmpsReal1(4)%PartRef(1:6) =  (/1, 4, 5, 2, 3, 6/)
        TreeAmpsReal1(5)%PartRef(1:6) =  (/1, 6, 5, 4, 3, 2/)

        do iTree=1,5
            call LinkTreeParticles2(TreeAmpsReal(iTree),ExtParticles(1:6))
            call LinkTreeParticles2(TreeAmpsReal1(iTree),ExtParticles1(1:6))
        enddo

      first_time = .false.
      endif


        if (TopDecays.ge.1) then
          call EvalPhasespace_TopDecay(p(1:4,1),yRnDk(1:4),.false., MomDK(1:4,1:3),PSWgt1)
          call EvalPhasespace_TopDecay(p(1:4,2),yRnDk(5:8),.false.,MomDK(1:4,4:6),PSWgt2)
        elseif(TopDecays.eq.0) then
          PSWgt1 = one
          PSWgt2 = one
        else
           print *, 'topdecays are not defined'
           pause
        endif

       if (TopDecays.ge.1) then
          Nmax = -1
        elseif(TopDecays.eq.0) then
           Nmax = 1
       endif


         do i1 =-1,Nmax,2
             hel(1) = i1
        do i2 = -1,Nmax,2
            hel(2) = i2
        do i3 = -1,1,2
            hel(3) = i3
        do i4 = -1,1,2
            hel(4) = i4
        do i5 = -1,1,2
            hel(5) = i5
        do i6 = -1,1,2
            hel(6) = i6


    call GenerateEventttqqqq(p,MomDK(1:4,1:6),hel,ExtParticles(1:6))
    call GenerateEventttqqqq(p,MomDK(1:4,1:6),hel,ExtParticles1(1:6))

             do i = 1,5
      call EvalTree2(TreeAmpsReal(i),A(i))
        Ac(i) = conjg(A(i))
             enddo

             do i=6,10
      call EvalTree2(TreeAmpsReal1(i-5),A(i))
        Ac(i) = conjg(A(i))
             enddo



       do i=1,10
          do j=1,10
           cres = cres + A(i)*Ac(j)*C(i,j)
       enddo
       enddo


      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      res = real(cres)

      end subroutine ttqqqq_id





SUBROUTINE InitProcess_TbTqbqqbq(ExtParticles) ! 80.25.2009
use ModProcess
use ModParameters
use ModMisc
implicit none
type(Particle) :: ExtParticles(:)

 ExtParticles(1)%PartType = -5
 ExtParticles(1)%ExtRef   = 1
 ExtParticles(1)%Mass = m_Top
 ExtParticles(1)%Mass2= ExtParticles(1)%Mass**2

 ExtParticles(2)%PartType = 5
 ExtParticles(2)%ExtRef   = 2
 ExtParticles(2)%Mass = m_Top
 ExtParticles(2)%Mass2= ExtParticles(2)%Mass**2

 ExtParticles(3)%PartType = -4
 ExtParticles(3)%ExtRef   = 3
 ExtParticles(3)%Mass = 0d0
 ExtParticles(3)%Mass2= 0d0

 ExtParticles(4)%PartType = 4
 ExtParticles(4)%ExtRef   = 4
 ExtParticles(4)%Mass = 0d0
 ExtParticles(4)%Mass2= 0d0

 ExtParticles(5)%PartType = -3
 ExtParticles(5)%ExtRef   = 5
 ExtParticles(5)%Mass = 0d0
 ExtParticles(5)%Mass2= 0d0

 ExtParticles(6)%PartType = 3
 ExtParticles(6)%ExtRef   = 6
 ExtParticles(6)%Mass = 0d0
 ExtParticles(6)%Mass2= 0d0

RETURN
END SUBROUTINE InitProcess_TbTqbqqbq




SUBROUTINE InitProcess_TbTqbqqbq1(ExtParticles) ! 80.25.2009
use ModProcess
use ModParameters
use ModMisc
implicit none
type(Particle) :: ExtParticles(:)

 ExtParticles(1)%PartType = -5
 ExtParticles(1)%ExtRef   = 1
 ExtParticles(1)%Mass = m_Top
 ExtParticles(1)%Mass2= ExtParticles(1)%Mass**2

 ExtParticles(2)%PartType = 5
 ExtParticles(2)%ExtRef   = 2
 ExtParticles(2)%Mass = m_Top
 ExtParticles(2)%Mass2= ExtParticles(2)%Mass**2

 ExtParticles(3)%PartType = -4
 ExtParticles(3)%ExtRef   = 3
 ExtParticles(3)%Mass = 0d0
 ExtParticles(3)%Mass2= 0d0

 ExtParticles(4)%PartType = 3
 ExtParticles(4)%ExtRef   = 4
 ExtParticles(4)%Mass = 0d0
 ExtParticles(4)%Mass2= 0d0

 ExtParticles(5)%PartType = -3
 ExtParticles(5)%ExtRef   = 5
 ExtParticles(5)%Mass = 0d0
 ExtParticles(5)%Mass2= 0d0

 ExtParticles(6)%PartType = 4
 ExtParticles(6)%ExtRef   = 6
 ExtParticles(6)%Mass = 0d0
 ExtParticles(6)%Mass2= 0d0

RETURN
END SUBROUTINE InitProcess_TbTqbqqbq1

SUBROUTINE InitTrees2(NumQuarks,NumGluons,NumTrees,TreeAmpsDip) ! NEW
use ModMisc
use ModProcess
use ModParameters
implicit none
integer :: NumQuarks,NumGluons,NumTrees
type(TreeProcess) :: TreeAmpsDip(:)
integer :: NumParticles,iTree

 NumParticles=NumQuarks+NumGluons
 do iTree=1,NumTrees
     TreeAmpsDip(iTree)%NumPart=NumParticles
     TreeAmpsDip(iTree)%NumQua=NumQuarks
     allocate( TreeAmpsDip(iTree)%NumGlu(0:NumQuarks)     )
     allocate( TreeAmpsDip(iTree)%PartRef(1:NumParticles) )
     allocate( TreeAmpsDip(iTree)%PartType(1:NumParticles))
     allocate( TreeAmpsDip(iTree)%Quarks(1:NumQuarks)     )
     allocate( TreeAmpsDip(iTree)%Gluons(1:NumGluons)     )
     TreeAmpsDip(iTree)%NumGlu(0) = NumGluons
 enddo

RETURN
END SUBROUTINE InitTrees2


SUBROUTINE LinkTreeParticles2(TheTreeAmp,TheParticles) ! NEW
use ModProcess
use ModParameters
implicit none
type(TreeProcess) :: TheTreeAmp
type(Particle),target :: TheParticles(:)
integer :: iPart,PartRef,PartType,ig,iq,NPart,counter,QuarkPos(1:6)


           counter = 0
           do NPart=1,TheTreeAmp%NumPart
                 TheTreeAmp%PartType(NPart) = TheParticles( TheTreeAmp%PartRef(NPart) )%PartType
                 if( TheTreeAmp%PartType(NPart).ne.10 ) then  ! not a gluon
!                      TheTreeAmp%NumQua = TheTreeAmp%NumQua + 1
                    counter = counter + 1
                    QuarkPos(counter) = NPart
                 endif
           enddo
!           set number of gluons between quark lines
           if( TheTreeAmp%PartType(1).ne.10 ) then ! not a gluon
             if( TheTreeAmp%NumQua .eq. 2 ) then
                   TheTreeAmp%NumGlu(1) = QuarkPos(2) - QuarkPos(1) - 1
                   TheTreeAmp%NumGlu(2) = TheTreeAmp%NumPart - QuarkPos(2)
             endif
             if( TheTreeAmp%NumQua .eq. 4 ) then
                   TheTreeAmp%NumGlu(1) = QuarkPos(2) - QuarkPos(1) - 1
                   TheTreeAmp%NumGlu(2) = QuarkPos(3) - QuarkPos(2) - 1
                   TheTreeAmp%NumGlu(3) = QuarkPos(4) - QuarkPos(3) - 1
                   TheTreeAmp%NumGlu(4) = TheTreeAmp%NumPart - QuarkPos(4)
             endif
             if( TheTreeAmp%NumQua .eq. 6 ) then
                   TheTreeAmp%NumGlu(1) = QuarkPos(2) - QuarkPos(1) - 1
                   TheTreeAmp%NumGlu(2) = QuarkPos(3) - QuarkPos(2) - 1
                   TheTreeAmp%NumGlu(3) = QuarkPos(4) - QuarkPos(3) - 1
                   TheTreeAmp%NumGlu(4) = QuarkPos(5) - QuarkPos(4) - 1
                   TheTreeAmp%NumGlu(5) = QuarkPos(6) - QuarkPos(5) - 1
                   TheTreeAmp%NumGlu(6) = TheTreeAmp%NumPart - QuarkPos(6)
             endif


           elseif( TheTreeAmp%PartType(1).eq.0 ) then
             if( TheTreeAmp%NumQua .eq. 2 ) then
                   TheTreeAmp%NumGlu(1) = QuarkPos(1) - 2
                   TheTreeAmp%NumGlu(2) = QuarkPos(2) - QuarkPos(1) - 1
                   TheTreeAmp%NumGlu(3) = TheTreeAmp%NumPart - QuarkPos(2)
             endif
             if( TheTreeAmp%NumQua .eq. 4 ) then
                   TheTreeAmp%NumGlu(1) = QuarkPos(1) - 2
                   TheTreeAmp%NumGlu(2) = QuarkPos(2) - QuarkPos(1) - 1
                   TheTreeAmp%NumGlu(3) = QuarkPos(3) - QuarkPos(2) - 1
                   TheTreeAmp%NumGlu(4) = QuarkPos(4) - QuarkPos(3) - 1
                   TheTreeAmp%NumGlu(5) = TheTreeAmp%NumPart - QuarkPos(4)
             endif
             if( TheTreeAmp%NumQua .eq. 6 ) then
                   TheTreeAmp%NumGlu(1) = QuarkPos(1) - 2
                   TheTreeAmp%NumGlu(2) = QuarkPos(2) - QuarkPos(1) - 1
                   TheTreeAmp%NumGlu(3) = QuarkPos(3) - QuarkPos(2) - 1
                   TheTreeAmp%NumGlu(4) = QuarkPos(4) - QuarkPos(3) - 1
                   TheTreeAmp%NumGlu(5) = QuarkPos(5) - QuarkPos(4) - 1
                   TheTreeAmp%NumGlu(6) = QuarkPos(6) - QuarkPos(5) - 1
                   TheTreeAmp%NumGlu(7) = TheTreeAmp%NumPart - QuarkPos(6)
             endif
           endif


  ig=0; iq=0;
  do iPart=1,TheTreeAmp%NumPart
     PartRef = TheTreeAmp%PartRef(iPart)
     PartType= TheParticles(PartRef)%PartType
     if( PartType.eq.10 ) then  ! PartType==Gluon
           ig=ig+1
           TheTreeAmp%Gluons(ig)%PartType => TheParticles(PartRef)%PartType
           TheTreeAmp%Gluons(ig)%ExtRef   => TheParticles(PartRef)%ExtRef
           TheTreeAmp%Gluons(ig)%Mass     => TheParticles(PartRef)%Mass
           TheTreeAmp%Gluons(ig)%Mass2    => TheParticles(PartRef)%Mass2
           TheTreeAmp%Gluons(ig)%Helicity => TheParticles(PartRef)%Helicity
           TheTreeAmp%Gluons(ig)%Mom      => TheParticles(PartRef)%Mom
           TheTreeAmp%Gluons(ig)%Pol      => TheParticles(PartRef)%Pol
           if( PartType.ne.TheParticles(PartRef)%PartType ) print *,"Error1 in LinkTreeParticles"
     else ! PartType==Quark
           iq=iq+1
           TheTreeAmp%Quarks(iq)%PartType => TheParticles(PartRef)%PartType
           TheTreeAmp%Quarks(iq)%ExtRef   => TheParticles(PartRef)%ExtRef
           TheTreeAmp%Quarks(iq)%Mass     => TheParticles(PartRef)%Mass
           TheTreeAmp%Quarks(iq)%Mass2    => TheParticles(PartRef)%Mass2
           TheTreeAmp%Quarks(iq)%Helicity => TheParticles(PartRef)%Helicity
           TheTreeAmp%Quarks(iq)%Mom      => TheParticles(PartRef)%Mom
           TheTreeAmp%Quarks(iq)%Pol      => TheParticles(PartRef)%Pol
           if( PartType.ne.TheParticles(PartRef)%PartType ) print *,"Error2 in LinkTreeParticles"
     endif
  enddo
  if( ig.ne.TheTreeAmp%NumGlu(0) .OR. iq.ne.TheTreeAmp%NumQua ) print *,"Error3 in LinkTreeParticles"

return
END SUBROUTINE


SUBROUTINE GenerateEventttqqqq(Mom,MomDK,hel,ExtParticles)
use ModMisc
use ModProcess
use ModParameters
implicit none
type(Particle) :: ExtParticles(:)
real(8) , intent(in) :: Mom(1:4,1:6),MomDK(1:4,1:6)
integer, intent(in) :: hel(1:6)
integer :: i1,i2,i3,i4,i5,i6

     i1 = hel(1)
     i2 = hel(2)
     i3 = hel(3)
     i4 = hel(4)
     i5 = hel(5)
     i6 = hel(6)


     ExtParticles(1)%Mom(1:4) = dcmplx(Mom(1:4,1))   ! HERE WAS A BUG: these two assignments were inside the TopDecay if-condition
     ExtParticles(2)%Mom(1:4) = dcmplx(Mom(1:4,2))
     if (TopDecays.ge.1) then
        call TopDecay(ExtParticles(1),DK_LO,MomDK(1:4,1:3))
        call TopDecay(ExtParticles(2),DK_LO,MomDK(1:4,4:6))
     else
        call vSpi(ExtParticles(1)%Mom(1:4),ExtParticles(1)%Mass,i1,ExtParticles(1)%Pol(1:4))
        call ubarSpi(ExtParticles(2)%Mom(1:4),ExtParticles(2)%Mass,i2,ExtParticles(2)%Pol(1:4))
    endif

    ExtParticles(3)%Mom(1:4) =dcmplx(Mom(1:4,3))
    call vSpi(ExtParticles(3)%Mom(1:4),ExtParticles(3)%Mass,i3,ExtParticles(3)%Pol(1:4))

    ExtParticles(4)%Mom(1:4) =dcmplx(Mom(1:4,4))
    call ubarSpi(ExtParticles(4)%Mom(1:4),ExtParticles(4)%Mass,i4,ExtParticles(4)%Pol(1:4))

    ExtParticles(5)%Mom(1:4) = dcmplx(Mom(1:4,5))
    call vSpi(ExtParticles(5)%Mom(1:4),ExtParticles(5)%Mass,i5,ExtParticles(5)%Pol(1:4))

    ExtParticles(6)%Mom(1:4) = dcmplx(Mom(1:4,6))
    call ubarSpi(ExtParticles(6)%Mom(1:4),ExtParticles(6)%Mass,i6,ExtParticles(6)%Pol(1:4))
RETURN
END SUBROUTINE GenerateEventttqqqq











!   subroutine sixquark_dip_int(p,z,in1,in2,res)
!   real(dp), intent(in) :: p(4,5)
!   integer, intent(in) :: in1, in2
!   real(dp), intent(out) :: res(1:3)
!   complex(dp) :: cres
!   type(Particle) :: ExtParticles1(1:5),ExtParticles2(1:5),ExtParticles3(1:5)
!   type(TreeProcess) :: TreeAmpsDip(1:4),TreeAmpsDipg(1:4),TreeAmpsDipq(1:6)
!   complex(dp) :: Bm(1:4),Bq(1:6)
!   complex(dp) :: AM(1:4,1:4)
!   integer :: i1,i2,i3,i4,hel(5),i,j,i5,i6,emi,n,iTree
!   real(dp) :: dipsoft_ii, dipfini_ii, dipplus_ii
!   real(dp) :: dipsoft_fi, dipfini_fi, dipplus_fi
!   real(dp) :: z, C(4,4)
!   real(dp) :: Cq(6,6)
!   real(dp) :: mtrs(2)
!   real(dp) :: q(4,5)
!   real(dp) :: p1(4,5)
!   character :: diptype*2
!
!
!        res = zero
!
!       dipsoft_ii =ii_qg(zero,zero,p,3,4,z,1)
!       dipfini_ii =ii_qg(zero,zero,p,3,4,z,2)
!       dipplus_ii =ii_qg(zero,zero,p,3,4,z,3)
!
!
!   if (abs(in1).ne.abs(in2)) then
!        if ((in1.gt.0).and.(in2.gt.0)) then
!           call qqq_ttqqqq_dip_int(p,mtrs)
!         !  res = res + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2))*fx1(in1)*fx2(in2) &
!         !  +mtrs(1)*(dipfini_ii+dipplus_ii)*fx1z(in1)*fx2(in2)/z  &
!         !  +mtrs(2)*(dipfini_ii+dipplus_ii)*fx1(in1)*fx2z(in2)/z
!           res(1) = res(1) + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2))
!           res(2) = res(2) + (dipfini_ii+dipplus_ii)*mtrs(1)
!           res(3) = res(3) + (dipfini_ii+dipplus_ii)*mtrs(2)
!
!        elseif ((in1.gt.0).and.(in2.lt.0)) then
!           call qqb_ttqqqq_dip_12_int(p,mtrs)
!         !  res = res + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2))*fx1(in1)*fx2(in2) &
!         !  +mtrs(1)*(dipfini_ii+dipplus_ii)*fx1z(in1)*fx2(in2)/z &
!         !  +mtrs(2)*(dipfini_ii+dipplus_ii)*fx1(in1)*fx2z(in2)/z
!           res(1) = res(1) + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2))
!           res(2) = res(2) + (dipfini_ii+dipplus_ii)*mtrs(1)
!           res(3) = res(3) + (dipfini_ii+dipplus_ii)*mtrs(2)
!
!        elseif ((in1.lt.0).and.(in2.gt.0)) then
!           p1 = p
!           p1(:,3) = p(:,4)
!           p1(:,4) = p(:,3)
!           call qqb_ttqqqq_dip_12_int(p1,mtrs)
!           !  res = res + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2))*fx1(in1)*fx2(in2) &
!           !  +mtrs(1)*(dipfini_ii+dipplus_ii)*fx1z(in1)*fx2(in2)/z &
!           !  +mtrs(2)*(dipfini_ii+dipplus_ii)*fx1(in1)*fx2z(in2)/z
!           res(1) = res(1) + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2))
!           res(2) = res(2) + (dipfini_ii+dipplus_ii)*mtrs(1)
!           res(3) = res(3) + (dipfini_ii+dipplus_ii)*mtrs(2)
!
!        elseif ((in1.lt.0).and.(in2.lt.0)) then
!           call qbb_ttqqqq_dip_int(p,mtrs)
!         !   res = res + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2))*fx1(in1)*fx2(in2) &
!         !     +mtrs(1)*(dipfini_ii+dipplus_ii)*fx1z(in1)*fx2(in2)/z &
!         !     +mtrs(2)*(dipfini_ii+dipplus_ii)*fx1(in1)*fx2z(in2)/z
!           res(1) = res(1) + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2))
!           res(2) = res(2) + (dipfini_ii+dipplus_ii)*mtrs(1)
!           res(3) = res(3) + (dipfini_ii+dipplus_ii)*mtrs(2)
!
!        endif
!   endif  ! for not same/ anti-same flavors
!
!
!
!   if (abs(in1).eq.abs(in2)) then
!        if ((in1.gt.0).and.(in2.gt.0)) then
!           call qqq_ttqqqq_dip_int(p,mtrs)
!         !  res = res + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2))*fx1(in1)*fx2(in2) &
!         !  +mtrs(1)*(dipfini_ii+dipplus_ii)*fx1z(in1)*fx2(in2)/z &
!         !  +mtrs(2)*(dipfini_ii+dipplus_ii)*fx1(in1)*fx2z(in2)/z
!           res(1) = res(1) + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2))
!           res(2) = res(2) + (dipfini_ii+dipplus_ii)*mtrs(1)
!           res(3) = res(3) + (dipfini_ii+dipplus_ii)*mtrs(2)
!           res(1:3) = two*res(1:3)   ! this factor of two is because there
!                                     ! are twice as many dipoles
!                                     ! but we also need the symmetry factor
!
!        elseif ((in1.gt.0).and.(in2.lt.0)) then
!           call qqb_ttqqqq_dip_12_int(p,mtrs)
!           !  res = res + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2))*fx1(in1)*fx2(in2) &
!           !  +mtrs(1)*(dipfini_ii+dipplus_ii)*fx1z(in1)*fx2(in2)/z &
!           !  +mtrs(2)*(dipfini_ii+dipplus_ii)*fx1(in1)*fx2z(in2)/z
!           res(1) = res(1) + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2))
!           res(2) = res(2) + (dipfini_ii+dipplus_ii)*mtrs(1)
!           res(3) = res(3) + (dipfini_ii+dipplus_ii)*mtrs(2)
!
!           call qqb_ttqqqq_dip_11_int(p,mtrs)
!           dipsoft_fi = fi_qg(zero,zero,p,3,5,z,1)
!           dipfini_fi = fi_qg(zero,zero,p,3,5,z,2)
!           dipplus_fi = fi_qg(zero,zero,p,3,5,z,3)
!         !  res = res + (dipsoft_fi-dipplus_fi)*(mtrs(1)+mtrs(2))*fx1(in1)*fx2(in2) &
!         !  +mtrs(1)*(dipfini_fi+dipplus_fi)*fx1z(in1)*fx2(in2)/z &
!         !  +mtrs(2)*(dipfini_fi+dipplus_fi)*fx1(in1)*fx2z(in2)/z
!           res(1) = res(1) + (dipsoft_fi-dipplus_fi)*(mtrs(1)+mtrs(2))
!           res(2) = res(2) + (dipfini_fi+dipplus_fi)*mtrs(1)
!           res(3) = res(3) + (dipfini_fi+dipplus_fi)*mtrs(2)
!
!        elseif ((in1.lt.0).and.(in2.gt.0)) then
!           p1 = p
!           p1(:,3) = p(:,4)
!           p1(:,4) = p(:,3)
!           call qqb_ttqqqq_dip_12_int(p1,mtrs)
!           !  res = res + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2))*fx1(in1)*fx2(in2) &
!           !  +mtrs(1)*(dipfini_ii+dipplus_ii)*fx1z(in1)*fx2(in2)/z &
!           !  +mtrs(2)*(dipfini_ii+dipplus_ii)*fx1(in1)*fx2z(in2)/z
!           res(1) = res(1) + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2))
!           res(2) = res(2) + (dipfini_ii+dipplus_ii)*mtrs(1)
!           res(3) = res(3) + (dipfini_ii+dipplus_ii)*mtrs(2)
!
!           dipsoft_fi = fi_qg(zero,zero,p,4,5,z,1)
!           dipfini_fi = fi_qg(zero,zero,p,4,5,z,2)
!           dipplus_fi = fi_qg(zero,zero,p,4,5,z,3)
!           call qqb_ttqqqq_dip_11_int(p1,mtrs)
!           !  res = res + (dipsoft_fi-dipplus_fi)*(mtrs(1)+mtrs(2))*fx1(in1)*fx2(in2) &
!           !  +mtrs(1)*(dipfini_fi+dipplus_fi)*fx1z(in1)*fx2(in2)/z  &
!           !  +mtrs(2)*(dipfini_fi+dipplus_fi)*fx1(in1)*fx2z(in2)/z
!           res(1) = res(1) + (dipsoft_fi-dipplus_fi)*(mtrs(1)+mtrs(2))
!           res(2) = res(2) + (dipfini_fi+dipplus_fi)*mtrs(1)
!           res(3) = res(3) + (dipfini_fi+dipplus_fi)*mtrs(2)
!
!        elseif ((in1.lt.0).and.(in2.lt.0)) then
!           call qbb_ttqqqq_dip_int(p,mtrs)
!           !   res = res + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2))*fx1(in1)*fx2(in2) &
!           !     +mtrs(1)*(dipfini_ii+dipplus_ii)*fx1z(in1)*fx2(in2)/z  &
!           !     +mtrs(2)*(dipfini_ii+dipplus_ii)*fx1(in1)*fx2z(in2)/z
!           res(1) = res(1) + (dipsoft_ii-dipplus_ii)*(mtrs(1)+mtrs(2))
!           res(2) = res(2) + (dipfini_ii+dipplus_ii)*mtrs(1)
!           res(3) = res(3) + (dipfini_ii+dipplus_ii)*mtrs(2)
!           res(1:3) = two*res(1:3)   ! this factor of two is because there
!                                     ! are twice as many dipoles
!                                     ! but we also need the symmetry factor
!        endif
!   endif  ! for not same/ anti-same flavors
!
! end subroutine sixquark_dip_int
!
!
!
!
!
!
!   subroutine qqb_ttqqqq_dip_12_int(p,res)
!   real(dp), intent(in) :: p(4,5)
!   real(dp), intent(out) :: res(2)
!   complex(dp) :: cres
!   type(Particle) :: ExtParticles1(1:5),ExtParticles2(1:5),ExtParticles3(1:5)
!   type(TreeProcess) :: TreeAmpsDip(1:4),TreeAmpsDipg(1:4),TreeAmpsDipq(1:6)
!   complex(dp) :: Bm(1:4),Bq(1:6)
!   complex(dp) :: AM(1:4,1:4)
!   integer :: i1,i2,i3,i4,hel(5),i,j,i5,i6,emi,n, in1, in2,iTree
!   real(dp) :: dipsoft, dipfini, dipplus
!   real(dp) :: z, C(4,4), fx1(-5:5), fx2(-5:5),fx1z(-5:5),fx2z(-5:5)
!   real(dp) :: Cq(6,6)
!   real(dp) :: mtrs
!   real(dp) :: q(4,5)
!   real(dp) :: pd1(4,5), pd2(4,5)
!   character :: diptype*2
!
!
!
!       res = zero
!       cres = (zero,zero)
!
! !---- there are two dipoles here
!        pd1 = p
!        pd1(:,3) = p(:,5)
!        pd1(:,5) = p(:,3)
!
!        pd2 = p
!        pd2(:,4) = p(:,5)
!        pd2(:,5) = p(:,4)
!
! !      init external particles
!       call InitProcess_TbTQbQG(ExtParticles1(1:5))
! !      init tree processes for 0-> tb t bq q g g
!        call InitTrees(4,1,4,TreeAmpsDip)
!
!       TreeAmpsDip(1)%PartRef(1:5) = (/1,2,3,4,5/)
!       TreeAmpsDip(2)%PartRef(1:5) = (/1,5,2,3,4/)
!       TreeAmpsDip(3)%PartRef(1:5) = (/1,2,5,3,4/)
!       TreeAmpsDip(4)%PartRef(1:5) = (/1,2,3,5,4/)
!
!
!       do iTree=1,4
!       call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles1(1:5))
!       enddo
!
!       do n = 1,2
!
!          if (n.eq.1) then
!             q = pd1
!          elseif(n.eq.2) then
!             q = pd2
!          endif
!
!
!       AM = (zero,zero) ! first, the ``amplitude matrix'' is set to zero
!
!       do i1=-1,1,2
!          do i2=-1,1,2
!             do i3 = -1,1,2
!                do i4 = -1,1,2
!                   do i5 = -1,1,2
!
!                  hel(1) = i1
!                  hel(2) = i2
!                  hel(3) = i3
!                  hel(4) = i4
!                  hel(5) = i5
!
!       call GenerateEventttqqg(q,hel,ExtParticles1(1:5))
!
!       do i6 = 1,4
!       call EvalTree2(TreeAmpsDip(i6),Bm(i6))
!       enddo
!
!        do i=1,4
!           do j=1,4
!        AM(i,j) = AM(i,j) + Bm(i)*conjg(Bm(j)) !color correlation matrix
!          enddo
!        enddo
!
!
!                   enddo
!                 enddo
!               enddo
!             enddo
!           enddo
!
!
!         call cc_ttqqg(C)
!
!         mtrs = zero
!
!       do i=1,4
!       do j=1,4
!     mtrs = mtrs + C(i,j)*real(AM(i,j),dp)
!       enddo
!       enddo
!
!
!       res(n) = mtrs
!
!        enddo
!
!
! end subroutine qqb_ttqqqq_dip_12_int
!
!
!
!   subroutine qqb_ttqqqq_dip_11_int(p,res)
!   real(dp), intent(in) :: p(4,5)
!   real(dp), intent(out) :: res(2)  !but in reality only one dip here
!   complex(dp) :: cres
!   type(Particle) :: ExtParticles1(1:5),ExtParticles2(1:5),ExtParticles3(1:5)
!   type(TreeProcess) :: TreeAmpsDip(1:4),TreeAmpsDipg(1:4),TreeAmpsDipq(1:6)
!   complex(dp) :: Bm(1:4),Bq(1:6)
!   complex(dp) :: AM(1:4,1:4)
!   integer :: i1,i2,i3,i4,hel(5),i,j,i5,i6,emi,n, in1, in2,iTree
!   real(dp) :: dipsoft, dipfini, dipplus
!   real(dp) :: z, C(4,4), fx1(-5:5), fx2(-5:5),fx1z(-5:5),fx2z(-5:5)
!   real(dp) :: Cq(6,6)
!   real(dp) :: mtrs
!   real(dp) :: q(4,5)
!   real(dp) :: pd1(4,5), pd2(4,5)
!   character :: diptype*2
!
!       res = zero
!       cres = (zero,zero)
!
! !---- there is only one dipole here
!
!        pd1 = p
!
! !      init external particles
!       call InitProcess_TbTQbQG(ExtParticles1(1:5))
! !      init tree processes for 0-> tb t bq q g g
!        call InitTrees(4,1,4,TreeAmpsDip)
!
!       TreeAmpsDip(1)%PartRef(1:5) = (/1,2,3,4,5/)
!       TreeAmpsDip(2)%PartRef(1:5) = (/1,5,2,3,4/)
!       TreeAmpsDip(3)%PartRef(1:5) = (/1,2,5,3,4/)
!       TreeAmpsDip(4)%PartRef(1:5) = (/1,2,3,5,4/)
!
!
!       do iTree=1,4
!       call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles1(1:5))
!       enddo
!
!       do n = 1,1
!
!          q = pd1
!
!       AM = (zero,zero) ! first, the ``amplitude matrix'' is set to zero
!
!       do i1=-1,1,2
!          do i2=-1,1,2
!             do i3 = -1,1,2
!                do i4 = -1,1,2
!                   do i5 = -1,1,2
!
!                  hel(1) = i1
!                  hel(2) = i2
!                  hel(3) = i3
!                  hel(4) = i4
!                  hel(5) = i5
!
!       call GenerateEventttqqg(q,hel,ExtParticles1(1:5))
!
!       do i6 = 1,4
!       call EvalTree2(TreeAmpsDip(i6),Bm(i6))
!       enddo
!
!        do i=1,4
!           do j=1,4
!        AM(i,j) = AM(i,j) + Bm(i)*conjg(Bm(j)) !color correlation matrix
!          enddo
!        enddo
!
!
!                   enddo
!                 enddo
!               enddo
!             enddo
!           enddo
!
!
!         call cc_ttqqg(C)
!
!         mtrs = zero
!
!       do i=1,4
!       do j=1,4
!     mtrs = mtrs + C(i,j)*real(AM(i,j),dp)
!       enddo
!       enddo
!
!
!       res(n) = mtrs
!
!        enddo
!
!
! end subroutine qqb_ttqqqq_dip_11_int
!
!
!
!
!   subroutine qqq_ttqqqq_dip_int(p,res)
!   real(dp), intent(in) :: p(4,5)
!   real(dp), intent(out) :: res(2)
!   complex(dp) :: cres
!   type(Particle) :: ExtParticles1(1:5),ExtParticles2(1:5),ExtParticles3(1:5)
!   type(TreeProcess) :: TreeAmpsDip(1:4),TreeAmpsDipg(1:4),TreeAmpsDipq(1:6)
!   complex(dp) :: Bm(1:4),Bq(1:6)
!   complex(dp) :: AM(1:4,1:4)
!   integer :: i1,i2,i3,i4,hel(5),i,j,i5,i6,emi,n, in1, in2,iTree
!   real(dp) :: dipsoft, dipfini, dipplus
!   real(dp) :: z, C(4,4), fx1(-5:5), fx2(-5:5),fx1z(-5:5),fx2z(-5:5)
!   real(dp) :: Cq(6,6)
!   real(dp) :: mtrs
!   real(dp) :: q(4,5)
!   real(dp) :: pd1(4,5), pd2(4,5)
!   character :: diptype*2
!
!
!
!       res = zero
!       cres = (zero,zero)
!
!
! !---- there are two dipoles here
!
!       pd1 = p
!       pd1(:,3) = p(:,4)
!       pd1(:,4) = p(:,5)
!       pd1(:,5) = p(:,3)
!
!
!       pd2 = p
!       pd2(:,4) = p(:,5)
!       pd2(:,5) = p(:,4)
!
!
!
! !      init external particles
!       call InitProcess_TbTQbQG(ExtParticles1(1:5))
! !      init tree processes for 0-> tb t bq q g g
!        call InitTrees(4,1,4,TreeAmpsDip)
!
!       TreeAmpsDip(1)%PartRef(1:5) = (/1,2,3,4,5/)
!       TreeAmpsDip(2)%PartRef(1:5) = (/1,5,2,3,4/)
!       TreeAmpsDip(3)%PartRef(1:5) = (/1,2,5,3,4/)
!       TreeAmpsDip(4)%PartRef(1:5) = (/1,2,3,5,4/)
!
!
!       do iTree=1,4
!       call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles1(1:5))
!       enddo
!
!
!       do n = 1,2
!
!          if (n.eq.1) then
!             q = pd1
!          elseif(n.eq.2) then
!             q = pd2
!          endif
!
!
!       AM = (zero,zero) ! first, the ``amplitude matrix'' is set to zero
!
!       do i1=-1,1,2
!          do i2=-1,1,2
!             do i3 = -1,1,2
!                do i4 = -1,1,2
!                   do i5 = -1,1,2
!
!                  hel(1) = i1
!                  hel(2) = i2
!                  hel(3) = i3
!                  hel(4) = i4
!                  hel(5) = i5
!
!       call GenerateEventttqqg(q,hel,ExtParticles1(1:5))
!
!       do i6 = 1,4
!       call EvalTree2(TreeAmpsDip(i6),Bm(i6))
!       enddo
!
!        do i=1,4
!           do j=1,4
!        AM(i,j) = AM(i,j) + Bm(i)*conjg(Bm(j)) !color correlation matrix
!          enddo
!        enddo
!
!
!                   enddo
!                 enddo
!               enddo
!             enddo
!           enddo
!
!
!         call cc_ttqqg(C)
!
!         mtrs = zero
!
!       do i=1,4
!       do j=1,4
!     mtrs = mtrs + C(i,j)*real(AM(i,j),dp)
!       enddo
!       enddo
!
!
!       res(n) = mtrs
!
!
!        enddo
!
!
! end subroutine qqq_ttqqqq_dip_int
!
!
!
! subroutine qbb_ttqqqq_dip_int(p,res)
!   real(dp), intent(in) :: p(4,5)
!   real(dp), intent(out) :: res(2)
!   complex(dp) :: cres
!   type(Particle) :: ExtParticles1(1:5),ExtParticles2(1:5),ExtParticles3(1:5)
!   type(TreeProcess) :: TreeAmpsDip(1:4),TreeAmpsDipg(1:4),TreeAmpsDipq(1:6)
!   complex(dp) :: Bm(1:4),Bq(1:6)
!   complex(dp) :: AM(1:4,1:4)
!   integer :: i1,i2,i3,i4,hel(5),i,j,i5,i6,emi,n, in1, in2,iTree
!   real(dp) :: dipsoft, dipfini, dipplus
!   real(dp) :: z, C(4,4), fx1(-5:5), fx2(-5:5),fx1z(-5:5),fx2z(-5:5)
!   real(dp) :: Cq(6,6)
!   real(dp) :: mtrs
!   real(dp) :: q(4,5)
!   real(dp) :: pd1(4,5), pd2(4,5)
!   character :: diptype*2
!
!       res = zero
!       cres = (zero,zero)
! !---- there are two dipoles here
!
!       pd1 = p
!
!       pd1(:,3) = p(:,5)
!       pd1(:,5) = p(:,3)
!
!
!       pd2 = p
!
!       pd2(:,3) = p(:,5)
!       pd2(:,4) = p(:,3)
!       pd2(:,5) = p(:,4)
!
! !      init external particles
!       call InitProcess_TbTQbQG(ExtParticles1(1:5))
! !      init tree processes for 0-> tb t bq q g g
!        call InitTrees(4,1,4,TreeAmpsDip)
!
!       TreeAmpsDip(1)%PartRef(1:5) = (/1,2,3,4,5/)
!       TreeAmpsDip(2)%PartRef(1:5) = (/1,5,2,3,4/)
!       TreeAmpsDip(3)%PartRef(1:5) = (/1,2,5,3,4/)
!       TreeAmpsDip(4)%PartRef(1:5) = (/1,2,3,5,4/)
!
!
!       do iTree=1,4
!       call LinkTreeParticles(TreeAmpsDip(iTree),ExtParticles1(1:5))
!       enddo
!
!
!       do n = 1,2
!
!          if (n.eq.1) then
!             q = pd1
!          elseif(n.eq.2) then
!             q = pd2
!          endif
!
!
!       AM = (zero,zero) ! first, the ``amplitude matrix'' is set to zero
!
!       do i1=-1,1,2
!          do i2=-1,1,2
!             do i3 = -1,1,2
!                do i4 = -1,1,2
!                   do i5 = -1,1,2
!
!                  hel(1) = i1
!                  hel(2) = i2
!                  hel(3) = i3
!                  hel(4) = i4
!                  hel(5) = i5
!
!       call GenerateEventttqqg(q,hel,ExtParticles1(1:5))
!
!       do i6 = 1,4
!       call EvalTree2(TreeAmpsDip(i6),Bm(i6))
!       enddo
!
!        do i=1,4
!           do j=1,4
!        AM(i,j) = AM(i,j) + Bm(i)*conjg(Bm(j)) !color correlation matrix
!          enddo
!        enddo
!
!
!                   enddo
!                 enddo
!               enddo
!             enddo
!           enddo
!
!
!         call cc_ttqqg(C)
!
!         mtrs = zero
!
!       do i=1,4
!       do j=1,4
!     mtrs = mtrs + C(i,j)*real(AM(i,j),dp)
!       enddo
!       enddo
!
!
!       res(n) = mtrs
!
!
!        enddo
!
!
! end subroutine qbb_ttqqqq_dip_int
!




end  module


