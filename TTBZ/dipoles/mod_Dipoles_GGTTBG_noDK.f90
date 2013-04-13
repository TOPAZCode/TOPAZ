#DEFINE _CHECK_DIPOLE_MOMMAP 0
#DEFINE _APPLY_CUTS 0

      MODULE ModDipoles_GGTTBG_noDK
      implicit none


      double precision, private, parameter :: NCol=3d0
      double precision, private, parameter :: TR=1d0
      double precision, private, parameter :: CA=2d0*TR*NCol
      double precision, private, parameter :: CF=TR*(NCol**2-1d0)/NCol

      type :: Dipole
         double precision DipoleValue
         double precision MomTd(0:3,1:4)
         double precision Mass2Td(1:4)
         double precision alphaCut
      end type

      type(Dipole),private :: TheDipoles(1:12)
      double precision, private :: yRndDK(1:8),r_sc
      integer, private :: NBin(1:20)

      contains







      FUNCTION applyDipoleCut(MomTd)
      implicit none
      logical :: applyDipoleCut
      double precision MomTd(0:3,1:4)



      applyDipoleCut=.false.
      return
      END FUNCTION




      SUBROUTINE EvalDipoles_GGTTBG_noDK(Mom,Mass2,yRnd,Wgt,ResultSum)
!     evaluate dipole contributions
      use ModKinematics
      use ModMisc
      use ModParameters
      implicit none
      double precision Mom(0:3,1:5),Mass2(1:5),Wgt
      double precision, parameter :: DipMinus=-1d0*0.5d0
      integer :: NHisto
      double precision ResultSum,yRnd(1:8)
      double precision RunFactor

      yRndDK(1:8)=yRnd(1:8)

!     evaluate all dipoles and do binning:
      ResultSum = (0d0)
      r_sc = 0d0
      TheDipoles(1)%alphaCut = alpha_ii
      call evalDipole1_152(Mom,Mass2,TheDipoles(1))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(1)%DipoleValue = DipMinus * TheDipoles(1)%DipoleValue * Wgt * RunFactor**3
! IF( ObsSet.EQ.8 ) THEN
!       call intoHisto(1,WhichBin(1,r_sc),r_sc*TheDipoles(1)%DipoleValue)
!       call intoHisto(2,WhichBin(2,r_sc),TheDipoles(1)%DipoleValue)
!       do NHisto=3,NumHistograms
!           call intoHisto(NHisto,NBin(NHisto),TheDipoles(1)%DipoleValue)
!       enddo
! ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(1)%DipoleValue)
      enddo
! ENDIF
      ResultSum = ResultSum + TheDipoles(1)%DipoleValue


      r_sc = 0d0
      TheDipoles(2)%alphaCut = alpha_if
      call evalDipole2_153(Mom,Mass2,TheDipoles(2))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(2)%DipoleValue = DipMinus * TheDipoles(2)%DipoleValue * Wgt * RunFactor**3
! IF( ObsSet.EQ.8 ) THEN
!       call intoHisto(1,WhichBin(1,r_sc),r_sc*TheDipoles(2)%DipoleValue)
!       call intoHisto(2,WhichBin(2,r_sc),TheDipoles(2)%DipoleValue)
!       do NHisto=3,NumHistograms
!           call intoHisto(NHisto,NBin(NHisto),TheDipoles(2)%DipoleValue)
!       enddo
! ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(2)%DipoleValue)
      enddo
! ENDIF
      ResultSum = ResultSum + TheDipoles(2)%DipoleValue


      r_sc = 0d0
      TheDipoles(3)%alphaCut = alpha_if
      call evalDipole3_154(Mom,Mass2,TheDipoles(3))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(3)%DipoleValue = DipMinus * TheDipoles(3)%DipoleValue * Wgt * RunFactor**3
! IF( ObsSet.EQ.8 ) THEN
!       call intoHisto(1,WhichBin(1,r_sc),r_sc*TheDipoles(3)%DipoleValue)
!       call intoHisto(2,WhichBin(2,r_sc),TheDipoles(3)%DipoleValue)
!       do NHisto=3,NumHistograms
!           call intoHisto(NHisto,NBin(NHisto),TheDipoles(3)%DipoleValue)
!       enddo
! ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(3)%DipoleValue)
      enddo
! ENDIF
      ResultSum = ResultSum + TheDipoles(3)%DipoleValue


      r_sc = 0d0
      TheDipoles(4)%alphaCut = alpha_ii
      call evalDipole4_251(Mom,Mass2,TheDipoles(4))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(4)%DipoleValue = DipMinus * TheDipoles(4)%DipoleValue * Wgt * RunFactor**3
! IF( ObsSet.EQ.8 ) THEN
!       call intoHisto(1,WhichBin(1,r_sc),r_sc*TheDipoles(4)%DipoleValue)
!       call intoHisto(2,WhichBin(2,r_sc),TheDipoles(4)%DipoleValue)
!       do NHisto=3,NumHistograms
!           call intoHisto(NHisto,NBin(NHisto),TheDipoles(4)%DipoleValue)
!       enddo
! ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(4)%DipoleValue)
      enddo
! ENDIF
      ResultSum = ResultSum + TheDipoles(4)%DipoleValue


      r_sc = 0d0
      TheDipoles(5)%alphaCut = alpha_if
      call evalDipole5_253(Mom,Mass2,TheDipoles(5))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(5)%DipoleValue = DipMinus * TheDipoles(5)%DipoleValue * Wgt * RunFactor**3
! IF( ObsSet.EQ.8 ) THEN
!       call intoHisto(1,WhichBin(1,r_sc),r_sc*TheDipoles(5)%DipoleValue)
!       call intoHisto(2,WhichBin(2,r_sc),TheDipoles(5)%DipoleValue)
!       do NHisto=3,NumHistograms
!           call intoHisto(NHisto,NBin(NHisto),TheDipoles(5)%DipoleValue)
!       enddo
! ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(5)%DipoleValue)
      enddo
! ENDIF
      ResultSum = ResultSum + TheDipoles(5)%DipoleValue


      r_sc = 0d0
      TheDipoles(6)%alphaCut = alpha_if
      call evalDipole6_254(Mom,Mass2,TheDipoles(6))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(6)%DipoleValue = DipMinus * TheDipoles(6)%DipoleValue * Wgt * RunFactor**3
! IF( ObsSet.EQ.8 ) THEN
!       call intoHisto(1,WhichBin(1,r_sc),r_sc*TheDipoles(6)%DipoleValue)
!       call intoHisto(2,WhichBin(2,r_sc),TheDipoles(6)%DipoleValue)
!       do NHisto=3,NumHistograms
!           call intoHisto(NHisto,NBin(NHisto),TheDipoles(6)%DipoleValue)
!       enddo
! ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(6)%DipoleValue)
      enddo
! ENDIF
      ResultSum = ResultSum + TheDipoles(6)%DipoleValue


      r_sc = 0d0
      TheDipoles(7)%alphaCut = alpha_fi
      call evalDipole7_351(Mom,Mass2,TheDipoles(7))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(7)%DipoleValue = DipMinus * TheDipoles(7)%DipoleValue * Wgt * RunFactor**3
! IF( ObsSet.EQ.8 ) THEN
!       call intoHisto(1,WhichBin(1,r_sc),r_sc*TheDipoles(7)%DipoleValue)
!       call intoHisto(2,WhichBin(2,r_sc),TheDipoles(7)%DipoleValue)
!       do NHisto=3,NumHistograms
!           call intoHisto(NHisto,NBin(NHisto),TheDipoles(7)%DipoleValue)
!       enddo
! ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(7)%DipoleValue)
      enddo
! ENDIF
      ResultSum = ResultSum + TheDipoles(7)%DipoleValue


      r_sc = 0d0
      TheDipoles(8)%alphaCut = alpha_fi
      call evalDipole8_352(Mom,Mass2,TheDipoles(8))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(8)%DipoleValue = DipMinus * TheDipoles(8)%DipoleValue * Wgt * RunFactor**3
! IF( ObsSet.EQ.8 ) THEN
!       call intoHisto(1,WhichBin(1,r_sc),r_sc*TheDipoles(8)%DipoleValue)
!       call intoHisto(2,WhichBin(2,r_sc),TheDipoles(8)%DipoleValue)
!       do NHisto=3,NumHistograms
!           call intoHisto(NHisto,NBin(NHisto),TheDipoles(8)%DipoleValue)
!       enddo
! ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(8)%DipoleValue)
      enddo
! ENDIF
      ResultSum = ResultSum + TheDipoles(8)%DipoleValue


      r_sc = 0d0
      TheDipoles(9)%alphaCut = alpha_ff
      call evalDipole9_354(Mom,Mass2,TheDipoles(9))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(9)%DipoleValue = DipMinus * TheDipoles(9)%DipoleValue * Wgt * RunFactor**3
! IF( ObsSet.EQ.8 ) THEN
!       call intoHisto(1,WhichBin(1,r_sc),r_sc*TheDipoles(9)%DipoleValue)
!       call intoHisto(2,WhichBin(2,r_sc),TheDipoles(9)%DipoleValue)
!       do NHisto=3,NumHistograms
!           call intoHisto(NHisto,NBin(NHisto),TheDipoles(9)%DipoleValue)
!       enddo
! ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(9)%DipoleValue)
      enddo
! ENDIF
      ResultSum = ResultSum + TheDipoles(9)%DipoleValue


      r_sc = 0d0
      TheDipoles(10)%alphaCut = alpha_fi
      call evalDipole10_451(Mom,Mass2,TheDipoles(10))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(10)%DipoleValue = DipMinus * TheDipoles(10)%DipoleValue * Wgt * RunFactor**3
! IF( ObsSet.EQ.8 ) THEN
!       call intoHisto(1,WhichBin(1,r_sc),r_sc*TheDipoles(10)%DipoleValue)
!       call intoHisto(2,WhichBin(2,r_sc),TheDipoles(10)%DipoleValue)
!       do NHisto=3,NumHistograms
!           call intoHisto(NHisto,NBin(NHisto),TheDipoles(10)%DipoleValue)
!       enddo
! ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(10)%DipoleValue)
      enddo
! ENDIF
      ResultSum = ResultSum + TheDipoles(10)%DipoleValue


      r_sc = 0d0
      TheDipoles(11)%alphaCut = alpha_fi
      call evalDipole11_452(Mom,Mass2,TheDipoles(11))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(11)%DipoleValue = DipMinus * TheDipoles(11)%DipoleValue * Wgt * RunFactor**3
! IF( ObsSet.EQ.8 ) THEN
!       call intoHisto(1,WhichBin(1,r_sc),r_sc*TheDipoles(11)%DipoleValue)
!       call intoHisto(2,WhichBin(2,r_sc),TheDipoles(11)%DipoleValue)
!       do NHisto=3,NumHistograms
!           call intoHisto(NHisto,NBin(NHisto),TheDipoles(11)%DipoleValue)
!       enddo
! ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(11)%DipoleValue)
      enddo
! ENDIF
      ResultSum = ResultSum + TheDipoles(11)%DipoleValue


      r_sc = 0d0
      TheDipoles(12)%alphaCut = alpha_ff
      call evalDipole12_453(Mom,Mass2,TheDipoles(12))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(12)%DipoleValue = DipMinus * TheDipoles(12)%DipoleValue * Wgt * RunFactor**3
! IF( ObsSet.EQ.8 ) THEN
!       call intoHisto(1,WhichBin(1,r_sc),r_sc*TheDipoles(12)%DipoleValue)
!       call intoHisto(2,WhichBin(2,r_sc),TheDipoles(12)%DipoleValue)
!       do NHisto=3,NumHistograms
!           call intoHisto(NHisto,NBin(NHisto),TheDipoles(12)%DipoleValue)
!       enddo
! ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(12)%DipoleValue)
      enddo
! ENDIF
      ResultSum = ResultSum + TheDipoles(12)%DipoleValue



      return
      END SUBROUTINE





      SUBROUTINE check_MomTd(e,i,s,MomTd,Mass2Td)
!     check on-shellness and momentum conservation for tilde momenta
      use modMisc
      implicit none
      integer e,i,s
      double precision MomTd(0:3,1:4),Mass2Td(1:4)
      double precision MomCons

      write(*,*) ''
      write(*,'(A,I1,I1,I1)') 'checking momenta mapping in dipole D_',e,i,s
      write(*,'(3X,A,1PE17.7)') "OSness: ", ((MomTd(0:3,1).dot.MomTd(0:3,1))-Mass2Td(1))/MomTd(0,1)**2
      write(*,'(3X,A,1PE17.7)') "OSness: ", ((MomTd(0:3,2).dot.MomTd(0:3,2))-Mass2Td(2))/MomTd(0,2)**2
      write(*,'(3X,A,1PE17.7)') "OSness: ", ((MomTd(0:3,3).dot.MomTd(0:3,3))-Mass2Td(3))/MomTd(0,3)**2
      write(*,'(3X,A,1PE17.7)') "OSness: ", ((MomTd(0:3,4).dot.MomTd(0:3,4))-Mass2Td(4))/MomTd(0,4)**2
      MomCons = (MomTd(0,1)+MomTd(0,2)-MomTd(0,3)-MomTd(0,4))/MomTd(0,1)
      write(*,'(3X,A,1PE17.7)') "EM-cons.: ", MomCons
      MomCons = (MomTd(1,1)+MomTd(1,2)-MomTd(1,3)-MomTd(1,4))/MomTd(0,1)
      write(*,'(3X,A,1PE17.7)') "EM-cons.: ", MomCons
      MomCons = (MomTd(2,1)+MomTd(2,2)-MomTd(2,3)-MomTd(2,4))/MomTd(0,1)
      write(*,'(3X,A,1PE17.7)') "EM-cons.: ", MomCons
      MomCons = (MomTd(3,1)+MomTd(3,2)-MomTd(3,3)-MomTd(3,4))/MomTd(0,1)
      write(*,'(3X,A,1PE17.7)') "EM-cons.: ", MomCons

      return
      END SUBROUTINE





      SUBROUTINE evalDipole1_152(Mom,Mass2,TheDipole)
!     IS-emitter-IS-spectater dipole 1
      use modMisc
      use modParameters
      implicit none
      type(Dipole) :: TheDipole
      double precision Mom(0:3,1:5),Mass2(1:5)
      double precision s15,s12,s25
      double precision LeadSing,SqAmp
      double precision x,v,Fac1,Fac2,K(0:3),KTd(0:3),KKTdSum(0:3)
      double precision Split_A,Split_B,Split_V(0:3)

!     kinematic variables
      s15 = 2d0*(Mom(0:3,1)).dot.(Mom(0:3,5))
      s12 = 2d0*(Mom(0:3,1)).dot.(Mom(0:3,2))
      s25 = 2d0*(Mom(0:3,5)).dot.(Mom(0:3,2))
      x = 1d0 - (s15+s25)/s12
      v = s15/s12
      if( TheDipole%AlphaCut.lt.v ) then
         TheDipole%DipoleValue = 0d0
         return
      endif
      LeadSing = -1d0/s15/x

!     momentum mapping
      TheDipole%MomTd(0:3,1) = x*Mom(0:3,1)
      K(0:3) = Mom(0:3,1) + Mom(0:3,2) - Mom(0:3,5)
      KTd(0:3) = TheDipole%MomTd(0:3,1) + Mom(0:3,2)
      KKTdSum(0:3) = K(0:3) + KTd(0:3)
      TheDipole%MomTd(0:3,2) = Mom(0:3,2)
      Fac1 = 2d0*((Mom(0:3,3)).dot.(KKTdSum))/(KKTdSum.dot.KKTdSum)
      Fac2 = 2d0*((Mom(0:3,3)).dot.(K))/(K.dot.K)
      TheDipole%MomTd(0:3,3) = Mom(0:3,3)-Fac1*KKTdSum(0:3)+Fac2*KTd(0:3)
      Fac1 = 2d0*((Mom(0:3,4)).dot.(KKTdSum))/(KKTdSum.dot.KKTdSum)
      Fac2 = 2d0*((Mom(0:3,4)).dot.(K))/(K.dot.K)
      TheDipole%MomTd(0:3,4) = Mom(0:3,4)-Fac1*KKTdSum(0:3)+Fac2*KTd(0:3)
      TheDipole%Mass2Td(1) = Mass2(1)
      TheDipole%Mass2Td(2) = Mass2(2)
      TheDipole%Mass2Td(3) = Mass2(3)
      TheDipole%Mass2Td(4) = Mass2(4)

#IF (_CHECK_DIPOLE_MOMMAP)
      call check_MomTd(1,5,2,TheDipole%MomTd,TheDipole%Mass2Td)
#ENDIF

      if( applyDipoleCut(TheDipole%MomTd).eq. .true. ) then
         TheDipole%DipoleValue = 0d0
         return
      endif

!     splitting 2: Gluon_I --> Gluon_I + Gluon_F
      Split_A =-4d0*alpha_s4Pi * ( x/(1d0-x)+x*(1d0-x) )
      Split_B = 4d0*alpha_s4Pi * ( (1d0/x-1d0)*(2d0*s12/s15/s25) )
      Split_V(0:3) = Mom(0:3,5) - s15/s12*Mom(0:3,2)

!     color and spin correlated tree process
      SqAmp = Tree_GG_TTb_152_2(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

      TheDipole%DipoleValue = LeadSing * SqAmp
      return
      END SUBROUTINE





      SUBROUTINE evalDipole2_153(Mom,Mass2,TheDipole)
!     IS-emitter-FS-spectater dipole 2
      use modMisc
      use modParameters
      implicit none
      type(Dipole) :: TheDipole
      double precision s15,s13,s35
      double precision Mom(0:3,1:5),Mass2(1:5)
      double precision x,u,LeadSing,SplitF,SqAmp
      double precision Split_A,Split_B,Split_V(0:3)

!     kinematic variables
      s15 = 2d0*(Mom(0:3,1)).dot.(Mom(0:3,5))
      s13 = 2d0*(Mom(0:3,1)).dot.(Mom(0:3,3))
      s35 = 2d0*(Mom(0:3,5)).dot.(Mom(0:3,3))

      x = 1d0-s35/(s13+s15)
      u = s15/(s15+s13)
      if( TheDipole%AlphaCut.lt.u ) then
         TheDipole%DipoleValue = 0d0
         return
      endif
      LeadSing = -1d0/s15/x

!     momentum mapping
      TheDipole%MomTd(0:3,1) = x*Mom(0:3,1)
      TheDipole%MomTd(0:3,2) = Mom(0:3,2)
      TheDipole%MomTd(0:3,3) = Mom(0:3,5) + Mom(0:3,3) - (1d0-x)*Mom(0:3,1)
      TheDipole%MomTd(0:3,4) = Mom(0:3,4)
      TheDipole%Mass2Td(1) = Mass2(1)
      TheDipole%Mass2Td(2) = Mass2(2)
      TheDipole%Mass2Td(3) = Mass2(3)
      TheDipole%Mass2Td(4) = Mass2(4)

#IF (_CHECK_DIPOLE_MOMMAP)
      call check_MomTd(1,5,3,TheDipole%MomTd,TheDipole%Mass2Td)
#ENDIF

      if( applyDipoleCut(TheDipole%MomTd).eq. .true. ) then
         TheDipole%DipoleValue = 0d0
         return
      endif

!     splitting 2: Gluon_I --> Gluon_I + Gluon_F
      Split_A =-4d0*alpha_s4Pi * ( 1d0/(1d0-x+u)-1d0+x*(1d0-x) )
      Split_B = 4d0*alpha_s4Pi * ( (1d0/x-1d0)*(2d0*u*(1d0-u)/s35))
      Split_V(0:3) = Mom(0:3,5)/u - Mom(0:3,3)/(1d0-u)

!     color and spin correlated tree process
      SqAmp = Tree_GG_TTb_153_2(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

      TheDipole%DipoleValue = LeadSing * SqAmp
      return
      END SUBROUTINE





      SUBROUTINE evalDipole3_154(Mom,Mass2,TheDipole)
!     IS-emitter-FS-spectater dipole 3
      use modMisc
      use modParameters
      implicit none
      type(Dipole) :: TheDipole
      double precision s15,s14,s45
      double precision Mom(0:3,1:5),Mass2(1:5)
      double precision x,u,LeadSing,SplitF,SqAmp
      double precision Split_A,Split_B,Split_V(0:3)

!     kinematic variables
      s15 = 2d0*(Mom(0:3,1)).dot.(Mom(0:3,5))
      s14 = 2d0*(Mom(0:3,1)).dot.(Mom(0:3,4))
      s45 = 2d0*(Mom(0:3,5)).dot.(Mom(0:3,4))

      x = 1d0-s45/(s14+s15)
      u = s15/(s15+s14)
      if( TheDipole%AlphaCut.lt.u ) then
         TheDipole%DipoleValue = 0d0
         return
      endif
      LeadSing = -1d0/s15/x

!     momentum mapping
      TheDipole%MomTd(0:3,1) = x*Mom(0:3,1)
      TheDipole%MomTd(0:3,2) = Mom(0:3,2)
      TheDipole%MomTd(0:3,3) = Mom(0:3,3)
      TheDipole%MomTd(0:3,4) = Mom(0:3,5) + Mom(0:3,4) - (1d0-x)*Mom(0:3,1)
      TheDipole%Mass2Td(1) = Mass2(1)
      TheDipole%Mass2Td(2) = Mass2(2)
      TheDipole%Mass2Td(3) = Mass2(3)
      TheDipole%Mass2Td(4) = Mass2(4)

#IF (_CHECK_DIPOLE_MOMMAP)
      call check_MomTd(1,5,4,TheDipole%MomTd,TheDipole%Mass2Td)
#ENDIF

      if( applyDipoleCut(TheDipole%MomTd).eq. .true. ) then
         TheDipole%DipoleValue = 0d0
         return
      endif

!     splitting 2: Gluon_I --> Gluon_I + Gluon_F
      Split_A =-4d0*alpha_s4Pi * ( 1d0/(1d0-x+u)-1d0+x*(1d0-x) )
      Split_B = 4d0*alpha_s4Pi * ( (1d0/x-1d0)*(2d0*u*(1d0-u)/s45))
      Split_V(0:3) = Mom(0:3,5)/u - Mom(0:3,4)/(1d0-u)

!     color and spin correlated tree process
      SqAmp = Tree_GG_TTb_154_2(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

      TheDipole%DipoleValue = LeadSing * SqAmp
      return
      END SUBROUTINE





      SUBROUTINE evalDipole4_251(Mom,Mass2,TheDipole)
!     IS-emitter-IS-spectater dipole 4
      use modMisc
      use modParameters
      implicit none
      type(Dipole) :: TheDipole
      double precision Mom(0:3,1:5),Mass2(1:5)
      double precision s25,s12,s15
      double precision LeadSing,SqAmp
      double precision x,v,Fac1,Fac2,K(0:3),KTd(0:3),KKTdSum(0:3)
      double precision Split_A,Split_B,Split_V(0:3)

!     kinematic variables
      s25 = 2d0*(Mom(0:3,2)).dot.(Mom(0:3,5))
      s12 = 2d0*(Mom(0:3,2)).dot.(Mom(0:3,1))
      s15 = 2d0*(Mom(0:3,5)).dot.(Mom(0:3,1))
      x = 1d0 - (s25+s15)/s12
      v = s25/s12
      if( TheDipole%AlphaCut.lt.v ) then
         TheDipole%DipoleValue = 0d0
         return
      endif
      LeadSing = -1d0/s25/x

!     momentum mapping
      TheDipole%MomTd(0:3,1) = Mom(0:3,1)
      TheDipole%MomTd(0:3,2) = x*Mom(0:3,2)
      K(0:3) = Mom(0:3,2) + Mom(0:3,1) - Mom(0:3,5)
      KTd(0:3) = TheDipole%MomTd(0:3,2) + Mom(0:3,1)
      KKTdSum(0:3) = K(0:3) + KTd(0:3)
      Fac1 = 2d0*((Mom(0:3,3)).dot.(KKTdSum))/(KKTdSum.dot.KKTdSum)
      Fac2 = 2d0*((Mom(0:3,3)).dot.(K))/(K.dot.K)
      TheDipole%MomTd(0:3,3) = Mom(0:3,3)-Fac1*KKTdSum(0:3)+Fac2*KTd(0:3)
      Fac1 = 2d0*((Mom(0:3,4)).dot.(KKTdSum))/(KKTdSum.dot.KKTdSum)
      Fac2 = 2d0*((Mom(0:3,4)).dot.(K))/(K.dot.K)
      TheDipole%MomTd(0:3,4) = Mom(0:3,4)-Fac1*KKTdSum(0:3)+Fac2*KTd(0:3)
      TheDipole%Mass2Td(1) = Mass2(1)
      TheDipole%Mass2Td(2) = Mass2(2)
      TheDipole%Mass2Td(3) = Mass2(3)
      TheDipole%Mass2Td(4) = Mass2(4)

#IF (_CHECK_DIPOLE_MOMMAP)
      call check_MomTd(2,5,1,TheDipole%MomTd,TheDipole%Mass2Td)
#ENDIF

      if( applyDipoleCut(TheDipole%MomTd).eq. .true. ) then
         TheDipole%DipoleValue = 0d0
         return
      endif

!     splitting 2: Gluon_I --> Gluon_I + Gluon_F
      Split_A =-4d0*alpha_s4Pi * ( x/(1d0-x)+x*(1d0-x) )
      Split_B = 4d0*alpha_s4Pi * ( (1d0/x-1d0)*(2d0*s12/s25/s15) )
      Split_V(0:3) = Mom(0:3,5) - s25/s12*Mom(0:3,1)

!     color and spin correlated tree process
      SqAmp = Tree_GG_TTb_251_2(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

      TheDipole%DipoleValue = LeadSing * SqAmp
      return
      END SUBROUTINE





      SUBROUTINE evalDipole5_253(Mom,Mass2,TheDipole)
!     IS-emitter-FS-spectater dipole 5
      use modMisc
      use modParameters
      implicit none
      type(Dipole) :: TheDipole
      double precision s25,s23,s35
      double precision Mom(0:3,1:5),Mass2(1:5)
      double precision x,u,LeadSing,SplitF,SqAmp
      double precision Split_A,Split_B,Split_V(0:3)

!     kinematic variables
      s25 = 2d0*(Mom(0:3,2)).dot.(Mom(0:3,5))
      s23 = 2d0*(Mom(0:3,2)).dot.(Mom(0:3,3))
      s35 = 2d0*(Mom(0:3,5)).dot.(Mom(0:3,3))

      x = 1d0-s35/(s23+s25)
      u = s25/(s25+s23)
      if( TheDipole%AlphaCut.lt.u ) then
         TheDipole%DipoleValue = 0d0
         return
      endif
      LeadSing = -1d0/s25/x

!     momentum mapping
      TheDipole%MomTd(0:3,1) = Mom(0:3,1)
      TheDipole%MomTd(0:3,2) = x*Mom(0:3,2)
      TheDipole%MomTd(0:3,3) = Mom(0:3,5) + Mom(0:3,3) - (1d0-x)*Mom(0:3,2)
      TheDipole%MomTd(0:3,4) = Mom(0:3,4)
      TheDipole%Mass2Td(1) = Mass2(1)
      TheDipole%Mass2Td(2) = Mass2(2)
      TheDipole%Mass2Td(3) = Mass2(3)
      TheDipole%Mass2Td(4) = Mass2(4)

#IF (_CHECK_DIPOLE_MOMMAP)
      call check_MomTd(2,5,3,TheDipole%MomTd,TheDipole%Mass2Td)
#ENDIF

      if( applyDipoleCut(TheDipole%MomTd).eq. .true. ) then
         TheDipole%DipoleValue = 0d0
         return
      endif

!     splitting 2: Gluon_I --> Gluon_I + Gluon_F
      Split_A =-4d0*alpha_s4Pi * ( 1d0/(1d0-x+u)-1d0+x*(1d0-x) )
      Split_B = 4d0*alpha_s4Pi * ( (1d0/x-1d0)*(2d0*u*(1d0-u)/s35))
      Split_V(0:3) = Mom(0:3,5)/u - Mom(0:3,3)/(1d0-u)

!     color and spin correlated tree process
      SqAmp = Tree_GG_TTb_253_2(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

      TheDipole%DipoleValue = LeadSing * SqAmp
      return
      END SUBROUTINE





      SUBROUTINE evalDipole6_254(Mom,Mass2,TheDipole)
!     IS-emitter-FS-spectater dipole 6
      use modMisc
      use modParameters
      implicit none
      type(Dipole) :: TheDipole
      double precision s25,s24,s45
      double precision Mom(0:3,1:5),Mass2(1:5)
      double precision x,u,LeadSing,SplitF,SqAmp
      double precision Split_A,Split_B,Split_V(0:3)

!     kinematic variables
      s25 = 2d0*(Mom(0:3,2)).dot.(Mom(0:3,5))
      s24 = 2d0*(Mom(0:3,2)).dot.(Mom(0:3,4))
      s45 = 2d0*(Mom(0:3,5)).dot.(Mom(0:3,4))

      x = 1d0-s45/(s24+s25)
      u = s25/(s25+s24)
      if( TheDipole%AlphaCut.lt.u ) then
         TheDipole%DipoleValue = 0d0
         return
      endif
      LeadSing = -1d0/s25/x

!     momentum mapping
      TheDipole%MomTd(0:3,1) = Mom(0:3,1)
      TheDipole%MomTd(0:3,2) = x*Mom(0:3,2)
      TheDipole%MomTd(0:3,3) = Mom(0:3,3)
      TheDipole%MomTd(0:3,4) = Mom(0:3,5) + Mom(0:3,4) - (1d0-x)*Mom(0:3,2)
      TheDipole%Mass2Td(1) = Mass2(1)
      TheDipole%Mass2Td(2) = Mass2(2)
      TheDipole%Mass2Td(3) = Mass2(3)
      TheDipole%Mass2Td(4) = Mass2(4)

#IF (_CHECK_DIPOLE_MOMMAP)
      call check_MomTd(2,5,4,TheDipole%MomTd,TheDipole%Mass2Td)
#ENDIF

      if( applyDipoleCut(TheDipole%MomTd).eq. .true. ) then
         TheDipole%DipoleValue = 0d0
         return
      endif

!     splitting 2: Gluon_I --> Gluon_I + Gluon_F
      Split_A =-4d0*alpha_s4Pi * ( 1d0/(1d0-x+u)-1d0+x*(1d0-x) )
      Split_B = 4d0*alpha_s4Pi * ( (1d0/x-1d0)*(2d0*u*(1d0-u)/s45))
      Split_V(0:3) = Mom(0:3,5)/u - Mom(0:3,4)/(1d0-u)

!     color and spin correlated tree process
      SqAmp = Tree_GG_TTb_254_2(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

      TheDipole%DipoleValue = LeadSing * SqAmp
      return
      END SUBROUTINE





      SUBROUTINE evalDipole7_351(Mom,Mass2,TheDipole)
!     FS-emitter-IS-spectater dipole 7
      use modMisc
      use modParameters
      implicit none
      type(Dipole) :: TheDipole
      double precision s35,s13,s15
      double precision Mom(0:3,1:5),Mass2(1:5)
      double precision x,z,LeadSing,SplitF,SqAmp
      double precision Split_A,Split_B,Split_V(0:3)

!     kinematic variables
      s35 = 2d0*(Mom(0:3,3)).dot.(Mom(0:3,5))
      s13 = 2d0*(Mom(0:3,3)).dot.(Mom(0:3,1))
      s15 = 2d0*(Mom(0:3,5)).dot.(Mom(0:3,1))
      x = 1d0-s35/(s13+s15)
      z = s13/(s15+s13)
      LeadSing = -1d0/s35/x

!     momentum mapping
      TheDipole%MomTd(0:3,1) = x*Mom(0:3,1)
      TheDipole%MomTd(0:3,2) = Mom(0:3,2)
      TheDipole%MomTd(0:3,3) = Mom(0:3,3) + Mom(0:3,5) - (1d0-x)*Mom(0:3,1)
      TheDipole%MomTd(0:3,4) = Mom(0:3,4)
      TheDipole%Mass2Td(1) = Mass2(1)
      TheDipole%Mass2Td(2) = Mass2(2)
      TheDipole%Mass2Td(3) = Mass2(3)
      TheDipole%Mass2Td(4) = Mass2(4)

#IF (_CHECK_DIPOLE_MOMMAP)
      call check_MomTd(3,5,1,TheDipole%MomTd,TheDipole%Mass2Td)
#ENDIF

      if( applyDipoleCut(TheDipole%MomTd).eq. .true. ) then
         TheDipole%DipoleValue = 0d0
         return
      endif

!     splitting 5: A/Quark_F --> A/Quark_F + Gluon_F
      Split_A = 2d0*alpha_s4Pi*(2d0/(2d0-z-x)-1d0-z-2d0*Mass2(3)/s35)

!     color and spin correlated tree process
      SqAmp = Tree_GG_TTb_351_5(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

      TheDipole%DipoleValue = LeadSing * SqAmp
      return
      END SUBROUTINE





      SUBROUTINE evalDipole8_352(Mom,Mass2,TheDipole)
!     FS-emitter-IS-spectater dipole 8
      use modMisc
      use modParameters
      implicit none
      type(Dipole) :: TheDipole
      double precision s35,s23,s25
      double precision Mom(0:3,1:5),Mass2(1:5)
      double precision x,z,LeadSing,SplitF,SqAmp
      double precision Split_A,Split_B,Split_V(0:3)

!     kinematic variables
      s35 = 2d0*(Mom(0:3,3)).dot.(Mom(0:3,5))
      s23 = 2d0*(Mom(0:3,3)).dot.(Mom(0:3,2))
      s25 = 2d0*(Mom(0:3,5)).dot.(Mom(0:3,2))
      x = 1d0-s35/(s23+s25)
      z = s23/(s25+s23)
      LeadSing = -1d0/s35/x

!     momentum mapping
      TheDipole%MomTd(0:3,1) = Mom(0:3,1)
      TheDipole%MomTd(0:3,2) = x*Mom(0:3,2)
      TheDipole%MomTd(0:3,3) = Mom(0:3,3) + Mom(0:3,5) - (1d0-x)*Mom(0:3,2)
      TheDipole%MomTd(0:3,4) = Mom(0:3,4)
      TheDipole%Mass2Td(1) = Mass2(1)
      TheDipole%Mass2Td(2) = Mass2(2)
      TheDipole%Mass2Td(3) = Mass2(3)
      TheDipole%Mass2Td(4) = Mass2(4)

#IF (_CHECK_DIPOLE_MOMMAP)
      call check_MomTd(3,5,2,TheDipole%MomTd,TheDipole%Mass2Td)
#ENDIF

      if( applyDipoleCut(TheDipole%MomTd).eq. .true. ) then
         TheDipole%DipoleValue = 0d0
         return
      endif

!     splitting 5: A/Quark_F --> A/Quark_F + Gluon_F
      Split_A = 2d0*alpha_s4Pi*(2d0/(2d0-z-x)-1d0-z-2d0*Mass2(3)/s35)

!     color and spin correlated tree process
      SqAmp = Tree_GG_TTb_352_5(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

      TheDipole%DipoleValue = LeadSing * SqAmp
      return
      END SUBROUTINE





      SUBROUTINE evalDipole9_354(Mom,Mass2,TheDipole)
!     FS-emitter-FS-spectater dipole 9
      use modMisc
      use modParameters
      implicit none
      type(Dipole) :: TheDipole
      double precision s35,s34,s45
      double precision Mom(0:3,1:5),Mass2(1:5)
      double precision y,z,LeadSing,SplitF,SqAmp
      double precision v,Q2,mu2_i,mu2_k,Q(0:3)
      double precision MomFac1,MomFac2,MomFac3
      double precision Split_A,Split_B,Split_V(0:3)

!     kinematic variables
      s35 = 2d0*(Mom(0:3,3)).dot.(Mom(0:3,5))
      s34 = 2d0*(Mom(0:3,3)).dot.(Mom(0:3,4))
      s45 = 2d0*(Mom(0:3,5)).dot.(Mom(0:3,4))
      y = s35/(s35+s45+s34)
      z = s34/(s45+s34)
      Q2 = Mass2(3)+Mass2(4) + s35+s34+s45
      MomFac1 = dsqrt(((s35+s34+s45)**2-4d0*Mass2(3)*Mass2(4))/((s34+s45)**2-4d0*Mass2(4)*(Mass2(3)+s35)))
      MomFac2 = (0.5d0*(s34+s45)+Mass2(4))/Q2
      MomFac3 = MomFac2 + 0.5d0*s35/Q2
      LeadSing = -1d0/s35

!     momentum mapping
      Q(0:3) = Mom(0:3,3)+Mom(0:3,5)+Mom(0:3,4)
      TheDipole%MomTd(0:3,4) = MomFac1*(Mom(0:3,4)-MomFac2*Q(0:3)) + MomFac3*Q(0:3)
      TheDipole%MomTd(0:3,1) = Mom(0:3,1)
      TheDipole%MomTd(0:3,2) = Mom(0:3,2)
      TheDipole%MomTd(0:3,3) = Q(0:3) - TheDipole%MomTd(0:3,4)
      TheDipole%Mass2Td(1) = Mass2(1)
      TheDipole%Mass2Td(2) = Mass2(2)
      TheDipole%Mass2Td(3) = Mass2(3)
      TheDipole%Mass2Td(4) = Mass2(4)

#IF (_CHECK_DIPOLE_MOMMAP)
      call check_MomTd(3,5,4,TheDipole%MomTd,TheDipole%Mass2Td)
#ENDIF

      if( applyDipoleCut(TheDipole%MomTd).eq. .true. ) then
        TheDipole%DipoleValue = 0d0
        return
      endif

!     splitting 5: A/Quark_F --> A/Quark_F + Gluon_F
      mu2_i = Mass2(3)/dabs(Q2)
      mu2_k = Mass2(4)/dabs(Q2)
      v = (1d0-y)*dsqrt((1d0-(mu2_i-mu2_k)**2-2d0*(mu2_i+mu2_k))/( (2d0*mu2_k+(1d0-mu2_i-mu2_k)*(1d0-y))**2-4d0*mu2_k))
      Split_A = 2d0*alpha_s4Pi*(2d0/(1d0-z+y*z)-v*(1d0+z+2d0*Mass2(3)/s35))

!     color and spin correlated tree process
      SqAmp = Tree_GG_TTb_354_5(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

      TheDipole%DipoleValue = LeadSing * SqAmp
      return
      END SUBROUTINE





      SUBROUTINE evalDipole10_451(Mom,Mass2,TheDipole)
!     FS-emitter-IS-spectater dipole 10
      use modMisc
      use modParameters
      implicit none
      type(Dipole) :: TheDipole
      double precision s45,s14,s15
      double precision Mom(0:3,1:5),Mass2(1:5)
      double precision x,z,LeadSing,SplitF,SqAmp
      double precision Split_A,Split_B,Split_V(0:3)

!     kinematic variables
      s45 = 2d0*(Mom(0:3,4)).dot.(Mom(0:3,5))
      s14 = 2d0*(Mom(0:3,4)).dot.(Mom(0:3,1))
      s15 = 2d0*(Mom(0:3,5)).dot.(Mom(0:3,1))
      x = 1d0-s45/(s14+s15)
      z = s14/(s15+s14)
      LeadSing = -1d0/s45/x

!     momentum mapping
      TheDipole%MomTd(0:3,1) = x*Mom(0:3,1)
      TheDipole%MomTd(0:3,2) = Mom(0:3,2)
      TheDipole%MomTd(0:3,3) = Mom(0:3,3)
      TheDipole%MomTd(0:3,4) = Mom(0:3,4) + Mom(0:3,5) - (1d0-x)*Mom(0:3,1)
      TheDipole%Mass2Td(1) = Mass2(1)
      TheDipole%Mass2Td(2) = Mass2(2)
      TheDipole%Mass2Td(3) = Mass2(3)
      TheDipole%Mass2Td(4) = Mass2(4)

#IF (_CHECK_DIPOLE_MOMMAP)
      call check_MomTd(4,5,1,TheDipole%MomTd,TheDipole%Mass2Td)
#ENDIF

      if( applyDipoleCut(TheDipole%MomTd).eq. .true. ) then
         TheDipole%DipoleValue = 0d0
         return
      endif

!     splitting 5: A/Quark_F --> A/Quark_F + Gluon_F
      Split_A = 2d0*alpha_s4Pi*(2d0/(2d0-z-x)-1d0-z-2d0*Mass2(4)/s45)

!     color and spin correlated tree process
      SqAmp = Tree_GG_TTb_451_5(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

      TheDipole%DipoleValue = LeadSing * SqAmp
      return
      END SUBROUTINE





      SUBROUTINE evalDipole11_452(Mom,Mass2,TheDipole)
!     FS-emitter-IS-spectater dipole 11
      use modMisc
      use modParameters
      implicit none
      type(Dipole) :: TheDipole
      double precision s45,s24,s25
      double precision Mom(0:3,1:5),Mass2(1:5)
      double precision x,z,LeadSing,SplitF,SqAmp
      double precision Split_A,Split_B,Split_V(0:3)

!     kinematic variables
      s45 = 2d0*(Mom(0:3,4)).dot.(Mom(0:3,5))
      s24 = 2d0*(Mom(0:3,4)).dot.(Mom(0:3,2))
      s25 = 2d0*(Mom(0:3,5)).dot.(Mom(0:3,2))
      x = 1d0-s45/(s24+s25)
      z = s24/(s25+s24)
      LeadSing = -1d0/s45/x

!     momentum mapping
      TheDipole%MomTd(0:3,1) = Mom(0:3,1)
      TheDipole%MomTd(0:3,2) = x*Mom(0:3,2)
      TheDipole%MomTd(0:3,3) = Mom(0:3,3)
      TheDipole%MomTd(0:3,4) = Mom(0:3,4) + Mom(0:3,5) - (1d0-x)*Mom(0:3,2)
      TheDipole%Mass2Td(1) = Mass2(1)
      TheDipole%Mass2Td(2) = Mass2(2)
      TheDipole%Mass2Td(3) = Mass2(3)
      TheDipole%Mass2Td(4) = Mass2(4)

#IF (_CHECK_DIPOLE_MOMMAP)
      call check_MomTd(4,5,2,TheDipole%MomTd,TheDipole%Mass2Td)
#ENDIF

      if( applyDipoleCut(TheDipole%MomTd).eq. .true. ) then
         TheDipole%DipoleValue = 0d0
         return
      endif

!     splitting 5: A/Quark_F --> A/Quark_F + Gluon_F
      Split_A = 2d0*alpha_s4Pi*(2d0/(2d0-z-x)-1d0-z-2d0*Mass2(4)/s45)

!     color and spin correlated tree process
      SqAmp = Tree_GG_TTb_452_5(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

      TheDipole%DipoleValue = LeadSing * SqAmp
      return
      END SUBROUTINE





      SUBROUTINE evalDipole12_453(Mom,Mass2,TheDipole)
!     FS-emitter-FS-spectater dipole 12
      use modMisc
      use modParameters
      implicit none
      type(Dipole) :: TheDipole
      double precision s45,s34,s35
      double precision Mom(0:3,1:5),Mass2(1:5)
      double precision y,z,LeadSing,SplitF,SqAmp
      double precision v,Q2,mu2_i,mu2_k,Q(0:3)
      double precision MomFac1,MomFac2,MomFac3
      double precision Split_A,Split_B,Split_V(0:3)

!     kinematic variables
      s45 = 2d0*(Mom(0:3,4)).dot.(Mom(0:3,5))
      s34 = 2d0*(Mom(0:3,4)).dot.(Mom(0:3,3))
      s35 = 2d0*(Mom(0:3,5)).dot.(Mom(0:3,3))
      y = s45/(s45+s35+s34)
      z = s34/(s35+s34)
      Q2 = Mass2(4)+Mass2(3) + s45+s34+s35
      MomFac1 = dsqrt(((s45+s34+s35)**2-4d0*Mass2(4)*Mass2(3))/((s34+s35)**2-4d0*Mass2(3)*(Mass2(4)+s45)))
      MomFac2 = (0.5d0*(s34+s35)+Mass2(3))/Q2
      MomFac3 = MomFac2 + 0.5d0*s45/Q2
      LeadSing = -1d0/s45

!     momentum mapping
      Q(0:3) = Mom(0:3,4)+Mom(0:3,5)+Mom(0:3,3)
      TheDipole%MomTd(0:3,3) = MomFac1*(Mom(0:3,3)-MomFac2*Q(0:3)) + MomFac3*Q(0:3)
      TheDipole%MomTd(0:3,1) = Mom(0:3,1)
      TheDipole%MomTd(0:3,2) = Mom(0:3,2)
      TheDipole%MomTd(0:3,4) = Q(0:3) - TheDipole%MomTd(0:3,3)
      TheDipole%Mass2Td(1) = Mass2(1)
      TheDipole%Mass2Td(2) = Mass2(2)
      TheDipole%Mass2Td(3) = Mass2(3)
      TheDipole%Mass2Td(4) = Mass2(4)

#IF (_CHECK_DIPOLE_MOMMAP)
      call check_MomTd(4,5,3,TheDipole%MomTd,TheDipole%Mass2Td)
#ENDIF

      if( applyDipoleCut(TheDipole%MomTd).eq. .true. ) then
        TheDipole%DipoleValue = 0d0
        return
      endif

!     splitting 5: A/Quark_F --> A/Quark_F + Gluon_F
      mu2_i = Mass2(4)/dabs(Q2)
      mu2_k = Mass2(3)/dabs(Q2)
      v = (1d0-y)*dsqrt((1d0-(mu2_i-mu2_k)**2-2d0*(mu2_i+mu2_k))/( (2d0*mu2_k+(1d0-mu2_i-mu2_k)*(1d0-y))**2-4d0*mu2_k))
      Split_A = 2d0*alpha_s4Pi*(2d0/(1d0-z+y*z)-v*(1d0+z+2d0*Mass2(4)/s45))

!     color and spin correlated tree process
      SqAmp = Tree_GG_TTb_453_5(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

      TheDipole%DipoleValue = LeadSing * SqAmp
      return
      END SUBROUTINE








      FUNCTION Tree_GG_TTb_152_2(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modKinematics
      use modMyRecurrence
      implicit none
      double precision Split_A,Split_B,Split_V(0:3),MomDK(0:3,1:6),PSWgt2,PSWgt3
      double complex SqAmp,Amp(2,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double precision :: Msq_T_BW,Msq_W_ENU,PSWgt,MomTmp(0:3)
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      real(8),parameter::ColCorr(1:2,1:2)=(/-72.D0,0.D0,0.D0,-72.D0/)
      integer,target :: HelList(1:4,1:3)
      logical :: applyPSCut


      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
      HelList(3,1:3)=(/0,1,0/)
      HelList(4,1:3)=(/1,1,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

      PSWgt = 1d0
IF( TOPDECAYS.LT.0 ) THEN
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)
ENDIF

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin)
      if( applyPSCut ) then
          SqAmp = 0d0
          return
      endif
! IF( OBSSET.EQ.8 ) THEN
!       if( Collider.eq.1 ) then
!         r_sc = calc_rgg(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       else
!         r_sc = calc_rqq(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       endif
! ENDIF

      MomTmp(0:3)  = MomTd(0:3,3)! swap back
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)


!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

!      set particles
      call pol_mless(MomTd(0:3,1),+1,PolV(0:3,plus,1))
      call pol_mless(MomTd(0:3,1),-1,PolV(0:3,minus,1))

      call pol_mless(MomTd(0:3,2),+1,PolV(0:3,plus,2))
      call pol_mless(MomTd(0:3,2),-1,PolV(0:3,minus,2))
      Gluons(2)%Mom => MomTd(:,2)
      Gluons(2)%PartType => Glu_
      Gluons(2)%ExtRef => ExtRef
      Gluons(2)%Mass2 => Mass2Td(2)
      Gluons(2)%Mass  => MassTd(2)

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => Top_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ATop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,4
      Gluons(2)%Pol => PolV(:,HelList(iHel,1),2)
      Gluons(2)%Helicity => HelList(iHel,1)
      Quarks(1)%Pol => PolV(:,HelList(iHel,2),3)
      Quarks(1)%Helicity => HelList(iHel,2)
      Quarks(2)%Pol => PolV(:,HelList(iHel,3),4)
      Quarks(2)%Helicity => HelList(iHel,3)

!      calc currents
      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,1,0,0/) )
      Amp(1,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(1,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,0,0,1/) )
      Amp(2,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(2,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      EpsDotP = PolV(0:3,plus,1).dot.Split_V(0:3)
      SpinCorr_plus_plus  =-Split_A + Split_B*cdabs(EpsDotP)**2
      SpinCorr_plus_minus = Split_B * EpsDotP**2
      do col =1,2
      do colP=1,2
        SqAmp = SqAmp + ColCorr(colP,col) *(       SpinCorr_plus_plus * ( dconjg(Amp(colP,plus)) *Amp(col,plus)     &
                                                                        + dconjg(Amp(colP,minus))*Amp(col,minus) )  &
                                           +        SpinCorr_plus_minus * dconjg(Amp(colP,plus)) *Amp(col,minus)    &
                                           + dconjg(SpinCorr_plus_minus)* dconjg(Amp(colP,minus))*Amp(col,plus)  )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * 2d0*SqAmp * PSWgt

      return
      END FUNCTION






      FUNCTION Tree_GG_TTb_153_2(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      use modKinematics
      implicit none
      double precision Split_A,Split_B,Split_V(0:3),MomDK(0:3,1:6)
      double complex SqAmp,Amp(2,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double precision :: Msq_T_BW,Msq_W_ENU,PSWgt,MomTmp(0:3),PSWgt3,PSWgt2
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      real(8),parameter::ColCorr(1:2,1:2)=(/-64.D0,8.D0,8.D0,8.D0/)
      integer,target :: HelList(1:4,1:3)
      logical :: applyPSCut

      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
      HelList(3,1:3)=(/0,1,0/)
      HelList(4,1:3)=(/1,1,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

      PSWgt = 1d0
IF( TOPDECAYS.LT.0 ) THEN
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)
ENDIF


      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin)
      if( applyPSCut ) then
          SqAmp = 0d0
          return
      endif
! IF( OBSSET.EQ.8 ) THEN
!       if( Collider.eq.1 ) then
!         r_sc = calc_rgg(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       else
!         r_sc = calc_rqq(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       endif
! ENDIF

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)


!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

!      set particles
      call pol_mless(MomTd(0:3,1),+1,PolV(0:3,plus,1))
      call pol_mless(MomTd(0:3,1),-1,PolV(0:3,minus,1))

      call pol_mless(MomTd(0:3,2),+1,PolV(0:3,plus,2))
      call pol_mless(MomTd(0:3,2),-1,PolV(0:3,minus,2))
      Gluons(2)%Mom => MomTd(:,2)
      Gluons(2)%PartType => Glu_
      Gluons(2)%ExtRef => ExtRef
      Gluons(2)%Mass2 => Mass2Td(2)
      Gluons(2)%Mass  => MassTd(2)

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => Top_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ATop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,4
      Gluons(2)%Pol => PolV(:,HelList(iHel,1),2)
      Gluons(2)%Helicity => HelList(iHel,1)
      Quarks(1)%Pol => PolV(:,HelList(iHel,2),3)
      Quarks(1)%Helicity => HelList(iHel,2)
      Quarks(2)%Pol => PolV(:,HelList(iHel,3),4)
      Quarks(2)%Helicity => HelList(iHel,3)

!      calc currents
      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,1,0,0/) )
      Amp(1,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(1,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,0,0,1/) )
      Amp(2,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(2,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      EpsDotP = PolV(0:3,plus,1).dot.Split_V(0:3)
      SpinCorr_plus_plus  =-Split_A + Split_B*cdabs(EpsDotP)**2
      SpinCorr_plus_minus = Split_B * EpsDotP**2
      do col =1,2
      do colP=1,2
        SqAmp = SqAmp + ColCorr(colP,col) *(       SpinCorr_plus_plus * ( dconjg(Amp(colP,plus)) *Amp(col,plus)     &
                                                                        + dconjg(Amp(colP,minus))*Amp(col,minus) )  &
                                           +        SpinCorr_plus_minus * dconjg(Amp(colP,plus)) *Amp(col,minus)    &
                                           + dconjg(SpinCorr_plus_minus)* dconjg(Amp(colP,minus))*Amp(col,plus)  )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * 2d0*SqAmp * PSWgt

      return
      END FUNCTION






      FUNCTION Tree_GG_TTb_154_2(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modKinematics
      use modMyRecurrence
      implicit none
      double precision Split_A,Split_B,Split_V(0:3),MomDK(0:3,1:6)
      double complex SqAmp,Amp(2,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double precision :: Msq_T_BW,Msq_W_ENU,PSWgt,MomTmp(0:3),PSWgt3,PSWgt2
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      real(8),parameter::ColCorr(1:2,1:2)=(/8.D0,8.D0,8.D0,-64.D0/)
      integer,target :: HelList(1:4,1:3)
      logical :: applyPSCut

      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
      HelList(3,1:3)=(/0,1,0/)
      HelList(4,1:3)=(/1,1,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

      PSWgt = 1d0
IF( TOPDECAYS.LT.0 ) THEN
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)
ENDIF


      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin)   ! top and anti-top are in wrong order for this call, but for pT-distr. it is okay
      if( applyPSCut ) then
          SqAmp = 0d0
          return
      endif
! IF( OBSSET.EQ.8 ) THEN
!       if( Collider.eq.1 ) then
!         r_sc = calc_rgg(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       else
!         r_sc = calc_rqq(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       endif
! ENDIF

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)

!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

!      set particles
      call pol_mless(MomTd(0:3,1),+1,PolV(0:3,plus,1))
      call pol_mless(MomTd(0:3,1),-1,PolV(0:3,minus,1))

      call pol_mless(MomTd(0:3,2),+1,PolV(0:3,plus,2))
      call pol_mless(MomTd(0:3,2),-1,PolV(0:3,minus,2))
      Gluons(2)%Mom => MomTd(:,2)
      Gluons(2)%PartType => Glu_
      Gluons(2)%ExtRef => ExtRef
      Gluons(2)%Mass2 => Mass2Td(2)
      Gluons(2)%Mass  => MassTd(2)

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => Top_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ATop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,4
      Gluons(2)%Pol => PolV(:,HelList(iHel,1),2)
      Gluons(2)%Helicity => HelList(iHel,1)
      Quarks(1)%Pol => PolV(:,HelList(iHel,2),3)
      Quarks(1)%Helicity => HelList(iHel,2)
      Quarks(2)%Pol => PolV(:,HelList(iHel,3),4)
      Quarks(2)%Helicity => HelList(iHel,3)

!      calc currents
      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,1,0,0/) )
      Amp(1,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(1,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,0,0,1/) )
      Amp(2,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(2,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      EpsDotP = PolV(0:3,plus,1).dot.Split_V(0:3)
      SpinCorr_plus_plus  =-Split_A + Split_B*cdabs(EpsDotP)**2
      SpinCorr_plus_minus = Split_B * EpsDotP**2
      do col =1,2
      do colP=1,2
        SqAmp = SqAmp + ColCorr(colP,col) *(       SpinCorr_plus_plus * ( dconjg(Amp(colP,plus)) *Amp(col,plus)     &
                                                                        + dconjg(Amp(colP,minus))*Amp(col,minus) )  &
                                           +        SpinCorr_plus_minus * dconjg(Amp(colP,plus)) *Amp(col,minus)    &
                                           + dconjg(SpinCorr_plus_minus)* dconjg(Amp(colP,minus))*Amp(col,plus)  )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * 2d0*SqAmp * PSWgt

      return
      END FUNCTION






      FUNCTION Tree_GG_TTb_251_2(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      use modKinematics
      implicit none
      double precision Split_A,Split_B,Split_V(0:3),MomDK(0:3,1:6)
      double complex SqAmp,Amp(2,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double precision :: Msq_T_BW,Msq_W_ENU,PSWgt,MomTmp(0:3),PSWgt3,PSWgt2
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      real(8),parameter::ColCorr(1:2,1:2)=(/-72.D0,0.D0,0.D0,-72.D0/)
      integer,target :: HelList(1:4,1:3)
      logical :: applyPSCut

      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
      HelList(3,1:3)=(/0,1,0/)
      HelList(4,1:3)=(/1,1,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))


      PSWgt = 1d0
IF( TOPDECAYS.LT.0 ) THEN
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)
ENDIF


      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin)   ! top and anti-top are in wrong order for this call, but for pT-distr. it is okay
      if( applyPSCut ) then
          SqAmp = 0d0
          return
      endif
! IF( OBSSET.EQ.8 ) THEN
!       if( Collider.eq.1 ) then
!         r_sc = calc_rgg(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       else
!         r_sc = calc_rqq(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       endif
! ENDIF

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)

!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

!      set particles
      call pol_mless(MomTd(0:3,1),+1,PolV(0:3,plus,1))
      call pol_mless(MomTd(0:3,1),-1,PolV(0:3,minus,1))
      Gluons(1)%Mom => MomTd(:,1)
      Gluons(1)%PartType => Glu_
      Gluons(1)%ExtRef => ExtRef
      Gluons(1)%Mass2 => Mass2Td(1)
      Gluons(1)%Mass  => MassTd(1)

      call pol_mless(MomTd(0:3,2),+1,PolV(0:3,plus,2))
      call pol_mless(MomTd(0:3,2),-1,PolV(0:3,minus,2))

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => Top_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ATop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,4
      Gluons(1)%Pol => PolV(:,HelList(iHel,1),1)
      Gluons(1)%Helicity => HelList(iHel,1)
      Quarks(1)%Pol => PolV(:,HelList(iHel,2),3)
      Quarks(1)%Helicity => HelList(iHel,2)
      Quarks(2)%Pol => PolV(:,HelList(iHel,3),4)
      Quarks(2)%Helicity => HelList(iHel,3)

!      calc currents
      Res(0:3) = cur_g_2f( (/Gluons(1)/),(/Quarks(2),Quarks(1)/),(/2,0,0,1/) )
      Amp(1,plus)  = Res(0:3).dot.PolV(0:3,plus,2)
      Amp(1,minus) = Res(0:3).dot.PolV(0:3,minus,2)

      Res(0:3) = cur_g_2f( (/Gluons(1)/),(/Quarks(2),Quarks(1)/),(/2,1,0,0/) )
      Amp(2,plus)  = Res(0:3).dot.PolV(0:3,plus,2)
      Amp(2,minus) = Res(0:3).dot.PolV(0:3,minus,2)

      EpsDotP = PolV(0:3,plus,2).dot.Split_V(0:3)
      SpinCorr_plus_plus  =-Split_A + Split_B*cdabs(EpsDotP)**2
      SpinCorr_plus_minus = Split_B * EpsDotP**2
      do col =1,2
      do colP=1,2
        SqAmp = SqAmp + ColCorr(colP,col) *(       SpinCorr_plus_plus * ( dconjg(Amp(colP,plus)) *Amp(col,plus)     &
                                                                        + dconjg(Amp(colP,minus))*Amp(col,minus) )  &
                                           +        SpinCorr_plus_minus * dconjg(Amp(colP,plus)) *Amp(col,minus)    &
                                           + dconjg(SpinCorr_plus_minus)* dconjg(Amp(colP,minus))*Amp(col,plus)  )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * 2d0*SqAmp * PSWgt

      return
      END FUNCTION






      FUNCTION Tree_GG_TTb_253_2(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modKinematics
      use modMyRecurrence
      implicit none
      double precision Split_A,Split_B,Split_V(0:3),MomDK(0:3,1:6)
      double complex SqAmp,Amp(2,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double precision :: Msq_T_BW,Msq_W_ENU,PSWgt,MomTmp(0:3),PSWgt3,PSWgt2
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      real(8),parameter::ColCorr(1:2,1:2)=(/8.D0,8.D0,8.D0,-64.D0/)
      integer,target :: HelList(1:4,1:3)
      logical :: applyPSCut

      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
      HelList(3,1:3)=(/0,1,0/)
      HelList(4,1:3)=(/1,1,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

      PSWgt = 1d0
IF( TOPDECAYS.LT.0 ) THEN
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)
ENDIF


      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin)   ! top and anti-top are in wrong order for this call, but for pT-distr. it is okay
      if( applyPSCut ) then
          SqAmp = 0d0
          return
      endif
! IF( OBSSET.EQ.8 ) THEN
!       if( Collider.eq.1 ) then
!         r_sc = calc_rgg(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       else
!         r_sc = calc_rqq(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       endif
! ENDIF
      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)

!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

!      set particles
      call pol_mless(MomTd(0:3,1),+1,PolV(0:3,plus,1))
      call pol_mless(MomTd(0:3,1),-1,PolV(0:3,minus,1))
      Gluons(1)%Mom => MomTd(:,1)
      Gluons(1)%PartType => Glu_
      Gluons(1)%ExtRef => ExtRef
      Gluons(1)%Mass2 => Mass2Td(1)
      Gluons(1)%Mass  => MassTd(1)

      call pol_mless(MomTd(0:3,2),+1,PolV(0:3,plus,2))
      call pol_mless(MomTd(0:3,2),-1,PolV(0:3,minus,2))

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => Top_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ATop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,4
      Gluons(1)%Pol => PolV(:,HelList(iHel,1),1)
      Gluons(1)%Helicity => HelList(iHel,1)
      Quarks(1)%Pol => PolV(:,HelList(iHel,2),3)
      Quarks(1)%Helicity => HelList(iHel,2)
      Quarks(2)%Pol => PolV(:,HelList(iHel,3),4)
      Quarks(2)%Helicity => HelList(iHel,3)

!      calc currents
      Res(0:3) = cur_g_2f( (/Gluons(1)/),(/Quarks(2),Quarks(1)/),(/2,0,0,1/) )
      Amp(1,plus)  = Res(0:3).dot.PolV(0:3,plus,2)
      Amp(1,minus) = Res(0:3).dot.PolV(0:3,minus,2)

      Res(0:3) = cur_g_2f( (/Gluons(1)/),(/Quarks(2),Quarks(1)/),(/2,1,0,0/) )
      Amp(2,plus)  = Res(0:3).dot.PolV(0:3,plus,2)
      Amp(2,minus) = Res(0:3).dot.PolV(0:3,minus,2)

      EpsDotP = PolV(0:3,plus,2).dot.Split_V(0:3)
      SpinCorr_plus_plus  =-Split_A + Split_B*cdabs(EpsDotP)**2
      SpinCorr_plus_minus = Split_B * EpsDotP**2
      do col =1,2
      do colP=1,2
        SqAmp = SqAmp + ColCorr(colP,col) *(       SpinCorr_plus_plus * ( dconjg(Amp(colP,plus)) *Amp(col,plus)     &
                                                                        + dconjg(Amp(colP,minus))*Amp(col,minus) )  &
                                           +        SpinCorr_plus_minus * dconjg(Amp(colP,plus)) *Amp(col,minus)    &
                                           + dconjg(SpinCorr_plus_minus)* dconjg(Amp(colP,minus))*Amp(col,plus)  )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * 2d0*SqAmp * PSWgt

      return
      END FUNCTION






      FUNCTION Tree_GG_TTb_254_2(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modKinematics
      use modMyRecurrence
      implicit none
      double precision Split_A,Split_B,Split_V(0:3),MomDK(0:3,1:6)
      double complex SqAmp,Amp(2,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double precision :: Msq_T_BW,Msq_W_ENU,PSWgt,MomTmp(0:3),PSWgt3,PSWgt2
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      real(8),parameter::ColCorr(1:2,1:2)=(/-64.D0,8.D0,8.D0,8.D0/)
      integer,target :: HelList(1:4,1:3)
      logical :: applyPSCut

      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
      HelList(3,1:3)=(/0,1,0/)
      HelList(4,1:3)=(/1,1,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

      PSWgt = 1d0
IF( TOPDECAYS.LT.0 ) THEN
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)
ENDIF


      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin)   ! top and anti-top are in wrong order for this call, but for pT-distr. it is okay
      if( applyPSCut ) then
          SqAmp = 0d0
          return
      endif
! IF( OBSSET.EQ.8 ) THEN
!       if( Collider.eq.1 ) then
!         r_sc = calc_rgg(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       else
!         r_sc = calc_rqq(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       endif
! ENDIF

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)

!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

!      set particles
      call pol_mless(MomTd(0:3,1),+1,PolV(0:3,plus,1))
      call pol_mless(MomTd(0:3,1),-1,PolV(0:3,minus,1))
      Gluons(1)%Mom => MomTd(:,1)
      Gluons(1)%PartType => Glu_
      Gluons(1)%ExtRef => ExtRef
      Gluons(1)%Mass2 => Mass2Td(1)
      Gluons(1)%Mass  => MassTd(1)

      call pol_mless(MomTd(0:3,2),+1,PolV(0:3,plus,2))
      call pol_mless(MomTd(0:3,2),-1,PolV(0:3,minus,2))

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => Top_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ATop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,4
      Gluons(1)%Pol => PolV(:,HelList(iHel,1),1)
      Gluons(1)%Helicity => HelList(iHel,1)
      Quarks(1)%Pol => PolV(:,HelList(iHel,2),3)
      Quarks(1)%Helicity => HelList(iHel,2)
      Quarks(2)%Pol => PolV(:,HelList(iHel,3),4)
      Quarks(2)%Helicity => HelList(iHel,3)

!      calc currents
      Res(0:3) = cur_g_2f( (/Gluons(1)/),(/Quarks(2),Quarks(1)/),(/2,0,0,1/) )
      Amp(1,plus)  = Res(0:3).dot.PolV(0:3,plus,2)
      Amp(1,minus) = Res(0:3).dot.PolV(0:3,minus,2)

      Res(0:3) = cur_g_2f( (/Gluons(1)/),(/Quarks(2),Quarks(1)/),(/2,1,0,0/) )
      Amp(2,plus)  = Res(0:3).dot.PolV(0:3,plus,2)
      Amp(2,minus) = Res(0:3).dot.PolV(0:3,minus,2)

      EpsDotP = PolV(0:3,plus,2).dot.Split_V(0:3)
      SpinCorr_plus_plus  =-Split_A + Split_B*cdabs(EpsDotP)**2
      SpinCorr_plus_minus = Split_B * EpsDotP**2
      do col =1,2
      do colP=1,2
        SqAmp = SqAmp + ColCorr(colP,col) *(       SpinCorr_plus_plus * ( dconjg(Amp(colP,plus)) *Amp(col,plus)     &
                                                                        + dconjg(Amp(colP,minus))*Amp(col,minus) )  &
                                           +        SpinCorr_plus_minus * dconjg(Amp(colP,plus)) *Amp(col,minus)    &
                                           + dconjg(SpinCorr_plus_minus)* dconjg(Amp(colP,minus))*Amp(col,plus)  )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * 2d0*SqAmp * PSWgt

      return
      END FUNCTION






      FUNCTION Tree_GG_TTb_351_5(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modKinematics
      use modParameters
      use modMyRecurrence
      implicit none
      double precision Split_A,Split_B,Split_V(0:3),MomDK(0:3,1:6)
      double complex SqAmp,Amp(2,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double precision :: Msq_T_BW,Msq_W_ENU,PSWgt,MomTmp(0:3),PSWgt3,PSWgt2
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      real(8),parameter::ColCorr(1:2,1:2)=(/-64.D0,8.D0,8.D0,8.D0/)
      integer,target :: HelList(1:4,1:3)
      logical :: applyPSCut

      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
      HelList(3,1:3)=(/0,1,0/)
      HelList(4,1:3)=(/1,1,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

      PSWgt = 1d0
IF( TOPDECAYS.LT.0 ) THEN
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)
ENDIF


      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin)   ! top and anti-top are in wrong order for this call, but for pT-distr. it is okay
      if( applyPSCut ) then
          SqAmp = 0d0
          return
      endif
! IF( OBSSET.EQ.8 ) THEN
!       if( Collider.eq.1 ) then
!         r_sc = calc_rgg(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       else
!         r_sc = calc_rqq(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       endif
! ENDIF

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)

!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

!      set particles
      call pol_mless(MomTd(0:3,1),+1,PolV(0:3,plus,1))
      call pol_mless(MomTd(0:3,1),-1,PolV(0:3,minus,1))

      call pol_mless(MomTd(0:3,2),+1,PolV(0:3,plus,2))
      call pol_mless(MomTd(0:3,2),-1,PolV(0:3,minus,2))
      Gluons(2)%Mom => MomTd(:,2)
      Gluons(2)%PartType => Glu_
      Gluons(2)%ExtRef => ExtRef
      Gluons(2)%Mass2 => Mass2Td(2)
      Gluons(2)%Mass  => MassTd(2)

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => Top_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ATop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,4
      Gluons(2)%Pol => PolV(:,HelList(iHel,1),2)
      Gluons(2)%Helicity => HelList(iHel,1)
      Quarks(1)%Pol => PolV(:,HelList(iHel,2),3)
      Quarks(1)%Helicity => HelList(iHel,2)
      Quarks(2)%Pol => PolV(:,HelList(iHel,3),4)
      Quarks(2)%Helicity => HelList(iHel,3)

!      calc currents
      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,1,0,0/) )
      Amp(1,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(1,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,0,0,1/) )
      Amp(2,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(2,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      do col =1,2
      do colP=1,2
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus)) *Amp(col,plus)     &
                                                    + dconjg(Amp(colP,minus))*Amp(col,minus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * 2d0*SqAmp * PSWgt

      return
      END FUNCTION






      FUNCTION Tree_GG_TTb_352_5(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modKinematics
      use modParameters
      use modMyRecurrence
      implicit none
      double precision Split_A,Split_B,Split_V(0:3),MomDK(0:3,1:6)
      double complex SqAmp,Amp(2,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double precision :: Msq_T_BW,Msq_W_ENU,PSWgt,MomTmp(0:3),PSWgt3,PSWgt2
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      real(8),parameter::ColCorr(1:2,1:2)=(/8.D0,8.D0,8.D0,-64.D0/)
      integer,target :: HelList(1:4,1:3)
      logical :: applyPSCut

      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
      HelList(3,1:3)=(/0,1,0/)
      HelList(4,1:3)=(/1,1,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

      PSWgt = 1d0
IF( TOPDECAYS.LT.0 ) THEN
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)
ENDIF


      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin)   ! top and anti-top are in wrong order for this call, but for pT-distr. it is okay
      if( applyPSCut ) then
          SqAmp = 0d0
          return
      endif
! IF( OBSSET.EQ.8 ) THEN
!       if( Collider.eq.1 ) then
!         r_sc = calc_rgg(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       else
!         r_sc = calc_rqq(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       endif
! ENDIF

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)

!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

!      set particles
      call pol_mless(MomTd(0:3,1),+1,PolV(0:3,plus,1))
      call pol_mless(MomTd(0:3,1),-1,PolV(0:3,minus,1))

      call pol_mless(MomTd(0:3,2),+1,PolV(0:3,plus,2))
      call pol_mless(MomTd(0:3,2),-1,PolV(0:3,minus,2))
      Gluons(2)%Mom => MomTd(:,2)
      Gluons(2)%PartType => Glu_
      Gluons(2)%ExtRef => ExtRef
      Gluons(2)%Mass2 => Mass2Td(2)
      Gluons(2)%Mass  => MassTd(2)

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => Top_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ATop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,4
      Gluons(2)%Pol => PolV(:,HelList(iHel,1),2)
      Gluons(2)%Helicity => HelList(iHel,1)
      Quarks(1)%Pol => PolV(:,HelList(iHel,2),3)
      Quarks(1)%Helicity => HelList(iHel,2)
      Quarks(2)%Pol => PolV(:,HelList(iHel,3),4)
      Quarks(2)%Helicity => HelList(iHel,3)

!      calc currents
      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,1,0,0/) )
      Amp(1,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(1,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,0,0,1/) )
      Amp(2,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(2,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      do col =1,2
      do colP=1,2
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus)) *Amp(col,plus)     &
                                                    + dconjg(Amp(colP,minus))*Amp(col,minus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * 2d0*SqAmp * PSWgt

      return
      END FUNCTION






      FUNCTION Tree_GG_TTb_354_5(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modKinematics
      use modMyRecurrence
      implicit none
      double precision Split_A,Split_B,Split_V(0:3),MomDK(0:3,1:6)
      double complex SqAmp,Amp(2,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double precision :: Msq_T_BW,Msq_W_ENU,PSWgt,MomTmp(0:3),PSWgt3,PSWgt2
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      real(8),parameter::ColCorr(1:2,1:2)=(/-8.D0/9.D0,-80.D0/9.D0,-80.D0/9.D0,-8.D0/9.D0/)
      integer,target :: HelList(1:4,1:3)
      logical :: applyPSCut

      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
      HelList(3,1:3)=(/0,1,0/)
      HelList(4,1:3)=(/1,1,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

      PSWgt = 1d0
IF( TOPDECAYS.LT.0 ) THEN
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)
ENDIF


      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin)   ! top and anti-top are in wrong order for this call, but for pT-distr. it is okay
      if( applyPSCut ) then
          SqAmp = 0d0
          return
      endif
! IF( OBSSET.EQ.8 ) THEN
!       if( Collider.eq.1 ) then
!         r_sc = calc_rgg(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       else
!         r_sc = calc_rqq(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       endif
! ENDIF

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)

!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

!      set particles
      call pol_mless(MomTd(0:3,1),+1,PolV(0:3,plus,1))
      call pol_mless(MomTd(0:3,1),-1,PolV(0:3,minus,1))

      call pol_mless(MomTd(0:3,2),+1,PolV(0:3,plus,2))
      call pol_mless(MomTd(0:3,2),-1,PolV(0:3,minus,2))
      Gluons(2)%Mom => MomTd(:,2)
      Gluons(2)%PartType => Glu_
      Gluons(2)%ExtRef => ExtRef
      Gluons(2)%Mass2 => Mass2Td(2)
      Gluons(2)%Mass  => MassTd(2)

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => Top_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ATop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,4
      Gluons(2)%Pol => PolV(:,HelList(iHel,1),2)
      Gluons(2)%Helicity => HelList(iHel,1)
      Quarks(1)%Pol => PolV(:,HelList(iHel,2),3)
      Quarks(1)%Helicity => HelList(iHel,2)
      Quarks(2)%Pol => PolV(:,HelList(iHel,3),4)
      Quarks(2)%Helicity => HelList(iHel,3)

!      calc currents
      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,1,0,0/) )
      Amp(1,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(1,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,0,0,1/) )
      Amp(2,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(2,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      do col =1,2
      do colP=1,2
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus)) *Amp(col,plus)     &
                                                    + dconjg(Amp(colP,minus))*Amp(col,minus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * 2d0*SqAmp * PSWgt

      return
      END FUNCTION






      FUNCTION Tree_GG_TTb_451_5(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modKinematics
      use modMyRecurrence
      implicit none
      double precision Split_A,Split_B,Split_V(0:3),MomDK(0:3,1:6)
      double complex SqAmp,Amp(2,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double precision :: Msq_T_BW,Msq_W_ENU,PSWgt,MomTmp(0:3),PSWgt3,PSWgt2
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      real(8),parameter::ColCorr(1:2,1:2)=(/8.D0,8.D0,8.D0,-64.D0/)
      integer,target :: HelList(1:4,1:3)
      logical :: applyPSCut

      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
      HelList(3,1:3)=(/0,1,0/)
      HelList(4,1:3)=(/1,1,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

      PSWgt = 1d0
IF( TOPDECAYS.LT.0 ) THEN
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)
ENDIF


      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin)   ! top and anti-top are in wrong order for this call, but for pT-distr. it is okay
      if( applyPSCut ) then
          SqAmp = 0d0
          return
      endif
! IF( OBSSET.EQ.8 ) THEN
!       if( Collider.eq.1 ) then
!         r_sc = calc_rgg(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       else
!         r_sc = calc_rqq(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       endif
! ENDIF

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)

!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

!      set particles
      call pol_mless(MomTd(0:3,1),+1,PolV(0:3,plus,1))
      call pol_mless(MomTd(0:3,1),-1,PolV(0:3,minus,1))

      call pol_mless(MomTd(0:3,2),+1,PolV(0:3,plus,2))
      call pol_mless(MomTd(0:3,2),-1,PolV(0:3,minus,2))
      Gluons(2)%Mom => MomTd(:,2)
      Gluons(2)%PartType => Glu_
      Gluons(2)%ExtRef => ExtRef
      Gluons(2)%Mass2 => Mass2Td(2)
      Gluons(2)%Mass  => MassTd(2)

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => Top_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ATop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,4
      Gluons(2)%Pol => PolV(:,HelList(iHel,1),2)
      Gluons(2)%Helicity => HelList(iHel,1)
      Quarks(1)%Pol => PolV(:,HelList(iHel,2),3)
      Quarks(1)%Helicity => HelList(iHel,2)
      Quarks(2)%Pol => PolV(:,HelList(iHel,3),4)
      Quarks(2)%Helicity => HelList(iHel,3)

!      calc currents
      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,1,0,0/) )
      Amp(1,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(1,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,0,0,1/) )
      Amp(2,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(2,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      do col =1,2
      do colP=1,2
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus)) *Amp(col,plus)     &
                                                    + dconjg(Amp(colP,minus))*Amp(col,minus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * 2d0*SqAmp * PSWgt

      return
      END FUNCTION






      FUNCTION Tree_GG_TTb_452_5(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modKinematics
      use modMyRecurrence
      implicit none
      double precision Split_A,Split_B,Split_V(0:3),MomDK(0:3,1:6)
      double complex SqAmp,Amp(2,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double precision :: Msq_T_BW,Msq_W_ENU,PSWgt,MomTmp(0:3),PSWgt3,PSWgt2
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      real(8),parameter::ColCorr(1:2,1:2)=(/-64.D0,8.D0,8.D0,8.D0/)
      integer,target :: HelList(1:4,1:3)
      logical :: applyPSCut

      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
      HelList(3,1:3)=(/0,1,0/)
      HelList(4,1:3)=(/1,1,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

      PSWgt = 1d0
IF( TOPDECAYS.LT.0 ) THEN
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)
ENDIF


      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin)   ! top and anti-top are in wrong order for this call, but for pT-distr. it is okay
      if( applyPSCut ) then
          SqAmp = 0d0
          return
      endif
! IF( OBSSET.EQ.8 ) THEN
!       if( Collider.eq.1 ) then
!         r_sc = calc_rgg(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       else
!         r_sc = calc_rqq(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       endif
! ENDIF

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)

!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

!      set particles
      call pol_mless(MomTd(0:3,1),+1,PolV(0:3,plus,1))
      call pol_mless(MomTd(0:3,1),-1,PolV(0:3,minus,1))

      call pol_mless(MomTd(0:3,2),+1,PolV(0:3,plus,2))
      call pol_mless(MomTd(0:3,2),-1,PolV(0:3,minus,2))
      Gluons(2)%Mom => MomTd(:,2)
      Gluons(2)%PartType => Glu_
      Gluons(2)%ExtRef => ExtRef
      Gluons(2)%Mass2 => Mass2Td(2)
      Gluons(2)%Mass  => MassTd(2)

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => Top_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ATop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,4
      Gluons(2)%Pol => PolV(:,HelList(iHel,1),2)
      Gluons(2)%Helicity => HelList(iHel,1)
      Quarks(1)%Pol => PolV(:,HelList(iHel,2),3)
      Quarks(1)%Helicity => HelList(iHel,2)
      Quarks(2)%Pol => PolV(:,HelList(iHel,3),4)
      Quarks(2)%Helicity => HelList(iHel,3)

!      calc currents
      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,1,0,0/) )
      Amp(1,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(1,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,0,0,1/) )
      Amp(2,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(2,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      do col =1,2
      do colP=1,2
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus)) *Amp(col,plus)     &
                                                    + dconjg(Amp(colP,minus))*Amp(col,minus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * 2d0*SqAmp * PSWgt

      return
      END FUNCTION






      FUNCTION Tree_GG_TTb_453_5(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modKinematics
      use modMyRecurrence
      implicit none
      double precision Split_A,Split_B,Split_V(0:3),MomDK(0:3,1:6)
      double complex SqAmp,Amp(2,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double precision :: Msq_T_BW,Msq_W_ENU,PSWgt,MomTmp(0:3),PSWgt3,PSWgt2
      double complex, target :: PolV(0:3,0:1,1:4)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      real(8),parameter::ColCorr(1:2,1:2)=(/-8.D0/9.D0,-80.D0/9.D0,-80.D0/9.D0,-8.D0/9.D0/)
      integer,target :: HelList(1:4,1:3)
      logical :: applyPSCut

      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
      HelList(3,1:3)=(/0,1,0/)
      HelList(4,1:3)=(/1,1,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

      PSWgt = 1d0
IF( TOPDECAYS.LT.0 ) THEN
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
   call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)
ENDIF


      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin)   ! top and anti-top are in wrong order for this call, but for pT-distr. it is okay
      if( applyPSCut ) then
          SqAmp = 0d0
          return
      endif
! IF( OBSSET.EQ.8 ) THEN
!       if( Collider.eq.1 ) then
!         r_sc = calc_rgg(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       else
!         r_sc = calc_rqq(dble(MomTd(1:4,1:4)),MomDK(1:4,1:6))
!       endif
! ENDIF

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)

!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

!      set particles
      call pol_mless(MomTd(0:3,1),+1,PolV(0:3,plus,1))
      call pol_mless(MomTd(0:3,1),-1,PolV(0:3,minus,1))

      call pol_mless(MomTd(0:3,2),+1,PolV(0:3,plus,2))
      call pol_mless(MomTd(0:3,2),-1,PolV(0:3,minus,2))
      Gluons(2)%Mom => MomTd(:,2)
      Gluons(2)%PartType => Glu_
      Gluons(2)%ExtRef => ExtRef
      Gluons(2)%Mass2 => Mass2Td(2)
      Gluons(2)%Mass  => MassTd(2)

      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),+1,PolV(0:3,plus,3))
      call ubarSpi(MomTd(0:3,3),dsqrt(Mass2Td(3)),-1,PolV(0:3,minus,3))
      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => Top_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),+1,PolV(0:3,plus,4))
      call vSpi(MomTd(0:3,4),dsqrt(Mass2Td(4)),-1,PolV(0:3,minus,4))
      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ATop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,4
      Gluons(2)%Pol => PolV(:,HelList(iHel,1),2)
      Gluons(2)%Helicity => HelList(iHel,1)
      Quarks(1)%Pol => PolV(:,HelList(iHel,2),3)
      Quarks(1)%Helicity => HelList(iHel,2)
      Quarks(2)%Pol => PolV(:,HelList(iHel,3),4)
      Quarks(2)%Helicity => HelList(iHel,3)

!      calc currents
      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,1,0,0/) )
      Amp(1,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(1,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      Res(0:3) = cur_g_2f( (/Gluons(2)/),(/Quarks(2),Quarks(1)/),(/2,0,0,1/) )
      Amp(2,plus)  = Res(0:3).dot.PolV(0:3,plus,1)
      Amp(2,minus) = Res(0:3).dot.PolV(0:3,minus,1)

      do col =1,2
      do colP=1,2
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus)) *Amp(col,plus)     &
                                                    + dconjg(Amp(colP,minus))*Amp(col,minus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * 2d0*SqAmp * PSWgt

      return
      END FUNCTION








      END MODULE
