#DEFINE _CHECK_DIPOLE_MOMMAP 0
#DEFINE _APPLY_CUTS 0

      MODULE ModDipoles_QQBTTBG
      use ModTopdecay
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

      double precision, private :: yRndDK(1:8),xFrac,r_sc
      integer, private :: NBin(1:20)

      contains









      SUBROUTINE EvalDipoles_QQBTTBG(Mom,Mass2,yRnd,Wgt,ResultSum)
!     evaluate dipole contributions
      use ModKinematics
      use ModMisc
      use ModParameters
      implicit none
      double precision Mom(0:3,1:5),Mass2(1:5),Wgt
      double precision, parameter :: DipMinus=-1d0, ColorConv=0.5d0
      integer NHisto
      double precision ResultSum,yRnd(1:9)
      double precision RunFactor

      yRndDK(1:8)=yRnd(1:8)
      xFrac = yRnd(9)

!     evaluate all dipoles and do binning:
!     ResultSum has to be added to the real emission amlpitude (minus is incoporated as DipMinus)
!     ColorConv accounts for factor 1/2 from different color convention in dipoles
      ResultSum = (0d0)
      r_sc=0d0
      TheDipoles(1)%alphaCut = alpha_ii
      call evalDipole1_152(Mom,Mass2,TheDipoles(1))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(1)%DipoleValue = DipMinus * ColorConv * Wgt * TheDipoles(1)%DipoleValue * RunFactor**3
IF( ObsSet.EQ.7 ) THEN
      do NHisto=1,3
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(1)%DipoleValue)
      enddo
      call intoHisto(4,1,MInv_LB*TheDipoles(1)%DipoleValue)
      call intoHisto(5,1,MInv_LB**2*TheDipoles(1)%DipoleValue)
      call intoHisto(6,1,xFrac*TheDipoles(1)%DipoleValue)
ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(1)%DipoleValue)
      enddo
ENDIF
      ResultSum = ResultSum + TheDipoles(1)%DipoleValue

      r_sc=0d0
      TheDipoles(2)%alphaCut = alpha_if
      call evalDipole2_153(Mom,Mass2,TheDipoles(2))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(2)%DipoleValue = DipMinus * ColorConv * Wgt * TheDipoles(2)%DipoleValue * RunFactor**3
IF( ObsSet.EQ.7 ) THEN
      do NHisto=1,3
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(2)%DipoleValue)
      enddo
      call intoHisto(4,1,MInv_LB*TheDipoles(2)%DipoleValue)
      call intoHisto(5,1,MInv_LB**2*TheDipoles(2)%DipoleValue)
      call intoHisto(6,1,xFrac*TheDipoles(2)%DipoleValue)
ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(2)%DipoleValue)
      enddo
ENDIF
      ResultSum = ResultSum + TheDipoles(2)%DipoleValue


      r_sc=0d0
      TheDipoles(3)%alphaCut = alpha_if
      call evalDipole3_154(Mom,Mass2,TheDipoles(3))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(3)%DipoleValue = DipMinus * ColorConv * Wgt * TheDipoles(3)%DipoleValue * RunFactor**3
IF( ObsSet.EQ.7 ) THEN
      do NHisto=1,3
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(3)%DipoleValue)
      enddo
      call intoHisto(4,1,MInv_LB*TheDipoles(3)%DipoleValue)
      call intoHisto(5,1,MInv_LB**2*TheDipoles(3)%DipoleValue)
      call intoHisto(6,1,xFrac*TheDipoles(3)%DipoleValue)
ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(3)%DipoleValue)
      enddo
ENDIF
      ResultSum = ResultSum + TheDipoles(3)%DipoleValue


      r_sc=0d0
      TheDipoles(4)%alphaCut = alpha_ii
      call evalDipole4_251(Mom,Mass2,TheDipoles(4))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(4)%DipoleValue = DipMinus * ColorConv * Wgt * TheDipoles(4)%DipoleValue * RunFactor**3
IF( ObsSet.EQ.7 ) THEN
      do NHisto=1,3
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(4)%DipoleValue)
      enddo
      call intoHisto(4,1,MInv_LB*TheDipoles(4)%DipoleValue)
      call intoHisto(5,1,MInv_LB**2*TheDipoles(4)%DipoleValue)
      call intoHisto(6,1,xFrac*TheDipoles(4)%DipoleValue)
ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(4)%DipoleValue)
      enddo
ENDIF
      ResultSum = ResultSum + TheDipoles(4)%DipoleValue


      r_sc=0d0
      TheDipoles(5)%alphaCut = alpha_if
      call evalDipole5_253(Mom,Mass2,TheDipoles(5))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(5)%DipoleValue = DipMinus * ColorConv * Wgt * TheDipoles(5)%DipoleValue * RunFactor**3
IF( ObsSet.EQ.7 ) THEN
      do NHisto=1,3
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(5)%DipoleValue)
      enddo
      call intoHisto(4,1,MInv_LB*TheDipoles(5)%DipoleValue)
      call intoHisto(5,1,MInv_LB**2*TheDipoles(5)%DipoleValue)
      call intoHisto(6,1,xFrac*TheDipoles(5)%DipoleValue)
ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(5)%DipoleValue)
      enddo
ENDIF
      ResultSum = ResultSum + TheDipoles(5)%DipoleValue


      r_sc=0d0
      TheDipoles(6)%alphaCut = alpha_if
      call evalDipole6_254(Mom,Mass2,TheDipoles(6))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(6)%DipoleValue = DipMinus * ColorConv * Wgt * TheDipoles(6)%DipoleValue * RunFactor**3
IF( ObsSet.EQ.7 ) THEN
      do NHisto=1,3
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(6)%DipoleValue)
      enddo
      call intoHisto(4,1,MInv_LB*TheDipoles(6)%DipoleValue)
      call intoHisto(5,1,MInv_LB**2*TheDipoles(6)%DipoleValue)
      call intoHisto(6,1,xFrac*TheDipoles(6)%DipoleValue)
ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(6)%DipoleValue)
      enddo
ENDIF
      ResultSum = ResultSum + TheDipoles(6)%DipoleValue


      r_sc=0d0
      TheDipoles(7)%alphaCut = alpha_fi
      call evalDipole7_351(Mom,Mass2,TheDipoles(7))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(7)%DipoleValue = DipMinus * ColorConv * Wgt * TheDipoles(7)%DipoleValue * RunFactor**3
IF( ObsSet.EQ.7 ) THEN
      do NHisto=1,3
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(7)%DipoleValue)
      enddo
      call intoHisto(4,1,MInv_LB*TheDipoles(7)%DipoleValue)
      call intoHisto(5,1,MInv_LB**2*TheDipoles(7)%DipoleValue)
      call intoHisto(6,1,xFrac*TheDipoles(7)%DipoleValue)
ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(7)%DipoleValue)
      enddo
ENDIF
      ResultSum = ResultSum + TheDipoles(7)%DipoleValue


      r_sc=0d0
      TheDipoles(8)%alphaCut = alpha_fi
      call evalDipole8_352(Mom,Mass2,TheDipoles(8))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(8)%DipoleValue = DipMinus * ColorConv * Wgt * TheDipoles(8)%DipoleValue * RunFactor**3
IF( ObsSet.EQ.7 ) THEN
      do NHisto=1,3
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(8)%DipoleValue)
      enddo
      call intoHisto(4,1,MInv_LB*TheDipoles(8)%DipoleValue)
      call intoHisto(5,1,MInv_LB**2*TheDipoles(8)%DipoleValue)
      call intoHisto(6,1,xFrac*TheDipoles(8)%DipoleValue)
ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(8)%DipoleValue)
      enddo
ENDIF
      ResultSum = ResultSum + TheDipoles(8)%DipoleValue

      r_sc=0d0
      TheDipoles(9)%alphaCut = 1d0 !alpha_ff not yet implemeted
      call evalDipole9_354(Mom,Mass2,TheDipoles(9))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(9)%DipoleValue = DipMinus * ColorConv * Wgt * TheDipoles(9)%DipoleValue * RunFactor**3
IF( ObsSet.EQ.7 ) THEN
      do NHisto=1,3
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(9)%DipoleValue)
      enddo
      call intoHisto(4,1,MInv_LB*TheDipoles(9)%DipoleValue)
      call intoHisto(5,1,MInv_LB**2*TheDipoles(9)%DipoleValue)
      call intoHisto(6,1,xFrac*TheDipoles(9)%DipoleValue)
ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(9)%DipoleValue)
      enddo
ENDIF
      ResultSum = ResultSum + TheDipoles(9)%DipoleValue


      r_sc=0d0
      TheDipoles(10)%alphaCut = alpha_fi
      call evalDipole10_451(Mom,Mass2,TheDipoles(10))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(10)%DipoleValue = DipMinus * ColorConv * Wgt * TheDipoles(10)%DipoleValue * RunFactor**3
IF( ObsSet.EQ.7 ) THEN
      do NHisto=1,3
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(10)%DipoleValue)
      enddo
      call intoHisto(4,1,MInv_LB*TheDipoles(10)%DipoleValue)
      call intoHisto(5,1,MInv_LB**2*TheDipoles(10)%DipoleValue)
      call intoHisto(6,1,xFrac*TheDipoles(10)%DipoleValue)
ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(10)%DipoleValue)
      enddo
ENDIF
      ResultSum = ResultSum + TheDipoles(10)%DipoleValue


      r_sc=0d0
      TheDipoles(11)%alphaCut = alpha_fi
      call evalDipole11_452(Mom,Mass2,TheDipoles(11))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(11)%DipoleValue = DipMinus * ColorConv * Wgt * TheDipoles(11)%DipoleValue * RunFactor**3
IF( ObsSet.EQ.7 ) THEN
      do NHisto=1,3
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(11)%DipoleValue)
      enddo
      call intoHisto(4,1,MInv_LB*TheDipoles(11)%DipoleValue)
      call intoHisto(5,1,MInv_LB**2*TheDipoles(11)%DipoleValue)
      call intoHisto(6,1,xFrac*TheDipoles(11)%DipoleValue)
ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(11)%DipoleValue)
      enddo
ENDIF
      ResultSum = ResultSum + TheDipoles(11)%DipoleValue


      r_sc=0d0
      TheDipoles(12)%alphaCut = 1d0 !alpha_ff not yet implemented
      call evalDipole12_453(Mom,Mass2,TheDipoles(12))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(12)%DipoleValue = DipMinus * ColorConv * Wgt * TheDipoles(12)%DipoleValue * RunFactor**3
IF( ObsSet.EQ.7 ) THEN
      do NHisto=1,3
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(12)%DipoleValue)
      enddo
      call intoHisto(4,1,MInv_LB*TheDipoles(12)%DipoleValue)
      call intoHisto(5,1,MInv_LB**2*TheDipoles(12)%DipoleValue)
      call intoHisto(6,1,xFrac*TheDipoles(12)%DipoleValue)
ELSE
      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),TheDipoles(12)%DipoleValue)
      enddo
ENDIF
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

!     splitting 1: A/Quark_I --> A/Quark_I + Gluon_F
      Split_A = 2d0*alpha_s4Pi*(1d0+x**2)/(1d0-x)

!     color and spin correlated tree process
!      SqAmp = Tree_UUb_TTb_152_1(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_UUb_TTb_ij(1,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)



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


!     splitting 1: A/Quark_I --> A/Quark_I + Gluon_F
      Split_A = 2d0*alpha_s4Pi*(2d0/(1d0-x+u)-1d0-x)

!     color and spin correlated tree process
!      SqAmp = Tree_UUb_TTb_153_1(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_UUb_TTb_ij(2,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

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

!     splitting 1: A/Quark_I --> A/Quark_I + Gluon_F
      Split_A = 2d0*alpha_s4Pi*(2d0/(1d0-x+u)-1d0-x)

!     color and spin correlated tree process
!      SqAmp = Tree_UUb_TTb_154_1(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_UUb_TTb_ij(3,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)


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

!     splitting 1: A/Quark_I --> A/Quark_I + Gluon_F
      Split_A = 2d0*alpha_s4Pi*(1d0+x**2)/(1d0-x)

!     color and spin correlated tree process
!      SqAmp = Tree_UUb_TTb_251_1(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_UUb_TTb_ij(4,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)


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


!     splitting 1: A/Quark_I --> A/Quark_I + Gluon_F
      Split_A = 2d0*alpha_s4Pi*(2d0/(1d0-x+u)-1d0-x)

!     color and spin correlated tree process
!      SqAmp = Tree_UUb_TTb_253_1(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_UUb_TTb_ij(5,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)



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


!     splitting 1: A/Quark_I --> A/Quark_I + Gluon_F
      Split_A = 2d0*alpha_s4Pi*(2d0/(1d0-x+u)-1d0-x)

!     color and spin correlated tree process
!      SqAmp = Tree_UUb_TTb_254_1(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_UUb_TTb_ij(6,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)



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
      if( TheDipole%AlphaCut.lt.1d0-x ) then  ! NEW
         TheDipole%DipoleValue = 0d0
         return
      endif
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

!     splitting 5: A/Quark_F --> A/Quark_F + Gluon_F
      Split_A = 2d0*alpha_s4Pi*(2d0/(2d0-z-x)-1d0-z-2d0*Mass2(3)/s35)

!     color and spin correlated tree process
!      SqAmp = Tree_UUb_TTb_351_5(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_UUb_TTb_ij(7,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)


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
      if( TheDipole%AlphaCut.lt.1d0-x ) then  ! NEW
         TheDipole%DipoleValue = 0d0
         return
      endif
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

!     splitting 5: A/Quark_F --> A/Quark_F + Gluon_F
      Split_A = 2d0*alpha_s4Pi*(2d0/(2d0-z-x)-1d0-z-2d0*Mass2(3)/s35)

!     color and spin correlated tree process
!      SqAmp = Tree_UUb_TTb_352_5(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_UUb_TTb_ij(8,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)


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


!     splitting 5: A/Quark_F --> A/Quark_F + Gluon_F
      mu2_i = Mass2(3)/dabs(Q2)
      mu2_k = Mass2(4)/dabs(Q2)
      v = (1d0-y)*dsqrt((1d0-(mu2_i-mu2_k)**2-2d0*(mu2_i+mu2_k))/( (2d0*mu2_k+(1d0-mu2_i-mu2_k)*(1d0-y))**2-4d0*mu2_k))
      Split_A = 2d0*alpha_s4Pi*(2d0/(1d0-z+y*z)-v*(1d0+z+2d0*Mass2(3)/s35))

!     color and spin correlated tree process
!      SqAmp = Tree_UUb_TTb_354_5(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_UUb_TTb_ij(9,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

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
      if( TheDipole%AlphaCut.lt.1d0-x ) then  ! NEW
         TheDipole%DipoleValue = 0d0
         return
      endif
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

!     splitting 5: A/Quark_F --> A/Quark_F + Gluon_F
      Split_A = 2d0*alpha_s4Pi*(2d0/(2d0-z-x)-1d0-z-2d0*Mass2(4)/s45)

!     color and spin correlated tree process
!      SqAmp = Tree_UUb_TTb_451_5(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_UUb_TTb_ij(10,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)


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
      if( TheDipole%AlphaCut.lt.1d0-x ) then  ! NEW
         TheDipole%DipoleValue = 0d0
         return
      endif
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

!     splitting 5: A/Quark_F --> A/Quark_F + Gluon_F
      Split_A = 2d0*alpha_s4Pi*(2d0/(2d0-z-x)-1d0-z-2d0*Mass2(4)/s45)

!     color and spin correlated tree process
!      SqAmp = Tree_UUb_TTb_452_5(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_UUb_TTb_ij(11,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)


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


!     splitting 5: A/Quark_F --> A/Quark_F + Gluon_F
      mu2_i = Mass2(4)/dabs(Q2)
      mu2_k = Mass2(3)/dabs(Q2)
      v = (1d0-y)*dsqrt((1d0-(mu2_i-mu2_k)**2-2d0*(mu2_i+mu2_k))/( (2d0*mu2_k+(1d0-mu2_i-mu2_k)*(1d0-y))**2-4d0*mu2_k))
      Split_A = 2d0*alpha_s4Pi*(2d0/(1d0-z+y*z)-v*(1d0+z+2d0*Mass2(4)/s45))

!     color and spin correlated tree process
!      SqAmp = Tree_UUb_TTb_453_5(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_UUb_TTb_ij(12,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)


      TheDipole%DipoleValue = LeadSing * SqAmp
      return
      END SUBROUTINE





      FUNCTION Tree_UUb_TTb_ij(icorr,MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      use modKinematics
      implicit none
      double precision Split_A,Split_B,Split_V(0:3)
      double complex SqAmp,Amp(1,0:1),Res(0:3)
      integer col,colP,iHel,icorr,ngl
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      double precision :: PSWgt2,PSWgt3,MomDK(0:3,1:6)
      type(Particle),target :: TopQuark(1:2)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:1),Quarks(1:4)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(9)
      real(8)::ColCorr(1:1,1:1)
      real(8),parameter::ColCorr1(1:1,1:1)=(/8.D0/3.D0/)
      real(8),parameter::ColCorr2(1:1,1:1)=(/-56.D0/3.D0/)
      real(8),parameter::ColCorr3(1:1,1:1)=(/-16.D0/3.D0/)
      real(8),parameter::ColCorr4(1:1,1:1)=(/8.D0/3.D0/)
      real(8),parameter::ColCorr5(1:1,1:1)=(/-16.D0/3.D0/)
      real(8),parameter::ColCorr6(1:1,1:1)=(/-56.D0/3.D0/)
      real(8),parameter::ColCorr7(1:1,1:1)=(/-56.D0/3.D0/)
      real(8),parameter::ColCorr8(1:1,1:1)=(/-16.D0/3.D0/)
      real(8),parameter::ColCorr9(1:1,1:1)=(/8.D0/3.D0/)
      real(8),parameter::ColCorr10(1:1,1:1)=(/-16.D0/3.D0/)
      real(8),parameter::ColCorr11(1:1,1:1)=(/-56.D0/3.D0/)
      real(8),parameter::ColCorr12(1:1,1:1)=(/8.D0/3.D0/)
      integer,target :: HelList(1:16,1:4)
      double complex :: MomTmp(0:3)
      logical :: applyPSCut


      if(icorr.eq.1) ColCorr(1:1,1:1)=ColCorr1(1:1,1:1)
      if(icorr.eq.2) ColCorr(1:1,1:1)=ColCorr2(1:1,1:1)
      if(icorr.eq.3) ColCorr(1:1,1:1)=ColCorr3(1:1,1:1)
      if(icorr.eq.4) ColCorr(1:1,1:1)=ColCorr4(1:1,1:1)
      if(icorr.eq.5) ColCorr(1:1,1:1)=ColCorr5(1:1,1:1)
      if(icorr.eq.6) ColCorr(1:1,1:1)=ColCorr6(1:1,1:1)
      if(icorr.eq.7) ColCorr(1:1,1:1)=ColCorr7(1:1,1:1)
      if(icorr.eq.8) ColCorr(1:1,1:1)=ColCorr8(1:1,1:1)
      if(icorr.eq.9) ColCorr(1:1,1:1)=ColCorr9(1:1,1:1)
      if(icorr.eq.10) ColCorr(1:1,1:1)=ColCorr10(1:1,1:1)
      if(icorr.eq.11) ColCorr(1:1,1:1)=ColCorr11(1:1,1:1)
      if(icorr.eq.12) ColCorr(1:1,1:1)=ColCorr12(1:1,1:1)


      HelList(1,1:4)=(/0,0,0,0/)
      HelList(2,1:4)=(/0,1,0,0/)


      MassTd(1:4) = dsqrt(Mass2Td(1:4))

      call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
      call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top momenta to be consistent with Kinematics conventions
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin,xJPsiFrag=xFrac)
      if( applyPSCut ) then
          SqAmp = 0d0
          return
      endif


      MomTmp(0:3)  = MomTd(0:3,3)! swap back
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)

!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

      TopQuark(1)%PartType = Top_
      TopQuark(1)%Mom(1:4) = MomTd(0:3,3)
      call TopDecay(TopQuark(1),DK_LO,MomDK(0:3,4:6))
      TopQuark(2)%PartType = ATop_
      TopQuark(2)%Mom(1:4) = MomTd(0:3,4)
      call TopDecay(TopQuark(2),DK_LO,MomDK(0:3,1:3))

      Quarks(3)%Pol => TopQuark(1)%Pol
      Quarks(4)%Pol => TopQuark(2)%Pol



!      set particles
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),+1,PolV(0:3,plus,1))
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),-1,PolV(0:3,minus,1))
      Quarks(1)%Mom => MomTd(:,1)
      Quarks(1)%PartType => AUp_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(1)
      Quarks(1)%Mass  => MassTd(1)

      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),+1,PolV(0:3,plus,2))
      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),-1,PolV(0:3,minus,2))
      Quarks(2)%Mom => MomTd(:,2)
      Quarks(2)%PartType => Up_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(2)
      Quarks(2)%Mass  => MassTd(2)

      Quarks(3)%Mom => MomTd(:,3)
      Quarks(3)%PartType => Top_
      Quarks(3)%ExtRef => ExtRef
      Quarks(3)%Mass2 => Mass2Td(3)
      Quarks(3)%Mass  => MassTd(3)

      Quarks(4)%Mom => MomTd(:,4)
      Quarks(4)%PartType => ATop_
      Quarks(4)%ExtRef => ExtRef
      Quarks(4)%Mass2 => Mass2Td(4)
      Quarks(4)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,2
      Quarks(2)%Pol => PolV(:,HelList(iHel,2),2)
      Quarks(2)%Helicity => HelList(iHel,2)

!      calc currents
      Res(0:3) = cur_f_4f((/Gluons(1)/),(/Quarks(2),Quarks(3),Quarks(4)/),Quarks(1)%PartType,(/0,0,0,0,0/),0)
      Amp(1,plus)  = psp1_(Res(0:3),PolV(0:3,plus,1))
      Amp(1,minus) = psp1_(Res(0:3),PolV(0:3,minus,1))

      do col =1,1
      do colP=1,1
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus))*Amp(col,plus) + dconjg(Amp(colP,minus))*Amp(col,minus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

      return
      END FUNCTION





      FUNCTION Tree_UUb_TTb_152_1(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      use modKinematics
      implicit none
      double precision Split_A,Split_B,Split_V(0:3)
      double complex SqAmp,Amp(1,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      double precision :: PSWgt2,PSWgt3,MomDK(0:3,1:6)
      type(Particle),target :: TopQuark(1:2)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:1),Quarks(1:4)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(9)
      real(8),parameter::ColCorr(1:1,1:1)=(/8.D0/3.D0/)
      integer,target :: HelList(1:16,1:4)
      double complex :: MomTmp(0:3)
      logical :: applyPSCut

      HelList(1,1:4)=(/0,0,0,0/)
      HelList(2,1:4)=(/0,1,0,0/)
!       HelList(1,1:4)=(/0,0,0,0/)
!       HelList(2,1:4)=(/1,0,0,0/)
!       HelList(3,1:4)=(/0,1,0,0/)
!       HelList(4,1:4)=(/1,1,0,0/)
!       HelList(5,1:4)=(/0,0,1,0/)
!       HelList(6,1:4)=(/1,0,1,0/)
!       HelList(7,1:4)=(/0,1,1,0/)
!       HelList(8,1:4)=(/1,1,1,0/)
!       HelList(9,1:4)=(/0,0,0,1/)
!       HelList(10,1:4)=(/1,0,0,1/)
!       HelList(11,1:4)=(/0,1,0,1/)
!       HelList(12,1:4)=(/1,1,0,1/)
!       HelList(13,1:4)=(/0,0,1,1/)
!       HelList(14,1:4)=(/1,0,1,1/)
!       HelList(15,1:4)=(/0,1,1,1/)
!       HelList(16,1:4)=(/1,1,1,1/)

      MassTd(1:4) = dsqrt(Mass2Td(1:4))


      call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
      call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top momenta to be consistent with Kinematics conventions
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin,xJPsiFrag=xFrac)
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

      TopQuark(1)%PartType = Top_
      TopQuark(1)%Mom(1:4) = MomTd(0:3,3)
      call TopDecay(TopQuark(1),DK_LO,MomDK(0:3,4:6))
      TopQuark(2)%PartType = ATop_
      TopQuark(2)%Mom(1:4) = MomTd(0:3,4)
      call TopDecay(TopQuark(2),DK_LO,MomDK(0:3,1:3))

      Quarks(3)%Pol => TopQuark(1)%Pol
      Quarks(4)%Pol => TopQuark(2)%Pol



!      set particles
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),+1,PolV(0:3,plus,1))
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),-1,PolV(0:3,minus,1))
      Quarks(1)%Mom => MomTd(:,1)
      Quarks(1)%PartType => AUp_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(1)
      Quarks(1)%Mass  => MassTd(1)

      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),+1,PolV(0:3,plus,2))
      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),-1,PolV(0:3,minus,2))
      Quarks(2)%Mom => MomTd(:,2)
      Quarks(2)%PartType => Up_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(2)
      Quarks(2)%Mass  => MassTd(2)

      Quarks(3)%Mom => MomTd(:,3)
      Quarks(3)%PartType => Top_
      Quarks(3)%ExtRef => ExtRef
      Quarks(3)%Mass2 => Mass2Td(3)
      Quarks(3)%Mass  => MassTd(3)

      Quarks(4)%Mom => MomTd(:,4)
      Quarks(4)%PartType => ATop_
      Quarks(4)%ExtRef => ExtRef
      Quarks(4)%Mass2 => Mass2Td(4)
      Quarks(4)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
!       do iHel=1,16
      do iHel=1,2
!       Quarks(1)%Pol => PolV(:,HelList(iHel,1),1)
!       Quarks(1)%Helicity => HelList(iHel,1)
      Quarks(2)%Pol => PolV(:,HelList(iHel,2),2)
      Quarks(2)%Helicity => HelList(iHel,2)
!       Quarks(3)%Pol => PolV(:,HelList(iHel,3),3)
!       Quarks(3)%Helicity => HelList(iHel,3)
!       Quarks(4)%Pol => PolV(:,HelList(iHel,4),4)
!       Quarks(4)%Helicity => HelList(iHel,4)

!      calc currents
      Res(0:3) = cur_f_4f((/Gluons(1)/),(/Quarks(2),Quarks(3),Quarks(4)/),Quarks(1)%PartType,(/0,0,0,0,0/),0)
      Amp(1,plus)  = psp1_(Res(0:3),PolV(0:3,plus,1))
      Amp(1,minus) = psp1_(Res(0:3),PolV(0:3,minus,1))

      do col =1,1
      do colP=1,1
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus))*Amp(col,plus) + dconjg(Amp(colP,minus))*Amp(col,minus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

      return
      END FUNCTION






      FUNCTION Tree_UUb_TTb_153_1(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      use modKinematics
      implicit none
      double precision Split_A,Split_B,Split_V(0:3)
      double complex SqAmp,Amp(1,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      double precision :: PSWgt2,PSWgt3,MomDK(0:3,1:6)
      type(Particle),target :: TopQuark(1:2)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:1),Quarks(1:4)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(9)
      real(8),parameter::ColCorr(1:1,1:1)=(/-56.D0/3.D0/)
      integer,target :: HelList(1:16,1:4)
      double complex :: MomTmp(0:3)
      logical :: applyPSCut

      HelList(1,1:4)=(/0,0,0,0/)
      HelList(2,1:4)=(/0,1,0,0/)
!       HelList(3,1:4)=(/0,1,0,0/)
!       HelList(4,1:4)=(/1,1,0,0/)
!       HelList(5,1:4)=(/0,0,1,0/)
!       HelList(6,1:4)=(/1,0,1,0/)
!       HelList(7,1:4)=(/0,1,1,0/)
!       HelList(8,1:4)=(/1,1,1,0/)
!       HelList(9,1:4)=(/0,0,0,1/)
!       HelList(10,1:4)=(/1,0,0,1/)
!       HelList(11,1:4)=(/0,1,0,1/)
!       HelList(12,1:4)=(/1,1,0,1/)
!       HelList(13,1:4)=(/0,0,1,1/)
!       HelList(14,1:4)=(/1,0,1,1/)
!       HelList(15,1:4)=(/0,1,1,1/)
!       HelList(16,1:4)=(/1,1,1,1/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))


      call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
      call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top momenta to be consistent with Kinematics conventions
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin,xJPsiFrag=xFrac)
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



      TopQuark(1)%PartType = Top_
      TopQuark(1)%Mom(1:4) = MomTd(0:3,3)
      call TopDecay(TopQuark(1),DK_LO,MomDK(0:3,4:6))
      TopQuark(2)%PartType = ATop_
      TopQuark(2)%Mom(1:4) = MomTd(0:3,4)
      call TopDecay(TopQuark(2),DK_LO,MomDK(0:3,1:3))

      Quarks(3)%Pol => TopQuark(1)%Pol
      Quarks(4)%Pol => TopQuark(2)%Pol


!      set particles
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),+1,PolV(0:3,plus,1))
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),-1,PolV(0:3,minus,1))
      Quarks(1)%Mom => MomTd(:,1)
      Quarks(1)%PartType => AUp_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(1)
      Quarks(1)%Mass  => MassTd(1)

      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),+1,PolV(0:3,plus,2))
      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),-1,PolV(0:3,minus,2))
      Quarks(2)%Mom => MomTd(:,2)
      Quarks(2)%PartType => Up_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(2)
      Quarks(2)%Mass  => MassTd(2)

      Quarks(3)%Mom => MomTd(:,3)
      Quarks(3)%PartType => Top_
      Quarks(3)%ExtRef => ExtRef
      Quarks(3)%Mass2 => Mass2Td(3)
      Quarks(3)%Mass  => MassTd(3)

      Quarks(4)%Mom => MomTd(:,4)
      Quarks(4)%PartType => ATop_
      Quarks(4)%ExtRef => ExtRef
      Quarks(4)%Mass2 => Mass2Td(4)
      Quarks(4)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
!       do iHel=1,16
      do iHel=1,2
!       Quarks(1)%Pol => PolV(:,HelList(iHel,1),1)
!       Quarks(1)%Helicity => HelList(iHel,1)
      Quarks(2)%Pol => PolV(:,HelList(iHel,2),2)
      Quarks(2)%Helicity => HelList(iHel,2)
!       Quarks(3)%Pol => PolV(:,HelList(iHel,3),3)
!       Quarks(3)%Helicity => HelList(iHel,3)
!       Quarks(4)%Pol => PolV(:,HelList(iHel,4),4)
!       Quarks(4)%Helicity => HelList(iHel,4)

!      calc currents
      Res(0:3) = cur_f_4f((/Gluons(1)/),(/Quarks(2),Quarks(3),Quarks(4)/),Quarks(1)%PartType,(/0,0,0,0,0/),0)
      Amp(1,plus)  = psp1_(Res(0:3),PolV(0:3,plus,1))
      Amp(1,minus) = psp1_(Res(0:3),PolV(0:3,minus,1))

      do col =1,1
      do colP=1,1
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus))*Amp(col,plus) + dconjg(Amp(colP,minus))*Amp(col,minus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

      return
      END FUNCTION






      FUNCTION Tree_UUb_TTb_154_1(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      use modKinematics
      implicit none
      double precision Split_A,Split_B,Split_V(0:3)
      double complex SqAmp,Amp(1,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      double precision :: PSWgt2,PSWgt3,MomDK(0:3,1:6)
      type(Particle),target :: TopQuark(1:2)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:1),Quarks(1:4)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(9)
      real(8),parameter::ColCorr(1:1,1:1)=(/-16.D0/3.D0/)
      integer,target :: HelList(1:16,1:4)
      double complex :: MomTmp(0:3)
      logical :: applyPSCut

      HelList(1,1:4)=(/0,0,0,0/)
      HelList(2,1:4)=(/0,1,0,0/)
!       HelList(3,1:4)=(/0,1,0,0/)
!       HelList(4,1:4)=(/1,1,0,0/)
!       HelList(5,1:4)=(/0,0,1,0/)
!       HelList(6,1:4)=(/1,0,1,0/)
!       HelList(7,1:4)=(/0,1,1,0/)
!       HelList(8,1:4)=(/1,1,1,0/)
!       HelList(9,1:4)=(/0,0,0,1/)
!       HelList(10,1:4)=(/1,0,0,1/)
!       HelList(11,1:4)=(/0,1,0,1/)
!       HelList(12,1:4)=(/1,1,0,1/)
!       HelList(13,1:4)=(/0,0,1,1/)
!       HelList(14,1:4)=(/1,0,1,1/)
!       HelList(15,1:4)=(/0,1,1,1/)
!       HelList(16,1:4)=(/1,1,1,1/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))


      call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
      call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top momenta to be consistent with Kinematics conventions
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin,xJPsiFrag=xFrac)
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


      TopQuark(1)%PartType = Top_
      TopQuark(1)%Mom(1:4) = MomTd(0:3,3)
      call TopDecay(TopQuark(1),DK_LO,MomDK(0:3,4:6))
      TopQuark(2)%PartType = ATop_
      TopQuark(2)%Mom(1:4) = MomTd(0:3,4)
      call TopDecay(TopQuark(2),DK_LO,MomDK(0:3,1:3))

      Quarks(3)%Pol => TopQuark(1)%Pol
      Quarks(4)%Pol => TopQuark(2)%Pol


!      set particles
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),+1,PolV(0:3,plus,1))
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),-1,PolV(0:3,minus,1))
      Quarks(1)%Mom => MomTd(:,1)
      Quarks(1)%PartType => AUp_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(1)
      Quarks(1)%Mass  => MassTd(1)

      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),+1,PolV(0:3,plus,2))
      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),-1,PolV(0:3,minus,2))
      Quarks(2)%Mom => MomTd(:,2)
      Quarks(2)%PartType => Up_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(2)
      Quarks(2)%Mass  => MassTd(2)

      Quarks(3)%Mom => MomTd(:,3)
      Quarks(3)%PartType => Top_
      Quarks(3)%ExtRef => ExtRef
      Quarks(3)%Mass2 => Mass2Td(3)
      Quarks(3)%Mass  => MassTd(3)

      Quarks(4)%Mom => MomTd(:,4)
      Quarks(4)%PartType => ATop_
      Quarks(4)%ExtRef => ExtRef
      Quarks(4)%Mass2 => Mass2Td(4)
      Quarks(4)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
!       do iHel=1,16
      do iHel=1,2
!       Quarks(1)%Pol => PolV(:,HelList(iHel,1),1)
!       Quarks(1)%Helicity => HelList(iHel,1)
      Quarks(2)%Pol => PolV(:,HelList(iHel,2),2)
      Quarks(2)%Helicity => HelList(iHel,2)
!       Quarks(3)%Pol => PolV(:,HelList(iHel,3),3)
!       Quarks(3)%Helicity => HelList(iHel,3)
!       Quarks(4)%Pol => PolV(:,HelList(iHel,4),4)
!       Quarks(4)%Helicity => HelList(iHel,4)

!      calc currents
      Res(0:3) = cur_f_4f((/Gluons(1)/),(/Quarks(2),Quarks(3),Quarks(4)/),Quarks(1)%PartType,(/0,0,0,0,0/),0)
      Amp(1,plus)  = psp1_(Res(0:3),PolV(0:3,plus,1))
      Amp(1,minus) = psp1_(Res(0:3),PolV(0:3,minus,1))

      do col =1,1
      do colP=1,1
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus))*Amp(col,plus) + dconjg(Amp(colP,minus))*Amp(col,minus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

      return
      END FUNCTION






      FUNCTION Tree_UUb_TTb_251_1(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      use modKinematics
      implicit none
      double precision Split_A,Split_B,Split_V(0:3)
      double complex SqAmp,Amp(1,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      double precision :: PSWgt2,PSWgt3,MomDK(0:3,1:6)
      type(Particle),target :: TopQuark(1:2)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:1),Quarks(1:4)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(9)
      real(8),parameter::ColCorr(1:1,1:1)=(/8.D0/3.D0/)
      integer,target :: HelList(1:16,1:4)
      double complex :: MomTmp(0:3)
      logical :: applyPSCut

      HelList(1,1:4)=(/0,0,0,0/)
      HelList(2,1:4)=(/0,1,0,0/)
!       HelList(3,1:4)=(/0,1,0,0/)
!       HelList(4,1:4)=(/1,1,0,0/)
!       HelList(5,1:4)=(/0,0,1,0/)
!       HelList(6,1:4)=(/1,0,1,0/)
!       HelList(7,1:4)=(/0,1,1,0/)
!       HelList(8,1:4)=(/1,1,1,0/)
!       HelList(9,1:4)=(/0,0,0,1/)
!       HelList(10,1:4)=(/1,0,0,1/)
!       HelList(11,1:4)=(/0,1,0,1/)
!       HelList(12,1:4)=(/1,1,0,1/)
!       HelList(13,1:4)=(/0,0,1,1/)
!       HelList(14,1:4)=(/1,0,1,1/)
!       HelList(15,1:4)=(/0,1,1,1/)
!       HelList(16,1:4)=(/1,1,1,1/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))


      call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
      call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top momenta to be consistent with Kinematics conventions
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin,xJPsiFrag=xFrac)
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


      TopQuark(1)%PartType = Top_
      TopQuark(1)%Mom(1:4) = MomTd(0:3,3)
      call TopDecay(TopQuark(1),DK_LO,MomDK(0:3,4:6))
      TopQuark(2)%PartType = ATop_
      TopQuark(2)%Mom(1:4) = MomTd(0:3,4)
      call TopDecay(TopQuark(2),DK_LO,MomDK(0:3,1:3))

      Quarks(3)%Pol => TopQuark(1)%Pol
      Quarks(4)%Pol => TopQuark(2)%Pol


!      set particles
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),+1,PolV(0:3,plus,1))
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),-1,PolV(0:3,minus,1))
      Quarks(1)%Mom => MomTd(:,1)
      Quarks(1)%PartType => AUp_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(1)
      Quarks(1)%Mass  => MassTd(1)

      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),+1,PolV(0:3,plus,2))
      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),-1,PolV(0:3,minus,2))
      Quarks(2)%Mom => MomTd(:,2)
      Quarks(2)%PartType => Up_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(2)
      Quarks(2)%Mass  => MassTd(2)

      Quarks(3)%Mom => MomTd(:,3)
      Quarks(3)%PartType => Top_
      Quarks(3)%ExtRef => ExtRef
      Quarks(3)%Mass2 => Mass2Td(3)
      Quarks(3)%Mass  => MassTd(3)

      Quarks(4)%Mom => MomTd(:,4)
      Quarks(4)%PartType => ATop_
      Quarks(4)%ExtRef => ExtRef
      Quarks(4)%Mass2 => Mass2Td(4)
      Quarks(4)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,2
!       Quarks(1)%Pol => PolV(:,HelList(iHel,1),1)
!       Quarks(1)%Helicity => HelList(iHel,1)
      Quarks(2)%Pol => PolV(:,HelList(iHel,2),2)
      Quarks(2)%Helicity => HelList(iHel,2)
!       Quarks(3)%Pol => PolV(:,HelList(iHel,3),3)
!       Quarks(3)%Helicity => HelList(iHel,3)
!       Quarks(4)%Pol => PolV(:,HelList(iHel,4),4)
!       Quarks(4)%Helicity => HelList(iHel,4)

!      calc currents
      Res(0:3) = cur_f_4f((/Gluons(1)/),(/Quarks(2),Quarks(3),Quarks(4)/),Quarks(1)%PartType,(/0,0,0,0,0/),0)
      Amp(1,plus)  = psp1_(Res(0:3),PolV(0:3,plus,1))
      Amp(1,minus) = psp1_(Res(0:3),PolV(0:3,minus,1))

      do col =1,1
      do colP=1,1
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus))*Amp(col,plus) + dconjg(Amp(colP,minus))*Amp(col,minus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

      return
      END FUNCTION






      FUNCTION Tree_UUb_TTb_253_1(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      use modKinematics
      implicit none
      double precision Split_A,Split_B,Split_V(0:3)
      double complex SqAmp,Amp(1,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      double precision :: PSWgt2,PSWgt3,MomDK(0:3,1:6)
      type(Particle),target :: TopQuark(1:2)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:1),Quarks(1:4)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(9)
      real(8),parameter::ColCorr(1:1,1:1)=(/-16.D0/3.D0/)
      integer,target :: HelList(1:16,1:4)
      double complex :: MomTmp(0:3)
      logical :: applyPSCut

      HelList(1,1:4)=(/0,0,0,0/)
      HelList(2,1:4)=(/0,1,0,0/)
!       HelList(3,1:4)=(/0,1,0,0/)
!       HelList(4,1:4)=(/1,1,0,0/)
!       HelList(5,1:4)=(/0,0,1,0/)
!       HelList(6,1:4)=(/1,0,1,0/)
!       HelList(7,1:4)=(/0,1,1,0/)
!       HelList(8,1:4)=(/1,1,1,0/)
!       HelList(9,1:4)=(/0,0,0,1/)
!       HelList(10,1:4)=(/1,0,0,1/)
!       HelList(11,1:4)=(/0,1,0,1/)
!       HelList(12,1:4)=(/1,1,0,1/)
!       HelList(13,1:4)=(/0,0,1,1/)
!       HelList(14,1:4)=(/1,0,1,1/)
!       HelList(15,1:4)=(/0,1,1,1/)
!       HelList(16,1:4)=(/1,1,1,1/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))


      call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
      call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top momenta to be consistent with Kinematics conventions
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin,xJPsiFrag=xFrac)
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


      TopQuark(1)%PartType = Top_
      TopQuark(1)%Mom(1:4) = MomTd(0:3,3)
      call TopDecay(TopQuark(1),DK_LO,MomDK(0:3,4:6))
      TopQuark(2)%PartType = ATop_
      TopQuark(2)%Mom(1:4) = MomTd(0:3,4)
      call TopDecay(TopQuark(2),DK_LO,MomDK(0:3,1:3))

      Quarks(3)%Pol => TopQuark(1)%Pol
      Quarks(4)%Pol => TopQuark(2)%Pol


!      set particles
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),+1,PolV(0:3,plus,1))
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),-1,PolV(0:3,minus,1))
      Quarks(1)%Mom => MomTd(:,1)
      Quarks(1)%PartType => AUp_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(1)
      Quarks(1)%Mass  => MassTd(1)

      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),+1,PolV(0:3,plus,2))
      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),-1,PolV(0:3,minus,2))
      Quarks(2)%Mom => MomTd(:,2)
      Quarks(2)%PartType => Up_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(2)
      Quarks(2)%Mass  => MassTd(2)

      Quarks(3)%Mom => MomTd(:,3)
      Quarks(3)%PartType => Top_
      Quarks(3)%ExtRef => ExtRef
      Quarks(3)%Mass2 => Mass2Td(3)
      Quarks(3)%Mass  => MassTd(3)

      Quarks(4)%Mom => MomTd(:,4)
      Quarks(4)%PartType => ATop_
      Quarks(4)%ExtRef => ExtRef
      Quarks(4)%Mass2 => Mass2Td(4)
      Quarks(4)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,2
!       Quarks(1)%Pol => PolV(:,HelList(iHel,1),1)
!       Quarks(1)%Helicity => HelList(iHel,1)
      Quarks(2)%Pol => PolV(:,HelList(iHel,2),2)
      Quarks(2)%Helicity => HelList(iHel,2)
!       Quarks(3)%Pol => PolV(:,HelList(iHel,3),3)
!       Quarks(3)%Helicity => HelList(iHel,3)
!       Quarks(4)%Pol => PolV(:,HelList(iHel,4),4)
!       Quarks(4)%Helicity => HelList(iHel,4)

!      calc currents
      Res(0:3) = cur_f_4f((/Gluons(1)/),(/Quarks(2),Quarks(3),Quarks(4)/),Quarks(1)%PartType,(/0,0,0,0,0/),0)
      Amp(1,plus)  = psp1_(Res(0:3),PolV(0:3,plus,1))
      Amp(1,minus) = psp1_(Res(0:3),PolV(0:3,minus,1))

      do col =1,1
      do colP=1,1
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus))*Amp(col,plus) + dconjg(Amp(colP,minus))*Amp(col,minus)  )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

      return
      END FUNCTION






      FUNCTION Tree_UUb_TTb_254_1(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      use modKinematics
      implicit none
      double precision Split_A,Split_B,Split_V(0:3)
      double complex SqAmp,Amp(1,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      double precision :: PSWgt2,PSWgt3,MomDK(0:3,1:6)
      type(Particle),target :: TopQuark(1:2)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:1),Quarks(1:4)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(9)
      real(8),parameter::ColCorr(1:1,1:1)=(/-56.D0/3.D0/)
      integer,target :: HelList(1:16,1:4)
      double complex :: MomTmp(0:3)
      logical :: applyPSCut

      HelList(1,1:4)=(/0,0,0,0/)
      HelList(2,1:4)=(/0,1,0,0/)
!       HelList(3,1:4)=(/0,1,0,0/)
!       HelList(4,1:4)=(/1,1,0,0/)
!       HelList(5,1:4)=(/0,0,1,0/)
!       HelList(6,1:4)=(/1,0,1,0/)
!       HelList(7,1:4)=(/0,1,1,0/)
!       HelList(8,1:4)=(/1,1,1,0/)
!       HelList(9,1:4)=(/0,0,0,1/)
!       HelList(10,1:4)=(/1,0,0,1/)
!       HelList(11,1:4)=(/0,1,0,1/)
!       HelList(12,1:4)=(/1,1,0,1/)
!       HelList(13,1:4)=(/0,0,1,1/)
!       HelList(14,1:4)=(/1,0,1,1/)
!       HelList(15,1:4)=(/0,1,1,1/)
!       HelList(16,1:4)=(/1,1,1,1/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))


      call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
      call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top momenta to be consistent with Kinematics conventions
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin,xJPsiFrag=xFrac)
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


      TopQuark(1)%PartType = Top_
      TopQuark(1)%Mom(1:4) = MomTd(0:3,3)
      call TopDecay(TopQuark(1),DK_LO,MomDK(0:3,4:6))
      TopQuark(2)%PartType = ATop_
      TopQuark(2)%Mom(1:4) = MomTd(0:3,4)
      call TopDecay(TopQuark(2),DK_LO,MomDK(0:3,1:3))

      Quarks(3)%Pol => TopQuark(1)%Pol
      Quarks(4)%Pol => TopQuark(2)%Pol

!      set particles
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),+1,PolV(0:3,plus,1))
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),-1,PolV(0:3,minus,1))
      Quarks(1)%Mom => MomTd(:,1)
      Quarks(1)%PartType => AUp_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(1)
      Quarks(1)%Mass  => MassTd(1)

      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),+1,PolV(0:3,plus,2))
      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),-1,PolV(0:3,minus,2))
      Quarks(2)%Mom => MomTd(:,2)
      Quarks(2)%PartType => Up_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(2)
      Quarks(2)%Mass  => MassTd(2)

      Quarks(3)%Mom => MomTd(:,3)
      Quarks(3)%PartType => Top_
      Quarks(3)%ExtRef => ExtRef
      Quarks(3)%Mass2 => Mass2Td(3)
      Quarks(3)%Mass  => MassTd(3)

      Quarks(4)%Mom => MomTd(:,4)
      Quarks(4)%PartType => ATop_
      Quarks(4)%ExtRef => ExtRef
      Quarks(4)%Mass2 => Mass2Td(4)
      Quarks(4)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,2
!       Quarks(1)%Pol => PolV(:,HelList(iHel,1),1)
!       Quarks(1)%Helicity => HelList(iHel,1)
      Quarks(2)%Pol => PolV(:,HelList(iHel,2),2)
      Quarks(2)%Helicity => HelList(iHel,2)
!       Quarks(3)%Pol => PolV(:,HelList(iHel,3),3)
!       Quarks(3)%Helicity => HelList(iHel,3)
!       Quarks(4)%Pol => PolV(:,HelList(iHel,4),4)
!       Quarks(4)%Helicity => HelList(iHel,4)

!      calc currents
      Res(0:3) = cur_f_4f((/Gluons(1)/),(/Quarks(2),Quarks(3),Quarks(4)/),Quarks(1)%PartType,(/0,0,0,0,0/),0)
      Amp(1,plus)  = psp1_(Res(0:3),PolV(0:3,plus,1))
      Amp(1,minus) = psp1_(Res(0:3),PolV(0:3,minus,1))

      do col =1,1
      do colP=1,1
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus))*Amp(col,plus) + dconjg(Amp(colP,minus))*Amp(col,minus)  )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

      return
      END FUNCTION






      FUNCTION Tree_UUb_TTb_351_5(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      use modKinematics
      implicit none
      double precision Split_A,Split_B,Split_V(0:3)
      double complex SqAmp,Amp(1,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      double precision :: PSWgt2,PSWgt3,MomDK(0:3,1:6)
      type(Particle),target :: TopQuark(1:2)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:1),Quarks(1:4)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(9)
      real(8),parameter::ColCorr(1:1,1:1)=(/-56.D0/3.D0/)
      integer,target :: HelList(1:16,1:4)
      double complex :: MomTmp(0:3)
      logical :: applyPSCut

      HelList(1,1:4)=(/0,0,0,0/)
      HelList(2,1:4)=(/0,1,0,0/)
!       HelList(3,1:4)=(/0,1,0,0/)
!       HelList(4,1:4)=(/1,1,0,0/)
!       HelList(5,1:4)=(/0,0,1,0/)
!       HelList(6,1:4)=(/1,0,1,0/)
!       HelList(7,1:4)=(/0,1,1,0/)
!       HelList(8,1:4)=(/1,1,1,0/)
!       HelList(9,1:4)=(/0,0,0,1/)
!       HelList(10,1:4)=(/1,0,0,1/)
!       HelList(11,1:4)=(/0,1,0,1/)
!       HelList(12,1:4)=(/1,1,0,1/)
!       HelList(13,1:4)=(/0,0,1,1/)
!       HelList(14,1:4)=(/1,0,1,1/)
!       HelList(15,1:4)=(/0,1,1,1/)
!       HelList(16,1:4)=(/1,1,1,1/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))


      call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
      call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top momenta to be consistent with Kinematics conventions
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin,xJPsiFrag=xFrac)
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


      TopQuark(1)%PartType = Top_
      TopQuark(1)%Mom(1:4) = MomTd(0:3,3)
      call TopDecay(TopQuark(1),DK_LO,MomDK(0:3,4:6))
      TopQuark(2)%PartType = ATop_
      TopQuark(2)%Mom(1:4) = MomTd(0:3,4)
      call TopDecay(TopQuark(2),DK_LO,MomDK(0:3,1:3))

      Quarks(3)%Pol => TopQuark(1)%Pol
      Quarks(4)%Pol => TopQuark(2)%Pol


!      set particles
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),+1,PolV(0:3,plus,1))
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),-1,PolV(0:3,minus,1))
      Quarks(1)%Mom => MomTd(:,1)
      Quarks(1)%PartType => AUp_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(1)
      Quarks(1)%Mass  => MassTd(1)

      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),+1,PolV(0:3,plus,2))
      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),-1,PolV(0:3,minus,2))
      Quarks(2)%Mom => MomTd(:,2)
      Quarks(2)%PartType => Up_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(2)
      Quarks(2)%Mass  => MassTd(2)

      Quarks(3)%Mom => MomTd(:,3)
      Quarks(3)%PartType => Top_
      Quarks(3)%ExtRef => ExtRef
      Quarks(3)%Mass2 => Mass2Td(3)
      Quarks(3)%Mass  => MassTd(3)

      Quarks(4)%Mom => MomTd(:,4)
      Quarks(4)%PartType => ATop_
      Quarks(4)%ExtRef => ExtRef
      Quarks(4)%Mass2 => Mass2Td(4)
      Quarks(4)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,2
!       Quarks(1)%Pol => PolV(:,HelList(iHel,1),1)
!       Quarks(1)%Helicity => HelList(iHel,1)
      Quarks(2)%Pol => PolV(:,HelList(iHel,2),2)
      Quarks(2)%Helicity => HelList(iHel,2)
!       Quarks(3)%Pol => PolV(:,HelList(iHel,3),3)
!       Quarks(3)%Helicity => HelList(iHel,3)
!       Quarks(4)%Pol => PolV(:,HelList(iHel,4),4)
!       Quarks(4)%Helicity => HelList(iHel,4)

!      calc currents
      Res(0:3) = cur_f_4f((/Gluons(1)/),(/Quarks(2),Quarks(3),Quarks(4)/),Quarks(1)%PartType,(/0,0,0,0,0/),0)
      Amp(1,plus)  = psp1_(Res(0:3),PolV(0:3,plus,1))
      Amp(1,minus) = psp1_(Res(0:3),PolV(0:3,minus,1))

      do col =1,1
      do colP=1,1
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus))*Amp(col,plus) + dconjg(Amp(colP,minus))*Amp(col,minus)  )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

      return
      END FUNCTION






      FUNCTION Tree_UUb_TTb_352_5(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      use modKinematics
      implicit none
      double precision Split_A,Split_B,Split_V(0:3)
      double complex SqAmp,Amp(1,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      double precision :: PSWgt2,PSWgt3,MomDK(0:3,1:6)
      type(Particle),target :: TopQuark(1:2)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:1),Quarks(1:4)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(9)
      real(8),parameter::ColCorr(1:1,1:1)=(/-16.D0/3.D0/)
      integer,target :: HelList(1:16,1:4)
      double complex :: MomTmp(0:3)
      logical :: applyPSCut

      HelList(1,1:4)=(/0,0,0,0/)
      HelList(2,1:4)=(/0,1,0,0/)
!       HelList(3,1:4)=(/0,1,0,0/)
!       HelList(4,1:4)=(/1,1,0,0/)
!       HelList(5,1:4)=(/0,0,1,0/)
!       HelList(6,1:4)=(/1,0,1,0/)
!       HelList(7,1:4)=(/0,1,1,0/)
!       HelList(8,1:4)=(/1,1,1,0/)
!       HelList(9,1:4)=(/0,0,0,1/)
!       HelList(10,1:4)=(/1,0,0,1/)
!       HelList(11,1:4)=(/0,1,0,1/)
!       HelList(12,1:4)=(/1,1,0,1/)
!       HelList(13,1:4)=(/0,0,1,1/)
!       HelList(14,1:4)=(/1,0,1,1/)
!       HelList(15,1:4)=(/0,1,1,1/)
!       HelList(16,1:4)=(/1,1,1,1/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

      call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
      call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top momenta to be consistent with Kinematics conventions
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin,xJPsiFrag=xFrac)
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


      TopQuark(1)%PartType = Top_
      TopQuark(1)%Mom(1:4) = MomTd(0:3,3)
      call TopDecay(TopQuark(1),DK_LO,MomDK(0:3,4:6))
      TopQuark(2)%PartType = ATop_
      TopQuark(2)%Mom(1:4) = MomTd(0:3,4)
      call TopDecay(TopQuark(2),DK_LO,MomDK(0:3,1:3))

      Quarks(3)%Pol => TopQuark(1)%Pol
      Quarks(4)%Pol => TopQuark(2)%Pol



!      set particles
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),+1,PolV(0:3,plus,1))
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),-1,PolV(0:3,minus,1))
      Quarks(1)%Mom => MomTd(:,1)
      Quarks(1)%PartType => AUp_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(1)
      Quarks(1)%Mass  => MassTd(1)

      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),+1,PolV(0:3,plus,2))
      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),-1,PolV(0:3,minus,2))
      Quarks(2)%Mom => MomTd(:,2)
      Quarks(2)%PartType => Up_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(2)
      Quarks(2)%Mass  => MassTd(2)

      Quarks(3)%Mom => MomTd(:,3)
      Quarks(3)%PartType => Top_
      Quarks(3)%ExtRef => ExtRef
      Quarks(3)%Mass2 => Mass2Td(3)
      Quarks(3)%Mass  => MassTd(3)

      Quarks(4)%Mom => MomTd(:,4)
      Quarks(4)%PartType => ATop_
      Quarks(4)%ExtRef => ExtRef
      Quarks(4)%Mass2 => Mass2Td(4)
      Quarks(4)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,2
!       Quarks(1)%Pol => PolV(:,HelList(iHel,1),1)
!       Quarks(1)%Helicity => HelList(iHel,1)
      Quarks(2)%Pol => PolV(:,HelList(iHel,2),2)
      Quarks(2)%Helicity => HelList(iHel,2)
!       Quarks(3)%Pol => PolV(:,HelList(iHel,3),3)
!       Quarks(3)%Helicity => HelList(iHel,3)
!       Quarks(4)%Pol => PolV(:,HelList(iHel,4),4)
!       Quarks(4)%Helicity => HelList(iHel,4)

!      calc currents
      Res(0:3) = cur_f_4f((/Gluons(1)/),(/Quarks(2),Quarks(3),Quarks(4)/),Quarks(1)%PartType,(/0,0,0,0,0/),0)
      Amp(1,plus)  = psp1_(Res(0:3),PolV(0:3,plus,1))
      Amp(1,minus) = psp1_(Res(0:3),PolV(0:3,minus,1))

      do col =1,1
      do colP=1,1
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus))*Amp(col,plus)  + dconjg(Amp(colP,minus))*Amp(col,minus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

      return
      END FUNCTION






      FUNCTION Tree_UUb_TTb_354_5(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      use modKinematics
      implicit none
      double precision Split_A,Split_B,Split_V(0:3)
      double complex SqAmp,Amp(1,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      double precision :: PSWgt2,PSWgt3,MomDK(0:3,1:6)
      type(Particle),target :: TopQuark(1:2)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:1),Quarks(1:4)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(9)
      real(8),parameter::ColCorr(1:1,1:1)=(/8.D0/3.D0/)
      integer,target :: HelList(1:16,1:4)
      double complex :: MomTmp(0:3)
      logical :: applyPSCut

      HelList(1,1:4)=(/0,0,0,0/)
      HelList(2,1:4)=(/0,1,0,0/)
!       HelList(3,1:4)=(/0,1,0,0/)
!       HelList(4,1:4)=(/1,1,0,0/)
!       HelList(5,1:4)=(/0,0,1,0/)
!       HelList(6,1:4)=(/1,0,1,0/)
!       HelList(7,1:4)=(/0,1,1,0/)
!       HelList(8,1:4)=(/1,1,1,0/)
!       HelList(9,1:4)=(/0,0,0,1/)
!       HelList(10,1:4)=(/1,0,0,1/)
!       HelList(11,1:4)=(/0,1,0,1/)
!       HelList(12,1:4)=(/1,1,0,1/)
!       HelList(13,1:4)=(/0,0,1,1/)
!       HelList(14,1:4)=(/1,0,1,1/)
!       HelList(15,1:4)=(/0,1,1,1/)
!       HelList(16,1:4)=(/1,1,1,1/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))


      call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
      call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top momenta to be consistent with Kinematics conventions
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin,xJPsiFrag=xFrac)
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


      TopQuark(1)%PartType = Top_
      TopQuark(1)%Mom(1:4) = MomTd(0:3,3)
      call TopDecay(TopQuark(1),DK_LO,MomDK(0:3,4:6))
      TopQuark(2)%PartType = ATop_
      TopQuark(2)%Mom(1:4) = MomTd(0:3,4)
      call TopDecay(TopQuark(2),DK_LO,MomDK(0:3,1:3))

      Quarks(3)%Pol => TopQuark(1)%Pol
      Quarks(4)%Pol => TopQuark(2)%Pol



!      set particles
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),+1,PolV(0:3,plus,1))
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),-1,PolV(0:3,minus,1))
      Quarks(1)%Mom => MomTd(:,1)
      Quarks(1)%PartType => AUp_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(1)
      Quarks(1)%Mass  => MassTd(1)

      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),+1,PolV(0:3,plus,2))
      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),-1,PolV(0:3,minus,2))
      Quarks(2)%Mom => MomTd(:,2)
      Quarks(2)%PartType => Up_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(2)
      Quarks(2)%Mass  => MassTd(2)

      Quarks(3)%Mom => MomTd(:,3)
      Quarks(3)%PartType => Top_
      Quarks(3)%ExtRef => ExtRef
      Quarks(3)%Mass2 => Mass2Td(3)
      Quarks(3)%Mass  => MassTd(3)

      Quarks(4)%Mom => MomTd(:,4)
      Quarks(4)%PartType => ATop_
      Quarks(4)%ExtRef => ExtRef
      Quarks(4)%Mass2 => Mass2Td(4)
      Quarks(4)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,2
!       Quarks(1)%Pol => PolV(:,HelList(iHel,1),1)
!       Quarks(1)%Helicity => HelList(iHel,1)
      Quarks(2)%Pol => PolV(:,HelList(iHel,2),2)
      Quarks(2)%Helicity => HelList(iHel,2)
!       Quarks(3)%Pol => PolV(:,HelList(iHel,3),3)
!       Quarks(3)%Helicity => HelList(iHel,3)
!       Quarks(4)%Pol => PolV(:,HelList(iHel,4),4)
!       Quarks(4)%Helicity => HelList(iHel,4)

!      calc currents
      Res(0:3) = cur_f_4f((/Gluons(1)/),(/Quarks(2),Quarks(3),Quarks(4)/),Quarks(1)%PartType,(/0,0,0,0,0/),0)
      Amp(1,plus)  = psp1_(Res(0:3),PolV(0:3,plus,1))
      Amp(1,minus) = psp1_(Res(0:3),PolV(0:3,minus,1))

      do col =1,1
      do colP=1,1
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus))*Amp(col,plus) + dconjg(Amp(colP,minus))*Amp(col,minus)  )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

      return
      END FUNCTION






      FUNCTION Tree_UUb_TTb_451_5(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      use modKinematics
      implicit none
      double precision Split_A,Split_B,Split_V(0:3)
      double complex SqAmp,Amp(1,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      double precision :: PSWgt2,PSWgt3,MomDK(0:3,1:6)
      type(Particle),target :: TopQuark(1:2)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:1),Quarks(1:4)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(9)
      real(8),parameter::ColCorr(1:1,1:1)=(/-16.D0/3.D0/)
      integer,target :: HelList(1:16,1:4)
      double complex :: MomTmp(0:3)
      logical :: applyPSCut

      HelList(1,1:4)=(/0,0,0,0/)
      HelList(2,1:4)=(/0,1,0,0/)
!       HelList(3,1:4)=(/0,1,0,0/)
!       HelList(4,1:4)=(/1,1,0,0/)
!       HelList(5,1:4)=(/0,0,1,0/)
!       HelList(6,1:4)=(/1,0,1,0/)
!       HelList(7,1:4)=(/0,1,1,0/)
!       HelList(8,1:4)=(/1,1,1,0/)
!       HelList(9,1:4)=(/0,0,0,1/)
!       HelList(10,1:4)=(/1,0,0,1/)
!       HelList(11,1:4)=(/0,1,0,1/)
!       HelList(12,1:4)=(/1,1,0,1/)
!       HelList(13,1:4)=(/0,0,1,1/)
!       HelList(14,1:4)=(/1,0,1,1/)
!       HelList(15,1:4)=(/0,1,1,1/)
!       HelList(16,1:4)=(/1,1,1,1/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))

      call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
      call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top momenta to be consistent with Kinematics conventions
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin,xJPsiFrag=xFrac)
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


      TopQuark(1)%PartType = Top_
      TopQuark(1)%Mom(1:4) = MomTd(0:3,3)
      call TopDecay(TopQuark(1),DK_LO,MomDK(0:3,4:6))
      TopQuark(2)%PartType = ATop_
      TopQuark(2)%Mom(1:4) = MomTd(0:3,4)
      call TopDecay(TopQuark(2),DK_LO,MomDK(0:3,1:3))

      Quarks(3)%Pol => TopQuark(1)%Pol
      Quarks(4)%Pol => TopQuark(2)%Pol




!      set particles
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),+1,PolV(0:3,plus,1))
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),-1,PolV(0:3,minus,1))
      Quarks(1)%Mom => MomTd(:,1)
      Quarks(1)%PartType => AUp_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(1)
      Quarks(1)%Mass  => MassTd(1)

      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),+1,PolV(0:3,plus,2))
      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),-1,PolV(0:3,minus,2))
      Quarks(2)%Mom => MomTd(:,2)
      Quarks(2)%PartType => Up_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(2)
      Quarks(2)%Mass  => MassTd(2)

      Quarks(3)%Mom => MomTd(:,3)
      Quarks(3)%PartType => Top_
      Quarks(3)%ExtRef => ExtRef
      Quarks(3)%Mass2 => Mass2Td(3)
      Quarks(3)%Mass  => MassTd(3)

      Quarks(4)%Mom => MomTd(:,4)
      Quarks(4)%PartType => ATop_
      Quarks(4)%ExtRef => ExtRef
      Quarks(4)%Mass2 => Mass2Td(4)
      Quarks(4)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,2
!       Quarks(1)%Pol => PolV(:,HelList(iHel,1),1)
!       Quarks(1)%Helicity => HelList(iHel,1)
      Quarks(2)%Pol => PolV(:,HelList(iHel,2),2)
      Quarks(2)%Helicity => HelList(iHel,2)
!       Quarks(3)%Pol => PolV(:,HelList(iHel,3),3)
!       Quarks(3)%Helicity => HelList(iHel,3)
!       Quarks(4)%Pol => PolV(:,HelList(iHel,4),4)
!       Quarks(4)%Helicity => HelList(iHel,4)

!      calc currents
      Res(0:3) = cur_f_4f((/Gluons(1)/),(/Quarks(2),Quarks(3),Quarks(4)/),Quarks(1)%PartType,(/0,0,0,0,0/),0)
      Amp(1,plus)  = psp1_(Res(0:3),PolV(0:3,plus,1))
      Amp(1,minus) = psp1_(Res(0:3),PolV(0:3,minus,1))

      do col =1,1
      do colP=1,1
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus))*Amp(col,plus) + dconjg(Amp(colP,minus))*Amp(col,minus)  )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

      return
      END FUNCTION






      FUNCTION Tree_UUb_TTb_452_5(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      use modKinematics
      implicit none
      double precision Split_A,Split_B,Split_V(0:3)
      double complex SqAmp,Amp(1,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      double precision :: PSWgt2,PSWgt3,MomDK(0:3,1:6)
      type(Particle),target :: TopQuark(1:2)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:1),Quarks(1:4)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(9)
      real(8),parameter::ColCorr(1:1,1:1)=(/-56.D0/3.D0/)
      integer,target :: HelList(1:16,1:4)
      double complex :: MomTmp(0:3)
      logical :: applyPSCut

      HelList(1,1:4)=(/0,0,0,0/)
      HelList(2,1:4)=(/0,1,0,0/)
!       HelList(3,1:4)=(/0,1,0,0/)
!       HelList(4,1:4)=(/1,1,0,0/)
!       HelList(5,1:4)=(/0,0,1,0/)
!       HelList(6,1:4)=(/1,0,1,0/)
!       HelList(7,1:4)=(/0,1,1,0/)
!       HelList(8,1:4)=(/1,1,1,0/)
!       HelList(9,1:4)=(/0,0,0,1/)
!       HelList(10,1:4)=(/1,0,0,1/)
!       HelList(11,1:4)=(/0,1,0,1/)
!       HelList(12,1:4)=(/1,1,0,1/)
!       HelList(13,1:4)=(/0,0,1,1/)
!       HelList(14,1:4)=(/1,0,1,1/)
!       HelList(15,1:4)=(/0,1,1,1/)
!       HelList(16,1:4)=(/1,1,1,1/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))


      call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
      call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top momenta to be consistent with Kinematics conventions
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin,xJPsiFrag=xFrac)
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


      TopQuark(1)%PartType = Top_
      TopQuark(1)%Mom(1:4) = MomTd(0:3,3)
      call TopDecay(TopQuark(1),DK_LO,MomDK(0:3,4:6))
      TopQuark(2)%PartType = ATop_
      TopQuark(2)%Mom(1:4) = MomTd(0:3,4)
      call TopDecay(TopQuark(2),DK_LO,MomDK(0:3,1:3))

      Quarks(3)%Pol => TopQuark(1)%Pol
      Quarks(4)%Pol => TopQuark(2)%Pol


!      set particles
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),+1,PolV(0:3,plus,1))
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),-1,PolV(0:3,minus,1))
      Quarks(1)%Mom => MomTd(:,1)
      Quarks(1)%PartType => AUp_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(1)
      Quarks(1)%Mass  => MassTd(1)

      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),+1,PolV(0:3,plus,2))
      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),-1,PolV(0:3,minus,2))
      Quarks(2)%Mom => MomTd(:,2)
      Quarks(2)%PartType => Up_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(2)
      Quarks(2)%Mass  => MassTd(2)

      Quarks(3)%Mom => MomTd(:,3)
      Quarks(3)%PartType => Top_
      Quarks(3)%ExtRef => ExtRef
      Quarks(3)%Mass2 => Mass2Td(3)
      Quarks(3)%Mass  => MassTd(3)

      Quarks(4)%Mom => MomTd(:,4)
      Quarks(4)%PartType => ATop_
      Quarks(4)%ExtRef => ExtRef
      Quarks(4)%Mass2 => Mass2Td(4)
      Quarks(4)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,2
!       Quarks(1)%Pol => PolV(:,HelList(iHel,1),1)
!       Quarks(1)%Helicity => HelList(iHel,1)
      Quarks(2)%Pol => PolV(:,HelList(iHel,2),2)
      Quarks(2)%Helicity => HelList(iHel,2)
!       Quarks(3)%Pol => PolV(:,HelList(iHel,3),3)
!       Quarks(3)%Helicity => HelList(iHel,3)
!       Quarks(4)%Pol => PolV(:,HelList(iHel,4),4)
!       Quarks(4)%Helicity => HelList(iHel,4)

!      calc currents
      Res(0:3) = cur_f_4f((/Gluons(1)/),(/Quarks(2),Quarks(3),Quarks(4)/),Quarks(1)%PartType,(/0,0,0,0,0/),0)
      Amp(1,plus)  = psp1_(Res(0:3),PolV(0:3,plus,1))
      Amp(1,minus) = psp1_(Res(0:3),PolV(0:3,minus,1))

      do col =1,1
      do colP=1,1
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus))*Amp(col,plus) + dconjg(Amp(colP,minus))*Amp(col,minus)  )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

      return
      END FUNCTION






      FUNCTION Tree_UUb_TTb_453_5(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      use modKinematics
      implicit none
      double precision Split_A,Split_B,Split_V(0:3)
      double complex SqAmp,Amp(1,0:1),Res(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      double precision :: PSWgt2,PSWgt3,MomDK(0:3,1:6)
      type(Particle),target :: TopQuark(1:2)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:1),Quarks(1:4)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(9)
      real(8),parameter::ColCorr(1:1,1:1)=(/8.D0/3.D0/)
      integer,target :: HelList(1:16,1:4)
      double complex :: MomTmp(0:3)
      logical :: applyPSCut

      HelList(1,1:4)=(/0,0,0,0/)
      HelList(2,1:4)=(/0,1,0,0/)
!       HelList(3,1:4)=(/0,1,0,0/)
!       HelList(4,1:4)=(/1,1,0,0/)
!       HelList(5,1:4)=(/0,0,1,0/)
!       HelList(6,1:4)=(/1,0,1,0/)
!       HelList(7,1:4)=(/0,1,1,0/)
!       HelList(8,1:4)=(/1,1,1,0/)
!       HelList(9,1:4)=(/0,0,0,1/)
!       HelList(10,1:4)=(/1,0,0,1/)
!       HelList(11,1:4)=(/0,1,0,1/)
!       HelList(12,1:4)=(/1,1,0,1/)
!       HelList(13,1:4)=(/0,0,1,1/)
!       HelList(14,1:4)=(/1,0,1,1/)
!       HelList(15,1:4)=(/0,1,1,1/)
!       HelList(16,1:4)=(/1,1,1,1/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))


      call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3)
      call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2)

      MomTmp(0:3)  = MomTd(0:3,3)! swap ATop and Top momenta to be consistent with Kinematics conventions
      MomTd(0:3,3) = MomTd(0:3,4)
      MomTd(0:3,4) = MomTmp(0:3)
      call Kinematics_TTBAR(.false.,dble(MomTd),MomDK,applyPSCut,NBin,xJPsiFrag=xFrac)
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


      TopQuark(1)%PartType = Top_
      TopQuark(1)%Mom(1:4) = MomTd(0:3,3)
      call TopDecay(TopQuark(1),DK_LO,MomDK(0:3,4:6))
      TopQuark(2)%PartType = ATop_
      TopQuark(2)%Mom(1:4) = MomTd(0:3,4)
      call TopDecay(TopQuark(2),DK_LO,MomDK(0:3,1:3))

      Quarks(3)%Pol => TopQuark(1)%Pol
      Quarks(4)%Pol => TopQuark(2)%Pol


!      set particles
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),+1,PolV(0:3,plus,1))
      call vSpi(MomTd(0:3,1),dsqrt(Mass2Td(1)),-1,PolV(0:3,minus,1))
      Quarks(1)%Mom => MomTd(:,1)
      Quarks(1)%PartType => AUp_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(1)
      Quarks(1)%Mass  => MassTd(1)

      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),+1,PolV(0:3,plus,2))
      call ubarSpi(MomTd(0:3,2),dsqrt(Mass2Td(2)),-1,PolV(0:3,minus,2))
      Quarks(2)%Mom => MomTd(:,2)
      Quarks(2)%PartType => Up_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(2)
      Quarks(2)%Mass  => MassTd(2)

      Quarks(3)%Mom => MomTd(:,3)
      Quarks(3)%PartType => Top_
      Quarks(3)%ExtRef => ExtRef
      Quarks(3)%Mass2 => Mass2Td(3)
      Quarks(3)%Mass  => MassTd(3)

      Quarks(4)%Mom => MomTd(:,4)
      Quarks(4)%PartType => ATop_
      Quarks(4)%ExtRef => ExtRef
      Quarks(4)%Mass2 => Mass2Td(4)
      Quarks(4)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,2
!       Quarks(1)%Pol => PolV(:,HelList(iHel,1),1)
!       Quarks(1)%Helicity => HelList(iHel,1)
      Quarks(2)%Pol => PolV(:,HelList(iHel,2),2)
      Quarks(2)%Helicity => HelList(iHel,2)
!       Quarks(3)%Pol => PolV(:,HelList(iHel,3),3)
!       Quarks(3)%Helicity => HelList(iHel,3)
!       Quarks(4)%Pol => PolV(:,HelList(iHel,4),4)
!       Quarks(4)%Helicity => HelList(iHel,4)

!      calc currents
      Res(0:3) = cur_f_4f((/Gluons(1)/),(/Quarks(2),Quarks(3),Quarks(4)/),Quarks(1)%PartType,(/0,0,0,0,0/),0)
      Amp(1,plus)  = psp1_(Res(0:3),PolV(0:3,plus,1))
      Amp(1,minus) = psp1_(Res(0:3),PolV(0:3,minus,1))

      do col =1,1
      do colP=1,1
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*( dconjg(Amp(colP,plus))*Amp(col,plus) + dconjg(Amp(colP,minus))*Amp(col,minus)  )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

      return
      END FUNCTION








      END MODULE
