#DEFINE _CHECK_DIPOLE_MOMMAP 0
#DEFINE _APPLY_CUTS 0

      MODULE ModDipoles_GGSTSTBG
      use ModParameters
      use ModTopdecay
      implicit none

      integer,private,parameter :: NumMaxHisto=45
!
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

      double precision, private :: yRndDK(8:19)
      integer, private :: NBin(1:NumMaxHisto)

      logical, parameter :: invert_alphaCut = .false.


      contains







      SUBROUTINE EvalDipoles_GGSTSTBG(Mom,Mass2,yRnd,Wgt,ResultSum)
!     evaluate dipole contributions
      use ModKinematics
      use ModMisc
      use ModParameters
      implicit none
      double precision Mom(0:3,1:5),Mass2(1:5),Wgt
      double precision, parameter :: DipMinus=-1d0*0.5d0
      integer NHisto
      double precision ResultSum,yRnd(8:19)
      double precision RunFactor

      yRndDK(8:19)=yRnd(8:19)


if(alpha_ff.ne.1d0) call Error("alpha_ff.ne.1d0 is not yet implemented")


!     evaluate all dipoles and do binning:
      ResultSum = (0d0)
      TheDipoles(1)%alphaCut = alpha_ii
      call evalDipole1_152(Mom,Mass2,TheDipoles(1))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(1)%DipoleValue = DipMinus * TheDipoles(1)%DipoleValue * Wgt * RunFactor**3
      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(1)%DipoleValue)
      enddo
      ResultSum = ResultSum + TheDipoles(1)%DipoleValue


      TheDipoles(2)%alphaCut = alpha_if
      call evalDipole2_153(Mom,Mass2,TheDipoles(2))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(2)%DipoleValue = DipMinus *  TheDipoles(2)%DipoleValue * Wgt * RunFactor**3
      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(2)%DipoleValue)
      enddo
      ResultSum = ResultSum + TheDipoles(2)%DipoleValue


      TheDipoles(3)%alphaCut = alpha_if
      call evalDipole3_154(Mom,Mass2,TheDipoles(3))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(3)%DipoleValue = TheDipoles(3)%DipoleValue * Wgt * DipMinus * RunFactor**3
      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(3)%DipoleValue)
      enddo
      ResultSum = ResultSum + TheDipoles(3)%DipoleValue


      TheDipoles(4)%alphaCut = alpha_ii
      call evalDipole4_251(Mom,Mass2,TheDipoles(4))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(4)%DipoleValue = TheDipoles(4)%DipoleValue * Wgt * DipMinus * RunFactor**3
      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(4)%DipoleValue)
      enddo
      ResultSum = ResultSum + TheDipoles(4)%DipoleValue


      TheDipoles(5)%alphaCut = alpha_if
      call evalDipole5_253(Mom,Mass2,TheDipoles(5))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(5)%DipoleValue = TheDipoles(5)%DipoleValue * Wgt * DipMinus * RunFactor**3
      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(5)%DipoleValue)
      enddo
      ResultSum = ResultSum + TheDipoles(5)%DipoleValue


      TheDipoles(6)%alphaCut = alpha_if
      call evalDipole6_254(Mom,Mass2,TheDipoles(6))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(6)%DipoleValue = TheDipoles(6)%DipoleValue * Wgt * DipMinus * RunFactor**3
      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(6)%DipoleValue)
      enddo
      ResultSum = ResultSum + TheDipoles(6)%DipoleValue


      TheDipoles(7)%alphaCut = alpha_fi
      call evalDipole7_351(Mom,Mass2,TheDipoles(7))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(7)%DipoleValue = TheDipoles(7)%DipoleValue * Wgt * DipMinus * RunFactor**3
      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(7)%DipoleValue)
      enddo
      ResultSum = ResultSum + TheDipoles(7)%DipoleValue


      TheDipoles(8)%alphaCut = alpha_fi
      call evalDipole8_352(Mom,Mass2,TheDipoles(8))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(8)%DipoleValue = TheDipoles(8)%DipoleValue * Wgt * DipMinus * RunFactor**3
      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(8)%DipoleValue)
      enddo
      ResultSum = ResultSum + TheDipoles(8)%DipoleValue


      TheDipoles(9)%alphaCut = alpha_ff
      call evalDipole9_354(Mom,Mass2,TheDipoles(9))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(9)%DipoleValue = TheDipoles(9)%DipoleValue * Wgt * DipMinus * RunFactor**3
      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(9)%DipoleValue)
      enddo
      ResultSum = ResultSum + TheDipoles(9)%DipoleValue


      TheDipoles(10)%alphaCut = alpha_fi
      call evalDipole10_451(Mom,Mass2,TheDipoles(10))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(10)%DipoleValue = TheDipoles(10)%DipoleValue * Wgt * DipMinus * RunFactor**3
      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(10)%DipoleValue)
      enddo
      ResultSum = ResultSum + TheDipoles(10)%DipoleValue


      TheDipoles(11)%alphaCut = alpha_fi
      call evalDipole11_452(Mom,Mass2,TheDipoles(11))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(11)%DipoleValue = TheDipoles(11)%DipoleValue * Wgt * DipMinus * RunFactor**3
      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(11)%DipoleValue)
      enddo
      ResultSum = ResultSum + TheDipoles(11)%DipoleValue


      TheDipoles(12)%alphaCut = alpha_ff
      call evalDipole12_453(Mom,Mass2,TheDipoles(12))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(12)%DipoleValue = TheDipoles(12)%DipoleValue * Wgt * DipMinus * RunFactor**3
      do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),TheDipoles(12)%DipoleValue)
      enddo
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
      if( .not. invert_alphacut ) then
          if( TheDipole%AlphaCut.lt.v ) then
            TheDipole%DipoleValue = 0d0
            return
          endif
      else
          if( TheDipole%AlphaCut.ge.v ) then
            TheDipole%DipoleValue = 0d0
            return
          endif
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


!     splitting 2: Gluon_I --> Gluon_I + Gluon_F
      Split_A =-4d0*alpha_s4Pi * ( x/(1d0-x)+x*(1d0-x) )
      Split_B = 4d0*alpha_s4Pi * ( (1d0/x-1d0)*(2d0*s12/s15/s25) )
      Split_V(0:3) = Mom(0:3,5) - s15/s12*Mom(0:3,2)

!     color and spin correlated tree process
!       SqAmp = Tree_GG_TTb_152_2(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_GG_TTb(1,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
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
      if( .not. invert_alphacut ) then
          if( TheDipole%AlphaCut.lt.u ) then
            TheDipole%DipoleValue = 0d0
            return
          endif
      else
          if( TheDipole%AlphaCut.ge.u ) then
            TheDipole%DipoleValue = 0d0
            return
          endif
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


!     splitting 2: Gluon_I --> Gluon_I + Gluon_F
      Split_A =-4d0*alpha_s4Pi * ( 1d0/(1d0-x+u)-1d0+x*(1d0-x) )
      Split_B = 4d0*alpha_s4Pi * ( (1d0/x-1d0)*(2d0*u*(1d0-u)/s35))
      Split_V(0:3) = Mom(0:3,5)/u - Mom(0:3,3)/(1d0-u)

!     color and spin correlated tree process
!       SqAmp = Tree_GG_TTb_153_2(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_GG_TTb(2,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

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
      if( .not. invert_alphacut ) then
          if( TheDipole%AlphaCut.lt.u ) then
            TheDipole%DipoleValue = 0d0
            return
          endif
      else
          if( TheDipole%AlphaCut.ge.u ) then
            TheDipole%DipoleValue = 0d0
            return
          endif
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


!     splitting 2: Gluon_I --> Gluon_I + Gluon_F
      Split_A =-4d0*alpha_s4Pi * ( 1d0/(1d0-x+u)-1d0+x*(1d0-x) )
      Split_B = 4d0*alpha_s4Pi * ( (1d0/x-1d0)*(2d0*u*(1d0-u)/s45))
      Split_V(0:3) = Mom(0:3,5)/u - Mom(0:3,4)/(1d0-u)

!     color and spin correlated tree process
!       SqAmp = Tree_GG_TTb_154_2(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_GG_TTb(3,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

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
      if( .not. invert_alphacut ) then
          if( TheDipole%AlphaCut.lt.v ) then
            TheDipole%DipoleValue = 0d0
            return
          endif
      else
          if( TheDipole%AlphaCut.ge.v ) then
            TheDipole%DipoleValue = 0d0
            return
          endif
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


!     splitting 2: Gluon_I --> Gluon_I + Gluon_F
      Split_A =-4d0*alpha_s4Pi * ( x/(1d0-x)+x*(1d0-x) )
      Split_B = 4d0*alpha_s4Pi * ( (1d0/x-1d0)*(2d0*s12/s25/s15) )
      Split_V(0:3) = Mom(0:3,5) - s25/s12*Mom(0:3,1)

!     color and spin correlated tree process
!       SqAmp = Tree_GG_TTb_251_2(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_GG_TTb(4,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

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
      if( .not. invert_alphacut ) then
          if( TheDipole%AlphaCut.lt.u ) then
            TheDipole%DipoleValue = 0d0
            return
          endif
      else
          if( TheDipole%AlphaCut.ge.u ) then
            TheDipole%DipoleValue = 0d0
            return
          endif
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


!     splitting 2: Gluon_I --> Gluon_I + Gluon_F
      Split_A =-4d0*alpha_s4Pi * ( 1d0/(1d0-x+u)-1d0+x*(1d0-x) )
      Split_B = 4d0*alpha_s4Pi * ( (1d0/x-1d0)*(2d0*u*(1d0-u)/s35))
      Split_V(0:3) = Mom(0:3,5)/u - Mom(0:3,3)/(1d0-u)

!     color and spin correlated tree process
!       SqAmp = Tree_GG_TTb_253_2(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_GG_TTb(5,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

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
      if( .not. invert_alphacut ) then
          if( TheDipole%AlphaCut.lt.u ) then
            TheDipole%DipoleValue = 0d0
            return
          endif
      else
          if( TheDipole%AlphaCut.ge.u ) then
            TheDipole%DipoleValue = 0d0
            return
          endif
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


!     splitting 2: Gluon_I --> Gluon_I + Gluon_F
      Split_A =-4d0*alpha_s4Pi * ( 1d0/(1d0-x+u)-1d0+x*(1d0-x) )
      Split_B = 4d0*alpha_s4Pi * ( (1d0/x-1d0)*(2d0*u*(1d0-u)/s45))
      Split_V(0:3) = Mom(0:3,5)/u - Mom(0:3,4)/(1d0-u)

!     color and spin correlated tree process
!       SqAmp = Tree_GG_TTb_254_2(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_GG_TTb(6,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

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
      if( .not. invert_alphacut ) then
          if( TheDipole%AlphaCut.lt.1d0-x ) then  ! NEW
            TheDipole%DipoleValue = 0d0
            return
          endif
      else
          if( TheDipole%AlphaCut.ge.1d0-x ) then  ! NEW
            TheDipole%DipoleValue = 0d0
            return
          endif
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
!       SqAmp = Tree_GG_TTb_351_5(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_GG_TTb(7,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

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
      if( .not. invert_alphacut ) then
          if( TheDipole%AlphaCut.lt.1d0-x   ) then  ! NEW
            TheDipole%DipoleValue = 0d0
            return
          endif
      else
          if( TheDipole%AlphaCut.ge.1d0-x   ) then  ! NEW
            TheDipole%DipoleValue = 0d0
            return
          endif
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
!       SqAmp = Tree_GG_TTb_352_5(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_GG_TTb(8,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

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
!       SqAmp = Tree_GG_TTb_354_5(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_GG_TTb(9,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

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
      if( .not. invert_alphacut ) then
          if( TheDipole%AlphaCut.lt.1d0-x ) then  ! NEW
            TheDipole%DipoleValue = 0d0
            return
          endif
      else
          if( TheDipole%AlphaCut.ge.1d0-x ) then  ! NEW
            TheDipole%DipoleValue = 0d0
            return
          endif
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
!       SqAmp = Tree_GG_TTb_451_5(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_GG_TTb(10,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

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
      if( .not. invert_alphacut ) then
          if( TheDipole%AlphaCut.lt.1d0-x ) then  ! NEW
            TheDipole%DipoleValue = 0d0
            return
          endif
      else
          if( TheDipole%AlphaCut.ge.1d0-x ) then  ! NEW
            TheDipole%DipoleValue = 0d0
            return
          endif
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
!       SqAmp = Tree_GG_TTb_452_5(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_GG_TTb(11,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

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
!       SqAmp = Tree_GG_TTb_453_5(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)
      SqAmp = Tree_GG_TTb(12,dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V)

      TheDipole%DipoleValue = LeadSing * SqAmp
      return
      END SUBROUTINE




      FUNCTION Tree_GG_TTb(icorr,MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      use modKinematics
      use ModExoticDecay
      implicit none
      double precision Split_A,Split_B,Split_V(0:3)
      double complex SqAmp,Amp(2,0:1),Res(0:3)
      integer col,colP,iHel,icorr,ngl,A0barHel,A0Hel,SecHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      double precision :: PSWgt2,PSWgt3,PSWgt4,PSWgt5,MomExt(1:4,1:16)
      type(Particle),target :: TopQuark(1:2)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      real(8)::ColCorr(1:2,1:2)
      real(8),parameter::ColCorr1(1:2,1:2)=(/-72.D0,0.D0,0.D0,-72.D0/)
      real(8),parameter::ColCorr2(1:2,1:2)=(/-64.D0,8.D0,8.D0,8.D0/)
      real(8),parameter::ColCorr3(1:2,1:2)=(/8.D0,8.D0,8.D0,-64.D0/)
      real(8),parameter::ColCorr4(1:2,1:2)=(/-72.D0,0.D0,0.D0,-72.D0/)
      real(8),parameter::ColCorr5(1:2,1:2)=(/8.D0,8.D0,8.D0,-64.D0/)
      real(8),parameter::ColCorr6(1:2,1:2)=(/-64.D0,8.D0,8.D0,8.D0/)
      real(8),parameter::ColCorr7(1:2,1:2)=(/-64.D0,8.D0,8.D0,8.D0/)
      real(8),parameter::ColCorr8(1:2,1:2)=(/8.D0,8.D0,8.D0,-64.D0/)
      real(8),parameter::ColCorr9(1:2,1:2)=(/-8.D0/9.D0,-80.D0/9.D0,-80.D0/9.D0,-8.D0/9.D0/)
      real(8),parameter::ColCorr10(1:2,1:2)=(/8.D0,8.D0,8.D0,-64.D0/)
      real(8),parameter::ColCorr11(1:2,1:2)=(/-64.D0,8.D0,8.D0,8.D0/)
      real(8),parameter::ColCorr12(1:2,1:2)=(/-8.D0/9.D0,-80.D0/9.D0,-80.D0/9.D0,-8.D0/9.D0/)
      integer,target :: HelList(1:4,1:3)
      double complex :: MomTmp(0:3)
      logical :: applyPSCut

      SqAmp = 0d0
      if(icorr.eq.1) ColCorr(1:2,1:2)=ColCorr1(1:2,1:2)
      if(icorr.eq.2) ColCorr(1:2,1:2)=ColCorr2(1:2,1:2)
      if(icorr.eq.3) ColCorr(1:2,1:2)=ColCorr3(1:2,1:2)
      if(icorr.eq.4) ColCorr(1:2,1:2)=ColCorr4(1:2,1:2)
      if(icorr.eq.5) ColCorr(1:2,1:2)=ColCorr5(1:2,1:2)
      if(icorr.eq.6) ColCorr(1:2,1:2)=ColCorr6(1:2,1:2)
      if(icorr.eq.7) ColCorr(1:2,1:2)=ColCorr7(1:2,1:2)
      if(icorr.eq.8) ColCorr(1:2,1:2)=ColCorr8(1:2,1:2)
      if(icorr.eq.9) ColCorr(1:2,1:2)=ColCorr9(1:2,1:2)
      if(icorr.eq.10) ColCorr(1:2,1:2)=ColCorr10(1:2,1:2)
      if(icorr.eq.11) ColCorr(1:2,1:2)=ColCorr11(1:2,1:2)
      if(icorr.eq.12) ColCorr(1:2,1:2)=ColCorr12(1:2,1:2)

      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))


      MomExt(1:4,1)=MomTd(0:3,1)
      MomExt(1:4,2)=MomTd(0:3,2)
      MomExt(1:4,3)=MomTd(0:3,4)
      MomExt(1:4,4)=MomTd(0:3,3)

      TopQuark(1)%PartType = STop_
      TopQuark(1)%Mom(1:4) = MomTd(0:3,3)
      TopQuark(2)%PartType = ASTop_
      TopQuark(2)%Mom(1:4) = MomTd(0:3,4)

          IF(XTOPDECAYS.EQ.3) THEN
              call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,3),yRndDK(8:9),MomExt(1:4,5:6),PSWgt2)!  Chi top
              call EvalPhasespace_StopDK(ST_Chi0_T,MomExt(1:4,4),yRndDK(10:11),MomExt(1:4,10:11),PSWgt3)
              call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,6),yRndDK(12:15),MomExt(1:4,7:9),PSWgt4)! bot lep neu
              call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,11),yRndDK(16:19),MomExt(1:4,12:14),PSWgt5)
          ENDIF

          call Kinematics_TTbarETmiss(.false.,MomExt,(/3,4,5,10,6,11,7,8,9,12,13,14,0/),applyPSCut,NBin)
          if( applyPSCut ) then
              SqAmp = 0d0
              return
          endif


!      momentum crossing
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

      if( icorr.eq.4 .or. icorr.eq.5 .or. icorr.eq.6 .or. icorr.eq.10 .or. icorr.eq.11.or. icorr.eq.12 ) then
        ngl=1
      else
        ngl=2
      endif
      if( icorr.eq.7 .or. icorr.eq.8 .or. icorr.eq.9 .or. icorr.eq.10 .or. icorr.eq.11.or. icorr.eq.12) then
          Split_A=-Split_A
      endif

!      set particles
      call pol_mless(MomTd(0:3,1),+1,PolV(0:3,plus,1))
      call pol_mless(MomTd(0:3,1),-1,PolV(0:3,minus,1))

      call pol_mless(MomTd(0:3,2),+1,PolV(0:3,plus,2))
      call pol_mless(MomTd(0:3,2),-1,PolV(0:3,minus,2))
      Gluons(ngl)%Mom => MomTd(:,ngl)
      Gluons(ngl)%PartType => Glu_
      Gluons(ngl)%ExtRef => ExtRef
      Gluons(ngl)%Mass2 => Mass2Td(ngl)
      Gluons(ngl)%Mass  => MassTd(ngl)

      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => STop_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ASTop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)


!      sum over helicities
      if( XTopDecays.eq.0 ) then
            SecHel=1
      else
            SecHel=-1
      endif

      SqAmp = (0d0,0d0)
      do A0barHel=1,SecHel,-2
      do A0Hel=1,SecHel,-2

            call STopDecay(TopQuark(2),DKX_STChi0_LO,A0barHel,MomExt(1:4,5:9))
            call STopDecay(TopQuark(1),DKX_STChi0_LO,A0Hel,MomExt(1:4,10:14))


      Quarks(1)%Pol => TopQuark(1)%Pol
      Quarks(2)%Pol => TopQuark(2)%Pol

      do iHel=1,2
      Gluons(ngl)%Pol => PolV(:,HelList(iHel,1),ngl)
      Gluons(ngl)%Helicity => HelList(iHel,1)

!      calc currents
      Res(0:3) = cur_g_2s( (/Gluons(ngl)/),(/Quarks(2),Quarks(1)/),(/2,1,0,0/) )
      Amp(3-ngl,plus)  = Res(0:3).dot.PolV(0:3,plus,3-ngl)
      Amp(3-ngl,minus) = Res(0:3).dot.PolV(0:3,minus,3-ngl)

      Res(0:3) = cur_g_2s( (/Gluons(ngl)/),(/Quarks(2),Quarks(1)/),(/2,0,0,1/) )
      Amp(ngl,plus)  = Res(0:3).dot.PolV(0:3,plus,3-ngl)
      Amp(ngl,minus) = Res(0:3).dot.PolV(0:3,minus,3-ngl)

      EpsDotP = PolV(0:3,plus,3-ngl).dot.Split_V(0:3)
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
      enddo
      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg *SqAmp



      return
      END FUNCTION







      END MODULE
