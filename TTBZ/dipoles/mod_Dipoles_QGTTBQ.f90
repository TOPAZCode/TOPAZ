#DEFINE _CHECK_DIPOLE_MOMMAP 0
#DEFINE _APPLY_CUTS 0

      MODULE ModDipoles_QGTTBQ
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

      type(Dipole),private :: TheDipoles(1:6)

      double precision, private :: yRndDK(1:8),xFrac,r_sc
      integer, private :: NBin(1:20)

      contains







      FUNCTION applyDipoleCut(MomTd)
      implicit none
      logical :: applyDipoleCut
      double precision MomTd(0:3,1:4)



      applyDipoleCut=.false.
      return
      END FUNCTION




      SUBROUTINE EvalDipoles_QGTTBQ(Mom,Mass2,yRnd,Wgt,ResultSum)
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
      r_sc = 0d0
      TheDipoles(1)%alphaCut = alpha_ii
      call evalDipole1_15x(Mom,Mass2,TheDipoles(1))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(1)%DipoleValue = DipMinus * Wgt * TheDipoles(1)%DipoleValue * RunFactor**3
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


      r_sc = 0d0
      TheDipoles(4)%alphaCut = alpha_ii
      call evalDipole4_25x(Mom,Mass2,TheDipoles(4))
      RunFactor = RunAlphaS(2,MuRen)
      TheDipoles(4)%DipoleValue = DipMinus * Wgt * TheDipoles(4)%DipoleValue * RunFactor**3
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


!       TheDipoles(2)%alphaCut = 1d0
!       call evalDipole2_153(Mom,Mass2,TheDipoles(2))
!       RunFactor = RunAlphaS(2,MuRen)
!       TheDipoles(2)%DipoleValue = DipMinus * ColorConv * Wgt * TheDipoles(2)%DipoleValue * RunFactor**3
!       pT_tilde = dsqrt( TheDipoles(2)%MomTd(1,3)**2 + TheDipoles(2)%MomTd(2,3)**2 )
!       Bin_PT = WhichBin(Histo_PT,pT_tilde)
!       call intoHisto(Histo_PT,Bin_PT,TheDipoles(2)%DipoleValue)
!       M_tilde = dsqrt( TheDipoles(2)%Mass2Td(3) + TheDipoles(2)%Mass2Td(4) +2d0*TheDipoles(2)%MomTd(0:3,3).dot.TheDipoles(2)%MomTd(0:3,4) )
!       Bin_M = WhichBin(Histo_M,M_tilde)
!       call intoHisto(Histo_M,Bin_M,TheDipoles(2)%DipoleValue)
!       ResultSum = ResultSum + TheDipoles(2)%DipoleValue
!
!
!       TheDipoles(3)%alphaCut = 1d0
!       call evalDipole3_154(Mom,Mass2,TheDipoles(3))
!       RunFactor = RunAlphaS(2,MuRen)
!       TheDipoles(3)%DipoleValue = DipMinus * ColorConv * Wgt * TheDipoles(3)%DipoleValue * RunFactor**3
!       pT_tilde = dsqrt( TheDipoles(3)%MomTd(1,3)**2 + TheDipoles(3)%MomTd(2,3)**2 )
!       Bin_PT = WhichBin(Histo_PT,pT_tilde)
!       call intoHisto(Histo_PT,Bin_PT,TheDipoles(3)%DipoleValue)
!       M_tilde = dsqrt( TheDipoles(3)%Mass2Td(3) + TheDipoles(3)%Mass2Td(4) +2d0*TheDipoles(3)%MomTd(0:3,3).dot.TheDipoles(3)%MomTd(0:3,4) )
!       Bin_M = WhichBin(Histo_M,M_tilde)
!       call intoHisto(Histo_M,Bin_M,TheDipoles(3)%DipoleValue)
!       ResultSum = ResultSum + TheDipoles(3)%DipoleValue
!
!
!       TheDipoles(4)%alphaCut = 1d0
!       call evalDipole4_251(Mom,Mass2,TheDipoles(4))
!       RunFactor = RunAlphaS(2,MuRen)
!       TheDipoles(4)%DipoleValue = DipMinus * ColorConv * Wgt * TheDipoles(4)%DipoleValue * RunFactor**3
!       pT_tilde = dsqrt( TheDipoles(4)%MomTd(1,3)**2 + TheDipoles(4)%MomTd(2,3)**2 )
!       Bin_PT = WhichBin(Histo_PT,pT_tilde)
!       call intoHisto(Histo_PT,Bin_PT,TheDipoles(4)%DipoleValue)
!       M_tilde = dsqrt( TheDipoles(4)%Mass2Td(3) + TheDipoles(4)%Mass2Td(4) +2d0*TheDipoles(4)%MomTd(0:3,3).dot.TheDipoles(4)%MomTd(0:3,4) )
!       Bin_M = WhichBin(Histo_M,M_tilde)
!       call intoHisto(Histo_M,Bin_M,TheDipoles(4)%DipoleValue)
!       ResultSum = ResultSum + TheDipoles(4)%DipoleValue
!
!
!       TheDipoles(5)%alphaCut = 1d0
!       call evalDipole5_253(Mom,Mass2,TheDipoles(5))
!       RunFactor = RunAlphaS(2,MuRen)
!       TheDipoles(5)%DipoleValue = DipMinus * ColorConv * Wgt * TheDipoles(5)%DipoleValue * RunFactor**3
!       pT_tilde = dsqrt( TheDipoles(5)%MomTd(1,3)**2 + TheDipoles(5)%MomTd(2,3)**2 )
!       Bin_PT = WhichBin(Histo_PT,pT_tilde)
!       call intoHisto(Histo_PT,Bin_PT,TheDipoles(5)%DipoleValue)
!       M_tilde = dsqrt( TheDipoles(5)%Mass2Td(3) + TheDipoles(5)%Mass2Td(4) +2d0*TheDipoles(5)%MomTd(0:3,3).dot.TheDipoles(5)%MomTd(0:3,4) )
!       Bin_M = WhichBin(Histo_M,M_tilde)
!       call intoHisto(Histo_M,Bin_M,TheDipoles(5)%DipoleValue)
!       ResultSum = ResultSum + TheDipoles(5)%DipoleValue
!
!
!       TheDipoles(6)%alphaCut = 1d0
!       call evalDipole6_254(Mom,Mass2,TheDipoles(6))
!       RunFactor = RunAlphaS(2,MuRen)
!       TheDipoles(6)%DipoleValue = DipMinus * ColorConv * Wgt * TheDipoles(6)%DipoleValue * RunFactor**3
!       pT_tilde = dsqrt( TheDipoles(6)%MomTd(1,3)**2 + TheDipoles(6)%MomTd(2,3)**2 )
!       Bin_PT = WhichBin(Histo_PT,pT_tilde)
!       call intoHisto(Histo_PT,Bin_PT,TheDipoles(6)%DipoleValue)
!       M_tilde = dsqrt( TheDipoles(6)%Mass2Td(3) + TheDipoles(6)%Mass2Td(4) +2d0*TheDipoles(6)%MomTd(0:3,3).dot.TheDipoles(6)%MomTd(0:3,4) )
!       Bin_M = WhichBin(Histo_M,M_tilde)
!       call intoHisto(Histo_M,Bin_M,TheDipoles(6)%DipoleValue)
!       ResultSum = ResultSum + TheDipoles(6)%DipoleValue



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





      SUBROUTINE evalDipole1_15x(Mom,Mass2,TheDipole)
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
      TheDipole%Mass2Td(1) = 0d0
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

!     splitting 3: A/Quark_I --> Gluon_I + A/Quark_F
      Split_A = 2d0*alpha_s4Pi * ( -x )
      Split_B = 2d0*alpha_s4Pi * ( (1d0/x-1d0)*(4d0*s12/s15/s25) )
      Split_V(0:3) = Mom(0:3,5) - s15/s12*Mom(0:3,2)

!     color and spin correlated tree process
      SqAmp = Tree_GG_TTb_000_1(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V) *(-4d0/3d0)


      TheDipole%DipoleValue = LeadSing * SqAmp
      return
      END SUBROUTINE




      SUBROUTINE evalDipole4_25x(Mom,Mass2,TheDipole)
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
      TheDipole%Mass2Td(2) = Mass2(5)
      TheDipole%Mass2Td(3) = Mass2(3)
      TheDipole%Mass2Td(4) = Mass2(4)

#IF (_CHECK_DIPOLE_MOMMAP)
      call check_MomTd(2,5,1,TheDipole%MomTd,TheDipole%Mass2Td)
#ENDIF

      if( applyDipoleCut(TheDipole%MomTd).eq. .true. ) then
         TheDipole%DipoleValue = 0d0
         return
      endif

!     splitting 4: Gluon_I --> A/Quark_I + A/Quark_F
      Split_A = 2d0*alpha_s4Pi*(1d0-2d0*x+2d0*x**2)

!     color and spin correlated tree process
      SqAmp = Tree_UUb_TTb_000_2(dcmplx(TheDipole%MomTd),TheDipole%Mass2Td,Split_A,Split_B,Split_V) * (-0.5d0)

      TheDipole%DipoleValue = LeadSing * SqAmp
      return
      END SUBROUTINE










      FUNCTION Tree_GG_TTb_000_1(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      use modKinematics
      implicit none
      double complex SqAmp,Amp(2,0:1),Res(0:3)
      double precision Split_A,Split_B,Split_V(0:3)
      integer col,colP,iHel
      double complex,target :: MomTd(0:3,1:4)
      double precision,target :: Mass2Td(1:4),MassTd(1:4)
      double complex, target :: PolV(0:3,0:1,1:4)
      double precision :: PSWgt2,PSWgt3,MomDK(0:3,1:6)
      type(Particle),target :: TopQuark(1:2)
      integer, parameter :: plus=1, minus=0
      integer,target :: ExtRef=-1
      double complex :: EpsDotP,SpinCorr_plus_plus,SpinCorr_plus_minus
      type(PtrToParticle) :: Gluons(1:2),Quarks(1:2)
      double precision,parameter :: SpinAvg=dble(4), ColAvg=dble(64)
      real(8),parameter::ColCorr(1:2,1:2)=(/64.D0/3.D0,-8.D0/3.D0,-8.D0/3.D0,64.D0/3.D0/)
      integer,target :: HelList(1:4,1:3)
      double complex :: MomTmp(0:3)
      logical :: applyPSCut

      HelList(1,1:3)=(/0,0,0/)
      HelList(2,1:3)=(/1,0,0/)
!       HelList(3,1:3)=(/0,1,0/)
!       HelList(4,1:3)=(/1,1,0/)
      MassTd(1:4) = dsqrt(Mass2Td(1:4))


      call EvalPhasespace_TopDecay(dble(MomTd(0:3,3)),yRndDK(5:8),.false.,MomDK(0:3,4:6),PSWgt3) ! Top
      call EvalPhasespace_TopDecay(dble(MomTd(0:3,4)),yRndDK(1:4),.false.,MomDK(0:3,1:3),PSWgt2) ! ATop


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

      Quarks(1)%Pol => TopQuark(1)%Pol
      Quarks(2)%Pol => TopQuark(2)%Pol


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

      Quarks(1)%Mom => MomTd(:,3)
      Quarks(1)%PartType => Top_
      Quarks(1)%ExtRef => ExtRef
      Quarks(1)%Mass2 => Mass2Td(3)
      Quarks(1)%Mass  => MassTd(3)

      Quarks(2)%Mom => MomTd(:,4)
      Quarks(2)%PartType => ATop_
      Quarks(2)%ExtRef => ExtRef
      Quarks(2)%Mass2 => Mass2Td(4)
      Quarks(2)%Mass  => MassTd(4)


!      sum over helicities
      SqAmp = (0d0,0d0)
      do iHel=1,2
      Gluons(2)%Pol => PolV(:,HelList(iHel,1),2)
      Gluons(2)%Helicity => HelList(iHel,1)
!       Quarks(1)%Pol => PolV(:,HelList(iHel,2),3)
!       Quarks(1)%Helicity => HelList(iHel,2)
!       Quarks(2)%Pol => PolV(:,HelList(iHel,3),4)
!       Quarks(2)%Helicity => HelList(iHel,3)

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
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

!      momentum crossing back
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

      return
      END FUNCTION






      FUNCTION Tree_UUb_TTb_000_2(MomTd,Mass2Td,Split_A,Split_B,Split_V) Result(SqAmp)
      use modProcess
      use modParameters
      use modMyRecurrence
      use modKinematics
      implicit none
      double complex SqAmp,Amp(1,0:1),Res(0:3)
      double precision Split_A,Split_B,Split_V(0:3)
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
      real(8),parameter::ColCorr(1:1,1:1)=(/8.D0/)
      integer,target :: HelList(1:8,1:4)
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
        SqAmp = SqAmp + ColCorr(colP,col) * Split_A*(dconjg(Amp(colP,plus))*Amp(col,plus) + dconjg(Amp(colP,minus))*Amp(col,minus) )
      enddo
      enddo

      enddo
      SqAmp = alpha_s4Pi**2 / SpinAvg / ColAvg * SqAmp

!      momentum crossing back
      MomTd(0:3,1) = -MomTd(0:3,1)
      MomTd(0:3,2) = -MomTd(0:3,2)

      return
      END FUNCTION








      END MODULE
