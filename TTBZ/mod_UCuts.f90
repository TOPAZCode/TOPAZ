MODULE ModUCuts
! this is the version before the matching routines have been moved out of mod_residues

public :: PentCut,QuadCut,TripCut,DoubCut,SingCut

contains


SUBROUTINE PentCut(ThePrimAmp)
use ModProcess
use ModNVBasis
use ModMisc
use ModResidues
implicit none
type(PrimitiveAmplitude),target :: ThePrimAmp
type(UnitarityCut),pointer :: PentCuts
type(Propagator),pointer :: IntPart(:)
type(TreeProcess),pointer :: TreeProcs(:)
integer :: CutNum
double precision :: KMom(1:4,1:4),VMom(1:4,1:4)
double precision :: V5(1:4),x1,x2,x3,x4
double complex :: SqrtPreF
double complex :: LoopMom(1:5),alpha_5
double complex :: Res(1:2)

   include 'misc/global_import'

   PentCuts => ThePrimAmp%UCuts(5)
   IntPart  => ThePrimAmp%IntPart(1:NumExtParticles)
   do CutNum = 1,PentCuts%NumCuts
      if (PentCuts%skip(CutNum)) then
         PentCuts%Coeff(CutNum,:)=(0d0,0d0)
         cycle
      endif
!DEC$ IF (_DebugShowCuts==1)
       print *, "PentCut: ",PentCuts%CutProp(CutNum,1:5)
!DEC$ ENDIF

      TreeProcs => ThePrimAmp%UCuts(5)%TreeProcess(CutNum,1:5)

! def. of k's as in paper
!       KMom(1,1:4) = ThePrimAmp%IntPart(PentCuts%CutProp(CutNum,2))%Mom(1:4) - ThePrimAmp%IntPart(PentCuts%CutProp(CutNum,1))%Mom(1:4)
!       KMom(2,1:4) = ThePrimAmp%IntPart(PentCuts%CutProp(CutNum,3))%Mom(1:4) - ThePrimAmp%IntPart(PentCuts%CutProp(CutNum,2))%Mom(1:4)
!       KMom(3,1:4) = ThePrimAmp%IntPart(PentCuts%CutProp(CutNum,4))%Mom(1:4) - ThePrimAmp%IntPart(PentCuts%CutProp(CutNum,3))%Mom(1:4)
!       KMom(4,1:4) = ThePrimAmp%IntPart(PentCuts%CutProp(CutNum,5))%Mom(1:4) - ThePrimAmp%IntPart(PentCuts%CutProp(CutNum,4))%Mom(1:4)
! def. of k's as in kirill's program
      KMom(1,1:4) = IntPart(PentCuts%CutProp(CutNum,2))%Mom(1:4) - IntPart(PentCuts%CutProp(CutNum,1))%Mom(1:4)
      KMom(2,1:4) = IntPart(PentCuts%CutProp(CutNum,3))%Mom(1:4) - IntPart(PentCuts%CutProp(CutNum,1))%Mom(1:4)
      KMom(3,1:4) = IntPart(PentCuts%CutProp(CutNum,4))%Mom(1:4) - IntPart(PentCuts%CutProp(CutNum,1))%Mom(1:4)
      KMom(4,1:4) = IntPart(PentCuts%CutProp(CutNum,5))%Mom(1:4) - IntPart(PentCuts%CutProp(CutNum,1))%Mom(1:4)


! construct NV-basis
      call NVBasis4(KMom(1:4,1:4),VMom(1:4,1:4))

      x1 = (KMom(1,1:4).dot.KMom(1,1:4)) + IntPart(PentCuts%CutProp(CutNum,1))%Mass2 - IntPart(PentCuts%CutProp(CutNum,2))%Mass2
      x2 = (KMom(2,1:4).dot.KMom(2,1:4)) + IntPart(PentCuts%CutProp(CutNum,1))%Mass2 - IntPart(PentCuts%CutProp(CutNum,3))%Mass2
      x3 = (KMom(3,1:4).dot.KMom(3,1:4)) + IntPart(PentCuts%CutProp(CutNum,1))%Mass2 - IntPart(PentCuts%CutProp(CutNum,4))%Mass2
      x4 = (KMom(4,1:4).dot.KMom(4,1:4)) + IntPart(PentCuts%CutProp(CutNum,1))%Mass2 - IntPart(PentCuts%CutProp(CutNum,5))%Mass2
      V5(1:4) = -0.5d0*( x1*VMom(1,1:4) + x2*VMom(2,1:4) + x3*VMom(3,1:4) + x4*VMom(4,1:4) )
      SqrtPreF = cdsqrt(dcmplx( IntPart(PentCuts%CutProp(CutNum,1))%Mass2 - (V5.dot.V5) ))

      LoopMom(1:4) = V5(1:4)
      alpha_5 = (0d0,+1d0)
      LoopMom(5)   = SqrtPreF * alpha_5
      call resid5(LoopMom(1:5),dcmplx(KMom(1,1:4)),dcmplx(KMom(2,1:4)),dcmplx(KMom(3,1:4)),dcmplx(KMom(4,1:4)),PentCuts%CutProp(CutNum,1:5),TreeProcs,Res(1))
      alpha_5 = (0d0,-1d0)
      LoopMom(5)   = SqrtPreF * alpha_5
      call resid5(LoopMom(1:5),dcmplx(KMom(1,1:4)),dcmplx(KMom(2,1:4)),dcmplx(KMom(3,1:4)),dcmplx(KMom(4,1:4)),PentCuts%CutProp(CutNum,1:5),TreeProcs,Res(2))

      PentCuts%Coeff(CutNum,0) = 0.5d0*(Res(1) + Res(2))/LoopMom(5)**2

!     set vars. for kirill's routines
      coeff5(CutNum,1)= PentCuts%Coeff(CutNum,0)
      mass5(CutNum,1) = IntPart( PentCuts%CutProp(CutNum,1) )%Mass
      mass5(CutNum,2) = IntPart( PentCuts%CutProp(CutNum,2) )%Mass
      mass5(CutNum,3) = IntPart( PentCuts%CutProp(CutNum,3) )%Mass
      mass5(CutNum,4) = IntPart( PentCuts%CutProp(CutNum,4) )%Mass
      mass5(CutNum,5) = IntPart( PentCuts%CutProp(CutNum,5) )%Mass
      propv5(CutNum,1:4)  = (0d0,0d0)
      propv5(CutNum,5:8)  = KMom(1,1:4)
      propv5(CutNum,9:12) = KMom(2,1:4)
      propv5(CutNum,13:16)= KMom(3,1:4)
      propv5(CutNum,17:20)= KMom(4,1:4)


!DEC$ IF (_DebugPrintPentCoeff==1)
       print *, coeff5(CutNum,1)
      enddo
      pause
!DEC$ ELSE
      enddo
!DEC$ ENDIF


END SUBROUTINE






SUBROUTINE QuadCut(ThePrimAmp)
use ModProcess
use ModNVBasis
use ModMisc
use ModResidues
use ifport
implicit none
type(PrimitiveAmplitude),target :: ThePrimAmp
type(UnitarityCut),pointer :: QuadCuts
type(Propagator) :: IntPart(NumExtParticles)
type(TreeProcess),pointer :: TreeProcs(:)
integer :: CutNum,k
double precision :: KMom(1:3,1:4),VMom(1:3,1:4)
double precision :: V4(1:4),x1,x2,x3
double complex :: SqrtPreF
double complex :: NMom(1,1:4),LoopMom(1:5),alpha_1,alpha_5
double complex :: Res(1:5),ResCheck
double complex :: LinSysEq(1:3,1:4),s1(1:5),se2(3:5)

   include 'misc/global_import'


   QuadCuts => ThePrimAmp%UCuts(4)
   IntPart(1:NumExtParticles) = ThePrimAmp%IntPart(1:NumExtParticles)

   do CutNum = 1,QuadCuts%NumCuts
      if (QuadCuts%skip(CutNum)) then
         QuadCuts%Coeff(CutNum,:)=(0d0,0d0)
         cycle
      endif

!DEC$ IF (_DebugShowCuts==1)
   print *, "QuadCut",CutNum,": ", QuadCuts%CutProp(CutNum,1:4)
!DEC$ ENDIF

      TreeProcs => ThePrimAmp%UCuts(4)%TreeProcess(CutNum,1:4)

! def. of k's as in paper
!       KMom(1,1:4) = IntPart(QuadCuts%CutProp(CutNum,2))%Mom(1:4) - IntPart(QuadCuts%CutProp(CutNum,1))%Mom(1:4)
!       KMom(2,1:4) = IntPart(QuadCuts%CutProp(CutNum,3))%Mom(1:4) - IntPart(QuadCuts%CutProp(CutNum,2))%Mom(1:4)
!       KMom(3,1:4) = IntPart(QuadCuts%CutProp(CutNum,4))%Mom(1:4) - IntPart(QuadCuts%CutProp(CutNum,3))%Mom(1:4)
! def. of k's as in kirill's program
      KMom(1,1:4) = IntPart(QuadCuts%CutProp(CutNum,2))%Mom(1:4) - IntPart(QuadCuts%CutProp(CutNum,1))%Mom(1:4)
      KMom(2,1:4) = IntPart(QuadCuts%CutProp(CutNum,3))%Mom(1:4) - IntPart(QuadCuts%CutProp(CutNum,1))%Mom(1:4)
      KMom(3,1:4) = IntPart(QuadCuts%CutProp(CutNum,4))%Mom(1:4) - IntPart(QuadCuts%CutProp(CutNum,1))%Mom(1:4)

! construct NV-basis
      call NVBasis3(KMom(1:3,1:4),VMom(1:3,1:4),NMom(1,1:4))
      x1 = (KMom(1,1:4).dot.KMom(1,1:4)) + IntPart(QuadCuts%CutProp(CutNum,1))%Mass2 - IntPart(QuadCuts%CutProp(CutNum,2))%Mass2
      x2 = (KMom(2,1:4).dot.KMom(2,1:4)) + IntPart(QuadCuts%CutProp(CutNum,1))%Mass2 - IntPart(QuadCuts%CutProp(CutNum,3))%Mass2
      x3 = (KMom(3,1:4).dot.KMom(3,1:4)) + IntPart(QuadCuts%CutProp(CutNum,1))%Mass2 - IntPart(QuadCuts%CutProp(CutNum,4))%Mass2
      V4(1:4) = -0.5d0*( x1*VMom(1,1:4) + x2*VMom(2,1:4) + x3*VMom(3,1:4) )
      SqrtPreF = cdsqrt(dcmplx( IntPart(QuadCuts%CutProp(CutNum,1))%Mass2 - (V4.dot.V4) ))


!---------------- 1 ----------------------
      alpha_1 = +1d0
!     alpha_5 = 0d0
      s1(1) = alpha_1*SqrtPreF
      LoopMom(1:4) = V4(1:4) + s1(1)*NMom(1,1:4)
      LoopMom(5)   = 0d0
      call resid4(LoopMom(1:5),dcmplx(KMom(1,1:4)),dcmplx(KMom(2,1:4)),dcmplx(KMom(3,1:4)),QuadCuts%CutProp(CutNum,1:4),TreeProcs,Res(1))
!----------------- 2 ---------------------
      alpha_1 = -1d0
!     alpha_5 = 0d0
      s1(2) = alpha_1*SqrtPreF
      LoopMom(1:4) = V4(1:4) + s1(2)*NMom(1,1:4)
      LoopMom(5)   = 0d0

      call resid4(LoopMom(1:5),dcmplx(KMom(1,1:4)),dcmplx(KMom(2,1:4)),dcmplx(KMom(3,1:4)),QuadCuts%CutProp(CutNum,1:4),TreeProcs,Res(2))
      QuadCuts%Coeff(CutNum,0) = 0.5d0*(Res(1) + Res(2))
      QuadCuts%Coeff(CutNum,1) = 0.5d0*(Res(1) - Res(2))/SqrtPreF
!----------------- 3 ---------------------


      alpha_1 = drand(1)
      do k=3,5
            alpha_1 = drand(0)
            alpha_5 = drand(0)
            s1(k)  =   alpha_1*SqrtPreF/cdsqrt(alpha_1**2+alpha_5**2)
            se2(k) = ( alpha_5*SqrtPreF/cdsqrt(alpha_1**2+alpha_5**2) )**2
            se2(k) = ( alpha_5*SqrtPreF )**2 /(alpha_1**2+alpha_5**2)
            LoopMom(1:4)  = V4(1:4) + s1(k)*NMom(1,1:4)
            LoopMom(5)    = alpha_5*SqrtPreF/cdsqrt(alpha_1**2+alpha_5**2) * (0d0,1d0)     !   n_5=(0,0,0,0,i), n_5.dot.n_5=-1
            LinSysEq(k-2,1) = se2(k)
            LinSysEq(k-2,2) = se2(k)*s1(k)
            LinSysEq(k-2,3) = se2(k)**2
            call resid4(LoopMom(1:5),dcmplx(KMom(1,1:4)),dcmplx(KMom(2,1:4)),dcmplx(KMom(3,1:4)),QuadCuts%CutProp(CutNum,1:4),TreeProcs,Res(k))
            LinSysEq(k-2,4) = Res(k) - QuadCuts%Coeff(CutNum,0) - QuadCuts%Coeff(CutNum,1) * s1(k)
      enddo
      QuadCuts%Coeff(CutNum,2:4) = go_Gauss(3,LinSysEq(1:3,1:4))

!     set vars. for kirill's routines
      coeff4(CutNum,1:5) = QuadCuts%Coeff(CutNum,0:4)
      refvect4(CutNum,1:4) = NMom(1,1:4)
      mass4(CutNum,1) = IntPart( QuadCuts%CutProp(CutNum,1) )%Mass
      mass4(CutNum,2) = IntPart( QuadCuts%CutProp(CutNum,2) )%Mass
      mass4(CutNum,3) = IntPart( QuadCuts%CutProp(CutNum,3) )%Mass
      mass4(CutNum,4) = IntPart( QuadCuts%CutProp(CutNum,4) )%Mass
      propv4(CutNum,1:4)  = (0d0,0d0)
      propv4(CutNum,5:8)  = KMom(1,1:4)
      propv4(CutNum,9:12) = KMom(2,1:4)
      propv4(CutNum,13:16)= KMom(3,1:4)




!DEC$ IF (_ResidueCheckQuad==1)
      ResCheck = Res(1) - (QuadCuts%Coeff(CutNum,0)+ QuadCuts%Coeff(CutNum,1)*s1(1))
      print *, "QuadRes 1",ResCheck
      ResCheck = Res(2) - (QuadCuts%Coeff(CutNum,0)+ QuadCuts%Coeff(CutNum,1)*s1(2))
      print *, "QuadRes 2",ResCheck
      ResCheck = Res(3) -  &
      (QuadCuts%Coeff(CutNum,0)+QuadCuts%Coeff(CutNum,1)*s1(3)+(QuadCuts%Coeff(CutNum,2)+QuadCuts%Coeff(CutNum,3)*s1(3))*se2(3)+QuadCuts%Coeff(CutNum,4)*se2(3)**2)
      print *, "QuadRes 3",ResCheck
      ResCheck = Res(4) -  &
      (QuadCuts%Coeff(CutNum,0)+QuadCuts%Coeff(CutNum,1)*s1(4)+(QuadCuts%Coeff(CutNum,2)+QuadCuts%Coeff(CutNum,3)*s1(4))*se2(4)+QuadCuts%Coeff(CutNum,4)*se2(4)**2)
      print *, "QuadRes 4",ResCheck
      ResCheck = Res(5) -  &
      (QuadCuts%Coeff(CutNum,0)+QuadCuts%Coeff(CutNum,1)*s1(5)+(QuadCuts%Coeff(CutNum,2)+QuadCuts%Coeff(CutNum,3)*s1(5))*se2(5)+QuadCuts%Coeff(CutNum,4)*se2(5)**2)
      print *, "QuadRes 5",ResCheck
!DEC$ ENDIF


!DEC$ IF (_DebugPrintQuadCoeff==1)
      print *, ""
      do k=1,5
         print *, k,coeff4(CutNum,k)
      enddo
      enddo
      pause
!DEC$ ELSE
      enddo
!DEC$ ENDIF


END SUBROUTINE









SUBROUTINE TripCut(ThePrimAmp)
use ModProcess
use ModParameters
use ModNVBasis
use ModMisc
use ModResidues
use ifport
implicit none
type(PrimitiveAmplitude),target :: ThePrimAmp
type(UnitarityCut),pointer :: TripCuts
type(Propagator) :: IntPart(NumExtParticles)
type(TreeProcess),pointer :: TreeProcs(:)
integer :: CutNum
double precision :: KMom(1:2,1:4),VMom(1:2,1:4)
double precision :: V3(1:4),x1,x2,x3,phi
double complex :: SqrtPreF,swap_tmp
double complex :: NMom(1:2,1:4),LoopMom(1:5),alpha_1,alpha_2,alpha_5
double complex :: Res8,Res9,Res10,ResCheck,s1(1:10),s2(1:10),se2(8:10)
double complex :: LinSysEq(1:8,1:8)
integer :: k,kk
double complex :: LHSEQ(1:4),LHSEQ_NEW(1:4),BCOEFF(0:3)
include "misc/global_import"


   TripCuts => ThePrimAmp%UCuts(3)
   IntPart(1:NumExtParticles) = ThePrimAmp%IntPart(1:NumExtParticles)

    do CutNum = 1,TripCuts%NumCuts
      if (TripCuts%skip(CutNum)) then
         TripCuts%Coeff(CutNum,:)=(0d0,0d0)
         cycle
      endif
!DEC$ IF (_DebugShowCuts==1)
   print *, "TripCut",CutNum,":",TripCuts%CutProp(CutNum,1:3)
!DEC$ ENDIF

      TreeProcs => ThePrimAmp%UCuts(3)%TreeProcess(CutNum,1:3)

! def. of k's as in paper
!       KMom(1,1:4) = IntPart(TripCuts%CutProp(CutNum,2))%Mom(1:4) - IntPart(TripCuts%CutProp(CutNum,1))%Mom(1:4)
!       KMom(2,1:4) = IntPart(TripCuts%CutProp(CutNum,3))%Mom(1:4) - IntPart(TripCuts%CutProp(CutNum,2))%Mom(1:4)
! def. of k's as in kirill's program
      KMom(1,1:4) = IntPart(TripCuts%CutProp(CutNum,2))%Mom(1:4) - IntPart(TripCuts%CutProp(CutNum,1))%Mom(1:4)
      KMom(2,1:4) = IntPart(TripCuts%CutProp(CutNum,3))%Mom(1:4) - IntPart(TripCuts%CutProp(CutNum,1))%Mom(1:4)

! construct NV-basis
      call NVBasis2(KMom(1:2,1:4),VMom(1:2,1:4),NMom(1:2,1:4))
      x1 = (KMom(1,1:4).dot.KMom(1,1:4)) + IntPart(TripCuts%CutProp(CutNum,1))%Mass2 - IntPart(TripCuts%CutProp(CutNum,2))%Mass2
      x2 = (KMom(2,1:4).dot.KMom(2,1:4)) + IntPart(TripCuts%CutProp(CutNum,1))%Mass2 - IntPart(TripCuts%CutProp(CutNum,3))%Mass2
      V3(1:4) = -0.5d0*( x1*VMom(1,1:4) + x2*VMom(2,1:4) )
      SqrtPreF = cdsqrt(dcmplx( IntPart(TripCuts%CutProp(CutNum,1))%Mass2 - (V3.dot.V3) ))


      alpha_1 = drand(12)
      if (cdabs(SqrtPreF).gt. 1d-1) then
         do k=1,7
            alpha_1 = drand(0)
            alpha_2 = drand(0)
            s1(k) = alpha_1*SqrtPreF/cdsqrt(alpha_1**2+alpha_2**2)
            s2(k) = alpha_2*SqrtPreF/cdsqrt(alpha_1**2+alpha_2**2)
            LoopMom(1:4) = V3(1:4) + s1(k)*NMom(1,1:4) + s2(k)*NMom(2,1:4)
            LoopMom(5)   = 0d0
!             LinSysEq(k,1) = 1d0
!             LinSysEq(k,2) = s1(k)
!             LinSysEq(k,3) = s2(k)
!             LinSysEq(k,4) = s1(k)**2 - s2(k)**2
!             LinSysEq(k,5) = s1(k)*s2(k)
!             LinSysEq(k,6) = s1(k)**2*s2(k)
!             LinSysEq(k,7) = s1(k)*s2(k)**2
            LinSysEq(k,1) = 1d0
            LinSysEq(k,2) = s1(k)
            LinSysEq(k,3) = s2(k)
            LinSysEq(k,4) = s1(k)*s2(k)
            LinSysEq(k,5) = s1(k)**2 - s2(k)**2
            LinSysEq(k,6) = s1(k)**2*s2(k)
            LinSysEq(k,7) = s1(k)*s2(k)**2
            call resid3(LoopMom(1:5),dcmplx(KMom(1,1:4)),dcmplx(KMom(2,1:4)),TripCuts%CutProp(CutNum,1:3),TreeProcs,LinSysEq(k,8))
         enddo

         TripCuts%Coeff(CutNum,0:6) = go_Gauss(7,LinSysEq(1:7,1:8))


         alpha_1 = +dsqrt(0.5d0)
         alpha_5 = +dsqrt(0.5d0)
         s1(8) = alpha_1*SqrtPreF
   !     s2(8) = 0d0
         se2(8) = (alpha_5*SqrtPreF)**2
         LoopMom(1:4) = V3(1:4) + s1(8)*NMom(1,1:4)
         LoopMom(5)   = alpha_5 * SqrtPreF * (0d0,1d0)
         LinSysEq(1,1) = se2(8)
         LinSysEq(1,2) = se2(8)*s1(8)
         LinSysEq(1,3) = 0d0
         call resid3(LoopMom(1:5),dcmplx(KMom(1,1:4)),dcmplx(KMom(2,1:4)),TripCuts%CutProp(CutNum,1:3),TreeProcs,Res8)
         LinSysEq(1,4) = Res8 - TripCuts%Coeff(CutNum,0)           &
                                - TripCuts%Coeff(CutNum,1) * s1(8)     &
!                                 - TripCuts%Coeff(CutNum,3) * s1(8)**2
                                - TripCuts%Coeff(CutNum,4) * s1(8)**2

         alpha_1 = -dsqrt(0.5d0)
         alpha_5 = +dsqrt(0.5d0)
         s1(9) = alpha_1*SqrtPreF
   !     s2(9) = 0d0
         se2(9) = (alpha_5*SqrtPreF)**2
         LoopMom(1:4) = V3(1:4) + s1(9)*NMom(1,1:4)
         LoopMom(5)   = alpha_5 * SqrtPreF * (0d0,1d0)
         LinSysEq(2,1) = se2(9)
         LinSysEq(2,2) = se2(9)*s1(9)
         LinSysEq(2,3) = 0d0
         call resid3(LoopMom(1:5),dcmplx(KMom(1,1:4)),dcmplx(KMom(2,1:4)),TripCuts%CutProp(CutNum,1:3),TreeProcs,Res9)
         LinSysEq(2,4) = Res9 - TripCuts%Coeff(CutNum,0)             &
                                - TripCuts%Coeff(CutNum,1) * s1(9)     &
!                                 - TripCuts%Coeff(CutNum,3) * s1(9)**2
                                - TripCuts%Coeff(CutNum,4) * s1(9)**2



         alpha_2 = +dsqrt(0.5d0)
         alpha_5 = +dsqrt(0.5d0)
   !     s1(10) = 0d0
         s2(10) = alpha_2*SqrtPreF
         se2(10) = (alpha_5*SqrtPreF)**2
         LoopMom(1:4) = V3(1:4) + s2(10)*NMom(2,1:4)
         LoopMom(5)   = alpha_5 * SqrtPreF * (0d0,1d0)
         LinSysEq(3,1) = se2(10)
         LinSysEq(3,2) = 0d0
         LinSysEq(3,3) = se2(10)*s2(10)
         call resid3(LoopMom(1:5),dcmplx(KMom(1,1:4)),dcmplx(KMom(2,1:4)),TripCuts%CutProp(CutNum,1:3),TreeProcs,Res10)
         LinSysEq(3,4) = Res10 - TripCuts%Coeff(CutNum,0)              &
                                 - TripCuts%Coeff(CutNum,2) * s2(10)     &
!                                  + TripCuts%Coeff(CutNum,3) * s2(10)**2
                                 + TripCuts%Coeff(CutNum,4) * s2(10)**2

         TripCuts%Coeff(CutNum,7:9) = go_Gauss(3,LinSysEq(1:3,1:4))

      else!if( cdabs(SqrtPreF).lt.1d-1 )

         alpha_1 = drand(7)
         do k=1,4
!             phi = 2d0*DblPi*(k-1d0)/4d0
!             alpha_1 = + cdexp( (0d0,1d0)*phi ) + (0d0,1d-12)
            alpha_1 = drand(0)
            alpha_2 = + cdsqrt(SqrtPreF**2-alpha_1**2)
            s1(k) = alpha_1
            s2(k) = alpha_2
            LoopMom(1:4) = V3(1:4) + s1(k)*NMom(1,1:4) + s2(k)*NMom(2,1:4)
            LoopMom(5)   = 0d0

            LinSysEq(k,1) = 1d0
            LinSysEq(k,2) = s1(k)
            LinSysEq(k,7) = s2(k)
            LinSysEq(k,4) = s1(k)*s2(k)
            LinSysEq(k,5) = s1(k)**2 - s2(k)**2
            LinSysEq(k,6) = s1(k)**2*s2(k)
            LinSysEq(k,3) = s1(k)*s2(k)**2
            call resid3(LoopMom(1:5),dcmplx(KMom(1,1:4)),dcmplx(KMom(2,1:4)),TripCuts%CutProp(CutNum,1:3),TreeProcs,LinSysEq(k,8))
         enddo



         do k=5,7
!             phi = 2d0*DblPi*(k-4d0-1d0)/4d0
!             alpha_1 = + cdexp( (0d0,1d0)*phi ) + (0d0,1d-12)
            alpha_1 = drand(0)
            alpha_2 = - cdsqrt(SqrtPreF**2-alpha_1**2)
            s1(k) = alpha_1
            s2(k) = alpha_2
            LoopMom(1:4) = V3(1:4) + s1(k)*NMom(1,1:4) + s2(k)*NMom(2,1:4)
            LoopMom(5)   = 0d0

            LinSysEq(k,1) = 1d0
            LinSysEq(k,2) = s1(k)
            LinSysEq(k,7) = s2(k)
            LinSysEq(k,4) = s1(k)*s2(k)
            LinSysEq(k,5) = s1(k)**2 - s2(k)**2
            LinSysEq(k,6) = s1(k)**2*s2(k)
            LinSysEq(k,3) = s1(k)*s2(k)**2
            call resid3(LoopMom(1:5),dcmplx(KMom(1,1:4)),dcmplx(KMom(2,1:4)),TripCuts%CutProp(CutNum,1:3),TreeProcs,LinSysEq(k,8))
         enddo

         TripCuts%Coeff(CutNum,0:6) = go_Gauss(7,LinSysEq(1:7,1:8))
!          print *, "check"
!          print *, check_MatSol(7,LinSysEq(1:7,1:7),TripCuts%Coeff(CutNum,0:6),LinSysEq(1:7,8))
!          print *, "check"
!        swap 3 <--> 7 due to part. pivoting
         swap_tmp = TripCuts%Coeff(CutNum,2)
         TripCuts%Coeff(CutNum,2) = TripCuts%Coeff(CutNum,6)
         TripCuts%Coeff(CutNum,6) = swap_tmp

!          do k=1,4
!             phi = 2d0*DblPi*(k-1d0)/4d0
!             alpha_1 = + cdexp( (0d0,1d0)*phi) + (0d0,1d-12)
!             alpha_2 = + cdsqrt(SqrtPreF**2-alpha_1**2)
!             s1(k) = alpha_1
!             s2(k) = alpha_2
!             LoopMom(1:4) = V3(1:4) + s1(k)*NMom(1,1:4) + s2(k)*NMom(2,1:4)
!             LoopMom(5)   = 0d0
!
!             LinSysEq(k,1) = 1d0
!             LinSysEq(k,2) = s1(k)
!             LinSysEq(k,7) = s2(k)
!             LinSysEq(k,4) = s1(k)*s2(k)
!             LinSysEq(k,5) = s1(k)**2 - s2(k)**2
!             LinSysEq(k,6) = s1(k)**2*s2(k)
!             LinSysEq(k,3) = s1(k)*s2(k)**2
!             call resid3(LoopMom(1:5),dcmplx(KMom(1,1:4)),dcmplx(KMom(2,1:4)),TripCuts%CutProp(CutNum,1:3),TreeProcs,LinSysEq(k,8))
!          enddo
!
!          do k=5,8
!             phi = 2d0*DblPi*(k-4d0-1d0)/4d0
!             alpha_1 = + cdexp( (0d0,1d0)*phi ) + (0d0,1d-12)
!             alpha_2 = - cdsqrt(SqrtPreF**2-alpha_1**2)
!             s1(k) = alpha_1
!             s2(k) = alpha_2
!             LoopMom(1:4) = V3(1:4) + s1(k)*NMom(1,1:4) + s2(k)*NMom(2,1:4)
!             LoopMom(5)   = 0d0
!
!             LinSysEq(k,1) = 1d0
!             LinSysEq(k,2) = s1(k)
!             LinSysEq(k,7) = s2(k)
!             LinSysEq(k,4) = s1(k)*s2(k)
!             LinSysEq(k,5) = s1(k)**2 - s2(k)**2
!             LinSysEq(k,6) = s1(k)**2*s2(k)
!             LinSysEq(k,3) = s1(k)*s2(k)**2
!             call resid3(LoopMom(1:5),dcmplx(KMom(1,1:4)),dcmplx(KMom(2,1:4)),TripCuts%CutProp(CutNum,1:3),TreeProcs,LinSysEq(k,8))
!
!             lhseq_new(k-4)=0.5d0/(-alpha_2)*(LinSysEq(k-4,8)-LinSysEq(k,8))
!          enddo
!          lhseq(1)=0.5d0*(LinSysEq(1,8)+LinSysEq(5,8))
!          lhseq(2)=0.5d0*(LinSysEq(2,8)+LinSysEq(6,8))
!          lhseq(3)=0.5d0*(LinSysEq(3,8)+LinSysEq(7,8))
!          lhseq(4)=0.5d0*(LinSysEq(4,8)+LinSysEq(8,8))
!          do k=1,4
!             bcoeff(k-1) = zero
!             do kk=1,4
!                 bcoeff(k-1)=bcoeff(k-1) + 0.25d0*lhseq(kk)*cdexp( -2d0*pi*dcmplx(0d0,1d0)/4d0*(kk-1)*(k-1) )
!             enddo
!          enddo
!         TripCuts%Coeff(CutNum,6)=-bcoeff(3)
!         TripCuts%Coeff(CutNum,4)=0.5d0*bcoeff(2)
!         TripCuts%Coeff(CutNum,1)=bcoeff(1) - TripCuts%Coeff(CutNum,6)*SqrtPreF**2
!         TripCuts%Coeff(CutNum,0)=bcoeff(0) + TripCuts%Coeff(CutNum,4)*SqrtPreF**2
!         TripCuts%Coeff(CutNum,3)=(lhseq_new(1)-lhseq_new(3))/2d0
!         TripCuts%Coeff(CutNum,2)=(lhseq_new(1)+lhseq_new(2)-dcmplx(1d0,1d0)*TripCuts%Coeff(CutNum,3))/2d0
!         TripCuts%Coeff(CutNum,5)=(lhseq_new(1)-lhseq_new(2)-dcmplx(1d0,-1d0)*TripCuts%Coeff(CutNum,3))/2d0






         alpha_1 = (1d0, -1d-12)   !!!!    THIS VIOLATES on-shellness of loop-mom by 10^-12   --> gauge cancellation only valid up to 10^-12
         alpha_5 = cdsqrt( SqrtPreF**2 - alpha_1**2 )
         s1(8) = alpha_1
         s2(8) = 0d0
         se2(8) = alpha_5**2
         LoopMom(1:4) = V3(1:4) + s1(8)*NMom(1,1:4)
         LoopMom(5)   = alpha_5 * (0d0,1d0)
         LinSysEq(1,1) = se2(8)
         LinSysEq(1,2) = se2(8)*s1(8)
         LinSysEq(1,3) = 0d0
         call resid3(LoopMom(1:5),dcmplx(KMom(1,1:4)),dcmplx(KMom(2,1:4)),TripCuts%CutProp(CutNum,1:3),TreeProcs,Res8)
         LinSysEq(1,4) = Res8 - TripCuts%Coeff(CutNum,0)             &
                              - TripCuts%Coeff(CutNum,1) * s1(8)     &
                              - TripCuts%Coeff(CutNum,4) * s1(8)**2

         alpha_1 = (-1d0, 0d0)
         alpha_5 = alpha_5
         s1(9) = alpha_1
         s2(9) = 0d0
         se2(9) = alpha_5**2
         LoopMom(1:4) = V3(1:4) + s1(9)*NMom(1,1:4)
         LoopMom(5)   = alpha_5 * (0d0,1d0)
         LinSysEq(2,1) = se2(9)
         LinSysEq(2,2) = se2(9)*s1(9)
         LinSysEq(2,3) = 0d0

         call resid3(LoopMom(1:5),dcmplx(KMom(1,1:4)),dcmplx(KMom(2,1:4)),TripCuts%CutProp(CutNum,1:3),TreeProcs,Res9)
         LinSysEq(2,4) = Res9 - TripCuts%Coeff(CutNum,0)             &
                                - TripCuts%Coeff(CutNum,1) * s1(9)     &
                                - TripCuts%Coeff(CutNum,4) * s1(9)**2


         alpha_2 = (1d0, 1d-12)
         alpha_5 = cdsqrt( SqrtPreF**2 - alpha_2**2 )
         s1(10) = 0d0
         s2(10) = alpha_2
         se2(10) = alpha_5**2
         LoopMom(1:4) = V3(1:4) + s2(10)*NMom(2,1:4)
         LoopMom(5)   = alpha_5 * (0d0,1d0)
         LinSysEq(3,1) = se2(10)
         LinSysEq(3,2) = 0d0
         LinSysEq(3,3) = se2(10)*s2(10)
         call resid3(LoopMom(1:5),dcmplx(KMom(1,1:4)),dcmplx(KMom(2,1:4)),TripCuts%CutProp(CutNum,1:3),TreeProcs,Res10)
         LinSysEq(3,4) = Res10 - TripCuts%Coeff(CutNum,0)             &
                                 - TripCuts%Coeff(CutNum,2) * s2(10)    &
                                 + TripCuts%Coeff(CutNum,4) * s2(10)**2

         TripCuts%Coeff(CutNum,7:9) = go_Gauss(3,LinSysEq(1:3,1:4))

!________________________________________________
      endif



!DEC$ IF (_ResidueCheckTrip==1)
      do k=1,7
      ResCheck = LinSysEq(k,8) - (TripCuts%Coeff(CutNum,0)+TripCuts%Coeff(CutNum,1)*s1(k)+TripCuts%Coeff(CutNum,2)*s2(k)+TripCuts%Coeff(CutNum,3)*s1(k)*s2(k)  &
                                 +TripCuts%Coeff(CutNum,4)*(s1(k)**2-s2(k)**2)+(TripCuts%Coeff(CutNum,5)*s1(k)+TripCuts%Coeff(CutNum,6)*s2(k))*s1(k)*s2(k) )
      print *, "TripRes",k,ResCheck
      enddo
      k=8
      ResCheck = Res8 - (TripCuts%Coeff(CutNum,0)+TripCuts%Coeff(CutNum,1)*s1(k)+TripCuts%Coeff(CutNum,2)*s2(k)+TripCuts%Coeff(CutNum,3)*s1(k)*s2(k)  &
                                 +TripCuts%Coeff(CutNum,4)*(s1(k)**2-s2(k)**2)+(TripCuts%Coeff(CutNum,5)*s1(k)+TripCuts%Coeff(CutNum,6)*s2(k))*s1(k)*s2(k)     &
                                 +TripCuts%Coeff(CutNum,7)*se2(k)+TripCuts%Coeff(CutNum,8)*se2(k)*s1(k)+TripCuts%Coeff(CutNum,9)*se2(k)*s2(k))
      print *, "TripRes",k,ResCheck
      k=9
      ResCheck = Res9 - (TripCuts%Coeff(CutNum,0)+TripCuts%Coeff(CutNum,1)*s1(k)+TripCuts%Coeff(CutNum,2)*s2(k)+TripCuts%Coeff(CutNum,3)*s1(k)*s2(k)  &
                                 +TripCuts%Coeff(CutNum,4)*(s1(k)**2-s2(k)**2)+(TripCuts%Coeff(CutNum,5)*s1(k)+TripCuts%Coeff(CutNum,6)*s2(k))*s1(k)*s2(k)     &
                                 +TripCuts%Coeff(CutNum,7)*se2(k)+TripCuts%Coeff(CutNum,8)*se2(k)*s1(k)+TripCuts%Coeff(CutNum,9)*se2(k)*s2(k))
      print *, "TripRes",k,ResCheck
      k=10
      ResCheck = Res10 - (TripCuts%Coeff(CutNum,0)+TripCuts%Coeff(CutNum,1)*s1(k)+TripCuts%Coeff(CutNum,2)*s2(k)+TripCuts%Coeff(CutNum,3)*s1(k)*s2(k)  &
                                 +TripCuts%Coeff(CutNum,4)*(s1(k)**2-s2(k)**2)+(TripCuts%Coeff(CutNum,5)*s1(k)+TripCuts%Coeff(CutNum,6)*s2(k))*s1(k)*s2(k)     &
                                 +TripCuts%Coeff(CutNum,7)*se2(k)+TripCuts%Coeff(CutNum,8)*se2(k)*s1(k)+TripCuts%Coeff(CutNum,9)*se2(k)*s2(k))
      print *, "TripRes",k,ResCheck
!DEC$ ENDIF



!     set vars. for kirill's routines
      coeff3(CutNum,1) = TripCuts%Coeff(CutNum,0)
      coeff3(CutNum,2) = TripCuts%Coeff(CutNum,1)
      coeff3(CutNum,3) = TripCuts%Coeff(CutNum,2)
      coeff3(CutNum,4) = TripCuts%Coeff(CutNum,3)
      coeff3(CutNum,5) = TripCuts%Coeff(CutNum,4)
      coeff3(CutNum,6) = TripCuts%Coeff(CutNum,5)
      coeff3(CutNum,7) = TripCuts%Coeff(CutNum,6)
      coeff3(CutNum,8) = TripCuts%Coeff(CutNum,7)
      coeff3(CutNum,9) = TripCuts%Coeff(CutNum,8)
      coeff3(CutNum,10)= TripCuts%Coeff(CutNum,9)

      refvect3(CutNum,1:4) = NMom(1,1:4)
      refvect3(CutNum,5:8) = NMom(2,1:4)
      mass3(CutNum,1) = IntPart( TripCuts%CutProp(CutNum,1) )%Mass
      mass3(CutNum,2) = IntPart( TripCuts%CutProp(CutNum,2) )%Mass
      mass3(CutNum,3) = IntPart( TripCuts%CutProp(CutNum,3) )%Mass
      propv3(CutNum,1:4)  = (0d0,0d0)
      propv3(CutNum,5:8)  = KMom(1,1:4)
      propv3(CutNum,9:12) = KMom(2,1:4)

!       print *, "# tripcut",CutNum
!       do k=1,10
!         write(*,"(1PE23.16,3X,1PE23.16)") dreal(coeff3(CutNum,k)),dimag(coeff3(CutNum,k))
!       enddo
!       do k=1,10
!       if( cdabs(coeff3(CutNum,k)).lt.1d-6 ) coeff3(CutNum,k) = (0d0,0d0)
!       enddo

!DEC$ IF (_DebugPrintTripCoeff==1)
      do k=1,10
         print *, coeff3(Cutnum,k)
      enddo
      enddo
      pause
!DEC$ ELSE
      enddo
!DEC$ ENDIF


END SUBROUTINE





SUBROUTINE DoubCut(ThePrimAmp)
use ModProcess
use ModParameters
use ModNVBasis
use ModMisc
use ModResidues
use ifport
implicit none
type(PrimitiveAmplitude),target :: ThePrimAmp
type(UnitarityCut),pointer :: DoubCuts
type(Propagator) :: IntPart(NumExtParticles)
type(TreeProcess),pointer :: TreeProcs(:)
integer :: CutNum
double precision :: KMom(1,1:4),VMom(1,1:4)
double complex :: NMom(1:3,1:4),SqrtPreF
double precision :: V2(1:4),x1,x2,phi
double complex :: x3,x4,r1,r2,r3,r4
double complex :: LoopMom(1:5),alpha_1,alpha_2,alpha_3,alpha_5,alpha_4,alpha_e
double complex :: Res(1:10),s1(1:10),s2(1:10),s3(1:10),se2
double complex :: ResCheck,s1CH,s2CH,s3CH,seCH,ResCH,xtr
double complex :: LinSysEq(1:10,1:11)
integer :: k

double complex :: lhseq1,lhseq47a,lhseq36a,lhseq3,lhseq47b,lhseq36b,lhseq1new,xe,SqrtPreF2a,SqrtPreF2b
   include 'misc/global_import'


   DoubCuts => ThePrimAmp%UCuts(2)
   IntPart(1:NumExtParticles) = ThePrimAmp%IntPart(1:NumExtParticles)

   do CutNum = 1,DoubCuts%NumCuts
!do CutNum = 1,1; print *, "only doub cut 1"
!DEC$ IF (_DebugShowCuts==1)
   print *, "DoubCut ",CutNum,":",DoubCuts%CutProp(CutNum,1:2)
!DEC$ ENDIF
   if (DoubCuts%skip(CutNum)) then
         DoubCuts%Coeff(CutNum,:)=(0d0,0d0)
         cycle
      endif

      TreeProcs => ThePrimAmp%UCuts(2)%TreeProcess(CutNum,1:2)

      KMom(1,1:4) = IntPart(DoubCuts%CutProp(CutNum,2))%Mom(1:4) - IntPart(DoubCuts%CutProp(CutNum,1))%Mom(1:4)
      if ( dabs(KMom(1,1:4).dot.KMom(1,1:4)) .ge. 1d-7 ) then    !  NON-light like vector
         call NVBasis1(KMom(1,1:4),VMom(1,1:4),NMom(1:3,1:4))
         x1 = (KMom(1,1:4).dot.KMom(1,1:4)) + IntPart(DoubCuts%CutProp(CutNum,1))%Mass2 - IntPart(DoubCuts%CutProp(CutNum,2))%Mass2
         V2(1:4) = -0.5d0*x1*VMom(1,1:4)
         SqrtPreF = cdsqrt(dcmplx( IntPart(DoubCuts%CutProp(CutNum,1))%Mass2 - (V2.dot.V2) ))

         if (cdabs(SqrtPreF).gt. 1d-1) then     !  NON-vanishing SqrtF
            alpha_1 = drand(1)
            do k=1,9
               alpha_1 = drand(0)
               alpha_2 = drand(0)
               alpha_3 = drand(0)

               s1(k) = alpha_1*SqrtPreF/cdsqrt(alpha_1**2+alpha_2**2+alpha_3**2)
               s2(k) = alpha_2*SqrtPreF/cdsqrt(alpha_1**2+alpha_2**2+alpha_3**2)
               s3(k) = alpha_3*SqrtPreF/cdsqrt(alpha_1**2+alpha_2**2+alpha_3**2)

               LoopMom(1:4) = V2(1:4) + s1(k)*NMom(1,1:4) + s2(k)*NMom(2,1:4)  + s3(k)*NMom(3,1:4)
               LoopMom(5)   = 0d0

               LinSysEq(k,1) = 1d0
               LinSysEq(k,2) = s1(k)
               LinSysEq(k,3) = s2(k)
               LinSysEq(k,4) = s3(k)
!                LinSysEq(k,5) = s1(k)**2-s3(k)**2
!                LinSysEq(k,6) = s2(k)**2-s3(k)**2
               LinSysEq(k,5) = s1(k)**2-s2(k)**2
               LinSysEq(k,6) = s1(k)**2+s2(k)**2-2d0*s3(k)**2
               LinSysEq(k,7) = s1(k)*s2(k)
               LinSysEq(k,8) = s1(k)*s3(k)
               LinSysEq(k,9) = s2(k)*s3(k)
               call resid2(LoopMom(1:5),dcmplx(KMom(1,1:4)),DoubCuts%CutProp(CutNum,1:2),TreeProcs,Res(k))
            enddo
            LinSysEq(1:9,10) = Res(1:9)
            DoubCuts%Coeff(CutNum,0:8) = go_Gauss(9,LinSysEq(1:9,1:10))


         else     !  vanishing  SqrtF


          alpha_1 = drand(1)
          do k=1,9
                alpha_1 = drand(0)
                alpha_2 = drand(0)
                alpha_3 = cdsqrt(SqrtPreF**2-alpha_1**2-alpha_2**2)

                s1(k) = alpha_1
                s2(k) = alpha_2
                s3(k) = alpha_3

                LoopMom(1:4) = V2(1:4) + s1(k)*NMom(1,1:4) + s2(k)*NMom(2,1:4) + s3(k)*NMom(3,1:4)
                LoopMom(5)   = 0d0

                LinSysEq(k,1) = 1d0
                LinSysEq(k,2) = s1(k)
                LinSysEq(k,3) = s2(k)
                LinSysEq(k,4) = s3(k)
!                 LinSysEq(k,5) = s1(k)**2-s3(k)**2
!                 LinSysEq(k,6) = s2(k)**2-s3(k)**2
                LinSysEq(k,5) = s1(k)**2-s2(k)**2
                LinSysEq(k,6) = s1(k)**2+s2(k)**2-2d0*s3(k)**2
                LinSysEq(k,7) = s1(k)*s2(k)
                LinSysEq(k,8) = s1(k)*s3(k)
                LinSysEq(k,9) = s2(k)*s3(k)
                call resid2(LoopMom(1:5),dcmplx(KMom(1,1:4)),DoubCuts%CutProp(CutNum,1:2),TreeProcs,Res(k))
          enddo

          LinSysEq(1,10) = Res(1)
          LinSysEq(2,10) = Res(2)
          LinSysEq(3,10) = Res(3)
          LinSysEq(4,10) = Res(4)
          LinSysEq(5,10) = Res(5)
          LinSysEq(6,10) = Res(6)
          LinSysEq(7,10) = Res(7)
          LinSysEq(8,10) = Res(8)
          LinSysEq(9,10) = Res(9)
          DoubCuts%Coeff(CutNum,0:8) = go_Gauss(9,LinSysEq(1:9,1:10))

         endif

         alpha_5 = 1d0
         alpha_3 = cdsqrt( SqrtPreF**2-alpha_5**2 )
         s3(10)  = alpha_3
         se2     = alpha_5
         LoopMom(1:4) = V2(1:4) + s3(10)*NMom(3,1:4)
         LoopMom(5)   = (0d0,1d0)
         call resid2(LoopMom(1:5),dcmplx(KMom(1,1:4)),DoubCuts%CutProp(CutNum,1:2),TreeProcs,Res(10))

         DoubCuts%Coeff(CutNum,9) = Res(10) - DoubCuts%Coeff(CutNum,0)  &
                                            - DoubCuts%Coeff(CutNum,3)*s3(10)  &
!                                             + DoubCuts%Coeff(CutNum,4)*s3(10)**2
                                            + DoubCuts%Coeff(CutNum,5)*s3(10)**2*2d0



!!!!!!!!!  RR  ADDED 24 SEPT 2013 !!!!!!!!!!!!
! Checks that the OPP equations are solved correctly,
! by reconstructing an arbitrary loop momentum l and using it to compute residue from currents,
! and compare that to residue from OPP coefficients.
! The global variable opp_err gives the relative difference between these two values. 
! Checked for ttbZ production. 
         
         if (cdabs(SqrtPreF).gt. 1d-1) then     !  NON-vanishing SqrtF
            alpha_1 = 17d0/19d0
            alpha_2 = -5d0/7d0
            alpha_3 = 23d0/29d0
            alpha_4 = -31d0/39d0
            xtr = SqrtPreF/sqrt(alpha_1**2+alpha_2**2+alpha_3**2+alpha_4**2)
            LoopMom(1:4) = V2(1:4) + xtr*(alpha_1*NMom(1,1:4) + alpha_2*NMom(2,1:4) + alpha_3*NMom(3,1:4))
            LoopMom(5)   = xtr*alpha_4*(0d0,1d0)
            
            s1CH = (LoopMom(1:4).dot.NMom(1,1:4))
            s2CH = (LoopMom(1:4).dot.NMom(2,1:4))
            s3CH = (LoopMom(1:4).dot.NMom(3,1:4))
            seCH = (xtr*alpha_4)**2

            call resid2(LoopMom(1:5),dcmplx(KMom(1,1:4)),DoubCuts%CutProp(CutNum,1:2),TreeProcs,ResCH)
            ResCheck = (DoubCuts%Coeff(CutNum,0)+DoubCuts%Coeff(CutNum,1)*s1CH+&
                 + DoubCuts%Coeff(CutNum,2)*s2CH+DoubCuts%Coeff(CutNum,3)*s3CH  &
                 + DoubCuts%Coeff(CutNum,4)*(s1CH**2-s2CH**2)+&
                 + DoubCuts%Coeff(CutNum,5)*(s1CH**2+s2CH**2-2d0*s3CH**2) &
                 + DoubCuts%Coeff(CutNum,6)*s1CH*s2CH &
                 + DoubCuts%Coeff(CutNum,7)*s1CH*s3CH &
                 + DoubCuts%Coeff(CutNum,8)*s2CH*s3CH &
                 + seCH*DoubCuts%Coeff(CutNum,9) )

         else

            alpha_1 = 11d0/19d0
            alpha_2 = -1d0/7d0
            alpha_3 = 3d0/29d0
            xtr = cdsqrt(SqrtPreF**2 - alpha_1**2-alpha_2**2-alpha_3**2)
            LoopMom(1:4) = V2(1:4) + (alpha_1*NMom(1,1:4) + alpha_2*NMom(2,1:4) + alpha_3*NMom(3,1:4))
            LoopMom(5)   = xtr*(0d0,1d0)
            s1CH = (LoopMom(1:4).dot.NMom(1,1:4))
            s2CH = (LoopMom(1:4).dot.NMom(2,1:4))
            s3CH = (LoopMom(1:4).dot.NMom(3,1:4))
            seCH = (xtr)**2
            
            call resid2(LoopMom(1:5),dcmplx(KMom(1,1:4)),DoubCuts%CutProp(CutNum,1:2),TreeProcs,ResCH)           
            ResCheck = (DoubCuts%Coeff(CutNum,0)+DoubCuts%Coeff(CutNum,1)*s1CH+&
                 + DoubCuts%Coeff(CutNum,2)*s2CH+DoubCuts%Coeff(CutNum,3)*s3CH  &
                 + DoubCuts%Coeff(CutNum,4)*(s1CH**2-s2CH**2)+&
                 + DoubCuts%Coeff(CutNum,5)*(s1CH**2+s2CH**2-2d0*s3CH**2) &
                 + DoubCuts%Coeff(CutNum,6)*s1CH*s2CH &
                 + DoubCuts%Coeff(CutNum,7)*s1CH*s3CH + &
                 + DoubCuts%Coeff(CutNum,8)*s2CH*s3CH + &
                 + seCH*DoubCuts%Coeff(CutNum,9) )
         endif
         
         if ( cdabs(ResCheck).gt.1d-4) then
            opp_err =abs( (ResCH-ResCheck)/ResCheck)
         else
            opp_err = abs( (ResCH-ResCheck) )
         endif
!!!!!! CHECK ENDS HERE !!!!!!!!!!!!!!


!     set vars. for kirill's routines
         tagdcut(CutNum,1)=666    ! needed for integrals
         DoubCuts%tagcuts(CutNum)=666

         coeff2(CutNum,1) = DoubCuts%Coeff(CutNum,0)
         coeff2(CutNum,2) = DoubCuts%Coeff(CutNum,1)
         coeff2(CutNum,3) = DoubCuts%Coeff(CutNum,2)
         coeff2(CutNum,4) = DoubCuts%Coeff(CutNum,3)
         coeff2(CutNum,5) = DoubCuts%Coeff(CutNum,4)
         coeff2(CutNum,6) = DoubCuts%Coeff(CutNum,5)
         coeff2(CutNum,7) = DoubCuts%Coeff(CutNum,6)
         coeff2(CutNum,8) = DoubCuts%Coeff(CutNum,7)
         coeff2(CutNum,9) = DoubCuts%Coeff(CutNum,8)
         coeff2(CutNum,10)= DoubCuts%Coeff(CutNum,9)

         refvect2(CutNum,1:4)  = NMom(1,1:4)
         refvect2(CutNum,5:8)  = NMom(2,1:4)
         refvect2(CutNum,9:12) = NMom(3,1:4)

         mass2(CutNum,1) = IntPart( DoubCuts%CutProp(CutNum,1) )%Mass
         mass2(CutNum,2) = IntPart( DoubCuts%CutProp(CutNum,2) )%Mass
         propv2(CutNum,1:4)  = (0d0,0d0)
         propv2(CutNum,5:8)  = KMom(1,1:4)


      else    !  light-like-vector
         

         call NVBasis1_Light(KMom(1,1:4),VMom(1,1:4),NMom(1:3,1:4))

! -------- #1
         r1 = KMom(1,1:4).dot.NMom(1,1:4)
         x1 = 0.5d0
         x2 = ( IntPart(DoubCuts%CutProp(CutNum,2))%Mass2 - IntPart(DoubCuts%CutProp(CutNum,1))%Mass2 )/(2d0*r1)
         x3 = 1d0
         SqrtPreF = cdsqrt(dcmplx( IntPart(DoubCuts%CutProp(CutNum,1))%Mass2*(1d0+x1) - x1*IntPart(DoubCuts%CutProp(CutNum,2))%Mass2 ))
         x4 = cdsqrt(SqrtPreF**2 - x3**2)

         LoopMom(1:4) = x1*VMom(1,1:4) + x2*NMom(1,1:4) + x3*NMom(2,1:4) + x4*NMom(3,1:4)
         LoopMom(5)   = 0d0
         call resid2(LoopMom(1:5),dcmplx(KMom(1,1:4)),DoubCuts%CutProp(CutNum,1:2),TreeProcs,Res(1))

! -------- #2
         x3 = -x3
         x4 = x4
         LoopMom(1:4) = x1*VMom(1,1:4) + x2*NMom(1,1:4) + x3*NMom(2,1:4) + x4*NMom(3,1:4)
         LoopMom(5)   = 0d0
         call resid2(LoopMom(1:5),dcmplx(KMom(1,1:4)),DoubCuts%CutProp(CutNum,1:2),TreeProcs,Res(2))


! -------- #3
         x3 = -x3
         x4 = -x4
         LoopMom(1:4) = x1*VMom(1,1:4) + x2*NMom(1,1:4) + x3*NMom(2,1:4) + x4*NMom(3,1:4)
         LoopMom(5)   = 0d0
         call resid2(LoopMom(1:5),dcmplx(KMom(1,1:4)),DoubCuts%CutProp(CutNum,1:2),TreeProcs,Res(3))


! -------- #4
         x3 = -x3
         x4 = x4
         LoopMom(1:4) = x1*VMom(1,1:4) + x2*NMom(1,1:4) + x3*NMom(2,1:4) + x4*NMom(3,1:4)
         LoopMom(5)   = 0d0
         call resid2(LoopMom(1:5),dcmplx(KMom(1,1:4)),DoubCuts%CutProp(CutNum,1:2),TreeProcs,Res(4))

         lhseq1   = 0.25d0*(Res(1)+Res(2)+Res(3)+Res(4))
         DoubCuts%Coeff(CutNum,8) = 0.25d0/(-x4)*(Res(1)-Res(3)-Res(2)+Res(4))
         lhseq47a = 0.25d0/(-x4)*(Res(1)-Res(3)+Res(2)-Res(4))
         lhseq36a = (Res(1)-Res(2))/2d0 - DoubCuts%Coeff(CutNum,8)*(-x4)


! -------- #5
         x3 = 0.5d0
         x4 = zsqrt(SqrtPreF**2 - x3**2)
         LoopMom(1:4) = x1*VMom(1,1:4) + x2*NMom(1,1:4) + x3*NMom(2,1:4) + x4*NMom(3,1:4)
         LoopMom(5)   = 0d0
         call resid2(LoopMom(1:5),dcmplx(KMom(1,1:4)),DoubCuts%CutProp(CutNum,1:2),TreeProcs,Res(5))


! -------- #6
         x1 = -0.5d0
         x3 = 1d0
         SqrtPreF = cdsqrt(dcmplx( IntPart(DoubCuts%CutProp(CutNum,1))%Mass2*(1d0+x1) - x1*IntPart(DoubCuts%CutProp(CutNum,2))%Mass2 ))
         x4 = cdsqrt(SqrtPreF**2 - x3**2)
         LoopMom(1:4) = x1*VMom(1,1:4) + x2*NMom(1,1:4) + x3*NMom(2,1:4) + x4*NMom(3,1:4)
         LoopMom(5)   = 0d0
         call resid2(LoopMom(1:5),dcmplx(KMom(1,1:4)),DoubCuts%CutProp(CutNum,1:2),TreeProcs,Res(1))


! -------- #7
         x3 = -x3
         x4 = x4
         LoopMom(1:4) = x1*VMom(1,1:4) + x2*NMom(1,1:4) + x3*NMom(2,1:4) + x4*NMom(3,1:4)
         LoopMom(5)   = 0d0
         call resid2(LoopMom(1:5),dcmplx(KMom(1,1:4)),DoubCuts%CutProp(CutNum,1:2),TreeProcs,Res(2))


! -------- #8
         x3 = -x3
         x4 = -x4
         LoopMom(1:4) = x1*VMom(1,1:4) + x2*NMom(1,1:4) + x3*NMom(2,1:4) + x4*NMom(3,1:4)
         LoopMom(5)   = 0d0
         call resid2(LoopMom(1:5),dcmplx(KMom(1,1:4)),DoubCuts%CutProp(CutNum,1:2),TreeProcs,Res(3))


! -------- #9
         x3 = -x3
         x4 =  x4
         LoopMom(1:4) = x1*VMom(1,1:4) + x2*NMom(1,1:4) + x3*NMom(2,1:4) + x4*NMom(3,1:4)
         LoopMom(5)   = 0d0
         call resid2(LoopMom(1:5),dcmplx(KMom(1,1:4)),DoubCuts%CutProp(CutNum,1:2),TreeProcs,Res(4))


         lhseq3=0.25d0*(Res(1)+Res(2)+Res(3)+Res(4))
         lhseq47b=0.25d0/(-x4)*(Res(1)-Res(3)+Res(2)-Res(4))
         lhseq36b=(Res(1)-Res(2))/2d0-DoubCuts%Coeff(CutNum,8)*(-x4)
         DoubCuts%Coeff(CutNum,3)=0.5d0*(lhseq47a+lhseq47b)
         DoubCuts%Coeff(CutNum,6)=(lhseq47a-lhseq47b)/r1
         DoubCuts%Coeff(CutNum,2)=0.5d0*(lhseq36a+lhseq36b)
         DoubCuts%Coeff(CutNum,5)=(lhseq36a-lhseq36b)/r1


         x1 = 0.5d0
         SqrtPreF2a = IntPart(DoubCuts%CutProp(CutNum,1))%Mass2*(1d0+x1) - x1*IntPart(DoubCuts%CutProp(CutNum,2))%Mass2
         x3 = dcmplx(0.5d0,0d0)
         x4 = cdsqrt(SqrtPreF2a - x3**2)
         LoopMom(1:4) = x1*VMom(1,1:4) + x2*NMom(1,1:4) + x3*NMom(2,1:4) + x4*NMom(3,1:4)
         LoopMom(5)   = 0d0


         alpha_1 = LoopMom(1:4).dot.NMom(1,1:4)
         alpha_2 = LoopMom(1:4).dot.NMom(2,1:4)
         alpha_3 = LoopMom(1:4).dot.NMom(3,1:4)
         Res(5)=Res(5) - DoubCuts%Coeff(CutNum,2)*alpha_2-DoubCuts%Coeff(CutNum,3)*alpha_3-DoubCuts%Coeff(CutNum,5)*alpha_1*alpha_2 -DoubCuts%Coeff(CutNum,6)*alpha_1*alpha_3 -DoubCuts%Coeff(CutNum,8)*alpha_2*alpha_3
         x1 = -0.5d0
         SqrtPreF2b = IntPart(DoubCuts%CutProp(CutNum,1))%Mass2*(1d0+x1) - x1*IntPart(DoubCuts%CutProp(CutNum,2))%Mass2
         DoubCuts%Coeff(CutNum,7)=2d0/3d0*(lhseq1-Res(5))
         DoubCuts%Coeff(CutNum,1)= ( (lhseq1-DoubCuts%Coeff(CutNum,7)*(2d0-SqrtPreF2a))- (lhseq3-DoubCuts%Coeff(CutNum,7)*(2d0-SqrtPreF2b)) )/r1

! -------- #10
         x1 = 0.25d0
         SqrtPreF = cdsqrt(dcmplx( IntPart(DoubCuts%CutProp(CutNum,1))%Mass2*(1d0+x1) - x1*IntPart(DoubCuts%CutProp(CutNum,2))%Mass2 ))
         x3 = 1d0
         x4 = cdsqrt(SqrtPreF**2 - x3**2)
         LoopMom(1:4) = x1*VMom(1,1:4) + x2*NMom(1,1:4) + x3*NMom(2,1:4) + x4*NMom(3,1:4)
         LoopMom(5)   = 0d0
         call resid2(LoopMom(1:5),dcmplx(KMom(1,1:4)),DoubCuts%CutProp(CutNum,1:2),TreeProcs,Res(6))


         alpha_1 = LoopMom(1:4).dot.NMom(1,1:4)
         alpha_2 = LoopMom(1:4).dot.NMom(2,1:4)
         alpha_3 = LoopMom(1:4).dot.NMom(3,1:4)
         Res(6)=Res(6)-DoubCuts%Coeff(CutNum,1)*alpha_1-DoubCuts%Coeff(CutNum,2)*alpha_2-DoubCuts%Coeff(CutNum,3)*alpha_3-DoubCuts%Coeff(CutNum,5)*alpha_1*alpha_2 -DoubCuts%Coeff(CutNum,6)*alpha_1*alpha_3 -DoubCuts%Coeff(CutNum,7)*(alpha_2**2-alpha_3**2)-DoubCuts%Coeff(CutNum,8)*alpha_2*alpha_3
         lhseq1new=lhseq1-DoubCuts%Coeff(CutNum,7)*(dcmplx(2d0,0d0)-SqrtPreF2a) -DoubCuts%Coeff(CutNum,1)/2d0*r1
         DoubCuts%Coeff(CutNum,0)=(4d0*Res(6) - lhseq1new)/3d0
         DoubCuts%Coeff(CutNum,4)=16d0/r1*(Res(6) - DoubCuts%Coeff(CutNum,0))


! -------- #11
         x1 = -0.5d0
         SqrtPreF = cdsqrt(dcmplx( IntPart(DoubCuts%CutProp(CutNum,1))%Mass2*(1d0+x1) - x1*IntPart(DoubCuts%CutProp(CutNum,2))%Mass2 ))
         x3 = (0d0,0d0)
         xe = 0.5768943d0
         x4 = cdsqrt(SqrtPreF**2 - xe**2 - x3**2)
         LoopMom(1:4) = x1*VMom(1,1:4) + x2*NMom(1,1:4) + x3*NMom(2,1:4) + x4*NMom(3,1:4)
         LoopMom(5)   = xe * (0d0,1d0)
         call resid2(LoopMom(1:5),dcmplx(KMom(1,1:4)),DoubCuts%CutProp(CutNum,1:2),TreeProcs,Res(1))


         DoubCuts%Coeff(CutNum,9) = 1d0/xe**2*(Res(1) - DoubCuts%Coeff(CutNum,0) - DoubCuts%Coeff(CutNum,1)*r1*x1 - DoubCuts%Coeff(CutNum,2)*x3-DoubCuts%Coeff(CutNum,3)*x4 -DoubCuts%Coeff(CutNum,4)*r1**2*x1**2 -DoubCuts%Coeff(CutNum,5)*r1*x1*x3 - DoubCuts%Coeff(CutNum,6)*r1*x1*x4-DoubCuts%Coeff(CutNum,7)*(x3**2 - x4**2) - DoubCuts%Coeff(CutNum,8)*x3*x4 )



!!!!!!!!!  RR  ADDED 24 SEPT 2013 !!!!!!!!!!!!
! Checks that the OPP equations are solved correctly,
! by reconstructing an arbitrary loop momentum l and using it to compute residue from currents,
! and compare that to residue from OPP coefficients.
! The global variable opp_err gives the relative difference between these two values. 
! N.B. NOT TESTED YET -- CAVEAT EMPTOR! 

         alpha_1 = 17d0/19d0
         alpha_2 = ( IntPart(DoubCuts%CutProp(CutNum,2))%Mass2 - IntPart(DoubCuts%CutProp(CutNum,1))%Mass2 )/(2d0*r1)
         alpha_3 = -1d0/7d0
         alpha_e = 3d0/29d0
         alpha_4 = cdsqrt(SqrtPreF**2 - alpha_e**2 - alpha_3**2)
         LoopMom(1:4) = alpha_1*VMom(1,1:4) + alpha_2*NMom(1,1:4) + alpha_3*NMom(2,1:4) + alpha_4*NMom(3,1:4)
         LoopMom(5)   = alpha_e * (0d0,1d0)
         s1CH = LoopMom(1:4) .dot. NMom(1,1:4)
         s2CH = LoopMom(1:4) .dot. NMom(2,1:4)
         s3CH = LoopMom(1:4) .dot. NMom(3,1:4)
         seCH = alpha_e**2

         call resid2(LoopMom(1:5),dcmplx(KMom(1,1:4)),DoubCuts%CutProp(CutNum,1:2),TreeProcs,ResCH)
         ResCheck = DoubCuts%Coeff(CutNum,0) + DoubCuts%Coeff(CutNum,1)*s1CH + DoubCuts%Coeff(CutNum,2)*s2CH + &
              & DoubCuts%Coeff(CutNum,3)*s3CH +  DoubCuts%Coeff(CutNum,4)*s1CH**2 + DoubCuts%Coeff(CutNum,5)*s1CH*s2CH &
              & + DoubCuts%Coeff(CutNum,6)*s1CH*s3CH + DoubCuts%Coeff(CutNum,7)*(s2CH**2 - s3CH**2) &
              & + DoubCuts%Coeff(CutNum,8)*s2CH*s3CH + DoubCuts%Coeff(CutNum,9)*seCH
!         print *, 'resCH',resCH
!         print *, 'ResCheck',ResCheck

          if ( cdabs(ResCheck).gt.1d-4) then
             opp_err = (ResCH-ResCheck)/ResCheck
          else
             opp_err = (ResCH-ResCheck)
          endif
!          print *, 'err', opp_err
!          pause

          ! end stability check in D dim

         

         tagdcut(CutNum,1)=999    ! needed for integrals
         DoubCuts%tagcuts(CutNum)=999


         coeff2(CutNum,1) = DoubCuts%Coeff(CutNum,0)
         coeff2(CutNum,2) = DoubCuts%Coeff(CutNum,1)
         coeff2(CutNum,3) = DoubCuts%Coeff(CutNum,2)
         coeff2(CutNum,4) = DoubCuts%Coeff(CutNum,3)
         coeff2(CutNum,5) = DoubCuts%Coeff(CutNum,4)
         coeff2(CutNum,6) = DoubCuts%Coeff(CutNum,5)
         coeff2(CutNum,7) = DoubCuts%Coeff(CutNum,6)
         coeff2(CutNum,8) = DoubCuts%Coeff(CutNum,7)
         coeff2(CutNum,9) = DoubCuts%Coeff(CutNum,8)
         coeff2(CutNum,10)= DoubCuts%Coeff(CutNum,9)

         refvect2(CutNum,1:4)  = NMom(1,1:4)
         refvect2(CutNum,5:8)  = NMom(2,1:4)
         refvect2(CutNum,9:12) = NMom(3,1:4)

         mass2(CutNum,1) = IntPart( DoubCuts%CutProp(CutNum,1) )%Mass
         mass2(CutNum,2) = IntPart( DoubCuts%CutProp(CutNum,2) )%Mass
         propv2(CutNum,1:4)  = (0d0,0d0)
         propv2(CutNum,5:8)  = KMom(1,1:4)

      endif

!       do k=1,10
!       if( cdabs(coeff2(CutNum,k)).lt.1d-6 ) coeff2(CutNum,k) = (0d0,0d0)
!       enddo
!       print *, "# doubcut",CutNum
!       do k=1,10
!       write(*,"(3X,1PE23.16,3X,1PE23.16)") dreal(coeff2(CutNum,k)),dimag(coeff2(CutNum,k))
!       enddo

!       if ( dabs(KMom(1,1:4).dot.KMom(1,1:4)) .ge. 1d-7 ) then    !  NON-light like vector
!       do k=1,10
!         write(*,"(I,3X,1PE23.16,3X,1PE23.16)") 1,dreal(coeff2(CutNum,k)),dimag(coeff2(CutNum,k))
!       enddo
!       else
!       if (cdabs(SqrtPreF).gt. 1d-1) then     !  NON-vanishing SqrtF
!       do k=1,10
!         write(*,"(I,3X,1PE23.16,3X,1PE23.16)") 2,dreal(coeff2(CutNum,k)),dimag(coeff2(CutNum,k))
!       enddo
!       else
!       do k=1,10
!         write(*,"(I,3X,1PE23.16,3X,1PE23.16)") 3,dreal(coeff2(CutNum,k)),dimag(coeff2(CutNum,k))
!       enddo
!       endif
!       endif

!DEC$ IF (_DebugPrintDoubCoeff==1)
      do k=1,10
         print *, cutnum,k,coeff2(Cutnum,k)
      enddo
      enddo
      pause
!DEC$ ELSE
      enddo
!DEC$ ENDIF


END SUBROUTINE





SUBROUTINE SingCut(ThePrimAmp)
use ModProcess
use ModParameters
use ModNVBasis
use ModMisc
use ModResidues
implicit none
type(PrimitiveAmplitude),target :: ThePrimAmp
type(UnitarityCut),pointer :: SingCuts
type(Propagator) :: IntPart(NumExtParticles)
type(TreeProcess),pointer :: TreeProcs(:)
integer :: CutNum
double complex :: NMom(1:4,1:4),SqrtPreF
double complex :: LoopMom(1:5)
double complex :: Res(1:2)

   include 'misc/global_import'


   SingCuts => ThePrimAmp%UCuts(1)
   IntPart(1:NumExtParticles) = ThePrimAmp%IntPart(1:NumExtParticles)

!  loop over all sing. cuts
   do CutNum = 1,SingCuts%NumCuts
!DEC$ IF (_DebugShowCuts==1)
      print *, "SingCut ",CutNum,":",SingCuts%CutProp(CutNum,1:1)
!DEC$ ENDIF
      if (SingCuts%skip(CutNum)) then
         SingCuts%Coeff(CutNum,:)=(0d0,0d0)
         cycle
      endif

       TreeProcs => ThePrimAmp%UCuts(1)%TreeProcess(CutNum,1:1)


       NMom(1,1) = (1d0,0d0)
       NMom(1,2) = (0d0,0d0)
       NMom(1,3) = (0d0,0d0)
       NMom(1,4) = (0d0,0d0)

       NMom(2,1) = (0d0,0d0)
       NMom(2,2) = (0d0,1d0)
       NMom(2,3) = (0d0,0d0)
       NMom(2,4) = (0d0,0d0)

       NMom(3,1) = (0d0,0d0)
       NMom(3,2) = (0d0,0d0)
       NMom(3,3) = (0d0,1d0)
       NMom(3,4) = (0d0,0d0)

       NMom(4,1) = (0d0,0d0)
       NMom(4,2) = (0d0,0d0)
       NMom(4,3) = (0d0,0d0)
       NMom(4,4) = (0d0,1d0)

       LoopMom(1:4) = IntPart(SingCuts%CutProp(CutNum,1))%Mass*(NMom(1,1:4)+NMom(2,1:4)+NMom(3,1:4)+NMom(4,1:4))/2d0
       LoopMom(5)   = 0d0
       call resid1(LoopMom(1:5),SingCuts%CutProp(CutNum,1:1),TreeProcs,Res(1))

       LoopMom(1:4) = - LoopMom(1:4)
       LoopMom(5)   = 0d0
       call resid1(LoopMom(1:5),SingCuts%CutProp(CutNum,1:1),TreeProcs,Res(2))
       SingCuts%Coeff(CutNum,0) = 0.5d0 * ( Res(1) + Res(2) )


!     set vars. for kirill's routines
      coeff1(CutNum,1) = SingCuts%Coeff(CutNum,0)
      mass1(CutNum,1)  = IntPart( SingCuts%CutProp(CutNum,1) )%Mass

!       print *, "# singcut",CutNum
!       write(*,"(1PE23.16,3X,1PE23.16)") dreal(coeff1(CutNum,1)),dimag(coeff1(CutNum,1))
!       if( cdabs(coeff1(CutNum,1)).lt.1d-6 ) coeff1(CutNum,1) = (0d0,0d0)


!DEC$ IF (_DebugPrintSingCoeff==1)
         print *, cutnum,coeff1(Cutnum,1)
      enddo
      pause
!DEC$ ELSE
      enddo
!DEC$ ENDIF


END SUBROUTINE

















END MODULE
