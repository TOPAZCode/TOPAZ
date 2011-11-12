MODULE ModIntegrand
implicit none




contains








FUNCTION Integrand_01(yRnd,VgsWgt)
use ModParameters
use ModKinematics
use ModMisc
implicit none
real(8) :: Integrand_01
real(8) :: yRnd(1:VegasMxDim),VgsWgt
include 'vegas_common.f'
real(8) :: Mom(1:4,1:5),PSWgt,TheDipole


  Integrand_01 = 0d0
  Mom(1:4,1) = (/m_Top,0d0,0d0,0d0/)
  call EvalPhasespace(4,m_Top,(/m_W,0d0,0d0,0d0/),yRnd(1:8),Mom(1:4,2:5),PSWgt)
!   TheDipole = Dipole_fi_qg(Mom(1:4,1:5))
  TheDipole = Dipole_fi_gg(Mom(1:4,1:5))

  Integrand_01 = PSWgt * TheDipole * 4d0*Pi

RETURN
END FUNCTION






FUNCTION Integrand_02(yRnd,VgsWgt)
use ModParameters
use ModKinematics
use ModMisc
implicit none
real(8) :: Integrand_02
real(8) :: yRnd(1:VegasMxDim),VgsWgt
include 'vegas_common.f'
real(8) :: Mom(1:4,1:4),PSWgt,TheIntDipole


  Integrand_02 = 0d0
  Mom(1:4,1) = (/m_Top,0d0,0d0,0d0/)
  call EvalPhasespace(3,m_Top,(/m_W,0d0,0d0/),yRnd(1:5),Mom(1:4,2:4),PSWgt)
!   TheIntDipole = IntDipole_fi_qg(Mom(1:4,1:4),1,3)
  TheIntDipole = IntDipole_fi_gg(Mom(1:4,1:4),1,4)


  Integrand_02 = PSWgt * TheIntDipole / 2d0/Pi


RETURN
END FUNCTION











FUNCTION Dipole_fi_qg(Mom)
use ModParameters
use ModKinematics
use ModMisc
use ModProcess
use ModAmplitudes
implicit none
real(8) :: Dipole_fi_qg
real(8) :: Mom(1:4,1:5),MomT(1:4,1:4),MomB(1:4),Mass,pbDpg,ptDpg,ptDpb,omz,rsq,z,y,ymax,x,paux(1:4),MomA(1:4),Pia(1:4),Pia2,Ria,r2,f
real(8) :: mt,SingDepth,MadGraph_tree,MadGraph_tree2,mw,DipColorCorr=0d0,MomK(1:4),MomKTilde(1:4)
real(8) :: dipole,dipole1,dipole2,dipole3,dipole4,dipole5,dipole6,dipole7
integer :: Pcol1,Pcol2,Steps,n,i,j,k,hel1,hel2,hel3,hel4,hel5
complex(8) :: Ampl1,Ampl2,AmplH(-1:1),SqMat,SqMatDip,diag,offdiag,xm,xp,HH(-1:1,-1:1)
integer :: NumQuarks,NumGluons,NumTrees,SingType=-99
type(TreeProcess),save :: TreeAmpsDip(1:2)
integer :: NumParticles,iTree
type(Particle),save :: ExtParticles(1:5)
integer,parameter :: ff_=1,fi_=2,if_=3
logical,save :: init=.true.

  mt=m_top; mw=m_w

  if(init) then
          NumQuarks=2
          NumGluons=1
          NumParticles=NumQuarks+NumGluons+1
          allocate( TreeAmpsDip(1)%NumGlu(0:NumQuarks)     )
          allocate( TreeAmpsDip(1)%PartRef(1:NumParticles) )
          allocate( TreeAmpsDip(1)%PartType(1:NumParticles))
          allocate( TreeAmpsDip(1)%Quarks(1:NumQuarks)     )
          allocate( TreeAmpsDip(1)%Gluons(1:NumGluons)     )
          TreeAmpsDip(1)%NumPart=NumParticles
          TreeAmpsDip(1)%NumQua=NumQuarks
          TreeAmpsDip(1)%NumGlu(0) = NumGluons
          ExtParticles(1)%PartType = 5! top
          ExtParticles(1)%ExtRef   = 1
          ExtParticles(1)%Mass = mt
          ExtParticles(1)%Mass2= ExtParticles(1)%Mass**2
          ExtParticles(2)%PartType = 4! str
          ExtParticles(2)%ExtRef   = 2
          ExtParticles(2)%Mass = 0d0
          ExtParticles(2)%Mass2= ExtParticles(2)%Mass**2
          ExtParticles(3)%PartType = 10! glu
          ExtParticles(3)%ExtRef   = 3
          ExtParticles(3)%Mass = 0d0
          ExtParticles(3)%Mass2= ExtParticles(2)%Mass**2
          ExtParticles(4)%PartType = 13! w+
          ExtParticles(4)%ExtRef   = 4
          ExtParticles(4)%Mass = mw
          ExtParticles(4)%Mass2= ExtParticles(4)%Mass**2
          TreeAmpsDip(1)%PartRef(1:4) = (/1,2,3,4/)
          call LinkTreeParticles(TreeAmpsDip(1),ExtParticles(1:4))
          init=.false.
  endif

  Dipole_fi_qg = 0d0


          i=3; j=5; k=1;
          z = 1d0-(    (Mom(1:4,k).dot.Mom(1:4,j))/((Mom(1:4,i).dot.Mom(1:4,k))+(Mom(1:4,k).dot.Mom(1:4,j))-(Mom(1:4,i).dot.Mom(1:4,j)))     )
          Ria = 2d0*(Mom(1:4,j).dot.Mom(1:4,k))/(1d0-z)! = mt^2 * (1-r^2)
          y = (Mom(1:4,i).dot.Mom(1:4,j))*2d0/m_Top**2 / (1d0-dsqrt(1d0-Ria/m_Top**2))**2
          ymax = (1d0+dsqrt(1d0-Ria/m_Top**2))**2/(z+(1d0-z)*(1d0-Ria/m_Top**2))*z*(1d0-z)

!        new mapping
         MomK(1:4) = Mom(1:4,2)+Mom(1:4,4)
         call DipoleTrafox(Mom(1:4,1),MomK(1:4),Mom(1:4,2),MomT(1:4,2))
         call DipoleTrafox(Mom(1:4,1),MomK(1:4),Mom(1:4,4),MomT(1:4,4))
         MomT(1:4,3) = Mom(1:4,1) - MomT(1:4,2) - MomT(1:4,4)

          if( MomT(1,4)/mt.lt.1d-2 .or. (MomT(1:4,3).dot.MomT(1:4,4))/mt**2.lt.1d-2) return


          ExtParticles(1)%Mom(1:4) = dcmplx(MomT(1:4,1))
          ExtParticles(2)%Mom(1:4) = dcmplx(MomT(1:4,3))
          ExtParticles(3)%Mom(1:4) = dcmplx(MomT(1:4,4))
          ExtParticles(4)%Mom(1:4) = dcmplx(MomT(1:4,2))
                  HH(+1,+1) = 1d0
                  HH(-1,-1) = 1d0
                  HH(-1,+1) = 0d0
                  HH(+1,-1) = 0d0
          SqMatDip = (0d0,0d0)
          do hel1=1,-1,-2
          do hel2=1,-1,-2
          do hel4=1,-1,-1
          do hel3=1,-1,-2
                  call uSpi(ExtParticles(1)%Mom(1:4),ExtParticles(1)%Mass,Hel1,ExtParticles(1)%Pol(1:4))
                  call ubarSpi(ExtParticles(2)%Mom(1:4),ExtParticles(2)%Mass,Hel2,ExtParticles(2)%Pol(1:4))
                  call pol_mless(ExtParticles(3)%Mom(1:4),Hel3,ExtParticles(3)%Pol(1:4))
                  call pol_massSR(ExtParticles(4)%Mom(1:4),mw,Hel4,ExtParticles(4)%Pol(1:4))

                  call EvalTree2(TreeAmpsDip(1),AmplH(hel3))
        enddo
                  SqMatDip = SqMatDip + HH(+1,+1)*AmplH(+1)*dconjg(AmplH(+1)) &
                                      + HH(-1,-1)*AmplH(-1)*dconjg(AmplH(-1)) &
                                      + HH(-1,+1)*AmplH(-1)*dconjg(AmplH(+1)) &
                                      + HH(+1,-1)*AmplH(+1)*dconjg(AmplH(-1))
        enddo
        enddo
        enddo
!                                        4pi*alpha/2/sw2              (4pi*alpha_s)
        SqMatDip = SqMatDip * 1d0/2d0/3d0* 0.2063277518095031333d0  *  (1.633628179866692484d0)


!             DipColorCorr = +2d0/3d0
!                           8*pi*alpha_s
!             Dipole_fi_qg = 3.26725635973d0 * SqMatDip * ( (1d0+z**2)/(1d0-z) - (Mom(1:4,i).dot.Mom(1:4,j))*mt**2/(Mom(1:4,k).dot.Mom(1:4,j))**2 ) /(Mom(1:4,i).dot.Mom(1:4,j))/2d0  * DipColorCorr
            Dipole_fi_qg = SqMatDip *( (1d0+z**2)/(1d0-z) - (Mom(1:4,i).dot.Mom(1:4,j))*mt**2/(Mom(1:4,k).dot.Mom(1:4,j))**2 ) /(Mom(1:4,i).dot.Mom(1:4,j))/2d0
            Dipole_fi_qg = Dipole_fi_qg * (-1d0) * StepFunc( 1d0-alpha_DKfi-z )*StepFunc( y-alpha_DKfi*ymax )
RETURN
END FUNCTION






FUNCTION Dipole_fi_gg(Mom)
use ModParameters
use ModKinematics
use ModMisc
use ModProcess
use ModAmplitudes
implicit none
real(8) :: Dipole_fi_gg
real(8) :: Mom(1:4,1:5),MomT(1:4,1:4),MomB(1:4),Mass,pbDpg,ptDpg,ptDpb,omz,rsq,z,y,ymax,x,paux(1:4),MomA(1:4),Pia(1:4),Pia2,Ria,r2,f
real(8) :: mt,SingDepth,MadGraph_tree,MadGraph_tree2,mw,DipColorCorr=0d0,MomK(1:4),MomKTilde(1:4)
real(8) :: dipole,dipole1,dipole2,dipole3,dipole4,dipole5,dipole6,dipole7
integer :: Pcol1,Pcol2,Steps,n,i,j,k,hel1,hel2,hel3,hel4,hel5
complex(8) :: Ampl1,Ampl2,AmplH(-1:1),SqMat,SqMatDip,diag,offdiag,xm,xp,HH(-1:1,-1:1)
integer :: NumQuarks,NumGluons,NumTrees,SingType=-99
type(TreeProcess),save :: TreeAmpsDip(1:2)
integer :: NumParticles,iTree
type(Particle),save :: ExtParticles(1:5)
integer,parameter :: ff_=1,fi_=2,if_=3
logical,save :: init=.true.


  mt=m_top; mw=m_w

  if(init) then
          NumQuarks=2
          NumGluons=1
          NumParticles=NumQuarks+NumGluons+1
          allocate( TreeAmpsDip(1)%NumGlu(0:NumQuarks)     )
          allocate( TreeAmpsDip(1)%PartRef(1:NumParticles) )
          allocate( TreeAmpsDip(1)%PartType(1:NumParticles))
          allocate( TreeAmpsDip(1)%Quarks(1:NumQuarks)     )
          allocate( TreeAmpsDip(1)%Gluons(1:NumGluons)     )
          TreeAmpsDip(1)%NumPart=NumParticles
          TreeAmpsDip(1)%NumQua=NumQuarks
          TreeAmpsDip(1)%NumGlu(0) = NumGluons
          ExtParticles(1)%PartType = 5! top
          ExtParticles(1)%ExtRef   = 1
          ExtParticles(1)%Mass = mt
          ExtParticles(1)%Mass2= ExtParticles(1)%Mass**2
          ExtParticles(2)%PartType = 4! str
          ExtParticles(2)%ExtRef   = 2
          ExtParticles(2)%Mass = 0d0
          ExtParticles(2)%Mass2= ExtParticles(2)%Mass**2
          ExtParticles(3)%PartType = 10! glu
          ExtParticles(3)%ExtRef   = 3
          ExtParticles(3)%Mass = 0d0
          ExtParticles(3)%Mass2= ExtParticles(2)%Mass**2
          ExtParticles(4)%PartType = 13! w+
          ExtParticles(4)%ExtRef   = 4
          ExtParticles(4)%Mass = mw
          ExtParticles(4)%Mass2= ExtParticles(4)%Mass**2
          TreeAmpsDip(1)%PartRef(1:4) = (/1,2,3,4/)
          call LinkTreeParticles(TreeAmpsDip(1),ExtParticles(1:4))
          init=.false.
  endif

  Dipole_fi_gg = 0d0

          i=4; j=5; k=1;
          z = 1d0-(    (Mom(1:4,k).dot.Mom(1:4,j))/((Mom(1:4,i).dot.Mom(1:4,k))+(Mom(1:4,k).dot.Mom(1:4,j))-(Mom(1:4,i).dot.Mom(1:4,j)))     )! C-E z

          Ria = 2d0*(Mom(1:4,j).dot.Mom(1:4,k))/(1d0-z)! = mt^2 * (1-r^2)
          y = (Mom(1:4,i).dot.Mom(1:4,j))*2d0/m_Top**2 / (1d0-dsqrt(1d0-Ria/m_Top**2))**2
          ymax = (1d0+dsqrt(1d0-Ria/m_Top**2))**2/(z+(1d0-z)*(1d0-Ria/m_Top**2))*z*(1d0-z)

          if( StepFunc( 1d0-alpha_DKfi-z )*StepFunc( y-alpha_DKfi*ymax ).eq.0d0 ) return

          Pia(1:4) = Mom(1:4,i)+Mom(1:4,j)-Mom(1:4,k)
          Pia2 = Pia(1:4).dot.Pia(1:4)
          MomT(1:4,1:4)=Mom(1:4,1:4)
          MomT(1:4,i) = dsqrt(lambda(Pia2,mt**2,0d0))/dsqrt(lambda(2d0*(Mom(1:4,i).dot.Mom(1:4,j)),Pia2,mt**2)) * (Mom(1:4,k)-(Pia(1:4).dot.Mom(1:4,k))/Pia2*Pia(1:4))  &
                      + (Pia2-mt**2)/2d0/Pia2*Pia(1:4)
          MomT(1:4,k) = MomT(1:4,i) - Pia(1:4)


!        new mapping
         MomK(1:4) = Mom(1:4,2)+Mom(1:4,3)
         call DipoleTrafox(Mom(1:4,1),MomK(1:4),Mom(1:4,2),MomT(1:4,2))
         call DipoleTrafox(Mom(1:4,1),MomK(1:4),Mom(1:4,3),MomT(1:4,3))
         MomT(1:4,4) = Mom(1:4,1) - MomT(1:4,2) - MomT(1:4,3)



          if( MomT(1,4)/mt.lt.1d-2 .or. (MomT(1:4,3).dot.MomT(1:4,4))/mt**2.lt.1d-2) return

          ExtParticles(1)%Mom(1:4) = dcmplx(MomT(1:4,1))
          ExtParticles(2)%Mom(1:4) = dcmplx(MomT(1:4,3))
          ExtParticles(3)%Mom(1:4) = dcmplx(MomT(1:4,4))
          ExtParticles(4)%Mom(1:4) = dcmplx(MomT(1:4,2))

                  diag= 2d0*( z/(1d0-z)   - mt**2/2d0*(Mom(1:4,i).dot.Mom(1:4,j))/(Mom(1:4,k).dot.Mom(1:4,j))**2  )  ! works with Dittmaier z and C-E z
!                   diag= 2d0*( z/(1d0-z)   + mt**2/2d0*(Mom(1:4,i).dot.Mom(1:4,j))/(Mom(1:4,k).dot.Mom(1:4,j))**2  )  ! works with Dittmaier z and C-E z
                  offdiag = 1d0/(Mom(1:4,i).dot.Mom(1:4,j))

                  MomA(1:4) = ( (Mom(1:4,k).dot.Mom(1:4,i))*Mom(1:4,i) - (Mom(1:4,k).dot.Mom(1:4,j))*Mom(1:4,j) )*2d0/Ria
                  Pia(1:4) = MomT(1:4,4)

                  f = 4d0/Ria**2 * ((Mom(1:4,k).dot.MomK(1:4))**2- mt**2*(mt**2-Ria) )
                  paux(1:4) = 1d0/(Mom(1:4,1).dot.Pia(1:4))*( (Mom(1:4,1).dot.Pia(1:4))*MomA(1:4)  - (MomA(1:4).dot.Pia(1:4))*Mom(1:4,1)  ) * f

                  call pol_mless(ExtParticles(3)%Mom(1:4),-1,ExtParticles(3)%Pol(1:4))
                  xm = dcmplx(paux(1:4)).dot.ExtParticles(3)%Pol(1:4)
                  call pol_mless(ExtParticles(3)%Mom(1:4),+1,ExtParticles(3)%Pol(1:4))
                  xp = dcmplx(paux(1:4)).dot.ExtParticles(3)%Pol(1:4)
                  HH(+1,+1) = diag + offdiag*conjg(xp)*xp
                  HH(-1,-1) = diag + offdiag*conjg(xm)*xm
                  HH(-1,+1) = offdiag*conjg(xm)*xp
                  HH(+1,-1) = offdiag*conjg(xp)*xm

          SqMatDip = (0d0,0d0)
          do hel1=1,-1,-2
          do hel2=1,-1,-2
          do hel4=1,-1,-1
          do hel3=1,-1,-2
                  call uSpi(ExtParticles(1)%Mom(1:4),ExtParticles(1)%Mass,Hel1,ExtParticles(1)%Pol(1:4))
                  call ubarSpi(ExtParticles(2)%Mom(1:4),ExtParticles(2)%Mass,Hel2,ExtParticles(2)%Pol(1:4))
                  call pol_mless(ExtParticles(3)%Mom(1:4),Hel3,ExtParticles(3)%Pol(1:4))
                  call pol_massSR(ExtParticles(4)%Mom(1:4),mw,Hel4,ExtParticles(4)%Pol(1:4))

                  call EvalTree2(TreeAmpsDip(1),AmplH(hel3))
        enddo
                  SqMatDip = SqMatDip + HH(+1,+1)*AmplH(+1)*dconjg(AmplH(+1)) &
                                      + HH(-1,-1)*AmplH(-1)*dconjg(AmplH(-1)) &
                                      + HH(-1,+1)*AmplH(-1)*dconjg(AmplH(+1)) &
                                      + HH(+1,-1)*AmplH(+1)*dconjg(AmplH(-1))
        enddo
        enddo
        enddo
!                                        4pi*alpha/2/sw2              (4pi*alpha_s)
        SqMatDip = SqMatDip * 1d0/2d0/3d0* 0.2063277518095031333d0  *  (1.633628179866692484d0)
!         SqMatDip = diag


!             DipColorCorr = -6d0
!             Dipole_fi_gg = 3.26725635973d0 * SqMatDip * 1d0/(Mom(1:4,i).dot.Mom(1:4,j))/2d0 *DipColorCorr !  CS f-i dipole for gg


            Dipole_fi_gg =SqMatDip * 1d0/(Mom(1:4,i).dot.Mom(1:4,j))/2d0
            Dipole_fi_gg = Dipole_fi_gg * (-1d0) * StepFunc( 1d0-alpha_DKfi-z )*StepFunc( y-alpha_DKfi*ymax )

RETURN
END FUNCTION










FUNCTION IntDipole_fi_qg(Mom,s,ij)
use ModParameters
use ModKinematics
use ModMisc
use ModProcess
use ModAmplitudes
implicit none
real(8) :: IntDipole_fi_qg
real(8) :: Mom(1:4,1:4),MomT(1:4,1:4),MomB(1:4),Mass,pbDpg,ptDpg,ptDpb,omz,z,y,ymax,x,paux(1:4),MomA(1:4),Pia(1:4),Pia2,Ria,r2,f
real(8) :: mt,SingDepth,MadGraph_tree,MadGraph_tree2,mw,DipColorCorr=0d0,MomK(1:4),MomKTilde(1:4),rsq,prec(1:4)
real(8) :: dipole,dipole1,dipole2,dipole3,dipole4,dipole5,dipole6,dipole7
integer :: Pcol1,Pcol2,Steps,n,i,j,k,hel1,hel2,hel3,hel4,hel5,s,ij
complex(8) :: Ampl1,Ampl2,AmplH(-1:1),SqMat,SqMatDip,diag,offdiag,xm,xp,HH(-1:1,-1:1)
integer :: NumQuarks,NumGluons,NumTrees,SingType=-99
type(TreeProcess),save :: TreeAmpsDip(1:2)
integer :: NumParticles,iTree
type(Particle),save :: ExtParticles(1:5)
integer,parameter :: ff_=1,fi_=2,if_=3
logical,save :: init=.true.

  mt=m_top; mw=m_w

  if(init) then
          NumQuarks=2
          NumGluons=1
          NumParticles=NumQuarks+NumGluons+1
          allocate( TreeAmpsDip(1)%NumGlu(0:NumQuarks)     )
          allocate( TreeAmpsDip(1)%PartRef(1:NumParticles) )
          allocate( TreeAmpsDip(1)%PartType(1:NumParticles))
          allocate( TreeAmpsDip(1)%Quarks(1:NumQuarks)     )
          allocate( TreeAmpsDip(1)%Gluons(1:NumGluons)     )
          TreeAmpsDip(1)%NumPart=NumParticles
          TreeAmpsDip(1)%NumQua=NumQuarks
          TreeAmpsDip(1)%NumGlu(0) = NumGluons
          ExtParticles(1)%PartType = 5! top
          ExtParticles(1)%ExtRef   = 1
          ExtParticles(1)%Mass = mt
          ExtParticles(1)%Mass2= ExtParticles(1)%Mass**2
          ExtParticles(2)%PartType = 4! str
          ExtParticles(2)%ExtRef   = 2
          ExtParticles(2)%Mass = 0d0
          ExtParticles(2)%Mass2= ExtParticles(2)%Mass**2
          ExtParticles(3)%PartType = 10! glu
          ExtParticles(3)%ExtRef   = 3
          ExtParticles(3)%Mass = 0d0
          ExtParticles(3)%Mass2= ExtParticles(2)%Mass**2
          ExtParticles(4)%PartType = 13! w+
          ExtParticles(4)%ExtRef   = 4
          ExtParticles(4)%Mass = mw
          ExtParticles(4)%Mass2= ExtParticles(4)%Mass**2
          TreeAmpsDip(1)%PartRef(1:4) = (/1,2,3,4/)
          call LinkTreeParticles(TreeAmpsDip(1),ExtParticles(1:4))
          init=.false.
  endif

  IntDipole_fi_qg = 0d0

   prec(:) = +Mom(:,s) -Mom(:,ij)  ! assumes all are outgoing
   rsq = scr(prec,prec)/mt**2
    MomT(1:4,1:4) = Mom(1:4,1:4)

          if( MomT(1,4)/mt.lt.1d-2 .or. (MomT(1:4,3).dot.MomT(1:4,4))/mt**2.lt.1d-2) return


          ExtParticles(1)%Mom(1:4) = dcmplx(MomT(1:4,1))
          ExtParticles(2)%Mom(1:4) = dcmplx(MomT(1:4,3))
          ExtParticles(3)%Mom(1:4) = dcmplx(MomT(1:4,4))
          ExtParticles(4)%Mom(1:4) = dcmplx(MomT(1:4,2))
                  HH(+1,+1) = 1d0
                  HH(-1,-1) = 1d0
                  HH(-1,+1) = 0d0
                  HH(+1,-1) = 0d0
          SqMatDip = (0d0,0d0)
          do hel1=1,-1,-2
          do hel2=1,-1,-2
          do hel4=1,-1,-1
          do hel3=1,-1,-2
                  call uSpi(ExtParticles(1)%Mom(1:4),ExtParticles(1)%Mass,Hel1,ExtParticles(1)%Pol(1:4))
                  call ubarSpi(ExtParticles(2)%Mom(1:4),ExtParticles(2)%Mass,Hel2,ExtParticles(2)%Pol(1:4))
                  call pol_mless(ExtParticles(3)%Mom(1:4),Hel3,ExtParticles(3)%Pol(1:4))
                  call pol_massSR(ExtParticles(4)%Mom(1:4),mw,Hel4,ExtParticles(4)%Pol(1:4))

                  call EvalTree2(TreeAmpsDip(1),AmplH(hel3))
        enddo
                  SqMatDip = SqMatDip + HH(+1,+1)*AmplH(+1)*dconjg(AmplH(+1)) &
                                      + HH(-1,-1)*AmplH(-1)*dconjg(AmplH(-1)) &
                                      + HH(-1,+1)*AmplH(-1)*dconjg(AmplH(+1)) &
                                      + HH(+1,-1)*AmplH(+1)*dconjg(AmplH(-1))
        enddo
        enddo
        enddo
!                                        4pi*alpha/2/sw2              (4pi*alpha_s)
        SqMatDip = SqMatDip * 1d0/2d0/3d0* 0.2063277518095031333d0  *  (1.633628179866692484d0)


!             DipColorCorr = +2d0/3d0
! !                           8*pi*alpha_s
!             Dipole_fi_qg = 3.26725635973d0 * SqMatDip * ( (1d0+z**2)/(1d0-z) - (Mom(1:4,i).dot.Mom(1:4,j))*mt**2/(Mom(1:4,k).dot.Mom(1:4,j))**2 ) /(Mom(1:4,i).dot.Mom(1:4,j))/2d0  * DipColorCorr


            IntDipole_fi_qg = -SqMatDip *  ( 2d0*dlog(alpha_DKfi)**2 -dlog(alpha_DKfi)*(-7d0/2d0+4d0*alpha_DKfi-alpha_DKfi**2/2d0) &
                               - 2d0*(1d0-alpha_DKfi)*rsq/(1d0-rsq)*dlog(rsq/(1d0-alpha_DKfi+rsq*alpha_DKfi)) )





RETURN
END FUNCTION






FUNCTION IntDipole_fi_gg(Mom,s,ij)
use ModParameters
use ModKinematics
use ModMisc
use ModProcess
use ModAmplitudes
implicit none
real(8) :: IntDipole_fi_gg
real(8) :: Mom(1:4,1:4),MomT(1:4,1:4),MomB(1:4),Mass,pbDpg,ptDpg,ptDpb,omz,z,y,ymax,x,paux(1:4),MomA(1:4),Pia(1:4),Pia2,Ria,r2,f
real(8) :: mt,SingDepth,MadGraph_tree,MadGraph_tree2,mw,DipColorCorr=0d0,MomK(1:4),MomKTilde(1:4),rsq,prec(1:4)
real(8) :: dipole,dipole1,dipole2,dipole3,dipole4,dipole5,dipole6,dipole7
integer :: Pcol1,Pcol2,Steps,n,i,j,k,hel1,hel2,hel3,hel4,hel5,s,ij
complex(8) :: Ampl1,Ampl2,AmplH(-1:1),SqMat,SqMatDip,diag,offdiag,xm,xp,HH(-1:1,-1:1)
integer :: NumQuarks,NumGluons,NumTrees,SingType=-99
type(TreeProcess),save :: TreeAmpsDip(1:2)
integer :: NumParticles,iTree
type(Particle),save :: ExtParticles(1:5)
integer,parameter :: ff_=1,fi_=2,if_=3
logical,save :: init=.true.
real(8) :: L,L2,epinv,epinv2
real(8):: r,r4,r6,r8,r10,dec_fi_gg
real(8) :: x1,x2,x3,x4,x5,x6,x7,x8,x9


  mt=m_top; mw=m_w

  if(init) then
          NumQuarks=2
          NumGluons=1
          NumParticles=NumQuarks+NumGluons+1
          allocate( TreeAmpsDip(1)%NumGlu(0:NumQuarks)     )
          allocate( TreeAmpsDip(1)%PartRef(1:NumParticles) )
          allocate( TreeAmpsDip(1)%PartType(1:NumParticles))
          allocate( TreeAmpsDip(1)%Quarks(1:NumQuarks)     )
          allocate( TreeAmpsDip(1)%Gluons(1:NumGluons)     )
          TreeAmpsDip(1)%NumPart=NumParticles
          TreeAmpsDip(1)%NumQua=NumQuarks
          TreeAmpsDip(1)%NumGlu(0) = NumGluons
          ExtParticles(1)%PartType = 5! top
          ExtParticles(1)%ExtRef   = 1
          ExtParticles(1)%Mass = mt
          ExtParticles(1)%Mass2= ExtParticles(1)%Mass**2
          ExtParticles(2)%PartType = 4! str
          ExtParticles(2)%ExtRef   = 2
          ExtParticles(2)%Mass = 0d0
          ExtParticles(2)%Mass2= ExtParticles(2)%Mass**2
          ExtParticles(3)%PartType = 10! glu
          ExtParticles(3)%ExtRef   = 3
          ExtParticles(3)%Mass = 0d0
          ExtParticles(3)%Mass2= ExtParticles(2)%Mass**2
          ExtParticles(4)%PartType = 13! w+
          ExtParticles(4)%ExtRef   = 4
          ExtParticles(4)%Mass = mw
          ExtParticles(4)%Mass2= ExtParticles(4)%Mass**2
          TreeAmpsDip(1)%PartRef(1:4) = (/1,2,3,4/)
          call LinkTreeParticles(TreeAmpsDip(1),ExtParticles(1:4))
          init=.false.
  endif

  IntDipole_fi_gg = 0d0

   prec(1:4) = +Mom(1:4,s) -Mom(1:4,ij)  ! assumes all are outgoing
   rsq = (prec.dot.prec)/mt**2
   r2 = rsq
   r = dsqrt(r2)
   MomT(1:4,1:4) = Mom(1:4,1:4)


          if( MomT(1,4)/mt.lt.1d-2 .or. (MomT(1:4,3).dot.MomT(1:4,4))/mt**2.lt.1d-2) return

          ExtParticles(1)%Mom(1:4) = dcmplx(MomT(1:4,1))
          ExtParticles(2)%Mom(1:4) = dcmplx(MomT(1:4,3))
          ExtParticles(3)%Mom(1:4) = dcmplx(MomT(1:4,4))
          ExtParticles(4)%Mom(1:4) = dcmplx(MomT(1:4,2))
                  HH(+1,+1) = 1d0
                  HH(-1,-1) = 1d0
                  HH(-1,+1) = 0d0
                  HH(+1,-1) = 0d0
          SqMatDip = (0d0,0d0)
          do hel1=1,-1,-2
          do hel2=1,-1,-2
          do hel4=1,-1,-1
          do hel3=1,-1,-2
                  call uSpi(ExtParticles(1)%Mom(1:4),ExtParticles(1)%Mass,Hel1,ExtParticles(1)%Pol(1:4))
                  call ubarSpi(ExtParticles(2)%Mom(1:4),ExtParticles(2)%Mass,Hel2,ExtParticles(2)%Pol(1:4))
                  call pol_mless(ExtParticles(3)%Mom(1:4),Hel3,ExtParticles(3)%Pol(1:4))
                  call pol_massSR(ExtParticles(4)%Mom(1:4),mw,Hel4,ExtParticles(4)%Pol(1:4))

                  call EvalTree2(TreeAmpsDip(1),AmplH(hel3))
        enddo
                  SqMatDip = SqMatDip + HH(+1,+1)*AmplH(+1)*dconjg(AmplH(+1)) &
                                      + HH(-1,-1)*AmplH(-1)*dconjg(AmplH(-1)) &
                                      + HH(-1,+1)*AmplH(-1)*dconjg(AmplH(+1)) &
                                      + HH(+1,-1)*AmplH(+1)*dconjg(AmplH(-1))
        enddo
        enddo
        enddo
!                                        4pi*alpha/2/sw2              (4pi*alpha_s)
        SqMatDip = SqMatDip * 1d0/2d0/3d0* 0.2063277518095031333d0  *  (1.633628179866692484d0)
!         SqMatDip = 1d0

    r4=r2**2
    r6=r4*r2
    r8=r6*r2
    r10=r8*r2
    x1=alpha_DKfi
    x2=x1*x1
    x3=x2*x1
    x4=x3*x1
    x5=x4*x1
    x6=x5*x1
    x7=x6*x1
    x8=x7*x1
    x9=x8*x1
    epinv = 0d0
    epinv2 = 0d0

      L = dlog(mt**2/MuRen**2)
      L2= L**2


         dec_fi_gg =  0.5d0*epinv*epinv2 &
            +(17d0/12d0-dlog(1d0-r2)-1d0/2d0*L)*epinv &!   maybe the log terms are wrong by a factor of 1/2
            +DLi2(1d0-r2)-5d0/12d0*Pi**2+1d0/6d0*r2*(-48d0*r6*x1 &
            +12d0*r8*x1+12d0*x1-36d0*x1*r2+60d0*r4*x1-6d0*r4*x3 &
            +15d0*x2*r4+65d0*r4-6d0*x3*r2+3d0*x2*r2-45d0*r2 &
            +11d0*r8-43d0*r6+12d0)/(-1d0+r)**5/(r+1d0)**5*dlog(r) &
            -dlog(x1)**2+1d0/12d0*(-1d0+x1)*(2d0*x2-x1+23d0)*dlog(x1) &
            +dlog(1d0-r2)**2+(L-17d0/6d0)*dlog(1d0-r2) &
            -1d0/4d0*r2*(-1d0+x1)*(3d0*r4*x1+23d0*r4-2d0*x2*r4 &
            -13d0*r2-2d0*x2*r2-x1*r2+4d0-16d0*r6 &
            +4d0*r8)/(-1d0+r)**5/(r+1d0)**5*dlog(1d0-x1+x1*r2) &
            -1d0/240d0*(-910d0+67d0*x2+1700d0*L*x1*r2+602d0*r6*x2-60d0*L2*r10*x1 &
            -393d0*r8*x3+1062d0*r6*x3+9520d0*r4*x1-263d0*r8*x2-9520d0*r6*x1 &
            +4630d0*r8*x1-4510d0*x1*r2+340d0*L*r10*x1+8d0*x9*r10+40d0*x9*r2 &
        -80d0*x9*r4+80d0*x9*r6-340d0*L*x1-910d0*r10*x1+60d0*L2*x1-40d0*x9*r8 &
        -300d0*L2*x1*r2-3400d0*L*r4*x1+600d0*L2*r4*x1+3400d0*L*r6*x1-600d0*L2*r6*x1 &
        -1700d0*L*r8*x1+300d0*L2*r8*x1+348d0*x8*r4-42d0*x8*r10+198d0*x8*r8 &
        -162d0*x8*r2+272d0*x2*r4-358d0*x2*r2-372d0*x8*r6+2040d0*L*r4 &
        +115d0*x7*r10+910d0*x1+363d0*x7*r2-802d0*x7*r4+898d0*x7*r6-8d0*x9+30d0*x8 &
        -67d0*x7+141d0*x6-241d0*x5+283d0*x4 &
        -507d0*x7*r8-5540d0*r4+1520d0*x6*r4-216d0*x6*r10+3680d0*r2+915d0*x6*r8 &
        -730d0*x6*r2-1630d0*x6*r6+240d0*L2*r2+340d0*L+305d0*x5*r10-2516d0*x5*r4 &
        +1219d0*x5*r2+3680d0*r6-910d0*r8-205d0*x3+2524d0*x5*r6 &
        -1291d0*x5*r8-360d0*L2*r4-1360d0*L*r2-1360d0*x4*r2-1360d0*L*r6+3050d0*r4*x4 &
        -278d0*x4*r10-2650d0*r6*x4+240d0*L2*r6+1195d0*r8*x4+40d0*r10*x2 &
        -2078d0*r4*x3+917d0*x3*r2+340d0*L*r8+97d0*x3*r10-60d0*L2*r8 &
        -60d0*L2)/(1d0-x1+x1*r2)/(-1d0+r)**4/(r+1d0)**4



            IntDipole_fi_gg = -SqMatDip *dec_fi_gg



!             IntDipole_fi_gg = -SqMatDip *  ( 2d0*r2*(-1d0+x1)/(-1d0+r)/(r+1d0)*dlog(r) +dlog(x1)**2  -r2*(-1d0+x1)/(-1d0+r)/(r+1d0)*dlog(1d0-x1+x1*r2)  )
!             IntDipole_fi_gg = -SqMatDip *  ( dlog(x1)**2 + (1d0-x1)*dlog(x1)  )

!             IntDipole_fi_gg = -SqMatDip *  ( 2d0*r2*(-1d0+x1)/(-1d0+r)/(r+1d0)*dlog(r) +dlog(x1)**2  -r2*(-1d0+x1)/(-1d0+r)/(r+1d0)*dlog(1d0-x1+x1*r2)  - (  dlog(x1)**2 + (1d0-x1)*dlog(x1) ) )

RETURN
END FUNCTION


















SUBROUTINE LinkTreeParticles(TheTreeAmp,TheParticles) ! NEW
use ModProcess
implicit none
type(TreeProcess) :: TheTreeAmp
type(Particle),target :: TheParticles(1:4)
integer :: iPart,PartRef,PartType,ig,iq,ib,NPart,counter,QuarkPos(1:6)
logical IsAQuark,IsABoson


! print *, size(TheParticles(:))
! pause

           counter = 0
           do NPart=1,TheTreeAmp%NumPart
                 TheTreeAmp%PartType(NPart) = TheParticles( TheTreeAmp%PartRef(NPart) )%PartType
                 if( IsAQuark( TheTreeAmp%PartType(NPart) ) ) then  ! not a gluon and not a boson
!                      TheTreeAmp%NumQua = TheTreeAmp%NumQua + 1
                    counter = counter + 1
                    QuarkPos(counter) = NPart
                 endif
           enddo

!           set number of gluons between quark lines
           if( IsAQuark( TheTreeAmp%PartType(1) ) ) then ! not a gluon and not a boson
             if( TheTreeAmp%NumQua .eq. 2 ) then
                   TheTreeAmp%NumGlu(1) = QuarkPos(2) - QuarkPos(1) - 1
                   TheTreeAmp%NumGlu(2) = TheTreeAmp%NumQua+TheTreeAmp%NumGlu(0) - QuarkPos(2)
             endif
             if( TheTreeAmp%NumQua .eq. 4 ) then
                   TheTreeAmp%NumGlu(1) = QuarkPos(2) - QuarkPos(1) - 1
                   TheTreeAmp%NumGlu(2) = QuarkPos(3) - QuarkPos(2) - 1
                   TheTreeAmp%NumGlu(3) = QuarkPos(4) - QuarkPos(3) - 1
                   TheTreeAmp%NumGlu(4) = TheTreeAmp%NumQua+TheTreeAmp%NumGlu(0) - QuarkPos(4)
             endif
             if( TheTreeAmp%NumQua .eq. 6 ) then
                   TheTreeAmp%NumGlu(1) = QuarkPos(2) - QuarkPos(1) - 1
                   TheTreeAmp%NumGlu(2) = QuarkPos(3) - QuarkPos(2) - 1
                   TheTreeAmp%NumGlu(3) = QuarkPos(4) - QuarkPos(3) - 1
                   TheTreeAmp%NumGlu(4) = QuarkPos(5) - QuarkPos(4) - 1
                   TheTreeAmp%NumGlu(5) = QuarkPos(6) - QuarkPos(5) - 1
                   TheTreeAmp%NumGlu(6) = TheTreeAmp%NumQua+TheTreeAmp%NumGlu(0) - QuarkPos(6)
             endif

           elseif( TheTreeAmp%PartType(1).eq.10 ) then ! is a gluon
             if( TheTreeAmp%NumQua .eq. 2 ) then
                   TheTreeAmp%NumGlu(1) = QuarkPos(1) - 2
                   TheTreeAmp%NumGlu(2) = QuarkPos(2) - QuarkPos(1) - 1
                   TheTreeAmp%NumGlu(3) = TheTreeAmp%NumQua+TheTreeAmp%NumGlu(0) - QuarkPos(2)
             endif
             if( TheTreeAmp%NumQua .eq. 4 ) then
                   TheTreeAmp%NumGlu(1) = QuarkPos(1) - 2
                   TheTreeAmp%NumGlu(2) = QuarkPos(2) - QuarkPos(1) - 1
                   TheTreeAmp%NumGlu(3) = QuarkPos(3) - QuarkPos(2) - 1
                   TheTreeAmp%NumGlu(4) = QuarkPos(4) - QuarkPos(3) - 1
                   TheTreeAmp%NumGlu(5) = TheTreeAmp%NumQua+TheTreeAmp%NumGlu(0) - QuarkPos(4)
             endif
             if( TheTreeAmp%NumQua .eq. 6 ) then
                   TheTreeAmp%NumGlu(1) = QuarkPos(1) - 2
                   TheTreeAmp%NumGlu(2) = QuarkPos(2) - QuarkPos(1) - 1
                   TheTreeAmp%NumGlu(3) = QuarkPos(3) - QuarkPos(2) - 1
                   TheTreeAmp%NumGlu(4) = QuarkPos(4) - QuarkPos(3) - 1
                   TheTreeAmp%NumGlu(5) = QuarkPos(5) - QuarkPos(4) - 1
                   TheTreeAmp%NumGlu(6) = QuarkPos(6) - QuarkPos(5) - 1
                   TheTreeAmp%NumGlu(7) = TheTreeAmp%NumQua+TheTreeAmp%NumGlu(0) - QuarkPos(6)
             endif
           endif


  ig=0; iq=0; ib=0;
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
     elseif( IsAQuark(PartType) ) then ! PartType==Quark
           iq=iq+1
           TheTreeAmp%Quarks(iq)%PartType => TheParticles(PartRef)%PartType
           TheTreeAmp%Quarks(iq)%ExtRef   => TheParticles(PartRef)%ExtRef
           TheTreeAmp%Quarks(iq)%Mass     => TheParticles(PartRef)%Mass
           TheTreeAmp%Quarks(iq)%Mass2    => TheParticles(PartRef)%Mass2
           TheTreeAmp%Quarks(iq)%Helicity => TheParticles(PartRef)%Helicity
           TheTreeAmp%Quarks(iq)%Mom      => TheParticles(PartRef)%Mom
           TheTreeAmp%Quarks(iq)%Pol      => TheParticles(PartRef)%Pol
           if( PartType.ne.TheParticles(PartRef)%PartType ) print *,"Error2 in LinkTreeParticles"
     elseif( IsABoson(PartType) ) then  ! PartType==Boson
           ib=ib+1
           TheTreeAmp%Boson%PartType => TheParticles(PartRef)%PartType
           TheTreeAmp%Boson%ExtRef   => TheParticles(PartRef)%ExtRef
           TheTreeAmp%Boson%Mass     => TheParticles(PartRef)%Mass
           TheTreeAmp%Boson%Mass2    => TheParticles(PartRef)%Mass2
           TheTreeAmp%Boson%Helicity => TheParticles(PartRef)%Helicity
           TheTreeAmp%Boson%Mom      => TheParticles(PartRef)%Mom
           TheTreeAmp%Boson%Pol      => TheParticles(PartRef)%Pol
           if( PartType.ne.TheParticles(PartRef)%PartType ) print *,"Error2 in LinkTreeParticles"
     endif
  enddo
  if( ig.ne.TheTreeAmp%NumGlu(0) .OR. iq.ne.TheTreeAmp%NumQua .OR. ib.ne.TheTreeAmp%NumPart-TheTreeAmp%NumGlu(0)-TheTreeAmp%NumQua) print *,"Error3 in LinkTreeParticles"

return
END SUBROUTINE




FUNCTION IsABoson(Type)
implicit none
logical IsABoson
integer :: Type

   if( abs(Type).eq.11 .or. abs(Type).eq.12 .or. abs(Type).eq.13 ) then
      IsABoson = .true.
   else
      IsABoson = .false.
   endif
END FUNCTION






FUNCTION IsAQuark(Type)
implicit none
logical IsAQuark
integer :: Type

   if( abs(Type).le.6 .and. Type.ne.0 ) then
      IsAQuark = .true.
   else
      IsAQuark = .false.
   endif
END FUNCTION



SUBROUTINE DipoleTrafox(MomT,MomK,MomIn,MomOut)
use ModMisc
implicit none
real(8),intent(in) :: MomT(1:4),MomK(1:4),MomIn(1:4)
real(8),intent(out) :: MomOut(1:4)
real(8) :: xh,sinhyp,coshyp,a,b
real(8) :: mt2,mK2,s_tP,s_KP,s_tK


    mt2  = MomT(1:4).dot.MomT(1:4)
    mK2  = MomK(1:4).dot.MomK(1:4)
    s_tP = MomT(1:4).dot.MomIn(1:4)
    s_KP = MomK(1:4).dot.MomIn(1:4)
    s_tK = MomT(1:4).dot.MomK(1:4)

    xh = s_tK**2 - mt2 * mK2
    sinhyp = 0.5d0/(mt2 * mK2)*( -(mt2-mK2)*s_tK + (mt2+mK2)*dsqrt(xh) )
    coshyp = 0.5d0/(mt2 * mK2)*( +(mt2+mK2)*s_tK - (mt2-mK2)*dsqrt(xh) )

    a = sinhyp/dsqrt(xh)
    b = (coshyp-1d0)/xh

    MomOut(1:4) = MomIn(1:4) - a*( s_tP*MomK(1:4) - s_KP*MomT(1:4) )  &
                  + b*( s_tK*( s_tP*MomK(1:4) + s_KP*MomT(1:4) )  -  mK2*s_tP*MomT(1:4) - mt2*s_KP*MomK(1:4))


RETURN
END SUBROUTINE





END MODULE
