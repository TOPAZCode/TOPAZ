program DKDipoles!     mat.el. for t->b w g g    and dipoles
use ModAmplitudes!     ff dipoles are taken 1:1 from CDST, if+fi dipoles are modified
use ModProcess
use modMisc
implicit none
real(8) :: Mom(1:4,1:5),MomT(1:4,1:4),MomB(1:4),Mass,pbDpg,ptDpg,ptDpb,omz,rsq,z,y,x,paux(1:4),MomA(1:4),Pia(1:4),Pia2,Ria,r2,f
real(8) :: SingDepth,MadGraph_tree,MadGraph_tree2
real(8),parameter :: mt=172d0,mw=80.419d0
real(8) :: DipColorCorr=0d0,MomK(1:4),MomKTilde(1:4),MomLep(1:4,2),MomLepT(1:4,2),PsWgt
real(8) :: dipole,dipole1,dipole2,dipole3,dipole4,dipole5,dipole6,dipole7
integer :: Pcol1,Pcol2,Steps,n,i,j,k,hel1,hel2,hel3,hel4,hel5
complex(8) :: Ampl1,Ampl2,AmplH(-1:1),SqMat,SqMatDip,diag,offdiag,xm,xp,HH(-1:1,-1:1),pol(4,8),spi(1:4),barSpi(1:4)
integer :: NumQuarks,NumGluons,NumTrees,SingType=-99
type(TreeProcess) :: TreeAmpsDip(1:2)
integer :: NumParticles,iTree
type(Particle) :: ExtParticles(1:5)
integer,parameter :: ff_=1,fi_=2,if_=3
integer,  parameter :: WDK_Lep=1
integer,  parameter :: WDK_LepPho=2
integer,  parameter :: WDK_Had=3
integer,  parameter :: WDK_HadPho=4
integer,  parameter :: WDK_LO=0
integer :: Wp_DKmode = WDK_Lep
integer :: Wm_DKmode = WDK_Lep
integer,parameter :: WPlus=+1, WMinus=-1
real(8), parameter :: GF = (1.16639d-5)
real(8),parameter :: g2_weak = 4d0*dsqrt(2d0)*mW**2*GF
real(8),parameter :: g_weak = dsqrt(g2_weak)
real(8), parameter :: Ga_W    = 2.14d0
real(8),parameter :: NWAFactor_W   = 1d0/dsqrt(2d0*Ga_W*mW)
complex(8),parameter :: WProp = (0d0,-1d0)*NWAFactor_W
logical :: WDecay=.true.





  print *, "cur_f_2fW in Mod_Amplitudes needs to be called in DIRAC convention!"; pause

  call coupsm(0)

  NumQuarks=2
  NumGluons=2
  NumParticles=NumQuarks+NumGluons+1
  TreeAmpsDip(1)%NumPart=NumParticles
  TreeAmpsDip(1)%NumQua=NumQuarks
  allocate( TreeAmpsDip(1)%NumGlu(0:NumQuarks)     )
  allocate( TreeAmpsDip(1)%PartRef(1:NumParticles) )
  allocate( TreeAmpsDip(1)%PartType(1:NumParticles))
  allocate( TreeAmpsDip(1)%Quarks(1:NumQuarks)     )
  allocate( TreeAmpsDip(1)%Gluons(1:NumGluons)     )
  TreeAmpsDip(1)%NumGlu(0) = NumGluons

  TreeAmpsDip(2)%NumPart=NumParticles
  TreeAmpsDip(2)%NumQua=NumQuarks
  allocate( TreeAmpsDip(2)%NumGlu(0:NumQuarks)     )
  allocate( TreeAmpsDip(2)%PartRef(1:NumParticles) )
  allocate( TreeAmpsDip(2)%PartType(1:NumParticles))
  allocate( TreeAmpsDip(2)%Quarks(1:NumQuarks)     )
  allocate( TreeAmpsDip(2)%Gluons(1:NumGluons)     )
  TreeAmpsDip(2)%NumGlu(0) = NumGluons

      do n=1,12

          NumQuarks=2
          NumGluons=2
          NumParticles=NumQuarks+NumGluons+1
          TreeAmpsDip(1)%NumPart=NumParticles
          TreeAmpsDip(1)%NumQua=NumQuarks
          TreeAmpsDip(1)%NumGlu(0) = NumGluons

          TreeAmpsDip(2)%NumPart=NumParticles
          TreeAmpsDip(2)%NumQua=NumQuarks
          TreeAmpsDip(2)%NumGlu(0) = NumGluons


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
          ExtParticles(3)%Mass2= ExtParticles(3)%Mass**2

          ExtParticles(4)%PartType = 10! glu
          ExtParticles(4)%ExtRef   = 4
          ExtParticles(4)%Mass = 0d0
          ExtParticles(4)%Mass2= ExtParticles(4)%Mass**2

          ExtParticles(5)%PartType = 13! w+
          ExtParticles(5)%ExtRef   = 5
          ExtParticles(5)%Mass = mw
          ExtParticles(5)%Mass2= ExtParticles(5)%Mass**2


          TreeAmpsDip(1)%PartRef(1:5) = (/1,2,3,4,5/)
          call LinkTreeParticles(TreeAmpsDip(1),ExtParticles(1:5))

          TreeAmpsDip(2)%PartRef(1:5) = (/1,2,4,3,5/)
          call LinkTreeParticles(TreeAmpsDip(2),ExtParticles(1:5))

           Dipole = 0d0;   Dipole1 = 0d0;   Dipole2 = 0d0;   Dipole3 = 0d0;   Dipole4 = 0d0;   Dipole5 = 0d0;  Dipole6 = 0d0;   Dipole7 = 0d0;



          Pcol1= 4 -1
          Pcol2= 6 -1
          SingDepth = 1d-8
          Steps = 10

          print *,n
          Mom(1:4,1)= (/mt,0d0,0d0,0d0/)
          call gensing(4,mt,(/mw,0d0,0d0,0d0/),Mom(1:4,2:5),Pcol1,Pcol2,SingDepth,Steps)
          if(WDecay) then
              call genps(2,mW,(/0.3d0,0.7d0/),(/0d0,0d0/),MomLep(1:4,1:2),PSWgt)
              call boost(MomLep(1:4,1),Mom(1:4,2),mW)
              call boost(MomLep(1:4,2),Mom(1:4,2),mW)
          endif

! Mom(1:4,1)= (/ 172.000000000000000d0,   0.000000000000000d0,   0.000000000000000d0,   0.000000000000000d0/)
! Mom(1:4,3)= (/   16.376837591417118d0,  -0.310086510680035d0,  -4.737621420368479d0, -15.673531801378175d0/)
! Mom(1:4,4)= (/   49.675220278715372d0,  41.167455334141593d0, -11.436293882299141d0,  25.339283993199938d0/)
! Mom(1:4,5)= (/   12.438147982055812d0,   3.826780070951152d0,  -0.424910214473981d0, -11.827203000782808d0/)
! MomLep(1:4,1)= (/   60.749128709182514d0, -37.325425205746200d0,  44.683866802785431d0,  17.338434753056529d0/)
! MomLep(1:4,2)= (/   32.760665438629182d0,  -7.358723688666508d0, -28.085041285643833d0, -15.176983944095486d0/)
! Mom(1:4,2)=MomLep(1:4,1)+MomLep(1:4,2)
! WPol (-1.21903215073211,-4.27390903321014) (3.82029649068264,9.24485929840147) (4.36789460340088,0.253540309470071) (-7.30373719250250,4.27390903321014)
! Sq.Mat (7.18977319378407,0.000000000000000E+000)



! Mom(1:4,1)= (/  172.000000000000000d0,   0.000000000000000d0,   0.000000000000000d0,   0.000000000000000d0/)
! Mom(1:4,3)= (/   16.275110865860512d0,   5.004917102368888d0,  -0.568229874651678d0, -15.476018651578958d0/)
! Mom(1:4,4)= (/   48.464930011006778d0,  39.408441349307118d0, -12.864129926036123d0,  25.106540037131111d0/)
! Mom(1:4,5)= (/   12.438147982055812d0,   3.826780070951152d0,  -0.424910214473981d0, -11.827203000782808d0/)
! MomLep(1:4,1)= (/   60.720895867051482d0, -39.172116318631623d0,  43.030050399813391d0,  17.348984425793461d0/)
! MomLep(1:4,2)= (/   34.100915274025397d0,  -9.068022203995540d0, -29.172780384651606d0, -15.152302810562805d0/)
! Mom(1:4,2)=MomLep(1:4,1)+MomLep(1:4,2)
! Sq.Mat (2365890.30922929,0.000000000000000E+000)


          if( Mom(1,4)/mt .lt. 1d-2 ) stop
          if( Mom(1,5)/mt .lt. 1d-2 ) stop
          if( (Mom(1:4,4).dot.Mom(1:4,5))/mt**2 .lt. 1d-2 ) stop
          if( (Mom(1:4,3).dot.Mom(1:4,4))/mt**2 .lt. 1d-2 ) stop

          print *, "Sing",(Mom(1:4,3).dot.Mom(1:4,5))/mt**2

          ExtParticles(1)%Mom(1:4) = dcmplx(Mom(1:4,1))
          ExtParticles(2)%Mom(1:4) = dcmplx(Mom(1:4,3))
          ExtParticles(3)%Mom(1:4) = dcmplx(Mom(1:4,4))
          ExtParticles(4)%Mom(1:4) = dcmplx(Mom(1:4,5))
          ExtParticles(5)%Mom(1:4) = dcmplx(Mom(1:4,2))
          SqMat = (0d0,0d0)

          do hel1=1,-1,-2
          do hel2=1,-1,-2
          do hel3=1,-1,-2
          do hel4=1,-1,-2
          do hel5=1,-1,-1
                  call uSpi(ExtParticles(1)%Mom(1:4),ExtParticles(1)%Mass,Hel1,ExtParticles(1)%Pol(1:4))
                  call ubarSpi(ExtParticles(2)%Mom(1:4),ExtParticles(2)%Mass,Hel2,ExtParticles(2)%Pol(1:4))
                  call pol_mless(ExtParticles(3)%Mom(1:4),Hel3,ExtParticles(3)%Pol(1:4))
                  call pol_mless(ExtParticles(4)%Mom(1:4),Hel4,ExtParticles(4)%Pol(1:4))
                  if(WDecay) then
                       if(hel5.ne.1) cycle
                       call WPolVec2(WPlus,Wp_DKmode,MomLep(1:4,1:2),WDK_LO,ExtParticles(5)%Pol(1:4))
                       ExtParticles(5)%Pol(1:4)=ExtParticles(5)%Pol(1:4)  *WProp*g_weak/dsqrt(2d0)
                  else
                      call pol_massSR(ExtParticles(5)%Mom(1:4),mw,Hel5,ExtParticles(5)%Pol(1:4))
                      ExtParticles(5)%Pol(1:4)=ExtParticles(5)%Pol(1:4)*g_weak/dsqrt(2d0)
                  endif

                  call EvalTree2(TreeAmpsDip(1),Ampl1)
                  call EvalTree2(TreeAmpsDip(2),Ampl2)

                  SqMat = SqMat + 64d0/3d0*( Ampl1*dconjg(Ampl1) + Ampl2*dconjg(Ampl2) )  -  8d0/3d0*( Ampl1*dconjg(Ampl2) + Ampl2*dconjg(Ampl1) )
        enddo
        enddo
        enddo
        enddo
        enddo
!                                           (4pi*alpha_s)^2           symm.fact
        SqMat = SqMat * 1d0/2d0/3d0 *  (1.633628179866692484d0)**2   * 0.5d0


!
!         call ST_BWPGG((/Mom(1:4,1),Mom(1:4,3),Mom(1:4,2),Mom(1:4,4),Mom(1:4,5)/),MadGraph_tree)
!         print *, "My Amp.:", SqMat
!         print *, "MadGraph2:", MadGraph_tree
!         print *, "ratio:", MadGraph_tree/SqMat
! pause;cycle
12 continue



! --------------------------------- dipole calculation -------------------------------------------------
! there are the following color correlations:
! T_t * T_b =  1/2*CA -CF = 1/6
! T_t * T_g = -1/2*CA = -3/2
! T_b * T_g = -1/2*CA = -3/2
          NumQuarks=2
          NumGluons=1
          NumParticles=NumQuarks+NumGluons+1
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





!        if( SingType.eq.ff_ ) then!    top w bot glu glu
          i=3; j=5; k=4;
          y = (Mom(1:4,i).dot.Mom(1:4,j))/((Mom(1:4,i).dot.Mom(1:4,j))+(Mom(1:4,i).dot.Mom(1:4,k))+(Mom(1:4,j).dot.Mom(1:4,k)))
          z = (Mom(1:4,i).dot.Mom(1:4,k))/((Mom(1:4,i).dot.Mom(1:4,k))+(Mom(1:4,k).dot.Mom(1:4,j)))
          MomT(1:4,1:4)=Mom(1:4,1:4)
          MomT(1:4,k) = 1d0/(1d0-y)*Mom(1:4,k)
          MomT(1:4,i) = Mom(1:4,i)+Mom(1:4,j)-y/(1d0-y)*Mom(1:4,k)
          if( MomT(1,4)/mt.lt.1d-2 .or. (MomT(1:4,3).dot.MomT(1:4,4))/mt**2.lt.1d-2) goto 13
!         endif

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
                  if(WDecay) then
                      call WPolVec2(WPlus,Wp_DKmode,MomLep(1:4,1:2),WDK_LO,ExtParticles(4)%Pol(1:4))
                      ExtParticles(4)%Pol(1:4)=ExtParticles(4)%Pol(1:4)*WProp*g_weak/dsqrt(2d0)
                      if( hel4.ne.1 ) ExtParticles(4)%Pol(1:4) = (0d0,0d0)
                  else
                      call pol_massSR(ExtParticles(4)%Mom(1:4),mw,Hel4,ExtParticles(4)%Pol(1:4))
                      ExtParticles(4)%Pol(1:4) = ExtParticles(4)%Pol(1:4)*g_weak/dsqrt(2d0)
                  endif


                  call EvalTree2(TreeAmpsDip(1),AmplH(hel3))
        enddo
                  SqMatDip = SqMatDip + HH(+1,+1)*AmplH(+1)*dconjg(AmplH(+1)) &
                                      + HH(-1,-1)*AmplH(-1)*dconjg(AmplH(-1)) &
                                      + HH(-1,+1)*AmplH(-1)*dconjg(AmplH(+1)) &
                                      + HH(+1,-1)*AmplH(+1)*dconjg(AmplH(-1))
        enddo
        enddo
        enddo
!                                            (4pi*alpha_s)          color
        SqMatDip = SqMatDip * 1d0/2d0/3d0  *  (1.633628179866692484d0) * 8d0
!         call ST_BWPG((/MomT(1:4,1),MomT(1:4,3),MomT(1:4,2),MomT(1:4,4)/),MadGraph_tree)
!         print *, "My Amp.:", SqMatDip
!         print *, "MadGraph2:", MadGraph_tree
!         print *, "MG ratio:", MadGraph_tree/SqMatDip


            DipColorCorr = -3d0/2d0
!                      8*pi*alpha_s
            dipole1 = 3.26725635973d0 * SqMatDip * ( 2d0/(1d0-z*(1d0-y)) -1d0-z ) /(Mom(1:4,i).dot.Mom(1:4,j))/4d0  * DipColorCorr !  CS f-f dipole for qg

! print *, "SqMatDip1",SqMatDip,dipole1/SqMatDip
! print *, "MomTd",MomT(1:4,3)
! print *, "MomTd",MomT(1:4,4)
! print *, "MomTd",MomT(1:4,2)


13 continue


!        if( SingType.eq.ff_ ) then!    top w bot glu glu
          i=3; j=4; k=5;
          y = (Mom(1:4,i).dot.Mom(1:4,j))/((Mom(1:4,i).dot.Mom(1:4,j))+(Mom(1:4,i).dot.Mom(1:4,k))+(Mom(1:4,j).dot.Mom(1:4,k)))
          z = (Mom(1:4,i).dot.Mom(1:4,k))/((Mom(1:4,i).dot.Mom(1:4,k))+(Mom(1:4,k).dot.Mom(1:4,j)))
          MomT(1:4,1:4)=Mom(1:4,1:4)
          MomT(1:4,k-1) = 1d0/(1d0-y)*Mom(1:4,k)
          MomT(1:4,i) = Mom(1:4,i)+Mom(1:4,j)-y/(1d0-y)*Mom(1:4,k)
          if( MomT(1,4)/mt.lt.1d-2 .or. (MomT(1:4,3).dot.MomT(1:4,4))/mt**2.lt.1d-2) goto 14
!         endif

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
                  if(WDecay) then
                      call WPolVec2(WPlus,Wp_DKmode,MomLep(1:4,1:2),WDK_LO,ExtParticles(4)%Pol(1:4))
                      ExtParticles(4)%Pol(1:4)=ExtParticles(4)%Pol(1:4)*WProp*g_weak/dsqrt(2d0)
		      if( hel4.ne.1 ) ExtParticles(4)%Pol(1:4) = (0d0,0d0)
                  else
                      call pol_massSR(ExtParticles(4)%Mom(1:4),mw,Hel4,ExtParticles(4)%Pol(1:4))
                      ExtParticles(4)%Pol(1:4) = ExtParticles(4)%Pol(1:4)*g_weak/dsqrt(2d0)
                  endif


                  call EvalTree2(TreeAmpsDip(1),AmplH(hel3))
        enddo
                  SqMatDip = SqMatDip + HH(+1,+1)*AmplH(+1)*dconjg(AmplH(+1)) &
                                      + HH(-1,-1)*AmplH(-1)*dconjg(AmplH(-1)) &
                                      + HH(-1,+1)*AmplH(-1)*dconjg(AmplH(+1)) &
                                      + HH(+1,-1)*AmplH(+1)*dconjg(AmplH(-1))
        enddo
        enddo
        enddo
!                                               (4pi*alpha_s)
        SqMatDip = SqMatDip * 1d0/2d0/3d0  *  (1.633628179866692484d0) * 8d0
!         call ST_BWPG((/MomT(1:4,1),MomT(1:4,3),MomT(1:4,2),MomT(1:4,4)/),MadGraph_tree)
!         print *, "My Amp.:", SqMatDip
!         print *, "MadGraph2:", MadGraph_tree
!         print *, "MG ratio:", MadGraph_tree/SqMatDip

            DipColorCorr = -3d0/2d0
            dipole2 = 3.26725635973d0 * SqMatDip * ( 2d0/(1d0-z*(1d0-y)) -1d0-z ) /(Mom(1:4,i).dot.Mom(1:4,j))/4d0  * DipColorCorr !  CS f-f dipole for qg



14 continue


!        if( SingType.eq.ff_ ) then!    top w bot glu glu
          i=4; j=5; k=3;
          y = (Mom(1:4,i).dot.Mom(1:4,j))/((Mom(1:4,i).dot.Mom(1:4,j))+(Mom(1:4,i).dot.Mom(1:4,k))+(Mom(1:4,j).dot.Mom(1:4,k)))
          z = (Mom(1:4,i).dot.Mom(1:4,k))/((Mom(1:4,i).dot.Mom(1:4,k))+(Mom(1:4,k).dot.Mom(1:4,j)))
          MomT(1:4,1:4)=Mom(1:4,1:4)
          MomT(1:4,k) = 1d0/(1d0-y)*Mom(1:4,k)
          MomT(1:4,i) = Mom(1:4,i)+Mom(1:4,j)-y/(1d0-y)*Mom(1:4,k)
          if( MomT(1,4)/mt.lt.1d-2 .or. (MomT(1:4,3).dot.MomT(1:4,4))/mt**2.lt.1d-2) goto 15
!         endif
          ExtParticles(1)%Mom(1:4) = dcmplx(MomT(1:4,1))
          ExtParticles(2)%Mom(1:4) = dcmplx(MomT(1:4,3))
          ExtParticles(3)%Mom(1:4) = dcmplx(MomT(1:4,4))
          ExtParticles(4)%Mom(1:4) = dcmplx(MomT(1:4,2))
                  diag= 2d0*(1d0/(1d0-z*(1d0-y)) + 1d0/(1d0-(1d0-z)*(1d0-y)) - 2d0 )
                  offdiag = 2d0/(Mom(1:4,i).dot.Mom(1:4,j))
                  paux(1:4) = z*Mom(1:4,i) - (1d0-z)*Mom(1:4,j)
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
                  if(WDecay) then
                      call WPolVec2(WPlus,Wp_DKmode,MomLep(1:4,1:2),WDK_LO,ExtParticles(4)%Pol(1:4))
                      ExtParticles(4)%Pol(1:4)=ExtParticles(4)%Pol(1:4)*WProp*g_weak/dsqrt(2d0)
		      if( hel4.ne.1 ) ExtParticles(4)%Pol(1:4) = (0d0,0d0)
                  else
                      call pol_massSR(ExtParticles(4)%Mom(1:4),mw,Hel4,ExtParticles(4)%Pol(1:4))
                      ExtParticles(4)%Pol(1:4) = ExtParticles(4)%Pol(1:4)*g_weak/dsqrt(2d0)
                  endif


                  call EvalTree2(TreeAmpsDip(1),AmplH(hel3))
        enddo
                  SqMatDip = SqMatDip + HH(+1,+1)*AmplH(+1)*dconjg(AmplH(+1)) &
                                      + HH(-1,-1)*AmplH(-1)*dconjg(AmplH(-1)) &
                                      + HH(-1,+1)*AmplH(-1)*dconjg(AmplH(+1)) &
                                      + HH(+1,-1)*AmplH(+1)*dconjg(AmplH(-1))
        enddo
        enddo
        enddo
!                                               (4pi*alpha_s)
        SqMatDip = SqMatDip * 1d0/2d0/3d0  *  (1.633628179866692484d0) * 8d0


            DipColorCorr = -3d0/2d0
            dipole3 = 3.26725635973d0 * SqMatDip /(Mom(1:4,i).dot.Mom(1:4,j))/4d0 *DipColorCorr !  CS f-f dipole for gg



15 continue









!        if( SingType.eq.fi_ ) then
          i=3; j=5; k=1;
!           x = 1d0 - (Mom(1:4,i).dot.Mom(1:4,j))/((Mom(1:4,k).dot.Mom(1:4,i))+(Mom(1:4,k).dot.Mom(1:4,j)))
!           z = (Mom(1:4,i).dot.Mom(1:4,k))/((Mom(1:4,i).dot.Mom(1:4,k))+(Mom(1:4,k).dot.Mom(1:4,j)))
          z = 1d0-(    (Mom(1:4,k).dot.Mom(1:4,j))/((Mom(1:4,i).dot.Mom(1:4,k))+(Mom(1:4,k).dot.Mom(1:4,j))-(Mom(1:4,i).dot.Mom(1:4,j)))     )


!           Pia(1:4) = Mom(1:4,i)+Mom(1:4,j)-Mom(1:4,k)
!           Pia2 = Pia(1:4).dot.Pia(1:4)
!           MomT(1:4,1:4)=Mom(1:4,1:4)
!           MomT(1:4,i) = dsqrt(lambda(Pia2,mt**2,0d0))/dsqrt(lambda(2d0*(Mom(1:4,i).dot.Mom(1:4,j)),Pia2,mt**2)) * (Mom(1:4,k)-(Pia(1:4).dot.Mom(1:4,k))/Pia2*Pia(1:4))  &
!                       + (Pia2-mt**2)/2d0/Pia2*Pia(1:4)
!           MomT(1:4,k) = MomT(1:4,i) - Pia(1:4)


!        new mapping
         MomK(1:4) = Mom(1:4,2)+Mom(1:4,4)
         call DipoleTrafox(Mom(1:4,1),MomK(1:4),Mom(1:4,2),MomT(1:4,2))
         call DipoleTrafox(Mom(1:4,1),MomK(1:4),MomLep(1:4,1),MomLepT(1:4,1))
         call DipoleTrafox(Mom(1:4,1),MomK(1:4),MomLep(1:4,2),MomLepT(1:4,2))
         call DipoleTrafox(Mom(1:4,1),MomK(1:4),Mom(1:4,4),MomT(1:4,4))
         MomT(1:4,3) = Mom(1:4,1) - MomT(1:4,2) - MomT(1:4,4)
!          print *, Mom(1:4,1).dot.Mom(1:4,1)
!          print *, Mom(1:4,2).dot.Mom(1:4,2)
!          print *, Mom(1:4,3).dot.Mom(1:4,3)
!          print *, Mom(1:4,4).dot.Mom(1:4,4)
!          print *, Mom(1:4,1) - MomT(1:4,2) - MomT(1:4,3) - MomT(1:4,4)



          if( MomT(1,4)/mt.lt.1d-2 .or. (MomT(1:4,3).dot.MomT(1:4,4))/mt**2.lt.1d-2) goto 16
!         endif
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
                  if(WDecay) then
                      call WPolVec2(WPlus,Wp_DKmode,MomLep(1:4,1:2),WDK_LO,ExtParticles(4)%Pol(1:4))
                      ExtParticles(4)%Pol(1:4)=ExtParticles(4)%Pol(1:4)*WProp*g_weak/dsqrt(2d0)
		      if( hel4.ne.1 ) ExtParticles(4)%Pol(1:4) = (0d0,0d0)
                  else
                      call pol_massSR(ExtParticles(4)%Mom(1:4),mw,Hel4,ExtParticles(4)%Pol(1:4))
                      ExtParticles(4)%Pol(1:4) = ExtParticles(4)%Pol(1:4)*g_weak/dsqrt(2d0)
                  endif

                  call EvalTree2(TreeAmpsDip(1),AmplH(hel3))
        enddo
                  SqMatDip = SqMatDip + HH(+1,+1)*AmplH(+1)*dconjg(AmplH(+1)) &
                                      + HH(-1,-1)*AmplH(-1)*dconjg(AmplH(-1)) &
                                      + HH(-1,+1)*AmplH(-1)*dconjg(AmplH(+1)) &
                                      + HH(+1,-1)*AmplH(+1)*dconjg(AmplH(-1))
        enddo
        enddo
        enddo
!                                         (4pi*alpha_s)
        SqMatDip = SqMatDip * 1d0/2d0/3d0  *  (1.633628179866692484d0) * 8d0
!         call ST_BWPG((/MomT(1:4,1),MomT(1:4,3),MomT(1:4,2),MomT(1:4,4)/),MadGraph_tree)
!         print *, "My Amp.:", SqMatDip
!         print *, "MadGraph2:", MadGraph_tree
!         print *, "MG ratio:", MadGraph_tree/SqMatDip


            DipColorCorr = +1d0/6d0

!           new dipole which contains i-f dipole contribution
            Dipole4 = 3.26725635973d0 * SqMatDip * ( (1d0+z**2)/(1d0-z) - (Mom(1:4,i).dot.Mom(1:4,j))*mt**2/(Mom(1:4,k).dot.Mom(1:4,j))**2 ) /(Mom(1:4,i).dot.Mom(1:4,j))/4d0  * DipColorCorr

! print *, "SqMatDip2",SqMatDip,Dipole4/SqMatDip
! print *, "MomTd",MomT(1:4,3)
! print *, "MomTd",MomT(1:4,4)
! print *, "MomTd",MomT(1:4,2)


16 continue


!        if( SingType.eq.fi_ ) then
          i=3; j=4; k=1;
!           x = 1d0 - (Mom(1:4,i).dot.Mom(1:4,j))/((Mom(1:4,k).dot.Mom(1:4,i))+(Mom(1:4,k).dot.Mom(1:4,j)))
!           z = (Mom(1:4,i).dot.Mom(1:4,k))/((Mom(1:4,i).dot.Mom(1:4,k))+(Mom(1:4,k).dot.Mom(1:4,j)))
          z = 1d0-(    (Mom(1:4,k).dot.Mom(1:4,j))/((Mom(1:4,i).dot.Mom(1:4,k))+(Mom(1:4,k).dot.Mom(1:4,j))-(Mom(1:4,i).dot.Mom(1:4,j)))     )

!           Pia(1:4) = Mom(1:4,i)+Mom(1:4,j)-Mom(1:4,k)
!           Pia2 = Pia(1:4).dot.Pia(1:4)
!           MomT(1:4,1:4)=Mom(1:4,1:4);  MomT(1:4,4)=Mom(1:4,5);
!           MomT(1:4,i) = dsqrt(lambda(Pia2,mt**2,0d0))/dsqrt(lambda(2d0*(Mom(1:4,i).dot.Mom(1:4,j)),Pia2,mt**2)) * (Mom(1:4,k)-(Pia(1:4).dot.Mom(1:4,k))/Pia2*Pia(1:4))  &
!                       + (Pia2-mt**2)/2d0/Pia2*Pia(1:4)
!           MomT(1:4,k) = MomT(1:4,i) - Pia(1:4)


!        new mapping
         MomK(1:4) = Mom(1:4,2)+Mom(1:4,5)
         call DipoleTrafox(Mom(1:4,1),MomK(1:4),Mom(1:4,2),MomT(1:4,2))
         call DipoleTrafox(Mom(1:4,1),MomK(1:4),MomLep(1:4,1),MomLepT(1:4,1))
         call DipoleTrafox(Mom(1:4,1),MomK(1:4),MomLep(1:4,2),MomLepT(1:4,2))
         call DipoleTrafox(Mom(1:4,1),MomK(1:4),Mom(1:4,5),MomT(1:4,4))
         MomT(1:4,3) = Mom(1:4,1) - MomT(1:4,2) - MomT(1:4,4)
!          print *, Mom(1:4,1).dot.Mom(1:4,1)
!          print *, Mom(1:4,2).dot.Mom(1:4,2)
!          print *, Mom(1:4,3).dot.Mom(1:4,3)
!          print *, Mom(1:4,4).dot.Mom(1:4,4)
!          print *, Mom(1:4,1) - MomT(1:4,2) - MomT(1:4,3) - MomT(1:4,4)


          if( MomT(1,4)/mt.lt.1d-2 .or. (MomT(1:4,3).dot.MomT(1:4,4))/mt**2.lt.1d-2) goto 17
!         endif
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
                  if(WDecay) then
                      call WPolVec2(WPlus,Wp_DKmode,MomLep(1:4,1:2),WDK_LO,ExtParticles(4)%Pol(1:4))
                      ExtParticles(4)%Pol(1:4)=ExtParticles(4)%Pol(1:4)*WProp*g_weak/dsqrt(2d0)
		      if( hel4.ne.1 ) ExtParticles(4)%Pol(1:4) = (0d0,0d0)
                  else
                      call pol_massSR(ExtParticles(4)%Mom(1:4),mw,Hel4,ExtParticles(4)%Pol(1:4))
                      ExtParticles(4)%Pol(1:4) = ExtParticles(4)%Pol(1:4)*g_weak/dsqrt(2d0)
                  endif


                  call EvalTree2(TreeAmpsDip(1),AmplH(hel3))
        enddo
                  SqMatDip = SqMatDip + HH(+1,+1)*AmplH(+1)*dconjg(AmplH(+1)) &
                                      + HH(-1,-1)*AmplH(-1)*dconjg(AmplH(-1)) &
                                      + HH(-1,+1)*AmplH(-1)*dconjg(AmplH(+1)) &
                                      + HH(+1,-1)*AmplH(+1)*dconjg(AmplH(-1))
        enddo
        enddo
        enddo
!                                            (4pi*alpha_s)
        SqMatDip = SqMatDip * 1d0/2d0/3d0  *  (1.633628179866692484d0) * 8d0
!         call ST_BWPG((/MomT(1:4,1),MomT(1:4,3),MomT(1:4,2),MomT(1:4,4)/),MadGraph_tree)
!         print *, "My Amp.:", SqMatDip
!         print *, "MadGraph2:", MadGraph_tree
!         print *, "MG ratio:", MadGraph_tree/SqMatDip


            DipColorCorr = +1d0/6d0

!           new dipole which contains i-f dipole contribution
            Dipole5 = 3.26725635973d0 * SqMatDip *( (1d0+z**2)/(1d0-z) - (Mom(1:4,i).dot.Mom(1:4,j))*mt**2/(Mom(1:4,k).dot.Mom(1:4,j))**2 ) /(Mom(1:4,i).dot.Mom(1:4,j))/4d0  * DipColorCorr



17 continue



!        if( SingType.eq.fi_ ) then
          i=4; j=5; k=1;
!           x = 1d0 - (Mom(1:4,i).dot.Mom(1:4,j))/((Mom(1:4,k).dot.Mom(1:4,i))+(Mom(1:4,k).dot.Mom(1:4,j)))
!           z = (Mom(1:4,i).dot.Mom(1:4,k))/((Mom(1:4,i).dot.Mom(1:4,k))+(Mom(1:4,k).dot.Mom(1:4,j)))!  Dittmaier z
          z = 1d0-(    (Mom(1:4,k).dot.Mom(1:4,j))/((Mom(1:4,i).dot.Mom(1:4,k))+(Mom(1:4,k).dot.Mom(1:4,j))-(Mom(1:4,i).dot.Mom(1:4,j)))     )! C-E z
          Ria = 2d0*(Mom(1:4,j).dot.Mom(1:4,k))/(1d0-z)! = mt^2 * (1-r^2)

! print *,"check", mt**2-Ria ,(Mom(1:4,k)-Mom(1:4,i)-Mom(1:4,j)).dot.(Mom(1:4,k)-Mom(1:4,i)-Mom(1:4,j))
! print *, "check2",(1d0-z)/z,(Mom(1:4,j).dot.Mom(1:4,k))/(Mom(1:4,i).dot.Mom(1:4,k))
! print *, "check3",z,2d0*(Mom(1:4,i).dot.Mom(1:4,k))/Ria

!           Pia(1:4) = Mom(1:4,i)+Mom(1:4,j)-Mom(1:4,k)
!           Pia2 = Pia(1:4).dot.Pia(1:4)
!           MomT(1:4,1:4)=Mom(1:4,1:4)
!           MomT(1:4,i) = dsqrt(lambda(Pia2,mt**2,0d0))/dsqrt(lambda(2d0*(Mom(1:4,i).dot.Mom(1:4,j)),Pia2,mt**2)) * (Mom(1:4,k)-(Pia(1:4).dot.Mom(1:4,k))/Pia2*Pia(1:4))  &
!                       + (Pia2-mt**2)/2d0/Pia2*Pia(1:4)
!           MomT(1:4,k) = MomT(1:4,i) - Pia(1:4)


!        new mapping
         MomK(1:4) = Mom(1:4,2)+Mom(1:4,3)
         call DipoleTrafox(Mom(1:4,1),MomK(1:4),Mom(1:4,2),MomT(1:4,2))
         call DipoleTrafox(Mom(1:4,1),MomK(1:4),MomLep(1:4,1),MomLepT(1:4,1))
         call DipoleTrafox(Mom(1:4,1),MomK(1:4),MomLep(1:4,2),MomLepT(1:4,2))
         call DipoleTrafox(Mom(1:4,1),MomK(1:4),Mom(1:4,3),MomT(1:4,3))
         MomT(1:4,4) = Mom(1:4,1) - MomT(1:4,2) - MomT(1:4,3)
!          print *, Mom(1:4,1).dot.Mom(1:4,1)
!          print *, Mom(1:4,2).dot.Mom(1:4,2)
!          print *, Mom(1:4,3).dot.Mom(1:4,3)
!          print *, Mom(1:4,4).dot.Mom(1:4,4)
!          print *, Mom(1:4,1) - MomT(1:4,2) - MomT(1:4,3) - MomT(1:4,4)



          if( MomT(1,4)/mt.lt.1d-2 .or. (MomT(1:4,3).dot.MomT(1:4,4))/mt**2.lt.1d-2) goto 18
!         endif
          ExtParticles(1)%Mom(1:4) = dcmplx(MomT(1:4,1))
          ExtParticles(2)%Mom(1:4) = dcmplx(MomT(1:4,3))
          ExtParticles(3)%Mom(1:4) = dcmplx(MomT(1:4,4))
          ExtParticles(4)%Mom(1:4) = dcmplx(MomT(1:4,2))

                  ! this is the new dipole which includes the i-f contribution
!                   diag= 2d0*( 0.5d0*(z/(1d0-z) + (1d0-z)/z)  - mt**2/2d0*(Mom(1:4,i).dot.Mom(1:4,j))/(Mom(1:4,k).dot.Mom(1:4,j))**2  )  ! works with Dittmaier z
                  diag= 2d0*( z/(1d0-z)   - mt**2/2d0*(Mom(1:4,i).dot.Mom(1:4,j))/(Mom(1:4,k).dot.Mom(1:4,j))**2  )  ! works with Dittmaier z and C-E z
                  offdiag = 1d0/(Mom(1:4,i).dot.Mom(1:4,j))

!                   paux(1:4) = z*Mom(1:4,i) - (1d0-z)*Mom(1:4,j)
!                   paux(1:4) = ( 2d0*(Mom(1:4,k).dot.Mom(1:4,i))*Mom(1:4,i) - 2d0*(Mom(1:4,k).dot.Mom(1:4,j))*Mom(1:4,j) )/Ria

                  MomA(1:4) = ( (Mom(1:4,k).dot.Mom(1:4,i))*Mom(1:4,i) - (Mom(1:4,k).dot.Mom(1:4,j))*Mom(1:4,j) )*2d0/Ria
!                   MomKTilde(1:4) = (mt**2-(MomK(1:4).dot.MomK(1:4)))/2d0/dsqrt( (MomK(1:4).dot.Mom(1:4,1))**2 - mt**2*(MomK(1:4).dot.MomK(1:4)) )  &
!                                       * ( MomK(1:4)-(MomK(1:4).dot.Mom(1:4,1))/mt**2*Mom(1:4,1) ) + (mt**2+(MomK(1:4).dot.MomK(1:4)))/2d0/mt**2*Mom(1:4,1)
!                   Pia(1:4) = Mom(1:4,1) - MomKTilde(1:4)
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
                  if(WDecay) then
                      call WPolVec2(WPlus,Wp_DKmode,MomLep(1:4,1:2),WDK_LO,ExtParticles(4)%Pol(1:4))
                      ExtParticles(4)%Pol(1:4)=ExtParticles(4)%Pol(1:4)*WProp*g_weak/dsqrt(2d0)
		      if( hel4.ne.1 ) ExtParticles(4)%Pol(1:4) = (0d0,0d0)
                  else
                      call pol_massSR(ExtParticles(4)%Mom(1:4),mw,Hel4,ExtParticles(4)%Pol(1:4))
                      ExtParticles(4)%Pol(1:4) = ExtParticles(4)%Pol(1:4)*g_weak/dsqrt(2d0)
                  endif



                  call EvalTree2(TreeAmpsDip(1),AmplH(hel3))
        enddo
                  SqMatDip = SqMatDip + HH(+1,+1)*AmplH(+1)*dconjg(AmplH(+1)) &
                                      + HH(-1,-1)*AmplH(-1)*dconjg(AmplH(-1)) &
                                      + HH(-1,+1)*AmplH(-1)*dconjg(AmplH(+1)) &
                                      + HH(+1,-1)*AmplH(+1)*dconjg(AmplH(-1))
        enddo
        enddo
        enddo
!                                           (4pi*alpha_s)
        SqMatDip = SqMatDip * 1d0/2d0/3d0  *  (1.633628179866692484d0) * 8d0


            DipColorCorr = -3d0/2d0
            dipole6 = 3.26725635973d0 * SqMatDip * 1d0/(Mom(1:4,i).dot.Mom(1:4,j))/4d0 *DipColorCorr !  CS f-i dipole for gg



18 continue



!        if( SingType.eq.fi_ ) then
          i=5; j=4; k=1;
!           x = 1d0 - (Mom(1:4,i).dot.Mom(1:4,j))/((Mom(1:4,k).dot.Mom(1:4,i))+(Mom(1:4,k).dot.Mom(1:4,j)))
!           z = (Mom(1:4,i).dot.Mom(1:4,k))/((Mom(1:4,i).dot.Mom(1:4,k))+(Mom(1:4,k).dot.Mom(1:4,j)))!  Dittmaier z
          z = 1d0-(    (Mom(1:4,k).dot.Mom(1:4,j))/((Mom(1:4,i).dot.Mom(1:4,k))+(Mom(1:4,k).dot.Mom(1:4,j))-(Mom(1:4,i).dot.Mom(1:4,j)))     )! C-E z
          Ria = 2d0*(Mom(1:4,j).dot.Mom(1:4,k))/(1d0-z)! = mt^2 * (1-r^2)

!           Pia(1:4) = Mom(1:4,i)+Mom(1:4,j)-Mom(1:4,k)
!           Pia2 = Pia(1:4).dot.Pia(1:4)
!           MomT(1:4,1:4)=Mom(1:4,1:4)
!           MomT(1:4,4) = dsqrt(lambda(Pia2,mt**2,0d0))/dsqrt(lambda(2d0*(Mom(1:4,i).dot.Mom(1:4,j)),Pia2,mt**2)) * (Mom(1:4,k)-(Pia(1:4).dot.Mom(1:4,k))/Pia2*Pia(1:4))  &
!                       + (Pia2-mt**2)/2d0/Pia2*Pia(1:4)
!           MomT(1:4,k) = MomT(1:4,4) - Pia(1:4)


!        new mapping
         MomK(1:4) = Mom(1:4,2)+Mom(1:4,3)
         call DipoleTrafox(Mom(1:4,1),MomK(1:4),Mom(1:4,2),MomT(1:4,2))
         call DipoleTrafox(Mom(1:4,1),MomK(1:4),MomLep(1:4,1),MomLepT(1:4,1))
         call DipoleTrafox(Mom(1:4,1),MomK(1:4),MomLep(1:4,2),MomLepT(1:4,2))
         call DipoleTrafox(Mom(1:4,1),MomK(1:4),Mom(1:4,3),MomT(1:4,3))
         MomT(1:4,4) = Mom(1:4,1) - MomT(1:4,2) - MomT(1:4,3)
!          print *, Mom(1:4,1).dot.Mom(1:4,1)
!          print *, Mom(1:4,2).dot.Mom(1:4,2)
!          print *, Mom(1:4,3).dot.Mom(1:4,3)
!          print *, Mom(1:4,4).dot.Mom(1:4,4)
!          print *, Mom(1:4,1) - MomT(1:4,2) - MomT(1:4,3) - MomT(1:4,4)



          if( MomT(1,4)/mt.lt.1d-2 .or. (MomT(1:4,3).dot.MomT(1:4,4))/mt**2.lt.1d-2) goto 19
!         endif
          ExtParticles(1)%Mom(1:4) = dcmplx(MomT(1:4,1))
          ExtParticles(2)%Mom(1:4) = dcmplx(MomT(1:4,3))
          ExtParticles(3)%Mom(1:4) = dcmplx(MomT(1:4,4))
          ExtParticles(4)%Mom(1:4) = dcmplx(MomT(1:4,2))

                  ! this is the new dipole which includes the i-f contribution
!                   diag= 2d0*( 0.5d0*(z/(1d0-z) + (1d0-z)/z)  - mt**2/2d0*(Mom(1:4,i).dot.Mom(1:4,j))/(Mom(1:4,k).dot.Mom(1:4,j))**2  )  ! works with Dittmaier z
                  diag= 2d0*( z/(1d0-z)  - mt**2/2d0*(Mom(1:4,i).dot.Mom(1:4,j))/(Mom(1:4,k).dot.Mom(1:4,j))**2  )  ! works with Dittmaier z and C-E z
                  offdiag = 1d0/(Mom(1:4,i).dot.Mom(1:4,j))

!                   paux(1:4) = z*Mom(1:4,i) - (1d0-z)*Mom(1:4,j)
!                   paux(1:4) = ( 2d0*(Mom(1:4,k).dot.Mom(1:4,i))*Mom(1:4,i) - 2d0*(Mom(1:4,k).dot.Mom(1:4,j))*Mom(1:4,j) )/Ria

                  MomA(1:4) = ( (Mom(1:4,k).dot.Mom(1:4,i))*Mom(1:4,i) - (Mom(1:4,k).dot.Mom(1:4,j))*Mom(1:4,j) )*2d0/Ria
!                   MomK(1:4) = Mom(1:4,2)+Mom(1:4,3)! P_W+P_B
!                   MomKTilde(1:4) = (mt**2-(MomK(1:4).dot.MomK(1:4)))/2d0/dsqrt( (MomK(1:4).dot.Mom(1:4,1))**2 - mt**2*(MomK(1:4).dot.MomK(1:4)) )  &
!                                       * ( MomK(1:4)-(MomK(1:4).dot.Mom(1:4,1))/mt**2*Mom(1:4,1) ) + (mt**2+(MomK(1:4).dot.MomK(1:4)))/2d0/mt**2*Mom(1:4,1)
!                   Pia(1:4) = Mom(1:4,1) - MomKTilde(1:4)
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

!print *, "z",z
!print *, HH(+1,+1)/(Mom(1:4,i).dot.Mom(1:4,j))/2d0
!print *, HH(-1,-1)/(Mom(1:4,i).dot.Mom(1:4,j))/2d0
!print *, HH(-1,+1)/(Mom(1:4,i).dot.Mom(1:4,j))/2d0
!print *, HH(+1,-1)/(Mom(1:4,i).dot.Mom(1:4,j))/2d0
!pause


          SqMatDip = (0d0,0d0)
          do hel1=1,-1,-2
          do hel2=1,-1,-2
          do hel4=1,-1,-1
          do hel3=1,-1,-2
                  call uSpi(ExtParticles(1)%Mom(1:4),ExtParticles(1)%Mass,Hel1,ExtParticles(1)%Pol(1:4))
                  call ubarSpi(ExtParticles(2)%Mom(1:4),ExtParticles(2)%Mass,Hel2,ExtParticles(2)%Pol(1:4))
                  call pol_mless(ExtParticles(3)%Mom(1:4),Hel3,ExtParticles(3)%Pol(1:4))
                  if(WDecay) then
                      call WPolVec2(WPlus,Wp_DKmode,MomLep(1:4,1:2),WDK_LO,ExtParticles(4)%Pol(1:4))
                      ExtParticles(4)%Pol(1:4)=ExtParticles(4)%Pol(1:4)*WProp*g_weak/dsqrt(2d0)
		      if( hel4.ne.1 ) ExtParticles(4)%Pol(1:4) = (0d0,0d0)
                  else
                      call pol_massSR(ExtParticles(4)%Mom(1:4),mw,Hel4,ExtParticles(4)%Pol(1:4))
                      ExtParticles(4)%Pol(1:4) = ExtParticles(4)%Pol(1:4)*g_weak/dsqrt(2d0)
                  endif


                  call EvalTree2(TreeAmpsDip(1),AmplH(hel3))
          enddo
                  SqMatDip = SqMatDip + HH(+1,+1)*AmplH(+1)*dconjg(AmplH(+1)) &
                                      + HH(-1,-1)*AmplH(-1)*dconjg(AmplH(-1)) &
                                      + HH(-1,+1)*AmplH(-1)*dconjg(AmplH(+1)) &
                                      + HH(+1,-1)*AmplH(+1)*dconjg(AmplH(-1))
          enddo
          enddo
          enddo
!                                               (4pi*alpha_s)
        SqMatDip = SqMatDip * 1d0/2d0/3d0  *  (1.633628179866692484d0) * 8d0


            DipColorCorr = -3d0/2d0
            dipole7 = 3.26725635973d0 * SqMatDip * 1d0/(Mom(1:4,i).dot.Mom(1:4,j))/4d0 *DipColorCorr !  CS f-i dipole for gg




19 continue


        Dipole = Dipole1 + Dipole2 + Dipole5 + Dipole3 + Dipole6 + Dipole7 + Dipole4


! i=1; j=5; k=3;
! print *, "CHECK", Dipole5 /(3.26725635973d0 * SqMatDip * ( (Mom(1:4,i).dot.Mom(1:4,k))/((Mom(1:4,i)+Mom(1:4,k)).dot.Mom(1:4,j)) - mt**2/2d0/(Mom(1:4,i).dot.Mom(1:4,j)) )/(Mom(1:4,i).dot.Mom(1:4,j))  *2d0/3d0  )

! print *, "diff",dble(SqMat)-Dipole


        print *, "SqMat",SqMat
        print *, "Dipoles",Dipole,dble(Dipole/SqMat)+1d0
!         print *, "-----------------------------------"

!
!         print *, "Dipole1",Dipole1,Dipole1/dble(SqMat)
!         print *, "Dipole2",Dipole2,Dipole2/dble(SqMat)
!         print *, "Dipole3",Dipole3,Dipole3/dble(SqMat)
!         print *, "Dipole4",Dipole4,Dipole4/dble(SqMat)
!         print *, "Dipole5",Dipole5,Dipole5/dble(SqMat)
!         print *, "Dipole6",Dipole6,Dipole6/dble(SqMat)
!         print *, "Dipole7",Dipole7,Dipole7/dble(SqMat)




      enddo







end program


! Lorentz transformation for the dipole matrix element in top decays with additional photon in w decay
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



SUBROUTINE WTransform(MomDK,MomDKTd,pbDpg,ptDpg,ptDpb)
use ModMisc
implicit none
real(8) :: MomDK(1:4,1:4),MomDKTd(1:4,1:3),pw(4),pt(4),lDt(3:4),lDw(3:4)
real(8) :: root,hsin,hcos,a,b
real(8) :: ptDpt,pwDpw,ptDpw,ptDpg,pbDpg,ptDpb,pWDpl,ptDpl,pWDpn,ptDpn


    pw(1:4) = MomDK(1:4,2) !+ MomDK(1:4,3)
    pt(1:4) = pw(1:4) + MomDK(1:4,1) + MomDK(1:4,4)

    pbDpg = dot(MomDK(1:4,1),MomDK(1:4,4))
    ptDpg = dot(pt(1:4),MomDK(1:4,4))
    ptDpb = dot(pt(1:4),MomDK(1:4,1))
    ptDpw = dot(pt(1:4),pw(1:4))
    ptDpt = dot(pt(1:4),pt(1:4))
    pwDpw = dot(pw(1:4),pw(1:4))
    pWDpl = dot(pw(1:4),MomDK(1:4,2))
    ptDpl = dot(pt(1:4),MomDK(1:4,2))
    pWDpn = dot(pw(1:4),MomDK(1:4,3))
    ptDpn = dot(pt(1:4),MomDK(1:4,3))

    root=dsqrt(ptDpw**2-ptDpt*pwDpw)
    hsin=0.5d0/(ptDpt*pwDpw)*(-(ptDpt-pwDpw)*ptDpw+(ptDpt+pwDpw)*root)
    hcos=0.5d0/(ptDpt*pwDpw)*(+(ptDpt+pwDpw)*ptDpw-(ptDpt-pwDpw)*root)
    a=hsin/root
    b=(hcos-1d0)/root**2

    MomDKTd(1:4,2) = MomDK(1:4,2) + a*( pt(1:4)*pwDpl-pw(1:4)*ptDpl )  &
                                  + b*( ptDpw*(pt(1:4)*pwDpl+pw(1:4)*ptDpl) - pt(1:4)*pwDpw*ptDpl - pw(1:4)*ptDpt*pwDpl )
    MomDKTd(1:4,3) = MomDK(1:4,3) + a*( pt(1:4)*pwDpn-pw(1:4)*ptDpn )  &
                                  + b*( ptDpw*(pt(1:4)*pwDpn+pw(1:4)*ptDpn) - pt(1:4)*pwDpw*ptDpn - pw(1:4)*ptDpt*pwDpn )
    MomDKTd(1:4,1) = pt(1:4) - MomDKTd(1:4,2) - MomDKTd(1:4,3)

RETURN
END SUBROUTINE


  subroutine Wpolvec2(Wpm,WDK,MomPol,order,Wpolvec)
    use ModParameters
    use ModMisc
    implicit none
    complex(8) :: Wpolvec(1:4)
    complex(8) ::  Wcurr1(1:4), Wcurr2(1:4), Wmom(1:4), wcurr(1:4), NeuSpi(1:4), LepSpi(1:4), wcurr3(1:4)
    integer :: WDK, Wpm,QW, order
    real(8) :: MomPol(1:4,1:2)
    real(8),parameter :: sq2 =  1.414213562373095D0
    real(8),parameter :: g2_weak = 4d0*dsqrt(2d0)*m_W**2*GF
    real(8),parameter :: g_weak = dsqrt(g2_weak)
    real(8) :: sw, EL, coup, Q_f1, Q_f2, check1, check2, check3, check4, check5
    complex(8) :: klep_dn(1:4), kneu_up(1:4), kpho(1:4), PropMom(1:4), EpsP(1:4)
    complex(8) :: SpiLep_dn(1:4), SpiNeu_up(1:4), SpiLep2(1:4),SpiLep1(1:4), coupl_sqrt
    integer :: Dv
    logical nan1, nan2, nan3


    QW=Wpm
! Resolving markus coupling structure, sin(Theta_w)
    sw = dsqrt(4.d0*DblPi*alpha/g2_weak)
    EL = dsqrt(4.d0*DblPi*alpha)
    Dv =4
    If (order .ne. 0 ) then
       Write(*,*) "Only, LO is implemented <--> order =0"
    endif
    If(QW .ne. -1 .and. QW .ne. 1) then
       write(*,*) "Sorry, we can only deliver W+/W- polarisation vectors"
       stop
    endif


    klep_dn(1:4) = dcmplx(MomPol(1:4,1),0.d0)
    kneu_up(1:4) = dcmplx(MomPol(1:4,2),0.d0)


    if(Q_top .eq. -4.d0/3.d0) QW=-QW

! W+/W- vertex is I *EL /sqrt(2)/sw  times ne * sq2
! absorbing Factor -ne/sq2 introduced by vbqq_Weyl
    coup = -ne*EL/sq2/sw  * (ne *sq2)
    ! No W decay



    ! ALL W+ cases
    If(QW .eq. 1) then
       call    vSpi_Weyl(klep_dn,+1,SpiLep_dn)     ! l+ or dn_bar
       call ubarSpi_Weyl(kneu_up,-1,SpiNeu_up)  ! nu or up
       ! Anti-Lepton + Neutrino
       if(WDK .eq. 1) then
          Wpolvec(1:4) = vbqq_Weyl(Dv,SpiNeu_up,SpiLep_dn)* coup
          return
       endif
       ! Anti-Down + Up
       if(WDK .eq. 3) then
          Wpolvec(1:4) = vbqq_Weyl(Dv,SpiNeu_up,SpiLep_dn)* coup * dsqrt(3.d0) * sq2 ! color-Factor dsqrt(3) *  Flavour  dsqrt(2)
          return
       endif

       ! Anti-Lepton + Neutrino + Photon
       ! Anti-Down + Up + Photon

    endif

    ! All W- cases
    If(QW .eq. -1) then
        call ubarSpi_Weyl(klep_dn,-1,SpiLep_dn)  ! l- or dn
        call    vSpi_Weyl(kneu_up,+1,SpiNeu_up)     ! nubar or up_bar

       !Lepton +  Anti-Neutrino
       if(WDK .eq. 1) then
          Wpolvec(1:4) = vbqq_Weyl(Dv,SpiLep_dn,SpiNeu_up)* coup
          return
       endif
       !Down +  Anti-Up
       if(WDK .eq. 3) then
          Wpolvec(1:4) =  vbqq_Weyl(Dv,SpiLep_dn,SpiNeu_up)* coup * dsqrt(3.d0) * sq2 ! color-Factor dsqrt(3) *  Flavour  dsqrt(2)
          return
       endif

       ! Lepton + Anti-Neutrino + Photon
       ! Down + Anti-Up + Photon



    endif

  end















