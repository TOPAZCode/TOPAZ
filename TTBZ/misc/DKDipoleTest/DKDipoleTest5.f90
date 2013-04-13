program DKDipoles!     mat.el. for t->b w g    and dipole    with new currents
use ModAmplitudes
use ModProcess
use modMisc
implicit none
real(8) :: Mom(1:4,1:5),MomT(1:4,1:4),pbDpg,ptDpg,ptDpb,omz,rsq,z,y,Q(1:4),Qsq
real(8) :: mt=172d0,SingDepth,MadGraph_tree,MadGraph_tree2,Dipole2,Dipole,mw=80.419d0
integer :: Pcol1,Pcol2,Steps,hel1,hel2,hel3,hel4,hel5,n,i,j,k,ij
complex(8) :: Ampl1,Ampl2,SqMat,SqMatDip
integer :: NumQuarks,NumGluons,NumTrees
type(TreeProcess) :: TreeAmpsDip(1:2)
integer :: NumParticles,iTree
type(Particle) :: ExtParticles(1:5)

 call coupsm(0)
print *, "HELAS INIT DISABLED"





  NumQuarks=2
  NumGluons=1
  NumParticles=NumQuarks+NumGluons+1
  TreeAmpsDip(1)%NumPart=NumParticles
  TreeAmpsDip(1)%NumQua=NumQuarks
  allocate( TreeAmpsDip(1)%NumGlu(0:NumQuarks)     )
  allocate( TreeAmpsDip(1)%PartRef(1:NumParticles) )
  allocate( TreeAmpsDip(1)%PartType(1:NumParticles))
  allocate( TreeAmpsDip(1)%Quarks(1:NumQuarks)     )
  allocate( TreeAmpsDip(1)%Gluons(1:NumGluons)     )
  TreeAmpsDip(1)%NumGlu(0) = NumGluons




  do n=1,12
          print *, n

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
          ExtParticles(4)%Mass2= ExtParticles(3)%Mass**2

          TreeAmpsDip(1)%PartRef(1:4) = (/1,2,3,4/)
          call LinkTreeParticles(TreeAmpsDip(1),ExtParticles(1:4))


          Pcol1= 4 -1
          Pcol2= 5 -1
          SingDepth = 1d-8
          Steps = 10

          Mom(1:4,1)= (/mt,0d0,0d0,0d0/)
          call gensing(3,mt,(/mw,0d0,0d0/),Mom(1:4,2:4),Pcol1,Pcol2,SingDepth,Steps)
!           call genps(3,mt,(/0.9d0,0.4d0,0.4d0,0.2d0,0.2d0/),(/mw,0d0,0d0/),Mom(1:4,2:4),SingDepth)

!           if( Mom(1,4)/mt .lt. 1d-2 ) stop
!           print *, "Sing ",(Mom(1:4,3).dot.Mom(1:4,4))/mt**2




          ExtParticles(1)%Mom(1:4) = dcmplx(Mom(1:4,1))
          ExtParticles(2)%Mom(1:4) = dcmplx(Mom(1:4,3))
          ExtParticles(3)%Mom(1:4) = dcmplx(Mom(1:4,4))
          ExtParticles(4)%Mom(1:4) = dcmplx(Mom(1:4,2))

          SqMat = (0d0,0d0)
          do hel1=1,-1,-2
          do hel2=1,-1,-2
          do hel3=1,-1,-2
          do hel4=1,-1,-1
                  call uSpi(ExtParticles(1)%Mom(1:4),ExtParticles(1)%Mass,Hel1,ExtParticles(1)%Pol(1:4))
                  call ubarSpi(ExtParticles(2)%Mom(1:4),ExtParticles(2)%Mass,Hel2,ExtParticles(2)%Pol(1:4))
                  call pol_mless(ExtParticles(3)%Mom(1:4),Hel3,ExtParticles(3)%Pol(1:4))
                  call pol_massSR(ExtParticles(4)%Mom(1:4),mw,Hel4,ExtParticles(4)%Pol(1:4))

                  call EvalTree2(TreeAmpsDip(1),Ampl1)

                  SqMat = SqMat + Ampl1*dconjg(Ampl1)

        enddo
        enddo
        enddo
        enddo
!                                        4pi*alpha/2/sw2              (4pi*alpha_s)        color.fact
        SqMat = SqMat * 1d0/2d0/3d0* 0.2063277518095031333d0  *  (1.633628179866692484d0)   * 8d0
        print *, "SqMat",SqMat

!         call ST_BWPG((/Mom(1:4,1),Mom(1:4,3),Mom(1:4,2),Mom(1:4,4)/),MadGraph_tree)
!         print *, "MadGraph2:", MadGraph_tree
!         print *, "ratio:", MadGraph_tree/SqMat


! --------------------------------- dipole calculation -------------------------------------------------




!           this is the f-f mapping of CS for a massive spectator
            i=3; j=4; k=2; ij=3;
            y = (Mom(1:4,i).dot.Mom(1:4,j))/((Mom(1:4,i).dot.Mom(1:4,j))+(Mom(1:4,j).dot.Mom(1:4,k))+(Mom(1:4,i).dot.Mom(1:4,k)))
            z = (Mom(1:4,i).dot.Mom(1:4,k))/((Mom(1:4,i).dot.Mom(1:4,k))+(Mom(1:4,j).dot.Mom(1:4,k)))
            Q(1:4) = Mom(1:4,i)+Mom(1:4,j)+Mom(1:4,k)
            Qsq = Q(1:4).dot.Q(1:4)
            MomT(1:4,1:4)=Mom(1:4,1:4)
            MomT(1:4,k) = SqrtLambda(Qsq,0d0,mw**2)/SqrtLambda(Qsq,2d0*(Mom(1:4,i).dot.Mom(1:4,j)),mw**2) * ( Mom(1:4,k) - (Q(1:4).dot.Mom(1:4,k))/Qsq*Q(1:4) ) + (Qsq+mw**2)/2d0/Qsq*Q(1:4)
            MomT(1:4,ij) = Q(1:4) - MomT(1:4,k)
! print *, "check", Mom(1:4,k).dot.Mom(1:4,k),MomT(1:4,k).dot.MomT(1:4,k)
! print *, "check", Mom(1:4,i).dot.Mom(1:4,i),MomT(1:4,ij).dot.MomT(1:4,ij)
! print *, "check", Mom(1:4,1)-Mom(1:4,2)-Mom(1:4,3)-Mom(1:4,4)
! print *, "check", MomT(1:4,1)-MomT(1:4,2)-MomT(1:4,3)


!           this is the i-f mapping for a massive emitter
            i=3; j=4; k=1; ij=3;
            z = dsqrt(mt**2/(mt**2-2d0*(Mom(1:4,k).dot.Mom(1:4,j)) ) )
            Q(1:4) = Mom(1:4,k)-Mom(1:4,j)
            y = 1d0+(z-1d0)*(Q(1:4).dot.Q(1:4))/(Mom(1:4,i).dot.Q(1:4))
            MomT(1:4,1:4)=Mom(1:4,1:4)
            MomT(1:4,k) = z * (Mom(1:4,k)-Mom(1:4,j))




            Q(1)   = Mom(1,2)
            Q(2:4) =-Mom(2:4,2)
            Qsq = Q(1:4).dot.Q(1:4)
            print *,qsq
pause
print *, MomT(1:4,k)-MomT(1:4,j)
            call boost(MomT(1:4,k),Q(1:4),Qsq)
            call boost(MomT(1:4,j),Q(1:4),Qsq)
print *, MomT(1:4,k)-MomT(1:4,j)
pause
            call boost(MomT(1:4,k),Q(1:4),mt)
            call boost(MomT(1:4,i),Q(1:4),mt)
            call boost(MomT(1:4,2),Q(1:4),mt)

!             MomT(1:4,ij)=  MomT(1:4,i) * y

! print *, "check", Mom(1:4,k).dot.Mom(1:4,k),MomT(1:4,k).dot.MomT(1:4,k)
! print *, "check", Mom(1:4,i).dot.Mom(1:4,i),MomT(1:4,ij).dot.MomT(1:4,ij)
! print *, "check", Mom(1:4,2).dot.Mom(1:4,2),MomT(1:4,2).dot.MomT(1:4,2)
! print *, "check", Mom(1:4,1)-Mom(1:4,2)-Mom(1:4,3)-Mom(1:4,4)
print *, "check", MomT(1:4,1)-MomT(1:4,2)-MomT(1:4,3)

          NumQuarks=2
          NumGluons=0
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

          ExtParticles(3)%PartType = 13! w+
          ExtParticles(3)%ExtRef   = 3
          ExtParticles(3)%Mass = mw
          ExtParticles(3)%Mass2= ExtParticles(3)%Mass**2

          TreeAmpsDip(1)%PartRef(1:3) = (/1,2,3/)
          call LinkTreeParticles(TreeAmpsDip(1),ExtParticles(1:3))



          ExtParticles(1)%Mom(1:4) = dcmplx(MomT(1:4,1))
          ExtParticles(2)%Mom(1:4) = dcmplx(MomT(1:4,3))
          ExtParticles(3)%Mom(1:4) = dcmplx(MomT(1:4,2))

          SqMatDip = (0d0,0d0)
          do hel1=1,-1,-2
          do hel2=1,-1,-2
          do hel3=1,-1,-1
                  call uSpi(ExtParticles(1)%Mom(1:4),ExtParticles(1)%Mass,Hel1,ExtParticles(1)%Pol(1:4))
                  call ubarSpi(ExtParticles(2)%Mom(1:4),ExtParticles(2)%Mass,Hel2,ExtParticles(2)%Pol(1:4))
                  call pol_massSR(ExtParticles(3)%Mom(1:4),mw,Hel3,ExtParticles(3)%Pol(1:4))

                  call EvalTree2(TreeAmpsDip(1),Ampl2)

                  SqMatDip = SqMatDip + Ampl2*dconjg(Ampl2)

          enddo
          enddo
          enddo
  !                                        4pi*alpha/2/sw2
          SqMatDip = SqMatDip * 1d0/2d0* 0.2063277518095031333d0
          print *, "SqMatDip",SqMatDip

!           call ST_BWP((/MomT(1:4,1),MomT(1:4,3),MomT(1:4,2)/),MadGraph_tree)
!           print *, "MadGraph2:", MadGraph_tree
!           print *, "ratio:", MadGraph_tree/SqMatDip


!                     8pi*alpha_s                   splitting fct
!           dipole = 3.26725635973d0 * SqMatDip * (1d0+z**2)/(1d0-z) /(Mom(1:4,i).dot.Mom(1:4,j))
          dipole = 3.26725635973d0 * SqMatDip *( 2d0/(1d0-z*(1d0-y))-1d0-z ) /(Mom(1:4,i).dot.Mom(1:4,j))   ! CS dipole
          print *, "Dipole", dipole / SqMat / 3d0  *2d0


!           Dipole = -3.26725635973d0 * 4d0/3d0 * ( 1d0/ptDpg*(ptDpb/(ptDpg+pbDpg) - mt**2/ptDpg/2d0) + 1d0/pbDpg*ptDpb/(ptDpg+pbDpg) )   ! eikonal approx.
!           print *, "Dipole", dipole / SqMat


    enddo




end program








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







SUBROUTINE LinkTreeParticles(TheTreeAmp,TheParticles) ! NEW
use ModProcess
implicit none
type(TreeProcess) :: TheTreeAmp
type(Particle),target :: TheParticles(1:5)
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







          subroutine uSpi(p,m,i,f)
          implicit none
          integer i
          real(8) m
          complex(8) p(4)
          complex(8) f(4),fc
          real(8)  p0,px,py,pz,fc2

          p0=dreal(p(1))
          px=dreal(p(2))
          py=dreal(p(3))
          pz=dreal(p(4))

          fc2=p0+m
          fc=cdsqrt( dcmplx(fc2))
!           fc=dsqrt(fc2)
          if (i.eq.1) then
            f(1)=fc
            f(2)=dcmplx(0d0,0d0)
            f(3)=pz*fc/fc2
            f(4)=(px+(0d0,1d0)*py)*fc/fc2
          elseif (i.eq.-1) then
            f(1)=dcmplx(0d0,0d0)
            f(2)=fc
            f(3)=(px-(0d0,1d0)*py)*fc/fc2
            f(4)=-pz*fc/fc2
          else
stop
          endif

          return
          end subroutine











