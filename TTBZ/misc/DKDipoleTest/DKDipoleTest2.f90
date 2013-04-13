program DKDipoles
implicit none
real(8) :: Mom(1:4,1:4),MomT(1:4,1:3),pbDpg,ptDpg,ptDpb,omz,rsq,z,y
real(8) :: mt=172d0,SingDepth,MadGraph_tree,MadGraph_tree2,dot,Dipole2,Dipole
integer :: Pcol1,Pcol2,Steps,i

     Pcol1= 4 -1
     Pcol2= 6 -1
     SingDepth = 1d-8
     Steps = 10
     Mom(1:4,1)= (/mt,0d0,0d0,0d0/)
     call coupsm(0)

      do i=1,20
          print *,i
          call gensing(4,mt,(/81.419d0,0d0,0d0,0d0/),Mom(1:4,2:5),Pcol1,Pcol2,SingDepth,Steps)
          call ST_BWPGG((/Mom(1:4,1),Mom(1:4,3),Mom(1:4,2),Mom(1:4,4),Mom(1:4,5)/),MadGraph_tree)

          if( Mom(1,4).lt.20d0 ) cycle
          if( dot(Mom(1:4,4),Mom(1:4,3))/mt**2.lt.1d-2 ) cycle
          if( dot(Mom(1:4,4),Mom(1:4,5))/mt**2.lt.1d-2 ) cycle

          print *, "Sing:", dot(Mom(1:4,3),Mom(1:4,5))/mt**2
          print *, "MadGraph:", MadGraph_tree


!           call WTransform((/Mom(1:4,3),Mom(1:4,2),Mom(1:4,5),Mom(1:4,4)/),MomT(1:4,1:3),pbDpg,ptDpg,ptDpb)
!           print *, dot(MomT(1:4,1),MomT(1:4,1))





          MomT(1:4,1)=Mom(1:4,3)+Mom(1:4,5)
!           MomT(1:4,2)=Mom(1:4,2)
          call ST_BWPG((/Mom(1:4,1),MomT(1:4,1),MomT(1:4,2),MomT(1:4,4)/),MadGraph_tree2)
          print *, "MadGraph2:", MadGraph_tree2


          pause
      enddo


end program




SUBROUTINE WTransform(MomDK,MomDKTd,pbDpg,ptDpg,ptDpb)
implicit none
real(8) :: MomDK(1:4,1:4),MomDKTd(1:4,1:3),pw(4),pt(4),lDt(3:4),lDw(3:4)
real(8) :: root,hsin,hcos,a,b
real(8) :: ptDpt,pwDpw,ptDpw,ptDpg,pbDpg,ptDpb,pWDpl,ptDpl,pWDpn,ptDpn,dot


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




FUNCTION dot(p1,p2)
implicit none
real(8), intent(in) :: p1(1:4),p2(1:4)
real(8)             :: dot

     dot = p1(1)*p2(1)  &
                      - p1(2)*p2(2)  &
                      - p1(3)*p2(3)  &
                       -p1(4)*p2(4)
END FUNCTION dot