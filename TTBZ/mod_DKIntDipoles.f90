MODULE ModDKIntDipoles
use ModParameters
use ModMisc
implicit none

public :: DKff_qq, DKff_qg, DKff_gg, DKfi_qq, DKfi_gg

real(8),parameter,private  :: pisqo6=pi**2/6d0
character,parameter :: scheme*(4)='tH-V'  ! this parameter should match the scheme of the virtual corrections, options: 'fdh', 'tH-V'

contains

! corrected bugs: DblPi -> DblPi**2


!  the global prefactor for all dipoles below has to be
!  alpha_s/(2*Pi) 1/Gamma(1-eps) * (4*Pi*mu^2/m2)^eps
!
! sijk = 2*p_ij.p_k in the all-outgoing convention

! References:
! NT:  hep-ph/9806317
! CET: hep-ph/0408158



FUNCTION  DKff_qq(xe,mu2,sijk)!  NT eq.(22)
implicit none
real(8) :: DKff_qq,DKff_qq_xe(-2:0),mu2,sijk,L,L2
integer :: xe


    L  = dlog(mu2/sijk)
    L2 = L**2

    DKff_qq_xe(-2) = 1d0
    DKff_qq_xe(-1) = 3d0/2d0
    DKff_qq_xe(0)  = 7d0/2d0 - 0.5d0*DblPi**2 -dlog(alpha_DKTff)**2 + 3d0/2d0*( alpha_DKTff - dlog(alpha_DKTff) )

!   expand (mu^2/sijk)^eps
    DKff_qq_xe(-1) = DKff_qq_xe(-1) + L*DKff_qq_xe(-2)
    DKff_qq_xe(0)  = DKff_qq_xe(0)  + L*DKff_qq_xe(-1) + 0.5d0*L2*DKff_qq_xe(-2)


    DKff_qq = DKff_qq_xe(xe)
RETURN
END FUNCTION




FUNCTION  DKff_qg(xe,mu2,sijk)!  NT eq.(23)
implicit none
real(8) :: DKff_qg,DKff_qg_xe(-2:0),mu2,sijk,L,L2
integer :: xe

    L  = dlog(mu2/sijk)
    L2 = L**2

    DKff_qg_xe(-2) = 0d0
    DKff_qg_xe(-1) = -2d0/3d0
    DKff_qg_xe(0)  = -10d0/9d0 - 2d0/3d0*( alpha_DKTff - dlog(alpha_DKTff) )

!   expand (mu^2/sijk)^eps
    DKff_qg_xe(-1) = DKff_qg_xe(-1) + L*DKff_qg_xe(-2)
    DKff_qg_xe(0)  = DKff_qg_xe(0)  + L*DKff_qg_xe(-1) + 0.5d0*L2*DKff_qg_xe(-2)

    DKff_qg = DKff_qg_xe(xe)
RETURN
END FUNCTION





FUNCTION  DKff_gg(xe,mu2,sijk)!  NT eq.(24)
implicit none
real(8) :: DKff_gg,DKff_gg_xe(-2:0),mu2,sijk,L,L2
integer :: xe

    L  = dlog(mu2/sijk)
    L2 = L**2

    DKff_gg_xe(-2) = 2d0
    DKff_gg_xe(-1) = 11d0/3d0
    DKff_gg_xe(0)  = 67d0/9d0 - DblPi**2 -2d0*dlog(alpha_DKTff)**2 + 11d0/3d0*( alpha_DKTff - dlog(alpha_DKTff) )

!   expand (mu^2/sijk)^eps
    DKff_gg_xe(-1) = DKff_gg_xe(-1) + L*DKff_gg_xe(-2)
    DKff_gg_xe(0)  = DKff_gg_xe(0)  + L*DKff_gg_xe(-1) + 0.5d0*L2*DKff_gg_xe(-2)

    DKff_gg = DKff_gg_xe(xe)
RETURN
END FUNCTION








!wrong!FUNCTION DKfi_qq(xe,m2,sijk)!  CET eq.(38)
!wrong!implicit none
!wrong!real(8) :: DKfi_qq,DKfi_qq_xe(-2:0),m2,sijk,L,L2,omrsq
!wrong!integer :: xe
!wrong!
!wrong!     L  = dlog(m2/sijk)!  note that 1-r2=sijk/m2
!wrong!     L2 = L**2
!wrong!    omrsq = -sijk/m2
!wrong!
!wrong!    DKfi_qq_xe(-2) = 1d0
!wrong!    DKfi_qq_xe(-1) = 5d0/2d0 - 2d0*dlog(omrsq)
!wrong!    DKfi_qq_xe(0)  = 25d0/4d0 + 0.5d0*(1d0/omrsq**2-8d0/omrsq+7d0)*dlog(1d0-omrsq) + 0.5d0/omrsq + 2d0*DLi2(omrsq) - 5d0*pisqo6  &
!wrong!                - 5d0*dlog(omrsq) + 2d0*dlog(omrsq)**2 &!    eta=0 (FHD scheme)
!wrong!                - ( 2d0*dlog(alpha_DKTfi)**2 -dlog(alpha_DKTfi)*(-7d0/2d0+4d0*alpha_DKTfi-alpha_DKTfi**2/2d0) &
!wrong!                - 2d0*(1d0-alpha_DKTfi)*(1d0-omrsq)/(omrsq)*dlog((1d0-omrsq)/(1d0-alpha_DKTfi*omrsq)) )
!wrong!
!wrong!
!wrong!    DKfi_qq  = DKfi_qq_xe(xe)
!wrong!RETURN
!wrong!END FUNCTION
!wrong!

! THIS ROUTINE WAS CHANGED BY SCHARF 08/16/2011 BECAUSE HE THINKS mu-dependence is wrong and eta =1, 't Hooft Veltman-scheme. ORIGINAL ABOVE  !!!!
FUNCTION DKfi_qq(xe,m2,mu2,sijk)!  CET eq.(38)
implicit none
real(8) :: DKfi_qq,DKfi_qq_xe(-2:0),m2,sijk,L,L2,omrsq, mu2
integer :: xe

    L  = dlog(mu2/m2)!  note that 1-r2=sijk/m2
    L2 = L**2
    omrsq = -sijk/m2

    DKfi_qq_xe(-2) = 1d0
    DKfi_qq_xe(-1) = 5d0/2d0 - 2d0*dlog(omrsq)
    DKfi_qq_xe(0)  = 25d0/4d0 + 0.5d0*(1d0/omrsq**2-8d0/omrsq+7d0)*dlog(1d0-omrsq) + 0.5d0/omrsq + 2d0*DLi2(omrsq) - 5d0*pisqo6  &
                - 5d0*dlog(omrsq) + 2d0*dlog(omrsq)**2 + 1.d0/2.d0 &!    eta=1 -> eta/2 = 1/2 (t'Hooft Veltman scheme, see ttb-photon paper)
                - ( 2d0*dlog(alpha_DKTfi)**2 -dlog(alpha_DKTfi)*(-7d0/2d0+4d0*alpha_DKTfi-alpha_DKTfi**2/2d0) &
                - 2d0*(1d0-alpha_DKTfi)*(1d0-omrsq)/(omrsq)*dlog((1d0-omrsq)/(1d0-alpha_DKTfi*omrsq)) )

!   factor out (mu^2/m2)^eps
    DKfi_qq_xe(-1) = DKfi_qq_xe(-1) + L*DKfi_qq_xe(-2)
    DKfi_qq_xe(0)  = DKfi_qq_xe(0) + 0.5d0*L2*DKfi_qq_xe(-2) + L*DKfi_qq_xe(-1)

    DKfi_qq  = DKfi_qq_xe(xe)
RETURN
END FUNCTION







!wrong!FUNCTION DKfi_gg(xe,m2,sijk)!  ask Kirill
!wrong!implicit none
!wrong!real(8) :: DKfi_gg,DKfi_gg_xe(-2:0),m2,sijk,L,L2
!wrong!real(8):: r,r2,r4,r6,r8,r10
!wrong!real(8) :: x1,x2,x3,x4,x5,x6,x7,x8,x9
!wrong!integer :: xe
!wrong!
!wrong!
!wrong!    r2 = 1d0 + sijk/m_Top**2
!wrong!    r = dsqrt(r2)
!wrong!
!wrong!    r4=r2**2
!wrong!    r6=r4*r2
!wrong!    r8=r6*r2
!wrong!    r10=r8*r2
!wrong!    x1=alpha_DKTfi
!wrong!    x2=x1*x1
!wrong!    x3=x2*x1
!wrong!    x4=x3*x1
!wrong!    x5=x4*x1
!wrong!    x6=x5*x1
!wrong!    x7=x6*x1
!wrong!    x8=x7*x1
!wrong!    x9=x8*x1
!wrong!
!wrong!    L = dlog(m_Top**2/MuRen**2)
!wrong!    L2= L**2
!wrong!
!wrong!
!wrong!
!wrong!!         DKfi_gg =  0.5d0*epinv*epinv2 &
!wrong!!             +(17d0/12d0-dlog(1d0-r2)-1d0/2d0*L)*epinv &
!wrong!        DKfi_gg_xe(-2) = 0.5d0
!wrong!        DKfi_gg_xe(-1) = 17d0/12d0  -dlog(1d0-r2)-1d0/2d0*L!   maybe the log terms are wrong by a factor of 1/2
!wrong!        DKfi_gg_xe(0)  = DLi2(1d0-r2)-5d0/12d0*Pi**2+1d0/6d0*r2*(-48d0*r6*x1 &
!wrong!                      +12d0*r8*x1+12d0*x1-36d0*x1*r2+60d0*r4*x1-6d0*r4*x3 &
!wrong!                      +15d0*x2*r4+65d0*r4-6d0*x3*r2+3d0*x2*r2-45d0*r2 &
!wrong!                      +11d0*r8-43d0*r6+12d0)/(-1d0+r)**5/(r+1d0)**5*dlog(r) &
!wrong!                      -dlog(x1)**2+1d0/12d0*(-1d0+x1)*(2d0*x2-x1+23d0)*dlog(x1) &
!wrong!                      +dlog(1d0-r2)**2+(L-17d0/6d0)*dlog(1d0-r2) &
!wrong!                      -1d0/4d0*r2*(-1d0+x1)*(3d0*r4*x1+23d0*r4-2d0*x2*r4 &
!wrong!                      -13d0*r2-2d0*x2*r2-x1*r2+4d0-16d0*r6 &
!wrong!                      +4d0*r8)/(-1d0+r)**5/(r+1d0)**5*dlog(1d0-x1+x1*r2) &
!wrong!                      -1d0/240d0*(-910d0+67d0*x2+1700d0*L*x1*r2+602d0*r6*x2-60d0*L2*r10*x1 &
!wrong!                      -393d0*r8*x3+1062d0*r6*x3+9520d0*r4*x1-263d0*r8*x2-9520d0*r6*x1 &
!wrong!                      +4630d0*r8*x1-4510d0*x1*r2+340d0*L*r10*x1+8d0*x9*r10+40d0*x9*r2 &
!wrong!                      -80d0*x9*r4+80d0*x9*r6-340d0*L*x1-910d0*r10*x1+60d0*L2*x1-40d0*x9*r8 &
!wrong!                      -300d0*L2*x1*r2-3400d0*L*r4*x1+600d0*L2*r4*x1+3400d0*L*r6*x1-600d0*L2*r6*x1 &
!wrong!                      -1700d0*L*r8*x1+300d0*L2*r8*x1+348d0*x8*r4-42d0*x8*r10+198d0*x8*r8 &
!wrong!                      -162d0*x8*r2+272d0*x2*r4-358d0*x2*r2-372d0*x8*r6+2040d0*L*r4 &
!wrong!                      +115d0*x7*r10+910d0*x1+363d0*x7*r2-802d0*x7*r4+898d0*x7*r6-8d0*x9+30d0*x8 &
!wrong!                      -67d0*x7+141d0*x6-241d0*x5+283d0*x4 &
!wrong!                      -507d0*x7*r8-5540d0*r4+1520d0*x6*r4-216d0*x6*r10+3680d0*r2+915d0*x6*r8 &
!wrong!                      -730d0*x6*r2-1630d0*x6*r6+240d0*L2*r2+340d0*L+305d0*x5*r10-2516d0*x5*r4 &
!wrong!                      +1219d0*x5*r2+3680d0*r6-910d0*r8-205d0*x3+2524d0*x5*r6 &
!wrong!                      -1291d0*x5*r8-360d0*L2*r4-1360d0*L*r2-1360d0*x4*r2-1360d0*L*r6+3050d0*r4*x4 &
!wrong!                      -278d0*x4*r10-2650d0*r6*x4+240d0*L2*r6+1195d0*r8*x4+40d0*r10*x2 &
!wrong!                      -2078d0*r4*x3+917d0*x3*r2+340d0*L*r8+97d0*x3*r10-60d0*L2*r8 &
!wrong!                      -60d0*L2)/(1d0-x1+x1*r2)/(-1d0+r)**4/(r+1d0)**4
!wrong!
!wrong!
!wrong!!       factor out (mu^2/m2)^eps
!wrong!        DKfi_gg_xe(-1) = DKfi_gg_xe(-1) - dlog(MuRen**2/m2)*DKfi_gg_xe(-2)
!wrong!        DKfi_gg_xe(0)  = DKfi_gg_xe(0) - 0.5d0*dlog(MuRen**2/m2)**2*DKfi_gg_xe(-2) - dlog(MuRen**2/m2)*DKfi_gg_xe(-1)
!wrong!
!wrong!        DKfi_gg = DKfi_gg_xe(xe)
!wrong!RETURN
!wrong!END FUNCTION



! THIS ROUTINE WAS CHANGED BY SCHARF 08/16/2011 BECAUSE HE THINKS mu-dependence is wrong. ORIGINAL ABOVE  !!!!
! CORRECTED HAS BEEN THE MU-DEPENDENCE of DKfi_gg_xe(-1)
! IS MU-DEPENDENCE of DKfi_gg_xe(0) CORRECT ???????????????????????? ASK KIRILL !!!!!!!
! WHATS ABOUT eta ?  FDH OR 't Hooft Veltman ???
FUNCTION DKfi_gg(xe,m2,mu2,sijk)!  ask Kirill
implicit none
real(8) :: DKfi_gg,DKfi_gg_xe(-2:0),m2,sijk,L,L2, mu2
real(8):: r,r2,r4,r6,r8,r10
real(8) :: x1,x2,x3,x4,x5,x6,x7,x8,x9
integer :: xe


    r2 = 1d0 + sijk/m_Top**2
    r = dsqrt(r2)

    r4=r2**2
    r6=r4*r2
    r8=r6*r2
    r10=r8*r2
    x1=alpha_DKTfi
    x2=x1*x1
    x3=x2*x1
    x4=x3*x1
    x5=x4*x1
    x6=x5*x1
    x7=x6*x1
    x8=x7*x1
    x9=x8*x1

    L = dlog(m2/mu2)
    L2= L**2



!         DKfi_gg =  0.5d0*epinv*epinv2 &
!             +(17d0/12d0-dlog(1d0-r2)-1d0/2d0*L)*epinv &
        DKfi_gg_xe(-2) = 0.5d0
        DKfi_gg_xe(-1) = 17d0/12d0  -dlog(1d0-r2)-1d0/2d0*L
        DKfi_gg_xe(0)  = DLi2(1d0-r2)-5d0/12d0*Pi**2+1d0/6d0*r2*(-48d0*r6*x1 &
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

! THIS ROUTINE WAS CHANGED BY SCHARF 08/16/2011 BECAUSE HE THINKS mu-dependence is wrong. ORIGINAL ABOVE  !!!!
!       factor out (mu^2/m2)^eps
        DKfi_gg_xe(-1) = DKfi_gg_xe(-1)
! CORRECTED HAS BEEN THE MU-DEPENDENCE of DKfi_gg_xe(-1)
! IS MU-DEPENDENCE of DKfi_gg_xe(0) CORRECT ???????????????????????? YES IT IS ALREADY, therefore no correction terms as above
        DKfi_gg_xe(0)  = DKfi_gg_xe(0)
! IS MU-DEPENDENCE of DKfi_gg_xe(0) CORRECT ????????????????????????

        DKfi_gg = DKfi_gg_xe(xe)
RETURN
END FUNCTION






SUBROUTINE DKDip_FFmapping(MomExt,i,j,k,y,z,MomExtTd)!  i has to be smaller than j
implicit none
real(8) :: MomExt(1:4,1:5)!   order: Bot,Lep,Lep,Glu,Glu
real(8) :: MomExtTd(1:4,1:4)! order: Bot,Lep,Lep,Glu
real(8) :: y,z
integer :: i,j,k

          MomExtTd(1:4,1:4)=MomExt(1:4,1:4)

          y = (MomExt(1:4,i).dot.MomExt(1:4,j))/((MomExt(1:4,i).dot.MomExt(1:4,j))+(MomExt(1:4,i).dot.MomExt(1:4,k))+(MomExt(1:4,j).dot.MomExt(1:4,k)))
          z = (MomExt(1:4,i).dot.MomExt(1:4,k))/((MomExt(1:4,i).dot.MomExt(1:4,k))+(MomExt(1:4,k).dot.MomExt(1:4,j)))
          if( k.eq.5) then! this means i<j<k
              MomExtTd(1:4,4) = 1d0/(1d0-y)*MomExt(1:4,k)
          else
              MomExtTd(1:4,k) = 1d0/(1d0-y)*MomExt(1:4,k)
          endif
          MomExtTd(1:4,i) = MomExt(1:4,i)+MomExt(1:4,j)-y/(1d0-y)*MomExt(1:4,k)


! print *, "FF Mapping",i,j,k
! print *, dsqrt((MomExtTd(1:4,5)+MomExtTd(1:4,6)+MomExtTd(1:4,7)+MomExtTd(1:4,8)).dot.(MomExtTd(1:4,5)+MomExtTd(1:4,6)+MomExtTd(1:4,7)+MomExtTd(1:4,8)))
! print *, dsqrt((MomExtTd(1:4,6)+MomExtTd(1:4,7)).dot.(MomExtTd(1:4,6)+MomExtTd(1:4,7)))
! print *, (MomExt(1:4,5)+MomExt(1:4,6)+MomExt(1:4,7)+MomExt(1:4,8)+MomExt(1:4,9))-(MomExtTd(1:4,5)+MomExtTd(1:4,6)+MomExtTd(1:4,7)+MomExtTd(1:4,8))
! print *, MomExtTd(1:4,5).dot.MomExtTd(1:4,5)
! print *, MomExtTd(1:4,6).dot.MomExtTd(1:4,6)
! print *, MomExtTd(1:4,7).dot.MomExtTd(1:4,7)
! print *, MomExtTd(1:4,8).dot.MomExtTd(1:4,8)
! print *, "----"
! pause

END SUBROUTINE



SUBROUTINE DKDip_FImapping(MomExt,i,j,z,MomExtTd)
implicit none
real(8) :: MomExt(1:4,1:5)!   order: Bot,Lep,Lep,Glu,Glu
real(8) :: MomExtTd(1:4,1:4)! order: Bot,Lep,Lep,Glu
real(8) :: MomTop(1:4),MomK(1:4)
real(8) :: z
integer :: i,j,n


         MomTop(1:4)=0d0; MomK(1:4) = 0d0
         do n=1,5
              MomTop(1:4) = MomTop(1:4) + MomExt(1:4,n)
              if( n.ne.i .and. n.ne.j ) MomK(1:4) = MomK(1:4) + MomExt(1:4,n)
         enddo
         z = 1d0-(    (MomTop(1:4).dot.MomExt(1:4,j))/((MomExt(1:4,i).dot.MomTop(1:4))+(MomTop(1:4).dot.MomExt(1:4,j))-(MomExt(1:4,i).dot.MomExt(1:4,j)))     )

         MomExtTd(1:4,min(i,j)) = MomTop(1:4)
         do n=1,5
            if( n.ne.i .and. n.ne.j ) then
                if( n.ne.5 ) then
                      call DipoleTrafo(MomTop(1:4),MomK(1:4),MomExt(1:4,n),MomExtTd(1:4,n))
                      MomExtTd(1:4,min(i,j)) = MomExtTd(1:4,min(i,j)) - MomExtTd(1:4,n)
                else
                      call DipoleTrafo(MomTop(1:4),MomK(1:4),MomExt(1:4,n),MomExtTd(1:4,4))
                      MomExtTd(1:4,min(i,j)) = MomExtTd(1:4,min(i,j)) - MomExtTd(1:4,4)
                endif
            endif
         enddo

! print *, "FI Mapping",i,j
! print *, dsqrt((MomExtTd(1:4,5)+MomExtTd(1:4,6)+MomExtTd(1:4,7)+MomExtTd(1:4,8)).dot.(MomExtTd(1:4,5)+MomExtTd(1:4,6)+MomExtTd(1:4,7)+MomExtTd(1:4,8)))
! print *, dsqrt((MomExtTd(1:4,6)+MomExtTd(1:4,7)).dot.(MomExtTd(1:4,6)+MomExtTd(1:4,7)))
! print *, (MomExt(1:4,5)+MomExt(1:4,6)+MomExt(1:4,7)+MomExt(1:4,8)+MomExt(1:4,9))-(MomExtTd(1:4,5)+MomExtTd(1:4,6)+MomExtTd(1:4,7)+MomExtTd(1:4,8))
! print *, MomExtTd(1:4,5).dot.MomExtTd(1:4,5)
! print *, MomExtTd(1:4,6).dot.MomExtTd(1:4,6)
! print *, MomExtTd(1:4,7).dot.MomExtTd(1:4,7)
! print *, MomExtTd(1:4,8).dot.MomExtTd(1:4,8)
! print *, "----"

END SUBROUTINE



END MODULE
