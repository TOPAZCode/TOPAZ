MODULE ModAmplitudes
implicit none

contains





SUBROUTINE InitTrees(NumQuarks,NumGluons,NumTrees,TreeAmpsDip,NumBoson) ! NEW
use ModMisc
use ModProcess
implicit none
integer :: NumQuarks,NumGluons,NumTrees
integer,optional :: NumBoson
type(TreeProcess) :: TreeAmpsDip(:)
integer :: NumParticles,iTree


  NumParticles=NumQuarks+NumGluons
  if( present(NumBoson) ) NumParticles=NumParticles+NumBoson
  do iTree=1,NumTrees
      TreeAmpsDip(iTree)%NumPart=NumParticles
      TreeAmpsDip(iTree)%NumQua=NumQuarks
      allocate( TreeAmpsDip(iTree)%NumGlu(0:NumQuarks)     )
      allocate( TreeAmpsDip(iTree)%PartRef(1:NumParticles) )
      allocate( TreeAmpsDip(iTree)%PartType(1:NumParticles))
      allocate( TreeAmpsDip(iTree)%Quarks(1:NumQuarks)     )
      allocate( TreeAmpsDip(iTree)%Gluons(1:NumGluons)     )
      TreeAmpsDip(iTree)%NumGlu(0) = NumGluons
  enddo

RETURN
END SUBROUTINE InitTrees





SUBROUTINE LinkTreeParticles(TheTreeAmp,TheParticles) ! NEW
use ModProcess
use ModMisc
use ModParameters
implicit none
type(TreeProcess) :: TheTreeAmp
type(Particle),target :: TheParticles(:)
integer :: iPart,PartRef,PartType,ig,iq,ib,NPart,counter,QuarkPos(1:6)



           counter = 0
           do NPart=1,TheTreeAmp%NumPart
                 TheTreeAmp%PartType(NPart) = TheParticles( TheTreeAmp%PartRef(NPart) )%PartType
                 if( IsAQuark( TheTreeAmp%PartType(NPart) ) ) then  ! not a gluon and not a boson
!                      TheTreeAmp%NumQua = TheTreeAmp%NumQua + 1    ! this is suppoed to be done outside this subroutine
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
           if( ib.ge.2 ) call Error("Too many bosons in LinkTreeParticles")
           if( abs(PartType).eq.abs(Wp_) ) TheTreeAmp%NumW = TheTreeAmp%NumW + 1
           if( abs(PartType).eq.abs(Z0_) ) TheTreeAmp%NumV = TheTreeAmp%NumV + 1
           if( abs(PartType).eq.abs(Pho_)) TheTreeAmp%NumV = TheTreeAmp%NumV + 1
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





SUBROUTINE new_calc_ampl(Dv,Ds,tag_f,TreeProc,Res)
use ModMyRecurrence
use ModMyWeylRecurrence
use ModProcess
use ModParameters
implicit none
integer :: Dv,Ds,tag_f,n
complex(8) :: Res(1:Ds)
type(TreeProcess) :: TreeProc
logical :: Boson
integer :: i,j,Order(1:6)


      call setDim(Dv,Ds)
      if ( TreeProc%NumQua+TreeProc%NumSca.eq.0 ) then!  no quarks and scalars
          Res(1:Dv) = cur_g(TreeProc%Gluons(2:TreeProc%NumGlu(0)),TreeProc%NumGlu(0))
!----------------------------------------
      elseif( TreeProc%NumQua.eq.2 .and. TreeProc%NumSca.eq.0 ) then!  2 quarks and no scalars
          if( TreeProc%NumQua+TreeProc%NumGlu(0).eq.TreeProc%NumPart ) then
                Boson=.false.
          else
                Boson=.true.
          endif
          if( TreeProc%PartType(1).eq.Glu_ ) then
              Res(1:Dv) = cur_g_2f( TreeProc%Gluons(2:TreeProc%NumGlu(0)),TreeProc%Quarks(1:TreeProc%NumQua),TreeProc%NumGlu(0:3) )
! print *, "CHECK POINT cur_g_2f"
! print *, sc_(Res(1:Dv),SumMom(TreeProc%Gluons(2:TreeProc%NumGlu(0)),1,TreeProc%NumGlu(0)-1) + TreeProc%Quarks(1)%Mom+ TreeProc%Quarks(2)%Mom)
! pause
          elseif( IsAQuark(TreeProc%PartType(1)) .and. .not.Boson ) then
             Res(1:Ds) = cur_f_2f( TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Quarks(2:2),TreeProc%Quarks(1)%PartType,TreeProc%NumGlu(0:2) ) 
!             Res(1:Ds) = cur_f_2f_new( TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Quarks(2:2),TreeProc%Quarks(1)%PartType,TreeProc%NumGlu(0:2) )
          elseif( IsAQuark(TreeProc%PartType(1)) .and. Boson ) then
             if( TreeProc%NumV.eq.1 ) then
                Res(1:Ds) = cur_f_2fV(TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Quarks(2:2),TreeProc%Quarks(1)%PartType,TreeProc%Boson,TreeProc%NumGlu(0:2))
             elseif( TreeProc%NumW.eq.1 ) then!  this should only be used in top quark decay (because it is a Weyl current)
!                 Res(1:Ds) = cur_f_2fW( TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Quarks(1:2),TreeProc%Boson,TreeProc%NumGlu(0:2) )
                Res(1:Ds) = cur_f_2fW_WEYL( TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Quarks(1:2),TreeProc%Boson,TreeProc%NumGlu(0:2) )! this should be default
             else
                call Error("requested current with a boson is not available")
             endif
          else
             call Error("requested current is not available 2q")
          endif
!----------------------------------------
      elseif( TreeProc%NumSca.eq.2 .and. TreeProc%NumQua.eq.0 ) then!  2 scalars and no quarks
          if( TreeProc%PartType(1).eq.Glu_ ) then
             Res(1:Dv) = cur_g_2s(TreeProc%Gluons(2:TreeProc%NumGlu(0)),TreeProc%Scalars(1:TreeProc%NumSca),TreeProc%NumGlu(0:3))
! print *, "CHECK POINT cur_g_2s"
! print *, sc_(Res(1:Dv),SumMom(TreeProc%Gluons(2:TreeProc%NumGlu(0)),1,TreeProc%NumGlu(0)-1) + TreeProc%Scalars(1)%Mom+ TreeProc%Scalars(2)%Mom)
          elseif( IsAScalar(TreeProc%PartType(1)) ) then
             Res(1) = cur_s_2s( TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Scalars(2:2),TreeProc%NumGlu(0:2) )
             Res(2:Ds) = 0d0
          else
             call Error("requested current is not available 2s")
          endif
!----------------------------------------
      elseif( TreeProc%NumQua.eq.2 .and. TreeProc%NumSca.eq.2 ) then!  2 quarks and 2 scalars
          j=1;
          do i=1,TreeProc%NumPart-1
              if( TreeProc%PartType(i).eq.Glu_ ) cycle
              Order(j)=abs( TreeProc%PartType(i) )
              j=j+1
          enddo

          if( TreeProc%PartType(1).eq.Glu_ ) then
                if( IsAQuark(Order(1)) .and. IsAQuark(Order(2)) .and. IsAScalar(Order(3)) ) then 
!                 elseif( TreeProc%Quarks(1)%ExtRef.lt.TreeProc%Scalars(1)%ExtRef .and. (TreeProc%Quarks(2)%ExtRef.lt.TreeProc%Scalars(2)%ExtRef .or. TreeProc%Scalars(2)%ExtRef.eq.-1) )then 
                      Res(1:Dv) = cur_g_ffss(TreeProc%Gluons(2:TreeProc%NumGlu(0)),TreeProc%Scalars(1:2),TreeProc%Quarks(1:2),TreeProc%NumGlu(0:5))
! print *, "CHECK POINT cur_g_ffss"
! print *, sc_(Res(1:Dv),SumMom(TreeProc%Gluons(2:TreeProc%NumGlu(0)),1,TreeProc%NumGlu(0)-1) + TreeProc%Quarks(1)%Mom+ TreeProc%Quarks(2)%Mom+ TreeProc%Scalars(1)%Mom+ TreeProc%Scalars(2)%Mom)
! pause
!                       print *, "calling cur_g_ffss"
!                       pause
                elseif( IsAScalar(Order(1)) .and. IsAScalar(Order(2)) .and. IsAQuark(Order(3)) ) then 
!                 elseif( TreeProc%Quarks(1)%ExtRef.gt.TreeProc%Scalars(1)%ExtRef .and. (TreeProc%Quarks(2)%ExtRef.gt.TreeProc%Scalars(2)%ExtRef .or. TreeProc%Quarks(2)%ExtRef.eq.-1) )then 
                      Res(1:Dv) = cur_g_ssff(TreeProc%Gluons(2:TreeProc%NumGlu(0)),TreeProc%Scalars(1:2),TreeProc%Quarks(1:2),TreeProc%NumGlu(0:5))
! print *, "CHECK POINT cur_g_ssff"
! print *, sc_(Res(1:Dv),SumMom(TreeProc%Gluons(2:TreeProc%NumGlu(0)),1,TreeProc%NumGlu(0)-1) + TreeProc%Quarks(1)%Mom+ TreeProc%Quarks(2)%Mom+ TreeProc%Scalars(1)%Mom+ TreeProc%Scalars(2)%Mom)
! pause
!                       print *, "calling cur_g_ssff"
!                       pause
                elseif( IsAScalar(Order(1)) .and. IsAQuark(Order(2)) .and. IsAQuark(Order(3)) ) then 
                      Res(1:Dv) = cur_g_sffs(TreeProc%Gluons(2:TreeProc%NumGlu(0)),TreeProc%Scalars(1:2),TreeProc%Quarks(1:2),TreeProc%NumGlu(0:5))
!                       print *, "calling cur_g_sffs"
!                       print *, "check gauge",sc_(Res(1:Dv),SumMom(TreeProc%Gluons(2:TreeProc%NumGlu(0)),1,TreeProc%NumGlu(0)-1) + TreeProc%Quarks(1)%Mom+ TreeProc%Quarks(2)%Mom+ TreeProc%Scalars(1)%Mom+ TreeProc%Scalars(2)%Mom)
                else
                      call Error("requested current is not available g_2q+2s")
                endif

          elseif( IsAScalar(TreeProc%PartType(1))  ) then
!                 if( TreeProc%Quarks(1)%ExtRef.lt.TreeProc%Scalars(2)%ExtRef .or. TreeProc%Scalars(2)%ExtRef.eq.-1 ) then 
                if( IsAQuark(Order(2)) ) then 
!                     print *, "calling cur_s_sffs"
                    Res(1) = cur_s_sffs( TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Scalars(2:2),TreeProc%Quarks(1:2),TreeProc%NumGlu(0:4))
                    Res(2:Ds) = 0d0
                else
                    Res(1) = cur_s_ssff( TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Scalars(2:2),TreeProc%Quarks(1:2),TreeProc%NumGlu(0:4))
!                     print *, "calling cur_s_ssff"
                    Res(2:Ds) = 0d0
                endif

          elseif( IsAQuark(TreeProc%PartType(1))  ) then
!                 if( TreeProc%Scalars(1)%ExtRef.lt.TreeProc%Quarks(2)%ExtRef .or. TreeProc%Quarks(2)%ExtRef.eq.-1 ) then 
                if( IsAScalar(Order(2)) ) then 
                    Res(1:Ds) = cur_f_fssf( TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Scalars(1:2),TreeProc%Quarks(2:2),TreeProc%NumGlu(0:4))
!                     print *, "calling cur_f_fssf"
                else
                    Res(1:Ds) = cur_f_ffss( TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Scalars(1:2),TreeProc%Quarks(2:2),TreeProc%NumGlu(0:4))
!                     print *, "calling cur_f_ffss"
                endif
          else
             call Error("requested current is not available 2q+2s")
          endif

!----------------------------------------
      elseif( TreeProc%NumQua.eq.4  .and. TreeProc%NumSca.eq.0) then!  4 quarks, no scalars
          if( TreeProc%PartType(1).eq.Glu_ ) then
              Res(1:Dv) = cur_g_4f( TreeProc%Gluons(2:TreeProc%NumGlu(0)),TreeProc%Quarks(1:TreeProc%NumQua),TreeProc%NumGlu(0:5) )
! print *, "CHECK POINT cur_g_4f"
! print *, sc_(Res(1:Dv),SumMom(TreeProc%Gluons(2:TreeProc%NumGlu(0)),1,TreeProc%NumGlu(0)-1) + TreeProc%Quarks(1)%Mom+ TreeProc%Quarks(2)%Mom+ TreeProc%Quarks(3)%Mom+ TreeProc%Quarks(4)%Mom)
! pause
!                       print *, "calling cur_g_4f"
          elseif( IsAQuark(TreeProc%PartType(1)) ) then
              Res(1:Ds) = cur_f_4f( TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Quarks(2:4),TreeProc%Quarks(1)%PartType,TreeProc%NumGlu(0:4),tag_f )
!                       print *, "calling cur_f_4f"
          else
             call Error("requested current is not available 4q")
          endif
!----------------------------------------
      elseif( TreeProc%NumQua.eq.0  .and. TreeProc%NumSca.eq.4) then!  0 quarks, 4 scalars
          if( IsAScalar(TreeProc%PartType(1)) ) then
!               print *, "calling cur_s_4s",tag_f
              Res(1) = cur_s_4s( TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Scalars(2:4),TreeProc%NumGlu(0:4),tag_f )
              Res(2:Ds) = 0d0
          else
             call Error("requested current is not available 4s")
          endif
!----------------------------------------
      elseif( TreeProc%NumQua.eq.4  .and. TreeProc%NumSca.eq.2) then!  4 quarks, 2 scalars
          j=1;
          do i=1,TreeProc%NumPart-1
              if( TreeProc%PartType(i).eq.Glu_ ) cycle
              Order(j)=abs( TreeProc%PartType(i) )
              j=j+1
          enddo

          if( IsAQuark(Order(1)) .and. IsAQuark(Order(2))  ) then
              Res(1:Ds) = cur_f_fffssf( TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Scalars(1:2),TreeProc%Quarks(2:4),TreeProc%NumGlu(0:6) )
!               print *, "calling cur_f_fffssf"
          elseif( IsAQuark(Order(1)) .and. IsAScalar(Order(2))  ) then
              Res(1:Ds) = cur_f_fssfff( TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Scalars(1:2),TreeProc%Quarks(2:4),TreeProc%NumGlu(0:6) )
!               print *, "calling cur_f_fssfff"
          else
             call Error("requested current is not available 4q+2s")
          endif
!----------------------------------------
      elseif( TreeProc%NumQua.eq.2  .and. TreeProc%NumSca.eq.4) then!  2 quarks, 4 scalars
          j=1;
          do i=1,TreeProc%NumPart-1
              if( TreeProc%PartType(i).eq.Glu_ ) cycle
              Order(j)=abs( TreeProc%PartType(i) )
              j=j+1
          enddo
          if( IsAScalar(Order(1)) .and. IsAScalar(Order(2)) .and. IsAQuark(Order(3))  ) then
!               print *, "calling cur_s_ssffss"
              Res(1) = cur_s_ssffss( TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Scalars(2:4),TreeProc%Quarks(1:2),TreeProc%NumGlu(0:6) )
              Res(2:Ds) = 0d0
          elseif( IsAScalar(Order(1)) .and. IsAQuark(Order(2)) .and. IsAQuark(Order(3))  ) then
              if( tag_f.ne.3 ) then 
!                   print *, "calling cur_s_sffsss"
                  Res(1) = cur_s_sffsss( TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Scalars(2:4),TreeProc%Quarks(1:2),TreeProc%NumGlu(0:6) )
                  Res(2:Ds) = 0d0
              else 
!                   print *, "calling cur_s_sffsss_CLOSEDLOOPCONTRIB"
                  Res(1) = cur_s_sffsss_CLOSEDLOOPCONTRIB( TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Scalars(2:4),TreeProc%Quarks(1:2),TreeProc%NumGlu(0:6) )
                  Res(2:Ds) = 0d0
              endif 
          elseif( IsAScalar(Order(1)) .and. IsAScalar(Order(2)) .and. IsAScalar(Order(3))  ) then
              if( tag_f.ne.3 ) then
                  Res(1) = cur_s_sssffs( TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Scalars(2:4),TreeProc%Quarks(1:2),TreeProc%NumGlu(0:6) )
!                   print *, "calling cur_s_sssffs"
                  Res(2:Ds) = 0d0
              else
!                   print *, "calling cur_s_sssffs_CLOSEDLOOPCONTRIB"
                  Res(1) = cur_s_sssffs_CLOSEDLOOPCONTRIB( TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Scalars(2:4),TreeProc%Quarks(1:2),TreeProc%NumGlu(0:6) )
                  Res(2:Ds) = 0d0
              endif
          else
             call Error("requested current is not available 2q+4s")
          endif
!----------------------------------------
      elseif( TreeProc%NumQua.eq.6  .and. TreeProc%NumSca.eq.0) then!  6 quarks, no scalars
          if( IsAQuark(TreeProc%PartType(1)) ) then
              Res(1:Ds) = cur_f_6f( TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Quarks(2:6),TreeProc%Quarks(1)%PartType,TreeProc%NumGlu(0:6),tag_f )
          else
              call Error("requested current is not available 6q")
          endif
      else
           call Error("requested current is not available xx")
      endif


return
END SUBROUTINE







SUBROUTINE new_ampl(Dv,Ds,Nj1,Nj2,POLI,BPOLF,tag_f,TreeProc,mur)
use ModProcess
use ModMisc
use ModParameters
implicit none
integer ::  Dv,Ds,Nj1,Nj2,tag_f
complex(8),target :: POLI(8,16),BPOLF(8,16)
complex(8) :: mur(10,10)
integer :: j2,j1,i
complex(8) :: vg(Dv), vs(Ds)
complex(8) :: Res(1:Ds)
type(TreeProcess) :: TreeProc

    if( TreeProc%PartType(TreeProc%NumPart).eq.Glu_ ) then
            do j2=1,Nj2
                   TreeProc%Gluons(TreeProc%NumGlu(0))%Pol => BPOLF(j2,1:Dv)
                   call new_calc_ampl(Dv,Ds,tag_f,TreeProc,Res(1:Dv))
                   if( TreeProc%PartType(1).eq.Glu_) then
                      do j1=1,Nj1
                        mur(j1,j2) = sc_(POLI(j1,1:Dv),Res(1:Dv))
                      enddo
                   elseif( IsAQuark(TreeProc%PartType(1)) ) then
                      do j1=1,Nj1
                        mur(j1,j2) = psp1_(POLI(j1,1:Ds),Res(1:Ds))
                      enddo
                   elseif( IsAScalar(TreeProc%PartType(1)) ) then
                         j1=1
                         mur(j1,j2) = POLI(j1,1)*Res(1)
                   else
                        call Error("new_ampl")
                   endif
            enddo
    elseif( IsAQuark(TreeProc%PartType(TreeProc%NumPart)) ) then
            do j2=1,Nj2
                   TreeProc%Quarks(TreeProc%NumQua)%Pol => BPOLF(j2,1:Ds)
                   call new_calc_ampl(Dv,Ds,tag_f,TreeProc,Res(1:Ds))
                   if( TreeProc%PartType(1).eq.Glu_) then
                      do j1=1,Nj1
                         mur(j1,j2) = sc_(POLI(j1,1:Dv),Res(1:Dv))
                      enddo
                   elseif( IsAQuark(TreeProc%PartType(1)) ) then
                      do j1=1,Nj1
                         mur(j1,j2) = psp1_(POLI(j1,1:Ds),Res(1:Ds))
                      enddo
                   elseif( IsAScalar(TreeProc%PartType(1)) ) then
                         j1=1
                         mur(j1,j2) = POLI(j1,1)*Res(1)
                   else
                        call Error("new_ampl")
                   endif
            enddo
    elseif( IsAScalar(TreeProc%PartType(TreeProc%NumPart)) ) then
            do j2=1,Nj2
                   TreeProc%Scalars(TreeProc%NumSca)%Pol => BPOLF(j2,1:Ds)
                   call new_calc_ampl(Dv,Ds,tag_f,TreeProc,Res(1:Ds))
                   if( TreeProc%PartType(1).eq.Glu_) then
                      do j1=1,Nj1
                         mur(j1,j2) = sc_(POLI(j1,1:Dv),Res(1:Dv))
                      enddo
                   elseif( IsAQuark(TreeProc%PartType(1)) ) then
                      do j1=1,Nj1
                         mur(j1,j2) = psp1_(POLI(j1,1:Ds),Res(1:Ds))
                      enddo
                   elseif( IsAScalar(TreeProc%PartType(1)) ) then
                         j1=1
                         mur(j1,j2) = POLI(j1,1)*Res(1)
                   else
                        call Error("new_ampl")
                   endif
            enddo
    endif



return
END SUBROUTINE







SUBROUTINE EvalTree(TheBornAmp)
use ModProcess
use ModMisc
use ModParameters
implicit none
type(BornAmplitude) :: TheBornAmp
complex(8) :: Res(1:4)

   call new_calc_ampl(4,4,0,TheBornAmp%TreeProc,Res(1:4))

   if( IsAQuark(TheBornAmp%TreeProc%PartType(1)) ) then
      TheBornAmp%Result = psp1_(Res(1:4),ExtParticle(TheBornAmp%ExtLine(1))%Pol(1:4))
   elseif( IsAScalar(TheBornAmp%TreeProc%PartType(1)) ) then
      TheBornAmp%Result = Res(1) * ExtParticle(TheBornAmp%ExtLine(1))%Pol(1)
   elseif( TheBornAmp%TreeProc%PartType(1).eq.Glu_ ) then
      TheBornAmp%Result = (Res(1:4)).dot.(ExtParticle(TheBornAmp%ExtLine(1))%Pol(1:4))
   else
      call Error("EvalTree")
   endif

return
END SUBROUTINE EvalTree



SUBROUTINE EvalTree2(TheTreeProcess,Res)
use ModProcess
use ModMisc
use ModParameters
implicit none
type(TreeProcess) :: TheTreeProcess
complex(8) :: Pol(1:4),Res

   call new_calc_ampl(4,4,0,TheTreeProcess,Pol(1:4))

   if( IsAQuark(TheTreeProcess%PartType(1)) ) then
      Res = psp1_(Pol(1:4),TheTreeProcess%Quarks(1)%Pol(1:4))  ! assumes that first particle is an Quark!!!!
   else
      call Error("EvalTree2")
   endif

return
END SUBROUTINE EvalTree2






SUBROUTINE EvalMassCTs(TheBornAmp,Mu2Ren,Res)
use ModProcess
! use ModRecurrence
use ModMyRecurrence
use ModParameters
use ModMisc
implicit none
type(BornAmplitude) :: TheBornAmp
integer :: NPart,iq,ig,ng1,i,pos1(1)
complex(8) :: p(1:4,1:NumExtParticles),sp(1:4,1:NumExtParticles)
character :: Fl(1:NumExtParticles)*3,Fl1(1:2)*3
complex(8) :: Res(-2:1),ResSpi(1:4),ResMCT
complex(8) :: pol_glu(1:4,1:NumExtParticles-2), p_glu(1:4,1:NumExtParticles-2)
complex(8) :: pol_q(1:4,1:2), p_q(1:4,1:2)
real(8) :: Mu2Ren

       if( TheBornAmp%TreeProc%NumQua.eq.2 .and. TheBornAmp%TreeProc%NumSca.eq.0 ) then
          ResSpi = cur_f_2f_massCT(TheBornAmp%TreeProc%Gluons,TheBornAmp%TreeProc%Quarks(2:2),TheBornAmp%TreeProc%Quarks(1)%PartType,TheBornAmp%TreeProc%NumGlu)
          ResMCT = psp1_(ResSpi,TheBornAmp%TreeProc%Quarks(1)%Pol(1:4)) * TheBornAmp%TreeProc%Quarks(1)%Mass
          res(-1) = res(-1) + ResMCT*3d0/2d0
          res(0)  = res(0)  + ResMCT*( 5d0/2d0-3d0*dlog(TheBornAmp%TreeProc%Quarks(1)%Mass2/Mu2Ren)*0.5d0 )

       elseif( TheBornAmp%TreeProc%NumQua.eq.4 .and. TheBornAmp%TreeProc%NumSca.eq.0 ) then
          ResSpi = cur_f_4f_massCT(TheBornAmp%TreeProc%Gluons,TheBornAmp%TreeProc%Quarks(2:4),TheBornAmp%TreeProc%Quarks(1)%PartType,TheBornAmp%TreeProc%NumGlu)
          ResMCT = psp1_(ResSpi,TheBornAmp%TreeProc%Quarks(1)%Pol(1:4)) * TheBornAmp%TreeProc%Quarks(1)%Mass
          res(-1) = res(-1) + ResMCT*3d0/2d0
          res(0)  = res(0)  + ResMCT*( 5d0/2d0-3d0*dlog(TheBornAmp%TreeProc%Quarks(1)%Mass2/Mu2Ren)*0.5d0 )

       elseif( TheBornAmp%TreeProc%NumQua.eq.0 .and. TheBornAmp%TreeProc%NumSca.eq.2 ) then
          ResSpi(1) = cur_s_2s_massCT(TheBornAmp%TreeProc%Gluons,TheBornAmp%TreeProc%Scalars(2:2),TheBornAmp%TreeProc%NumGlu)
          ResMCT = ResSpi(1) * TheBornAmp%TreeProc%Scalars(1)%Pol(1) * TheBornAmp%TreeProc%Scalars(1)%Mass2
          res(-1) = res(-1) + ResMCT*(3d0/2d0 - 1d0/2d0 )
          res(0)  = res(0)  + ResMCT*(7d0/2d0 - 3d0*dlog(TheBornAmp%TreeProc%Scalars(1)%Mass2/Mu2Ren)*0.5d0     -1d0/2d0+dlog(TheBornAmp%TreeProc%Scalars(1)%Mass2/Mu2Ren)*0.5d0 )

!           res(-1) = res(-1) + ResMCT*(3d0/2d0  )
!           res(0)  = res(0)  + ResMCT*(7d0/2d0 - 3d0*dlog(TheBornAmp%TreeProc%Scalars(1)%Mass2/Mu2Ren)*0.5d0    )

       endif

return
END SUBROUTINE EvalMassCTs




SUBROUTINE  RenormalizeUV(ThePrimAmp,TheBornAmp,Mu2)
use ModProcess
use ModParameters
use ModMisc
implicit none
type(PrimitiveAmplitude) :: ThePrimAmp
type(BornAmplitude) :: TheBornAmp
real(8) :: Mu2,beta0,massCT
integer :: q


   if( ThePrimAmp%AmpType.eq.1 .and. ThePrimAmp%FermionLines.eq.0 .and. ThePrimAmp%ScalarLines.eq.0) then
   elseif( ThePrimAmp%AmpType.eq.1 .and. (ThePrimAmp%FermionLines.eq.1 .or. ThePrimAmp%ScalarLines.eq.1)) then
      call EvalMassCTs(TheBornAmp,Mu2,ThePrimAmp%Result(-2:1))
   elseif( ThePrimAmp%AmpType.eq.1 .and. ThePrimAmp%FermionLines.eq.2) then
      call EvalMassCTs(TheBornAmp,Mu2,ThePrimAmp%Result(-2:1))
   elseif( ThePrimAmp%AmpType.eq.3 ) then
      call EvalMassCTs(TheBornAmp,Mu2,ThePrimAmp%Result(-2:1))
      if( ThePrimAmp%ScalarLines.ne.0 ) call Error("check again if CT is correct here")
   endif


return
END SUBROUTINE RenormalizeUV





SUBROUTINE  OneLoopDiv(ThePrimAmp,mu2,Ng,rdiv2,rdiv1)
use ModProcess
use ModParameters
use ModMisc
implicit none
real(8) :: mu2,Mass_i,Mass_j
real(8),parameter :: TR=1d0
complex(8) :: rdiv2, rdiv1
complex(8) :: s_ij,beta0,dZ2_top
integer :: NPart,i,j,q,Ng
type(PrimitiveAmplitude) :: ThePrimAmp
complex(8) :: MomV(1:4,NumExtParticles)
integer :: FermionLoop


   do NPart=1,NumExtParticles
      MomV(1:4,NPart) = ExtParticle( ThePrimAmp%ExtLine(NPart) )%Mom(1:4)
   enddo

   rdiv2 = (0d0,0d0)
   rdiv1 = (0d0,0d0)
   if( ThePrimAmp%AmpType.eq.1 .and. ThePrimAmp%FermionLines.eq.0 .and. ThePrimAmp%ScalarLines.eq.0) then
      FermionLoop = 0
      do NPart=1,NumExtParticles
              i = NPart
              j = NPart+1
              if( NPart.eq.NumExtParticles ) then
                i = 1
                j = NPart
              endif
              s_ij = 2d0*sc_(MomV(1:4,i),MomV(1:4,j)) + (0d0,1d-15)
              Mass_i = ExtParticle(ThePrimAmp%ExtLine(i))%Mass
              Mass_j = ExtParticle(ThePrimAmp%ExtLine(j))%Mass
              call J_sing(s_ij,Mass_i,Mass_j,mu2,rdiv2,rdiv1)
      enddo
      do NPart=1,NumExtParticles
            call Gamma_sing(ExtParticle(ThePrimAmp%ExtLine(NPart))%PartType,FermionLoop,rdiv1)
      enddo
      q = NumExtParticles-2
      beta0 = q/2d0*11d0/3d0
      dZ2_top= 0d0

!--------------------------------------------------------------------------
   elseif( ThePrimAmp%AmpType.eq.1 .and. ThePrimAmp%FermionLines.eq.1 .and. ThePrimAmp%ScalarLines.eq.0) then
      FermionLoop = 0
      do NPart=ThePrimAmp%FermLine1Out,NumExtParticles
              i = NPart
              j = NPart+1
              if( NPart.eq.NumExtParticles ) then
                i = 1
                j = NPart
              endif
              s_ij = 2d0*sc_(MomV(1:4,i),MomV(1:4,j)) + (0d0,1d-15)
              Mass_i = ExtParticle(ThePrimAmp%ExtLine(i))%Mass
              Mass_j = ExtParticle(ThePrimAmp%ExtLine(j))%Mass
              call J_sing(s_ij,Mass_i,Mass_j,mu2,rdiv2,rdiv1)
      enddo
      do NPart=1,NumExtParticles
            call Gamma_sing(ExtParticle(ThePrimAmp%ExtLine(NPart))%PartType,FermionLoop,rdiv1)
      enddo

      q = NumExtParticles-2
      if( Ng.eq.2 ) then
         beta0 = -q/2d0*11d0/3d0   +2d0/3d0!- 4d0/3d0*TR*Nf_light
      elseif( Ng.eq.3 ) then
         beta0 = -q/2d0*11d0/3d0   +3d0/3d0!- 4d0/3d0*TR*Nf_light
      endif
      dZ2_top= 3d0/2d0
!      beta0 = 0d0
!      dZ2_top=0d0


!--------------------------------------------------------------------------
   elseif( ThePrimAmp%AmpType.eq.1 .and. ThePrimAmp%FermionLines.eq.2) then
      FermionLoop = 0
      do NPart=ThePrimAmp%FermLine1Out,ThePrimAmp%FermLine2In-1
              i = NPart
              j = NPart+1
              s_ij = 2d0*sc_(MomV(1:4,i),MomV(1:4,j)) + (0d0,1d-15)
              Mass_i = ExtParticle(ThePrimAmp%ExtLine(i))%Mass
              Mass_j = ExtParticle(ThePrimAmp%ExtLine(j))%Mass
              call J_sing(s_ij,Mass_i,Mass_j,mu2,rdiv2,rdiv1)
      enddo
      do NPart=ThePrimAmp%FermLine2Out,NumExtParticles
              i = NPart
              j = NPart+1
              if( NPart.eq.NumExtParticles ) then
                i = 1
                j = NPart
              endif
              s_ij = 2d0*sc_(MomV(1:4,i),MomV(1:4,j)) + (0d0,1d-15)
              Mass_i = ExtParticle(ThePrimAmp%ExtLine(i))%Mass
              Mass_j = ExtParticle(ThePrimAmp%ExtLine(j))%Mass
              call J_sing(s_ij,Mass_i,Mass_j,mu2,rdiv2,rdiv1)
      enddo
      do NPart=1,NumExtParticles
            call Gamma_sing(ExtParticle(ThePrimAmp%ExtLine(NPart))%PartType,FermionLoop,rdiv1)
      enddo
      q = NumExtParticles-2
      if( Ng.eq.2 ) then
        beta0 = q/2d0*11d0/3d0
      elseif( Ng.eq.3 ) then
         beta0 = q/2d0*11d0/3d0 - 10d0/3d0
      endif
      dZ2_top= 3d0/2d0

!--------------------------------------------------------------------------
   elseif( ThePrimAmp%AmpType.eq.2  .and. ThePrimAmp%ScalarLines.eq.0) then
      if( GetMass(ThePrimAmp%FermLoopPart).eq.0d0 ) then
          FermionLoop = 1
          q = NumExtParticles-2
          beta0 = -q/2d0*4d0/3d0/2d0   !minus??
      else
          FermionLoop = 1
          q = NumExtParticles-2
          beta0 = -q/2d0*4d0/3d0/2d0   !minus??
      endif
      do NPart=1,NumExtParticles
            call Gamma_sing(ExtParticle(ThePrimAmp%ExtLine(NPart))%PartType,FermionLoop,rdiv1)
      enddo
      dZ2_top= 0d0


!--------------------------------------------------------------------------
   elseif( ThePrimAmp%AmpType.eq.3  .and. ThePrimAmp%ScalarLines.eq.0) then
      FermionLoop = 0
      do NPart=ThePrimAmp%FermLine1Out,NumExtParticles
              i = NPart
              j = NPart + 1
              if( j.gt.NumExtParticles ) then
                j = j-NumExtParticles
              endif
              s_ij = 2d0*sc_(MomV(1:4,i),MomV(1:4,j)) + (0d0,1d-15)
              Mass_i = ExtParticle(ThePrimAmp%ExtLine(i))%Mass
              Mass_j = ExtParticle(ThePrimAmp%ExtLine(j))%Mass
              call J_sing(s_ij,Mass_i,Mass_j,mu2,rdiv2,rdiv1)
      enddo
      do NPart=0,NumExtParticles-ThePrimAmp%FermLine1Out+1
            i = ThePrimAmp%FermLine1Out + NPart
            if( i.gt.NumExtParticles ) then
                i = i-NumExtParticles
            endif
            call Gamma_sing(ExtParticle(ThePrimAmp%ExtLine(i))%PartType,FermionLoop,rdiv1)
      enddo

      q = NumExtParticles-ThePrimAmp%FermLine1Out
      if( Ng.eq.2 ) then
         beta0 = q/2d0*11d0/3d0
      elseif( Ng.eq.3 ) then
         beta0 = q/2d0*11d0/3d0 - q*10d0/3d0
      endif
      dZ2_top= 3d0/2d0


!--------------------------------------------------------------------------
   elseif( ThePrimAmp%AmpType.eq.4  .and. ThePrimAmp%ScalarLines.eq.0) then
      FermionLoop = 0
      do NPart=ThePrimAmp%FermLine2In,ThePrimAmp%FermLine2Out-1
              i = NPart
              j = NPart+1
              s_ij = 2d0*sc_(MomV(1:4,i),MomV(1:4,j)) + (0d0,1d-15)
              Mass_i = ExtParticle(ThePrimAmp%ExtLine(i))%Mass
              Mass_j = ExtParticle(ThePrimAmp%ExtLine(j))%Mass
              call J_sing(s_ij,Mass_i,Mass_j,mu2,rdiv2,rdiv1)
      enddo
      do NPart=ThePrimAmp%FermLine2In,ThePrimAmp%FermLine2Out
            call Gamma_sing(ExtParticle(ThePrimAmp%ExtLine(NPart))%PartType,FermionLoop,rdiv1)
      enddo

      q = (ThePrimAmp%FermLine2Out-ThePrimAmp%FermLine2In-1)
      if( Ng.eq.2 ) then
           beta0 = q/2d0*11d0/3d0
      elseif( Ng.eq.3 ) then
           beta0 = q/2d0*11d0/3d0 - q*10d0/3d0
      endif
      dZ2_top= 0d0


!--------------------------------------------------------------------------
   elseif( ThePrimAmp%AmpType.eq.1 .and. ThePrimAmp%FermionLines.eq.0 .and. ThePrimAmp%ScalarLines.eq.1) then
      FermionLoop = 0
      do NPart=ThePrimAmp%ScaLine1Out,NumExtParticles
              i = NPart
              j = NPart+1
              if( NPart.eq.NumExtParticles ) then
                i = 1
                j = NPart
              endif
              s_ij = 2d0*sc_(MomV(1:4,i),MomV(1:4,j)) + (0d0,1d-15)
              Mass_i = ExtParticle(ThePrimAmp%ExtLine(i))%Mass
              Mass_j = ExtParticle(ThePrimAmp%ExtLine(j))%Mass
              call J_sing(s_ij,Mass_i,Mass_j,mu2,rdiv2,rdiv1)
      enddo
      do NPart=1,NumExtParticles
            call Gamma_sing(ExtParticle(ThePrimAmp%ExtLine(NPart))%PartType,FermionLoop,rdiv1)
      enddo

      q = NumExtParticles-2
      if( Ng.eq.2 ) then
         beta0 = -q/2d0*11d0/3d0   +2d0/3d0!- 4d0/3d0*TR*Nf_light
      elseif( Ng.eq.3 ) then
         beta0 = -q/2d0*11d0/3d0   +3d0/3d0!- 4d0/3d0*TR*Nf_light
      endif
      dZ2_top= 0d0
!      beta0 = 0d0
!      dZ2_top=0d0



!--------------------------------needs modification------------------------------------------
   elseif( ThePrimAmp%AmpType.eq.1 .and. ThePrimAmp%FermionLines.eq.1 .and. ThePrimAmp%ScalarLines.eq.1) then
      FermionLoop = 0
      do NPart=ThePrimAmp%ScaLine1Out,ThePrimAmp%FermLine1In-1
              i = NPart
              j = NPart+1
              s_ij = 2d0*sc_(MomV(1:4,i),MomV(1:4,j)) + (0d0,1d-15)
              Mass_i = ExtParticle(ThePrimAmp%ExtLine(i))%Mass
              Mass_j = ExtParticle(ThePrimAmp%ExtLine(j))%Mass
              call J_sing(s_ij,Mass_i,Mass_j,mu2,rdiv2,rdiv1)
      enddo
      do NPart=ThePrimAmp%FermLine1Out,NumExtParticles
              i = NPart
              j = NPart+1
              if( NPart.eq.NumExtParticles ) then
                i = 1
                j = NPart
              endif
              s_ij = 2d0*sc_(MomV(1:4,i),MomV(1:4,j)) + (0d0,1d-15)
              Mass_i = ExtParticle(ThePrimAmp%ExtLine(i))%Mass
              Mass_j = ExtParticle(ThePrimAmp%ExtLine(j))%Mass
              call J_sing(s_ij,Mass_i,Mass_j,mu2,rdiv2,rdiv1)
      enddo
      do NPart=1,NumExtParticles
            call Gamma_sing(ExtParticle(ThePrimAmp%ExtLine(NPart))%PartType,FermionLoop,rdiv1)
      enddo
      q = NumExtParticles-2
      if( Ng.eq.2 ) then
        beta0 = q/2d0*11d0/3d0
      elseif( Ng.eq.3 ) then
         beta0 = q/2d0*11d0/3d0 - 10d0/3d0
      endif
      dZ2_top= 0d0


!--------------------------------needs modification-------------------------------------------
   elseif( ThePrimAmp%AmpType.eq.4  .and. ThePrimAmp%ScalarLines.eq.1) then
      FermionLoop = 0
      do NPart=ThePrimAmp%FermLine1In,ThePrimAmp%FermLine1Out-1
              i = NPart
              j = NPart+1
              s_ij = 2d0*sc_(MomV(1:4,i),MomV(1:4,j)) + (0d0,1d-15)
              Mass_i = ExtParticle(ThePrimAmp%ExtLine(i))%Mass
              Mass_j = ExtParticle(ThePrimAmp%ExtLine(j))%Mass
              call J_sing(s_ij,Mass_i,Mass_j,mu2,rdiv2,rdiv1)
      enddo
      do NPart=ThePrimAmp%FermLine1In,ThePrimAmp%FermLine1Out
            call Gamma_sing(ExtParticle(ThePrimAmp%ExtLine(NPart))%PartType,FermionLoop,rdiv1)
      enddo

      q = (ThePrimAmp%FermLine1Out-ThePrimAmp%FermLine1In-1)
      if( Ng.eq.2 ) then
           beta0 = q/2d0*11d0/3d0
      elseif( Ng.eq.3 ) then
           beta0 = q/2d0*11d0/3d0 - q*10d0/3d0
      endif
      dZ2_top= 0d0


!------------------------------needs modification---------------------------------------------
   elseif( ThePrimAmp%AmpType.eq.3  .and. ThePrimAmp%ScalarLines.eq.1) then
      FermionLoop = 0
      do NPart=ThePrimAmp%ScaLine1Out,NumExtParticles
              i = NPart
              j = NPart + 1
              if( j.gt.NumExtParticles ) then
                j = j-NumExtParticles
              endif
              s_ij = 2d0*sc_(MomV(1:4,i),MomV(1:4,j)) + (0d0,1d-15)
              Mass_i = ExtParticle(ThePrimAmp%ExtLine(i))%Mass
              Mass_j = ExtParticle(ThePrimAmp%ExtLine(j))%Mass
              call J_sing(s_ij,Mass_i,Mass_j,mu2,rdiv2,rdiv1)
      enddo
      do NPart=0,NumExtParticles-ThePrimAmp%ScaLine1Out+1
            i = ThePrimAmp%ScaLine1Out + NPart
            if( i.gt.NumExtParticles ) then
                i = i-NumExtParticles
            endif
            call Gamma_sing(ExtParticle(ThePrimAmp%ExtLine(i))%PartType,FermionLoop,rdiv1)
      enddo

      q = NumExtParticles-ThePrimAmp%ScaLine1Out

      if( Ng.eq.2 ) then
         beta0 = q/2d0*11d0/3d0
      elseif( Ng.eq.3 ) then
         beta0 = q/2d0*11d0/3d0 - q*10d0/3d0
      endif
      dZ2_top= 0d0
   endif




   rdiv1 = rdiv1 + beta0 + dZ2_top ! needs to be checked


return
END SUBROUTINE OneLoopDiv




! J_sing = (mu**2/|s_jk|)^eps * [ V^cc_jk(s_jk; m_j,m_k;eps) + 1/v_jk*(I*pi/eps - pi^2/2)*Theta[s_jk] ]
! refers to hep-ph/0011222, Eq.(9), line 2
SUBROUTINE J_sing(s_jk,Mass_j,Mass_k,mu2,rdiv2,rdiv1)
use ModParameters
implicit none
complex(8) :: s_jk,v_ij,rdiv2,rdiv1
real(8) :: Mass_j,Mass_k,mu2
real(8) :: mu2j,mu2k,rho,vij

  if( Mass_j.eq.0d0 .AND. Mass_k.eq.0d0 ) then
      rdiv2 = rdiv2 - 1d0
      rdiv1 = rdiv1 - cdlog(dcmplx(mu2)/cdabs(s_jk))
      if( dreal(s_jk).gt.0d0 ) rdiv1 = rdiv1 - dcmplx(0d0,DblPi)
  elseif( Mass_j.eq.0d0 .AND. Mass_k.ne.0d0 ) then
      rdiv2 = rdiv2 - 0.5d0
      rdiv1 = rdiv1 - 0.5d0*(cdlog(dcmplx(mu2)/cdabs(s_jk))+cdlog(dcmplx(Mass_k**2)/cdabs(s_jk)))
      if( dreal(s_jk).gt.0d0 ) rdiv1 = rdiv1 - dcmplx(0d0,DblPi)
  elseif( Mass_j.ne.0d0 .AND. Mass_k.eq.0d0 ) then
      rdiv2 = rdiv2 - 0.5d0
      rdiv1 = rdiv1 - 0.5d0*(cdlog(dcmplx(mu2)/cdabs(s_jk))+cdlog(dcmplx(Mass_j**2)/cdabs(s_jk)))
      if( dreal(s_jk).gt.0d0 ) rdiv1 = rdiv1 - dcmplx(0d0,DblPi)
  elseif( Mass_j.ne.0d0 .AND. Mass_k.ne.0d0 ) then
      v_ij = cdsqrt( 1d0-(2d0*Mass_j*Mass_k/s_jk)**2 )
      rdiv1 = rdiv1 - 0.5d0/v_ij*cdlog((1d0-v_ij)/(1d0+v_ij))
      if( dreal(s_jk).gt.0d0 ) rdiv1 = rdiv1 - dcmplx(0d0,DblPi)/v_ij
  endif

return
END SUBROUTINE


SUBROUTINE Gamma_sing(PartType,NfTerms,rdiv1)
use ModParameters
use ModMisc
implicit none
integer :: PartType
complex(8) :: rdiv1
real(8),parameter :: TR=1d0
integer :: NfTerms

    if( PartType.eq.Glu_.and. NfTerms.eq.0) then
        rdiv1 = rdiv1 - 11d0/6d0 +2d0/3d0*TR*Nf_light
    elseif( PartType.eq.Glu_ .and. NfTerms.eq.1) then
        rdiv1 = rdiv1 + 2d0/3d0/2d0
    elseif( IsAQuark(PartType) .and. NfTerms.eq.0) then
        if(GetMass(PartType).eq.0d0) then
            rdiv1 = rdiv1 - 3d0/2d0/2d0
        else
            rdiv1 = rdiv1 - 1d0/2d0
        endif
    elseif( IsAScalar(PartType) .and. NfTerms.eq.0) then
        if(GetMass(PartType).ne.0d0) then
            rdiv1 = rdiv1 - 1d0/2d0
        endif
    endif
return
END SUBROUTINE







FUNCTION CheckPoles(ThePrimAmp,TheBornAmp,rdiv)
use ModProcess
implicit none
type(PrimitiveAmplitude) :: ThePrimAmp
type(BornAmplitude) :: TheBornAmp
complex(8) :: rdiv(1:2)
real(8) :: CheckPoles

      if(cdabs(rdiv(1)).eq.0d0 .or.cdabs(TheBornAmp%Result).lt.1d-13 ) then
        CheckPoles = 0d0
        return
      endif
      CheckPoles = cdabs( ThePrimAmp%Result(-1)/TheBornAmp%Result/rdiv(1)-1d0 )

END FUNCTION






SUBROUTINE WritePrimAmpResult(ThePrimAmp,TheBornAmp,rdiv,other)
use ModProcess
implicit none
type(PrimitiveAmplitude) :: ThePrimAmp
type(BornAmplitude) :: TheBornAmp
complex(8) :: rdiv(1:2)
real(8),optional  :: other(:)

!            write(*,"(A,I2,A,10I10)") "Primitive Amplitude (type ",ThePrimAmp%AmpType,")",ThePrimAmp%ExtLine(:)
!            print *, "tree amp:",TheBornAmp%Result
!           print *, "1-Loop, res(-2):",ThePrimAmp%Result(-2)
!           print *, "1-Loop, res(-1):",ThePrimAmp%Result(-1)
!           print *, "1-Loop, res(0): ",ThePrimAmp%Result(0)
!           print *, "1-Loop, res(1): ",ThePrimAmp%Result(1)
!            print *, "norm. div. res(-2):", ThePrimAmp%Result(-2)/TheBornAmp%Result
!            print *, "check div. res(-2):", rdiv(2)
!            print *, "norm. div. res(-1):", ThePrimAmp%Result(-1)/TheBornAmp%Result
!            print *, "check div. res(-1):", rdiv(1)

!           print *, "1-Loop, rel cut: ",ThePrimAmp%Result(0)/TheBornAmp%Result
!           print *, "1-Loop, rel rat: ",ThePrimAmp%Result(1)/TheBornAmp%Result
!           print *, "1-Loop, rel fin: ",(ThePrimAmp%Result(0)+ThePrimAmp%Result(1))/TheBornAmp%Result
!---------
!           write(*,"(A,I2,A,10I10)") "Primitive Amplitude (type ",ThePrimAmp%AmpType,")",ThePrimAmp%ExtLine(:)
!           if(rdiv(2).ne.(0d0,0d0) ) then
!             print *, "precision res(-2):",ThePrimAmp%Result(-2)/TheBornAmp%Result/rdiv(2) - 1d0
!           else
!             print *, "precision res(-2): vanishing"
!           endif
!           if(rdiv(1).ne.(0d0,0d0) ) then
!             print *, "precision res(-1):",ThePrimAmp%Result(-1)/TheBornAmp%Result/rdiv(1) - 1d0
!           else
!             print *, "precision res(-1): vanishing"
!           endif
!---------
!           if(rdiv(2).ne.(0d0,0d0) ) then
!             if(cdabs(ThePrimAmp%Result(-2)/TheBornAmp%Result/rdiv(2) - 1d0).gt.1d-9) then
              print *, ""
              write(*,"(A,10I10)") "Primitive Amplitude ",ThePrimAmp%ExtLine(1:NumExtParticles)
!               print *, "precision res(-2):",ThePrimAmp%Result(-2)/TheBornAmp%Result/rdiv(2) - 1d0
              print *, "norm. div. res(-2):", ThePrimAmp%Result(-2)/TheBornAmp%Result
              print *, "check div. res(-2):", rdiv(2)
!             endif
!           endif

!           if( rdiv(1).ne.(0d0,0d0) ) then
!             if(cdabs(ThePrimAmp%Result(-1)/TheBornAmp%Result/rdiv(1) - 1d0).gt.1d0) then
              print *, ""
              write(*,"(A,10I10)") "Primitive Amplitude ",ThePrimAmp%ExtLine(1:NumExtParticles)
!               print *, "rel. dev(-1):",cdabs(ThePrimAmp%Result(-1)/TheBornAmp%Result/rdiv(1) - 1d0)
              print *, "norm. div. res(-1):", ThePrimAmp%Result(-1)/TheBornAmp%Result
              print *, "check div. res(-1):", rdiv(1)
!               print *, "norm. div. res(0+1):", (ThePrimAmp%Result(0)+ThePrimAmp%Result(1))/TheBornAmp%Result
!               print *, "tree", TheBornAmp%Result
!               print *, other(:)
!            endif
!           endif

END SUBROUTINE





SUBROUTINE WriteLatexOutput(iPrimAmp,PrimAmps,BornAmps)
use ModProcess
use ModParameters
implicit none
type(PrimitiveAmplitude) :: PrimAmps
type(BornAmplitude) :: BornAmps
integer :: iPrimAmp,i,j
character :: csign(1:5)*(5)
character :: ctype(1:5)*(10),iPart(1:5)*(1),sign1*(1),sign2*(1)
character :: tex(1:5)*(100)


  do i=1,5
    if(ExtParticle(PrimAmps%ExtLine(i))%Helicity.eq.+1) csign(i)="+"
    if(ExtParticle(PrimAmps%ExtLine(i))%Helicity.eq.-1) csign(i)="-"

    if(ExtParticle(PrimAmps%ExtLine(i))%PartType.eq.Glu_) ctype(i)="g"
    if(ExtParticle(PrimAmps%ExtLine(i))%PartType.eq.Top_)  ctype(i)="t"
    if(ExtParticle(PrimAmps%ExtLine(i))%PartType.eq.ATop_)  ctype(i)="\bar{t}"
    if(ExtParticle(PrimAmps%ExtLine(i))%PartType.eq.Str_)  ctype(i)="q"
    if(ExtParticle(PrimAmps%ExtLine(i))%PartType.eq.AStr_)  ctype(i)="\bar{q}"
    write(iPart(i),"(I1)") PrimAmps%ExtLine(i)
    tex(i) = trim(iPart(i))//"_{"//trim(ctype(i))//"}^{"//trim(csign(i))//"}"
  enddo

  sign1=""
  sign2=""
  if(dimag(PrimAmps%Result(-1)/BornAmps%Result).ge.0d0) sign1="+"
  if(dimag((PrimAmps%Result(0)+PrimAmps%Result(1))/BornAmps%Result).ge.0d0) sign2="+"
  write(*,"(A,F9.6,A,F13.8,A,F13.8,A,F13.8,A,F13.8,A)") "$A("//trim(tex(1))//","//trim(tex(2))//","//trim(tex(3))//","//trim(tex(4))//","//trim(tex(5))//")"// &
                 "$ & $",dreal(PrimAmps%Result(-2)/BornAmps%Result), &
                 "$ & $",dreal(PrimAmps%Result(-1)/BornAmps%Result),sign1,dimag(PrimAmps%Result(-1)/BornAmps%Result)," \mathrm{i} $ & $", &
                      dreal((PrimAmps%Result(0)+PrimAmps%Result(1))/BornAmps%Result),sign2,dimag((PrimAmps%Result(0)+PrimAmps%Result(1))/BornAmps%Result)," \mathrm{i}$ \\"

!   print *, "---",iPrimAmp
!   sign1=""
!   if(dimag(BornAmps%Result).ge.0d0) sign1="+"
!   write(*,"(F13.8,A,F13.8,A)") dreal(BornAmps%Result),sign1,dimag(BornAmps%Result),"\mathrm{i}"

END SUBROUTINE








!       SUBROUTINE calc_ampl(TreeProc,Dv,Ds,ng,nq,nw,sp,p,fl,tag_f,vg,vs)
!       use ModRecurrence
!       use ModMyRecurrence
!       use ModProcess
!       use ModParameters
!       implicit none
!       complex(8), intent(in) :: sp(:,:), p(:,:)
!       character :: fl(:)*3
!       ! RESTORE INTENT IN
!       integer :: Dv,Ds,ng,nq,nw,tag_f
!       integer :: npa
!       integer :: i, ig,iq
!       complex(8), intent(out) :: vg(Dv),vs(Ds)
!       complex(8) :: pol_glu(Dv,1:ng),vgnew(Dv),vsnew(Ds),vsold(Ds),check
!       complex(8) :: p_glu(Dv,1:ng)
!       complex(8) :: pol_q(Ds,1:nq),tmp(1:ds)
!       complex(8) :: p_q(Dv,1:nq)
!       complex(8) :: p_W(Dv), pol_W(Dv)
!       character :: fl1(1:nq)*3
!       integer :: pos1(1),pos2(2),pos3(4),pos5(5)
!       integer :: ng1,ng2, ng3, ng4,ng5, iW, sw
!       integer :: calc,flcheck
!       type(TreeProcess) :: TreeProc
!       type(PtrToParticle) :: TmpQuark
!
!       include 'global_import'
!
!       calc = 0
!       npa=size(p,dim=2)
!
!       call SetDim(Dv,Ds)
!
! !________________________________________________________________
!
!       if (fl(1) == 'glu'.and.nq == 0) then
!           pol_glu(:,:) = sp(1:Dv,:)
!           p_glu = p
!           vg=g(pol_glu(:,2:ng),p_glu(:,2:ng))
!           calc = 1
!
! !DEC$ IF (_DebugCheckMyImpl1==1)
!           vgnew = cur_g(TreeProc%Gluons(2:TreeProc%NumGlu(0)),TreeProc%NumGlu(0))
!           do iq=1,Dv
!              check = vg(iq)-vgnew(iq)
!           enddo
!           if(abs(check).gt.1d-10) then
!             print *, "error in cur_g"
!           endif
! !DEC$ ENDIF
!
!       endif
! !________________________________________________________________
!
!       if (fl(1) == 'glu'.and.nq == 2) then
!            iq = 0
!            ig = 0
!            do i=2,npa
!                if (fl(i) == 'glu') then
!                   ig = ig + 1
!                   pol_glu(:,ig) = sp(1:Dv,i)
!                   p_glu(:,ig) = p(:,i)
!                endif
!                if (fl(i) == 'top'.or.fl(i) == 'bot'.or.fl(i)=='str') then
!                   iq = iq + 1
!                   pol_q(:,iq) = sp(:,i)
!                   p_q(:,iq) = p(:,i)
!                   pos2(iq) = i
!                   fl1(iq) = fl(i)
!                endif
!             enddo
!
!             ng1 = pos2(1) -2
!             ng2 = pos2(2) - pos2(1) - 1
!             pol_glu(:,ng) = sp(1:Dv,1)
!
!             if (fl1(1) == fl1(2)) then
!                vg=g_fbf( pol_glu(:,1:ng-1),p_glu(:,1:ng-1),pol_q(:,1),p_q(:,1),fl1(1),pol_q(:,2),p_q(:,2),fl1(2),ng1,ng2 )
!
! ! CHECK
! !       if ( ng.ge.5 .and.nq.eq.2 ) then
! !         print *, "HERE: ",nq,ng,fl
! !         ng1=1
! !         ng2=0
! !         ng3=4-ng1-ng2
! !
! !         ng1=3
! !         ng2=0
! !         ng3=1
! !         ng=1+ng1+ng2+ng3
! !
! !
! !         vg=g_fbf( pol_glu(:,1:ng-1),p_glu(:,1:ng-1),pol_q(:,1),p_q(:,1),fl1(1),pol_q(:,2),p_q(:,2),fl1(2),ng1,ng2 )
! !         vgnew=cur_g_2f(TreeProc%Gluons(2:ng),TreeProc%Quarks(1:2),(/ng,ng1,ng2,ng3/))
! !
! !         do iq=1,Dv
! !            check = vg(iq)-vgnew(iq)
! !         enddo
! !         print *, "check: ",ng,check
! !         print *, ""
! !         print *, vg(:)
! !         print *, ".."
! !         print *, vgnew(:)
! !
! !         stop
! !       endif
! ! CHECK
!
!
! !DEC$ IF (_DebugCheckMyImpl1==1)
!                vgnew=cur_g_2f(TreeProc%Gluons(2:TreeProc%NumGlu(0)),TreeProc%Quarks(1:2),TreeProc%NumGlu(0:3))
!                do iq=1,Dv
!                   check = vg(iq)-vgnew(iq)
!                enddo
!                if(abs(check).gt.1d-10) print *, "error in cur_g_2fer"
! !DEC$ ENDIF
!             else
!                vg = czero
!             endif
!             calc = 1
!       endif
! !________________________________________________________________
!
!
!
!       if (((fl(1) == 'top').or.(fl(1) == 'bot').or.(fl(1) == 'str')).and.nq == 2) then
!            iq = 0
!            ig = 0
!           do i=2,npa
!                if (fl(i) == 'glu') then
!                   ig = ig + 1
!                   pol_glu(:,ig) = sp(1:Dv,i)
!                   p_glu(:,ig) = p(:,i)
!                endif
!                if (fl(i) == 'top'.or.fl(i) == 'bot'.or.fl(i)=='str') then
!                   iq = iq + 1
!                   pol_q(:,iq) = sp(:,i)
!                   p_q(:,iq) = p(:,i)
!                   fl1(iq) = fl(i)
!                   pos1(iq) = i
!                endif
!           enddo
!           ng1 = pos1(1) -2
!
!          if (fl(1) == fl1(1)) then
!             vs=f(pol_glu,p_glu,pol_q(:,1),p_q(:,1),fl1(1),fl(1),ng1)
! !DEC$ IF (_DebugCheckMyImpl1==1)
!             vsnew=cur_f_2f(TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Quarks(2:2),TreeProc%Quarks(1)%PartType,TreeProc%NumGlu(0:2))
! !             print *, "old: ",vs
! !             print *,""
! !             print *,"new: ",vsnew
!             do iq=1,Ds
!               check = vs(iq)-vsnew(iq)
!             enddo
!            if(abs(check).gt.1d-10)  print *, "error in cur_f_2f",ng1,TreeProc%NumGlu
! !DEC$ ENDIF
!
!          else
!             vs = czero
!          endif
!          calc = 1
!       endif
!
! !________________________________________________________________
!
!
!       if ((fl(1) == 'top'.or.fl(1) == 'bot'.or.fl(1) == 'str').and.nq == 4) then
!
!            iq = 1
!            ig = 0
!            pol_q(:,1) = sp(:,1)
!            p_q(:,1) = p(:,1)
!            do i=2,npa
!
!                if (fl(i) == 'glu') then
!                   ig = ig + 1
!                   pol_glu(:,ig) = sp(1:Dv,i)
!                   p_glu(:,ig) = p(:,i)
!                endif
!
!                if (fl(i) == 'top'.or.fl(i) == 'bot'.or.fl(i) == 'str') then
!                   pos3(iq) = i
!                   iq = iq + 1
!                   pol_q(:,iq) = sp(:,i)
!                   p_q(:,iq) = p(:,i)
!                   fl1(iq) = fl(i)
!                endif
!
!             enddo
!             ng1 = pos3(1) - 2
!             ng2 = pos3(2) - pos3(1) - 1
!             ng3 = pos3(3) - pos3(2) - 1
!
!
!             tmp(1:ds) = pol_q(1:ds,2)
!             pol_q(1:ds,2) = pol_q(1:ds,3)
!             pol_q(1:ds,3) = tmp(1:ds)
!             tmp(1:dv) = p_q(1:dv,2)
!             p_q(1:dv,2) = p_q(1:dv,3)
!             p_q(1:dv,3) = tmp(1:dv)
!
!             vsold=f_bffbf(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),fl1(2:nq),fl(1),tag_f,ng1,ng2,ng3)
!             vsold = -vsold
!
!             vs=cur_f_4f(TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Quarks(2:4),TreeProc%Quarks(1)%PartType,TreeProc%NumGlu(0:4),tag_f)
!
!
! !DEC$ IF (_DebugCheckMyImpl1==1)
!             do iq=1,Ds
!                   check = vsold(iq)-vs(iq)
!            enddo
!
!            if(abs(check).gt.1d-10) then
! !               print *, "error in cur_f_4f ",treeproc%parttype
! !               print *, tag_f
! !               print *, "old: ",vsold(1:Ds)
! !               print *,""
! !               print *, "new: ",vs(1:Ds)
! !               stop
!            endif
! !DEC$ ENDIF
!
!             tmp(1:ds) = pol_q(1:ds,2)
!             pol_q(1:ds,2) = pol_q(1:ds,3)
!             pol_q(1:ds,3) = tmp(1:ds)
!             tmp(1:dv) = p_q(1:dv,2)
!             p_q(1:dv,2) = p_q(1:dv,3)
!             p_q(1:dv,3) = tmp(1:dv)
!
! !             TmpQuark = TreeProc%Quarks(3)
! !             TreeProc%Quarks(3) = TreeProc%Quarks(4)
! !             TreeProc%Quarks(4) = TmpQuark
!             calc = 1
!       endif
!
! !________________________________________________________________
!
!
!        if (fl(1) == 'glu'.and.nq == 4) then
!            iq = 0
!            ig = 0
!            do i=2,npa
!                if (fl(i) == 'glu') then
!                   ig = ig + 1
!                   pol_glu(:,ig) = sp(1:Dv,i)
!                   p_glu(:,ig) = p(:,i)
!                endif
!                if (fl(i) == 'top'.or.fl(i) == 'bot'.or.fl(i)=='str') then
!                   iq = iq + 1
!                   pol_q(:,iq) = sp(:,i)
!                   p_q(:,iq) = p(:,i)
!                   pos3(iq) = i
!                   fl1(iq) = fl(i)
!                endif
!            enddo
!          ng1 = pos3(1) - 2
!          ng2 = pos3(2) - pos3(1) - 1
!          ng3 = pos3(3) - pos3(2) - 1
!          ng4 = pos3(4) - pos3(3) - 1
!          ng5 = ig-ng1-ng2-ng3-ng4
!
!          pol_glu(:,ng) = sp(1:Dv,1)
!
!          print *, "current g_sbsfbf with ",fl(1:npa)
! !          vg=g_sbsfbf(pol_glu(:,1:ng-1),p_glu(:,1:ng-1),pol_q,p_q,fl1,tag_f,ng1,ng2,ng3,ng4)
! !          print *, "old ",vg
! !          vg=new_g_fbffbf(pol_glu(:,1:ng-1),p_glu(:,1:ng-1),pol_q,p_q,fl1,tag_f,ng1,ng2,ng3,ng4,ng5)
! !          print *, "new ",vg
!
!          !stop
!
!          calc = 1
!       endif
! !________________________________________________________________
!
!
!
!       if (fl(1) == 'bot'.and.nq == 6) then
!          print *, "fl(1) == 'bot'.and.nq == 6"
!           iq = 1
!           ig = 0
!           pol_q(:,1) = sp(:,1)
!           p_q(:,1) = p(:,1)
!           do i=2,npa
!                if (fl(i) == 'glu') then
!                   ig = ig + 1
!                   pol_glu(:,ig) = sp(1:Dv,i)
!                   p_glu(:,ig) = p(:,i)
!                endif
!
!                if (fl(i) == 'top'.or.fl(i) == 'bot'.or.fl(i) == 'str') then
!                   pos5(iq) = i
!                   iq = iq + 1
!                   pol_q(:,iq) = sp(:,i)
!                   p_q(:,iq) = p(:,i)
!                   fl1(iq) = fl(i)
!                endif
!           enddo
!          ng1 = pos5(1) -2
!          ng2 = pos5(2) - pos5(1) - 1
!          ng3 = pos5(3) - pos5(2) - 1
!          ng4 = pos5(4) - pos5(3) - 1
!          vs=f_bffbffbf(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),fl1(2:nq),fl(1),tag_f,ng1,ng2,ng3,ng4)
!          calc = 1
!       endif
!
!
!       if (calc == 0) then
!          print *, 'UNCALCULATED: ', fl(1:npa)
!          !stop
!       endif
!
!
! !________________________________________________________________
!
!
! !DEC$ IF (_DebugCheckMyImpl1==1)
!           do i=1,npa
!             if (fl(i) == 'glu') flcheck=10
!             if (fl(i) == 'top') flcheck=5
!             if (fl(i) == 'bot') flcheck=6
!             if (fl(i) == 'str') flcheck=4
!             if (flcheck-abs(TreeProc%PartType(i)).ne.0) then
!                print *, "Error in flavor check"
!             endif
!             if(ng-TreeProc%NumGlu(0).ne.0 .or. nq-TreeProc%NumQua.ne.0) then
!                     print *, "Error in nq,ng check"
!             endif
!           enddo
!
!
!           if (fl(1) .ne. 'glu'.and.nq .eq. 2) then
!              if(ng1-TreeProc%NumGlu(1).ne.0..or.ng-ng1-TreeProc%NumGlu(2).ne.0) print *,"NumGlu check"
!
!              if(cdabs(p_q(1,1)-treeproc%quarks(2)%mom(1)+p_q(2,1)-treeproc%quarks(2)%mom(2) &
!                      +p_q(3,1)-treeproc%quarks(2)%mom(3)+p_q(4,1)-treeproc%quarks(2)%mom(4)).gt.1d-10) print *, "wrong momentum q"
!              do i=1,ng
!                 if(cdabs(p_glu(1,i)-treeproc%gluons(i)%mom(1)+p_glu(2,i)-treeproc%gluons(i)%mom(2) &
!                         +p_glu(3,i)-treeproc%gluons(i)%mom(3)+p_glu(4,i)-treeproc%gluons(i)%mom(4)).gt.1d-10) print *, "wrong momentum g"
!              enddo
!           endif
!
!
!
!           if (fl(1) .ne. 'glu'.and.nq .eq. 4) then
!              if(ng1-TreeProc%NumGlu(1).ne.0.or.ng2-TreeProc%NumGlu(2).ne.0.or.ng3-TreeProc%NumGlu(3).ne.0.or.ng-ng1-ng2-ng3-TreeProc%NumGlu(4).ne.0) print *,"NumGlu check"
!
!              do i=2,nq
!                 if(cdabs(p_q(1,i)-treeproc%quarks(i)%mom(1)+p_q(2,i)-treeproc%quarks(i)%mom(2) &
!                         +p_q(3,i)-treeproc%quarks(i)%mom(3)+p_q(4,i)-treeproc%quarks(i)%mom(4)).gt.1d-10) print *, "wrong momentum q",i
!              enddo
!              do i=1,ng
!                 if(cdabs(p_glu(1,i)-treeproc%gluons(i)%mom(1)+p_glu(2,i)-treeproc%gluons(i)%mom(2) &
!                         +p_glu(3,i)-treeproc%gluons(i)%mom(3)+p_glu(4,i)-treeproc%gluons(i)%mom(4)).gt.1d-10) print *, "wrong momentum g",i
!              enddo
!           endif
!
!
!
!           if (fl(1) .eq. 'glu'.and.nq .eq. 2) then
!              if(ng1-TreeProc%NumGlu(1).ne.0.or.ng2-TreeProc%NumGlu(2).ne.0.or.ng-1-ng1-ng2-TreeProc%NumGlu(3).ne.0) print *,"NumGlu check"
!              do i=2,ng
!                 if(cdabs(p_glu(1,i-1)-treeproc%gluons(i)%mom(1)+p_glu(2,i-1)-treeproc%gluons(i)%mom(2) &
!                         +p_glu(3,i-1)-treeproc%gluons(i)%mom(3)+p_glu(4,i-1)-treeproc%gluons(i)%mom(4)).gt.1d-10) print *, "wrong momentum g",i
!                 if(cdabs(pol_glu(1,i-1)-treeproc%gluons(i)%pol(1)+pol_glu(2,i-1)-treeproc%gluons(i)%pol(2) &
!                         +pol_glu(3,i-1)-treeproc%gluons(i)%pol(3)+pol_glu(4,i-1)-treeproc%gluons(i)%pol(4)).gt.1d-10) print *, "wrong pol g",i
!              enddo
!              do i=1,nq
!                 if(cdabs(p_q(1,i)-treeproc%quarks(i)%mom(1)+p_q(2,i)-treeproc%quarks(i)%mom(2) &
!                         +p_q(3,i)-treeproc%quarks(i)%mom(3)+p_q(4,i)-treeproc%quarks(i)%mom(4)).gt.1d-10) print *, "wrong momentum q",i
!                 if(cdabs(pol_q(1,i)-treeproc%quarks(i)%pol(1)+pol_q(2,i)-treeproc%quarks(i)%pol(2) &
!                         +pol_q(3,i)-treeproc%quarks(i)%pol(3)+pol_q(4,i)-treeproc%quarks(i)%pol(4)).gt.1d-10) print *, "wrong pol q",i
!              enddo
!
!           endif
!
!           if (fl(1) .eq. 'glu'.and.nq .eq. 4) then
!             if(ng1-TreeProc%NumGlu(1).ne.0.or.ng2-TreeProc%NumGlu(2).ne.0.or.ng3-TreeProc%NumGlu(3).ne.0 &
!                .or.ng4-TreeProc%NumGlu(4).ne.0.or.ng-1-ng1-ng2-ng3-ng4-TreeProc%NumGlu(5).ne.0) print *,"NumGlu check"
!           endif
! !DEC$ ENDIF
!
!
!       END SUBROUTINE calc_ampl









!       SUBROUTINE ampl(TreeProc,Dv,Ds,Nj1,Nj2,POLI,q1,lab1,tg,tag_f,ia,ib,BPOLF,q2,lab2,mur)
!       use ModProcess
!       use ModMisc
!       use ModParameters
!       implicit none
!       include 'global_import'
!       integer, intent(in) ::  Dv,Ds,Nj1,Nj2,tg,ia,ib,tag_f
!       character, intent(in) ::  lab1*3,lab2*3
!       complex(8), intent(in),target :: POLI(8,16),BPOLF(8,16)
!       complex(8), intent(in) :: q1(5),q2(5)
!       complex(8), intent(out) :: mur(10,10)
!       integer :: i,npa,i1,i2,j,ip,j2,j1,comp,qqs
!       character :: Lab_ex1(ib-ia+3)*3
!       character :: fl(ib-ia+3)*3,fl_aux(ib-ia+3)*3
!       character :: fl1(ib-ia+3)*3
!       complex(8) :: res
!       complex(8) :: pol1(Ds),pol2(Ds)
!       complex(8) :: po(ib-ia+3,Ds),pe(ib-ia+3,Dv)
!       complex(8),target :: p(Dv,ib-ia+3)
!       complex(8) :: sp(Ds,ib-ia+3),res1
!       complex(8) :: Jv(Dv),vpol1(Ds),Jsp(Ds)
!       complex(8) :: Jsp1(Ds)
!       complex(8) :: sp_aux(Ds,ib-ia+3),p_aux(Dv,ib-ia+3)
!       complex(8) :: sp1(Ds,ib-ia+3),p1(Dv,ib-ia+3)
!       complex(8) :: vg(Dv), vs(Ds),pol_glu(Dv)
!       integer :: ipo
!       integer :: ng,nq
!       type(TreeProcess) :: TreeProc
!
!
!           i1=tg
!           npa = ib-ia+3
!
!
!            do i=1,ib-ia+1
!             i2=i1 + (i-1)
!             if (i2 <= npoint) then
!                i2=i2
!             else
!                i2=i2 - Npoint
!             endif
!             do j=1,4
!                po(i+1,j)=hel(i2,j)
!                pe(i+1,j)=mom(i2,j)
!             enddo
!
!             do j=5,Ds
!                po(i+1,j)=czero
!             enddo
!             do j=5,Dv
!                pe(i+1,j)=czero
!             enddo
!             Lab_ex1(i+1)=Lab_ex(i2)
!           enddo
!
!           ng=0
!           nq=0
!
!           fl(1) =lab1
!           fl(npa)=lab2
!
!           do i=2,npa-1
!             fl(i)=Lab_ex1(i)
!           enddo
!
!           do i=1,npa
!             if (fl(i) == 'glu') ng=ng+1
!             if ((fl(i) == 'top').or.(fl(i) == 'bot').or.(fl(i) =='str')) nq=nq+1
!           enddo
!
!           if(Dv.eq.4) then
!             do j=1,4
!                p(j,1)  =q1(j)
!                p(j,npa)=-q2(j)
!             enddo
!           else
!              do j=1,5
!                p(j,1)  =q1(j)
!                p(j,npa)=-q2(j)
!              enddo
!              do j=6,Dv
!                p(j,1)  =czero
!                p(j,npa)=czero
!              enddo
!            endif
!
!           do j=1,Dv
!           do i=2,npa-1
!             p(j,i)=pe(i,j)
!           enddo
!           enddo
!
!           do j=1,Ds
!           do i=2,npa-1
!             sp(j,i)=po(i,j)
!           enddo
!           enddo
!
! !DEC$ IF (_DebugCheckMyImpl1==1)
!            if( TreeProc%PartType(TreeProc%NumPart).eq.Glu_ ) then
!                 if( abs(TreeProc%Gluons(TreeProc%NumGlu(0))%Mom(1) - p(1,npa) ).gt.1d-10) print *,"error in last mom"
!                 if( abs(TreeProc%Gluons(TreeProc%NumGlu(0))%Mom(2) - p(2,npa) ).gt.1d-10) print *,"error in last mom"
!                 if( abs(TreeProc%Gluons(TreeProc%NumGlu(0))%Mom(3) - p(3,npa) ).gt.1d-10) print *,"error in last mom"
!                 if( abs(TreeProc%Gluons(TreeProc%NumGlu(0))%Mom(4) - p(4,npa) ).gt.1d-10) print *,"error in last mom"
!            elseif( IsAQuark(TreeProc%PartType(TreeProc%NumPart)) ) then
!                 if( abs(TreeProc%Quarks(TreeProc%NumQua)%Mom(1) - p(1,npa) ).gt.1d-10) print *,"error in last mom"
!                 if( abs(TreeProc%Quarks(TreeProc%NumQua)%Mom(2) - p(2,npa) ).gt.1d-10) print *,"error in last mom"
!                 if( abs(TreeProc%Quarks(TreeProc%NumQua)%Mom(3) - p(3,npa) ).gt.1d-10) print *,"error in last mom"
!                 if( abs(TreeProc%Quarks(TreeProc%NumQua)%Mom(4) - p(4,npa) ).gt.1d-10) print *,"error in last mom"
!             endif
! !DEC$ ENDIF
!
!
!             do j2=1,Nj2
!                 do j=1,Ds
!                   sp(j,npa)=BPOLF(j2,j)
!                 enddo
!
!                 if( TreeProc%PartType(TreeProc%NumPart).eq.Glu_ ) then
!                    TreeProc%Gluons(TreeProc%NumGlu(0))%Pol => BPOLF(j2,1:Ds)
!                 elseif( IsAQuark(TreeProc%PartType(TreeProc%NumPart)) ) then
!                    TreeProc%Quarks(TreeProc%NumQua)%Pol => BPOLF(j2,1:Ds)
!                 endif
!
!                 sp(:,1) = czero
!                 call calc_ampl(TreeProc,Dv,Ds,ng,nq,0,sp,p,fl,tag_f,vg,vs)
!                 do j1=1,Nj1
!                   do j=1,Ds
!                      sp(j,1)=POLI(j1,j)     ! first particle is loop particle (off-shell in rec.rel.)
!                   enddo
!                   if (fl(1).eq.'glu') then
!                      pol_glu= sp(1:Dv,1)
!                      res = sc_(pol_glu,vg)
!                   endif
!                   if ((fl(1).eq.'top').or.(fl(1).eq.'bot').or.(fl(1).eq.'str')) then
!                      res = psp1_(vs,sp(:,1))
!                   endif
!                   mur(j1,j2)=res
!                 enddo
!             enddo
!
!         return
!         END SUBROUTINE ampl













END MODULE ModAmplitudes

