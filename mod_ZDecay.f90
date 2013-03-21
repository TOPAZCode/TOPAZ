MODULE ModZDecay
implicit none


private
public :: ZDecay
public :: ZGamPolVec,ZGamQcoupl,ZGamLcoupl

contains 


SUBROUTINE ZDecay(ZBoson,Topol,MomDK)
use ModParameters
use ModMisc
use ModProcess
implicit none
type(Particle) :: ZBoson
integer :: Topol
integer :: LepHel
real(8) :: MomDK(1:4,1:2),PropPhoton
complex(8) :: propZ
real(8) :: zeros(1:4)
complex(8) :: couplZFF_right,couplZFF_left,couplZFF


!DEC$ IF(_CheckMomenta .EQ.1)
   zeros(:) = 0d0
   zeros(1:4) = dble(ZBoson%Mom(1:4))
   zeros(1:4) = zeros(1:4) - MomDK(1:4,1) - MomDK(1:4,2)
   if( ZDecays.lt.10 .and. any(abs(zeros(1:4)/dble(ZBoson%Mom(1))).gt.1d-8) ) then! only check if Z boson is on-shell
      print *, "ERROR: energy-momentum violation in SUBROUTINE ZDecay(): ",zeros(1:4)
      print *, "momentum dump:"
      print *, MomDK(:,:)
   endif

   if( ZDecays.lt.10 ) zeros(1) = dble(ZBoson%Mom(1:4).dot.ZBoson%Mom(1:4)) - m_Z**2! only check if Z boson is on-shell
   zeros(2)=  MomDK(1:4,1).dot.MomDK(1:4,1)
   zeros(3)=  MomDK(1:4,2).dot.MomDK(1:4,2)
   if( any(abs(zeros(1:3)/dble(ZBoson%Mom(1))**2).gt.1d-5) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE ZDecay(): ",zeros(1:3)
      print *, "momentum dump:"
      print *, MomDK(:,:)
   endif
!DEC$ ENDIF


   if( ZDecays.lt.10 ) then! on-shell Z boson
       PropZ = (1d0,0d0)/dsqrt(2d0*Ga_Zexp*m_Z)
   elseif( ZDecays.gt.10 ) then! off-shell Z boson
       PropZ = (1d0,0d0)/( sc_(ZBoson%Mom(1:4),ZBoson%Mom(1:4))-m_Z**2 + ci*Ga_Zexp*m_Z )
   endif
   PropZ = PropZ * sc_(ZBoson%Mom(1:4),ZBoson%Mom(1:4))! this factor will cancel against a 1/kV^2 term in the decay matrix element

   PropPhoton = (1d0,0d0)! the photon propagator is included as 1/kV^2 term in the decay matrix element

    if( ZDecays.eq.1 .or. ZDecays.eq.11 ) then 
        couplZFF_right = couplZEE_right 
        couplZFF_left  = couplZEE_left
    elseif( ZDecays.eq.2 .or. ZDecays.eq.12 ) then 
        couplZFF_right = couplZNN_right
        couplZFF_left  = couplZNN_left
    else
        call Error("ZDecay not yet implemented",ZDecays)
    endif
    LepHel=ZBoson%Helicity
    
    if( LepHel.eq.+1 ) couplZFF = couplZFF_right 
    if( LepHel.eq.-1 ) couplZFF = couplZFF_left  


   if( ZDecays.lt.10  ) then  ! Z is on-shell
      couplZTT_left_dyn  = couplZTT_left *PropZ *couplZFF
      couplZTT_right_dyn = couplZTT_right*PropZ *couplZFF
   elseif( ZDecays.gt.10 ) then  ! Z is off-shell
      couplZTT_left_dyn  = couplZTT_left *PropZ *couplZFF + Q_Top*PropPhoton*Q_el
      couplZTT_right_dyn = couplZTT_right*PropZ *couplZFF + Q_Top*PropPhoton*Q_el
   endif


    if( Topol.eq.DK_LO ) then 
        if( abs(ZBoson%Helicity).ne.1 ) call Error("invalid Z boson polarization in ZDecay",ZBoson%Helicity)

        ZBoson%Pol(1:4)  = ZGamPolVec(dcmplx(MomDK(1:4,1)),dcmplx(MomDK(1:4,2)),ZBoson%Helicity)
        ZBoson%Pol(5:16) = 0d0
        
    else
        call Error("Wrong topology in ZDecay",Topol)
    endif


END SUBROUTINE




  function ZGamPolVec(pl,pa,LepHel)
! Polarization vector for Z/gamma decay to leptons
! pl, pa are momenta of lepton and antilepton
! LepType = 1 for charged leptons, 2 for neutrinos; 
! LepHel is lepton helicity
    use ModParameters
    use ModMisc
    implicit none
    complex(8)  :: pl(4), pa(4)
    integer     :: LepHel
    complex(8)  :: ZGamPolVec(4)
    complex(8)  :: pZ(4),leppol(4),aleppol(4)
    real(8) :: couplZFF_right,couplZFF_left
    
    pZ=pl+pa
    
    call ubarSpi_Weyl(pl,LepHel,leppol)
    call vSpi_Weyl(pa,-LepHel,aleppol)
    ZGamPolVec(1:4)=vbqq_Weyl(4,leppol,aleppol)*sqrt2

! NB: note factor 1/pZ^2 in here:
    ZGamPolVec=ZGamPolVec/(sc_(pZ,pZ))
    ZGamPolVec=ZGamPolVec*dsqrt(alpha4Pi)

! RR -- these are all set in the ZDecay routine
!    if( ZDecays.eq.1 .or. ZDecays.eq.11 ) then 
!        couplZFF_right = couplZEE_right 
!        couplZFF_left  = couplZEE_left
!    elseif( ZDecays.eq.2 .or. ZDecays.eq.12 ) then 
!        couplZFF_right = couplZNN_right
!        couplZFF_left  = couplZNN_left
!    else
!        call Error("ZDecay not yet implemented",ZDecays)
!    endif
!    
!    if( -LepHel.eq.+1 ) ZGamPolVec = ZGamPolVec * couplZFF_right * dsqrt(alpha4Pi)
!    if( -LepHel.eq.-1 ) ZGamPolVec = ZGamPolVec * couplZFF_left  * dsqrt(alpha4Pi)


  end function ZGamPolVec




  subroutine ZGamQcoupl(QType,QHel,couplZQQ,couplGQQ)
! QType = (u,d,c,s,t,b), QHel is quark helicity
! Can modify it to allow BSM ttZ couplings, F_V-F_A for LH top, F_V+F_A for RH
    use ModParameters
    integer, intent(in)     :: QType, QHel
    real(8), intent(out)     :: couplZQQ, couplGQQ

! Decide on quark charge and isospin    
    if ( mod(QType,2) .eq. 1) then      ! up-type quark
       if (QHel .eq. -1) then
          couplZQQ=couplZUU_left
       elseif (QHel .eq. 1) then
          couplZQQ=couplZUU_right
       else
          stop "Error in ZDecay: QHel undefined"
       endif
       couplGQQ=Q_up
    elseif ( mod(QType,2) .eq. 0) then      ! down-type quark
       if (QHel .eq. -1) then
          couplZQQ=couplZDD_left
       elseif (QHel .eq. 1) then
          couplZQQ=couplZDD_right
       else
          stop "Error in ZDecay: QHel undefined"
       endif    
       couplGQQ=Q_dn
    else
       stop "Error in ZDecay: QType undefined"
    endif

  end subroutine ZGamQcoupl
    

  subroutine ZGamLCoupl(LepType,LepHel,couplZLL,couplGLL)
! LepType = 1 for charged leptons, 2 for neutrinos; 
! LepHel is lepton helicity
    use ModParameters
    use ModMisc
    implicit none
    integer, intent(in)     :: LepType,LepHel
    real(8),intent(out)     :: couplZLL,couplGLL
    

 ! Decide on lepton couplings
    if (LepType .eq. 1) then        ! charged lepton
       if (LepHel .eq. -1) then
          couplZLL=couplZEE_left
       elseif (LepHel .eq. 1) then
          couplZLL=couplZEE_right
       else
          stop "Error in ZDecay: LepHel undefined"
       endif
       couplGLL=Q_el
    elseif (LepType .eq. 2) then        ! neutrino
       if (LepHel .eq. -1) then
          couplZLL=couplZNN_left
       elseif (LepHel .eq. 1) then
          couplZLL=couplZNN_right
       else
          stop "Error in ZDecay: LepHel undefined"
       endif
       couplGLL=Q_nu       ! 0
    else
       stop "Error in ZDecay: LepType undefined"
    endif

  end subroutine ZGamLCoupl

end MODULE ModZDecay
