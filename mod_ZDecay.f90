MODULE ModZDecay
implicit none

contains 

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
    
    pZ=pl+pa
    
    call ubarSpi_Weyl(pl,LepHel,leppol)
    call vSpi_Weyl(pa,-LepHel,aleppol)
    ZGamPolVec(1:4)=vbqq_Weyl(4,leppol,aleppol)

! NB: note factor 1/pZ^2 in here:
    ZGamPolVec=ZGamPolVec/(sc_(pZ,pZ))


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
