MODULE ModZDecay
implicit none

contains 

  function ZGamPolVec(pl,pa,QType,LepType,QHel,LepHel)
! Polarization vector for Z/gamma decay to leptons
! pl, pa are momenta of lepton and antilepton
! QType = (u,d,c,s,t,b), QHel is quark helicity
! LepType = 1 for charged leptons, 2 for neutrinos; 
! LepHel is lepton helicity
! Can modify it to allow BSM ttZ couplings, F_V-F_A for LH top, F_V+F_A for RH
    use ModParameters
    use ModMisc
    implicit none
    complex(8)  :: pl(4), pa(4)
    integer     :: QType, LepType,QHel, LepHel
    complex(8)  :: ZGamPolVec(4)
    complex(8)  :: pZ(4),leppol(4),aleppol(4),propZ,propG
    real(8)     :: couplZQQ, couplGQQ,couplZLL,couplGLL
    
    pZ=pl+pa
    propZ=1d0/(sc_(pZ,pZ)-m_Z**2+ci*Ga_ZExp*m_Z)
    propG=1d0/sc_(pZ,pZ)

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
    
    call ubarSpi_Weyl(pl,LepHel,leppol)
    call vSpi_Weyl(pa,-LepHel,aleppol)
    ZGamPolVec(1:4)=vbqq_Weyl(4,aleppol,leppol)

! now add couplings
    ZGamPolVec=ZGamPolVec*(couplZQQ*couplZLL*propZ+couplGQQ*couplGLL*propG)

  end function ZGamPolVec

end MODULE ModZDecay
