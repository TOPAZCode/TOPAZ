module modIntDipoles_ZprimeTTB
use ModIntDipoles
implicit none

integer, parameter, private :: dp = selected_real_kind(15)

double precision, private :: MomDK(0:3,1:6)

contains

  Subroutine IntDip_qg_ZprimeInt_ttb(p,z,IDip)
    use ModParameters
    use ModMisc
    implicit none
    real(8) :: p(4,4),pflip(4,4)
    real(8) :: IDip(1:3),APsoft,APfini,APplus,z,sHat
    real(8) :: dipsoft,dipfini,dipplus,epcorr!, epinv, epinv2
    real(8) :: CF 
    integer :: n,emi

    print *, 'IntDip_qg_ZprimeInt_ttb not yet implemented'

    stop

  END Subroutine IntDip_qg_ZprimeInt_ttb

  Subroutine IntDip_qbg_ZprimeInt_ttb(p,z,IDip)
    use ModParameters
    use ModMisc
    implicit none
    real(8) :: p(4,4),pflip(4,4)
    real(8) :: IDip(1:3),APsoft,APfini,APplus,z,sHat
    real(8) :: dipsoft,dipfini,dipplus,epcorr!, epinv, epinv2
    real(8) :: CF 
    integer :: n,emi

    print *, 'IntDip_qbg_ZprimeInt_ttb not yet implemented'

    stop

  END Subroutine IntDip_qbg_ZprimeInt_ttb


  SUBROUTINE IntDip_qqb_ZprimeInt_ttb(p,z,IDip)
    use ModParameters
    use ModMisc
    implicit none
    real(8) :: p(4,4),pflip(4,4)
    real(8) :: IDip(1:3),APsoft,APfini,APplus,z,sHat
    real(8) :: dipsoft,dipfini,dipplus,epcorr!, epinv, epinv2
    real(8) :: CF 
    integer :: n,emi

!    epinv2=0d0
!    epinv =0d0

   ! all-outgoing
    pflip(:,1) = -p(:,1)
    pflip(:,2) = -p(:,2)
    pflip(:,3) = p(:,3)
    pflip(:,4) = p(:,4)

    IDip(1:3) = 0d0
    do n=1,8
       if(n.eq.1) then
          dipsoft =-if_qq(0d0,m_Top,pflip,1,3,z,1)
          dipfini =-if_qq(0d0,m_Top,pflip,1,3,z,2) 
          dipplus =-if_qq(0d0,m_Top,pflip,1,3,z,3) 
          emi = 1
       elseif(n.eq.2) then
          dipsoft =if_qq(0d0,m_Top,pflip,2,3,z,1)
          dipfini =if_qq(0d0,m_Top,pflip,2,3,z,2)
          dipplus =if_qq(0d0,m_Top,pflip,2,3,z,3)
          emi = 2
       elseif(n.eq.3) then
          dipsoft =if_qq(0d0,m_Top,pflip,1,4,z,1)
          dipfini =if_qq(0d0,m_Top,pflip,1,4,z,2) 
          dipplus =if_qq(0d0,m_Top,pflip,1,4,z,3) 
          emi = 1
       elseif(n.eq.4) then
          dipsoft =-if_qq(0d0,m_Top,pflip,2,4,z,1)
          dipfini =-if_qq(0d0,m_Top,pflip,2,4,z,2)
          dipplus =-if_qq(0d0,m_Top,pflip,2,4,z,3)
          emi = 2
       elseif(n.eq.5) then
          dipsoft =-fi_qq(m_Top,0d0,pflip,3,1,z,1)
          dipfini =-fi_qq(m_Top,0d0,pflip,3,1,z,2)
          dipplus =-fi_qq(m_Top,0d0,pflip,3,1,z,3)
          emi = 1
       elseif(n.eq.6) then
          dipsoft =fi_qq(m_Top,0d0,pflip,3,2,z,1)
          dipfini =fi_qq(m_Top,0d0,pflip,3,2,z,2)
          dipplus =fi_qq(m_Top,0d0,pflip,3,2,z,3)
          emi = 2
       elseif(n.eq.7) then
          dipsoft =fi_qq(m_Top,0d0,pflip,4,1,z,1)
          dipfini =fi_qq(m_Top,0d0,pflip,4,1,z,2)
          dipplus =fi_qq(m_Top,0d0,pflip,4,1,z,3)
          emi = 1
       elseif(n.eq.8) then
          dipsoft =-fi_qq(m_Top,0d0,pflip,4,2,z,1)
          dipfini =-fi_qq(m_Top,0d0,pflip,4,2,z,2)
          dipplus =-fi_qq(m_Top,0d0,pflip,4,2,z,3)
          emi = 2
       endif

!       print *, 'dipoles set to zero for 1/ep'
!       dipplus = 0d0; print *, 'int dip ep check'
!       dipfini = 0d0; print *, 'int dip ep check'

       if(emi.eq.1) then
          IDip(1) = IDip(1) + (dipsoft-dipplus)
          IDip(2) = IDip(2) + (dipfini+dipplus)
       elseif(emi.eq.2) then
          IDip(1) = IDip(1) + (dipsoft-dipplus)
          IDip(3) = IDip(3) + (dipfini+dipplus)
       elseif(emi.eq.3) then
          IDip(1) = IDip(1) + (dipsoft-dipplus)
          IDip(2) = IDip(2) + (dipfini+dipplus)*0.5d0
          IDip(3) = IDip(3) + (dipfini+dipplus)*0.5d0
       endif

    enddo

    ! !        epcorr=epinv+2d0*dlog(renscale/facscale)
 !   epcorr=epinv

!    AP not needed for interference
!    APsoft= 3d0/2d0*CF     * epcorr
!    APfini= (-1d0-z)*CF    * epcorr
!    APplus= 2d0*CF/(1d0-z) * epcorr

    IDip(1) = IDip(1) !+ (APsoft - APplus)*2d0
    IDip(2) = IDip(2) !+ (APfini + APplus)
    IDip(3) = IDip(3) !+ (APfini + APplus)

!     print *, "AP",(APsoft - APplus),(APfini + APplus)
!     print *, "sum",IDip(2:3)
!     pause

  END SUBROUTINE IntDip_qqb_ZprimeInt_ttb




  SUBROUTINE IntDip_qqb_Zprime_ttb(p,z,IDip)
    use ModParameters
    use ModMisc
    implicit none
    real(8) :: p(4,4)
    real(8) :: IDip(1:3),APsoft,APfini,APplus,z,sHat
    real(8) :: dipsoft,dipfini,dipplus,epcorr!, epinv, epinv2
    real(8) :: CF 
    integer :: n,emi


    CF = 4d0/3d0

!    epinv2=0d0
!    epinv =0d0

    IDip(1:3) = 0d0
    do n=1,4
       if(n.eq.1) then
          dipsoft =ii_qq(0d0,0d0,p,1,2,z,1) * CF
          dipfini =ii_qq(0d0,0d0,p,1,2,z,2) * CF
          dipplus =ii_qq(0d0,0d0,p,1,2,z,3) * CF
          emi = 1
       elseif(n.eq.2) then
          dipsoft =ii_qq(0d0,0d0,p,2,1,z,1) * CF
          dipfini =ii_qq(0d0,0d0,p,2,1,z,2) * CF
          dipplus =ii_qq(0d0,0d0,p,2,1,z,3) * CF
          emi = 2
       elseif(n.eq.3) then
          dipsoft =ff_qq(m_Top,m_Top,p,3,4,z,1) * CF
          dipfini =ff_qq(m_Top,m_Top,p,3,4,z,2) * CF
          dipplus =ff_qq(m_Top,m_Top,p,3,4,z,3) * CF
          emi = 3
       elseif(n.eq.4) then
          dipsoft =ff_qq(m_Top,m_Top,p,4,3,z,1) * CF
          dipfini =ff_qq(m_Top,m_Top,p,4,3,z,2) * CF
          dipplus =ff_qq(m_Top,m_Top,p,4,3,z,3) * CF
          emi = 3
       endif

       ! this is for check against virtual amplitude
!       dipplus = 0d0; print *, 'dipoles check'
!       dipfini = 0d0; print *, 'dipoles check'

       if(emi.eq.1) then
          IDip(1) = IDip(1) + (dipsoft-dipplus)
          IDip(2) = IDip(2) + (dipfini+dipplus)
       elseif(emi.eq.2) then
          IDip(1) = IDip(1) + (dipsoft-dipplus)
          IDip(3) = IDip(3) + (dipfini+dipplus)
       elseif(emi.eq.3) then
          IDip(1) = IDip(1) + (dipsoft-dipplus)
          IDip(2) = IDip(2) + (dipfini+dipplus)*0.5d0
          IDip(3) = IDip(3) + (dipfini+dipplus)*0.5d0
       endif
    enddo

    
     !print *, epinv
     !print *, "IntDip",IDip(2:3)

    ! !        epcorr=epinv+2d0*dlog(renscale/facscale)
    epcorr=epinv

    !      this is for qqb-->Zpr-->ttb+g process
    APsoft= 3d0/2d0*CF     * epcorr
    APfini= (-1d0-z)*CF    * epcorr
    APplus= 2d0*CF/(1d0-z) * epcorr

    ! this is for check against virtual amplitude
!    APplus = 0d0; print *, 'dipoles check'
!    APfini = 0d0; print *, 'dipoles check'

    IDip(1) = IDip(1) + (APsoft - APplus)*2d0
    IDip(2) = IDip(2) + (APfini + APplus)
    IDip(3) = IDip(3) + (APfini + APplus)

     !print *, "AP",(APsoft - APplus),(APfini + APplus)
     !print *, "sum",IDip(2:3)
     !pause

  END SUBROUTINE IntDip_qqb_Zprime_ttb


  

  SUBROUTINE IntDip_gqb_Zprime_ttb(p,z,IDip)
    use ModParameters
    use ModMisc
    implicit none
    real(8) :: p(4,4)
    real(8) :: IDip(1:3),APsoft,APfini,APplus,z,sij,sHat
    real(8) :: dipsoft,dipfini,dipplus,epcorr
    integer :: n,emi
    real(8) :: TR

    TR = 1d0/2d0

!    epinv2=0d0
!    epinv =0d0

    IDip(1:3) = 0d0

    dipsoft =ii_gq(0d0,0d0,p,1,2,z,1) !* TR Already taken care of in the dipole subroutine
    dipfini =ii_gq(0d0,0d0,p,1,2,z,2) !* TR
    dipplus =ii_gq(0d0,0d0,p,1,2,z,3) !* TR
    emi = 1


    if(emi.eq.1) then
       IDip(1) = IDip(1) + (dipsoft-dipplus)
       IDip(2) = IDip(2) + (dipfini+dipplus)
    elseif(emi.eq.2) then! never called
       IDip(1) = IDip(1) + (dipsoft-dipplus)
       IDip(3) = IDip(3) + (dipfini+dipplus)
    elseif(emi.eq.3) then! never called
       IDip(1) = IDip(1) + (dipsoft-dipplus)
       IDip(2) = IDip(2) + (dipfini+dipplus)*0.5d0
       IDip(3) = IDip(3) + (dipfini+dipplus)*0.5d0
    endif


     !print *, epinv
     !print *, "IntDip",IDip(2:3)

    ! !        epcorr=epinv+2d0*dlog(renscale/facscale)
    epcorr=epinv

    !      this is for gqb-->Zpr-->ttb+qb process
    APsoft= 0d0
    APfini= TR*(z**2+(1d0-z)**2) * epcorr
    APplus= 0d0
    IDip(1) = IDip(1) + (APsoft - APplus)
    IDip(2) = IDip(2) + (APfini + APplus)
    
    ! print *, "AP",(APsoft - APplus),(APfini + APplus)
    ! print *, "sum",IDip(1:3)
    ! pause
    
  END SUBROUTINE IntDip_gqb_Zprime_ttb




  SUBROUTINE IntDip_qg_Zprime_ttb(p,z,IDip)
    use ModParameters
    use ModMisc
    implicit none
    real(8) :: p(4,4)
    real(8) :: IDip(1:3),APsoft,APfini,APplus,z,sij,sHat
    real(8) :: dipsoft,dipfini,dipplus,epcorr
    integer :: n,emi
    real(8) :: TR

    TR = 1d0/2d0

    !epinv2=0d0
    !epinv =0d0

    IDip(1:3) = 0d0

    dipsoft =ii_gq(0d0,0d0,p,2,1,z,1) !* TR Already taken care of in the dipole subroutine
    dipfini =ii_gq(0d0,0d0,p,2,1,z,2) !* TR
    dipplus =ii_gq(0d0,0d0,p,2,1,z,3) !* TR
    emi = 2


    if(emi.eq.1) then! never called
       IDip(1) = IDip(1) + (dipsoft-dipplus)
       IDip(2) = IDip(2) + (dipfini+dipplus)
    elseif(emi.eq.2) then
       IDip(1) = IDip(1) + (dipsoft-dipplus)
       IDip(3) = IDip(3) + (dipfini+dipplus)
    elseif(emi.eq.3) then! never called
       IDip(1) = IDip(1) + (dipsoft-dipplus)
       !         IDip(2) = IDip(2) + (dipfini+dipplus)*0.5d0
       !         IDip(3) = IDip(3) + (dipfini+dipplus)*0.5d0
    endif
    

     !print *, epinv
     !print *, "IntDip",IDip(2:3)
    
    ! !        epcorr=epinv+2d0*dlog(renscale/facscale)
    epcorr=epinv

    !      this is for gqb-->Zpr-->ttb+qb process
    APsoft= 0d0
    APfini= TR*(z**2+(1d0-z)**2) * epcorr
    APplus= 0d0
    IDip(1) = IDip(1) + (APsoft - APplus)
    IDip(3) = IDip(3) + (APfini + APplus)
    
     !print *, "AP",(APsoft - APplus),(APfini + APplus)
     !print *, "sum",IDip(1:3)
     !pause
    


  END SUBROUTINE IntDip_qg_Zprime_ttb







END module modIntDipoles_ZprimeTTB
