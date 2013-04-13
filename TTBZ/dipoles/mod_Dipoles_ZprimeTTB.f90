#DEFINE _CHECK_DIPOLE_MOMMAP 0
#DEFINE _APPLY_CUTS 0

  MODULE ModDipoles_Zprime
    use ModTopdecay
    implicit none


    double precision, private, parameter :: NCol=3d0
!--F check TR
    double precision, private, parameter :: TR=1d0/2d0
    double precision, private, parameter :: CA=2d0*TR*NCol
    double precision, private, parameter :: CF=TR*(NCol**2-1d0)/NCol

    type :: Dipole
       double precision DipoleValue
       double precision MomTd(0:3,1:4)
       double precision Mass2Td(1:4)
       double precision alphaCut
    end type Dipole

    type(Dipole),private :: TheDipoles(1:4)

    double precision, private :: yRndDK(1:8),xFrac,r_sc
    integer, private :: NBin(1:20)
    
    logical, parameter :: invert_alphaCut = .false.

  contains



    
    SUBROUTINE Dipoles_qqb_Zprime_ttb(nDipole,MomExt,MomExtTd,Dipole)! global norm:   4d0*Pi*alpha_s
      use ModParameters
      use ModKinematics
      use ModMisc
      implicit none
      integer :: nDipole,a,i,b,j,k
      real(8) :: MomExt(1:4,1:5),MomExtTd(1:4,1:4),Q(1:4),QTd(1:4),KSum(1:4)
      real(8) :: sab,sai,sbi,sij,sik,skj,x,v,y,yp,z,Q2,mu2,mu,MomFac1,MomFac2,MomFac3
      real(8) :: Dipole
      
      
      if(nDipole.eq.1) then
         a=1; i=3; b=2!   initial-initial
         
         sab = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,b))
         sai = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,i))
         sbi = 2d0*(MomExt(1:4,b).dot.MomExt(1:4,i))
         x = 1d0 - (sai+sbi)/sab
         v = sai/sab
         if( .not. invert_alphacut ) then
            if( alpha_ii.lt.v ) then
                Dipole = (0d0,0d0)
                return
            endif
         else
            if( alpha_ii.gt.v ) then
                Dipole = (0d0,0d0)
                return
            endif
         endif
         
         MomExtTd(1:4,a) = x*MomExt(1:4,a)
         MomExtTd(1:4,b) = MomExt(1:4,b)
         
         Q(1:4)   = MomExt(1:4,a)+MomExt(1:4,b)-MomExt(1:4,i)
         QTd(1:4) = MomExtTd(1:4,a)+MomExtTd(1:4,b)
         KSum(1:4) = Q(1:4)+QTd(1:4)
         MomExtTd(1:4,3) = MomExt(1:4,4) - 2d0*(MomExt(1:4,4).dot.KSum(1:4))/(KSum(1:4).dot.KSum(1:4))*KSum(1:4) + 2d0*(MomExt(1:4,4).dot.Q(1:4))/(Q(1:4).dot.Q(1:4))*QTd(1:4)
         MomExtTd(1:4,4) = MomExt(1:4,5) - 2d0*(MomExt(1:4,5).dot.KSum(1:4))/(KSum(1:4).dot.KSum(1:4))*KSum(1:4) + 2d0*(MomExt(1:4,5).dot.Q(1:4))/(Q(1:4).dot.Q(1:4))*QTd(1:4)
         
         Dipole = -1d0/sai/x * (2d0/(1d0-x)-1d0-x)

      elseif(nDipole.eq.2) then
         a=2; i=3; b=1!   initial-initial
         
         sab = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,b))
         sai = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,i))
         sbi = 2d0*(MomExt(1:4,b).dot.MomExt(1:4,i))
         x = 1d0 - (sai+sbi)/sab
         v = sai/sab
         if( .not. invert_alphacut ) then
              if( alpha_ii.lt.v ) then
                  Dipole = (0d0,0d0)
                  return
              endif
         else
              if( alpha_ii.gt.v ) then
                  Dipole = (0d0,0d0)
                  return
              endif
         endif

         MomExtTd(1:4,a) = x*MomExt(1:4,a)
         MomExtTd(1:4,b) = MomExt(1:4,b)
         
         Q(1:4)   = MomExt(1:4,a)+MomExt(1:4,b)-MomExt(1:4,i)
         QTd(1:4) = MomExtTd(1:4,a)+MomExtTd(1:4,b)
         KSum(1:4) = Q(1:4)+QTd(1:4)
         MomExtTd(1:4,3) = MomExt(1:4,4) - 2d0*(MomExt(1:4,4).dot.KSum(1:4))/(KSum(1:4).dot.KSum(1:4))*KSum(1:4) + 2d0*(MomExt(1:4,4).dot.Q(1:4))/(Q(1:4).dot.Q(1:4))*QTd(1:4)
         MomExtTd(1:4,4) = MomExt(1:4,5) - 2d0*(MomExt(1:4,5).dot.KSum(1:4))/(KSum(1:4).dot.KSum(1:4))*KSum(1:4) + 2d0*(MomExt(1:4,5).dot.Q(1:4))/(Q(1:4).dot.Q(1:4))*QTd(1:4)
         
         Dipole = -1d0/sai/x * (2d0/(1d0-x)-1d0-x)




      elseif(nDipole.eq.3) then
         i=4; j=3; k=5! final-final
         
         sij = 2d0*(MomExt(1:4,i)).dot.(MomExt(1:4,j))
         sik = 2d0*(MomExt(1:4,i)).dot.(MomExt(1:4,k))
         skj = 2d0*(MomExt(1:4,k)).dot.(MomExt(1:4,j))
         Q2 = 2d0*m_Top**2+ sij+sik+skj
         mu2 = m_Top**2/dabs(Q2)
         mu = dsqrt(mu2)
         y = sij/(sij+skj+sik)
         z = sik/(skj+sik)
         v = (1d0-y)*dsqrt((1d0-4d0*mu2)/( (2d0*mu2+(1d0-2d0*mu2)*(1d0-y))**2-4d0*mu2))
         yp = 1d0-2d0*mu*(1d0-mu)/(1d0-2d0*mu2)
         if( .not. invert_alphacut ) then
              if( y.gt.alpha_ff*yp) then
                  Dipole = (0d0,0d0)
                  return
              endif
         else
              if( y.lt.alpha_ff*yp) then
                  Dipole = (0d0,0d0)
                  return
              endif
         endif

         MomFac1 = dsqrt(((sij+sik+skj)**2-4d0*m_Top**4)/((sik+skj)**2-4d0*m_Top**2*(m_Top**2+sij)))
         MomFac2 = (0.5d0*(sik+skj)+m_Top**2)/Q2
         MomFac3 = MomFac2 + 0.5d0*sij/Q2
         
         Q(1:4) = MomExt(1:4,i)+MomExt(1:4,j)+MomExt(1:4,k)
         MomExtTd(1:4,k-1) = MomFac1*(MomExt(1:4,k)-MomFac2*Q(1:4)) + MomFac3*Q(1:4)
         MomExtTd(1:4,i-1) = Q(1:4) - MomExtTd(1:4,k-1)
         MomExtTd(1:4,1:2) = MomExt(1:4,1:2)
         
         Dipole = -1d0/sij * (2d0/(1d0-z+y*z)-v*(1d0+z+2d0*m_Top**2/sij))


      elseif(nDipole.eq.4) then
         i=5; j=3; k=4! final-final

         sij = 2d0*(MomExt(1:4,i)).dot.(MomExt(1:4,j))
         sik = 2d0*(MomExt(1:4,i)).dot.(MomExt(1:4,k))
         skj = 2d0*(MomExt(1:4,k)).dot.(MomExt(1:4,j))
         Q2 = 2d0*m_Top**2 + sij+sik+skj
         mu2 = m_Top**2/dabs(Q2)
         mu = dsqrt(mu2)
         y = sij/(sij+skj+sik)
         z = sik/(skj+sik)
         v = (1d0-y)*dsqrt((1d0-4d0*mu2)/( (2d0*mu2+(1d0-2d0*mu2)*(1d0-y))**2-4d0*mu2))
         yp = 1d0-2d0*mu*(1d0-mu)/(1d0-2d0*mu2)
         if( .not. invert_alphacut ) then
              if( y.gt.alpha_ff*yp) then
                  Dipole = (0d0,0d0)
                  return
              endif
         else
              if( y.lt.alpha_ff*yp) then
                  Dipole = (0d0,0d0)
                  return
              endif
         endif
         MomFac1 = dsqrt(((sij+sik+skj)**2-4d0*m_Top**4)/((sik+skj)**2-4d0*m_Top**2*(m_Top**2+sij)))
         MomFac2 = (0.5d0*(sik+skj)+m_Top**2)/Q2
         MomFac3 = MomFac2 + 0.5d0*sij/Q2
         
         Q(1:4) = MomExt(1:4,i)+MomExt(1:4,j)+MomExt(1:4,k)
         MomExtTd(1:4,k-1) = MomFac1*(MomExt(1:4,k)-MomFac2*Q(1:4)) + MomFac3*Q(1:4)
         MomExtTd(1:4,i-1) = Q(1:4) - MomExtTd(1:4,k-1)
         MomExtTd(1:4,1:2) = MomExt(1:4,1:2)
         

         Dipole = -1d0/sij * (2d0/(1d0-z+y*z)-v*(1d0+z+2d0*m_Top**2/sij))

      endif

    END SUBROUTINE Dipoles_qqb_Zprime_ttb





    SUBROUTINE Dipoles_gqb_Zprime_ttbqb(nDipole,MomExt,MomExtTd,Dipole)! global norm:   4d0*Pi*alpha_s
      use ModParameters
      use ModKinematics
      use ModMisc
      implicit none
      integer :: nDipole,a,i,b,j,k
      real(8) :: MomExt(1:4,1:5),MomExtTd(1:4,1:4),Q(1:4),QTd(1:4),KSum(1:4)
      real(8) :: sab,sai,sbi,sij,sik,skj,x,v,y,yp,z,Q2,mu2,mu,MomFac1,MomFac2,MomFac3
      real(8) :: Dipole
      
      a=1; i=3; b=2 ! initial-initial
      
      sab = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,b))
      sai = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,i))
      sbi = 2d0*(MomExt(1:4,b).dot.MomExt(1:4,i))
      x = 1d0 - (sai+sbi)/sab
      v = sai/sab
      
      if( .not. invert_alphacut ) then
            if( alpha_ii.lt.v ) then
              Dipole = (0d0,0d0)
              return
            endif
      else
            if( alpha_ii.gt.v ) then
              Dipole = (0d0,0d0)
              return
            endif
      endif

      MomExtTd(1:4,a) = x*MomExt(1:4,a)
      MomExtTd(1:4,b) = MomExt(1:4,b)

      Q(1:4)   = MomExt(1:4,a)+MomExt(1:4,b)-MomExt(1:4,i)
      QTd(1:4) = MomExtTd(1:4,a)+MomExtTd(1:4,b)
      KSum(1:4) = Q(1:4)+QTd(1:4)
      MomExtTd(1:4,3) = MomExt(1:4,4) - 2d0*(MomExt(1:4,4).dot.KSum(1:4))/(KSum(1:4).dot.KSum(1:4))*KSum(1:4) &
           + 2d0*(MomExt(1:4,4).dot.Q(1:4))/(Q(1:4).dot.Q(1:4))*QTd(1:4)
      MomExtTd(1:4,4) = MomExt(1:4,5) - 2d0*(MomExt(1:4,5).dot.KSum(1:4))/(KSum(1:4).dot.KSum(1:4))*KSum(1:4) &
           + 2d0*(MomExt(1:4,5).dot.Q(1:4))/(Q(1:4).dot.Q(1:4))*QTd(1:4)
      Dipole = -1d0/sai/x * 2d0*TR * (1d0-2d0*x*(1d0-x))
            
    END SUBROUTINE Dipoles_gqb_Zprime_ttbqb



    SUBROUTINE Dipoles_qg_Zprime_ttbq(nDipole,MomExt,MomExtTd,Dipole)! global norm:   4d0*Pi*alpha_s
      use ModParameters
      use ModKinematics
      use ModMisc
      implicit none
      integer :: nDipole,a,i,b,j,k
      real(8) :: MomExt(1:4,1:5),MomExtTd(1:4,1:4),Q(1:4),QTd(1:4),KSum(1:4)
      real(8) :: sab,sai,sbi,sij,sik,skj,x,v,y,yp,z,Q2,mu2,mu,MomFac1,MomFac2,MomFac3
      real(8) :: Dipole
      
      a=2; i=3; b=1!   initial-initial

      sab = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,b))
      sai = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,i))
      sbi = 2d0*(MomExt(1:4,b).dot.MomExt(1:4,i))
      x = 1d0 - (sai+sbi)/sab
      v = sai/sab
      if( .not. invert_alphacut ) then
            if( alpha_ii.lt.v ) then
              Dipole = (0d0,0d0)
              return
            endif
      else
            if( alpha_ii.gt.v ) then
              Dipole = (0d0,0d0)
              return
            endif
      endif

      MomExtTd(1:4,a) = x*MomExt(1:4,a)
      MomExtTd(1:4,b) = MomExt(1:4,b)
      
      Q(1:4)   = MomExt(1:4,a)+MomExt(1:4,b)-MomExt(1:4,i)
      QTd(1:4) = MomExtTd(1:4,a)+MomExtTd(1:4,b)
      KSum(1:4) = Q(1:4)+QTd(1:4)
      MomExtTd(1:4,3) = MomExt(1:4,4) - 2d0*(MomExt(1:4,4).dot.KSum(1:4))/(KSum(1:4).dot.KSum(1:4))*KSum(1:4) &
           + 2d0*(MomExt(1:4,4).dot.Q(1:4))/(Q(1:4).dot.Q(1:4))*QTd(1:4)
      MomExtTd(1:4,4) = MomExt(1:4,5) - 2d0*(MomExt(1:4,5).dot.KSum(1:4))/(KSum(1:4).dot.KSum(1:4))*KSum(1:4) &
           + 2d0*(MomExt(1:4,5).dot.Q(1:4))/(Q(1:4).dot.Q(1:4))*QTd(1:4)

      Dipole = -1d0/sai/x * 2d0*TR * (1d0-2d0*x*(1d0-x))

    END SUBROUTINE Dipoles_qg_Zprime_ttbq





    SUBROUTINE Dipoles_qqb_Zprime_interf(nDipole,MomExt,MomExtTd,Dipole)! global norm:   4d0*Pi*alpha_s
      use ModParameters
      use ModKinematics
      use ModMisc
      implicit none
      integer :: nDipole,a,i,b,j,k
      real(8) :: MomExt(1:4,1:5),MomExtTd(1:4,1:4),Q(1:4),QTd(1:4),KSum(1:4)
      real(8) :: sab,sai,saj,sbj,sij,x,u,v,y,yp,z,Q2,mu2,mu,MomFac1,MomFac2,MomFac3
      real(8) :: Dipole
      
      Dipole = 0d0
      
      if (nDipole.eq.1) then

         a=1; i=3; j=4; !initial-final, q tb
         
         sai = 2d0*(MomExt(1:4,a)).dot.(MomExt(1:4,i)) 
         saj = 2d0*(MomExt(1:4,a)).dot.(MomExt(1:4,j)) 
         sij = 2d0*(MomExt(1:4,i)).dot.(MomExt(1:4,j)) 
         
         x = 1d0-sij/(sai+saj)
         u = sai/(sai+saj)
         if( alpha_if.lt.u ) then
!         if( alpha_if.gt.u ) then !-- inverted alpha check
            Dipole = 0d0
            return
         endif

         Dipole = -1d0/sai/x
         
         MomExtTd(1:4,1) = x*MomExt(1:4,1) ! out q
         MomExtTd(1:4,2) = MomExt(1:4,2) ! out qb
         MomExtTd(1:4,3) = MomExt(1:4,4) + MomExt(1:4,3) - (1d0-x)*MomExt(1:4,1) ! tb
         MomExtTd(1:4,4) = MomExt(1:4,5) ! t

         !     splitting 1: A/Quark_I --> A/Quark_I + Gluon_F
         Dipole = Dipole * (2d0/(1d0-x+u)-1d0-x)

      elseif (nDipole.eq.2) then

         a=2; i=3; j=4; !initial-final, qb tb
         
         sai = 2d0*(MomExt(1:4,a)).dot.(MomExt(1:4,i)) 
         saj = 2d0*(MomExt(1:4,a)).dot.(MomExt(1:4,j)) 
         sij = 2d0*(MomExt(1:4,i)).dot.(MomExt(1:4,j)) 
         
         x = 1d0-sij/(sai+saj)
         u = sai/(sai+saj)
         if( alpha_if.lt.u ) then
!         if( alpha_if.gt.u ) then !-- inverted alpha check
            Dipole = 0d0
            return
         endif
         Dipole = -1d0/sai/x
         
         MomExtTd(1:4,1) = MomExt(1:4,1) ! out q
         MomExtTd(1:4,2) = x*MomExt(1:4,2) ! out qb
         MomExtTd(1:4,3) = MomExt(1:4,4) + MomExt(1:4,3) - (1d0-x)*MomExt(1:4,2) ! tb
         MomExtTd(1:4,4) = MomExt(1:4,5) ! t

         !     splitting 1: A/Quark_I --> A/Quark_I + Gluon_F
         Dipole = Dipole * (2d0/(1d0-x+u)-1d0-x)

         Dipole = - Dipole ! from color matrix

      elseif (nDipole.eq.3) then

         a=1; i=3; j=5; !initial-final, q t
         
         sai = 2d0*(MomExt(1:4,a)).dot.(MomExt(1:4,i)) 
         saj = 2d0*(MomExt(1:4,a)).dot.(MomExt(1:4,j)) 
         sij = 2d0*(MomExt(1:4,i)).dot.(MomExt(1:4,j)) 
         
         x = 1d0-sij/(sai+saj)
         u = sai/(sai+saj)
         if( alpha_if.lt.u ) then
!         if( alpha_if.gt.u ) then !-- inverted alpha check
            Dipole = 0d0
            return
         endif
         Dipole = -1d0/sai/x
         
         MomExtTd(1:4,1) = x*MomExt(1:4,1) ! out q
         MomExtTd(1:4,2) = MomExt(1:4,2) ! out qb
         MomExtTd(1:4,3) = MomExt(1:4,4) ! tb
         MomExtTd(1:4,4) = MomExt(1:4,5) + MomExt(1:4,3) - (1d0-x)*MomExt(1:4,1) ! t

         !     splitting 1: A/Quark_I --> A/Quark_I + Gluon_F
         Dipole = Dipole * (2d0/(1d0-x+u)-1d0-x)

         Dipole = - Dipole ! from color matrix

      elseif (nDipole.eq.4) then

         a=2; i=3; j=5; !initial-final, qb t
         
         sai = 2d0*(MomExt(1:4,a)).dot.(MomExt(1:4,i)) 
         saj = 2d0*(MomExt(1:4,a)).dot.(MomExt(1:4,j)) 
         sij = 2d0*(MomExt(1:4,i)).dot.(MomExt(1:4,j)) 
         
         x = 1d0-sij/(sai+saj)
         u = sai/(sai+saj)
         if( alpha_if.lt.u ) then
!         if( alpha_if.gt.u ) then !-- inverted alpha check
            Dipole = 0d0
            return
         endif
         Dipole = -1d0/sai/x
         
         MomExtTd(1:4,1) = MomExt(1:4,1) ! out q
         MomExtTd(1:4,2) = x*MomExt(1:4,2) ! out qb
         MomExtTd(1:4,3) = MomExt(1:4,4) ! tb
         MomExtTd(1:4,4) = MomExt(1:4,5) + MomExt(1:4,3) - (1d0-x)*MomExt(1:4,2) ! t

         !     splitting 1: A/Quark_I --> A/Quark_I + Gluon_F
         Dipole = Dipole * (2d0/(1d0-x+u)-1d0-x)

      elseif (nDipole.eq.5) then

         j=4; i=3; a=1 !final-initial, tb q

         sai = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,i))
         saj = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,j))
         sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))

         x = 1d0-sij/(sai+saj)
         z = saj/(sai+saj)
         if( alpha_fi.lt.1d0-x ) then  ! NEW
!         if( alpha_fi.gt.1d0-x ) then  ! NEW !-- inverted alpha check
            Dipole = 0d0
            return
         endif
         Dipole = -1d0/sij/x

         MomExtTd(1:4,1) = x*MomExt(1:4,1) ! out q
         MomExtTd(1:4,2) = MomExt(1:4,2) ! out qb
         MomExtTd(1:4,3) = MomExt(1:4,4) + MomExt(1:4,3) - (1d0-x)*MomExt(1:4,1) ! tb
         MomExtTd(1:4,4) = MomExt(1:4,5) ! out t

         !     splitting 5: A/Quark_F --> A/Quark_F + Gluon_F
         Dipole = Dipole * (2d0/(2d0-z-x)-1d0-z-2d0*m_Top**2/sij)

      elseif (nDipole.eq.6) then

         j=4; i=3; a=2 !final-initial, tb qb

         sai = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,i))
         saj = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,j))
         sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))

         x = 1d0-sij/(sai+saj)
         z = saj/(sai+saj)
         if( alpha_fi.lt.1d0-x ) then  ! NEW
!         if( alpha_fi.gt.1d0-x ) then  ! NEW !-- inverted alpha check
            Dipole = 0d0
            return
         endif
         Dipole = -1d0/sij/x

         MomExtTd(1:4,1) = MomExt(1:4,1) ! out q
         MomExtTd(1:4,2) = x*MomExt(1:4,2) ! out qb
         MomExtTd(1:4,3) = MomExt(1:4,4) + MomExt(1:4,3) - (1d0-x)*MomExt(1:4,2) ! tb
         MomExtTd(1:4,4) = MomExt(1:4,5) ! out t

         !     splitting 5: A/Quark_F --> A/Quark_F + Gluon_F
         Dipole = Dipole * (2d0/(2d0-z-x)-1d0-z-2d0*m_Top**2/sij)

         Dipole = - Dipole ! from color matrix

      elseif (nDipole.eq.7) then

         j=5; i=3; a=1 !final-initial, t q

         sai = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,i))
         saj = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,j))
         sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))

         x = 1d0-sij/(sai+saj)
         z = saj/(sai+saj)
         if( alpha_fi.lt.1d0-x ) then  ! NEW
!         if( alpha_fi.gt.1d0-x ) then  ! NEW !-- inverted alpha check
            Dipole = 0d0
            return
         endif
         Dipole = -1d0/sij/x

         MomExtTd(1:4,1) = x*MomExt(1:4,1) ! out q
         MomExtTd(1:4,2) = MomExt(1:4,2) ! out qb
         MomExtTd(1:4,3) = MomExt(1:4,4) ! out tb
         MomExtTd(1:4,4) = MomExt(1:4,5) + MomExt(1:4,3) - (1d0-x)*MomExt(1:4,1) ! t

         !     splitting 5: A/Quark_F --> A/Quark_F + Gluon_F
         Dipole = Dipole * (2d0/(2d0-z-x)-1d0-z-2d0*m_Top**2/sij)

         Dipole = - Dipole ! from color matrix

      elseif (nDipole.eq.8) then

         j=5; i=3; a=2 !final-initial, t qb

         sai = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,i))
         saj = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,j))
         sij = 2d0*(MomExt(1:4,i).dot.MomExt(1:4,j))

         x = 1d0-sij/(sai+saj)
         z = saj/(sai+saj)
         if( alpha_fi.lt.1d0-x ) then  ! NEW
!         if( alpha_fi.gt.1d0-x ) then  ! NEW !-- inverted alpha check
            Dipole = 0d0
            return
         endif
         Dipole = -1d0/sij/x

         MomExtTd(1:4,1) = MomExt(1:4,1) ! out q
         MomExtTd(1:4,2) = x*MomExt(1:4,2) ! out qb
         MomExtTd(1:4,3) = MomExt(1:4,4) ! out tb
         MomExtTd(1:4,4) = MomExt(1:4,5) + MomExt(1:4,3) - (1d0-x)*MomExt(1:4,2) ! t

         !     splitting 5: A/Quark_F --> A/Quark_F + Gluon_F
         Dipole = Dipole * (2d0/(2d0-z-x)-1d0-z-2d0*m_Top**2/sij)

      endif
      
      return

    END SUBROUTINE Dipoles_qqb_Zprime_interf


  END MODULE ModDipoles_Zprime
