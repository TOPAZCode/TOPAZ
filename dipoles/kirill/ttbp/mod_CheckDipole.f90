module ModCheckDipole
use ModParameters
use ModMisc
use ModKinematics
!--------------------------------
use Mod_TGAMTBGGG
use Mod_TGAMTQQG
use ModDipoles_TGAMTBGGG
use ModDipoles_QQ_TGAMTBG
use ModDipoles_QG_TGAMTBQ
use ModDipoles_QBG_TGAMTBQB
!---------------------------------
use ModIntDipoles_GG_TTBGAMG
use ModIntDipoles_QQ_TTBGAMG
!---------------------------------
use ModIntDipoles_QG_TTBGAMQ
use ModIntDipoles_QBG_TTBGAMQB
!---------------------------------
implicit none 
private 

      integer, parameter  :: dp = selected_real_kind(15)
      public :: checkdipoles


contains 

subroutine checkdipoles()
implicit none
integer :: iTree
integer, parameter :: Nout=4
integer :: hel(6), i 
real(dp) :: Mom(1:4,Nout+2),MomA(1:4,Nout+2),Mass(1:Nout+2)
real(dp) :: Mom1(1:4,Nout+2)
real(dp) :: Momr(1:4,5), Mom5(1:4,5)
real(dp) :: MomDK(1:4,1:6)
real(dp) :: Ehat, xRndPS(1:3*Nout-4),MassA(1:Nout+2), PSWgt
integer :: Pcol1, Pcol2, Steps
real(dp) :: SingDepth, PSWgt1,PSWgt2
real(dp) :: pjetout(4,6), weight, pdip(4,5), res(2), resdip(3),sum_dip(2)
real(dp) :: resdip3(1:2,1:3)
integer  :: Njet, j1, j2
real(dp) :: yRnDk(1:8), Wgt(2), z
!type(Particle) :: ExtParticles(1:6)
!type(TreeProcess) :: TreeAmpsReal(1:24)
!complex(8)        :: ResAmpsReal(1:24)
!type(TreeProcess) :: TreeAmpsDip(1:6)
!complex(8)        :: ResAmpsDip(1:6)


!----------  Random numbers for decays 


    yRnDk(1) = 0.15d0
    yRnDk(2) = 0.23d0
    yRnDk(3) = 0.2d0
    yRnDk(4) = 0.27d0
    yRnDk(5) = 0.40d0
    yRnDk(6) = 0.33d0
    yRnDk(7) = 0.76d0
    yRnDk(8) = 0.55d0    

    Wgt = 1d0

!----------------------------------------------------
!    Mom(lorentz index, 1-4, particle number)

    Ehat = 1960.0_dp
    MassA = zero
    MassA(3) = m_TOP
    MassA(4) = m_TOP

    xRndPS(1) = 0.1
    xRndPS(2) = 0.21
    xRndPS(3) = 0.38
    xRndPS(4) = 0.23
    xRndPS(5) = 0.47
    xRndPS(6) = 0.35
    xRndPS(7) = 0.76
    xRndPS(8) = 0.55    


!--- how to generate phase-space points in singular kinematics -----------------
!    1 and 2 are on the z-axis; ! for SOFT singularity, request Pcol1 = Pcol2
    do i=1,4
    Pcol1 =  5 - 1  ! particle 1  ; -1 is a convention 
    Pcol2 =  2 - 1  ! particle 3  ; -1 is a convention 
    Steps = 4
    SingDepth = 1d-7
    PSWgt = 1d0
    call gensing(Nout,Ehat,MassA(3:Nout+2),MomA(1:4,3:Nout+2), & 
    Pcol1,Pcol2,SingDepth,Steps)
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------
!    call genps(Nout,Ehat,xRndPS(1:3*Nout-4),MassA(3:Nout+2) &
! ,MomA(1:4,3:Nout+2),PSWgt)
!    PSWgt = PSWgt * (2d0*pi)**(4-Nout*3) * (4d0*pi)**(Nout-1)

    MomA(1,1) =   -1.0_dp*Ehat/2.0_dp
    MomA(2,1) =   zero
    MomA(3,1) =   zero
    MomA(4,1) =  -1.0_dp*Ehat/2.0_dp
    MomA(1,2) =  -1.0_dp*Ehat/2.0_dp
    MomA(2,2) =   zero
    MomA(3,2) =   zero
    MomA(4,2) =   1.0_dp*Ehat/2.0_dp


!--- now rewriting everything in the following order: (outgoing convention) 
!    tbar(1), t(2), g_initial(3), g_initial(4), g_final(5), ....

    Mass = zero 
    Mom(:,1) = MomA(:,3)
    Mass(1) = MassA(3) 
    Mom(:,2) = MomA(:,4)
    Mass(2) = MassA(4) 
    Mom(:,3) = MomA(:,1) 
    Mom(:,4) = MomA(:,2) 

     do i = 5,Nout+2
    Mom(:,i) = MomA(:,i)
     enddo


     print *, 'mass',Mass 
     print *, '----------------------'


!    Mom(1,1) = 893.329069311565d0       
!    Mom(2,1) = -380.916001324909D0   
!    Mom(3,1) = 673.706912170134D0     
!    Mom(4,1) = -410.406898706505d0     

!   Mom(1,2) = 787.238346369073d0   
!   Mom(2,2) = 455.644260470862d0  
!   Mom(3,2) = -476.574873698717d0     
!   Mom(4,2) = 392.917181671818d0     

!    Mom(1,3) = -980.000000000000d0    
!    Mom(2,3) = 0.000000000000000d0
!    Mom(3,3) = 0.000000000000000d0
!    Mom(4,3) = -980.000000000000d0     

!    Mom(1,4) = -980.000000000000d0   
!    Mom(2,4) = 0.000000000000000d0
!    Mom(3,4) = 0.000000000000000d0
!    Mom(4,4) = 980.000000000000d0
     
!    Mom(1,5) = 65.3791347507327d0      
!    Mom(2,5) = -0.944639100049881d0       
!    Mom(3,5) = -1.97830330233341d0     
!    Mom(4,5) = 65.3423693614572d0 

!---------------------------------------------- photon     
!    Mom(1,6) =  214.053449568629d0       
!    Mom(2,6) = -73.7836200459032d0       
!    Mom(3,6) = -195.153735169083d0     
!    Mom(4,6) = -47.8526523267703d0   


!   Mom(1,1) = 3.23802785116962       
!   Mom(2,1) =  -1.79079039615013      
!   Mom(3,1) =  -0.803579324606970     
!   Mom(4,1) =  -1.88935290273437     

!   Mom(1,2) = 3.36788686977249        
!   Mom(2,2) =  2.76685475072411      
!   Mom(3,2) =  -0.328113686517814     
!   Mom(4,2) =  -0.719039751822576     

!   Mom(1,3) = -9.80000000000000       
!   Mom(2,3) = 0.000000000000000E+000  
!   Mom(3,3) = 0.000000000000000E+000
!   Mom(4,3) = -9.80000000000000     

!   Mom(1,4) = -9.80000000000000     
!   Mom(2,4) = 0.000000000000000E+000  
!   Mom(3,4) = 0.000000000000000E+000
!   Mom(4,4) =  9.80000000000000     

!   Mom(1,5) =  7.71820733963581      
!   Mom(2,5) = -0.574576444471603      
!   Mom(3,5) = -0.214533319145605     
!   Mom(4,5) = 7.69380022498913     

!-------------------------------------   this is the photon
!   Mom(1,6) = 5.27587793942208      
!   Mom(2,6) = -0.401487910102377        
!   Mom(3,6) = 1.34622633027039     
!   Mom(4,6) = -5.08540757043218     


!-----------------------------------    to check int. dipoles 
 Mom(1,1) = 846.481298737869d0
 Mom(2,1) =  561.907802275647d0 
 Mom(3,1)= -69.8412892696248d0      
 Mom(4,1) = -604.390110109038d0 
     
  Mom(1,2) =  972.382375341243d0 
  Mom(2,2) =  -601.284201964586d0 
  Mom(3,2) = 89.8360707785765d0 
  Mom(4,2) = 738.437047261593d0 
    
  Mom(1,3) = -980.000000000000d0
  Mom(2,3)=0.000000000000000d0
  Mom(3,3) = 0.000000000000000d0 
  Mom(4,3) = -980.000000000000d0 
    
  Mom(1,4) = -980.000000000000d0 
  Mom(2,4) =  0.000000000000000d0 
  Mom(3,4) =  0.000000000000000d0
  Mom(4,4) =  980.000000000000d0 
     
  Mom(:,5) = 0d0 

  Mom(1,6) = 141.135223691194d0 
  Mom(2,6) = 39.3762137374185d0 
  Mom(3,6) = -19.9948163876085d0      
  Mom(4,6) = -134.048023023547d0     


     print *, Mom(:,1) 
     print *, Mom(:,2) 
     print *, Mom(:,3) 
     print *, Mom(:,4) 
     print *, Mom(:,5) 
     print *, Mom(:,6) 
     pause

!     print *, scr(Mom(:,1),Mom(:,1)) 
!     print *, scr(Mom(:,2),Mom(:,2)) 
!     print *, scr(Mom(:,3),Mom(:,3))
!     print *, scr(Mom(:,4),Mom(:,4)) 
!     print *, scr(Mom(:,5),Mom(:,5)) 
!     print *, scr(Mom(:,6),Mom(:,6)) 
!     print *, scr(Mom(:,7),Mom(:,7)) 
!     print *, scr(Mom(:,8),Mom(:,8)) 
!     print *, scr(Mom(:,9),Mom(:,9))
!     pause


!---- now call the amplitude for a particular subprocess 
      res = 0d0 
      sum_dip = 0d0

!-------various things for gg -> t tbar + gamma + glue
!      call TGAMTBGGG(Mom,yRnDk,Wgt(1),res(1))
!      print *, 'amplitude', res(1)
!      pause
!------ This result seems to be 32 times larger than the madgraph
!      print *, 'result_amplitude', & 
!      (4d0*pi*0.13d0)**3*(2d0/3d0)**2*4d0*pi/128d0*res/4d0

!      pause
!      call TGAMTBGG(Mom5,yRnDk,Wgt(1),res)
!      print *, 'res_ampl_approx', res 
!      pause
!---- now call the dipoles for gg -> tbartg+gamma ; 
!      here p is ordered as bar t gam t g3(in) g4(in) g(fin)
!      call EvalDipoles_TGAMTBGGG((/Mom(:,1),Mom(:,6),Mom(:,2), & 
!      Mom(:,3),Mom(:,4),Mom(:,5)/),yRnDk,Wgt(1),sum_dip(1))
!      print *, 'dipole', sum_dip(1) 

!-------- for q qbar -> tbar t + gamma + glue channel-------------
!      call  TGAMTQQG(Mom,yRnDk,Wgt(1),res)
!      print *, 'result_amplitude', res(1), res(2)
!      pause

!     call EvalDipoles_QQ_TGAMTBG((/Mom(:,1),Mom(:,6),Mom(:,2), & 
!     Mom(:,3),Mom(:,4),Mom(:,5)/),yRnDk,Wgt,sum_dip)
!      print *, 'result_dipole', sum_dip(1),sum_dip(2)
!      pause
!------------------------------------------------------------------


!------- for q glue -> tbar t + gamma + glue channel --------------
! we take barq_final(position3) and gluon(position5) as the inital momenta
!       call TGAMTQQG((/Mom(:,1),Mom(:,2),Mom(:,3),Mom(:,5), & 
!       Mom(:,4),Mom(:,6)/),yRnDk,Wgt(1),res)
!      print *, 'result_amplitude_qg', res(1), res(2)
!       pause
!--- input here tbar, gamma, t , qbar, q, g
!       call EvalDipoles_QG_TGAMTBQ((/Mom(:,1),Mom(:,6),Mom(:,2), & 
!      Mom(:,3),Mom(:,5),Mom(:,4)/),yRnDk,Wgt,sum_dip)
!           print *, 'result_dip', sum_dip(1), sum_dip(2) 
!       pause
!-----------------------------------------------------------------


!------- for barq glue -> tbar t + gamma + glue channel --------------
! we take q_final(position4) and gluon(position5) as the inital momenta
!       call TGAMTQQG((/Mom(:,1),Mom(:,2),Mom(:,5),Mom(:,3), & 
!       Mom(:,4),Mom(:,6)/),yRnDk,Wgt(1),res)
!      print *, 'result_amplitude_qbarg', res(1), res(2)
!       pause
!--- input here tbar, gamma, t , qbar, q, g
!       call EvalDipoles_QBG_TGAMTBQB((/Mom(:,1),Mom(:,6),Mom(:,2), & 
!      Mom(:,5),Mom(:,3),Mom(:,4)/),yRnDk,Wgt,sum_dip)
!           print *, 'result_dip', sum_dip(1), sum_dip(2) 
!       pause
!-----------------------------------------------------------------
       

!--- to check that the integrated dipoles actually produce a number 

      z = 0.3d0

       if (TopDecays.ge.1) then
         call EvalPhasespace_TopDecay(Mom(1:4,1),yRnDk(1:4),.false., &
  MomDK(1:4,1:3),PSWgt1)
         call EvalPhasespace_TopDecay(Mom(1:4,2),yRnDk(5:8),.false., &
  MomDK(1:4,4:6),PSWgt2)
       elseif(TopDecays.eq.0) then
         MomDK = 0d0 
         PSWgt1 = one
         PSWgt2 = one
       else
          print *, 'topdecays are not defined'
          pause
       endif


!-------------- GG/integrated
!      call EvalIntDipoles_GG_TTBGAMG((/Mom(:,1),Mom(:,6),Mom(:,2),Mom(:,3), & 
!      Mom(:,4)/),MomDK,z,resdip)
!      print *, resdip(1), resdip(2), resdip(3) 

!------------------------------------------------------------------


!-------------- QQB/integrated
!      call EvalIntDipoles_QQ_TTBGAMG((/Mom(:,1),Mom(:,6),Mom(:,2),Mom(:,3), & 
!      Mom(:,4)/),MomDK,z,resdip3)
!      print *, resdip3(1,1), resdip3(1,2), resdip3(1,3) 
!      print *, resdip3(2,1), resdip3(2,2), resdip3(2,3) 
!------------------------------------------------------------------


!-------------- QG/integrated  / input tbar gamma t q glue
!      call EvalIntDipoles_QG_TTBGAMQ((/Mom(:,1),Mom(:,6),Mom(:,2),Mom(:,3), & 
!      Mom(:,4)/),MomDK,z,resdip3)
!      print *, 'A', resdip3(1,1), resdip3(1,2), resdip3(1,3) 
!      print *, 'A', resdip3(2,1), resdip3(2,2), resdip3(2,3) 
!------------------------------------------------------------------


!-------------- QBG/integrated  / input tbar, gamma, t, qbar, glue
      call EvalIntDipoles_QBG_TTBGAMQB((/Mom(:,1),Mom(:,6),Mom(:,2),Mom(:,3), & 
      Mom(:,4)/),MomDK,z,resdip3)
      print *, 'B', resdip3(1,1), resdip3(1,2), resdip3(1,3) 
      print *, 'B', resdip3(2,1), resdip3(2,2), resdip3(2,3) 
!------------------------------------------------------------------





end subroutine checkdipoles


end module ModCheckDipole
