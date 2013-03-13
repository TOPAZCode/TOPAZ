MODULE ModMyWeylRecurrence
use ModProcess
use ModMisc
use ModMyRecurrence
implicit none



public  :: f_Weyl,bf_Weyl,cur_f_2fW_Weyl
integer,parameter,private :: Dv=4,Ds=4
real(8),parameter :: PropCut = 1.0d-10


CONTAINS



! ---------------------------- currents with a W boson coupling to a top-bot-line -----------------------------------------



FUNCTION cur_f_2fW_Weyl(Gluons,Quarks,Boson,NumGlu) result(Res)           ! Quarks(:) DOES include the OFF-shell quark, in contrast to all other routines!
implicit none
complex(8) :: Res(1:Ds)
integer :: NumGlu(0:2),i,rIn,rOut,Quark1PartType
type(PtrToParticle) :: Gluons(1:),Boson,Quarks(1:2)      ! off-shell quark is not included in Quarks(:)
complex(8) :: GluMom(1:Dv,NumGlu(0)), QuarkMom(1:Dv)
complex(8) :: GluPol(1:Dv,NumGlu(0)), QuarkPol(1:Ds)
character :: FerFla1*3,FerFla2*3
integer :: PartKey,HelKey,CurrKey,Hel_Tmp
real(8) :: Quark1Mass



   do i=1,NumGlu(0)
    GluMom(1:Dv,i) = Gluons(i)%Mom(1:Dv)
    GluPol(1:Dv,i) = Gluons(i)%Pol(1:Dv)
   enddo
   QuarkMom(1:Dv) = Quarks(2)%Mom(1:Dv)
   QuarkPol(1:Ds) = Quarks(2)%Pol(1:Ds)

   if( abs(Quarks(1)%PartType).eq.5 ) then!  5=Top_
      FerFla1="top"
   else
      call Error("This Flavor is not allowed in cur_f_2fW_Weyl",Quark1PartType)
   endif

   if( abs(Quarks(2)%PartType).eq.4 ) then!  4=Str_
     FerFla2="str"
   else
      call Error("This Flavor is not allowed in cur_f_2fW_Weyl",Quarks(2)%PartType)
   endif

   if( NumGlu(0)-NumGlu(1)-NumGlu(2).ne.0 ) call Error("Wrong NumGlu in cur_f_2fV",NumGlu(0)-NumGlu(1)-NumGlu(2))


   rIn =1
   rOut=NumGlu(0)
   if( Quarks(2)%PartType .gt.0 ) then      !    X----->----
      Res(:) = fW_Weyl(GluPol(1:Dv,rIn:rOut),GluMom(1:Dv,rIn:rOut),QuarkPol(1:Ds),QuarkMom(1:Dv),Quarks(1)%Mass,FerFla2,FerFla1,Boson%Pol(1:Dv),Boson%Mom(1:Dv),NumGlu(1))
   else                                     !    X-----<----
      Res(:) = bfW_Weyl(GluPol(1:Dv,rIn:rOut),GluMom(1:Dv,rIn:rOut),QuarkPol(1:Ds),QuarkMom(1:Dv),Quarks(1)%Mass,FerFla2,FerFla1,Boson%Pol(1:Dv),Boson%Mom(1:Dv),NumGlu(1))
   endif


return
END FUNCTION







      recursive function f_Weyl(e,k,sp,p,mass,flout,flin,ms) result(res)
      implicit none
      complex(8), intent(in) :: e(:,:), k(:,:)
      complex(8), intent(in) :: sp(:), p(:)
      character,  intent(in) :: flin*3  ! flavor of off-shell f-line
      character,  intent(in)  :: flout*3 ! flavor of on-shell f-line
      integer, intent(in) ::  ms
      integer             :: ms1,m,ng1, ng2
      integer :: ngluon
      complex(8)             :: res(size(sp))
      complex(8)             :: tmp(size(sp))
      complex(8)             :: k1(size(p))
      complex(8)             :: k2(size(p))
      complex(8)             :: sp2(size(sp))
      complex(8)             :: sp3(size(sp))
      complex(8)             :: e1(size(e,dim=1))
      complex(8)             :: e2(size(e,dim=1))
      complex(8)  :: k1sq,k2sq,k3sq
      real(8) :: mass

      ngluon = size(e,dim=2)
      ng1 = ms   !#gluons to the left of a f-line
      ng2 = ngluon - ms  !#gluons to the right of the f-line

      if (flout.ne.flin) then
         res = (0d0,0d0)
      else

      if (ng2 < 0) write(*,*) 'WRONG DEFINITION OF CURRENT A'

      if (ngluon == 0) then
         res = sp

      else

         res = (0d0,0d0)

       do m=0,ng2-1
           k1 = sum(k(:,ng1+1+m:ngluon),dim=2)
           e1=g(e(:,ng1+1+m:ngluon),k(:,ng1+1+m:ngluon))
           k1sq=sc_(k1,k1)

           k2 = sum(k(:,1:ng1+m),dim=2)
           k2 = k2 + p
           k2sq = sc_(k2,k2)-mass**2
           sp2 = f_Weyl(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,mass,flout,flin,ng1)
           if (ng1 >0.or.m>0) sp2 = spb2_Weyl(sp2,k2)+mass*sp2
           tmp = vqg_Weyl(sp2,e1)

           if (m < ng2-1)  then
              if(abs(k1sq) > propcut) then
                  tmp = -(0d0,1d0)/k1sq*tmp
              else
                  tmp = (0d0,0d0)
              endif
           endif

           if (ng1>0.or.m>0) then
              if (abs(k2sq) > propcut) then
                  tmp =  (0d0,1d0)/k2sq*tmp
              else
                  tmp = (0d0,0d0)
               endif
           endif
           res = res + tmp
        enddo

        do m=1,ng1

           k1 = sum(k(:,1:m),dim=2)
           e1=g(e(:,1:m),k(:,1:m))
           k1sq = sc_(k1,k1)

           k2 = sum(k(:,m+1:ngluon),dim=2)
           k2 = k2 + p
           k2sq = sc_(k2,k2) - mass**2
           ms1 = ng1 - m
           sp2=f_Weyl(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,mass,flout,flin,ms1)

           if (ng2 > 0.or.m < ng1) sp2 = spb2_Weyl(sp2,k2)+mass*sp2
           tmp = vgq_Weyl(e1,sp2)
           if (m > 1) then
              if (abs(k1sq) > propcut) then
                  tmp=-(0d0,1d0)/k1sq*tmp
              else
                  tmp = (0d0,0d0)
              endif
           endif

           if (ng2 > 0.or. m < ng1) then
            if (abs(k2sq) > propcut) then
              tmp=(0d0,1d0)/k2sq*tmp
            else
              tmp = (0d0,0d0)
            endif
           endif

           res = res + tmp
           enddo

         endif
         endif ! endif for flavor consistency condition

      end function f_Weyl



      recursive function bf_Weyl(e,k,sp,p,mass,flout,flin,ms) result(res)
      implicit none
      complex(8), intent(in) :: e(:,:), k(:,:)
      complex(8), intent(in) :: sp(:), p(:)
      integer, intent(in) ::  ms
      character, intent(in) :: flout*3  ! flavor of on-shell f-line
      character, intent(in) :: flin*3   ! flavor of off-shell f-line
      integer             :: ms1,m,ng1, ng2
      integer :: ngluon
      complex(8)             :: res(size(sp))
      complex(8)             :: tmp(size(sp))
      complex(8)             :: k1(size(p))
      complex(8)             :: k2(size(p))
      complex(8)             :: sp2(size(sp))
      complex(8)             :: sp3(size(sp))
      complex(8)             :: e1(size(e,dim=1))
      complex(8)             :: e2(size(e,dim=1))
      complex(8)  :: k1sq,k2sq,k3sq
      real(8) :: mass


       if (flout.ne.flin) then
          res = (0d0,0d0)
       else

       ngluon = size(e,dim=2)
      ng1 = ms   !#gluons to the left of a f-line
      ng2 = ngluon - ms  !#gluons to the right of the f-line
      if (ng2 < 0) write(*,*) 'WRONG DEFINITION OF CURRENT B'
      if (ngluon == 0) then
         res = sp
      else

       res = (0d0,0d0)
       do m=0,ng2-1
           k1 = sum(k(:,ng1+1+m:ngluon),dim=2)
           e1=g(e(:,ng1+1+m:ngluon),k(:,ng1+1+m:ngluon))
           k1sq=sc_(k1,k1)

           k2 = sum(k(:,1:ng1+m),dim=2)
           k2 = -k2 - p
           k2sq = sc_(k2,k2)-mass**2
           sp2 = bf_Weyl(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,mass,flout,flin,ng1)
           if (ng1>0.or.m>0) sp2 = spi2_Weyl(k2,sp2)+mass*sp2

            tmp = vbqg_Weyl(sp2,e1)

           if (m < ng2-1) then
            if (abs(k1sq) > propcut) then
           tmp = -(0d0,1d0)/k1sq*tmp
            else
           tmp = (0d0,0d0)
            endif
            endif
           if (ng1>0.or.m>0) then
            if (abs(k2sq) > propcut) then
           tmp =  (0d0,1d0)/k2sq*tmp
            else
            tmp = (0d0,0d0)
             endif
             endif

           res = res + tmp


        enddo


        do m=1,ng1

           k1 = sum(k(:,1:m),dim=2)
           e1=g(e(:,1:m),k(:,1:m))
           k1sq = sc_(k1,k1)

           k2 = sum(k(:,m+1:ngluon),dim=2)
           k2 = -k2 - p
           k2sq = sc_(k2,k2) - mass**2
           ms1 = ng1 - m
           sp2=bf_Weyl(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,mass,flout,flin,ms1)

           if (ng2 > 0.or.m < ng1) sp2 = spi2_Weyl(k2,sp2)+mass*sp2

           tmp = vgbq_Weyl(e1,sp2)

           if (m > 1) then
               if (abs(k1sq) > propcut) then
              tmp=-(0d0,1d0)/k1sq*tmp
              else
              tmp = (0d0,0d0)
              endif
           endif

           if (ng2 > 0.or. m < ng1) then
              if (abs(k2sq) > propcut) then
              tmp=(0d0,1d0)/k2sq*tmp
              else
              tmp = (0d0,0d0)
              endif
           endif

           res = res + tmp

           enddo

          endif
          endif ! endif for flavor consisency condition
      end function bf_Weyl



















      recursive function fW_Weyl(e,k,sp,p,mass,flout,flin,eW,kW,ms) result(res)
      implicit none
      complex(8), intent(in) :: e(:,:), k(:,:)
      complex(8), intent(in) :: sp(:), p(:)
      complex(8), intent(in) :: eW(:), kW(:)
      character,  intent(in) :: flin*3  ! flavor of off-shell f-line
      character,  intent(in) :: flout*3 ! flavor of on-shell f-line
      integer, intent(in) ::  ms
      integer             :: ms1,m,ng1, ng2
      integer :: ngluon
      complex(8)             :: res(size(sp))
      complex(8)             :: tmp(size(sp))
      complex(8)             :: k1(size(p))
      complex(8)             :: k2(size(p))
      complex(8)             :: sp2(size(sp))
      complex(8)             :: sp3(size(sp))
      complex(8)             :: e1(size(e,dim=1))
      complex(8)             :: e2(size(e,dim=1))
      complex(8)  :: k1sq,k2sq,k3sq
      real(8) :: mass

      ngluon = size(e,dim=2)
      ng1 = ms   !#gluons to the left of a f-line
      ng2 = ngluon - ms  !#gluons to the right of the f-line


      if (ng2 < 0) write(*,*) 'WRONG DEFINITION OF CURRENT fW_Weyl'

if (ngluon == 0) then
         res = vbqW_Weyl(sp,eW)
else

       res = (0d0,0d0)
       do m=0,ng2-1
           k1 = sum(k(:,ng1+1+m:ngluon),dim=2)
           e1=g(e(:,ng1+1+m:ngluon),k(:,ng1+1+m:ngluon))
           k1sq=sc_(k1,k1)

           k2 = sum(k(:,1:ng1+m),dim=2)
           k2 = k2 + p + kW
           k2sq = sc_(k2,k2)-mass**2
           sp2 = fW_Weyl(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,mass,flout,flin,eW,kW,ng1)
           sp2 = spb2_Weyl(sp2,k2)+mass*sp2

           tmp = vqg_Weyl(sp2,e1)
           if (m < ng2-1)  then
              if(abs(k1sq) > propcut) then
                  tmp = -(0d0,1d0)/k1sq*tmp
              else
                  tmp = (0d0,0d0)
              endif
           endif

           if (abs(k2sq) > propcut) then
                  tmp =  (0d0,1d0)/k2sq*tmp
           else
                  tmp = (0d0,0d0)
           endif
           res = res + tmp
        enddo




        do m=1,ng1
           k1 = sum(k(:,1:m),dim=2)
           e1=g(e(:,1:m),k(:,1:m))
           k1sq = sc_(k1,k1)

           k2 = sum(k(:,m+1:ngluon),dim=2)
           k2 = k2 + p + kW
           k2sq = sc_(k2,k2) - mass**2
           ms1 = ng1 - m
           sp2=fW_Weyl(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,mass,flout,flin,eW,kW,ms1)

           sp2 = spb2_Weyl(sp2,k2)+mass*sp2

           tmp = vgq_Weyl(e1,sp2)
           if (m > 1) then
              if (abs(k1sq) > propcut) then
                  tmp=-(0d0,1d0)/k1sq*tmp
              else
                  tmp = (0d0,0d0)
              endif
           endif

           if (abs(k2sq) > propcut) then
              tmp=(0d0,1d0)/k2sq*tmp
           else
              tmp = (0d0,0d0)
           endif

           res = res + tmp
        enddo



        sp2 = f_Weyl(e,k,sp,p,0d0,flout,flout,ms) !  mass set to zero because this is a bottom quark
        k2 = sum(k(:,1:ngluon),dim=2)
        k2 = k2 + p
        k2sq = sc_(k2,k2)   !-mass**2

        sp2 = spb2_Weyl(sp2,k2) !+mass*sp2

        tmp = vbqW_Weyl(sp2,eW)
        if (abs(k2sq) > propcut) then
               tmp = (0d0,1d0)/k2sq*tmp
        else
                tmp = (0d0,0d0)
        endif
        res = res + tmp

endif

end function fW_Weyl




      recursive function bfW_Weyl(e,k,sp,p,mass,flout,flin,eW,kW,ms) result(res)
      implicit none
      complex(8), intent(in) :: e(:,:), k(:,:)
      complex(8), intent(in) :: sp(:), p(:)
      complex(8), intent(in) :: eW(:), kW(:)
      integer, intent(in) ::  ms
      character, intent(in) :: flout*3  ! flavor of on-shell f-line
      character, intent(in) :: flin*3   ! flavor of off-shell f-line
      integer             :: ms1,m,ng1, ng2
      integer :: ngluon
      complex(8)             :: res(size(sp))
      complex(8)             :: tmp(size(sp))
      complex(8)             :: k1(size(p))
      complex(8)             :: k2(size(p))
      complex(8)             :: sp2(size(sp))
      complex(8)             :: sp3(size(sp))
      complex(8)             :: e1(size(e,dim=1))
      complex(8)             :: e2(size(e,dim=1))
      complex(8)  :: k1sq,k2sq,k3sq
      real(8) :: mass


      ngluon = size(e,dim=2)
      ng1 = ms   !#gluons to the left of a f-line
      ng2 = ngluon - ms  !#gluons to the right of the f-line
      if (ng2 < 0) write(*,*) 'WRONG DEFINITION OF CURRENT fbW'

if (ngluon == 0) then
         res = vWq_Weyl(eW,sp)
else

       res = (0d0,0d0)
       do m=0,ng2-1
           k1 = sum(k(:,ng1+1+m:ngluon),dim=2)
           e1=g(e(:,ng1+1+m:ngluon),k(:,ng1+1+m:ngluon))
           k1sq=sc_(k1,k1)

           k2 = sum(k(:,1:ng1+m),dim=2)
           k2 = -k2 - p - kW
           k2sq = sc_(k2,k2)-mass**2

           sp2 = bfW_Weyl(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,mass,flout,flin,eW,kW,ng1)
           sp2 = spi2_Weyl(k2,sp2)+mass*sp2

           tmp = vbqg_Weyl(sp2,e1)

           if (m < ng2-1) then
            if (abs(k1sq) > propcut) then
                tmp = -(0d0,1d0)/k1sq*tmp
            else
                tmp = (0d0,0d0)
            endif
           endif

           if (abs(k2sq) > propcut) then
              tmp =  (0d0,1d0)/k2sq*tmp
           else
               tmp = (0d0,0d0)
           endif

           res = res + tmp
        enddo


        do m=1,ng1

           k1 = sum(k(:,1:m),dim=2)
           e1=g(e(:,1:m),k(:,1:m))
           k1sq = sc_(k1,k1)

           k2 = sum(k(:,m+1:ngluon),dim=2)
           k2 = -k2 - p - kW
           k2sq = sc_(k2,k2) - mass**2
           ms1 = ng1 - m
           sp2=bfW_Weyl(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,mass,flout,flin,eW,kW,ms1)

           sp2 = spi2_Weyl(k2,sp2)+mass*sp2

           tmp = vgbq_Weyl(e1,sp2)

           if (m > 1) then
              if (abs(k1sq) > propcut) then
                  tmp=-(0d0,1d0)/k1sq*tmp
              else
                  tmp = (0d0,0d0)
              endif
           endif

           if (abs(k2sq) > propcut) then
                 tmp=(0d0,1d0)/k2sq*tmp
              else
                 tmp = (0d0,0d0)
           endif

           res = res + tmp

        enddo

        sp2 = bf_Weyl(e,k,sp,p,0d0,flout,flout,ms) !  mass set to zero because this is a bottom quark
        k2 = sum(k(:,1:ngluon),dim=2)
        k2 = -k2 - p
        k2sq = sc_(k2,k2)   !-mass**2

        sp2 = spi2_Weyl(k2,sp2) !+mass*sp2

        tmp = vWq_Weyl(eW,sp2)
        if (abs(k2sq) > propcut) then
               tmp = (0d0,1d0)/k2sq*tmp
        else
               tmp = (0d0,0d0)
        endif
        res = res + tmp

endif





      end function bfW_Weyl








END MODULE ModMyWeylRecurrence



