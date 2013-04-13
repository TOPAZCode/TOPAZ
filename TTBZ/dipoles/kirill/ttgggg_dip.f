! this is the file to subtract dipoles
      module ttgggg_dip
%      use
      implicit none
      private 

      public :: gg_ttgggg_dip

      contains

--- so, we assume that there is a list of dipoles that we need to go through
--- I will label things as 0-> bar t(p1)+t(p2) + g(p3)+g(p4)+g(p5) + g(p6) 
--- and I will assume that gluons in the initial state have momenta p3 and p4
--- so only p5 or p6 can be emitted

-- after that, the list of dipoles looks like that 

      subroutine gg_ttgggg_dip(res)
      real(dp), intent(out), res
      integer, parameter :: ndip = 40
      integer, parameter :: in1 = 3
      integer, parameter :: in2 = 4
      integer ::  dip(20,3)
      data dip(1,1)/5/ , dip(1,2)/3/, dip(1,3)/1/
      data dip(2,1)/5/ , dip(2,2)/3/, dip(2,3)/2/
      data dip(3,1)/5/ , dip(3,2)/3/, dip(3,3)/4/
      data dip(4,1)/5/ , dip(4,2)/3/, dip(4,3)/6/
      data dip(5,1)/5/ , dip(5,2)/4/, dip(5,3)/1/
      data dip(6,1)/5/ , dip(6,2)/4/, dip(6,3)/2/
      data dip(7,1)/5/ , dip(7,2)/4/, dip(7,3)/3/
      data dip(8,1)/5/ , dip(8,2)/4/, dip(8,3)/6/
      data  dip(9,1)/5/ ,  dip(9,2)/1/, dip(9,3)/2/
      data dip(10,1)/5/ , dip(10,2)/1/, dip(10,3)/3/
      data dip(11,1)/5/ , dip(11,2)/1/, dip(11,3)/4/
      data dip(12,1)/5/ , dip(12,2)/1/, dip(12,3)/6/
      data dip(13,1)/5/ , dip(13,2)/2/, dip(13,3)/1/
      data dip(14,1)/5/ , dip(14,2)/2/, dip(14,3)/3/
      data dip(15,1)/5/ , dip(15,2)/2/, dip(15,3)/4/
      data dip(16,1)/5/ , dip(16,2)/2/, dip(16,3)/6/
      data dip(17,1)/5/ , dip(17,2)/6/, dip(17,3)/1/
      data dip(18,1)/5/ , dip(18,2)/6/, dip(18,3)/2/
      data dip(19,1)/5/ , dip(19,2)/6/, dip(19,3)/3/
      data dip(20,1)/5/ , dip(20,2)/6/, dip(20,3)/4/
      data dip(21,1)/6/ , dip(21,2)/3/, dip(21,3)/1/
      data dip(22,1)/6/ , dip(22,2)/3/, dip(22,3)/2/
      data dip(23,1)/6/ , dip(23,2)/3/, dip(23,3)/4/
      data dip(24,1)/6/ , dip(24,2)/3/, dip(24,3)/5/
      data dip(25,1)/6/ , dip(25,2)/4/, dip(25,3)/1/
      data dip(26,1)/6/ , dip(26,2)/4/, dip(26,3)/2/
      data dip(27,1)/6/ , dip(27,2)/4/, dip(27,3)/3/
      data dip(28,1)/6/ , dip(28,2)/4/, dip(28,3)/5/
      data dip(29,1)/6/, dip(29,2)/1/, dip(29,3)/2/
      data dip(30,1)/6/, dip(30,2)/1/, dip(30,3)/3/
      data dip(31,1)/6/, dip(31,2)/1/, dip(31,3)/4/
      data dip(32,1)/6/, dip(32,2)/1/, dip(32,3)/5/
      data dip(33,1)/6/ , dip(33,2)/2/, dip(33,3)/1/
      data dip(34,1)/6/ , dip(34,2)/2/, dip(34,3)/3/
      data dip(35,1)/6/ , dip(35,2)/2/, dip(35,3)/4/
      data dip(36,1)/6/ , dip(36,2)/2/, dip(36,3)/5/
      data dip(37,1)/6/ , dip(37,2)/5/, dip(37,3)/1/
      data dip(38,1)/6/ , dip(38,2)/5/, dip(38,3)/2/
      data dip(39,1)/6/ , dip(39,2)/5/, dip(39,3)/3/
      data dip(40,1)/6/ , dip(40,2)/5/, dip(40,3)/4/
      real(dp), mass(6)

!-----flavor list qm -- massive quark, qu -- massless quark
!----- gl -- gluon
      character :: fl(6)*3
      integer :: n,ndip,i,j,k

      fl =    'glu'
      fl(1) = 'qm'
      fl(2) = 'qm'

!     mass assignment 

      mass = zero
      mass(1) = mt 
      mass(2) = mt 


            do n=1,ndip

         i=dip(n,1)  ! emitted
         j=dip(n,2)  ! emittor  
         k=dip(n,3)  ! spectator

      if (j.ne.in1.and.j.ne.in2.and.k.ne.in1.and.k.ne.in2) then 
      call dipff(n,i,j,k,mass(i),mass(j),mass(k),fl(i),fl(j),p,rrr)
      endif
         
         if ((j.eq.in1.or.j.eq.in2).and.(k.eq.in1.or.k.eq.in2)) then 
      call dipfi(n,i,j,k,mass(i),mass(j),mass(k),fl(i),fl(j),p,rrr)
         endif 

         if ( (j.eq.in1.and.k.ne.in2).or.(j.eq.in2.and.k.ne.in1)) then 
      call dipif(n,i,j,k,mass(i),mass(j),mass(k),fl(i),fl(j),p,rrr)
         endif    


         if ( (j.ne.in1.and.k.eq.in2).or.(j.ne.in2.and.k.eq.in1)) then 
      call dipii(n,i,j,k,mass(i),mass(j),mass(k),fl(i),fl(j),p,rrr)
         endif    
             
      enddo
      
      end subroutine

       double precision function vl(x1,x2,x3
       real(dp), x1,x2,x3
       
       vl = x1**2+x2**2+x3**2-2.0_dp*x1*x2-2.0_dp*x1*x3
     , - 2.0_dp*x2*x3
       end function vl


c-------dipole subroutines

       subroutine dipff(n,i,j,k,mi,mj,mk,fl1,fl2,p)
       implicit none
       integer, intent(in) :: n,i,j,k
       real(dp), intent(in) :: mi,mj,mk
       character, intent(in) :: fl1*3,fl2*3
       real(dp), intent(in) :: p(5,4)
       integer :: oh(6,5)
       real(dp) :: mij
       real(dp) :: C(6,6)
       data oh(1,1)/1/, oh(1,2)/2/, oh(1,3)/3/, oh(1,4)/4/, oh(1,5)/5/
       data oh(2,1)/1/, oh(2,2)/2/, oh(2,3)/3/, oh(2,4)/5/, oh(2,5)/4/
       data oh(3,1)/1/, oh(3,2)/2/, oh(3,3)/4/, oh(3,4)/3/, oh(3,5)/5/
       data oh(4,1)/1/, oh(4,2)/2/, oh(4,3)/4/, oh(4,4)/5/, oh(4,5)/3/
       data oh(5,1)/1/, oh(5,2)/2/, oh(5,3)/5/, oh(5,4)/3/, oh(5,5)/4/
       data oh(6,1)/1/, oh(6,2)/2/, oh(6,3)/5/, oh(6,4)/4/, oh(6,5)/3/

!       momentum mapping 
        pi=p(i,:)
        pj=p(j,:) 
        pk=p(k,:)

        qq = pi+pj+pk
        qq2 = sc(qq,qq)
        qqij = pi+pj
        qqij2 = sc(qqij,qqij)
        
        if (mi.eq.zero) mij = mj
        
        pkt = sqrt(vl(qq2,mij**2,mk**2))/sqrt(vl(qq2,pij2,mk**2))*(
     ,pk - sc(qq,pk)/q2*qq) + (q2+mk**2-mij**2)/2/q2*qq
        pij = qq - pkt

        zi = sc(pi,pk)/(sc(pi,pk)+sc(pj,pk))
        zj = 1-zi
        yijk = sc(pi,pj)/(sc(pi,pj)+sc(pi,pk)+sc(pj,pk))

        vijk = sqrt(vl(qq2,qij2,mk**2))/(qq2-qij2-mk**2)

       zim = zi - one/two*(one-vijk)
       zjm = zj - one/two*(one-vijk)

       paux = zim*pi - zjm*pj

          q(1,:) = p(1,:)
          fl(1,:) = 'qm'
          q(2,:) = p(2,:)
          fl(2,:) = 'qm'
          q(3,:) = p(3,:)
          fl(3,:) = 'gl'
          q(4,:) = p(4,:)
          fl(4,:) = 'gl'
          q(5,:) = p(6,:)         
          fl(5,:) = 'gl'

       if (i.eq.6) then 
          q(j,:) = pij(:)
          q(k,:) = pkt(:) 
          pos = j
       endif 

       if (i.eq.5) then 

            if (j.eq.6) then 
            q(5,:) = pij(:)
            pos = 5
            else 
            q(j,:) = pij(:)
            pos = j
            endif 

            if(k.eq.6) then 
            q(5,:) = pkt(:)
             else 
            q(k,:) = pkt(:) 
            endif 
       endif      

        call colorcor_ttgggg(n,C)
       
!--- after momentum mapping -- sum over colors and polarizations


       do i1=-1,1,2

          do i2 = -1,1,2
             do i3 = -1,1,2
                do i4 = -1,1,2
                   do i5=-1,1,2

           hel(1) = i1
           hel(2) = i2
           hel(3) = i3
           hel(4) = i4
           hel(5) = i5

           if (pos.eq.3) then 
              hel(3) = i1
              hel(1) = i3
           endif

           if (pos.eq.4) then 
              hel(4) = i1
              hel(1) = i4
           endif 

           if (pos.eq.5) then 
              hel(5) = i1
              hel(1) =i5
           endif    

      call GenerateEvent5(q,hel,ExtParticles(1:5))

       if (pos.eq.3) then 
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsReal(1)%gLUONS(1)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsReal(1)%gLUONS(1)%pOL(1:4)       
      endif 

      if (pos.eq.4) then 
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsReal(1)%gLUONS(2)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsReal(1)%gLUONS(2)%pOL(1:4)       
      endif 

      if (pos.eq.5) then 
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsReal(1)%gLUONS(3)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsReal(1)%gLUONS(3)%pOL(1:4)       
      endif

      do i = 1,6
      call EvalTree(TreeAmpsReal(i),ResAmpsReal(i1,i))
      enddo

      enddo
      enddo
      enddo
      enddo
      enddo
     
!     now build the helicity matrix 

      if ( (fl1.eq.'gl').and.(fl2.eq.'gl') ) then 

      diag= one/(one-zi*(one-yijk))+one/(one-zj*(one-yijk))
     , -(two - kap*zp*zm)/vijk 
       offdiag = 1/vijk/sc(pi,pj)      

       xm = sc(POL1(-1,:),paux)  
       xp = sc(POL1(1,:),paux)  

      HH(-1,1)=offdiag*cnjg(xm)*xp
      HH(1,-1)=offdiag*cnjg(xp)*xm
      HH(-1,-1)=diag + offdiag*cnjg(xm)*xp
      HH(1,1) = diag + offdiag*cnjg(xp)*xp

      endif    


      if ( (ff1.eq.'gl').and.(fl2.eq.'qm') ) then 
       diag = two/(one-zj*(one-yijk))
     ,-tvijk/vijk*(one + zj+mj**2/sc(pi,pj))
       HH(-1,1) = zero
       HH(1,-1) = zero 
       HH(1,1) = diag
       HH(-1,-1) = diag
      endif 


      do i1 = -1,1,2
         do i2 = -1,1,2
        do i=1,6
           do j=1,6

      res = res + C(i,j)
     ,*ResAmpsReal(i1,i)*cnjg(ResAmpsReal(i1,j))*HH(i1,i2)

            enddo 
          enddo

          if ( (fl1.eq.'gl').and.(fl2.eq.'gl')) then 

         endif 



         enddo 
         enddo 
         enddo 
         enddo 
         enddo 


       end subroutine dipff



       subroutine dipfi(n,i,j,a,mi,mi,ma,fl1,fl2,p)
       implicit none
       integer, intent(in) :: n,i,j,k
       real(dp), intent(in) :: mi,mj,mk
       character, intent(in) :: fl1*3,fl2*3
       real(dp), intent(in) :: p(5,4)
       integer :: oh(6,5)
       real(dp) :: mij
       real(dp) :: C(6,6)
       data oh(1,1)/1/, oh(1,2)/2/, oh(1,3)/3/, oh(1,4)/4/, oh(1,5)/5/
       data oh(2,1)/1/, oh(2,2)/2/, oh(2,3)/3/, oh(2,4)/5/, oh(2,5)/4/
       data oh(3,1)/1/, oh(3,2)/2/, oh(3,3)/4/, oh(3,4)/3/, oh(3,5)/5/
       data oh(4,1)/1/, oh(4,2)/2/, oh(4,3)/4/, oh(4,4)/5/, oh(4,5)/3/
       data oh(5,1)/1/, oh(5,2)/2/, oh(5,3)/5/, oh(5,4)/3/, oh(5,5)/4/
       data oh(6,1)/1/, oh(6,2)/2/, oh(6,3)/5/, oh(6,4)/4/, oh(6,5)/3/

!       momentum mapping 
        pi=p(i,:)
        pj=p(j,:) 
        pa=(-1.0_dp)*p(a,:)  ! sign change because ``a'' is the inital

        if (mi.eq.zero) mij = mj

        xija = ( sc(pa,pi)+sc(pa,pj)-sc(pi,pj)
     ,+one/two*(mij**2-mi**2-mj**2))/(sc(pa,pi)+sc(pa,pj))
        zi = sc(pa,pi)/(sc(pa,pi)+sc(pa,pj))
        zj = sc(pa,pj)/(sc(pa,pi)+sc(pa,pj))

        pij = pi+pj-(1-xija)*pa
        pat = xija*pa



        qq = pi+pj+pk
        qq2 = sc(qq,qq)
        qqij = pi+pj
        qqij2 = sc(qqij,qqij)
        

        pkt = sqrt(vl(qq2,mij**2,mk**2))/sqrt(vl(qq2,pij2,mk**2))*(
     ,pk - sc(qq,pk)/q2*qq) + (q2+mk**2-mij**2)/2/q2*qq
        pij = qq - pkt

        zi = sc(pi,pk)/(sc(pi,pk)+sc(pj,pk))
        zj = 1-zi
        yijk = sc(pi,pj)/(sc(pi,pj)+sc(pi,pk)+sc(pj,pk))

        vijk = sqrt(vl(qq2,qij2,mk**2))/(qq2-qij2-mk**2)

       zim = zi - one/two*(one-vijk)
       zjm = zj - one/two*(one-vijk)

       paux = zim*pi - zjm*pj

          q(1,:) = p(1,:)
          fl(1,:) = 'qm'
          q(2,:) = p(2,:)
          fl(2,:) = 'qm'
          q(3,:) = p(3,:)
          fl(3,:) = 'gl'
          q(4,:) = p(4,:)
          fl(4,:) = 'gl'
          q(5,:) = p(6,:)         
          fl(5,:) = 'gl'

       if (i.eq.6) then 
          q(j,:) = pij(:)
          q(k,:) = pkt(:) 
          pos = j
       endif 

       if (i.eq.5) then 

            if (j.eq.6) then 
            q(5,:) = pij(:)
            pos = 5
            else 
            q(j,:) = pij(:)
            pos = j
            endif 

            if(k.eq.6) then 
            q(5,:) = pkt(:)
             else 
            q(k,:) = pkt(:) 
            endif 
       endif      

        call colorcor_ttgggg(n,C)
       
!--- after momentum mapping -- sum over colors and polarizations


       do i1=-1,1,2

          do i2 = -1,1,2
             do i3 = -1,1,2
                do i4 = -1,1,2
                   do i5=-1,1,2

           hel(1) = i1
           hel(2) = i2
           hel(3) = i3
           hel(4) = i4
           hel(5) = i5

           if (pos.eq.3) then 
              hel(3) = i1
              hel(1) = i3
           endif

           if (pos.eq.4) then 
              hel(4) = i1
              hel(1) = i4
           endif 

           if (pos.eq.5) then 
              hel(5) = i1
              hel(1) =i5
           endif    

      call GenerateEvent(Mom,hel,ExtParticles(1:5))

       if (pos.eq.3) then 
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsReal(iTree)%gLUONS(1)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsReal(iTree)%gLUONS(1)%pOL(1:4)       
      endif 

      if (pos.eq.4) then 
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsReal(iTree)%gLUONS(2)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsReal(iTree)%gLUONS(2)%pOL(1:4)       
      endif 

      if (pos.eq.5) then 
      if (i1.eq.-1) POL1(i1,:)= TreeAmpsReal(iTree)%gLUONS(3)%pOL(1:4)
      if (i1.eq.1)  POL1(i1,:)= TreeAmpsReal(iTree)%gLUONS(3)%pOL(1:4)       
      endif

      do i = 1,6
      call EvalTree(TreeAmpsReal(i),ResAmpsReal(i1,i))
      enddo

      enddo
      enddo
      enddo
      enddo
      enddo
     

!     now build the helicity matrix 

      if ( (fl1.eq.'gl').and.(fl2.eq.'gl') ) then 

      diag= one/(one-zi*(one-yijk))+one/(one-zj*(one-yijk))
     , -(two - kap*zp*zm)/vijk 
       offdiag = 1/vijk/sc(pi,pj)      

       xm = sc(POL1(-1,:),paux)  
       xp = sc(POL1(1,:),paux)  

      HH(-1,1)=offdiag*cnjg(xm)*xp
      HH(1,-1)=offdiag*cnjg(xp)*xm
      HH(-1,-1)=diag + offdiag*cnjg(xm)*xp
      HH(1,1) = diag + offdiag*cnjg(xp)*xp

      endif    


      if ( (ff1.eq.'gl').and.(fl2.eq.'qm') ) then 
       diag = two/(one-zj*(one-yijk))
     ,-tvijk/vijk*(one + zj+mj**2/sc(pi,pj))
       HH(-1,1) = zero
       HH(1,-1) = zero 
       HH(1,1) = diag
       HH(-1,-1) = diag
      endif 


      do i1 = -1,1,2
         do i2 = -1,1,2
        do i=1,6
           do j=1,6

      res = res + C(i,j)
     ,*ResAmpsReal(i1,i)*cnjg(ResAmpsReal(i1,j))*HH(i1,i2)

            enddo 
          enddo

          if ( (fl1.eq.'gl').and.(fl2.eq.'gl')) then 

         endif 



         enddo 
         enddo 
         enddo 
         enddo 
         enddo 


       end subroutine dipfi



