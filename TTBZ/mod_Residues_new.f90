       MODULE ModResidues_new
! this is the version before matching routines have been removed from here to mod_UCuts
! not 100% sure about precision when swiching to 128bit
       public  :: resid5_new,resid4_new,resid3_new,resid2_new,resid1_new
       public  :: resid4Subtr,resid3Subtr,resid2Subtr,resid1Subtr
       private :: givepol


       contains



       SUBROUTINE resid5_new(lv,k1,k2,k3,k4,l5c,TreeProcs,res)
       use ModAmplitudes
       use ModProcess
       use ModMisc
       use ModParameters
       implicit none
       include 'misc/global_import'
       integer Ds,i,tag_pol
       integer j,j1,j2,tag_f,tag_Z
       integer i1,i2,i3,i4,ia,ib,j3,j4,j5
       double complex lv(5)
       double complex k1(4),k2(4),k3(4),k4(4)
       complex(8) q1(1:8),q2(1:8),q3(1:8),q4(1:8),q5(1:8)
       complex(8),target :: q1neg(1:8),q2neg(1:8),q3neg(1:8),q4neg(1:8),q5neg(1:8)
       complex(8) BPOL1(8,16),POL1(8,16)
       complex(8) BPOL2(8,16),POL2(8,16)
       complex(8) BPOL3(8,16),POL3(8,16)
       complex(8) BPOL4(8,16),POL4(8,16)
       complex(8) BPOL5(8,16),POL5(8,16)
       complex(8) mur1(10,10),mur2(10,10)
       complex(8) mur3(10,10),mur4(10,10)
       complex(8) mur5(10,10)
       double complex res,resaux,res6,res8
       integer Nj1,Nj2,Nj3,Nj4,Nj5
       character lab1*3, lab2*3
       integer l5c(5)
       type(TreeProcess) :: TreeProcs(:)


          tag_pol = 0


          if( lv(5).eq.0d0 ) then
               Ds = 4
          else
               Ds = 5
          endif

          do i=1,4
          q1(i)=lv(i)
          q2(i)=lv(i)+k1(i)
          q3(i)=lv(i)+k2(i)
          q4(i)=lv(i)+k3(i)
          q5(i)=lv(i)+k4(i)
          enddo
!     I assume here that Ds = 5, always, to have a pentuple cut
          q1(5)=lv(5)
          q2(5)=lv(5)
          q3(5)=lv(5)
          q4(5)=lv(5)
          q5(5)=lv(5)
          q1(6:8)=(0d0,0d0)
          q2(6:8)=(0d0,0d0)
          q3(6:8)=(0d0,0d0)
          q4(6:8)=(0d0,0d0)
          q5(6:8)=(0d0,0d0)
          q1neg(1:8) = -q1(1:8)
          q2neg(1:8) = -q2(1:8)
          q3neg(1:8) = -q3(1:8)
          q4neg(1:8) = -q4(1:8)
          q5neg(1:8) = -q5(1:8)

!         set momentum vector for last particle
          if( TreeProcs(1)%PartType(TreeProcs(1)%NumPart).eq.Glu_ ) then
                TreeProcs(1)%Gluons(TreeProcs(1)%NumGlu(0))%Mom => q2neg(:)
          elseif( IsAQuark(TreeProcs(1)%PartType(TreeProcs(1)%NumPart)) ) then
                TreeProcs(1)%Quarks(TreeProcs(1)%NumQua)%Mom => q2neg(:)
          elseif( IsAScalar(TreeProcs(1)%PartType(TreeProcs(1)%NumPart)) ) then
                TreeProcs(1)%Scalars(TreeProcs(1)%NumSca)%Mom => q2neg(:)
          endif
          if( TreeProcs(2)%PartType(TreeProcs(2)%NumPart).eq.Glu_ ) then
                TreeProcs(2)%Gluons(TreeProcs(2)%NumGlu(0))%Mom => q3neg(:)
          elseif( IsAQuark(TreeProcs(2)%PartType(TreeProcs(2)%NumPart)) ) then
                TreeProcs(2)%Quarks(TreeProcs(2)%NumQua)%Mom => q3neg(:)
          elseif( IsAScalar(TreeProcs(2)%PartType(TreeProcs(2)%NumPart)) ) then
                TreeProcs(2)%Scalars(TreeProcs(2)%NumSca)%Mom => q3neg(:)
          endif
          if( TreeProcs(3)%PartType(TreeProcs(3)%NumPart).eq.Glu_ ) then
                TreeProcs(3)%Gluons(TreeProcs(3)%NumGlu(0))%Mom => q4neg(:)
          elseif( IsAQuark(TreeProcs(3)%PartType(TreeProcs(3)%NumPart)) ) then
                TreeProcs(3)%Quarks(TreeProcs(3)%NumQua)%Mom => q4neg(:)
          elseif( IsAScalar(TreeProcs(3)%PartType(TreeProcs(3)%NumPart)) ) then
                TreeProcs(3)%Scalars(TreeProcs(3)%NumSca)%Mom => q4neg(:)
          endif
          if( TreeProcs(4)%PartType(TreeProcs(4)%NumPart).eq.Glu_ ) then
                TreeProcs(4)%Gluons(TreeProcs(4)%NumGlu(0))%Mom => q5neg(:)
          elseif( IsAQuark(TreeProcs(4)%PartType(TreeProcs(4)%NumPart)) ) then
                TreeProcs(4)%Quarks(TreeProcs(4)%NumQua)%Mom => q5neg(:)
          elseif( IsAScalar(TreeProcs(4)%PartType(TreeProcs(4)%NumPart)) ) then
                TreeProcs(4)%Scalars(TreeProcs(4)%NumSca)%Mom => q5neg(:)
          endif
          if( TreeProcs(5)%PartType(TreeProcs(5)%NumPart).eq.Glu_ ) then
                TreeProcs(5)%Gluons(TreeProcs(5)%NumGlu(0))%Mom => q1neg(:)
          elseif( IsAQuark(TreeProcs(5)%PartType(TreeProcs(5)%NumPart)) ) then
                TreeProcs(5)%Quarks(TreeProcs(5)%NumQua)%Mom => q1neg(:)
          elseif( IsAScalar(TreeProcs(5)%PartType(TreeProcs(5)%NumPart)) ) then
                TreeProcs(5)%Scalars(TreeProcs(5)%NumSca)%Mom => q1neg(:)
          endif


          if (Ds.eq.5) then

          call givepol(TreeProcs(1)%PartType(1),q1(1:5),6,8,tag_pol,Nj1,BPOL1,POL1)
          call givepol(TreeProcs(2)%PartType(1),q2(1:5),6,8,tag_pol,Nj2,BPOL2,POL2)
          call givepol(TreeProcs(3)%PartType(1),q3(1:5),6,8,tag_pol,Nj3,BPOL3,POL3)
          call givepol(TreeProcs(4)%PartType(1),q4(1:5),6,8,tag_pol,Nj4,BPOL4,POL4)
          call givepol(TreeProcs(5)%PartType(1),q5(1:5),6,8,tag_pol,Nj5,BPOL5,POL5)

          res=dcmplx(0d0,0d0)

!DEC$ IF (_DebugUseMyAmps==0)
          if( any(Lab_in(l5c(1:5)).eq.'glu') ) then
            tag_f = 0
          elseif( all(Lab_in(l5c(1:5)).eq.'top') .or .all(Lab_in(l5c(1:5)).eq.'sto') .or. all(Lab_in(l5c(1:5)).eq.'str') ) then
            tag_f = 1
          elseif( all(Lab_in(l5c(1:5)).eq.'bot') .or. all(Lab_in(l5c(1:5)).eq.'chm') ) then
            tag_f = 2
          elseif( all(Lab_in(l5c(1:5)).eq.'sbo')  ) then
            tag_f = 3
          else
            tag_f = 99
          endif
          tag_Z=0

!------------ next loop
!!          corr. to arguments of amps in eq.(18) for ext. momenta
            i1=l5c(2)-l5c(1)
            i2=l5c(3)-l5c(2)
            i3=l5c(4)-l5c(3)
            i4=l5c(5)-l5c(4)
            lab1=Lab_in(l5c(1))
            lab2=Lab_in(l5c(2))
            call ampl(TreeProcs(1),6,8,Nj1,Nj2,POL1,q1(1:5),lab1,l5c(1),tag_f,1,i1,BPOL2,q2(1:5),lab2,mur1)

!----------next loop
            ia =i1+1
            ib = i1+i2
            lab1=Lab_in(l5c(2))
            lab2=Lab_in(l5c(3))
            call ampl(TreeProcs(2),6,8,Nj2,Nj3,POL2,q2(1:5),lab1,l5c(2),tag_f,ia,ib,BPOL3,q3(1:5),lab2,mur2)

!----------next loop
            ia=i1+i2+1
            ib=i1+i2+i3
            lab1=Lab_in(l5c(3))
            lab2=Lab_in(l5c(4))
            call ampl(TreeProcs(3),6,8,Nj3,Nj4,POL3,q3(1:5),lab1,l5c(3),tag_f,ia,ib,BPOL4,q4(1:5),lab2,mur3)

!----------next loop
            ia=i1+i2+i3+1
            ib=i1+i2+i3+i4
            lab1=Lab_in(l5c(4))
            lab2=Lab_in(l5c(5))
            call ampl(TreeProcs(4),6,8,Nj4,Nj5,POL4,q4(1:5),lab1,l5c(4),tag_f,ia,ib,BPOL5,q5(1:5),lab2,mur4)

!----------  next loop
            lab1=Lab_in(l5c(5))
            lab2=Lab_in(l5c(1))
            ia=i1+i2+i3+i4+1
            ib=Npoint
            call ampl(TreeProcs(5),6,8,Nj5,Nj1,POL5,q5(1:5),lab1,l5c(5),tag_f,ia,ib,BPOL1,q1(1:5),lab2,mur5)
!DEC$ ELSE

          if( any(Lab_in(l5c(1:5)).eq.'glu') ) then
            tag_f = 0
          elseif( all(Lab_in(l5c(1:5)).eq.'top') .or. all(Lab_in(l5c(1:5)).eq.'sto') .or. all(Lab_in(l5c(1:5)).eq.'str') ) then
            tag_f = 1
          elseif( all(Lab_in(l5c(1:5)).eq.'bot') .or. all(Lab_in(l5c(1:5)).eq.'chm') ) then
            tag_f = 2
          elseif( all(Lab_in(l5c(1:5)).eq.'sbo') ) then
            tag_f = 3
          else
            tag_f = 99
          endif

            call new_ampl(6,8,Nj1,Nj2,POL1,BPOL2,tag_f,tag_Z,TreeProcs(1),mur1)
            call new_ampl(6,8,Nj2,Nj3,POL2,BPOL3,tag_f,tag_Z,TreeProcs(2),mur2)
            call new_ampl(6,8,Nj3,Nj4,POL3,BPOL4,tag_f,tag_Z,TreeProcs(3),mur3)
            call new_ampl(6,8,Nj4,Nj5,POL4,BPOL5,tag_f,tag_Z,TreeProcs(4),mur4)
            call new_ampl(6,8,Nj5,Nj1,POL5,BPOL1,tag_f,tag_Z,TreeProcs(5),mur5)
!DEC$ ENDIF


!!       sum over internal polarizations (~1024 steps)
          do j1=1,Nj1
            do j2=1,Nj2
              do j3=1,Nj3
                do j4=1,Nj4
                   do j5=1,Nj5
                        res=res+ mur1(j1,j2)*mur2(j2,j3)*mur3(j3,j4)*mur4(j4,j5)*mur5(j5,j1)
                enddo
              enddo
            enddo
           enddo
          enddo

!         to multiply by proper power of I:
          do j=1,5
              if (Lab_in(l5c(j)).eq.'top' .or. Lab_in(l5c(j)).eq.'bot' .or. Lab_in(l5c(j)).eq.'sto' .or. Lab_in(l5c(j)).eq.'sbo' & 
             .or. Lab_in(l5c(j)).eq.'str' .or. Lab_in(l5c(j)).eq.'chm') then
                res = dcmplx(0d0,1d0)*res
              else
                 res = dcmplx(0d0,-1d0)*res
              endif
          enddo

!         intermediate 6-dim result
          res6 = res




!           8-dim calculation
          call givepol(TreeProcs(1)%PartType(1),q1(1:5),8,16,tag_pol,Nj1,BPOL1,POL1)
          call givepol(TreeProcs(2)%PartType(1),q2(1:5),8,16,tag_pol,Nj2,BPOL2,POL2)
          call givepol(TreeProcs(3)%PartType(1),q3(1:5),8,16,tag_pol,Nj3,BPOL3,POL3)
          call givepol(TreeProcs(4)%PartType(1),q4(1:5),8,16,tag_pol,Nj4,BPOL4,POL4)
          call givepol(TreeProcs(5)%PartType(1),q5(1:5),8,16,tag_pol,Nj5,BPOL5,POL5)
          res=dcmplx(0d0,0d0)



!DEC$ IF (_DebugUseMyAmps==0)
!------------ next loop
            lab1=Lab_in(l5c(1))
            lab2=Lab_in(l5c(2))
            call ampl(TreeProcs(1),8,16,Nj1,Nj2,POL1,q1(1:5),lab1,l5c(1),tag_f,1,i1,BPOL2,q2(1:5),lab2,mur1)

!------------ next loop
            ia =i1+1
            ib = i1+i2
            lab1=Lab_in(l5c(2))
            lab2=Lab_in(l5c(3))
            call ampl(TreeProcs(2),8,16,Nj2,Nj3,POL2,q2(1:5),lab1,l5c(2),tag_f,ia,ib,BPOL3,q3(1:5),lab2,mur2)

!----------- next loop
            ia=i1+i2+1
            ib=i1+i2+i3
            lab1=Lab_in(l5c(3))
            lab2=Lab_in(l5c(4))
            call ampl(TreeProcs(3),8,16,Nj3,Nj4,POL3,q3(1:5),lab1,l5c(3),tag_f,ia,ib,BPOL4,q4(1:5),lab2,mur3)

!--------- next loop
            ia=i1+i2+i3+1
            ib=i1+i2+i3+i4
            lab1=Lab_in(l5c(4))
            lab2=Lab_in(l5c(5))
            call ampl(TreeProcs(4),8,16,Nj4,Nj5,POL4,q4(1:5),lab1,l5c(4),tag_f,ia,ib,BPOL5,q5(1:5),lab2,mur4)

!------- next loop
            ia=i1+i2+i3+i4+1
            ib=Npoint
            lab1=Lab_in(l5c(5))
            lab2=Lab_in(l5c(1))
            call ampl(TreeProcs(5),8,16,Nj5,Nj1,POL5,q5(1:5),lab1,l5c(5),tag_f,ia,ib,BPOL1,q1(1:5),lab2,mur5)

!DEC$ ELSE
            call new_ampl(8,16,Nj1,Nj2,POL1,BPOL2,tag_f,tag_Z,TreeProcs(1),mur1)
            call new_ampl(8,16,Nj2,Nj3,POL2,BPOL3,tag_f,tag_Z,TreeProcs(2),mur2)
            call new_ampl(8,16,Nj3,Nj4,POL3,BPOL4,tag_f,tag_Z,TreeProcs(3),mur3)
            call new_ampl(8,16,Nj4,Nj5,POL4,BPOL5,tag_f,tag_Z,TreeProcs(4),mur4)
            call new_ampl(8,16,Nj5,Nj1,POL5,BPOL1,tag_f,tag_Z,TreeProcs(5),mur5)
!DEC$ ENDIF



          do j1=1,Nj1
            do j2=1,Nj2
              do j3=1,Nj3
                do j4=1,Nj4
                   do j5=1,Nj5
                     res=res + mur1(j1,j2)*mur2(j2,j3)*mur3(j3,j4)*mur4(j4,j5)*mur5(j5,j1)
                enddo
              enddo
            enddo
           enddo
          enddo

!           to multiply by proper power of I:
           do j=1,5
              if (Lab_in(l5c(j)).eq.'top' .or. Lab_in(l5c(j)).eq.'bot' .or. Lab_in(l5c(j)).eq.'sto' .or. Lab_in(l5c(j)).eq.'sbo' .or. & 
                  Lab_in(l5c(j)).eq.'str' .or. Lab_in(l5c(j)).eq.'chm') then
                res = dcmplx(0d0,1d0)*res
              else
                 res = dcmplx(0d0,-1d0)*res
              endif
           enddo

!          intermediate 8-dim result
           res8 = res
!           final result                          ! res=1/(D2-D1)*{ (D2-4)*(2/Nh6)*res6 - (D1-4)*(2/Nh8)*res8 }     with D1=6, D2=8
            if( tag_f.eq.0 ) then                 ! at least one gluon in loop,   Nh6=4, Nh8=2
               res = res6 -res8
            elseif (tag_f.eq.2) then              ! closed fermion loop,          Nh6=4, Nh8=8
               res = res6 - res8/4d0
            elseif (tag_f.eq.3) then              ! closed scalar loop,           Nh6=2, Nh8=2
               res = 2d0*res6 - res8
            else                                  ! no gluon in the loop,         Nh6=2, Nh8=2
               res = 2d0*res6 - res8
            endif
! endif for Ds = 5
         endif
END SUBROUTINE





       subroutine resid4_new(lv,k1,k2,k3,l4c,TreeProcs,res)
       use ModAmplitudes
       use ModProcess
       use ModMisc
       use ModParameters
       implicit none
       include 'misc/global_import'
       integer Ds,i,tag_pol
       integer j,j1,j2, tag_f,tag_Z
       integer i1,i2,i3,ia,ib,j3,j4
       double complex lv(5)
       double complex vpol1(16),vpol2(16)
       double complex k1(4),k2(4),k3(4)
       complex(8) :: q1(1:8),q2(1:8),q3(1:8),q4(1:8)
       complex(8),target :: q1neg(1:8),q2neg(1:8),q3neg(1:8),q4neg(1:8)
       complex(8) BPOL1(8,16),POL1(8,16)
       complex(8) BPOL2(8,16),POL2(8,16)
       complex(8) BPOL3(8,16),POL3(8,16)
       complex(8) BPOL4(8,16),POL4(8,16)
       complex(8) mur1(10,10),mur2(10,10)
       complex(8) mur3(10,10),mur4(10,10)
       double complex res,res0,res6,res8,propX
       double complex krefa(4),lvt(5)
       double complex vprop(5)
       integer Nj1,Nj2,Nj3,Nj4
       character lab1*3, lab2*3
       integer l4c(4),n45,lmatch45(50),im,pos,pos1
       type(TreeProcess) :: TreeProcs(:)

          tag_pol = 0


          if( lv(5).eq.0d0 ) then
               Ds = 4
          else
               Ds = 5
          endif


          do i=1,4
          q1(i)=lv(i)
          q2(i)=lv(i)+k1(i)
          q3(i)=lv(i)+k2(i)
          q4(i)=lv(i)+k3(i)
          enddo

          if (Ds.eq.4) then
          q1(5)=dcmplx(0d0,0d0)
          q2(5)=dcmplx(0d0,0d0)
          q3(5)=dcmplx(0d0,0d0)
          q4(5)=dcmplx(0d0,0d0)
          else
          q1(5)=lv(5)
          q2(5)=lv(5)
          q3(5)=lv(5)
          q4(5)=lv(5)
          endif
          q1(6:8)=(0d0,0d0)
          q2(6:8)=(0d0,0d0)
          q3(6:8)=(0d0,0d0)
          q4(6:8)=(0d0,0d0)
          q1neg(1:8)=-q1(1:8)
          q2neg(1:8)=-q2(1:8)
          q3neg(1:8)=-q3(1:8)
          q4neg(1:8)=-q4(1:8)

          if( any(Lab_in(l4c(1:4)).eq.'glu') ) then
            tag_f = 0
          elseif( all(Lab_in(l4c(1:4)).eq.'top') .or. all(Lab_in(l4c(1:4)).eq.'sto') .or. all(Lab_in(l4c(1:4)).eq.'str') ) then
            tag_f = 1
          elseif( all(Lab_in(l4c(1:4)).eq.'bot') .or. all(Lab_in(l4c(1:4)).eq.'chm')) then
            tag_f = 2
          elseif( all(Lab_in(l4c(1:4)).eq.'sbo') .or. all(Lab_in(l4c(1:4)).eq.'sch')) then
            tag_f = 3
          else
            tag_f = 99
          endif
          tag_Z=0


!         set momentum vector for last particle
          if( TreeProcs(1)%PartType(TreeProcs(1)%NumPart).eq.Glu_ ) then
                TreeProcs(1)%Gluons(TreeProcs(1)%NumGlu(0))%Mom => q2neg(:)
          elseif( IsAQuark(TreeProcs(1)%PartType(TreeProcs(1)%NumPart)) ) then
                TreeProcs(1)%Quarks(TreeProcs(1)%NumQua)%Mom => q2neg(:)
          elseif( IsAScalar(TreeProcs(1)%PartType(TreeProcs(1)%NumPart)) ) then
                TreeProcs(1)%Scalars(TreeProcs(1)%NumSca)%Mom => q2neg(:)
          endif
          if( TreeProcs(2)%PartType(TreeProcs(2)%NumPart).eq.Glu_ ) then
                TreeProcs(2)%Gluons(TreeProcs(2)%NumGlu(0))%Mom => q3neg(:)
          elseif( IsAQuark(TreeProcs(2)%PartType(TreeProcs(2)%NumPart)) ) then
                TreeProcs(2)%Quarks(TreeProcs(2)%NumQua)%Mom => q3neg(:)
          elseif( IsAScalar(TreeProcs(2)%PartType(TreeProcs(2)%NumPart)) ) then
                TreeProcs(2)%Scalars(TreeProcs(2)%NumSca)%Mom => q3neg(:)
          endif
          if( TreeProcs(3)%PartType(TreeProcs(3)%NumPart).eq.Glu_ ) then
                TreeProcs(3)%Gluons(TreeProcs(3)%NumGlu(0))%Mom => q4neg(:)
          elseif( IsAQuark(TreeProcs(3)%PartType(TreeProcs(3)%NumPart)) ) then
                TreeProcs(3)%Quarks(TreeProcs(3)%NumQua)%Mom => q4neg(:)
          elseif( IsAScalar(TreeProcs(3)%PartType(TreeProcs(3)%NumPart)) ) then
                TreeProcs(3)%Scalars(TreeProcs(3)%NumSca)%Mom => q4neg(:)
          endif
          if( TreeProcs(4)%PartType(TreeProcs(4)%NumPart).eq.Glu_ ) then
                TreeProcs(4)%Gluons(TreeProcs(4)%NumGlu(0))%Mom => q1neg(:)
          elseif( IsAQuark(TreeProcs(4)%PartType(TreeProcs(4)%NumPart)) ) then
                TreeProcs(4)%Quarks(TreeProcs(4)%NumQua)%Mom => q1neg(:)
          elseif( IsAScalar(TreeProcs(4)%PartType(TreeProcs(4)%NumPart)) ) then
                TreeProcs(4)%Scalars(TreeProcs(4)%NumSca)%Mom => q1neg(:)
          endif


          if (Ds.eq.4) then
          call givepol(TreeProcs(1)%PartType(1),q1(1:5),4,4,tag_pol,Nj1,BPOL1,POL1)        ! POL1/BPOL1=u-/ubar-spinor for fermions,   POL1==BPOL1 for bosons
          call givepol(TreeProcs(2)%PartType(1),q2(1:5),4,4,tag_pol,Nj2,BPOL2,POL2)
          call givepol(TreeProcs(3)%PartType(1),q3(1:5),4,4,tag_pol,Nj3,BPOL3,POL3)
          call givepol(TreeProcs(4)%PartType(1),q4(1:5),4,4,tag_pol,Nj4,BPOL4,POL4)

!         calculate product of amplitudes
!         first  organize 4 lists that will be arguments
!         for amplitude procedure
!         start with the first label on the cut list ll
!         Lab_ex and Lab_in

            res=dcmplx(0d0,0d0)


!DEC$ IF (_DebugUseMyAmps==0)
            i1=l4c(2)-l4c(1)
            i2=l4c(3)-l4c(2)
            i3=l4c(4)-l4c(3)
!----------
            lab1=Lab_in(l4c(1))
            lab2=Lab_in(l4c(2))
            call ampl(TreeProcs(1),4,4,Nj1,Nj2,POL1,q1(1:5),lab1,l4c(1),tag_f,1,i1,BPOL2,q2(1:5),lab2,mur1)
!-----------
            ia =i1+1
            ib = i1+i2
            lab1=Lab_in(l4c(2))
            lab2=Lab_in(l4c(3))
            call ampl(TreeProcs(2),4,4,Nj2,Nj3,POL2,q2(1:5),lab1,l4c(2),tag_f,ia,ib,BPOL3,q3(1:5),lab2,mur2)
!----------
            ia=i1+i2+1
            ib=i1+i2+i3
            lab1=Lab_in(l4c(3))
            lab2=Lab_in(l4c(4))
            call ampl(TreeProcs(3),4,4,Nj3,Nj4,POL3,q3(1:5),lab1,l4c(3),tag_f,ia,ib,BPOL4,q4(1:5),lab2,mur3)
!-----------
            ia=i1+i2+i3+1
            ib=Npoint
            lab1=Lab_in(l4c(4))
            lab2=Lab_in(l4c(1))
            call ampl(TreeProcs(4),4,4,Nj4,Nj1,POL4,q4(1:5),lab1,l4c(4),tag_f,ia,ib,BPOL1,q1(1:5),lab2,mur4)
!DEC$ ELSE
            call new_ampl(4,4,Nj1,Nj2,POL1,BPOL2,tag_f,tag_Z,TreeProcs(1),mur1)
            call new_ampl(4,4,Nj2,Nj3,POL2,BPOL3,tag_f,tag_Z,TreeProcs(2),mur2)
            call new_ampl(4,4,Nj3,Nj4,POL3,BPOL4,tag_f,tag_Z,TreeProcs(3),mur3)
            call new_ampl(4,4,Nj4,Nj1,POL4,BPOL1,tag_f,tag_Z,TreeProcs(4),mur4)
!DEC$ ENDIF


            do j1=1,Nj1
               do j2=1,Nj2
               do j3=1,Nj3
                  do j4=1,Nj4
!                     print *, j1,j2,j3,j4
!                     print *,mur1(j1,j2),mur2(j2,j3),mur3(j3,j4),mur4(j4,j1)
!                     print *, 'add',mur1(j1,j2)*mur2(j2,j3)*mur3(j3,j4)*mur4(j4,j1)
                        res=res+ mur1(j1,j2)*mur2(j2,j3)*mur3(j3,j4)*mur4(j4,j1)
                  enddo
               enddo
               enddo
            enddo

!           to multiply by proper power of I:
           do j=1,4
              if (Lab_in(l4c(j)).eq.'top' .or. Lab_in(l4c(j)).eq.'bot' .or. Lab_in(l4c(j)).eq.'sto' .or. Lab_in(l4c(j)).eq.'sbo' .or. Lab_in(l4c(j)).eq.'str' .or. Lab_in(l4c(j)).eq.'chm') then
                 res = dcmplx(0d0,1d0)*res
              else
                 res = dcmplx(0d0,-1d0)*res
              endif
           enddo


!cccccccccccccccccccc endif for Ds = 4
         endif

          if (Ds.eq.5) then

          call givepol(TreeProcs(1)%PartType(1),q1(1:5),6,8,tag_pol,Nj1,BPOL1,POL1)
          call givepol(TreeProcs(2)%PartType(1),q2(1:5),6,8,tag_pol,Nj2,BPOL2,POL2)
          call givepol(TreeProcs(3)%PartType(1),q3(1:5),6,8,tag_pol,Nj3,BPOL3,POL3)
          call givepol(TreeProcs(4)%PartType(1),q4(1:5),6,8,tag_pol,Nj4,BPOL4,POL4)


!         calculate product of amplitudes
!         first  organize 4 lists that will be arguments
!         for amplitude procedure
!         start with the first label on the cut list ll
!         Lab_ex and Lab_in

          res=dcmplx(0d0,0d0)

            i1=l4c(2)-l4c(1)
            i2=l4c(3)-l4c(2)
            i3=l4c(4)-l4c(3)

!DEC$ IF (_DebugUseMyAmps==0)
!---------
            lab1=Lab_in(l4c(1))
            lab2=Lab_in(l4c(2))
            call ampl(TreeProcs(1),6,8,Nj1,Nj2,POL1,q1(1:5),lab1,l4c(1),tag_f,1,i1,BPOL2,q2(1:5),lab2,mur1)
!-----------
            ia =i1+1
            ib = i1+i2
            lab1=Lab_in(l4c(2))
            lab2=Lab_in(l4c(3))
            call ampl(TreeProcs(2),6,8,Nj2,Nj3,POL2,q2(1:5),lab1,l4c(2),tag_f,ia,ib,BPOL3,q3(1:5),lab2,mur2)
!----------
            ia=i1+i2+1
            ib=i1+i2+i3
            lab1=Lab_in(l4c(3))
            lab2=Lab_in(l4c(4))
            call ampl(TreeProcs(3),6,8,Nj3,Nj4,POL3,q3(1:5),lab1,l4c(3),tag_f,ia,ib,BPOL4,q4(1:5),lab2,mur3)
!----------
            ia=i1+i2+i3+1
            ib=Npoint
            lab1=Lab_in(l4c(4))
            lab2=Lab_in(l4c(1))
            call ampl(TreeProcs(4),6,8,Nj4,Nj1,POL4,q4(1:5),lab1,l4c(4),tag_f,ia,ib,BPOL1,q1(1:5),lab2,mur4)
!DEC$ ELSE
            call new_ampl(6,8,Nj1,Nj2,POL1,BPOL2,tag_f,tag_Z,TreeProcs(1),mur1)
            call new_ampl(6,8,Nj2,Nj3,POL2,BPOL3,tag_f,tag_Z,TreeProcs(2),mur2)
            call new_ampl(6,8,Nj3,Nj4,POL3,BPOL4,tag_f,tag_Z,TreeProcs(3),mur3)
            call new_ampl(6,8,Nj4,Nj1,POL4,BPOL1,tag_f,tag_Z,TreeProcs(4),mur4)
!DEC$ ENDIF


            do j1=1,Nj1
               do j2=1,Nj2
               do j3=1,Nj3
                  do j4=1,Nj4
                        res=res+ mur1(j1,j2)*mur2(j2,j3)*mur3(j3,j4)*mur4(j4,j1)
                  enddo
               enddo
               enddo
            enddo

!           to multiply by proper power of I:
           do j=1,4
              if (Lab_in(l4c(j)).eq.'top' .or. Lab_in(l4c(j)).eq.'bot' .or. Lab_in(l4c(j)).eq.'sto' .or. Lab_in(l4c(j)).eq.'sbo' .or. Lab_in(l4c(j)).eq.'str' .or. Lab_in(l4c(j)).eq.'chm') then
                res = dcmplx(0d0,1d0)*res
              else
                 res = dcmplx(0d0,-1d0)*res
              endif
           enddo

!           intermediate 6-dim result
            res6 = res

!         8-dim calculation
          call givepol(TreeProcs(1)%PartType(1),q1(1:5),8,16,tag_pol,Nj1,BPOL1,POL1)
          call givepol(TreeProcs(2)%PartType(1),q2(1:5),8,16,tag_pol,Nj2,BPOL2,POL2)
          call givepol(TreeProcs(3)%PartType(1),q3(1:5),8,16,tag_pol,Nj3,BPOL3,POL3)
          call givepol(TreeProcs(4)%PartType(1),q4(1:5),8,16,tag_pol,Nj4,BPOL4,POL4)


!         calculate product of amplitudes
!         first  organize 4 lists that will be arguments
!         for amplitude procedure
!         start with the first label on the cut list ll
!         Lab_ex and Lab_in

          res=dcmplx(0d0,0d0)

!DEC$ IF (_DebugUseMyAmps==0)
          lab1=Lab_in(l4c(1))
          lab2=Lab_in(l4c(2))
          call ampl(TreeProcs(1),8,16,Nj1,Nj2,POL1,q1(1:5),lab1,l4c(1),tag_f,1,i1,BPOL2,q2(1:5),lab2,mur1)
!---------
            ia =i1+1
            ib = i1+i2
            lab1=Lab_in(l4c(2))
            lab2=Lab_in(l4c(3))
            call ampl(TreeProcs(2),8,16,Nj2,Nj3,POL2,q2(1:5),lab1,l4c(2),tag_f,ia,ib,BPOL3,q3(1:5),lab2,mur2)
!----------
            ia=i1+i2+1
            ib=i1+i2+i3
            lab1=Lab_in(l4c(3))
            lab2=Lab_in(l4c(4))
            call ampl(TreeProcs(3),8,16,Nj3,Nj4,POL3,q3(1:5),lab1,l4c(3),tag_f,ia,ib,BPOL4,q4(1:5),lab2,mur3)
!---------
            ia=i1+i2+i3+1
            ib=Npoint
            lab1=Lab_in(l4c(4))
            lab2=Lab_in(l4c(1))
            call ampl(TreeProcs(4),8,16,Nj4,Nj1,POL4,q4(1:5),lab1,l4c(4),tag_f,ia,ib,BPOL1,q1(1:5),lab2,mur4)
!DEC$ ELSE
            call new_ampl(8,16,Nj1,Nj2,POL1,BPOL2,tag_f,tag_Z,TreeProcs(1),mur1)
            call new_ampl(8,16,Nj2,Nj3,POL2,BPOL3,tag_f,tag_Z,TreeProcs(2),mur2)
            call new_ampl(8,16,Nj3,Nj4,POL3,BPOL4,tag_f,tag_Z,TreeProcs(3),mur3)
            call new_ampl(8,16,Nj4,Nj1,POL4,BPOL1,tag_f,tag_Z,TreeProcs(4),mur4)
!DEC$ ENDIF

            do j1=1,Nj1
               do j2=1,Nj2
               do j3=1,Nj3
                  do j4=1,Nj4
                        res=res+ mur1(j1,j2)*mur2(j2,j3)*mur3(j3,j4)*mur4(j4,j1)
                  enddo
               enddo
               enddo
            enddo

!           to multiply by proper power of I:
           do j=1,4
              if (Lab_in(l4c(j)).eq.'top' .or. Lab_in(l4c(j)).eq.'bot' .or. Lab_in(l4c(j)).eq.'sto' .or. Lab_in(l4c(j)).eq.'sbo' .or. Lab_in(l4c(j)).eq.'str' .or. Lab_in(l4c(j)).eq.'chm') then
                res = dcmplx(0d0,1d0)*res
              else
                 res = dcmplx(0d0,-1d0)*res
              endif
           enddo
           res8 = res

!           final result
            if( tag_f.eq.0 ) then
               res = res6 -res8
            elseif (tag_f.eq.2) then
               res = res6 - res8/4d0
            elseif (tag_f.eq.3) then
               res = 2d0*res6 - res8
            else
               res = 2d0*res6 - res8
            endif

!cccccccccccccccccccc endif for Ds = 5
         endif


       return
       end subroutine








       subroutine resid3_new(lv,k1,k2,l3c,TreeProcs,res)
       use ModAmplitudes
       use ModProcess
       use ModMisc
       use ModParameters
       implicit none
       include 'misc/global_import'
       integer Ds,i,tag_pol,pos,pos1,pos2
       integer j,j1,j2,tag_f,tag_Z
       integer i1,i2,ia,ib,j3,n34,lmatch34(50)
       integer lmatch35(50)
       double complex lv(5),lvt(5)
       double complex v45(4),vprop(5),propX,r1,r2,krefa(4)
       double complex vpol1(16),vpol2(16)
       double complex k1(4),k2(4)
       complex(8), target :: q1(1:8),q2(1:8),q3(1:8),q1neg(1:8),q2neg(1:8),q3neg(1:8)
       complex(8) BPOL1(8,16),POL1(8,16)
       complex(8) BPOL2(8,16),POL2(8,16)
       complex(8) BPOL3(8,16),POL3(8,16)
       complex(8) mur1(10,10),mur2(10,10)
       complex(8) mur3(10,10)
       double complex e(1),f(5)
       double complex vne(5)
       double complex res,res0,res6,res8
       double complex tpol1(16),tpol2(16)
       double complex prop1,prop2
       integer Nj1,Nj2,Nj3,n35,im,lpos1(2)
       character lab1*3, lab2*3
       integer l3c(3)
       type(TreeProcess) :: TreeProcs(:)

          tag_pol=0

          if( lv(5).eq.0d0 ) then
               Ds = 4
          else
               Ds = 5
          endif

          do i=1,4
            q1(i)=lv(i)
            q2(i)=lv(i)+k1(i)
            q3(i)=lv(i)+k2(i)
          enddo

          if (Ds.eq.4) then
            q1(5)=dcmplx(0d0,0d0)
            q2(5)=dcmplx(0d0,0d0)
            q3(5)=dcmplx(0d0,0d0)
          else
            q1(5)=lv(5)
            q2(5)=lv(5)
            q3(5)=lv(5)
          endif
          q1(6:8)=(0d0,0d0)
          q2(6:8)=(0d0,0d0)
          q3(6:8)=(0d0,0d0)
          q1neg(1:8) = -q1(1:8)
          q2neg(1:8) = -q2(1:8)
          q3neg(1:8) = -q3(1:8)

          if( any(Lab_in(l3c(1:3)).eq.'glu') ) then
            tag_f = 0
          elseif( all(Lab_in(l3c(1:3)).eq.'top') .or. all(Lab_in(l3c(1:3)).eq.'sto') .or. all(Lab_in(l3c(1:3)).eq.'str') ) then
            tag_f = 1
          elseif( all(Lab_in(l3c(1:3)).eq.'bot') .or. all(Lab_in(l3c(1:3)).eq.'chm') ) then
            tag_f = 2
          elseif( all(Lab_in(l3c(1:3)).eq.'sbo')  ) then
            tag_f = 3
          else
            tag_f = 99
          endif
          tag_Z=0

!     amplitudes for different space dimensions
!         set momentum vector for last particle
          if( TreeProcs(1)%PartType(TreeProcs(1)%NumPart).eq.Glu_ ) then
                TreeProcs(1)%Gluons(TreeProcs(1)%NumGlu(0))%Mom => q2neg(:)
          elseif( IsAQuark(TreeProcs(1)%PartType(TreeProcs(1)%NumPart)) ) then
                TreeProcs(1)%Quarks(TreeProcs(1)%NumQua)%Mom => q2neg(:)
          elseif( IsAScalar(TreeProcs(1)%PartType(TreeProcs(1)%NumPart)) ) then
                TreeProcs(1)%Scalars(TreeProcs(1)%NumSca)%Mom => q2neg(:)
          endif
          if( TreeProcs(2)%PartType(TreeProcs(2)%NumPart).eq.Glu_ ) then
                TreeProcs(2)%Gluons(TreeProcs(2)%NumGlu(0))%Mom => q3neg(:)
          elseif( IsAQuark(TreeProcs(2)%PartType(TreeProcs(2)%NumPart)) ) then
                TreeProcs(2)%Quarks(TreeProcs(2)%NumQua)%Mom => q3neg(:)
          elseif( IsAScalar(TreeProcs(2)%PartType(TreeProcs(2)%NumPart)) ) then
                TreeProcs(2)%Scalars(TreeProcs(2)%NumSca)%Mom => q3neg(:)
          endif
          if( TreeProcs(3)%PartType(TreeProcs(3)%NumPart).eq.Glu_ ) then
                TreeProcs(3)%Gluons(TreeProcs(3)%NumGlu(0))%Mom => q1neg(:)
          elseif( IsAQuark(TreeProcs(3)%PartType(TreeProcs(3)%NumPart)) ) then
                TreeProcs(3)%Quarks(TreeProcs(3)%NumQua)%Mom => q1neg(:)
          elseif( IsAScalar(TreeProcs(3)%PartType(TreeProcs(3)%NumPart)) ) then
                TreeProcs(3)%Scalars(TreeProcs(3)%NumSca)%Mom => q1neg(:)
          endif


          if (Ds.eq.4) then
            call givepol(TreeProcs(1)%PartType(1),q1(1:5),4,4,tag_pol,Nj1,BPOL1,POL1)
            call givepol(TreeProcs(2)%PartType(1),q2(1:5),4,4,tag_pol,Nj2,BPOL2,POL2)
            call givepol(TreeProcs(3)%PartType(1),q3(1:5),4,4,tag_pol,Nj3,BPOL3,POL3)
            res=dcmplx(0d0,0d0)

!DEC$ IF (_DebugUseMyAmps==0)
            i1=l3c(2)-l3c(1)
            i2=l3c(3)-l3c(2)
            lab1=Lab_in(l3c(1))
            lab2=Lab_in(l3c(2))
            call ampl(TreeProcs(1),4,4,Nj1,Nj2,POL1,q1(1:5),lab1,l3c(1),tag_f,1,i1,BPOL2,q2(1:5),lab2,mur1)
!----------
            lab1=Lab_in(l3c(2))
            lab2=Lab_in(l3c(3))
            ia =i1+1
            ib = i1+i2
            call ampl(TreeProcs(2),4,4,Nj2,Nj3,POL2,q2(1:5),lab1,l3c(2),tag_f,ia,ib,BPOL3,q3(1:5),lab2,mur2)
!------------
            lab1=Lab_in(l3c(3))
            lab2=Lab_in(l3c(1))
            ia=i1+i2+1
            ib=Npoint
            call ampl(TreeProcs(3),4,4,Nj3,Nj1,POL3,q3(1:5),lab1,l3c(3),tag_f,ia,ib,BPOL1,q1(1:5),lab2,mur3)
!DEC$ ELSE
            call new_ampl(4,4,Nj1,Nj2,POL1,BPOL2,tag_f,tag_Z,TreeProcs(1),mur1)
            call new_ampl(4,4,Nj2,Nj3,POL2,BPOL3,tag_f,tag_Z,TreeProcs(2),mur2)
            call new_ampl(4,4,Nj3,Nj1,POL3,BPOL1,tag_f,tag_Z,TreeProcs(3),mur3)
!DEC$ ENDIF

            do j1=1,Nj1
            do j2=1,Nj2
            do j3=1,Nj3
                    res=res+ mur1(j1,j2)*mur2(j2,j3)*mur3(j3,j1)
            enddo
            enddo
            enddo

!           to multiply by proper power of I:
           do j=1,3
              if (Lab_in(l3c(j)).eq.'top' .or. Lab_in(l3c(j)).eq.'bot' .or. Lab_in(l3c(j)).eq.'sto' .or. Lab_in(l3c(j)).eq.'sbo' .or. Lab_in(l3c(j)).eq.'str' .or. Lab_in(l3c(j)).eq.'chm') then
                res = dcmplx(0d0,1d0)*res
              else
                 res = dcmplx(0d0,-1d0)*res
              endif
           enddo


!cccccccccccccccccccc endif for Ds = 4
         endif



          if (Ds.eq.5) then

          call givepol(TreeProcs(1)%PartType(1),q1(1:5),6,8,tag_pol,Nj1,BPOL1,POL1)
          call givepol(TreeProcs(2)%PartType(1),q2(1:5),6,8,tag_pol,Nj2,BPOL2,POL2)
          call givepol(TreeProcs(3)%PartType(1),q3(1:5),6,8,tag_pol,Nj3,BPOL3,POL3)

!         calculate product of amplitudes
!         first  organize 4 lists that will be arguments
!         for amplitude procedure
!         start with the first label on the cut list ll
!         Lab_ex and Lab_in

          res=dcmplx(0d0,0d0)

!DEC$ IF (_DebugUseMyAmps==0)
            i1=l3c(2)-l3c(1)
            i2=l3c(3)-l3c(2)
            lab1=Lab_in(l3c(1))
            lab2=Lab_in(l3c(2))
            call ampl(TreeProcs(1),6,8,Nj1,Nj2,POL1,q1(1:5),lab1,l3c(1),tag_f,1,i1,BPOL2,q2(1:5),lab2,mur1)
!-----------
            ia =i1+1
            ib = i1+i2
            lab1=Lab_in(l3c(2))
            lab2=Lab_in(l3c(3))
            call ampl(TreeProcs(2),6,8,Nj2,Nj3,POL2,q2(1:5),lab1,l3c(2),tag_f,ia,ib,BPOL3,q3(1:5),lab2,mur2)
!----------
            lab1=Lab_in(l3c(3))
            lab2=Lab_in(l3c(1))
            ia=i1+i2+1
            ib=Npoint
            call ampl(TreeProcs(3),6,8,Nj3,Nj1,POL3,q3(1:5),lab1,l3c(3),tag_f,ia,ib,BPOL1,q1(1:5),lab2,mur3)
!DEC$ ELSE
            call new_ampl(6,8,Nj1,Nj2,POL1,BPOL2,tag_f,tag_Z,TreeProcs(1),mur1)
            call new_ampl(6,8,Nj2,Nj3,POL2,BPOL3,tag_f,tag_Z,TreeProcs(2),mur2)
            call new_ampl(6,8,Nj3,Nj1,POL3,BPOL1,tag_f,tag_Z,TreeProcs(3),mur3)
!DEC$ ENDIF

            do j1=1,Nj1
            do j2=1,Nj2
            do j3=1,Nj3
                  res=res+ mur1(j1,j2)*mur2(j2,j3)*mur3(j3,j1)
            enddo
            enddo
            enddo


!           to multiply by proper power of I:
           do j=1,3
              if (Lab_in(l3c(j)).eq.'top' .or. Lab_in(l3c(j)).eq.'bot' .or. Lab_in(l3c(j)).eq.'sto' .or. Lab_in(l3c(j)).eq.'sbo' .or. Lab_in(l3c(j)).eq.'str' .or. Lab_in(l3c(j)).eq.'chm') then
                 res = dcmplx(0d0,1d0)*res
                else
                 res = dcmplx(0d0,-1d0)*res
              endif
           enddo

!           intermediate 6-dim result
            res6 = res



!           8-dim calculation
          call givepol(TreeProcs(1)%PartType(1),q1(1:5),8,16,tag_pol,Nj1,BPOL1,POL1)
          call givepol(TreeProcs(2)%PartType(1),q2(1:5),8,16,tag_pol,Nj2,BPOL2,POL2)
          call givepol(TreeProcs(3)%PartType(1),q3(1:5),8,16,tag_pol,Nj3,BPOL3,POL3)


!         calculate product of amplitudes
!         first  organize 4 lists that will be arguments
!         for amplitude procedure
!         start with the first label on the cut list ll
!         Lab_ex and Lab_in

          res=dcmplx(0d0,0d0)

!DEC$ IF (_DebugUseMyAmps==0)
            i1=l3c(2)-l3c(1)
            i2=l3c(3)-l3c(2)
            lab1=Lab_in(l3c(1))
            lab2=Lab_in(l3c(2))
            call ampl(TreeProcs(1),8,16,Nj1,Nj2,POL1,q1(1:5),lab1,l3c(1),tag_f,1,i1,BPOL2,q2(1:5),lab2,mur1)
!------------new loop
            ia =i1+1
            ib = i1+i2
            lab1=Lab_in(l3c(2))
            lab2=Lab_in(l3c(3))
            call ampl(TreeProcs(2),8,16,Nj2,Nj3,POL2,q2(1:5),lab1,l3c(2),tag_f,ia,ib,BPOL3,q3(1:5),lab2,mur2)
!----------
            ia=i1+i2+1
            ib=Npoint
            lab1=Lab_in(l3c(3))
            lab2=Lab_in(l3c(1))
            call ampl(TreeProcs(3),8,16,Nj3,Nj1,POL3,q3(1:5),lab1,l3c(3),tag_f,ia,ib,BPOL1,q1(1:5),lab2,mur3)
!DEC$ ELSE
            call new_ampl(8,16,Nj1,Nj2,POL1,BPOL2,tag_f,tag_Z,TreeProcs(1),mur1)
            call new_ampl(8,16,Nj2,Nj3,POL2,BPOL3,tag_f,tag_Z,TreeProcs(2),mur2)
            call new_ampl(8,16,Nj3,Nj1,POL3,BPOL1,tag_f,tag_Z,TreeProcs(3),mur3)
!DEC$ ENDIF


            do j1=1,Nj1
            do j2=1,Nj2
            do j3=1,Nj3
                  res=res+ mur1(j1,j2)*mur2(j2,j3)*mur3(j3,j1)
            enddo
            enddo
            enddo


!           to multiply by proper power of I:
           do j=1,3
              if (Lab_in(l3c(j)).eq.'top' .or. Lab_in(l3c(j)).eq.'bot' .or. Lab_in(l3c(j)).eq.'sto' .or. Lab_in(l3c(j)).eq.'sbo' .or. Lab_in(l3c(j)).eq.'str' .or. Lab_in(l3c(j)).eq.'chm') then
                res = dcmplx(0d0,1d0)*res
              else
                res = dcmplx(0d0,-1d0)*res
              endif
           enddo
           res8 = res

!           final result
            if( tag_f.eq.0 ) then
               res = res6 - res8
            elseif (tag_f.eq.2) then
               res = res6 - res8/4d0
            elseif (tag_f.eq.3) then
               res = 2d0*res6 - res8
            else
               res = 2d0*res6 - res8
            endif

!cccccccccccccccccccc endif for Ds = 5
         endif



       return
       end subroutine





       subroutine resid2_new(lv,k1,l2c,TreeProcs,res)
       use ModAmplitudes
       use ModProcess
       use ModMisc
       use ModParameters
       implicit none
       include 'misc/global_import'
       integer Ds,i,tag_pol,pos,pos1,pos2,n25
       integer j,j1,j2,lpos(2),lpos2(3),pos3
       integer i1,ia,ib,n24,n23
       integer lmatch25(50),lmatch24(50),lmatch23(50)
       integer tag_f,tag_Z
       double complex lv(5),lvt(5),e(1)
       double complex v45(4),vprop(5),propX,r1,r2,krefa(4)
       double complex d(5),c(10),re
       double complex prop1,prop2
       double complex vpol1(16),vpol2(16)
       double complex k1(4)
       double complex r3,r4
       double complex vne(5),trikoeff
       double complex v3(4),v4(4)
       complex(8) q1(1:8),q2(1:8)
       complex(8),target :: q1neg(1:8),q2neg(1:8)
       complex(8) BPOL1(8,16),POL1(8,16)
       complex(8) BPOL2(8,16),POL2(8,16)
       complex(8) mur1(10,10),mur2(10,10)
       double complex res,res0,res6,res8
       double complex q12(5),prop3,r22,r32,r42
       integer Nj1,Nj2
       character lab1*3, lab2*3
       integer l2c(2),im
       type(TreeProcess) :: TreeProcs(:)
!        complex(16) :: res_Impr


          if( lv(5).eq.0d0 ) then
               Ds = 4
          else
               Ds = 5
          endif

          do i=1,4
            q1(i)=lv(i)
            q2(i)=lv(i)+k1(i)
          enddo

          if (Ds.eq.4) then
            q1(5)=dcmplx(0d0,0d0)
            q2(5)=dcmplx(0d0,0d0)
          else
            q1(5)=lv(5)
            q2(5)=lv(5)
          endif
          q1(6:8)=(0d0,0d0)
          q2(6:8)=(0d0,0d0)
          q1neg(1:8)=-q1(1:8)
          q2neg(1:8)=-q2(1:8)

          q12(1:5) = q1(1:5)+q2(1:5)

          call sc(5,q12,q12,r1)
          if(zabs(r1-dcmplx(ExtParticle(1)%Mass2)).lt.1D-6) then   !  r1=(q1+q2)^2 = q1^2+q2^2+2*q1.q2 = 0+mt**2+0    because lv^2=lv.k1=0
              tag_pol=1                                    !  this checks if there is a self-energy contribution
          else                                             !  if yes (i.e. tag_pol=1) then use full gluon polarization sum
              tag_pol=0                                    !  note: tag_pol does not affect fermion loops!
          endif

          if( any(Lab_in(l2c(1:2)).eq.'glu') ) then
            tag_f = 0
          elseif( all(Lab_in(l2c(1:2)).eq.'top') .or.  all(Lab_in(l2c(1:2)).eq.'sto') .or. all(Lab_in(l2c(1:2)).eq.'str') ) then
            tag_f = 1
          elseif( all(Lab_in(l2c(1:2)).eq.'bot') .or. all(Lab_in(l2c(1:2)).eq.'chm') ) then
            tag_f = 2
          elseif( all(Lab_in(l2c(1:2)).eq.'sbo') ) then
            tag_f = 3
          else
            tag_f = 99
          endif

!          if ( (tag_f .eq. 1 .or. tag_f .eq. 2) .and. (Lab_ex(l2c(1)) .eq. 'zee') .and. (Lab_ex(l2c(2)-1) .eq. 'zee') ) then
!         tag_Z=1 removes diagrams of a Z coupl. to a closed fermion loop coupl. to the rest 
          if ( (TreeProcs(1)%NumPart .eq. 3) .and. ( abs(TreeProcs(1)%PartType(1)) .eq. Bot_ .or. abs(TreeProcs(1)%PartType(1)) .eq. Chm_ )&
               & .and. ( abs(TreeProcs(1)%PartType(3)) .eq. Bot_ .or. abs(TreeProcs(1)%PartType(3)) .eq. Chm_) &
               & .and. ( TreeProcs(1)%PartType(2) .eq. Z0_) ) then
             tag_Z=1
!              print *, "setting tag_Z=1 (2-cut)"
          elseif ( (TreeProcs(2)%NumPart .eq. 3) .and. ( abs(TreeProcs(2)%PartType(1)) .eq. Bot_ .or. abs(TreeProcs(2)%PartType(1)) .eq. Chm_ )&
               & .and. ( abs(TreeProcs(2)%PartType(3)) .eq. Bot_ .or. abs(TreeProcs(2)%PartType(3)) .eq. Chm_) &
               & .and. ( TreeProcs(2)%PartType(2) .eq. Z0_) ) then
             tag_Z = 1
!              print *, "setting tag_Z=1 (2-cut) 2"
          else
             tag_Z=0
          endif


!         set momentum vector for last particle
          if( TreeProcs(1)%PartType(TreeProcs(1)%NumPart).eq.Glu_ ) then
                TreeProcs(1)%Gluons(TreeProcs(1)%NumGlu(0))%Mom => q2neg(:)
          elseif( IsAQuark(TreeProcs(1)%PartType(TreeProcs(1)%NumPart)) ) then
                TreeProcs(1)%Quarks(TreeProcs(1)%NumQua)%Mom => q2neg(:)
          elseif( IsAScalar(TreeProcs(1)%PartType(TreeProcs(1)%NumPart)) ) then
                TreeProcs(1)%Scalars(TreeProcs(1)%NumSca)%Mom => q2neg(:)
          endif
          if( TreeProcs(2)%PartType(TreeProcs(2)%NumPart).eq.Glu_ ) then
                TreeProcs(2)%Gluons(TreeProcs(2)%NumGlu(0))%Mom => q1neg(:)
          elseif( IsAQuark(TreeProcs(2)%PartType(TreeProcs(2)%NumPart)) ) then
                TreeProcs(2)%Quarks(TreeProcs(2)%NumQua)%Mom => q1neg(:)
          elseif( IsAScalar(TreeProcs(2)%PartType(TreeProcs(2)%NumPart)) ) then
                TreeProcs(2)%Scalars(TreeProcs(2)%NumSca)%Mom => q1neg(:)
          endif

          if (Ds.eq.4) then

          call givepol(TreeProcs(1)%PartType(1),q1(1:5),4,4,tag_pol,Nj1,BPOL1,POL1)
          call givepol(TreeProcs(2)%PartType(1),q2(1:5),4,4,tag_pol,Nj2,BPOL2,POL2)

!         calculate product of amplitudes
!         first  organize 4 lists that will be arguments
!         for amplitude procedure
!         start with the first label on the cut list ll
!         Lab_ex and Lab_in

          res=dcmplx(0d0,0d0)



!DEC$ IF (_DebugUseMyAmps==0)
            i1=l2c(2)-l2c(1)
            lab1=Lab_in(l2c(1))
            lab2=Lab_in(l2c(2))
            call ampl(TreeProcs(1),4,4,Nj1,Nj2,POL1,q1(1:5),lab1,l2c(1),tag_f,1,i1,BPOL2,q2(1:5),lab2,mur1)
!-----------new loop
            ia =i1+1
            ib=Npoint
            lab1=Lab_in(l2c(2))
            lab2=Lab_in(l2c(1))
            call ampl(TreeProcs(2),4,4,Nj2,Nj1,POL2,q2(1:5),lab1,l2c(2),tag_f,ia,ib,BPOL1,q1(1:5),lab2,mur2)
!DEC$ ELSE
            call new_ampl(4,4,Nj1,Nj2,POL1,BPOL2,tag_f,tag_Z,TreeProcs(1),mur1)
            call new_ampl(4,4,Nj2,Nj1,POL2,BPOL1,tag_f,tag_Z,TreeProcs(2),mur2)
!DEC$ ENDIF

            do j1=1,Nj1
            do j2=1,Nj2
               res=res+ mur1(j1,j2)*mur2(j2,j1)
            enddo
            enddo

!--------- to multiply by proper power of I:
           do j=1,2
              if (Lab_in(l2c(j)).eq.'top' .or. Lab_in(l2c(j)).eq.'bot' .or. Lab_in(l2c(j)).eq.'sto' .or. Lab_in(l2c(j)).eq.'sbo' .or. Lab_in(l2c(j)).eq.'str' .or. Lab_in(l2c(j)).eq.'chm') then
                res = dcmplx(0d0,1d0)*res
              else
                res = dcmplx(0d0,-1d0)*res
              endif
           enddo

!cccccccccccccccccccc endif for Ds = 4
         endif





          if (Ds.eq.5) then

          call givepol(TreeProcs(1)%PartType(1),q1(1:5),6,8,tag_pol,Nj1,BPOL1,POL1)
          call givepol(TreeProcs(2)%PartType(1),q2(1:5),6,8,tag_pol,Nj2,BPOL2,POL2)

!         calculate product of amplitudes
!         first  organize 4 lists that will be arguments
!         for amplitude procedure
!         start with the first label on the cut list ll
!         Lab_ex and Lab_in

          res=dcmplx(0d0,0d0)



!DEC$ IF (_DebugUseMyAmps==0)
            i1=l2c(2)-l2c(1)
            lab1=Lab_in(l2c(1))
            lab2=Lab_in(l2c(2))
            call ampl(TreeProcs(1),6,8,Nj1,Nj2,POL1,q1(1:5),lab1,l2c(1),tag_f,1,i1,BPOL2,q2(1:5),lab2,mur1)
!------------new loop
            ia =i1+1
            ib = Npoint
            lab1=Lab_in(l2c(2))
            lab2=Lab_in(l2c(1))
            call ampl(TreeProcs(2),6,8,Nj2,Nj1,POL2,q2(1:5),lab1,l2c(2),tag_f,ia,ib,BPOL1,q1(1:5),lab2,mur2)
!DEC$ ELSE
            call new_ampl(6,8,Nj1,Nj2,POL1,BPOL2,tag_f,tag_Z,TreeProcs(1),mur1)
            call new_ampl(6,8,Nj2,Nj1,POL2,BPOL1,tag_f,tag_Z,TreeProcs(2),mur2)
!DEC$ ENDIF

            do j1=1,Nj1
            do j2=1,Nj2
               res=res+ mur1(j1,j2)*mur2(j2,j1)
            enddo
            enddo


!           to multiply by proper power of I:
           do j=1,2
              if (Lab_in(l2c(j)).eq.'top' .or. Lab_in(l2c(j)).eq.'bot' .or. Lab_in(l2c(j)).eq.'sto' .or. Lab_in(l2c(j)).eq.'sbo' .or. Lab_in(l2c(j)).eq.'str' .or. Lab_in(l2c(j)).eq.'chm') then
                res = dcmplx(0d0,1d0)*res
              else
                 res = dcmplx(0d0,-1d0)*res
              endif
           enddo

!----------- intermediate 6-dim result
            res6 = res


!----------8-dim calculation

          call givepol(TreeProcs(1)%PartType(1),q1(1:5),8,16,tag_pol,Nj1,BPOL1,POL1)
          call givepol(TreeProcs(2)%PartType(1),q2(1:5),8,16,tag_pol,Nj2,BPOL2,POL2)

!         calculate product of amplitudes
!         first  organize 4 lists that will be arguments
!         for amplitude procedure
!         start with the first label on the cut list ll
!         Lab_ex and Lab_in

          res=dcmplx(0d0,0d0)


!DEC$ IF (_DebugUseMyAmps==0)
            lab1=Lab_in(l2c(1))
            lab2=Lab_in(l2c(2))
            call ampl(TreeProcs(1),8,16,Nj1,Nj2,POL1,q1(1:5),lab1,l2c(1),tag_f,1,i1,BPOL2,q2(1:5),lab2,mur1)
!------------new loop
            ia =i1+1
            ib = Npoint
            lab1=Lab_in(l2c(2))
            lab2=Lab_in(l2c(1))
            call ampl(TreeProcs(2),8,16,Nj2,Nj1,POL2,q2(1:5),lab1,l2c(2),tag_f,ia,ib,BPOL1,q1(1:5),lab2,mur2)
!DEC$ ELSE
            call new_ampl(8,16,Nj1,Nj2,POL1,BPOL2,tag_f,tag_Z,TreeProcs(1),mur1)
            call new_ampl(8,16,Nj2,Nj1,POL2,BPOL1,tag_f,tag_Z,TreeProcs(2),mur2)
!DEC$ ENDIF

            do j1=1,Nj1
            do j2=1,Nj2
               res=res+ mur1(j1,j2)*mur2(j2,j1)
            enddo
            enddo


!---------- multiply by proper power of I:
           do j=1,2
              if (Lab_in(l2c(j)).eq.'top' .or. Lab_in(l2c(j)).eq.'bot' .or. Lab_in(l2c(j)).eq.'sto' .or. Lab_in(l2c(j)).eq.'sbo' .or. Lab_in(l2c(j)).eq.'str' .or. Lab_in(l2c(j)).eq.'chm') then
                 res = dcmplx(0d0,1d0)*res
              else
                 res = dcmplx(0d0,-1d0)*res
              endif
           enddo

           res8 = res

!           final result
            if( tag_f.eq.0 ) then
               res = res6 -res8
            elseif (tag_f.eq.2) then
               res = res6 - res8/4d0
            elseif (tag_f.eq.3) then
               res = 2d0*res6 -res8
            else
               res = 2d0*res6 - res8
            endif


!cccccccccccccccccccc endif for Ds = 5
         endif

 

 
       return
       end subroutine





       subroutine resid1_new(lv,l1c,TreeProcs,res)
       use ModAmplitudes
       use ModProcess
       use ModMisc
       use ModParameters
       implicit none
       include 'misc/global_import'
       integer Nj1
       integer Ds,i,tag_pol,pos,pos1,pos2,pos3,pos4
       integer j,j1,lpos(3),lpos1(2),lpos3(4),tag_f,tag_Z
       integer ib,n14,n13,n12,n15
       integer l1c(1),im
       integer lmatch14(50),lmatch13(50)
       integer lmatch12(50),lmatch15(50)
       character lab1*3, lab2*3
       double complex lv(5),lvt(5)
       double complex v45(4),vprop(5),propX,r1,r2,krefa(4)
       double complex e(1),d(5),c(10),b(10)
       double complex prop1,prop2,prop3,prop4
       double complex vpol1(16),vpol2(16)
       double complex r3,r4
       double complex vne(5),trikoeff
       double complex v2(4),v3(4),v4(4)
       complex(8) BPOL1(8,16),POL1(8,16)
       complex(8) q1(1:8),mur1(10,10)
       complex(8),target :: q1neg(1:8)
       double complex res,res0,re
       double complex bkoeff
       type(TreeProcess) :: TreeProcs(:)




          if( Lab_in(l1c(1)).eq.'glu' ) then
            tag_f = 0
          elseif( Lab_in(l1c(1)).eq.'top' .or. Lab_in(l1c(1)).eq.'sto' .or. Lab_in(l1c(1)).eq.'str') then
            tag_f = 1
          elseif( Lab_in(l1c(1)).eq.'bot' .or. Lab_in(l1c(1)).eq.'chm' ) then
            tag_f = 2
          elseif( Lab_in(l1c(1)).eq.'sbo'  ) then
            tag_f = 3
          else
            tag_f = 99
          endif

          if ( (tag_f .eq. 1 .or. tag_f .eq. 2) .and.  ( any(Lab_ex(1:NumExtParticles) .eq. 'zee')) ) then
             tag_Z=1
!             print *, "setting tag_Z=1 (1-cut)"
          else
             tag_Z=0
          endif


          if( lv(5).eq.0d0 ) then
               Ds = 4
          else
               Ds = 5
          endif

          do i=1,4
          q1(i)=lv(i)
          enddo
          if (Ds.eq.4) then
          q1(5)=dcmplx(0d0,0d0)
          else
          q1(5)=lv(5)
          endif
          q1(6:8)=(0d0,0d0)
          q1neg(1:8)=-q1(1:8)

          tag_pol = 1


!         set momentum vector for last particle
          if( TreeProcs(1)%PartType(TreeProcs(1)%NumPart).eq.Glu_ ) then
                TreeProcs(1)%Gluons(TreeProcs(1)%NumGlu(0))%Mom => q1neg(:)
          elseif( IsAQuark(TreeProcs(1)%PartType(TreeProcs(1)%NumPart)) ) then
                TreeProcs(1)%Quarks(TreeProcs(1)%NumQua)%Mom => q1neg(:)
          elseif( IsAScalar(TreeProcs(1)%PartType(TreeProcs(1)%NumPart)) ) then
                TreeProcs(1)%Scalars(TreeProcs(1)%NumSca)%Mom => q1neg(:)
          endif

          if (Ds.eq.4) then
          call givepol(TreeProcs(1)%PartType(1),q1(1:5),4,4,tag_pol,Nj1,BPOL1,POL1)


          res=dcmplx(0d0,0d0)
!DEC$ IF (_DebugUseMyAmps==0)
            ib = Npoint
            lab1=Lab_in(l1c(1))
            lab2=Lab_in(l1c(1))
            call ampl(TreeProcs(1),4,4,Nj1,Nj1,POL1,q1(1:5),lab1,l1c(1),tag_f,1,ib,BPOL1,q1(1:5),lab2,mur1)
!DEC$ ELSE
            call new_ampl(4,4,Nj1,Nj1,POL1,BPOL1,tag_f,tag_Z,TreeProcs(1),mur1)
!DEC$ ENDIF


            do j1=1,Nj1
               res = res + mur1(j1,j1)
            enddo


              if (Lab_in(l1c(1)).eq.'top' .or. Lab_in(l1c(1)).eq.'bot' .or. Lab_in(l1c(1)).eq.'sto' .or. Lab_in(l1c(1)).eq.'sbo' .or. Lab_in(l1c(1)).eq.'str' .or. Lab_in(l1c(1)).eq.'chm') then
                res = dcmplx(0d0,1d0)*res
              else
                res = dcmplx(0d0,-1d0)*res
              endif

!ccccccccccccccccccccccccc Ds=4 end
          endif


 



 
          return
          end subroutine





       subroutine resid4Subtr(lv,ThePrimAmps,iPrimAmp,CutNum,res)
       use ModProcess
       use ModMisc
       implicit none
       double complex lv(5),res
       type(PrimitiveAmplitude),target :: ThePrimAmps(:)
       type(PrimitiveAmplitude),pointer :: ThePrimAmp
       integer i,CutNum,iPrimAmp,PrimSelect,k
       double complex lvt(5),propX
       double complex vprop(5)
       integer MatchCutNum,FirstCutMatch
       integer :: NumMatch,MatchHiCuts(1:10),FirstHiProp(1:10),MissHiProp(1:10,1:1)


!        subtracting the 5-cut

       ThePrimAmp => ThePrimAmps(iPrimAmp)

       do k=0,ThePrimAmp%NumSisters
         if(k.eq.0) then
             PrimSelect = iPrimAmp
         else
             PrimSelect = ThePrimAmp%Sisters(k)
         endif
         do i = 1,ThePrimAmp%UCuts(4)%Match(CutNum)%Subt(5)%NumMatch(k)
            MatchCutNum = ThePrimAmp%UCuts(4)%Match(CutNum)%Subt(5)%MatchHiCuts(k,i)
            FirstCutMatch =   ThePrimAmp%UCuts(4)%Match(CutNum)%Subt(5)%FirstHiProp(k,i)
            MissHiProp(i,1) = ThePrimAmp%UCuts(4)%Match(CutNum)%Subt(5)%MissHiProp(k,i,1)


            vprop(1:4) = lv(1:4) + dcmplx( -ThePrimAmps(PrimSelect)%IntPart( ThePrimAmp%UCuts(5)%CutProp(MatchCutNum,FirstCutMatch) )%Mom(1:4) &
                       + ThePrimAmps(PrimSelect)%IntPart(MissHiProp(i,1))%Mom(1:4) )
            vprop(5:5) = lv(5:5)
            call sc(5,vprop,vprop,propX)
            propX = propX-dcmplx(ThePrimAmps(PrimSelect)%IntPart( MissHiProp(i,1) )%Mass2,0d0)
!            res=res - ThePrimAmp%UCuts(5)%Coeff(MatchCutNum,0)*lv(5)**2/propX
            res=res - ThePrimAmps(PrimSelect)%UCuts(5)%Coeff(MatchCutNum,0)*lv(5)**2/propX
         enddo
      enddo

       return
       end subroutine




       subroutine resid3Subtr(lv,ThePrimAmps,iPrimAmp,CutNum,res)
       use ModProcess
       use ModMisc
       implicit none
       double complex v45(4),vprop(5),propX,r1,r2,lv(5),res,lvt(5)
       integer CutNum,i,iPrimAmp,k,PrimSelect
       double complex prop1,prop2
       double complex vne(5)
       type(PrimitiveAmplitude),target :: ThePrimAmps(:)
       type(PrimitiveAmplitude),pointer :: ThePrimAmp
       integer MatchCutNum,FirstCutMatch
       integer :: NumMatch,MatchHiCuts(1:10),FirstHiProp(1:10),MissHiProp(1:10,1:2)


       ThePrimAmp => ThePrimAmps(iPrimAmp)

       do k=0,ThePrimAmp%NumSisters

         if(k.eq.0) then
             PrimSelect = iPrimAmp
         else
             PrimSelect = ThePrimAmp%Sisters(k)
         endif


!        subtracting the 5-cut
         do i = 1,ThePrimAmp%UCuts(3)%Match(CutNum)%Subt(5)%NumMatch(k)
            MatchCutNum = ThePrimAmp%UCuts(3)%Match(CutNum)%Subt(5)%MatchHiCuts(k,i)
            FirstCutMatch =   ThePrimAmp%UCuts(3)%Match(CutNum)%Subt(5)%FirstHiProp(k,i)
            MissHiProp(i,1) = ThePrimAmp%UCuts(3)%Match(CutNum)%Subt(5)%MissHiProp(k,i,1)
            MissHiProp(i,2) = ThePrimAmp%UCuts(3)%Match(CutNum)%Subt(5)%MissHiProp(k,i,2)

            vprop(1:4) = lv(1:4) + dcmplx( -ThePrimAmps(PrimSelect)%IntPart(ThePrimAmp%UCuts(5)%CutProp(MatchCutNum,FirstCutMatch))%Mom(1:4) &
                       + ThePrimAmps(PrimSelect)%IntPart(MissHiProp(i,1))%Mom(1:4) )
            vprop(5:5) = lv(5:5)
            call sc(5,vprop,vprop,propX)
            prop1 = propX-dcmplx(ThePrimAmps(PrimSelect)%IntPart( MissHiProp(i,1) )%Mass2,0d0)

            vprop(1:4) = lv(1:4) + dcmplx( -ThePrimAmps(PrimSelect)%IntPart(ThePrimAmp%UCuts(5)%CutProp(MatchCutNum,FirstCutMatch))%Mom(1:4) &
                       + ThePrimAmps(PrimSelect)%IntPart(MissHiProp(i,2))%Mom(1:4) )
            vprop(5:5) = lv(5:5)
            call sc(5,vprop,vprop,propX)
            prop2 = propX-dcmplx(ThePrimAmps(PrimSelect)%IntPart( MissHiProp(i,2) )%Mass2,0d0)

!            res=res - ThePrimAmp%UCuts(5)%Coeff(MatchCutNum,0)*lv(5)**2/prop1/prop2
            res=res - ThePrimAmps(PrimSelect)%UCuts(5)%Coeff(MatchCutNum,0)*lv(5)**2/prop1/prop2
         enddo




!        subtracting the 4-cut
         do i = 1,ThePrimAmp%UCuts(3)%Match(CutNum)%Subt(4)%NumMatch(k)
            MatchCutNum = ThePrimAmp%UCuts(3)%Match(CutNum)%Subt(4)%MatchHiCuts(k,i)
            FirstCutMatch =   ThePrimAmp%UCuts(3)%Match(CutNum)%Subt(4)%FirstHiProp(k,i)
            MissHiProp(i,1) = ThePrimAmp%UCuts(3)%Match(CutNum)%Subt(4)%MissHiProp(k,i,1)

            vprop(1:4) = lv(1:4) + dcmplx( -ThePrimAmps(PrimSelect)%IntPart(ThePrimAmp%UCuts(4)%CutProp(MatchCutNum,FirstCutMatch))%Mom(1:4) &
                       + ThePrimAmps(PrimSelect)%IntPart(MissHiProp(i,1))%Mom(1:4) )
            vprop(5:5) = lv(5:5)
            call sc(5,vprop,vprop,propX)
            propX = propX-dcmplx(ThePrimAmps(PrimSelect)%IntPart( MissHiProp(i,1) )%Mass2,0d0)

            lvt(1:5) = lv(1:5)
            if( FirstCutMatch.gt.1 ) lvt(1:4) = lvt(1:4) - ThePrimAmps(PrimSelect)%UCuts(4)%KMom(MatchCutNum,FirstCutMatch-1,1:4)
            v45(1:4) = ThePrimAmps(PrimSelect)%UCuts(4)%NMom(MatchCutNum,1,1:4)
            call sc(4,v45,lvt,r1)
            r2 = -(0d0,1d0)*lvt(5) != (0,0,0,0,i).dot.lvt(1:5)

!            res=res-(ThePrimAmp%UCuts(4)%Coeff(MatchCutNum,0)     &
!                    +ThePrimAmp%UCuts(4)%Coeff(MatchCutNum,1)*r1   &
!                    +(r2**2)*(ThePrimAmp%UCuts(4)%Coeff(MatchCutNum,2)    &
!                    +ThePrimAmp%UCuts(4)%Coeff(MatchCutNum,3)*r1    &
!                    +ThePrimAmp%UCuts(4)%Coeff(MatchCutNum,4)*(r2**2))    &
!                    )/propX

            res=res-(ThePrimAmps(PrimSelect)%UCuts(4)%Coeff(MatchCutNum,0)     &
                    +ThePrimAmps(PrimSelect)%UCuts(4)%Coeff(MatchCutNum,1)*r1   &
                    +(r2**2)*(ThePrimAmps(PrimSelect)%UCuts(4)%Coeff(MatchCutNum,2)    &
                    +ThePrimAmps(PrimSelect)%UCuts(4)%Coeff(MatchCutNum,3)*r1    &
                    +ThePrimAmps(PrimSelect)%UCuts(4)%Coeff(MatchCutNum,4)*(r2**2))    &
                    )/propX
       enddo


      enddo

       return
       end subroutine






       subroutine resid2Subtr(lv,ThePrimAmps,iPrimAmp,CutNum,res)
       use ModProcess
       use ModMisc
       implicit none
       integer CutNum,iPrimAmp
       double complex lv(1:5),res
       double complex v45(4),vprop(5),propX,r1,r2,r22,lvt(5)
       type(PrimitiveAmplitude),target :: ThePrimAmps(:)
       type(PrimitiveAmplitude),pointer :: ThePrimAmp
       double complex prop1,prop2,prop3
       double complex r3,r4,re
       double complex vne(5)
       double complex v3(4),v4(4),r32,r42
       integer i,MatchCutNum,FirstCutMatch,k,PrimSelect
       integer :: NumMatch,MatchHiCuts(1:10),FirstHiProp(1:10),MissHiProp(1:10,1:3)


       ThePrimAmp => ThePrimAmps(iPrimAmp)

       do k=0,ThePrimAmp%NumSisters

         if(k.eq.0) then
             PrimSelect = iPrimAmp
         else
             PrimSelect = ThePrimAmp%Sisters(k)
         endif



!        subtracting the 5-cut
         do i = 1,ThePrimAmp%UCuts(2)%Match(CutNum)%Subt(5)%NumMatch(k)
            MatchCutNum=ThePrimAmp%UCuts(2)%Match(CutNum)%Subt(5)%MatchHiCuts(k,i)
            FirstCutMatch=ThePrimAmp%UCuts(2)%Match(CutNum)%Subt(5)%FirstHiProp(k,i)
            MissHiProp(i,1) = ThePrimAmp%UCuts(2)%Match(CutNum)%Subt(5)%MissHiProp(k,i,1)
            MissHiProp(i,2) = ThePrimAmp%UCuts(2)%Match(CutNum)%Subt(5)%MissHiProp(k,i,2)
            MissHiProp(i,3) = ThePrimAmp%UCuts(2)%Match(CutNum)%Subt(5)%MissHiProp(k,i,3)

            vprop(1:4) = lv(1:4) + dcmplx( -ThePrimAmps(PrimSelect)%IntPart(ThePrimAmp%UCuts(5)%CutProp(MatchCutNum,FirstCutMatch))%Mom(1:4) &
                       + ThePrimAmps(PrimSelect)%IntPart(MissHiProp(i,1))%Mom(1:4) )
            vprop(5:5) = lv(5:5)
            call sc(5,vprop,vprop,propX)
            prop1 = propX-dcmplx(ThePrimAmps(PrimSelect)%IntPart( MissHiProp(i,1) )%Mass2,0d0)

            vprop(1:4) = lv(1:4) + dcmplx( -ThePrimAmps(PrimSelect)%IntPart(ThePrimAmp%UCuts(5)%CutProp(MatchCutNum,FirstCutMatch))%Mom(1:4) &
                       + ThePrimAmps(PrimSelect)%IntPart(MissHiProp(i,2))%Mom(1:4) )
            vprop(5:5) = lv(5:5)
            call sc(5,vprop,vprop,propX)
            prop2 = propX-dcmplx(ThePrimAmps(PrimSelect)%IntPart( MissHiProp(i,2) )%Mass2,0d0)

            vprop(1:4) = lv(1:4) + dcmplx( -ThePrimAmps(PrimSelect)%IntPart(ThePrimAmp%UCuts(5)%CutProp(MatchCutNum,FirstCutMatch))%Mom(1:4) &
                       + ThePrimAmps(PrimSelect)%IntPart(MissHiProp(i,3))%Mom(1:4) )
            vprop(5:5) = lv(5:5)
            call sc(5,vprop,vprop,propX)
            prop3 = propX-dcmplx(ThePrimAmps(PrimSelect)%IntPart( MissHiProp(i,3) )%Mass2,0d0)

!            res=res - ThePrimAmp%UCuts(5)%Coeff(MatchCutNum,0)*lv(5)**2/prop1/prop2/prop3
            res=res - ThePrimAmps(PrimSelect)%UCuts(5)%Coeff(MatchCutNum,0)*lv(5)**2/prop1/prop2/prop3
         enddo





!        subtracting the 4-cut
         do i = 1,ThePrimAmp%UCuts(2)%Match(CutNum)%Subt(4)%NumMatch(k)
            MatchCutNum = ThePrimAmp%UCuts(2)%Match(CutNum)%Subt(4)%MatchHiCuts(k,i)
            FirstCutMatch =   ThePrimAmp%UCuts(2)%Match(CutNum)%Subt(4)%FirstHiProp(k,i)
            MissHiProp(i,1) = ThePrimAmp%UCuts(2)%Match(CutNum)%Subt(4)%MissHiProp(k,i,1)
            MissHiProp(i,2) = ThePrimAmp%UCuts(2)%Match(CutNum)%Subt(4)%MissHiProp(k,i,2)

            vprop(1:4) = lv(1:4) + dcmplx( -ThePrimAmps(PrimSelect)%IntPart(ThePrimAmp%UCuts(4)%CutProp(MatchCutNum,FirstCutMatch))%Mom(1:4) &
                       + ThePrimAmps(PrimSelect)%IntPart(MissHiProp(i,1))%Mom(1:4) )
            vprop(5:5) = lv(5:5)
            call sc(5,vprop,vprop,propX)
            prop1 = propX-dcmplx(ThePrimAmps(PrimSelect)%IntPart( MissHiProp(i,1) )%Mass2,0d0)

            vprop(1:4) = lv(1:4) + dcmplx( -ThePrimAmps(PrimSelect)%IntPart(ThePrimAmp%UCuts(4)%CutProp(MatchCutNum,FirstCutMatch))%Mom(1:4) &
                       + ThePrimAmps(PrimSelect)%IntPart(MissHiProp(i,2))%Mom(1:4) )
            vprop(5:5) = lv(5:5)
            call sc(5,vprop,vprop,propX)
            prop2 = propX-dcmplx(ThePrimAmps(PrimSelect)%IntPart( MissHiProp(i,2) )%Mass2,0d0)

            lvt(1:5) = lv(1:5)
            if( FirstCutMatch.gt.1 ) lvt(1:4) = lvt(1:4) - ThePrimAmps(PrimSelect)%UCuts(4)%KMom(MatchCutNum,FirstCutMatch-1,1:4)
            v45(1:4) = ThePrimAmps(PrimSelect)%UCuts(4)%NMom(MatchCutNum,1,1:4)
            call sc(4,v45,lvt,r1)
            r2 = -(0d0,1d0)*lvt(5) != (0,0,0,0,i).dot.lvt(1:5)

!            res=res-(ThePrimAmp%UCuts(4)%Coeff(MatchCutNum,0)     &
!                    +ThePrimAmp%UCuts(4)%Coeff(MatchCutNum,1)*r1   &
!                    +(r2**2)*(ThePrimAmp%UCuts(4)%Coeff(MatchCutNum,2)    &
!                    +ThePrimAmp%UCuts(4)%Coeff(MatchCutNum,3)*r1    &
!                    +ThePrimAmp%UCuts(4)%Coeff(MatchCutNum,4)*(r2**2))    &
!                    )/prop1/prop2

            res=res-(ThePrimAmps(PrimSelect)%UCuts(4)%Coeff(MatchCutNum,0)     &
                    +ThePrimAmps(PrimSelect)%UCuts(4)%Coeff(MatchCutNum,1)*r1   &
                    +(r2**2)*(ThePrimAmps(PrimSelect)%UCuts(4)%Coeff(MatchCutNum,2)    &
                    +ThePrimAmps(PrimSelect)%UCuts(4)%Coeff(MatchCutNum,3)*r1    &
                    +ThePrimAmps(PrimSelect)%UCuts(4)%Coeff(MatchCutNum,4)*(r2**2))    &
                    )/prop1/prop2
       enddo






!      subtracting the 3-cut
       do i = 1,ThePrimAmp%UCuts(2)%Match(CutNum)%Subt(3)%NumMatch(k)
            MatchCutNum = ThePrimAmp%UCuts(2)%Match(CutNum)%Subt(3)%MatchHiCuts(k,i)
            FirstCutMatch =   ThePrimAmp%UCuts(2)%Match(CutNum)%Subt(3)%FirstHiProp(k,i)
            MissHiProp(i,1) = ThePrimAmp%UCuts(2)%Match(CutNum)%Subt(3)%MissHiProp(k,i,1)

            vprop(1:4) = lv(1:4) + dcmplx( -ThePrimAmps(PrimSelect)%IntPart(ThePrimAmp%UCuts(3)%CutProp(MatchCutNum,FirstCutMatch))%Mom(1:4) &
                       + ThePrimAmps(PrimSelect)%IntPart(MissHiProp(i,1))%Mom(1:4) )
            vprop(5:5) = lv(5:5)
            call sc(5,vprop,vprop,propX)
            propX = propX-dcmplx(ThePrimAmps(PrimSelect)%IntPart( MissHiProp(i,1) )%Mass2,0d0)

            lvt(1:5) = lv(1:5)
            if( FirstCutMatch.gt.1 ) lvt(1:4) = lvt(1:4) - ThePrimAmps(PrimSelect)%UCuts(3)%KMom(MatchCutNum,FirstCutMatch-1,1:4)
            v3(1:4) = ThePrimAmps(PrimSelect)%UCuts(3)%NMom(MatchCutNum,1,1:4)
            v4(1:4) = ThePrimAmps(PrimSelect)%UCuts(3)%NMom(MatchCutNum,2,1:4)
            call sc(4,v3,lvt,r3)
            call sc(4,v4,lvt,r4)
            re = -(0d0,1d0)*lvt(5) != (0,0,0,0,i).dot.lvt(1:5)

!            res=res - ( ThePrimAmp%UCuts(3)%Coeff(MatchCutNum,0)    &
!                      + ThePrimAmp%UCuts(3)%Coeff(MatchCutNum,1)*r3    &
!                      + ThePrimAmp%UCuts(3)%Coeff(MatchCutNum,2)*r4    &
!                      + ThePrimAmp%UCuts(3)%Coeff(MatchCutNum,3)*r3*r4    &
!                      + ThePrimAmp%UCuts(3)%Coeff(MatchCutNum,4)*(r3**2-r4**2)    &
!                      + ThePrimAmp%UCuts(3)%Coeff(MatchCutNum,5)*r3**2*r4    &
!                      + ThePrimAmp%UCuts(3)%Coeff(MatchCutNum,6)*r3*r4**2    &
!                      + re**2*( ThePrimAmp%UCuts(3)%Coeff(MatchCutNum,7)+ThePrimAmp%UCuts(3)%Coeff(MatchCutNum,8)*r3    &
!                              + ThePrimAmp%UCuts(3)%Coeff(MatchCutNum,9)*r4)     &
!                      )/propX

            res=res - ( ThePrimAmps(PrimSelect)%UCuts(3)%Coeff(MatchCutNum,0)    &
                      + ThePrimAmps(PrimSelect)%UCuts(3)%Coeff(MatchCutNum,1)*r3    &
                      + ThePrimAmps(PrimSelect)%UCuts(3)%Coeff(MatchCutNum,2)*r4    &
                      + ThePrimAmps(PrimSelect)%UCuts(3)%Coeff(MatchCutNum,3)*r3*r4    &
                      + ThePrimAmps(PrimSelect)%UCuts(3)%Coeff(MatchCutNum,4)*(r3**2-r4**2)    &
                      + ThePrimAmps(PrimSelect)%UCuts(3)%Coeff(MatchCutNum,5)*r3**2*r4    &
                      + ThePrimAmps(PrimSelect)%UCuts(3)%Coeff(MatchCutNum,6)*r3*r4**2    &
                      + re**2*( ThePrimAmps(PrimSelect)%UCuts(3)%Coeff(MatchCutNum,7)+ThePrimAmps(PrimSelect)%UCuts(3)%Coeff(MatchCutNum,8)*r3    &
                              + ThePrimAmps(PrimSelect)%UCuts(3)%Coeff(MatchCutNum,9)*r4)     &
                      )/propX
       enddo


       enddo

       return
       end subroutine







       subroutine resid1Subtr(lv,ThePrimAmps,iPrimAmp,CutNum,res)
       use ModProcess
       use ModMisc
       implicit none
       include 'misc/global_import'!   this is only needed for tagcut(:)
       integer CutNum,l1c(1:1),iPrimAmp
       double complex lv(1:5),res,lvt(5),re
       double complex v45(4),vprop(5),propX,r1,r2
       double complex prop1,prop2,prop3,prop4
       double complex r3,r4,vne(5)
       double complex v2(4),v3(4),v4(4)
       type(PrimitiveAmplitude),target :: ThePrimAmps(:)
       type(PrimitiveAmplitude),pointer :: ThePrimAmp
       integer MatchCutNum,FirstCutMatch,k,PrimSelect
       integer :: i,NumMatch,MatchHiCuts(1:10),FirstHiProp(1:10),MissHiProp(1:10,1:4)


       ThePrimAmp => ThePrimAmps(iPrimAmp)


       do k=0,ThePrimAmp%NumSisters

         if(k.eq.0) then
             PrimSelect = iPrimAmp
         else
             PrimSelect = ThePrimAmp%Sisters(k)
         endif


!        subtracting the 5-cut
         do i = 1,ThePrimAmp%UCuts(1)%Match(CutNum)%Subt(5)%NumMatch(k)
            MatchCutNum=ThePrimAmp%UCuts(1)%Match(CutNum)%Subt(5)%MatchHiCuts(k,i)
            FirstCutMatch=ThePrimAmp%UCuts(1)%Match(CutNum)%Subt(5)%FirstHiProp(k,i)
            MissHiProp(i,1) = ThePrimAmp%UCuts(1)%Match(CutNum)%Subt(5)%MissHiProp(k,i,1)
            MissHiProp(i,2) = ThePrimAmp%UCuts(1)%Match(CutNum)%Subt(5)%MissHiProp(k,i,2)
            MissHiProp(i,3) = ThePrimAmp%UCuts(1)%Match(CutNum)%Subt(5)%MissHiProp(k,i,3)
            MissHiProp(i,4) = ThePrimAmp%UCuts(1)%Match(CutNum)%Subt(5)%MissHiProp(k,i,4)

            vprop(1:4) = lv(1:4) + dcmplx( -ThePrimAmps(PrimSelect)%IntPart(ThePrimAmp%UCuts(5)%CutProp(MatchCutNum,FirstCutMatch))%Mom(1:4) &
                       + ThePrimAmps(PrimSelect)%IntPart(MissHiProp(i,1))%Mom(1:4) )
            vprop(5:5) = lv(5:5)
            call sc(5,vprop,vprop,propX)
            prop1 = propX-dcmplx(ThePrimAmps(PrimSelect)%IntPart( MissHiProp(i,1) )%Mass2,0d0)

            vprop(1:4) = lv(1:4) + dcmplx( -ThePrimAmps(PrimSelect)%IntPart(ThePrimAmp%UCuts(5)%CutProp(MatchCutNum,FirstCutMatch))%Mom(1:4) &
                       + ThePrimAmps(PrimSelect)%IntPart(MissHiProp(i,2))%Mom(1:4) )
            vprop(5:5) = lv(5:5)
            call sc(5,vprop,vprop,propX)
            prop2 = propX-dcmplx(ThePrimAmps(PrimSelect)%IntPart( MissHiProp(i,2) )%Mass2,0d0)

            vprop(1:4) = lv(1:4) + dcmplx( -ThePrimAmps(PrimSelect)%IntPart(ThePrimAmp%UCuts(5)%CutProp(MatchCutNum,FirstCutMatch))%Mom(1:4) &
                       + ThePrimAmps(PrimSelect)%IntPart(MissHiProp(i,3))%Mom(1:4) )
            vprop(5:5) = lv(5:5)
            call sc(5,vprop,vprop,propX)
            prop3 = propX-dcmplx(ThePrimAmps(PrimSelect)%IntPart( MissHiProp(i,3) )%Mass2,0d0)

            vprop(1:4) = lv(1:4) + dcmplx( -ThePrimAmps(PrimSelect)%IntPart(ThePrimAmp%UCuts(5)%CutProp(MatchCutNum,FirstCutMatch))%Mom(1:4) &
                       + ThePrimAmps(PrimSelect)%IntPart(MissHiProp(i,4))%Mom(1:4) )
            vprop(5:5) = lv(5:5)
            call sc(5,vprop,vprop,propX)
            prop4 = propX-dcmplx(ThePrimAmps(PrimSelect)%IntPart( MissHiProp(i,4) )%Mass2,0d0)

!            res=res - ThePrimAmp%UCuts(5)%Coeff(MatchCutNum,0)*lv(5)**2/prop1/prop2/prop3/prop4
            res=res - ThePrimAmps(PrimSelect)%UCuts(5)%Coeff(MatchCutNum,0)*lv(5)**2/prop1/prop2/prop3/prop4
         enddo




!        subtracting the 4-cut
         do i = 1,ThePrimAmp%UCuts(1)%Match(CutNum)%Subt(4)%NumMatch(k)
            MatchCutNum = ThePrimAmp%UCuts(1)%Match(CutNum)%Subt(4)%MatchHiCuts(k,i)
            FirstCutMatch =   ThePrimAmp%UCuts(1)%Match(CutNum)%Subt(4)%FirstHiProp(k,i)
            MissHiProp(i,1) = ThePrimAmp%UCuts(1)%Match(CutNum)%Subt(4)%MissHiProp(k,i,1)
            MissHiProp(i,2) = ThePrimAmp%UCuts(1)%Match(CutNum)%Subt(4)%MissHiProp(k,i,2)
            MissHiProp(i,3) = ThePrimAmp%UCuts(1)%Match(CutNum)%Subt(4)%MissHiProp(k,i,3)

            vprop(1:4) = lv(1:4) + dcmplx( -ThePrimAmps(PrimSelect)%IntPart(ThePrimAmp%UCuts(4)%CutProp(MatchCutNum,FirstCutMatch))%Mom(1:4) &
                       + ThePrimAmps(PrimSelect)%IntPart(MissHiProp(i,1))%Mom(1:4) )
            vprop(5:5) = lv(5:5)
            call sc(5,vprop,vprop,propX)
            prop1 = propX-dcmplx(ThePrimAmps(PrimSelect)%IntPart( MissHiProp(i,1) )%Mass2,0d0)

            vprop(1:4) = lv(1:4) + dcmplx( -ThePrimAmps(PrimSelect)%IntPart(ThePrimAmp%UCuts(4)%CutProp(MatchCutNum,FirstCutMatch))%Mom(1:4) &
                       + ThePrimAmps(PrimSelect)%IntPart(MissHiProp(i,2))%Mom(1:4) )
            vprop(5:5) = lv(5:5)
            call sc(5,vprop,vprop,propX)
            prop2 = propX-dcmplx(ThePrimAmps(PrimSelect)%IntPart( MissHiProp(i,2) )%Mass2,0d0)

            vprop(1:4) = lv(1:4) + dcmplx( -ThePrimAmps(PrimSelect)%IntPart(ThePrimAmp%UCuts(4)%CutProp(MatchCutNum,FirstCutMatch))%Mom(1:4) &
                       + ThePrimAmps(PrimSelect)%IntPart(MissHiProp(i,3))%Mom(1:4) )
            vprop(5:5) = lv(5:5)
            call sc(5,vprop,vprop,propX)
            prop3 = propX-dcmplx(ThePrimAmps(PrimSelect)%IntPart( MissHiProp(i,3) )%Mass2,0d0)

            lvt(1:5) = lv(1:5)
            if( FirstCutMatch.gt.1 ) lvt(1:4) = lvt(1:4) - ThePrimAmps(PrimSelect)%UCuts(4)%KMom(MatchCutNum,FirstCutMatch-1,1:4)
            v45(1:4) = ThePrimAmps(PrimSelect)%UCuts(4)%NMom(MatchCutNum,1,1:4)
            call sc(4,v45,lvt,r1)
            r2 = -(0d0,1d0)*lvt(5) != (0,0,0,0,i).dot.lvt(1:5)

!            res=res-(ThePrimAmp%UCuts(4)%Coeff(MatchCutNum,0)     &
!                    +ThePrimAmp%UCuts(4)%Coeff(MatchCutNum,1)*r1   &
!                    +(r2**2)*(ThePrimAmp%UCuts(4)%Coeff(MatchCutNum,2)    &
!                    +ThePrimAmp%UCuts(4)%Coeff(MatchCutNum,3)*r1    &
!                    +ThePrimAmp%UCuts(4)%Coeff(MatchCutNum,4)*(r2**2))    &
!                    )/prop1/prop2/prop3

            res=res-(ThePrimAmps(PrimSelect)%UCuts(4)%Coeff(MatchCutNum,0)     &
                    +ThePrimAmps(PrimSelect)%UCuts(4)%Coeff(MatchCutNum,1)*r1   &
                    +(r2**2)*(ThePrimAmps(PrimSelect)%UCuts(4)%Coeff(MatchCutNum,2)    &
                    +ThePrimAmps(PrimSelect)%UCuts(4)%Coeff(MatchCutNum,3)*r1    &
                    +ThePrimAmps(PrimSelect)%UCuts(4)%Coeff(MatchCutNum,4)*(r2**2))    &
                    )/prop1/prop2/prop3
         enddo




!        subtracting the 3-cut
         do i = 1,ThePrimAmp%UCuts(1)%Match(CutNum)%Subt(3)%NumMatch(k)
            MatchCutNum = ThePrimAmp%UCuts(1)%Match(CutNum)%Subt(3)%MatchHiCuts(k,i)
            FirstCutMatch =   ThePrimAmp%UCuts(1)%Match(CutNum)%Subt(3)%FirstHiProp(k,i)
            MissHiProp(i,1) = ThePrimAmp%UCuts(1)%Match(CutNum)%Subt(3)%MissHiProp(k,i,1)
            MissHiProp(i,2) = ThePrimAmp%UCuts(1)%Match(CutNum)%Subt(3)%MissHiProp(k,i,2)

            vprop(1:4) = lv(1:4) + dcmplx( -ThePrimAmps(PrimSelect)%IntPart(ThePrimAmp%UCuts(3)%CutProp(MatchCutNum,FirstCutMatch))%Mom(1:4) &
                       + ThePrimAmps(PrimSelect)%IntPart(MissHiProp(i,1))%Mom(1:4) )
            vprop(5:5) = lv(5:5)
            call sc(5,vprop,vprop,propX)
            prop1 = propX-dcmplx(ThePrimAmps(PrimSelect)%IntPart( MissHiProp(i,1) )%Mass2,0d0)

            vprop(1:4) = lv(1:4) + dcmplx( -ThePrimAmps(PrimSelect)%IntPart(ThePrimAmp%UCuts(3)%CutProp(MatchCutNum,FirstCutMatch))%Mom(1:4) &
                       + ThePrimAmps(PrimSelect)%IntPart(MissHiProp(i,2))%Mom(1:4) )
            vprop(5:5) = lv(5:5)
            call sc(5,vprop,vprop,propX)
            prop2 = propX-dcmplx(ThePrimAmps(PrimSelect)%IntPart( MissHiProp(i,2) )%Mass2,0d0)

            lvt(1:5) = lv(1:5)
            if( FirstCutMatch.gt.1 ) lvt(1:4) = lvt(1:4) - ThePrimAmps(PrimSelect)%UCuts(3)%KMom(MatchCutNum,FirstCutMatch-1,1:4)
            v3(1:4) = ThePrimAmps(PrimSelect)%UCuts(3)%NMom(MatchCutNum,1,1:4)
            v4(1:4) = ThePrimAmps(PrimSelect)%UCuts(3)%NMom(MatchCutNum,2,1:4)
            call sc(4,v3,lvt,r3)
            call sc(4,v4,lvt,r4)
            re = -(0d0,1d0)*lvt(5) != (0,0,0,0,i).dot.lvt(1:5)

!            res=res - ( ThePrimAmp%UCuts(3)%Coeff(MatchCutNum,0)    &
!                      + ThePrimAmp%UCuts(3)%Coeff(MatchCutNum,1)*r3    &
!                      + ThePrimAmp%UCuts(3)%Coeff(MatchCutNum,2)*r4    &
!                      + ThePrimAmp%UCuts(3)%Coeff(MatchCutNum,3)*r3*r4    &
!                      + ThePrimAmp%UCuts(3)%Coeff(MatchCutNum,4)*(r3**2-r4**2)    &
!                      + ThePrimAmp%UCuts(3)%Coeff(MatchCutNum,5)*r3**2*r4    &
!                      + ThePrimAmp%UCuts(3)%Coeff(MatchCutNum,6)*r3*r4**2    &
!                      + re**2*( ThePrimAmp%UCuts(3)%Coeff(MatchCutNum,7)+ThePrimAmp%UCuts(3)%Coeff(MatchCutNum,8)*r3    &
!                              + ThePrimAmp%UCuts(3)%Coeff(MatchCutNum,9)*r4)     &
!                      )/prop1/prop2

            res=res - ( ThePrimAmps(PrimSelect)%UCuts(3)%Coeff(MatchCutNum,0)    &
                      + ThePrimAmps(PrimSelect)%UCuts(3)%Coeff(MatchCutNum,1)*r3    &
                      + ThePrimAmps(PrimSelect)%UCuts(3)%Coeff(MatchCutNum,2)*r4    &
                      + ThePrimAmps(PrimSelect)%UCuts(3)%Coeff(MatchCutNum,3)*r3*r4    &
                      + ThePrimAmps(PrimSelect)%UCuts(3)%Coeff(MatchCutNum,4)*(r3**2-r4**2)    &
                      + ThePrimAmps(PrimSelect)%UCuts(3)%Coeff(MatchCutNum,5)*r3**2*r4    &
                      + ThePrimAmps(PrimSelect)%UCuts(3)%Coeff(MatchCutNum,6)*r3*r4**2    &
                      + re**2*( ThePrimAmps(PrimSelect)%UCuts(3)%Coeff(MatchCutNum,7)+ThePrimAmps(PrimSelect)%UCuts(3)%Coeff(MatchCutNum,8)*r3    &
                              + ThePrimAmps(PrimSelect)%UCuts(3)%Coeff(MatchCutNum,9)*r4)     &
                      )/prop1/prop2
        enddo




!        subtracting the 2-cut
         do i = 1,ThePrimAmp%UCuts(1)%Match(CutNum)%Subt(2)%NumMatch(k)
            MatchCutNum = ThePrimAmp%UCuts(1)%Match(CutNum)%Subt(2)%MatchHiCuts(k,i)
            FirstCutMatch =   ThePrimAmp%UCuts(1)%Match(CutNum)%Subt(2)%FirstHiProp(k,i)
            MissHiProp(i,1) = ThePrimAmp%UCuts(1)%Match(CutNum)%Subt(2)%MissHiProp(k,i,1)

            vprop(1:4) = lv(1:4) + dcmplx( -ThePrimAmps(PrimSelect)%IntPart(ThePrimAmp%UCuts(2)%CutProp(MatchCutNum,FirstCutMatch))%Mom(1:4) &
                       + ThePrimAmps(PrimSelect)%IntPart(MissHiProp(i,1))%Mom(1:4) )
            vprop(5:5) = lv(5:5)
            call sc(5,vprop,vprop,propX)
            prop1 = propX-dcmplx(ThePrimAmps(PrimSelect)%IntPart( MissHiProp(i,1) )%Mass2,0d0)

            lvt(1:5) = lv(1:5)
            if( FirstCutMatch.gt.1 ) lvt(1:4) = lvt(1:4) - ThePrimAmps(PrimSelect)%UCuts(2)%KMom(MatchCutNum,FirstCutMatch-1,1:4)
            v2(1:4) = ThePrimAmps(PrimSelect)%UCuts(2)%NMom(MatchCutNum,1,1:4)
            v3(1:4) = ThePrimAmps(PrimSelect)%UCuts(2)%NMom(MatchCutNum,2,1:4)
            v4(1:4) = ThePrimAmps(PrimSelect)%UCuts(2)%NMom(MatchCutNum,3,1:4)
            call sc(4,v2,lvt,r2)
            call sc(4,v3,lvt,r3)
            call sc(4,v4,lvt,r4)
            re = -(0d0,1d0)*lvt(5) != (0,0,0,0,i).dot.lvt(1:5)

!            if (tagdcut(MatchCutNum,1).eq.666) then
            if (ThePrimAmps(PrimSelect)%UCuts(2)%tagcuts(MatchCutNum).eq.666) then
!              res = res - ( ThePrimAmp%UCuts(2)%Coeff(MatchCutNum,0)     &
!                          + ThePrimAmp%UCuts(2)%Coeff(MatchCutNum,1)*r2     &
!                          + ThePrimAmp%UCuts(2)%Coeff(MatchCutNum,2)*r3     &
!                          + ThePrimAmp%UCuts(2)%Coeff(MatchCutNum,3)*r4     &
!                          + ThePrimAmp%UCuts(2)%Coeff(MatchCutNum,4)*(r2**2-r3**2)     &
!                          + ThePrimAmp%UCuts(2)%Coeff(MatchCutNum,5)*(r2**2+r3**2-2*r4**2)     &
!                          + ThePrimAmp%UCuts(2)%Coeff(MatchCutNum,6)*r2*r3     &
!                          + ThePrimAmp%UCuts(2)%Coeff(MatchCutNum,7)*r2*r4           &
!                          + ThePrimAmp%UCuts(2)%Coeff(MatchCutNum,8)*r3*r4     &
!                          + ThePrimAmp%UCuts(2)%Coeff(MatchCutNum,9)*re**2 )/prop1
!
              res = res - ( ThePrimAmps(PrimSelect)%UCuts(2)%Coeff(MatchCutNum,0)     &
                          + ThePrimAmps(PrimSelect)%UCuts(2)%Coeff(MatchCutNum,1)*r2     &
                          + ThePrimAmps(PrimSelect)%UCuts(2)%Coeff(MatchCutNum,2)*r3     &
                          + ThePrimAmps(PrimSelect)%UCuts(2)%Coeff(MatchCutNum,3)*r4     &
                          + ThePrimAmps(PrimSelect)%UCuts(2)%Coeff(MatchCutNum,4)*(r2**2-r3**2)     &
                          + ThePrimAmps(PrimSelect)%UCuts(2)%Coeff(MatchCutNum,5)*(r2**2+r3**2-2*r4**2)     &
                          + ThePrimAmps(PrimSelect)%UCuts(2)%Coeff(MatchCutNum,6)*r2*r3     &
                          + ThePrimAmps(PrimSelect)%UCuts(2)%Coeff(MatchCutNum,7)*r2*r4           &
                          + ThePrimAmps(PrimSelect)%UCuts(2)%Coeff(MatchCutNum,8)*r3*r4     &
                          + ThePrimAmps(PrimSelect)%UCuts(2)%Coeff(MatchCutNum,9)*re**2 )/prop1

!            elseif (tagdcut(MatchCutNum,1).eq.999) then
           elseif (ThePrimAmps(PrimSelect)%UCuts(2)%tagcuts(MatchCutNum).eq.999) then
!              res = res - ( ThePrimAmp%UCuts(2)%Coeff(MatchCutNum,0)     &
!                          + ThePrimAmp%UCuts(2)%Coeff(MatchCutNum,1)*r2     &
!                          + ThePrimAmp%UCuts(2)%Coeff(MatchCutNum,2)*r3     &
!                          + ThePrimAmp%UCuts(2)%Coeff(MatchCutNum,3)*r4     &
!                          + ThePrimAmp%UCuts(2)%Coeff(MatchCutNum,4)*r2**2        &
!                          + ThePrimAmp%UCuts(2)%Coeff(MatchCutNum,5)*r2*r3              &
!                          + ThePrimAmp%UCuts(2)%Coeff(MatchCutNum,6)*r2*r4     &
!                          + ThePrimAmp%UCuts(2)%Coeff(MatchCutNum,7)*(r3**2-r4**2)     &
!                          + ThePrimAmp%UCuts(2)%Coeff(MatchCutNum,8)*r3*r4     &
!                          + ThePrimAmp%UCuts(2)%Coeff(MatchCutNum,9)*re**2 )/prop1

              res = res - ( ThePrimAmps(PrimSelect)%UCuts(2)%Coeff(MatchCutNum,0)     &
                          + ThePrimAmps(PrimSelect)%UCuts(2)%Coeff(MatchCutNum,1)*r2     &
                          + ThePrimAmps(PrimSelect)%UCuts(2)%Coeff(MatchCutNum,2)*r3     &
                          + ThePrimAmps(PrimSelect)%UCuts(2)%Coeff(MatchCutNum,3)*r4     &
                          + ThePrimAmps(PrimSelect)%UCuts(2)%Coeff(MatchCutNum,4)*r2**2        &
                          + ThePrimAmps(PrimSelect)%UCuts(2)%Coeff(MatchCutNum,5)*r2*r3              &
                          + ThePrimAmps(PrimSelect)%UCuts(2)%Coeff(MatchCutNum,6)*r2*r4     &
                          + ThePrimAmps(PrimSelect)%UCuts(2)%Coeff(MatchCutNum,7)*(r3**2-r4**2)     &
                          + ThePrimAmps(PrimSelect)%UCuts(2)%Coeff(MatchCutNum,8)*r3*r4     &
                          + ThePrimAmps(PrimSelect)%UCuts(2)%Coeff(MatchCutNum,9)*re**2 )/prop1
            endif
        enddo

       enddo

       return
       end subroutine









!       procedure that returns polarization vectors for fermions and bosons
        SUBROUTINE givepol(PartType,q1,Dv,Ds,tag_pol,Nj,BPOL,POL)
!!      q1: momentum
!!      Nj: number of helicity states
!!      tag_pol: 0=phys. gluon pol.sum, 1=full gluon pol.sum (needed for doubcuts+single cuts)
!!      POL: first spinors
!!      BPOL: last spinors
        use ModNVBasis
        use ModParameters
        use ModMisc
        implicit none
        integer tag_pol,Ds,Dv,Nhf,Nhg,i,j,Nj
        complex(8) kuks, r12, r14
        complex(8) q1(5),p(Dv), q14(4)
        complex(8) BPOL(8,16),POL(8,16)
        complex(8) v(4,4)
        complex(8) ulist(Ds,8),barulist(Ds,8)
        complex(8) ulist_CHECK(Ds,8),barulist_CHECK(Ds,8)
        integer PartType


!       initializing polarization lists
        POL(1:8,1:16) =  (0d0,0d0)
        BPOL(1:8,1:16) = (0d0,0d0)

        if (tag_pol.eq.0) then
!!            number of hel.states: fermions, gluons
              if (Dv.eq.4) then
                Nhf=2
                Nhg=2
              elseif (Dv.eq.6) then
                Nhf=4
                Nhg=4
              elseif (Dv.eq.8) then
                Nhf=8
                Nhg=2
              endif
        elseif (tag_pol.eq.1) then
              if (Dv.eq.4) then
                Nhf=2
                Nhg=4
              elseif (Dv.eq.6) then
                Nhf=4
                Nhg=6
              elseif (Dv.eq.8) then
                Nhf=8
                Nhg=2
              endif
         endif




         if( IsAQuark(PartType) ) then

                if( PartType.lt.0 ) then
                    if( Dv.eq.4 ) then
                        p(1:4)=-q1(1:4)
                    else
                        p(1:5)=-q1(1:5)
                        p(6:Dv)= (0d0,0d0)
                    endif
                    call    give_usp(Nhf,Dv,Ds,p,GetMass(PartType),ulist)
                    call give_barusp(Nhf,Dv,Ds,p,GetMass(PartType),barulist)

!                       if(Dv.eq.4 .and. GetMass(PartType).eq.0d0) then   ! check completeness relation for j=1...4
!                         j=4
!                         q14(1:4)=(0d0,0d0)
!                         q14(j)=(1d0,0d0)
!                         call spi2(4,4,p,q14(1:4),q1(1:4))
!                         print *,ulist(1:4,1)*barulist(j,1)+ulist(1:4,2)*barulist(j,2)-q1(1:4)
!                         pause
!                       endif

                    do j=1,Ds
                    do i=1,Nhf
                        BPOL(i,j)=barulist(j,i)
                        POL(i,j)=ulist(j,i)
                    enddo
                    enddo
                else
                    if( Dv.eq.4 ) then
                        p(1:4)=+q1(1:4)
                    else
                        p(1:5)=+q1(1:5)
                        p(6:Dv)= (0d0,0d0)
                    endif
                    call    give_usp(Nhf,Dv,Ds,p,GetMass(PartType),ulist)
                    call give_barusp(Nhf,Dv,Ds,p,GetMass(PartType),barulist)
                    do j=1,Ds
                    do i=1,Nhf
                        POL(i,j)=barulist(j,i)
                        BPOL(i,j)=ulist(j,i)
                    enddo
                    enddo
                endif

                Nj=Nhf

         elseif( PartType.eq.Glu_ ) then

            if( Dv.eq.4 ) then
                  call give1to4vect_light(q1,v)
                  POL(1,1:4) =v(3,1:4)
                  POL(2,1:4) =v(4,1:4)
                  BPOL(1,1:4)=v(3,1:4)
                  BPOL(2,1:4)=v(4,1:4)

                  call sc(4,v(1,1:4),v(2,1:4),r12)
                  kuks=1d0/zsqrt(2d0*r12)

                  POL(3,1:4)=kuks*(v(1,1:4)+v(2,1:4))
                  BPOL(3,1:4)=POL(3,1:4)

                  POL(4,1:4)=dcmplx(0d0,1d0)*kuks*(v(1,1:4)-v(2,1:4))
                  BPOL(4,1:4)=POL(4,1:4)

            elseif (Dv.eq.6) then
                  q14(1:4)=q1(1:4)

                  call sc(4,q14,q14,r14)
                  r14=zsqrt(1d0*r14)

                  call give1to4vect(q14,v)
                  POL(1,1:4)=v(2,1:4)
                  BPOL(1,1:4)=v(2,1:4)
                  POL(2,1:4)=v(3,1:4)
                  BPOL(2,1:4)=v(3,1:4)
                  POL(3,1:4)=v(4,1:4)
                  BPOL(3,1:4)=v(4,1:4)
                  POL(4,1:4)=zero
                  BPOL(4,1:4)=zero
                  POL(5,1:4)=zero
                  BPOL(5,1:4)=zero
                  POL(6,1:4)=q14(1:4)/r14
                  BPOL(6,1:4)=q14(1:4)/r14
                  POL(4,5)=zero
                  BPOL(4,5)=zero
                  POL(4,6)=dcmplx(0d0,1d0)
                  BPOL(4,6)=dcmplx(0d0,1d0)
                  POL(5,5)=dcmplx(0d0,1d0)
                  BPOL(5,5)=dcmplx(0d0,1d0)
                  POL(6,5)=zero
                  BPOL(6,5)=zero
                  POL(6,6)=zero
                  BPOL(6,6)=zero
            elseif (Dv.eq.8) then
                  POL(1,1:6)=zero
                  BPOL(1,1:6)=zero
                  POL(2,1:6)=zero
                  BPOL(2,1:6)=zero
                  POL(1,7)=dcmplx(0d0,1d0)
                  BPOL(1,7)=dcmplx(0d0,1d0)
                  POL(2,7)=zero
                  BPOL(2,7)=zero
                  POL(2,8)=dcmplx(0d0,1d0)
                  BPOL(2,8)=dcmplx(0d0,1d0)
             endif
            Nj=Nhg

         elseif( IsAScalar(PartType) ) then
              Nj=1
              POL(1,1)=(1d0,0d0)
              BPOL(1,1)=(1d0,0d0)
         else
              call Error("givepol")
         endif

          END SUBROUTINE








END MODULE ModResidues_new
